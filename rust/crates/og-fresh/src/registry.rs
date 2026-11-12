//! Container registries, through the v2 API every registry speaks.

use crate::http::{Body, Fetch, Request, Want, join, origin};
use crate::{Lookup, LookupError, array, failed, need_str, opt_field, string};
use regex::Regex;
use std::sync::LazyLock;

/// An image reference, split.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Image {
    pub host: String,
    pub repo: String,
    /// `None` for a digest-only pin.
    pub tag: Option<String>,
    pub digest: Option<String>,
}

/// `quay.io/biocontainers/x:1--h_0` is (quay.io, biocontainers/x,
/// 1--h_0); `python:3.11.17@sha256:ab` is (docker.io, library/python,
/// 3.11.17, sha256:ab): a tag pinned to a digest keeps its tag; `a/b@sha256:ab`
/// has no tag; a reference with neither is `latest`.
pub fn split_image(image: &str) -> Image {
    let (name, digest) = match image.split_once('@') {
        Some((n, d)) => (n, Some(d.to_owned())),
        None => (image, None),
    };
    let last = name.rsplit('/').next().unwrap_or(name);
    let (name, tag) = if last.contains(':') {
        let (n, t) = name.rsplit_once(':').expect("a colon");
        (n, Some(t.to_owned()))
    } else {
        (
            name,
            if digest.is_some() {
                None
            } else {
                Some("latest".to_owned())
            },
        )
    };
    let parts: Vec<&str> = name.split('/').collect();
    let (host, mut repo) = if parts.len() > 1 && (parts[0].contains('.') || parts[0].contains(':'))
    {
        (parts[0].to_owned(), parts[1..].join("/"))
    } else {
        ("docker.io".to_owned(), name.to_owned())
    };
    if host == "docker.io" && !repo.contains('/') {
        repo = format!("library/{repo}");
    }
    Image {
        host,
        repo,
        tag,
        digest,
    }
}

/// The tag of an image reference (`a/b:1.2@sha256:..` is `1.2`); the
/// whole reference when it has no tag.
pub fn image_tag(reference: &str) -> &str {
    let r = reference.split('@').next().unwrap_or("");
    let last = r.rsplit('/').next().unwrap_or(r);
    if last.contains(':') {
        r.rsplit_once(':').map_or(r, |(_, t)| t)
    } else {
        r
    }
}

/// What the images section asks of a registry.
pub trait Registry: Sync {
    fn list_tags(&self, host: &str, repo: &str) -> Lookup<Vec<String>>;
    fn manifest_digest(&self, host: &str, repo: &str, reference: &str) -> Lookup<String>;
}

/// The registries on the network.
pub struct Online<'a>(pub &'a dyn Fetch);

fn registry_host(host: &str) -> &str {
    if host == "docker.io" {
        "registry-1.docker.io"
    } else {
        host
    }
}

impl Online<'_> {
    fn auth(&self, host: &str, repo: &str) -> Lookup<Option<String>> {
        let url = match host {
            "docker.io" => format!(
                "https://auth.docker.io/token?service=registry.docker.io&scope=repository:{repo}:pull"
            ),
            "ghcr.io" => format!("https://ghcr.io/token?scope=repository:{repo}:pull"),
            _ => return Ok(None),
        };
        let tok = need_str(
            &self.0.json(&url)?,
            "token",
            &format!("registry token for {host}/{repo}"),
        )?;
        Ok(Some(format!("Bearer {tok}")))
    }

    /// Every tag, page by page. The token goes only to the registry that
    /// issued it: a `next` page on another host is refused, and nothing is
    /// sent there.
    pub fn list_tags_with(
        client: &dyn Fetch,
        host: &str,
        repo: &str,
        auth: Option<String>,
    ) -> Lookup<Vec<String>> {
        static NEXT: LazyLock<Regex> =
            LazyLock::new(|| Regex::new(r#"<([^>]+)>;\s*rel="next""#).expect("link"));
        let mut url = format!("https://{}/v2/{repo}/tags/list?n=1000", registry_host(host));
        let home = {
            let (s, h) = origin(&url);
            (s.to_owned(), h.to_owned())
        };
        let mut tags = vec![];
        for _ in 0..200 {
            let mut req = Request::get(&url, Want::Json);
            if let Some(a) = &auth {
                req = req.header("Authorization", a);
            }
            let (body, headers) = client.request(&req)?;
            let Body::Json(body) = body else {
                return Err(LookupError::Shape(format!("{url}: no JSON body")));
            };
            if let Some(t) = opt_field(&body, "tags")? {
                for x in array(t, "tags")? {
                    tags.push(string(x, "a tag")?.to_owned());
                }
            }
            let Some(next) = headers.get("Link").and_then(|l| NEXT.captures(l)) else {
                break;
            };
            url = join(&url, &next[1]);
            let (s, h) = origin(&url);
            if (s, h) != (home.0.as_str(), home.1.as_str()) {
                return failed(format!(
                    "tag list of {host}/{repo} points its next page at another host: {url}"
                ));
            }
        }
        if tags.is_empty() {
            return failed(format!("tag list of {host}/{repo}: empty"));
        }
        Ok(tags)
    }
}

impl Registry for Online<'_> {
    fn list_tags(&self, host: &str, repo: &str) -> Lookup<Vec<String>> {
        let auth = self.auth(host, repo)?;
        Online::list_tags_with(self.0, host, repo, auth)
    }

    fn manifest_digest(&self, host: &str, repo: &str, reference: &str) -> Lookup<String> {
        let url = format!(
            "https://{}/v2/{repo}/manifests/{reference}",
            registry_host(host)
        );
        let mut req = Request::get(&url, Want::Head).header(
            "Accept",
            "application/vnd.oci.image.index.v1+json, \
             application/vnd.docker.distribution.manifest.list.v2+json, \
             application/vnd.oci.image.manifest.v1+json, \
             application/vnd.docker.distribution.manifest.v2+json",
        );
        if let Some(a) = self.auth(host, repo)? {
            req = req.header("Authorization", &a);
        }
        let h = self.0.head_with(req)?;
        match h.get("Docker-Content-Digest") {
            Some(d) if !d.is_empty() => Ok(d.to_owned()),
            _ => Err(LookupError::Failed(format!(
                "digest of {host}/{repo}:{reference}: empty"
            ))),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn img(h: &str, r: &str, t: Option<&str>, d: Option<&str>) -> Image {
        Image {
            host: h.into(),
            repo: r.into(),
            tag: t.map(Into::into),
            digest: d.map(Into::into),
        }
    }

    #[test]
    fn references() {
        assert_eq!(
            split_image("python:3.11"),
            img("docker.io", "library/python", Some("3.11"), None)
        );
        assert_eq!(
            split_image("quay.io/biocontainers/t1k:1.0.9--h5ca1c30_0"),
            img(
                "quay.io",
                "biocontainers/t1k",
                Some("1.0.9--h5ca1c30_0"),
                None
            )
        );
        assert_eq!(
            split_image("a/b@sha256:00"),
            img("docker.io", "a/b", None, Some("sha256:00"))
        );
        assert_eq!(
            split_image("python:3.11.17@sha256:ab"),
            img(
                "docker.io",
                "library/python",
                Some("3.11.17"),
                Some("sha256:ab")
            )
        );
        assert_eq!(
            split_image("localhost:5000/x"),
            img("localhost:5000", "x", Some("latest"), None)
        );
        assert_eq!(image_tag("pgkb/pharmcat:3.4.0@sha256:aa"), "3.4.0");
        assert_eq!(image_tag("host:5000/x"), "host:5000/x");
    }
}
