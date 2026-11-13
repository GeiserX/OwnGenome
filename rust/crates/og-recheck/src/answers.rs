//! Where the upstream answers come from: the network (through `og-fresh`'s
//! typed lookups), an offline JSON file, and a memo so a check taken again
//! after `--rerun` asks nobody a second time.

use crate::{date8, value_text};
use indexmap::IndexMap;
use og_fresh::http::{Client, Fetch, Request, Want};
use og_fresh::lookups::{self, Ctx, http_date};
use og_fresh::pins::{Pins, hla_release_in, plain_pins};
use og_fresh::registry::Online as Registries;
use og_fresh::{Lookup, Upstream, failed, field, need_str, string};
use pycompat::json::{Value, loads};
use std::cell::RefCell;
use std::collections::HashMap;
use std::path::{Path, PathBuf};

pub const PGS_API: &str = "https://www.pgscatalog.org/rest/score";

/// The newest pipeline release: its tag and its pins (`hla_database` from
/// its `scripts/lib/common.sh`).
pub type Release = (String, IndexMap<String, String>);

/// The answers one check needs.
pub trait Answers {
    /// The newest ClinVar file's date, `YYYYMMDD`.
    fn clinvar(&self) -> Lookup<Option<String>>;
    /// The newest upstream release of a looked-up pin (`VEP_CACHE_RELEASE`,
    /// `PCGR_DATA_BUNDLE`, `PYPGX_BUNDLE_VERSION`, `PHARMCAT_IMAGE`,
    /// `hla_database`).
    fn pinned(&self, key: &str) -> Lookup<Option<String>>;
    /// The date of a PGS score's harmonised file upstream, `YYYYMMDD`.
    fn pgs(&self, id: &str) -> Lookup<Option<String>>;
    fn release(&self) -> Lookup<Release>;
    /// The hosts asked, in order.
    fn hosts(&self) -> Vec<String>;
    fn offline(&self) -> bool;
}

/// Answers from a JSON file: `{item: {"upstream": value} | {"error":
/// text}, "release": {"tag", "pins"}}`.
pub struct Offline(IndexMap<String, Value>);

impl Offline {
    pub fn read(path: &Path) -> Result<Offline, String> {
        let text = std::fs::read_to_string(path)
            .map_err(|e| format!("cannot read {}: {e}", path.display()))?;
        match loads(&text) {
            Ok(Value::Object(m)) => Ok(Offline(m)),
            Ok(_) => Err(format!(
                "{}: the offline answers are not an object",
                path.display()
            )),
            Err(e) => Err(format!("{}: {e}", path.display())),
        }
    }

    fn entry(&self, key: &str) -> Lookup<&IndexMap<String, Value>> {
        let Some(Value::Object(e)) = self.0.get(key) else {
            return failed(format!("{key} is not in the offline answers"));
        };
        match e.get("error") {
            Some(v) if pycompat::value::truthy(v) => failed(value_text(v).unwrap_or_default()),
            _ => Ok(e),
        }
    }

    fn upstream(&self, key: &str) -> Lookup<Option<String>> {
        Ok(self.entry(key)?.get("upstream").and_then(value_text))
    }
}

impl Answers for Offline {
    fn clinvar(&self) -> Lookup<Option<String>> {
        Ok(self.upstream("clinvar_file_date")?.and_then(|u| date8(&u)))
    }
    fn pinned(&self, key: &str) -> Lookup<Option<String>> {
        self.upstream(key)
    }
    fn pgs(&self, id: &str) -> Lookup<Option<String>> {
        Ok(self.upstream(&format!("pgs:{id}"))?.and_then(|u| date8(&u)))
    }
    fn release(&self) -> Lookup<Release> {
        let e = self.entry("release")?;
        let tag = match e.get("tag") {
            Some(Value::Str(t)) if !t.is_empty() => t.clone(),
            _ => return failed("the offline answers' release needs a tag and a pins object"),
        };
        let pins = match e.get("pins") {
            None | Some(Value::Null) => IndexMap::new(),
            Some(Value::Object(p)) => p
                .iter()
                .filter_map(|(k, v)| value_text(v).map(|t| (k.clone(), t)))
                .collect(),
            Some(_) => return failed("the offline answers' release needs a tag and a pins object"),
        };
        Ok((tag, pins))
    }
    fn hosts(&self) -> Vec<String> {
        vec![]
    }
    fn offline(&self) -> bool {
        true
    }
}

/// Answers from the network: `og-fresh`'s lookups for ClinVar and the pins
/// (on the run's own pins), one PGS Catalog lookup per score, and the
/// pipeline's newest release.
pub struct Online {
    pub client: Client,
    /// The run's manifest pins (`version` rows).
    pub pins: Pins,
    pub repo: String,
    pub pgs_api: String,
    pub raw: String,
    pub root: PathBuf,
}

impl Online {
    fn ctx<'a>(&'a self, reg: &'a Registries<'a>, paths: &'a [PathBuf; 3]) -> Ctx<'a> {
        Ctx {
            fetch: &self.client,
            registry: reg,
            pins: &self.pins,
            versions: &paths[0],
            setup: &paths[1],
            common: &paths[2],
            now: jiff::Timestamp::now(),
        }
    }

    fn found(&self, f: lookups::LookupFn) -> Lookup<Upstream> {
        let reg = Registries(&self.client);
        let paths = [
            self.root.join("versions.env"),
            self.root.join("scripts/setup.sh"),
            self.root.join("scripts/lib/common.sh"),
        ];
        Ok(f(&self.ctx(&reg, &paths))?.upstream)
    }
}

impl Answers for Online {
    fn clinvar(&self) -> Lookup<Option<String>> {
        match self.found(lookups::clinvar)? {
            Upstream::Date(d) => Ok(Some(d)),
            other => failed(format!("no date in the ClinVar answer: {other:?}")),
        }
    }

    fn pinned(&self, key: &str) -> Lookup<Option<String>> {
        let f: lookups::LookupFn = match key {
            "VEP_CACHE_RELEASE" => lookups::ensembl,
            "PCGR_DATA_BUNDLE" => lookups::pcgr,
            "PYPGX_BUNDLE_VERSION" => lookups::pypgx,
            "PHARMCAT_IMAGE" => lookups::pharmcat,
            "hla_database" => lookups::imgt,
            _ => return failed(format!("no upstream lookup for {key}")),
        };
        Ok(match self.found(f)? {
            Upstream::Release(r) => Some(r),
            // No pypgx past the held one has both an image and a bundle
            // tag: upstream offers nothing newer than the run's pin.
            Upstream::CanMoveTo(to) => {
                to.or_else(|| self.pins.get("PYPGX_BUNDLE_VERSION").map(str::to_owned))
            }
            Upstream::Date(d) => Some(d),
            Upstream::ReportOnly => None,
        })
    }

    fn pgs(&self, id: &str) -> Lookup<Option<String>> {
        let body = self.client.json(&format!("{}/{id}", self.pgs_api))?;
        let same = matches!(&body, Value::Object(m) if m.get("id") == Some(&Value::Str(id.into())));
        if !same {
            return failed(format!("the PGS Catalog has no score {id}"));
        }
        let files = field(&body, "ftp_harmonized_scoring_files")?;
        let url = string(
            field(field(files, "GRCh38")?, "positions")?,
            "the harmonised file URL",
        )?;
        let h = self.client.head_with(Request::get(url, Want::Head))?;
        let lm = match h.get("Last-Modified") {
            Some(v) if !v.is_empty() => v,
            _ => return failed(format!("Last-Modified of {url}: empty")),
        };
        match http_date(lm) {
            Some(t) => Ok(Some(
                t.to_zoned(jiff::tz::TimeZone::UTC)
                    .date()
                    .strftime("%Y%m%d")
                    .to_string(),
            )),
            None => failed(format!("cannot read the Last-Modified of {url}: {lm:?}")),
        }
    }

    fn release(&self) -> Lookup<Release> {
        let tag = need_str(
            &self
                .client
                .github(&format!("/repos/{}/releases/latest", self.repo))?,
            "tag_name",
            &format!("latest release of {}", self.repo),
        )?;
        let text = self
            .client
            .text(&format!("{}/{}/{tag}/versions.env", self.raw, self.repo))?;
        let common = self.client.text(&format!(
            "{}/{}/{tag}/scripts/lib/common.sh",
            self.raw, self.repo
        ))?;
        let mut pins = plain_pins(&text);
        if let Some(h) = hla_release_in(&common) {
            pins.insert("hla_database".into(), h);
        }
        Ok((tag, pins))
    }

    fn hosts(&self) -> Vec<String> {
        self.client.hosts()
    }

    fn offline(&self) -> bool {
        false
    }
}

/// The answers of one check, kept: the check taken again after `--rerun`
/// asks nobody a second time.
pub struct Memo<A> {
    inner: A,
    seen: RefCell<HashMap<String, Lookup<Option<String>>>>,
    release: RefCell<Option<Lookup<Release>>>,
}

impl<A: Answers> Memo<A> {
    pub fn new(inner: A) -> Self {
        Memo {
            inner,
            seen: RefCell::new(HashMap::new()),
            release: RefCell::new(None),
        }
    }

    fn ask(
        &self,
        key: String,
        f: impl FnOnce() -> Lookup<Option<String>>,
    ) -> Lookup<Option<String>> {
        if let Some(v) = self.seen.borrow().get(&key) {
            return v.clone();
        }
        let v = f();
        self.seen.borrow_mut().insert(key, v.clone());
        v
    }
}

impl<A: Answers> Answers for Memo<A> {
    fn clinvar(&self) -> Lookup<Option<String>> {
        self.ask("clinvar".into(), || self.inner.clinvar())
    }
    fn pinned(&self, key: &str) -> Lookup<Option<String>> {
        self.ask(format!("pinned:{key}"), || self.inner.pinned(key))
    }
    fn pgs(&self, id: &str) -> Lookup<Option<String>> {
        self.ask(format!("pgs:{id}"), || self.inner.pgs(id))
    }
    fn release(&self) -> Lookup<Release> {
        if let Some(r) = self.release.borrow().as_ref() {
            return r.clone();
        }
        let r = self.inner.release();
        *self.release.borrow_mut() = Some(r.clone());
        r
    }
    fn hosts(&self) -> Vec<String> {
        self.inner.hosts()
    }
    fn offline(&self) -> bool {
        self.inner.offline()
    }
}

/// Online answers for `og recheck` on this machine.
pub fn online(pins: Pins, repo: &str, pgs_api: &str, raw: &str, root: &Path) -> Online {
    let retries = std::env::var("OWNGENOME_RECHECK_RETRIES")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(og_fresh::http::RETRIES);
    let mut client = Client::new(retries);
    client.announce = Some("recheck");
    Online {
        client,
        pins,
        repo: repo.into(),
        pgs_api: pgs_api.into(),
        raw: raw.into(),
        root: root.to_path_buf(),
    }
}

impl Answers for Box<dyn Answers> {
    fn clinvar(&self) -> Lookup<Option<String>> {
        (**self).clinvar()
    }
    fn pinned(&self, key: &str) -> Lookup<Option<String>> {
        (**self).pinned(key)
    }
    fn pgs(&self, id: &str) -> Lookup<Option<String>> {
        (**self).pgs(id)
    }
    fn release(&self) -> Lookup<Release> {
        (**self).release()
    }
    fn hosts(&self) -> Vec<String> {
        (**self).hosts()
    }
    fn offline(&self) -> bool {
        (**self).offline()
    }
}
