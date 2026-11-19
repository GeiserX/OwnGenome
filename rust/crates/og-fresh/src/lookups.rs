//! One lookup per database, bundle or tool Renovate cannot see, and the
//! image and digest sections. Each lookup returns its report row and a
//! typed answer ([`Found`]); a lookup that cannot answer is an error.

use crate::http::Fetch;
use crate::pins::{Pins, hla_db_release, setting, vep_cache_template};
use crate::registry::{Registry, split_image};
use crate::version::{change_kind, max_by, newer_tags, ver};
use crate::{
    DbRow, Found, Lookup, LookupError, Upstream, array, failed, field, need_str, opt_field, string,
};
use indexmap::IndexMap;
use pycompat::json::Value;
use regex::Regex;
use std::collections::BTreeSet;
use std::path::Path;
use std::sync::LazyLock;

pub const CLINVAR_MAX_AGE_DAYS: i64 = 35;
/// pypgx 0.27.0 pulled pandas 3.0 and broke every gene (the hold note in
/// versions.env), so it never counts as a version the held pin could move to.
pub const PYPGX_BROKEN: &[&str] = &["0.27.0"];
pub const RAW: &str = "https://raw.githubusercontent.com";

/// What the lookups read besides the network.
pub struct Ctx<'a> {
    pub fetch: &'a dyn Fetch,
    pub registry: &'a dyn Registry,
    pub pins: &'a Pins,
    /// `versions.env` (the gnomAD hold note).
    pub versions: &'a Path,
    /// `scripts/setup.sh` (the ClinVar URL, Parascopy's Zenodo record).
    pub setup: &'a Path,
    /// `scripts/lib/common.sh` (the VEP cache URL, HLA_DB_RELEASE).
    pub common: &'a Path,
    /// Now, for the ClinVar file's age.
    pub now: jiff::Timestamp,
}

fn re(p: &str) -> Regex {
    Regex::new(p).expect("a fixed pattern")
}

fn row(
    item: &str,
    pinned: String,
    upstream: String,
    behind: bool,
    todo: impl FnOnce() -> String,
) -> DbRow {
    DbRow {
        item: item.into(),
        pinned,
        upstream,
        state: if behind { "behind" } else { "current" }.into(),
        todo: if behind { todo() } else { String::new() },
    }
}

fn tag_of(pins: &Pins, var: &str) -> Lookup<String> {
    split_image(pins.need(var)?).tag.ok_or_else(|| {
        LookupError::Failed(format!("{var} is pinned by digest only: no tag to compare"))
    })
}

/// `text` without its leading `v`s (Python's `lstrip("v")`).
fn unv(text: &str) -> &str {
    text.trim_start_matches('v')
}

/// A `Last-Modified` header as a time.
pub fn http_date(value: &str) -> Option<jiff::Timestamp> {
    jiff::fmt::rfc2822::DateTimeParser::new()
        .parse_timestamp(value.trim())
        .ok()
}

/// The ClinVar file `setup.sh` downloads: its `Last-Modified`.
pub fn clinvar(c: &Ctx<'_>) -> Lookup<Found> {
    static URL: LazyLock<Regex> = LazyLock::new(|| re(r#"CLINVAR_URL="(https?://[^"]+)""#));
    let url = setting(c.setup, &URL, "CLINVAR_URL")?;
    let h = c.fetch.head(&url)?;
    let lm = match h.get("Last-Modified") {
        Some(v) if !v.is_empty() => v.to_owned(),
        _ => return failed(format!("Last-Modified of {url}: empty")),
    };
    let when = http_date(&lm).ok_or_else(|| {
        LookupError::Shape(format!("cannot read the Last-Modified of {url}: {lm:?}"))
    })?;
    let age = (c.now.as_second() - when.as_second()).div_euclid(86_400);
    let current = age <= CLINVAR_MAX_AGE_DAYS;
    let date = when.to_zoned(jiff::tz::TimeZone::UTC).date();
    let name = url.rsplit('/').next().unwrap_or(&url);
    let mut todo = format!(
        "NCBI refreshes it monthly. Refresh a local copy older than {CLINVAR_MAX_AGE_DAYS} days \
         (its `##fileDate`): delete `$GENOME_DIR/clinvar/` and rerun `./scripts/setup.sh`."
    );
    if !current {
        todo += &format!(
            " Upstream itself is older than {CLINVAR_MAX_AGE_DAYS} days: check that the URL still is the live file."
        );
    }
    Ok(Found {
        row: DbRow {
            item: "ClinVar (steps 6, 23)".into(),
            pinned: "downloaded once by setup.sh".into(),
            upstream: format!(
                "{name} dated {} ({age} days ago)",
                date.strftime("%Y-%m-%d")
            ),
            state: if current { "current" } else { "behind" }.into(),
            todo,
        },
        upstream: Upstream::Date(date.strftime("%Y%m%d").to_string()),
    })
}

/// An integer of a JSON answer (Ensembl's releases are numbers).
fn int_of(v: &Value) -> Lookup<i128> {
    match v {
        Value::Int(i) => Ok(*i),
        Value::Str(s) => s
            .trim()
            .parse()
            .map_err(|_| LookupError::Shape(format!("{s:?} is not a number"))),
        _ => Err(LookupError::Shape("a release is not a number".into())),
    }
}

pub fn ensembl(c: &Ctx<'_>) -> Lookup<Found> {
    let cache = c.pins.need("VEP_CACHE_RELEASE")?;
    let image_tag = tag_of(c.pins, "VEP_IMAGE")?;
    let body = c
        .fetch
        .json("https://rest.ensembl.org/info/data/?content-type=application/json")?;
    let releases = match opt_field(&body, "releases")? {
        Some(v) if !matches!(v, Value::Array(a) if a.is_empty()) => array(v, "Ensembl releases")?,
        _ => return failed("Ensembl releases: empty"),
    };
    let current = releases.iter().map(int_of).collect::<Lookup<Vec<_>>>()?;
    let current = current.into_iter().max().expect("non-empty");
    c.fetch
        .head(&vep_cache_template(c.common)?.replace("%s", cache))?;
    let pinned: i128 = cache
        .parse()
        .map_err(|_| LookupError::Shape(format!("VEP_CACHE_RELEASE {cache:?} is not a number")))?;
    Ok(Found {
        row: row(
            "Ensembl / VEP cache (step 13)",
            format!("VEP_CACHE_RELEASE={cache}, VEP_IMAGE {image_tag}"),
            format!("Ensembl release {current}"),
            current > pinned,
            || {
                format!(
                    "Move VEP_IMAGE to release_{current}.x and VEP_CACHE_RELEASE to {current} in one change \
                     (a major: Renovate waits for approval); rerun steps 13, 30, 23 and 31."
                )
            },
        ),
        upstream: Upstream::Release(current.to_string()),
    })
}

pub fn pcgr(c: &Ctx<'_>) -> Lookup<Found> {
    let tag = tag_of(c.pins, "PCGR_IMAGE")?;
    let (bundle, vep) = (
        c.pins.need("PCGR_DATA_BUNDLE")?,
        c.pins.need("PCGR_VEP_CACHE_RELEASE")?,
    );
    let latest = need_str(
        &c.fetch.github("/repos/sigven/pcgr/releases/latest")?,
        "tag_name",
        "PCGR latest release",
    )?;
    let version = unv(&latest);
    let src = c
        .fetch
        .text(&format!("{RAW}/sigven/pcgr/{latest}/pcgr/pcgr_vars.py"))?;
    static DB: LazyLock<Regex> = LazyLock::new(|| re(r"(?m)^DB_VERSION\s*=\s*'(\d+)'"));
    static VV: LazyLock<Regex> = LazyLock::new(|| re(r"(?m)^VEP_VERSION\s*=\s*'(\d+)'"));
    let (Some(db), Some(vv)) = (DB.captures(&src), VV.captures(&src)) else {
        return failed(format!(
            "no DB_VERSION or VEP_VERSION in pcgr_vars.py at {latest}"
        ));
    };
    let (db, vv) = (&db[1], &vv[1]);
    let mut on = vec![];
    for (name, host) in [
        ("Docker Hub sigven/pcgr", "docker.io"),
        ("ghcr.io/sigven/pcgr", "ghcr.io"),
    ] {
        if c.registry
            .list_tags(host, "sigven/pcgr")?
            .iter()
            .any(|t| t == version)
        {
            on.push(name);
        }
    }
    let behind = ver(version) > ver(&tag) || db != bundle || vv != vep;
    let where_ = if on.is_empty() {
        "no registry yet".to_owned()
    } else {
        on.join(", ")
    };
    Ok(Found {
        row: row(
            "PCGR / CPSR (step 17)",
            format!("{tag}, bundle {bundle}, VEP {vep}"),
            format!("{version}, bundle {db}, VEP {vv}; image on {where_}"),
            behind,
            || {
                format!(
                    "PCGR_IMAGE, PCGR_DATA_BUNDLE and PCGR_VEP_CACHE_RELEASE move together; rerun step 17.{}",
                    if on.is_empty() {
                        " Wait for a published image."
                    } else {
                        ""
                    }
                )
            },
        ),
        upstream: Upstream::Release(db.to_owned()),
    })
}

/// Verify a pinned cache URL and its entry in Ensembl's small CHECKSUMS file.
/// Large archives are only read with HEAD, never downloaded by freshness.
fn cache_download(c: &Ctx<'_>, pin: &str, label: &str) -> Lookup<Found> {
    let release = c.pins.need(pin)?;
    let url = vep_cache_template(c.common)?.replace("%s", release);
    c.fetch.head(&url)?;
    let (base, name) = url.rsplit_once('/').expect("cache URL has a basename");
    let checksums = c.fetch.text(&format!("{base}/CHECKSUMS"))?;
    let listed = checksums.lines().any(|line| {
        let fields: Vec<_> = line.split_whitespace().collect();
        fields.len() == 3
            && fields[2] == name
            && fields[0].parse::<u16>().is_ok()
            && fields[1].parse::<u64>().is_ok_and(|n| n > 0)
    });
    if !listed {
        return failed(format!("{name} has no usable entry in {base}/CHECKSUMS"));
    }
    Ok(Found {
        row: row(
            label,
            release.into(),
            "URL live; CHECKSUMS entry present".into(),
            false,
            String::new,
        ),
        upstream: Upstream::ReportOnly,
    })
}

/// The release used by steps 13, 30, 23 and 31.
pub fn vep_download(c: &Ctx<'_>) -> Lookup<Found> {
    cache_download(c, "VEP_CACHE_RELEASE", "VEP cache download")
}

/// The separate release used inside CPSR's image.
pub fn cpsr_download(c: &Ctx<'_>) -> Lookup<Found> {
    cache_download(c, "PCGR_VEP_CACHE_RELEASE", "CPSR cache download")
}

/// The pinned PCGR bundle URL and size. SHA-256 is checked by setup, which
/// downloads the file; freshness deliberately reads only its headers.
pub fn pcgr_download(c: &Ctx<'_>) -> Lookup<Found> {
    let bundle = c.pins.need("PCGR_DATA_BUNDLE")?;
    let expected = c
        .pins
        .need("PCGR_DATA_BYTES")?
        .parse::<u64>()
        .map_err(|_| LookupError::Shape("PCGR_DATA_BYTES is not a byte count".into()))?;
    let url = format!("https://insilico.hpc.uio.no/pcgr/pcgr_ref_data.{bundle}.grch38.tgz");
    let headers = c.fetch.head(&url)?;
    let size = headers
        .get("Content-Length")
        .and_then(|s| s.parse::<u64>().ok());
    if expected == 0 || size != Some(expected) {
        return failed(format!(
            "{url}: Content-Length {size:?}, expected {expected}"
        ));
    }
    Ok(Found {
        row: row(
            "PCGR bundle download",
            bundle.into(),
            format!("URL live; {expected} bytes"),
            false,
            String::new,
        ),
        upstream: Upstream::ReportOnly,
    })
}

pub fn pypgx(c: &Ctx<'_>) -> Lookup<Found> {
    let image_tag = tag_of(c.pins, "PYPGX_IMAGE")?;
    let held = image_tag.split("--").next().unwrap_or("").to_owned();
    let bundle = c.pins.need("PYPGX_BUNDLE_VERSION")?;
    let tags = c
        .fetch
        .github("/repos/sbslee/pypgx-bundle/tags?per_page=100")?;
    let tags = array(&tags, "pypgx-bundle tags")?;
    if tags.is_empty() {
        return failed("pypgx-bundle tags: empty");
    }
    static TAG: LazyLock<Regex> = LazyLock::new(|| re(r"^v?\d+\.\d+\.\d+$"));
    let mut bundles = BTreeSet::new();
    for t in tags {
        let name = string(field(t, "name")?, "a tag name")?;
        if TAG.is_match(name) {
            bundles.insert(unv(name).to_owned());
        }
    }
    if !bundles.contains(bundle) {
        return failed(format!(
            "PYPGX_BUNDLE_VERSION {bundle} is not a pypgx-bundle tag"
        ));
    }
    let pypi = c.fetch.json("https://pypi.org/pypi/pypgx/json")?;
    let pypi = match opt_field(&pypi, "info")? {
        Some(info) => need_str(info, "version", "pypgx on PyPI")?,
        None => return failed("pypgx on PyPI: empty"),
    };
    let images: BTreeSet<String> = c
        .registry
        .list_tags("quay.io", "biocontainers/pypgx")?
        .iter()
        .filter_map(|t| crate::version::bioconda_version(t).map(str::to_owned))
        .collect();
    if images.is_empty() {
        return failed("pypgx biocontainer versions: empty");
    }
    let newest_bundle = max_by(bundles.iter(), |v| ver(v)).expect("non-empty");
    let newest_image = max_by(images.iter(), |v| ver(v)).expect("non-empty");
    let mut movable: Vec<&String> = images
        .intersection(&bundles)
        .filter(|v| ver(v) > ver(&held) && !PYPGX_BROKEN.contains(&v.as_str()))
        .collect();
    movable.sort_by_key(|v| ver(v));
    let to = movable.last().map(|v| (*v).clone());
    Ok(Found {
        row: DbRow {
            item: "pypgx and pypgx-bundle (step 32)".into(),
            pinned: format!("pypgx {held} (held), bundle {bundle}"),
            upstream: format!(
                "PyPI {pypi}, newest biocontainer {newest_image}, newest bundle tag {newest_bundle}"
            ),
            state: to.as_ref().map_or_else(
                || "held, cannot move yet".to_owned(),
                |v| format!("can move to {v}"),
            ),
            todo: match &to {
                Some(v) => format!(
                    "pypgx {v} has an image and a bundle tag: move PYPGX_IMAGE and PYPGX_BUNDLE_VERSION together and rerun step 32."
                ),
                None => format!(
                    "The hold waits for a pypgx release past {} with a biocontainer and a pypgx-bundle tag of the same version.",
                    PYPGX_BROKEN.join(", ")
                ),
            },
        },
        upstream: Upstream::CanMoveTo(to),
    })
}

pub fn pharmcat(c: &Ctx<'_>) -> Lookup<Found> {
    let tag = tag_of(c.pins, "PHARMCAT_IMAGE")?;
    let latest = need_str(
        &c.fetch.github("/repos/PharmGKB/PharmCAT/releases/latest")?,
        "tag_name",
        "PharmCAT latest release",
    )?;
    let latest = unv(&latest).to_owned();
    let has_image = c
        .registry
        .list_tags("docker.io", "pgkb/pharmcat")?
        .contains(&latest);
    Ok(Found {
        row: row(
            "PharmCAT (steps 7, 27)",
            tag.clone(),
            format!(
                "{latest}{}",
                if has_image {
                    ""
                } else {
                    " (no pgkb/pharmcat image yet)"
                }
            ),
            ver(&latest) > ver(&tag),
            || "Revalidate steps 7 and 27 on a known sample before moving PHARMCAT_IMAGE.".into(),
        ),
        upstream: Upstream::Release(latest),
    })
}

/// (version, date) from the `release_version.txt` of an ANHIG
/// repository's `Latest` branch.
pub fn release_version(
    f: &dyn Fetch,
    repo: &str,
    product: &str,
) -> Lookup<(String, Option<String>)> {
    let text = f.text(&format!("{RAW}/ANHIG/{repo}/Latest/release_version.txt"))?;
    let v = re(&format!(
        r"#\s*version:\s*{}\s+(\d+(?:\.\d+)+)",
        regex::escape(product)
    ));
    static D: LazyLock<Regex> = LazyLock::new(|| re(r"#\s*date:\s*(\S+)"));
    let Some(v) = v.captures(&text) else {
        return failed(format!(
            "no '{product}' version line in {repo} release_version.txt"
        ));
    };
    Ok((v[1].to_owned(), D.captures(&text).map(|d| d[1].to_owned())))
}

fn dated(v: &str, d: &Option<String>) -> String {
    match d {
        Some(d) => format!("{v} ({d})"),
        None => v.to_owned(),
    }
}

pub fn imgt(c: &Ctx<'_>) -> Lookup<Found> {
    let pinned = hla_db_release(c.common)?;
    let (latest, date) = release_version(c.fetch, "IMGTHLA", "IPD-IMGT/HLA")?;
    Ok(Found {
        row: row(
            "IPD-IMGT/HLA (step 8)",
            format!("HLA_DB_RELEASE={pinned} (scripts/lib/common.sh)"),
            dated(&latest, &date),
            ver(&latest) > ver(&pinned),
            || {
                "Renovate opens this bump (ANHIG/IMGTHLA). After it merges, rerun `setup.sh` and step 8, \
                 which builds a new T1K index for the release."
                    .into()
            },
        ),
        upstream: Upstream::Release(latest),
    })
}

pub fn kir(c: &Ctx<'_>) -> Lookup<Found> {
    let pinned = c.pins.need("KIR_DB_RELEASE")?;
    let (latest, date) = release_version(c.fetch, "IPDKIR", "IPD-KIR")?;
    Ok(Found {
        row: row(
            "IPD-KIR (step 8 with KIR=true)",
            pinned.into(),
            dated(&latest, &date),
            ver(&latest) > ver(pinned),
            || {
                format!(
                    "Move KIR_DB_RELEASE to {latest} (`setup.sh --kir-data` reads branch {} of ANHIG/IPDKIR), \
                     then rerun step 8 with KIR=true.",
                    latest.replace('.', "")
                )
            },
        ),
        upstream: Upstream::Release(latest),
    })
}

pub fn pgsc_calc(c: &Ctx<'_>) -> Lookup<Found> {
    let pinned = c.pins.need("PGSC_CALC_VERSION")?;
    let latest = need_str(
        &c.fetch
            .github("/repos/PGScatalog/pgsc_calc/releases/latest")?,
        "tag_name",
        "pgsc_calc latest release",
    )?;
    Ok(Found {
        row: row(
            "pgsc_calc (steps 25, 26)",
            pinned.into(),
            latest.clone(),
            ver(unv(&latest)) > ver(unv(pinned)),
            || {
                "Renovate lists it on the dependency dashboard. Read its conf/modules.config and nextflow.config \
                 and move the PGSC_*_IMAGE lines, PGSC_CALC_SHA256 and PGSC_CALC_NF_SCHEMA with it (docs/25-prs.md)."
                    .into()
            },
        ),
        upstream: Upstream::Release(latest),
    })
}

pub fn dbsnp(c: &Ctx<'_>) -> Lookup<Found> {
    let pinned = c.pins.need("CONSUMER_DBSNP_BUILD")?;
    let text = c.fetch.text("https://ftp.ncbi.nih.gov/snp/archive/")?;
    static B: LazyLock<Regex> = LazyLock::new(|| re(r#"href="b(\d+)/""#));
    let builds: BTreeSet<u64> = B
        .captures_iter(&text)
        .filter_map(|m| m[1].parse().ok())
        .collect();
    let Some(&newest) = builds.last() else {
        return failed("dbSNP builds in the archive listing: empty");
    };
    let Some(have) = pinned
        .strip_prefix('b')
        .and_then(|n| n.parse::<u64>().ok())
        .filter(|_| pinned.len() > 1 && pinned[1..].bytes().all(|b| b.is_ascii_digit()))
    else {
        return failed(format!("CONSUMER_DBSNP_BUILD {pinned:?} is not b<number>"));
    };
    Ok(Found {
        row: row(
            "dbSNP (step 40 site table)",
            pinned.into(),
            format!("b{newest}"),
            newest > have,
            || {
                format!(
                    "Move CONSUMER_DBSNP_BUILD to b{newest} with its md5, rebuild the table once \
                     (`setup.sh --consumer-sites`, about 30 GB) and record CONSUMER_SITES_TABLE_SHA256 (docs/updating.md)."
                )
            },
        ),
        upstream: Upstream::Release(format!("b{newest}")),
    })
}

/// A JSON value as Python's `str()` prints it in a report cell.
fn cell(v: Option<&Value>) -> String {
    match v {
        None | Some(Value::Null) => "None".into(),
        Some(Value::Str(s)) => s.clone(),
        Some(Value::Int(i)) => i.to_string(),
        Some(other) => pycompat::value::py_str(other).unwrap_or_default(),
    }
}

pub fn parascopy_data(c: &Ctx<'_>) -> Lookup<Found> {
    let pinned = c.pins.need("PARASCOPY_DATA_VERSION")?;
    static REC: LazyLock<Regex> = LazyLock::new(|| re(r"(?m)^PARASCOPY_DATA_RECORD=(\d+)"));
    let record = setting(c.setup, &REC, "PARASCOPY_DATA_RECORD")?;
    let body = c.fetch.json(&format!(
        "https://zenodo.org/api/records/{record}/versions?size=25&sort=version"
    ))?;
    let what = format!("versions of Zenodo record {record}");
    let hits = match opt_field(&body, "hits")? {
        Some(h) => opt_field(h, "hits")?,
        None => None,
    };
    let hits = match hits {
        Some(Value::Array(a)) if !a.is_empty() => a,
        Some(Value::Array(_)) | None => return failed(format!("{what}: empty")),
        Some(_) => return Err(LookupError::Shape(format!("{what} is not a list"))),
    };
    static FILE: LazyLock<Regex> = LazyLock::new(|| re(r"^GRCh38_v(\d+(?:\.\d+)*)\.tar\.gz$"));
    let mut found: indexmap::IndexMap<String, String> = indexmap::IndexMap::new();
    for h in hits {
        let id = cell(opt_field(h, "id")?);
        if let Some(files) = opt_field(h, "files")? {
            for f in array(files, "files")? {
                let key = match opt_field(f, "key")? {
                    Some(k) => string(k, "a file key")?,
                    None => "",
                };
                if let Some(m) = FILE.captures(key) {
                    found.insert(m[1].to_owned(), id.clone());
                }
            }
        }
    }
    let Some(newest) = max_by(found.keys(), |v| ver(v)).cloned() else {
        return failed(format!("a GRCh38_v<version>.tar.gz in the {what}: empty"));
    };
    Ok(Found {
        row: row(
            "Parascopy data (step 35)",
            format!("{pinned} (Zenodo record {record})"),
            format!("{newest} (Zenodo record {})", found[&newest]),
            ver(&newest) > ver(pinned),
            || {
                "Check that PARASCOPY_IMAGE reads it, then move PARASCOPY_DATA_VERSION, PARASCOPY_DATA_RECORD \
                 and the two md5s in scripts/setup.sh together and rerun step 35."
                    .into()
            },
        ),
        upstream: Upstream::Release(newest),
    })
}

pub fn gnomad(c: &Ctx<'_>) -> Lookup<Found> {
    static HOLD: LazyLock<Regex> = LazyLock::new(|| re(r"gnomad_v(\d+(?:\.\d+)*)_constraint"));
    let text = std::fs::read_to_string(c.versions).unwrap_or_default();
    let Some(pinned) = HOLD.captures(&text).map(|m| m[1].to_owned()) else {
        return failed("no gnomad_v<version>_constraint hold note in versions.env");
    };
    let api = "https://storage.googleapis.com/storage/v1/b/gcp-public-data--gnomad/o";
    let body = c
        .fetch
        .json(&format!("{api}?prefix=release/&delimiter=/"))?;
    let prefixes = match opt_field(&body, "prefixes")? {
        Some(Value::Array(a)) if !a.is_empty() => a,
        Some(Value::Array(_)) | None => return failed("gnomAD releases: empty"),
        Some(_) => return Err(LookupError::Shape("gnomAD releases is not a list".into())),
    };
    static REL: LazyLock<Regex> = LazyLock::new(|| re(r"^\d+(\.\d+)*$"));
    let mut found = vec![];
    for p in prefixes {
        let rel = string(p, "a prefix")?
            .split('/')
            .nth(1)
            .ok_or_else(|| LookupError::Shape("a prefix with no release".into()))?
            .to_owned();
        if REL.is_match(&rel) && ver(&rel) > ver(&pinned) {
            let items = c.fetch.json(&format!(
                "{api}?prefix=release/{rel}/constraint/&maxResults=1"
            ))?;
            if opt_field(&items, "items")?.is_some_and(pycompat::value::truthy) {
                found.push(rel);
            }
        }
    }
    let newest = max_by(found.iter(), |v| ver(v)).cloned();
    Ok(Found {
        row: DbRow {
            item: "gnomAD constraint (steps 23, 31)".into(),
            pinned: format!("v{pinned} (held)"),
            upstream: match &newest {
                Some(n) => format!("newest release with constraint data: {n}"),
                None => format!("no release past v{pinned}"),
            },
            state: "report only".into(),
            todo: if newest.is_some() {
                "Held: the 4.1.2 table has a different layout; read it and adapt steps 23 and 31 before moving.".into()
            } else {
                String::new()
            },
        },
        upstream: Upstream::ReportOnly,
    })
}

/// The first two parts of a dotted version (`25.10.4` is `25.10`).
fn line_of(v: &str) -> String {
    v.split('.').take(2).collect::<Vec<_>>().join(".")
}

pub fn nextflow(c: &Ctx<'_>) -> Lookup<Found> {
    let pinned = c.pins.need("NEXTFLOW_VERSION")?;
    let rels = c
        .fetch
        .github("/repos/nextflow-io/nextflow/releases?per_page=100")?;
    let rels = array(&rels, "Nextflow releases")?;
    if rels.is_empty() {
        return failed("Nextflow releases: empty");
    }
    static STABLE: LazyLock<Regex> = LazyLock::new(|| re(r"^v?\d+\.\d+\.\d+$"));
    let mut stable = vec![];
    for r in rels {
        let tag = string(field(r, "tag_name")?, "a release tag")?;
        let pre = opt_field(r, "prerelease")?.is_some_and(pycompat::value::truthy);
        if !pre && STABLE.is_match(tag) {
            stable.push(unv(tag).to_owned());
        }
    }
    if stable.is_empty() {
        return failed("stable Nextflow releases: empty");
    }
    let line = line_of(pinned);
    let in_line = max_by(
        stable.iter().filter(|s| s.starts_with(&format!("{line}."))),
        |v| ver(v),
    )
    .cloned();
    let mut newer: indexmap::IndexMap<String, String> = indexmap::IndexMap::new();
    for s in &stable {
        let ln = line_of(s);
        let best = newer.get(&ln).map_or_else(|| ver("0"), |b| ver(b));
        if ver(&ln) > ver(&line) && ver(s) > best {
            newer.insert(ln, s.clone());
        }
    }
    let newest_line = max_by(newer.values(), |v| ver(v)).cloned();
    let behind = in_line.as_ref().is_some_and(|v| ver(v) > ver(pinned));
    let mut todo = String::new();
    if behind {
        todo += "Move NEXTFLOW_VERSION within the line. ";
    }
    if let Some(n) = &newest_line {
        todo += &format!("Moving to {n} is a new line: run the stub and e2e jobs on it first.");
    }
    Ok(Found {
        row: DbRow {
            item: "Nextflow".into(),
            pinned: pinned.into(),
            upstream: format!(
                "{line} line: {}; newest stable line: {}",
                in_line.as_deref().unwrap_or("none"),
                newest_line.as_deref().unwrap_or(&line)
            ),
            state: if behind {
                "behind"
            } else if newest_line.is_some() {
                "newer line"
            } else {
                "current"
            }
            .into(),
            todo,
        },
        upstream: Upstream::Release(in_line.unwrap_or_else(|| pinned.to_owned())),
    })
}

pub fn cyrius(c: &Ctx<'_>) -> Lookup<Found> {
    let pinned = c.pins.need("CYRIUS_VERSION")?;
    let body = c.fetch.json("https://pypi.org/pypi/cyrius/json")?;
    let latest = match opt_field(&body, "info")? {
        Some(info) => need_str(info, "version", "Cyrius on PyPI")?,
        None => return failed("Cyrius on PyPI: empty"),
    };
    Ok(Found {
        row: row(
            "Cyrius (step 21, legacy)",
            pinned.into(),
            format!("PyPI {latest}"),
            ver(&latest) > ver(pinned),
            || {
                "Check scripts/cyrius-constraints.txt still resolves before moving CYRIUS_VERSION."
                    .into()
            },
        ),
        upstream: Upstream::Release(latest),
    })
}

// The curated PGx tables. Each row of assets/pgx_contraindications.tsv,
// assets/pgx_inhibitors.tsv and assets/pgx_genes.tsv was copied by hand from
// a source that can change under it: a CPIC level moves, a label is
// withdrawn, the FDA reclassifies an inhibitor, pypgx gains a copy-number
// model. These lookups compare every row with its source now. A row that no
// longer agrees is `drift` (a warning in the issue, not an error); a source
// that does not answer, or answers empty, is an error as everywhere here.

/// The CPIC API's pairs of one gene (the `Pair_source` of a row).
pub const CPIC_PAIRS: &str = "https://api.cpicpgx.org/v1/pair_view?genesymbol=eq.";
/// One DailyMed label's version history, by its set id.
pub const DAILYMED_HISTORY: &str = "https://dailymed.nlm.nih.gov/dailymed/services/v2/spls/";
/// The FDA's table of clinical CYP inhibitors (the `Class_source` of the
/// inhibitor rows).
pub const FDA_CYP_TABLE: &str = "https://www.fda.gov/drugs/drug-interactions-labeling/healthcare-professionals-fdas-examples-drugs-interact-cyp-enzymes-and-transporter-systems";

/// The rows of a curated table under `assets/` beside `versions.env`, by
/// column name. `#` lines are comments; the first other line is the header.
pub fn curated(c: &Ctx<'_>, name: &str, need: &[&str]) -> Lookup<Vec<IndexMap<String, String>>> {
    let path = c
        .versions
        .parent()
        .unwrap_or_else(|| Path::new(""))
        .join("assets")
        .join(name);
    let text = std::fs::read_to_string(&path)
        .map_err(|e| LookupError::Failed(format!("cannot read {}: {e}", path.display())))?;
    let mut lines = text
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty());
    let Some(header) = lines.next() else {
        return failed(format!("{name}: no header"));
    };
    let header: Vec<&str> = header.split('\t').collect();
    if let Some(miss) = need.iter().find(|n| !header.contains(n)) {
        return Err(LookupError::Shape(format!("{name} has no {miss} column")));
    }
    let rows: Vec<IndexMap<String, String>> = lines
        .map(|l| {
            header
                .iter()
                .zip(l.split('\t'))
                .map(|(h, v)| ((*h).to_owned(), v.trim().to_owned()))
                .collect()
        })
        .collect();
    if rows.is_empty() {
        return failed(format!("{name}: no rows"));
    }
    Ok(rows)
}

fn col<'a>(row: &'a IndexMap<String, String>, name: &str) -> &'a str {
    row.get(name).map_or("", String::as_str)
}

/// The report row of a table check: `current` when nothing drifted, else
/// `drift` with each drifted row named.
fn table_row(item: &str, pinned: String, agree: String, drift: &[String], todo: &str) -> Found {
    let upstream = if drift.is_empty() {
        agree
    } else {
        format!("{} difference(s): {}", drift.len(), drift.join("; "))
    };
    Found {
        row: DbRow {
            item: item.into(),
            pinned,
            upstream,
            state: if drift.is_empty() { "current" } else { "drift" }.into(),
            todo: if drift.is_empty() {
                String::new()
            } else {
                todo.into()
            },
        },
        upstream: Upstream::ReportOnly,
    }
}

/// `assets/pgx_contraindications.tsv` against the CPIC API (the pair, its
/// level, still provisional) and DailyMed (the label the note quotes is
/// still published).
pub fn pgx_contraindications(c: &Ctx<'_>) -> Lookup<Found> {
    let rows = curated(
        c,
        "pgx_contraindications.tsv",
        &["Gene", "Drug", "CPIC_level", "Label_source"],
    )?;
    let mut drift = vec![];
    let mut pairs: IndexMap<String, Vec<Value>> = IndexMap::new();
    let mut labels: IndexMap<String, bool> = IndexMap::new();
    static SETID: LazyLock<Regex> = LazyLock::new(|| re(r"setid=([0-9A-Fa-f-]{36})\b"));
    for row in &rows {
        let (gene, drug, level) = (col(row, "Gene"), col(row, "Drug"), col(row, "CPIC_level"));
        if !pairs.contains_key(gene) {
            let body = c.fetch.json(&format!("{CPIC_PAIRS}{gene}"))?;
            let list = array(&body, &format!("CPIC pairs of {gene}"))?;
            // An empty answer cannot tell a dropped gene from an outage.
            if list.is_empty() {
                return failed(format!("CPIC pairs of {gene}: empty"));
            }
            pairs.insert(gene.to_owned(), list.to_vec());
        }
        let found = pairs[gene].iter().find(|p| {
            matches!(opt_field(p, "drugname"), Ok(Some(Value::Str(d))) if d.eq_ignore_ascii_case(drug))
        });
        let Some(pair) = found else {
            drift.push(format!("{gene} {drug}: CPIC no longer lists the pair"));
            continue;
        };
        let up = need_str(pair, "cpiclevel", &format!("CPIC level of {gene} {drug}"))?;
        if up != level {
            drift.push(format!(
                "{gene} {drug}: CPIC level {up}, the table says {level}"
            ));
        }
        if !matches!(opt_field(pair, "provisional")?, Some(Value::Bool(true))) {
            drift.push(format!(
                "{gene} {drug}: no longer provisional in CPIC (a guideline: PharmCAT may call the gene)"
            ));
        }
        let Some(setid) = SETID
            .captures(col(row, "Label_source"))
            .map(|m| m[1].to_owned())
        else {
            return Err(LookupError::Shape(format!(
                "{gene} {drug}: Label_source has no DailyMed setid"
            )));
        };
        let live = match labels.get(&setid) {
            Some(l) => *l,
            None => {
                let body = c
                    .fetch
                    .json(&format!("{DAILYMED_HISTORY}{setid}/history.json"))?;
                let history = match opt_field(&body, "data")? {
                    Some(d) => opt_field(d, "history")?,
                    None => return failed(format!("DailyMed label {setid}: empty")),
                };
                let live = matches!(history, Some(Value::Array(h)) if !h.is_empty());
                labels.insert(setid.clone(), live);
                live
            }
        };
        if !live {
            drift.push(format!(
                "{gene} {drug}: DailyMed no longer publishes label {setid}"
            ));
        }
    }
    Ok(table_row(
        "Drug warnings for disease genes (assets/pgx_contraindications.tsv, step 27)",
        format!("{} rows, CPIC levels and FDA labels", rows.len()),
        format!(
            "CPIC lists every pair at the same level, provisional; {} label(s) published",
            labels.len()
        ),
        &drift,
        "Read the CPIC pair and the label again, then fix, move or drop the row in \
         assets/pgx_contraindications.tsv (its date and Note too) and rerun step 27.",
    ))
}

/// `name` of the FDA table as a medication list names it: lower case, no
/// footnote marks, no markup.
fn fda_cell(raw: &str) -> String {
    static SUP: LazyLock<Regex> = LazyLock::new(|| re(r"(?s)<sup>.*?</sup>"));
    static TAG: LazyLock<Regex> = LazyLock::new(|| re(r"(?s)<[^>]*>"));
    let t = SUP.replace_all(raw, "");
    let t = TAG.replace_all(&t, " ");
    let t = t.replace("&nbsp;", " ").replace("&amp;", "&");
    t.split_whitespace()
        .collect::<Vec<_>>()
        .join(" ")
        .to_lowercase()
}

/// The CYP enzymes of an FDA table cell (`2C19; 2D6 strong inhibitor`).
fn fda_genes(cell: &str) -> Vec<String> {
    static G: LazyLock<Regex> = LazyLock::new(|| re(r"(?i)\b(1A2|2B6|2C8|2C9|2C19|2D6|3A)\b"));
    G.captures_iter(cell).map(|m| m[1].to_uppercase()).collect()
}

/// The FDA table's enzyme name for a table gene: CYP3A4 and CYP3A5 are its
/// `3A`.
fn fda_gene(gene: &str) -> String {
    let g = gene.trim_start_matches("CYP");
    if g.starts_with("3A") {
        "3A".into()
    } else {
        g.to_owned()
    }
}

/// One substance of the FDA table with its strong and moderate CYP
/// inhibitions as `(enzyme, strength)`.
pub type FdaRow = (String, Vec<(String, &'static str)>);

/// The FDA table: (its date, each substance with its strong and moderate
/// CYP inhibitions as `(enzyme, strength)`).
pub fn fda_inhibitors(html: &str) -> Lookup<(Option<String>, Vec<FdaRow>)> {
    static TABLE: LazyLock<Regex> = LazyLock::new(|| re(r#"(?s)summary="Table 1\."(.*?)</table>"#));
    static TD: LazyLock<Regex> = LazyLock::new(|| re(r"(?s)<td[^>]*>(.*?)</td>"));
    static DATE: LazyLock<Regex> =
        LazyLock::new(|| re(r#"(?s)Content current as of.{0,300}?datetime="(\d{4}-\d{2}-\d{2})"#));
    let Some(table) = TABLE.captures(html) else {
        return Err(LookupError::Shape(
            "the FDA page has no Table 1 of substrates, inhibitors and inducers".into(),
        ));
    };
    let mut out = vec![];
    for tr in table[1].split("<tr") {
        let cells: Vec<String> = TD.captures_iter(tr).map(|m| fda_cell(&m[1])).collect();
        if cells.len() < 3 || cells[0].is_empty() {
            continue;
        }
        let mut inh = vec![];
        for (cell, strength) in [(&cells[1], "strong"), (&cells[2], "moderate")] {
            inh.extend(fda_genes(cell).into_iter().map(|g| (g, strength)));
        }
        out.push((cells[0].clone(), inh));
    }
    if out.is_empty() {
        return failed("the FDA table of CYP inhibitors: empty");
    }
    Ok((DATE.captures(html).map(|m| m[1].to_owned()), out))
}

/// The parts of a substance of the FDA table: `lopinavir and ritonavir` is
/// lopinavir and ritonavir.
fn fda_parts(name: &str) -> Vec<&str> {
    name.split(" and ")
        .map(|p| p.trim_matches(|ch: char| !ch.is_alphanumeric()))
        .filter(|p| !p.is_empty())
        .collect()
}

/// `assets/pgx_inhibitors.tsv` against the FDA's table of CYP inhibitors:
/// every row's class, and every strong or moderate inhibitor of a gene the
/// table covers (strong only for CYP3A) that has no row.
pub fn pgx_inhibitors(c: &Ctx<'_>) -> Lookup<Found> {
    let rows = curated(c, "pgx_inhibitors.tsv", &["Drug", "Gene", "Strength"])?;
    let html = c.fetch.text(FDA_CYP_TABLE)?;
    let (date, fda) = fda_inhibitors(&html)?;
    let strength_of = |names: &[&FdaRow], g: &str| {
        let s: Vec<&str> = names
            .iter()
            .flat_map(|(_, inh)| inh.iter().filter(|(e, _)| e == g).map(|(_, s)| *s))
            .collect();
        if s.contains(&"strong") {
            Some("strong")
        } else if s.contains(&"moderate") {
            Some("moderate")
        } else {
            None
        }
    };
    let mut drift = vec![];
    for row in &rows {
        let (drug, gene, ours) = (col(row, "Drug"), col(row, "Gene"), col(row, "Strength"));
        // Its own rows and the ritonavir-boosted combinations it is part
        // of: a partner row takes the combination's class (indinavir alone
        // is no strong inhibitor; indinavir and ritonavir is).
        let named: Vec<_> = fda
            .iter()
            .filter(|(n, _)| fda_parts(n).contains(&drug))
            .collect();
        if named.is_empty() {
            drift.push(format!("{drug} {gene}: not in the FDA table"));
            continue;
        }
        match strength_of(&named, &fda_gene(gene)) {
            Some(s) if s == ours => {}
            Some(s) => drift.push(format!(
                "{drug} {gene}: the FDA says {s}, the table says {ours}"
            )),
            None => drift.push(format!(
                "{drug} {gene}: the FDA no longer lists it as a strong or moderate inhibitor"
            )),
        }
    }
    // What the table covers: strong and moderate inhibitors of CYP2D6,
    // CYP2C19, CYP2C9 and CYP2B6, strong inhibitors of CYP3A.
    let has = |drug: &str, g: &str| {
        rows.iter()
            .any(|r| col(r, "Drug") == drug && fda_gene(col(r, "Gene")) == g)
    };
    for (name, inh) in &fda {
        let first = fda_parts(name).first().copied().unwrap_or(name.as_str());
        for (g, s) in inh {
            let covered = matches!(g.as_str(), "2D6" | "2C19" | "2C9" | "2B6")
                || (g == "3A" && *s == "strong");
            if covered && !has(name, g) && !has(first, g) {
                drift.push(format!(
                    "{name}: the FDA lists a {s} CYP{g} inhibitor the table has no row for"
                ));
            }
        }
    }
    Ok(table_row(
        "Drug interactions (assets/pgx_inhibitors.tsv, step 27)",
        format!("{} rows, FDA classes read 2026-11-03", rows.len()),
        format!(
            "the FDA table{} classes every row the same and has no inhibitor the table lacks",
            date.map(|d| format!(" (current as of {d})"))
                .unwrap_or_default()
        ),
        &drift,
        "Read the FDA table again, then fix, add or drop the rows in assets/pgx_inhibitors.tsv \
         (and the date in its header) and rerun step 27 with a medication list.",
    ))
}

/// The fields of one CSV line; a quoted field may hold commas.
fn csv_fields(line: &str) -> Vec<String> {
    let (mut out, mut cur, mut quoted) = (vec![], String::new(), false);
    for ch in line.chars() {
        match ch {
            '"' => quoted = !quoted,
            ',' if !quoted => out.push(std::mem::take(&mut cur)),
            _ => cur.push(ch),
        }
    }
    out.push(cur);
    out
}

/// `assets/pgx_genes.tsv` against the gene table of pypgx's newest release:
/// every row is a pypgx gene, and the bam rows are exactly the genes with a
/// copy-number model (SV=TRUE). tests/test_pypgx_genes.sh checks the same
/// at the pinned version; this says when a newer pypgx moves it.
pub fn pgx_genes(c: &Ctx<'_>) -> Lookup<Found> {
    let rows = curated(c, "pgx_genes.tsv", &["gene", "mode"])?;
    let latest = need_str(
        &c.fetch.github("/repos/sbslee/pypgx/releases/latest")?,
        "tag_name",
        "pypgx latest release",
    )?;
    let csv = c.fetch.text(&format!(
        "{RAW}/sbslee/pypgx/{latest}/pypgx/api/data/gene-table.csv"
    ))?;
    let mut lines = csv.lines().filter(|l| !l.trim().is_empty());
    let header = csv_fields(lines.next().unwrap_or(""));
    let (Some(gi), Some(si)) = (
        header.iter().position(|h| h == "Gene"),
        header.iter().position(|h| h == "SV"),
    ) else {
        return Err(LookupError::Shape(format!(
            "pypgx {latest} gene-table.csv has no Gene or SV column"
        )));
    };
    let mut sv: indexmap::IndexMap<String, bool> = indexmap::IndexMap::new();
    for l in lines {
        let f = csv_fields(l);
        if let (Some(g), Some(s)) = (f.get(gi), f.get(si)) {
            sv.insert(g.clone(), s == "TRUE");
        }
    }
    if sv.is_empty() {
        return failed(format!("pypgx {latest} gene-table.csv: empty"));
    }
    let mut drift = vec![];
    for row in &rows {
        let (gene, mode) = (col(row, "gene"), col(row, "mode"));
        match sv.get(gene) {
            None => drift.push(format!("{gene}: not in pypgx {latest}'s gene table")),
            Some(false) if mode == "bam" => drift.push(format!(
                "{gene}: a bam row, but pypgx {latest} has no copy-number model for it"
            )),
            _ => {}
        }
    }
    for (gene, has_sv) in &sv {
        if *has_sv
            && !rows
                .iter()
                .any(|r| col(r, "gene") == gene && col(r, "mode") == "bam")
        {
            drift.push(format!(
                "{gene}: pypgx {latest} has a copy-number model, but it is not a bam row"
            ));
        }
    }
    let bam = rows.iter().filter(|r| col(r, "mode") == "bam").count();
    Ok(table_row(
        "pypgx genes (assets/pgx_genes.tsv, step 32)",
        format!(
            "{} genes, {bam} bam (pypgx {})",
            rows.len(),
            c.pins.get("PYPGX_BUNDLE_VERSION").unwrap_or("?")
        ),
        format!("pypgx {latest}'s gene table agrees"),
        &drift,
        "When PYPGX_IMAGE moves to this pypgx, change the rows of assets/pgx_genes.tsv \
         with it (tests/test_pypgx_genes.sh checks them at the pin) and rerun step 32.",
    ))
}

/// A lookup by its report name.
pub type LookupFn = fn(&Ctx<'_>) -> Lookup<Found>;

/// The database lookups in report order, by the name an error row carries.
pub const DATABASES: &[(&str, LookupFn)] = &[
    ("Ensembl / VEP cache", ensembl),
    ("PCGR / CPSR", pcgr),
    ("VEP cache download", vep_download),
    ("CPSR cache download", cpsr_download),
    ("PCGR bundle download", pcgr_download),
    ("pypgx", pypgx),
    ("PharmCAT", pharmcat),
    ("IPD-IMGT/HLA", imgt),
    ("IPD-KIR", kir),
    ("Nextflow", nextflow),
    ("pgsc_calc", pgsc_calc),
    ("dbSNP", dbsnp),
    ("Cyrius", cyrius),
    ("Parascopy data", parascopy_data),
    ("gnomAD constraint", gnomad),
    ("PGx contraindications table", pgx_contraindications),
    ("PGx inhibitors table", pgx_inhibitors),
    ("PGx genes table", pgx_genes),
];

/// One image of the images section.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ImageState {
    /// Nothing newer of the pin's shape.
    Current,
    /// Only a new build of the pinned version.
    Rebuild,
    /// A `N.N` tag that follows its own patch releases; the newest line.
    Line(String),
    /// Behind: (newest, kind), then the newest on the pin's major when it
    /// differs.
    Behind(Vec<(String, &'static str)>),
}

/// The images section's answer for one `*_IMAGE` pin.
pub fn image(reg: &dyn Registry, reference: &str) -> Option<Lookup<(String, ImageState)>> {
    let img = split_image(reference);
    let tag = img.tag.clone()?; // digest only: the digests section reports it
    Some((|| {
        let tags = reg.list_tags(&img.host, &img.repo)?;
        if !tags.contains(&tag) {
            return failed(format!(
                "{reference}: the pinned tag is not in the registry's tag list"
            ));
        }
        let (newest, same) = newer_tags(&tag, &tags)?;
        let Some(newest) = newest else {
            return Ok((tag, ImageState::Current));
        };
        // python:3.11 follows 3.11.x by itself; 3.14 is a new line, a
        // choice rather than a missed update, so it is reported apart.
        static NN: LazyLock<Regex> = LazyLock::new(|| re(r"^\d+\.\d+$"));
        if NN.is_match(&tag) {
            let patch = re(&format!(r"^{}\.\d+$", regex::escape(&tag)));
            if tags.iter().any(|t| patch.is_match(t)) {
                return Ok((tag, ImageState::Line(newest)));
            }
        }
        if change_kind(&tag, &newest) == "build" {
            return Ok((tag, ImageState::Rebuild));
        }
        let mut rows = vec![];
        if let Some(s) = same.filter(|s| *s != newest) {
            rows.push((s.clone(), change_kind(&tag, &s)));
        }
        let kind = change_kind(&tag, &newest);
        rows.push((newest, kind));
        Ok((tag, ImageState::Behind(rows)))
    })())
}

/// The digests section's answer for a pin with a digest: (pinned, what the
/// tag points at now, state).
pub fn digest(reg: &dyn Registry, reference: &str) -> Option<Lookup<(String, String, String)>> {
    let img = split_image(reference);
    let pinned_digest = img.digest.clone()?;
    Some((|| {
        let now =
            reg.manifest_digest(&img.host, &img.repo, img.tag.as_deref().unwrap_or("latest"))?;
        let versioned: Vec<String> = if img.tag.is_some() {
            vec![]
        } else {
            reg.list_tags(&img.host, &img.repo)?
                .into_iter()
                .filter(|t| t.bytes().any(|b| b.is_ascii_digit()))
                .collect()
        };
        let mut state = match &img.tag {
            Some(tag) if now == pinned_digest => {
                format!("`{tag}` still points at the pinned digest")
            }
            Some(tag) => {
                format!("`{tag}` was rebuilt: move the digest after the image test passes on it")
            }
            None if now == pinned_digest => "same as publisher's latest".into(),
            None => "publisher's latest moved".into(),
        };
        if !versioned.is_empty() {
            let mut v = versioned;
            v.sort_by_key(|t| ver(t));
            state += &format!(
                "; versioned tags now exist: {}",
                v[v.len().saturating_sub(3)..].join(", ")
            );
        }
        let short = |d: &str| d.chars().take(19).collect::<String>();
        let pinned = match &img.tag {
            Some(t) => format!("{}:{t}@{}", img.repo, short(&pinned_digest)),
            None => format!("{}@{}", img.repo, short(&pinned_digest)),
        };
        let at = format!(
            "{}: {}",
            img.tag.as_deref().unwrap_or("latest"),
            short(&now)
        );
        Ok((pinned, at, state))
    })())
}

/// The ClinVar file's date from `upstream`, when the lookup gave one.
pub fn date_of(found: &Found) -> Option<&str> {
    match &found.upstream {
        Upstream::Date(d) => Some(d),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn http_dates() {
        let t = http_date("Tue, 06 Oct 2026 12:00:00 GMT").unwrap();
        assert_eq!(t.as_second(), 1_791_288_000);
        assert!(http_date("not a date").is_none());
        assert_eq!(line_of("25.10.4"), "25.10");
    }
}
