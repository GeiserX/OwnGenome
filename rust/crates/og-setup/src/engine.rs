//! setup's container runs, as scripts/lib/common.sh's run_in builds them:
//! `ENGINE run --rm` with
//!
//! - `--network none`, unless the run asks for the network (`net`);
//! - `--user UID:GID -e HOME=/tmp`, or `--user 0:0` for an image that
//!   cannot run unprivileged (`root`);
//! - the install folder read-only at `/genome`, the downloads folder
//!   read-only at `/downloads` when it lies outside it, and each `rw`
//!   folder (inside one of them) writable at its container path ([`cpath`]);
//! - SAMPLE's folder writable at `/genome/SAMPLE` when SAMPLE is set;
//!
//! then the caller's own engine options, the image and its command. The
//! engine is CONTAINER_ENGINE (docker by default); the files og writes
//! itself keep umask 077, which `og setup install` sets.

use std::path::Path;
use std::process::Command;

/// The container path of HOST inside the install folder GENOME (or the
/// downloads folder DOWNLOADS), or common.sh's message.
pub fn cpath(genome: &str, downloads: Option<&str>, host: &str) -> Result<String, String> {
    if host == genome {
        return Ok("/genome".to_owned());
    }
    if let Some(rest) = host.strip_prefix(genome).filter(|r| r.starts_with('/')) {
        return Ok(format!("/genome{rest}"));
    }
    let dl = downloads
        .filter(|d| !d.is_empty())
        .unwrap_or("/nonexistent");
    if let Some(rest) = host.strip_prefix(dl).filter(|r| r.starts_with('/')) {
        return Ok(format!("/downloads{rest}"));
    }
    Err(format!(
        "ERROR: {host} is outside GENOME_DIR ({genome}); containers only see GENOME_DIR (and OG_DOWNLOADS)."
    ))
}

/// One container run.
#[derive(Debug, Clone, Default)]
pub struct RunIn {
    /// `--net`: the run downloads something.
    pub net: bool,
    /// `--root`: the container runs as root.
    pub root: bool,
    /// `--rw DIR`: folders made writable, in order.
    pub rw: Vec<String>,
}

/// What a run is placed in: the install and downloads folders, SAMPLE, and
/// the caller's user and group.
#[derive(Debug, Clone)]
pub struct Place {
    pub genome: String,
    pub downloads: Option<String>,
    pub sample: Option<String>,
    pub uid: u32,
    pub gid: u32,
}

impl Place {
    /// From GENOME_DIR, OG_DOWNLOADS and SAMPLE of VARS, and this process's
    /// user and group.
    pub fn from_vars(v: &crate::install::Vars) -> Place {
        Place {
            genome: v.get("GENOME_DIR").unwrap_or("").to_owned(),
            downloads: v.nonempty("OG_DOWNLOADS").map(str::to_owned),
            sample: v.nonempty("SAMPLE").map(str::to_owned),
            uid: rustix::process::getuid().as_raw(),
            gid: rustix::process::getgid().as_raw(),
        }
    }
}

impl RunIn {
    /// The engine's arguments, before the caller's own (REST: engine
    /// options, the image, its command), as run_in builds them; the folders
    /// it makes writable are created.
    pub fn argv(&self, p: &Place, rest: &[String]) -> Result<Vec<String>, String> {
        let mut extra = Vec::new();
        for d in &self.rw {
            let c = cpath(&p.genome, p.downloads.as_deref(), d)?;
            std::fs::create_dir_all(d).map_err(|e| format!("ERROR: cannot make {d}: {e}"))?;
            extra.extend(["-v".to_owned(), format!("{d}:{c}")]);
        }
        let mut a: Vec<String> = vec!["run".into(), "--rm".into()];
        if !self.net {
            a.extend(["--network".into(), "none".into()]);
        }
        if self.root {
            a.extend(["--user".into(), "0:0".into()]);
        } else {
            a.extend([
                "--user".into(),
                format!("{}:{}", p.uid, p.gid),
                "-e".into(),
                "HOME=/tmp".into(),
            ]);
        }
        a.extend(["-v".into(), format!("{}:/genome:ro", p.genome)]);
        if let Some(dl) = &p.downloads
            && !dl.starts_with(&format!("{}/", p.genome))
            && Path::new(dl).is_dir()
        {
            a.extend(["-v".into(), format!("{dl}:/downloads:ro")]);
        }
        if let Some(s) = &p.sample {
            let host = format!("{}/{s}", p.genome);
            std::fs::create_dir_all(&host)
                .map_err(|e| format!("ERROR: cannot make {host}: {e}"))?;
            a.extend(["-v".into(), format!("{host}:/genome/{s}")]);
        }
        a.extend(extra);
        a.extend(rest.iter().cloned());
        Ok(a)
    }

    /// Runs it with ENGINE; the engine's exit code (1 when it could not
    /// start), its output on the terminal.
    pub fn run(&self, engine: &str, p: &Place, rest: &[String]) -> Result<i32, String> {
        use std::io::Write;
        let a = self.argv(p, rest)?;
        let _ = std::io::stdout().flush();
        Ok(
            crate::child::status(Command::new(engine).args(&a))
                .map_or(1, |s| s.code().unwrap_or(1)),
        )
    }
}
