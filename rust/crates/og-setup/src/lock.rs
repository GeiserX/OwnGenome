//! The installers' lock, as scripts/lib/common.sh's lock_acquire makes it:
//! a folder holding its owner's process id in `pid`.
//!
//! A lock is stale when that process is gone (a killed run), or when it has
//! no pid a minute after it was made (a run killed between making the
//! folder and writing the pid). A stale lock is renamed before it is
//! removed, so of two waiters that judged it stale only one takes it over.
//! A live one is waited for, five seconds at a time. The lock goes when its
//! [`Lock`] is dropped.

use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::{Duration, SystemTime};

/// A held lock; dropping it removes the folder (lock_release).
#[derive(Debug)]
pub struct Lock {
    path: PathBuf,
}

impl Drop for Lock {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.path);
    }
}

/// Whether a process with this pid text is alive, as `kill -0 PID` says:
/// a pid that is not a number, a process that is gone, or one this user may
/// not signal all count as gone.
pub fn alive(pid: &str) -> bool {
    match pid.parse::<i32>() {
        Ok(0) => true,
        Ok(n) if n > 0 => rustix::process::Pid::from_raw(n)
            .is_some_and(|p| rustix::process::test_kill_process(p).is_ok()),
        _ => false,
    }
}

/// Whether a lock folder with no pid is older than a minute
/// (`find LOCK -maxdepth 0 -mmin +1`).
pub fn older_than_a_minute(lock: &Path, now: SystemTime) -> bool {
    std::fs::metadata(lock)
        .and_then(|m| m.modified())
        .ok()
        .and_then(|t| now.duration_since(t).ok())
        .is_some_and(|age| age > Duration::from_secs(60))
}

/// Why a lock is stale, or `None` when its owner may still be working.
pub fn stale(lock: &Path, now: SystemTime) -> Option<String> {
    let pid = std::fs::read_to_string(lock.join("pid"))
        .map(|p| p.trim_end_matches('\n').to_owned())
        .unwrap_or_default();
    if !pid.is_empty() {
        return (!alive(&pid)).then(|| format!("process {pid} is gone"));
    }
    older_than_a_minute(lock, now).then(|| "no owner after a minute".to_owned())
}

/// lock_acquire LOCK: wait for LOCK, take it over when it is stale, and hold
/// it with this process's id. Messages go to OUT, as the shell's echo did.
pub fn acquire(lock: &Path, out: &mut dyn Write) -> std::io::Result<Lock> {
    loop {
        if std::fs::create_dir(lock).is_ok() {
            break;
        }
        if let Some(why) = stale(lock, SystemTime::now()) {
            let aside = PathBuf::from(format!("{}.stale.{}", lock.display(), std::process::id()));
            if std::fs::rename(lock, &aside).is_ok() {
                let _ = writeln!(out, "  Removing stale lock {} ({why}).", lock.display());
                let _ = std::fs::remove_dir_all(&aside);
            }
            continue;
        }
        // A lock that vanished between the two reads is taken next round.
        if !lock.exists() {
            continue;
        }
        let pid = std::fs::read_to_string(lock.join("pid"))
            .map(|p| p.trim_end_matches('\n').to_owned())
            .unwrap_or_default();
        let held = if pid.is_empty() {
            "starting"
        } else {
            pid.as_str()
        };
        let _ = writeln!(
            out,
            "  Waiting for {} (held by process {held})...",
            lock.display()
        );
        let ran = std::process::Command::new("sleep")
            .arg("5")
            .stdin(std::process::Stdio::null())
            .status()
            .is_ok_and(|s| s.success());
        if !ran {
            std::thread::sleep(Duration::from_secs(5));
        }
    }
    let held = Lock {
        path: lock.to_path_buf(),
    };
    std::fs::write(lock.join("pid"), format!("{}\n", std::process::id()))?;
    Ok(held)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tmp(name: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!("og-lock-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    #[test]
    fn a_dead_pid_is_stale_and_taken_over() {
        let d = tmp("dead");
        let lock = d.join("x.lock");
        std::fs::create_dir(&lock).unwrap();
        // A pid far above any kernel's default pid_max.
        std::fs::write(lock.join("pid"), "2147483646\n").unwrap();
        let mut out = Vec::new();
        let held = acquire(&lock, &mut out).unwrap();
        let text = String::from_utf8(out).unwrap();
        assert_eq!(
            text,
            format!(
                "  Removing stale lock {} (process 2147483646 is gone).\n",
                lock.display()
            )
        );
        assert_eq!(
            std::fs::read_to_string(lock.join("pid")).unwrap(),
            format!("{}\n", std::process::id())
        );
        drop(held);
        assert!(!lock.exists(), "the lock stays after it is released");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn a_live_pid_is_not_stale() {
        let d = tmp("live");
        let lock = d.join("x.lock");
        std::fs::create_dir(&lock).unwrap();
        std::fs::write(lock.join("pid"), format!("{}\n", std::process::id())).unwrap();
        assert_eq!(stale(&lock, SystemTime::now()), None);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn a_pidless_lock_is_stale_only_after_a_minute() {
        let d = tmp("pidless");
        let lock = d.join("x.lock");
        std::fs::create_dir(&lock).unwrap();
        let now = SystemTime::now();
        assert_eq!(
            stale(&lock, now),
            None,
            "a new pid-less lock is a run starting"
        );
        assert_eq!(stale(&lock, now + Duration::from_secs(59)), None);
        assert_eq!(
            stale(&lock, now + Duration::from_secs(62)).as_deref(),
            Some("no owner after a minute")
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn a_pid_that_is_no_number_is_gone() {
        assert!(!alive("abc"));
        assert!(!alive("-3"));
        assert!(alive(&std::process::id().to_string()));
    }
}
