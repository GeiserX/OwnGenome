//! The installers' lock, as scripts/lib/common.sh's lock_acquire makes it:
//! a folder holding its owner's process id in `pid`.
//!
//! A lock is stale when that process is gone (a killed run), or when it has
//! no pid a minute after it was made (a run killed between making the
//! folder and writing the pid). A stale lock is renamed before it is
//! removed, so of two waiters that judged it stale only one takes it over.
//! A live one is waited for, five seconds at a time. The lock goes when its
//! [`Lock`] is dropped, and when SIGINT, SIGTERM or SIGHUP ends og (Ctrl-C,
//! a kill, a closed terminal): the shell's `trap 'lock_release ...' EXIT`
//! removed it then too. A signal ignored when og started (nohup's SIGHUP) stays
//! ignored, as bash keeps it. A lock left by a run killed with SIGKILL is
//! taken over as stale.

use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::{Mutex, MutexGuard, OnceLock};
use std::time::{Duration, SystemTime};

/// A held lock; dropping it removes the folder (lock_release).
#[derive(Debug)]
pub struct Lock {
    path: PathBuf,
}

impl Drop for Lock {
    fn drop(&mut self) {
        let mut h = held();
        let _ = std::fs::remove_dir_all(&self.path);
        h.retain(|p| p != &self.path);
    }
}

/// The locks this process holds, for [`on_signals`].
static HELD: Mutex<Vec<PathBuf>> = Mutex::new(Vec::new());

fn held() -> MutexGuard<'static, Vec<PathBuf>> {
    HELD.lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner)
}

/// Once per process, before the first lock is made: on SIGINT, SIGTERM or
/// SIGHUP, remove every held lock and exit 128 + the signal (130, 143,
/// 129), the status a shell gives a script those signals end. A signal
/// this process ignored when it started stays ignored, as bash keeps it:
/// nohup's SIGHUP, the SIGINT of a job a script started with `&`.
fn on_signals() -> std::io::Result<()> {
    use signal_hook::consts::{SIGHUP, SIGINT, SIGTERM};
    static ARMED: OnceLock<Result<(), String>> = OnceLock::new();
    ARMED
        .get_or_init(|| {
            let ignored = ignored_at_start();
            let watch: Vec<i32> = [SIGINT, SIGTERM, SIGHUP]
                .into_iter()
                .filter(|s| !ignored.contains(s))
                .collect();
            if watch.is_empty() {
                return Ok(());
            }
            let mut signals =
                signal_hook::iterator::Signals::new(watch).map_err(|e| e.to_string())?;
            std::thread::Builder::new()
                .name("og-lock-signals".to_owned())
                .spawn(move || {
                    if let Some(sig) = signals.forever().next() {
                        // Held to the end, so no lock is made or dropped
                        // after the list is emptied.
                        let mut h = held();
                        for p in h.drain(..) {
                            let _ = std::fs::remove_dir_all(&p);
                        }
                        std::process::exit(128 + sig);
                    }
                })
                .map(drop)
                .map_err(|e| e.to_string())
        })
        .clone()
        .map_err(|e| std::io::Error::other(format!("cannot watch for signals: {e}")))
}

/// The signals this process ignores (SIG_IGN), read before og sets any
/// handler, so they are the ones it started with. Linux lists them in
/// /proc/self/status (`SigIgn:`, a hex mask, bit N-1 for signal N).
/// Elsewhere a child bash inherits them and `trap -p` prints each as
/// `trap -- '' SIGHUP` (bash 4.2 or later, and setup needs 4.4; the 3.2 of
/// macOS's /bin/bash, or no bash, lists none).
fn ignored_at_start() -> Vec<i32> {
    if let Ok(status) = std::fs::read_to_string("/proc/self/status") {
        return status
            .lines()
            .find_map(|l| l.strip_prefix("SigIgn:"))
            .and_then(|m| u64::from_str_radix(m.trim(), 16).ok())
            .map(|mask| (1..=64).filter(|n| mask >> (n - 1) & 1 == 1).collect())
            .unwrap_or_default();
    }
    let Ok(o) = std::process::Command::new("bash")
        .args(["-c", "trap -p"])
        .stdin(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .output()
    else {
        return Vec::new();
    };
    String::from_utf8_lossy(&o.stdout)
        .lines()
        .filter_map(|l| l.strip_prefix("trap -- '' "))
        .filter_map(|name| {
            [
                signal_hook::consts::SIGINT,
                signal_hook::consts::SIGTERM,
                signal_hook::consts::SIGHUP,
            ]
            .into_iter()
            .find(|&s| signal_hook::low_level::signal_name(s) == Some(name.trim()))
        })
        .collect()
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
///
/// Only a lock that is there already is waited for. Any other reason the
/// folder cannot be made (its parent is not writable, the disk is
/// read-only) is an error at once, `cannot create the lock LOCK: ...`.
/// lock_acquire waits forever in that case, printing that the lock is held.
pub fn acquire(lock: &Path, out: &mut dyn Write) -> std::io::Result<Lock> {
    on_signals()?;
    loop {
        // Made and listed under the list's mutex, so a signal removes
        // every lock this process made and none it did not.
        let mut h = held();
        match std::fs::create_dir(lock) {
            Ok(()) => {
                h.push(lock.to_path_buf());
                break;
            }
            Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => drop(h),
            Err(e) => {
                return Err(std::io::Error::new(
                    e.kind(),
                    format!("cannot create the lock {}: {e}", lock.display()),
                ));
            }
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
        let by = if pid.is_empty() {
            "starting"
        } else {
            pid.as_str()
        };
        let _ = writeln!(
            out,
            "  Waiting for {} (held by process {by})...",
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
    let guard = Lock {
        path: lock.to_path_buf(),
    };
    std::fs::write(lock.join("pid"), format!("{}\n", std::process::id())).map_err(|e| {
        std::io::Error::new(
            e.kind(),
            format!("cannot write {}/pid: {e}", lock.display()),
        )
    })?;
    Ok(guard)
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

    /// A folder this user may not write into is no contention: acquire
    /// says so at once instead of trying again (it spun with no message).
    #[test]
    fn a_lock_that_cannot_be_made_is_an_error() {
        use std::os::unix::fs::PermissionsExt;
        if rustix::process::geteuid().is_root() {
            return; // root writes into a read-only folder
        }
        let d = tmp("readonly");
        let ro = d.join("ro");
        std::fs::create_dir(&ro).unwrap();
        std::fs::set_permissions(&ro, std::fs::Permissions::from_mode(0o555)).unwrap();
        let lock = ro.join("x.lock");
        let (tx, rx) = std::sync::mpsc::channel();
        let l = lock.clone();
        std::thread::spawn(move || {
            let _ = tx.send(acquire(&l, &mut std::io::sink()).map(drop));
        });
        let got = rx.recv_timeout(Duration::from_secs(20));
        std::fs::set_permissions(&ro, std::fs::Permissions::from_mode(0o755)).unwrap();
        let _ = std::fs::remove_dir_all(&d);
        let e = got
            .expect("acquire did not return for a folder it may not write into")
            .expect_err("acquire took a lock it could not make");
        assert_eq!(e.kind(), std::io::ErrorKind::PermissionDenied, "{e}");
        assert_eq!(
            e.to_string(),
            format!(
                "cannot create the lock {}: Permission denied (os error 13)",
                lock.display()
            )
        );
    }

    /// Ctrl-C, a kill or a closed terminal removes the lock and exits
    /// 128 + the signal, as setup.sh's `trap 'lock_release ...' EXIT` does.
    /// This test runs itself as the child that holds the lock.
    #[test]
    fn a_signal_removes_the_held_lock() {
        use std::io::BufRead;
        use std::process::{Command, Stdio};
        const CHILD: &str = "OG_LOCK_TEST_HOLD";
        if let Some(lock) = std::env::var_os(CHILD) {
            let _held = acquire(Path::new(&lock), &mut std::io::sink()).unwrap();
            println!("lock-held");
            let _ = std::io::stdout().flush();
            std::thread::sleep(Duration::from_secs(60));
            std::process::exit(9);
        }
        let d = tmp("signal");
        // The child, started by bash after IGNORE (a `trap '' SIG` that
        // ignores signals, as nohup or a script's `&` does); it says when
        // it holds LOCK.
        let hold = |lock: &Path, ignore: &str| {
            let mut child = Command::new("bash")
                .args(["-c", &format!("{ignore}exec \"$0\" \"$@\"")])
                .arg(std::env::current_exe().unwrap())
                .args([
                    "--exact",
                    "lock::tests::a_signal_removes_the_held_lock",
                    "--nocapture",
                    "--test-threads=1",
                ])
                .env(CHILD, lock)
                .stdin(Stdio::null())
                .stdout(Stdio::piped())
                .stderr(Stdio::null())
                .spawn()
                .unwrap();
            let mut lines = std::io::BufReader::new(child.stdout.take().unwrap()).lines();
            assert!(
                lines.any(|l| l.is_ok_and(|l| l.contains("lock-held"))),
                "the child never held {}",
                lock.display()
            );
            assert!(lock.join("pid").exists());
            let pid = rustix::process::Pid::from_raw(child.id() as i32).unwrap();
            (child, pid)
        };
        for (sig, code) in [
            (rustix::process::Signal::INT, 130),
            (rustix::process::Signal::TERM, 143),
            (rustix::process::Signal::HUP, 129),
        ] {
            // A signal this test itself started with ignored (under nohup)
            // is ignored in the child too, which is the case below.
            if ignored_at_start().contains(&sig.as_raw()) {
                continue;
            }
            let lock = d.join(format!("{code}.lock"));
            let (mut child, pid) = hold(&lock, "");
            rustix::process::kill_process(pid, sig).unwrap();
            let status = child.wait().unwrap();
            assert_eq!(status.code(), Some(code), "after {sig:?}: {status}");
            assert!(!lock.exists(), "{sig:?} left {}", lock.display());
        }
        // Ignored when og started: still ignored, and the lock stays held.
        // Without /proc, only a bash that lists inherited ignores (4.2 or
        // later; macOS's /bin/bash 3.2 does not) tells og about them.
        let told = Path::new("/proc/self/status").exists()
            || Command::new("bash")
                .args(["-c", "trap '' HUP; bash -c 'trap -p'"])
                .output()
                .is_ok_and(|o| String::from_utf8_lossy(&o.stdout).contains("SIGHUP"));
        if !told {
            let _ = std::fs::remove_dir_all(&d);
            return;
        }
        let lock = d.join("nohup.lock");
        let (mut child, pid) = hold(&lock, "trap '' HUP INT; ");
        for sig in [rustix::process::Signal::HUP, rustix::process::Signal::INT] {
            rustix::process::kill_process(pid, sig).unwrap();
        }
        std::thread::sleep(Duration::from_millis(500));
        let gone = child.try_wait().unwrap();
        assert!(
            gone.is_none(),
            "an ignored signal ended the child: {gone:?}"
        );
        assert!(
            lock.join("pid").exists(),
            "an ignored signal removed the lock"
        );
        rustix::process::kill_process(pid, rustix::process::Signal::TERM).unwrap();
        assert_eq!(child.wait().unwrap().code(), Some(143));
        assert!(!lock.exists());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn a_pid_that_is_no_number_is_gone() {
        assert!(!alive("abc"));
        assert!(!alive("-3"));
        assert!(alive(&std::process::id().to_string()));
    }
}
