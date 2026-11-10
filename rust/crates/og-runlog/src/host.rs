//! The machine `og run` runs on: its cores, memory, host name and file
//! systems. Behind a trait, so a test can stand in a box of its own.

use std::io;
use std::path::Path;

/// What `statvfs(2)` says of a file system.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct FsStat {
    pub frsize: u64,
    pub bsize: u64,
    pub blocks: u64,
    pub bfree: u64,
    pub bavail: u64,
    pub files: u64,
    pub ffree: u64,
}

/// The machine.
pub trait Host {
    /// `statvfs(2)` of the file system holding `path`.
    fn statvfs(&self, path: &Path) -> io::Result<FsStat>;
    /// The device id of `path`.
    fn dev(&self, path: &Path) -> io::Result<u64>;
    /// Whether `path` exists (a broken symlink does not).
    fn exists(&self, path: &Path) -> bool;
    /// The RAM this system sees, in bytes (inside WSL2 or a VM: the VM's).
    fn physical_memory(&self) -> Option<u128>;
}

/// The real machine.
pub struct System;

impl Host for System {
    fn statvfs(&self, path: &Path) -> io::Result<FsStat> {
        let s = rustix::fs::statvfs(path)?;
        Ok(FsStat {
            frsize: s.f_frsize,
            bsize: s.f_bsize,
            blocks: s.f_blocks,
            bfree: s.f_bfree,
            bavail: s.f_bavail,
            files: s.f_files,
            ffree: s.f_ffree,
        })
    }

    fn dev(&self, path: &Path) -> io::Result<u64> {
        use std::os::unix::fs::MetadataExt;
        Ok(std::fs::metadata(path)?.dev())
    }

    fn exists(&self, path: &Path) -> bool {
        path.exists()
    }

    /// `MemTotal` of `/proc/meminfo`, the `sysinfo(2)` total RAM that
    /// `sysconf(_SC_PHYS_PAGES)` counts. `None` where there is no such
    /// file (og runs on Linux; elsewhere the memory line says unknown).
    fn physical_memory(&self) -> Option<u128> {
        let text = std::fs::read_to_string("/proc/meminfo").ok()?;
        let line = text.lines().find(|l| l.starts_with("MemTotal:"))?;
        let kb: u128 = line.split_whitespace().nth(1)?.parse().ok()?;
        Some(kb * 1024)
    }
}

/// The CPUs this system has online: `sysconf(_SC_NPROCESSORS_ONLN)`, what
/// Python's `os.cpu_count()` returns. On Linux that is the list in
/// `/sys/devices/system/cpu/online` (a container's CPU quota does not
/// change it); elsewhere the standard library's count.
pub fn cpu_count() -> usize {
    if let Ok(list) = std::fs::read_to_string("/sys/devices/system/cpu/online")
        && let Some(n) = count_cpu_list(list.trim())
    {
        return n;
    }
    std::thread::available_parallelism().map_or(1, usize::from)
}

/// The number of CPUs in a list such as `0-3,6,8-9`.
fn count_cpu_list(list: &str) -> Option<usize> {
    let mut n = 0;
    for part in list.split(',') {
        n += match part.split_once('-') {
            Some((a, b)) => {
                b.parse::<usize>()
                    .ok()?
                    .checked_sub(a.parse::<usize>().ok()?)?
                    + 1
            }
            None => {
                part.parse::<usize>().ok()?;
                1
            }
        };
    }
    (n > 0).then_some(n)
}

/// The short host name: the system's node name up to its first `.`.
pub fn short_hostname() -> String {
    let name = rustix::system::uname();
    let node = name.nodename().to_string_lossy();
    node.split('.').next().unwrap_or("").to_owned()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cpu_lists() {
        assert_eq!(count_cpu_list("0-3"), Some(4));
        assert_eq!(count_cpu_list("0-3,6,8-9"), Some(7));
        assert_eq!(count_cpu_list("0"), Some(1));
        assert_eq!(count_cpu_list(""), None);
        assert!(cpu_count() >= 1);
    }
}
