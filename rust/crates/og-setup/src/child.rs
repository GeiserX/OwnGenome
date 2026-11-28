//! Every process og's setup starts goes through here, so a signal that ends
//! og can stop them before an installer lock goes ([`crate::lock`]).
//!
//! The shell runs its `trap 'lock_release ...' EXIT` only after the
//! command it is waiting for has exited, so a curl killed with setup never
//! writes into a `.part` file a new run has locked. og does the same: on
//! SIGINT, SIGTERM or SIGHUP no new child starts, each running child gets
//! the signal, and the locks go once every one of them has been reaped (a
//! child still there after 30 seconds gets SIGKILL). A thread whose child
//! ended while og is stopping waits there for og to exit, so it cannot
//! print a failure or exit with another status first.

use std::io;
use std::process::{Child, Command, ExitStatus, Output, Stdio};
use std::sync::{Condvar, Mutex, MutexGuard, PoisonError};
use std::time::{Duration, Instant};

struct Running {
    /// A signal is ending og: nothing starts any more.
    stopping: bool,
    /// The children started and not yet reaped.
    pids: Vec<u32>,
}

static RUNNING: Mutex<Running> = Mutex::new(Running {
    stopping: false,
    pids: Vec::new(),
});
/// Notified each time a child leaves the list.
static LEFT: Condvar = Condvar::new();

fn running() -> MutexGuard<'static, Running> {
    RUNNING.lock().unwrap_or_else(PoisonError::into_inner)
}

/// Wait here until og exits (the signal thread ends the process).
fn hold_forever() -> ! {
    loop {
        std::thread::park();
    }
}

/// A child in the list, until [`Ticket::done`] after its wait.
#[derive(Debug)]
pub struct Ticket(u32);

impl Ticket {
    /// The child has been reaped: take it off the list. When a signal is
    /// ending og, this thread goes no further.
    pub fn done(self) {
        drop(self);
        if running().stopping {
            hold_forever();
        }
    }
}

impl Drop for Ticket {
    fn drop(&mut self) {
        let mut r = running();
        r.pids.retain(|&p| p != self.0);
        LEFT.notify_all();
    }
}

/// `cmd.spawn()`, the child listed until its [`Ticket`] is done. Once a
/// signal is ending og, no child starts: the caller waits here for og to exit.
pub fn spawn(cmd: &mut Command) -> io::Result<(Child, Ticket)> {
    let mut r = running();
    if r.stopping {
        drop(r);
        hold_forever();
    }
    let child = cmd.spawn()?;
    let pid = child.id();
    r.pids.push(pid);
    Ok((child, Ticket(pid)))
}

/// `cmd.status()`: stdin, stdout and stderr as the caller set them, else
/// og's own.
pub fn status(cmd: &mut Command) -> io::Result<ExitStatus> {
    let (mut child, ticket) = spawn(cmd)?;
    let status = child.wait();
    ticket.done();
    status
}

/// `cmd.output()` for a caller that sets stdin and stderr: stdout is captured.
pub fn output(cmd: &mut Command) -> io::Result<Output> {
    cmd.stdout(Stdio::piped());
    let (child, ticket) = spawn(cmd)?;
    let output = child.wait_with_output();
    ticket.done();
    output
}

/// For the signal thread: start nothing more, send SIG to every running
/// child, and return once each has been reaped. One still there after 30
/// seconds gets SIGKILL, and is waited for 10 seconds more.
pub(crate) fn stop_all(sig: i32) {
    let mut r = running();
    r.stopping = true;
    send(&r.pids, sig);
    r = wait_for_none(r, Duration::from_secs(30));
    if !r.pids.is_empty() {
        send(&r.pids, signal_hook::consts::SIGKILL);
        drop(wait_for_none(r, Duration::from_secs(10)));
    }
}

fn send(pids: &[u32], sig: i32) {
    let Some(sig) = rustix::process::Signal::from_named_raw(sig) else {
        return;
    };
    for &p in pids {
        if let Some(pid) = i32::try_from(p)
            .ok()
            .and_then(rustix::process::Pid::from_raw)
        {
            let _ = rustix::process::kill_process(pid, sig);
        }
    }
}

fn wait_for_none(
    mut r: MutexGuard<'static, Running>,
    limit: Duration,
) -> MutexGuard<'static, Running> {
    let end = Instant::now() + limit;
    while !r.pids.is_empty() {
        let left = end.saturating_duration_since(Instant::now());
        if left.is_zero() {
            break;
        }
        r = LEFT
            .wait_timeout(r, left)
            .unwrap_or_else(PoisonError::into_inner)
            .0;
    }
    r
}
