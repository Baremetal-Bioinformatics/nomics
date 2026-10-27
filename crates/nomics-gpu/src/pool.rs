//! The process-wide host worker pool: CPU work around the GPU calls (packing
//! here, and alignment's per-block work in nomics-align) runs on persistent
//! threads instead of threads spawned per call.
//!
//! [`for_each`] keeps the scheduling of the per-call threads it replaces: at
//! most `threads` workers per call drain one shared queue, and the pool is
//! large enough (twice the CPUs by default, `NOMICS_HOST_POOL_THREADS`) that
//! concurrent calls (two batches in flight, two lanes per batch) all get
//! their workers at once and share the cores through the OS scheduler as
//! before. Jobs never wait on the GPU or on channels, only on the queue, so
//! an oversubscribed pool just queues workers.
use std::sync::{Mutex, OnceLock};
use std::time::Instant;

fn pool() -> &'static rayon::ThreadPool {
    static POOL: OnceLock<rayon::ThreadPool> = OnceLock::new();
    POOL.get_or_init(|| {
        let cpus = std::thread::available_parallelism().map_or(8, |n| n.get());
        let threads = std::env::var("NOMICS_HOST_POOL_THREADS")
            .ok()
            .and_then(|v| v.parse::<usize>().ok())
            .filter(|&n| n > 0)
            .unwrap_or(2 * cpus);
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .thread_name(|i| format!("nomics-host-{i}"))
            .build()
            .expect("host worker pool")
    })
}

/// Run `f` over `items` with at most `threads` workers on the host pool.
/// Items are independent, so the result does not depend on scheduling.
/// Returns the workers' summed busy time in seconds. With one worker (or at
/// most one item) everything runs inline on the caller.
pub fn for_each<T: Send>(threads: usize, items: Vec<T>, f: impl Fn(T) + Sync) -> f64 {
    let workers = threads.max(1).min(items.len());
    let started = Instant::now();
    if workers <= 1 {
        items.into_iter().for_each(f);
        return started.elapsed().as_secs_f64();
    }
    let queue = Mutex::new(items);
    let busy = Mutex::new(0f64);
    let (queue_ref, busy_ref, f) = (&queue, &busy, &f);
    // The caller blocks here like a thread-scope join; it never runs jobs.
    pool().in_place_scope(|s| {
        for _ in 0..workers {
            s.spawn(move |_| {
                let t = Instant::now();
                loop {
                    let item = queue_ref.lock().unwrap_or_else(|e| e.into_inner()).pop();
                    match item {
                        Some(item) => f(item),
                        None => break,
                    }
                }
                *busy_ref.lock().unwrap_or_else(|e| e.into_inner()) += t.elapsed().as_secs_f64();
            });
        }
    });
    busy.into_inner().unwrap_or_else(|e| e.into_inner())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicUsize, Ordering};

    #[test]
    fn every_item_runs_once_whatever_the_worker_count() {
        for threads in [0, 1, 3, 64] {
            for n in [0, 1, 2, 1000] {
                let seen: Vec<AtomicUsize> = (0..n).map(|_| AtomicUsize::new(0)).collect();
                let busy = for_each(threads, (0..n).collect(), |i: usize| {
                    seen[i].fetch_add(1, Ordering::Relaxed);
                });
                assert!(busy >= 0.0);
                assert!(seen.iter().all(|s| s.load(Ordering::Relaxed) == 1));
            }
        }
    }

    #[test]
    fn nested_and_concurrent_calls_complete() {
        let total = AtomicUsize::new(0);
        std::thread::scope(|s| {
            for _ in 0..4 {
                s.spawn(|| {
                    for_each(8, (0..50).collect(), |_: usize| {
                        for_each(3, (0..10).collect(), |_: usize| {
                            total.fetch_add(1, Ordering::Relaxed);
                        });
                    });
                });
            }
        });
        assert_eq!(total.load(Ordering::Relaxed), 4 * 50 * 10);
    }
}
