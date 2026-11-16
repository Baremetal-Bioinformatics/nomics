//! Exclusive wall-time counters. Diagnostics only; never used by scientific code.
use std::{
    cell::RefCell,
    marker::PhantomData,
    rc::Rc,
    time::{Duration, Instant},
};
pub const PHASES: [&str; 8] = [
    "setup",
    "reading_filtering",
    "activity_regions",
    "assembly",
    "likelihoods",
    "genotyping_annotation",
    "f1r2",
    "writing",
];
thread_local! { static CLOCK: RefCell<Clock> = RefCell::new(Clock::default()); }
#[derive(Default)]
struct Clock {
    totals: [Duration; 8],
    stack: Vec<(usize, Instant, Duration)>,
}
pub(crate) struct Span(PhantomData<Rc<()>>);
pub(crate) fn phase(index: usize) -> Span {
    CLOCK.with(|c| {
        c.borrow_mut()
            .stack
            .push((index, Instant::now(), Duration::ZERO))
    });
    Span(PhantomData)
}
impl Drop for Span {
    fn drop(&mut self) {
        CLOCK.with(|c| {
            let mut c = c.borrow_mut();
            let (i, start, children) = c.stack.pop().expect("timing scope");
            let elapsed = start.elapsed();
            c.totals[i] += elapsed.saturating_sub(children);
            if let Some(parent) = c.stack.last_mut() {
                parent.2 += elapsed;
            }
        });
    }
}
pub(crate) fn take() -> [Duration; 8] {
    CLOCK.with(|c| std::mem::take(&mut c.borrow_mut().totals))
}
