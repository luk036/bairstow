use std::collections::HashMap;
use std::*;

use __future__::print_function;
use bairstow::autocorr::{initial_autocorr, pbairstow_autocorr};
use bairstow::rootfinding::{initial_guess, pbairstow_even};
fn run_autocorr<RT>() -> RT {
    "[summary]

    Returns:
        [type]: [description]
    ";
    let p = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
    let vr0s = initial_autocorr(p);
    let (_, niter, _) = pbairstow_autocorr(p, vr0s);
    return niter;
}
fn run_pbairstow<RT>() -> RT {
    "[summary]

    Returns:
        [type]: [description]
    ";
    let p = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
    let vr0s = initial_guess(p);
    let (_, niter, _) = pbairstow_even(p, vr0s);
    return niter;
}
fn test_autocorr<T0>(benchmark: T0) {
    "[summary]

    Arguments:
        benchmark ([type]): [description]
    ";
    let result = benchmark(run_autocorr);
    assert!(result <= 13);
}
fn test_pbairstow<T0>(benchmark: T0) {
    "[summary]

    Arguments:
        benchmark ([type]): [description]
    ";
    let result = benchmark(run_pbairstow);
    assert!(result <= 12);
}
