use std::collections::HashMap;
use std::*;

use bairstow::autocorr::{extract_autocorr, initial_autocorr, pbairstow_autocorr};
use bairstow::rootfinding::{find_rootq, Options};
fn test_autocorr() {
    let h = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
    let vr0s = initial_autocorr(h);
    let opts = Options();
    let (vrs, niter, found) = pbairstow_autocorr(h, vr0s, opts);
    println!("{:?} ", vec![niter, found]);
    for vr in vrs {
        vr = extract_autocorr(vr);
        println!("{:?} ", find_rootq(vr));
    }
    assert!(niter <= 13);
}
