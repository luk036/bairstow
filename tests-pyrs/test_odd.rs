use std::collections::HashMap;
use std::*;

use bairstow::rootfinding::{find_rootq, initial_guess, pbairstow_even};
fn test_odd2() {
    let h = vec![5.0, 2.5, 9.2, 6.9, 2.6, 0.2, 0];
    let vr0s = initial_guess(h);
    let (vrs, niter, found) = pbairstow_even(h, vr0s);
    println!("{:?} ", vec![niter, found]);
    println!(
        "{:?} ",
        vrs.iter().map(|vr| find_rootq(vr)).collect::<Vec<_>>()
    );
    assert!(niter <= 114);
}
