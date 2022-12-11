use std::collections::HashMap;
use std::*;

use bairstow::rootfinding::{find_rootq, initial_guess, pbairstow_even};
fn test_rootfind() {
    let h = vec![5.0, 2.0, 9.0, 6.0, 2.0];
    let vr0s = initial_guess(h);
    let (_, niter, found) = pbairstow_even(h, vr0s);
    println!("{:?} ", vec![niter, found]);
    assert!(niter <= 5);
}
fn test_rootfind2() {
    let h = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
    let vr0s = initial_guess(h);
    let (vrs, niter, found) = pbairstow_even(h, vr0s);
    println!("{:?} ", vec![niter, found]);
    println!(
        "{:?} ",
        vrs.iter().map(|vr| find_rootq(vr)).collect::<Vec<_>>()
    );
    assert!(niter <= 12);
}
