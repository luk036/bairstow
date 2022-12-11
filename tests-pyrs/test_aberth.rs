use std::collections::HashMap;
use std::*;

use bairstow::aberth::{
    aberth, aberth_autocorr, initial_aberth, initial_aberth_autocorr, initial_aberth_autocorr_orig,
    initial_aberth_orig,
};
use bairstow::rootfinding::Options;
fn test_aberth1() {
    let h = vec![5.0, 2.0, 9.0, 6.0, 2.0];
    let z0s = initial_aberth(h);
    let (_, niter, found) = aberth(h, z0s);
    println!("{:?} ", vec![niter, found]);
    assert!(niter <= 8);
}
fn test_aberth2() {
    let h = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
    let z0s = initial_aberth(h);
    let (zs, niter, found) = aberth(h, z0s);
    println!("{:?} ", vec![niter, found]);
    println!("{:?} ", zs.iter().map(|z| z).collect::<Vec<_>>());
    assert!(niter <= 8);
}
static r: Vec<f64> = vec![
    -0.00196191,
    -0.00094597,
    -0.00023823,
    0.00134667,
    0.00380494,
    0.00681596,
    0.0097864,
    0.01186197,
    0.0121238,
    0.00985211,
    0.00474894,
    -0.00281751,
    -0.01173923,
    -0.0201885,
    -0.02590168,
    -0.02658216,
    -0.02035729,
    -0.00628271,
    0.01534627,
    0.04279982,
    0.0732094,
    0.10275561,
    0.12753013,
    0.14399228,
    0.15265722,
    0.14399228,
    0.12753013,
    0.10275561,
    0.0732094,
    0.04279982,
    0.01534627,
    -0.00628271,
    -0.02035729,
    -0.02658216,
    -0.02590168,
    -0.0201885,
    -0.01173923,
    -0.00281751,
    0.00474894,
    0.00985211,
    0.0121238,
    0.01186197,
    0.0097864,
    0.00681596,
    0.00380494,
    0.00134667,
    -0.00023823,
    -0.00094597,
    -0.00196191,
];
fn test_aberth_fir() {
    let z0s = initial_aberth(r);
    let opt = Options();
    opt.tol = 1e-08;
    let (zs, niter, found) = aberth(r, z0s, opt);
    println!("{:?} ", vec![niter, found]);
    for z in zs {
        println!("{:?} ", z);
    }
    assert!(niter <= 11);
}
fn test_aberth_autocorr_fir() {
    let z0s = initial_aberth_autocorr(r);
    let opt = Options();
    opt.tol = 1e-14;
    let (zs, niter, found) = aberth_autocorr(r, z0s, opt);
    println!("{:?} ", vec![niter, found]);
    for z in zs {
        println!("{:?} ", z);
    }
    assert!(niter <= 12);
}
fn test_aberth_fir_orig() {
    let z0s = initial_aberth_orig(r);
    let opt = Options();
    opt.tol = 1e-08;
    let (zs, niter, found) = aberth(r, z0s, opt);
    println!("{:?} ", vec![niter, found]);
    for z in zs {
        println!("{:?} ", z);
    }
    assert!(niter <= 13);
}
fn test_aberth_autocorr_fir_orig() {
    let z0s = initial_aberth_autocorr_orig(r);
    let opt = Options();
    opt.tol = 1e-14;
    let (zs, niter, found) = aberth_autocorr(r, z0s, opt);
    println!("{:?} ", vec![niter, found]);
    for z in zs {
        println!("{:?} ", z);
    }
    assert!(niter <= 12);
}
