use std::collections::HashMap;
use std::*;

use bairstow::autocorr::{extract_autocorr, initial_autocorr, pbairstow_autocorr};
use bairstow::rootfinding::{find_rootq, initial_guess, pbairstow_even, Options};
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
fn test_fir_even() {
    let vr0s = initial_guess(r);
    let opts = Options();
    opts.tol = 1e-07;
    let (vrs, niter, found) = pbairstow_even(r, vr0s, opts);
    println!("{:?} ", vec![niter, found]);
    for vr in vrs {
        println!("{:?} ", find_rootq(vr));
    }
    assert!(niter <= 119);
}
fn test_fir_auto() {
    let vr0s = initial_autocorr(r);
    println!("{:?} ", "vrs: {}".format(vr0s.len()));
    let opts = Options();
    opts.tol = 1e-07;
    let (vrs, niter, found) = pbairstow_autocorr(r, vr0s, opts);
    println!("{:?} ", vec![niter, found]);
    for vr in vrs {
        vr = extract_autocorr(vr);
        println!("{:?} ", find_rootq(vr));
    }
    assert!(niter <= 16);
}
