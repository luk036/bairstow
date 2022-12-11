use std::collections::HashMap;
use std::*;

use bairstow::aberth::{horner_backward, horner_eval};
use pytest::approx;
fn test_backward2() {
    let p = vec![1.0, -6.798, 2.9948, -0.043686, 8.9248e-05];
    let n = (p.len() - 1);
    let alpha = 6.3256;
    let P = horner_backward(p, n, alpha);
    assert!((-(P) * alpha.pow(5)) == approx(-0.013355264987140483));
    assert!(p[3] == approx(0.006920331351966613));
}
fn test_backward1() {
    let p = vec![1.0, -6.798, 2.9948, -0.043686, 8.9248e-05];
    let n = (p.len() - 1);
    let P = horner_eval(p, n, 6.3256);
    assert!(P == approx(-0.012701469838522064));
    assert!(p[3] == approx(-0.0020220560640132265));
}
