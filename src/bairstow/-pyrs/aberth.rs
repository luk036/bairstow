use std::collections::HashMap;
use std::*;

use cmath::exp;
use math::pi;

use lds::Vdcorput;
use robin::Robin;
use rootfinding::{horner_eval, Options};
const PI: _ = pi;
fn horner_backward(pb: List, n: i32, alpha: complex) -> complex {
    // "[summary]
    //
    //     Args:
    //         pa (List[float]): [description]
    //         r (float): [description]
    //
    //     Returns:
    //         float: [description]
    //
    //     Examples:
    //         >>> p = [1.0, -6.7980, 2.9948, -0.043686, 0.000089248]
    //         >>> n = len(p) - 1
    //         >>> alpha = 6.3256
    //         >>> P = horner_backward(p, n, alpha)
    //         >>> -P * (alpha ** 5)
    //         -0.013355264987140483
    //         >>> p[3]
    //         0.006920331351966613
    //     "
    for i in (2..(n + 2)) {
        pb[-(i)] -= pb[-(i - 1)];
        pb[-(i)] /= -(alpha);
    }
    return pb[-(n + 1)];
}
fn initial_aberth(pa: Vec<complex>) -> Vec<complex> {
    // "[summary]
    //
    //     Args:
    //         pa (List): [description]
    //
    //     Returns:
    //         List: [description]
    //
    //     Examples:
    //         >>> h = [5.0, 2.0, 9.0, 6.0, 2.0]
    //         >>> z0s = initial_aberth(h)
    //     "
    let N: int = (pa.len() - 1);
    let c: complex = (-(pa[1]) / (N * pa[0]));
    let Pc: complex = horner_eval(pa.copy(), N, c);
    let re: complex = -(Pc).pow((1.0 / N));
    let z0s: Vec<complex> = vec![];
    let vgen = Vdcorput(2);
    vgen.reseed(1);
    for i in (0..N) {
        let vdc = ((2 * PI) * vgen.pop());
        z0s += vec![(c + (re*exp((vdc*1j))))];
    }
    return z0s;
}
fn initial_aberth_orig(pa: Vec<complex>) -> Vec<complex> {
    // "[summary]
    //
    //     Args:
    //         pa (List): [description]
    //
    //     Returns:
    //         List: [description]
    //
    //     Examples:
    //         >>> h = [5.0, 2.0, 9.0, 6.0, 2.0]
    //         >>> z0s = initial_aberth_orig(h)
    //     ";
    let N: int = (pa.len() - 1);
    let c: complex = (-(pa[1]) / (N * pa[0]));
    let Pc: complex = horner_eval(pa.copy(), N, c);
    let re: complex = -(Pc).pow((1.0 / N));
    let k = ((2 * PI) / N);
    let z0s: Vec<complex> = vec![];
    for i in (0..N) {
        let theta = (k * (0.25 + i));
        z0s += vec![(c + (re*exp((theta*1j))))];
    }
    return z0s;
}
fn aberth<T0>(pa: Vec<complex>, zs: Vec<complex>, options: T0) -> (Vec<complex>, i32, bool) {
    // "Aberth's method for polynomial root-finding
    //     Args:
    //         pa (List): [description]
    //         zs (List): [description]
    //         options (Options, optional): [description]. Defaults to Options().
    //
    //     Returns:
    //         [type]: [description]
    //
    //     Examples:
    //         >>> h = [5.0, 2.0, 9.0, 6.0, 2.0]
    //         >>> z0s = initial_aberth(h)
    //         >>> opt = Options()
    //         >>> opt.tol = 1e-8
    //         >>> zs, niter, found = aberth(h, z0s, opt)
    //     ";
    let M = zs.len();
    let N = (pa.len() - 1);
    let converged = (vec![false] * M);
    let robin = Robin(M);
    for niter in (1..options.max_iter) {
        let mut tol = 0.0;
        for i in (0..M).into_iter().filter(|i| !converged[i]) {
            let pb = pa.copy();
            let P = horner_eval(pb, N, zs[i]);
            let tol_i = abs(P);
            if tol_i < options.tol_ind {
                converged[i] = true;
                continue;
            }
            let mut P1 = horner_eval(pb, (N - 1), zs[i]);
            tol = tol_i.iter().max().unwrap();
            for j in robin.exclude(i) {
                P1 -= (P / (zs[i] - zs[j]));
            }
            zs[i] -= (P / P1);
        }
        if tol < options.tol {
            return (zs, niter, true);
        }
    }
    return (zs, options.max_iter, false);
}
fn initial_aberth_autocorr(pa: Vec<f32>) -> Vec<complex> {
    // "[summary]
    //
    //     Args:
    //         pa (List): [description]
    //
    //     Returns:
    //         List: [description]
    //
    //     Examples:
    //         >>> h = [5.0, 2.0, 9.0, 6.0, 2.0]
    //         >>> z0s = initial_aberth_autocorr(h)
    //     ";
    let N: int = (pa.len() - 1);
    let re: float = abs(pa[-1]).pow((1.0 / N));
    if abs(re) > 1 {
        let mut re = (1 / re);
    }
    N /= 2;
    let mut z0s = vec![];
    let vgen = Vdcorput(2);
    vgen.reseed(1);
    for i in (0..N) {
        let vdc = ((2 * PI) * vgen.pop());
        z0s += vec![(re*exp((vdc*1j)))];
    }
    return z0s;
}
fn initial_aberth_autocorr_orig(pa: Vec<f32>) -> Vec<complex> {
    // "[summary]
    //
    //     Args:
    //         pa (List): [description]
    //
    //     Returns:
    //         List: [description]
    //
    //     Examples:
    //         >>> h = [5.0, 2.0, 9.0, 6.0, 2.0]
    //         >>> z0s = initial_aberth_autocorr_orig(h)
    //     ";
    let N: int = (pa.len() - 1);
    let re: float = abs(pa[-1]).pow((1.0 / N));
    if abs(re) > 1 {
        let mut re = (1 / re);
    }
    N /= 2;
    let k = ((2 * PI) / N);
    let mut z0s = vec![];
    for i in (0..N) {
        let theta = (k * (0.25 + i));
        z0s += vec![(re*exp((theta*1j)))];
    }
    return z0s;
}
fn aberth_autocorr<T0>(pa: Vec<f32>, zs: Vec<complex>, options: T0) -> (Vec<complex>, i32, bool) {
    // "[summary]
    //
    //     Args:
    //         pa (List): [description]
    //         zs (List): [description]
    //         options (Options, optional): [description]. Defaults to Options().
    //
    //     Returns:
    //         [type]: [description]
    //
    //     Examples:
    //         >>> h = [5.0, 2.0, 9.0, 6.0, 2.0]
    //         >>> z0s = initial_aberth_autocorr(h)
    //         >>> zs, niter, found = aberth_autocorr(h, z0s)
    //         >>> opt = Options()
    //         >>> opt.tol = 1e-8
    //         >>> zs, niter, found = aberth_autocorr(h, z0s, opt)
    //     ";
    let M: int = zs.len();
    let N: int = (pa.len() - 1);
    let converged: Vec<bool> = (vec![false] * M);
    let robin = Robin(M);
    for niter in (1..options.max_iter) {
        let tol: float = 0.0;
        for i in (0..M).into_iter().filter(|i| !converged[i]) {
            let pb = pa.copy();
            let P = horner_eval(pb, N, zs[i]);
            let tol_i = abs(P);
            if tol_i < options.tol_ind {
                converged[i] = true;
                continue;
            }
            let mut P1 = horner_eval(pb, (N - 1), zs[i]);
            let mut tol = tol_i.iter().max().unwrap();
            for j in robin.exclude(i) {
                P1 -= (P / (zs[i] - zs[j]));
                let zsn = (1.0 / zs[j]);
                P1 -= (P / (zs[i] - zsn));
            }
            zs[i] -= (P / P1);
        }
        if tol < options.tol {
            return (zs, niter, true);
        }
    }
    return (zs, options.max_iter, false);
}
