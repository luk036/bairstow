use std::collections::HashMap;
use std::*;

use math::{acos, cos, sqrt};

use lds::Vdcorput;
use matrix2::Matrix2;
use robin::Robin;
use vector2::Vector2;
const PI: _ = acos(-1.0);
struct Options {
    max_iter: int,
    tol: f64,
    tol_ind: f64,
}

// impl Options {
// let max_iter: int = 2000;
// let tol: f64 = 1e-12;
// let tol_ind: f64 = 1e-15;
// }
fn delta1(vA: Vector2, vr: Vector2, vp: Vector2) -> Vector2 {
    // "for ri - rj
    //
    //     Args:
    //         vA (Vector2): [description]
    //         vr (Vector2): [description]
    //         vp (Vector2): [description]
    //
    //     Returns:
    //         Vector2: [description]
    //
    //     Examples:
    //         >>> d = delta1(Vector2(1, 2), Vector2(-2, 0), Vector2(4, -5))
    //         >>> print(d)
    //         <0.2, 0.4>
    //     ";
    let (r, q) = (vr.x, vr.y);
    let (p, s) = (vp.x, vp.y);
    let mp = Matrix2(Vector2(-(s), -(p)), Vector2((p * q), ((p * r) - s)));
    return (mp.mdot(vA) / mp.det());
}
fn delta2(vA: Vector2, vr: Vector2, vA1: Vector2) -> Vector2 {
    // "for A1
    //
    //     Args:
    //         vA (Vector2): [description]
    //         vr (Vector2): [description]
    //         vp (Vector2): [description]
    //
    //     Returns:
    //         Vector2: [description]
    //
    //     Examples:
    //         >>> d = delta2(Vector2(1, 2), Vector2(-2, 0), Vector2(4, 5))
    //         >>> print(d)
    //         <0.2, -0.4>
    //     ";
    let (r, q) = (vr.x, vr.y);
    let (p, s) = (vA1.x, vA1.y);
    let mp = Matrix2(Vector2(-(s), p), Vector2((p * q), ((p * r) + s)));
    return (mp.mdot(vA) / mp.det());
}
fn delta(vA: Vector2, vr: Vector2, vp: Vector2) -> Vector2 {
    // "for -vA1
    //
    //     Args:
    //         vA (Vector2): [description]
    //         vr (Vector2): [description]
    //         vp (Vector2): [description]
    //
    //     Returns:
    //         Vector2: [description]
    //
    //     Examples:
    //         >>> d = delta(Vector2(1, 2), Vector2(-2, 0), Vector2(4, -5))
    //         >>> print(d)
    //         <0.2, -0.4>
    //     ";
    let (r, q) = (vr.x, vr.y);
    let (p, s) = (vp.x, vp.y);
    let mp = Matrix2(Vector2(s, p), Vector2((p * q), ((p * r) - s)));
    return (mp.mdot(vA) / mp.det());
}
fn suppress_old(vA: Vector2, vA1: Vector2, vri: Vector2, vrj: Vector2) {
    // "[summary]
    //
    //     Args:
    //         vA (Vector2): [description]
    //         vr (Vector2): [description]
    //         vp (Vector2): [description]
    //
    //     Returns:
    //         Vector2: [description]
    //
    //     Examples:
    //         >>> vA = Vector2(3, 3)
    //         >>> vA1 = Vector2(1, 2)
    //         >>> vri = Vector2(-2, 0)
    //         >>> vrj = Vector2(4, -5)
    //         >>> suppress_old(vA, vA1, vri, vrj)
    //         >>> dr = delta(vA, vri, Vector2(vA1._x, -vA1._y))
    //         >>> print(dr)
    //         <-16.780821917808325, -1.4383561643835612>
    //     ";
    let (A, B) = (vA.x, vA.y);
    let (A1, B1) = (vA1.x, vA1.y);
    let vp = (vri - vrj);
    let (r, q) = (vri.x, -(vri.y));
    let (p, s) = (vp.x, -(vp.y));
    let f = ((r * p) + s);
    let qp = (q * p);
    let e = ((f * s) - (qp * p));
    let a = (((A * s) - (B * p)) / e);
    let b = (((B * f) - (A * qp)) / e);
    let c = (A1 - a);
    let d = ((B1 - b) - (a * p));
    vA._x = a;
    vA._y = b;
    vA1._x = (((c * s) - (d * p)) / e);
    vA1._y = (((d * f) - (c * qp)) / e);
}
fn suppress(vA: Vector2, vA1: Vector2, vri: Vector2, vrj: Vector2) {
    // "[summary]
    //
    //     Args:
    //         vA (Vector2): [description]
    //         vr (Vector2): [description]
    //         vp (Vector2): [description]
    //
    //     Returns:
    //         Vector2: [description]
    //
    //     Examples:
    //         >>> vA = Vector2(3, 3)
    //         >>> vA1 = Vector2(1, 2)
    //         >>> vri = Vector2(-2, 0)
    //         >>> vrj = Vector2(4, -5)
    //         >>> suppress(vA, vA1, vri, vrj)
    //         >>> dr = delta2(vA, vri, vA1)
    //         >>> print(dr)
    //         <-16.780821917808325, -1.4383561643835612>
    //     ";
    let vp = (vri - vrj);
    let vAnew = delta1(vA, vri, vp);
    vA._x = vAnew._x;
    vA._y = vAnew._y;
    vA1._x -= vA._x;
    vA1._y -= ((vA._x * vp._x) + vA._y);
    let vA1new = delta1(vA1, vri, vp);
    vA1._x = vA1new._x;
    vA1._y = vA1new._y;
}
fn horner_eval<T0, RT>(coeffs: List, degree: i32, zval: T0) -> RT {
    // "Polynomial evaluation using Horner's scheme
    //
    //     Note: coeffs becomes the quotient after calling this function
    //
    //     Args:
    //         coeffs (List): List of coefficients of polynomial
    //         val (f64 or complex): value to be evaluated
    //
    //     Returns:
    //         f64 or complex
    //
    //     Examples:
    //         >>> coeffs = [1, -8, -72, 382, 727, -2310]
    //         >>> horner_eval(coeffs, 5, 3)
    //         960
    //         >>> coeffs
    //         [1, -5, -87, 121, 1090, 960]
    //     ";
    for i in (0..degree) {
        coeffs[(i + 1)] += (coeffs[i] * zval);
    }
    return coeffs[degree];
}
fn horner_backward<T0, RT>(coeffs: List, degree: i32, val: T0) -> RT {
    // "Polynomial evaluation using Horner's scheme
    //
    //     Note: coeffs becomes the quotient after calling this function
    //
    //     Args:
    //         coeffs (List): List of coefficients of polynomial
    //         val (f64 or complex): value to be evaluated
    //
    //     Returns:
    //         f64 or complex
    //
    //     Examples:
    //         >>> coeffs = [1.0, -6.7980, 2.9948, -0.043686, 0.000089248]
    //         >>> n = len(coeffs) - 1
    //         >>> alpha = 6.3256
    //         >>> P = horner_backward(coeffs, 4, alpha)
    //         >>> -P * (alpha ** 5)
    //         -0.013355264987140483
    //         >>> coeffs[3]
    //         0.006920331351966613
    //     ";
    for i in (2..(degree + 2)) {
        coeffs[-(i)] -= coeffs[-(i - 1)];
        coeffs[-(i)] /= -(val);
    }
    return coeffs[-(degree + 1)];
}
fn horner(coeffs: Vec<f32>, degree: i32, vr: Vector2) -> Vector2 {
    // "[summary]
    //
    //     Note: pb becomes the quotient after calling this function
    //
    //     Args:
    //         coeffs (List[f64]): [description]
    //         vr (Vector2): [description]
    //
    //     Returns:
    //         Vector2: [description]
    //
    //     Examples:
    //         >>> coeffs = [1, -8, -72, 382, 727, -2310]
    //         >>> vp = horner(coeffs, 5, Vector2(-1, -6))  # x^2 + x - 6
    //         >>> coeffs
    //         [1, -9, -57, 385, 0, 0]
    //         >>> coeffs = [1, -8, -72, 382, 727, -2310]
    //         >>> vp = horner(coeffs, 5, Vector2(2, -3))  # x^2 - 2x - 3
    //         >>> coeffs
    //         [1, -6, -81, 202, 888, -1704]
    //     ";
    coeffs[1] += (coeffs[0] * vr.x);
    for i in (2..degree) {
        coeffs[i] += ((coeffs[(i - 1)] * vr.x) - (coeffs[(i - 2)] * vr.y));
    }
    coeffs[degree] -= (coeffs[(degree - 2)] * vr.y);
    return Vector2(coeffs[(degree - 1)], coeffs[degree]);
}
fn initial_guess_orig(coeffs: Vec<f32>) -> Vec<Vector2> {
    // "[summary]
    //
    //     Args:
    //         pa (List[f64]): [description]
    //
    //     Returns:
    //         List[Vector2]: [description]
    //
    //     Examples:
    //         >>> h = [10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0]
    //         >>> vr0s = initial_guess(h)
    //     ";
    let mut degree = (coeffs.len() - 1);
    let center = (-(coeffs[1]) / (degree * coeffs[0]));
    let Pc = horner_eval(coeffs.copy(), degree, center);
    let reff = abs(Pc).pow((1 / degree));
    let m = ((center * center) + (reff * reff));
    let mut vr0s = vec![];
    degree /= 2;
    degree *= 2;
    let k = (PI / degree);
    for i in (1..degree).step_by(2) {
        let temp = (reff * cos((k * i)));
        let r0 = (2 * (center + temp));
        let t0 = (m + ((2 * center) * temp));
        vr0s += vec![Vector2(r0, t0)];
    }
    return vr0s;
}
fn initial_guess(coeffs: Vec<f32>) -> Vec<Vector2> {
    // "[summary]
    //
    //     Args:
    //         pa (List[f64]): [description]
    //
    //     Returns:
    //         List[Vector2]: [description]
    //
    //     Examples:
    //         >>> h = [10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0]
    //         >>> vr0s = initial_guess(h)
    //     ";
    let mut degree = (coeffs.len() - 1);
    let center = (-(coeffs[1]) / (degree * coeffs[0]));
    let Pc = horner_eval(coeffs.copy(), degree, center);
    let reff = abs(Pc).pow((1 / degree));
    let m = ((center * center) + (reff * reff));
    let mut vr0s = vec![];
    degree /= 2;
    degree *= 2;
    let vgen = Vdcorput(2);
    vgen.reseed(1);
    for i in (1..degree).step_by(2) {
        let temp = (reff * cos((PI * vgen.pop())));
        let r0 = (2 * (center + temp));
        let t0 = (m + ((2 * center) * temp));
        vr0s += vec![Vector2(r0, t0)];
    }
    return vr0s;
}
fn pbairstow_even<T0>(pa: Vec<f32>, vrs: Vec<Vector2>, options: T0) -> (Vec<Vector2>, i32, bool) {
    // "Parallel Bairstow's method
    //
    //     Args:
    //         pa (List[f64]): [description]
    //         vrs (List[Vector2]): [description]
    //         options (Options, optional): [description]. Defaults to Options().
    //
    //     Returns:
    //         [type]: [description]
    //
    //     Examples:
    //         >>> h = [10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0]
    //         >>> vr0s = initial_guess(h)
    //         >>> vrs, niter, found = pbairstow_even(h, vr0s)
    //     ";
    let M = vrs.len();
    let N = (pa.len() - 1);
    let converged = (vec![false] * M);
    let robin = Robin(M);
    for niter in (1..options.max_iter) {
        let mut tol = 0.0;
        for i in (0..M).into_iter().filter(|i| converged[i] == false) {
            let pb = pa.copy();
            let vA = horner(pb, N, vrs[i]);
            let tol_i = abs(vA.x).iter().max().unwrap();
            if tol_i < options.tol_ind {
                converged[i] = true;
                continue;
            }
            let vA1 = horner(pb, (N - 2), vrs[i]);
            tol = tol_i.iter().max().unwrap();
            for j in robin.exclude(i) {
                suppress_old(vA, vA1, vrs[i], vrs[j]);
            }
            vrs[i] -= delta2(vA, vrs[i], vA1);
        }
        if tol < options.tol {
            return (vrs, niter, true);
        }
    }
    return (vrs, options.max_iter, false);
}
fn find_rootq(vr: Vector2) -> (f32, f32) {
    // "Solve x^2 - r*x + q = 0
    //
    //     (x - x1)(x - x2) = x^2 - (x1 + x2) x + x1 * x2
    //
    //     Args:
    //         vr (Vector2): [description]
    //
    //     Returns:
    //         Tuple[f64, f64]: [description]
    //
    //     Examples:
    //         >>> vr = find_rootq(Vector2(5, 6))
    //         >>> print(vr)
    //         (3.0, 2.0)
    //     ";
    let hr = (vr.x / 2);
    let d = ((hr * hr) - vr.y);
    if d < 0 {
        // x1 = (hr + (sqrt(-(d))*1j));
    } else {
        x1 = (hr + if hr >= 0 { sqrt(d) } else { -sqrt(d) });
    }
    let x2 = (vr.y / x1);
    return (x1, x2);
}
