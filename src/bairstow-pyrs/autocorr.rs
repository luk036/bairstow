use std::collections::HashMap;
use std::*;

use math::{acos, cos, sqrt};

use lds::Vdcorput;
use robin::Robin;
use rootfinding::{delta1, delta2, horner, suppress_old, Options};
use vector2::Vector2;
const PI: _ = acos(-1.0);
fn initial_autocorr_new(pa: Vec<f32>) -> Vec<Vector2> {
    "[summary]

    Args:
        pa (List[float]): [description]

    Returns:
        List[Vector2]: [description]

    Examples:
        >>> h = [10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0]
        >>> vr0s = initial_autocorr(h)
    ";
    let mut N = (pa.len() - 1);
    let mut re = abs(pa[-1]).pow((1.0 / N));
    if re > 1 {
        re = (1 / re);
    }
    N /= 2;
    let m = (re * re);
    let vgen = Vdcorput(2);
    vgen.reseed(1);
    let vr0s = (1..N)
        .step_by(2)
        .iter()
        .map(|_| Vector2(((2 * re) * cos((PI * vgen.pop()))), m))
        .collect::<Vec<_>>();
    return vr0s;
}
fn initial_autocorr(pa: Vec<f32>) -> Vec<Vector2> {
    "[summary]

    Args:
        pa (List[float]): [description]

    Returns:
        List[Vector2]: [description]
    ";
    let mut N = (pa.len() - 1);
    let mut re = abs(pa[-1]).pow((1.0 / N));
    if re < 1 {
        re = (1 / re);
    }
    N /= 2;
    let k = (PI / N);
    let m = (re * re);
    let vr0s = (1..N)
        .step_by(2)
        .iter()
        .map(|i| Vector2(((2 * re) * cos((k * i))), m))
        .collect::<Vec<_>>();
    return vr0s;
}
fn pbairstow_autocorr<RT>(pa: Vec<f32>, vrs: Vec<Vector2>, options: Options) -> RT {
    "[summary]

    Args:
        pa (List[float]): [description]
        vrs (List[Vector2]): [description]
        options (Options, optional): [description]. Defaults to Options().

    Returns:
        [type]: [description]

    Examples:
        >>> h = [10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0]
        >>> vr0s = initial_autocorr(h)
        >>> vrs, niter, found = pbairstow_autocorr(h, vr0s)
    ";
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
            tol = tol.iter().max().unwrap();
            let vA1 = horner(pb, (N - 2), vrs[i]);
            for j in robin.exclude(i) {
                suppress_old(vA, vA1, vrs[i], vrs[j]);
                let vrn = (Vector2(vrs[j].x, 1.0) / vrs[j].y);
                suppress_old(vA, vA1, vrs[i], vrn);
            }
            vrs[i] -= delta2(vA, vrs[i], vA1);
        }
        if tol < options.tol {
            return (vrs, niter, true);
        }
    }
    return (vrs, options.max_iter, false);
}
fn extract_autocorr(vr: Vector2) -> Vector2 {
    "Extract the quadratic function where its roots are within a unit circle

    x^2 - r*x + q  or (1/q) - (r/q) * x + x^2
    (x - a1)(x - a2) = x^2 - (a1 + a2) x + a1 * a2

    Args:
        vr (Vector2): [description]

    Returns:
        Vector2: [description]

    Examples:
        >>> vr = extract_autocorr(Vector2(1, 4))
        >>> print(vr)
        <0.25, 0.25>
    ";
    let (r, q) = (vr.x, vr.y);
    let hr = (r / 2.0);
    let d = ((hr * hr) - q);
    if d < 0.0 {
        if q > 1.0 {
            vr = (Vector2(r, 1.0) / q);
        }
    } else {
        let mut a1 = (hr + if hr >= 0.0 { sqrt(d) } else { -sqrt(d) });
        let mut a2 = (q / a1);
        if abs(a1) > 1.0 {
            if abs(a2) > 1.0 {
                a2 = (1.0 / a2);
            }
            a1 = (1.0 / a1);
            vr = Vector2((a1 + a2), (a1 * a2));
        } else {
            if abs(a2) > 1.0 {
                a2 = (1.0 / a2);
                vr = Vector2((a1 + a2), (a1 * a2));
            }
        }
    }
    return vr;
}
