from math import acos, cos, pow, sqrt

from .rootfinding import Options, delta, horner, makeadjoint
from .vector2 import vector2

PI = acos(-1.0)


def initial_autocorr(pa):
    """[summary]

    Args:
        pa ([type]): [description]

    Returns:
        [type]: [description]
    """
    N = len(pa) - 1
    re = pow(abs(pa[-1]), 1.0 / N)
    N //= 2
    k = PI / N
    m = re * re
    vr0s = [vector2(2 * re * cos(k * i), m) for i in range(1, N, 2)]
    return vr0s


def pbairstow_autocorr(pa, vrs, options=Options()):
    """[summary]

    Args:
        pa ([type]): [description]
        vrs ([type]): [description]
        options ([type], optional): [description]. Defaults to Options().

    Returns:
        [type]: [description]
    """
    M = len(vrs)  # assume polynomial of h is even
    found = False
    converged = [False] * M
    for niter in range(options.max_iter):
        tol = 0.0
        for i in filter(lambda i: converged[i] is False, range(M)):  # exclude converged
            vA, pb = horner(pa, vrs[i])
            tol_i = max(abs(vA.x), abs(vA.y))
            if tol_i < options.tol:
                converged[i] = True
                continue
            tol = max(tol, tol_i)
            vA1, _ = horner(pb, vrs[i])
            for j in filter(lambda j: j != i, range(M)):  # exclude i
                vA1 -= delta(vA, vrs[j], vrs[i] - vrs[j])
            for j in range(M):
                vrn = vector2(vrs[j].x, 1.0) / vrs[j].y
                vA1 -= delta(vA, vrn, vrs[i] - vrn)
            vrs[i] -= delta(vA, vrs[i], vA1)
        if tol < options.tol:
            found = True
            break
    return vrs, niter + 1, found


def extract_autocorr(vr):
    """Extract the quadratic function where its roots are within a unit circle

    x^2 - r*x + t  or x^2 - (r/t) * x + (1/t)

    (x - x1)(x - x2) = x^2 - (x1 + x2) x + x1 * x2

    determinant r/2 + q

    Args:
        vr ([type]): [description]

    Returns:
        [type]: [description]
    """
    r, t = vr.x, vr.y
    hr = r / 2.0
    d = hr * hr - t
    if d < 0.0:  # complex conjugate root
        if t > 0.0:
            vr = vector2(r, 1.0) / t
    else:
        # two real roots
        x1 = hr + (sqrt(d) if hr >= 0.0 else -sqrt(d))
        x2 = t / x1
        if abs(x1) > 1.0:
            if abs(x2) > 1.0:
                x2 = 1.0 / x2
            x1 = 1.0 / x1
            vr = vector2(x1 + x2, x1 * x2)
        elif abs(x2) > 1.0:
            x2 = 1.0 / x2
            vr = vector2(x1 + x2, x1 * x2)

    return vr
