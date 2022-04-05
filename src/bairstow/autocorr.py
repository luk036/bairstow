from math import acos, cos, sqrt
from typing import List

from .aberth import Options
from .rootfinding import delta, horner
from .vector2 import vector2

PI = acos(-1.0)


def initial_autocorr(pa: List[float]) -> List[vector2]:
    """[summary]

    Args:
        pa (List[float]): [description]

    Returns:
        List[vector2]: [description]
    """
    N = len(pa) - 1
    re = abs(pa[-1]) ** (1.0 / N)
    if re > 1:
        re = 1 / re
    N //= 2
    k = PI / N
    m = re * re
    vr0s = [vector2(-2 * re * cos(k * i), m) for i in range(1, N, 2)]
    return vr0s


def pbairstow_autocorr(
    pa: List[float], vrs: List[vector2], options: Options = Options()
):
    """[summary]

    Args:
        pa (List[float]): [description]
        vrs (List[vector2]): [description]
        options (Options, optional): [description]. Defaults to Options().

    Returns:
        [type]: [description]
    """
    M = len(vrs)  # assume polynomial of h is even
    N = len(pa) - 1
    converged = [False] * M
    for niter in range(1, options.max_iter):
        tol = 0.0
        # found = True  # initial
        for i in filter(lambda i: converged[i] is False, range(M)):  # exclude converged
            pb = pa.copy()
            vA = horner(pb, N, vrs[i])
            tol_i = max(abs(vA.x), abs(vA.y))
            if tol_i < options.tol_ind:
                converged[i] = True
                continue
            tol = max(tol, tol_i)
            found = False
            vA1 = horner(pb, N - 2, vrs[i])
            for j in filter(lambda j: j != i, range(M)):  # exclude i
                vA1 -= delta(vA, vrs[j], vrs[i] - vrs[j])
            for j in range(M):
                vrn = vector2(vrs[j].x, 1.0) / vrs[j].y
                vA1 -= delta(vA, vrn, vrs[i] - vrn)
            vrs[i] -= delta(vA, vrs[i], vA1)
            if vrs[i].y > 1:
                vrs[i] = vector2(vrs[i].x, 1.0) / vrs[i].y
        if tol < options.tol:
            return vrs, niter, True
    return vrs, options.max_iter, False


def extract_autocorr(vr: vector2) -> vector2:
    """Extract the quadratic function where its roots are within a unit circle

    x^2 + r*x + t  or x^2 + (r/t) * x + (1/t)

    (x + a1)(x + a2) = x^2 + (a1 + a2) x + a1 * a2

    Args:
        vr (vector2): [description]

    Returns:
        vector2: [description]
    """
    r, t = vr.x, vr.y
    hr = r / 2.0
    d = hr * hr - t
    if d < 0.0:  # complex conjugate root
        if t > 1.0:
            vr = vector2(r, 1.0) / t
    else:
        # two real roots
        a1 = hr + (sqrt(d) if hr >= 0.0 else -sqrt(d))
        a2 = t / a1
        if abs(a1) > 1.0:
            if abs(a2) > 1.0:
                a2 = 1.0 / a2
            a1 = 1.0 / a1
            vr = vector2(a1 + a2, a1 * a2)
        elif abs(a2) > 1.0:
            a2 = 1.0 / a2
            vr = vector2(a1 + a2, a1 * a2)
    return vr
