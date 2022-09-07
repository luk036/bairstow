from math import acos, cos, sqrt
from typing import List

from .rootfinding import Options
from .rootfinding import delta, horner
from .vector2 import vector2
from .lds import Vdcorput

PI = acos(-1.0)


def initial_autocorr(pa: List[float]) -> List[vector2]:
    """[summary]

    Args:
        pa (List[float]): [description]

    Returns:
        List[vector2]: [description]

    Examples:
        >>> h = [10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0]
        >>> vr0s = initial_autocorr(h)
    """
    N = len(pa) - 1
    re = abs(pa[-1]) ** (1.0 / N)
    if re > 1:
        re = 1 / re
    N //= 2
    # k = PI / N
    m = re * re
    # vr0s = [vector2(-2 * re * cos(k * i), m) for i in range(1, N, 2)]
    vgen = Vdcorput(2)
    vgen.reseed(1)
    vr0s = [vector2(-2 * re * cos(PI * vgen.pop()), m) for _ in range(1, N, 2)]
    return vr0s


def initial_autocorr_bad(pa: List[float]) -> List[vector2]:
    """[summary]

    Args:
        pa (List[float]): [description]

    Returns:
        List[vector2]: [description]
    """
    N = len(pa) - 1
    re = abs(pa[-1]) ** (1.0 / N)
    if re < 1:  # use those outside the unit circle
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

    Examples:
        >>> h = [10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0]
        >>> vr0s = initial_autocorr(h)
        >>> vrs, niter, found = pbairstow_autocorr(h, vr0s)
        >>> print(vrs[0])
        <-0.1711207835281031, 0.5573808087014712>
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
            vA1 = horner(pb, N - 2, vrs[i])
            for j in filter(lambda j: j != i, range(M)):  # exclude i
                vA1 -= delta(vA, vrs[j], vrs[i] - vrs[j])
            for j in range(M):
                vrn = vector2(vrs[j].x, 1.0) / vrs[j].y
                vA1 -= delta(vA, vrn, vrs[i] - vrn)
            vrs[i] -= delta(vA, vrs[i], vA1)
            vrs[i] = extract_autocorr(vrs[i])
        # if vrs[i].y > 1.0:
        #     vrs[i] = vector2(vrs[i].x, 1.0) / vrs[i].y
        # for i in range(M):  # exclude converged
        if tol < options.tol:
            return vrs, niter, True
    return vrs, options.max_iter, False


def pbairstow_autocorr_bad(
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
        for i in filter(lambda i: converged[i] is False, range(M)):  # exclude converged
            pb = pa.copy()
            vA = horner(pb, N, vrs[i])
            tol_i = max(abs(vA.x), abs(vA.y))
            if tol_i < options.tol_ind:
                converged[i] = True
                continue
            tol = max(tol, tol_i)
            vA1 = horner(pb, N - 2, vrs[i])
            for j in filter(lambda j: j != i, range(M)):  # exclude i
                vA1 -= delta(vA, vrs[j], vrs[i] - vrs[j])
            for j in range(M):
                vrn = vector2(vrs[j].x, 1.0) / vrs[j].y
                vA1 -= delta(vA, vrn, vrs[i] - vrn)
            vrs[i] -= delta(vA, vrs[i], vA1)
        # if vrs[i].y > 1.0:
        #     vrs[i] = vector2(vrs[i].x, 1.0) / vrs[i].y
        for i in range(M):  # exclude converged
            vrs[i] = extract_autocorr(vrs[i])
        if tol < options.tol:
            return vrs, niter, True
    return vrs, options.max_iter, False


def extract_autocorr(vr: vector2) -> vector2:
    """Extract the quadratic function where its roots are within a unit circle

    x^2 + r*x + t  or (1/t) + (r/t) * x + x^2
    (x + a1)(x + a2) = x^2 + (a1 + a2) x + a1 * a2

    Args:
        vr (vector2): [description]

    Returns:
        vector2: [description]

    Examples:
        >>> vr = extract_autocorr(vector2(-1, 4)) 
        >>> print(vr)
        <-0.25, 0.25>
    """
    r, t = vr.x, vr.y
    hr = r / 2.0
    d = hr * hr - t
    if d < 0.0:  # complex conjugate root
        if t > 1.0:
            vr = vector2(r / t, 1.0 / t)
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
