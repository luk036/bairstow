from math import acos, cos, sqrt, pi
from typing import List

from .lds import VdCorput
from .robin import Robin
from .rootfinding import Options, delta, horner, suppress_old
from .vector2 import Vector2

PI = pi


def initial_autocorr_new(coeffs: List[float]) -> List[Vector2]:
    """[summary]

    Args:
        coeffs (List[float]): [description]

    Returns:
        List[Vector2]: [description]

    Examples:
        >>> h = [10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0]
        >>> vr0s = initial_autocorr(h)
    """
    degree = len(coeffs) - 1
    re = pow(abs(coeffs[-1]), 1.0 / degree)
    if re > 1:
        re = 1 / re
    degree //= 2
    # k = PI / degree
    m = re * re
    # vr0s = [Vector2(2 * re * cos(k * i), -m) for i in range(1, degree, 2)]
    vgen = VdCorput(2)
    vgen.reseed(1)
    vr0s = [Vector2(2 * re * cos(PI * vgen.pop()), -m) for _ in range(1, degree, 2)]
    return vr0s


def initial_autocorr(coeffs: List[float]) -> List[Vector2]:
    """[summary]

    Args:
        coeffs (List[float]): [description]

    Returns:
        List[Vector2]: [description]
    """
    degree = len(coeffs) - 1
    re = pow(abs(coeffs[-1]), 1.0 / degree)
    if re < 1:  # use those outside the unit circle
        re = 1 / re
    degree //= 2
    k = PI / degree
    m = re * re
    vr0s = [Vector2(2 * re * cos(k * i), -m) for i in range(1, degree, 2)]
    return vr0s


def pbairstow_autocorr(
    coeffs: List[float], vrs: List[Vector2], options: Options = Options()
):
    """[summary]

    Args:
        coeffs (List[float]): [description]
        vrs (List[Vector2]): [description]
        options (Options, optional): [description]. Defaults to Options().

    Returns:
        [type]: [description]

    Examples:
        >>> h = [10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0]
        >>> vr0s = initial_autocorr(h)
        >>> vrs, niter, found = pbairstow_autocorr(h, vr0s)
    """
    M = len(vrs)  # assume polynomial of h is even
    degree = len(coeffs) - 1
    converged = [False] * M
    robin = Robin(M)
    for niter in range(options.max_iters):
        tol = 0.0
        # found = True  # initial
        for i in filter(lambda i: converged[i] is False, range(M)):
            pb = coeffs.copy()
            vA = horner(pb, degree, vrs[i])
            tol_i = max(abs(vA.x), abs(vA.y))
            if tol_i < options.tol_ind:
                converged[i] = True
                continue
            tol = max(tol, tol_i)
            vA1 = horner(pb, degree - 2, vrs[i])
            # for j in filter(lambda j: j != i, range(M)):  # exclude i
            # for j in robin.exclude(i):
            #     suppress_old(vA, vA1, vrs[i], vrs[j])
            #     # vA1 -= delta1(vA, vrs[j], vrs[i] - vrs[j])
            # vrs[i] -= delta2(vA, vrs[i], vA1)
            for j in robin.exclude(i):
                suppress_old(vA, vA1, vrs[i], vrs[j])
                # for j in range(M):
                vrn = Vector2(-vrs[j].x, 1.0) / vrs[j].y
                suppress_old(vA, vA1, vrs[i], vrn)
            vrs[i] -= delta(vA, vrs[i], vA1)
            # vrs[i] = extract_autocorr(vrs[i])
        # if vrs[i].y > 1.0:
        #     vrs[i] = Vector2(vrs[i].x, 1.0) / vrs[i].y
        # for i in range(M):  # exclude converged
        if tol < options.tol:
            return vrs, niter, True
    return vrs, options.max_iters, False


def extract_autocorr(vr: Vector2) -> Vector2:
    """Extract the quadratic function where its roots are within a unit circle

    x^2 - r*x - q  or (-1/q) + (r/q) * x + x^2
    (x - a1)(x - a2) = x^2 - (a1 + a2) x + a1 * a2

    Args:
        vr (Vector2): [description]

    Returns:
        Vector2: [description]

    Examples:
        >>> vr = extract_autocorr(Vector2(1, -4))
        >>> print(vr)
        <0.25, -0.25>
    """
    r, q = vr.x, vr.y
    hr = r / 2.0
    d = hr * hr + q
    if d < 0.0:  # complex conjugate root
        if q < -1.0:
            vr = Vector2(-r, 1.0) / q
    else:
        # two real roots
        a1 = hr + (sqrt(d) if hr >= 0.0 else -sqrt(d))
        a2 = -q / a1
        if abs(a1) > 1.0:
            if abs(a2) > 1.0:
                a2 = 1.0 / a2
            a1 = 1.0 / a1
            vr = Vector2(a1 + a2, -a1 * a2)
        elif abs(a2) > 1.0:
            a2 = 1.0 / a2
            vr = Vector2(a1 + a2, -a1 * a2)
    return vr
