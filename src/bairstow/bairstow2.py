from .vector2 import Vector2
from .matrix2 import Matrix2
from .rootfinding import Options
from .robin import Robin
from typing import List, Tuple


def horner2(coeffs: List, degree: int, vr: Vector2) -> Vector2:
    """[summary]

                       ⎛ 2            ⎞
        P(x) = P (x) ⋅ ⎝x  - r ⋅ x - q⎠ + A ⋅ x + B
                1

    Note: pb becomes the quotient after calling this function

    Args:
        coeffs (List[float]): [description]
        vr (Vector2): [description]

    Returns:
        Vector2: [description]

    Examples:
        >>> coeffs = [1, -8, -72, 382, 727, -2310]
        >>> vp = horner2(coeffs, 5, Vector2(2, 3))  # x^2 - 2 x - 3
        >>> coeffs
        [1, -6, -81, 202, 888, 72]
        >>> vp = horner2(coeffs, 3, Vector2(2, 3))
        >>> coeffs[:4]
        [1, -4, -86, 18]
    """
    b0 = 0
    b1 = coeffs[0]
    for i in range(degree):
        coeffs[i + 1] += vr.x * b1 + vr.y * b0
        b0 = b1
        b1 = coeffs[i + 1]
    return Vector2(b0, b1)


def bairstow2(
    pa: List[float], vr: Vector2, options=Options()
) -> Tuple[Vector2, int, bool]:
    """Original Bairstow's method

    Args:
        pa (List[float]): [description]
        vr (Vector2): [description]
        options (Options, optional): [description]. Defaults to Options().

    Returns:
        [type]: [description]

    Examples:
        >>> coeffs = [1, -8, -72, 382, 727, -2310]
        >>> vr0 = Vector2(2, 3)
        >>> vr, niter, found = bairstow2(coeffs, vr0)
        >>> found
        True
        >>> niter
        6
        >>> print(vr)
        <2.0, 15.0>
        >>> coeffs = [1, -6, -69, 154]
        >>> vr0 = Vector2(2, 3)
        >>> vr, niter, found = bairstow2(coeffs, vr0)
        >>> found
        True
        >>> niter
        7
        >>> print(vr)
        <4.0, 77.0>
    """
    n = len(pa) - 1
    for niter in range(options.max_iters):
        pb = pa.copy()
        vb = horner2(pb, n, vr)
        tol = max(abs(vb.x), abs(vb.y))
        if tol < options.tol:
            return vr, niter, True
        vc = horner2(pb, n - 2, vr)
        cb = vb.x + vr.x * vc.y + vr.y * vc.x
        mat_c = Matrix2(Vector2(vc.y, -vc.x), Vector2(-cb, vc.y))
        vr -= mat_c.mdot(vb) / mat_c.det()
    return vr, options.max_iters, False


def suppress(vb: Vector2, vc: Vector2, vri: Vector2, vrj: Vector2):
    """[summary]

    Args:
        vA (Vector2): [description]
        vr (Vector2): [description]
        vp (Vector2): [description]

    Returns:
        Vector2: [description]

    Examples:
        >>> vb = Vector2(3, 3)
        >>> vc = Vector2(1, 2)
        >>> vri = Vector2(-2, 0)
        >>> vrj = Vector2(4, 5)
        >>> vb, vc = suppress(vb, vc, vri, vrj)
    """
    vp = vri - vrj
    r, q = vri.x, vri.y
    p, s = vp.x, vp.y
    m_adjoint = Matrix2(Vector2(r * p + s, -p), Vector2(-q * p, s))
    e = m_adjoint.det()
    vbb = m_adjoint.mdot(vb)
    vt = vc * e - vbb
    vt._y -= vbb._x * p
    vbb *= e
    vcc = m_adjoint.mdot(vt)
    return vbb, vcc


def pbairstow2_even(
    pa: List[float], vrs: List[Vector2], options=Options()
) -> Tuple[List[Vector2], int, bool]:
    """Parallel Bairstow's method

    Args:
        pa (List[float]): [description]
        vrs (List[Vector2]): [description]
        options (Options, optional): [description]. Defaults to Options().

    Returns:
        [type]: [description]
    """
    M = len(vrs)
    N = len(pa) - 1
    converged = [False] * M
    robin = Robin(M)
    for niter in range(options.max_iters):
        tol = 0.0
        # exclude converged
        for i in filter(lambda i: converged[i] is False, range(M)):
            # for i in range(M):
            pb = pa.copy()
            vb = horner2(pb, N, vrs[i])
            tol_i = max(abs(vb.x), abs(vb.y))
            if tol_i < options.tol_ind:
                converged[i] = True
                continue
            vc = horner2(pb, N - 2, vrs[i])
            tol = max(tol_i, tol)
            for j in robin.exclude(i):
                vb, vc = suppress(vb, vc, vrs[i], vrs[j])

            cb = vb.x + vrs[i].x * vc.y + vrs[i].y * vc.x
            mat_c = Matrix2(Vector2(vc.y, -vc.x), Vector2(-cb, vc.y))
            vrs[i] -= mat_c.mdot(vb) / mat_c.det()
        if tol < options.tol:
            return vrs, niter, True
    return vrs, options.max_iters, False
