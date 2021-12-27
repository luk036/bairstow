from math import acos, cos, pow, sqrt
from typing import List, Tuple

from .matrix2 import matrix2
from .vector2 import vector2

PI = acos(-1.0)


def delta(vA: vector2, vr: vector2, vp: vector2) -> vector2:
    """[summary]

    r * p - m   -p
    q * p       -m

    Args:
        vA (vector2): [description]
        vr (vector2): [description]
        vp (vector2): [description]

    Returns:
        vector2: [description]
    """
    r, q = vr.x, vr.y
    p, m = vp.x, vp.y
    mp = matrix2(vector2(-m, p), vector2(-p * q, p * r - m))
    return mp.mdot(vA) / mp.det()  # 6 mul's + 2 div's


def horner_eval(pa: List[float], r: float) -> float:
    """[summary]

    Args:
        pa (List[float]): [description]
        r (float): [description]

    Returns:
        float: [description]
    """
    pb = pa.copy()
    for i in range(len(pa)):
        pb[i] += pb[i - 1] * r
    return pb[-1]


def horner(pa: List[float], vr: vector2) -> Tuple[vector2, List[float]]:
    """[summary]

    Args:
        pa (List[float]): [description]
        vr (vector2): [description]

    Returns:
        vector2: [description]
    """
    r, q = vr.x, vr.y
    n = len(pa) - 1
    pb = pa.copy()
    pb[1] -= pb[0] * r
    for i in range(2, n):
        pb[i] -= pb[i - 1] * r + pb[i - 2] * q
    pb[n] -= pb[n - 2] * q
    return vector2(pb[n - 1], pb[n]), pb[:-2]


class Options:
    max_iter: int = 2000
    tol: float = 1e-12


def initial_guess(pa: List[float]) -> List[vector2]:
    """[summary]

    Args:
        pa (List[float]): [description]

    Returns:
        List[vector2]: [description]
    """
    N = len(pa) - 1
    c = -pa[1] / (N * pa[0])
    # P = np.poly1d(pa)
    Pc = horner_eval(pa, c)
    re = pow(abs(Pc), 1.0 / N)
    m = c * c + re * re
    vr0s = []
    N //= 2
    N *= 2  # make even
    k = PI / N
    for i in range(1, N, 2):
        temp = re * cos(k * i)
        r0 = -2 * (c + temp)
        t0 = m + 2 * c * temp
        vr0s += [vector2(r0, t0)]
    return vr0s


def pbairstow_even(pa: List[float], vrs: List[vector2], options: Options = Options()):
    """[summary]

    Args:
        pa (List[float]): [description]
        vrs (List[vector2]): [description]
        options (Options, optional): [description]. Defaults to Options().

    Returns:
        [type]: [description]
    """
    M = len(vrs)
    found = False
    converged = [False] * M
    for niter in range(options.max_iter):
        tol = 0
        # found = True
        for i in filter(lambda i: converged[i] is False, range(M)):  # exclude converged
            # for i in range(M):
            vA, pb = horner(pa, vrs[i])
            tol_i = max(abs(vA.x), abs(vA.y))
            if tol_i < 1e-15:
                converged[i] = True
                continue
            # found = False
            vA1, _ = horner(pb, vrs[i])
            tol = max(tol_i, tol)
            for j in filter(lambda j: j != i, range(M)):  # exclude i
                vA1 -= delta(vA, vrs[j], vrs[i] - vrs[j])
            vrs[i] -= delta(vA, vrs[i], vA1)
        if tol < options.tol:
            found = True
            break
    return vrs, niter + 1, found


def find_rootq(vr: vector2) -> Tuple[float, float]:
    """[summary]

    x^2 + r*x + t

    (x - x1)(x - x2) = x^2 - (x1 + x2) x + x1 * x2

    Args:
        vr (vector2): [description]

    Returns:
        Tuple[float, float]: [description]
    """
    r, t = vr.x, vr.y
    hr = r / 2.0
    d = hr * hr - t
    if d < 0.0:
        x1 = -hr + sqrt(-d) * 1j
    else:
        x1 = -hr + (sqrt(d) if hr <= 0.0 else -sqrt(d))
    x2 = t / x1
    return x1, x2
