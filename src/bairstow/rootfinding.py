from math import acos, cos, sqrt
from typing import List, Tuple

from .matrix2 import Matrix2
from .vector2 import Vector2

PI = acos(-1.0)

class Options:
    max_iter: int = 2000
    tol: float = 1e-12
    tol_ind: float = 1e-15
    # tol_suppress: float = 1e-1


# def horner_eval(pb: List[float], z):

def delta(vA: Vector2, vr: Vector2, vp: Vector2) -> Vector2:
    """[summary]

    r * p - m   -p
    q * p       -m

    Args:
        vA (Vector2): [description]
        vr (Vector2): [description]
        vp (Vector2): [description]

    Returns:
        Vector2: [description]

    Examples:
        >>> d = delta(Vector2(1, 2), Vector2(2, 0), Vector2(4, 5))
        >>> print(d)
        <-0.2, -0.4>
    """
    r, q = vr.x, vr.y
    p, m = vp.x, vp.y
    mp = Matrix2(Vector2(-m, p), Vector2(-p * q, p * r - m))
    return mp.mdot(vA) / mp.det()  # 6 mul's + 2 div's


def horner_eval(pb: List[float], n: int, z: float) -> float:
    """[summary]

    Args:
        pa (List[float]): [description]
        r (float): [description]

    Returns:
        float: [description]

    Examples:
        >>> p = [10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0]
        >>> n = len(p) - 1
        >>> P = horner_eval(p, n, 2.0)
        >>> P
        18250.0
        >>> p[3]
        460.0
    """
    for i in range(n):
        pb[i + 1] += pb[i] * z
    return pb[n]


def horner(pb: List[float], n: int, vr: Vector2) -> Vector2:
    """[summary]

    Args:
        pa (List[float]): [description]
        vr (Vector2): [description]

    Returns:
        Vector2: [description]

    Examples:
        >>> p = [10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0]
        >>> n = len(p) - 1
        >>> P = horner(p, n, Vector2(-1.0, -2.0))
        >>> P.x * 2.0 + P.y
        18250.0
    """
    r, q = vr.x, vr.y
    pb[1] -= pb[0] * r
    for i in range(2, n):
        pb[i] -= pb[i - 1] * r + pb[i - 2] * q
    pb[n] -= pb[n - 2] * q
    return Vector2(pb[n - 1], pb[n])


def initial_guess(pa: List[float]) -> List[Vector2]:
    """[summary]

    Args:
        pa (List[float]): [description]

    Returns:
        List[Vector2]: [description]

    Examples:
        >>> h = [10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0]
        >>> vr0s = initial_guess(h)
        >>> print(vr0s[0])
        <-1.4537520283010816, 0.7559947795353074>
    """
    N = len(pa) - 1
    c = -pa[1] / (N * pa[0])
    # P = np.poly1d(pa)
    Pc = horner_eval(pa.copy(), N, c)
    reff = abs(Pc) ** (1.0 / N)
    m = c * c + reff * reff
    vr0s = []
    N //= 2
    N *= 2  # make even
    k = PI / N
    for i in range(1, N, 2):
        temp = reff * cos(k * i)
        r0 = -2 * (c + temp)
        t0 = m + 2 * c * temp
        vr0s += [Vector2(r0, t0)]
    return vr0s


def pbairstow_even(pa: List[float], vrs: List[Vector2], options = Options()) \
-> Tuple[List[Vector2], int, bool]:
    """[summary]

    Args:
        pa (List[float]): [description]
        vrs (List[Vector2]): [description]
        options (Options, optional): [description]. Defaults to Options().

    Returns:
        [type]: [description]

    Examples:
        >>> h = [10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0]
        >>> vr0s = initial_guess(h)
        >>> vrs, niter, found = pbairstow_even(h, vr0s)
        >>> print(vrs[0])
        <-0.1711207835281031, 0.5573808087014712>
    """
    M = len(vrs)
    N = len(pa) - 1
    converged = [False] * M
    for niter in range(1, options.max_iter):
        tol = 0.0
        for i in filter(lambda i: converged[i] is False, range(M)):  # exclude converged
            # for i in range(M):
            pb = pa.copy()
            vA = horner(pb, N, vrs[i])
            tol_i = max(abs(vA.x), abs(vA.y))
            if tol_i < options.tol_ind:
                converged[i] = True
                continue
            vA1 = horner(pb, N - 2, vrs[i])
            tol = max(tol_i, tol)
            for j in filter(lambda j: j != i, range(M)):  # exclude i
                vA1 -= delta(vA, vrs[j], vrs[i] - vrs[j])
            vrs[i] -= delta(vA, vrs[i], vA1)
        if tol < options.tol:
            return vrs, niter, True
    return vrs, options.max_iter, False


def find_rootq(vr: Vector2) -> Tuple[float, float]:
    """Solve x^2 + r*x + t = 0

    (x - x1)(x - x2) = x^2 - (x1 + x2) x + x1 * x2

    Args:
        vr (Vector2): [description]

    Returns:
        Tuple[float, float]: [description]

    Examples:
        >>> vr = find_rootq(Vector2(-5, 6)) 
        >>> print(vr)
        (3.0, 2.0)
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
