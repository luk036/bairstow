from math import acos, cos, pow, sqrt

from .matrix2 import matrix2
from .vector2 import vector2

PI = acos(-1.0)


def makeadjoint(vr, vp):
    """[summary]

    Args:
        vr ([type]): [description]
        vp ([type]): [description]

    Returns:
        [type]: [description]
    """
    r, q = vr.x, vr.y
    p, s = vp.x, vp.y
    return matrix2(vector2(s, -p), vector2(-p * q, p * r + s))


def suppress_orig(vA1, vA, vr, vrj):
    """[summary]

    Args:
        vA ([type]): [description]
        vA1 ([type]): [description]
        vr ([type]): [description]
        vrj ([type]): [description]

    Returns:
        [type]: [description]
    """
    vp = vr - vrj
    p = vp.x
    Ap = makeadjoint(vr, vp)  # 2 mul's
    e = Ap.det()  # 2 mul's
    va = Ap.dot(vA)  # 4 mul's
    vA = va  # 2 mul's
    a, b = va.x, va.y
    vc = vA1 - vector2(a, a * p + b) / e  # 3 mul
    vA1 = Ap.dot(vc)  # 4 mul
    return vA, vA1


def suppress(vA1, vA, vr, vrj):
    """[summary]

    Args:
        vA ([type]): [description]
        vA1 ([type]): [description]
        vr ([type]): [description]
        vrj ([type]): [description]

    Returns:
        [type]: [description]
    """
    vp = vr - vrj
    mp = makeadjoint(vrj, vp)  # 2 mul's
    vA1 -= mp.mdot(vA) / mp.det()  # 6 mul's + 2 div's
    return vA1


def delta(vA, vr, vp):
    """[summary]

    Args:
        vA ([type]): [description]
        vr ([type]): [description]
        vrj ([type]): [description]

    Returns:
        [type]: [description]
    """
    mp = makeadjoint(vr, vp)  # 2 mul's
    return mp.mdot(vA) / mp.det()  # 6 mul's + 2 div's


def check_newton(vA, vA1, vr):
    """[summary]

    Args:
        vA ([type]): [description]
        vA1 ([type]): [description]
        vr ([type]): [description]

    Returns:
        [type]: [description]
    """
    mA1 = makeadjoint(vr, vA1)  # 2 mul's
    return mA1.mdot(vA) / mA1.det()  # 6 mul's + 2 div's


def horner_eval(pa, r):
    """[summary]

    Args:
        pa ([type]): [description]
        r ([type]): [description]

    Returns:
        [type]: [description]
    """
    pb = pa.copy()
    for i in range(len(pa)):
        pb[i] += pb[i - 1] * r
    return pb[-1]


def horner(pa, vr):
    """[summary]

    Args:
        pa ([type]): [description]
        vr ([type]): [description]

    Returns:
        [type]: [description]
    """
    r, q = vr.x, vr.y
    n = len(pa) - 1
    pb = pa.copy()
    pb[1] += pb[0] * r
    for i in range(2, n):
        pb[i] += pb[i - 2] * q + pb[i - 1] * r
    pb[n] += pb[n - 2] * q
    return vector2(pb[n - 1], pb[n]), pb[:-2]


class Options:
    max_iter: int = 2000
    tol: float = 1e-12


def initial_guess(pa):
    """[summary]

    Args:
        pa ([type]): [description]

    Returns:
        [type]: [description]
    """
    N = len(pa) - 1
    c = -pa[1] / (N * pa[0])
    # P = np.poly1d(pa)
    Pc = horner_eval(pa, c)
    re = pow(abs(Pc), 1.0 / N)
    k = PI / N
    m = c * c + re * re
    vr0s = []
    for i in range(1, N, 2):
        temp = re * cos(k * i)
        r0 = 2 * (c + temp)
        t0 = m + 2 * c * temp
        vr0s += [vector2(r0, -t0)]
    return vr0s


def pbairstow_even(pa, vrs, options=Options()):
    """[summary]

    Args:
        pa ([type]): [description]
        vrs ([type]): [description]
        options ([type], optional): [description]. Defaults to Options().

    Returns:
        [type]: [description]
    """
    M = len(vrs)
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
            vrs[i] -= delta(vA, vrs[i], vA1)
            
        if tol < options.tol:
            found = True
            break
    return vrs, niter + 1, found


def find_rootq(r, t):
    """[summary]

    x^2 - r*x + t  or x^2 - (r/t) * x + (1/t)

    (x - x1)(x - x2) = x^2 - (x1 + x2) x + x1 * x2

    determinant r/2 + q

    Args:
        r ([type]): [description]
        t ([type]): [description]

    Returns:
        [type]: [description]
    """
    hr = r / 2.0
    d = hr * hr - t
    if d < 0.0:
        x1 = hr + sqrt(-d) * 1j
    else:
        x1 = hr + (sqrt(d) if hr >= 0.0 else -sqrt(d))
    x2 = t / x1
    return x1, x2
