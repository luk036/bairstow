from cmath import exp
from math import pi
from typing import List
# from pytest import approx
from .rootfinding import Options
from .lds import Vdcorput


# PI = acos(-1.0)
PI = pi

#     """[summary]
#
#     Args:
#         pa (List[float]): [description]
#         r (float): [description]
#
#     Returns:
#         float: [description]
#     """
#     ans = pb[0]
#     for i in range(1, len(pb)):
#         ans = ans * z + pb[i]
#     return ans


def horner_eval(pb: List[float], n: int, alpha: float) -> float:
    """[summary]

    Args:
        pb (List[float]): [description]
        n (int): [description]
        alpha (float): [description]

    Returns:
        float: [description]

    Examples:
        >>> p = [1.0, -6.7980, 2.9948, -0.043686, 0.000089248]
        >>> n = len(p) - 1
        >>> alpha = 6.3256
        >>> P = horner_eval(p, n, alpha)
        >>> P
        -0.012701469838522064
        >>> p[3]
        -0.0020220560640132265
    """
    for i in range(n):
        pb[i + 1] += pb[i] * alpha
    return pb[n]


def horner_backward(pb: List[float], n: int, alpha: float) -> float:
    """[summary]

    Args:
        pa (List[float]): [description]
        r (float): [description]

    Returns:
        float: [description]

    Examples:
        >>> p = [1.0, -6.7980, 2.9948, -0.043686, 0.000089248]
        >>> n = len(p) - 1
        >>> alpha = 6.3256
        >>> P = horner_backward(p, n, alpha)
        >>> -P * (alpha ** 5)
        -0.013355264987140483
        >>> p[3]
        0.006920331351966613
    """
    for i in range(2, n + 2):
        pb[-i] -= pb[-(i - 1)]
        pb[-i] /= -alpha
    return pb[-(n + 1)]


# def initial_aberth_lds(pa: List) -> List:
#     """[summary]

#     Args:
#         pa (List): [description]

#     Returns:
#         List: [description]
#     """
#     N = len(pa) - 1
#     c = -pa[1] / (N * pa[0])
#     Pc = horner_eval(pa.copy(), N, c)
#     re = (-Pc) ** (1.0 / N)
#     z0s = []
#     two_PI = 2 * PI
#     vdc_gen = vdcorput()

#     for i in range(N):
#         theta = two_PI * vdc_gen() + 0.25
#         z0s += [c + re * exp(theta * 1j)]
#     return z0s


def initial_aberth(pa: List) -> List:
    """[summary]

    Args:
        pa (List): [description]

    Returns:
        List: [description]

    Examples:
        >>> h = [5.0, 2.0, 9.0, 6.0, 2.0]
        >>> z0s = initial_aberth(h)
    """
    N = len(pa) - 1
    c = -pa[1] / (N * pa[0])
    Pc = horner_eval(pa.copy(), N, c)
    re = (-Pc) ** (1.0 / N)
    # k = 2 * PI / N
    z0s = []
    vgen = Vdcorput(2)
    vgen.reseed(1)
    for i in range(N):
        vdc = 2 * PI * vgen.pop()
        z0s += [c + re * exp(vdc * 1j)]
    return z0s


def aberth(pa: List, zs: List, options: Options = Options()):
    """[summary]

    Args:
        pa (List): [description]
        zs (List): [description]
        options (Options, optional): [description]. Defaults to Options().

    Returns:
        [type]: [description]

    Examples:
        >>> h = [5.0, 2.0, 9.0, 6.0, 2.0]
        >>> z0s = initial_aberth(h)
        >>> opt = Options()
        >>> opt.tol = 1e-8
        >>> zs, niter, found = aberth(h, z0s, opt)
    """
    M = len(zs)
    N = len(pa) - 1
    converged = [False] * M
    for niter in range(1, options.max_iter):
        tol = 0.0
        for i in filter(lambda i: not converged[i], range(M)):  # exclude converged
            pb = pa.copy()
            P = horner_eval(pb, N, zs[i])
            tol_i = abs(P)
            if tol_i < options.tol_ind:
                converged[i] = True
                continue
            P1 = horner_eval(pb, N - 1, zs[i])
            tol = max(tol_i, tol)
            for j in filter(lambda j: j != i, range(M)):  # exclude i
                P1 -= P / (zs[i] - zs[j])
            zs[i] -= P / P1
        if tol < options.tol:
            return zs, niter, True
    return zs, options.max_iter, False


def initial_aberth_autocorr(pa: List) -> List:
    """[summary]

    Args:
        pa (List): [description]

    Returns:
        List: [description]

    Examples:
        >>> h = [5.0, 2.0, 9.0, 6.0, 2.0]
        >>> z0s = initial_aberth_autocorr(h)
    """
    N = len(pa) - 1
    re = abs(pa[-1]) ** (1.0 / N)
    # c = -pa[1] / (N * pa[0])
    # Pc = horner_eval(pa.copy(), N, c)
    # re = (-Pc) ** (1.0 / N)
    if abs(re) > 1:
        re = 1 / re
    N //= 2
    # k = 2 * PI / N
    z0s = []
    vgen = Vdcorput(2)
    vgen.reseed(1)
    for i in range(N):
        vdc = 2 * PI * vgen.pop()
        z0s += [re * exp(vdc * 1j)]
    return z0s


def aberth_autocorr(pa: List, zs: List, options: Options = Options()):
    """[summary]

    Args:
        pa (List): [description]
        zs (List): [description]
        options (Options, optional): [description]. Defaults to Options().

    Returns:
        [type]: [description]

    Examples:
        >>> h = [5.0, 2.0, 9.0, 6.0, 2.0]
        >>> z0s = initial_aberth_autocorr(h)
        >>> zs, niter, found = aberth_autocorr(h, z0s)
        >>> zs[0]
        (-0.35350437336258744+0.3130287231135712j)
        >>> opt = Options()
        >>> opt.tol = 1e-8
        >>> zs, niter, found = aberth_autocorr(h, z0s, opt)
        >>> zs[0]
        (-0.35350437336258744+0.3130287231135712j)
    """
    M = len(zs)
    N = len(pa) - 1
    converged = [False] * M
    for niter in range(1, options.max_iter):
        tol = 0.0
        for i in filter(lambda i: not converged[i], range(M)):  # exclude converged
            pb = pa.copy()
            P = horner_eval(pb, N, zs[i])
            tol_i = abs(P)
            if tol_i < options.tol_ind:
                converged[i] = True
                continue
            P1 = horner_eval(pb, N - 1, zs[i])
            tol = max(tol_i, tol)
            for j in filter(lambda j: j != i, range(M)):  # exclude i
                P1 -= P / (zs[i] - zs[j])
            for j in range(M):  # exclude i
                zsn = 1.0 / zs[j]
                P1 -= P / (zs[i] - zsn)
            zs[i] -= P / P1
            if abs(zs[i]) > 1.0:  # pick those inside the unit circle
                zs[i] = 1.0 / zs[i]
        if tol < options.tol:
            return zs, niter, True
    return zs, options.max_iter, False
