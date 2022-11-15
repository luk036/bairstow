from cmath import exp
from math import pi
from typing import List, Tuple, Union

from .lds import Vdcorput
from .robin import Robin

# from pytest import approx
from .rootfinding import Options, horner_eval

FoC = Union[float, complex]
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


def horner_backward(pb: List, n: int, alpha: FoC) -> FoC:
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


def initial_aberth(pa: List[FoC]) -> List[complex]:
    """[summary]

    Args:
        pa (List): [description]

    Returns:
        List: [description]

    Examples:
        >>> h = [5.0, 2.0, 9.0, 6.0, 2.0]
        >>> z0s = initial_aberth(h)
    """
    N: int = len(pa) - 1
    c: FoC = -pa[1] / (N * pa[0])
    Pc: FoC = horner_eval(pa.copy(), N, c)
    re: FoC = (-Pc) ** (1.0 / N)
    # k = 2 * PI / N
    z0s: List[complex] = []
    vgen = Vdcorput(2)
    vgen.reseed(1)
    for i in range(N):
        vdc = 2 * PI * vgen.pop()
        z0s += [c + re * exp(vdc * 1j)]
    return z0s


def initial_aberth_orig(pa: List[FoC]) -> List[complex]:
    """[summary]

    Args:
        pa (List): [description]

    Returns:
        List: [description]

    Examples:
        >>> h = [5.0, 2.0, 9.0, 6.0, 2.0]
        >>> z0s = initial_aberth_orig(h)
    """
    N: int = len(pa) - 1
    c: FoC = -pa[1] / (N * pa[0])
    Pc: FoC = horner_eval(pa.copy(), N, c)
    re: FoC = (-Pc) ** (1.0 / N)
    k = 2 * PI / N
    z0s: List[complex] = []
    for i in range(N):
        theta = k * (0.25 + i)
        z0s += [c + re * exp(theta * 1j)]
    return z0s


def aberth(
    pa: List[FoC], zs: List[complex], options=Options()
) -> Tuple[List[complex], int, bool]:
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
    robin = Robin(M)
    for niter in range(1, options.max_iter):
        tol = 0.0
        for i in filter(lambda i: not converged[i], range(M)):
            pb = pa.copy()
            P = horner_eval(pb, N, zs[i])
            tol_i = abs(P)
            if tol_i < options.tol_ind:
                converged[i] = True
                continue
            P1 = horner_eval(pb, N - 1, zs[i])
            tol = max(tol_i, tol)
            # for j in filter(lambda j: j != i, range(M)):  # exclude i
            for j in robin.exclude(i):
                P1 -= P / (zs[i] - zs[j])
            zs[i] -= P / P1
        if tol < options.tol:
            return zs, niter, True
    return zs, options.max_iter, False


def initial_aberth_autocorr(pa: List[float]) -> List[complex]:
    """[summary]

    Args:
        pa (List): [description]

    Returns:
        List: [description]

    Examples:
        >>> h = [5.0, 2.0, 9.0, 6.0, 2.0]
        >>> z0s = initial_aberth_autocorr(h)
    """
    N: int = len(pa) - 1
    re: float = abs(pa[-1]) ** (1.0 / N)
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


def initial_aberth_autocorr_orig(pa: List[float]) -> List[complex]:
    """[summary]

    Args:
        pa (List): [description]

    Returns:
        List: [description]

    Examples:
        >>> h = [5.0, 2.0, 9.0, 6.0, 2.0]
        >>> z0s = initial_aberth_autocorr_orig(h)
    """
    N: int = len(pa) - 1
    re: float = abs(pa[-1]) ** (1.0 / N)
    # c = -pa[1] / (N * pa[0])
    # Pc = horner_eval(pa.copy(), N, c)
    # re = (-Pc) ** (1.0 / N)
    if abs(re) > 1:
        re = 1 / re
    N //= 2
    k = 2 * PI / N
    z0s = []
    for i in range(N):
        theta = k * (0.25 + i)
        z0s += [re * exp(theta * 1j)]
    return z0s


def aberth_autocorr(
    pa: List[float], zs: List[complex], options=Options()
) -> Tuple[List[complex], int, bool]:
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
        >>> opt = Options()
        >>> opt.tol = 1e-8
        >>> zs, niter, found = aberth_autocorr(h, z0s, opt)
    """
    M: int = len(zs)
    N: int = len(pa) - 1
    converged: List[bool] = [False] * M
    robin = Robin(M)
    for niter in range(1, options.max_iter):
        tol: float = 0.0
        # exclude converged
        for i in filter(lambda i: not converged[i], range(M)):
            pb = pa.copy()
            P = horner_eval(pb, N, zs[i])
            tol_i = abs(P)
            if tol_i < options.tol_ind:
                converged[i] = True
                continue
            P1 = horner_eval(pb, N - 1, zs[i])
            tol = max(tol_i, tol)
            # for j in filter(lambda j: j != i, range(M)):  # exclude i
            for j in robin.exclude(i):
                P1 -= P / (zs[i] - zs[j])
                # for j in range(M):  # exclude i
                zsn = 1.0 / zs[j]
                P1 -= P / (zs[i] - zsn)
            zs[i] -= P / P1
            # if abs(zs[i]) > 1.0:  # pick those inside the unit circle
            #     zs[i] = 1.0 / zs[i]
        if tol < options.tol:
            return zs, niter, True
    return zs, options.max_iter, False
