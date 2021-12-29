from cmath import exp
from math import acos
from typing import List

PI = acos(-1.0)


def horner_eval(pb: List[float], n: int, z: float) -> float:
    """[summary]

    Args:
        pa (List[float]): [description]
        r (float): [description]

    Returns:
        float: [description]
    """
    for i in range(n):
        pb[i + 1] += pb[i] * z
    return pb[n]


class Options:
    max_iter: int = 2000
    tol: float = 1e-12
    tol_ind: float = 1e-15


def initial_aberth(pa: List) -> List:
    """[summary]

    Args:
        pa (List): [description]

    Returns:
        List: [description]
    """
    N = len(pa) - 1
    c = -pa[1] / (N * pa[0])
    Pc = horner_eval(pa.copy(), N, c)
    re = (-Pc) ** (1.0 / N)
    z0s = []
    k = 2 * PI / N
    for i in range(N):
        z0s += [c + re * exp(k * (i + 0.25) * 1j)]
    return z0s


def paberth(pa: List, zs: List, options: Options = Options()):
    """[summary]

    Args:
        pa (List): [description]
        zs (List): [description]
        options (Options, optional): [description]. Defaults to Options().

    Returns:
        [type]: [description]
    """
    M = len(zs)
    N = len(pa) - 1
    found = False
    converged = [False] * M
    for niter in range(options.max_iter):
        tol = 0
        for i in filter(lambda i: converged[i] is False, range(M)):  # exclude converged
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
            found = True
            break
    return zs, niter + 1, found
