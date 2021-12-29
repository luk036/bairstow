from cmath import exp
from math import acos
from typing import List

# from pylds.low_discr_seq import vdcorput

PI = acos(-1.0)


class Options:
    max_iter: int = 2000
    tol: float = 1e-12
    tol_ind: float = 1e-15


def horner_eval(pb: List[float], z):
    """[summary]

    Args:
        pa (List[float]): [description]
        r (float): [description]

    Returns:
        float: [description]
    """
    ans = pb[0]
    for i in range(1, len(pb)):
        ans = ans * z + pb[i]
    return ans


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
    """
    N = len(pa) - 1
    c = -pa[1] / (N * pa[0])
    Pc = horner_eval(pa, c)
    re = (-Pc) ** (1.0 / N)
    k = 2 * PI / N
    z0s = []
    for i in range(N):
        z0s += [c + re * exp(k * (i + 0.25) * 1j)]
    return z0s


def aberth(pa: List, zs: List, options: Options = Options()):
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
    pb = [(N - i) * p for i, p in enumerate(pa[:-1])]
    for niter in range(options.max_iter):
        tol = 0
        for i in filter(lambda i: converged[i] is False, range(M)):  # exclude converged
            P = horner_eval(pa, zs[i])
            tol_i = abs(P)
            if tol_i < options.tol_ind:
                converged[i] = True
                continue
            P1 = horner_eval(pb, zs[i])
            tol = max(tol_i, tol)
            for j in filter(lambda j: j != i, range(M)):  # exclude i
                P1 -= P / (zs[i] - zs[j])
            zs[i] -= P / P1
        if tol < options.tol:
            found = True
            break
    return zs, niter + 1, found
