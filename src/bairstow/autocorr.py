from math import acos, cos, pow, sqrt
from .rootfinding import makeadjoint, horner, Options, delta
from .vector2 import vector2

PI = acos(-1.0)


def initial_autocorr(pa):
    """[summary]

    Args:
        pa ([type]): [description]

    Returns:
        [type]: [description]
    """
    N = len(pa) - 1
    re = pow(abs(pa[-1]), 1.0 / N)
    N //= 2
    k = PI / N
    m = re * re
    vr0s = [vector2(2*re*cos(k*i), m) for i in range(1, N, 2)]
    return vr0s


def pbairstow_autocorr(pa, vrs, options=Options()):
    """[summary]

    Args:
        pa ([type]): [description]
        vrs ([type]): [description]
        options ([type], optional): [description]. Defaults to Options().

    Returns:
        [type]: [description]
    """
    M = len(vrs)  # assume polynomial of h is even
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
                # mp = makeadjoint(vrs[j], vrs[i] - vrs[j])  # 2 mul's
                # vA1 -= mp.mdot(vA) / mp.det()  # 6 mul's + 2 div's
                vA1 -= delta(vA, vrs[i], vrs[j])

            for j in range(M):
                vrn = vector2(-vrs[j].x, 1.0) / vrs[j].y
                # mp = makeadjoint(vrn, vrs[i] - vrn)  # 2 mul's
                # vA1 -= mp.mdot(vA) / mp.det()  # 6 mul's + 2 div's
                vA1 -= delta(vA, vrs[i], vrn)

            mA1 = makeadjoint(vrs[i], vA1)  # 2 mul's
            vrs[i] -= mA1.mdot(vA) / mA1.det()  # Gauss-Seidel fashion
        if tol < options.tol:
            found = True
            break
    return vrs, niter + 1, found