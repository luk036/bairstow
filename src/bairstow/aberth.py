from cmath import exp
from math import pi
from typing import List, Tuple, Union

from .lds import VdCorput
from .robin import Robin

# from pytest import approx
from .rootfinding import Options, horner_eval

FoC = Union[float, complex]
PI = pi


def horner_backward(pb: List, degree: int, alpha: FoC) -> FoC:
    """[summary]

    Args:
        pb (List): _description_
        degree (int): _description_
        alpha (FoC): _description_

    Returns:
        FoC: _description_

    Examples:
        >>> coeffs = [1.0, -6.7980, 2.9948, -0.043686, 0.000089248]
        >>> degree = len(coeffs) - 1
        >>> alpha = 6.3256
        >>> p_eval = horner_backward(coeffs, degree, alpha)
        >>> -p_eval * pow(alpha, 5)
        -0.013355264987140483
        >>> coeffs[3]
        0.006920331351966613
    """
    for i in range(2, degree + 2):
        pb[-i] -= pb[-(i - 1)]
        pb[-i] /= -alpha
    return pb[-(degree + 1)]


def initial_aberth(coeffs: List[FoC]) -> List[complex]:
    """[summary]

    Args:
        coeffs (List): [description]

    Returns:
        List: [description]

    Examples:
        >>> h = [5.0, 2.0, 9.0, 6.0, 2.0]
        >>> z0s = initial_aberth(h)
    """
    degree: int = len(coeffs) - 1
    center: FoC = -coeffs[1] / (degree * coeffs[0])
    p_center: FoC = horner_eval(coeffs.copy(), degree, center)
    re: FoC = pow(-p_center, 1.0 / degree)
    # k = 2 * PI / degree
    z0s: List[complex] = []
    vgen = VdCorput(2)
    vgen.reseed(1)
    for _ in range(degree):
        vdc = 2 * PI * vgen.pop()
        z0s += [center + re * exp(vdc * 1j)]
    return z0s


def initial_aberth_orig(coeffs: List[FoC]) -> List[complex]:
    """[summary]

    Args:
        coeffs (List): [description]

    Returns:
        List: [description]

    Examples:
        >>> h = [5.0, 2.0, 9.0, 6.0, 2.0]
        >>> z0s = initial_aberth_orig(h)
    """
    degree: int = len(coeffs) - 1
    center: FoC = -coeffs[1] / (degree * coeffs[0])
    p_center: FoC = horner_eval(coeffs.copy(), degree, center)
    re: FoC = pow(-p_center, 1.0 / degree)
    k = 2 * PI / degree
    z0s: List[complex] = []
    for i in range(degree):
        theta = k * (0.25 + i)
        z0s += [center + re * exp(theta * 1j)]
    return z0s


def aberth(
    coeffs: List[FoC], zs: List[complex], options: Options = Options()
) -> Tuple[List[complex], int, bool]:
    """Aberth's method for polynomial root-finding

                    P ⎛z ⎞
         new          ⎝ i⎠
        z    = z  - ───────
         i      i   P' ⎛z ⎞
                       ⎝ i⎠
    where
                              degree
                            _____
                            ╲
                             ╲    P ⎛z ⎞
                              ╲     ⎝ i⎠
        P' ⎛z ⎞ = P  ⎛z ⎞ -   ╱   ───────
           ⎝ i⎠    1 ⎝ i⎠    ╱    z  - z
                            ╱      i    j
                            ‾‾‾‾‾
                            j ≠ i

    Args:
        coeffs (List): [description]
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
    degree = len(coeffs) - 1
    converged = [False] * M
    robin = Robin(M)
    for niter in range(options.max_iters):
        tol = 0.0
        for i in filter(lambda i: not converged[i], range(M)):
            pb = coeffs.copy()
            p_eval = horner_eval(pb, degree, zs[i])
            tol_i = abs(p_eval)
            if tol_i < options.tol_ind:
                converged[i] = True
                continue
            p1_eval = horner_eval(pb, degree - 1, zs[i])
            tol = max(tol_i, tol)
            # for j in filter(lambda j: j != i, range(M)):  # exclude i
            for j in robin.exclude(i):
                p1_eval -= p_eval / (zs[i] - zs[j])
            zs[i] -= p_eval / p1_eval
        if tol < options.tol:
            return zs, niter, True
    return zs, options.max_iters, False


def initial_aberth_autocorr(coeffs: List[float]) -> List[complex]:
    """[summary]

    Args:
        coeffs (List): [description]

    Returns:
        List: [description]

    Examples:
        >>> h = [5.0, 2.0, 9.0, 6.0, 2.0]
        >>> z0s = initial_aberth_autocorr(h)
    """
    degree: int = len(coeffs) - 1
    re: float = pow(abs(coeffs[-1]), 1.0 / degree)
    # center = -coeffs[1] / (degree * coeffs[0])
    # p_center = horner_eval(coeffs.copy(), degree, center)
    # re = (-p_center) ** (1.0 / degree)
    if abs(re) > 1:
        re = 1 / re
    degree //= 2
    # k = 2 * PI / degree
    z0s = []
    vgen = VdCorput(2)
    vgen.reseed(1)
    for _ in range(degree):
        vdc = 2 * PI * vgen.pop()
        z0s += [re * exp(vdc * 1j)]
    return z0s


def initial_aberth_autocorr_orig(coeffs: List[float]) -> List[complex]:
    """[summary]

    Args:
        coeffs (List): [description]

    Returns:
        List: [description]

    Examples:
        >>> h = [5.0, 2.0, 9.0, 6.0, 2.0]
        >>> z0s = initial_aberth_autocorr_orig(h)
    """
    degree: int = len(coeffs) - 1
    re: float = pow(abs(coeffs[-1]), 1.0 / degree)
    # center = -coeffs[1] / (degree * coeffs[0])
    # p_center = horner_eval(coeffs.copy(), degree, center)
    # re = (-p_center) ** (1.0 / degree)
    if abs(re) > 1:
        re = 1 / re
    degree //= 2
    k = 2 * PI / degree
    z0s = []
    for i in range(degree):
        theta = k * (0.25 + i)
        z0s += [re * exp(theta * 1j)]
    return z0s


def aberth_autocorr(
    coeffs: List[float], zs: List[complex], options=Options()
) -> Tuple[List[complex], int, bool]:
    """[summary]

    Args:
        coeffs (List): [description]
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
    degree: int = len(coeffs) - 1
    converged: List[bool] = [False] * M
    robin = Robin(M)
    for niter in range(options.max_iters):
        tol: float = 0.0
        # exclude converged
        for i in filter(lambda i: not converged[i], range(M)):
            pb = coeffs.copy()
            p_eval = horner_eval(pb, degree, zs[i])
            tol_i = abs(p_eval)
            if tol_i < options.tol_ind:
                converged[i] = True
                continue
            p1_eval = horner_eval(pb, degree - 1, zs[i])
            tol = max(tol_i, tol)
            # for j in filter(lambda j: j != i, range(M)):  # exclude i
            for j in robin.exclude(i):
                p1_eval -= p_eval / (zs[i] - zs[j])
                # for j in range(M):  # exclude i
                zsn = 1.0 / zs[j]
                p1_eval -= p_eval / (zs[i] - zsn)
            zs[i] -= p_eval / p1_eval
            # if abs(zs[i]) > 1.0:  # pick those inside the unit circle
            #     zs[i] = 1.0 / zs[i]
        if tol < options.tol:
            return zs, niter, True
    return zs, options.max_iters, False


# def test_aberth():
#     h = [5.0, 2.0, 9.0, 6.0, 2.0]
#     z0s = initial_aberth(h)
#     zs, niter, found = aberth(h, z0s)
#     assert (niter == 2)
#     assert (found)
#     zs, niter, found = aberth(h, z0s, Options(tol=1e-10))
#     assert (niter == 2)
#     assert (found)
#     zs, niter, found = aberth(h, z0s, Options(max_iters=1))
#     assert (niter == 1)
#     assert (found)
#     zs, niter, found = aberth(h, z0s, Options(max_iters=1, tol=1e-10))
#     assert (niter == 1)
#     assert (found)
#     zs, niter, found = aberth(h, z0s, Options(max_iters=1, tol=1e-11))
#     assert (niter == 0)