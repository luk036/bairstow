# -*- coding: utf-8 -*-
from __future__ import print_function

from bairstow.autocorr import initial_autocorr, pbairstow_autocorr
from bairstow.rootfinding import initial_guess, pbairstow_even


def run_autocorr():
    """[summary]

    Returns:
        [type]: [description]
    """
    p = [10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0]
    vr0s = initial_autocorr(p)
    _, niter, _ = pbairstow_autocorr(p, vr0s)
    return niter


def run_pbairstow():
    """[summary]

    Returns:
        [type]: [description]
    """
    p = [10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0]
    vr0s = initial_guess(p)
    _, niter, _ = pbairstow_even(p, vr0s)
    return niter


def test_autocorr(benchmark):
    """[summary]

    Arguments:
        benchmark ([type]): [description]
    """
    result = benchmark(run_autocorr)
    assert result <= 12


def test_pbairstow(benchmark):
    """[summary]

    Arguments:
        benchmark ([type]): [description]
    """
    result = benchmark(run_pbairstow)
    assert result <= 11
