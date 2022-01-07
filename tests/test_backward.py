from pytest import approx

from bairstow.aberth import horner_backward, horner_eval


def test_backward2():
    p = [1.0, -6.7980, 2.9948, -0.043686, 0.000089248]
    n = len(p) - 1
    alpha = 6.3256
    P = horner_backward(p, n, alpha)
    assert -P * (alpha ** 5) == approx(-0.013355264987140483)
    assert p[3] == approx(0.006920331351966613)


def test_backward1():
    p = [1.0, -6.7980, 2.9948, -0.043686, 0.000089248]
    n = len(p) - 1
    P = horner_eval(p, n, 6.3256)
    assert P == approx(-0.012701469838522064)
    assert p[3] == approx(-0.0020220560640132265)
