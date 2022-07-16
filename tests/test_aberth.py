from bairstow.aberth import (
    aberth,
    aberth_autocorr,
    initial_aberth,
    initial_aberth_autocorr,
)
from bairstow.rootfinding import Options


def test_aberth1():
    h = [5.0, 2.0, 9.0, 6.0, 2.0]
    z0s = initial_aberth(h)
    _, niter, found = aberth(h, z0s)
    print([niter, found])
    assert niter <= 8


def test_aberth2():
    h = [10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0]
    z0s = initial_aberth(h)
    zs, niter, found = aberth(h, z0s)
    print([niter, found])
    print([z for z in zs])
    assert niter <= 8


r = [
    -0.00196191,
    -0.00094597,
    -0.00023823,
    0.00134667,
    0.00380494,
    0.00681596,
    0.0097864,
    0.01186197,
    0.0121238,
    0.00985211,
    0.00474894,
    -0.00281751,
    -0.01173923,
    -0.0201885,
    -0.02590168,
    -0.02658216,
    -0.02035729,
    -0.00628271,
    0.01534627,
    0.04279982,
    0.0732094,
    0.10275561,
    0.12753013,
    0.14399228,
    0.15265722,
    0.14399228,
    0.12753013,
    0.10275561,
    0.0732094,
    0.04279982,
    0.01534627,
    -0.00628271,
    -0.02035729,
    -0.02658216,
    -0.02590168,
    -0.0201885,
    -0.01173923,
    -0.00281751,
    0.00474894,
    0.00985211,
    0.0121238,
    0.01186197,
    0.0097864,
    0.00681596,
    0.00380494,
    0.00134667,
    -0.00023823,
    -0.00094597,
    -0.00196191,
]


def test_aberth_fir():
    z0s = initial_aberth(r)
    opt = Options()
    opt.tol = 1e-8
    zs, niter, found = aberth(r, z0s, opt)
    print([niter, found])
    for z in zs:
        print(z)
    assert niter <= 11


def test_aberth_autocorr_fir():
    z0s = initial_aberth_autocorr(r)
    opt = Options()
    opt.tol = 1e-14
    zs, niter, found = aberth_autocorr(r, z0s, opt)
    print([niter, found])
    for z in zs:
        print(z)
    assert niter <= 9


# def test_aberth_fir_lds():
#     r = [
#         -0.00196191,
#         -0.00094597,
#         -0.00023823,
#         0.00134667,
#         0.00380494,
#         0.00681596,
#         0.0097864,
#         0.01186197,
#         0.0121238,
#         0.00985211,
#         0.00474894,
#         -0.00281751,
#         -0.01173923,
#         -0.0201885,
#         -0.02590168,
#         -0.02658216,
#         -0.02035729,
#         -0.00628271,
#         0.01534627,
#         0.04279982,
#         0.0732094,
#         0.10275561,
#         0.12753013,
#         0.14399228,
#         0.15265722,
#         0.14399228,
#         0.12753013,
#         0.10275561,
#         0.0732094,
#         0.04279982,
#         0.01534627,
#         -0.00628271,
#         -0.02035729,
#         -0.02658216,
#         -0.02590168,
#         -0.0201885,
#         -0.01173923,
#         -0.00281751,
#         0.00474894,
#         0.00985211,
#         0.0121238,
#         0.01186197,
#         0.0097864,
#         0.00681596,
#         0.00380494,
#         0.00134667,
#         -0.00023823,
#         -0.00094597,
#         -0.00196191,
#     ]
#     z0s = initial_aberth_lds(r)
#     opt = Options()
#     opt.tol = 1e-8
#     zs, niter, found = aberth(r, z0s, opt)
#     print([niter, found])
#     print([z for z in zs])
#     assert niter <= 13
