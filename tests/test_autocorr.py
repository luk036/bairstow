from bairstow.autocorr import initial_autocorr, pbairstow_autocorr


def test_autocorr():
    # vA = vector2(0.1, 1.2)
    # vA1 = vector2(2.3, 3.4)
    # vr = vector2(4.5, 5.6)
    # vrj = vector2(6.7, 7.8)
    # vAorig, vA1orig = suppress_orig(vA, vA1, vr, vrj)
    # print(check_newton(vAorig, vA1orig, vr))
    # vA1 = suppress(vA, vA1, vr, vrj)
    # print(check_newton(vA, vA1, vr))
    h = [10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0]
    vr0s = initial_autocorr(h)
    # print(vr0s[1])
    # vA, pb = horner(h, vr0s[1])
    # print(vA)
    # print(pb)
    # vA1, _ = horner(pb, vr0s[1])
    # print(vA1)
    _, niter, found = pbairstow_autocorr(h, vr0s)
    print([niter, found])
    # print([find_rootq(-r.x, -r.y) for r in vrs])
    assert niter <= 31
