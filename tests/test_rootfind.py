from bairstow.rootfinding import find_rootq, initial_guess, pbairstow_even


def test_rootfind():
    h = [5.0, 2.0, 9.0, 6.0, 2.0]
    vr0s = initial_guess(h)
    _, niter, found = pbairstow_even(h, vr0s)
    print([niter, found])
    assert niter <= 4


def test_rootfind2():
    h = [10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0]
    vr0s = initial_guess(h)
    vrs, niter, found = pbairstow_even(h, vr0s)
    print([niter, found])
    print(find_rootq(vr) for vr in vrs)
    assert niter <= 11
