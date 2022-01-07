from bairstow.autocorr import extract_autocorr, initial_autocorr, pbairstow_autocorr
from bairstow.rootfinding import find_rootq


def test_autocorr():
    h = [10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0]
    vr0s = initial_autocorr(h)
    vrs, niter, found = pbairstow_autocorr(h, vr0s)
    print([niter, found])

    for vr in vrs:
        vr = extract_autocorr(vr)

    print(find_rootq(vr) for vr in vrs)

    assert niter <= 7
