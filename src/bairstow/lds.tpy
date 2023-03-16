from math import cos, pi, sin
from typing import List

TWO_PI = 2.0 * pi


def vdc(k: int, base: int) -> float:
    """_summary_

    Args:
        k (int): _description_
        base (int): _description_

    Returns:
        float: _description_
    """
    vdc = 0.0
    denom = 1.0
    while k != 0:
        denom *= base
        remainder = k % base
        k //= base
        vdc += remainder / denom
    return vdc


class Vdcorput:
    """Van der Corput sequence generator

    Examples:
        >>> vgen = Vdcorput(2)
        >>> vgen.reseed(0)
        >>> for _ in range(10):
        ...     print(vgen.pop())
        ...
        0.5
        0.25
        0.75
        0.125
        0.625
        0.375
        0.875
        0.0625
        0.5625
        0.3125
    """

    count: int
    base: int

    def __init__(self, base: int = 2):
        """_summary_

        Args:
            base (int, optional): _description_. Defaults to 2.
        """
        self.count = 0
        self.base = base

    def pop(self) -> float:
        """_summary_

        Returns:
            float: _description_
        """
        self.count += 1
        return vdc(self.count, self.base)

    # [allow(dead_code)]
    def reseed(self, seed: int):
        """_summary_

        Args:
            seed (int): _description_
        """
        self.count = seed


class Circle:
    """Circle sequence generator

    Examples:
        >>> cgen = Circle(2)
        >>> cgen.reseed(0)
        >>> for _ in range(10):
        ...     print(cgen.pop())
        ...
        [1.2246467991473532e-16, -1.0]
        [1.0, 6.123233995736766e-17]
        [-1.0, -1.8369701987210297e-16]
        [0.7071067811865475, 0.7071067811865476]
        [-0.7071067811865475, -0.7071067811865477]
        [0.7071067811865476, -0.7071067811865475]
        [-0.7071067811865477, 0.7071067811865474]
        [0.3826834323650898, 0.9238795325112867]
        [-0.38268343236508967, -0.9238795325112868]
        [0.9238795325112867, -0.3826834323650897]
    """

    vdc: Vdcorput

    def __init__(self, base: int):
        """_summary_

        Args:
            base (int): _description_
        """
        self.vdc = Vdcorput(base)

    def pop(self) -> List[float]:
        """_summary_

        Returns:
            List[float]: _description_
        """
        theta = self.vdc.pop() * TWO_PI  # map to [0, 2*pi]
        return [sin(theta), cos(theta)]

    # [allow(dead_code)]
    def reseed(self, seed: int):
        """_summary_

        Args:
            seed (int): _description_
        """
        self.vdc.reseed(seed)


if __name__ == "__main__":
    base = [2, 3, 5, 7]

    vgen = Vdcorput(2)
    for _ in range(10):
        print(vgen.pop())

    cgen = Circle(2)
    for _ in range(10):
        print(cgen.pop())
