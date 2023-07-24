class Vector2:
    __slots__ = ("_x", "_y")
    _x: float
    _y: float

    def __init__(self, x, y):
        """[summary]

        Args:
            x ([type]): [description]
            y ([type]): [description]
        """
        self._x = x
        self._y = y

    @property
    def x(self):
        """[summary]

        Returns:
            [type]: [description]
        """
        return self._x

    @property
    def y(self):
        """[summary]

        Returns:
            [type]: [description]
        """
        return self._y

    def dot(self, rhs):
        """[summary]

        Args:
            rhs ([type]): [description]

        Returns:
            [type]: [description]

        Examples:
            >>> v1 = Vector2(1, 2)
            >>> v2 = Vector2(3, 4)
            >>> v1.dot(v2)
            11
        """
        return self._x * rhs._x + self._y * rhs._y

    def __iadd__(self, rhs):
        """[summary]

        Args:
            rhs ([type]): [description]

        Returns:
            [type]: [description]

        Examples:
            >>> v1 = Vector2(1, 2)
            >>> v2 = Vector2(3, 4)
            >>> v1 += v2
            >>> print(v1)
            <4, 6>
        """
        self._x += rhs.x
        self._y += rhs.y
        return self

    def __add__(self, rhs):
        """[summary]

        Args:
            rhs ([type]): [description]

        Returns:
            [type]: [description]

        Examples:
            >>> v1 = Vector2(1, 2)
            >>> v2 = Vector2(3, 4)
            >>> print(v1 + v2)
            <4, 6>
            >>> print(v1)
            <1, 2>
        """
        return Vector2(self.x + rhs.x, self.y + rhs.y)

    def __isub__(self, rhs):
        """[summary]

        Args:
            rhs ([type]): [description]

        Returns:
            [type]: [description]

        Examples:
            >>> v1 = Vector2(1, 2)
            >>> v2 = Vector2(3, 4)
            >>> v1 -= v2
            >>> print(v1)
            <-2, -2>
            >>> print(v2)
            <3, 4>
        """
        self._x -= rhs.x
        self._y -= rhs.y
        return self

    def __sub__(self, rhs):
        """[summary]

        Args:
            rhs ([type]): [description]

        Returns:
            [type]: [description]

        Examples:
            >>> v1 = Vector2(1, 2)
            >>> v2 = Vector2(3, 4)
            >>> print(v1 - v2)
            <-2, -2>
            >>> print(v1)
            <1, 2>
            >>> print(v2)
            <3, 4>
        """
        return Vector2(self.x - rhs.x, self.y - rhs.y)

    def __imul__(self, alpha: float):
        """[summary]

        Args:
            alpha (float): scalar

        Returns:
            [type]: [description]

        Examples:
            >>> v1 = Vector2(1, 2)
            >>> v1 *= 2
            >>> print(v1)
            <2, 4>
        """
        self._x *= alpha
        self._y *= alpha
        return self

    def __mul__(self, alpha: float):
        """[summary]

        Args:
            alpha (float): scalar

        Returns:
            [type]: [description]

        Examples:
            >>> v1 = Vector2(1, 2)
            >>> print(v1 * 2)
            <2, 4>
            >>> print(v1)
            <1, 2>
        """
        return Vector2(self.x * alpha, self.y * alpha)

    def __itruediv__(self, alpha: float):
        """[summary]

        Args:
            alpha (float): scalar

        Returns:
            [type]: [description]

        Examples:
            >>> v1 = Vector2(1, 2)
            >>> v1 /= 2
            >>> print(v1)
            <0.5, 1.0>
            >>> print(v1)
            <0.5, 1.0>
            >>> v1 /= 0
            Traceback (most recent call last):
                ...
            ZeroDivisionError: float division by zero
            >>> print(v1)
            <0.5, 1.0>
            >>> print(v1)
            <0.5, 1.0>
            >>> v1 /= 1
            >>> print(v1)
            <0.5, 1.0>
            >>> print(v1)
            <0.5, 1.0>
            >>> v1 /= 2
            >>> print(v1)
            <0.25, 0.5>
            >>> print(v1)
            <0.25, 0.5>
            >>> v1 /= 2
            >>> print(v1)
            <0.125, 0.25>
        """
        self._x /= alpha
        self._y /= alpha
        return self

    def __truediv__(self, alpha: float):
        """[summary]

        Args:
            alpha (float): scalar

        Returns:
            [type]: [description]

        Examples:
            >>> v1 = Vector2(1, 2)
            >>> print(v1 / 2)
            <0.5, 1.0>
            >>> print(v1)
            <1, 2>
            >>> v1 /= 1
            >>> print(v1)
            <1.0, 2.0>
            >>> v1 /= 2
            >>> print(v1)
            <0.5, 1.0>
            >>> print(v1)
            <0.5, 1.0>
            >>> v1 /= 2
            >>> print(v1)
            <0.25, 0.5>
        """
        return Vector2(self.x / alpha, self.y / alpha)

    def __str__(self):
        """[summary]

        Returns:
            [type]: [description]

        Examples:
            >>> v1 = Vector2(1, 2)
            >>> print(v1)
            <1, 2>
            >>> v2 = Vector2(3, 4)
            >>> print(v2)
            <3, 4>
            >>> v3 = Vector2(5, 6)
            >>> print(v3)
            <5, 6>
        """
        # return "<{self.x}, {self.y}>".format(self=self)
        return f"<{self.x}, {self.y}>"


if __name__ == "__main__":
    v = Vector2(3.0, 4.0)
    w = Vector2(5.0, 6.0)
    print(v.dot(w))
