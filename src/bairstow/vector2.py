class vector2:
    __slots__ = ("_x", "_y")

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
        """
        return self._x * rhs._x + rhs._y * self._y

    def __iadd__(self, rhs):
        """[summary]

        Args:
            rhs ([type]): [description]

        Returns:
            [type]: [description]
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
        """
        return vector2(self.x + rhs.x, self.y + rhs.y)

    def __isub__(self, rhs):
        """[summary]

        Args:
            rhs ([type]): [description]

        Returns:
            [type]: [description]
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
        """
        return vector2(self.x - rhs.x, self.y - rhs.y)

    def __imul__(self, alpha: float):
        """[summary]

        Args:
            alpha (float): scalar

        Returns:
            [type]: [description]
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
        """
        return vector2(self.x * alpha, self.y * alpha)

    def __itruediv__(self, alpha: float):
        """[summary]

        Args:
            alpha (float): scalar

        Returns:
            [type]: [description]
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
        """
        return vector2(self.x / alpha, self.y / alpha)

    def __str__(self):
        """[summary]

        Returns:
            [type]: [description]
        """
        return "<{self.x}, {self.y}>".format(self=self)
