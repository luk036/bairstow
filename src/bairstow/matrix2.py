from .vector2 import Vector2


class Matrix2:
    _x: Vector2
    _y: Vector2

    def __init__(self, x: Vector2, y: Vector2):
        """[summary]

        Args:
            x (Vector2): [description]
            y (Vector2): [description]
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

    def mdot(self, rhs: Vector2) -> Vector2:
        """matrix-vector product

        Args:
            rhs (Vector2): [description]

        Returns:
            Vector2: [description]
        """
        return Vector2(self._x.dot(rhs), self._y.dot(rhs))

    def det(self) -> float:
        """determinant

        Returns:
            float: [description]
        """
        a11, a12 = self.x.x, self.x.y
        a21, a22 = self.y.x, self.y.y
        return a11 * a22 - a12 * a21


if __name__ == "__main__":
    v = Vector2(3.0, 4.0)
    w = Vector2(5.0, 6.0)
    m = Matrix2(v, w)
    print(m.det())
