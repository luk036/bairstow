from .vector2 import Vector2


class Matrix2(Vector2):
    def __init__(self, x: Vector2, y: Vector2):
        """[summary]

        Args:
            x (Vector2): [description]
            y (Vector2): [description]
        """
        Vector2.__init__(self, x, y)

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
