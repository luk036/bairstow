from .vector2 import vector2


class matrix2(vector2):
    def __init__(self, x: vector2, y: vector2):
        """[summary]

        Args:
            x (vector2): [description]
            y (vector2): [description]
        """
        vector2.__init__(self, x, y)

    def mdot(self, rhs: vector2) -> vector2:
        """matrix-vector product

        Args:
            rhs (vector2): [description]

        Returns:
            vector2: [description]
        """
        return vector2(self._x.dot(rhs), self._y.dot(rhs))

    def det(self) -> float:
        """determinant

        Returns:
            float: [description]
        """
        a11, a12 = self.x.x, self.x.y
        a21, a22 = self.y.x, self.y.y
        return a11 * a22 - a12 * a21
