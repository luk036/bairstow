use std::collections::HashMap;
use std::*;

use vector2::Vector2;
struct Matrix2 {}

impl Matrix2 {
    fn __init__(&self, x: Vector2, y: Vector2) {
        "[summary]

        Args:
            x (Vector2): [description]
            y (Vector2): [description]
        ";
        Vector2::__init__(self, x, y);
    }
    fn mdot(&self, rhs: Vector2) -> Vector2 {
        "matrix-vector product

        Args:
            rhs (Vector2): [description]

        Returns:
            Vector2: [description]
        ";
        return Vector2(self._x.dot(rhs), self._y.dot(rhs));
    }
    fn det(&self) -> f32 {
        "determinant

        Returns:
            float: [description]
        ";
        let (a11, a12) = (self.x.x, self.x.y);
        let (a21, a22) = (self.y.x, self.y.y);
        return ((a11 * a22) - (a12 * a21));
    }
}
