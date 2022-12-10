use std::collections::HashMap;
use std::*;

struct Vector2 {
    _x: ST0,
    _y: ST1,
}

impl Vector2 {
    const __slots__: _ = ("_x", "_y");
    fn __init__<T0, T1>(&self, x: T0, y: T1) {
        "[summary]

        Args:
            x ([type]): [description]
            y ([type]): [description]
        ";
        self._x = x;
        self._y = y;
    }
    fn x<RT>(&self) -> RT {
        "[summary]

        Returns:
            [type]: [description]
        ";
        return self._x;
    }
    fn y<RT>(&self) -> RT {
        "[summary]

        Returns:
            [type]: [description]
        ";
        return self._y;
    }
    fn dot<T0, RT>(&self, rhs: T0) -> RT {
        "[summary]

        Args:
            rhs ([type]): [description]

        Returns:
            [type]: [description]
        ";
        return ((self._x * rhs._x) + (rhs._y * self._y));
    }
    fn __iadd__<T0, RT>(&self, rhs: T0) -> RT {
        "[summary]

        Args:
            rhs ([type]): [description]

        Returns:
            [type]: [description]
        ";
        self._x += rhs.x;
        self._y += rhs.y;
        return self;
    }
    fn __add__<T0, RT>(&self, rhs: T0) -> RT {
        "[summary]

        Args:
            rhs ([type]): [description]

        Returns:
            [type]: [description]
        ";
        return Vector2((self.x + rhs.x), (self.y + rhs.y));
    }
    fn __isub__<T0, RT>(&self, rhs: T0) -> RT {
        "[summary]

        Args:
            rhs ([type]): [description]

        Returns:
            [type]: [description]
        ";
        self._x -= rhs.x;
        self._y -= rhs.y;
        return self;
    }
    fn __sub__<T0, RT>(&self, rhs: T0) -> RT {
        "[summary]

        Args:
            rhs ([type]): [description]

        Returns:
            [type]: [description]
        ";
        return Vector2((self.x - rhs.x), (self.y - rhs.y));
    }
    fn __imul__<RT>(&self, alpha: f32) -> RT {
        "[summary]

        Args:
            alpha (float): scalar

        Returns:
            [type]: [description]
        ";
        self._x *= alpha;
        self._y *= alpha;
        return self;
    }
    fn __mul__<RT>(&self, alpha: f32) -> RT {
        "[summary]

        Args:
            alpha (float): scalar

        Returns:
            [type]: [description]
        ";
        return Vector2((self.x * alpha), (self.y * alpha));
    }
    fn __itruediv__<RT>(&self, alpha: f32) -> RT {
        "[summary]

        Args:
            alpha (float): scalar

        Returns:
            [type]: [description]
        ";
        self._x /= alpha;
        self._y /= alpha;
        return self;
    }
    fn __truediv__<RT>(&self, alpha: f32) -> RT {
        "[summary]

        Args:
            alpha (float): scalar

        Returns:
            [type]: [description]
        ";
        return Vector2((self.x / alpha), (self.y / alpha));
    }
    fn __str__<RT>(&self) -> RT {
        "[summary]

        Returns:
            [type]: [description]
        ";
        return "<{self.x}, {self.y}>".format(self);
    }
}
