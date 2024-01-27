
[PYTHON_CODE]
"""
The file 'vector2.py' contains a class 'Vector2'. The code below demonstrates how to use the Vector2 class by creating several instances of it and using the various methods to perform calculations with them.
"""

#!/usr/bin/env python3
import math
from typing import Tuple, Any

class Vector2:
    def __init__(self, x: float = 0.0, y: float = 0.0):
        """
        The function initializes the Vector2 object with a specified x and y value. If no argument is provided, it initializes the x and y values to zero.

        :param x: The parameter "x" is a scalar value that represents the x component of the vector.
        :type x: float
        :param y: The parameter "y" is a scalar value that represents the y component of the vector.
        :type y: float
        """
        self._x = x
        self._y = y
    
    @property
    def x(self) -> float:
        return self._x
    
    @x.setter
    def x(self, value: float):
        if not isinstance(value, (float, int)):
            raise TypeError("The 'x' attribute must be a float or an integer.")
        self._x = value
    
    @property
    def y(self) -> float:
        return self._y
    
    @y.setter
    def y(self, value: float):
        if not isinstance(value, (float, int)):
            raise TypeError("The 'y' attribute must be a float or an integer.")
        self._y = value

    @staticmethod
    def distance(v1: Any, v2: Any) -> float:
        """
        The function calculates the distance between two vectors. If no argument is provided, it assumes that both arguments are instances of Vector2 and calculates the distance between them.

        :param v1: The parameter "v1" represents one of the vectors to calculate the distance between.
        :type v1: Vector2
        :param v2: The parameter "v2" represents another vector to calculate the distance between.
        :type v2: Vector2
        :return: The function returns a scalar value representing the distance between two vectors. If no argument is provided, it returns the distance between both vectors.
        :rtype: float
        """
        if not isinstance(v1, Vector2) or not isinstance(v2, Vector2):
            raise TypeError("The 'distance' function expects both arguments to be instances of the Vector2 class.")
        return math.sqrt((v2.x - v1.x) ** 2 + (v2.y - v1.y) ** 2)
    
    @staticmethod
    def magnitude(v: Any) -> float:
        """
        The function calculates the magnitude of a given vector. If no argument is provided, it assumes that the argument is an instance of Vector2 and calculates its magnitude.

        :param v: The parameter "v" represents the vector whose magnitude is to be calculated.
        :type v: float
        :return: The function returns a scalar value representing the magnitude of a given vector. If no argument is provided, it returns the magnitude of the Vector2 object.
        :rtype: float
        """
        if not isinstance(v, Vector2):
            raise TypeError("The 'magnitude' function expects an argument to be an instance of the Vector2 class.")
        return math.sqrt(v.x ** 2 + v.y ** 2)
    
    @staticmethod
    def angle(v1: Any, v2: Any) -> float:
        """
        The function calculates the angle between two vectors in radians. If no argument is provided, it assumes that both arguments are instances of Vector2 and calculates the angle between them.

        :param v1: The parameter "v1" represents one of the vectors to calculate the angle between.
        :type v1: Vector2
        :param v2: The parameter "v2" represents another vector to calculate the angle between.
        :type v2: Vector2
        :return: The function returns a scalar value representing the angle between two vectors in radians. If no argument is provided, it returns the angle between both vectors in radians.
        :rtype: float
        """
        if not isinstance(v1, Vector2) or not isinstance(v2, Vector2):
            raise TypeError("The 'angle' function expects both arguments to be instances of the Vector2 class.")
        return math.acos((v1.x * v2.x + v1.y * v2.y) / (Vector2.magnitude(v1) * Vector2.magnitude(v2)))
    
    def add(self, other: Any) -> Tuple[float, float]:
        """
        The function adds two vectors together and returns their sum as a tuple. If no argument is provided, it assumes that the first argument is an instance of Vector2 and calculates its addition with a specified scalar value.

        :param other: The parameter "other" represents the second vector to add together or the scalar value to add to the vector.
        :type other: float or Vector2
        :return: The function returns a tuple representing the sum of two vectors (or the result of adding a scalar value to a vector). If no argument is provided, it returns the sum of both vectors (i.e., the addition of the first vector with a specified scalar value).
        :rtype: Tuple[float, float]
        """
        if not isinstance(other, Vector2) and not isinstance(other, (float, int)):
            raise TypeError("The 'add' function expects either one argument to be an instance of the Vector2 class or another argument to be a scalar value.")
        
        if isinstance(other, Vector2):
            return self.x + other.x, self.y + other.y
        
        elif isinstance(other, (float, int)):
            return self.x + other, self.y + other
    
    def sub(self, other: Any) -> Tuple[float, float]:
        """
        The function subtracts a specified vector or scalar value from another vector and returns their difference as a tuple. If no argument is provided, it assumes that the first argument is an instance of Vector2 and calculates its subtraction with a specified scalar value.

        :param other: The parameter "other" represents the second vector to subtract together or the scalar value to subtract from the vector.
        :type other: float or Vector2
        :return: The function returns a tuple representing the difference of two vectors (or the result of subtracting a scalar value from a vector). If no argument is provided, it returns the difference of both vectors (i.e., the subtraction of the first vector with a specified scalar value).
        :rtype: Tuple[float, float]
        """
        if not isinstance(other, Vector2) and not isinstance(other, (float, int)):
            raise TypeError("The 'sub' function expects either one argument to be an instance of the Vector2 class or another argument to be a scalar value.")
        
        elif isinstance(other, Vector2):
            return self.x - other.x, self.y - other.y
        
        elif isinstance(other, (float, int)):
            return self.x - other, self.y - other
    
    def mul(self, other: Any) -> Tuple[float, float]:
        """
        The function multiplies two vectors together and returns their product as a tuple. If no argument is provided, it assumes that the first argument is an instance of Vector2 and calculates its multiplication with a specified scalar value.

        :param other: The parameter "other" represents the second vector to multiply together or the scalar value to multiply from the vector.
        :type other: float or Vector2
        :return: The function returns a tuple representing the product of two vectors (or the result of multiplying a scalar value with a vector). If no argument is provided, it returns the product of both vectors (i.e., the multiplication of the first vector with a specified scalar value).
        :rtype: Tuple[float, float]
        """
        if not isinstance(other, Vector2) and not isinstance(other, (float, int)):
            raise TypeError("The 'mul' function expects either one argument to be an instance of the Vector2 class or another argument to be a scalar value.")
        
        elif isinstance(other, Vector2):
            return self.x * other.x, self.y * other.y
        
        elif isinstance(other, (float, int)):
            return self.x * other, self.y * other
    
    def div(self, other: Any) -> Tuple[float, float]:
        """
        The function divides one vector by another and returns their ratio as a tuple. If no argument is provided, it assumes that the first argument is an instance of Vector2 and calculates its division with a specified scalar value.

        :param other: The parameter "other" represents the second vector to divide or the scalar value to divide from the vector.
        :type other: float or Vector2
        :return: The function returns a tuple representing the ratio of two vectors (or the result of dividing a vector by a specified scalar value). If no argument is provided, it returns the ratio of both vectors (i.e., the division of the first vector by a specified scalar value).
        :rtype: Tuple[float, float]
        """
        if not isinstance(other, Vector2) and not isinstance(other, (float, int)):
            raise TypeError("The 'div' function expects either one argument to be an instance of the Vector2 class or another argument to be a scalar value.")
        
        elif isinstance(other, Vector2):
            return self.x / other.x, self.y / other.y
        
        elif isinstance(other, (float, int)):
            return self.x / other, self.y / other
    
    def normalize(self) -> Tuple[float, float]:
        """
        The function normalizes a vector and returns its magnitude as a scalar value. If no argument is provided, it assumes that the first argument is an instance of Vector2 and calculates its normalization.

        :return: The function returns a scalar value representing the magnitude of a normalized vector. If no argument is provided, it returns the magnitude of the Vector2 object (i.e., the result of normalizing it).
        :rtype: float
        """
        if not isinstance(self, Vector2):
            raise TypeError("The 'normalize' function expects one argument to be an instance of the Vector2 class.")
        
        return (self.x / self.mag), (self.y / self.mag)
    
    def rotate(self, angle: float = 0) -> Tuple[float, float]:
        """
        The function rotates a vector by an angle in degrees and returns its new coordinates as a tuple. If no argument is provided, it assumes that the first argument is an instance of Vector2 and calculates its rotation with no specified angle.

        :param angle: The parameter "angle" represents the angle (in degrees) to rotate the vector by. Defaults to 0.
        :type angle: float
        :return: The function returns a tuple representing the new coordinates of the rotated vector. If no argument is provided, it returns the coordinates of the rotated Vector2 object.
        :rtype: Tuple[float, float]
        """
        if not isinstance(self, Vector2):
            raise TypeError("The 'rotate' function expects one argument to be an instance of the Vector2 class.")
        
        cos = math.cos(math.radians(angle))
        sin = math.sin(math.radians(angle))
        x_new = self.x * cos - self.y * sin
        y_new = self.x * sin + self.y * cos
        
        return (x_new, y_new)
    
    def distance(self, other: Vector2) -> float:
        """
        The function calculates the Euclidean distance between two vectors and returns it as a scalar value. If no argument is provided, it assumes that the first argument is an instance of Vector2 and calculates its distance from another specified vector.

        :param other: The parameter "other" represents the second vector to calculate the distance between.
        :type other: Vector2
        :return: The function returns a scalar value representing the Euclidean distance between two vectors. If no argument is provided, it returns the distance from another specified vector.
        :rtype: float
        """
        if not isinstance(self, Vector2) or not isinstance(other, Vector2):
            raise TypeError("The 'distance' function expects one or two arguments to be instances of the Vector2 class.")
        
        return ((self.x - other.x) ** 2 + (self.y - other.y) ** 2) ** 0.5
    
    def dot(self, other: Vector2) -> float:
        """
        The function calculates the dot product between two vectors and returns it as a scalar value. If no argument is provided, it assumes that the first argument is an instance of Vector2 and calculates its dot product with another specified vector.

        :param other: The parameter "other" represents the second vector to calculate the dot product between.
        :type other: Vector2
        :return: The function returns a scalar value representing the dot product between two vectors. If no argument is provided, it returns the dot product of another specified vector.
        :rtype: float
        """
        if not isinstance(self, Vector2) or not isinstance(other, Vector2):
            raise TypeError("The 'dot' function expects one or two arguments to be instances of the Vector2 class.")
        
        return self.x * other.x + self.y * other.y
    
    def cross(self, other: Vector2) -> float:
        """
        The function calculates the cross product between two vectors and returns it as a scalar value. If no argument is provided, it assumes that the first argument is an instance of Vector2 and calculates its cross product with another specified vector.

        :param other: The parameter "other" represents the second vector to calculate the cross product between.
        :type other: Vector2
        :return: The function returns a scalar value representing the cross product between two vectors. If no argument is provided, it returns the cross product of another specified vector.
        :rtype: float
        """
        if not isinstance(self, Vector2) or not isinstance(other, Vector2):
            raise TypeError("The 'cross' function expects one or two arguments to be instances of the Vector2 class.")
        
        return self.x * other.y - self.y * other.x
    
    def project(self, normal: Vector2) -> float:
        """
        The function projects a vector onto another vector and returns it as a scalar value. If no argument is provided, it assumes that the first argument is an instance of Vector2 and calculates its projection onto another specified vector.

        :param normal: The parameter "normal" represents the vector to project onto.
        :type normal: Vector2
        :return: The function returns a scalar value representing the projection of one vector onto another. If no argument is provided, it returns the projection onto another specified vector.
        :rtype: float
        """
        if not isinstance(self, Vector2) or not isinstance(normal, Vector2):
            raise TypeError("The 'project' function expects one or two arguments to be instances of the Vector2 class.")
        
        return normal.dot(self) / normal.mag
    
    def reflect(self, normal: Vector2) -> float:
        """
        The function reflects a vector around another vector and returns it as a scalar value. If no argument is provided, it assumes that the first argument is an instance of Vector2 and calculates its reflection around another specified vector.

        :param normal: The parameter "normal" represents the vector to reflect around.
        :type normal: Vector2
        :return: The function returns a scalar value representing the reflected vector. If no argument is provided, it returns the reflection around another specified vector.
        :rtype: float
        """
        if not isinstance(self, Vector2) or not isinstance(normal, Vector2):
            raise TypeError("The 'reflect' function expects one or two arguments to be instances of the Vector2 class.")
        
        return self - 2 * normal.project(self) * normal
    
    def __add__(self, other: Vector2) -> Vector2:
        """
        The function adds another vector to this one and returns it as a new instance of Vector2. If no argument is provided, it assumes that the first argument is an instance of Vector2 and calculates its addition with another specified vector.

        :param other: The parameter "other" represents the vector to add to this one.
        :type other: Vector2
        :return: The function returns a new instance of Vector2 representing the sum of two vectors. If no argument is provided, it returns the addition with another specified vector.
        :rtype: Vector2
        """
        if not isinstance(self, Vector2) or not isinstance(other, Vector2):
            raise TypeError("The '+' operator expects one or two arguments to be instances of the Vector2 class.")
        
        return Vector2((self.x + other.x), (self.y + other.y))
    
    def __sub__(self, other: Vector2) -> Vector2:
        """
        The function subtracts another vector from this one and returns it as a new instance of Vector2. If no argument is provided, it assumes that the first argument is an instance of Vector2 and calculates its subtraction with another specified vector.

        :param other: The parameter "other" represents the vector to subtract from this one.
        :type other: Vector2
        :return: The function returns a new instance of Vector2 representing the difference between two vectors. If no argument is provided, it returns the subtraction with another specified vector.
        :rtype: Vector2
        """
        if not isinstance(self, Vector2) or not isinstance(other, Vector2):
            raise TypeError("The '-' operator expects one or two arguments to be instances of the Vector2 class.")
        
        return Vector2((self.x - other.x), (self.y - other.y))
    
    def __mul__(self, other: float) -> Vector2:
        """
        The function multiplies this vector by a scalar and returns it as a new instance of Vector2. If no argument is provided, it assumes that the first argument is an instance of Vector2 and calculates its multiplication with another specified scalar.

        :param other: The parameter "other" represents the scalar to multiply this vector by.
        :type other: float
        :return: The function returns a new instance of Vector2 representing the product of this vector and another scalar. If no argument is provided, it returns the multiplication with another specified scalar.
        :rtype: Vector2
        """
        if not isinstance(self, Vector2) or not isinstance(other, float):
            raise TypeError("The '*' operator expects one or two arguments to be instances of the Vector2 class and a scalar.")
        
        return Vector2((self.x * other), (self.y * other))
    
    def __truediv__(self, other: float) -> Vector2:
        """
        The function divides this vector by a scalar and returns it as a new instance of Vector2. If no argument is provided, it assumes that the first argument is an instance of Vector2 and calculates its division by another specified scalar.

        :param other: The parameter "other" represents the scalar to divide this vector by.
        :type other: float
        :return: The function returns a new instance of Vector2 representing the quotient of this vector and another scalar. If no argument is provided, it returns the division by another specified scalar.
        :rtype: Vector2
        """
        if not isinstance(self, Vector2) or not isinstance(other, float):
            raise TypeError("The '/' operator expects one or two arguments to be instances of the Vector2 class and a scalar.")
        
        return Vector2((self.x / other), (self.y / other))
    
    def __str__(self) -> str:
        """
        The function returns this vector as a string representation.

        :return: The function returns the string representation of this vector.
        :rtype: str
        """
        return f"{self.x}, {self.y}"

