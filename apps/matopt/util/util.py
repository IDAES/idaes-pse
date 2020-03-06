import numpy as np


def isZero(x, atol):
    """Determine if a floating point number is equal to zero.

    Args:
    x (float): a number
    atol (float): absolute tolerance from zero

    Returns:
    (bool) true if the number is sufficiently close to zero
    """
    return atol > x > -atol


def areEqual(x, y, atol):
    """Determine if two floating point numbers are equal.

    Args:
    x (float): a number
    y (float): a number
    atol (float): absolute tolerance for equality

    Returns:
    (bool) true if the two numbers are sufficiently close
    """
    return isZero(x - y, atol)


def myArrayEq(x, y, atol):
    """Determine if two numpy arrays of floating point numbers are equal.

    Args:
    x (numpy.ndarray): a numpy array 
    y (numpy.ndarray): a numpy array
    atol (float): absolute tolerance for equality

    Returns:
    (bool) true if the two arrays are equal
    """
    return (np.abs(x - y) < atol).all()


myPointEq = myArrayEq


def myPointsEq(x, y, atol):
    """Determine if two lists of numpy arrays are equal.

    Args:
    x (list<numpy.ndarray>): a list of numpy arrays
    y (list<numpy.ndarray>): a list of numpy arrays
    atol (float): absolute tolerance for equality

    Returns:
    (bool) true if the two arrays are equal
    """
    if len(x) != len(y):
        return False
    for i in range(len(x)):
        if not myArrayEq(x[i], y[i], atol):
            return False
    return True


def ListHasPoint(L, P, atol):
    """Determine if a list of numpy arrays contains a specific point.

    Args:
    x (list<numpy.ndarray>): a list of numpy arrays
    y (list<numpy.ndarray>): a list of numpy arrays
    atol (float): absolute tolerance for equality

    Returns:
    (bool) true if the two arrays are equal
    """
    for l in L:
        if myArrayEq(l, P, atol):
            return True
    return False
