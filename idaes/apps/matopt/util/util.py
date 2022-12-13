#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################


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
    return (
        (x[0] - y[0] < atol)
        and (x[0] - y[0] > -atol)
        and (x[1] - y[1] < atol)
        and (x[1] - y[1] > -atol)
        and (x[2] - y[2] < atol)
        and (x[2] - y[2] > -atol)
    )


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
