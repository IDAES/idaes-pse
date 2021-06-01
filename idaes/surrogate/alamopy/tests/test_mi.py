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
# This file applies the ALAMOpython module to the six hump camel problem
# More information on this function can be found at :
#            https://www.sfu.ca/~ssurjano/camel6.html
# This problem utilizes ALAMO's sampling features
from idaes.surrogate import alamopy

# Import additional python modules for creating the synthetic data
import math
import numpy as np
import pytest


def sixcamel(*x):
    x1, x2 = x
    t1 = np.multiply(
        4.0 - 2.1 * np.power(x1, 2) + np.divide(np.power(x1, 4), 3.0), np.power(x1, 2)
    )
    t2 = np.multiply(4 * np.power(x2, 2) - 4, np.power(x2, 2))
    z = t1 + np.multiply(x1, x2) + t2
    return z

@pytest.mark.unit
def test_multiple_input():

    if not alamopy.has_alamo():
        return False

    # Specify number of points to be used in the training set
    # Validation data can be provided optionally
    ndata = 10
    x = np.random.uniform([-2, -1], [2, 1], (ndata, 2))
    z = np.zeros((ndata, 2))
    # specify simulator as examples.sixcamel
    sim = sixcamel
    for i in range(ndata):
        z[i,0] = sim(x[i][0], x[i][1])
        z[i,1] = sim(x[i][0], x[i][1])

    # # Use alamopy's python function wrapper to avoid using ALAMO's I/O format
    almsim = alamopy.wrapwriter(sim)


    # NON GENERIC
    alamo_settings = {'almname': "cam6",
                'monomialpower':(1, 2, 3, 4, 5, 6),
                'multi2power':(1, 2),
                'simulator':almsim,
                'expandoutput':True}

    res = alamopy.doalamo(x, z)

@pytest.mark.unit
def test_single_input_CV():

    if not alamopy.has_alamo():
        return False

    # Specify number of points to be used in the training set
    # Validation data can be provided optionally
    ndata = 10
    x = np.random.uniform([-2, -1], [2, 1], (ndata, 2))
    z = np.zeros((ndata, 1))
    # specify simulator as examples.sixcamel
    sim = sixcamel
    for i in range(ndata):
        z[i,0] = sim(x[i][0], x[i][1])
        # z[i,1] = sim(x[i][0], x[i][1])

    # # Use alamopy's python function wrapper to avoid using ALAMO's I/O format
    almsim = alamopy.wrapwriter(sim)


    # NON GENERIC
    alamo_settings = {'almname': "cam6",
                'monomialpower':(1, 2, 3, 4, 5, 6),
                'multi2power':(1, 2),
                'simulator':almsim,
                'maxiter':20,
                'cvfun':True}

    res = alamopy.doalamo(x, z, lmo = 3)
