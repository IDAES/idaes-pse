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
"""
Alamopy tests with sixhumpcamel examples
"""
import pytest
from idaes.surrogate import alamopy
from idaes.surrogate.alamopy import alamo, almconfidence, almplot, wrapwriter
from idaes.surrogate.alamopy.multos import deletefile
import numpy as np

has_alamo_flag = alamopy.multos.has_alamo()


def sixcamel(*x):
    x1, x2 = x
    t1 = np.multiply(
        4.0 - 2.1 * np.power(x1, 2) + np.divide(np.power(x1, 4), 3.0), np.power(x1, 2)
    )
    t2 = np.multiply(4 * np.power(x2, 2) - 4, np.power(x2, 2))
    z = t1 + np.multiply(x1, x2) + t2
    return z


@pytest.mark.skipif(not has_alamo_flag, reason="alamo executable not found")
@pytest.mark.unit
def test_basic():

    if has_alamo_flag:
        ndata=10
        x = np.random.uniform([-2,-1],[2,1],(ndata,2))
        z = [0]*ndata
        # specify simulator as examples.sixcamel
        sim = sixcamel
        for i in range(ndata):
            z[i]=sim(x[i][0],x[i][1])

        # Use alamopy's python function wrapper to avoid using ALAMO's I/O format
        almsim = wrapwriter(sim)
        # Call alamo through the alamopy wrapper
        res = alamo(x, z, almname='cam6',monomialpower=(1,2,3,4,5,6), multi2power=(1,2), simulator=almsim, expandoutput=True,maxiter=20)#,cvfun=True)
        #conf_inv = almconfidence(res)

        #print('Model: {}'.format(res['model']))
        #print('Confidence Intervals : {}'.format(conf_inv['conf_inv']))
        almplot(res, show=False)
