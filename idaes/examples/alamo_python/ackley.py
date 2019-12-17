##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
# This file applies the ALAMOpython module to the ackley problem
# More information on this problemc an be found at
# https://www.sfu.ca/~ssurjano/ackley.html
# The function is corrupted by normally distributed noise of varaince 0.1


from idaes.surrogate import alamopy
from idaes.examples.alamo_python import examples

#Import additional python modules for creating the synthetic data
import math
import numpy as np

def main():
    # Specify number of poitns to be used in the training set
    # Validation data can be provided optionally
    ndata= 50
    nval = 500
    lb = [-15,-15]
    ub = [15,15]
    x = np.random.uniform(lb,ub,(ndata,2))
    xval = np.random.uniform(lb,ub,(nval,2))
    z = [0]*ndata
    zval = [0]*nval
    # specify simulator as examples.sixcamel
    sim = examples.ackley
    for i in range(ndata):
        z[i]=sim(x[i][0],x[i][1])
    for i in range(nval):
        zval[i] = sim(xval[i][0],xval[i][1])

    # Use alamopy's python function wrapper to avoid using ALAMO's I/O format
    almsim = alamopy.wrapwriter(sim)
    # Call alamo through the alamopy wrapper
    res = alamopy.doalamo(x,z,xval=xval,zval=zval,almname='ackley',monomialpower=(1,2,3,4,5,6),expfcns=1,multi2power=(1,2),expandoutput=True)
    # Calculate confidence intervals
    conf_inv = alamopy.almconfidence(res)

    print('Model: {}'.format(res['model']))
    print('Confidence Intervals : {}'.format(conf_inv['conf_inv']))

if __name__ == '__main__':
    np.random.seed(100)
    main()
