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
# This file applies the ALAMOpython module to the six hump camel problem
# More information on this function can be found at :
#            https://www.sfu.ca/~ssurjano/camel6.html
# This problem utilizes ALAMO's sampling features
from idaes.surrogate import alamopy
from idaes.examples.alamo_python import examples

#Import additional python modules for creating the synthetic data
import math
import numpy as np

def main():
    # Specify number of poitns to be used in the training set
    # Validation data can be provided optionally
    ndata=10
    x = np.random.uniform([-2,-1],[2,1],(ndata,2))
    z = [0]*ndata
    # specify simulator as examples.sixcamel
    sim = examples.sixcamel
    for i in range(ndata):
        z[i]=sim(x[i][0],x[i][1])

    # Use alamopy's python function wrapper to avoid using ALAMO's I/O format
    almsim = alamopy.wrapwriter(sim)
    # Call alamo through the alamopy wrapper
    res = alamopy.doalamo(x,z,almname='cam6',monomialpower=(1,2,3,4,5,6),multi2power=(1,2),simulator=almsim, expandoutput=True,maxiter=20,cvfun=True)
#    print res
    print('Model: {}'.format(res['model']))

if __name__ == '__main__':
    np.random.seed(100)
    main()
