##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################


def mapminmax(x, *args):
    import numpy as np
    if args == ():
        xmax = np.amax(x, 0)
        xmin = np.amin(x, 0)
    else:
        xmax = args[0]
        xmin = args[1]
    try:
        ndata, ninputs = np.shape(x)
        y = np.ones([ndata, ninputs])
    except Exception:
        ndata = np.size(x)
        ninputs = 1
        y = np.ones([ndata])
    if ninputs > 1:
        # check for constant column
        for j in range(ninputs):
            if (xmin[j] != xmax[j]):
                y[:, j] = 2 * (x[:, j] - xmin[j]) / float(xmax[j] - xmin[j]) - 1
            else:
                y[:, j] = 1.0
    else:
        if (xmin != xmax):
            y = 2 * (x - xmin) / float(xmax - xmin) - 1
        else:
            y[:] = 1.0
    return y, xmax, xmin
