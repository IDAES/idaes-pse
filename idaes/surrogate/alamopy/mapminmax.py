###############################################################################
# Copyright
# =========
#
# Institute for the Design of Advanced Energy Systems Process Systems Engineering
# Framework (IDAES PSE Framework) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021 by the
# software owners: The Regents of the University of California, through Lawrence
# Berkeley National Laboratory,  National Technology & Engineering Solutions of
# Sandia, LLC, Carnegie Mellon University, West Virginia University Research
# Corporation, et al.  All rights reserved.
#
# NOTICE.  This Software was developed under funding from the U.S. Department of
# Energy and the U.S. Government consequently retains certain rights. As such, the
# U.S. Government has been granted for itself and others acting on its behalf a
# paid-up, nonexclusive, irrevocable, worldwide license in the Software to
# reproduce, distribute copies to the public, prepare derivative works, and
# perform publicly and display publicly, and to permit other to do so. Copyright
# (C) 2018-2019 IDAES - All Rights Reserved
#
###############################################################################


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
