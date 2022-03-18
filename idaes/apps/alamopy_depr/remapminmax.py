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


def remapminmax(x, xmax, xmin):
    import numpy as np

    try:
        ndata, ninputs = np.shape(x)
        y = np.ones([ndata, ninputs])
    except Exception:
        ndata = np.size(x)
        ninputs = 1
        y = np.ones([ndata])
    if ninputs > 1:
        for j in range(ninputs):
            y[:, j] = (xmax[j] - xmin[j]) * (x[:, j] + 1) / (2.0) + xmin[j]
    else:
        y = (xmax - xmin) * (x + 1) / (2.0) + xmin
    return y
