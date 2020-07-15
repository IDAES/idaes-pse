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
