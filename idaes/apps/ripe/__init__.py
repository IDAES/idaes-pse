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
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes".
"""
__all__ = [
    "ripemodel",
    "ems",
    "rspace",
    "sharedata",
    "debug",
    "powerlawp5",
    "powerlaw2",
    "powerlaw3",
    "powerlaw4",
    "avrami2",
    "avrami3",
    "avrami4",
    "avrami5",
    "randomnuc",
    "ptompkins",
    "jander",
    "antijander",
    "valensi",
    "parabolic",
    "gb3d",
    "zlt",
    "grain",
    # PYLINT-TODO-FIX: this seems to be a genuine error since "massact" is not imported from .mechs
    "massact",  # pylint: disable=undefined-all-variable
    "massactm",
    "getmechs",
]

from .main import ripemodel, ripewrite, print_results  # noqa: F401
from .shared import rspace, sharedata, debug  # noqa: F401
from .atermconstruct import (
    makeaterm,
    formatinputs,
    checkargs,
    normalizefeatures,
)  # noqa: F401
from .kinforms import lin, linjac, arr, arrjac, refarr, refarrjac  # noqa: F401
from .mechs import (
    powerlawp5,
    powerlaw2,
    powerlaw3,
    powerlaw4,
    avrami2,
    avrami3,
    avrami4,
    avrami5,
    randomnuc,
    ptompkins,
    jander,
    antijander,
    valensi,
    parabolic,
    gb3d,
    zlt,
    grain,
    getmechs,
    massactm,
)  # noqa: F401
from .genpyomo import ripeomo  # noqa: F401
from .targets import (
    doalamo,
    dopwalamo,
    gentargets,
    sstargets,
    dynamictargets,
)  # noqa: F401
from .confinv import confinv  # noqa: F401
from .emsampling import constructmodel, ems  # noqa: F401
from .checkoptions import checkoptions  # noqa: F401
from .bounds import stoich_cons, count_neg, get_bounds  # noqa: F401
