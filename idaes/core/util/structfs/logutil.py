#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Utility functions for logging
"""

import logging
import warnings

g_quiet = {}


def quiet(roots=("idaes", "pyomo"), level=logging.CRITICAL):
    """Be very quiet. I'm hunting wabbits.

    Ignore warnings and set all loggers starting with one of
    the values in 'roots' to the given level (default=CRITICAL).
    """
    warnings.filterwarnings("ignore")
    all_loggers = [logging.getLogger()] + [
        logging.getLogger(name) for name in logging.root.manager.loggerDict
    ]
    for lg in all_loggers:
        for root in roots:
            if lg.name.startswith(root + "."):
                g_quiet[lg.name] = lg.level
                lg.setLevel(level)


def unquiet():
    """Reverse previous quiet()"""
    for k in list(g_quiet.keys()):
        v = g_quiet[k]
        lg = logging.getLogger(k)
        lg.setLevel(v)
        del g_quiet[k]
