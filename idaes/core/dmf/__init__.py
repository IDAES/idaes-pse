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
IDAES Data Management Framework (DMF)

The DMF lets you save, search, and retrieve provenance related
to your models.
"""
__author__ = "Dan Gunter"

import logging

from .dmfbase import DMF, DMFConfig  # noqa: F401
from .dmfbase import create_configuration  # noqa: F401
from .getver import get_version_info  # noqa: F401
from .userapi import get_workspace  # noqa: F401
from . import resource  # noqa: F401

# DMF version is the same as IDAES version
from idaes import __version__  # noqa


# default log format
h = logging.StreamHandler()
h.setFormatter(
    logging.Formatter("%(asctime)s [%(levelname)s] " "%(name)s: %(message)s")
)
logging.getLogger("idaes.core.dmf").addHandler(h)
logging.getLogger("idaes.core.dmf").propagate = False
