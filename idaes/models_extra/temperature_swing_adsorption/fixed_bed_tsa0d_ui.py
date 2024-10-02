#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################

"""
UI exports for 0D Fixed Bed TSA unit model.
"""
from .fixed_bed_tsa0d import FixedBedTSA0D as Model
# from idaes.core.ui import fsapi as w -- the world according to dang
from watertap.ui import fsapi as w


def export_to_ui() -> w.FlowsheetInterface:
    # TO-DO: The wrapper
    return None