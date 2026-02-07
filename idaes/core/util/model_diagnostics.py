# -*- coding: utf-8 -*-
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
This module was the original container for the diagnostics tools.

As these have grown significantly, the code has been broken out into sub-modules in the
diagnostics_tools subpackage. This module is retained for backward compatibility,
but all tools are now imported from the diagnostics_tools subpackage.
"""

__author__ = "Alexander Dowling, Douglas Allan, Andrew Lee, Robby Parker, Ben Knueven"


# Imports for backward compatibility
from idaes.core.util.diagnostics_tools import (
    ConstraintTermAnalysisVisitor,
    compute_ill_conditioning_certificate,
    # TODO: Rename and redirect once old DegeneracyHunter is removed
    DegeneracyHunter as DegeneracyHunter2,
    DHCONFIG,
    DiagnosticsToolbox,
    # IpoptConvergenceAnalysis,
    SVDToolbox,
    SVDCONFIG,
    get_valid_range_of_component,
    list_components_with_values_outside_valid_range,
    set_bounds_from_valid_range,
    ipopt_halt_on_error,
)
from idaes.core.util.diagnostics_tools.deprecated.degeneracy_hunter_legacy import DegeneracyHunter
