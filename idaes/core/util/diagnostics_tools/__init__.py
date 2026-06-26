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
Centralized imports for model diagnostics tools.
"""

from .bounds import (
    get_valid_range_of_component,
    list_components_with_values_outside_valid_range,
    set_bounds_from_valid_range,
)
from .constraint_term_analysis import (
    ConstraintTermAnalysisVisitor,
)
from .degeneracy_hunter import (
    DegeneracyHunter,
    DHCONFIG,
)
from .ill_conditioning import (
    compute_ill_conditioning_certificate,
)
from .ipopt_halt_on_error import (
    ipopt_solve_halt_on_error,
)
from .svd_toolbox import (
    SVDToolbox,
    SVDCONFIG,
    svd_dense,
    svd_sparse,
)

# These imports need to go last to avoid circular imports
from .diagnostics_toolbox import (
    DiagnosticsToolbox,
)
from .convergence_analysis import (
    CACONFIG,
    IpoptConvergenceAnalysis,
)
