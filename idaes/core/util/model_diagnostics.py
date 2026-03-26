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


from pyomo.common.deprecation import relocated_module_attribute

relocated_module_attribute(
    "ConstraintTermAnalysisVisitor",
    "idaes.core.util.diagnostics_tools.constraint_term_analysis.ConstraintTermAnalysisVisitor",
    version="2.12.0",
)

relocated_module_attribute(
    "compute_ill_conditioning_certificate",
    "idaes.core.util.diagnostics_tools.ill_conditioning.compute_ill_conditioning_certificate",
    version="2.12.0",
)

relocated_module_attribute(
    "DegeneracyHunter",
    "idaes.core.util.diagnostics_tools.deprecated.degeneracy_hunter_legacy.DegeneracyHunter",
    version="2.12.0",
)
relocated_module_attribute(
    "DegeneracyHunter2",
    "idaes.core.util.diagnostics_tools.degeneracy_hunter.DegeneracyHunter",
    version="2.12.0",
)
relocated_module_attribute(
    "DHCONFIG",
    "idaes.core.util.diagnostics_tools.degeneracy_hunter.DHCONFIG",
    version="2.12.0",
)

relocated_module_attribute(
    "DiagnosticsToolbox",
    "idaes.core.util.diagnostics_tools.diagnostics_toolbox.DiagnosticsToolbox",
    version="2.12.0",
)
relocated_module_attribute(
    "CONFIG",
    "idaes.core.util.diagnostics_tools.diagnostics_toolbox.CONFIG",
    version="2.12.0",
)


relocated_module_attribute(
    "IpoptConvergenceAnalysis",
    "idaes.core.util.diagnostics_tools.convergence_analysis.IpoptConvergenceAnalysis",
    version="2.12.0",
)
relocated_module_attribute(
    "psweep_runner_validator",
    "idaes.core.util.diagnostics_tools.convergence_analysis.psweep_runner_validator",
    version="2.12.0",
)
relocated_module_attribute(
    "CACONFIG",
    "idaes.core.util.diagnostics_tools.convergence_analysis.CACONFIG",
    version="2.12.0",
)

relocated_module_attribute(
    "ipopt_solve_halt_on_error",
    "idaes.core.util.diagnostics_tools.ipopt_halt_on_error.ipopt_solve_halt_on_error",
    version="2.12.0",
)

relocated_module_attribute(
    "get_valid_range_of_component",
    "idaes.core.util.diagnostics_tools.bounds.get_valid_range_of_component",
    version="2.12.0",
)
relocated_module_attribute(
    "list_components_with_values_outside_valid_range",
    "idaes.core.util.diagnostics_tools.bounds.list_components_with_values_outside_valid_range",
    version="2.12.0",
)
relocated_module_attribute(
    "set_bounds_from_valid_range",
    "idaes.core.util.diagnostics_tools.bounds.set_bounds_from_valid_range",
    version="2.12.0",
)

relocated_module_attribute(
    "SVDToolbox",
    "idaes.core.util.diagnostics_tools.svd_toolbox.SVDToolbox",
    version="2.12.0",
)
relocated_module_attribute(
    "SVDCONFIG",
    "idaes.core.util.diagnostics_tools.svd_toolbox.SVDCONFIG",
    version="2.12.0",
)
relocated_module_attribute(
    "svd_callback_validator",
    "idaes.core.util.diagnostics_tools.svd_toolbox.svd_callback_validator",
    version="2.12.0",
)
relocated_module_attribute(
    "svd_dense",
    "idaes.core.util.diagnostics_tools.svd_toolbox.svd_dense",
    version="2.12.0",
)
relocated_module_attribute(
    "svd_sparse",
    "idaes.core.util.diagnostics_tools.svd_toolbox.svd_sparse",
    version="2.12.0",
)

relocated_module_attribute(
    "check_parallel_jacobian",
    "idaes.core.util.diagnostics_tools.utils.check_parallel_jacobian",
    version="2.12.0",
)
