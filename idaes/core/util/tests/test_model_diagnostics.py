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
Tests for backward compatibility of model_diagnostics module.

This module tests the deprecated imports that have been relocated to the
diagnostics_tools subpackage to ensure backward compatibility is maintained.
"""

import pytest


@pytest.mark.unit
class TestDeprecationPaths:
    """
    Test that all deprecated import paths in model_diagnostics.py still work.

    These imports should work but may emit deprecation warnings.
    """

    def test_constraint_term_analysis_visitor_import(self):
        """Test ConstraintTermAnalysisVisitor can be imported from old location."""
        from idaes.core.util.model_diagnostics import ConstraintTermAnalysisVisitor
        from idaes.core.util.diagnostics_tools.constraint_term_analysis import (
            ConstraintTermAnalysisVisitor as NewConstraintTermAnalysisVisitor,
        )

        assert ConstraintTermAnalysisVisitor is NewConstraintTermAnalysisVisitor

    def test_compute_ill_conditioning_certificate_import(self):
        """Test compute_ill_conditioning_certificate can be imported from old location."""
        from idaes.core.util.model_diagnostics import (
            compute_ill_conditioning_certificate,
        )
        from idaes.core.util.diagnostics_tools.ill_conditioning import (
            compute_ill_conditioning_certificate as new_compute_ill_conditioning_certificate,
        )

        assert (
            compute_ill_conditioning_certificate
            is new_compute_ill_conditioning_certificate
        )

    def test_degeneracy_hunter_import(self):
        """Test DegeneracyHunter can be imported from old location."""
        from idaes.core.util.model_diagnostics import DegeneracyHunter
        from idaes.core.util.diagnostics_tools.deprecated.degeneracy_hunter_legacy import (
            DegeneracyHunter as NewDegeneracyHunter,
        )

        assert DegeneracyHunter is NewDegeneracyHunter

    def test_degeneracy_hunter2_import(self):
        """Test DegeneracyHunter2 can be imported from old location."""
        from idaes.core.util.model_diagnostics import DegeneracyHunter2
        from idaes.core.util.diagnostics_tools.degeneracy_hunter import (
            DegeneracyHunter as NewDegeneracyHunter,
        )

        assert DegeneracyHunter2 is NewDegeneracyHunter

    def test_dhconfig_import(self):
        """Test DHCONFIG can be imported from old location."""
        from idaes.core.util.model_diagnostics import DHCONFIG
        from idaes.core.util.diagnostics_tools.degeneracy_hunter import (
            DHCONFIG as NewDHCONFIG,
        )

        assert DHCONFIG is NewDHCONFIG

    def test_diagnostics_toolbox_import(self):
        """Test DiagnosticsToolbox can be imported from old location."""
        from idaes.core.util.model_diagnostics import DiagnosticsToolbox
        from idaes.core.util.diagnostics_tools.diagnostics_toolbox import (
            DiagnosticsToolbox as NewDiagnosticsToolbox,
        )

        assert DiagnosticsToolbox is NewDiagnosticsToolbox

    def test_config_import(self):
        """Test CONFIG can be imported from old location."""
        from idaes.core.util.model_diagnostics import CONFIG
        from idaes.core.util.diagnostics_tools.diagnostics_toolbox import (
            CONFIG as NewCONFIG,
        )

        assert CONFIG is NewCONFIG

    def test_ipopt_convergence_analysis_import(self):
        """Test IpoptConvergenceAnalysis can be imported from old location."""
        from idaes.core.util.model_diagnostics import IpoptConvergenceAnalysis
        from idaes.core.util.diagnostics_tools.convergence_analysis import (
            IpoptConvergenceAnalysis as NewIpoptConvergenceAnalysis,
        )

        assert IpoptConvergenceAnalysis is NewIpoptConvergenceAnalysis

    def test_psweep_runner_validator_import(self):
        """Test psweep_runner_validator can be imported from old location."""
        from idaes.core.util.model_diagnostics import psweep_runner_validator
        from idaes.core.util.diagnostics_tools.convergence_analysis import (
            psweep_runner_validator as new_psweep_runner_validator,
        )

        assert psweep_runner_validator is new_psweep_runner_validator

    def test_caconfig_import(self):
        """Test CACONFIG can be imported from old location."""
        from idaes.core.util.model_diagnostics import CACONFIG
        from idaes.core.util.diagnostics_tools.convergence_analysis import (
            CACONFIG as NewCACONFIG,
        )

        assert CACONFIG is NewCACONFIG

    def test_ipopt_solve_halt_on_error_import(self):
        """Test ipopt_solve_halt_on_error can be imported from old location."""
        from idaes.core.util.model_diagnostics import ipopt_solve_halt_on_error
        from idaes.core.util.diagnostics_tools.ipopt_halt_on_error import (
            ipopt_solve_halt_on_error as new_ipopt_solve_halt_on_error,
        )

        assert ipopt_solve_halt_on_error is new_ipopt_solve_halt_on_error

    def test_get_valid_range_of_component_import(self):
        """Test get_valid_range_of_component can be imported from old location."""
        from idaes.core.util.model_diagnostics import get_valid_range_of_component
        from idaes.core.util.diagnostics_tools.bounds import (
            get_valid_range_of_component as new_get_valid_range_of_component,
        )

        assert get_valid_range_of_component is new_get_valid_range_of_component

    def test_list_components_with_values_outside_valid_range_import(self):
        """Test list_components_with_values_outside_valid_range can be imported from old location."""
        from idaes.core.util.model_diagnostics import (
            list_components_with_values_outside_valid_range,
        )
        from idaes.core.util.diagnostics_tools.bounds import (
            list_components_with_values_outside_valid_range as new_list_components_with_values_outside_valid_range,
        )

        assert (
            list_components_with_values_outside_valid_range
            is new_list_components_with_values_outside_valid_range
        )

    def test_set_bounds_from_valid_range_import(self):
        """Test set_bounds_from_valid_range can be imported from old location."""
        from idaes.core.util.model_diagnostics import set_bounds_from_valid_range
        from idaes.core.util.diagnostics_tools.bounds import (
            set_bounds_from_valid_range as new_set_bounds_from_valid_range,
        )

        assert set_bounds_from_valid_range is new_set_bounds_from_valid_range

    def test_svd_toolbox_import(self):
        """Test SVDToolbox can be imported from old location."""
        from idaes.core.util.model_diagnostics import SVDToolbox
        from idaes.core.util.diagnostics_tools.svd_toolbox import (
            SVDToolbox as NewSVDToolbox,
        )

        assert SVDToolbox is NewSVDToolbox

    def test_svdconfig_import(self):
        """Test SVDCONFIG can be imported from old location."""
        from idaes.core.util.model_diagnostics import SVDCONFIG
        from idaes.core.util.diagnostics_tools.svd_toolbox import (
            SVDCONFIG as NewSVDCONFIG,
        )

        assert SVDCONFIG is NewSVDCONFIG

    def test_svd_callback_validator_import(self):
        """Test svd_callback_validator can be imported from old location."""
        from idaes.core.util.model_diagnostics import svd_callback_validator
        from idaes.core.util.diagnostics_tools.svd_toolbox import (
            svd_callback_validator as new_svd_callback_validator,
        )

        assert svd_callback_validator is new_svd_callback_validator

    def test_svd_dense_import(self):
        """Test svd_dense can be imported from old location."""
        from idaes.core.util.model_diagnostics import svd_dense
        from idaes.core.util.diagnostics_tools.svd_toolbox import (
            svd_dense as new_svd_dense,
        )

        assert svd_dense is new_svd_dense

    def test_svd_sparse_import(self):
        """Test svd_sparse can be imported from old location."""
        from idaes.core.util.model_diagnostics import svd_sparse
        from idaes.core.util.diagnostics_tools.svd_toolbox import (
            svd_sparse as new_svd_sparse,
        )

        assert svd_sparse is new_svd_sparse

    def test_check_parallel_jacobian_import(self):
        """Test check_parallel_jacobian can be imported from old location."""
        from idaes.core.util.model_diagnostics import check_parallel_jacobian
        from idaes.core.util.diagnostics_tools.utils import (
            check_parallel_jacobian as new_check_parallel_jacobian,
        )

        assert check_parallel_jacobian is new_check_parallel_jacobian


@pytest.mark.unit
class TestAllDeprecatedImportsTogether:
    """
    Test importing all deprecated items in a single test to ensure no conflicts.
    """

    def test_all_imports_work(self):
        """Test that all deprecated imports can be imported together."""
        from idaes.core.util.model_diagnostics import (
            CACONFIG,
            CONFIG,
            DHCONFIG,
            SVDCONFIG,
            ConstraintTermAnalysisVisitor,
            DegeneracyHunter,
            DegeneracyHunter2,
            DiagnosticsToolbox,
            IpoptConvergenceAnalysis,
            SVDToolbox,
            check_parallel_jacobian,
            compute_ill_conditioning_certificate,
            get_valid_range_of_component,
            ipopt_solve_halt_on_error,
            list_components_with_values_outside_valid_range,
            psweep_runner_validator,
            set_bounds_from_valid_range,
            svd_callback_validator,
            svd_dense,
            svd_sparse,
        )

        # Verify all imports are not None
        assert ConstraintTermAnalysisVisitor is not None
        assert compute_ill_conditioning_certificate is not None
        assert DegeneracyHunter is not None
        assert DegeneracyHunter2 is not None
        assert DHCONFIG is not None
        assert DiagnosticsToolbox is not None
        assert CONFIG is not None
        assert IpoptConvergenceAnalysis is not None
        assert psweep_runner_validator is not None
        assert CACONFIG is not None
        assert ipopt_solve_halt_on_error is not None
        assert get_valid_range_of_component is not None
        assert list_components_with_values_outside_valid_range is not None
        assert set_bounds_from_valid_range is not None
        assert SVDToolbox is not None
        assert SVDCONFIG is not None
        assert svd_callback_validator is not None
        assert svd_dense is not None
        assert svd_sparse is not None
        assert check_parallel_jacobian is not None
