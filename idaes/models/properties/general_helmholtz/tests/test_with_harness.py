#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################

__author__ = "John Eslick"

import pytest
from idaes.models.properties.general_helmholtz import (
    HelmholtzParameterBlock,
    PhaseType,
    StateVars,
    AmountBasis,
)
from idaes.models.properties.tests.test_harness import PropertyTestHarness
from idaes.models.properties.general_helmholtz import helmholtz_available


@pytest.mark.skipif(not helmholtz_available(), reason="General Helmholtz not available")
class TestBasicMix(PropertyTestHarness):
    def configure(self):
        self.prop_pack = HelmholtzParameterBlock
        self.param_args = {"pure_component": "h2o", "phase_presentation": PhaseType.MIX}
        self.prop_args = {}
        self.has_density_terms = True
        # Helmholtz package initialization has no solver calls, so can't fail
        self.skip_initialization_raises_exception_test = True


@pytest.mark.skipif(not helmholtz_available(), reason="General Helmholtz not available")
class TestBasicLV(PropertyTestHarness):
    def configure(self):
        self.prop_pack = HelmholtzParameterBlock
        self.param_args = {"pure_component": "h2o", "phase_presentation": PhaseType.LG}
        self.prop_args = {}
        self.has_density_terms = True
        # Helmholtz package initialization has no solver calls, so can't fail
        self.skip_initialization_raises_exception_test = True


@pytest.mark.skipif(not helmholtz_available(), reason="General Helmholtz not available")
class TestBasicL(PropertyTestHarness):
    def configure(self):
        self.prop_pack = HelmholtzParameterBlock
        self.param_args = {"pure_component": "h2o", "phase_presentation": PhaseType.L}
        self.prop_args = {}
        self.has_density_terms = True
        # Helmholtz package initialization has no solver calls, so can't fail
        self.skip_initialization_raises_exception_test = True


@pytest.mark.skipif(not helmholtz_available(), reason="General Helmholtz not available")
class TestBasicV(PropertyTestHarness):
    def configure(self):
        self.prop_pack = HelmholtzParameterBlock
        self.param_args = {"pure_component": "h2o", "phase_presentation": PhaseType.G}
        self.prop_args = {}
        self.has_density_terms = True
        # Helmholtz package initialization has no solver calls, so can't fail
        self.skip_initialization_raises_exception_test = True


@pytest.mark.skipif(not helmholtz_available(), reason="General Helmholtz not available")
class TestBasicMixSP(PropertyTestHarness):
    def configure(self):
        self.prop_pack = HelmholtzParameterBlock
        self.param_args = {
            "pure_component": "h2o",
            "phase_presentation": PhaseType.MIX,
            "state_vars": StateVars.PS,
        }
        self.prop_args = {}
        self.has_density_terms = True
        # Helmholtz package initialization has no solver calls, so can't fail
        self.skip_initialization_raises_exception_test = True


@pytest.mark.skipif(not helmholtz_available(), reason="General Helmholtz not available")
class TestBasicLVSP(PropertyTestHarness):
    def configure(self):
        self.prop_pack = HelmholtzParameterBlock
        self.param_args = {
            "pure_component": "h2o",
            "phase_presentation": PhaseType.LG,
            "state_vars": StateVars.PS,
        }
        self.prop_args = {}
        self.has_density_terms = True
        # Helmholtz package initialization has no solver calls, so can't fail
        self.skip_initialization_raises_exception_test = True


@pytest.mark.skipif(not helmholtz_available(), reason="General Helmholtz not available")
class TestBasicLSP(PropertyTestHarness):
    def configure(self):
        self.prop_pack = HelmholtzParameterBlock
        self.param_args = {
            "pure_component": "h2o",
            "phase_presentation": PhaseType.L,
            "state_vars": StateVars.PS,
        }
        self.prop_args = {}
        self.has_density_terms = True
        # Helmholtz package initialization has no solver calls, so can't fail
        self.skip_initialization_raises_exception_test = True


@pytest.mark.skipif(not helmholtz_available(), reason="General Helmholtz not available")
class TestBasicVSP(PropertyTestHarness):
    def configure(self):
        self.prop_pack = HelmholtzParameterBlock
        self.param_args = {
            "pure_component": "h2o",
            "phase_presentation": PhaseType.G,
            "state_vars": StateVars.PS,
        }
        self.prop_args = {}
        self.has_density_terms = True
        # Helmholtz package initialization has no solver calls, so can't fail
        self.skip_initialization_raises_exception_test = True


@pytest.mark.skipif(not helmholtz_available(), reason="General Helmholtz not available")
class TestBasicMixUP(PropertyTestHarness):
    def configure(self):
        self.prop_pack = HelmholtzParameterBlock
        self.param_args = {
            "pure_component": "h2o",
            "phase_presentation": PhaseType.MIX,
            "state_vars": StateVars.PU,
        }
        self.prop_args = {}
        self.has_density_terms = True
        # Helmholtz package initialization has no solver calls, so can't fail
        self.skip_initialization_raises_exception_test = True


@pytest.mark.skipif(not helmholtz_available(), reason="General Helmholtz not available")
class TestBasicLVUP(PropertyTestHarness):
    def configure(self):
        self.prop_pack = HelmholtzParameterBlock
        self.param_args = {
            "pure_component": "h2o",
            "phase_presentation": PhaseType.LG,
            "state_vars": StateVars.PU,
        }
        self.prop_args = {}
        self.has_density_terms = True
        # Helmholtz package initialization has no solver calls, so can't fail
        self.skip_initialization_raises_exception_test = True


@pytest.mark.skipif(not helmholtz_available(), reason="General Helmholtz not available")
class TestBasicLUP(PropertyTestHarness):
    def configure(self):
        self.prop_pack = HelmholtzParameterBlock
        self.param_args = {
            "pure_component": "h2o",
            "phase_presentation": PhaseType.L,
            "state_vars": StateVars.PU,
        }
        self.prop_args = {}
        self.has_density_terms = True
        # Helmholtz package initialization has no solver calls, so can't fail
        self.skip_initialization_raises_exception_test = True


@pytest.mark.skipif(not helmholtz_available(), reason="General Helmholtz not available")
class TestBasicVUP(PropertyTestHarness):
    def configure(self):
        self.prop_pack = HelmholtzParameterBlock
        self.param_args = {
            "pure_component": "h2o",
            "phase_presentation": PhaseType.G,
            "state_vars": StateVars.PU,
        }
        self.prop_args = {}
        self.has_density_terms = True
        # Helmholtz package initialization has no solver calls, so can't fail
        self.skip_initialization_raises_exception_test = True


@pytest.mark.skipif(not helmholtz_available(), reason="General Helmholtz not available")
class TestBasicMixMass(PropertyTestHarness):
    def configure(self):
        self.prop_pack = HelmholtzParameterBlock
        self.param_args = {
            "pure_component": "h2o",
            "phase_presentation": PhaseType.MIX,
            "amount_basis": AmountBasis.MASS,
        }
        self.prop_args = {}
        self.has_density_terms = True
        # Helmholtz package initialization has no solver calls, so can't fail
        self.skip_initialization_raises_exception_test = True


@pytest.mark.skipif(not helmholtz_available(), reason="General Helmholtz not available")
class TestBasicLVMass(PropertyTestHarness):
    def configure(self):
        self.prop_pack = HelmholtzParameterBlock
        self.param_args = {
            "pure_component": "h2o",
            "phase_presentation": PhaseType.LG,
            "amount_basis": AmountBasis.MASS,
        }
        self.prop_args = {}
        self.has_density_terms = True
        # Helmholtz package initialization has no solver calls, so can't fail
        self.skip_initialization_raises_exception_test = True


@pytest.mark.skipif(not helmholtz_available(), reason="General Helmholtz not available")
class TestBasicLMass(PropertyTestHarness):
    def configure(self):
        self.prop_pack = HelmholtzParameterBlock
        self.param_args = {
            "pure_component": "h2o",
            "phase_presentation": PhaseType.L,
            "amount_basis": AmountBasis.MASS,
        }
        self.prop_args = {}
        self.has_density_terms = True
        # Helmholtz package initialization has no solver calls, so can't fail
        self.skip_initialization_raises_exception_test = True


@pytest.mark.skipif(not helmholtz_available(), reason="General Helmholtz not available")
class TestBasicVMass(PropertyTestHarness):
    def configure(self):
        self.prop_pack = HelmholtzParameterBlock
        self.param_args = {
            "pure_component": "h2o",
            "phase_presentation": PhaseType.G,
            "amount_basis": AmountBasis.MASS,
        }
        self.prop_args = {}
        self.has_density_terms = True
        # Helmholtz package initialization has no solver calls, so can't fail
        self.skip_initialization_raises_exception_test = True


@pytest.mark.skipif(not helmholtz_available(), reason="General Helmholtz not available")
class TestBasicMixSPMass(PropertyTestHarness):
    def configure(self):
        self.prop_pack = HelmholtzParameterBlock
        self.param_args = {
            "pure_component": "h2o",
            "phase_presentation": PhaseType.MIX,
            "state_vars": StateVars.PS,
            "amount_basis": AmountBasis.MASS,
        }
        self.prop_args = {}
        self.has_density_terms = True
        # Helmholtz package initialization has no solver calls, so can't fail
        self.skip_initialization_raises_exception_test = True


@pytest.mark.skipif(not helmholtz_available(), reason="General Helmholtz not available")
class TestBasicLVSPMass(PropertyTestHarness):
    def configure(self):
        self.prop_pack = HelmholtzParameterBlock
        self.param_args = {
            "pure_component": "h2o",
            "phase_presentation": PhaseType.LG,
            "state_vars": StateVars.PS,
            "amount_basis": AmountBasis.MASS,
        }
        self.prop_args = {}
        self.has_density_terms = True
        # Helmholtz package initialization has no solver calls, so can't fail
        self.skip_initialization_raises_exception_test = True


@pytest.mark.skipif(not helmholtz_available(), reason="General Helmholtz not available")
class TestBasicLSPMass(PropertyTestHarness):
    def configure(self):
        self.prop_pack = HelmholtzParameterBlock
        self.param_args = {
            "pure_component": "h2o",
            "phase_presentation": PhaseType.L,
            "state_vars": StateVars.PS,
            "amount_basis": AmountBasis.MASS,
        }
        self.prop_args = {}
        self.has_density_terms = True
        # Helmholtz package initialization has no solver calls, so can't fail
        self.skip_initialization_raises_exception_test = True


@pytest.mark.skipif(not helmholtz_available(), reason="General Helmholtz not available")
class TestBasicVSPMass(PropertyTestHarness):
    def configure(self):
        self.prop_pack = HelmholtzParameterBlock
        self.param_args = {
            "pure_component": "h2o",
            "phase_presentation": PhaseType.G,
            "state_vars": StateVars.PS,
            "amount_basis": AmountBasis.MASS,
        }
        self.prop_args = {}
        self.has_density_terms = True
        # Helmholtz package initialization has no solver calls, so can't fail
        self.skip_initialization_raises_exception_test = True


@pytest.mark.skipif(not helmholtz_available(), reason="General Helmholtz not available")
class TestBasicMixUPMass(PropertyTestHarness):
    def configure(self):
        self.prop_pack = HelmholtzParameterBlock
        self.param_args = {
            "pure_component": "h2o",
            "phase_presentation": PhaseType.MIX,
            "state_vars": StateVars.PU,
            "amount_basis": AmountBasis.MASS,
        }
        self.prop_args = {}
        self.has_density_terms = True
        # Helmholtz package initialization has no solver calls, so can't fail
        self.skip_initialization_raises_exception_test = True


@pytest.mark.skipif(not helmholtz_available(), reason="General Helmholtz not available")
class TestBasicLVUPMass(PropertyTestHarness):
    def configure(self):
        self.prop_pack = HelmholtzParameterBlock
        self.param_args = {
            "pure_component": "h2o",
            "phase_presentation": PhaseType.LG,
            "state_vars": StateVars.PU,
            "amount_basis": AmountBasis.MASS,
        }
        self.prop_args = {}
        self.has_density_terms = True
        # Helmholtz package initialization has no solver calls, so can't fail
        self.skip_initialization_raises_exception_test = True


@pytest.mark.skipif(not helmholtz_available(), reason="General Helmholtz not available")
class TestBasicLUPMass(PropertyTestHarness):
    def configure(self):
        self.prop_pack = HelmholtzParameterBlock
        self.param_args = {
            "pure_component": "h2o",
            "phase_presentation": PhaseType.L,
            "state_vars": StateVars.PU,
            "amount_basis": AmountBasis.MASS,
        }
        self.prop_args = {}
        self.has_density_terms = True
        # Helmholtz package initialization has no solver calls, so can't fail
        self.skip_initialization_raises_exception_test = True


@pytest.mark.skipif(not helmholtz_available(), reason="General Helmholtz not available")
class TestBasicVUPMass(PropertyTestHarness):
    def configure(self):
        self.prop_pack = HelmholtzParameterBlock
        self.param_args = {
            "pure_component": "h2o",
            "phase_presentation": PhaseType.G,
            "state_vars": StateVars.PU,
            "amount_basis": AmountBasis.MASS,
        }
        self.prop_args = {}
        self.has_density_terms = True
        # Helmholtz package initialization has no solver calls, so can't fail
        self.skip_initialization_raises_exception_test = True


class TestBasicMixTPXMass(PropertyTestHarness):
    def configure(self):
        self.prop_pack = HelmholtzParameterBlock
        self.param_args = {
            "pure_component": "h2o",
            "phase_presentation": PhaseType.MIX,
            "state_vars": StateVars.TPX,
            "amount_basis": AmountBasis.MASS,
        }
        self.prop_args = {}
        self.has_density_terms = True
        # Helmholtz package initialization has no solver calls, so can't fail
        self.skip_initialization_raises_exception_test = True


class TestBasicLVTPXMass(PropertyTestHarness):
    def configure(self):
        self.prop_pack = HelmholtzParameterBlock
        self.param_args = {
            "pure_component": "h2o",
            "phase_presentation": PhaseType.LG,
            "state_vars": StateVars.TPX,
            "amount_basis": AmountBasis.MASS,
        }
        self.prop_args = {}
        self.has_density_terms = True
        # Helmholtz package initialization has no solver calls, so can't fail
        self.skip_initialization_raises_exception_test = True


class TestBasicLTPXMass(PropertyTestHarness):
    def configure(self):
        self.prop_pack = HelmholtzParameterBlock
        self.param_args = {
            "pure_component": "h2o",
            "phase_presentation": PhaseType.L,
            "state_vars": StateVars.TPX,
            "amount_basis": AmountBasis.MASS,
        }
        self.prop_args = {}
        self.has_density_terms = True
        # Helmholtz package initialization has no solver calls, so can't fail
        self.skip_initialization_raises_exception_test = True


class TestBasicVTPXMass(PropertyTestHarness):
    def configure(self):
        self.prop_pack = HelmholtzParameterBlock
        self.param_args = {
            "pure_component": "h2o",
            "phase_presentation": PhaseType.G,
            "state_vars": StateVars.TPX,
            "amount_basis": AmountBasis.MASS,
        }
        self.prop_args = {}
        self.has_density_terms = True
        # Helmholtz package initialization has no solver calls, so can't fail
        self.skip_initialization_raises_exception_test = True
