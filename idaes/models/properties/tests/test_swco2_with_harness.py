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
import idaes.models.properties.swco2 as swco2
from idaes.models.properties.tests.test_harness import PropertyTestHarness
from idaes.models.properties.general_helmholtz import helmholtz_available


@pytest.mark.skipif(not helmholtz_available(), reason="General Helmholtz not available")
class TestBasicMix(PropertyTestHarness):
    def configure(self):
        self.prop_pack = swco2.SWCO2ParameterBlock
        self.param_args = {"phase_presentation": swco2.PhaseType.MIX}
        self.prop_args = {}
        self.has_density_terms = True
        # Helmholtz package initialization has no solver calls, so can't fail
        self.skip_initialization_raises_exception_test = True


@pytest.mark.skipif(not helmholtz_available(), reason="General Helmholtz not available")
class TestBasicLV(PropertyTestHarness):
    def configure(self):
        self.prop_pack = swco2.SWCO2ParameterBlock
        self.param_args = {"phase_presentation": swco2.PhaseType.LG}
        self.prop_args = {}
        self.has_density_terms = True
        # Helmholtz package initialization has no solver calls, so can't fail
        self.skip_initialization_raises_exception_test = True


@pytest.mark.skipif(not helmholtz_available(), reason="General Helmholtz not available")
class TestBasicL(PropertyTestHarness):
    def configure(self):
        self.prop_pack = swco2.SWCO2ParameterBlock
        self.param_args = {"phase_presentation": swco2.PhaseType.L}
        self.prop_args = {}
        self.has_density_terms = True
        # Helmholtz package initialization has no solver calls, so can't fail
        self.skip_initialization_raises_exception_test = True


@pytest.mark.skipif(not helmholtz_available(), reason="General Helmholtz not available")
class TestBasicV(PropertyTestHarness):
    def configure(self):
        self.prop_pack = swco2.SWCO2ParameterBlock
        self.param_args = {"phase_presentation": swco2.PhaseType.G}
        self.prop_args = {}
        self.has_density_terms = True
        # Helmholtz package initialization has no solver calls, so can't fail
        self.skip_initialization_raises_exception_test = True
