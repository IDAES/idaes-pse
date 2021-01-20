##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################

__author__ = "John Eslick"

import pytest
import idaes.generic_models.properties.iapws95 as iapws95
from idaes.generic_models.properties.tests.test_harness import \
    PropertyTestHarness
from idaes.generic_models.properties.iapws95 import \
    iapws95_available as prop_available
import pyomo.environ as pyo

if pyo.SolverFactory('ipopt').available():
    solver = pyo.SolverFactory('ipopt')
    solver.options = {'tol': 1e-6}
else:
    solver = None


@pytest.mark.unit
@pytest.mark.skipif(not prop_available(), reason="Property lib not available")
class TestBasicMix(PropertyTestHarness):
    def configure(self):
        self.prop_pack = iapws95.Iapws95ParameterBlock
        self.param_args = {"phase_presentation": iapws95.PhaseType.MIX}
        self.prop_args = {}
        self.has_density_terms = True


@pytest.mark.unit
@pytest.mark.skipif(not prop_available(), reason="Property lib not available")
class TestBasicLV(PropertyTestHarness):
    def configure(self):
        self.prop_pack = iapws95.Iapws95ParameterBlock
        self.param_args = {"phase_presentation": iapws95.PhaseType.LG}
        self.prop_args = {}
        self.has_density_terms = True


@pytest.mark.unit
@pytest.mark.skipif(not prop_available(), reason="Property lib not available")
class TestBasicL(PropertyTestHarness):
    def configure(self):
        self.prop_pack = iapws95.Iapws95ParameterBlock
        self.param_args = {"phase_presentation": iapws95.PhaseType.L}
        self.prop_args = {}
        self.has_density_terms = True


@pytest.mark.unit
@pytest.mark.skipif(not prop_available(), reason="Property lib not available")
class TestBasicV(PropertyTestHarness):
    def configure(self):
        self.prop_pack = iapws95.Iapws95ParameterBlock
        self.param_args = {"phase_presentation": iapws95.PhaseType.G}
        self.prop_args = {}
        self.has_density_terms = True
