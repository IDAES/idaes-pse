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
import idaes.generic_models.properties.swco2 as swco2
from idaes.generic_models.properties.tests.test_harness import \
    PropertyTestHarness
import pyomo.environ as pyo


@pytest.mark.unit
class TestBasicMix(PropertyTestHarness):
    def configure(self):
        self.prop_pack = swco2.SWCO2ParameterBlock
        self.param_args = {"phase_presentation": swco2.PhaseType.MIX}
        self.prop_args = {}
        self.has_density_terms = True


@pytest.mark.unit
class TestBasicLV(PropertyTestHarness):
    def configure(self):
        self.prop_pack = swco2.SWCO2ParameterBlock
        self.param_args = {"phase_presentation": swco2.PhaseType.LG}
        self.prop_args = {}
        self.has_density_terms = True


@pytest.mark.unit
class TestBasicL(PropertyTestHarness):
    def configure(self):
        self.prop_pack = swco2.SWCO2ParameterBlock
        self.param_args = {"phase_presentation": swco2.PhaseType.L}
        self.prop_args = {}
        self.has_density_terms = True


@pytest.mark.unit
class TestBasicV(PropertyTestHarness):
    def configure(self):
        self.prop_pack = swco2.SWCO2ParameterBlock
        self.param_args = {"phase_presentation": swco2.PhaseType.G}
        self.prop_args = {}
        self.has_density_terms = True
