###############################################################################
# ** Copyright Notice **
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021 by the
# software owners: The Regents of the University of California, through Lawrence
# Berkeley National Laboratory,  National Technology & Engineering Solutions of
# Sandia, LLC, Carnegie Mellon University, West Virginia University Research
# Corporation, et al.  All rights reserved.
#
# NOTICE.  This Software was developed under funding from the U.S. Department of
# Energy and the U.S. Government consequently retains certain rights. As such, the
# U.S. Government has been granted for itself and others acting on its behalf a
# paid-up, nonexclusive, irrevocable, worldwide license in the Software to
# reproduce, distribute copies to the public, prepare derivative works, and
# perform publicly and display publicly, and to permit other to do so.
###############################################################################

__author__ = "John Eslick"

import pytest
import idaes.generic_models.properties.iapws95 as iapws95
from idaes.generic_models.properties.tests.test_harness import \
    PropertyTestHarness
import pyomo.environ as pyo


@pytest.mark.unit
class TestBasicMix(PropertyTestHarness):
    def configure(self):
        self.prop_pack = iapws95.Iapws95ParameterBlock
        self.param_args = {"phase_presentation": iapws95.PhaseType.MIX}
        self.prop_args = {}
        self.has_density_terms = True


@pytest.mark.unit
class TestBasicLV(PropertyTestHarness):
    def configure(self):
        self.prop_pack = iapws95.Iapws95ParameterBlock
        self.param_args = {"phase_presentation": iapws95.PhaseType.LG}
        self.prop_args = {}
        self.has_density_terms = True


@pytest.mark.unit
class TestBasicL(PropertyTestHarness):
    def configure(self):
        self.prop_pack = iapws95.Iapws95ParameterBlock
        self.param_args = {"phase_presentation": iapws95.PhaseType.L}
        self.prop_args = {}
        self.has_density_terms = True


@pytest.mark.unit
class TestBasicV(PropertyTestHarness):
    def configure(self):
        self.prop_pack = iapws95.Iapws95ParameterBlock
        self.param_args = {"phase_presentation": iapws95.PhaseType.G}
        self.prop_args = {}
        self.has_density_terms = True
