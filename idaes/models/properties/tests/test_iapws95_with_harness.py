#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################

__author__ = "John Eslick"

import pytest
import idaes.models.properties.iapws95 as iapws95
from idaes.models.properties.tests.test_harness import PropertyTestHarness


class TestBasicMix(PropertyTestHarness):
    def configure(self):
        self.prop_pack = iapws95.Iapws95ParameterBlock
        self.param_args = {"phase_presentation": iapws95.PhaseType.MIX}
        self.prop_args = {}
        self.has_density_terms = True
        # Helmholtz package initialization has no solver calls, so can't fail
        self.skip_initialization_raises_exception_test = True


class TestBasicLV(PropertyTestHarness):
    def configure(self):
        self.prop_pack = iapws95.Iapws95ParameterBlock
        self.param_args = {"phase_presentation": iapws95.PhaseType.LG}
        self.prop_args = {}
        self.has_density_terms = True
        # Helmholtz package initialization has no solver calls, so can't fail
        self.skip_initialization_raises_exception_test = True


class TestBasicL(PropertyTestHarness):
    def configure(self):
        self.prop_pack = iapws95.Iapws95ParameterBlock
        self.param_args = {"phase_presentation": iapws95.PhaseType.L}
        self.prop_args = {}
        self.has_density_terms = True
        # Helmholtz package initialization has no solver calls, so can't fail
        self.skip_initialization_raises_exception_test = True


class TestBasicV(PropertyTestHarness):
    def configure(self):
        self.prop_pack = iapws95.Iapws95ParameterBlock
        self.param_args = {"phase_presentation": iapws95.PhaseType.G}
        self.prop_args = {}
        self.has_density_terms = True
        # Helmholtz package initialization has no solver calls, so can't fail
        self.skip_initialization_raises_exception_test = True
