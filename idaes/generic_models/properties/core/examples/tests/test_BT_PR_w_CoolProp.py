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
"""
Author: Andrew Lee
"""

import pytest

try:
    from CoolProp.CoolProp import PropsSI
    coolprop_present = True
except ModuleNotFoundError:
    coolprop_present = False

from idaes.core import FlowsheetBlock
from idaes.generic_models.properties.core.eos.ceos import \
    cubic_roots_available
from idaes.generic_models.properties.core.examples.BT_PR_w_CoolProp import \
    configuration
from idaes.generic_models.properties.core.generic.generic_property import (
        GenericParameterBlock)
from pyomo.util.check_units import assert_units_consistent

from pyomo.environ import (ConcreteModel, Param, value, Var)

from idaes.core.util import get_solver

import idaes.logger as idaeslog
SOUT = idaeslog.INFO

# Set module level pyest marker
pytestmark = pytest.mark.cubic_root
prop_available = cubic_roots_available()


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


# -----------------------------------------------------------------------------
# Test loading parameters
@pytest.mark.skipif(not coolprop_present, reason="CoolProp not installed")
class TestCoolPropParameters(object):
    @pytest.fixture()
    def m(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(default={'dynamic': False})

        m.fs.props = GenericParameterBlock(default=configuration)

        return m

    @pytest.mark.unit
    def test_parameter_block(self, m):
        # Benzene parameters
        assert isinstance(m.fs.props.benzene.temperature_crit, Var)
        assert m.fs.props.benzene.temperature_crit.fixed
        assert value(m.fs.props.benzene.temperature_crit) == PropsSI(
            "TCRIT", "Benzene")

        assert isinstance(m.fs.props.benzene.pressure_crit, Var)
        assert m.fs.props.benzene.pressure_crit.fixed
        assert value(m.fs.props.benzene.pressure_crit) == PropsSI(
            "PCRIT", "Benzene")

        assert isinstance(m.fs.props.benzene.mw, Param)
        assert value(m.fs.props.benzene.mw) == PropsSI(
            "molarmass", "Benzene")

        assert isinstance(m.fs.props.benzene.omega, Var)
        assert m.fs.props.benzene.omega.fixed
        assert value(m.fs.props.benzene.omega) == PropsSI(
            "acentric", "Benzene")

        # Toluene parameters
        assert isinstance(m.fs.props.toluene.temperature_crit, Var)
        assert m.fs.props.toluene.temperature_crit.fixed
        assert value(m.fs.props.toluene.temperature_crit) == PropsSI(
            "TCRIT", "toluene")

        assert isinstance(m.fs.props.toluene.pressure_crit, Var)
        assert m.fs.props.toluene.pressure_crit.fixed
        assert value(m.fs.props.toluene.pressure_crit) == PropsSI(
            "PCRIT", "toluene")

        assert isinstance(m.fs.props.toluene.mw, Param)
        assert value(m.fs.props.toluene.mw) == PropsSI(
            "molarmass", "toluene")

        assert isinstance(m.fs.props.toluene.omega, Var)
        assert m.fs.props.toluene.omega.fixed
        assert value(m.fs.props.toluene.omega) == PropsSI(
            "acentric", "toluene")
