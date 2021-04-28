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
"""
Tests for methods from electrolyte library

Authors: Andrew Lee
"""

import pytest
import types

from pyomo.environ import \
    ConcreteModel, Block, Expression, value, Var, units as pyunits
from pyomo.common.config import ConfigBlock
from pyomo.util.check_units import assert_units_equivalent

from idaes.generic_models.properties.core.pure.electrolyte import \
    relative_permittivity_constant
from idaes.core.util.misc import add_object_reference
from idaes.core.property_meta import PropertyClassMetadata


@pytest.fixture()
def frame():
    m = ConcreteModel()

    # Create a dummy parameter block
    m.params = Block()

    m.params.config = ConfigBlock(implicit=True)
    m.params.config.parameter_data = {
        "relative_permittivity_liq_comp": 101}

    m.meta_object = PropertyClassMetadata()
    m.meta_object.default_units["temperature"] = pyunits.K
    m.meta_object.default_units["mass"] = pyunits.kg
    m.meta_object.default_units["length"] = pyunits.m
    m.meta_object.default_units["time"] = pyunits.s
    m.meta_object.default_units["amount"] = pyunits.mol

    def get_metadata(self):
        return m.meta_object
    m.get_metadata = types.MethodType(get_metadata, m)
    m.params.get_metadata = types.MethodType(get_metadata, m.params)

    # Add necessary parameters to parameter block
    m.params.temperature_ref = Var(initialize=298.15, units=pyunits.K)
    m.params.pressure_ref = Var(initialize=1e5, units=pyunits.Pa)

    # Create a dummy state block
    m.props = Block([1])
    add_object_reference(m.props[1], "params", m.params)

    m.props[1].temperature = Var(initialize=500, units=pyunits.K)
    m.props[1].pressure = Var(initialize=101325, units=pyunits.Pa)

    return m


@pytest.mark.unit
def test_relative_permitivity_constant(frame):
    relative_permittivity_constant.build_parameters(frame.params)

    assert isinstance(frame.params.relative_permittivity_liq_comp, Var)
    assert value(frame.params.relative_permittivity_liq_comp) == 101

    expr = relative_permittivity_constant.return_expression(
        frame.props[1], frame.params, frame.props[1].temperature)
    assert value(expr) == 101

    frame.props[1].temperature.value = 600
    assert value(expr) == 101

    assert_units_equivalent(expr, pyunits.dimensionless)
