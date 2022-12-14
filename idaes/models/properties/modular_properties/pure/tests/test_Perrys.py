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
Tests for methods from Perry's

All methods and parameters from:

Perry's Chemical Engineers' Handbook, 7th Edition
Perry, Green, Maloney, 1997, McGraw-Hill

All parameter indicies based on conventions used by the source

Authors: Andrew Lee
"""

import pytest
import types

from pyomo.environ import ConcreteModel, Block, value, Var, Param, units as pyunits
from pyomo.common.config import ConfigBlock
from pyomo.util.check_units import assert_units_equivalent

from idaes.models.properties.modular_properties.pure.Perrys import *
from idaes.core.util.misc import add_object_reference
from idaes.core.base.property_meta import PropertyClassMetadata, UnitSet
from idaes.core.util.exceptions import ConfigurationError


@pytest.fixture()
def frame():
    m = ConcreteModel()

    # Create a dummy parameter block
    m.params = Block()

    m.params.config = ConfigBlock(implicit=True)
    m.params.config.parameter_data = {
        "dens_mol_liq_comp_coeff": {
            "eqn_type": 1,
            "1": 5.459,
            "2": 0.30542,
            "3": 647.13,
            "4": 0.081,
        },
        "cp_mol_liq_comp_coeff": {
            "1": 2.7637e05,
            "2": -2.0901e03,
            "3": 8.1250e00,
            "4": -1.4116e-2,
            "5": 9.3701e-06,
        },
        "enth_mol_form_liq_comp_ref": -285.83e3,
        "entr_mol_form_liq_comp_ref": 69.95,
    }
    m.params.config.include_enthalpy_of_formation = True

    # Also need to dummy configblock on the model for the test
    m.config = ConfigBlock(implicit=True)
    m.config.include_enthalpy_of_formation = True

    m.meta_object = PropertyClassMetadata()
    m.meta_object._default_units = UnitSet(
        temperature=pyunits.K,
        mass=pyunits.kg,
        length=pyunits.m,
        time=pyunits.s,
        amount=pyunits.mol,
    )

    def get_metadata(self):
        return m.meta_object

    m.get_metadata = types.MethodType(get_metadata, m)
    m.params.get_metadata = types.MethodType(get_metadata, m.params)

    # Add necessary parameters to parameter block
    m.params.temperature_ref = Var(initialize=273.16, units=pyunits.K)

    # Create a dummy state block
    m.props = Block([1])
    add_object_reference(m.props[1], "params", m.params)

    m.props[1].temperature = Var(initialize=273.16, units=pyunits.K)

    return m


@pytest.mark.unit
def test_cp_mol_liq_comp(frame):
    Perrys.cp_mol_liq_comp.build_parameters(frame.params)

    assert isinstance(frame.params.cp_mol_liq_comp_coeff_1, Var)
    assert value(frame.params.cp_mol_liq_comp_coeff_1) == 2.7637e05
    assert isinstance(frame.params.cp_mol_liq_comp_coeff_2, Var)
    assert value(frame.params.cp_mol_liq_comp_coeff_2) == -2.0901e03
    assert isinstance(frame.params.cp_mol_liq_comp_coeff_3, Var)
    assert value(frame.params.cp_mol_liq_comp_coeff_3) == 8.1250e00
    assert isinstance(frame.params.cp_mol_liq_comp_coeff_4, Var)
    assert value(frame.params.cp_mol_liq_comp_coeff_4) == -1.4116e-2
    assert isinstance(frame.params.cp_mol_liq_comp_coeff_5, Var)
    assert value(frame.params.cp_mol_liq_comp_coeff_5) == 9.3701e-06

    expr = Perrys.cp_mol_liq_comp.return_expression(
        frame.props[1], frame.params, frame.props[1].temperature
    )
    assert value(expr) == pytest.approx(76.150, rel=1e-3)

    frame.props[1].temperature.value = 533.15
    assert value(expr) == pytest.approx(89.390, rel=1e-3)

    assert_units_equivalent(expr, pyunits.J / pyunits.mol / pyunits.K)


@pytest.mark.unit
def test_enth_mol_liq_comp(frame):
    Perrys.enth_mol_liq_comp.build_parameters(frame.params)

    assert isinstance(frame.params.enth_mol_form_liq_comp_ref, Var)
    assert value(frame.params.enth_mol_form_liq_comp_ref) == -285.83e3

    expr = Perrys.enth_mol_liq_comp.return_expression(
        frame.props[1], frame.params, frame.props[1].temperature
    )
    assert value(expr) == value(frame.params.enth_mol_form_liq_comp_ref)

    frame.props[1].temperature.value = 533.15
    assert value(expr) == pytest.approx(-265423, rel=1e-3)

    assert_units_equivalent(expr, pyunits.J / pyunits.mol)


@pytest.mark.unit
def test_enth_mol_liq_comp_no_formation(frame):
    frame.config.include_enthalpy_of_formation = False
    frame.params.config.include_enthalpy_of_formation = False

    Perrys.enth_mol_liq_comp.build_parameters(frame.params)

    assert not hasattr(frame.params, "enth_mol_form_liq_comp_ref")

    expr = Perrys.enth_mol_liq_comp.return_expression(
        frame.props[1], frame.params, frame.props[1].temperature
    )
    assert value(expr) == 0

    frame.props[1].temperature.value = 533.15
    assert value(expr) == pytest.approx(-265423 + 285830, rel=1e-3)

    assert_units_equivalent(expr, pyunits.J / pyunits.mol)


@pytest.mark.unit
def test_entr_mol_liq_comp(frame):
    Perrys.entr_mol_liq_comp.build_parameters(frame.params)

    assert isinstance(frame.params.entr_mol_form_liq_comp_ref, Var)
    assert value(frame.params.entr_mol_form_liq_comp_ref) == 69.95

    expr = Perrys.entr_mol_liq_comp.return_expression(
        frame.props[1], frame.params, frame.props[1].temperature
    )
    assert value(expr) == value(frame.params.entr_mol_form_liq_comp_ref)

    frame.props[1].temperature.value = 533.15
    assert value(expr) == pytest.approx(122.05, rel=1e-3)

    assert_units_equivalent(expr, pyunits.J / pyunits.mol / pyunits.K)


@pytest.mark.unit
def test_dens_mol_liq_comp(frame):
    Perrys.dens_mol_liq_comp.build_parameters(frame.params)

    assert isinstance(frame.params.dens_mol_liq_comp_coeff_eqn_type, Param)
    assert value(frame.params.dens_mol_liq_comp_coeff_eqn_type) == 1

    expr = Perrys.dens_mol_liq_comp.return_expression(
        frame.props[1], frame.params, frame.props[1].temperature
    )
    assert value(expr) == pytest.approx(55.583e3, rel=1e-4)

    frame.props[1].temperature.value = 333.15
    assert value(expr) == pytest.approx(54.703e3, rel=1e-4)

    assert_units_equivalent(expr, pyunits.mol / pyunits.m**3)


@pytest.mark.unit
def test_dens_mol_liq_comp_eqn_1(frame):
    Perrys.dens_mol_liq_comp_eqn_1.build_parameters(frame.params)

    assert isinstance(frame.params.dens_mol_liq_comp_coeff_1, Var)
    assert value(frame.params.dens_mol_liq_comp_coeff_1) == 5.459
    assert isinstance(frame.params.dens_mol_liq_comp_coeff_2, Var)
    assert value(frame.params.dens_mol_liq_comp_coeff_2) == 0.30542
    assert isinstance(frame.params.dens_mol_liq_comp_coeff_3, Var)
    assert value(frame.params.dens_mol_liq_comp_coeff_3) == 647.13
    assert isinstance(frame.params.dens_mol_liq_comp_coeff_4, Var)
    assert value(frame.params.dens_mol_liq_comp_coeff_4) == 0.081

    expr = Perrys.dens_mol_liq_comp_eqn_1.return_expression(
        frame.props[1], frame.params, frame.props[1].temperature
    )
    assert value(expr) == pytest.approx(55.583e3, rel=1e-4)

    frame.props[1].temperature.value = 333.15
    assert value(expr) == pytest.approx(54.703e3, rel=1e-4)

    assert_units_equivalent(expr, pyunits.mol / pyunits.m**3)


@pytest.mark.unit
def test_dens_mol_liq_comp_eqn_2(frame):
    Perrys.dens_mol_liq_comp_eqn_2.build_parameters(frame.params)

    assert isinstance(frame.params.dens_mol_liq_comp_coeff_1, Var)
    assert value(frame.params.dens_mol_liq_comp_coeff_1) == 5.459
    assert isinstance(frame.params.dens_mol_liq_comp_coeff_2, Var)
    assert value(frame.params.dens_mol_liq_comp_coeff_2) == 0.30542
    assert isinstance(frame.params.dens_mol_liq_comp_coeff_3, Var)
    assert value(frame.params.dens_mol_liq_comp_coeff_3) == 647.13
    assert isinstance(frame.params.dens_mol_liq_comp_coeff_4, Var)
    assert value(frame.params.dens_mol_liq_comp_coeff_4) == 0.081

    expr = Perrys.dens_mol_liq_comp_eqn_2.return_expression(
        frame.props[1], frame.params, frame.props[1].temperature
    )
    assert value(expr) == pytest.approx(4.99375e10, rel=1e-4)

    frame.props[1].temperature.value = 333.15
    assert value(expr) == pytest.approx(7.48194e10, rel=1e-4)

    assert_units_equivalent(expr, pyunits.mol / pyunits.m**3)


@pytest.mark.unit
def test_dens_mol_liq_comp_exceptions(frame):
    frame.params.config.parameter_data["dens_mol_liq_comp_coeff"]["eqn_type"] = 3
    with pytest.raises(ConfigurationError):
        Perrys.dens_mol_liq_comp.build_parameters(frame.params)

    assert isinstance(frame.params.dens_mol_liq_comp_coeff_eqn_type, Param)
    assert value(frame.params.dens_mol_liq_comp_coeff_eqn_type) == 3

    with pytest.raises(ConfigurationError):
        Perrys.dens_mol_liq_comp.return_expression(
            frame.props[1], frame.params, frame.props[1].temperature
        )
