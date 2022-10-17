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
Tests for methods from Reid, Prausnitz and Poling

All methods and parameters from:

The Properties of Gases & Liquids, 5th Edition
Reid, Prausnitz and Polling, 2001, McGraw-Hill

All parameter indicies based on conventions used by the source

Authors: Andrew Lee, Alejandro Garciadiego
"""

import pytest
import types

from pyomo.environ import ConcreteModel, Block, value, Var, units as pyunits
from pyomo.common.config import ConfigBlock
from pyomo.util.check_units import assert_units_equivalent

from idaes.models.properties.modular_properties.pure.RPP5 import *
from idaes.core.util.misc import add_object_reference
from idaes.core.base.property_meta import PropertyClassMetadata, UnitSet


@pytest.fixture()
def frame():
    m = ConcreteModel()

    # Create a dummy parameter block
    m.params = Block()

    m.params.config = ConfigBlock(implicit=True)
    m.params.config.parameter_data = {
        "cp_mol_ig_comp_coeff": {
            "a0": 4.395,
            "a1": -4.186e-3,
            "a2": 1.405e-5,
            "a3": -1.564e-8,
            "a4": 0.632e-11,
        },
        "enth_mol_form_vap_comp_ref": -241.81e3,
        "entr_mol_form_vap_comp_ref": 188.84,
        "pressure_sat_comp_coeff": {"A": 5.11564, "B": 1687.537, "C": 230.14},
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
    m.params.temperature_ref = Var(initialize=273.15, units=pyunits.K)
    m.params.pressure_ref = Var(initialize=1e5, units=pyunits.Pa)

    m.params.temperature_crit = Var(initialize=647.15, units=pyunits.K)
    m.params.pressure_crit = Var(initialize=220.64e5, units=pyunits.Pa)

    # Create a dummy state block
    m.props = Block([1])
    add_object_reference(m.props[1], "params", m.params)

    m.props[1].temperature = Var(initialize=298.15, units=pyunits.K)
    m.props[1].pressure = Var(initialize=101325, units=pyunits.Pa)

    return m


@pytest.mark.unit
def test_cp_mol_ig_comp(frame):
    RPP5.cp_mol_ig_comp.build_parameters(frame.params)

    assert isinstance(frame.params.cp_mol_ig_comp_coeff_a0, Var)
    assert value(frame.params.cp_mol_ig_comp_coeff_a0) == 4.395
    assert isinstance(frame.params.cp_mol_ig_comp_coeff_a1, Var)
    assert value(frame.params.cp_mol_ig_comp_coeff_a1) == -4.186e-3
    assert isinstance(frame.params.cp_mol_ig_comp_coeff_a2, Var)
    assert value(frame.params.cp_mol_ig_comp_coeff_a2) == 1.405e-5
    assert isinstance(frame.params.cp_mol_ig_comp_coeff_a3, Var)
    assert value(frame.params.cp_mol_ig_comp_coeff_a3) == -1.564e-8
    assert isinstance(frame.params.cp_mol_ig_comp_coeff_a4, Var)
    assert value(frame.params.cp_mol_ig_comp_coeff_a4) == 0.632e-11

    expr = RPP5.cp_mol_ig_comp.return_expression(
        frame.props[1], frame.params, frame.props[1].temperature
    )
    assert value(expr) == pytest.approx(33.518, abs=1e-3)

    frame.props[1].temperature.value = 400
    assert value(expr) == pytest.approx(34.334, abs=1e-3)

    assert_units_equivalent(expr, pyunits.J / pyunits.mol / pyunits.K)


@pytest.mark.unit
def test_enth_mol_ig_comp(frame):
    RPP5.enth_mol_ig_comp.build_parameters(frame.params)

    assert isinstance(frame.params.enth_mol_form_vap_comp_ref, Var)
    assert value(frame.params.enth_mol_form_vap_comp_ref) == -241.81e3

    expr = RPP5.enth_mol_ig_comp.return_expression(
        frame.props[1], frame.params, frame.props[1].temperature
    )
    assert value(expr) == pytest.approx(-240973.683, abs=1e-3)

    frame.props[1].temperature.value = 400
    assert value(expr) == pytest.approx(-237521.691, abs=1e-3)

    assert_units_equivalent(expr, pyunits.J / pyunits.mol)


@pytest.mark.unit
def test_entr_mol_ig_comp(frame):
    RPP5.entr_mol_ig_comp.build_parameters(frame.params)

    assert isinstance(frame.params.entr_mol_form_vap_comp_ref, Var)
    assert value(frame.params.entr_mol_form_vap_comp_ref) == 188.84

    expr = RPP5.entr_mol_ig_comp.return_expression(
        frame.props[1], frame.params, frame.props[1].temperature
    )
    assert value(expr) == pytest.approx(191.769, abs=1e-3)

    frame.props[1].temperature.value = 400
    assert value(expr) == pytest.approx(201.723, abs=1e-3)

    assert_units_equivalent(expr, pyunits.J / pyunits.mol / pyunits.K)


@pytest.mark.unit
def test_pressure_sat_comp(frame):
    RPP5.pressure_sat_comp.build_parameters(frame.params)

    assert isinstance(frame.params.pressure_sat_comp_coeff_A, Var)
    assert value(frame.params.pressure_sat_comp_coeff_A) == 5.11564
    assert isinstance(frame.params.pressure_sat_comp_coeff_B, Var)
    assert value(frame.params.pressure_sat_comp_coeff_B) == 1687.537
    assert isinstance(frame.params.pressure_sat_comp_coeff_C, Var)
    assert value(frame.params.pressure_sat_comp_coeff_C) == 230.14

    expr = RPP5.pressure_sat_comp.return_expression(
        frame.props[1], frame.params, frame.props[1].temperature
    )
    assert value(expr) == pytest.approx(3173.066, abs=1e-3)

    frame.props[1].temperature.value = 373.15
    assert value(expr) == pytest.approx(100939.248, rel=1e-3)

    assert_units_equivalent(expr, pyunits.Pa)


@pytest.mark.unit
def test_pressure_sat_comp_dT(frame):
    RPP5.pressure_sat_comp.build_parameters(frame.params)

    expr = RPP5.pressure_sat_comp.dT_expression(
        frame.props[1], frame.params, frame.props[1].temperature
    )

    delta = 1e-4 * pyunits.K
    val = RPP5.pressure_sat_comp.return_expression(
        frame.props[1], frame.params, frame.props[1].temperature
    )
    val_p = RPP5.pressure_sat_comp.return_expression(
        frame.props[1], frame.params, frame.props[1].temperature + delta
    )

    dPdT = value((val - val_p) / -delta)

    assert value(expr) == pytest.approx(dPdT, 1e-4)

    frame.props[1].temperature.value = 373.15

    val = RPP5.pressure_sat_comp.return_expression(
        frame.props[1], frame.params, frame.props[1].temperature
    )
    val_p = RPP5.pressure_sat_comp.return_expression(
        frame.props[1], frame.params, frame.props[1].temperature + delta
    )

    dPdT = value((val - val_p) / -delta)

    assert value(expr) == pytest.approx(dPdT, 1e-4)

    assert_units_equivalent(expr, pyunits.Pa / pyunits.degK)
