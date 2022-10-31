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

The Properties of Gases & Liquids, 3rd Edition
Reid, Prausnitz and Polling, 1977, McGraw-Hill

All parameter indicies based on conventions used by the source

Authors: Andrew Lee, Alejandro Garciadiego
"""

import pytest
import types

from pyomo.environ import ConcreteModel, Block, value, Var, units as pyunits
from pyomo.common.config import ConfigBlock
from pyomo.util.check_units import assert_units_equivalent

from idaes.models.properties.modular_properties.pure.RPP3 import *
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
            "A": 7.701,
            "B": 4.595e-4,
            "C": 2.521e-6,
            "D": -0.859e-9,
        },
        "enth_mol_form_vap_comp_ref": (-57.797e3, pyunits.cal / pyunits.mol),
        "entr_mol_form_vap_comp_ref": (45.13, pyunits.cal / pyunits.mol / pyunits.K),
        "pressure_sat_comp_coeff": {"A": 18.3036, "B": 3816.44, "C": -46.13},
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

    m.params.temperature_crit = Var(initialize=647.3, units=pyunits.K)
    m.params.pressure_crit = Var(initialize=217.6e5, units=pyunits.Pa)

    # Create a dummy state block
    m.props = Block([1])
    add_object_reference(m.props[1], "params", m.params)

    m.props[1].temperature = Var(initialize=298.15, units=pyunits.K)
    m.props[1].pressure = Var(initialize=101325, units=pyunits.Pa)

    return m


@pytest.mark.unit
def test_cp_mol_ig_comp(frame):
    RPP3.cp_mol_ig_comp.build_parameters(frame.params)

    assert isinstance(frame.params.cp_mol_ig_comp_coeff_A, Var)
    assert value(frame.params.cp_mol_ig_comp_coeff_A) == 7.701
    assert isinstance(frame.params.cp_mol_ig_comp_coeff_B, Var)
    assert value(frame.params.cp_mol_ig_comp_coeff_B) == 4.595e-4
    assert isinstance(frame.params.cp_mol_ig_comp_coeff_C, Var)
    assert value(frame.params.cp_mol_ig_comp_coeff_C) == 2.521e-6
    assert isinstance(frame.params.cp_mol_ig_comp_coeff_D, Var)
    assert value(frame.params.cp_mol_ig_comp_coeff_D) == -0.859e-9

    expr = RPP3.cp_mol_ig_comp.return_expression(
        frame.props[1], frame.params, frame.props[1].temperature
    )
    assert value(expr) == pytest.approx(33.636, abs=1e-3)

    frame.props[1].temperature.value = 400
    assert value(expr) == pytest.approx(34.447, abs=1e-3)

    assert_units_equivalent(expr, pyunits.J / pyunits.mol / pyunits.K)


@pytest.mark.unit
def test_enth_mol_ig_comp(frame):
    RPP3.enth_mol_ig_comp.build_parameters(frame.params)

    assert isinstance(frame.params.enth_mol_form_vap_comp_ref, Var)
    assert value(frame.params.enth_mol_form_vap_comp_ref) == (
        pytest.approx(-241822.6, abs=1e-1)
    )

    expr = RPP3.enth_mol_ig_comp.return_expression(
        frame.props[1], frame.params, frame.props[1].temperature
    )
    assert value(expr) == pytest.approx(-240983.962, abs=1e-3)

    frame.props[1].temperature.value = 400
    assert value(expr) == pytest.approx(-237517.968, abs=1e-3)

    assert_units_equivalent(expr, pyunits.J / pyunits.mol)


@pytest.mark.unit
def test_enth_mol_ig_comp_no_form(frame):
    frame.config.include_enthalpy_of_formation = False
    frame.params.config.include_enthalpy_of_formation = False
    RPP3.enth_mol_ig_comp.build_parameters(frame.params)

    assert not hasattr(frame.params, "enth_mol_form_vap_comp_ref")

    expr = RPP3.enth_mol_ig_comp.return_expression(
        frame.props[1], frame.params, frame.props[1].temperature
    )
    assert value(expr) == pytest.approx(838.686, abs=1e-3)

    frame.props[1].temperature.value = 400
    assert value(expr) == pytest.approx(4304.680, abs=1e-3)

    assert_units_equivalent(expr, pyunits.J / pyunits.mol)


@pytest.mark.unit
def test_entr_mol_ig_comp(frame):
    RPP3.entr_mol_ig_comp.build_parameters(frame.params)

    assert isinstance(frame.params.entr_mol_form_vap_comp_ref, Var)
    assert value(frame.params.entr_mol_form_vap_comp_ref) == (
        pytest.approx(188.8, abs=1e-1)
    )

    expr = RPP3.entr_mol_ig_comp.return_expression(
        frame.props[1], frame.params, frame.props[1].temperature
    )
    assert value(expr) == pytest.approx(191.761, abs=1e-3)

    frame.props[1].temperature.value = 400
    assert value(expr) == pytest.approx(201.756, abs=1e-3)

    assert_units_equivalent(expr, pyunits.J / pyunits.mol / pyunits.K)


@pytest.mark.unit
def test_pressure_sat_comp(frame):
    RPP3.pressure_sat_comp.build_parameters(frame.params)

    assert isinstance(frame.params.pressure_sat_comp_coeff_A, Var)
    assert value(frame.params.pressure_sat_comp_coeff_A) == 18.3036
    assert isinstance(frame.params.pressure_sat_comp_coeff_B, Var)
    assert value(frame.params.pressure_sat_comp_coeff_B) == 3816.44
    assert isinstance(frame.params.pressure_sat_comp_coeff_C, Var)
    assert value(frame.params.pressure_sat_comp_coeff_C) == -46.13

    expr = RPP3.pressure_sat_comp.return_expression(
        frame.props[1], frame.params, frame.props[1].temperature
    )
    assert value(expr) == pytest.approx(3143.1125, abs=1e-3)

    frame.props[1].temperature.value = 373.15
    assert value(expr) == pytest.approx(101317, rel=1e-4)

    assert_units_equivalent(expr, pyunits.Pa)


@pytest.mark.unit
def test_pressure_sat_comp_dT(frame):
    RPP3.pressure_sat_comp.build_parameters(frame.params)

    expr = RPP3.pressure_sat_comp.dT_expression(
        frame.props[1], frame.params, frame.props[1].temperature
    )

    delta = 1e-4 * pyunits.K
    val = RPP3.pressure_sat_comp.return_expression(
        frame.props[1], frame.params, frame.props[1].temperature
    )
    val_p = RPP3.pressure_sat_comp.return_expression(
        frame.props[1], frame.params, frame.props[1].temperature + delta
    )

    dPdT = value((val - val_p) / -delta)

    assert value(expr) == pytest.approx(dPdT, 1e-4)

    frame.props[1].temperature.value = 373.15

    val = RPP3.pressure_sat_comp.return_expression(
        frame.props[1], frame.params, frame.props[1].temperature
    )
    val_p = RPP3.pressure_sat_comp.return_expression(
        frame.props[1], frame.params, frame.props[1].temperature + delta
    )

    dPdT = value((val - val_p) / -delta)

    assert value(expr) == pytest.approx(dPdT, 1e-4)

    assert_units_equivalent(expr, pyunits.Pa / pyunits.K)
