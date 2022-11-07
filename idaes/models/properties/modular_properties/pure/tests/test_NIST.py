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
Tests for methods from NIST

All methods and parameters from:

https://webbook.nist.gov, retrieved 26 November 2019

All parameter indicies based on conventions used by the source

Authors: Andrew Lee
"""

import pytest
import types

from pyomo.environ import ConcreteModel, Block, Expression, value, Var, units as pyunits
from pyomo.common.config import ConfigBlock
from pyomo.util.check_units import assert_units_equivalent

from idaes.models.properties.modular_properties.pure.NIST import *
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
            "A": 30.09200,  # parameters for water
            "B": 6.832514,
            "C": 6.793435,
            "D": -2.534480,
            "E": 0.082139,
            "F": -250.8810,
            "G": 223.3967,
            "H": -241.8264,
        },
        "pressure_sat_comp_coeff": {
            "A": 3.55959,  # units bar, K
            "B": 643.748,
            "C": -198.043,
        },
    }
    m.params.config.include_enthalpy_of_formation = True

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
    m.params.temperature_ref = Var(initialize=298.15, units=pyunits.K)
    m.params.pressure_ref = Var(initialize=1e5, units=pyunits.Pa)

    # Create a dummy state block
    m.props = Block([1])
    add_object_reference(m.props[1], "params", m.params)

    m.props[1].temperature = Var(initialize=500, units=pyunits.K)
    m.props[1].pressure = Var(initialize=101325, units=pyunits.Pa)

    return m


@pytest.mark.unit
def test_cp_mol_ig_comp(frame):
    NIST.cp_mol_ig_comp.build_parameters(frame.params)

    assert isinstance(frame.params.cp_mol_ig_comp_coeff_A, Var)
    assert value(frame.params.cp_mol_ig_comp_coeff_A) == 30.09200
    assert isinstance(frame.params.cp_mol_ig_comp_coeff_B, Var)
    assert value(frame.params.cp_mol_ig_comp_coeff_B) == 6.832514
    assert isinstance(frame.params.cp_mol_ig_comp_coeff_C, Var)
    assert value(frame.params.cp_mol_ig_comp_coeff_C) == 6.793435
    assert isinstance(frame.params.cp_mol_ig_comp_coeff_D, Var)
    assert value(frame.params.cp_mol_ig_comp_coeff_D) == -2.534480
    assert isinstance(frame.params.cp_mol_ig_comp_coeff_E, Var)
    assert value(frame.params.cp_mol_ig_comp_coeff_E) == 0.082139
    assert isinstance(frame.params.cp_mol_ig_comp_coeff_F, Var)
    assert value(frame.params.cp_mol_ig_comp_coeff_F) == -250.8810
    assert isinstance(frame.params.cp_mol_ig_comp_coeff_G, Var)
    assert value(frame.params.cp_mol_ig_comp_coeff_G) == 223.3967
    assert isinstance(frame.params.cp_mol_ig_comp_coeff_H, Var)
    assert value(frame.params.cp_mol_ig_comp_coeff_H) == -241.8264

    assert isinstance(frame.params.enth_mol_form_vap_comp_ref, Expression)
    assert value(frame.params.enth_mol_form_vap_comp_ref) == -241.8264 * 1e3

    expr = NIST.cp_mol_ig_comp.return_expression(
        frame.props[1], frame.params, frame.props[1].temperature
    )
    assert value(expr) == pytest.approx(35.22, abs=1e-2)  # value from NIST

    frame.props[1].temperature.value = 600
    assert value(expr) == pytest.approx(36.32, abs=1e-2)  # value from NIST

    assert_units_equivalent(expr, pyunits.J / pyunits.mol / pyunits.K)


@pytest.mark.unit
def test_enth_mol_ig_comp(frame):
    NIST.enth_mol_ig_comp.build_parameters(frame.params)

    expr = NIST.enth_mol_ig_comp.return_expression(
        frame.props[1], frame.params, frame.props[1].temperature
    )
    assert value(expr) == pytest.approx(6.92e3 - 241826.4, rel=1e-3)  # value from NIST

    frame.props[1].temperature.value = 600
    assert value(expr) == pytest.approx(10.5e3 - 241826.4, rel=1e-3)  # value from NIST

    assert_units_equivalent(expr, pyunits.J / pyunits.mol)


@pytest.mark.unit
def test_enth_mol_ig_comp_no_formation(frame):
    NIST.enth_mol_ig_comp.build_parameters(frame.params)
    frame.params.config.include_enthalpy_of_formation = False

    expr = NIST.enth_mol_ig_comp.return_expression(
        frame.props[1], frame.params, frame.props[1].temperature
    )
    assert value(expr) == pytest.approx(6.92e3, rel=1e-3)  # value from NIST

    frame.props[1].temperature.value = 600
    assert value(expr) == pytest.approx(10.5e3, rel=1e-3)  # value from NIST

    assert_units_equivalent(expr, pyunits.J / pyunits.mol)


@pytest.mark.unit
def test_entr_mol_ig_comp_no_formation(frame):
    NIST.entr_mol_ig_comp.build_parameters(frame.params)
    frame.params.config.include_entropy_of_formation = False

    expr = NIST.entr_mol_ig_comp.return_expression(
        frame.props[1], frame.params, frame.props[1].temperature
    )
    assert value(expr) == pytest.approx(206.5, rel=1e-3)

    frame.props[1].temperature.value = 600
    assert value(expr) == pytest.approx(213.1, rel=1e-3)

    assert_units_equivalent(expr, pyunits.J / pyunits.mol / pyunits.K)


@pytest.mark.unit
def test_pressure_sat_comp(frame):
    NIST.pressure_sat_comp.build_parameters(frame.params)

    assert isinstance(frame.params.pressure_sat_comp_coeff_A, Var)
    assert value(frame.params.pressure_sat_comp_coeff_A) == 3.55959
    assert isinstance(frame.params.pressure_sat_comp_coeff_B, Var)
    assert value(frame.params.pressure_sat_comp_coeff_B) == 643.748
    assert isinstance(frame.params.pressure_sat_comp_coeff_C, Var)
    assert value(frame.params.pressure_sat_comp_coeff_C) == -198.043

    expr = NIST.pressure_sat_comp.return_expression(
        frame.props[1], frame.params, frame.props[1].temperature
    )
    assert value(expr) == pytest.approx(2677137, rel=1e-4)

    frame.props[1].temperature.value = 379
    assert value(expr) == pytest.approx(100490, rel=1e-4)

    assert_units_equivalent(expr, pyunits.Pa)


@pytest.mark.unit
def test_pressure_sat_comp_dT(frame):
    NIST.pressure_sat_comp.build_parameters(frame.params)

    expr = NIST.pressure_sat_comp.dT_expression(
        frame.props[1], frame.params, frame.props[1].temperature
    )

    delta = 1e-4 * pyunits.K
    val = NIST.pressure_sat_comp.return_expression(
        frame.props[1], frame.params, frame.props[1].temperature
    )
    val_p = NIST.pressure_sat_comp.return_expression(
        frame.props[1], frame.params, frame.props[1].temperature + delta
    )

    dPdT = value((val - val_p) / -delta)

    assert value(expr) == pytest.approx(dPdT, 1e-4)

    frame.props[1].temperature.value = 373.15

    val = NIST.pressure_sat_comp.return_expression(
        frame.props[1], frame.params, frame.props[1].temperature
    )
    val_p = NIST.pressure_sat_comp.return_expression(
        frame.props[1], frame.params, frame.props[1].temperature + delta
    )

    dPdT = value((val - val_p) / -delta)

    assert value(expr) == pytest.approx(dPdT, 1e-4)

    assert_units_equivalent(expr, pyunits.Pa / pyunits.K)
