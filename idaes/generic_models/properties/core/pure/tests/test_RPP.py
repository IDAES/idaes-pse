##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
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
Tests for methods from Reid, Prausnitz and Poling

All methods and parameters from:

The Properties of Gases & Liquids, 4th Edition
Reid, Prausnitz and Polling, 1987, McGraw-Hill

All parameter indicies based on conventions used by the source

Authors: Andrew Lee
"""

import pytest

from pyomo.environ import ConcreteModel, Block, value, Var

from idaes.generic_models.properties.core.pure.RPP import *
from idaes.core.util.misc import add_object_reference


@pytest.fixture()
def frame():
    m = ConcreteModel()

    # Create a dummy parameter block
    m.params = Block()

    # Add necessary parameters to parameter block
    m.params.temperature_ref = Var(initialize=273.15)
    m.params.pressure_ref = Var(initialize=1e5)

    m.params.pressure_sat_comp_coeff = Var(["H2O"], ["A", "B", "C", "D"])
    m.params.temperature_crit_comp = Var(["H2O"], initialize=647.3)
    m.params.pressure_crit_comp = Var(["H2O"], initialize=221.2e5)
    m.params.cp_mol_ig_comp_coeff = Var(["H2O"], ["A", "B", "C", "D"])
    m.params.enth_mol_form_phase_comp_ref = Var(["Vap"], ["H2O"])
    m.params.entr_mol_phase_comp_ref = Var(["Vap"], ["H2O"])

    m.params.cp_mol_ig_comp_coeff["H2O", "A"].value = 3.224e1
    m.params.cp_mol_ig_comp_coeff["H2O", "B"].value = 1.924e-3
    m.params.cp_mol_ig_comp_coeff["H2O", "C"].value = 1.055e-5
    m.params.cp_mol_ig_comp_coeff["H2O", "D"].value = -3.596e-9

    m.params.pressure_sat_comp_coeff["H2O", "A"].value = -7.76451
    m.params.pressure_sat_comp_coeff["H2O", "B"].value = 1.45838
    m.params.pressure_sat_comp_coeff["H2O", "C"].value = -2.77580
    m.params.pressure_sat_comp_coeff["H2O", "D"].value = -1.23303

    m.params.enth_mol_form_phase_comp_ref["Vap", "H2O"].value = -241.83e3
    m.params.entr_mol_phase_comp_ref["Vap", "H2O"].value = 188.84

    # Create a dummy state block
    m.props = Block([1])
    add_object_reference(m.props[1], "params", m.params)

    m.props[1].temperature = Var(initialize=298.15)
    m.props[1].pressure = Var(initialize=101325)

    return m


def test_cp_mol_ig_comp(frame):
    expr = cp_mol_ig_comp(frame.props[1], "H2O", frame.props[1].temperature)
    assert value(expr) == pytest.approx(33.656, abs=1e-3)

    frame.props[1].temperature.value = 400
    assert value(expr) == pytest.approx(34.467, abs=1e-3)


def test_enth_mol_ig_comp(frame):
    expr = enth_mol_ig_comp(frame.props[1], "H2O", frame.props[1].temperature)
    assert value(expr) == pytest.approx(-240990.825, abs=1e-3)

    frame.props[1].temperature.value = 400
    assert value(expr) == pytest.approx(-237522.824, abs=1e-3)


def test_entr_mol_ig_comp(frame):
    expr = entr_mol_ig_comp(frame.props[1], "H2O", frame.props[1].temperature)
    assert value(expr) == pytest.approx(373.541, abs=1e-3)

    frame.props[1].temperature.value = 400
    assert value(expr) == pytest.approx(383.542, abs=1e-3)


def test_pressure_sat_comp(frame):
    expr = pressure_sat_comp(frame.props[1], "H2O", frame.props[1].temperature)
    assert value(expr) == pytest.approx(3171.4391, abs=1e-3)

    frame.props[1].temperature.value = 373.15
    assert value(expr) == pytest.approx(101378, rel=1e-4)


def test_pressure_sat_comp_dT(frame):
    expr = pressure_sat_comp_dT(
            frame.props[1], "H2O", frame.props[1].temperature)

    delta = 1e-4
    val = pressure_sat_comp(frame.props[1], "H2O", frame.props[1].temperature)
    val_p = pressure_sat_comp(frame.props[1],
                              "H2O",
                              frame.props[1].temperature+delta)

    dPdT = value((val-val_p)/-delta)

    assert value(expr) == pytest.approx(dPdT, 1e-4)

    frame.props[1].temperature.value = 373.15

    val = pressure_sat_comp(frame.props[1], "H2O", frame.props[1].temperature)
    val_p = pressure_sat_comp(frame.props[1],
                              "H2O",
                              frame.props[1].temperature+delta)

    dPdT = value((val-val_p)/-delta)

    assert value(expr) == pytest.approx(dPdT, 1e-4)
