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

from idaes.property_models.core.pure.RPP import *
from idaes.core.util.misc import add_object_reference


@pytest.fixture()
def frame():
    m = ConcreteModel()

    # Create a dummy parameter block
    m.params = Block()

    # Add necessary parameters to parameter block
    m.params.temperature_ref = Var(initialize=273.15)
    m.params.pressure_ref = Var(initialize=1e5)

    m.params.pressure_sat_coeff = Var(["H2O"], ["A", "B", "C", "D"])
    m.params.temperature_crit = Var(["H2O"], initialize=647.3)
    m.params.pressure_crit = Var(["H2O"], initialize=221.2e5)
    m.params.cp_ig_coeff = Var(["H2O"], ["A", "B", "C", "D"])

    m.params.cp_ig_coeff["H2O", "A"].value = 3.224e1
    m.params.cp_ig_coeff["H2O", "B"].value = 1.924e-3
    m.params.cp_ig_coeff["H2O", "C"].value = 1.055e-5
    m.params.cp_ig_coeff["H2O", "D"].value = -3.596e-9

    m.params.pressure_sat_coeff["H2O", "A"].value = -7.76451
    m.params.pressure_sat_coeff["H2O", "B"].value = 1.45838
    m.params.pressure_sat_coeff["H2O", "C"].value = -2.77580
    m.params.pressure_sat_coeff["H2O", "D"].value = -1.23303

    # Create a dummy state block
    m.props = Block([1])
    add_object_reference(m.props[1], "_params", m.params)

    m.props[1].temperature = Var(initialize=298.15)
    m.props[1].pressure = Var(initialize=101325)

    return m


def test_cp_mol_ig(frame):
    expr = cp_mol_ig(frame.props[1], "H2O", frame.props[1].temperature)
    assert value(expr) == pytest.approx(33.656, abs=1e-3)

    frame.props[1].temperature.value = 400
    assert value(expr) == pytest.approx(34.467, abs=1e-3)


def test_enth_mol_ig(frame):
    expr = enth_mol_ig(frame.props[1], "H2O", frame.props[1].temperature)
    assert value(expr) == pytest.approx(839.175, abs=1e-3)

    frame.props[1].temperature.value = 400
    assert value(expr) == pytest.approx(4307.176, abs=1e-3)


def test_entr_mol_ig(frame):
    expr = entr_mol_ig(frame.props[1], "H2O", frame.props[1].temperature)
    assert value(expr) == pytest.approx(184.701, abs=1e-3)

    frame.props[1].temperature.value = 400
    assert value(expr) == pytest.approx(194.702, abs=1e-3)


def test_pressure_sat(frame):
    expr = pressure_sat(frame.props[1], "H2O", frame.props[1].temperature)
    assert value(expr) == pytest.approx(3171.4391, abs=1e-3)

    frame.props[1].temperature.value = 373.15
    assert value(expr) == pytest.approx(101378, rel=1e-4)


def test_pressure_sat_dT(frame):
    expr = pressure_sat_dT(frame.props[1], "H2O", frame.props[1].temperature)

    delta = 1e-4
    val = pressure_sat(frame.props[1], "H2O", frame.props[1].temperature)
    val_p = pressure_sat(frame.props[1],
                         "H2O",
                         frame.props[1].temperature+delta)

    dPdT = value((val-val_p)/-delta)

    assert value(expr) == pytest.approx(dPdT, 1e-4)

    frame.props[1].temperature.value = 373.15

    val = pressure_sat(frame.props[1], "H2O", frame.props[1].temperature)
    val_p = pressure_sat(frame.props[1],
                         "H2O",
                         frame.props[1].temperature+delta)

    dPdT = value((val-val_p)/-delta)

    assert value(expr) == pytest.approx(dPdT, 1e-4)
