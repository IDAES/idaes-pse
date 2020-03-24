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
Tests for methods from Perry's

All methods and parameters from:

Perry's Chemical Engineers' Handbook, 7th Edition
Perry, Green, Maloney, 1997, McGraw-Hill

All parameter indicies based on conventions used by the source

Authors: Andrew Lee
"""

import pytest

from pyomo.environ import ConcreteModel, Block, value, Var

from idaes.generic_models.properties.core.pure.Perrys import *
from idaes.core.util.misc import add_object_reference


@pytest.fixture()
def frame():
    m = ConcreteModel()

    # Create a dummy parameter block
    m.params = Block()

    # Add necessary parameters to parameter block
    m.params.temperature_ref = Var(initialize=273.16)

    m.params.dens_mol_liq_comp_coeff = Var(["H2O"], ["1", "2", "3", "4"])
    m.params.cp_mol_liq_comp_coeff = Var(["H2O"], ["1", "2", "3", "4", "5"])
    m.params.enth_mol_form_phase_comp_ref = Var(["Liq"], ["H2O"])
    m.params.entr_mol_phase_comp_ref = Var(["Liq"], ["H2O"])

    m.params.cp_mol_liq_comp_coeff["H2O", "1"].value = 2.7637e+05
    m.params.cp_mol_liq_comp_coeff["H2O", "2"].value = -2.0901e+03
    m.params.cp_mol_liq_comp_coeff["H2O", "3"].value = 8.1250e+00
    m.params.cp_mol_liq_comp_coeff["H2O", "4"].value = -1.4116e-2
    m.params.cp_mol_liq_comp_coeff["H2O", "5"].value = 9.3701e-06

    m.params.dens_mol_liq_comp_coeff["H2O", "1"].value = 5.459e3  # Factor 1e3 for unit conversion
    m.params.dens_mol_liq_comp_coeff["H2O", "2"].value = 0.30542
    m.params.dens_mol_liq_comp_coeff["H2O", "3"].value = 647.13
    m.params.dens_mol_liq_comp_coeff["H2O", "4"].value = 0.081

    m.params.enth_mol_form_phase_comp_ref["Liq", "H2O"].value = -285.83e3
    m.params.entr_mol_phase_comp_ref["Liq", "H2O"].value = 69.95

    # Create a dummy state block
    m.props = Block([1])
    add_object_reference(m.props[1], "params", m.params)

    m.props[1].temperature = Var(initialize=273.16)

    return m


def test_cp_mol_liq_comp(frame):
    expr = cp_mol_liq_comp(frame.props[1], "H2O", frame.props[1].temperature)
    assert value(expr) == pytest.approx(76.150, rel=1e-3)

    frame.props[1].temperature.value = 533.15
    assert value(expr) == pytest.approx(89.390, rel=1e-3)


def test_enth_mol_liq_comp(frame):
    expr = enth_mol_liq_comp(frame.props[1], "H2O", frame.props[1].temperature)
    assert value(expr) == value(
            frame.params.enth_mol_form_phase_comp_ref["Liq", "H2O"])

    frame.props[1].temperature.value = 533.15
    assert value(expr) == pytest.approx(-265423, rel=1e-3)


def test_entr_mol_liq_comp(frame):
    expr = entr_mol_liq_comp(frame.props[1], "H2O", frame.props[1].temperature)
    assert value(expr) == pytest.approx(1270, rel=1e-3)

    frame.props[1].temperature.value = 533.15
    assert value(expr) == pytest.approx(1322, rel=1e-3)


def test_dens_mol_liq_comp(frame):
    expr = dens_mol_liq_comp(frame.props[1], "H2O", frame.props[1].temperature)
    assert value(expr) == pytest.approx(55.583e3, rel=1e-4)

    frame.props[1].temperature.value = 333.15
    assert value(expr) == pytest.approx(54.703e3, rel=1e-4)
