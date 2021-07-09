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
Make sure the supercritical steam cycle example solves.
"""

__author__ = "John Eslick"

import pytest
import pyomo.environ as pyo
from pyomo.util.check_units import assert_units_consistent

from idaes.power_generation.flowsheets.supercritical_steam_cycle.supercritical_steam_cycle import main
from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              activated_equalities_generator)
from idaes.generic_models.properties import iapws95


@pytest.fixture(scope="module")
def model():
    m, solver = main()
    m.solver = solver

    return m


def gross_power_mw(model):
    # pyo.value(m.fs.turb.power[0]) is the power consumed in Watts
    return -pyo.value(model.fs.turb.power[0])/1e6


@pytest.mark.integration
def test_init(model):
    # check that the model solved properly and has 0 degrees of freedom
    assert(degrees_of_freedom(model) == 0)
    for c in activated_equalities_generator(model):
        assert(abs(c.body() - c.lower) < 5e-4)


# @pytest.mark.integration
# def test_unit_consistency(model):
#     assert_units_consistent(model)


@pytest.mark.integration
def test_init_value(model):
    assert gross_power_mw(model) == pytest.approx(622.38, abs=1e-2)


@pytest.mark.integration
def test_valve_change(model):
    model.fs.turb.throttle_valve[1].valve_opening[:].value = 0.25
    model.solver.solve(model, tee=True)
    assert gross_power_mw(model) == pytest.approx(594.66, abs=1e-2)
