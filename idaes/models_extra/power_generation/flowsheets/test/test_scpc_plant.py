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

__author__ = "Miguel Zamarripa"

import pytest
from pyomo.environ import check_optimal_termination, value
from pyomo.util.check_units import assert_units_consistent

import idaes.models_extra.power_generation.flowsheets.supercritical_power_plant.boiler_subflowsheet_build as blr
import idaes.models_extra.power_generation.flowsheets.supercritical_power_plant.SCPC_full_plant as SCPC
from idaes.core.util.model_statistics import degrees_of_freedom


@pytest.fixture(scope="module")
def boiler():
    m, solver = blr.main()

    m.solver = solver

    return m


@pytest.mark.integration
def test_init(boiler):
    # initialize each unit at the time
    blr.initialize(boiler)
    blr.unfix_inlets(boiler)
    # check that the model solved properly and has 0 degrees of freedom
    assert degrees_of_freedom(boiler) == 0


@pytest.mark.integration
def test_unit_consistency(boiler):
    assert_units_consistent(boiler)


@pytest.mark.integration
def test_boiler(boiler):
    # unfix inlets to build arcs at the flowsheet level
    boiler.fs.ATMP1.outlet.enth_mol[0].fix(62710.01)
    boiler.fs.ATMP1.SprayWater.flow_mol[0].unfix()
    result = boiler.solver.solve(boiler, tee=False)
    assert check_optimal_termination(result)
    assert value(
        boiler.fs.ECON.cold_side.properties_out[0].temperature
    ) == pytest.approx(521.009, 1)


@pytest.mark.integration
def test_power_plant():
    # SCPC.main imports and solves the SCPC Power Plant Flowsheet
    m, results = SCPC.main()
    assert check_optimal_termination(results)
