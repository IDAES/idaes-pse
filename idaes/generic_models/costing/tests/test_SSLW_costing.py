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
Tests for costing package based on methods from:

    Process and Product Design Principles: Synthesis, Analysis, and
    Evaluation
    Seider, Seader, Lewin, Windagdo, 3rd Ed. John Wiley and Sons
    Chapter 22. Cost Accounting and Capital Cost Estimation
    22.2 Cost Indexes and Capital Investment
"""
import pytest

from pyomo.environ import (Block,
                           check_optimal_termination,
                           ConcreteModel,
                           Constraint,
                           Param,
                           units as pyunits,
                           value,
                           Var)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.generic_models.costing import FlowsheetCostingBlock
from idaes.core.util import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom

from idaes.generic_models.costing.SSLW import SSLWCosting


# Some more information about this module
__author__ = "Andrew Lee"


solver = get_solver()


@pytest.fixture
def model():
    m = ConcreteModel()

    m.fs = FlowsheetBlock()

    m.fs.costing = FlowsheetCostingBlock(
        default={"costing_package": SSLWCosting})

    # Add a placeholder to represent a unit model
    m.fs.unit = Block()

    return m


@pytest.mark.unit
def test_global_definitions(model):
    assert "USD_500" in pyunits.pint_registry

    CEI = {"USD2010": 550.8,
           "USD2011": 585.7,
           "USD2012": 584.6,
           "USD2013": 567.3,
           "USD2014": 576.1,
           "USD2015": 556.8,
           "USD2016": 541.7,
           "USD2017": 567.5,
           "USD2018": 671.1,
           "USD2019": 680.0}

    for c, conv in CEI.items():
        assert c in pyunits.pint_registry

        assert pytest.approx(conv/500, rel=1e-10) == pyunits.convert_value(
            1, pyunits.USD_500, getattr(pyunits, c))


@pytest.mark.component
def test_cost_heat_exchanger(model):
    model.fs.unit.area = Param(initialize=1000,
                               units=pyunits.m**2)
    model.fs.unit.tube = Block()
    model.fs.unit.tube.properties_in = Block(model.fs.time)
    model.fs.unit.tube.properties_in[0].pressure = Param(
        initialize=2, units=pyunits.atm)

    model.fs.costing.cost_unit(model.fs.unit, SSLWCosting.hx_costing)

    assert isinstance(model.fs.unit_costing["unit"].base_cost_per_unit, Var)
    assert isinstance(model.fs.unit_costing["unit"].capital_cost, Var)
    assert isinstance(model.fs.unit_costing["unit"].number_of_units, Var)
    assert isinstance(model.fs.unit_costing["unit"].pressure_factor, Var)
    assert isinstance(model.fs.unit_costing["unit"].material_factor, Var)

    assert isinstance(model.fs.unit_costing["unit"].capital_cost_constraint,
                      Constraint)
    assert isinstance(model.fs.unit_costing["unit"].hx_material_eqn,
                      Constraint)
    assert isinstance(model.fs.unit_costing["unit"].p_factor_eq,
                      Constraint)
    assert isinstance(model.fs.unit_costing["unit"].capital_cost_constraint,
                      Constraint)

    assert degrees_of_freedom(model) == 0
    assert_units_consistent(model.fs.unit_costing)

    res = solver.solve(model)

    assert check_optimal_termination(res)

    assert pytest.approx(87704.6, 1e-5) == value(
        model.fs.unit_costing["unit"].base_cost_per_unit)
    assert pytest.approx(0.982982, 1e-5) == value(
        model.fs.unit_costing["unit"].pressure_factor)
    assert pytest.approx(4.08752, 1e-5) == value(
        model.fs.unit_costing["unit"].material_factor)

    assert pytest.approx(529737, 1e-5) == value(pyunits.convert(
        model.fs.unit_costing["unit"].capital_cost, to_units=pyunits.USD2018))
