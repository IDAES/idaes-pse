##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
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
Test for custom unit model.
"""

import pytest
from pyomo.environ import ConcreteModel, Constraint, Var, Set, value, \
    SolverStatus, TerminationCondition
from idaes.core import FlowsheetBlock
from idaes.generic_models.unit_models import CustomModel
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util import get_solver
from idaes.core.util.exceptions import ConfigurationError

__author__ = "Jaffer Ghouse"


# -----------------------------------------------------------------------------
# Get default solver for testing
opt = get_solver()
# -----------------------------------------------------------------------------

m = ConcreteModel()
m.fs = FlowsheetBlock(default={"dynamic": False})

m.fs.surrogate = CustomModel()
m.fs.surrogate.comp_list = Set(initialize=["c1", "c2"])

# input vars for surrogate
m.fs.surrogate.flow_comp_in = \
    Var(m.fs.time, m.fs.surrogate.comp_list, initialize=1.0)
m.fs.surrogate.temperature_in = Var(m.fs.time, initialize=298.15)
m.fs.surrogate.pressure_in = Var(m.fs.time, initialize=101)

# output vars for surrogate
m.fs.surrogate.flow_comp_out = \
    Var(m.fs.time, m.fs.surrogate.comp_list, initialize=1.0)
m.fs.surrogate.temperature_out = Var(m.fs.time, initialize=298.15)
m.fs.surrogate.pressure_out = Var(m.fs.time, initialize=101)

# Surrogate model equations
# Flow equation


def rule_flow(m, t, i):
    return m.flow_comp_out[t, i] == m.flow_comp_in[t, i]


m.fs.surrogate.eq_flow = \
    Constraint(m.fs.time, m.fs.surrogate.comp_list, rule=rule_flow)

m.fs.surrogate.eq_temperature = Constraint(
    expr=m.fs.surrogate.temperature_out[0] ==
    m.fs.surrogate.temperature_in[0] + 10)
m.fs.surrogate.eq_pressure = Constraint(
    expr=m.fs.surrogate.pressure_out[0] ==
    m.fs.surrogate.pressure_in[0] - 10)

inlet_dict = {"flow_mol_comp": m.fs.surrogate.flow_comp_in,
              "temperature": m.fs.surrogate.temperature_in,
              "pressure": m.fs.surrogate.pressure_in}
outlet_dict = {"flow_mol_comp": m.fs.surrogate.flow_comp_out,
               "temperature": m.fs.surrogate.temperature_out,
               "pressure": m.fs.surrogate.pressure_out}

m.fs.surrogate.add_ports(name="inlet", member_dict=inlet_dict)
m.fs.surrogate.add_ports(name="outlet", member_dict=outlet_dict)


def my_initialize():
    # Callback for user provided initialization sequence
    m.fs.surrogate.c2.deactivate()
    m.fs.surrogate.x_1.fix(1)
    m.fs.surrogate.x_2.fix(2)

    opt.solve(m)
    m.fs.surrogate.c2.activate()

    results = opt.solve(m)

    assert results.solver.termination_condition == \
        TerminationCondition.optimal
    assert results.solver.status == SolverStatus.ok


def test_ports():
    # Check inlet port and port members
    assert hasattr(m.fs.surrogate, "inlet")
    assert len(m.fs.surrogate.inlet.vars) == 3
    assert hasattr(m.fs.surrogate.inlet, "flow_mol_comp")
    assert hasattr(m.fs.surrogate.inlet, "temperature")
    assert hasattr(m.fs.surrogate.inlet, "pressure")

    # Check outlet port and port members
    assert hasattr(m.fs.surrogate, "outlet")
    assert len(m.fs.surrogate.outlet.vars) == 3
    assert hasattr(m.fs.surrogate.outlet, "flow_mol_comp")
    assert hasattr(m.fs.surrogate.outlet, "temperature")
    assert hasattr(m.fs.surrogate.outlet, "pressure")


def test_build():
    # Check build of model
    assert hasattr(m.fs.surrogate, "flow_comp_in")
    assert hasattr(m.fs.surrogate, "temperature_in")
    assert hasattr(m.fs.surrogate, "pressure_in")

    assert hasattr(m.fs.surrogate, "eq_flow")

    assert hasattr(m.fs.surrogate, "flow_comp_out")
    assert hasattr(m.fs.surrogate, "temperature_out")
    assert hasattr(m.fs.surrogate, "pressure_out")


def test_default_initialize():
    # Check default initialize method raises ConfigurationError

    with pytest.raises(ConfigurationError):
        m.fs.surrogate.initialize()

    m.fs.surrogate.inlet.flow_mol_comp[0, "c1"].fix(2)
    m.fs.surrogate.inlet.flow_mol_comp[0, "c2"].fix(2)
    m.fs.surrogate.inlet.temperature.fix(325)
    m.fs.surrogate.inlet.pressure.fix(200)

    assert degrees_of_freedom(m) == 0

    m.fs.surrogate.initialize()

    results = opt.solve(m)

    assert results.solver.termination_condition == \
        TerminationCondition.optimal
    assert results.solver.status == SolverStatus.ok

    assert value(m.fs.surrogate.outlet.flow_mol_comp[0, "c1"]) == \
        pytest.approx(2, abs=1e-3)
    assert value(m.fs.surrogate.outlet.flow_mol_comp[0, "c2"]) == \
        pytest.approx(2, abs=1e-3)

    assert value(m.fs.surrogate.outlet.temperature[0]) == \
        pytest.approx(335, abs=1e-3)
    assert value(m.fs.surrogate.outlet.pressure[0]) == \
        pytest.approx(190, abs=1e-3)


# def test_custom_initialize():
#     # Check custom initialize method as a callback from user
#     m.fs.surrogate.initialize(custom_initialize=my_initialize())

#     results = opt.solve(m)

#     assert results.solver.termination_condition == \
#         TerminationCondition.optimal
#     assert results.solver.status == SolverStatus.ok

#     assert value(m.fs.surrogate.x_1) == pytest.approx(1, abs=1e-3)
#     assert value(m.fs.surrogate.x_2) == pytest.approx(2, abs=1e-3)

#     assert value(m.fs.surrogate.y_1) == pytest.approx(1, abs=1e-3)
#     assert value(m.fs.surrogate.y_2) == pytest.approx(4, abs=1e-3)
