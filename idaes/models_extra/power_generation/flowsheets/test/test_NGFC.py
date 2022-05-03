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
Tests to make sure the NGFC example builds and solves correctly.
"""

__author__ = "Alex Noring"

import os

import pytest
import pyomo.environ as pyo
from pyomo.environ import check_optimal_termination
from pyomo.common.fileutils import this_file_dir

from idaes.models_extra.power_generation.flowsheets.NGFC.NGFC_flowsheet \
    import (
        build_power_island,
        build_reformer,
        set_power_island_inputs,
        set_reformer_inputs,
        scale_flowsheet,
        initialize_power_island,
        initialize_reformer,
        connect_reformer_to_power_island,
        SOFC_ROM_setup,
        add_SOFC_energy_balance,
        add_result_constraints,
        )

from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.solvers import get_solver
from idaes.core.util.scaling import (extreme_jacobian_columns,
                                     extreme_jacobian_rows)
from idaes.core.util import model_serializer as ms

solver_available = pyo.SolverFactory("ipopt").available()
solver = get_solver()


@pytest.fixture(scope="module")
def m():
    m = pyo.ConcreteModel(name="NGFC without carbon capture")
    m.fs = FlowsheetBlock(default={"dynamic": False})

    build_power_island(m)
    build_reformer(m)
    scale_flowsheet(m)
    set_power_island_inputs(m)
    set_reformer_inputs(m)

    return m


@pytest.mark.component
def test_build(m):
    """Build NGFC and check for unit ops"""

    assert degrees_of_freedom(m) == 0

    assert hasattr(m.fs, "anode")
    assert hasattr(m.fs.anode, "heat_duty")
    assert hasattr(m.fs, "cathode_heat")
    assert hasattr(m.fs.cathode_heat, "heat_duty")
    assert hasattr(m.fs, "reformer")
    assert hasattr(m.fs.reformer, "heat_duty")
    assert hasattr(m.fs.reformer, "deltaP")


@pytest.mark.unit
def test_scaling(m):
    scale_flowsheet(m)

    # check that less than 10% of model variables are badly scaled pre-solve
    badly_scaled_vars = extreme_jacobian_columns(m)
    all_var = list(m.component_data_objects(pyo.Var, descend_into=True))
    assert len(badly_scaled_vars)/len(all_var) < 0.1

    # check that less than 10% of model constraints are badly scaled pre-solve
    badly_scaled_cons = extreme_jacobian_rows(m)
    all_con = list(m.component_data_objects(pyo.Constraint, descend_into=True))
    assert len(badly_scaled_cons)/len(all_con) < 0.1


@pytest.mark.integration
def test_initialize(m):
    solver.options = {
        "max_iter": 500,
        "tol": 1e-7,
        "bound_push": 1e-12,
        "linear_solver": "ma57",
        "ma57_pivtol": 1e-3,
          }
    initialize_power_island(m)
    initialize_reformer(m)

    assert pyo.value(m.fs.reformer.lagrange_mult[(0, "H")]) == \
        pytest.approx(83053, 1e-3)
    assert pyo.value(m.fs.prereformer.lagrange_mult[(0, "H")]) == \
        pytest.approx(66902, 1e-3)
    assert pyo.value(m.fs.anode.lagrange_mult[(0, "H")]) == \
        pytest.approx(75542, 1e-3)
    assert pyo.value(m.fs.bypass_rejoin.outlet.temperature[0]) == \
        pytest.approx(1019, 1e-3)
    assert pyo.value(m.fs.anode_hx.shell_outlet.temperature[0]) == \
        pytest.approx(869, 1e-3)
    assert pyo.value(m.fs.cathode_hx.shell_outlet.temperature[0]) == \
        pytest.approx(475.8, 1e-3)


@pytest.mark.integration
def test_ROM(m):
    connect_reformer_to_power_island(m)
    SOFC_ROM_setup(m)
    add_SOFC_energy_balance(m)
    add_result_constraints(m)

    assert degrees_of_freedom(m) == 0

    assert not m.fs.anode_mix.feed.flow_mol[0].fixed

    assert hasattr(m.fs, "ROM_fuel_inlet_temperature")
    assert hasattr(m.fs, "ROM_internal_reformation")
    assert m.fs.SOFC.current_density.fixed
    assert not m.fs.air_blower.inlet.flow_mol[0].fixed
    assert m.fs.SOFC.deltaT_cell.fixed
    assert hasattr(m.fs, "stack_power")
    assert hasattr(m.fs, "SOFC_energy_balance")
    assert not m.fs.anode.outlet.temperature[0].fixed
    assert hasattr(m.fs, "CO2_emissions")


@pytest.mark.integration
def test_solve(m):
    solve_iteration = 0
    for i in range(1, 10):  # keep looping until condition is met
        solve_iteration += 1
        print('Solve # ', solve_iteration)
        res = solver.solve(m.fs.air_compressor_s2, tee=True)
        if 'Optimal Solution Found' in res.solver.message:
            break

    assert solve_iteration == 1  # this value can be increased if needed
    assert check_optimal_termination(res)


@pytest.mark.integration
def test_json_load(m):
    fname = os.path.join(
        os.path.join(os.path.dirname(this_file_dir()), "NGFC"),
        "NGFC_flowsheet_solution.json.gz",
    )

    ms.from_json(m, fname=fname)

    assert pyo.value(m.fs.cathode.ion_outlet.flow_mol[0]) == \
        pytest.approx(1670, 1e-5)
    assert pyo.value(m.fs.reformer_recuperator.area) == \
        pytest.approx(2328, 1e-3)
    assert pyo.value(m.fs.anode.heat_duty[0]) == \
        pytest.approx(-871373774, 1e-5)
    assert pyo.value(m.fs.CO2_emissions) == \
        pytest.approx(300, 1e-5)
    assert pyo.value(m.fs.net_power) == \
        pytest.approx(660, 1e-5)
