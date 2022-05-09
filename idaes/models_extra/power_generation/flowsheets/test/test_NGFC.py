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
Tests to make sure the NGFC example builds and solves correctly.
"""

__author__ = "Alex Noring, Brandon Paul"

import os
from collections import OrderedDict

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
        make_stream_dict,
        pfd_result,
        )

from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.solvers import get_solver
from idaes.core.util import scaling as iscale
from idaes.core.util.scaling import (extreme_jacobian_columns,
                                     extreme_jacobian_rows)
from idaes.core.util import model_serializer as ms
from idaes.core.util.tables import create_stream_table_dataframe

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

    # check that less than 10% of model variables are badly scaled pre-solve
    badly_scaled_vars = extreme_jacobian_columns(m)
    all_var = list(m.component_data_objects(pyo.Var, descend_into=True))
    assert len(badly_scaled_vars)/len(all_var) < 0.1

    # check that less than 10% of model constraints are badly scaled pre-solve
    badly_scaled_cons = extreme_jacobian_rows(m)
    all_con = list(m.component_data_objects(pyo.Constraint, descend_into=True))
    assert len(badly_scaled_cons)/len(all_con) < 0.1

    # some specific scaling checks

    # heat exchanger areas and overall heat transfer coefficiencts
    assert iscale.get_scaling_factor(
        m.fs.anode_hx.area) == \
        pytest.approx(1e-4, rel=1e-5)
    assert iscale.get_scaling_factor(
        m.fs.anode_hx.overall_heat_transfer_coefficient) == \
        pytest.approx(1, rel=1e-5)
    assert iscale.get_scaling_factor(
        m.fs.cathode_hx.area) == \
        pytest.approx(1e-4, rel=1e-5)
    assert iscale.get_scaling_factor(
        m.fs.cathode_hx.overall_heat_transfer_coefficient) == \
        pytest.approx(1, rel=1e-5)
    assert iscale.get_scaling_factor(
        m.fs.reformer_recuperator.area) == \
        pytest.approx(1e-4, rel=1e-5)
    assert iscale.get_scaling_factor(
        m.fs.reformer_recuperator.overall_heat_transfer_coefficient) == \
        pytest.approx(1, rel=1e-5)

    # control volume heats
    assert iscale.get_scaling_factor(
        m.fs.anode_hx.tube.heat) == \
        pytest.approx(1e-7, rel=1e-5)
    assert iscale.get_scaling_factor(
        m.fs.anode_hx.shell.heat) == \
        pytest.approx(1e-7, rel=1e-5)
    assert iscale.get_scaling_factor(
        m.fs.anode.control_volume.heat) == \
        pytest.approx(1e-8, rel=1e-5)
    assert iscale.get_scaling_factor(
        m.fs.cathode_hx.tube.heat) == \
        pytest.approx(1e-8, rel=1e-5)
    assert iscale.get_scaling_factor(
        m.fs.cathode_hx.shell.heat) == \
        pytest.approx(1e-8, rel=1e-5)
    assert iscale.get_scaling_factor(
        m.fs.cathode_heat.control_volume.heat) == \
        pytest.approx(1e-8, rel=1e-5)
    assert iscale.get_scaling_factor(
        m.fs.cathode_HRSG.control_volume.heat) == \
        pytest.approx(1e-6, rel=1e-5)
    assert iscale.get_scaling_factor(
        m.fs.intercooler_s1.control_volume.heat) == \
        pytest.approx(1e-6, rel=1e-5)
    assert iscale.get_scaling_factor(
        m.fs.intercooler_s2.control_volume.heat) == \
        pytest.approx(1e-6, rel=1e-5)
    assert iscale.get_scaling_factor(
        m.fs.anode_HRSG.control_volume.heat) == \
        pytest.approx(1e-8, rel=1e-5)
    assert iscale.get_scaling_factor(
        m.fs.prereformer.control_volume.heat) == \
        pytest.approx(1e-6, rel=1e-5)
    assert iscale.get_scaling_factor(
        m.fs.reformer.control_volume.heat) == \
        pytest.approx(1e-6, rel=1e-5)
    assert iscale.get_scaling_factor(
        m.fs.reformer_recuperator.shell.heat) == \
        pytest.approx(1e-6, rel=1e-5)
    assert iscale.get_scaling_factor(
        m.fs.reformer_recuperator.tube.heat) == \
        pytest.approx(1e-6, rel=1e-5)

    # work
    assert iscale.get_scaling_factor(
        m.fs.anode_blower.control_volume.work) == \
        pytest.approx(1e-5, rel=1e-5)
    assert iscale.get_scaling_factor(
        m.fs.air_blower.control_volume.work) == \
        pytest.approx(1e-6, rel=1e-5)
    assert iscale.get_scaling_factor(
        m.fs.cathode_blower.control_volume.work) == \
        pytest.approx(1e-5, rel=1e-5)
    assert iscale.get_scaling_factor(
        m.fs.air_compressor_s1.control_volume.work) == \
        pytest.approx(1e-6, rel=1e-5)
    assert iscale.get_scaling_factor(
        m.fs.air_compressor_s2.control_volume.work) == \
        pytest.approx(1e-6, rel=1e-5)
    assert iscale.get_scaling_factor(
        m.fs.cathode_expander.control_volume.work) == \
        pytest.approx(1e-6, rel=1e-5)
    assert iscale.get_scaling_factor(
        m.fs.combustor_expander.control_volume.work) == \
        pytest.approx(1e-6, rel=1e-5)
    assert iscale.get_scaling_factor(
        m.fs.NG_expander.control_volume.work) == \
        pytest.approx(1e-6, rel=1e-5)

    # reaction extents
    assert iscale.get_scaling_factor(
        m.fs.combustor.control_volume.rate_reaction_extent[0, "h2_cmb"]) == \
        pytest.approx(1e2, rel=1e-5)
    assert iscale.get_scaling_factor(
        m.fs.combustor.control_volume.rate_reaction_extent[0, "co_cmb"]) == \
        pytest.approx(1e2, rel=1e-5)
    assert iscale.get_scaling_factor(
        m.fs.combustor.control_volume.rate_reaction_extent[0, "ch4_cmb"]) == \
        pytest.approx(1e5, rel=1e-5)


@pytest.mark.integration
def test_initialize_power_island(m):
    solver.options = {
        "max_iter": 500,
        "tol": 1e-7,
        "bound_push": 1e-12,
        "linear_solver": "ma57",
        "ma57_pivtol": 1e-3,
          }
    initialize_power_island(m)

    assert pyo.value(m.fs.cathode.inlet.flow_mol[0]) == \
        pytest.approx(34168, rel=1e-3)
    assert pyo.value(m.fs.cathode.inlet.temperature[0]) == \
        pytest.approx(892, rel=1e-3)
    assert pyo.value(m.fs.cathode.inlet.pressure[0]) == \
        pytest.approx(105490, rel=1e-3)
    assert pyo.value(m.fs.cathode.inlet.mole_frac_comp[0, "H2O"]) == \
        pytest.approx(0.0109, rel=1e-3)
    assert pyo.value(m.fs.cathode.inlet.mole_frac_comp[0, "CO2"]) == \
        pytest.approx(0.0003073, rel=1e-3)
    assert pyo.value(m.fs.cathode.inlet.mole_frac_comp[0, "N2"]) == \
        pytest.approx(0.8099, rel=1e-3)
    assert pyo.value(m.fs.cathode.inlet.mole_frac_comp[0, "O2"]) == \
        pytest.approx(0.1690, rel=1e-3)
    assert pyo.value(m.fs.cathode.inlet.mole_frac_comp[0, "Ar"]) == \
        pytest.approx(0.009880, rel=1e-3)

    assert pyo.value(m.fs.anode.inlet.flow_mol[0]) == \
        pytest.approx(5898, rel=1e-3)
    assert pyo.value(m.fs.anode.inlet.temperature[0]) == \
        pytest.approx(1413, rel=1e-3)
    assert pyo.value(m.fs.anode.inlet.pressure[0]) == \
        pytest.approx(186611, rel=1e-3)
    assert pyo.value(m.fs.anode.inlet.mole_frac_comp[0, "CH4"]) == \
        pytest.approx(0.2382, rel=1e-3)
    assert pyo.value(m.fs.anode.inlet.mole_frac_comp[0, "CO"]) == \
        pytest.approx(0.3772, rel=1e-3)
    assert pyo.value(m.fs.anode.inlet.mole_frac_comp[0, "CO2"]) == \
        pytest.approx(0.2585, rel=1e-3)
    assert pyo.value(m.fs.anode.inlet.mole_frac_comp[0, "H2"]) == \
        pytest.approx(0.3433, rel=1e-3)
    assert pyo.value(m.fs.anode.inlet.mole_frac_comp[0, "H2O"]) == \
        pytest.approx(0.2422, rel=1e-3)
    assert pyo.value(m.fs.anode.inlet.mole_frac_comp[0, "N2"]) == \
        pytest.approx(0.2880, rel=1e-3)
    assert pyo.value(m.fs.anode.inlet.mole_frac_comp[0, "O2"]) == \
        pytest.approx(0.3086, rel=1e-3)
    assert pyo.value(m.fs.anode.inlet.mole_frac_comp[0, "Ar"]) == \
        pytest.approx(0.2871, rel=1e-3)

    assert pyo.value(m.fs.anode.lagrange_mult[0, "C"]) == \
        pytest.approx(26137, rel=1e-3)
    assert pyo.value(m.fs.anode.lagrange_mult[0, "H"]) == \
        pytest.approx(75542, rel=1e-3)
    assert pyo.value(m.fs.anode.lagrange_mult[0, "O"]) == \
        pytest.approx(304676, rel=1e-3)

    assert pyo.value(m.fs.anode.outlet.mole_frac_comp[0, "O2"]) == \
        pytest.approx(8.8565e-12, rel=1e-3)

    assert pyo.value(m.fs.prereformer.gibbs_scaling) == \
        pytest.approx(1e-4, rel=1e-3)

    assert pyo.value(m.fs.prereformer.lagrange_mult[0, "C"]) == \
        pytest.approx(7764, rel=1e-3)
    assert pyo.value(m.fs.prereformer.lagrange_mult[0, "H"]) == \
        pytest.approx(66902, rel=1e-3)
    assert pyo.value(m.fs.prereformer.lagrange_mult[0, "O"]) == \
        pytest.approx(301919, rel=1e-3)

    assert pyo.value(m.fs.prereformer.outlet.mole_frac_comp[0, "O2"]) == \
        pytest.approx(0, 1e-3)
    assert pyo.value(m.fs.prereformer.outlet.mole_frac_comp[0, "Ar"]) == \
        pytest.approx(0.07344, rel=1e-3)
    assert pyo.value(m.fs.prereformer.outlet.mole_frac_comp[0, "C2H6"]) == \
        pytest.approx(7.0975e-7, rel=1e-3)
    assert pyo.value(m.fs.prereformer.outlet.mole_frac_comp[0, "C3H8"]) == \
        pytest.approx(0, 1e-3)
    assert pyo.value(m.fs.prereformer.outlet.mole_frac_comp[0, "C4H10"]) == \
        pytest.approx(0, 1e-3)


@pytest.mark.integration
def test_initialize_reformer(m):
    solver.options = {
        "max_iter": 500,
        "tol": 1e-7,
        "bound_push": 1e-12,
        "linear_solver": "ma57",
        "ma57_pivtol": 1e-3,
          }
    initialize_reformer(m)

    assert pyo.value(m.fs.reformer.inlet.flow_mol[0]) == \
        pytest.approx(3879, rel=1e-3)  # mol/s
    assert pyo.value(m.fs.reformer.inlet.temperature[0]) == \
        pytest.approx(1251, rel=1e-3)  # K
    assert pyo.value(m.fs.reformer.inlet.pressure[0]) == \
        pytest.approx(5555795, rel=1e-3)  # Pa
    assert pyo.value(m.fs.reformer.inlet.mole_frac_comp[0, "CH4"]) == \
        pytest.approx(0.47992, rel=1e-3)
    assert pyo.value(m.fs.reformer.inlet.mole_frac_comp[0, "C2H6"]) == \
        pytest.approx(0.46836, rel=1e-3)
    assert pyo.value(m.fs.reformer.inlet.mole_frac_comp[0, "C3H8"]) == \
        pytest.approx(0.45690, rel=1e-3)
    assert pyo.value(m.fs.reformer.inlet.mole_frac_comp[0, "C4H10"]) == \
        pytest.approx(0.44557, rel=1e-3)
    assert pyo.value(m.fs.reformer.inlet.mole_frac_comp[0, "H2"]) == \
        pytest.approx(0.49154, rel=1e-3)
    assert pyo.value(m.fs.reformer.inlet.mole_frac_comp[0, "CO"]) == \
        pytest.approx(0.47632, rel=1e-3)
    assert pyo.value(m.fs.reformer.inlet.mole_frac_comp[0, "CO2"]) == \
        pytest.approx(0.45510, rel=1e-3)
    assert pyo.value(m.fs.reformer.inlet.mole_frac_comp[0, "H2O"]) == \
        pytest.approx(0.47011, rel=1e-3)
    assert pyo.value(m.fs.reformer.inlet.mole_frac_comp[0, "N2"]) == \
        pytest.approx(0.52541, rel=1e-3)
    assert pyo.value(m.fs.reformer.inlet.mole_frac_comp[0, "O2"]) == \
        pytest.approx(0.45772, rel=1e-3)
    assert pyo.value(m.fs.reformer.inlet.mole_frac_comp[0, "Ar"]) == \
        pytest.approx(0.52542, rel=1e-3)

    assert pyo.value(m.fs.reformer.lagrange_mult[0, "C"]) == \
        pytest.approx(-25493, rel=1e-3)
    assert pyo.value(m.fs.reformer.lagrange_mult[0, "H"]) == \
        pytest.approx(83053, rel=1e-3)
    assert pyo.value(m.fs.reformer.lagrange_mult[0, "O"]) == \
        pytest.approx(370367, rel=1e-3)

    assert pyo.value(m.fs.reformer.outlet.mole_frac_comp[0, "O2"]) == \
        pytest.approx(0, rel=1e-3)
    assert pyo.value(m.fs.reformer.outlet.mole_frac_comp[0, "Ar"]) == \
        pytest.approx(0.062046, rel=1e-3)
    assert pyo.value(m.fs.reformer.outlet.mole_frac_comp[0, "CH4"]) == \
        pytest.approx(0.32244, rel=1e-3)
    assert pyo.value(m.fs.reformer.outlet.mole_frac_comp[0, "C2H6"]) == \
        pytest.approx(0.000140979, rel=1e-3)
    assert pyo.value(m.fs.reformer.outlet.mole_frac_comp[0, "C3H8"]) == \
        pytest.approx(2.2842e-7, rel=1e-3)
    assert pyo.value(m.fs.reformer.outlet.mole_frac_comp[0, "C4H10"]) == \
        pytest.approx(3.460668e-10, rel=1e-3)

    assert pyo.value(m.fs.reformer.lagrange_mult[(0, "H")]) == \
        pytest.approx(83053, rel=1e-3)
    assert pyo.value(m.fs.prereformer.lagrange_mult[(0, "H")]) == \
        pytest.approx(66902, rel=1e-3)
    assert pyo.value(m.fs.anode.lagrange_mult[(0, "H")]) == \
        pytest.approx(75542, rel=1e-3)
    assert pyo.value(m.fs.bypass_rejoin.outlet.temperature[0]) == \
        pytest.approx(1019, rel=1e-3)
    assert pyo.value(m.fs.anode_hx.shell_outlet.temperature[0]) == \
        pytest.approx(869, rel=1e-3)
    assert pyo.value(m.fs.cathode_hx.shell_outlet.temperature[0]) == \
        pytest.approx(475.8, rel=1e-3)


@pytest.mark.unit
def test_connect_reformer_to_power_island(m):
    connect_reformer_to_power_island(m)

    assert degrees_of_freedom(m) == 0

    assert not m.fs.anode_mix.feed.flow_mol[0].fixed


@pytest.mark.unit
def test_ROM(m):
    SOFC_ROM_setup(m)

    assert hasattr(m.fs, "ROM_fuel_inlet_temperature")
    assert hasattr(m.fs, "ROM_internal_reformation")
    assert m.fs.SOFC.current_density.fixed
    assert not m.fs.air_blower.inlet.flow_mol[0].fixed
    assert m.fs.SOFC.deltaT_cell.fixed
    assert hasattr(m.fs, "stack_power")


@pytest.mark.unit
def test_SOFC_energy_balance(m):
    add_SOFC_energy_balance(m)

    assert hasattr(m.fs, "SOFC_energy_balance")
    assert not m.fs.anode.outlet.temperature[0].fixed
    assert not m.fs.cathode_heat.outlet.temperature[0].fixed


@pytest.mark.unit
def test_add_result_constraints(m):
    add_result_constraints(m)

    assert hasattr(m.fs, "HRSG_heat_duty_constraint")
    assert hasattr(m.fs, "reformer_steam_heat_constraint")
    assert hasattr(m.fs, "steam_cycle_heat_constraint")
    assert hasattr(m.fs, "steam_cycle_power_constraint")
    assert hasattr(m.fs, "gross_power_constraint")
    assert hasattr(m.fs, "auxiliary_load_constraint")
    assert hasattr(m.fs, "net_power_constraint")
    assert hasattr(m.fs, "efficiency_rule")
    assert hasattr(m.fs, "CO2_emission_constraint")

    assert degrees_of_freedom(m) == 0


@pytest.mark.integration
def test_solve(m):
    assert solver_available
    solver.options = {
        "max_iter": 100,
        "tol": 1e-7,
        "bound_push": 1e-12,
        "linear_solver": "ma57",
        "ma57_pivtol": 1e-3,
          }
    solve_iteration = 0
    for i in range(1, 10):  # keep looping until condition is met
        solve_iteration += 1
        print('Solve # ', solve_iteration)
        res = solver.solve(m, tee=True)
        if 'Optimal Solution Found' in res.solver.message:
            break

    assert check_optimal_termination(res)
    assert solve_iteration == 1  # this value may be updated as needed


@pytest.mark.unit
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


@pytest.mark.unit
def test_make_stream_dict(m):
    make_stream_dict(m)

    assert hasattr(m, "_streams")

    ordered_dict = OrderedDict(
        [
            ("FUEL_IN", m.fs.reformer_recuperator.tube_inlet),
            ("HOT_FUEL", m.fs.reformer_recuperator.tube_outlet),
            ("REF_IN", m.fs.reformer_mix.gas_inlet),
            ("AIR_IN", m.fs.reformer_mix.oxygen_inlet),
            ("STEAM_IN", m.fs.reformer_mix.steam_inlet),
            ("REF_OUT", m.fs.reformer.outlet),
            ("SYN_IN", m.fs.anode_mix.feed),
            ("ANO_HX_CI", m.fs.ANODE_MIXER),
            ("ANODE_IN", m.fs.fuel_cell_mix.fuel_inlet),
            ("ANODE_OUT", m.fs.anode.outlet),
            ("ANODE_REC", m.fs.ANODE_BLOWER),
            ("ANO_HX_HI", m.fs.ANODE_RECYCLE_HX),
            ("ANO_HX_HO", m.fs.anode_hx.shell_outlet),
            ("AIR", m.fs.air_blower.inlet),
            ("CATH_HX_CI", m.fs.AIR_BLOWER),
            ("CATH_HX_CO", m.fs.CATHODE_HX_COLD),
            ("CATH_IN", m.fs.CATHODE_MIXER),
            ("CATH_OUT", m.fs.CATHODE_HEAT),
            ("CATH_REC", m.fs.CATHODE_BLOWER),
            ("CATH_HX_HI", m.fs.CATHODE_RECYCLE_HX),
            ("CATH_HX_HO", m.fs.cathode_hx.shell_outlet),
            ("CATH_HRSG_IN", m.fs.cathode_HRSG.inlet),
            ("CATH_EXH", m.fs.cathode_HRSG.outlet),
            ("COMB_AIR", m.fs.combustor_mix.cathode_inlet),
            ("COMB_OUT", m.fs.combustor.outlet),
            ("ANOD_EXH", m.fs.anode_HRSG.outlet),
        ]
    )

    assert m._streams == ordered_dict


@pytest.mark.unit
def test_pfd_result(m):
    df = create_stream_table_dataframe(streams=m._streams, orient="index")
    pfd_result("results.svg", m, df)

    # check that a results file was created in the test directory
    # assert os.path.exists(os.path.join(this_file_dir(), "results.svg")) - file not created on CI client, how can we do this?
    # remove the file to keep the test directory clean
    # os.remove(os.path.join(this_file_dir(), "results.svg"))
    # assert not os.path.exists(os.path.join(this_file_dir(), "results.svg"))
