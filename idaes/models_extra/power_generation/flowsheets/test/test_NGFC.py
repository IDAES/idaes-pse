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

from io import StringIO
import sys
import os
from collections import OrderedDict

import pytest
import pyomo.environ as pyo
from pyomo.environ import check_optimal_termination
from pyomo.common.fileutils import this_file_dir
from pyomo.common.tempfiles import TempfileManager

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
        main,
        )

from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.solvers import get_solver
from idaes.core.util import scaling as iscale
from idaes.core.util.scaling import (extreme_jacobian_columns,
                                     extreme_jacobian_rows)
from idaes.core.util import model_serializer as ms

from idaes.core.util.tables import create_stream_table_dataframe


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


@pytest.mark.component
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


@pytest.mark.component
def test_initialize_power_island(m):
    initialize_power_island(m)

    # a few checks to make sure this method worked as expected
    assert pyo.value(m.fs.cathode.inlet.flow_mol[0]) == \
        pytest.approx(34168, rel=1e-3)
    assert pyo.value(m.fs.cathode.inlet.temperature[0]) == \
        pytest.approx(892, rel=1e-3)
    assert pyo.value(m.fs.cathode.inlet.pressure[0]) == \
        pytest.approx(105490, rel=1e-3)


@pytest.mark.component
def test_initialize_reformer(m):
    initialize_reformer(m)

    # a few checks to make sure this method worked as expected
    assert pyo.value(m.fs.reformer.inlet.flow_mol[0]) == \
        pytest.approx(3880, rel=1e-3)  # mol/s
    assert pyo.value(m.fs.reformer.inlet.temperature[0]) == \
        pytest.approx(1387, rel=1e-3)  # K
    assert pyo.value(m.fs.reformer.inlet.pressure[0]) == \
        pytest.approx(50025000, rel=1e-3)  # Pa

    assert pyo.value(m.fs.reformer.lagrange_mult[(0, "H")]) == \
        pytest.approx(83053, rel=1e-3)
    assert pyo.value(m.fs.prereformer.lagrange_mult[(0, "H")]) == \
        pytest.approx(69351, rel=1e-3)
    assert pyo.value(m.fs.anode.lagrange_mult[(0, "H")]) == \
        pytest.approx(74261, rel=1e-3)
    assert pyo.value(m.fs.bypass_rejoin.outlet.temperature[0]) == \
        pytest.approx(1019, rel=1e-3)
    assert pyo.value(m.fs.anode_hx.shell_outlet.temperature[0]) == \
        pytest.approx(1361, rel=1e-3)
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


@pytest.mark.unit
def test_json_load_init(m):
    fname = os.path.join(this_file_dir(), "NGFC_flowsheet_init.json.gz")

    ms.from_json(m, fname=fname)

    # check values
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
        pytest.approx(7256, rel=1e-3)
    assert pyo.value(m.fs.anode.inlet.temperature[0]) == \
        pytest.approx(1387, rel=1e-3)
    assert pyo.value(m.fs.anode.inlet.pressure[0]) == \
        pytest.approx(50025000, rel=1e-3)
    assert pyo.value(m.fs.anode.inlet.mole_frac_comp[0, "CH4"]) == \
        pytest.approx(0.5882, rel=1e-3)
    assert pyo.value(m.fs.anode.inlet.mole_frac_comp[0, "CO"]) == \
        pytest.approx(0.5110, rel=1e-3)
    assert pyo.value(m.fs.anode.inlet.mole_frac_comp[0, "CO2"]) == \
        pytest.approx(0.4126, rel=1e-3)
    assert pyo.value(m.fs.anode.inlet.mole_frac_comp[0, "H2"]) == \
        pytest.approx(0.4899, rel=1e-3)
    assert pyo.value(m.fs.anode.inlet.mole_frac_comp[0, "H2O"]) == \
        pytest.approx(0.3938, rel=1e-3)
    assert pyo.value(m.fs.anode.inlet.mole_frac_comp[0, "N2"]) == \
        pytest.approx(0.4898, rel=1e-3)
    assert pyo.value(m.fs.anode.inlet.mole_frac_comp[0, "O2"]) == \
        pytest.approx(0.3234, rel=1e-3)
    assert pyo.value(m.fs.anode.inlet.mole_frac_comp[0, "Ar"]) == \
        pytest.approx(0.4898, rel=1e-3)

    assert pyo.value(m.fs.anode.lagrange_mult[0, "C"]) == \
        pytest.approx(16759, rel=1e-3)
    assert pyo.value(m.fs.anode.lagrange_mult[0, "H"]) == \
        pytest.approx(74261, rel=1e-3)
    assert pyo.value(m.fs.anode.lagrange_mult[0, "O"]) == \
        pytest.approx(311619, rel=1e-3)

    assert pyo.value(m.fs.anode.outlet.mole_frac_comp[0, "O2"]) == \
        pytest.approx(1.4486e-22, rel=1e-3)

    assert pyo.value(m.fs.prereformer.gibbs_scaling) == \
        pytest.approx(1e-4, rel=1e-3)

    assert pyo.value(m.fs.prereformer.lagrange_mult[0, "C"]) == \
        pytest.approx(5615, rel=1e-3)
    assert pyo.value(m.fs.prereformer.lagrange_mult[0, "H"]) == \
        pytest.approx(69351, rel=1e-3)
    assert pyo.value(m.fs.prereformer.lagrange_mult[0, "O"]) == \
        pytest.approx(308901, rel=1e-3)

    assert pyo.value(m.fs.reformer.lagrange_mult[0, "C"]) == \
        pytest.approx(-25493, rel=1e-3)
    assert pyo.value(m.fs.reformer.lagrange_mult[0, "H"]) == \
        pytest.approx(83053, rel=1e-3)
    assert pyo.value(m.fs.reformer.lagrange_mult[0, "O"]) == \
        pytest.approx(370367, rel=1e-3)


@pytest.mark.unit
def test_json_load_solution(m):
    fname = os.path.join(this_file_dir(), "NGFC_flowsheet_solution.json.gz")

    ms.from_json(m, fname=fname)

    assert pyo.value(m.fs.cathode.ion_outlet.flow_mol[0]) == \
        pytest.approx(1670, rel=1e-3)
    assert pyo.value(m.fs.reformer_recuperator.area) == \
        pytest.approx(4511, rel=1e-3)
    assert pyo.value(m.fs.anode.heat_duty[0]) == \
        pytest.approx(-672918765, rel=1e-3)
    assert pyo.value(m.fs.CO2_emissions) == \
        pytest.approx(291, rel=1e-3)
    assert pyo.value(m.fs.net_power) == \
        pytest.approx(660, rel=1e-3)


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

    with TempfileManager.new_context() as tf:
        dname = tf.mkdtemp()
        assert os.path.isdir(dname)
        pfd_result(os.path.join(dname, "results.svg"), m, df)

        # check that a results file was created in the results directory
        assert os.path.isfile(os.path.join(dname, "results.svg"))
    # check that the temporary directory and file are cleaned up
    assert not os.path.isdir(dname)


@pytest.mark.integration
def test_main_init_and_solve(m):
    # remove both json, have main() reinitialize and resolve/serialize

    jsontestdir = this_file_dir()
    assert os.path.exists(os.path.join(jsontestdir,
                                       "NGFC_flowsheet_init.json.gz"))
    assert os.path.exists(os.path.join(jsontestdir,
                                       "NGFC_flowsheet_solution.json.gz"))

    # remove both json so we can reinitialize/resolve and regenerate them
    os.remove(os.path.join(jsontestdir, "NGFC_flowsheet_init.json.gz"))
    os.remove(os.path.join(jsontestdir, "NGFC_flowsheet_solution.json.gz"))
    assert not os.path.exists(os.path.join(jsontestdir,
                                           "NGFC_flowsheet_init.json.gz"))
    assert not os.path.exists(os.path.join(jsontestdir,
                                           "NGFC_flowsheet_solution.json.gz"))

    with TempfileManager.new_context() as tf:
        dname = tf.mkdtemp()
        assert os.path.isdir(dname)

        stream = StringIO()
        sys.stdout = stream
        m = main(dname, jsontestdir)
        sys.stdout = sys.__stdout__

        output1 = """Scaling flowsheet variables
overwriting mole_frac lower bound, set to 0 to remove warnings
Scaling flowsheet constraints
Calculating scaling factors
"""

        output2 = """Starting ROM initialization
ROM initialization completed
"""
        output3 = """EXIT: Optimal Solution Found.
PFD Results Created
"""

        assert isinstance(m, pyo.ConcreteModel)
        assert output1 in stream.getvalue()  # check pre-init printed output
        assert output2 in stream.getvalue()  # check the ROM solves after init
        # split output1 and output2 because initialization prints statuses
        assert output3 in stream.getvalue()  # check that final solve is
        # optimal and results are exported post-solve
        # checking together so internal statuses from init are not flagged

        # check that a results file was created in the results directory
        assert os.path.isfile(os.path.join(dname, "NGFC_results.svg"))
    # check that the temporary directory and file are cleaned up
    assert not os.path.isdir(dname)


@pytest.mark.integration
def test_main_load_and_solve(m):
    # remove solution json, have main() load init json and resolve/serialize

    jsontestdir = this_file_dir()
    assert os.path.exists(os.path.join(jsontestdir,
                                       "NGFC_flowsheet_init.json.gz"))
    assert os.path.exists(os.path.join(jsontestdir,
                                       "NGFC_flowsheet_solution.json.gz"))

    # remove solution json so we can resolve and regenerate it
    os.remove(os.path.join(jsontestdir, "NGFC_flowsheet_solution.json.gz"))
    assert not os.path.exists(os.path.join(jsontestdir,
                                           "NGFC_flowsheet_solution.json.gz"))

    with TempfileManager.new_context() as tf:
        dname = tf.mkdtemp()
        assert os.path.isdir(dname)

        stream = StringIO()
        sys.stdout = stream
        m = main(dname, jsontestdir)
        sys.stdout = sys.__stdout__

        output1 = """Scaling flowsheet variables
overwriting mole_frac lower bound, set to 0 to remove warnings
Scaling flowsheet constraints
Calculating scaling factors

Starting ROM initialization
ROM initialization completed
Loading initialized model
"""
        output2 = """EXIT: Optimal Solution Found.
PFD Results Created
"""

        assert isinstance(m, pyo.ConcreteModel)
        assert output1 in stream.getvalue()  # check pre-solve printed output
        # check entire string because initialization is skipped here
        assert output2 in stream.getvalue()  # check that final solve is
        # optimal and results are exported post-solve
        # checking together so internal statuses from init are not flagged

        # check that a results file was created in the results directory
        assert os.path.isfile(os.path.join(dname, "NGFC_results.svg"))
    # check that the temporary directory and file are cleaned up
    assert not os.path.isdir(dname)


@pytest.mark.unit
def test_main_load_json():
    # load solution json, do not call solver

    jsontestdir = this_file_dir()
    assert os.path.exists(os.path.join(jsontestdir,
                                       "NGFC_flowsheet_init.json.gz"))
    assert os.path.exists(os.path.join(jsontestdir,
                                       "NGFC_flowsheet_solution.json.gz"))

    with TempfileManager.new_context() as tf:
        dname = tf.mkdtemp()
        assert os.path.isdir(dname)

        stream = StringIO()
        sys.stdout = stream
        m = main(dname, jsontestdir)
        sys.stdout = sys.__stdout__

        output = """Scaling flowsheet variables
overwriting mole_frac lower bound, set to 0 to remove warnings
Scaling flowsheet constraints
Calculating scaling factors

Starting ROM initialization
ROM initialization completed
Loading solved model
PFD Results Created
"""

        assert isinstance(m, pyo.ConcreteModel)
        assert output in stream.getvalue()  # check printed output
        # a single string because no solver steps are called

        # check that a results file was created in the results directory
        assert os.path.isfile(os.path.join(dname, "NGFC_results.svg"))
    # check that the temporary directory and file are cleaned up
    assert not os.path.isdir(dname)
