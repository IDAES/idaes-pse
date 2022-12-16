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
# Import Pyomo libraries
import pyomo.environ as pyo
from idaes.core.util.model_statistics import degrees_of_freedom

import idaes.models_extra.power_generation.flowsheets.subcritical_power_plant.subcritical_power_plant as subcrit_plant
import idaes.models_extra.power_generation.flowsheets.subcritical_power_plant.steam_cycle_flowsheet as steam_cycle
import idaes.models_extra.power_generation.flowsheets.subcritical_power_plant.subcritical_boiler_flowsheet as blr
import idaes.models_extra.power_generation.flowsheets.subcritical_power_plant.subcritical_boiler as recyrc
import pytest

__author__ = "Boiler Subsystem Team (J. Ma, M. Zamarripa)"


@pytest.mark.component
def test_subcritical_boiler_ss_build():
    m = blr.get_model(dynamic=False, init=False)
    assert degrees_of_freedom(m) == 12


@pytest.mark.component
def test_subcritical_boiler_dynamic_build():
    m = blr.get_model(dynamic=True, init=False)
    assert degrees_of_freedom(m) == 223


@pytest.mark.integration
def test_subcritical_boiler():
    m = blr.main_steady_state()
    assert degrees_of_freedom(m) == 0
    # mass balance
    assert pytest.approx(0, abs=1e-2) == pyo.value(
        m.fs_main.fs_blr.aRH1.tube_inlet.flow_mol[0]
        + m.fs_main.fs_blr.aECON.tube_inlet.flow_mol[0]
        - m.fs_main.fs_blr.aPlaten.outlet.flow_mol[0]
        - m.fs_main.fs_blr.aRH2.tube_outlet.flow_mol[0]
        - m.fs_main.fs_blr.blowdown_split.FW_Blowdown.flow_mol[0]
        + m.fs_main.fs_blr.Attemp.Water_inlet.flow_mol[0]
    )
    # FEGT temperature
    assert pytest.approx(1399.9583, rel=1e-5) == pyo.value(
        m.fs_main.fs_blr.aBoiler.flue_gas_outlet.temperature[0]
    )
    assert pytest.approx(329148469.4485, rel=1e-5) == pyo.value(
        m.fs_main.fs_blr.aBoiler.heat_total_ww[0]
    )
    assert pytest.approx(418839867.42785, rel=1e-5) == pyo.value(
        m.fs_main.fs_blr.aBoiler.heat_total[0]
    )


@pytest.mark.integration
def test_subcritical_boiler_dynamic():
    m = blr.main_dynamic()
    assert degrees_of_freedom(m) == 0
    assert pytest.approx(1399.9583, rel=1e-5) == pyo.value(
        m.fs_main.fs_blr.aBoiler.flue_gas_outlet.temperature[0]
    )
    assert pytest.approx(329148469.4568, rel=1e-5) == pyo.value(
        m.fs_main.fs_blr.aBoiler.heat_total_ww[0]
    )
    assert pytest.approx(418839867.42785, rel=1e-5) == pyo.value(
        m.fs_main.fs_blr.aBoiler.heat_total[0]
    )
    assert pytest.approx(1408.9636, rel=1e-5) == pyo.value(
        m.fs_main.fs_blr.aBoiler.flue_gas_outlet.temperature[60]
    )
    assert pytest.approx(334512549.5052, rel=1e-5) == pyo.value(
        m.fs_main.fs_blr.aBoiler.heat_total_ww[60]
    )
    assert pytest.approx(426313190.1620, rel=1e-5) == pyo.value(
        m.fs_main.fs_blr.aBoiler.heat_total[60]
    )


@pytest.mark.integration
def test_steam_cycle():
    m = steam_cycle.main_steady_state()
    assert degrees_of_freedom(m) == 0
    assert pytest.approx(61958, rel=1e-5) == pyo.value(
        m.fs_main.fs_stc.turb.inlet_split.inlet.enth_mol[0]
    )
    assert pytest.approx(12473.27146, rel=1e-5) == pyo.value(
        m.fs_main.fs_stc.turb.inlet_split.inlet.flow_mol[0]
    )
    assert pytest.approx(266.8806, rel=1e-5) == pyo.value(
        m.fs_main.fs_stc.power_output[0]
    )
    assert pytest.approx(0, abs=1e-2) == pyo.value(
        m.fs_main.fs_stc.turb.inlet_split.inlet.flow_mol[0]  # turbine inlet
        - m.fs_main.fs_stc.turb.hp_split[14].outlet_1.flow_mol[0]  # out to reheat
        + m.fs_main.fs_stc.turb.ip_stages[1].inlet.flow_mol[0]  # in from reheat
        - m.fs_main.fs_stc.spray_valve.outlet.flow_mol[0]  # out to attemperator
        - m.fs_main.fs_stc.fwh6.desuperheat.cold_side_outlet.flow_mol[
            0
        ]  # out to economizer
        + m.fs_main.fs_stc.makeup_valve.inlet.flow_mol[0]  # in from makeup
    )


@pytest.mark.integration
def test_subc_power_plant():
    m = subcrit_plant.main_steady_state()
    assert degrees_of_freedom(m) == 0
    assert pytest.approx(61634.3740, rel=1e-5) == pyo.value(
        m.fs_main.fs_stc.turb.inlet_split.inlet.enth_mol[0]
    )
    assert pytest.approx(14908.39189, rel=1e-5) == pyo.value(
        m.fs_main.fs_stc.turb.inlet_split.inlet.flow_mol[0]
    )
    assert pytest.approx(320, rel=1e-5) == pyo.value(m.fs_main.fs_stc.power_output[0])
    assert pytest.approx(0, abs=1e-2) == pyo.value(
        m.fs_main.fs_stc.turb.inlet_split.inlet.flow_mol[0]  # turbine inlet
        - m.fs_main.fs_stc.turb.hp_split[14].outlet_1.flow_mol[0]  # out to reheat
        + m.fs_main.fs_stc.turb.ip_stages[1].inlet.flow_mol[0]  # in from reheat
        - m.fs_main.fs_stc.spray_valve.outlet.flow_mol[0]  # out to attemperator
        - m.fs_main.fs_stc.fwh6.desuperheat.cold_side_outlet.flow_mol[
            0
        ]  # out to economizer
        + m.fs_main.fs_stc.makeup_valve.inlet.flow_mol[0]  # in from makeup
    )


@pytest.mark.component
def test_dynamic_power_plant_build():
    # constructing and initializing dynamic power plant
    # not solving due to simulation time >20 min
    m = subcrit_plant.get_model(dynamic=True, init=False)
    assert m.dynamic is True
    assert degrees_of_freedom(m) == 168


@pytest.mark.component
def test_steadystate_power_plant_build():
    # constructing and initializing dynamic power plant
    # not solving due to simulation time >20 min
    m = subcrit_plant.get_model(dynamic=False, init=False)
    assert m.dynamic is False
    assert degrees_of_freedom(m) == -5


@pytest.mark.component
def test_dynamic_steam_cycle():
    # constructing and initializing dynamic steam cycle flowsheet
    m = steam_cycle.get_model(dynamic=True)
    assert m.dynamic is True
    assert degrees_of_freedom(m) == 7


@pytest.mark.component
def test_subcritical_recirculation_system():
    m = recyrc.main()
    assert degrees_of_freedom(m) == 0
