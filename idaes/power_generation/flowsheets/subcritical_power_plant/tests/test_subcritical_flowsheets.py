# Import Pyomo libraries
import pyomo.environ as pyo
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.generic_models.properties import iapws95
# import idaes.logger as idaeslog
from idaes.core.util.dyn_utils import (
    copy_values_at_time,
    copy_non_time_indexed_values)

import idaes.power_generation.flowsheets.subcritical_power_plant.subcritical_power_plant as subcrit_plant
import idaes.power_generation.flowsheets.subcritical_power_plant.steam_cycle_flowsheet as steam_cycle
import idaes.power_generation.flowsheets.subcritical_power_plant.subcritical_boiler_flowsheet as blr

import pytest

__author__ = "Boiler Subsystem Team (J. Ma, M. Zamarripa)"

solver_available = pyo.SolverFactory('ipopt').available()
prop_available = iapws95.iapws95_available()


@pytest.mark.integration
@pytest.mark.solver
@pytest.mark.skipif(not prop_available, reason="IAPWS not available")
@pytest.mark.skipif(not solver_available, reason="Solver not available")
def test_subcritical_boiler():
    m = blr.main_steady_state()
    assert degrees_of_freedom(m) == 0
    # mass balance
    assert (pytest.approx(0, abs=1e-3) ==
            pyo.value(m.fs_main.fs_blr.aRH1.tube_inlet.flow_mol[0]
                      + m.fs_main.fs_blr.aECON.tube_inlet.flow_mol[0]
                      - m.fs_main.fs_blr.aPlaten.outlet.flow_mol[0]
                      - m.fs_main.fs_blr.aRH2.tube_outlet.flow_mol[0]
                      - m.fs_main.fs_blr.blowdown_split.FW_Blowdown.flow_mol[0]
                      + m.fs_main.fs_blr.Attemp.Water_inlet.flow_mol[0]
                      ))
    # FEGT temperature
    assert (pytest.approx(1399.9583, abs=1e-3) ==
            pyo.value(m.fs_main.fs_blr.aBoiler.flue_gas_outlet.temperature[0]))
    assert (pytest.approx(329148487.5973, abs=1e-3) ==
            pyo.value(m.fs_main.fs_blr.aBoiler.heat_total_ww[0]))
    assert (pytest.approx(418839856.6386, abs=1e-3) ==
            pyo.value(m.fs_main.fs_blr.aBoiler.heat_total[0]))


@pytest.mark.integration
@pytest.mark.solver
@pytest.mark.skipif(not prop_available, reason="IAPWS not available")
@pytest.mark.skipif(not solver_available, reason="Solver not available")
def test_subcritical_boiler_dynamic():
    m = blr.main_dynamic()
    assert degrees_of_freedom(m) == 0
    assert (pytest.approx(1399.9583, abs=1e-3) ==
            pyo.value(m.fs_main.fs_blr.aBoiler.flue_gas_outlet.temperature[0]))
    assert (pytest.approx(329148487.53436375, abs=1e-3) ==
            pyo.value(m.fs_main.fs_blr.aBoiler.heat_total_ww[0]))
    assert (pytest.approx(418839856.6744472, abs=1e-3) ==
            pyo.value(m.fs_main.fs_blr.aBoiler.heat_total[0]))
    assert (pytest.approx(1408.9636, abs=1e-3) ==
            pyo.value(m.fs_main.fs_blr.aBoiler.flue_gas_outlet.
                      temperature[60]))
    assert (pytest.approx(334512567.44442254, abs=1e-3) ==
            pyo.value(m.fs_main.fs_blr.aBoiler.heat_total_ww[60]))
    assert (pytest.approx(426313180.661969, abs=1e-3) ==
            pyo.value(m.fs_main.fs_blr.aBoiler.heat_total[60]))


@pytest.mark.integration
@pytest.mark.solver
@pytest.mark.skipif(not prop_available, reason="IAPWS not available")
@pytest.mark.skipif(not solver_available, reason="Solver not available")
def test_steam_cycle():
    m = steam_cycle.main_steady_state()
    assert degrees_of_freedom(m) == 0
    assert (pytest.approx(61958, abs=1e-3) ==
            pyo.value(m.fs_main.fs_stc.turb.inlet_split.inlet.enth_mol[0]))
    assert (pytest.approx(12473.27146, abs=1e-3) ==
            pyo.value(m.fs_main.fs_stc.turb.inlet_split.inlet.flow_mol[0]))
    assert (pytest.approx(266.8806, abs=1e-3) ==
            pyo.value(m.fs_main.fs_stc.power_output[0]))
    assert (pytest.approx(0, abs=1e-3) ==
            pyo.value(m.fs_main.fs_stc.turb.inlet_split.inlet.flow_mol[0] # turbine inlet
                      - m.fs_main.fs_stc.turb.hp_split[14].outlet_1.flow_mol[0] # out to reheat
                      + m.fs_main.fs_stc.turb.ip_stages[1].inlet.flow_mol[0] # in from reheat
                      - m.fs_main.fs_stc.spray_valve.outlet.flow_mol[0] # out to attemperator
                      - m.fs_main.fs_stc.fwh6.desuperheat.outlet_2.flow_mol[0] # out to economizer
                      + m.fs_main.fs_stc.makeup_valve.inlet.flow_mol[0] # in from makeup
                      ))
