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
Tests for IAPWS using a heater block.  This is an easy way to test the
different calculation method options.

Author: John Eslick
"""
import pytest

from pyomo.environ import ConcreteModel, value

from idaes.core import FlowsheetBlock
from idaes.generic_models.unit_models import Heater
from idaes.generic_models.properties import iapws95
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core import MaterialBalanceType
from idaes.core.util import get_solver

# Set module level pyest marker
pytestmark = pytest.mark.iapws
prop_available = iapws95.iapws95_available()

# -----------------------------------------------------------------------------
# set up solver
solver = get_solver()


@pytest.mark.unit
def test_heater_ph_mixed_byphase():
    """Test mixed phase form with P-H state vars and phase mass balances
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = iapws95.Iapws95ParameterBlock(
        default={})
    m.fs.heater = Heater(
        default={"property_package": m.fs.properties})
    m.fs.heater.inlet.enth_mol.fix(4000)
    m.fs.heater.inlet.flow_mol.fix(100)
    m.fs.heater.inlet.pressure.fix(101325)
    m.fs.heater.heat_duty[0].fix(100*20000)
    m.fs.heater.initialize()
    assert degrees_of_freedom(m) == 0
    solver.solve(m)
    prop_in = m.fs.heater.control_volume.properties_in[0]
    prop_out = m.fs.heater.control_volume.properties_out[0]
    assert value(prop_in.temperature) == pytest.approx(326.166707, rel=1e-5)
    assert value(prop_out.temperature) == pytest.approx(373.12429, rel=1e-5)
    assert value(prop_in.phase_frac["Liq"]) == pytest.approx(1, rel=1e-6)
    assert value(prop_out.phase_frac["Liq"]) == pytest.approx(
        0.5953219, rel=1e-5)
    assert value(prop_in.phase_frac["Vap"]) == pytest.approx(0, abs=1e-6)
    assert value(prop_out.phase_frac["Vap"]) == pytest.approx(
        0.4046781, rel=1e-4)


@pytest.mark.unit
def test_heater_phmixed_mixed_total():
    """Test mixed phase form with P-H state vars and total mass balances
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = iapws95.Iapws95ParameterBlock(
        default={})
    m.fs.heater = Heater(
        default={"property_package": m.fs.properties,
                 "material_balance_type": MaterialBalanceType.componentTotal})
    m.fs.heater.inlet.enth_mol.fix(4000)
    m.fs.heater.inlet.flow_mol.fix(100)
    m.fs.heater.inlet.pressure.fix(101325)
    m.fs.heater.heat_duty[0].fix(100*20000)
    m.fs.heater.initialize()
    assert degrees_of_freedom(m) == 0
    solver.solve(m)
    prop_in = m.fs.heater.control_volume.properties_in[0]
    prop_out = m.fs.heater.control_volume.properties_out[0]
    assert value(prop_in.temperature) == pytest.approx(326.166707, rel=1e-5)
    assert value(prop_out.temperature) == pytest.approx(373.12429, rel=1e-5)
    assert value(prop_in.phase_frac["Liq"]) == pytest.approx(1, rel=1e-6)
    assert value(prop_out.phase_frac["Liq"]) == pytest.approx(
        0.5953219, rel=1e-5)
    assert value(prop_in.phase_frac["Vap"]) == pytest.approx(0, abs=1e-5)
    assert value(prop_out.phase_frac["Vap"]) == pytest.approx(
        0.4046781, rel=1e-4)


@pytest.mark.unit
def test_heater_ph_lg_total():
    """Test liquid/vapor form with P-H state vars and total mass balances
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = iapws95.Iapws95ParameterBlock(
        default={"phase_presentation": iapws95.PhaseType.LG})
    m.fs.heater = Heater(
        default={"property_package": m.fs.properties,
                 "material_balance_type": MaterialBalanceType.componentTotal})
    m.fs.heater.inlet.enth_mol.fix(4000)
    m.fs.heater.inlet.flow_mol.fix(100)
    m.fs.heater.inlet.pressure.fix(101325)
    m.fs.heater.heat_duty[0].fix(100*20000)
    m.fs.heater.initialize()
    assert degrees_of_freedom(m) == 0
    solver.solve(m)
    prop_in = m.fs.heater.control_volume.properties_in[0]
    prop_out = m.fs.heater.control_volume.properties_out[0]
    assert value(prop_in.temperature) == pytest.approx(326.166707, rel=1e-5)
    assert value(prop_out.temperature) == pytest.approx(373.12429, rel=1e-5)
    assert value(prop_in.phase_frac["Liq"]) == pytest.approx(1, rel=1e-5)
    assert value(prop_out.phase_frac["Liq"]) == pytest.approx(
        0.5953219, rel=1e-5)
    assert value(prop_in.phase_frac["Vap"]) == pytest.approx(0, abs=1e-5)
    assert value(prop_out.phase_frac["Vap"]) == pytest.approx(
        0.4046781, rel=1e-4)


@pytest.mark.unit
def test_heater_ph_lg_phase():
    """Test liquid/vapor form with P-H state vars and phase mass balances

    This is known not to work so this is a place holder to demostrate the extra
    equation and reserve a place for a test once the problem is fixed. The
    issue in case I forget is that the phase based mass balance writes a mass
    balnace for each phase, but the P-H formulation always calculates a phase
    fraction making one equation redundant.
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = iapws95.Iapws95ParameterBlock(
        default={"phase_presentation": iapws95.PhaseType.LG})
    m.fs.heater = Heater(
        default={"property_package": m.fs.properties,
                 "material_balance_type": MaterialBalanceType.componentPhase})
    m.fs.heater.inlet.enth_mol.fix(4000)
    m.fs.heater.inlet.flow_mol.fix(100)
    m.fs.heater.inlet.pressure.fix(101325)
    m.fs.heater.heat_duty[0].fix(100*20000)
    assert degrees_of_freedom(m) == -1


@pytest.mark.unit
def test_heater_ph_l_phase_two():
    """Test liquid phase only form with P-H state vars and phase mass balances
    where the result is a two phase mixture.  For the P-H vars if there are
    two phases, using the liquid only properties will case the phase fraction
    to report all liquid even though really it is not.  Only use this option
    if you are sure you do not have a vapor phase.
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = iapws95.Iapws95ParameterBlock(
        default={"phase_presentation": iapws95.PhaseType.L})
    m.fs.heater = Heater(
        default={"property_package": m.fs.properties,
                 "material_balance_type": MaterialBalanceType.componentPhase})
    m.fs.heater.inlet.enth_mol.fix(3000)
    m.fs.heater.inlet.flow_mol.fix(100)
    m.fs.heater.inlet.pressure.fix(101325)
    m.fs.heater.heat_duty[0].fix(100*1000)
    m.fs.heater.initialize()
    assert degrees_of_freedom(m) == 0
    solver.solve(m)
    prop_in = m.fs.heater.control_volume.properties_in[0]
    prop_out = m.fs.heater.control_volume.properties_out[0]
    assert abs(value(prop_in.temperature) - 312.88896252921666) <= 1e-4
    assert abs(value(prop_out.temperature) - 326.166707507874) <= 1e-4
    assert abs(value(prop_in.phase_frac["Liq"]) - 1) <= 1e-6
    assert abs(value(prop_out.phase_frac["Liq"]) - 1) <= 1e-6
    assert abs(value(prop_in.phase_frac["Vap"]) - 0) <= 1e-6
    assert abs(value(prop_out.phase_frac["Vap"]) - 0) <= 1e-6


@pytest.mark.unit
def test_heater_ph_l_phase():
    """Test liquid phase only form with P-H state vars and phase mass balances
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = iapws95.Iapws95ParameterBlock(
        default={"phase_presentation": iapws95.PhaseType.L})
    m.fs.heater = Heater(
        default={"property_package": m.fs.properties,
                 "material_balance_type": MaterialBalanceType.componentPhase})
    m.fs.heater.inlet.enth_mol.fix(3000)
    m.fs.heater.inlet.flow_mol.fix(100)
    m.fs.heater.inlet.pressure.fix(101325)
    m.fs.heater.heat_duty[0].fix(100*1000)
    m.fs.heater.initialize()
    assert degrees_of_freedom(m) == 0
    solver.solve(m)
    prop_in = m.fs.heater.control_volume.properties_in[0]
    prop_out = m.fs.heater.control_volume.properties_out[0]
    assert abs(value(prop_in.temperature) - 312.88896252921666) <= 1e-4
    assert abs(value(prop_out.temperature) - 326.166707507874) <= 1e-4
    assert abs(value(prop_in.phase_frac["Liq"]) - 1) <= 1e-6
    assert abs(value(prop_out.phase_frac["Liq"]) - 1) <= 1e-6
    assert abs(value(prop_in.phase_frac["Vap"]) - 0) <= 1e-6
    assert abs(value(prop_out.phase_frac["Vap"]) - 0) <= 1e-6


@pytest.mark.unit
def test_heater_ph_g_phase():
    """Test vapor phase only form with P-H state vars and phase mass balances
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = iapws95.Iapws95ParameterBlock(
        default={"phase_presentation": iapws95.PhaseType.G})
    m.fs.heater = Heater(
        default={"property_package": m.fs.properties,
                 "material_balance_type": MaterialBalanceType.componentPhase})
    m.fs.heater.inlet.enth_mol.fix(50000)
    m.fs.heater.inlet.flow_mol.fix(100)
    m.fs.heater.inlet.pressure.fix(101325)
    m.fs.heater.heat_duty[0].fix(100*10000)
    m.fs.heater.initialize()
    assert degrees_of_freedom(m) == 0
    solver.solve(m)
    prop_in = m.fs.heater.control_volume.properties_in[0]
    prop_out = m.fs.heater.control_volume.properties_out[0]
    assert abs(value(prop_in.temperature) - 422.60419933276177) <= 1e-4
    assert abs(value(prop_out.temperature) - 698.1604861702295) <= 1e-4
    assert abs(value(prop_in.phase_frac["Liq"]) - 0) <= 1e-6
    assert abs(value(prop_out.phase_frac["Liq"]) - 0) <= 1e-6
    assert abs(value(prop_in.phase_frac["Vap"]) - 1) <= 1e-6
    assert abs(value(prop_out.phase_frac["Vap"]) - 1) <= 1e-6


@pytest.mark.unit
def test_heater_tpx_g_phase():
    """Test vapor phase only form with T-P-x state vars and phase mass balances
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = iapws95.Iapws95ParameterBlock(
        default={"phase_presentation": iapws95.PhaseType.G,
                 "state_vars": iapws95.StateVars.TPX})
    m.fs.heater = Heater(
        default={"property_package": m.fs.properties,
                 "material_balance_type": MaterialBalanceType.componentPhase})
    m.fs.heater.inlet.flow_mol.fix(100)
    m.fs.heater.inlet.temperature.fix(422.60419933276177)
    m.fs.heater.inlet.pressure.fix(101325)
    m.fs.heater.heat_duty[0].fix(100*10000)
    prop_in = m.fs.heater.control_volume.properties_in[0]
    prop_out = m.fs.heater.control_volume.properties_out[0]
    prop_out.temperature = 550
    prop_out.pressure = 101325
    prop_out.flow_mol = 100
    m.fs.heater.initialize(outlvl=5)
    assert degrees_of_freedom(m) == 0
    solver.solve(m)
    assert abs(value(prop_out.temperature) - 698.1604861702295) <= 1e-4
    assert abs(value(prop_in.phase_frac["Liq"]) - 0) <= 1e-6
    assert abs(value(prop_out.phase_frac["Liq"]) - 0) <= 1e-6
    assert abs(value(prop_in.phase_frac["Vap"]) - 1) <= 1e-6
    assert abs(value(prop_out.phase_frac["Vap"]) - 1) <= 1e-6


@pytest.mark.unit
def test_heater_tpx_lg_total():
    """Test liquid/vapor form with T-P-x state vars and total mass balances. In
    this case you end up with two-phases at the end.
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = iapws95.Iapws95ParameterBlock(
        default={"phase_presentation": iapws95.PhaseType.LG,
                 "state_vars": iapws95.StateVars.TPX})
    m.fs.heater = Heater(
        default={"property_package": m.fs.properties,
                 "material_balance_type": MaterialBalanceType.componentTotal})
    m.fs.heater.inlet.temperature.fix(326.1667075078748)
    m.fs.heater.inlet.vapor_frac.fix(0.0)
    m.fs.heater.inlet.flow_mol.fix(100)
    m.fs.heater.inlet.pressure.fix(101325)
    m.fs.heater.heat_duty[0].fix(100*20000)
    prop_in = m.fs.heater.control_volume.properties_in[0]
    prop_out = m.fs.heater.control_volume.properties_out[0]
    prop_out.temperature = 400
    prop_out.vapor_frac = 0.9
    m.fs.heater.initialize(outlvl=5)
    assert degrees_of_freedom(m) == 0
    solver.solve(m, tee=True)
    assert value(prop_in.temperature) == pytest.approx(326.166707, rel=1e-5)
    assert value(prop_out.temperature) == pytest.approx(373.12429, rel=1e-5)
    assert value(prop_in.phase_frac["Liq"]) == pytest.approx(1, rel=1e-5)
    assert value(prop_out.phase_frac["Liq"]) == pytest.approx(
        0.5953219, rel=1e-5)
    assert value(prop_in.phase_frac["Vap"]) == pytest.approx(0, abs=1e-5)
    assert value(prop_out.phase_frac["Vap"]) == pytest.approx(
        0.4046781, rel=1e-4)


@pytest.mark.unit
def test_heater_tpx_lg_total_2():
    """Test liquid/vapor form with T-P-x state vars and total mass balances. In
    this case you end up with all vapor at the end.
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = iapws95.Iapws95ParameterBlock(
        default={"phase_presentation": iapws95.PhaseType.LG,
                 "state_vars": iapws95.StateVars.TPX})
    m.fs.heater = Heater(
        default={"property_package": m.fs.properties,
                 "material_balance_type": MaterialBalanceType.componentTotal})
    m.fs.heater.inlet.temperature.fix(326.1667075078748)
    m.fs.heater.inlet.vapor_frac.fix(0.0)
    m.fs.heater.inlet.flow_mol.fix(100)
    m.fs.heater.inlet.pressure.fix(101325)
    m.fs.heater.heat_duty[0].fix(100*50000)
    prop_in = m.fs.heater.control_volume.properties_in[0]
    prop_out = m.fs.heater.control_volume.properties_out[0]
    prop_out.temperature = 400
    prop_out.vapor_frac = 0.9
    m.fs.heater.initialize(outlvl=5)
    assert degrees_of_freedom(m) == 0
    solver.solve(m, tee=True)
    assert abs(value(prop_out.temperature) - 534.6889772922356) <= 2e-4
    assert abs(value(prop_in.phase_frac["Liq"]) - 1) <= 1e-2
    assert abs(value(prop_out.phase_frac["Liq"]) - 0) <= 1e-2
    assert abs(value(prop_in.phase_frac["Vap"]) - 0) <= 1e-2
    assert abs(value(prop_out.phase_frac["Vap"]) - 1) <= 1e-2


@pytest.mark.unit
def test_heater_tpx_lg_phase():
    """Test liquid/vapor form with T-P-x state vars and phase mass balances.

    This again has the problem where the property package calculates a vapor
    fraction, but you have a mass balance for each phase.
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = iapws95.Iapws95ParameterBlock(
        default={"phase_presentation": iapws95.PhaseType.LG,
                 "state_vars": iapws95.StateVars.TPX})
    m.fs.heater = Heater(
        default={"property_package": m.fs.properties,
                 "material_balance_type": MaterialBalanceType.componentPhase})
    m.fs.heater.inlet.temperature.fix(326.1667075078748)
    m.fs.heater.inlet.vapor_frac.fix(0.0)
    m.fs.heater.inlet.flow_mol.fix(100)
    m.fs.heater.inlet.pressure.fix(101325)
    m.fs.heater.heat_duty[0].fix(100*50000)

    assert degrees_of_freedom(m) == -1
