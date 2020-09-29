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
Test for costing correlations of flash tank, fired heater,
and vessels (including platform and ladders, and distillation columns)

Created on Sept 25, 2020 by M. Zamarripa
"""
import pytest
# Import Pyomo libraries
import pyomo.environ as pyo
# Import IDAES core
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
# Import Unit Model Modules
from idaes.generic_models.properties import iapws95
from idaes.core.util.testing import get_default_solver
import idaes.core.util.unit_costing as cs
from idaes.generic_models.unit_models.flash import Flash, EnergySplittingType
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from idaes.power_generation.properties import FlueGasParameterBlock
from idaes.generic_models.unit_models.pressure_changer import (
    PressureChanger,
    ThermodynamicAssumption,
)
# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_default_solver()

# -----------------------------------------------------------------------------


@pytest.mark.unit
def test_costing_FH_build():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.get_costing()
    m.fs.costing.CE_index = 550  # for testing only
    m.fs.unit = pyo.Block()
    m.fs.unit.heat_duty = pyo.Var(initialize=1e6)
    m.fs.unit.pressure = pyo.Var(initialize=1e5)
    m.fs.unit.costing = pyo.Block()
    m.fs.unit.heat_duty.fix(18390000)  # Btu/hr
    m.fs.unit.pressure.fix(700)  # psig

    cs.fired_heater_costing(m.fs.unit.costing,
                            fired_type='fuel',
                            Mat_factor='stain_steel',
                            ref_parameter_pressure=m.fs.unit.pressure,
                            ref_parameter_heat_duty=m.fs.unit.heat_duty)

    assert degrees_of_freedom(m) == 0
    # Check unit config arguments
    assert isinstance(m.fs.unit.costing.purchase_cost, pyo.Var)
    assert isinstance(m.fs.unit.costing.base_cost, pyo.Var)


@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_costing_FH_solve():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.get_costing()
    m.fs.costing.CE_index = 550  # for testing only
    m.fs.unit = pyo.Block()
    m.fs.unit.heat_duty = pyo.Var(initialize=1e6)
    m.fs.unit.pressure = pyo.Var(initialize=1e5)
    m.fs.unit.costing = pyo.Block()
    m.fs.unit.heat_duty.fix(18390000)  # Btu/hr
    m.fs.unit.pressure.fix(700)  # psig

    cs.fired_heater_costing(m.fs.unit.costing,
                            fired_type='fuel',
                            Mat_factor='stain_steel',
                            ref_parameter_pressure=m.fs.unit.pressure,
                            ref_parameter_heat_duty=m.fs.unit.heat_duty)

    assert degrees_of_freedom(m) == 0
    # Check unit config arguments
    assert isinstance(m.fs.unit.costing.purchase_cost, pyo.Var)
    assert isinstance(m.fs.unit.costing.base_cost, pyo.Var)

    results = solver.solve(m, tee=False)
    # Check for optimal solution
    assert results.solver.termination_condition == \
        pyo.TerminationCondition.optimal
    assert results.solver.status == pyo.SolverStatus.ok
    assert (pytest.approx(pyo.value(m.fs.unit.costing.purchase_cost),
                          abs=1e-2) == 962795.521)  # Example 22.1 Ref Book


@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_costing_flash():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = iapws95.Iapws95ParameterBlock(default={
        "phase_presentation": iapws95.PhaseType.LG})
    m.fs.unit = Flash(default={"property_package": m.fs.properties,
                               "ideal_separation": False,
                               "energy_split_basis":
                                   EnergySplittingType.enthalpy_split})
    m.fs.unit.inlet.flow_mol.fix(10000)
    m.fs.unit.inlet.enth_mol.fix(24000)
    m.fs.unit.inlet.pressure.fix(101325)
    m.fs.unit.heat_duty.fix(0)
    m.fs.unit.deltaP.fix(0)
    # calculate flash diameter based on heuristics for process vessels
    # Analysis, Synthesis, and Design of Chemical Processes
    # Table 11.6 Heuristics for Process Vessels
    # as a function of density, mass flowrate, and gas velocity
    # dens mass phase = kg/m3, flow_mass in kg/s
    # gas velocity, u = k(rhoL/rhoVap - 1)^0.5 in m/s
    m.fs.get_costing()  # create CE_index in costing global parameters
    # create variables to calculate vessel cost
    m.fs.unit.diameter = pyo.Var(initialize=10,
                                 domain=pyo.NonNegativeReals,
                                 doc='unit diameter in m')
    m.fs.unit.lenght = pyo.Var(initialize=10,
                               domain=pyo.NonNegativeReals,
                               doc='unit lenght in m')

    m.fs.unit.diameter.fix(0.155238)  # in ft
    m.fs.unit.lenght.fix(0.388095)   # in ft

    m.fs.unit.costing = pyo.Block()
    cs.vessel_costing(m.fs.unit.costing,
                      alignment='vertical',
                      weight_limit='option1',
                      L_D_range='option1',
                      PL=True,
                      plates=False)

    assert degrees_of_freedom(m) == 0
    # Check unit config arguments
    assert isinstance(m.fs.unit.costing.purchase_cost, pyo.Var)
    assert isinstance(m.fs.unit.costing.base_cost, pyo.Var)

    results = solver.solve(m, tee=False)
    # Check for optimal solution
    assert results.solver.termination_condition == \
        pyo.TerminationCondition.optimal
    assert results.solver.status == pyo.SolverStatus.ok
    assert isinstance(m.fs.unit.costing.purchase_cost, pyo.Var)
    assert isinstance(m.fs.unit.costing.base_cost, pyo.Var)
    assert isinstance(m.fs.unit.costing.base_cost_platf_ladders, pyo.Var)
    assert isinstance(m.fs.unit.costing.purchase_cost_trays, pyo.Var)
    assert isinstance(m.fs.unit.costing.base_cost, pyo.Var)


@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_costing_distillation_solve():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.get_costing()
    m.fs.costing.CE_index = 550
    # create a unit model and variables
    m.fs.unit = pyo.Block()
    m.fs.unit.heat_duty = pyo.Var(initialize=1e6)
    m.fs.unit.pressure = pyo.Var(initialize=1e5)
    m.fs.unit.diameter = pyo.Var(initialize=10,
                                 domain=pyo.NonNegativeReals,
                                 doc='unit diameter in m')
    m.fs.unit.lenght = pyo.Var(initialize=10,
                               domain=pyo.NonNegativeReals,
                               doc='unit lenght in m')
    # create costing block
    m.fs.unit.costing = pyo.Block()

    cs.vessel_costing(m.fs.unit.costing,
                      alignment='vertical',
                      weight_limit='option2',
                      L_D_range='option2',
                      PL=True,
                      plates=True,
                      number_tray=100,
                      ref_parameter_diameter=m.fs.unit.diameter,
                      ref_parameter_lenght=m.fs.unit.lenght)
    # pressure design and shell thickness from Example 22.13 Product and
    # Process Design Principless
    m.fs.unit.heat_duty.fix(18390000)  # Btu/hr
    m.fs.unit.pressure.fix(123)  # psig
    # pressure design minimum thickness tp = 0.582 in
    # vessel is vertical + quite tall the tower is subject to wind load,
    # and earthquake. Assume wall thickness of 1.25 in.
    # The additional wall thickness at the bottom of the tower is 0.889 in
    # average thickness is 1.027, plus corrosion allowance of 1/8
    # 1.152 in, therefore steel plate thickness is 1.250 (ts)
    m.fs.unit.costing.shell_thickness = 1.250  # inches
    m.fs.unit.diameter.fix(10)  # ft
    m.fs.unit.lenght.fix(212)   # ft
    m.fs.unit.costing.number_trays = 100

    assert degrees_of_freedom(m) == 0
    # Check unit config arguments
    assert isinstance(m.fs.unit.costing.purchase_cost, pyo.Var)
    assert isinstance(m.fs.unit.costing.base_cost, pyo.Var)

    results = solver.solve(m, tee=False)
    # Check for optimal solution
    assert results.solver.termination_condition == \
        pyo.TerminationCondition.optimal
    assert results.solver.status == pyo.SolverStatus.ok
    assert (pytest.approx(pyo.value(m.fs.unit.costing.base_cost),
                          abs=1e-2) == 636959.6929)  # Example 22.13 Ref Book
    assert (pytest.approx(pyo.value(m.fs.unit.costing.base_cost_platf_ladders),
                          abs=1e-2) == 97542.9005)  # Example 22.13 Ref Book
    assert (pytest.approx(pyo.value(m.fs.unit.costing.purchase_cost_trays),
                          abs=1e-2) == 293006.086)  # Example 22.13 Ref Book
    assert (pytest.approx(pyo.value(m.fs.unit.costing.purchase_cost),
                          abs=1e-2) == 1100958.9396)  # Example 22.13 Ref Book


def test_blower_build_and_solve():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.properties = iapws95.Iapws95ParameterBlock()
    # Add property packages to flowsheet library
    m.fs.prop_fluegas = FlueGasParameterBlock()
    m.fs.unit = PressureChanger(default={
            "property_package": m.fs.prop_fluegas,
            "thermodynamic_assumption": ThermodynamicAssumption.isentropic,
            "compressor": True})

    # # FLUE GAS Inlet from Primary Superheater
    FGrate = 8516.27996  # mol/s
    # # Use FG molar composition to set component flow rates (baseline report)
    m.fs.unit.control_volume.properties_in[0].flow_mol_comp["H2O"].\
        fix(FGrate*8.69/100)
    m.fs.unit.control_volume.properties_in[0].flow_mol_comp["CO2"].\
        fix(FGrate*14.49/100)
    m.fs.unit.control_volume.properties_in[0].flow_mol_comp["N2"].\
        fix(FGrate*74.34/100)
    m.fs.unit.control_volume.properties_in[0].flow_mol_comp["O2"].\
        fix(FGrate*2.47/100)
    m.fs.unit.control_volume.properties_in[0].flow_mol_comp["NO"].\
        fix(FGrate*0.0006)
    m.fs.unit.control_volume.properties_in[0].flow_mol_comp["SO2"]\
        .fix(FGrate*0.002)
    m.fs.unit.control_volume.properties_in[0].temperature.fix(200.335)
    m.fs.unit.control_volume.properties_in[0].pressure.fix(99973.98)

    m.fs.unit.deltaP.fix(144790-99973.98)
    m.fs.unit.efficiency_isentropic.fix(0.9)
    m.fs.unit.initialize()
    m.fs.unit.get_costing(mover_type='fan')

    calculate_variable_from_constraint(
        m.fs.unit.costing.purchase_cost,
        m.fs.unit.costing.cp_cost_eq)
    assert degrees_of_freedom(m) == 0
    # Check unit config arguments
    assert isinstance(m.fs.unit.costing.purchase_cost, pyo.Var)
    assert isinstance(m.fs.unit.costing.base_cost, pyo.Var)
    results = solver.solve(m, tee=True)
    assert results.solver.termination_condition == \
        pyo.TerminationCondition.optimal
    assert results.solver.status == pyo.SolverStatus.ok
    assert (pytest.approx(pyo.value(m.fs.unit.costing.base_cost),
                          abs=1e-2) == 56026.447)
    assert (pytest.approx(pyo.value(m.fs.unit.costing.purchase_cost),
                          abs=1e-2) == 272595.280)


def test_compressor_fan():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.properties = iapws95.Iapws95ParameterBlock()
    # Add property packages to flowsheet library
    m.fs.prop_fluegas = FlueGasParameterBlock()
    m.fs.unit = PressureChanger(default={
            "property_package": m.fs.prop_fluegas,
            "thermodynamic_assumption": ThermodynamicAssumption.isentropic,
            "compressor": True})

    # # FLUE GAS Inlet from Primary Superheater
    FGrate = 425.813998  # mol/s
    # # Use FG molar composition to set component flow rates (baseline report)
    m.fs.unit.control_volume.properties_in[0].flow_mol_comp["H2O"].\
        fix(FGrate*8.69/100)
    m.fs.unit.control_volume.properties_in[0].flow_mol_comp["CO2"].\
        fix(FGrate*14.49/100)
    m.fs.unit.control_volume.properties_in[0].flow_mol_comp["N2"].\
        fix(FGrate*74.34/100)
    m.fs.unit.control_volume.properties_in[0].flow_mol_comp["O2"].\
        fix(FGrate*2.47/100)
    m.fs.unit.control_volume.properties_in[0].flow_mol_comp["NO"].\
        fix(FGrate*0.0006)
    m.fs.unit.control_volume.properties_in[0].flow_mol_comp["SO2"]\
        .fix(FGrate*0.002)
    m.fs.unit.control_volume.properties_in[0].temperature.fix(200.335)
    m.fs.unit.control_volume.properties_in[0].pressure.fix(98658.6)

    m.fs.unit.deltaP.fix(101325-98658.6)
    m.fs.unit.efficiency_isentropic.fix(0.9)
    m.fs.unit.initialize()
    m.fs.unit.get_costing(mover_type='fan')

    calculate_variable_from_constraint(
            m.fs.unit.costing.purchase_cost,
            m.fs.unit.costing.cp_cost_eq)
    assert degrees_of_freedom(m) == 0
    # Check unit config arguments
    assert isinstance(m.fs.unit.costing.purchase_cost, pyo.Var)
    assert isinstance(m.fs.unit.costing.base_cost, pyo.Var)
    results = solver.solve(m, tee=True)
    assert results.solver.termination_condition == \
        pyo.TerminationCondition.optimal
    assert results.solver.status == pyo.SolverStatus.ok
    assert (pytest.approx(pyo.value(m.fs.unit.costing.base_cost),
                          abs=1e-2) == 4543.6428)
    assert (pytest.approx(pyo.value(m.fs.unit.costing.purchase_cost),
                          abs=1e-2) == 22106.9807)
