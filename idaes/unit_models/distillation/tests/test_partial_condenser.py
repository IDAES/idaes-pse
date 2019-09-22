##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
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
Tests for partial condenser unit model.

Author: Jaffer Ghouse
"""
import pytest
from pyomo.environ import (ConcreteModel, TerminationCondition,
                           SolverStatus, value)

from idaes.core import (FlowsheetBlock, MaterialBalanceType, EnergyBalanceType,
                        MomentumBalanceType, useDefault)
from idaes.unit_models.distillation import Condenser
from idaes.unit_models.distillation.condenser import CondenserType
from idaes.property_models.activity_coeff_models.BTX_activity_coeff_VLE \
    import BTXParameterBlock
from idaes.core.util.model_statistics import degrees_of_freedom, \
    number_variables, number_total_constraints
from idaes.core.util.testing import get_default_solver


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_default_solver()

m = ConcreteModel()
m.fs = FlowsheetBlock(default={"dynamic": False})
m.fs.properties = BTXParameterBlock(default={"valid_phase":
                                             ('Liq', 'Vap'),
                                             "activity_coeff_model":
                                             "Ideal"})
m.fs.properties_2 = BTXParameterBlock(default={"valid_phase":
                                               ('Liq', 'Vap'),
                                               "activity_coeff_model":
                                               "Ideal",
                                               "state_vars":
                                               "FcTP"})

###############################################################################


def test_build_partial_condenser():
    m.fs.C101_partial = Condenser(
        default={"property_package": m.fs.properties,
                 "condenser_type": CondenserType.partialCondenser})

    assert len(m.fs.C101_partial.config) == 9
    assert m.fs.C101_partial.config.condenser_type ==\
        CondenserType.partialCondenser
    assert m.fs.C101_partial.config.material_balance_type == \
        MaterialBalanceType.componentPhase
    assert m.fs.C101_partial.config.energy_balance_type == \
        EnergyBalanceType.enthalpyTotal
    assert m.fs.C101_partial.config.momentum_balance_type == \
        MomentumBalanceType.pressureTotal
    assert m.fs.C101_partial.config.has_pressure_change

    assert hasattr(m.fs.C101_partial. "reflux_ratio")
    assert not hasattr(m.fs.C101_partial. "eq_total_cond_spec")

    assert hasattr(m.fs.C101_partial. "inlet")

    assert hasattr(m.fs.C101_partial.inlet, "flow_mol")
    assert hasattr(m.fs.C101_partial.inlet, "mole_frac")
    assert hasattr(m.fs.C101_partial.inlet, "temperature")
    assert hasattr(m.fs.C101_partial.inlet, "pressure")

    assert hasattr(m.fs.C101_partial. "reflux")

    assert hasattr(m.fs.C101_partial.reflux, "flow_mol")
    assert hasattr(m.fs.C101_partial.reflux, "mole_frac")
    assert hasattr(m.fs.C101_partial.reflux, "temperature")
    assert hasattr(m.fs.C101_partial.reflux, "pressure")

    assert hasattr(m.fs.C101_partial. "distillate")
    assert hasattr(m.fs.C101_partial.distillate, "flow_mol")
    assert hasattr(m.fs.C101_partial.distillate, "mole_frac")
    assert hasattr(m.fs.C101_partial.distillate, "temperature")
    assert hasattr(m.fs.C101_partial.distillate, "pressure")

    assert hasattr(m.fs.C101_partial. "vapor_outlet")
    assert hasattr(m.fs.C101_partial.vapor_outlet, "flow_mol")
    assert hasattr(m.fs.C101_partial.vapor_outlet, "mole_frac")
    assert hasattr(m.fs.C101_partial.vapor_outlet, "temperature")
    assert hasattr(m.fs.C101_partial.vapor_outlet, "pressure")

    m.fs.C101_total_FcTP = Condenser(
        default={"property_package": m.fs.properties_2,
                 "condenser_type": CondenserType.totalCondenser})

    assert len(m.fs.C101_total_FcTP.config) == 9
    assert m.fs.C101_total_FcTP.config.condenser_type == \
        CondenserType.partialCondenser
    assert m.fs.C101_partial_FcTP.config.material_balance_type == \
        MaterialBalanceType.componentPhase
    assert m.fs.C101_partial_FcTP.config.energy_balance_type == \
        EnergyBalanceType.enthalpyTotal
    assert m.fs.C101_partial_FcTP.config.momentum_balance_type == \
        MomentumBalanceType.pressureTotal
    assert m.fs.C101_partial_FcTP.config.has_pressure_change

    assert hasattr(m.fs.C101_partial_FcTP, "reflux_ratio")
    assert hasattr(m.fs.C101_partial_FcTP, "eq_total_cond_spec")

    assert hasattr(m.fs.C101_partial_FcTP, "inlet")

    assert hasattr(m.fs.C101_partial_FcTP.inlet, "flow_mol_comp")
    assert hasattr(m.fs.C101_partial_FcTP.inlet, "temperature")
    assert hasattr(m.fs.C101_partial_FcTP.inlet, "pressure")

    assert hasattr(m.fs.C101_partial_FcTP, "reflux")

    assert hasattr(m.fs.C101_partial_FcTP.reflux, "flow_mol_comp")
    assert hasattr(m.fs.C101_partial_FcTP.reflux, "temperature")
    assert hasattr(m.fs.C101_partial_FcTP.reflux, "pressure")

    assert hasattr(m.fs.C101_partial_FcTP, "distillate")
    assert hasattr(m.fs.C101_partial_FcTP.distillate, "flow_mol_comp")
    assert hasattr(m.fs.C101_partial_FcTP.distillate, "temperature")
    assert hasattr(m.fs.C101_partial_FcTP.distillate, "pressure")

def test_set_inputs():

    # Check variables and constraints when using FTPz
    assert number_variables(m.fs.C101_total) == 51
    assert number_total_constraints(m.fs.C101_total) == 44

    # Fix the partial condenser variables
    m.fs.C101_partial.reflux_ratio.fix(1)
    m.fs.C101_partial.deltaP.fix(0)

    # Fix the inputs (typically this will be the outlet vapor from the top tray)
    m.fs.C101_partial.inlet.flow_mol.fix(1)
    m.fs.C101_partial.inlet.temperature.fix(375)
    m.fs.C101_partial.inlet.pressure.fix(101325)
    m.fs.C101_partial.inlet.mole_frac[0, "benzene"].fix(0.5)
    m.fs.C101_partial.inlet.mole_frac[0, "toluene"].fix(0.5)

    assert degrees_of_freedom(m.fs.C101_partial) == 0

    # Check variables and constraints when using FcTP
    assert number_variables(m.fs.C101_partial_FcTP) == 53
    assert number_total_constraints(m.fs.C101_partial_FcTP) == 47

    # Fix the partial condenser variables
    m.fs.C101_partial_FcTP.reflux_ratio.fix(1)
    m.fs.C101_partial_FcTP.deltaP.fix(0)

    # Fix the inputs (typically this will be the outlet vapor from the top tray)
    m.fs.C101_partial_FcTP.inlet.flow_mol_comp[0, "benzene"].fix(0.5)
    m.fs.C101_partial_FcTP.inlet.flow_mol_comp[0, "toluene"].fix(0.5)
    m.fs.C101_partial_FcTP.inlet.temperature.fix(375)
    m.fs.C101_partial_FcTP.inlet.pressure.fix(101325)

    assert degrees_of_freedom(m.fs.C101_partial_FcTP) == 0

def test_solve():
    # Test the solve to optimality for FTPz state vars
    m.fs.C101_partial.initialize(solver=solver, outlvl=1)

    solve_status = solver.solve(m.fs.C101_partial)

    assert solve_status.solver.termination_condition == \
        TerminationCondition.optimal
    assert solve_status.solver.status == SolverStatus.ok

    # Test the solve to optimality for FcTP state vars
    m.fs.C101_partial_FcTP.initialize(solver=solver, outlvl=1)

    solve_status = solver.solve(m.fs.C101_partial_FcTP)

    assert solve_status.solver.termination_condition == \
        TerminationCondition.optimal
    assert solve_status.solver.status == SolverStatus.ok

def test_solution():

    # Total condenser when using FTPz
    # Reflux port
    assert (pytest.approx(0.4999, abs=1e-3) ==
            value(m.fs.C101_partial.reflux.flow_mol[0]))
    assert (pytest.approx(0.5, abs=1e-3) ==
            value(m.fs.C101_partial.reflux.mole_frac[0, "benzene"]))
    assert (pytest.approx(0.5, abs=1e-3) ==
            value(m.fs.C101_partial.reflux.mole_frac[0, "toluene"]))
    assert (pytest.approx(365.347, abs=1e-3) ==
            value(m.fs.C101_partial.reflux.temperature[0]))
    assert (pytest.approx(101325, abs=1e-3) ==
            value(m.fs.C101_partial.reflux.pressure[0]))

    # Distillate port
    assert (pytest.approx(0.4999, abs=1e-3) ==
            value(m.fs.C101_partial.distillate.flow_mol[0]))
    assert (pytest.approx(0.5, abs=1e-3) ==
            value(m.fs.C101_partial.distillate.mole_frac[0, "benzene"]))
    assert (pytest.approx(0.5, abs=1e-3) ==
            value(m.fs.C101_partial.distillate.mole_frac[0, "toluene"]))
    assert (pytest.approx(365.347, abs=1e-3) ==
            value(m.fs.C101_partial.distillate.temperature[0]))
    assert (pytest.approx(101325, abs=1e-3) ==
            value(m.fs.C101_partial.distillate.pressure[0]))

    # Unit level
    assert (pytest.approx(-33727.350, abs=1e-3) ==
            value(m.fs.C101_partial.heat_duty[0]))

    # Total condenser when using FcTP
    # Reflux port
    assert (pytest.approx(0.25, abs=1e-3) ==
            value(m.fs.C101_partial_FcTP.reflux.flow_mol_comp[0, "benzene"]))
    assert (pytest.approx(0.25, abs=1e-3) ==
            value(m.fs.C101_partial_FcTP.reflux.flow_mol_comp[0, "toluene"]))
    assert (pytest.approx(365.347, abs=1e-3) ==
            value(m.fs.C101_partial_FcTP.reflux.temperature[0]))
    assert (pytest.approx(101325, abs=1e-3) ==
            value(m.fs.C101_partial_FcTP.reflux.pressure[0]))

    # Distillate port
    assert (pytest.approx(0.25, abs=1e-3) ==
            value(m.fs.C101_partial_FcTP.distillate.flow_mol_comp[0, "benzene"]))
    assert (pytest.approx(0.25, abs=1e-3) ==
            value(m.fs.C101_partial_FcTP.distillate.flow_mol_comp[0, "toluene"]))
    assert (pytest.approx(365.347, abs=1e-3) ==
            value(m.fs.C101_partial_FcTP.distillate.temperature[0]))
    assert (pytest.approx(101325, abs=1e-3) ==
            value(m.fs.C101_partial_FcTP.distillate.pressure[0]))

    # Unit level
    assert (pytest.approx(-33727.333, abs=1e-3) ==
            value(m.fs.C101_partial_FcTP.heat_duty[0]))
 print("-----------------------------------------------------------------------")

 ###############################################################################
 # Partial Condenser with FTPz
 m.fs.C101_partial = Condenser(default={"property_package": m.fs.properties,
                               "condenser_type": CondenserType.partialCondenser})

 # Fix the partial condenser variables
 m.fs.C101_partial.reflux_ratio.fix(0.5)
 m.fs.C101_partial.deltaP.fix(0)
 m.fs.C101_partial.vapor_outlet.temperature.fix(367.5)

 # Fix the inputs (typically this will be the outlet vapor from the top tray)
 m.fs.C101_partial.inlet.flow_mol.fix(1)
 m.fs.C101_partial.inlet.temperature.fix(368)
 m.fs.C101_partial.inlet.pressure.fix(101325)
 m.fs.C101_partial.inlet.mole_frac[0, "benzene"].fix(0.5)
 m.fs.C101_partial.inlet.mole_frac[0, "toluene"].fix(0.5)

 print("Partial Condenser - FTPz")
 print("The degrees of freedom is ", degrees_of_freedom(m))

 m.fs.C101_partial.initialize(solver=solver, outlvl=2)

 m.fs.C101_partial.reflux.display()
 m.fs.C101_partial.distillate.display()
 m.fs.C101_partial.vapor_outlet.display()
 m.fs.C101_partial.heat_duty.display()

 ###############################################################################
 # total condenser with FcTP
 m.fs.C101_total_FcTP = Condenser(default={"property_package": m.fs.properties_2,
                                           "condenser_type":
                                           CondenserType.totalCondenser})

 # Fix the partial condenser variables
 m.fs.C101_total_FcTP.reflux_ratio.fix(1)
 m.fs.C101_total_FcTP.deltaP.fix(0)

 # Fix the inputs (typically this will be the outlet vapor from the top tray)
 m.fs.C101_total_FcTP.inlet.flow_mol_comp[0, "benzene"].fix(0.5)
 m.fs.C101_total_FcTP.inlet.flow_mol_comp[0, "toluene"].fix(0.5)
 m.fs.C101_total_FcTP.inlet.temperature.fix(368)
 m.fs.C101_total_FcTP.inlet.pressure.fix(101325)


 print("-----------------------------------------------------------------------")
 print("Total Condenser - FcTP")
 print("The degrees of freedom is ", degrees_of_freedom(m.fs.C101_total_FcTP))


 solver = get_default_solver()
 m.fs.C101_total_FcTP.initialize(solver=solver, outlvl=2)

 m.fs.C101_total_FcTP.reflux.display()
 m.fs.C101_total_FcTP.distillate.display()
 m.fs.C101_total_FcTP.heat_duty.display()

 print("-----------------------------------------------------------------------")

 # Partial Condenser with FcTP
 m.fs.C101_partial_FcTP = Condenser(default={"property_package": m.fs.properties_2,
                                             "condenser_type":
                                             CondenserType.partialCondenser})

 # Fix the partial condenser variables
 m.fs.C101_partial_FcTP.reflux_ratio.fix(0.5)
 m.fs.C101_partial_FcTP.deltaP.fix(0)
 m.fs.C101_partial_FcTP.vapor_outlet.temperature.fix(367.5)

 # Fix the inputs (typically this will be the outlet vapor from the top tray)
 m.fs.C101_partial_FcTP.inlet.flow_mol_comp[0, "benzene"].fix(0.5)
 m.fs.C101_partial_FcTP.inlet.flow_mol_comp[0, "toluene"].fix(0.5)
 m.fs.C101_partial_FcTP.inlet.temperature.fix(368)
 m.fs.C101_partial_FcTP.inlet.pressure.fix(101325)


 print("Partial Condenser - FcTP")
 print("The degrees of freedom is ", degrees_of_freedom(m.fs.C101_partial_FcTP))

 m.fs.C101_partial_FcTP.initialize(solver=solver, outlvl=2)

 m.fs.C101_partial_FcTP.reflux.display()
 m.fs.C101_partial_FcTP.distillate.display()
 m.fs.C101_partial_FcTP.vapor_outlet.display()
 m.fs.C101_partial_FcTP.heat_duty.display()
