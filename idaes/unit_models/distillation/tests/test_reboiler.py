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
Tests for Reboiler unit model.

Author: Jaffer Ghouse
"""
import pytest
from pyomo.environ import (ConcreteModel, TerminationCondition,
                           SolverStatus, value)

from idaes.core import FlowsheetBlock, MaterialBalanceType, EnergyBalanceType, \
    MomentumBalanceType
from idaes.unit_models.distillation import Reboiler
from idaes.property_models.activity_coeff_models.BTX_activity_coeff_VLE \
    import BTXParameterBlock
from idaes.core.util.model_statistics import degrees_of_freedom, \
    number_variables, number_total_constraints
from idaes.core.util.testing import get_default_solver

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
# total reboiler with FTPz


def test_build():
    m.fs.R101 = Reboiler(default={"property_package": m.fs.properties})

    assert len(m.fs.R101.config) == 9

    assert m.fs.R101.config.has_boilup_ratio
    assert m.fs.R101.config.material_balance_type == \
        MaterialBalanceType.componentPhase
    assert m.fs.R101.config.energy_balance_type == \
        EnergyBalanceType.enthalpyTotal
    assert m.fs.R101.config.momentum_balance_type == \
        MomentumBalanceType.pressureTotal
    assert m.fs.R101.config.has_pressure_change

    assert hasattr(m.fs.R101, "inlet")
    assert hasattr(m.fs.R101.inlet, "flow_mol")
    assert hasattr(m.fs.R101.inlet, "mole_frac_comp")
    assert hasattr(m.fs.R101.inlet, "temperature")
    assert hasattr(m.fs.R101.inlet, "pressure")

    assert hasattr(m.fs.R101, "bottoms")
    assert hasattr(m.fs.R101.bottoms, "flow_mol")
    assert hasattr(m.fs.R101.bottoms, "mole_frac_comp")
    assert hasattr(m.fs.R101.bottoms, "temperature")
    assert hasattr(m.fs.R101.bottoms, "pressure")

    assert hasattr(m.fs.R101, "vapor_reboil")
    assert hasattr(m.fs.R101.vapor_reboil, "flow_mol")
    assert hasattr(m.fs.R101.vapor_reboil, "mole_frac_comp")
    assert hasattr(m.fs.R101.vapor_reboil, "temperature")
    assert hasattr(m.fs.R101.vapor_reboil, "pressure")

    m.fs.R101_FcTP = Reboiler(default={"property_package": m.fs.properties_2})

    assert len(m.fs.R101_FcTP.config) == 9

    assert m.fs.R101_FcTP.config.has_boilup_ratio
    assert m.fs.R101_FcTP.config.material_balance_type == \
        MaterialBalanceType.componentPhase
    assert m.fs.R101_FcTP.config.energy_balance_type == \
        EnergyBalanceType.enthalpyTotal
    assert m.fs.R101_FcTP.config.momentum_balance_type == \
        MomentumBalanceType.pressureTotal
    assert m.fs.R101_FcTP.config.has_pressure_change

    assert hasattr(m.fs.R101_FcTP, "inlet")
    assert hasattr(m.fs.R101_FcTP.inlet, "flow_mol_comp")
    assert hasattr(m.fs.R101_FcTP.inlet, "temperature")
    assert hasattr(m.fs.R101_FcTP.inlet, "pressure")

    assert hasattr(m.fs.R101_FcTP, "bottoms")
    assert hasattr(m.fs.R101_FcTP.bottoms, "flow_mol_comp")
    assert hasattr(m.fs.R101_FcTP.bottoms, "temperature")
    assert hasattr(m.fs.R101_FcTP.bottoms, "pressure")

    assert hasattr(m.fs.R101_FcTP, "vapor_reboil")
    assert hasattr(m.fs.R101_FcTP.vapor_reboil, "flow_mol_comp")
    assert hasattr(m.fs.R101_FcTP.vapor_reboil, "temperature")
    assert hasattr(m.fs.R101_FcTP.vapor_reboil, "pressure")


def test_set_inputs():

    assert number_variables(m.fs.R101) == 51
    assert number_total_constraints(m.fs.R101) == 44

    # Fix the reboiler variables
    m.fs.R101.boilup_ratio.fix(1)
    m.fs.R101.deltaP.fix(0)

    # Fix the inputs (typically this will be the outlet liquid from the
    # bottom tray)
    m.fs.R101.inlet.flow_mol.fix(1)
    m.fs.R101.inlet.temperature.fix(363)
    m.fs.R101.inlet.pressure.fix(101325)
    m.fs.R101.inlet.mole_frac_comp[0, "benzene"].fix(0.5)
    m.fs.R101.inlet.mole_frac_comp[0, "toluene"].fix(0.5)

    assert degrees_of_freedom(m.fs.R101) == 0

    assert number_variables(m.fs.R101_FcTP) == 53
    assert number_total_constraints(m.fs.R101_FcTP) == 47

    # Fix the reboiler variables
    m.fs.R101_FcTP.boilup_ratio.fix(1)
    m.fs.R101_FcTP.deltaP.fix(0)

    # Fix the inputs (typically this will be the outlet liquid from the
    # bottom tray)
    m.fs.R101_FcTP.inlet.flow_mol_comp[0, "benzene"].fix(0.5)
    m.fs.R101_FcTP.inlet.flow_mol_comp[0, "toluene"].fix(0.5)
    m.fs.R101_FcTP.inlet.temperature.fix(363)
    m.fs.R101_FcTP.inlet.pressure.fix(101325)

    assert degrees_of_freedom(m.fs.R101_FcTP) == 0


def test_solve():

    solver = get_default_solver()
    m.fs.R101.initialize(solver=solver, outlvl=2)

    solve_status = solver.solve(m.fs.R101)

    assert solve_status.solver.termination_condition == \
        TerminationCondition.optimal
    assert solve_status.solver.status == SolverStatus.ok

    m.fs.R101_FcTP.initialize(solver=solver, outlvl=2)

    solve_status = solver.solve(m.fs.R101_FcTP)

    assert solve_status.solver.termination_condition == \
        TerminationCondition.optimal
    assert solve_status.solver.status == SolverStatus.ok


def test_solution():

    # Reboiler when using FTPz

    # Bottoms port
    assert (pytest.approx(0.5, abs=1e-3) ==
            value(m.fs.R101.bottoms.flow_mol[0]))
    assert (pytest.approx(0.3891, abs=1e-3) ==
            value(m.fs.R101.bottoms.mole_frac_comp[0, "benzene"]))
    assert (pytest.approx(0.6109, abs=1e-3) ==
            value(m.fs.R101.bottoms.mole_frac_comp[0, "toluene"]))
    assert (pytest.approx(368.728, abs=1e-3) ==
            value(m.fs.R101.bottoms.temperature[0]))
    assert (pytest.approx(101325, abs=1e-3) ==
            value(m.fs.R101.bottoms.pressure[0]))

    # Vapor reboil port
    assert (pytest.approx(0.5, abs=1e-3) ==
            value(m.fs.R101.vapor_reboil.flow_mol[0]))
    assert (pytest.approx(0.6108, abs=1e-3) ==
            value(m.fs.R101.vapor_reboil.mole_frac_comp[0, "benzene"]))
    assert (pytest.approx(0.3892, abs=1e-3) ==
            value(m.fs.R101.vapor_reboil.mole_frac_comp[0, "toluene"]))
    assert (pytest.approx(368.728, abs=1e-3) ==
            value(m.fs.R101.vapor_reboil.temperature[0]))
    assert (pytest.approx(101325, abs=1e-3) ==
            value(m.fs.R101.bottoms.pressure[0]))

    # Unit level
    assert (pytest.approx(16926.526, abs=1e-3) ==
            value(m.fs.R101.heat_duty[0]))

    # Reboiler when using FcTP

    # Bottoms port
    assert (pytest.approx(0.19455, abs=1e-3) ==
            value(m.fs.R101_FcTP.bottoms.flow_mol_comp[0, "benzene"]))
    assert (pytest.approx(0.30545, abs=1e-3) ==
            value(m.fs.R101_FcTP.bottoms.flow_mol_comp[0, "toluene"]))
    assert (pytest.approx(368.728, abs=1e-3) ==
            value(m.fs.R101_FcTP.bottoms.temperature[0]))
    assert (pytest.approx(101325, abs=1e-3) ==
            value(m.fs.R101_FcTP.bottoms.pressure[0]))

    # Vapor reboil port
    assert (pytest.approx(0.3054, abs=1e-3) ==
            value(m.fs.R101_FcTP.vapor_reboil.flow_mol_comp[0, "benzene"]))
    assert (pytest.approx(0.1946, abs=1e-3) ==
            value(m.fs.R101_FcTP.vapor_reboil.flow_mol_comp[0, "toluene"]))
    assert (pytest.approx(368.728, abs=1e-3) ==
            value(m.fs.R101_FcTP.vapor_reboil.temperature[0]))
    assert (pytest.approx(101325, abs=1e-3) ==
            value(m.fs.R101_FcTP.bottoms.pressure[0]))

    # Unit level
    assert (pytest.approx(16926.370, abs=1e-3) ==
            value(m.fs.R101_FcTP.heat_duty[0]))
