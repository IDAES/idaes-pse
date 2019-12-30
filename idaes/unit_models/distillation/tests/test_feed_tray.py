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
Tests for feed tray unit model.

Author: Jaffer Ghouse
"""
import pytest
from pyomo.environ import (ConcreteModel, TerminationCondition,
                           SolverStatus, value)

from idaes.core import (FlowsheetBlock, MaterialBalanceType, EnergyBalanceType,
                        MomentumBalanceType)
from idaes.unit_models.distillation import Tray
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

###############################################################################


def test_build():
    m.fs.tray = Tray(default={"property_package": m.fs.properties,
                              "is_feed_tray": True,
                              "has_heat_transfer": True,
                              "has_pressure_change": True})

    assert len(m.fs.tray.config) == 9

    assert m.fs.tray.config.is_feed_tray

    assert hasattr(m.fs.tray, "feed")
    assert hasattr(m.fs.tray.feed, "flow_mol")
    assert hasattr(m.fs.tray.feed, "mole_frac_comp")
    assert hasattr(m.fs.tray.feed, "temperature")
    assert hasattr(m.fs.tray.feed, "pressure")

    assert hasattr(m.fs.tray, "liq_in")
    assert hasattr(m.fs.tray.liq_in, "flow_mol")
    assert hasattr(m.fs.tray.liq_in, "mole_frac_comp")
    assert hasattr(m.fs.tray.liq_in, "temperature")
    assert hasattr(m.fs.tray.liq_in, "pressure")

    assert hasattr(m.fs.tray, "vap_in")
    assert hasattr(m.fs.tray.vap_in, "flow_mol")
    assert hasattr(m.fs.tray.vap_in, "mole_frac_comp")
    assert hasattr(m.fs.tray.vap_in, "temperature")
    assert hasattr(m.fs.tray.vap_in, "pressure")

    assert hasattr(m.fs.tray, "liq_out")
    assert hasattr(m.fs.tray.liq_out, "flow_mol")
    assert hasattr(m.fs.tray.liq_out, "mole_frac_comp")
    assert hasattr(m.fs.tray.liq_out, "temperature")
    assert hasattr(m.fs.tray.liq_out, "pressure")

    assert hasattr(m.fs.tray, "vap_out")
    assert hasattr(m.fs.tray.vap_out, "flow_mol")
    assert hasattr(m.fs.tray.vap_out, "mole_frac_comp")
    assert hasattr(m.fs.tray.vap_out, "temperature")
    assert hasattr(m.fs.tray.vap_out, "pressure")

    assert m.fs.tray.config.has_heat_transfer
    assert hasattr(m.fs.tray, "heat_duty")

    assert m.fs.tray.config.has_pressure_change
    assert hasattr(m.fs.tray, "deltaP")

    assert not m.fs.tray.config.has_liquid_side_draw
    assert not hasattr(m.fs.tray, "liq_side_sf")
    assert not m.fs.tray.config.has_vapor_side_draw
    assert not hasattr(m.fs.tray, "vap_side_sf")


def test_set_inputs():

    # Check variables and constraints when using FTPz
    assert number_variables(m.fs.tray) == 94
    assert number_total_constraints(m.fs.tray) == 77

    # Set inputs
    m.fs.tray.feed.flow_mol.fix(1)
    m.fs.tray.feed.temperature.fix(369)
    m.fs.tray.feed.pressure.fix(101325)
    m.fs.tray.feed.mole_frac_comp[0, "benzene"].fix(0.5)
    m.fs.tray.feed.mole_frac_comp[0, "toluene"].fix(0.5)

    m.fs.tray.liq_in.flow_mol.fix(1)
    m.fs.tray.liq_in.temperature.fix(369)
    m.fs.tray.liq_in.pressure.fix(101325)
    m.fs.tray.liq_in.mole_frac_comp[0, "benzene"].fix(0.5)
    m.fs.tray.liq_in.mole_frac_comp[0, "toluene"].fix(0.5)

    m.fs.tray.vap_in.flow_mol.fix(1)
    m.fs.tray.vap_in.temperature.fix(372)
    m.fs.tray.vap_in.pressure.fix(101325)
    m.fs.tray.vap_in.mole_frac_comp[0, "benzene"].fix(0.5)
    m.fs.tray.vap_in.mole_frac_comp[0, "toluene"].fix(0.5)

    m.fs.tray.deltaP.fix(0)
    m.fs.tray.heat_duty.fix(0)

    assert degrees_of_freedom(m.fs.tray) == 0


def test_initialize():
    m.fs.tray.initialize(solver=solver)

    assert degrees_of_freedom(m.fs.tray) == 0


def test_solve():
    solve_status = solver.solve(m)
    assert solve_status.solver.termination_condition == \
        TerminationCondition.optimal
    assert solve_status.solver.status == SolverStatus.ok


def test_solution():

    # liq_out port
    assert (pytest.approx(0.92409, abs=1e-3) ==
            value(m.fs.tray.liq_out.flow_mol[0]))
    assert (pytest.approx(0.34840, abs=1e-3) ==
            value(m.fs.tray.liq_out.mole_frac_comp[0, "benzene"]))
    assert (pytest.approx(0.65159, abs=1e-3) ==
            value(m.fs.tray.liq_out.mole_frac_comp[0, "toluene"]))
    assert (pytest.approx(370.056, abs=1e-3) ==
            value(m.fs.tray.liq_out.temperature[0]))
    assert (pytest.approx(101325, abs=1e-3) ==
            value(m.fs.tray.liq_out.pressure[0]))

    # vap_out port
    assert (pytest.approx(2.0759, abs=1e-3) ==
            value(m.fs.tray.vap_out.flow_mol[0]))
    assert (pytest.approx(0.56748, abs=1e-3) ==
            value(m.fs.tray.vap_out.mole_frac_comp[0, "benzene"]))
    assert (pytest.approx(0.43252, abs=1e-3) ==
            value(m.fs.tray.vap_out.mole_frac_comp[0, "toluene"]))
    assert (pytest.approx(370.056, abs=1e-3) ==
            value(m.fs.tray.vap_out.temperature[0]))
    assert (pytest.approx(101325, abs=1e-3) ==
            value(m.fs.tray.vap_out.pressure[0]))
