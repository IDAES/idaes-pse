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
Tests for IAPWS using a mixer/separator blocks.  This tests whether the
different options work with mixers and separators, since they are not based
no control volumes.

Author: John Eslick
"""
import pytest

from pyomo.environ import ConcreteModel, SolverFactory, value

from idaes.core import FlowsheetBlock
from idaes.unit_models import (Separator, Mixer, SplittingType,
                               EnergySplittingType, MomentumMixingType)
from idaes.property_models import iapws95
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core import MaterialBalanceType
prop_available = iapws95.iapws95_available()

# -----------------------------------------------------------------------------
# See if ipopt is available and set up solver
if SolverFactory('ipopt').available():
    solver = SolverFactory('ipopt')
    solver.options = {'tol': 1e-6}
else:
    solver = None

@pytest.mark.skipif(not prop_available, reason="IAPWS not available")
@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_sep_ph_mixed_default():
    """Test mixed phase form with P-H state vars and phase mass balances
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = iapws95.Iapws95ParameterBlock(default={})
    m.fs.sep = Separator(default={
        "property_package":m.fs.properties,
        "split_basis":SplittingType.totalFlow,
        "ideal_separation":False,
        "num_outlets":2,
        "energy_split_basis":EnergySplittingType.equal_molar_enthalpy})

    m.fs.sep.inlet.enth_mol.fix(24000)
    m.fs.sep.inlet.flow_mol.fix(100)
    m.fs.sep.inlet.pressure.fix(101325)
    m.fs.sep.split_fraction[0,"outlet_2"].fix(0.10)
    m.fs.sep.initialize()
    assert degrees_of_freedom(m) == 0
    solver.solve(m)
    prop_in = m.fs.sep.mixed_state[0]
    prop_out = m.fs.sep.outlet_1_state[0]
    prop_out2 =  m.fs.sep.outlet_2_state[0]
    assert abs(value(prop_out.temperature) - 373.12429584768876) <= 1e-4
    assert abs(value(prop_out.phase_frac["Liq"]) - 0.5953218682380845) <= 1e-6
    assert abs(value(prop_out.phase_frac["Vap"]) - 0.40467813176191547) <= 1e-6
    assert abs(value(prop_out.flow_mol - 90.0) <= 1e-6)
    assert abs(value(prop_out2.flow_mol - 10.0) <= 1e-6)

@pytest.mark.skipif(not prop_available, reason="IAPWS not available")
@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_sep_ph_lg_default():
    """Test mixed phase form with P-H state vars and phase mass balances
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = iapws95.Iapws95ParameterBlock(default={
        "phase_presentation":iapws95.PhaseType.LG})
    m.fs.sep = Separator(default={
        "property_package":m.fs.properties,
        "split_basis":SplittingType.totalFlow,
        "ideal_separation":False,
        "num_outlets":2,
        "energy_split_basis":EnergySplittingType.equal_molar_enthalpy})

    m.fs.sep.inlet.enth_mol.fix(24000)
    m.fs.sep.inlet.flow_mol.fix(100)
    m.fs.sep.inlet.pressure.fix(101325)
    m.fs.sep.split_fraction[0,"outlet_2"].fix(0.10)
    m.fs.sep.initialize()
    # by default the separator writes phase mass balances so -2 degress of
    # freedom is right.  for IAPWS need total mass balance.
    assert degrees_of_freedom(m) == -2

@pytest.mark.skipif(not prop_available, reason="IAPWS not available")
@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_sep_ph_lg_total_mb():
    """Test mixed phase form with P-H state vars and phase mass balances
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = iapws95.Iapws95ParameterBlock(default={
        "phase_presentation":iapws95.PhaseType.LG})
    m.fs.sep = Separator(default={
        "property_package":m.fs.properties,
        "split_basis":SplittingType.totalFlow,
        "ideal_separation":False,
        "num_outlets":2,
        "material_balance_type":MaterialBalanceType.total,
        "energy_split_basis":EnergySplittingType.equal_molar_enthalpy})

    m.fs.sep.inlet.enth_mol.fix(24000)
    m.fs.sep.inlet.flow_mol.fix(100)
    m.fs.sep.inlet.pressure.fix(101325)
    m.fs.sep.split_fraction[0,"outlet_2"].fix(0.10)
    m.fs.sep.initialize()
    assert degrees_of_freedom(m) == 0
    solver.solve(m)
    prop_in = m.fs.sep.mixed_state[0]
    prop_out = m.fs.sep.outlet_1_state[0]
    prop_out2 =  m.fs.sep.outlet_2_state[0]
    assert abs(value(prop_out.temperature) - 373.12429584768876) <= 1e-4
    assert abs(value(prop_out.phase_frac["Liq"]) - 0.5953218682380845) <= 1e-6
    assert abs(value(prop_out.phase_frac["Vap"]) - 0.40467813176191547) <= 1e-6
    assert abs(value(prop_out.flow_mol - 90.0) <= 1e-6)
    assert abs(value(prop_out2.flow_mol - 10.0) <= 1e-6)
