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
Tests for ControlVolume0D scaling.

Author: John Eslick
"""
import pytest
import pyomo.environ as pyo
from idaes.core import (
    ControlVolume0DBlock,
    ControlVolumeBlockData,
    FlowDirection,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    FlowsheetBlock,
)
from idaes.core.util.testing import (
    PhysicalParameterTestBlock,
    ReactionParameterTestBlock,
)
from idaes.core.util import scaling as iscale
import idaes.logger as idaeslog


# -----------------------------------------------------------------------------
# Basic tests
@pytest.mark.unit
def test_base_build():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.cv = ControlVolume0DBlock(default={
            "property_package": m.fs.pp
        }
    )
    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_material_balances(
        balance_type=MaterialBalanceType.componentTotal,
        has_phase_equilibrium=False)
    m.fs.cv.add_energy_balances(
        balance_type=EnergyBalanceType.enthalpyTotal,
        has_heat_transfer=True,
        has_work_transfer=True)
    # add momentum balance
    m.fs.cv.add_momentum_balances(
        balance_type=MomentumBalanceType.pressureTotal,
        has_pressure_change=True)

    # The scaling factors used for this test were selected to be easy values to
    # test, they do not represent typical scaling factors.
    iscale.set_scaling_factor(m.fs.cv.heat, 11)
    iscale.set_scaling_factor(m.fs.cv.work, 12)
    iscale.calculate_scaling_factors(m)
    # Make sure the heat and work scaling factors are set and not overwitten
    # by the defaults in calculate_scaling_factors
    assert iscale.get_scaling_factor(m.fs.cv.heat) == 11
    assert iscale.get_scaling_factor(m.fs.cv.work) == 12

    # Didn't specify a deltaP scaling factor, so by default pressure in scaling
    # factor * 10 is used.
    for v in m.fs.cv.deltaP.values(): #deltaP is time indexed
        assert iscale.get_scaling_factor(v) == 1040

    # check scaling on mass, energy, and pressure balances.
    for c in m.fs.cv.material_balances.values():
        # this uses the minmum material flow term scale
        assert iscale.get_constraint_transform_applied_scaling_factor(c) == 112
    for c in m.fs.cv.enthalpy_balances.values():
        # this uses the minmum enthalpy flow term scale
        assert iscale.get_constraint_transform_applied_scaling_factor(c) == 110
    for c in m.fs.cv.pressure_balance.values():
        # This uses the inlet pressure scale
        assert iscale.get_constraint_transform_applied_scaling_factor(c) == 104
