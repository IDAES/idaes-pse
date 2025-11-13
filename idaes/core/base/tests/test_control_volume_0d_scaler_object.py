#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Tests for ControlVolume0D scaling.

Author: Douglas Allan
"""
import pytest
import pyomo.environ as pyo
from pyomo.common.collections import ComponentMap
from idaes.core import (
    ControlVolume0DBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    FlowsheetBlock,
)
from idaes.core.scaling import get_scaling_factor, set_scaling_factor
from idaes.core.scaling.util import list_unscaled_variables, list_unscaled_constraints
from idaes.core.util.testing import (
    PhysicalParameterTestBlock,
    ReactionParameterTestBlock,
)
from idaes.core.base.control_volume0d import ControlVolume0DScaler


def approx(x):
    return pytest.approx(x, rel=1e-15, abs=0)


@pytest.mark.unit
def test_basic_scaling():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    # Set flag to include inherent reactions
    m.fs.pp._has_inherent_reactions = True

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)

    m.fs.cv.add_material_balances(
        balance_type=MaterialBalanceType.componentTotal, has_phase_equilibrium=False
    )

    m.fs.cv.add_energy_balances(balance_type=EnergyBalanceType.enthalpyTotal)

    m.fs.cv.add_momentum_balances(
        balance_type=MomentumBalanceType.pressureTotal, has_pressure_change=True
    )

    assert m.fs.cv.default_scaler is ControlVolume0DScaler
    cv_scaler = ControlVolume0DScaler()
    cv_scaler.default_scaling_factors["volume"] = 1 / 83
    cv_scaler.scale_model(m.fs.cv)

    # Check scaling on select variables
    assert get_scaling_factor(m.fs.cv.volume[0]) == 1 / 83
    assert get_scaling_factor(m.fs.cv.deltaP[0]) == 1 / 13

    # Check scaling on mass, energy, and pressure balances.
    for c in m.fs.cv.material_balances.values():
        # This uses the maximum nominal value in the material balance to scale
        assert get_scaling_factor(c) == 1 / 43
    for c in m.fs.cv.enthalpy_balances.values():
        # This uses the maximum nominal value in the enthalpy balance to scale
        # Note that we're expecting the user to set a scaling value or hint
        # for enthalpy_flow_terms. Otherwise this scaling may give bad values
        assert get_scaling_factor(c) == 1 / 37
    for c in m.fs.cv.pressure_balance.values():
        assert get_scaling_factor(c) == 1 / 13

    for c in m.fs.cv.inherent_reaction_stoichiometry_constraint.values():
        assert get_scaling_factor(c) == 1 / 43

    assert len(list_unscaled_variables(m, include_fixed=True)) == 0
    assert len(list_unscaled_constraints(m)) == 0


@pytest.mark.unit
def test_user_set_scaling():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp)
    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_material_balances(
        balance_type=MaterialBalanceType.componentTotal, has_phase_equilibrium=False
    )
    m.fs.cv.add_energy_balances(
        balance_type=EnergyBalanceType.enthalpyTotal,
        has_heat_transfer=True,
        has_work_transfer=True,
    )
    # add momentum balance
    m.fs.cv.add_momentum_balances(
        balance_type=MomentumBalanceType.pressureTotal, has_pressure_change=True
    )

    # The scaling factors used for this test were selected to be easy values to
    # test, they do not represent typical scaling factors.
    set_scaling_factor(m.fs.cv.heat[0], 1 / 73)
    set_scaling_factor(m.fs.cv.work[0], 1 / 79)

    cv_scaler = ControlVolume0DScaler()
    cv_scaler.default_scaling_factors["volume"] = 1 / 83

    cv_scaler.scale_model(m.fs.cv)

    # Make sure the heat and work scaling factors are set and not overwritten
    # by the defaults in calculate_scaling_factors
    assert get_scaling_factor(m.fs.cv.heat[0]) == 1 / 73
    assert get_scaling_factor(m.fs.cv.work[0]) == 1 / 79

    for c in m.fs.cv.material_balances.values():
        assert get_scaling_factor(c) == 1 / 43
    for c in m.fs.cv.enthalpy_balances.values():
        assert get_scaling_factor(c) == 1 / 37
    for c in m.fs.cv.pressure_balance.values():
        assert get_scaling_factor(c) == 1 / 13

    assert len(list_unscaled_variables(m, include_fixed=True)) == 0
    assert len(list_unscaled_constraints(m)) == 0


@pytest.mark.unit
def test_full_auto_scaling_dynamic():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyo.units.s)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)
    m.fs.cv = ControlVolume0DBlock(
        property_package=m.fs.pp, reaction_package=m.fs.rp, dynamic=True
    )
    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)
    m.fs.cv.add_reaction_blocks(has_equilibrium=True)

    m.fs.cv.add_material_balances(
        balance_type=MaterialBalanceType.componentTotal,
        has_rate_reactions=True,
        has_equilibrium_reactions=True,
        has_phase_equilibrium=True,
        has_mass_transfer=True,
    )

    m.fs.cv.add_energy_balances(
        balance_type=EnergyBalanceType.enthalpyTotal,
        has_heat_of_reaction=True,
        has_heat_transfer=True,
        has_work_transfer=True,
        has_enthalpy_transfer=True,
    )

    m.fs.cv.add_momentum_balances(
        balance_type=MomentumBalanceType.pressureTotal, has_pressure_change=True
    )

    m.discretizer = pyo.TransformationFactory("dae.finite_difference")
    m.discretizer.apply_to(m, nfe=3, wrt=m.fs.time, scheme="BACKWARD")

    cv_scaler = ControlVolume0DScaler()
    cv_scaler.default_scaling_factors["volume"] = 1 / 71

    cv_scaler.scale_model(m.fs.cv)

    cv = m.fs.cv

    for v in cv.volume.values():
        assert get_scaling_factor(v) == 1 / 71

    for v in cv.material_holdup.values():
        assert get_scaling_factor(v) == approx(10 / (47 * 71))

    for v in cv.rate_reaction_generation.values():
        assert get_scaling_factor(v) == 1 / 43

    for v in cv.rate_reaction_extent.values():
        assert get_scaling_factor(v) == 1 / 43

    for v in cv.equilibrium_reaction_generation.values():
        assert get_scaling_factor(v) == 1 / 43

    for v in cv.equilibrium_reaction_extent.values():
        assert get_scaling_factor(v) == 1 / 43

    for v in cv.mass_transfer_term.values():
        assert get_scaling_factor(v) == 1 / 43

    for v in cv.phase_fraction.values():
        assert get_scaling_factor(v) == 10

    for v in cv.energy_holdup.values():
        assert get_scaling_factor(v) == approx(10 / (41 * 71))

    for v in cv.heat.values():
        assert get_scaling_factor(v) == 1 / 37

    for v in cv.work.values():
        assert get_scaling_factor(v) == 1 / 37

    for v in cv.enthalpy_transfer.values():
        assert get_scaling_factor(v) == 1 / 37

    for v in cv.deltaP.values():
        assert get_scaling_factor(v) == 1 / 13

    for c in cv.sum_of_phase_fractions.values():
        assert get_scaling_factor(c) == 1

    for c in cv.material_holdup_calculation.values():
        assert get_scaling_factor(c) == approx(10 / (47 * 71))

    for c in cv.rate_reaction_stoichiometry_constraint.values():
        assert get_scaling_factor(c) == 1 / 43

    for c in cv.equilibrium_reaction_stoichiometry_constraint.values():
        assert get_scaling_factor(c) == 1 / 43

    for c in cv.material_balances.values():
        assert get_scaling_factor(c) == 1 / 43

    for c in cv.enthalpy_balances.values():
        assert get_scaling_factor(c) == 1 / 37

    for c in cv.energy_holdup_calculation.values():
        assert get_scaling_factor(c) == approx(10 / (41 * 71))

    for c in cv.pressure_balance.values():
        assert get_scaling_factor(c) == 1 / 13

    # Unscaled variables are DerivativeVars and unscaled constraints
    # are time discretization equations. Both should be scaled by
    # some global method
    assert len(list_unscaled_variables(m, include_fixed=True)) == 24
    assert len(list_unscaled_constraints(m)) == 18


@pytest.mark.unit
def test_full_auto_scaling_mbtype_phase():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)
    m.fs.cv = ControlVolume0DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
    )
    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)
    m.fs.cv.add_reaction_blocks(has_equilibrium=True)

    m.fs.cv.add_material_balances(
        balance_type=MaterialBalanceType.componentPhase,
        has_rate_reactions=True,
        has_equilibrium_reactions=True,
        has_phase_equilibrium=True,
        has_mass_transfer=True,
    )

    m.fs.cv.add_energy_balances(
        balance_type=EnergyBalanceType.enthalpyTotal,
        has_heat_of_reaction=True,
        has_heat_transfer=True,
        has_work_transfer=True,
        has_enthalpy_transfer=True,
    )

    m.fs.cv.add_momentum_balances(
        balance_type=MomentumBalanceType.pressureTotal, has_pressure_change=True
    )

    scaler_obj = ControlVolume0DScaler()
    scaler_obj.default_scaling_factors["volume"] = 1 / 71

    scaler_obj.scale_model(m.fs.cv)

    assert len(list_unscaled_variables(m, include_fixed=True)) == 0
    assert len(list_unscaled_constraints(m)) == 0

    cv = m.fs.cv

    # Variables
    for v in cv.volume.values():
        assert get_scaling_factor(v) == 1 / 71

    for v in cv.rate_reaction_generation.values():
        assert get_scaling_factor(v) == 1 / 43

    for v in cv.rate_reaction_extent.values():
        assert get_scaling_factor(v) == 1 / 43

    for v in cv.equilibrium_reaction_generation.values():
        assert get_scaling_factor(v) == 1 / 43

    for v in cv.equilibrium_reaction_extent.values():
        assert get_scaling_factor(v) == 1 / 43

    for v in cv.mass_transfer_term.values():
        assert get_scaling_factor(v) == 1 / 43

    for v in cv.heat.values():
        assert get_scaling_factor(v) == 1 / 37

    for v in cv.work.values():
        assert get_scaling_factor(v) == 1 / 37

    for v in cv.enthalpy_transfer.values():
        assert get_scaling_factor(v) == 1 / 37

    for v in cv.deltaP.values():
        assert get_scaling_factor(v) == 1 / 13

    # Constraints
    for c in cv.rate_reaction_stoichiometry_constraint.values():
        assert get_scaling_factor(c) == 1 / 43

    for c in cv.equilibrium_reaction_stoichiometry_constraint.values():
        assert get_scaling_factor(c) == 1 / 43

    for c in cv.material_balances.values():
        assert get_scaling_factor(c) == 1 / 43

    for c in cv.enthalpy_balances.values():
        assert get_scaling_factor(c) == 1 / 37

    for c in cv.pressure_balance.values():
        assert get_scaling_factor(c) == 1 / 13

    for c in cv.phase_equilibrium_generation.values():
        assert get_scaling_factor(c) == 1 / 43


# TODO test element balance when implemented
# @pytest.mark.unit
# def test_full_auto_scaling_mbtype_element():
#     m = pyo.ConcreteModel()
#     m.fs = FlowsheetBlock(dynamic=True, time_units=pyo.units.s)
#     m.fs.pp = PhysicalParameterTestBlock()
#     m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)
#     m.fs.cv = ControlVolume0DBlock(
#         property_package=m.fs.pp, reaction_package=m.fs.rp, dynamic=True
#     )
#     m.fs.cv.add_geometry()
#     m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
#     m.fs.cv.add_reaction_blocks(has_equilibrium=False)

#     m.fs.cv.add_total_element_balances(has_mass_transfer=True)

#     m.discretizer = pyo.TransformationFactory("dae.finite_difference")
#     m.discretizer.apply_to(m, nfe=3, wrt=m.fs.time, scheme="BACKWARD")

#     iscale.calculate_scaling_factors(m)

#     # check that all variables have scaling factors
#     unscaled_var_list = list(iscale.unscaled_variables_generator(m))
#     # Unscaled variables are:
#     # cp at inlet and outlet (4 time points)
#     assert len(unscaled_var_list) == 8
#     # check that all constraints have been scaled
#     unscaled_constraint_list = list(iscale.unscaled_constraints_generator(m))
#     assert len(unscaled_constraint_list) == 0
