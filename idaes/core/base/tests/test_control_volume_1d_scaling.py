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
"""
Tests for ControlVolume0D scaling.

Author: John Eslick
"""
import pytest
import pyomo.environ as pyo
from idaes.core import (
    ControlVolume1DBlock,
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


# -----------------------------------------------------------------------------
# Basic tests
@pytest.mark.unit
def test_basic_scaling():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    # Set flag to include inherent reactions
    m.fs.pp._has_inherent_reactions = True

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)

    m.fs.cv.add_material_balances(
        balance_type=MaterialBalanceType.componentTotal, has_phase_equilibrium=False
    )

    m.fs.cv.add_energy_balances(balance_type=EnergyBalanceType.enthalpyTotal)

    m.fs.cv.add_momentum_balances(
        balance_type=MomentumBalanceType.pressureTotal, has_pressure_change=True
    )

    m.fs.cv.apply_transformation()

    iscale.calculate_scaling_factors(m)

    # check scaling on select variables
    assert iscale.get_scaling_factor(m.fs.cv.area) == 1
    for (t, x), v in m.fs.cv.deltaP.items():
        assert (
            iscale.get_scaling_factor(v) == 1040
        )  # 10x the properties pressure scaling factor
    for t in m.fs.time:
        for x in m.fs.cv.length_domain:
            assert iscale.get_scaling_factor(m.fs.cv.properties[t, x].flow_vol) == 100

    # check scaling on mass, energy, and pressure balances.
    for c in m.fs.cv.material_balances.values():
        # this uses the minimum material flow term scale
        assert iscale.get_constraint_transform_applied_scaling_factor(c) == 112
    for c in m.fs.cv.enthalpy_balances.values():
        # this uses the minimum enthalpy flow term scale
        assert iscale.get_constraint_transform_applied_scaling_factor(c) == 110
    for c in m.fs.cv.pressure_balance.values():
        # This uses the inlet pressure scale
        assert iscale.get_constraint_transform_applied_scaling_factor(c) == 104


@pytest.mark.unit
def test_user_set_scaling():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )
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

    m.fs.cv.apply_transformation()

    # The scaling factors used for this test were selected to be easy values to
    # test, they do not represent typical scaling factors.
    iscale.set_scaling_factor(m.fs.cv.heat, 11)
    iscale.set_scaling_factor(m.fs.cv.work, 12)
    iscale.set_scaling_factor(m.fs.cv.heat[0, 0], 17)
    iscale.calculate_scaling_factors(m)
    # Make sure the heat and work scaling factors are set and not overwritten
    # by the defaults in calculate_scaling_factors
    assert iscale.get_scaling_factor(m.fs.cv.heat) == 11
    assert iscale.get_scaling_factor(m.fs.cv.work) == 12
    # Make sure scaling factor is propagated but does not overwrite index specific factor
    for t in m.fs.time:
        for x in m.fs.cv.length_domain:
            if t == 0 and x == 0:
                assert iscale.get_scaling_factor(m.fs.cv.heat[t, x]) == 17
            else:
                assert iscale.get_scaling_factor(m.fs.cv.heat[t, x]) == 11
            assert iscale.get_scaling_factor(m.fs.cv.work[t, x]) == 12


@pytest.mark.unit
def test_full_auto_scaling():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)
    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
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

    m.fs.cv.apply_transformation()

    iscale.calculate_scaling_factors(m)

    # check that all variables have scaling factors
    unscaled_var_list = list(iscale.unscaled_variables_generator(m))
    # Unscaled variables are:
    # rate_reaction_extent (2 reactions, 11 spatial points)
    # equilibrium_reaction_extent  (2 reactions, 11 spatial points)
    # cp  (11 spatial points)
    assert len(unscaled_var_list) == 55
    # check that all constraints have been scaled
    unscaled_constraint_list = list(iscale.unscaled_constraints_generator(m))
    assert len(unscaled_constraint_list) == 0


@pytest.mark.unit
def test_full_auto_scaling_dynamic():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyo.units.s)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)
    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
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

    m.fs.cv.apply_transformation()
    m.discretizer = pyo.TransformationFactory("dae.finite_difference")
    m.discretizer.apply_to(m, nfe=3, wrt=m.fs.time, scheme="BACKWARD")

    iscale.calculate_scaling_factors(m)

    # check that all variables have scaling factors
    unscaled_var_list = list(iscale.unscaled_variables_generator(m))
    # Unscaled variables are:
    # rate_reaction_extent (2 reactions, 44 time & space points)
    # equilibrium_reaction_extent  (2 reactions, 44 time & space points)
    # cp  (44 time & space points)
    assert len(unscaled_var_list) == 220
    # check that all constraints have been scaled
    unscaled_constraint_list = list(iscale.unscaled_constraints_generator(m))
    assert len(unscaled_constraint_list) == 0


@pytest.mark.unit
def test_full_auto_scaling_mbtype_phase():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)
    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
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

    m.fs.cv.apply_transformation()

    iscale.calculate_scaling_factors(m)

    # check that all variables have scaling factors
    unscaled_var_list = list(iscale.unscaled_variables_generator(m))
    # Unscaled variables are:
    # rate_reaction_extent (2 reactions, 11 spatial points)
    # equilibrium_reaction_extent  (2 reactions, 11 spatial points)
    # phase_equilibrium_generation  (2 reactions, 11 spatial points)
    # cp  (11 spatial points)
    assert len(unscaled_var_list) == 77
    # check that all constraints have been scaled
    unscaled_constraint_list = list(iscale.unscaled_constraints_generator(m))
    assert len(unscaled_constraint_list) == 0


@pytest.mark.unit
def test_full_auto_scaling_mbtype_element():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyo.units.s)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)
    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )
    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    m.fs.cv.add_total_element_balances(has_mass_transfer=True)

    m.fs.cv.apply_transformation()
    m.discretizer = pyo.TransformationFactory("dae.finite_difference")
    m.discretizer.apply_to(m, nfe=3, wrt=m.fs.time, scheme="BACKWARD")

    iscale.calculate_scaling_factors(m)

    # check that all variables have scaling factors
    unscaled_var_list = list(iscale.unscaled_variables_generator(m))
    # Unscaled variables are:
    # cp  (44 space and time points)
    assert len(unscaled_var_list) == 44
    # check that all constraints have been scaled
    unscaled_constraint_list = list(iscale.unscaled_constraints_generator(m))
    assert len(unscaled_constraint_list) == 0
