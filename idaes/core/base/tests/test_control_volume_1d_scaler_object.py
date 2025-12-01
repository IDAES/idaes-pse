#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Tests for ControlVolume1D scaler object.

Author: Douglas Allan
"""
import pytest
import re
import pyomo.environ as pyo
from idaes.core import (
    ControlVolume1DBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    FlowsheetBlock,
    FlowDirection,
)
from idaes.core.base.control_volume1d import ControlVolume1DScaler
from idaes.core.util.testing import (
    PhysicalParameterTestBlock,
    ReactionParameterTestBlock,
)
from idaes.core.util.exceptions import BurntToast
from idaes.core.scaling import get_scaling_factor
from idaes.core.scaling.util import list_unscaled_variables, list_unscaled_constraints


def approx(x):
    return pytest.approx(x, rel=1e-15, abs=0)


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

    assert m.fs.cv.default_scaler is ControlVolume1DScaler

    scaler_obj = ControlVolume1DScaler()
    scaler_obj.default_scaling_factors["area"] = 1 / 83
    scaler_obj.default_scaling_factors["length"] = 1 / 599
    scaler_obj.scale_model(m.fs.cv)

    # check scaling on select variables
    assert get_scaling_factor(m.fs.cv.area) == 1 / 83

    for v in m.fs.cv._flow_terms.values():
        assert get_scaling_factor(v) == approx(1 / 43)

    for v in m.fs.cv.material_flow_dx.values():
        assert get_scaling_factor(v) == approx(1 / 43)

    for v in m.fs.cv._enthalpy_flow.values():
        assert get_scaling_factor(v) == approx(1 / 37)

    for v in m.fs.cv.enthalpy_flow_dx.values():
        assert get_scaling_factor(v) == approx(1 / 37)

    for v in m.fs.cv.deltaP.values():
        assert get_scaling_factor(v) == approx(599 / 13)
    for v in m.fs.cv.pressure_dx.values():
        get_scaling_factor(v) == approx(1 / 13)

    for t in m.fs.time:
        for x in m.fs.cv.length_domain:
            assert get_scaling_factor(m.fs.cv.properties[t, x].flow_vol) == 1 / 3

    # check scaling on mass, energy, and pressure balances.
    for c in m.fs.cv.material_balances.values():
        assert get_scaling_factor(c) == 1 / 43
    for c in m.fs.cv.enthalpy_balances.values():
        # this uses the minimum enthalpy flow term scale
        assert get_scaling_factor(c) == 1 / 37
    for c in m.fs.cv.pressure_balance.values():
        # This uses the inlet pressure scale
        assert get_scaling_factor(c) == 1 / 13


@pytest.mark.unit
def test_user_set_scaling_FD_forward():
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

    scaler_obj = ControlVolume1DScaler()

    # Scaling factors should get propagated from inlet state block to the other
    # state blocks
    scaler_obj.set_component_scaling_factor(
        m.fs.cv.properties[0, 0].temperature, 1 / 1009
    )
    scaler_obj.set_component_scaling_factor(m.fs.cv.properties[0, 0].pressure, 1 / 1301)
    scaler_obj.set_component_scaling_factor(
        m.fs.cv.properties[0, 0].flow_mol_phase_comp["p1", "c1"], 1 / 1117
    )
    scaler_obj.set_component_scaling_factor(
        m.fs.cv.properties[0, 0].flow_mol_phase_comp["p2", "c1"], 1 / 1931
    )
    scaler_obj.set_component_scaling_factor(
        m.fs.cv.properties[0, 0].flow_mol_phase_comp["p1", "c2"], 1 / 1549
    )
    scaler_obj.set_component_scaling_factor(
        m.fs.cv.properties[0, 0].flow_mol_phase_comp["p2", "c2"], 1 / 1999
    )

    scaler_obj.set_component_scaling_factor(m.fs.cv.properties[0, 1].pressure, 1 / 2251)

    scaler_obj.default_scaling_factors["area"] = 1 / 83
    scaler_obj.default_scaling_factors["length"] = 1 / 599

    # The scaling factors used for this test were selected to be easy values to
    # test, they do not represent typical scaling factors.
    for v in m.fs.cv.heat.values():
        scaler_obj.set_component_scaling_factor(v, 1 / 73)
    for v in m.fs.cv.work.values():
        scaler_obj.set_component_scaling_factor(v, 1 / 79)
    scaler_obj.set_component_scaling_factor(m.fs.cv.heat[0, 0], 17, overwrite=True)

    scaler_obj.scale_model(m.fs.cv)
    for t in m.fs.time:
        for x in m.fs.cv.length_domain:
            if t == 0 and x == 0:
                assert get_scaling_factor(m.fs.cv.heat[t, x]) == 17
            else:
                assert get_scaling_factor(m.fs.cv.heat[t, x]) == 1 / 73
            assert get_scaling_factor(m.fs.cv.work[t, x]) == 1 / 79

            assert get_scaling_factor(m.fs.cv.properties[t, x].temperature) == 1 / 1009
            if t == 0 and x == 1:
                # Propagation from inlet state block shouldn't overwrite existing scaling factor
                assert get_scaling_factor(m.fs.cv.properties[t, x].pressure) == 1 / 2251
            else:
                assert get_scaling_factor(m.fs.cv.properties[t, x].pressure) == 1 / 1301
            assert (
                get_scaling_factor(
                    m.fs.cv.properties[t, x].flow_mol_phase_comp["p1", "c1"]
                )
                == 1 / 1117
            )
            assert (
                get_scaling_factor(
                    m.fs.cv.properties[t, x].flow_mol_phase_comp["p2", "c1"]
                )
                == 1 / 1931
            )
            assert (
                get_scaling_factor(
                    m.fs.cv.properties[t, x].flow_mol_phase_comp["p1", "c2"]
                )
                == 1 / 1549
            )
            assert (
                get_scaling_factor(
                    m.fs.cv.properties[t, x].flow_mol_phase_comp["p2", "c2"]
                )
                == 1 / 1999
            )


@pytest.mark.unit
def test_user_set_scaling_FD_backward():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )
    m.fs.cv.add_geometry(flow_direction=FlowDirection.backward)
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

    scaler_obj = ControlVolume1DScaler()

    scaler_obj.set_component_scaling_factor(
        m.fs.cv.properties[0, 0].temperature, 1 / 1009
    )
    scaler_obj.set_component_scaling_factor(m.fs.cv.properties[0, 0].pressure, 1 / 1301)
    scaler_obj.set_component_scaling_factor(
        m.fs.cv.properties[0, 0].flow_mol_phase_comp["p1", "c1"], 1 / 1117
    )
    scaler_obj.set_component_scaling_factor(
        m.fs.cv.properties[0, 0].flow_mol_phase_comp["p2", "c1"], 1 / 1931
    )
    scaler_obj.set_component_scaling_factor(
        m.fs.cv.properties[0, 0].flow_mol_phase_comp["p1", "c2"], 1 / 1549
    )
    scaler_obj.set_component_scaling_factor(
        m.fs.cv.properties[0, 0].flow_mol_phase_comp["p2", "c2"], 1 / 1999
    )

    scaler_obj.set_component_scaling_factor(m.fs.cv.properties[0, 1].pressure, 1 / 2251)

    scaler_obj.default_scaling_factors["area"] = 1 / 83
    scaler_obj.default_scaling_factors["length"] = 1 / 599

    # The scaling factors used for this test were selected to be easy values to
    # test, they do not represent typical scaling factors.
    for v in m.fs.cv.heat.values():
        scaler_obj.set_component_scaling_factor(v, 1 / 73)
    for v in m.fs.cv.work.values():
        scaler_obj.set_component_scaling_factor(v, 1 / 79)
    scaler_obj.set_component_scaling_factor(m.fs.cv.heat[0, 0], 17, overwrite=True)

    scaler_obj.scale_model(m.fs.cv)
    # Since the inlet state block is [0, 1] and the only scaling factor set on
    # that block is for pressure, only that scaling factor is propagated to the
    # other state blocks. All other state variables use default values
    for t in m.fs.time:
        for x in m.fs.cv.length_domain:
            if t == 0 and x == 0:
                assert get_scaling_factor(m.fs.cv.heat[t, x]) == 17
            else:
                assert get_scaling_factor(m.fs.cv.heat[t, x]) == 1 / 73
            assert get_scaling_factor(m.fs.cv.work[t, x]) == 1 / 79

            if t == 0 and x == 0:
                # Scaling factors were specifically set at [0, 0], so those
                # shouldn't get overwritten
                assert (
                    get_scaling_factor(m.fs.cv.properties[t, x].temperature) == 1 / 1009
                )
                assert get_scaling_factor(m.fs.cv.properties[t, x].pressure) == 1 / 1301
                assert (
                    get_scaling_factor(
                        m.fs.cv.properties[t, x].flow_mol_phase_comp["p1", "c1"]
                    )
                    == 1 / 1117
                )
                assert (
                    get_scaling_factor(
                        m.fs.cv.properties[t, x].flow_mol_phase_comp["p2", "c1"]
                    )
                    == 1 / 1931
                )
                assert (
                    get_scaling_factor(
                        m.fs.cv.properties[t, x].flow_mol_phase_comp["p1", "c2"]
                    )
                    == 1 / 1549
                )
                assert (
                    get_scaling_factor(
                        m.fs.cv.properties[t, x].flow_mol_phase_comp["p2", "c2"]
                    )
                    == 1 / 1999
                )
            else:
                assert (
                    get_scaling_factor(m.fs.cv.properties[t, x].temperature) == 1 / 17
                )
                assert get_scaling_factor(m.fs.cv.properties[t, x].pressure) == 1 / 2251
                assert (
                    get_scaling_factor(
                        m.fs.cv.properties[t, x].flow_mol_phase_comp["p1", "c1"]
                    )
                    == 1 / 7
                )
                assert (
                    get_scaling_factor(
                        m.fs.cv.properties[t, x].flow_mol_phase_comp["p2", "c1"]
                    )
                    == 1 / 7
                )
                assert (
                    get_scaling_factor(
                        m.fs.cv.properties[t, x].flow_mol_phase_comp["p1", "c2"]
                    )
                    == 1 / 7
                )
                assert (
                    get_scaling_factor(
                        m.fs.cv.properties[t, x].flow_mol_phase_comp["p2", "c2"]
                    )
                    == 1 / 7
                )


@pytest.mark.unit
def test_user_set_scaling_FD_not_set():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    scaler_obj = ControlVolume1DScaler()

    with pytest.raises(
        RuntimeError,
        match=re.escape(
            "Scaler called on ControlVolume1D without a flow "
            "direction set. The unit model containing the ControlVolume1D "
            "should use the add_geometry method to specify a flow direction "
            "as part of model construction."
        ),
    ):
        scaler_obj.scale_model(m.fs.cv)


@pytest.mark.unit
def test_user_set_scaling_FD_burnt_toast():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )
    m.fs.cv._flow_direction = "foo"

    scaler_obj = ControlVolume1DScaler()

    with pytest.raises(
        BurntToast,
        match=re.escape(
            "Unknown flow direction specified. This indicates "
            "a new flow direction was added without support being "
            "extended to the scaler. Please contact the IDAES "
            "development team with this error."
        ),
    ):
        scaler_obj.scale_model(m.fs.cv)


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

    scaler_obj = ControlVolume1DScaler()
    scaler_obj.default_scaling_factors["area"] = 1 / 83
    scaler_obj.default_scaling_factors["length"] = 1 / 599
    scaler_obj.scale_model(m.fs.cv)

    cv = m.fs.cv

    # Variables
    assert get_scaling_factor(cv.length) == 1 / 599
    for v in cv.area.values():
        assert get_scaling_factor(v) == 1 / 83

    for v in cv._flow_terms.values():
        assert get_scaling_factor(v) == 1 / 43

    for v in cv._enthalpy_flow.values():
        assert get_scaling_factor(v) == 1 / 37

    for v in cv.material_flow_dx.values():
        assert get_scaling_factor(v) == 1 / 43

    for v in cv.enthalpy_flow_dx.values():
        assert get_scaling_factor(v) == 1 / 37

    for v in cv.pressure_dx.values():
        assert get_scaling_factor(v) == 1 / 13

    for v in cv.rate_reaction_generation.values():
        assert get_scaling_factor(v) == approx(599 / 43)

    for v in cv.rate_reaction_extent.values():
        assert get_scaling_factor(v) == approx(599 / 43)

    for v in cv.equilibrium_reaction_generation.values():
        assert get_scaling_factor(v) == approx(599 / 43)

    for v in cv.equilibrium_reaction_extent.values():
        assert get_scaling_factor(v) == approx(599 / 43)

    for v in cv.mass_transfer_term.values():
        assert get_scaling_factor(v) == approx(599 / 43)

    for v in cv.heat.values():
        assert get_scaling_factor(v) == approx(599 / 37)

    for v in cv.work.values():
        assert get_scaling_factor(v) == approx(599 / 37)

    for v in cv.enthalpy_transfer.values():
        assert get_scaling_factor(v) == approx(599 / 37)

    for v in cv.deltaP.values():
        assert get_scaling_factor(v) == approx(599 / 13)

    # Constraints
    for v in cv.material_flow_dx_disc_eq.values():
        assert get_scaling_factor(v) == 1 / 43

    for v in cv.enthalpy_flow_dx_disc_eq.values():
        assert get_scaling_factor(v) == 1 / 37

    for v in cv.pressure_dx_disc_eq.values():
        assert get_scaling_factor(v) == 1 / 13

    for c in cv.rate_reaction_stoichiometry_constraint.values():
        assert get_scaling_factor(c) == approx(599 / 43)

    for c in cv.equilibrium_reaction_stoichiometry_constraint.values():
        assert get_scaling_factor(c) == approx(599 / 43)

    for c in cv.material_balances.values():
        assert get_scaling_factor(c) == approx(1 / 43)

    for c in cv.enthalpy_balances.values():
        assert get_scaling_factor(c) == approx(1 / 37)

    for c in cv.pressure_balance.values():
        assert get_scaling_factor(c) == approx(1 / 13)

    assert len(list_unscaled_variables(m, include_fixed=True)) == 0
    assert len(list_unscaled_constraints(m)) == 0


@pytest.mark.unit
def test_full_auto_scaling_dynamic():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_set=[0, 1], time_units=pyo.units.s)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)
    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.collocation",
        transformation_scheme="LAGRANGE-LEGENDRE",
        finite_elements=4,
        collocation_points=3,
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

    scaler_obj = ControlVolume1DScaler()
    scaler_obj.default_scaling_factors["area"] = 1 / 83
    scaler_obj.default_scaling_factors["length"] = 1 / 599

    scaler_obj.scale_model(m.fs.cv)

    # Four finite elements with three collocation points on the element's
    # interior plus two boundary points (which overlap with other finite elements)
    # That results in 4 points per finite element plus one point at the end of the
    # final finite element.
    n_points_colloc = 4 * (3 + 1) + 1
    # Unscaled variables are material accumulation and energy accumulation DerivativeVars
    # Two phases with two components. Energy accumulation is also indexed by phase

    assert len(list_unscaled_variables(m, include_fixed=True)) == 4 * (
        2 * 2 * n_points_colloc + 2 * n_points_colloc
    )
    # There are variables for all 4 timepoints but only constraints for 3 timepoints
    assert len(list_unscaled_constraints(m)) == 3 * (
        2 * 2 * n_points_colloc + 2 * n_points_colloc
    )


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

    scaler_obj = ControlVolume1DScaler()
    scaler_obj.default_scaling_factors["area"] = 1 / 83
    scaler_obj.default_scaling_factors["length"] = 1 / 599

    scaler_obj.scale_model(m.fs.cv)

    assert len(list_unscaled_variables(m, include_fixed=True)) == 0
    assert len(list_unscaled_constraints(m)) == 0

    cv = m.fs.cv

    # Variables
    assert get_scaling_factor(cv.length) == 1 / 599
    for v in cv.area.values():
        assert get_scaling_factor(v) == 1 / 83

    for v in cv._flow_terms.values():
        assert get_scaling_factor(v) == 1 / 43

    for v in cv._enthalpy_flow.values():
        assert get_scaling_factor(v) == 1 / 37

    for v in cv.material_flow_dx.values():
        assert get_scaling_factor(v) == 1 / 43

    for v in cv.enthalpy_flow_dx.values():
        assert get_scaling_factor(v) == 1 / 37

    for v in cv.pressure_dx.values():
        assert get_scaling_factor(v) == 1 / 13

    for v in cv.rate_reaction_generation.values():
        assert get_scaling_factor(v) == approx(599 / 43)

    for v in cv.rate_reaction_extent.values():
        assert get_scaling_factor(v) == approx(599 / 43)

    for v in cv.equilibrium_reaction_generation.values():
        assert get_scaling_factor(v) == approx(599 / 43)

    for v in cv.equilibrium_reaction_extent.values():
        assert get_scaling_factor(v) == approx(599 / 43)

    for v in cv.mass_transfer_term.values():
        assert get_scaling_factor(v) == approx(599 / 43)

    for v in cv.heat.values():
        assert get_scaling_factor(v) == approx(599 / 37)

    for v in cv.work.values():
        assert get_scaling_factor(v) == approx(599 / 37)

    for v in cv.enthalpy_transfer.values():
        assert get_scaling_factor(v) == approx(599 / 37)

    for v in cv.deltaP.values():
        assert get_scaling_factor(v) == approx(599 / 13)

    # Constraints
    for v in cv.material_flow_dx_disc_eq.values():
        assert get_scaling_factor(v) == 1 / 43

    for v in cv.enthalpy_flow_dx_disc_eq.values():
        assert get_scaling_factor(v) == 1 / 37

    for v in cv.pressure_dx_disc_eq.values():
        assert get_scaling_factor(v) == 1 / 13

    for c in cv.rate_reaction_stoichiometry_constraint.values():
        assert get_scaling_factor(c) == approx(599 / 43)

    for c in cv.equilibrium_reaction_stoichiometry_constraint.values():
        assert get_scaling_factor(c) == approx(599 / 43)

    for c in cv.material_balances.values():
        assert get_scaling_factor(c) == approx(1 / 43)

    for c in cv.enthalpy_balances.values():
        assert get_scaling_factor(c) == approx(1 / 37)

    for c in cv.pressure_balance.values():
        assert get_scaling_factor(c) == approx(1 / 13)

    for c in cv.phase_equilibrium_generation.values():
        assert get_scaling_factor(c) == approx(1 / 43)


@pytest.mark.unit
def test_auto_scaling_just_material_balance():
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

    m.fs.cv.apply_transformation()

    scaler_obj = ControlVolume1DScaler()
    scaler_obj.default_scaling_factors["area"] = 1 / 83
    scaler_obj.default_scaling_factors["length"] = 1 / 599
    scaler_obj.scale_model(m.fs.cv)

    cv = m.fs.cv

    # Variables
    assert get_scaling_factor(cv.length) == 1 / 599
    for v in cv.area.values():
        assert get_scaling_factor(v) == 1 / 83

    for v in cv._flow_terms.values():
        assert get_scaling_factor(v) == 1 / 43

    for v in cv.material_flow_dx.values():
        assert get_scaling_factor(v) == 1 / 43

    for v in cv.rate_reaction_generation.values():
        assert get_scaling_factor(v) == approx(599 / 43)

    for v in cv.rate_reaction_extent.values():
        assert get_scaling_factor(v) == approx(599 / 43)

    for v in cv.equilibrium_reaction_generation.values():
        assert get_scaling_factor(v) == approx(599 / 43)

    for v in cv.equilibrium_reaction_extent.values():
        assert get_scaling_factor(v) == approx(599 / 43)

    for v in cv.mass_transfer_term.values():
        assert get_scaling_factor(v) == approx(599 / 43)

    # Constraints
    for v in cv.material_flow_dx_disc_eq.values():
        assert get_scaling_factor(v) == 1 / 43
    for c in cv.rate_reaction_stoichiometry_constraint.values():
        assert get_scaling_factor(c) == approx(599 / 43)

    for c in cv.equilibrium_reaction_stoichiometry_constraint.values():
        assert get_scaling_factor(c) == approx(599 / 43)

    for c in cv.material_balances.values():
        assert get_scaling_factor(c) == approx(1 / 43)

    assert len(list_unscaled_variables(m, include_fixed=True)) == 0
    assert len(list_unscaled_constraints(m)) == 0


@pytest.mark.unit
def test_auto_scaling_enthalpy_and_pressure_balance():
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

    m.fs.cv.add_energy_balances(
        balance_type=EnergyBalanceType.enthalpyTotal,
        has_heat_of_reaction=False,
        has_heat_transfer=True,
        has_work_transfer=True,
        has_enthalpy_transfer=True,
    )

    m.fs.cv.add_momentum_balances(
        balance_type=MomentumBalanceType.pressureTotal, has_pressure_change=True
    )

    m.fs.cv.apply_transformation()

    scaler_obj = ControlVolume1DScaler()
    scaler_obj.default_scaling_factors["area"] = 1 / 83
    scaler_obj.default_scaling_factors["length"] = 1 / 599
    scaler_obj.scale_model(m.fs.cv)

    cv = m.fs.cv

    # Variables
    assert get_scaling_factor(cv.length) == 1 / 599
    for v in cv.area.values():
        assert get_scaling_factor(v) == 1 / 83

    for v in cv._enthalpy_flow.values():
        assert get_scaling_factor(v) == 1 / 37

    for v in cv.enthalpy_flow_dx.values():
        assert get_scaling_factor(v) == 1 / 37

    for v in cv.pressure_dx.values():
        assert get_scaling_factor(v) == 1 / 13

    for v in cv.heat.values():
        assert get_scaling_factor(v) == approx(599 / 37)

    for v in cv.work.values():
        assert get_scaling_factor(v) == approx(599 / 37)

    for v in cv.enthalpy_transfer.values():
        assert get_scaling_factor(v) == approx(599 / 37)

    for v in cv.deltaP.values():
        assert get_scaling_factor(v) == approx(599 / 13)

    # Constraints
    for v in cv.enthalpy_flow_dx_disc_eq.values():
        assert get_scaling_factor(v) == 1 / 37

    for v in cv.pressure_dx_disc_eq.values():
        assert get_scaling_factor(v) == 1 / 13

    for c in cv.enthalpy_balances.values():
        assert get_scaling_factor(c) == approx(1 / 37)

    for c in cv.pressure_balance.values():
        assert get_scaling_factor(c) == approx(1 / 13)

    assert len(list_unscaled_variables(m, include_fixed=True)) == 0
    assert len(list_unscaled_constraints(m)) == 0


# TODO element balance
# @pytest.mark.unit
# def test_full_auto_scaling_mbtype_element():
#     m = pyo.ConcreteModel()
#     m.fs = FlowsheetBlock(dynamic=True, time_units=pyo.units.s)
#     m.fs.pp = PhysicalParameterTestBlock()
#     m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)
#     m.fs.cv = ControlVolume1DBlock(
#         property_package=m.fs.pp,
#         reaction_package=m.fs.rp,
#         transformation_method="dae.finite_difference",
#         transformation_scheme="BACKWARD",
#         finite_elements=10,
#     )
#     m.fs.cv.add_geometry()
#     m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
#     m.fs.cv.add_reaction_blocks(has_equilibrium=False)

#     m.fs.cv.add_total_element_balances(has_mass_transfer=True)

#     m.fs.cv.apply_transformation()
#     m.discretizer = pyo.TransformationFactory("dae.finite_difference")
#     m.discretizer.apply_to(m, nfe=3, wrt=m.fs.time, scheme="BACKWARD")

#     iscale.calculate_scaling_factors(m)

#     # check that all variables have scaling factors
#     unscaled_var_list = list(iscale.unscaled_variables_generator(m))
#     # Unscaled variables are:
#     # cp  (44 space and time points)
#     assert len(unscaled_var_list) == 44
#     # check that all constraints have been scaled
#     unscaled_constraint_list = list(iscale.unscaled_constraints_generator(m))
#     assert len(unscaled_constraint_list) == 0
