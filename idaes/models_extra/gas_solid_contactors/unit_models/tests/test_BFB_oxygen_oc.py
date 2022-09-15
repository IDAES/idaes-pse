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
Test for ControlVolumeBlockData, and for initializing the
bubbling fluidized bed module

Author: Chinedu Okoli
"""

import pytest

from pyomo.environ import (
    ConcreteModel,
    check_optimal_termination,
    SolverStatus,
    value,
    Var,
    Constraint,
)
from pyomo.util.check_units import assert_units_consistent
from pyomo.common.config import ConfigBlock
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
)
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
    unused_variables_set,
)
from idaes.core.util.testing import initialization_tester
from idaes.core.solvers import get_solver
from idaes.core.util import scaling as iscale
from idaes.core.util.exceptions import InitializationError

from idaes.models_extra.gas_solid_contactors.unit_models.bubbling_fluidized_bed import (
    BubblingFluidizedBed,
)
from idaes.models_extra.gas_solid_contactors.properties.oxygen_iron_OC_oxidation.gas_phase_thermo import (
    GasPhaseParameterBlock,
)
from idaes.models_extra.gas_solid_contactors.properties.oxygen_iron_OC_oxidation.solid_phase_thermo import (
    SolidPhaseParameterBlock,
)
from idaes.models_extra.gas_solid_contactors.properties.oxygen_iron_OC_oxidation.hetero_reactions import (
    HeteroReactionParameterBlock,
)

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


# -----------------------------------------------------------------------------
@pytest.mark.unit
def test_config():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    # Set up thermo props and reaction props
    m.fs.gas_properties = GasPhaseParameterBlock()
    m.fs.solid_properties = SolidPhaseParameterBlock()
    m.fs.hetero_reactions = HeteroReactionParameterBlock(
        solid_property_package=m.fs.solid_properties,
        gas_property_package=m.fs.gas_properties,
    )

    m.fs.unit = BubblingFluidizedBed(
        flow_type="co_current",
        finite_elements=5,
        transformation_method="dae.finite_difference",
        gas_phase_config={"property_package": m.fs.gas_properties},
        solid_phase_config={
            "property_package": m.fs.solid_properties,
            "reaction_package": m.fs.hetero_reactions,
        },
    )

    # Check unit config arguments
    assert len(m.fs.unit.config) == 14
    assert isinstance(m.fs.unit.config.gas_phase_config, ConfigBlock)
    assert isinstance(m.fs.unit.config.solid_phase_config, ConfigBlock)

    assert m.fs.unit.config.finite_elements == 5
    assert m.fs.unit.config.length_domain_set == [0.0, 1.0]
    assert m.fs.unit.config.transformation_method == "dae.finite_difference"
    assert m.fs.unit.config.transformation_scheme == "BACKWARD"
    assert m.fs.unit.config.collocation_points == 3
    assert m.fs.unit.config.flow_type == "co_current"
    assert m.fs.unit.config.material_balance_type == MaterialBalanceType.componentTotal
    assert m.fs.unit.config.energy_balance_type == EnergyBalanceType.enthalpyTotal
    assert m.fs.unit.config.momentum_balance_type == MomentumBalanceType.pressureTotal
    assert m.fs.unit.config.has_pressure_change is True

    # Check gas phase config arguments
    assert len(m.fs.unit.config.gas_phase_config) == 7
    assert m.fs.unit.config.gas_phase_config.has_equilibrium_reactions is False
    assert m.fs.unit.config.gas_phase_config.property_package is m.fs.gas_properties
    assert m.fs.unit.config.gas_phase_config.reaction_package is None

    # Check solid phase config arguments
    assert len(m.fs.unit.config.solid_phase_config) == 7
    assert m.fs.unit.config.solid_phase_config.has_equilibrium_reactions is False
    assert m.fs.unit.config.solid_phase_config.property_package is m.fs.solid_properties
    assert m.fs.unit.config.solid_phase_config.reaction_package is m.fs.hetero_reactions


# -----------------------------------------------------------------------------
class TestIronOC(object):
    @pytest.fixture(scope="class")
    def iron_oc(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        # Set up thermo props and reaction props
        m.fs.gas_properties = GasPhaseParameterBlock()
        m.fs.solid_properties = SolidPhaseParameterBlock()
        m.fs.hetero_reactions = HeteroReactionParameterBlock(
            solid_property_package=m.fs.solid_properties,
            gas_property_package=m.fs.gas_properties,
        )

        m.fs.unit = BubblingFluidizedBed(
            flow_type="co_current",
            finite_elements=5,
            transformation_method="dae.finite_difference",
            gas_phase_config={"property_package": m.fs.gas_properties},
            solid_phase_config={
                "property_package": m.fs.solid_properties,
                "reaction_package": m.fs.hetero_reactions,
            },
        )

        # Fix geometry variables
        m.fs.unit.number_orifice.fix(2500)  # [-]
        m.fs.unit.bed_diameter.fix(6.5)  # m
        m.fs.unit.bed_height.fix(5)  # m

        # Fix inlet port variables for gas and solid
        m.fs.unit.gas_inlet.flow_mol[0].fix(1567.79)  # mol/s
        m.fs.unit.gas_inlet.temperature[0].fix(373)  # K
        m.fs.unit.gas_inlet.pressure[0].fix(1.86e5)  # Pa = 1E5 bar
        m.fs.unit.gas_inlet.mole_frac_comp[0, "O2"].fix(0.2095)
        m.fs.unit.gas_inlet.mole_frac_comp[0, "N2"].fix(0.7808)
        m.fs.unit.gas_inlet.mole_frac_comp[0, "CO2"].fix(0.0004)
        m.fs.unit.gas_inlet.mole_frac_comp[0, "H2O"].fix(0.0093)

        m.fs.unit.solid_inlet.flow_mass[0].fix(1230.865)  # kg/s
        m.fs.unit.solid_inlet.particle_porosity[0].fix(0.27)  # (-)
        m.fs.unit.solid_inlet.temperature[0].fix(1173.9)  # K
        m.fs.unit.solid_inlet.mass_frac_comp[0, "Fe2O3"].fix(0.244162011502)
        m.fs.unit.solid_inlet.mass_frac_comp[0, "Fe3O4"].fix(0.201998299487)
        m.fs.unit.solid_inlet.mass_frac_comp[0, "Al2O3"].fix(0.553839689011)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, iron_oc):
        assert hasattr(iron_oc.fs.unit, "gas_inlet")
        assert len(iron_oc.fs.unit.gas_inlet.vars) == 4
        assert isinstance(iron_oc.fs.unit.gas_inlet.flow_mol, Var)
        assert isinstance(iron_oc.fs.unit.gas_inlet.mole_frac_comp, Var)
        assert isinstance(iron_oc.fs.unit.gas_inlet.temperature, Var)
        assert isinstance(iron_oc.fs.unit.gas_inlet.pressure, Var)

        assert hasattr(iron_oc.fs.unit, "solid_inlet")
        assert len(iron_oc.fs.unit.solid_inlet.vars) == 4
        assert isinstance(iron_oc.fs.unit.solid_inlet.flow_mass, Var)
        assert isinstance(iron_oc.fs.unit.solid_inlet.particle_porosity, Var)
        assert isinstance(iron_oc.fs.unit.solid_inlet.mass_frac_comp, Var)
        assert isinstance(iron_oc.fs.unit.solid_inlet.temperature, Var)

        assert hasattr(iron_oc.fs.unit, "gas_outlet")
        assert len(iron_oc.fs.unit.gas_outlet.vars) == 4
        assert isinstance(iron_oc.fs.unit.gas_outlet.flow_mol, Var)
        assert isinstance(iron_oc.fs.unit.gas_outlet.mole_frac_comp, Var)
        assert isinstance(iron_oc.fs.unit.gas_outlet.temperature, Var)
        assert isinstance(iron_oc.fs.unit.gas_outlet.pressure, Var)

        assert hasattr(iron_oc.fs.unit, "solid_outlet")
        assert len(iron_oc.fs.unit.solid_outlet.vars) == 4
        assert isinstance(iron_oc.fs.unit.solid_outlet.flow_mass, Var)
        assert isinstance(iron_oc.fs.unit.solid_outlet.particle_porosity, Var)
        assert isinstance(iron_oc.fs.unit.solid_outlet.mass_frac_comp, Var)
        assert isinstance(iron_oc.fs.unit.solid_outlet.temperature, Var)

        assert isinstance(iron_oc.fs.unit.orifice_area, Constraint)
        assert isinstance(iron_oc.fs.unit.bed_area_eqn, Constraint)
        assert isinstance(iron_oc.fs.unit.bubble_area, Constraint)
        assert isinstance(iron_oc.fs.unit.gas_emulsion_area, Constraint)
        assert isinstance(iron_oc.fs.unit.solid_emulsion_area, Constraint)
        assert isinstance(iron_oc.fs.unit.bubble_cloud_heat_trans_coeff, Constraint)
        assert isinstance(iron_oc.fs.unit.convective_heat_trans_coeff, Constraint)
        assert isinstance(iron_oc.fs.unit.convective_heat_transfer, Constraint)
        assert isinstance(iron_oc.fs.unit.bubble_cloud_bulk_heat_trans, Constraint)
        assert isinstance(iron_oc.fs.unit.bubble_mass_transfer, Constraint)
        assert isinstance(iron_oc.fs.unit.gas_emulsion_mass_transfer, Constraint)
        assert isinstance(iron_oc.fs.unit.bubble_heat_transfer, Constraint)
        assert isinstance(iron_oc.fs.unit.gas_emulsion_heat_transfer, Constraint)
        assert isinstance(iron_oc.fs.unit.solid_emulsion_heat_transfer, Constraint)

        assert number_variables(iron_oc) == 915
        assert number_total_constraints(iron_oc) == 867
        assert number_unused_variables(iron_oc) == 21

    @pytest.mark.unit
    def test_dof(self, iron_oc):
        assert degrees_of_freedom(iron_oc) == 0

    @pytest.fixture(scope="class")
    def iron_oc_unscaled(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        # Set up thermo props and reaction props
        m.fs.gas_properties = GasPhaseParameterBlock()
        m.fs.solid_properties = SolidPhaseParameterBlock()
        m.fs.hetero_reactions = HeteroReactionParameterBlock(
            solid_property_package=m.fs.solid_properties,
            gas_property_package=m.fs.gas_properties,
        )

        m.fs.unit = BubblingFluidizedBed(
            flow_type="co_current",
            finite_elements=5,
            transformation_method="dae.finite_difference",
            gas_phase_config={"property_package": m.fs.gas_properties},
            solid_phase_config={
                "property_package": m.fs.solid_properties,
                "reaction_package": m.fs.hetero_reactions,
            },
        )

        # Fix geometry variables
        m.fs.unit.number_orifice.fix(2500)  # [-]
        m.fs.unit.bed_diameter.fix(6.5)  # m
        m.fs.unit.bed_height.fix(5)  # m

        # Fix inlet port variables for gas and solid
        m.fs.unit.gas_inlet.flow_mol[0].fix(1567.79)  # mol/s
        m.fs.unit.gas_inlet.temperature[0].fix(373)  # K
        m.fs.unit.gas_inlet.pressure[0].fix(1.86e5)  # Pa = 1E5 bar
        m.fs.unit.gas_inlet.mole_frac_comp[0, "O2"].fix(0.2095)
        m.fs.unit.gas_inlet.mole_frac_comp[0, "N2"].fix(0.7808)
        m.fs.unit.gas_inlet.mole_frac_comp[0, "CO2"].fix(0.0004)
        m.fs.unit.gas_inlet.mole_frac_comp[0, "H2O"].fix(0.0093)

        m.fs.unit.solid_inlet.flow_mass[0].fix(1230.865)  # kg/s
        m.fs.unit.solid_inlet.particle_porosity[0].fix(0.27)  # (-)
        m.fs.unit.solid_inlet.temperature[0].fix(1173.9)  # K
        m.fs.unit.solid_inlet.mass_frac_comp[0, "Fe2O3"].fix(0.244162011502)
        m.fs.unit.solid_inlet.mass_frac_comp[0, "Fe3O4"].fix(0.201998299487)
        m.fs.unit.solid_inlet.mass_frac_comp[0, "Al2O3"].fix(0.553839689011)

        return m

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize_unscaled(self, iron_oc_unscaled):
        initialization_tester(
            iron_oc_unscaled,
            optarg={"tol": 1e-6},
            gas_phase_state_args={
                "flow_mol": 1567.79,
                "temperature": 1173.9,
                "pressure": 1.86e5,
            },
            solid_phase_state_args={"flow_mass": 1230.865, "temperature": 1173.9},
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve_unscaled(self, iron_oc_unscaled):
        results = solver.solve(iron_oc_unscaled)

        # Check for optimal solution
        assert check_optimal_termination(results)
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.component
    def test_scaling(self, iron_oc):
        BFB = iron_oc.fs.unit  # alias to keep test lines short
        iscale.calculate_scaling_factors(BFB)

        assert pytest.approx(0.01538, abs=1e-5) == iscale.get_scaling_factor(
            BFB.bed_diameter
        )
        assert pytest.approx(0.00301, abs=1e-5) == iscale.get_scaling_factor(
            BFB.bed_area
        )
        for (t, x), v in BFB.delta.items():
            assert pytest.approx(1.00000, abs=1e-5) == iscale.get_scaling_factor(v)
        for (t, x), v in BFB.bubble_diameter_max.items():
            assert pytest.approx(0.15385, abs=1e-5) == iscale.get_scaling_factor(v)
        for (t, x, j), v in BFB.Kgbulk_c.items():
            assert pytest.approx(3.014e-5, abs=1e-7) == iscale.get_scaling_factor(v)
        for (t, x), v in BFB._reform_var_1.items():
            assert pytest.approx(65.00000, abs=1e-5) == iscale.get_scaling_factor(v)
        for (t, x, j), v in BFB._reform_var_2.items():
            assert pytest.approx(53.53089, abs=1e-5) == iscale.get_scaling_factor(v)
        for (t, x), v in BFB._reform_var_3.items():
            assert pytest.approx(2.83941, abs=1e-5) == iscale.get_scaling_factor(v)
        for (t, x), v in BFB.Hbe.items():
            assert pytest.approx(1.000e-9, abs=1e-11) == iscale.get_scaling_factor(v)
        for (t, x), v in BFB.Hgbulk.items():
            assert pytest.approx(1.000e-9, abs=1e-8) == iscale.get_scaling_factor(v)
        for (t, x), v in BFB._reform_var_4.items():
            assert pytest.approx(1.093e-9, abs=1e-8) == iscale.get_scaling_factor(v)
        for (t, x), v in BFB._reform_var_5.items():
            assert pytest.approx(6.667e-6, abs=1e-8) == iscale.get_scaling_factor(v)
        for (t, x), v in BFB.htc_conv.items():
            assert pytest.approx(6.667e-6, abs=1e-5) == iscale.get_scaling_factor(v)
        for (t, x), v in BFB.ht_conv.items():
            assert pytest.approx(6.667e-6, abs=1e-5) == iscale.get_scaling_factor(v)
        for (t, x), v in BFB.bubble.heat.items():
            assert pytest.approx(1e-9, abs=1e-11) == iscale.get_scaling_factor(v)
        for (t, x), v in BFB.gas_emulsion.heat.items():
            assert pytest.approx(1e-9, abs=1e-11) == iscale.get_scaling_factor(v)
        for (t, x), v in BFB.solid_emulsion.heat.items():
            assert pytest.approx(1e-9, abs=1e-11) == iscale.get_scaling_factor(v)

        for c in BFB.orifice_area.values():
            assert pytest.approx(
                0.00301, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for c in BFB.bed_area_eqn.values():
            assert pytest.approx(
                0.00301, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB.bubble_area.items():
            assert pytest.approx(
                0.00301, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB.gas_emulsion_area.items():
            assert pytest.approx(
                0.00301, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB.solid_emulsion_area.items():
            assert pytest.approx(
                0.00301, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB.emulsion_vol_frac.items():
            assert pytest.approx(
                1.0000, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB.average_voidage.items():
            assert pytest.approx(
                1.0000, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB.bubble_growth_coefficient.items():
            assert pytest.approx(
                0.01538, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB.bubble_diameter_maximum.items():
            assert pytest.approx(
                8.619e-10, abs=1e-12
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB._reformulation_eqn_1.items():
            assert pytest.approx(
                6.50000, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB.bubble_diameter_eqn.items():
            assert pytest.approx(
                0.01538, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB.ddia_bubbledx_disc_eq.items():
            assert pytest.approx(
                0.01538, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB.bubble_velocity_rise.items():
            assert pytest.approx(
                0.15385, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB.emulsion_voidage.items():
            assert pytest.approx(
                1.00000, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB.bubble_velocity.items():
            assert pytest.approx(
                0.00100, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB.average_gas_density_eqn.items():
            assert pytest.approx(
                1e-6, abs=1e-8
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB.velocity_gas_superficial.items():
            assert pytest.approx(
                0.00100, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB.bubble_vol_frac_eqn.items():
            assert pytest.approx(
                0.00100, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB.solid_super_vel.items():
            assert pytest.approx(
                0.00001, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB.gas_emulsion_pressure_drop.items():
            assert pytest.approx(
                1e-5, abs=1e-7
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x, j), c in BFB._reformulation_eqn_2.items():
            assert pytest.approx(
                0.53531, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB._reformulation_eqn_3.items():
            assert pytest.approx(
                2.83941, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x, j), c in BFB.bubble_cloud_mass_trans_coeff.items():
            assert pytest.approx(
                2.83941, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x, j), c in BFB.bubble_cloud_bulk_mass_trans.items():
            assert pytest.approx(
                3.014e-5, abs=1e-7
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x, j), c in BFB.bubble_mass_transfer.items():
            assert pytest.approx(
                3.014e-5, abs=1e-7
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x, j), c in BFB.gas_emulsion_mass_transfer.items():
            assert pytest.approx(
                3.014e-5, abs=1e-7
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x, r), c in BFB.solid_emulsion_rxn_ext_constraint.items():
            assert pytest.approx(
                30.13585, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x, j), c in BFB.gas_emulsion_hetero_rxn_eqn.items():
            assert pytest.approx(
                30.13585, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB.bubble_gas_flowrate.items():
            assert pytest.approx(
                0.00100, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB.emulsion_gas_flowrate.items():
            assert pytest.approx(
                0.00100, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for t, c in BFB.gas_emulsion_pressure_in.items():
            assert pytest.approx(
                1e-5, abs=1e-7
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for t, c in BFB.gas_mole_flow_in.items():
            assert pytest.approx(
                0.00100, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for t, c in BFB.particle_porosity_in.items():
            assert pytest.approx(
                100, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for t, c in BFB.emulsion_gas_velocity_in.items():
            assert pytest.approx(
                25.23722, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, j), c in BFB.bubble_mole_frac_in.items():
            assert pytest.approx(
                10.00000, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, j), c in BFB.gas_emulsion_mole_frac_in.items():
            assert pytest.approx(
                10.00000, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for t, c in BFB.solid_emulsion_mass_flow_in.items():
            assert pytest.approx(
                0.00100, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, j), c in BFB.solid_emulsion_mass_frac_in.items():
            assert pytest.approx(
                10.00000, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for t, c in BFB.gas_pressure_out.items():
            assert pytest.approx(
                1e-5, abs=1e-7
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, p, j), c in BFB.gas_material_balance_out.items():
            assert pytest.approx(
                0.00100, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, p, j), c in BFB.solid_material_balance_out.items():
            assert pytest.approx(
                0.00100, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for t, c in BFB.solid_particle_porosity_out.items():
            assert pytest.approx(
                100.00000, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB._reformulation_eqn_4.items():
            assert pytest.approx(
                0, abs=1e-10
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB._reformulation_eqn_5.items():
            assert pytest.approx(
                0.66666, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB.bubble_cloud_heat_trans_coeff.items():
            assert pytest.approx(
                1e-11, abs=1e-11
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB.convective_heat_trans_coeff.items():
            assert pytest.approx(
                6.667e-6, abs=1e-8
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB.convective_heat_transfer.items():
            assert pytest.approx(
                6.667e-8, abs=1e-10
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB.bubble_cloud_bulk_heat_trans.items():
            assert pytest.approx(
                1e-9, abs=1e-11
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB.bubble_heat_transfer.items():
            assert pytest.approx(
                1e-11, abs=1e-11
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB.gas_emulsion_heat_transfer.items():
            if x == 0:
                assert pytest.approx(
                    1e-11, abs=1e-11
                ) == iscale.get_constraint_transform_applied_scaling_factor(c)
            else:
                assert pytest.approx(
                    1e-9, abs=1e-11
                ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB.solid_emulsion_heat_transfer.items():
            if x == 0:
                assert pytest.approx(
                    1e-11, abs=1e-11
                ) == iscale.get_constraint_transform_applied_scaling_factor(c)
            else:
                assert pytest.approx(
                    1e-9, abs=1e-11
                ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, p), c in BFB.gas_energy_balance_in.items():
            assert pytest.approx(
                1e-9, abs=1e-11
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for t, c in BFB.gas_emulsion_temperature_in.items():
            assert pytest.approx(
                0.01000, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for t, c in BFB.solid_energy_balance_in.items():
            assert pytest.approx(
                1e-9, abs=1e-11
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for t, c in BFB.gas_energy_balance_out.items():
            assert pytest.approx(
                1e-9, abs=1e-11
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for t, c in BFB.solid_energy_balance_out.items():
            assert pytest.approx(
                1e-9, abs=1e-11
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, iron_oc):
        BFB = iron_oc.fs.unit  # alias to keep test lines short
        iscale.calculate_scaling_factors(BFB)
        initialization_tester(
            iron_oc,
            optarg={"tol": 1e-6},
            gas_phase_state_args={
                "flow_mol": 1567.79,
                "temperature": 1173.9,
                "pressure": 1.86e5,
            },
            solid_phase_state_args={"flow_mass": 1230.865, "temperature": 1173.9},
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, iron_oc):
        results = solver.solve(iron_oc, tee=True)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, iron_oc):
        assert (
            pytest.approx(0.80244, abs=1e-5)
            == iron_oc.fs.unit.velocity_superficial_gas[0, 0].value
        )
        assert (
            pytest.approx(3.08621, abs=1e-5)
            == iron_oc.fs.unit.velocity_superficial_gas[0, 1].value
        )
        assert (
            pytest.approx(0.03431, abs=1e-5)
            == iron_oc.fs.unit.bubble_diameter[0, 0].value
        )
        assert (
            pytest.approx(1.76158, abs=1e-5)
            == iron_oc.fs.unit.bubble_diameter[0, 1].value
        )
        assert (
            pytest.approx(1.17521, abs=1e-5)
            == iron_oc.fs.unit.velocity_bubble[0, 0].value
        )
        assert (
            pytest.approx(6.00175, abs=1e-5)
            == iron_oc.fs.unit.velocity_bubble[0, 1].value
        )
        assert pytest.approx(0.64909, abs=1e-5) == iron_oc.fs.unit.delta[0, 0].value
        assert pytest.approx(0.50762, abs=1e-5) == iron_oc.fs.unit.delta[0, 1].value
        assert (
            pytest.approx(143561.507565, abs=1e-5)
            == iron_oc.fs.unit.gas_outlet.pressure[0].value
        )
        assert (
            pytest.approx(42438.49244, abs=1e-5)
            == iron_oc.fs.unit.gas_inlet.pressure[0].value
            - iron_oc.fs.unit.gas_outlet.pressure[0].value
        )

    @pytest.mark.component
    def test_units_consistent(self, iron_oc):
        assert_units_consistent(iron_oc)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, iron_oc):
        # Conservation of material check
        calculate_variable_from_constraint(
            iron_oc.fs.unit.gas_inlet_block[0].mw,
            iron_oc.fs.unit.gas_inlet_block[0].mw_eqn,
        )
        calculate_variable_from_constraint(
            iron_oc.fs.unit.gas_outlet_block[0].mw,
            iron_oc.fs.unit.gas_outlet_block[0].mw_eqn,
        )
        mbal_gas = value(
            (
                iron_oc.fs.unit.gas_inlet.flow_mol[0]
                * iron_oc.fs.unit.gas_inlet_block[0].mw
            )
            - (
                iron_oc.fs.unit.gas_outlet.flow_mol[0]
                * iron_oc.fs.unit.gas_outlet_block[0].mw
            )
        )
        mbal_solid = value(
            iron_oc.fs.unit.solid_inlet.flow_mass[0]
            - iron_oc.fs.unit.solid_outlet.flow_mass[0]
        )
        mbal_tol = mbal_gas + mbal_solid
        assert abs(mbal_tol) <= 1e-2

        # Reaction stoichiometric ratio check
        # Overall reducer reactions for iron oxidation:
        # O2 + 4Fe3O4 => 6Fe2O3
        mole_gas_reacted = value(
            iron_oc.fs.unit.gas_inlet.flow_mol[0]
            * iron_oc.fs.unit.gas_inlet.mole_frac_comp[0, "O2"]
            - iron_oc.fs.unit.gas_outlet.flow_mol[0]
            * iron_oc.fs.unit.gas_outlet.mole_frac_comp[0, "O2"]
        )
        mole_solid_reacted = value(
            (
                iron_oc.fs.unit.solid_inlet.flow_mass[0]
                * iron_oc.fs.unit.solid_inlet.mass_frac_comp[0, "Fe3O4"]
                / iron_oc.fs.unit.solid_inlet_block[0]._params.mw_comp["Fe3O4"]
            )
            - (
                iron_oc.fs.unit.solid_outlet.flow_mass[0]
                * iron_oc.fs.unit.solid_outlet.mass_frac_comp[0, "Fe3O4"]
                / iron_oc.fs.unit.solid_outlet_block[0]._params.mw_comp["Fe3O4"]
            )
        )
        stoichiometric_ratio = mole_solid_reacted / mole_gas_reacted
        assert pytest.approx(4, abs=1e-6) == stoichiometric_ratio

        # Conservation of energy check
        ebal_gas = value(
            (
                iron_oc.fs.unit.gas_inlet.flow_mol[0]
                * iron_oc.fs.unit.gas_inlet_block[0].enth_mol
            )
            - (
                iron_oc.fs.unit.gas_outlet.flow_mol[0]
                * iron_oc.fs.unit.gas_outlet_block[0].enth_mol
            )
        )
        ebal_solid = value(
            (
                iron_oc.fs.unit.solid_inlet.flow_mass[0]
                * iron_oc.fs.unit.solid_inlet_block[0].enth_mass
            )
            - (
                iron_oc.fs.unit.solid_outlet.flow_mass[0]
                * iron_oc.fs.unit.solid_outlet_block[0].enth_mass
            )
        )
        e_reaction = value(
            mole_gas_reacted
            * iron_oc.fs.unit.solid_emulsion.reactions[0, 0]._params.dh_rxn["R1"]
        )
        ebal_tol = ebal_gas + ebal_solid - e_reaction
        assert abs(ebal_tol) <= 1e-2

    @pytest.mark.ui
    @pytest.mark.unit
    def test_report(self, iron_oc):
        iron_oc.fs.unit.report()

    @pytest.mark.component
    def test_initialization_error(self, iron_oc):
        iron_oc.fs.unit.gas_outlet.flow_mol[0].fix(1)

        with pytest.raises(InitializationError):
            iron_oc.fs.unit.initialize()


# -----------------------------------------------------------------------------
class TestIronOC_EnergyBalanceType(object):
    @pytest.fixture(scope="class")
    def iron_oc(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        # Set up thermo props and reaction props
        m.fs.gas_properties = GasPhaseParameterBlock()
        m.fs.solid_properties = SolidPhaseParameterBlock()
        m.fs.hetero_reactions = HeteroReactionParameterBlock(
            solid_property_package=m.fs.solid_properties,
            gas_property_package=m.fs.gas_properties,
        )

        m.fs.unit = BubblingFluidizedBed(
            energy_balance_type=EnergyBalanceType.none,
            flow_type="co_current",
            finite_elements=5,
            transformation_method="dae.finite_difference",
            gas_phase_config={"property_package": m.fs.gas_properties},
            solid_phase_config={
                "property_package": m.fs.solid_properties,
                "reaction_package": m.fs.hetero_reactions,
            },
        )

        # Fix geometry variables
        m.fs.unit.number_orifice.fix(2500)  # [-]
        m.fs.unit.bed_diameter.fix(6.5)  # m
        m.fs.unit.bed_height.fix(5)  # m

        # Fix inlet port variables for gas and solid
        m.fs.unit.gas_inlet.flow_mol[0].fix(1567.79)  # mol/s
        m.fs.unit.gas_inlet.temperature[0].fix(373)  # K
        m.fs.unit.gas_inlet.pressure[0].fix(1.86e5)  # Pa = 1E5 bar
        m.fs.unit.gas_inlet.mole_frac_comp[0, "O2"].fix(0.2095)
        m.fs.unit.gas_inlet.mole_frac_comp[0, "N2"].fix(0.7808)
        m.fs.unit.gas_inlet.mole_frac_comp[0, "CO2"].fix(0.0004)
        m.fs.unit.gas_inlet.mole_frac_comp[0, "H2O"].fix(0.0093)

        m.fs.unit.solid_inlet.flow_mass[0].fix(1230.865)  # kg/s
        m.fs.unit.solid_inlet.particle_porosity[0].fix(0.27)  # (-)
        m.fs.unit.solid_inlet.temperature[0].fix(1173.9)  # K
        m.fs.unit.solid_inlet.mass_frac_comp[0, "Fe2O3"].fix(0.244162011502)
        m.fs.unit.solid_inlet.mass_frac_comp[0, "Fe3O4"].fix(0.201998299487)
        m.fs.unit.solid_inlet.mass_frac_comp[0, "Al2O3"].fix(0.553839689011)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, iron_oc):
        assert hasattr(iron_oc.fs.unit, "gas_inlet")
        assert len(iron_oc.fs.unit.gas_inlet.vars) == 4
        assert isinstance(iron_oc.fs.unit.gas_inlet.flow_mol, Var)
        assert isinstance(iron_oc.fs.unit.gas_inlet.mole_frac_comp, Var)
        assert isinstance(iron_oc.fs.unit.gas_inlet.temperature, Var)
        assert isinstance(iron_oc.fs.unit.gas_inlet.pressure, Var)

        assert hasattr(iron_oc.fs.unit, "solid_inlet")
        assert len(iron_oc.fs.unit.solid_inlet.vars) == 4
        assert isinstance(iron_oc.fs.unit.solid_inlet.flow_mass, Var)
        assert isinstance(iron_oc.fs.unit.solid_inlet.particle_porosity, Var)
        assert isinstance(iron_oc.fs.unit.solid_inlet.mass_frac_comp, Var)
        assert isinstance(iron_oc.fs.unit.solid_inlet.temperature, Var)

        assert hasattr(iron_oc.fs.unit, "gas_outlet")
        assert len(iron_oc.fs.unit.gas_outlet.vars) == 4
        assert isinstance(iron_oc.fs.unit.gas_outlet.flow_mol, Var)
        assert isinstance(iron_oc.fs.unit.gas_outlet.mole_frac_comp, Var)
        assert isinstance(iron_oc.fs.unit.gas_outlet.temperature, Var)
        assert isinstance(iron_oc.fs.unit.gas_outlet.pressure, Var)

        assert hasattr(iron_oc.fs.unit, "solid_outlet")
        assert len(iron_oc.fs.unit.solid_outlet.vars) == 4
        assert isinstance(iron_oc.fs.unit.solid_outlet.flow_mass, Var)
        assert isinstance(iron_oc.fs.unit.solid_outlet.particle_porosity, Var)
        assert isinstance(iron_oc.fs.unit.solid_outlet.mass_frac_comp, Var)
        assert isinstance(iron_oc.fs.unit.solid_outlet.temperature, Var)

        assert isinstance(iron_oc.fs.unit.gas_energy_balance_out, Constraint)
        assert isinstance(iron_oc.fs.unit.solid_energy_balance_out, Constraint)
        assert isinstance(iron_oc.fs.unit.isothermal_gas_emulsion, Constraint)
        assert isinstance(iron_oc.fs.unit.isothermal_solid_emulsion, Constraint)
        assert isinstance(iron_oc.fs.unit.isothermal_bubble, Constraint)

        assert number_variables(iron_oc) == 741
        assert number_total_constraints(iron_oc) == 660
        assert number_unused_variables(iron_oc) == 55
        print(unused_variables_set(iron_oc))

    @pytest.mark.unit
    def test_dof(self, iron_oc):
        assert degrees_of_freedom(iron_oc) == 0

    @pytest.fixture(scope="class")
    def iron_oc_unscaled(self, iron_oc):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        # Set up thermo props and reaction props
        m.fs.gas_properties = GasPhaseParameterBlock()
        m.fs.solid_properties = SolidPhaseParameterBlock()
        m.fs.hetero_reactions = HeteroReactionParameterBlock(
            solid_property_package=m.fs.solid_properties,
            gas_property_package=m.fs.gas_properties,
        )

        m.fs.unit = BubblingFluidizedBed(
            energy_balance_type=EnergyBalanceType.none,
            flow_type="co_current",
            finite_elements=5,
            transformation_method="dae.finite_difference",
            gas_phase_config={"property_package": m.fs.gas_properties},
            solid_phase_config={
                "property_package": m.fs.solid_properties,
                "reaction_package": m.fs.hetero_reactions,
            },
        )

        # Fix geometry variables
        m.fs.unit.number_orifice.fix(2500)  # [-]
        m.fs.unit.bed_diameter.fix(6.5)  # m
        m.fs.unit.bed_height.fix(5)  # m

        # Fix inlet port variables for gas and solid
        m.fs.unit.gas_inlet.flow_mol[0].fix(1567.79)  # mol/s
        m.fs.unit.gas_inlet.temperature[0].fix(373)  # K
        m.fs.unit.gas_inlet.pressure[0].fix(1.86e5)  # Pa = 1E5 bar
        m.fs.unit.gas_inlet.mole_frac_comp[0, "O2"].fix(0.2095)
        m.fs.unit.gas_inlet.mole_frac_comp[0, "N2"].fix(0.7808)
        m.fs.unit.gas_inlet.mole_frac_comp[0, "CO2"].fix(0.0004)
        m.fs.unit.gas_inlet.mole_frac_comp[0, "H2O"].fix(0.0093)

        m.fs.unit.solid_inlet.flow_mass[0].fix(1230.865)  # kg/s
        m.fs.unit.solid_inlet.particle_porosity[0].fix(0.27)  # (-)
        m.fs.unit.solid_inlet.temperature[0].fix(1173.9)  # K
        m.fs.unit.solid_inlet.mass_frac_comp[0, "Fe2O3"].fix(0.244162011502)
        m.fs.unit.solid_inlet.mass_frac_comp[0, "Fe3O4"].fix(0.201998299487)
        m.fs.unit.solid_inlet.mass_frac_comp[0, "Al2O3"].fix(0.553839689011)

        return m

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize_unscaled(self, iron_oc_unscaled):
        # needs a little more help converging in under 100 iterations (default)
        # so we provide a better set of state argument flags for initialization
        # alternatively, could raise maximum iterations to 105 in optarg
        initialization_tester(
            iron_oc_unscaled,
            optarg={"tol": 1e-6},
            gas_phase_state_args={
                "flow_mol": 1567.79,
                "temperature": 1173.9,
                "pressure": 1.86e5,
                "mole_frac_comp": {
                    "O2": 0.2095,
                    "N2": 0.7808,
                    "CO2": 0.0004,
                    "H2O": 0.0093,
                },
            },
            solid_phase_state_args={
                "flow_mass": 1230.865,
                "temperature": 1173.9,
                "mole_frac_comp": {
                    "Fe2O3": 0.244162011502,
                    "Fe3O4": 0.201998299487,
                    "Al2O3": 0.553839689011,
                },
            },
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve_unscaled(self, iron_oc_unscaled):
        results = solver.solve(iron_oc_unscaled)

        # Check for optimal solution
        assert check_optimal_termination(results)
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.component
    def test_scaling(self, iron_oc):
        BFB = iron_oc.fs.unit  # alias to keep test lines short
        iscale.calculate_scaling_factors(BFB)

        assert pytest.approx(0.01538, abs=1e-5) == iscale.get_scaling_factor(
            BFB.bed_diameter
        )
        assert pytest.approx(0.00301, abs=1e-5) == iscale.get_scaling_factor(
            BFB.bed_area
        )
        for (t, x), v in BFB.delta.items():
            assert pytest.approx(1.00000, abs=1e-5) == iscale.get_scaling_factor(v)
        for (t, x), v in BFB.bubble_diameter_max.items():
            assert pytest.approx(0.15385, abs=1e-5) == iscale.get_scaling_factor(v)
        for (t, x, j), v in BFB.Kgbulk_c.items():
            assert pytest.approx(3.014e-5, abs=1e-5) == iscale.get_scaling_factor(v)
        for (t, x), v in BFB._reform_var_1.items():
            assert pytest.approx(65.00000, abs=1e-5) == iscale.get_scaling_factor(v)
        for (t, x, j), v in BFB._reform_var_2.items():
            assert pytest.approx(53.53089, abs=1e-5) == iscale.get_scaling_factor(v)
        for (t, x), v in BFB._reform_var_3.items():
            assert pytest.approx(2.83941, abs=1e-5) == iscale.get_scaling_factor(v)

        for c in BFB.orifice_area.values():
            assert pytest.approx(
                3.014e-3, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for c in BFB.bed_area_eqn.values():
            assert pytest.approx(
                3.014e-3, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB.bubble_area.items():
            assert pytest.approx(
                3.014e-3, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB.gas_emulsion_area.items():
            assert pytest.approx(
                3.014e-3, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB.solid_emulsion_area.items():
            assert pytest.approx(
                3.014e-3, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB.emulsion_vol_frac.items():
            assert pytest.approx(
                1.00000, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB.average_voidage.items():
            assert pytest.approx(
                1.00000, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB.bubble_growth_coefficient.items():
            assert pytest.approx(
                0.015385, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB.bubble_diameter_maximum.items():
            assert pytest.approx(
                8.619e-10, abs=1e-12
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB._reformulation_eqn_1.items():
            assert pytest.approx(
                6.50000, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB.bubble_diameter_eqn.items():
            assert pytest.approx(
                0.01538, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB.ddia_bubbledx_disc_eq.items():
            assert pytest.approx(
                0.01538, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB.bubble_velocity_rise.items():
            assert pytest.approx(
                0.15385, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB.emulsion_voidage.items():
            assert pytest.approx(
                1.00000, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB.bubble_velocity.items():
            assert pytest.approx(
                0.00100, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB.average_gas_density_eqn.items():
            assert pytest.approx(
                1e-6, abs=1e-8
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB.velocity_gas_superficial.items():
            assert pytest.approx(
                0.00100, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB.bubble_vol_frac_eqn.items():
            assert pytest.approx(
                0.00100, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB.solid_super_vel.items():
            assert pytest.approx(
                0.00001, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB.gas_emulsion_pressure_drop.items():
            assert pytest.approx(
                1e-5, abs=1e-7
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x, j), c in BFB._reformulation_eqn_2.items():
            assert pytest.approx(
                0.53531, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB._reformulation_eqn_3.items():
            assert pytest.approx(
                2.83941, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x, j), c in BFB.bubble_cloud_mass_trans_coeff.items():
            assert pytest.approx(
                2.83941, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x, j), c in BFB.bubble_cloud_bulk_mass_trans.items():
            assert pytest.approx(
                3.014e-5, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x, j), c in BFB.bubble_mass_transfer.items():
            assert pytest.approx(
                3.014e-5, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x, j), c in BFB.gas_emulsion_mass_transfer.items():
            assert pytest.approx(
                3.014e-5, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x, r), c in BFB.solid_emulsion_rxn_ext_constraint.items():
            assert pytest.approx(
                30.13585, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x, j), c in BFB.gas_emulsion_hetero_rxn_eqn.items():
            assert pytest.approx(
                30.13585, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB.bubble_gas_flowrate.items():
            assert pytest.approx(
                0.00100, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in BFB.emulsion_gas_flowrate.items():
            assert pytest.approx(
                0.00100, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for t, c in BFB.gas_emulsion_pressure_in.items():
            assert pytest.approx(
                1e-5, abs=1e-7
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for t, c in BFB.gas_mole_flow_in.items():
            assert pytest.approx(
                0.00100, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for t, c in BFB.particle_porosity_in.items():
            assert pytest.approx(
                100.00000, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for t, c in BFB.emulsion_gas_velocity_in.items():
            assert pytest.approx(
                25.23722, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, j), c in BFB.bubble_mole_frac_in.items():
            assert pytest.approx(
                10.00000, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, j), c in BFB.gas_emulsion_mole_frac_in.items():
            assert pytest.approx(
                10.00000, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for t, c in BFB.solid_emulsion_mass_flow_in.items():
            assert pytest.approx(
                0.00100, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, j), c in BFB.solid_emulsion_mass_frac_in.items():
            assert pytest.approx(
                10.00000, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for t, c in BFB.gas_pressure_out.items():
            assert pytest.approx(
                1e-5, abs=1e-7
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, p, j), c in BFB.gas_material_balance_out.items():
            assert pytest.approx(
                0.00100, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, p, j), c in BFB.solid_material_balance_out.items():
            assert pytest.approx(
                0.00100, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for t, c in BFB.solid_particle_porosity_out.items():
            assert pytest.approx(
                100.00000, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for t, c in BFB.gas_energy_balance_out.items():
            assert pytest.approx(
                0.01, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for t, c in BFB.solid_energy_balance_out.items():
            assert pytest.approx(
                0.01, abs=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, iron_oc):
        BFB = iron_oc.fs.unit  # alias to keep test lines short
        iscale.calculate_scaling_factors(BFB)
        initialization_tester(
            iron_oc,
            optarg={"tol": 1e-6},
            gas_phase_state_args={
                "flow_mol": 1567.79,
                "temperature": 1173.9,
                "pressure": 1.86e5,
            },
            solid_phase_state_args={"flow_mass": 1230.865, "temperature": 1173.9},
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, iron_oc):
        results = solver.solve(iron_oc, tee=True)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, iron_oc):
        assert (
            pytest.approx(0.80244, abs=1e-5)
            == iron_oc.fs.unit.velocity_superficial_gas[0, 0].value
        )
        assert (
            pytest.approx(1.10353, abs=1e-5)
            == iron_oc.fs.unit.velocity_superficial_gas[0, 1].value
        )
        assert (
            pytest.approx(0.03431, abs=1e-5)
            == iron_oc.fs.unit.bubble_diameter[0, 0].value
        )
        assert (
            pytest.approx(1.09929, abs=1e-5)
            == iron_oc.fs.unit.bubble_diameter[0, 1].value
        )
        assert (
            pytest.approx(1.17521, abs=1e-5)
            == iron_oc.fs.unit.velocity_bubble[0, 0].value
        )
        assert (
            pytest.approx(3.39836, abs=1e-5)
            == iron_oc.fs.unit.velocity_bubble[0, 1].value
        )
        assert pytest.approx(0.64909, abs=1e-5) == iron_oc.fs.unit.delta[0, 0].value
        assert pytest.approx(0.31306, abs=1e-5) == iron_oc.fs.unit.delta[0, 1].value
        assert (
            pytest.approx(124978.80283, abs=1e-5)
            == iron_oc.fs.unit.gas_outlet.pressure[0].value
        )
        assert (
            pytest.approx(61021.19717, abs=1e-5)
            == iron_oc.fs.unit.gas_inlet.pressure[0].value
            - iron_oc.fs.unit.gas_outlet.pressure[0].value
        )

    @pytest.mark.component
    def test_units_consistent(self, iron_oc):
        assert_units_consistent(iron_oc)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, iron_oc):
        # Conservation of material check
        calculate_variable_from_constraint(
            iron_oc.fs.unit.gas_inlet_block[0].mw,
            iron_oc.fs.unit.gas_inlet_block[0].mw_eqn,
        )
        calculate_variable_from_constraint(
            iron_oc.fs.unit.gas_outlet_block[0].mw,
            iron_oc.fs.unit.gas_outlet_block[0].mw_eqn,
        )
        mbal_gas = value(
            (
                iron_oc.fs.unit.gas_inlet.flow_mol[0]
                * iron_oc.fs.unit.gas_inlet_block[0].mw
            )
            - (
                iron_oc.fs.unit.gas_outlet.flow_mol[0]
                * iron_oc.fs.unit.gas_outlet_block[0].mw
            )
        )
        mbal_solid = value(
            iron_oc.fs.unit.solid_inlet.flow_mass[0]
            - iron_oc.fs.unit.solid_outlet.flow_mass[0]
        )
        mbal_tol = mbal_gas + mbal_solid
        assert abs(mbal_tol) <= 1e-2

        # Reaction stoichiometric ratio check
        # Overall reducer reactions for methane combustion:
        # CH4 + 12Fe2O3 => 8Fe3O4 + CO2 + 2H2O
        mole_gas_reacted = value(
            iron_oc.fs.unit.gas_inlet.flow_mol[0]
            * iron_oc.fs.unit.gas_inlet.mole_frac_comp[0, "O2"]
            - iron_oc.fs.unit.gas_outlet.flow_mol[0]
            * iron_oc.fs.unit.gas_outlet.mole_frac_comp[0, "O2"]
        )
        mole_solid_reacted = value(
            (
                iron_oc.fs.unit.solid_inlet.flow_mass[0]
                * iron_oc.fs.unit.solid_inlet.mass_frac_comp[0, "Fe3O4"]
                / iron_oc.fs.unit.solid_inlet_block[0]._params.mw_comp["Fe3O4"]
            )
            - (
                iron_oc.fs.unit.solid_outlet.flow_mass[0]
                * iron_oc.fs.unit.solid_outlet.mass_frac_comp[0, "Fe3O4"]
                / iron_oc.fs.unit.solid_outlet_block[0]._params.mw_comp["Fe3O4"]
            )
        )
        stoichiometric_ratio = mole_solid_reacted / mole_gas_reacted
        assert pytest.approx(4, abs=1e-6) == stoichiometric_ratio

    @pytest.mark.ui
    @pytest.mark.unit
    def test_report(self, iron_oc):
        iron_oc.fs.unit.report()

    @pytest.mark.component
    def test_initialization_error(self, iron_oc):
        iron_oc.fs.unit.gas_outlet.flow_mol[0].fix(1)

        with pytest.raises(InitializationError):
            iron_oc.fs.unit.initialize()


# -----------------------------------------------------------------------------
class TestIronOC_TransformationMethod(object):
    @pytest.fixture(scope="class")
    def iron_oc(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        # Set up thermo props and reaction props
        m.fs.gas_properties = GasPhaseParameterBlock()
        m.fs.solid_properties = SolidPhaseParameterBlock()
        m.fs.hetero_reactions = HeteroReactionParameterBlock(
            solid_property_package=m.fs.solid_properties,
            gas_property_package=m.fs.gas_properties,
        )

        m.fs.unit = BubblingFluidizedBed(
            transformation_method="dae.collocation",
            gas_phase_config={"property_package": m.fs.gas_properties},
            solid_phase_config={
                "property_package": m.fs.solid_properties,
                "reaction_package": m.fs.hetero_reactions,
            },
        )

        # Fix geometry variables
        m.fs.unit.number_orifice.fix(2500)  # [-]
        m.fs.unit.bed_diameter.fix(6.5)  # m
        m.fs.unit.bed_height.fix(5)  # m

        # Fix inlet port variables for gas and solid
        m.fs.unit.gas_inlet.flow_mol[0].fix(1567.79)  # mol/s
        m.fs.unit.gas_inlet.temperature[0].fix(373)  # K
        m.fs.unit.gas_inlet.pressure[0].fix(1.86e5)  # Pa = 1E5 bar
        m.fs.unit.gas_inlet.mole_frac_comp[0, "O2"].fix(0.2095)
        m.fs.unit.gas_inlet.mole_frac_comp[0, "N2"].fix(0.7808)
        m.fs.unit.gas_inlet.mole_frac_comp[0, "CO2"].fix(0.0004)
        m.fs.unit.gas_inlet.mole_frac_comp[0, "H2O"].fix(0.0093)

        m.fs.unit.solid_inlet.flow_mass[0].fix(1230.865)  # kg/s
        m.fs.unit.solid_inlet.particle_porosity[0].fix(0.27)  # (-)
        m.fs.unit.solid_inlet.temperature[0].fix(1173.9)  # K
        m.fs.unit.solid_inlet.mass_frac_comp[0, "Fe2O3"].fix(0.244162011502)
        m.fs.unit.solid_inlet.mass_frac_comp[0, "Fe3O4"].fix(0.201998299487)
        m.fs.unit.solid_inlet.mass_frac_comp[0, "Al2O3"].fix(0.553839689011)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, iron_oc):
        assert iron_oc.fs.unit.config.transformation_method == "dae.collocation"
        assert hasattr(iron_oc.fs.unit, "gas_inlet")
        assert len(iron_oc.fs.unit.gas_inlet.vars) == 4
        assert isinstance(iron_oc.fs.unit.gas_inlet.flow_mol, Var)
        assert isinstance(iron_oc.fs.unit.gas_inlet.mole_frac_comp, Var)
        assert isinstance(iron_oc.fs.unit.gas_inlet.temperature, Var)
        assert isinstance(iron_oc.fs.unit.gas_inlet.pressure, Var)

        assert hasattr(iron_oc.fs.unit, "solid_inlet")
        assert len(iron_oc.fs.unit.solid_inlet.vars) == 4
        assert isinstance(iron_oc.fs.unit.solid_inlet.flow_mass, Var)
        assert isinstance(iron_oc.fs.unit.solid_inlet.particle_porosity, Var)
        assert isinstance(iron_oc.fs.unit.solid_inlet.mass_frac_comp, Var)
        assert isinstance(iron_oc.fs.unit.solid_inlet.temperature, Var)

        assert hasattr(iron_oc.fs.unit, "gas_outlet")
        assert len(iron_oc.fs.unit.gas_outlet.vars) == 4
        assert isinstance(iron_oc.fs.unit.gas_outlet.flow_mol, Var)
        assert isinstance(iron_oc.fs.unit.gas_outlet.mole_frac_comp, Var)
        assert isinstance(iron_oc.fs.unit.gas_outlet.temperature, Var)
        assert isinstance(iron_oc.fs.unit.gas_outlet.pressure, Var)

        assert hasattr(iron_oc.fs.unit, "solid_outlet")
        assert len(iron_oc.fs.unit.solid_outlet.vars) == 4
        assert isinstance(iron_oc.fs.unit.solid_outlet.flow_mass, Var)
        assert isinstance(iron_oc.fs.unit.solid_outlet.particle_porosity, Var)
        assert isinstance(iron_oc.fs.unit.solid_outlet.mass_frac_comp, Var)
        assert isinstance(iron_oc.fs.unit.solid_outlet.temperature, Var)

        assert isinstance(iron_oc.fs.unit.orifice_area, Constraint)
        assert isinstance(iron_oc.fs.unit.bed_area_eqn, Constraint)
        assert isinstance(iron_oc.fs.unit.bubble_area, Constraint)
        assert isinstance(iron_oc.fs.unit.gas_emulsion_area, Constraint)
        assert isinstance(iron_oc.fs.unit.solid_emulsion_area, Constraint)
        assert isinstance(iron_oc.fs.unit.bubble_cloud_heat_trans_coeff, Constraint)
        assert isinstance(iron_oc.fs.unit.convective_heat_trans_coeff, Constraint)
        assert isinstance(iron_oc.fs.unit.convective_heat_transfer, Constraint)
        assert isinstance(iron_oc.fs.unit.bubble_cloud_bulk_heat_trans, Constraint)
        assert isinstance(iron_oc.fs.unit.bubble_mass_transfer, Constraint)
        assert isinstance(iron_oc.fs.unit.gas_emulsion_mass_transfer, Constraint)
        assert isinstance(iron_oc.fs.unit.bubble_heat_transfer, Constraint)
        assert isinstance(iron_oc.fs.unit.gas_emulsion_heat_transfer, Constraint)
        assert isinstance(iron_oc.fs.unit.solid_emulsion_heat_transfer, Constraint)

        assert number_variables(iron_oc) == 2335
        assert number_total_constraints(iron_oc) == 2287
        assert number_unused_variables(iron_oc) == 21

    @pytest.mark.unit
    def test_dof(self, iron_oc):
        assert degrees_of_freedom(iron_oc) == 0

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, iron_oc):
        BFB = iron_oc.fs.unit  # alias to keep test lines short
        iscale.calculate_scaling_factors(BFB)
        initialization_tester(
            iron_oc,
            optarg={"tol": 1e-6},
            gas_phase_state_args={
                "flow_mol": 1567.79,
                "temperature": 1173.9,
                "pressure": 1.86e5,
            },
            solid_phase_state_args={"flow_mass": 1230.865, "temperature": 1173.9},
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, iron_oc):
        results = solver.solve(iron_oc, tee=True)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, iron_oc):
        assert (
            pytest.approx(0.80244, abs=1e-5)
            == iron_oc.fs.unit.velocity_superficial_gas[0, 0].value
        )
        assert (
            pytest.approx(3.02457, abs=1e-5)
            == iron_oc.fs.unit.velocity_superficial_gas[0, 1].value
        )
        assert (
            pytest.approx(0.03431, abs=1e-5)
            == iron_oc.fs.unit.bubble_diameter[0, 0].value
        )
        assert (
            pytest.approx(1.78570, abs=1e-5)
            == iron_oc.fs.unit.bubble_diameter[0, 1].value
        )
        assert (
            pytest.approx(1.17521, abs=1e-5)
            == iron_oc.fs.unit.velocity_bubble[0, 0].value
        )
        assert (
            pytest.approx(5.96027, abs=1e-5)
            == iron_oc.fs.unit.velocity_bubble[0, 1].value
        )
        assert pytest.approx(0.64909, abs=1e-5) == iron_oc.fs.unit.delta[0, 0].value
        assert pytest.approx(0.50080, abs=1e-5) == iron_oc.fs.unit.delta[0, 1].value
        assert (
            pytest.approx(145217.21159, abs=1e-5)
            == iron_oc.fs.unit.gas_outlet.pressure[0].value
        )
        assert (
            pytest.approx(40782.78841, abs=1e-5)
            == iron_oc.fs.unit.gas_inlet.pressure[0].value
            - iron_oc.fs.unit.gas_outlet.pressure[0].value
        )

    @pytest.mark.component
    def test_units_consistent(self, iron_oc):
        assert_units_consistent(iron_oc)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, iron_oc):
        # Conservation of material check
        calculate_variable_from_constraint(
            iron_oc.fs.unit.gas_inlet_block[0].mw,
            iron_oc.fs.unit.gas_inlet_block[0].mw_eqn,
        )
        calculate_variable_from_constraint(
            iron_oc.fs.unit.gas_outlet_block[0].mw,
            iron_oc.fs.unit.gas_outlet_block[0].mw_eqn,
        )
        mbal_gas = value(
            (
                iron_oc.fs.unit.gas_inlet.flow_mol[0]
                * iron_oc.fs.unit.gas_inlet_block[0].mw
            )
            - (
                iron_oc.fs.unit.gas_outlet.flow_mol[0]
                * iron_oc.fs.unit.gas_outlet_block[0].mw
            )
        )
        mbal_solid = value(
            iron_oc.fs.unit.solid_inlet.flow_mass[0]
            - iron_oc.fs.unit.solid_outlet.flow_mass[0]
        )
        mbal_tol = mbal_gas + mbal_solid
        assert abs(mbal_tol) <= 1e-2

        # Reaction stoichiometric ratio check
        # Overall reducer reactions for methane combustion:
        # CH4 + 12Fe2O3 => 8Fe3O4 + CO2 + 2H2O
        mole_gas_reacted = value(
            iron_oc.fs.unit.gas_inlet.flow_mol[0]
            * iron_oc.fs.unit.gas_inlet.mole_frac_comp[0, "O2"]
            - iron_oc.fs.unit.gas_outlet.flow_mol[0]
            * iron_oc.fs.unit.gas_outlet.mole_frac_comp[0, "O2"]
        )
        mole_solid_reacted = value(
            (
                iron_oc.fs.unit.solid_inlet.flow_mass[0]
                * iron_oc.fs.unit.solid_inlet.mass_frac_comp[0, "Fe3O4"]
                / iron_oc.fs.unit.solid_inlet_block[0]._params.mw_comp["Fe3O4"]
            )
            - (
                iron_oc.fs.unit.solid_outlet.flow_mass[0]
                * iron_oc.fs.unit.solid_outlet.mass_frac_comp[0, "Fe3O4"]
                / iron_oc.fs.unit.solid_outlet_block[0]._params.mw_comp["Fe3O4"]
            )
        )
        stoichiometric_ratio = mole_solid_reacted / mole_gas_reacted
        assert pytest.approx(4, abs=1e-6) == stoichiometric_ratio

        # Conservation of energy check
        ebal_gas = value(
            (
                iron_oc.fs.unit.gas_inlet.flow_mol[0]
                * iron_oc.fs.unit.gas_inlet_block[0].enth_mol
            )
            - (
                iron_oc.fs.unit.gas_outlet.flow_mol[0]
                * iron_oc.fs.unit.gas_outlet_block[0].enth_mol
            )
        )
        ebal_solid = value(
            (
                iron_oc.fs.unit.solid_inlet.flow_mass[0]
                * iron_oc.fs.unit.solid_inlet_block[0].enth_mass
            )
            - (
                iron_oc.fs.unit.solid_outlet.flow_mass[0]
                * iron_oc.fs.unit.solid_outlet_block[0].enth_mass
            )
        )
        e_reaction = value(
            mole_gas_reacted
            * iron_oc.fs.unit.solid_emulsion.reactions[0, 0]._params.dh_rxn["R1"]
        )
        ebal_tol = ebal_gas + ebal_solid - e_reaction
        assert abs(ebal_tol) <= 1e-2

    @pytest.mark.ui
    @pytest.mark.unit
    def test_report(self, iron_oc):
        iron_oc.fs.unit.report()

    @pytest.mark.component
    def test_initialization_error(self, iron_oc):
        iron_oc.fs.unit.gas_outlet.flow_mol[0].fix(1)

        with pytest.raises(InitializationError):
            iron_oc.fs.unit.initialize()
