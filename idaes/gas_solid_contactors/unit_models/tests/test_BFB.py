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
Test for ControlVolumeBlockData, and for initializing the
bubbling fluidized bed module

Author: Chinedu Okoli
"""

import pytest

from pyomo.environ import (ConcreteModel,
                           TerminationCondition,
                           SolverStatus,
                           value,
                           Var,
                           Constraint)
from pyomo.common.config import ConfigBlock
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from idaes.core import (FlowsheetBlock,
                        MaterialBalanceType,
                        EnergyBalanceType,
                        MomentumBalanceType)
from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              number_variables,
                                              number_total_constraints,
                                              number_unused_variables,
                                              unused_variables_set)
from idaes.core.util.testing import initialization_tester
from idaes.core.util import get_solver

from idaes.gas_solid_contactors.unit_models. \
    bubbling_fluidized_bed import BubblingFluidizedBed
from idaes.gas_solid_contactors.properties.methane_iron_OC_reduction. \
    gas_phase_thermo import GasPhaseParameterBlock
from idaes.gas_solid_contactors.properties.methane_iron_OC_reduction. \
    solid_phase_thermo import SolidPhaseParameterBlock
from idaes.gas_solid_contactors.properties.methane_iron_OC_reduction. \
    hetero_reactions import HeteroReactionParameterBlock

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


# -----------------------------------------------------------------------------
@pytest.mark.unit
def test_config():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    # Set up thermo props and reaction props
    m.fs.gas_properties = GasPhaseParameterBlock()
    m.fs.solid_properties = SolidPhaseParameterBlock()
    m.fs.hetero_reactions = HeteroReactionParameterBlock(
            default={"solid_property_package": m.fs.solid_properties,
                     "gas_property_package": m.fs.gas_properties})

    m.fs.unit = BubblingFluidizedBed(
            default={
                    "gas_phase_config":
                    {"property_package": m.fs.gas_properties},
                    "solid_phase_config":
                    {"property_package": m.fs.solid_properties,
                     "reaction_package": m.fs.hetero_reactions
                     }})

    # Check unit config arguments
    assert len(m.fs.unit.config) == 14
    assert isinstance(m.fs.unit.config.gas_phase_config, ConfigBlock)
    assert isinstance(m.fs.unit.config.solid_phase_config, ConfigBlock)

    assert m.fs.unit.config.finite_elements == 10
    assert m.fs.unit.config.length_domain_set == [0.0, 1.0]
    assert m.fs.unit.config.transformation_method == "dae.finite_difference"
    assert m.fs.unit.config.transformation_scheme == 'BACKWARD'
    assert m.fs.unit.config.collocation_points == 3
    assert m.fs.unit.config.flow_type == "co_current"
    assert m.fs.unit.config.material_balance_type == \
        MaterialBalanceType.componentTotal
    assert m.fs.unit.config.energy_balance_type == \
        EnergyBalanceType.enthalpyTotal
    assert m.fs.unit.config.momentum_balance_type == \
        MomentumBalanceType.pressureTotal
    assert m.fs.unit.config.has_pressure_change is True

    # Check gas phase config arguments
    assert len(m.fs.unit.config.gas_phase_config) == 7
    assert m.fs.unit.config.gas_phase_config.has_equilibrium_reactions is False
    assert m.fs.unit.config.gas_phase_config.property_package is \
        m.fs.gas_properties
    assert m.fs.unit.config.gas_phase_config.reaction_package is None

    # Check solid phase config arguments
    assert len(m.fs.unit.config.solid_phase_config) == 7
    assert m.fs.unit.config.solid_phase_config.has_equilibrium_reactions is \
        False
    assert m.fs.unit.config.solid_phase_config.property_package is \
        m.fs.solid_properties
    assert m.fs.unit.config.solid_phase_config.reaction_package is \
        m.fs.hetero_reactions


# -----------------------------------------------------------------------------
class TestIronOC(object):
    @pytest.fixture(scope="class")
    def iron_oc(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        # Set up thermo props and reaction props
        m.fs.gas_properties = GasPhaseParameterBlock()
        m.fs.solid_properties = SolidPhaseParameterBlock()
        m.fs.hetero_reactions = HeteroReactionParameterBlock(
                default={"solid_property_package": m.fs.solid_properties,
                         "gas_property_package": m.fs.gas_properties})

        m.fs.unit = BubblingFluidizedBed(
                default={
                        "gas_phase_config":
                        {"property_package": m.fs.gas_properties},
                        "solid_phase_config":
                        {"property_package": m.fs.solid_properties,
                         "reaction_package": m.fs.hetero_reactions
                         }})

        # Fix geometry variables
        m.fs.unit.number_orifice.fix(2500)  # [-]
        m.fs.unit.bed_diameter.fix(6.5)  # m
        m.fs.unit.bed_height.fix(5)  # m

        # Fix inlet port variables for gas and solid
        m.fs.unit.gas_inlet.flow_mol[0].fix(272.81)  # mol/s
        m.fs.unit.gas_inlet.temperature[0].fix(373)  # K
        m.fs.unit.gas_inlet.pressure[0].fix(1.86)  # bar
        m.fs.unit.gas_inlet.mole_frac_comp[0, "CO2"].fix(0.4772)
        m.fs.unit.gas_inlet.mole_frac_comp[0, "H2O"].fix(0.0646)
        m.fs.unit.gas_inlet.mole_frac_comp[0, "CH4"].fix(0.4582)

        m.fs.unit.solid_inlet.flow_mass[0].fix(1422)  # kg/s
        m.fs.unit.solid_inlet.particle_porosity[0].fix(0.27)  # (-)
        m.fs.unit.solid_inlet.temperature[0].fix(1186)  # K
        m.fs.unit.solid_inlet.mass_frac_comp[0, "Fe2O3"].fix(0.45)
        m.fs.unit.solid_inlet.mass_frac_comp[0, "Fe3O4"].fix(1e-9)
        m.fs.unit.solid_inlet.mass_frac_comp[0, "Al2O3"].fix(0.55)

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
        assert isinstance(iron_oc.fs.unit.bubble_cloud_heat_trans_coeff,
                          Constraint)
        assert isinstance(iron_oc.fs.unit.convective_heat_trans_coeff,
                          Constraint)
        assert isinstance(iron_oc.fs.unit.convective_heat_transfer,
                          Constraint)
        assert isinstance(iron_oc.fs.unit.bubble_cloud_bulk_heat_trans,
                          Constraint)
        assert isinstance(iron_oc.fs.unit.bubble_mass_transfer,
                          Constraint)
        assert isinstance(iron_oc.fs.unit.gas_emulsion_mass_transfer,
                          Constraint)
        assert isinstance(iron_oc.fs.unit.bubble_heat_transfer, Constraint)
        assert isinstance(iron_oc.fs.unit.gas_emulsion_heat_transfer,
                          Constraint)
        assert isinstance(iron_oc.fs.unit.solid_emulsion_heat_transfer,
                          Constraint)

        assert number_variables(iron_oc) == 1436
        assert number_total_constraints(iron_oc) == 1392
        assert number_unused_variables(iron_oc) == 18

    @pytest.mark.unit
    def test_dof(self, iron_oc):
        assert degrees_of_freedom(iron_oc) == 0

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, iron_oc):
        initialization_tester(
                iron_oc,
                optarg={'tol': 1e-6},
                gas_phase_state_args={"flow_mol": 272.81,
                                      "temperature": 1186,
                                      "pressure": 1.86},
                solid_phase_state_args={"flow_mass": 1230,
                                        "temperature": 1186})

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, iron_oc):
        results = solver.solve(iron_oc)

        # Check for optimal solution
        assert results.solver.termination_condition == \
            TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, iron_oc):
        assert (pytest.approx(0.14, abs=1e-2) ==
                iron_oc.fs.unit.velocity_superficial_gas[0, 0].value)
        assert (pytest.approx(1.052, abs=1e-2) ==
                iron_oc.fs.unit.velocity_superficial_gas[0, 1].value)
        assert (pytest.approx(0.015, abs=1e-2) ==
                iron_oc.fs.unit.bubble_diameter[0, 0].value)
        assert (pytest.approx(1.04, abs=1e-2) ==
                iron_oc.fs.unit.bubble_diameter[0, 1].value)
        assert (pytest.approx(0.375, abs=1e-2) ==
                iron_oc.fs.unit.velocity_bubble[0, 0].value)
        assert (pytest.approx(3.290, abs=1e-2) ==
                iron_oc.fs.unit.velocity_bubble[0, 1].value)
        assert (pytest.approx(0.267, abs=1e-2) ==
                iron_oc.fs.unit.delta[0, 0].value)
        assert (pytest.approx(0.307, abs=1e-2) ==
                iron_oc.fs.unit.delta[0, 1].value)
        assert (pytest.approx(1.24, abs=1e-2) ==
                iron_oc.fs.unit.gas_outlet.pressure[0].value)
        # Check that pressure drop occurs across the bed
        assert value(
                iron_oc.fs.unit.gas_inlet.pressure[0] -
                iron_oc.fs.unit.gas_outlet.pressure[0]) >= 0

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, iron_oc):
        # Conservation of material check
        calculate_variable_from_constraint(
                    iron_oc.fs.unit.gas_inlet_block[0].mw,
                    iron_oc.fs.unit.gas_inlet_block[0].mw_eqn)
        calculate_variable_from_constraint(
                    iron_oc.fs.unit.gas_outlet_block[0].mw,
                    iron_oc.fs.unit.gas_outlet_block[0].mw_eqn)
        mbal_gas = value(
                (iron_oc.fs.unit.gas_inlet.flow_mol[0] *
                 iron_oc.fs.unit.gas_inlet_block[0].mw) -
                (iron_oc.fs.unit.gas_outlet.flow_mol[0] *
                 iron_oc.fs.unit.gas_outlet_block[0].mw))
        mbal_solid = value(
                iron_oc.fs.unit.solid_inlet.flow_mass[0] -
                iron_oc.fs.unit.solid_outlet.flow_mass[0])
        mbal_tol = mbal_gas + mbal_solid
        assert abs(mbal_tol) <= 1e-2

        # Reaction stoichiometric ratio check
        # Overall reducer reactions for methane combustion:
        # CH4 + 12Fe2O3 => 8Fe3O4 + CO2 + 2H2O
        mole_gas_reacted = value(
                iron_oc.fs.unit.gas_inlet.flow_mol[0] *
                iron_oc.fs.unit.gas_inlet.mole_frac_comp[0, 'CH4'] -
                iron_oc.fs.unit.gas_outlet.flow_mol[0] *
                iron_oc.fs.unit.gas_outlet.mole_frac_comp[0, 'CH4'])
        mole_solid_reacted = value(
            (iron_oc.fs.unit.solid_inlet.flow_mass[0] *
             iron_oc.fs.unit.solid_inlet.mass_frac_comp[0, 'Fe2O3'] /
             iron_oc.fs.unit.solid_inlet_block[0]._params.mw_comp['Fe2O3']) -
            (iron_oc.fs.unit.solid_outlet.flow_mass[0] *
             iron_oc.fs.unit.solid_outlet.mass_frac_comp[0, 'Fe2O3'] /
             iron_oc.fs.unit.solid_outlet_block[0]._params.mw_comp['Fe2O3']))
        stoichiometric_ratio = mole_solid_reacted/mole_gas_reacted
        assert (pytest.approx(12, abs=1e-6) == stoichiometric_ratio)

        # Conservation of energy check
        ebal_gas = value(
            (iron_oc.fs.unit.gas_inlet.flow_mol[0] *
             iron_oc.fs.unit.gas_inlet_block[0].enth_mol) -
            (iron_oc.fs.unit.gas_outlet.flow_mol[0] *
             iron_oc.fs.unit.gas_outlet_block[0].enth_mol))
        ebal_solid = value(
            (iron_oc.fs.unit.solid_inlet.flow_mass[0] *
             iron_oc.fs.unit.solid_inlet_block[0].enth_mass) -
            (iron_oc.fs.unit.solid_outlet.flow_mass[0] *
             iron_oc.fs.unit.solid_outlet_block[0].enth_mass))
        e_reaction = value(
                mole_gas_reacted *
                iron_oc.fs.unit.solid_emulsion.reactions[0, 0].
                _params.dh_rxn["R1"])
        ebal_tol = ebal_gas + ebal_solid - e_reaction
        assert abs(ebal_tol) <= 1e-2

    @pytest.mark.ui
    @pytest.mark.unit
    def test_report(self, iron_oc):
        iron_oc.fs.unit.report()


# -----------------------------------------------------------------------------
class TestIronOC_EnergyBalanceType(object):
    @pytest.fixture(scope="class")
    def iron_oc(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        # Set up thermo props and reaction props
        m.fs.gas_properties = GasPhaseParameterBlock()
        m.fs.solid_properties = SolidPhaseParameterBlock()
        m.fs.hetero_reactions = HeteroReactionParameterBlock(
                default={"solid_property_package": m.fs.solid_properties,
                         "gas_property_package": m.fs.gas_properties})

        m.fs.unit = BubblingFluidizedBed(
                default={
                        "energy_balance_type": EnergyBalanceType.none,
                        "gas_phase_config":
                        {"property_package": m.fs.gas_properties},
                        "solid_phase_config":
                        {"property_package": m.fs.solid_properties,
                         "reaction_package": m.fs.hetero_reactions
                         }})

        # # Fix geometry variables
        m.fs.unit.number_orifice.fix(2500)  # [-]
        m.fs.unit.bed_diameter.fix(6.5)  # m
        m.fs.unit.bed_height.fix(5)  # m

        # # Fix inlet port variables for gas and solid
        m.fs.unit.gas_inlet.flow_mol[0].fix(272.81)  # mol/s
        m.fs.unit.gas_inlet.temperature[0].fix(1186)  # K
        m.fs.unit.gas_inlet.pressure[0].fix(1.86)  # bar
        m.fs.unit.gas_inlet.mole_frac_comp[0, "CO2"].fix(0.4772)
        m.fs.unit.gas_inlet.mole_frac_comp[0, "H2O"].fix(0.0646)
        m.fs.unit.gas_inlet.mole_frac_comp[0, "CH4"].fix(0.4582)

        m.fs.unit.solid_inlet.flow_mass[0].fix(1422)  # kg/
        m.fs.unit.solid_inlet.particle_porosity.fix(0.27)  # (-)
        m.fs.unit.solid_inlet.temperature[0].fix(1186)  # K
        m.fs.unit.solid_inlet.mass_frac_comp[0, "Fe2O3"].fix(0.45)
        m.fs.unit.solid_inlet.mass_frac_comp[0, "Fe3O4"].fix(1e-9)
        m.fs.unit.solid_inlet.mass_frac_comp[0, "Al2O3"].fix(0.55)

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
        assert isinstance(iron_oc.fs.unit.isothermal_solid_emulsion,
                          Constraint)
        assert isinstance(iron_oc.fs.unit.isothermal_bubble, Constraint)

        assert number_variables(iron_oc) == 1156
        assert number_total_constraints(iron_oc) == 1049
        assert number_unused_variables(iron_oc) == 82
        print(unused_variables_set(iron_oc))

    @pytest.mark.unit
    def test_dof(self, iron_oc):
        assert degrees_of_freedom(iron_oc) == 0

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, iron_oc):
        initialization_tester(
                iron_oc,
                optarg={'tol': 1e-6},
                gas_phase_state_args={"flow_mol": 272.81,
                                      "temperature": 1186,
                                      "pressure": 1.86},
                solid_phase_state_args={"flow_mass": 1422,
                                        "temperature": 1186})

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, iron_oc):
        results = solver.solve(iron_oc)

        # Check for optimal solution
        assert results.solver.termination_condition == \
            TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, iron_oc):
        assert (pytest.approx(0.44, abs=1e-2) ==
                iron_oc.fs.unit.velocity_superficial_gas[0, 0].value)
        assert (pytest.approx(1.05, abs=1e-2) ==
                iron_oc.fs.unit.velocity_superficial_gas[0, 1].value)
        assert (pytest.approx(0.03, abs=1e-2) ==
                iron_oc.fs.unit.bubble_diameter[0, 0].value)
        assert (pytest.approx(1.05, abs=1e-2) ==
                iron_oc.fs.unit.bubble_diameter[0, 1].value)
        assert (pytest.approx(0.77, abs=1e-2) ==
                iron_oc.fs.unit.velocity_bubble[0, 0].value)
        assert (pytest.approx(3.30, abs=1e-2) ==
                iron_oc.fs.unit.velocity_bubble[0, 1].value)
        assert (pytest.approx(0.53, abs=1e-2) ==
                iron_oc.fs.unit.delta[0, 0].value)
        assert (pytest.approx(0.307, abs=1e-2) ==
                iron_oc.fs.unit.delta[0, 1].value)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, iron_oc):
        # Conservation of material check
        calculate_variable_from_constraint(
                    iron_oc.fs.unit.gas_inlet_block[0].mw,
                    iron_oc.fs.unit.gas_inlet_block[0].mw_eqn)
        calculate_variable_from_constraint(
                    iron_oc.fs.unit.gas_outlet_block[0].mw,
                    iron_oc.fs.unit.gas_outlet_block[0].mw_eqn)
        mbal_gas = value(
                (iron_oc.fs.unit.gas_inlet.flow_mol[0] *
                 iron_oc.fs.unit.gas_inlet_block[0].mw) -
                (iron_oc.fs.unit.gas_outlet.flow_mol[0] *
                 iron_oc.fs.unit.gas_outlet_block[0].mw))
        mbal_solid = value(
                iron_oc.fs.unit.solid_inlet.flow_mass[0] -
                iron_oc.fs.unit.solid_outlet.flow_mass[0])
        mbal_tol = mbal_gas + mbal_solid
        assert abs(mbal_tol) <= 1e-2

        # Reaction stoichiometric ratio check
        # Overall reducer reactions for methane combustion:
        # CH4 + 12Fe2O3 => 8Fe3O4 + CO2 + 2H2O
        mole_gas_reacted = value(
                iron_oc.fs.unit.gas_inlet.flow_mol[0] *
                iron_oc.fs.unit.gas_inlet.mole_frac_comp[0, 'CH4'] -
                iron_oc.fs.unit.gas_outlet.flow_mol[0] *
                iron_oc.fs.unit.gas_outlet.mole_frac_comp[0, 'CH4'])
        mole_solid_reacted = value(
            (iron_oc.fs.unit.solid_inlet.flow_mass[0] *
             iron_oc.fs.unit.solid_inlet.mass_frac_comp[0, 'Fe2O3'] /
             iron_oc.fs.unit.solid_inlet_block[0]._params.mw_comp['Fe2O3']) -
            (iron_oc.fs.unit.solid_outlet.flow_mass[0] *
             iron_oc.fs.unit.solid_outlet.mass_frac_comp[0, 'Fe2O3'] /
             iron_oc.fs.unit.solid_outlet_block[0]._params.mw_comp['Fe2O3']))
        stoichiometric_ratio = mole_solid_reacted/mole_gas_reacted
        assert (pytest.approx(12, abs=1e-6) == stoichiometric_ratio)

    @pytest.mark.ui
    @pytest.mark.unit
    def test_report(self, iron_oc):
        iron_oc.fs.unit.report()
