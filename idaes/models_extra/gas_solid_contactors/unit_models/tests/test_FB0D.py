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
Tests for ControlVolumeBlockData, and for initializing the 0D fixed bed module

Author: Brandon Paul
"""

import pytest

from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    TransformationFactory,
    value,
    units as pyunits,
    Constraint,
)
from pyomo.util.check_units import assert_units_consistent
from idaes.core import FlowsheetBlock, EnergyBalanceType
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
    unused_variables_set,
)
from idaes.core.util.testing import initialization_tester
from idaes.core.util.initialization import initialize_by_time_element
from idaes.core.util import scaling as iscale
from idaes.core.solvers import get_solver
from idaes.core.util.exceptions import InitializationError

# Import FixedBed0D unit model
from idaes.models_extra.gas_solid_contactors.unit_models.fixed_bed_0D import FixedBed0D

# Import property packages
from idaes.models_extra.gas_solid_contactors.properties.methane_iron_OC_reduction.gas_phase_thermo import (
    GasPhaseParameterBlock,
)
from idaes.models_extra.gas_solid_contactors.properties.methane_iron_OC_reduction.solid_phase_thermo import (
    SolidPhaseParameterBlock,
)
from idaes.models_extra.gas_solid_contactors.properties.methane_iron_OC_reduction.hetero_reactions import (
    HeteroReactionParameterBlock,
)

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


# -----------------------------------------------------------------------------
@pytest.mark.unit
def test_config():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_set=[0, 3600], time_units=pyunits.s)

    # Set up thermo props and reaction props
    m.fs.gas_props = GasPhaseParameterBlock()
    m.fs.solid_props = SolidPhaseParameterBlock()
    m.fs.solid_rxns = HeteroReactionParameterBlock(
        solid_property_package=m.fs.solid_props, gas_property_package=m.fs.gas_props
    )

    m.fs.unit = FixedBed0D(
        gas_property_package=m.fs.gas_props,
        solid_property_package=m.fs.solid_props,
        reaction_package=m.fs.solid_rxns,
    )

    # Check unit config arguments
    assert len(m.fs.unit.config) == 9

    assert m.fs.unit.config.dynamic is True
    assert m.fs.unit.config.has_holdup is True
    assert m.fs.unit.config.energy_balance_type == EnergyBalanceType.enthalpyTotal
    assert m.fs.unit.config.gas_property_package is m.fs.gas_props
    assert m.fs.unit.config.gas_property_package_args is None
    assert m.fs.unit.config.solid_property_package is m.fs.solid_props
    assert m.fs.unit.config.solid_property_package_args is None
    assert m.fs.unit.config.reaction_package is m.fs.solid_rxns


# -----------------------------------------------------------------------------
class TestIronOC(object):
    @pytest.fixture(scope="class")
    def iron_oc(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=True, time_set=[0, 3600], time_units=pyunits.s)

        m.fs.gas_props = GasPhaseParameterBlock()
        m.fs.solid_props = SolidPhaseParameterBlock()
        m.fs.solid_rxns = HeteroReactionParameterBlock(
            solid_property_package=m.fs.solid_props, gas_property_package=m.fs.gas_props
        )

        m.fs.unit = FixedBed0D(
            gas_property_package=m.fs.gas_props,
            solid_property_package=m.fs.solid_props,
            reaction_package=m.fs.solid_rxns,
        )

        # Discretize time domain
        m.discretizer = TransformationFactory("dae.finite_difference")
        m.discretizer.apply_to(m, nfe=100, wrt=m.fs.time, scheme="BACKWARD")

        # Set reactor design conditions
        m.fs.unit.bed_diameter.fix(1)  # diameter of the TGA reactor [m]
        m.fs.unit.bed_height.fix(1)  # height of solids in the TGA reactor [m]

        # Set initial conditions of the solid phase
        m.fs.unit.solids[0].temperature.fix(1273.15)
        m.fs.unit.solids[0].particle_porosity.fix(0.20)
        m.fs.unit.solids[0].mass_frac_comp["Fe2O3"].fix(0.45)
        m.fs.unit.solids[0].mass_frac_comp["Fe3O4"].fix(0)
        m.fs.unit.solids[0].mass_frac_comp["Al2O3"].fix(0.55)

        # Set conditions of the gas phase (this is all fixed as gas side
        # assumption is excess gas flowrate which means all state variables
        # remain unchanged)
        for t in m.fs.time:
            m.fs.unit.gas[t].temperature.fix(1273.15)
            m.fs.unit.gas[t].pressure.fix(1.01325e5)  # 1atm
            m.fs.unit.gas[t].mole_frac_comp["CO2"].fix(0.4)
            m.fs.unit.gas[t].mole_frac_comp["H2O"].fix(0.5)
            m.fs.unit.gas[t].mole_frac_comp["CH4"].fix(0.1)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, iron_oc):
        assert isinstance(iron_oc.fs.unit.volume_bed_constraint, Constraint)
        assert isinstance(iron_oc.fs.unit.volume_solid_constraint, Constraint)
        assert isinstance(
            iron_oc.fs.unit.solids_material_holdup_constraints, Constraint
        )
        assert isinstance(
            iron_oc.fs.unit.solids_material_accumulation_constraints, Constraint
        )
        assert isinstance(iron_oc.fs.unit.mass_solids_constraint, Constraint)
        assert isinstance(iron_oc.fs.unit.sum_component_constraint, Constraint)
        assert isinstance(iron_oc.fs.unit.solids_energy_holdup_constraints, Constraint)
        assert isinstance(
            iron_oc.fs.unit.solids_energy_accumulation_constraints, Constraint
        )

        assert number_variables(iron_oc) == 3650
        assert number_total_constraints(iron_oc) == 2925
        assert number_unused_variables(iron_oc) == 206

    @pytest.mark.unit
    def test_dof(self, iron_oc):
        assert degrees_of_freedom(iron_oc) == 0

    @pytest.fixture(scope="class")
    def iron_oc_unscaled(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=True, time_set=[0, 360], time_units=pyunits.s)

        m.fs.gas_props = GasPhaseParameterBlock()
        m.fs.solid_props = SolidPhaseParameterBlock()
        m.fs.solid_rxns = HeteroReactionParameterBlock(
            solid_property_package=m.fs.solid_props, gas_property_package=m.fs.gas_props
        )

        m.fs.unit = FixedBed0D(
            gas_property_package=m.fs.gas_props,
            solid_property_package=m.fs.solid_props,
            reaction_package=m.fs.solid_rxns,
        )

        # Discretize time domain
        m.discretizer = TransformationFactory("dae.finite_difference")
        m.discretizer.apply_to(m, nfe=10, wrt=m.fs.time, scheme="BACKWARD")

        # Set reactor design conditions
        m.fs.unit.bed_diameter.fix(1)  # diameter of the TGA reactor [m]
        m.fs.unit.bed_height.fix(1)  # height of solids in the TGA reactor [m]

        # Set initial conditions of the solid phase
        m.fs.unit.solids[0].temperature.fix(1273.15)
        m.fs.unit.solids[0].particle_porosity.fix(0.20)
        m.fs.unit.solids[0].mass_frac_comp["Fe2O3"].fix(0.45)
        m.fs.unit.solids[0].mass_frac_comp["Fe3O4"].fix(0)
        m.fs.unit.solids[0].mass_frac_comp["Al2O3"].fix(0.55)

        # Set conditions of the gas phase (this is all fixed as gas side
        # assumption is excess gas flowrate which means all state variables
        # remain unchanged)
        for t in m.fs.time:
            m.fs.unit.gas[t].temperature.fix(1273.15)
            m.fs.unit.gas[t].pressure.fix(1.01325e5)  # 1atm
            m.fs.unit.gas[t].mole_frac_comp["CO2"].fix(0.4)
            m.fs.unit.gas[t].mole_frac_comp["H2O"].fix(0.5)
            m.fs.unit.gas[t].mole_frac_comp["CH4"].fix(0.1)

        return m

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize_unscaled(self, iron_oc_unscaled):
        optarg = {"tol": 1e-6}

        initialization_tester(iron_oc_unscaled, optarg=optarg)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize_by_time_unscaled(self, iron_oc_unscaled):

        optarg = {"tol": 1e-6}
        solver = get_solver("ipopt", optarg)  # create solver

        initialize_by_time_element(
            iron_oc_unscaled.fs, iron_oc_unscaled.fs.time, solver=solver
        )

        assert degrees_of_freedom(iron_oc_unscaled) == 0

        # Assert that model is still fixed and deactivated as expected
        assert iron_oc_unscaled.fs.unit.solids[
            iron_oc_unscaled.fs.time.first()
        ].particle_porosity.fixed

        for t in iron_oc_unscaled.fs.time:
            if t != iron_oc_unscaled.fs.time.first():
                assert not iron_oc_unscaled.fs.unit.solids[t].particle_porosity.fixed
            assert iron_oc_unscaled.fs.unit.solids[t].density_particle_constraint.active
            assert iron_oc_unscaled.fs.unit.solids[t].density_skeletal_constraint.active

        # Assert that constraints are feasible after initialization
        for con in iron_oc_unscaled.fs.component_data_objects(Constraint, active=True):
            assert value(con.body) - value(con.upper) < 1e-4
            assert value(con.lower) - value(con.body) < 1e-4

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve_unscaled(self, iron_oc_unscaled):
        results = solver.solve(iron_oc_unscaled)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.component
    def test_scaling(self, iron_oc):

        # Set scaling gas phase for state variables
        iron_oc.fs.gas_props.set_default_scaling("flow_mol", 1e-3)
        iron_oc.fs.gas_props.set_default_scaling("pressure", 1e-5)
        iron_oc.fs.gas_props.set_default_scaling("temperature", 1e-2)
        for comp in iron_oc.fs.gas_props.component_list:
            iron_oc.fs.gas_props.set_default_scaling("mole_frac_comp", 1e1, index=comp)
        # Set scaling for gas phase thermophysical and transport properties
        iron_oc.fs.gas_props.set_default_scaling("enth_mol", 1e-6)
        iron_oc.fs.gas_props.set_default_scaling("enth_mol_comp", 1e-6)
        iron_oc.fs.gas_props.set_default_scaling("cp_mol", 1e-6)
        iron_oc.fs.gas_props.set_default_scaling("cp_mol_comp", 1e-6)
        iron_oc.fs.gas_props.set_default_scaling("cp_mass", 1e-6)
        iron_oc.fs.gas_props.set_default_scaling("entr_mol", 1e-2)
        iron_oc.fs.gas_props.set_default_scaling("entr_mol_phase", 1e-2)
        iron_oc.fs.gas_props.set_default_scaling("dens_mol", 1)
        iron_oc.fs.gas_props.set_default_scaling("dens_mol_comp", 1)
        iron_oc.fs.gas_props.set_default_scaling("dens_mass", 1e2)
        iron_oc.fs.gas_props.set_default_scaling("visc_d_comp", 1e4)
        iron_oc.fs.gas_props.set_default_scaling("diffusion_comp", 1e5)
        iron_oc.fs.gas_props.set_default_scaling("therm_cond_comp", 1e2)
        iron_oc.fs.gas_props.set_default_scaling("visc_d", 1e5)
        iron_oc.fs.gas_props.set_default_scaling("therm_cond", 1e0)
        iron_oc.fs.gas_props.set_default_scaling("mw", 1e2)

        # Set scaling for solid phase state variables
        iron_oc.fs.solid_props.set_default_scaling("flow_mass", 1e-3)
        iron_oc.fs.solid_props.set_default_scaling("particle_porosity", 1e2)
        iron_oc.fs.solid_props.set_default_scaling("temperature", 1e-2)
        for comp in iron_oc.fs.solid_props.component_list:
            iron_oc.fs.solid_props.set_default_scaling(
                "mass_frac_comp", 1e1, index=comp
            )

        # Set scaling for solid phase thermophysical and transport properties
        iron_oc.fs.solid_props.set_default_scaling("enth_mass", 1e-6)
        iron_oc.fs.solid_props.set_default_scaling("enth_mol_comp", 1e-6)
        iron_oc.fs.solid_props.set_default_scaling("cp_mol_comp", 1e-6)
        iron_oc.fs.solid_props.set_default_scaling("cp_mass", 1e-6)
        iron_oc.fs.solid_props.set_default_scaling("dens_mass_particle", 1e-2)
        iron_oc.fs.solid_props.set_default_scaling("dens_mass_skeletal", 1e-2)

        FB0D = iron_oc.fs.unit  # alias to keep test lines short

        # Calculate scaling factors
        iscale.calculate_scaling_factors(FB0D)

        assert pytest.approx(0.10000, rel=1e-5) == iscale.get_scaling_factor(
            FB0D.bed_diameter
        )
        assert pytest.approx(0.10000, rel=1e-5) == iscale.get_scaling_factor(
            FB0D.bed_height
        )
        assert pytest.approx(0.127324, rel=1e-5) == iscale.get_scaling_factor(
            FB0D.volume_bed
        )

        for c in FB0D.volume_bed_constraint.values():
            assert pytest.approx(
                0.127324, rel=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for t, c in FB0D.volume_solid_constraint.items():
            assert pytest.approx(
                0.127324, rel=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, j), c in FB0D.solids_material_holdup_constraints.items():
            assert pytest.approx(
                0.00127324, rel=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, j), c in FB0D.solids_material_accumulation_constraints.items():
            assert pytest.approx(
                0.127324, rel=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for t, c in FB0D.mass_solids_constraint.items():
            assert pytest.approx(
                0.00127324, rel=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for t, c in FB0D.sum_component_constraint.items():
            assert pytest.approx(
                10.00000, rel=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for t, c in FB0D.solids_energy_holdup_constraints.items():
            assert pytest.approx(
                1.27324e-7, rel=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for t, c in FB0D.solids_energy_accumulation_constraints.items():
            assert pytest.approx(
                1.27324e-7, rel=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, iron_oc):
        optarg = {"tol": 1e-6}

        initialization_tester(iron_oc, optarg=optarg)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize_by_time(self, iron_oc):

        optarg = {"tol": 1e-6}
        solver = get_solver("ipopt", optarg)  # create solver

        initialize_by_time_element(iron_oc.fs, iron_oc.fs.time, solver=solver)

        assert degrees_of_freedom(iron_oc) == 0

        # Assert that model is still fixed and deactivated as expected
        assert iron_oc.fs.unit.solids[iron_oc.fs.time.first()].particle_porosity.fixed

        for t in iron_oc.fs.time:
            if t != iron_oc.fs.time.first():
                assert not iron_oc.fs.unit.solids[t].particle_porosity.fixed
            assert iron_oc.fs.unit.solids[t].density_particle_constraint.active
            assert iron_oc.fs.unit.solids[t].density_skeletal_constraint.active

        # Assert that constraints are feasible after initialization
        for con in iron_oc.fs.component_data_objects(Constraint, active=True):
            assert value(con.body) - value(con.upper) < 1e-4
            assert value(con.lower) - value(con.body) < 1e-4

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, iron_oc):
        results = solver.solve(iron_oc)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, iron_oc):
        assert (
            pytest.approx(1798.8532, rel=1e-5)
            == iron_oc.fs.unit.mass_solids[3600].value
        )
        assert (
            pytest.approx(0.195476, rel=1e-5)
            == iron_oc.fs.unit.solids[3600].particle_porosity.value
        )
        assert (
            pytest.approx(0.558299, rel=1e-5)
            == iron_oc.fs.unit.solids[3600].mass_frac_comp["Al2O3"].value
        )
        assert (
            pytest.approx(0.00511114, rel=1e-5)
            == iron_oc.fs.unit.solids[3600].mass_frac_comp["Fe2O3"].value
        )
        assert (
            pytest.approx(0.43659, rel=1e-5)
            == iron_oc.fs.unit.solids[3600].mass_frac_comp["Fe3O4"].value
        )
        assert (
            pytest.approx(1255.59, rel=1e-5)
            == iron_oc.fs.unit.solids[3600].temperature.value
        )

    @pytest.mark.component
    def test_units_consistent(self, iron_oc):
        assert_units_consistent(iron_oc)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, iron_oc):
        # Conservation of material check
        dt, T = 36, 3600  # size of time interval and time domain
        # first, check if Final Mass - Initial Mass = Total Amount Reacted
        mbal_solid_total = value(
            (iron_oc.fs.unit.mass_solids[T] - iron_oc.fs.unit.mass_solids[0])
            - (
                sum(
                    iron_oc.fs.unit.reactions[t].reaction_rate[r]
                    * iron_oc.fs.unit.config.reaction_package.rate_reaction_stoichiometry[
                        r, "Sol", j
                    ]
                    * iron_oc.fs.unit.solids[t]._params.mw_comp[j]
                    * iron_oc.fs.unit.volume_solid[t]
                    * dt
                    for t in list(iron_oc.fs.time)[1:]
                    for j in ("Al2O3", "Fe2O3", "Fe3O4")
                    for r in iron_oc.fs.unit.config.reaction_package.rate_reaction_idx
                )
            )
        )
        assert abs(mbal_solid_total) <= 1e-2

        # second, check if Change in Mass = Amount Reacted at each time t
        mbal_solid = []
        for t in list(iron_oc.fs.time)[1:]:  # skip t = 0
            mbal_solid.append(
                value(
                    (
                        iron_oc.fs.unit.mass_solids[t]
                        - iron_oc.fs.unit.mass_solids[t - dt]
                    )
                    - (
                        sum(
                            iron_oc.fs.unit.reactions[t].reaction_rate[r]
                            * iron_oc.fs.unit.config.reaction_package.rate_reaction_stoichiometry[
                                r, "Sol", j
                            ]
                            * iron_oc.fs.unit.solids[t]._params.mw_comp[j]
                            * iron_oc.fs.unit.volume_solid[t]
                            * dt
                            for j in ("Al2O3", "Fe2O3", "Fe3O4")
                            for r in iron_oc.fs.unit.config.reaction_package.rate_reaction_idx
                        )
                    )
                )
            )
        for val in mbal_solid:
            assert abs(val) <= 1e-2

        # Reaction stoichiometric ratio (gas is fixed, can only check solids)
        # Overall reducer reactions for methane combustion:
        # CH4 + 12Fe2O3 => 8Fe3O4 + CO2 + 2H2O
        fe2o3_reacted = value(
            (
                iron_oc.fs.unit.mass_solids[0]
                * iron_oc.fs.unit.solids[0].mass_frac_comp["Fe2O3"]
                - iron_oc.fs.unit.mass_solids[T]
                * iron_oc.fs.unit.solids[T].mass_frac_comp["Fe2O3"]
            )
            / iron_oc.fs.unit.solids[0]._params.mw_comp["Fe2O3"]
        )
        fe3o4_produced = value(
            (
                iron_oc.fs.unit.mass_solids[T]
                * iron_oc.fs.unit.solids[T].mass_frac_comp["Fe3O4"]
                - iron_oc.fs.unit.mass_solids[0]
                * iron_oc.fs.unit.solids[0].mass_frac_comp["Fe3O4"]
            )
            / iron_oc.fs.unit.solids[0]._params.mw_comp["Fe3O4"]
        )
        stoichiometric_ratio = fe2o3_reacted / fe3o4_produced
        assert pytest.approx(1.5, rel=1e-5) == stoichiometric_ratio

        # Conservation of energy check
        # first, check if Final Energy - Initial Energy = Reaction Enthalpy
        ebal_solid_total = value(
            (
                iron_oc.fs.unit.mass_solids[T] * iron_oc.fs.unit.solids[T].enth_mass
                - iron_oc.fs.unit.mass_solids[0] * iron_oc.fs.unit.solids[0].enth_mass
            )
            - (
                -sum(
                    iron_oc.fs.unit.reactions[t].reaction_rate[r]
                    * iron_oc.fs.unit.reactions[t].dh_rxn[r]
                    * iron_oc.fs.unit.volume_solid[t]
                    * dt
                    for t in list(iron_oc.fs.time)[1:]
                    for r in iron_oc.fs.unit.config.reaction_package.rate_reaction_idx
                )
            )
        )
        assert abs(ebal_solid_total) <= 1e-2

        # second, check if Change in Energy = Change in Holdup at each time t
        ebal_solid = []
        for t in list(iron_oc.fs.time)[1:]:  # skip t = 0
            ebal_solid.append(
                value(
                    (
                        iron_oc.fs.unit.mass_solids[t]
                        * iron_oc.fs.unit.solids[t].enth_mass
                        - iron_oc.fs.unit.mass_solids[t - dt]
                        * iron_oc.fs.unit.solids[t - dt].enth_mass
                    )
                    - (
                        -sum(
                            iron_oc.fs.unit.reactions[t].reaction_rate[r]
                            * iron_oc.fs.unit.reactions[t].dh_rxn[r]
                            * iron_oc.fs.unit.volume_solid[t]
                            * dt
                            for r in iron_oc.fs.unit.config.reaction_package.rate_reaction_idx
                        )
                    )
                )
            )
        for val in ebal_solid:
            assert abs(val) <= 1e-2

    @pytest.mark.ui
    @pytest.mark.unit
    def test_report(self, iron_oc):
        # no inlet/outlet ports (no stream table), test block reports instead
        iron_oc.fs.unit.gas.report()
        iron_oc.fs.unit.solids.report()

    @pytest.mark.component
    def test_initialization_error(self, iron_oc):
        iron_oc.fs.unit.volume_bed.fix(1)

        with pytest.raises(InitializationError):
            optarg = {"tol": 1e-6}
            iron_oc.fs.unit.initialize(optarg=optarg)


# -----------------------------------------------------------------------------
class TestIronOC_EnergyBalanceType(object):
    @pytest.fixture(scope="class")
    def iron_oc(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=True, time_set=[0, 3600], time_units=pyunits.s)

        m.fs.gas_props = GasPhaseParameterBlock()
        m.fs.solid_props = SolidPhaseParameterBlock()
        m.fs.solid_rxns = HeteroReactionParameterBlock(
            solid_property_package=m.fs.solid_props, gas_property_package=m.fs.gas_props
        )

        m.fs.unit = FixedBed0D(
            energy_balance_type=EnergyBalanceType.none,
            gas_property_package=m.fs.gas_props,
            solid_property_package=m.fs.solid_props,
            reaction_package=m.fs.solid_rxns,
        )

        # Discretize time domain
        m.discretizer = TransformationFactory("dae.finite_difference")
        m.discretizer.apply_to(m, nfe=100, wrt=m.fs.time, scheme="BACKWARD")

        # Set reactor design conditions
        m.fs.unit.bed_diameter.fix(1)  # diameter of the TGA reactor [m]
        m.fs.unit.bed_height.fix(1)  # height of solids in the TGA reactor [m]

        # Set initial conditions of the solid phase
        m.fs.unit.solids[0].temperature.fix(1273.15)
        m.fs.unit.solids[0].particle_porosity.fix(0.20)
        m.fs.unit.solids[0].mass_frac_comp["Fe2O3"].fix(0.45)
        m.fs.unit.solids[0].mass_frac_comp["Fe3O4"].fix(0)
        m.fs.unit.solids[0].mass_frac_comp["Al2O3"].fix(0.55)

        # Set conditions of the gas phase (this is all fixed as gas side
        # assumption is excess gas flowrate which means all state variables
        # remain unchanged)
        for t in m.fs.time:
            m.fs.unit.gas[t].temperature.fix(1273.15)
            m.fs.unit.gas[t].pressure.fix(1.01325e5)  # 1atm
            m.fs.unit.gas[t].mole_frac_comp["CO2"].fix(0.4)
            m.fs.unit.gas[t].mole_frac_comp["H2O"].fix(0.5)
            m.fs.unit.gas[t].mole_frac_comp["CH4"].fix(0.1)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, iron_oc):
        assert isinstance(iron_oc.fs.unit.volume_bed_constraint, Constraint)
        assert isinstance(iron_oc.fs.unit.volume_solid_constraint, Constraint)
        assert isinstance(
            iron_oc.fs.unit.solids_material_holdup_constraints, Constraint
        )
        assert isinstance(
            iron_oc.fs.unit.solids_material_accumulation_constraints, Constraint
        )
        assert isinstance(iron_oc.fs.unit.mass_solids_constraint, Constraint)
        assert isinstance(iron_oc.fs.unit.sum_component_constraint, Constraint)

        assert number_variables(iron_oc) == 3044
        assert number_total_constraints(iron_oc) == 2319
        assert number_unused_variables(iron_oc) == 206
        print(unused_variables_set(iron_oc))

    @pytest.mark.unit
    def test_dof(self, iron_oc):
        assert degrees_of_freedom(iron_oc) == 0

    @pytest.fixture(scope="class")
    def iron_oc_unscaled(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=True, time_set=[0, 3600], time_units=pyunits.s)

        m.fs.gas_props = GasPhaseParameterBlock()
        m.fs.solid_props = SolidPhaseParameterBlock()
        m.fs.solid_rxns = HeteroReactionParameterBlock(
            solid_property_package=m.fs.solid_props, gas_property_package=m.fs.gas_props
        )

        m.fs.unit = FixedBed0D(
            energy_balance_type=EnergyBalanceType.none,
            gas_property_package=m.fs.gas_props,
            solid_property_package=m.fs.solid_props,
            reaction_package=m.fs.solid_rxns,
        )

        # Discretize time domain
        m.discretizer = TransformationFactory("dae.finite_difference")
        m.discretizer.apply_to(m, nfe=100, wrt=m.fs.time, scheme="BACKWARD")

        # Set reactor design conditions
        m.fs.unit.bed_diameter.fix(1)  # diameter of the TGA reactor [m]
        m.fs.unit.bed_height.fix(1)  # height of solids in the TGA reactor [m]

        # Set initial conditions of the solid phase
        m.fs.unit.solids[0].temperature.fix(1273.15)
        m.fs.unit.solids[0].particle_porosity.fix(0.20)
        m.fs.unit.solids[0].mass_frac_comp["Fe2O3"].fix(0.45)
        m.fs.unit.solids[0].mass_frac_comp["Fe3O4"].fix(0)
        m.fs.unit.solids[0].mass_frac_comp["Al2O3"].fix(0.55)

        # Set conditions of the gas phase (this is all fixed as gas side
        # assumption is excess gas flowrate which means all state variables
        # remain unchanged)
        for t in m.fs.time:
            m.fs.unit.gas[t].temperature.fix(1273.15)
            m.fs.unit.gas[t].pressure.fix(1.01325e5)  # 1atm
            m.fs.unit.gas[t].mole_frac_comp["CO2"].fix(0.4)
            m.fs.unit.gas[t].mole_frac_comp["H2O"].fix(0.5)
            m.fs.unit.gas[t].mole_frac_comp["CH4"].fix(0.1)

        return m

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize_unscaled(self, iron_oc_unscaled):
        optarg = {"tol": 1e-6}

        initialization_tester(iron_oc_unscaled, optarg=optarg)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize_by_time_unscaled(self, iron_oc_unscaled):

        optarg = {"tol": 1e-6}
        solver = get_solver("ipopt", optarg)  # create solver

        initialize_by_time_element(
            iron_oc_unscaled.fs, iron_oc_unscaled.fs.time, solver=solver
        )

        assert degrees_of_freedom(iron_oc_unscaled) == 0

        # Assert that model is still fixed and deactivated as expected
        assert iron_oc_unscaled.fs.unit.solids[
            iron_oc_unscaled.fs.time.first()
        ].particle_porosity.fixed

        for t in iron_oc_unscaled.fs.time:
            if t != iron_oc_unscaled.fs.time.first():
                assert not iron_oc_unscaled.fs.unit.solids[t].particle_porosity.fixed
            assert iron_oc_unscaled.fs.unit.solids[t].density_particle_constraint.active
            assert iron_oc_unscaled.fs.unit.solids[t].density_skeletal_constraint.active

        # Assert that constraints are feasible after initialization
        for con in iron_oc_unscaled.fs.component_data_objects(Constraint, active=True):
            assert value(con.body) - value(con.upper) < 1e-4
            assert value(con.lower) - value(con.body) < 1e-4

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve_unscaled(self, iron_oc_unscaled):
        results = solver.solve(iron_oc_unscaled)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.component
    def test_scaling(self, iron_oc):

        # Set scaling gas phase for state variables
        iron_oc.fs.gas_props.set_default_scaling("flow_mol", 1e-3)
        iron_oc.fs.gas_props.set_default_scaling("pressure", 1e-5)
        iron_oc.fs.gas_props.set_default_scaling("temperature", 1e-2)
        for comp in iron_oc.fs.gas_props.component_list:
            iron_oc.fs.gas_props.set_default_scaling("mole_frac_comp", 1e1, index=comp)
        # Set scaling for gas phase thermophysical and transport properties
        iron_oc.fs.gas_props.set_default_scaling("enth_mol", 1e-6)
        iron_oc.fs.gas_props.set_default_scaling("enth_mol_comp", 1e-6)
        iron_oc.fs.gas_props.set_default_scaling("cp_mol", 1e-6)
        iron_oc.fs.gas_props.set_default_scaling("cp_mol_comp", 1e-6)
        iron_oc.fs.gas_props.set_default_scaling("cp_mass", 1e-6)
        iron_oc.fs.gas_props.set_default_scaling("entr_mol", 1e-2)
        iron_oc.fs.gas_props.set_default_scaling("entr_mol_phase", 1e-2)
        iron_oc.fs.gas_props.set_default_scaling("dens_mol", 1)
        iron_oc.fs.gas_props.set_default_scaling("dens_mol_comp", 1)
        iron_oc.fs.gas_props.set_default_scaling("dens_mass", 1e2)
        iron_oc.fs.gas_props.set_default_scaling("visc_d_comp", 1e4)
        iron_oc.fs.gas_props.set_default_scaling("diffusion_comp", 1e5)
        iron_oc.fs.gas_props.set_default_scaling("therm_cond_comp", 1e2)
        iron_oc.fs.gas_props.set_default_scaling("visc_d", 1e5)
        iron_oc.fs.gas_props.set_default_scaling("therm_cond", 1e0)
        iron_oc.fs.gas_props.set_default_scaling("mw", 1e2)

        # Set scaling for solid phase state variables
        iron_oc.fs.solid_props.set_default_scaling("flow_mass", 1e-3)
        iron_oc.fs.solid_props.set_default_scaling("particle_porosity", 1e2)
        iron_oc.fs.solid_props.set_default_scaling("temperature", 1e-2)
        for comp in iron_oc.fs.solid_props.component_list:
            iron_oc.fs.solid_props.set_default_scaling(
                "mass_frac_comp", 1e1, index=comp
            )

        # Set scaling for solid phase thermophysical and transport properties
        iron_oc.fs.solid_props.set_default_scaling("enth_mass", 1e-6)
        iron_oc.fs.solid_props.set_default_scaling("enth_mol_comp", 1e-6)
        iron_oc.fs.solid_props.set_default_scaling("cp_mol_comp", 1e-6)
        iron_oc.fs.solid_props.set_default_scaling("cp_mass", 1e-6)
        iron_oc.fs.solid_props.set_default_scaling("dens_mass_particle", 1e-2)
        iron_oc.fs.solid_props.set_default_scaling("dens_mass_skeletal", 1e-2)

        FB0D = iron_oc.fs.unit  # alias to keep test lines short

        # Calculate scaling factors
        iscale.calculate_scaling_factors(FB0D)

        assert pytest.approx(0.10000, rel=1e-5) == iscale.get_scaling_factor(
            FB0D.bed_diameter
        )
        assert pytest.approx(0.10000, rel=1e-5) == iscale.get_scaling_factor(
            FB0D.bed_height
        )
        assert pytest.approx(0.127324, rel=1e-5) == iscale.get_scaling_factor(
            FB0D.volume_bed
        )

        for c in FB0D.volume_bed_constraint.values():
            assert pytest.approx(
                0.127324, rel=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for t, c in FB0D.volume_solid_constraint.items():
            assert pytest.approx(
                0.127324, rel=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, j), c in FB0D.solids_material_holdup_constraints.items():
            assert pytest.approx(
                0.00127324, rel=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, j), c in FB0D.solids_material_accumulation_constraints.items():
            assert pytest.approx(
                0.127324, rel=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for t, c in FB0D.mass_solids_constraint.items():
            assert pytest.approx(
                0.00127324, rel=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for t, c in FB0D.sum_component_constraint.items():
            assert pytest.approx(
                10.00000, rel=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for t, c in FB0D.isothermal_solid_phase.items():
            assert pytest.approx(
                0.01000, rel=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, iron_oc):
        optarg = {"tol": 1e-6}

        initialization_tester(iron_oc, optarg=optarg)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize_by_time(self, iron_oc):

        optarg = {"tol": 1e-6}
        solver = get_solver("ipopt", optarg)  # create solver

        initialize_by_time_element(iron_oc.fs, iron_oc.fs.time, solver=solver)

        assert degrees_of_freedom(iron_oc) == 0

        # Assert that model is still fixed and deactivated as expected
        assert iron_oc.fs.unit.solids[iron_oc.fs.time.first()].particle_porosity.fixed

        for t in iron_oc.fs.time:
            if t != iron_oc.fs.time.first():
                assert not iron_oc.fs.unit.solids[t].particle_porosity.fixed
            assert iron_oc.fs.unit.solids[t].density_particle_constraint.active
            assert iron_oc.fs.unit.solids[t].density_skeletal_constraint.active

        # Assert that constraints are feasible after initialization
        for con in iron_oc.fs.component_data_objects(Constraint, active=True):
            assert value(con.body) - value(con.upper) < 1e-4
            assert value(con.lower) - value(con.body) < 1e-4

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, iron_oc):
        results = solver.solve(iron_oc)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, iron_oc):  # need to update these values
        assert (
            pytest.approx(1798.8281, rel=1e-5)
            == iron_oc.fs.unit.mass_solids[3600].value
        )
        assert (
            pytest.approx(0.195472, rel=1e-5)
            == iron_oc.fs.unit.solids[3600].particle_porosity.value
        )
        assert (
            pytest.approx(0.558307, rel=1e-5)
            == iron_oc.fs.unit.solids[3600].mass_frac_comp["Al2O3"].value
        )
        assert (
            pytest.approx(0.00469410, rel=1e-5)
            == iron_oc.fs.unit.solids[3600].mass_frac_comp["Fe2O3"].value
        )
        assert (
            pytest.approx(0.4370, rel=1e-5)
            == iron_oc.fs.unit.solids[3600].mass_frac_comp["Fe3O4"].value
        )
        assert (
            pytest.approx(1273.15, rel=1e-5)
            == iron_oc.fs.unit.solids[3600].temperature.value
        )

    @pytest.mark.component
    def test_units_consistent(self, iron_oc):
        assert_units_consistent(iron_oc)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, iron_oc):
        # Conservation of material check
        dt, T = 36, 3600  # size of time interval and time domain
        # first, check if Final Mass - Initial Mass = Total Amount Reacted
        mbal_solid_total = value(
            (iron_oc.fs.unit.mass_solids[T] - iron_oc.fs.unit.mass_solids[0])
            - (
                sum(
                    iron_oc.fs.unit.reactions[t].reaction_rate[r]
                    * iron_oc.fs.unit.config.reaction_package.rate_reaction_stoichiometry[
                        r, "Sol", j
                    ]
                    * iron_oc.fs.unit.solids[t]._params.mw_comp[j]
                    * iron_oc.fs.unit.volume_solid[t]
                    * dt
                    for t in list(iron_oc.fs.time)[1:]
                    for j in ("Al2O3", "Fe2O3", "Fe3O4")
                    for r in iron_oc.fs.unit.config.reaction_package.rate_reaction_idx
                )
            )
        )
        assert abs(mbal_solid_total) <= 1e-2

        # second, check if Change in Mass = Amount Reacted at each time t
        mbal_solid = []
        for t in list(iron_oc.fs.time)[1:]:  # skip t = 0
            mbal_solid.append(
                value(
                    (
                        iron_oc.fs.unit.mass_solids[t]
                        - iron_oc.fs.unit.mass_solids[t - dt]
                    )
                    - (
                        sum(
                            iron_oc.fs.unit.reactions[t].reaction_rate[r]
                            * iron_oc.fs.unit.config.reaction_package.rate_reaction_stoichiometry[
                                r, "Sol", j
                            ]
                            * iron_oc.fs.unit.solids[t]._params.mw_comp[j]
                            * iron_oc.fs.unit.volume_solid[t]
                            * dt
                            for j in ("Al2O3", "Fe2O3", "Fe3O4")
                            for r in iron_oc.fs.unit.config.reaction_package.rate_reaction_idx
                        )
                    )
                )
            )
        for val in mbal_solid:
            assert abs(val) <= 1e-2

        # Reaction stoichiometric ratio (gas is fixed, can only check solids)
        # Overall reducer reactions for methane combustion:
        # CH4 + 12Fe2O3 => 8Fe3O4 + CO2 + 2H2O
        fe2o3_reacted = value(
            (
                iron_oc.fs.unit.mass_solids[0]
                * iron_oc.fs.unit.solids[0].mass_frac_comp["Fe2O3"]
                - iron_oc.fs.unit.mass_solids[T]
                * iron_oc.fs.unit.solids[T].mass_frac_comp["Fe2O3"]
            )
            / iron_oc.fs.unit.solids[0]._params.mw_comp["Fe2O3"]
        )
        fe3o4_produced = value(
            (
                iron_oc.fs.unit.mass_solids[T]
                * iron_oc.fs.unit.solids[T].mass_frac_comp["Fe3O4"]
                - iron_oc.fs.unit.mass_solids[0]
                * iron_oc.fs.unit.solids[0].mass_frac_comp["Fe3O4"]
            )
            / iron_oc.fs.unit.solids[0]._params.mw_comp["Fe3O4"]
        )
        stoichiometric_ratio = fe2o3_reacted / fe3o4_produced
        assert pytest.approx(1.5, rel=1e-5) == stoichiometric_ratio

    @pytest.mark.ui
    @pytest.mark.unit
    def test_report(self, iron_oc):
        # no inlet/outlet ports (no stream table), test block reports instead
        iron_oc.fs.unit.gas.report()
        iron_oc.fs.unit.solids.report()

    @pytest.mark.component
    def test_initialization_error(self, iron_oc):
        iron_oc.fs.unit.volume_bed.fix(1)

        with pytest.raises(InitializationError):
            optarg = {"tol": 1e-6}
            iron_oc.fs.unit.initialize(optarg=optarg)
