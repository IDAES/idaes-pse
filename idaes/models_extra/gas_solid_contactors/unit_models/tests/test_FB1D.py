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
Tests for the 1D fixed bed module

Author: Chinedu Okoli
"""

import pytest

from pyomo.environ import (
    assert_optimal_termination,
    ConcreteModel,
    TransformationFactory,
    value,
    units as pyunits,
    Constraint,
    Var,
)
from pyomo.dae import ContinuousSet, Integral
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from pyomo.util.check_units import assert_units_consistent
from pyomo.common.config import ConfigBlock
from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
)
import pyomo.common.unittest as unittest

import idaes
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.testing import initialization_tester
from idaes.core.util.initialization import initialize_by_time_element
from idaes.core.util import scaling as iscale
from idaes.core.solvers import get_solver
from idaes.core.util.exceptions import ConfigurationError, InitializationError
from idaes.core.util.performance import PerformanceBaseClass

# Import FixedBed1D unit model
from idaes.models_extra.gas_solid_contactors.unit_models.fixed_bed_1D import FixedBed1D

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
    m.fs = FlowsheetBlock(dynamic=True, time_set=[0, 360], time_units=pyunits.s)

    # Set up thermo props and reaction props
    m.fs.gas_props = GasPhaseParameterBlock()
    m.fs.solid_props = SolidPhaseParameterBlock()
    m.fs.solid_rxns = HeteroReactionParameterBlock(
        solid_property_package=m.fs.solid_props, gas_property_package=m.fs.gas_props
    )

    m.fs.unit = FixedBed1D(
        gas_phase_config={"property_package": m.fs.gas_props},
        solid_phase_config={
            "property_package": m.fs.solid_props,
            "reaction_package": m.fs.solid_rxns,
        },
    )

    # Check unit config arguments
    assert len(m.fs.unit.config) == 15

    assert m.fs.unit.config.dynamic is True
    assert m.fs.unit.config.has_holdup is True
    assert isinstance(m.fs.unit.config.gas_phase_config, ConfigBlock)
    assert isinstance(m.fs.unit.config.solid_phase_config, ConfigBlock)

    assert m.fs.unit.config.finite_elements == 10
    assert m.fs.unit.config.length_domain_set == [0.0, 1.0]
    assert m.fs.unit.config.transformation_method == "dae.finite_difference"
    assert m.fs.unit.config.transformation_scheme == "BACKWARD"
    assert m.fs.unit.config.collocation_points == 3
    assert m.fs.unit.config.flow_type == "forward_flow"
    assert m.fs.unit.config.material_balance_type == MaterialBalanceType.componentTotal
    assert m.fs.unit.config.energy_balance_type == EnergyBalanceType.enthalpyTotal
    assert m.fs.unit.config.momentum_balance_type == MomentumBalanceType.pressureTotal
    assert m.fs.unit.config.has_pressure_change is True

    # Check gas phase config arguments
    assert len(m.fs.unit.config.gas_phase_config) == 7
    assert m.fs.unit.config.gas_phase_config.has_equilibrium_reactions is False
    assert m.fs.unit.config.gas_phase_config.property_package is m.fs.gas_props
    assert bool(m.fs.unit.config.gas_phase_config.property_package_args) is False
    assert m.fs.unit.config.gas_phase_config.reaction_package is None

    # Check solid phase config arguments
    assert len(m.fs.unit.config.solid_phase_config) == 7
    assert m.fs.unit.config.solid_phase_config.has_equilibrium_reactions is False
    assert m.fs.unit.config.solid_phase_config.property_package is m.fs.solid_props
    assert bool(m.fs.unit.config.solid_phase_config.property_package_args) is False
    assert m.fs.unit.config.solid_phase_config.reaction_package is m.fs.solid_rxns


@pytest.mark.unit
def test_config_validation():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_set=[0, 360], time_units=pyunits.s)

    # Set up thermo props and reaction props
    m.fs.gas_props = GasPhaseParameterBlock()
    m.fs.solid_props = SolidPhaseParameterBlock()
    m.fs.solid_rxns = HeteroReactionParameterBlock(
        solid_property_package=m.fs.solid_props, gas_property_package=m.fs.gas_props
    )

    m.fs.unit = FixedBed1D(
        gas_phase_config={"property_package": m.fs.gas_props},
        solid_phase_config={
            "property_package": m.fs.solid_props,
            "reaction_package": m.fs.solid_rxns,
        },
    )

    with pytest.raises(ConfigurationError):
        m.fs.unit = FixedBed1D(
            flow_type="reverse_flow",
            transformation_method="dae.collocation",
            gas_phase_config={"property_package": m.fs.gas_props},
            solid_phase_config={
                "property_package": m.fs.solid_props,
                "reaction_package": m.fs.solid_rxns,
            },
        )

    with pytest.raises(ConfigurationError):
        m.fs.unit = FixedBed1D(
            flow_type="reverse_flow",
            transformation_method="dae.finite_difference",
            transformation_scheme="BACKWARD",
            gas_phase_config={"property_package": m.fs.gas_props},
            solid_phase_config={
                "property_package": m.fs.solid_props,
                "reaction_package": m.fs.solid_rxns,
            },
        )

    with pytest.raises(ConfigurationError):
        m.fs.unit = FixedBed1D(
            flow_type="forward_flow",
            transformation_method="dae.finite_difference",
            transformation_scheme="FORWARD",
            gas_phase_config={"property_package": m.fs.gas_props},
            solid_phase_config={
                "property_package": m.fs.solid_props,
                "reaction_package": m.fs.solid_rxns,
            },
        )

    with pytest.raises(ConfigurationError):
        m.fs.unit = FixedBed1D(
            transformation_method="dae.collocation",
            transformation_scheme="BACKWARD",
            gas_phase_config={"property_package": m.fs.gas_props},
            solid_phase_config={
                "property_package": m.fs.solid_props,
                "reaction_package": m.fs.solid_rxns,
            },
        )

    with pytest.raises(ConfigurationError):
        m.fs.unit = FixedBed1D(
            transformation_method="dae.collocation",
            transformation_scheme="FORWARD",
            gas_phase_config={"property_package": m.fs.gas_props},
            solid_phase_config={
                "property_package": m.fs.solid_props,
                "reaction_package": m.fs.solid_rxns,
            },
        )


# -----------------------------------------------------------------------------
def build_model():
    m = ConcreteModel()
    horizon = 360
    m.fs = FlowsheetBlock(dynamic=True, time_set=[0, horizon], time_units=pyunits.s)

    m.fs.gas_props = GasPhaseParameterBlock()
    m.fs.solid_props = SolidPhaseParameterBlock()
    m.fs.solid_rxns = HeteroReactionParameterBlock(
        solid_property_package=m.fs.solid_props, gas_property_package=m.fs.gas_props
    )

    m.fs.unit = FixedBed1D(
        gas_phase_config={"property_package": m.fs.gas_props},
        solid_phase_config={
            "property_package": m.fs.solid_props,
            "reaction_package": m.fs.solid_rxns,
        },
    )

    # Discretize time domain
    t_element_size = 120  # s
    ntfe = int(horizon / t_element_size)
    m.discretizer = TransformationFactory("dae.finite_difference")
    m.discretizer.apply_to(m, nfe=ntfe, wrt=m.fs.time, scheme="BACKWARD")

    # Set reactor design conditions
    m.fs.unit.bed_diameter.fix(9)  # diameter of the reactor [m]
    m.fs.unit.bed_height.fix(5)  # height of solids in the reactor [m]

    # Fix inlet port variables (boundary conditions) for gas
    for t in m.fs.time:
        m.fs.unit.gas_inlet.flow_mol[t].fix(20)  # mol/s
        m.fs.unit.gas_inlet.temperature[t].fix(723.15)  # K
        m.fs.unit.gas_inlet.pressure[t].fix(1.565e5)  # Pa
        m.fs.unit.gas_inlet.mole_frac_comp[t, "CO2"].fix(0.02499)
        m.fs.unit.gas_inlet.mole_frac_comp[t, "H2O"].fix(0.00001)
        m.fs.unit.gas_inlet.mole_frac_comp[t, "CH4"].fix(0.975)

    # Specify gas phase initial conditions
    t0 = m.fs.time.first()
    for x in m.fs.unit.length_domain:
        if x != m.fs.unit.length_domain.first():  # skip inlet as it's already fixed
            m.fs.unit.gas_phase.properties[t0, x].mole_frac_comp["CO2"].fix(0.02499)
            m.fs.unit.gas_phase.properties[t0, x].mole_frac_comp["H2O"].fix(0.00001)
            m.fs.unit.gas_phase.properties[t0, x].mole_frac_comp["CH4"].fix(0.975)
            m.fs.unit.gas_phase.properties[t0, x].temperature.fix(1273.15)  # K
            m.fs.unit.gas_phase.properties[t0, x].flow_mol.fix(20)  # mol/s

    # Specify solid phase initial conditions
    t0 = m.fs.time.first()
    for x in m.fs.unit.length_domain:
        m.fs.unit.solid_properties[t0, x].temperature.fix(1273.15)
        m.fs.unit.solid_properties[t0, x].particle_porosity.fix(0.20)
        m.fs.unit.solid_properties[t0, x].mass_frac_comp["Fe2O3"].fix(0.45)
        m.fs.unit.solid_properties[t0, x].mass_frac_comp["Fe3O4"].fix(0)
        m.fs.unit.solid_properties[t0, x].mass_frac_comp["Al2O3"].fix(0.55)

    return m


def set_scaling(iron_oc):
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
    iron_oc.fs.solid_props.set_default_scaling("particle_porosity", 1e2)
    iron_oc.fs.solid_props.set_default_scaling("temperature", 1e-2)
    for comp in iron_oc.fs.solid_props.component_list:
        iron_oc.fs.solid_props.set_default_scaling("mass_frac_comp", 1e1, index=comp)

    # Set scaling for solid phase thermophysical and transport properties
    iron_oc.fs.solid_props.set_default_scaling("enth_mass", 1e-6)
    iron_oc.fs.solid_props.set_default_scaling("enth_mol_comp", 1e-6)
    iron_oc.fs.solid_props.set_default_scaling("cp_mol_comp", 1e-6)
    iron_oc.fs.solid_props.set_default_scaling("cp_mass", 1e-6)
    iron_oc.fs.solid_props.set_default_scaling("dens_mass_particle", 1e-2)
    iron_oc.fs.solid_props.set_default_scaling("dens_mass_skeletal", 1e-2)

    FB1D = iron_oc.fs.unit  # alias to keep test lines short

    # Calculate scaling factors
    iscale.calculate_scaling_factors(FB1D)


@pytest.mark.performance
class Test_FixedBed1D_Performance(PerformanceBaseClass, unittest.TestCase):
    # TODO: Remove this once Pyomo bug in DAE units is fixed
    TEST_UNITS = False

    def build_model(self):
        model = build_model()
        set_scaling(model)
        return model

    def initialize_model(self, model):
        with idaes.temporary_config_ctx():
            optarg = {
                "tol": 1e-5,
                "bound_push": 1e-22,
                "nlp_scaling_method": "user-scaling",
            }
            model.fs.unit.initialize(optarg=optarg)

    def solve_model(self, model):
        with idaes.temporary_config_ctx():
            solver.options = {
                "tol": 1e-5,
                "bound_push": 1e-22,
                "nlp_scaling_method": "user-scaling",
            }
            results = solver.solve(model)

            # Check for optimal solution
            assert_optimal_termination(results)


class TestIronOC(object):
    @pytest.fixture(scope="class")
    def iron_oc(self):
        return build_model()

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, iron_oc):
        assert hasattr(iron_oc.fs.unit, "gas_inlet")
        assert len(iron_oc.fs.unit.gas_inlet.vars) == 4
        assert isinstance(iron_oc.fs.unit.gas_inlet.flow_mol, Var)
        assert isinstance(iron_oc.fs.unit.gas_inlet.mole_frac_comp, Var)
        assert isinstance(iron_oc.fs.unit.gas_inlet.temperature, Var)
        assert isinstance(iron_oc.fs.unit.gas_inlet.pressure, Var)

        assert hasattr(iron_oc.fs.unit, "gas_outlet")
        assert len(iron_oc.fs.unit.gas_outlet.vars) == 4
        assert isinstance(iron_oc.fs.unit.gas_outlet.flow_mol, Var)
        assert isinstance(iron_oc.fs.unit.gas_outlet.mole_frac_comp, Var)
        assert isinstance(iron_oc.fs.unit.gas_outlet.temperature, Var)
        assert isinstance(iron_oc.fs.unit.gas_outlet.pressure, Var)

        assert isinstance(iron_oc.fs.unit.bed_area_eqn, Constraint)
        assert isinstance(iron_oc.fs.unit.gas_phase_area_constraint, Constraint)
        assert isinstance(iron_oc.fs.unit.solid_phase_area_constraint, Constraint)
        assert isinstance(iron_oc.fs.unit.gas_super_vel, Constraint)
        assert isinstance(iron_oc.fs.unit.gas_phase_config_pressure_drop, Constraint)
        assert isinstance(iron_oc.fs.unit.gas_comp_hetero_rxn, Constraint)
        assert isinstance(iron_oc.fs.unit.solid_material_holdup_calculation, Constraint)
        assert isinstance(iron_oc.fs.unit.solid_material_balances, Constraint)
        assert isinstance(iron_oc.fs.unit.solid_sum_component_eqn, Constraint)
        assert isinstance(iron_oc.fs.unit.solid_phase_heat_transfer, Constraint)
        assert isinstance(iron_oc.fs.unit.solid_energy_holdup_calculation, Constraint)
        assert isinstance(iron_oc.fs.unit.solid_enthalpy_balances, Constraint)

        assert isinstance(iron_oc.fs.unit.reynolds_number_particle, Constraint)
        assert isinstance(iron_oc.fs.unit.prandtl_number, Constraint)
        assert isinstance(iron_oc.fs.unit.nusselt_number_particle, Constraint)
        assert isinstance(iron_oc.fs.unit.gas_solid_htc_eqn, Constraint)
        assert isinstance(iron_oc.fs.unit.gas_phase_heat_transfer, Constraint)

        assert number_variables(iron_oc) == 3402
        assert number_total_constraints(iron_oc) == 3202
        assert number_unused_variables(iron_oc) == 71

    @pytest.mark.unit
    def test_dof(self, iron_oc):
        assert degrees_of_freedom(iron_oc) == 0

    @pytest.mark.component
    def test_scaling(self, iron_oc):
        set_scaling(iron_oc)

        FB1D = iron_oc.fs.unit  # alias to keep test lines short

        assert pytest.approx(0.01111, rel=1e-3) == iscale.get_scaling_factor(
            FB1D.bed_diameter
        )
        assert pytest.approx(1, rel=1e-3) == iscale.get_scaling_factor(FB1D.bed_height)

        for c in FB1D.bed_area_eqn.values():
            assert pytest.approx(
                0.0015719, rel=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for t, c in FB1D.solid_phase_area_constraint.items():
            assert pytest.approx(
                0.0015719, rel=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x, j), c in FB1D.solid_material_holdup_calculation.items():
            assert pytest.approx(
                1.5719e-05, rel=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x, j), c in FB1D.solid_material_accumulation_disc_eq.items():
            assert pytest.approx(
                1e-07, rel=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x, j), c in FB1D.solid_material_balances.items():
            assert pytest.approx(
                0.00015719, rel=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in FB1D.solid_sum_component_eqn.items():
            assert pytest.approx(
                10.00000, rel=1e-5
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in FB1D.solid_energy_holdup_calculation.items():
            assert pytest.approx(
                1.5719e-11, rel=1e-8
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)
        for (t, x), c in FB1D.solid_energy_accumulation_disc_eq.items():
            assert pytest.approx(
                1e-11, rel=1e-8
            ) == iscale.get_constraint_transform_applied_scaling_factor(c)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, iron_oc):
        with idaes.temporary_config_ctx():
            optarg = {
                "tol": 1e-5,
                "bound_push": 1e-22,
                "nlp_scaling_method": "user-scaling",
            }

            initialization_tester(iron_oc, optarg=optarg)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_block_triangularization_initialization(self, iron_oc):
        iron_oc.fs.unit.block_triangularization_initialize(calc_var_kwds={"eps": 1e-5})

        # Calculate scaling factors
        iscale.calculate_scaling_factors(iron_oc)

        assert degrees_of_freedom(iron_oc) == 0

        # Assert that model is still fixed and deactivated as expected
        t0 = iron_oc.fs.time.first()
        x0 = iron_oc.fs.unit.length_domain.first()
        assert iron_oc.fs.unit.solid_properties[t0, x0].particle_porosity.fixed

        for t in iron_oc.fs.time:
            if t != iron_oc.fs.time.first():
                assert not iron_oc.fs.unit.solid_properties[
                    t, x0
                ].particle_porosity.fixed
            assert iron_oc.fs.unit.solid_properties[
                t, x0
            ].density_particle_constraint.active
            assert iron_oc.fs.unit.solid_properties[
                t, x0
            ].density_skeletal_constraint.active

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize_by_time(self, iron_oc):
        with idaes.temporary_config_ctx():
            optarg = {
                "tol": 1e-5,
                "bound_push": 1e-22,
                "nlp_scaling_method": "user-scaling",
            }
            solver = get_solver("ipopt", optarg)  # create solver

            initialize_by_time_element(iron_oc.fs, iron_oc.fs.time, solver=solver)

        assert degrees_of_freedom(iron_oc) == 0

        # Assert that model is still fixed and deactivated as expected
        t0 = iron_oc.fs.time.first()
        x0 = iron_oc.fs.unit.length_domain.first()
        assert iron_oc.fs.unit.solid_properties[t0, x0].particle_porosity.fixed

        for t in iron_oc.fs.time:
            if t != iron_oc.fs.time.first():
                assert not iron_oc.fs.unit.solid_properties[
                    t, x0
                ].particle_porosity.fixed
            assert iron_oc.fs.unit.solid_properties[
                t, x0
            ].density_particle_constraint.active
            assert iron_oc.fs.unit.solid_properties[
                t, x0
            ].density_skeletal_constraint.active

        # Assert that constraints are feasible after initialization
        for con in iron_oc.fs.component_data_objects(Constraint, active=True):
            assert value(con.body) - value(con.upper) < 1e-4
            assert value(con.lower) - value(con.body) < 1e-4

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, iron_oc):
        with idaes.temporary_config_ctx():
            solver.options = {
                "tol": 1e-5,
                "bound_push": 1e-22,
                "nlp_scaling_method": "user-scaling",
            }
            results = solver.solve(iron_oc)

            # Check for optimal solution
            assert_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, iron_oc):
        t0 = iron_oc.fs.time.first()
        x0 = iron_oc.fs.unit.length_domain.first()
        tf = iron_oc.fs.time.last()
        xf = iron_oc.fs.unit.length_domain.last()
        assert (
            pytest.approx(0.01209, abs=1e-2)
            == iron_oc.fs.unit.velocity_superficial_gas[t0, x0].value
        )
        assert (
            pytest.approx(0.06714, abs=1e-2)
            == iron_oc.fs.unit.velocity_superficial_gas[tf, xf].value
        )
        # Check the pressure drop that occurs across the bed
        assert (
            pytest.approx(146126.06009, abs=1e-2)
            == iron_oc.fs.unit.gas_outlet.pressure[tf].value
        )
        assert (
            pytest.approx(10373.9409, abs=1e-2)
            == iron_oc.fs.unit.gas_inlet.pressure[tf].value
            - iron_oc.fs.unit.gas_outlet.pressure[tf].value
        )
        assert (
            pytest.approx(0.19996, rel=1e-5)
            == iron_oc.fs.unit.solid_properties[tf, xf].particle_porosity.value
        )
        assert (
            pytest.approx(0.55007, rel=1e-5)
            == iron_oc.fs.unit.solid_properties[tf, xf].mass_frac_comp["Al2O3"].value
        )
        assert (
            pytest.approx(0.44599, rel=1e-5)
            == iron_oc.fs.unit.solid_properties[tf, xf].mass_frac_comp["Fe2O3"].value
        )
        assert (
            pytest.approx(0.00393951, rel=1e-5)
            == iron_oc.fs.unit.solid_properties[tf, xf].mass_frac_comp["Fe3O4"].value
        )
        assert (
            pytest.approx(1272.86689, rel=1e-5)
            == iron_oc.fs.unit.solid_properties[tf, xf].temperature.value
        )

    @pytest.mark.component
    @pytest.mark.xfail(reason="known DAE unit consistency issue, see Pyomo/pyomo#1790")
    def test_units_consistent(self, iron_oc):
        assert_units_consistent(iron_oc)

    @pytest.mark.ui
    @pytest.mark.unit
    def test_report(self, iron_oc):
        iron_oc.fs.unit.report()
        # no solid inlet/outlet ports (no stream table), test block reports instead
        t0 = iron_oc.fs.time.first()
        xf = iron_oc.fs.unit.length_domain.last()
        iron_oc.fs.unit.solid_properties[t0, xf].report()

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_performance_contents(self, iron_oc):
        perf_dict = iron_oc.fs.unit._get_performance_contents()

        assert perf_dict == {
            "vars": {
                "Bed Height": iron_oc.fs.unit.bed_height,
                "Bed Area": iron_oc.fs.unit.bed_area,
                "Gas Inlet Velocity": iron_oc.fs.unit.velocity_superficial_gas[0, 0],
                "Gas Outlet Velocity": iron_oc.fs.unit.velocity_superficial_gas[0, 1],
            }
        }

    @pytest.mark.component
    def test_initialization_error(self, iron_oc):
        tf = iron_oc.fs.time.last()
        iron_oc.fs.unit.gas_outlet.flow_mol[tf].fix(0)

        with pytest.raises(InitializationError):
            iron_oc.fs.unit.initialize()


# -----------------------------------------------------------------------------
class TestIronOC_reverse_flow(object):
    @pytest.fixture(scope="class")
    def iron_oc(self):
        m = ConcreteModel()
        horizon = 360
        m.fs = FlowsheetBlock(dynamic=True, time_set=[0, horizon], time_units=pyunits.s)

        m.fs.gas_props = GasPhaseParameterBlock()
        m.fs.solid_props = SolidPhaseParameterBlock()
        m.fs.solid_rxns = HeteroReactionParameterBlock(
            solid_property_package=m.fs.solid_props, gas_property_package=m.fs.gas_props
        )

        m.fs.unit = FixedBed1D(
            flow_type="reverse_flow",
            gas_phase_config={"property_package": m.fs.gas_props},
            solid_phase_config={
                "property_package": m.fs.solid_props,
                "reaction_package": m.fs.solid_rxns,
            },
        )

        # Discretize time domain
        t_element_size = 120  # s
        ntfe = int(horizon / t_element_size)
        m.discretizer = TransformationFactory("dae.finite_difference")
        m.discretizer.apply_to(m, nfe=ntfe, wrt=m.fs.time, scheme="BACKWARD")

        # Set reactor design conditions
        m.fs.unit.bed_diameter.fix(9)  # diameter of the reactor [m]
        m.fs.unit.bed_height.fix(5)  # height of solids in the reactor [m]

        # Fix inlet port variables (boundary conditions) for gas
        for t in m.fs.time:
            m.fs.unit.gas_inlet.flow_mol[t].fix(20)  # mol/s
            m.fs.unit.gas_inlet.temperature[t].fix(723.15)  # K
            m.fs.unit.gas_inlet.pressure[t].fix(1.565e5)  # Pa
            m.fs.unit.gas_inlet.mole_frac_comp[t, "CO2"].fix(0.02499)
            m.fs.unit.gas_inlet.mole_frac_comp[t, "H2O"].fix(0.00001)
            m.fs.unit.gas_inlet.mole_frac_comp[t, "CH4"].fix(0.975)

        # Specify gas phase initial conditions
        t0 = m.fs.time.first()
        for x in m.fs.unit.length_domain:
            if x != m.fs.unit.length_domain.last():  # skip inlet as it's already fixed
                m.fs.unit.gas_phase.properties[t0, x].mole_frac_comp["CO2"].fix(0.02499)
                m.fs.unit.gas_phase.properties[t0, x].mole_frac_comp["H2O"].fix(0.00001)
                m.fs.unit.gas_phase.properties[t0, x].mole_frac_comp["CH4"].fix(0.975)
                m.fs.unit.gas_phase.properties[t0, x].temperature.fix(1273.15)  # K
                m.fs.unit.gas_phase.properties[t0, x].flow_mol.fix(20)  # mol/s

        # Specify solid phase initial conditions
        t0 = m.fs.time.first()
        for x in m.fs.unit.length_domain:
            m.fs.unit.solid_properties[t0, x].temperature.fix(1273.15)
            m.fs.unit.solid_properties[t0, x].particle_porosity.fix(0.20)
            m.fs.unit.solid_properties[t0, x].mass_frac_comp["Fe2O3"].fix(0.45)
            m.fs.unit.solid_properties[t0, x].mass_frac_comp["Fe3O4"].fix(0)
            m.fs.unit.solid_properties[t0, x].mass_frac_comp["Al2O3"].fix(0.55)

        return m

    @pytest.mark.unit
    def test_dof(self, iron_oc):
        assert degrees_of_freedom(iron_oc) == 0

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, iron_oc):
        # Calculate scaling factors
        iscale.calculate_scaling_factors(iron_oc)

        optarg = {
            "tol": 1e-5,
            "bound_push": 1e-22,
            "nlp_scaling_method": "user-scaling",
        }

        initialization_tester(iron_oc, optarg=optarg)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_block_triangularization_initialization(self, iron_oc):
        iron_oc.fs.unit.block_triangularization_initialize(calc_var_kwds={"eps": 1e-5})

        # Calculate scaling factors
        iscale.calculate_scaling_factors(iron_oc)

        assert degrees_of_freedom(iron_oc) == 0

        # Assert that model is still fixed and deactivated as expected
        t0 = iron_oc.fs.time.first()
        x0 = iron_oc.fs.unit.length_domain.first()
        assert iron_oc.fs.unit.solid_properties[t0, x0].particle_porosity.fixed

        for t in iron_oc.fs.time:
            if t != iron_oc.fs.time.first():
                assert not iron_oc.fs.unit.solid_properties[
                    t, x0
                ].particle_porosity.fixed
            assert iron_oc.fs.unit.solid_properties[
                t, x0
            ].density_particle_constraint.active
            assert iron_oc.fs.unit.solid_properties[
                t, x0
            ].density_skeletal_constraint.active

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize_by_time(self, iron_oc):

        optarg = {
            "tol": 1e-5,
            "bound_push": 1e-22,
            "nlp_scaling_method": "user-scaling",
        }
        solver = get_solver("ipopt", optarg)  # create solver

        initialize_by_time_element(iron_oc.fs, iron_oc.fs.time, solver=solver)

        assert degrees_of_freedom(iron_oc) == 0

        # Assert that model is still fixed and deactivated as expected
        t0 = iron_oc.fs.time.first()
        x0 = iron_oc.fs.unit.length_domain.first()
        assert iron_oc.fs.unit.solid_properties[t0, x0].particle_porosity.fixed

        for t in iron_oc.fs.time:
            if t != iron_oc.fs.time.first():
                assert not iron_oc.fs.unit.solid_properties[
                    t, x0
                ].particle_porosity.fixed
            assert iron_oc.fs.unit.solid_properties[
                t, x0
            ].density_particle_constraint.active
            assert iron_oc.fs.unit.solid_properties[
                t, x0
            ].density_skeletal_constraint.active

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
        assert_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, iron_oc):
        t0 = iron_oc.fs.time.first()
        x0 = iron_oc.fs.unit.length_domain.first()
        tf = iron_oc.fs.time.last()
        xf = iron_oc.fs.unit.length_domain.last()
        assert (
            pytest.approx(0.01209, abs=1e-2)
            == iron_oc.fs.unit.velocity_superficial_gas[t0, x0].value
        )
        assert (
            pytest.approx(0.012078, abs=1e-2)
            == iron_oc.fs.unit.velocity_superficial_gas[tf, xf].value
        )
        # Check the pressure drop that occurs across the bed
        assert (
            pytest.approx(146126.06009, abs=1e-2)
            == iron_oc.fs.unit.gas_outlet.pressure[tf].value
        )
        assert (
            pytest.approx(10373.9409, abs=1e-2)
            == iron_oc.fs.unit.gas_inlet.pressure[tf].value
            - iron_oc.fs.unit.gas_outlet.pressure[tf].value
        )
        assert (
            pytest.approx(0.196012, rel=1e-5)
            == iron_oc.fs.unit.solid_properties[tf, xf].particle_porosity.value
        )
        assert (
            pytest.approx(0.55730, rel=1e-5)
            == iron_oc.fs.unit.solid_properties[tf, xf].mass_frac_comp["Al2O3"].value
        )
        assert (
            pytest.approx(0.05858366, rel=1e-5)
            == iron_oc.fs.unit.solid_properties[tf, xf].mass_frac_comp["Fe2O3"].value
        )
        assert (
            pytest.approx(0.38411633, rel=1e-5)
            == iron_oc.fs.unit.solid_properties[tf, xf].mass_frac_comp["Fe3O4"].value
        )
        assert (
            pytest.approx(723.266955, rel=1e-5)
            == iron_oc.fs.unit.solid_properties[tf, xf].temperature.value
        )

    @pytest.mark.component
    @pytest.mark.xfail(reason="known DAE unit consistency issue, see Pyomo/pyomo#1790")
    def test_units_consistent(self, iron_oc):
        assert_units_consistent(iron_oc)

    @pytest.mark.ui
    @pytest.mark.unit
    def test_report(self, iron_oc):
        iron_oc.fs.unit.report()
        # no solid inlet/outlet ports (no stream table), test block reports instead
        t0 = iron_oc.fs.time.first()
        xf = iron_oc.fs.unit.length_domain.last()
        iron_oc.fs.unit.solid_properties[t0, xf].report()

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_performance_contents(self, iron_oc):
        perf_dict = iron_oc.fs.unit._get_performance_contents()

        assert perf_dict == {
            "vars": {
                "Bed Height": iron_oc.fs.unit.bed_height,
                "Bed Area": iron_oc.fs.unit.bed_area,
                "Gas Inlet Velocity": iron_oc.fs.unit.velocity_superficial_gas[0, 0],
                "Gas Outlet Velocity": iron_oc.fs.unit.velocity_superficial_gas[0, 1],
            }
        }

    @pytest.mark.component
    def test_initialization_error(self, iron_oc):
        tf = iron_oc.fs.time.last()
        iron_oc.fs.unit.gas_outlet.flow_mol[tf].fix(0)

        with pytest.raises(InitializationError):
            iron_oc.fs.unit.initialize()


# -----------------------------------------------------------------------------
class TestIronOC_NoReaction(object):
    @pytest.fixture(scope="class")
    def iron_oc(self):
        m = ConcreteModel()
        horizon = 360
        m.fs = FlowsheetBlock(dynamic=True, time_set=[0, horizon], time_units=pyunits.s)

        m.fs.gas_props = GasPhaseParameterBlock()
        m.fs.solid_props = SolidPhaseParameterBlock()

        m.fs.unit = FixedBed1D(
            gas_phase_config={"property_package": m.fs.gas_props},
            solid_phase_config={"property_package": m.fs.solid_props},
        )

        # Discretize time domain
        t_element_size = 120  # s
        ntfe = int(horizon / t_element_size)
        m.discretizer = TransformationFactory("dae.finite_difference")
        m.discretizer.apply_to(m, nfe=ntfe, wrt=m.fs.time, scheme="BACKWARD")

        # Set reactor design conditions
        m.fs.unit.bed_diameter.fix(9)  # diameter of the reactor [m]
        m.fs.unit.bed_height.fix(5)  # height of solids in the reactor [m]

        # Fix inlet port variables (boundary conditions) for gas
        for t in m.fs.time:
            m.fs.unit.gas_inlet.flow_mol[t].fix(20)  # mol/s
            m.fs.unit.gas_inlet.temperature[t].fix(1273.15)  # K
            m.fs.unit.gas_inlet.pressure[t].fix(1.565e5)  # Pa
            m.fs.unit.gas_inlet.mole_frac_comp[t, "CO2"].fix(0.02499)
            m.fs.unit.gas_inlet.mole_frac_comp[t, "H2O"].fix(0.00001)
            m.fs.unit.gas_inlet.mole_frac_comp[t, "CH4"].fix(0.975)

        # Specify gas phase initial conditions
        t0 = m.fs.time.first()
        for x in m.fs.unit.length_domain:
            if x != m.fs.unit.length_domain.first():  # skip inlet as it's already fixed
                m.fs.unit.gas_phase.properties[t0, x].mole_frac_comp["CO2"].fix(0.02499)
                m.fs.unit.gas_phase.properties[t0, x].mole_frac_comp["H2O"].fix(0.00001)
                m.fs.unit.gas_phase.properties[t0, x].mole_frac_comp["CH4"].fix(0.975)
                m.fs.unit.gas_phase.properties[t0, x].temperature.fix(1273.15)  # K
                m.fs.unit.gas_phase.properties[t0, x].flow_mol.fix(20)  # mol/s

        # Specify solid phase initial conditions
        t0 = m.fs.time.first()
        for x in m.fs.unit.length_domain:
            m.fs.unit.solid_properties[t0, x].temperature.fix(1273.15)
            m.fs.unit.solid_properties[t0, x].particle_porosity.fix(0.20)
            m.fs.unit.solid_properties[t0, x].mass_frac_comp["Fe2O3"].fix(0.45)
            m.fs.unit.solid_properties[t0, x].mass_frac_comp["Fe3O4"].fix(0)
            m.fs.unit.solid_properties[t0, x].mass_frac_comp["Al2O3"].fix(0.55)

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

        assert hasattr(iron_oc.fs.unit, "gas_outlet")
        assert len(iron_oc.fs.unit.gas_outlet.vars) == 4
        assert isinstance(iron_oc.fs.unit.gas_outlet.flow_mol, Var)
        assert isinstance(iron_oc.fs.unit.gas_outlet.mole_frac_comp, Var)
        assert isinstance(iron_oc.fs.unit.gas_outlet.temperature, Var)
        assert isinstance(iron_oc.fs.unit.gas_outlet.pressure, Var)

        assert isinstance(iron_oc.fs.unit.bed_area_eqn, Constraint)
        assert isinstance(iron_oc.fs.unit.gas_phase_area_constraint, Constraint)
        assert isinstance(iron_oc.fs.unit.solid_phase_area_constraint, Constraint)
        assert isinstance(iron_oc.fs.unit.gas_super_vel, Constraint)
        assert isinstance(iron_oc.fs.unit.gas_phase_config_pressure_drop, Constraint)
        assert isinstance(iron_oc.fs.unit.solid_material_holdup_calculation, Constraint)
        assert isinstance(iron_oc.fs.unit.solid_material_balances, Constraint)
        assert isinstance(iron_oc.fs.unit.solid_sum_component_eqn, Constraint)
        assert isinstance(iron_oc.fs.unit.solid_phase_heat_transfer, Constraint)
        assert isinstance(iron_oc.fs.unit.solid_energy_holdup_calculation, Constraint)
        assert isinstance(iron_oc.fs.unit.solid_enthalpy_balances, Constraint)

        assert isinstance(iron_oc.fs.unit.reynolds_number_particle, Constraint)
        assert isinstance(iron_oc.fs.unit.prandtl_number, Constraint)
        assert isinstance(iron_oc.fs.unit.nusselt_number_particle, Constraint)
        assert isinstance(iron_oc.fs.unit.gas_solid_htc_eqn, Constraint)
        assert isinstance(iron_oc.fs.unit.gas_phase_heat_transfer, Constraint)

        assert number_variables(iron_oc) == 3088
        assert number_total_constraints(iron_oc) == 2894
        assert number_unused_variables(iron_oc) == 71

    @pytest.mark.unit
    def test_dof(self, iron_oc):
        assert degrees_of_freedom(iron_oc) == 0

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, iron_oc):
        optarg = {
            "tol": 1e-5,
            "bound_push": 1e-22,
            "nlp_scaling_method": "user-scaling",
        }
        # Calculate scaling factors
        iscale.calculate_scaling_factors(iron_oc)

        initialization_tester(iron_oc, optarg=optarg)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_block_triangularization_initialization(self, iron_oc):
        iron_oc.fs.unit.block_triangularization_initialize(calc_var_kwds={"eps": 1e-5})

        # Calculate scaling factors
        iscale.calculate_scaling_factors(iron_oc)

        assert degrees_of_freedom(iron_oc) == 0

        # Assert that model is still fixed and deactivated as expected
        t0 = iron_oc.fs.time.first()
        x0 = iron_oc.fs.unit.length_domain.first()
        assert iron_oc.fs.unit.solid_properties[t0, x0].particle_porosity.fixed

        for t in iron_oc.fs.time:
            if t != iron_oc.fs.time.first():
                assert not iron_oc.fs.unit.solid_properties[
                    t, x0
                ].particle_porosity.fixed
            assert iron_oc.fs.unit.solid_properties[
                t, x0
            ].density_particle_constraint.active
            assert iron_oc.fs.unit.solid_properties[
                t, x0
            ].density_skeletal_constraint.active

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize_by_time(self, iron_oc):

        optarg = {
            "tol": 1e-5,
            "bound_push": 1e-22,
            "nlp_scaling_method": "user-scaling",
        }
        solver = get_solver("ipopt", optarg)  # create solver

        initialize_by_time_element(iron_oc.fs, iron_oc.fs.time, solver=solver)

        assert degrees_of_freedom(iron_oc) == 0

        # Assert that model is still fixed and deactivated as expected
        t0 = iron_oc.fs.time.first()
        x0 = iron_oc.fs.unit.length_domain.first()
        assert iron_oc.fs.unit.solid_properties[t0, x0].particle_porosity.fixed

        for t in iron_oc.fs.time:
            if t != iron_oc.fs.time.first():
                assert not iron_oc.fs.unit.solid_properties[
                    t, x0
                ].particle_porosity.fixed
            assert iron_oc.fs.unit.solid_properties[
                t, x0
            ].density_particle_constraint.active
            assert iron_oc.fs.unit.solid_properties[
                t, x0
            ].density_skeletal_constraint.active

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
        assert_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, iron_oc):
        t0 = iron_oc.fs.time.first()
        x0 = iron_oc.fs.unit.length_domain.first()
        tf = iron_oc.fs.time.last()
        xf = iron_oc.fs.unit.length_domain.last()
        assert (
            pytest.approx(0.01209, abs=1e-2)
            == iron_oc.fs.unit.velocity_superficial_gas[t0, x0].value
        )
        assert (
            pytest.approx(0.02159, abs=1e-2)
            == iron_oc.fs.unit.velocity_superficial_gas[tf, xf].value
        )
        # Check the pressure drop that occurs across the bed
        assert (
            pytest.approx(154128.50977, abs=1e-2)
            == iron_oc.fs.unit.gas_outlet.pressure[tf].value
        )
        assert (
            pytest.approx(2371.49023, abs=1e-2)
            == iron_oc.fs.unit.gas_inlet.pressure[tf].value
            - iron_oc.fs.unit.gas_outlet.pressure[tf].value
        )
        assert (
            pytest.approx(0.2, rel=1e-5)
            == iron_oc.fs.unit.solid_properties[tf, xf].particle_porosity.value
        )
        assert (
            pytest.approx(0.55, rel=1e-5)
            == iron_oc.fs.unit.solid_properties[tf, xf].mass_frac_comp["Al2O3"].value
        )
        assert (
            pytest.approx(0.45, rel=1e-5)
            == iron_oc.fs.unit.solid_properties[tf, xf].mass_frac_comp["Fe2O3"].value
        )
        assert (
            pytest.approx(0.00, rel=1e-5)
            == iron_oc.fs.unit.solid_properties[tf, xf].mass_frac_comp["Fe3O4"].value
        )
        assert (
            pytest.approx(1273.15, rel=1e-5)
            == iron_oc.fs.unit.solid_properties[tf, xf].temperature.value
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, iron_oc):

        iron_oc.fs.unit.time_set = ContinuousSet(
            bounds=(iron_oc.fs.time.first(), iron_oc.fs.time.last()),
            initialize=iron_oc.fs.time.get_finite_elements(),
            doc="time domain",
        )

        ###########################################################################
        #  Mass conservation check
        ###########################################################################
        def gas_input_total(b, t):
            calculate_variable_from_constraint(
                b.gas_phase.properties[t, 0].mw,
                b.gas_phase.properties[t, 0].mw_eqn,
            )
            return b.gas_phase.properties[t, 0].mw * b.gas_inlet.flow_mol[t]

        iron_oc.fs.unit.gas_input_total = Integral(
            iron_oc.fs.unit.time_set,
            wrt=iron_oc.fs.unit.time_set,
            rule=gas_input_total,
            doc="Gas input total, kg",
        )

        def gas_output_total(b, t):
            calculate_variable_from_constraint(
                b.gas_phase.properties[t, 1].mw,
                b.gas_phase.properties[t, 1].mw_eqn,
            )
            return b.gas_phase.properties[t, 1].mw * b.gas_outlet.flow_mol[t]

        iron_oc.fs.unit.gas_output_total = Integral(
            iron_oc.fs.unit.time_set,
            wrt=iron_oc.fs.unit.time_set,
            rule=gas_output_total,
            doc="Gas output total, kg",
        )

        def gas_holdup(b, x):
            tf = b.flowsheet().time.last()
            t0 = b.flowsheet().time.first()
            calculate_variable_from_constraint(
                b.gas_phase.properties[tf, x].mw,
                b.gas_phase.properties[tf, x].mw_eqn,
            )
            calculate_variable_from_constraint(
                b.gas_phase.properties[t0, x].mw,
                b.gas_phase.properties[t0, x].mw_eqn,
            )
            mass_final_per_length = (
                b.gas_phase.properties[tf, x].dens_mol
                * b.gas_phase.properties[tf, 1].mw
                * b.gas_phase.area[tf, x]
            )
            mass_initial_per_length = (
                b.gas_phase.properties[t0, x].dens_mol
                * b.gas_phase.properties[t0, 1].mw
                * b.gas_phase.area[t0, x]
            )
            return (mass_final_per_length - mass_initial_per_length) * b.bed_height

        iron_oc.fs.unit.gas_holdup = Integral(
            iron_oc.fs.unit.length_domain,
            wrt=iron_oc.fs.unit.length_domain,
            rule=gas_holdup,
            doc="Gas mass holdup, kg",
        )

        def oc_at_initial_time(b, x):
            t0 = b.flowsheet().time.first()
            calculate_variable_from_constraint(
                b.solid_properties[t0, x].dens_mass_particle,
                b.solid_properties[t0, x].density_particle_constraint,
            )
            mass_per_length = (
                b.solid_properties[t0, x].dens_mass_particle * b.solid_phase_area[t0, x]
            )
            return mass_per_length * b.bed_height

        iron_oc.fs.unit.oc_at_initial_time = Integral(
            iron_oc.fs.unit.length_domain,
            wrt=iron_oc.fs.unit.length_domain,
            rule=oc_at_initial_time,
            doc="OC at initial reaction time, kg",
        )

        def oc_at_final_time(b, x):
            tf = b.flowsheet().time.last()
            calculate_variable_from_constraint(
                b.solid_properties[tf, x].dens_mass_particle,
                b.solid_properties[tf, x].density_particle_constraint,
            )
            mass_per_length = (
                b.solid_properties[tf, x].dens_mass_particle * b.solid_phase_area[tf, x]
            )
            return mass_per_length * b.bed_height

        iron_oc.fs.unit.oc_at_final_time = Integral(
            iron_oc.fs.unit.length_domain,
            wrt=iron_oc.fs.unit.length_domain,
            rule=oc_at_final_time,
            doc="OC at final reaction time, kg",
        )

        mass_change_gas = value(
            iron_oc.fs.unit.gas_input_total
            - iron_oc.fs.unit.gas_output_total
            - iron_oc.fs.unit.gas_holdup
        )
        mass_change_solid = value(
            iron_oc.fs.unit.oc_at_initial_time - iron_oc.fs.unit.oc_at_final_time
        )

        assert abs(mass_change_gas) <= 1e-2
        assert abs(mass_change_solid) <= 1e-2

        ###########################################################################
        #  Energy conservation check
        ###########################################################################
        def ch4_reacted_mol(b, t):
            return (
                b.gas_inlet.flow_mol[t] * b.gas_inlet.mole_frac_comp[t, "CH4"]
                - b.gas_outlet.flow_mol[t] * b.gas_outlet.mole_frac_comp[t, "CH4"]
            )

        iron_oc.fs.unit.ch4_reacted_mol = Integral(
            iron_oc.fs.unit.time_set,
            wrt=iron_oc.fs.unit.time_set,
            rule=ch4_reacted_mol,
            doc="CH4 reacted total, mol",
        )

        def fe2o3_at_initial_time_mol(b, x):
            t0 = b.flowsheet().time.first()
            calculate_variable_from_constraint(
                b.solid_properties[t0, x].dens_mass_particle,
                b.solid_properties[t0, x].density_particle_constraint,
            )
            mol_fe2o3_per_length = (
                b.solid_properties[t0, x].dens_mass_particle
                * b.solid_properties[t0, x].mass_frac_comp["Fe2O3"]
                * b.solid_phase_area[t0, x]
                / b.solid_properties[t0, x]._params.mw_comp["Fe2O3"]
            )
            return mol_fe2o3_per_length * b.bed_height

        iron_oc.fs.unit.fe2o3_at_initial_time_mol = Integral(
            iron_oc.fs.unit.length_domain,
            wrt=iron_oc.fs.unit.length_domain,
            rule=fe2o3_at_initial_time_mol,
            doc="Fe2O3 at initial reaction time, mol",
        )

        def fe2o3_at_final_time_mol(b, x):
            tf = b.flowsheet().time.last()
            calculate_variable_from_constraint(
                b.solid_properties[tf, x].dens_mass_particle,
                b.solid_properties[tf, x].density_particle_constraint,
            )
            mol_fe2o3_per_length = (
                b.solid_properties[tf, x].dens_mass_particle
                * b.solid_properties[tf, x].mass_frac_comp["Fe2O3"]
                * b.solid_phase_area[tf, x]
                / b.solid_properties[tf, x]._params.mw_comp["Fe2O3"]
            )
            return mol_fe2o3_per_length * b.bed_height

        iron_oc.fs.unit.fe2o3_at_final_time_mol = Integral(
            iron_oc.fs.unit.length_domain,
            wrt=iron_oc.fs.unit.length_domain,
            rule=fe2o3_at_final_time_mol,
            doc="Fe2O3 at final reaction time, mol",
        )

        def gas_input_energy_total(b, t):
            return b.gas_phase.properties[t, 0].enth_mol * b.gas_inlet.flow_mol[t]

        iron_oc.fs.unit.gas_input_energy_total = Integral(
            iron_oc.fs.unit.time_set,
            wrt=iron_oc.fs.unit.time_set,
            rule=gas_input_energy_total,
            doc="Gas input energy total, J",
        )

        def gas_output_energy_total(b, t):
            return b.gas_phase.properties[t, 1].enth_mol * b.gas_outlet.flow_mol[t]

        iron_oc.fs.unit.gas_output_energy_total = Integral(
            iron_oc.fs.unit.time_set,
            wrt=iron_oc.fs.unit.time_set,
            rule=gas_output_energy_total,
            doc="Gas output energy total, J",
        )

        def gas_energy_holdup(b, x):
            tf = b.flowsheet().time.last()
            t0 = b.flowsheet().time.first()
            energy_final_per_length = (
                b.gas_phase.properties[tf, x].dens_mol
                * b.gas_phase.properties[tf, x].enth_mol
                * b.gas_phase.area[tf, x]
            )
            energy_initial_per_length = (
                b.gas_phase.properties[t0, x].dens_mol
                * b.gas_phase.properties[t0, x].enth_mol
                * b.gas_phase.area[t0, x]
            )
            return (energy_final_per_length - energy_initial_per_length) * b.bed_height

        iron_oc.fs.unit.gas_energy_holdup = Integral(
            iron_oc.fs.unit.length_domain,
            wrt=iron_oc.fs.unit.length_domain,
            rule=gas_energy_holdup,
            doc="Gas energy holdup, J",
        )

        def oc_energy_at_initial_time(b, x):
            t0 = b.flowsheet().time.first()
            energy_per_length = (
                b.solid_properties[t0, x].dens_mass_particle
                * b.solid_properties[t0, x].enth_mass
                * b.solid_phase_area[t0, x]
            )
            return energy_per_length * b.bed_height

        iron_oc.fs.unit.oc_energy_at_initial_time = Integral(
            iron_oc.fs.unit.length_domain,
            wrt=iron_oc.fs.unit.length_domain,
            rule=oc_energy_at_initial_time,
            doc="OC energy at initial reaction time, J",
        )

        def oc_energy_at_final_time(b, x):
            tf = b.flowsheet().time.last()
            energy_per_length = (
                b.solid_properties[tf, x].dens_mass_particle
                * b.solid_properties[tf, x].enth_mass
                * b.solid_phase_area[tf, x]
            )
            return energy_per_length * b.bed_height

        iron_oc.fs.unit.oc_energy_at_final_time = Integral(
            iron_oc.fs.unit.length_domain,
            wrt=iron_oc.fs.unit.length_domain,
            rule=oc_energy_at_final_time,
            doc="OC at final reaction time, J",
        )

        if (
            iron_oc.fs.unit.config.energy_balance_type is not EnergyBalanceType.none
            and iron_oc.fs.unit.config.solid_phase_config.reaction_package is not None
        ):
            e_reaction = (
                value(
                    (
                        iron_oc.fs.unit.fe2o3_at_final_time_mol
                        - iron_oc.fs.unit.fe2o3_at_initial_time_mol
                    )
                    * iron_oc.fs.unit.solid_reactions[0, 1]._params.dh_rxn["R1"]
                )
                / 12
            )
        else:
            e_reaction = 0

        enthalpy_change_gas = value(  # in - out - holdup/accumulation
            iron_oc.fs.unit.gas_input_energy_total
            - iron_oc.fs.unit.gas_output_energy_total
            - iron_oc.fs.unit.gas_energy_holdup
        )

        enthalpy_change_solid = value(
            iron_oc.fs.unit.oc_energy_at_final_time
            - iron_oc.fs.unit.oc_energy_at_initial_time
            - e_reaction
        )

        ebal_abs_tol = enthalpy_change_solid - enthalpy_change_gas

        assert abs(ebal_abs_tol) <= 1e-8

    @pytest.mark.component
    @pytest.mark.xfail(reason="known DAE unit consistency issue, see Pyomo/pyomo#1790")
    def test_units_consistent(self, iron_oc):
        assert_units_consistent(iron_oc)

    @pytest.mark.ui
    @pytest.mark.unit
    def test_report(self, iron_oc):
        iron_oc.fs.unit.report()
        # no solid inlet/outlet ports (no stream table), test block reports instead
        t0 = iron_oc.fs.time.first()
        xf = iron_oc.fs.unit.length_domain.last()
        iron_oc.fs.unit.solid_properties[t0, xf].report()

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_performance_contents(self, iron_oc):
        perf_dict = iron_oc.fs.unit._get_performance_contents()

        assert perf_dict == {
            "vars": {
                "Bed Height": iron_oc.fs.unit.bed_height,
                "Bed Area": iron_oc.fs.unit.bed_area,
                "Gas Inlet Velocity": iron_oc.fs.unit.velocity_superficial_gas[0, 0],
                "Gas Outlet Velocity": iron_oc.fs.unit.velocity_superficial_gas[0, 1],
            }
        }

    @pytest.mark.component
    def test_initialization_error(self, iron_oc):
        tf = iron_oc.fs.time.last()
        iron_oc.fs.unit.gas_outlet.flow_mol[tf].fix(0)

        with pytest.raises(InitializationError):
            iron_oc.fs.unit.initialize()


class TestIronOC_conservation_with_reaction(object):
    @pytest.fixture(scope="class")
    def iron_oc(self):
        m = ConcreteModel()
        horizon = 3600
        m.fs = FlowsheetBlock(dynamic=True, time_set=[0, horizon], time_units=pyunits.s)

        m.fs.gas_props = GasPhaseParameterBlock()
        m.fs.solid_props = SolidPhaseParameterBlock()
        m.fs.solid_rxns = HeteroReactionParameterBlock(
            solid_property_package=m.fs.solid_props, gas_property_package=m.fs.gas_props
        )

        m.fs.unit = FixedBed1D(
            finite_elements=20,
            gas_phase_config={"property_package": m.fs.gas_props},
            solid_phase_config={
                "property_package": m.fs.solid_props,
                "reaction_package": m.fs.solid_rxns,
            },
        )

        # Discretize time domain
        t_element_size = 120  # s
        ntfe = int(horizon / t_element_size)
        m.discretizer = TransformationFactory("dae.finite_difference")
        m.discretizer.apply_to(m, nfe=ntfe, wrt=m.fs.time, scheme="BACKWARD")

        # Set reactor design conditions
        m.fs.unit.bed_diameter.fix(9)  # diameter of the reactor [m]
        m.fs.unit.bed_height.fix(5)  # height of solids in the reactor [m]

        # Fix inlet port variables (boundary conditions) for gas
        for t in m.fs.time:
            m.fs.unit.gas_inlet.flow_mol[t].fix(20)  # mol/s
            m.fs.unit.gas_inlet.temperature[t].fix(1273.15)  # K
            m.fs.unit.gas_inlet.pressure[t].fix(1.565e5)  # Pa
            m.fs.unit.gas_inlet.mole_frac_comp[t, "CO2"].fix(0.02499)
            m.fs.unit.gas_inlet.mole_frac_comp[t, "H2O"].fix(0.00001)
            m.fs.unit.gas_inlet.mole_frac_comp[t, "CH4"].fix(0.975)

        # Specify gas phase initial conditions
        t0 = m.fs.time.first()
        for x in m.fs.unit.length_domain:
            if x != m.fs.unit.length_domain.first():  # skip inlet as it's already fixed
                m.fs.unit.gas_phase.properties[t0, x].mole_frac_comp["CO2"].fix(0.02499)
                m.fs.unit.gas_phase.properties[t0, x].mole_frac_comp["H2O"].fix(0.00001)
                m.fs.unit.gas_phase.properties[t0, x].mole_frac_comp["CH4"].fix(0.975)
                m.fs.unit.gas_phase.properties[t0, x].temperature.fix(1273.15)  # K
                m.fs.unit.gas_phase.properties[t0, x].flow_mol.fix(20)  # mol/s

        # Specify solid phase initial conditions
        t0 = m.fs.time.first()
        for x in m.fs.unit.length_domain:
            m.fs.unit.solid_properties[t0, x].temperature.fix(1273.15)
            m.fs.unit.solid_properties[t0, x].particle_porosity.fix(0.20)
            m.fs.unit.solid_properties[t0, x].mass_frac_comp["Fe2O3"].fix(0.45)
            m.fs.unit.solid_properties[t0, x].mass_frac_comp["Fe3O4"].fix(0)
            m.fs.unit.solid_properties[t0, x].mass_frac_comp["Al2O3"].fix(0.55)

        return m

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, iron_oc):
        # Calculate scaling factors
        iscale.calculate_scaling_factors(iron_oc)

        # Initialize model
        iron_oc.fs.unit.block_triangularization_initialize(calc_var_kwds={"eps": 1e-5})

        # Solve model
        initialize_by_time_element(iron_oc.fs, iron_oc.fs.time, solver=solver)

        iron_oc.fs.unit.time_set = ContinuousSet(
            bounds=(iron_oc.fs.time.first(), iron_oc.fs.time.last()),
            initialize=iron_oc.fs.time.get_finite_elements(),
            doc="time domain",
        )

        ###########################################################################
        #  Mass conservation check
        ###########################################################################
        def gas_input_total(b, t):
            calculate_variable_from_constraint(
                b.gas_phase.properties[t, 0].mw,
                b.gas_phase.properties[t, 0].mw_eqn,
            )
            return b.gas_phase.properties[t, 0].mw * b.gas_inlet.flow_mol[t]

        iron_oc.fs.unit.gas_input_total = Integral(
            iron_oc.fs.unit.time_set,
            wrt=iron_oc.fs.unit.time_set,
            rule=gas_input_total,
            doc="Gas input total, kg",
        )

        def gas_output_total(b, t):
            calculate_variable_from_constraint(
                b.gas_phase.properties[t, 1].mw,
                b.gas_phase.properties[t, 1].mw_eqn,
            )
            return b.gas_phase.properties[t, 1].mw * b.gas_outlet.flow_mol[t]

        iron_oc.fs.unit.gas_output_total = Integral(
            iron_oc.fs.unit.time_set,
            wrt=iron_oc.fs.unit.time_set,
            rule=gas_output_total,
            doc="Gas output total, kg",
        )

        def gas_holdup(b, x):
            tf = b.flowsheet().time.last()
            t0 = b.flowsheet().time.first()
            calculate_variable_from_constraint(
                b.gas_phase.properties[tf, x].mw,
                b.gas_phase.properties[tf, x].mw_eqn,
            )
            calculate_variable_from_constraint(
                b.gas_phase.properties[t0, x].mw,
                b.gas_phase.properties[t0, x].mw_eqn,
            )
            mass_final_per_length = (
                b.gas_phase.properties[tf, x].dens_mol
                * b.gas_phase.properties[tf, 1].mw
                * b.gas_phase.area[tf, x]
            )
            mass_initial_per_length = (
                b.gas_phase.properties[t0, x].dens_mol
                * b.gas_phase.properties[t0, 1].mw
                * b.gas_phase.area[t0, x]
            )
            return (mass_final_per_length - mass_initial_per_length) * b.bed_height

        iron_oc.fs.unit.gas_holdup = Integral(
            iron_oc.fs.unit.length_domain,
            wrt=iron_oc.fs.unit.length_domain,
            rule=gas_holdup,
            doc="Gas mass holdup, kg",
        )

        def oc_at_initial_time(b, x):
            t0 = b.flowsheet().time.first()
            calculate_variable_from_constraint(
                b.solid_properties[t0, x].dens_mass_particle,
                b.solid_properties[t0, x].density_particle_constraint,
            )
            mass_per_length = (
                b.solid_properties[t0, x].dens_mass_particle * b.solid_phase_area[t0, x]
            )
            return mass_per_length * b.bed_height

        iron_oc.fs.unit.oc_at_initial_time = Integral(
            iron_oc.fs.unit.length_domain,
            wrt=iron_oc.fs.unit.length_domain,
            rule=oc_at_initial_time,
            doc="OC at initial reaction time, kg",
        )

        def oc_at_final_time(b, x):
            tf = b.flowsheet().time.last()
            calculate_variable_from_constraint(
                b.solid_properties[tf, x].dens_mass_particle,
                b.solid_properties[tf, x].density_particle_constraint,
            )
            mass_per_length = (
                b.solid_properties[tf, x].dens_mass_particle * b.solid_phase_area[tf, x]
            )
            return mass_per_length * b.bed_height

        iron_oc.fs.unit.oc_at_final_time = Integral(
            iron_oc.fs.unit.length_domain,
            wrt=iron_oc.fs.unit.length_domain,
            rule=oc_at_final_time,
            doc="OC at final reaction time, kg",
        )

        mass_change_gas = value(
            iron_oc.fs.unit.gas_input_total
            - iron_oc.fs.unit.gas_output_total
            - iron_oc.fs.unit.gas_holdup
        )
        mass_change_solid = value(
            iron_oc.fs.unit.oc_at_initial_time - iron_oc.fs.unit.oc_at_final_time
        )
        mbal_abs_tol = abs(mass_change_gas + mass_change_solid)
        mbal_rel_tol = mbal_abs_tol / mass_change_solid
        assert abs(mbal_rel_tol) <= 0.15

        ###########################################################################
        #  Energy conservation check
        ###########################################################################
        def ch4_reacted_mol(b, t):
            return (
                b.gas_inlet.flow_mol[t] * b.gas_inlet.mole_frac_comp[t, "CH4"]
                - b.gas_outlet.flow_mol[t] * b.gas_outlet.mole_frac_comp[t, "CH4"]
            )

        iron_oc.fs.unit.ch4_reacted_mol = Integral(
            iron_oc.fs.unit.time_set,
            wrt=iron_oc.fs.unit.time_set,
            rule=ch4_reacted_mol,
            doc="CH4 reacted total, mol",
        )

        def fe2o3_at_initial_time_mol(b, x):
            t0 = b.flowsheet().time.first()
            calculate_variable_from_constraint(
                b.solid_properties[t0, x].dens_mass_particle,
                b.solid_properties[t0, x].density_particle_constraint,
            )
            mol_fe2o3_per_length = (
                b.solid_properties[t0, x].dens_mass_particle
                * b.solid_properties[t0, x].mass_frac_comp["Fe2O3"]
                * b.solid_phase_area[t0, x]
                / b.solid_properties[t0, x]._params.mw_comp["Fe2O3"]
            )
            return mol_fe2o3_per_length * b.bed_height

        iron_oc.fs.unit.fe2o3_at_initial_time_mol = Integral(
            iron_oc.fs.unit.length_domain,
            wrt=iron_oc.fs.unit.length_domain,
            rule=fe2o3_at_initial_time_mol,
            doc="Fe2O3 at initial reaction time, mol",
        )

        def fe2o3_at_final_time_mol(b, x):
            tf = b.flowsheet().time.last()
            calculate_variable_from_constraint(
                b.solid_properties[tf, x].dens_mass_particle,
                b.solid_properties[tf, x].density_particle_constraint,
            )
            mol_fe2o3_per_length = (
                b.solid_properties[tf, x].dens_mass_particle
                * b.solid_properties[tf, x].mass_frac_comp["Fe2O3"]
                * b.solid_phase_area[tf, x]
                / b.solid_properties[tf, x]._params.mw_comp["Fe2O3"]
            )
            return mol_fe2o3_per_length * b.bed_height

        iron_oc.fs.unit.fe2o3_at_final_time_mol = Integral(
            iron_oc.fs.unit.length_domain,
            wrt=iron_oc.fs.unit.length_domain,
            rule=fe2o3_at_final_time_mol,
            doc="Fe2O3 at final reaction time, mol",
        )

        def gas_input_energy_total(b, t):
            return b.gas_phase.properties[t, 0].enth_mol * b.gas_inlet.flow_mol[t]

        iron_oc.fs.unit.gas_input_energy_total = Integral(
            iron_oc.fs.unit.time_set,
            wrt=iron_oc.fs.unit.time_set,
            rule=gas_input_energy_total,
            doc="Gas input energy total, J",
        )

        def gas_output_energy_total(b, t):
            return b.gas_phase.properties[t, 1].enth_mol * b.gas_outlet.flow_mol[t]

        iron_oc.fs.unit.gas_output_energy_total = Integral(
            iron_oc.fs.unit.time_set,
            wrt=iron_oc.fs.unit.time_set,
            rule=gas_output_energy_total,
            doc="Gas output energy total, J",
        )

        def gas_energy_holdup(b, x):
            tf = b.flowsheet().time.last()
            t0 = b.flowsheet().time.first()
            energy_final_per_length = (
                b.gas_phase.properties[tf, x].dens_mol
                * b.gas_phase.properties[tf, x].enth_mol
                * b.gas_phase.area[tf, x]
            )
            energy_initial_per_length = (
                b.gas_phase.properties[t0, x].dens_mol
                * b.gas_phase.properties[t0, x].enth_mol
                * b.gas_phase.area[t0, x]
            )
            return (energy_final_per_length - energy_initial_per_length) * b.bed_height

        iron_oc.fs.unit.gas_energy_holdup = Integral(
            iron_oc.fs.unit.length_domain,
            wrt=iron_oc.fs.unit.length_domain,
            rule=gas_energy_holdup,
            doc="Gas energy holdup, J",
        )

        def oc_energy_at_initial_time(b, x):
            t0 = b.flowsheet().time.first()
            energy_per_length = (
                b.solid_properties[t0, x].dens_mass_particle
                * b.solid_properties[t0, x].enth_mass
                * b.solid_phase_area[t0, x]
            )
            return energy_per_length * b.bed_height

        iron_oc.fs.unit.oc_energy_at_initial_time = Integral(
            iron_oc.fs.unit.length_domain,
            wrt=iron_oc.fs.unit.length_domain,
            rule=oc_energy_at_initial_time,
            doc="OC energy at initial reaction time, J",
        )

        def oc_energy_at_final_time(b, x):
            tf = b.flowsheet().time.last()
            energy_per_length = (
                b.solid_properties[tf, x].dens_mass_particle
                * b.solid_properties[tf, x].enth_mass
                * b.solid_phase_area[tf, x]
            )
            return energy_per_length * b.bed_height

        iron_oc.fs.unit.oc_energy_at_final_time = Integral(
            iron_oc.fs.unit.length_domain,
            wrt=iron_oc.fs.unit.length_domain,
            rule=oc_energy_at_final_time,
            doc="OC at final reaction time, J",
        )

        if iron_oc.fs.unit.config.solid_phase_config.reaction_package is not None:
            e_reaction = (
                value(
                    (
                        iron_oc.fs.unit.fe2o3_at_final_time_mol
                        - iron_oc.fs.unit.fe2o3_at_initial_time_mol
                    )
                    * iron_oc.fs.unit.solid_reactions[0, 1]._params.dh_rxn["R1"]
                )
                / 12
            )
        else:
            e_reaction = 0

        enthalpy_change_gas = value(  # in - out - holdup/accumulation
            iron_oc.fs.unit.gas_input_energy_total
            - iron_oc.fs.unit.gas_output_energy_total
            - iron_oc.fs.unit.gas_energy_holdup
        )

        enthalpy_change_solid = value(
            iron_oc.fs.unit.oc_energy_at_final_time
            - iron_oc.fs.unit.oc_energy_at_initial_time
            - e_reaction
        )

        ebal_abs_tol = enthalpy_change_solid - enthalpy_change_gas
        ebal_rel_tol = ebal_abs_tol / enthalpy_change_solid
        assert abs(ebal_rel_tol) <= 0.15
