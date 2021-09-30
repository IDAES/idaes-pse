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

from pyomo.environ import (ConcreteModel,
                           TransformationFactory,
                           TerminationCondition,
                           SolverStatus,
                           value,
                           Constraint)
from idaes.core import (FlowsheetBlock,
                        EnergyBalanceType)
from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              number_variables,
                                              number_total_constraints,
                                              number_unused_variables,
                                              unused_variables_set)
from idaes.core.util.testing import initialization_tester
from idaes.core.util import get_solver

# Import FixedBed0D unit model
from idaes.gas_solid_contactors.unit_models.fixed_bed_0D import FixedBed0D

# Import property packages
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
    m.fs = FlowsheetBlock(default={"dynamic": True, "time_set": [0, 3600]})

    # Set up thermo props and reaction props
    m.fs.gas_props = GasPhaseParameterBlock()
    m.fs.solid_props = SolidPhaseParameterBlock()
    m.fs.solid_rxns = HeteroReactionParameterBlock(
            default={"solid_property_package": m.fs.solid_props,
                     "gas_property_package": m.fs.gas_props})

    m.fs.unit = FixedBed0D(
        default={"gas_property_package": m.fs.gas_props,
                 "solid_property_package": m.fs.solid_props,
                 "reaction_package": m.fs.solid_rxns})

    # Check unit config arguments
    assert len(m.fs.unit.config) == 9

    assert m.fs.unit.config.dynamic is True
    assert m.fs.unit.config.has_holdup is True
    assert m.fs.unit.config.energy_balance_type == \
        EnergyBalanceType.enthalpyTotal
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
        m.fs = FlowsheetBlock(default={"dynamic": True, "time_set": [0, 3600]})

        # Set up thermo props and reaction props
        m.fs.gas_props = GasPhaseParameterBlock()
        m.fs.solid_props = SolidPhaseParameterBlock()
        m.fs.solid_rxns = HeteroReactionParameterBlock(
                default={"solid_property_package": m.fs.solid_props,
                         "gas_property_package": m.fs.gas_props})

        m.fs.unit = FixedBed0D(
            default={"gas_property_package": m.fs.gas_props,
                     "solid_property_package": m.fs.solid_props,
                     "reaction_package": m.fs.solid_rxns})

        # Discretize time domain
        m.discretizer = TransformationFactory('dae.finite_difference')
        m.discretizer.apply_to(m,
                               nfe=100,
                               wrt=m.fs.time,
                               scheme="BACKWARD")

        # Set initial conditions of the solid phase
        m.fs.unit.mass_solids[0].fix(1)
        m.fs.unit.solids[0].particle_porosity.fix(0.20)
        m.fs.unit.solids[0].mass_frac_comp['Fe2O3'].fix(0.45)
        m.fs.unit.solids[0].mass_frac_comp['Fe3O4'].fix(0)
        m.fs.unit.solids[0].mass_frac_comp['Al2O3'].fix(0.55)
        m.fs.unit.solids[0].temperature.fix(1273.15)

        # Set conditions of the gas phase (this is all fixed as gas side
        # assumption is excess gas flowrate which means all state variables
        # remain unchanged)
        for t in m.fs.time:
            m.fs.unit.gas[t].temperature.fix(1273.15)
            m.fs.unit.gas[t].pressure.fix(1.01325)  # 1atm
            m.fs.unit.gas[t].mole_frac_comp['CO2'].fix(0.4)
            m.fs.unit.gas[t].mole_frac_comp['H2O'].fix(0.5)
            m.fs.unit.gas[t].mole_frac_comp['CH4'].fix(0.1)
            if m.fs.unit.config.energy_balance_type == EnergyBalanceType.none:
                m.fs.unit.solids[t].temperature.fix(1273.15)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, iron_oc):
        assert isinstance(iron_oc.fs.unit.solids_material_holdup_constraints,
                          Constraint)
        assert isinstance(iron_oc.fs.unit
                          .solids_material_accumulation_constraints,
                          Constraint)
        assert isinstance(iron_oc.fs.unit.mass_solids_eq, Constraint)
        assert isinstance(iron_oc.fs.unit.sum_component_eqn, Constraint)
        assert isinstance(iron_oc.fs.unit.solids_energy_holdup_constraints,
                          Constraint)
        assert isinstance(iron_oc.fs.unit
                          .solids_energy_accumulation_constraints, Constraint)

        assert number_variables(iron_oc) == 3546
        assert number_total_constraints(iron_oc) == 2823
        assert number_unused_variables(iron_oc) == 206

    @pytest.mark.unit
    def test_dof(self, iron_oc):
        assert degrees_of_freedom(iron_oc) == 0

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, iron_oc):
        initialization_tester(iron_oc)

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
    def test_solution(self, iron_oc):  # need to update these values
        assert (pytest.approx(0.9851, abs=1e-2) ==
                iron_oc.fs.unit.mass_solids[3600].value)
        assert (pytest.approx(0.1955, abs=1e-2) ==
                iron_oc.fs.unit.solids[3600].particle_porosity.value)
        assert (pytest.approx(0.5583, abs=1e-2) ==
                iron_oc.fs.unit.solids[3600].mass_frac_comp['Al2O3'].value)
        assert (pytest.approx(0.0051, abs=1e-2) ==
                iron_oc.fs.unit.solids[3600].mass_frac_comp['Fe2O3'].value)
        assert (pytest.approx(0.4366, abs=1e-2) ==
                iron_oc.fs.unit.solids[3600].mass_frac_comp['Fe3O4'].value)
        assert (pytest.approx(1255.59, abs=1e-2) ==
                iron_oc.fs.unit.solids[3600].temperature.value)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, iron_oc):
        # Conservation of material check
        # first, check if Final Mass - Initial Mass = Total Holdup Change
        mbal_solid_total = value(
                (iron_oc.fs.unit.mass_solids[3600] -
                 iron_oc.fs.unit.mass_solids[0]) -
                (sum(iron_oc.fs.unit.solids_material_holdup[3600, j]
                     for j in ('Al2O3', 'Fe2O3', 'Fe3O4')) -
                 sum(iron_oc.fs.unit.solids_material_holdup[0, j]
                     for j in ('Al2O3', 'Fe2O3', 'Fe3O4'))))
        assert abs(mbal_solid_total) <= 1e-2

        # second, check if Change in Mass = Change in Holdup at each time t
        dt = 36  # time step in seconds, used in time point discretization
        mbal_solid = []
        for t in iron_oc.fs.time:
            if t == 0:
                continue
            else:
                mbal_solid.append(value(
                    (iron_oc.fs.unit.mass_solids[t] -
                     iron_oc.fs.unit.mass_solids[t - dt]) -
                    (sum(iron_oc.fs.unit.solids_material_holdup[t, j]
                         for j in ('Al2O3', 'Fe2O3', 'Fe3O4')) -
                     sum(iron_oc.fs.unit.solids_material_holdup[t - dt, j]
                         for j in ('Al2O3', 'Fe2O3', 'Fe3O4')))))
        for val in mbal_solid:
            assert abs(val) <= 1e-2

        # Reaction stoichiometric ratio (gas is fixed, can only check solids)
        # Overall reducer reactions for methane combustion:
        # CH4 + 12Fe2O3 => 8Fe3O4 + CO2 + 2H2O
        fe2o3_reacted = value(
            (iron_oc.fs.unit.mass_solids[0] *
             iron_oc.fs.unit.solids[0].mass_frac_comp['Fe2O3'] -
             iron_oc.fs.unit.mass_solids[3600] *
             iron_oc.fs.unit.solids[3600].mass_frac_comp['Fe2O3']) /
            iron_oc.fs.unit.solids[0]._params.mw_comp['Fe2O3'])
        fe3o4_produced = value(
            (iron_oc.fs.unit.mass_solids[3600] *
             iron_oc.fs.unit.solids[3600].mass_frac_comp['Fe3O4'] -
             iron_oc.fs.unit.mass_solids[0] *
             iron_oc.fs.unit.solids[0].mass_frac_comp['Fe3O4']) /
            iron_oc.fs.unit.solids[0]._params.mw_comp['Fe3O4'])
        stoichiometric_ratio = fe2o3_reacted/fe3o4_produced
        assert (pytest.approx(1.5, abs=1e-6) == stoichiometric_ratio)

        # Conservation of energy check
        # first, check if Final Energy - Initial Energy = Total Holdup Change
        ebal_solid_total = value(
                (iron_oc.fs.unit.mass_solids[3600] *
                 iron_oc.fs.unit.solids[3600].enth_mass -
                 iron_oc.fs.unit.mass_solids[0] *
                 iron_oc.fs.unit.solids[0].enth_mass) -
                (iron_oc.fs.unit.solids_energy_holdup[3600] -
                 iron_oc.fs.unit.solids_energy_holdup[0]))
        assert abs(ebal_solid_total) <= 1e-2

        # second, check if Change in Energy = Change in Holdup at each time t
        dt = 36  # seconds, used in time point discretization
        ebal_solid = []
        for t in iron_oc.fs.time:
            if t == 0:
                continue
            else:
                ebal_solid.append(value(
                    (iron_oc.fs.unit.mass_solids[t] *
                     iron_oc.fs.unit.solids[t].enth_mass -
                     iron_oc.fs.unit.mass_solids[t - dt] *
                     iron_oc.fs.unit.solids[t - dt].enth_mass) -
                    (iron_oc.fs.unit.solids_energy_holdup[t] -
                     iron_oc.fs.unit.solids_energy_holdup[t - dt])))
        for val in ebal_solid:
            assert abs(val) <= 1e-2

    @pytest.mark.ui
    @pytest.mark.unit
    def test_report(self, iron_oc):
        # no inlet/outlet ports (no stream table), test block reports instead
        iron_oc.fs.unit.gas.report()
        iron_oc.fs.unit.solids.report()


# -----------------------------------------------------------------------------
class TestIronOC_EnergyBalanceType(object):
    @pytest.fixture(scope="class")
    def iron_oc(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": True, "time_set": [0, 3600]})

        # Set up thermo props and reaction props
        m.fs.gas_properties = GasPhaseParameterBlock()
        m.fs.solid_properties = SolidPhaseParameterBlock()
        m.fs.hetero_reactions = HeteroReactionParameterBlock(
                default={"solid_property_package": m.fs.solid_properties,
                         "gas_property_package": m.fs.gas_properties})

        m.fs.unit = FixedBed0D(
                default={"energy_balance_type": EnergyBalanceType.none,
                         "gas_property_package": m.fs.gas_properties,
                         "solid_property_package": m.fs.solid_properties,
                         "reaction_package": m.fs.hetero_reactions})

        # Discretize time domain
        m.discretizer = TransformationFactory('dae.finite_difference')
        m.discretizer.apply_to(m,
                               nfe=100,
                               wrt=m.fs.time,
                               scheme="BACKWARD")

        # Set initial conditions of the solid phase
        m.fs.unit.mass_solids[0].fix(1)
        m.fs.unit.solids[0].particle_porosity.fix(0.20)
        m.fs.unit.solids[0].mass_frac_comp['Fe2O3'].fix(0.45)
        m.fs.unit.solids[0].mass_frac_comp['Fe3O4'].fix(0)
        m.fs.unit.solids[0].mass_frac_comp['Al2O3'].fix(0.55)
        m.fs.unit.solids[0].temperature.fix(1273.15)

        # Set conditions of the gas phase (this is all fixed as gas side
        # assumption is excess gas flowrate which means all state variables
        # remain unchanged)
        for t in m.fs.time:
            m.fs.unit.gas[t].temperature.fix(1273.15)
            m.fs.unit.gas[t].pressure.fix(1.01325)  # 1atm
            m.fs.unit.gas[t].mole_frac_comp['CO2'].fix(0.4)
            m.fs.unit.gas[t].mole_frac_comp['H2O'].fix(0.5)
            m.fs.unit.gas[t].mole_frac_comp['CH4'].fix(0.1)
            if m.fs.unit.config.energy_balance_type == EnergyBalanceType.none:
                m.fs.unit.solids[t].temperature.fix(1273.15)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, iron_oc):
        assert isinstance(iron_oc.fs.unit.solids_material_holdup_constraints,
                          Constraint)
        assert isinstance(iron_oc.fs.unit
                          .solids_material_accumulation_constraints,
                          Constraint)
        assert isinstance(iron_oc.fs.unit.mass_solids_eq, Constraint)
        assert isinstance(iron_oc.fs.unit.sum_component_eqn, Constraint)

        assert number_variables(iron_oc) == 2940
        assert number_total_constraints(iron_oc) == 2117
        assert number_unused_variables(iron_oc) == 206
        print(unused_variables_set(iron_oc))

    @pytest.mark.unit
    def test_dof(self, iron_oc):
        assert degrees_of_freedom(iron_oc) == 0

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, iron_oc):
        initialization_tester(iron_oc)

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
    def test_solution(self, iron_oc):  # need to update these values
        assert (pytest.approx(0.9851, abs=1e-2) ==
                iron_oc.fs.unit.mass_solids[3600].value)
        assert (pytest.approx(0.1955, abs=1e-2) ==
                iron_oc.fs.unit.solids[3600].particle_porosity.value)
        assert (pytest.approx(0.5583, abs=1e-2) ==
                iron_oc.fs.unit.solids[3600].mass_frac_comp['Al2O3'].value)
        assert (pytest.approx(0.0047, abs=1e-2) ==
                iron_oc.fs.unit.solids[3600].mass_frac_comp['Fe2O3'].value)
        assert (pytest.approx(0.4370, abs=1e-2) ==
                iron_oc.fs.unit.solids[3600].mass_frac_comp['Fe3O4'].value)
        assert (pytest.approx(1273.15, abs=1e-2) ==
                iron_oc.fs.unit.solids[3600].temperature.value)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, iron_oc):
        # Conservation of material check
        # first, check if Final Mass - Initial Mass = Total Holdup Change
        mbal_solid_total = value(
                (iron_oc.fs.unit.mass_solids[3600] -
                 iron_oc.fs.unit.mass_solids[0]) -
                (sum(iron_oc.fs.unit.solids_material_holdup[3600, j]
                     for j in ('Al2O3', 'Fe2O3', 'Fe3O4')) -
                 sum(iron_oc.fs.unit.solids_material_holdup[0, j]
                     for j in ('Al2O3', 'Fe2O3', 'Fe3O4'))))
        assert abs(mbal_solid_total) <= 1e-2

        # second, check if Change in Mass = Change in Holdup at each time t
        dt = 36  # time step in seconds, used in time point discretization
        mbal_solid = []
        for t in iron_oc.fs.time:
            if t == 0:
                continue
            else:
                mbal_solid.append(value(
                    (iron_oc.fs.unit.mass_solids[t] -
                     iron_oc.fs.unit.mass_solids[t - dt]) -
                    (sum(iron_oc.fs.unit.solids_material_holdup[t, j]
                         for j in ('Al2O3', 'Fe2O3', 'Fe3O4')) -
                     sum(iron_oc.fs.unit.solids_material_holdup[t - dt, j]
                         for j in ('Al2O3', 'Fe2O3', 'Fe3O4')))))
        for val in mbal_solid:
            assert abs(val) <= 1e-2

        # Reaction stoichiometric ratio (gas is fixed, can only check solids)
        # Overall reducer reactions for methane combustion:
        # CH4 + 12Fe2O3 => 8Fe3O4 + CO2 + 2H2O
        fe2o3_reacted = value(
            (iron_oc.fs.unit.mass_solids[0] *
             iron_oc.fs.unit.solids[0].mass_frac_comp['Fe2O3'] -
             iron_oc.fs.unit.mass_solids[3600] *
             iron_oc.fs.unit.solids[3600].mass_frac_comp['Fe2O3']) /
            iron_oc.fs.unit.solids[0]._params.mw_comp['Fe2O3'])
        fe3o4_produced = value(
            (iron_oc.fs.unit.mass_solids[3600] *
             iron_oc.fs.unit.solids[3600].mass_frac_comp['Fe3O4'] -
             iron_oc.fs.unit.mass_solids[0] *
             iron_oc.fs.unit.solids[0].mass_frac_comp['Fe3O4']) /
            iron_oc.fs.unit.solids[0]._params.mw_comp['Fe3O4'])
        stoichiometric_ratio = fe2o3_reacted/fe3o4_produced
        assert (pytest.approx(1.5, abs=1e-6) == stoichiometric_ratio)

    @pytest.mark.ui
    @pytest.mark.unit
    def test_report(self, iron_oc):
        # no inlet/outlet ports (no stream table), test block reports instead
        iron_oc.fs.unit.gas.report()
        iron_oc.fs.unit.solids.report()
