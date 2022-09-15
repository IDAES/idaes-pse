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
Tests for CLC heterogeneous reaction block; tests for construction and solve
Author: Chinedu Okoli
"""

import pytest

from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    Var,
    Constraint,
    value,
    units as pyunits,
)
from pyomo.util.check_units import assert_units_consistent
import pyomo.common.unittest as unittest
from pyomo.common.collections import ComponentMap
from pyomo.util.subsystems import ParamSweeper
from pyomo.contrib.incidence_analysis import (
    solve_strongly_connected_components,
)

from idaes.core import FlowsheetBlock

from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_large_residuals,
    fixed_variables_set,
    activated_constraints_set,
)

from idaes.core.solvers import get_solver

from idaes.models_extra.gas_solid_contactors.properties.oxygen_iron_OC_oxidation.gas_phase_thermo import (
    GasPhaseParameterBlock,
)
from idaes.models_extra.gas_solid_contactors.properties.oxygen_iron_OC_oxidation.solid_phase_thermo import (
    SolidPhaseParameterBlock,
)
from idaes.models_extra.gas_solid_contactors.properties.oxygen_iron_OC_oxidation.hetero_reactions import (
    HeteroReactionParameterBlock,
)


# Get default solver for testing
solver = get_solver()


# -----------------------------------------------------------------------------
@pytest.fixture(scope="class")
def rxn_prop():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    # Set up thermo props and reaction props
    m.fs.solid_properties = SolidPhaseParameterBlock()
    m.fs.solid_state_block = m.fs.solid_properties.build_state_block(
        [0], parameters=m.fs.solid_properties, defined_state=True
    )

    m.fs.gas_properties = GasPhaseParameterBlock()
    m.fs.gas_state_block = m.fs.gas_properties.build_state_block(
        [0], parameters=m.fs.gas_properties, defined_state=True
    )

    m.fs.reactions = HeteroReactionParameterBlock(
        solid_property_package=m.fs.solid_properties,
        gas_property_package=m.fs.gas_properties,
    )
    m.fs.unit = m.fs.reactions.reaction_block_class(
        [0],
        parameters=m.fs.reactions,
        solid_state_block=m.fs.solid_state_block,
        gas_state_block=m.fs.gas_state_block,
        has_equilibrium=False,
    )

    # Fix required variables to make reaction model square
    # (gas mixture and component densities,
    # solid particle porosity, density and component fractions)
    m.fs.gas_state_block[0].dens_mol.fix(10)
    m.fs.gas_state_block[0].dens_mol_comp.fix(10)
    m.fs.solid_state_block[0].particle_porosity.fix(0.23)
    m.fs.solid_state_block[0].mass_frac_comp["Fe2O3"].fix(0.244)
    m.fs.solid_state_block[0].mass_frac_comp["Fe3O4"].fix(0.202)
    m.fs.solid_state_block[0].mass_frac_comp["Al2O3"].fix(0.554)
    m.fs.solid_state_block[0].dens_mass_skeletal.fix(3125)
    m.fs.solid_state_block[0].temperature.fix(1183.15)  # K

    return m


@pytest.mark.unit
def test_build_reaction_block(rxn_prop):
    assert isinstance(rxn_prop.fs.unit[0].k_rxn, Var)
    assert isinstance(rxn_prop.fs.unit[0].OC_conv, Var)
    assert isinstance(rxn_prop.fs.unit[0].OC_conv_temp, Var)
    assert isinstance(rxn_prop.fs.unit[0].reaction_rate, Var)


@pytest.mark.unit
def test_setInputs_reaction_block(rxn_prop):
    assert degrees_of_freedom(rxn_prop.fs.unit[0]) == 0


@pytest.mark.solver
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_initialize(rxn_prop):
    orig_fixed_vars = fixed_variables_set(rxn_prop)
    orig_act_consts = activated_constraints_set(rxn_prop)

    rxn_prop.fs.unit.initialize()

    assert degrees_of_freedom(rxn_prop) == 0

    fin_fixed_vars = fixed_variables_set(rxn_prop)
    fin_act_consts = activated_constraints_set(rxn_prop)

    assert len(fin_act_consts) == len(orig_act_consts)
    assert len(fin_fixed_vars) == len(orig_fixed_vars)

    for c in fin_act_consts:
        assert c in orig_act_consts
    for v in fin_fixed_vars:
        assert v in orig_fixed_vars


@pytest.mark.solver
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_solve(rxn_prop):

    assert hasattr(rxn_prop.fs.unit[0], "k_rxn")
    assert hasattr(rxn_prop.fs.unit[0], "OC_conv")
    assert hasattr(rxn_prop.fs.unit[0], "reaction_rate")

    results = solver.solve(rxn_prop.fs.unit[0])

    # Check for optimal solution
    assert check_optimal_termination(results)


@pytest.mark.solver
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_solution(rxn_prop):
    assert pytest.approx(1, rel=1e-5) == rxn_prop.fs.unit[0].k_rxn["R1"].value
    assert pytest.approx(0, rel=1e-5) == rxn_prop.fs.unit[0].OC_conv.value
    assert pytest.approx(0, rel=1e-5) == rxn_prop.fs.unit[0].reaction_rate["R1"].value


@pytest.mark.component
def test_units_consistent(rxn_prop):
    assert_units_consistent(rxn_prop)


@pytest.mark.unit
def test_state_vars():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.solid_properties = SolidPhaseParameterBlock()
    m.fs.solid_state = m.fs.solid_properties.build_state_block(
        [0], parameters=m.fs.solid_properties, defined_state=True
    )

    m.fs.gas_properties = GasPhaseParameterBlock()
    m.fs.gas_state = m.fs.gas_properties.build_state_block(
        [0], parameters=m.fs.gas_properties, defined_state=True
    )

    m.fs.reaction_properties = HeteroReactionParameterBlock(
        solid_property_package=m.fs.solid_properties,
        gas_property_package=m.fs.gas_properties,
    )
    m.fs.reaction_block = m.fs.reaction_properties.reaction_block_class(
        [0],
        parameters=m.fs.reaction_properties,
        solid_state_block=m.fs.solid_state,
        gas_state_block=m.fs.gas_state,
        has_equilibrium=False,
    )

    for var in m.fs.gas_state[0].define_state_vars().values():
        var.fix()
    for var in m.fs.solid_state[0].define_state_vars().values():
        var.fix()

    # Note that these checks are necessary to trigger the construction of
    # the reaction block variables
    assert isinstance(m.fs.reaction_block[0].k_rxn, Var)
    assert isinstance(m.fs.reaction_block[0].OC_conv, Var)
    assert isinstance(m.fs.reaction_block[0].OC_conv_temp, Var)
    assert isinstance(m.fs.reaction_block[0].reaction_rate, Var)

    assert degrees_of_freedom(m) == 0

    rxn_vars = list(m.fs.reaction_block[0].component_data_objects(Var))
    rxn_cons = list(m.fs.reaction_block[0].component_data_objects(Constraint))
    assert len(rxn_vars) == 4
    assert len(rxn_cons) == 4


@pytest.mark.unit
class TestProperties(unittest.TestCase):
    def _make_model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.solid_properties = SolidPhaseParameterBlock()
        m.fs.solid_state = m.fs.solid_properties.build_state_block(
            parameters=m.fs.solid_properties, defined_state=True
        )

        m.fs.gas_properties = GasPhaseParameterBlock()
        m.fs.gas_state = m.fs.gas_properties.build_state_block(
            parameters=m.fs.gas_properties, defined_state=True
        )

        m.fs.reaction_properties = HeteroReactionParameterBlock(
            solid_property_package=m.fs.solid_properties,
            gas_property_package=m.fs.gas_properties,
        )
        m.fs.reaction_block = m.fs.reaction_properties.reaction_block_class(
            parameters=m.fs.reaction_properties,
            solid_state_block=m.fs.solid_state,
            gas_state_block=m.fs.gas_state,
            has_equilibrium=False,
        )

        for var in m.fs.gas_state.define_state_vars().values():
            var.fix()
        for var in m.fs.solid_state.define_state_vars().values():
            var.fix()
        return m

    def test_k_rxn(self):
        m = self._make_model()
        rxn_block = m.fs.reaction_block

        n_scen = 4
        mols = pyunits.mol / pyunits.s
        kgs = pyunits.kg / pyunits.s
        K = pyunits.K
        bar = pyunits.bar
        state_values = {
            "gas_state.flow_mol": [1.0 * mols] * n_scen,
            "solid_state.flow_mass": [1.0 * kgs] * n_scen,
            "gas_state.temperature": [1000.0 * K, 1100.0 * K, 1200.0 * K, 1300.0 * K],
            "solid_state.temperature": [1000.0 * K, 1100.0 * K, 1200.0 * K, 1300.0 * K],
            "gas_state.pressure": [1.0 * bar] * n_scen,
            "solid_state.particle_porosity": [0.27] * n_scen,
            "gas_state.mole_frac_comp[O2]": [0.25] * n_scen,
            "gas_state.mole_frac_comp[N2]": [0.25] * n_scen,
            "gas_state.mole_frac_comp[H2O]": [0.25] * n_scen,
            "gas_state.mole_frac_comp[CO2]": [0.25] * n_scen,
            "solid_state.mass_frac_comp[Fe2O3]": [1.0 / 3.0] * n_scen,
            "solid_state.mass_frac_comp[Fe3O4]": [1.0 / 3.0] * n_scen,
            "solid_state.mass_frac_comp[Al2O3]": [1.0 / 3.0] * n_scen,
        }
        state_values = ComponentMap(
            (m.fs.find_component(name), values) for name, values in state_values.items()
        )

        # Units of k_rxn are "non-physical" si units, chosen to be
        # consistent with the reaction rate rule.
        k_rxn_units = pyunits.m / pyunits.s
        target_values = {
            "reaction_block.k_rxn[R1]": [
                5.7556e-5 * k_rxn_units,
                6.7076e-5 * k_rxn_units,
                7.6203e-5 * k_rxn_units,
                8.4888e-5 * k_rxn_units,
            ],
        }
        target_values = ComponentMap(
            (m.fs.find_component(name), values)
            for name, values in target_values.items()
        )

        assert degrees_of_freedom(m.fs) == 0

        param_sweeper = ParamSweeper(n_scen, state_values, output_values=target_values)
        with param_sweeper:
            for inputs, outputs in param_sweeper:
                solve_strongly_connected_components(m.fs)

                # Make sure property equalites have been converged
                assert number_large_residuals(m.fs, tol=1e-8) == 0

                # Sanity checks that inputs are properly set
                for var, val in inputs.items():
                    val = value(pyunits.convert(val, var.get_units()))
                    assert var.value == pytest.approx(value(val), rel=1e-3)

                # Make sure properties have been calculated as expected
                for var, val in outputs.items():
                    val = value(pyunits.convert(val, var.get_units()))
                    assert var.value == pytest.approx(value(val), rel=1e-3)

    def test_oc_conv(self):
        m = self._make_model()
        rxn_block = m.fs.reaction_block

        n_scen = 3
        mols = pyunits.mol / pyunits.s
        kgs = pyunits.kg / pyunits.s
        K = pyunits.K
        bar = pyunits.bar
        state_values = {
            "gas_state.flow_mol": [1.0 * mols] * n_scen,
            "solid_state.flow_mass": [1.0 * kgs] * n_scen,
            "gas_state.temperature": [1273.0 * K] * n_scen,
            "solid_state.temperature": [1273.0 * K] * n_scen,
            "gas_state.pressure": [1.0 * bar] * n_scen,
            "solid_state.particle_porosity": [0.27] * n_scen,
            "gas_state.mole_frac_comp[O2]": [0.25] * n_scen,
            "gas_state.mole_frac_comp[N2]": [0.25] * n_scen,
            "gas_state.mole_frac_comp[H2O]": [0.25] * n_scen,
            "gas_state.mole_frac_comp[CO2]": [0.25] * n_scen,
            "solid_state.mass_frac_comp[Fe2O3]": [2 / 3, 0.0, 1 / 3],
            "solid_state.mass_frac_comp[Fe3O4]": [0.0, 2 / 3, 1 / 3],
            "solid_state.mass_frac_comp[Al2O3]": [1 / 3] * n_scen,
        }
        state_values = ComponentMap(
            (m.fs.find_component(name), values) for name, values in state_values.items()
        )

        target_values = {
            "reaction_block.OC_conv": [
                1.0,
                0.0,
                0.4915,
            ],
            "reaction_block.OC_conv_temp": [
                0.000,
                1.0,
                0.6371,
            ],
        }
        target_values = ComponentMap(
            (m.fs.find_component(name), values)
            for name, values in target_values.items()
        )

        assert degrees_of_freedom(m.fs) == 0

        param_sweeper = ParamSweeper(n_scen, state_values, output_values=target_values)
        with param_sweeper:
            for inputs, outputs in param_sweeper:
                solve_strongly_connected_components(m.fs)

                # Make sure property equalites have been converged
                assert number_large_residuals(m.fs, tol=1e-8) == 0

                # Sanity checks that inputs are properly set
                for var, val in inputs.items():
                    val = value(pyunits.convert(val, var.get_units()))
                    assert var.value == pytest.approx(value(val), rel=1e-3)

                # Make sure properties have been calculated as expected
                for var, val in outputs.items():
                    val = value(pyunits.convert(val, var.get_units()))
                    assert var.value == pytest.approx(value(val), rel=1e-3)

    def test_reaction_rate(self):
        m = self._make_model()
        rxn_block = m.fs.reaction_block

        n_scen = 9
        mols = pyunits.mol / pyunits.s
        kgs = pyunits.kg / pyunits.s
        K = pyunits.K
        bar = pyunits.bar
        state_values = {
            "gas_state.flow_mol": [1.0 * mols] * n_scen,
            "solid_state.flow_mass": [1.0 * kgs] * n_scen,
            "gas_state.temperature": [1273.0 * K] * n_scen,
            "solid_state.temperature": [
                1273.0 * K,
                1273.0 * K,
                1273.0 * K,
                1273.0 * K,
                1273.0 * K,
                1273.0 * K,
                1100.0 * K,
                1200.0 * K,
                1300.0 * K,
            ],
            "gas_state.pressure": [1.0 * bar] * n_scen,
            "solid_state.particle_porosity": [0.27] * n_scen,
            "gas_state.mole_frac_comp[O2]": [
                1.0,
                0.7,
                0.0,
                0.25,
                0.25,
                0.25,
                0.25,
                0.25,
                0.25,
            ],
            "gas_state.mole_frac_comp[N2]": [
                0.0,
                0.1,
                1 / 3,
                0.25,
                0.25,
                0.25,
                0.25,
                0.25,
                0.25,
            ],
            "gas_state.mole_frac_comp[H2O]": [
                0.0,
                0.1,
                1 / 3,
                0.25,
                0.25,
                0.25,
                0.25,
                0.25,
                0.25,
            ],
            "gas_state.mole_frac_comp[CO2]": [
                0.0,
                0.1,
                1 / 3,
                0.25,
                0.25,
                0.25,
                0.25,
                0.25,
                0.25,
            ],
            "solid_state.mass_frac_comp[Fe2O3]": [
                1 / 3,
                1 / 3,
                1 / 3,
                2 / 3,
                0.0,
                1 / 3,
                1 / 3,
                1 / 3,
                1 / 3,
            ],
            "solid_state.mass_frac_comp[Fe3O4]": [
                1 / 3,
                1 / 3,
                1 / 3,
                0.0,
                2 / 3,
                1 / 3,
                1 / 3,
                1 / 3,
                1 / 3,
            ],
            "solid_state.mass_frac_comp[Al2O3]": [1 / 3] * n_scen,
        }
        state_values = ComponentMap(
            (m.fs.find_component(name), values) for name, values in state_values.items()
        )

        molm3s = pyunits.mol / pyunits.m**3 / pyunits.s
        target_values = {
            "reaction_block.reaction_rate[R1]": [
                351.367 * molm3s,
                245.957 * molm3s,
                3.7189815316196473e-07 * molm3s,
                0.000 * molm3s,
                271.731 * molm3s,
                87.842 * molm3s,
                71.344 * molm3s,
                81.051 * molm3s,
                90.288 * molm3s,
            ],
        }
        target_values = ComponentMap(
            (m.fs.find_component(name), values)
            for name, values in target_values.items()
        )

        assert degrees_of_freedom(m.fs) == 0

        param_sweeper = ParamSweeper(n_scen, state_values, output_values=target_values)
        with param_sweeper:
            for inputs, outputs in param_sweeper:
                solve_strongly_connected_components(m.fs)

                # Make sure property equalites have been converged
                assert number_large_residuals(m.fs, tol=1e-8) == 0

                # Sanity checks that inputs are properly set
                for var, val in inputs.items():
                    val = value(pyunits.convert(val, var.get_units()))
                    assert var.value == pytest.approx(value(val), rel=1e-3)

                # Make sure properties have been calculated as expected
                for var, val in outputs.items():
                    if value(val) != 0:
                        # To get around Pyomo issue #1627
                        val = value(pyunits.convert(val, var.get_units()))
                    assert var.value == pytest.approx(value(val), rel=1e-3)


if __name__ == "__main__":
    unittest.main()
