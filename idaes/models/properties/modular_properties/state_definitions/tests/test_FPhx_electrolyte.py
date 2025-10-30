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
Tests for constructing and using component lists in electrolyte systems
"""
# Import Python libraries
from copy import deepcopy
import pytest

# Import Pyomo units
from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    Constraint,
    Set,
    value,
    Var,
    units as pyunits,
)

# Import IDAES cores
from idaes.models.properties.modular_properties.state_definitions import FPhx
from idaes.models.properties.modular_properties.state_definitions.tests.electrolyte_testing_config import (
    get_config_no_inherent_reactions,
    get_config_with_inherent_reactions,
)

from idaes.core import FlowsheetBlock
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
    StateIndex,
)

from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.solvers import get_solver

solver = get_solver("ipopt_v2")


def dummy_method(b, *args, **kwargs):
    return 42


state_bounds = {
    "flow_mol": (0, 100, 1000, pyunits.mol / pyunits.s),
    "enth_mol": (0, 10000, 1e8, pyunits.J / pyunits.mol),
    "temperature": (273.15, 300, 500, pyunits.K),
    "pressure": (5e4, 1e5, 1e6, pyunits.Pa),
}


# -----------------------------------------------------------------------------
class TestApparentSpeciesBasisNoInherent:
    config = get_config_no_inherent_reactions()
    config["state_components"] = StateIndex.apparent
    config["state_definition"] = FPhx
    config["state_bounds"] = deepcopy(state_bounds)

    @pytest.fixture(scope="class")
    def frame(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.props = GenericParameterBlock(**TestApparentSpeciesBasisNoInherent.config)

        m.fs.state = m.fs.props.build_state_block([1], defined_state=True)

        return m

    @pytest.mark.unit
    def test_vars_and_constraints(self, frame):
        m = frame

        assert m.fs.state[1].component_list is m.fs.props.apparent_species_set

        assert (
            m.fs.state[1].phase_component_set is m.fs.props.apparent_phase_component_set
        )

        assert isinstance(m.fs.state[1].flow_mol, Var)
        assert len(m.fs.state[1].flow_mol) == 1
        assert isinstance(m.fs.state[1].pressure, Var)
        assert len(m.fs.state[1].pressure) == 1
        assert isinstance(m.fs.state[1].temperature, Var)
        assert len(m.fs.state[1].temperature) == 1

        assert isinstance(m.fs.state[1].flow_mol_phase, Var)
        assert len(m.fs.state[1].flow_mol_phase) == 2
        for p in m.fs.state[1].flow_mol_phase:
            assert p in ["Liq", "Vap"]
        assert isinstance(m.fs.state[1].phase_frac, Var)
        assert len(m.fs.state[1].phase_frac) == 2
        for p in m.fs.state[1].phase_frac:
            assert p in ["Liq", "Vap"]

        assert isinstance(m.fs.state[1].mole_frac_comp, Var)
        assert len(m.fs.state[1].mole_frac_comp) == 4
        for j in m.fs.state[1].mole_frac_comp:
            assert j in ["KHCO3", "H2O", "CO2", "N2"]

        assert isinstance(m.fs.state[1].total_flow_balance, Constraint)
        assert len(m.fs.state[1].total_flow_balance) == 1
        assert isinstance(m.fs.state[1].phase_fraction_constraint, Constraint)
        assert len(m.fs.state[1].phase_fraction_constraint) == 2
        assert isinstance(m.fs.state[1].component_flow_balances, Constraint)
        assert len(m.fs.state[1].component_flow_balances) == 4
        for j in m.fs.state[1].component_flow_balances:
            assert j in ["KHCO3", "H2O", "CO2", "N2"]

        assert isinstance(m.fs.state[1].mole_frac_phase_comp, Var)
        assert len(m.fs.state[1].mole_frac_phase_comp) == 7
        for j in m.fs.state[1].mole_frac_phase_comp:
            assert j in [
                ("Liq", "H2O"),
                ("Liq", "CO2"),
                ("Liq", "KHCO3"),
                ("Vap", "H2O"),
                ("Vap", "CO2"),
                ("Vap", "KHCO3"),
                ("Vap", "N2"),
            ]

        # Check references to base state variables
        assert m.fs.state[1].flow_mol_apparent is m.fs.state[1].flow_mol
        for i in m.fs.state[1].flow_mol_phase_apparent:
            assert (
                m.fs.state[1].flow_mol_phase_apparent[i]
                is m.fs.state[1].flow_mol_phase[i]
            )
            assert i in m.fs.props.phase_list
        for i in m.fs.state[1].flow_mol_phase_comp_apparent:
            assert (
                m.fs.state[1].flow_mol_phase_comp_apparent[i]
                is m.fs.state[1].flow_mol_phase_comp[i]
            )
            assert i in m.fs.props.apparent_phase_component_set
        for i in m.fs.state[1].mole_frac_phase_comp_apparent:
            assert (
                m.fs.state[1].mole_frac_phase_comp_apparent[i]
                is m.fs.state[1].mole_frac_phase_comp[i]
            )
            assert i in m.fs.props.apparent_phase_component_set

        # Check for true species components
        assert isinstance(m.fs.state[1].flow_mol_phase_comp_true, Var)
        assert len(m.fs.state[1].flow_mol_phase_comp_true) == 8
        assert isinstance(m.fs.state[1].appr_to_true_species, Constraint)
        assert len(m.fs.state[1].appr_to_true_species) == 8

    @pytest.mark.component
    def test_solve_for_true_species(self, frame):
        m = frame

        m.fs.state[1].flow_mol.fix(2)
        m.fs.state[1].phase_frac["Liq"].fix(0.5)
        m.fs.state[1].temperature.fix(300)
        m.fs.state[1].pressure.fix(1e5)

        m.fs.state[1].mole_frac_comp["H2O"].fix(0.25)
        m.fs.state[1].mole_frac_comp["CO2"].fix(0.25)
        m.fs.state[1].mole_frac_comp["KHCO3"].fix(0.25)
        m.fs.state[1].mole_frac_comp["N2"].fix(0.25)

        m.fs.state[1].mole_frac_phase_comp["Vap", "KHCO3"].fix(1 / 6)
        m.fs.state[1].mole_frac_phase_comp["Vap", "CO2"].fix(1 / 6)

        assert degrees_of_freedom(m.fs) == 0

        res = solver.solve(m.fs, tee=True)

        # Check for optimal solution
        assert check_optimal_termination(res)

        # Check true species flowrates
        assert value(
            m.fs.state[1].flow_mol_phase_comp_true["Vap", "H2O"]
        ) == pytest.approx(1 / 6, rel=1e-5)
        assert value(
            m.fs.state[1].flow_mol_phase_comp_true["Vap", "CO2"]
        ) == pytest.approx(1 / 6, rel=1e-5)
        assert value(
            m.fs.state[1].flow_mol_phase_comp_true["Vap", "KHCO3"]
        ) == pytest.approx(1 / 6, rel=1e-5)
        assert value(
            m.fs.state[1].flow_mol_phase_comp_true["Vap", "N2"]
        ) == pytest.approx(0.5, rel=1e-5)
        assert value(
            m.fs.state[1].flow_mol_phase_comp_true["Liq", "H2O"]
        ) == pytest.approx(1 / 3, rel=1e-5)
        assert value(
            m.fs.state[1].flow_mol_phase_comp_true["Liq", "CO2"]
        ) == pytest.approx(1 / 3, rel=1e-5)
        assert value(
            m.fs.state[1].flow_mol_phase_comp_true["Liq", "K+"]
        ) == pytest.approx(1 / 3, rel=1e-5)
        assert value(
            m.fs.state[1].flow_mol_phase_comp_true["Liq", "HCO3-"]
        ) == pytest.approx(1 / 3, rel=1e-5)


# -----------------------------------------------------------------------------
def dens_mol_H2O(*args, **kwargs):
    return 55e3


class TestApparentSpeciesBasisInherent:
    config = get_config_with_inherent_reactions()
    config["state_components"] = StateIndex.apparent
    config["state_definition"] = FPhx
    config["state_bounds"] = deepcopy(state_bounds)

    @pytest.fixture(scope="class")
    def frame(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.props = GenericParameterBlock(**TestApparentSpeciesBasisInherent.config)

        m.fs.state = m.fs.props.build_state_block([1], defined_state=True)

        m.fs.state[1].calculate_scaling_factors()

        return m

    @pytest.mark.unit
    def test_vars_and_constraints(self, frame):
        m = frame

        assert m.fs.state[1].component_list is m.fs.props.apparent_species_set

        assert (
            m.fs.state[1].phase_component_set is m.fs.props.apparent_phase_component_set
        )

        assert isinstance(m.fs.state[1].flow_mol, Var)
        assert len(m.fs.state[1].flow_mol) == 1
        assert isinstance(m.fs.state[1].pressure, Var)
        assert len(m.fs.state[1].pressure) == 1
        assert isinstance(m.fs.state[1].temperature, Var)
        assert len(m.fs.state[1].temperature) == 1

        assert isinstance(m.fs.state[1].flow_mol_phase, Var)
        assert len(m.fs.state[1].flow_mol_phase) == 1
        for p in m.fs.state[1].flow_mol_phase:
            assert p in ["Liq"]
        assert isinstance(m.fs.state[1].phase_frac, Var)
        assert len(m.fs.state[1].phase_frac) == 1
        for p in m.fs.state[1].phase_frac:
            assert p in ["Liq"]

        assert isinstance(m.fs.state[1].mole_frac_comp, Var)
        assert len(m.fs.state[1].mole_frac_comp) == 4
        for j in m.fs.state[1].mole_frac_comp:
            assert j in ["KHCO3", "H2O", "K2CO3", "KOH"]

        assert isinstance(m.fs.state[1].total_flow_balance, Constraint)
        assert len(m.fs.state[1].total_flow_balance) == 1
        assert isinstance(m.fs.state[1].phase_fraction_constraint, Constraint)
        assert len(m.fs.state[1].phase_fraction_constraint) == 1
        assert isinstance(m.fs.state[1].component_flow_balances, Constraint)
        assert len(m.fs.state[1].component_flow_balances) == 4
        for j in m.fs.state[1].component_flow_balances:
            assert j in ["KHCO3", "H2O", "K2CO3", "KOH"]

        assert isinstance(m.fs.state[1].mole_frac_phase_comp, Var)
        assert len(m.fs.state[1].mole_frac_phase_comp) == 4
        for j in m.fs.state[1].mole_frac_phase_comp:
            assert j in [
                ("Liq", "H2O"),
                ("Liq", "K2CO3"),
                ("Liq", "KHCO3"),
                ("Liq", "KOH"),
            ]

        # Check references to base state variables
        assert m.fs.state[1].flow_mol_apparent is m.fs.state[1].flow_mol
        for i in m.fs.state[1].flow_mol_phase_apparent:
            assert (
                m.fs.state[1].flow_mol_phase_apparent[i]
                is m.fs.state[1].flow_mol_phase[i]
            )
            assert i in m.fs.props.phase_list
        for i in m.fs.state[1].flow_mol_phase_comp_apparent:
            assert (
                m.fs.state[1].flow_mol_phase_comp_apparent[i]
                is m.fs.state[1].flow_mol_phase_comp[i]
            )
            assert i in m.fs.props.apparent_phase_component_set
        for i in m.fs.state[1].mole_frac_phase_comp_apparent:
            assert (
                m.fs.state[1].mole_frac_phase_comp_apparent[i]
                is m.fs.state[1].mole_frac_phase_comp[i]
            )
            assert i in m.fs.props.apparent_phase_component_set

        # Check for true species components
        assert isinstance(m.fs.state[1].flow_mol_phase_comp_true, Var)
        assert len(m.fs.state[1].flow_mol_phase_comp_true) == 6
        assert isinstance(m.fs.state[1].mole_frac_phase_comp_true, Var)
        assert len(m.fs.state[1].mole_frac_phase_comp_true) == 6
        assert isinstance(m.fs.state[1].appr_to_true_species, Constraint)
        assert len(m.fs.state[1].appr_to_true_species) == 6
        assert isinstance(m.fs.state[1].true_mole_frac_constraint, Constraint)
        assert len(m.fs.state[1].true_mole_frac_constraint) == 6

        # Check for inherent reactions
        assert m.fs.state[1].has_inherent_reactions
        assert not m.fs.state[1].include_inherent_reactions
        assert isinstance(m.fs.state[1].params.inherent_reaction_idx, Set)
        assert len(m.fs.state[1].params.inherent_reaction_idx) == 2
        for i in m.fs.state[1].params.inherent_reaction_idx:
            assert i in ["h2o_si", "co3_hco3"]

        assert isinstance(m.fs.state[1].apparent_inherent_reaction_extent, Var)
        assert len(m.fs.state[1].apparent_inherent_reaction_extent) == 2

    @pytest.mark.component
    def test_solve_for_true_species(self, frame):
        m = frame

        m.fs.state[1].flow_mol.fix(2)
        m.fs.state[1].enth_mol.fix(35000)
        m.fs.state[1].temperature.set_value(350)
        m.fs.state[1].pressure.fix(1e5)

        m.fs.state[1].mole_frac_comp["H2O"].fix(0.8)
        m.fs.state[1].mole_frac_comp["K2CO3"].fix(0.2)
        m.fs.state[1].mole_frac_comp["KHCO3"].fix(0)
        m.fs.state[1].mole_frac_comp["KOH"].fix(0)

        assert degrees_of_freedom(m.fs) == 0

        m.fs.state.initialize()

        res = solver.solve(m.fs)

        # Check for optimal solution
        assert check_optimal_termination(res)

        # Check apparent species flowrates
        for j in m.fs.state[1].mole_frac_comp:
            assert value(
                m.fs.state[1].flow_mol_phase_comp_apparent["Liq", j]
            ) == pytest.approx(
                value(m.fs.state[1].flow_mol * m.fs.state[1].mole_frac_comp[j]),
                rel=1e-5,
            )

        # Check element balances
        assert value(
            m.fs.state[1].flow_mol_phase_comp_apparent["Liq", "K2CO3"] * 2
            + m.fs.state[1].flow_mol_phase_comp_apparent["Liq", "KHCO3"]
            + m.fs.state[1].flow_mol_phase_comp_apparent["Liq", "KOH"]
        ) == pytest.approx(
            value(m.fs.state[1].flow_mol_phase_comp_true["Liq", "K+"]), rel=1e-5
        )
        assert value(
            m.fs.state[1].flow_mol_phase_comp_apparent["Liq", "K2CO3"]
            + m.fs.state[1].flow_mol_phase_comp_apparent["Liq", "KHCO3"]
        ) == pytest.approx(
            value(
                m.fs.state[1].flow_mol_phase_comp_true["Liq", "CO3--"]
                + m.fs.state[1].flow_mol_phase_comp_true["Liq", "HCO3-"]
            ),
            rel=1e-5,
        )
        assert value(
            m.fs.state[1].flow_mol_phase_comp_apparent["Liq", "K2CO3"] * 3
            + m.fs.state[1].flow_mol_phase_comp_apparent["Liq", "KHCO3"] * 3
            + m.fs.state[1].flow_mol_phase_comp_apparent["Liq", "KOH"]
            + m.fs.state[1].flow_mol_phase_comp_apparent["Liq", "H2O"]
        ) == pytest.approx(
            value(
                m.fs.state[1].flow_mol_phase_comp_true["Liq", "CO3--"] * 3
                + m.fs.state[1].flow_mol_phase_comp_true["Liq", "HCO3-"] * 3
                + m.fs.state[1].flow_mol_phase_comp_true["Liq", "OH-"]
                + m.fs.state[1].flow_mol_phase_comp_true["Liq", "H2O"]
            ),
            rel=1e-5,
        )
        assert value(
            m.fs.state[1].flow_mol_phase_comp_apparent["Liq", "KHCO3"]
            + m.fs.state[1].flow_mol_phase_comp_apparent["Liq", "KOH"]
            + m.fs.state[1].flow_mol_phase_comp_apparent["Liq", "H2O"] * 2
        ) == pytest.approx(
            value(
                m.fs.state[1].flow_mol_phase_comp_true["Liq", "HCO3-"]
                + m.fs.state[1].flow_mol_phase_comp_true["Liq", "OH-"]
                + m.fs.state[1].flow_mol_phase_comp_true["Liq", "H+"]
                + m.fs.state[1].flow_mol_phase_comp_true["Liq", "H2O"] * 2
            ),
            rel=1e-5,
        )

        # Check true species mole fractions
        assert value(
            m.fs.state[1].mole_frac_phase_comp_true["Liq", "CO3--"]
        ) == pytest.approx(0.142857, rel=1e-5)
        assert value(
            m.fs.state[1].mole_frac_phase_comp_true["Liq", "H+"]
        ) == pytest.approx(2.90081e-16, rel=1e-5)
        assert value(
            m.fs.state[1].mole_frac_phase_comp_true["Liq", "H2O"]
        ) == pytest.approx(0.571429, rel=1e-5)
        assert value(
            m.fs.state[1].mole_frac_phase_comp_true["Liq", "HCO3-"]
        ) == pytest.approx(1.13961e-08, rel=1e-5)
        assert value(
            m.fs.state[1].mole_frac_phase_comp_true["Liq", "K+"]
        ) == pytest.approx(0.285714, rel=1e-5)
        assert value(
            m.fs.state[1].mole_frac_phase_comp_true["Liq", "OH-"]
        ) == pytest.approx(1.139606e-08, rel=1e-5)


# -----------------------------------------------------------------------------
class TestTrueSpeciesBasisNoInherent:
    config = get_config_no_inherent_reactions()
    config["state_components"] = StateIndex.true
    config["state_definition"] = FPhx
    config["state_bounds"] = deepcopy(state_bounds)

    @pytest.fixture(scope="class")
    def frame(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.props = GenericParameterBlock(**TestTrueSpeciesBasisNoInherent.config)

        m.fs.state = m.fs.props.build_state_block([1], defined_state=True)

        return m

    @pytest.mark.unit
    def test_vars_and_constraints(self, frame):
        m = frame

        assert m.fs.state[1].component_list is m.fs.props.true_species_set

        assert m.fs.state[1].phase_component_set is m.fs.props.true_phase_component_set

        assert isinstance(m.fs.state[1].flow_mol, Var)
        assert len(m.fs.state[1].flow_mol) == 1
        assert isinstance(m.fs.state[1].pressure, Var)
        assert len(m.fs.state[1].pressure) == 1
        assert isinstance(m.fs.state[1].temperature, Var)
        assert len(m.fs.state[1].temperature) == 1

        assert isinstance(m.fs.state[1].flow_mol_phase, Var)
        assert len(m.fs.state[1].flow_mol_phase) == 2
        for p in m.fs.state[1].flow_mol_phase:
            assert p in ["Liq", "Vap"]
        assert isinstance(m.fs.state[1].phase_frac, Var)
        assert len(m.fs.state[1].phase_frac) == 2
        for p in m.fs.state[1].phase_frac:
            assert p in ["Liq", "Vap"]

        assert isinstance(m.fs.state[1].mole_frac_comp, Var)
        assert len(m.fs.state[1].mole_frac_comp) == 5
        for j in m.fs.state[1].mole_frac_comp:
            assert j in ["K+", "HCO3-", "H2O", "CO2", "N2"]

        assert isinstance(m.fs.state[1].total_flow_balance, Constraint)
        assert len(m.fs.state[1].total_flow_balance) == 1
        assert isinstance(m.fs.state[1].phase_fraction_constraint, Constraint)
        assert len(m.fs.state[1].phase_fraction_constraint) == 2
        assert isinstance(m.fs.state[1].component_flow_balances, Constraint)
        assert len(m.fs.state[1].component_flow_balances) == 5
        for j in m.fs.state[1].component_flow_balances:
            assert j in ["K+", "HCO3-", "H2O", "CO2", "N2"]

        assert isinstance(m.fs.state[1].mole_frac_phase_comp, Var)
        assert len(m.fs.state[1].mole_frac_phase_comp) == 8
        for j in m.fs.state[1].mole_frac_phase_comp:
            assert j in [
                ("Liq", "H2O"),
                ("Liq", "CO2"),
                ("Liq", "K+"),
                ("Liq", "HCO3-"),
                ("Vap", "H2O"),
                ("Vap", "CO2"),
                ("Vap", "KHCO3"),
                ("Vap", "N2"),
            ]

        # Check references to base state variables
        assert m.fs.state[1].flow_mol_true is m.fs.state[1].flow_mol
        for i in m.fs.state[1].flow_mol_phase_true:
            assert (
                m.fs.state[1].flow_mol_phase_true[i] is m.fs.state[1].flow_mol_phase[i]
            )
            assert i in m.fs.props.phase_list
        for i in m.fs.state[1].flow_mol_phase_comp_true:
            assert (
                m.fs.state[1].flow_mol_phase_comp_true[i]
                is m.fs.state[1].flow_mol_phase_comp[i]
            )
            assert i in m.fs.props.true_phase_component_set
        for i in m.fs.state[1].mole_frac_phase_comp_true:
            assert (
                m.fs.state[1].mole_frac_phase_comp_true[i]
                is m.fs.state[1].mole_frac_phase_comp[i]
            )
            assert i in m.fs.props.true_phase_component_set

        # Check for apparent species components
        assert isinstance(m.fs.state[1].flow_mol_phase_comp_apparent, Var)
        assert len(m.fs.state[1].flow_mol_phase_comp_apparent) == 7
        assert isinstance(m.fs.state[1].true_to_appr_species, Constraint)
        assert len(m.fs.state[1].true_to_appr_species) == 7


# -----------------------------------------------------------------------------
class TestTrueSpeciesBasisInherent:
    config = get_config_with_inherent_reactions()
    config["state_components"] = StateIndex.true
    config["state_definition"] = FPhx
    config["state_bounds"] = deepcopy(state_bounds)

    @pytest.fixture(scope="class")
    def frame(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.props = GenericParameterBlock(**TestTrueSpeciesBasisInherent.config)

        m.fs.state = m.fs.props.build_state_block([1], defined_state=True)

        m.fs.state[1].calculate_scaling_factors()

        return m

    @pytest.mark.unit
    def test_vars_and_constraints(self, frame):
        m = frame

        assert m.fs.state[1].component_list is m.fs.props.true_species_set

        assert m.fs.state[1].phase_component_set is m.fs.props.true_phase_component_set

        assert isinstance(m.fs.state[1].flow_mol, Var)
        assert len(m.fs.state[1].flow_mol) == 1
        assert isinstance(m.fs.state[1].pressure, Var)
        assert len(m.fs.state[1].pressure) == 1
        assert isinstance(m.fs.state[1].temperature, Var)
        assert len(m.fs.state[1].temperature) == 1

        assert isinstance(m.fs.state[1].flow_mol_phase, Var)
        assert len(m.fs.state[1].flow_mol_phase) == 1
        for p in m.fs.state[1].flow_mol_phase:
            assert p in ["Liq"]
        assert isinstance(m.fs.state[1].phase_frac, Var)
        assert len(m.fs.state[1].phase_frac) == 1
        for p in m.fs.state[1].phase_frac:
            assert p in ["Liq"]

        assert isinstance(m.fs.state[1].mole_frac_comp, Var)
        assert len(m.fs.state[1].mole_frac_comp) == 6
        for j in m.fs.state[1].mole_frac_comp:
            assert j in ["H+", "K+", "HCO3-", "CO3--", "OH-", "H2O"]

        assert isinstance(m.fs.state[1].total_flow_balance, Constraint)
        assert len(m.fs.state[1].total_flow_balance) == 1
        assert isinstance(m.fs.state[1].phase_fraction_constraint, Constraint)
        assert len(m.fs.state[1].phase_fraction_constraint) == 1
        assert isinstance(m.fs.state[1].component_flow_balances, Constraint)
        assert len(m.fs.state[1].component_flow_balances) == 6
        for j in m.fs.state[1].component_flow_balances:
            assert j in ["H+", "K+", "HCO3-", "CO3--", "OH-", "H2O"]

        assert isinstance(m.fs.state[1].mole_frac_phase_comp, Var)
        assert len(m.fs.state[1].mole_frac_phase_comp) == 6
        for j in m.fs.state[1].mole_frac_phase_comp:
            assert j in [
                ("Liq", "H+"),
                ("Liq", "K+"),
                ("Liq", "HCO3-"),
                ("Liq", "CO3--"),
                ("Liq", "OH-"),
                ("Liq", "H2O"),
            ]

        # Check references to base state variables
        assert m.fs.state[1].flow_mol_true is m.fs.state[1].flow_mol
        for i in m.fs.state[1].flow_mol_phase_true:
            assert (
                m.fs.state[1].flow_mol_phase_true[i] is m.fs.state[1].flow_mol_phase[i]
            )
            assert i in m.fs.props.phase_list
        for i in m.fs.state[1].flow_mol_phase_comp_true:
            assert (
                m.fs.state[1].flow_mol_phase_comp_true[i]
                is m.fs.state[1].flow_mol_phase_comp[i]
            )
            assert i in m.fs.props.true_phase_component_set
        for i in m.fs.state[1].mole_frac_phase_comp_true:
            assert (
                m.fs.state[1].mole_frac_phase_comp_true[i]
                is m.fs.state[1].mole_frac_phase_comp[i]
            )
            assert i in m.fs.props.true_phase_component_set

        # Check for apparent species components
        assert isinstance(m.fs.state[1].flow_mol_phase_comp_apparent, Var)
        assert len(m.fs.state[1].flow_mol_phase_comp_apparent) == 4
        assert isinstance(m.fs.state[1].mole_frac_phase_comp_apparent, Var)
        assert len(m.fs.state[1].mole_frac_phase_comp_apparent) == 4
        assert isinstance(m.fs.state[1].true_to_appr_species, Constraint)
        assert len(m.fs.state[1].true_to_appr_species) == 4
        assert isinstance(m.fs.state[1].appr_mole_frac_constraint, Constraint)
        assert len(m.fs.state[1].appr_mole_frac_constraint) == 4

        # Check for inherent reactions
        assert m.fs.state[1].has_inherent_reactions
        assert m.fs.state[1].include_inherent_reactions
        assert isinstance(m.fs.state[1].params.inherent_reaction_idx, Set)
        assert len(m.fs.state[1].params.inherent_reaction_idx) == 2
        for i in m.fs.state[1].params.inherent_reaction_idx:
            assert i in ["h2o_si", "co3_hco3"]

        assert not hasattr(m.fs.state[1], "apparent_inherent_reaction_extent")

    @pytest.mark.component
    def test_solve_for_true_species(self, frame):
        m = frame

        m.fs.state[1].flow_mol.fix(2)
        m.fs.state[1].enth_mol.fix(35000)
        m.fs.state[1].temperature.set_value(350)
        m.fs.state[1].pressure.fix(1e5)

        m.fs.state[1].mole_frac_comp["H2O"].fix(0.5716)
        m.fs.state[1].mole_frac_comp["K+"].fix(0.2858)
        m.fs.state[1].mole_frac_comp["H+"].fix(2.9e-16)
        m.fs.state[1].mole_frac_comp["HCO3-"].fix(1.14e-8)
        m.fs.state[1].mole_frac_comp["CO3--"].fix(0.1429)
        m.fs.state[1].mole_frac_comp["OH-"].fix(1.14e-8)

        assert degrees_of_freedom(m.fs) == 0

        m.fs.state.initialize()

        res = solver.solve(m.fs)

        # Check for optimal solution
        assert check_optimal_termination(res)

        assert m.fs.state[1].temperature.value == 350

        # Check true species flowrates
        for j in m.fs.state[1].mole_frac_comp:
            assert value(
                m.fs.state[1].flow_mol_phase_comp_true["Liq", j]
            ) == pytest.approx(
                value(m.fs.state[1].flow_mol * m.fs.state[1].mole_frac_comp[j]),
                rel=1e-5,
            )

        # Check element balances
        assert value(
            m.fs.state[1].flow_mol_phase_comp_apparent["Liq", "K2CO3"] * 2
            + m.fs.state[1].flow_mol_phase_comp_apparent["Liq", "KHCO3"]
            + m.fs.state[1].flow_mol_phase_comp_apparent["Liq", "KOH"]
        ) == pytest.approx(
            value(m.fs.state[1].flow_mol_phase_comp_true["Liq", "K+"]), rel=1e-5
        )
        assert value(
            m.fs.state[1].flow_mol_phase_comp_apparent["Liq", "K2CO3"]
            + m.fs.state[1].flow_mol_phase_comp_apparent["Liq", "KHCO3"]
        ) == pytest.approx(
            value(
                m.fs.state[1].flow_mol_phase_comp_true["Liq", "CO3--"]
                + m.fs.state[1].flow_mol_phase_comp_true["Liq", "HCO3-"]
            ),
            rel=1e-5,
        )
        assert value(
            m.fs.state[1].flow_mol_phase_comp_apparent["Liq", "K2CO3"] * 3
            + m.fs.state[1].flow_mol_phase_comp_apparent["Liq", "KHCO3"] * 3
            + m.fs.state[1].flow_mol_phase_comp_apparent["Liq", "KOH"]
            + m.fs.state[1].flow_mol_phase_comp_apparent["Liq", "H2O"]
        ) == pytest.approx(
            value(
                m.fs.state[1].flow_mol_phase_comp_true["Liq", "CO3--"] * 3
                + m.fs.state[1].flow_mol_phase_comp_true["Liq", "HCO3-"] * 3
                + m.fs.state[1].flow_mol_phase_comp_true["Liq", "OH-"]
                + m.fs.state[1].flow_mol_phase_comp_true["Liq", "H2O"]
            ),
            rel=1e-5,
        )
        assert value(
            m.fs.state[1].flow_mol_phase_comp_apparent["Liq", "KHCO3"]
            + m.fs.state[1].flow_mol_phase_comp_apparent["Liq", "KOH"]
            + m.fs.state[1].flow_mol_phase_comp_apparent["Liq", "H2O"] * 2
        ) == pytest.approx(
            value(
                m.fs.state[1].flow_mol_phase_comp_true["Liq", "HCO3-"]
                + m.fs.state[1].flow_mol_phase_comp_true["Liq", "OH-"]
                + m.fs.state[1].flow_mol_phase_comp_true["Liq", "H+"]
                + m.fs.state[1].flow_mol_phase_comp_true["Liq", "H2O"] * 2
            ),
            rel=1e-5,
        )

        # Check apparent species mole fractions
        assert value(
            m.fs.state[1].mole_frac_phase_comp_apparent["Liq", "K2CO3"]
        ) == pytest.approx(0.2, rel=1e-5)
        assert value(
            m.fs.state[1].mole_frac_phase_comp_apparent["Liq", "H2O"]
        ) == pytest.approx(0.8, rel=1e-5)
        assert value(
            m.fs.state[1].mole_frac_phase_comp_apparent["Liq", "KHCO3"]
        ) == pytest.approx(1.6e-8, abs=1e-9)
        assert value(
            m.fs.state[1].mole_frac_phase_comp_apparent["Liq", "KHCO3"]
        ) == pytest.approx(1.6e-8, abs=1e-9)
