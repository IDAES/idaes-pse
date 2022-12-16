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
Tests for constructing and using component lists in electrolyte systems
"""
# Import Python libraries
import pytest

# Import Pyomo units
from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    Constraint,
    Expression,
    Set,
    value,
    Var,
    units as pyunits,
)

# Import IDAES cores
from idaes.core import AqueousPhase, VaporPhase
from idaes.core.base.components import *

from idaes.models.properties.modular_properties.state_definitions import FcTP
from idaes.models.properties.modular_properties.base.tests.dummy_eos import DummyEoS
from idaes.models.properties.modular_properties.eos.ideal import Ideal
from idaes.models.properties.modular_properties.reactions.dh_rxn import constant_dh_rxn
from idaes.models.properties.modular_properties.reactions.equilibrium_constant import (
    van_t_hoff,
)
from idaes.models.properties.modular_properties.reactions.equilibrium_forms import (
    power_law_equil,
)
from idaes.models.properties.modular_properties.base.utility import ConcentrationForm

from idaes.core import FlowsheetBlock
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
    StateIndex,
)

from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.solvers import get_solver


# -----------------------------------------------------------------------------
class TestApparentSpeciesBasisNoInherent:
    config = {
        # Specifying components
        "components": {
            "H2O": {
                "type": Solvent,
                "parameter_data": {"mw": (18e-3, pyunits.kg / pyunits.mol)},
            },
            "CO2": {
                "type": Solute,
                "parameter_data": {"mw": (44e-3, pyunits.kg / pyunits.mol)},
            },
            "KHCO3": {
                "type": Apparent,
                "dissociation_species": {"K+": 1, "HCO3-": 1},
                "parameter_data": {"mw": (100.1e-3, pyunits.kg / pyunits.mol)},
            },
            "K+": {
                "type": Cation,
                "charge": +1,
                "parameter_data": {"mw": (39.1e-3, pyunits.kg / pyunits.mol)},
            },
            "HCO3-": {
                "type": Anion,
                "charge": -1,
                "parameter_data": {"mw": (61e-3, pyunits.kg / pyunits.mol)},
            },
            "N2": {
                "type": Component,
                "parameter_data": {"mw": (28e-3, pyunits.kg / pyunits.mol)},
            },
        },
        # Specifying phases
        "phases": {
            "Liq": {
                "type": AqueousPhase,
                "equation_of_state": DummyEoS,
                "equation_of_state_options": {"pH_range": "basic"},
            },
            "Vap": {"type": VaporPhase, "equation_of_state": Ideal},
        },
        # Set base units of measurement
        "base_units": {
            "time": pyunits.s,
            "length": pyunits.m,
            "mass": pyunits.kg,
            "amount": pyunits.mol,
            "temperature": pyunits.K,
        },
        # Specifying state definition
        "state_definition": FcTP,
        "state_bounds": {
            "flow_mol_comp": (0, 100, 1000, pyunits.mol / pyunits.s),
            "temperature": (273.15, 300, 500, pyunits.K),
            "pressure": (5e4, 1e5, 1e6, pyunits.Pa),
        },
        "state_components": StateIndex.apparent,
        "pressure_ref": (101325, pyunits.Pa),
        "temperature_ref": (298.15, pyunits.K),
    }

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

        assert isinstance(m.fs.state[1].flow_mol, Expression)
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

        m.fs.state[1].phase_frac["Liq"].fix(0.5)
        m.fs.state[1].temperature.fix(300)
        m.fs.state[1].pressure.fix(1e5)

        m.fs.state[1].flow_mol_comp["H2O"].fix(0.5)
        m.fs.state[1].flow_mol_comp["CO2"].fix(0.5)
        m.fs.state[1].flow_mol_comp["KHCO3"].fix(0.5)
        m.fs.state[1].flow_mol_comp["N2"].fix(0.5)

        m.fs.state[1].mole_frac_phase_comp["Vap", "KHCO3"].fix(1 / 6)
        m.fs.state[1].mole_frac_phase_comp["Vap", "CO2"].fix(1 / 6)

        assert degrees_of_freedom(m.fs) == 0

        solver = get_solver()
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
    config = {
        # Specifying components
        "components": {
            "H2O": {
                "type": Solvent,
                "dens_mol_liq_comp": dens_mol_H2O,
                "parameter_data": {"mw": (18e-3, pyunits.kg / pyunits.mol)},
            },
            "KHCO3": {
                "type": Apparent,
                "dissociation_species": {"K+": 1, "HCO3-": 1},
                "parameter_data": {"mw": (100.1e-3, pyunits.kg / pyunits.mol)},
            },
            "K2CO3": {
                "type": Apparent,
                "dissociation_species": {"K+": 2, "CO3--": 1},
                "parameter_data": {"mw": (138.2e-3, pyunits.kg / pyunits.mol)},
            },
            "KOH": {
                "type": Apparent,
                "dissociation_species": {"K+": 1, "OH-": 1},
                "parameter_data": {"mw": (56.1e-3, pyunits.kg / pyunits.mol)},
            },
            "H+": {
                "type": Cation,
                "charge": +1,
                "dens_mol_liq_comp": dens_mol_H2O,
                "parameter_data": {"mw": (1e-3, pyunits.kg / pyunits.mol)},
            },
            "K+": {
                "type": Cation,
                "charge": +1,
                "dens_mol_liq_comp": dens_mol_H2O,
                "parameter_data": {"mw": (39.1e-3, pyunits.kg / pyunits.mol)},
            },
            "OH-": {
                "type": Anion,
                "charge": -1,
                "dens_mol_liq_comp": dens_mol_H2O,
                "parameter_data": {"mw": (17e-3, pyunits.kg / pyunits.mol)},
            },
            "HCO3-": {
                "type": Anion,
                "charge": -1,
                "dens_mol_liq_comp": dens_mol_H2O,
                "parameter_data": {"mw": (61e-3, pyunits.kg / pyunits.mol)},
            },
            "CO3--": {
                "type": Anion,
                "charge": -2,
                "dens_mol_liq_comp": dens_mol_H2O,
                "parameter_data": {"mw": (60e-3, pyunits.kg / pyunits.mol)},
            },
        },
        # Specifying phases
        "phases": {
            "Liq": {
                "type": AqueousPhase,
                "equation_of_state": DummyEoS,
                "equation_of_state_options": {"pH_range": "basic"},
            }
        },
        # Set base units of measurement
        "base_units": {
            "time": pyunits.s,
            "length": pyunits.m,
            "mass": pyunits.kg,
            "amount": pyunits.mol,
            "temperature": pyunits.K,
        },
        # Specifying state definition
        "state_definition": FcTP,
        "state_bounds": {
            "flow_mol_comp": (0, 100, 1000, pyunits.mol / pyunits.s),
            "temperature": (273.15, 300, 500, pyunits.K),
            "pressure": (5e4, 1e5, 1e6, pyunits.Pa),
        },
        "state_components": StateIndex.apparent,
        "pressure_ref": (101325, pyunits.Pa),
        "temperature_ref": (298.15, pyunits.K),
        "inherent_reactions": {
            "H2O_si": {
                "stoichiometry": {
                    ("Liq", "H2O"): -1,
                    ("Liq", "H+"): 1,
                    ("Liq", "OH-"): 1,
                },
                "heat_of_reaction": constant_dh_rxn,
                "equilibrium_constant": van_t_hoff,
                "equilibrium_form": power_law_equil,
                "concentration_form": ConcentrationForm.molarity,
                "parameter_data": {
                    "reaction_order": {("Liq", "H+"): 1, ("Liq", "OH-"): 1},
                    "dh_rxn_ref": 1,
                    "k_eq_ref": 1e-14,
                    "T_eq_ref": 350,
                },
            },
            "co3_hco3": {
                "stoichiometry": {
                    ("Liq", "CO3--"): -1,
                    ("Liq", "H2O"): -1,
                    ("Liq", "HCO3-"): 1,
                    ("Liq", "OH-"): 1,
                },
                "heat_of_reaction": constant_dh_rxn,
                "equilibrium_constant": van_t_hoff,
                "equilibrium_form": power_law_equil,
                "concentration_form": ConcentrationForm.molarity,
                "parameter_data": {
                    "reaction_order": {
                        ("Liq", "CO3--"): -1,
                        ("Liq", "HCO3-"): 1,
                        ("Liq", "OH-"): 1,
                    },
                    "dh_rxn_ref": 1,
                    "k_eq_ref": 5e-11,
                    "T_eq_ref": 350,
                },
            },
        },
    }

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

        assert isinstance(m.fs.state[1].flow_mol, Expression)
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
            assert i in ["H2O_si", "co3_hco3"]

        assert isinstance(m.fs.state[1].apparent_inherent_reaction_extent, Var)
        assert len(m.fs.state[1].apparent_inherent_reaction_extent) == 2

    @pytest.mark.component
    def test_solve_for_true_species(self, frame):
        m = frame

        m.fs.state[1].temperature.fix(350)
        m.fs.state[1].pressure.fix(1e5)

        m.fs.state[1].flow_mol_comp["H2O"].fix(0.8 * 2)
        m.fs.state[1].flow_mol_comp["K2CO3"].fix(0.2 * 2)
        m.fs.state[1].flow_mol_comp["KHCO3"].fix(0)
        m.fs.state[1].flow_mol_comp["KOH"].fix(0)

        assert degrees_of_freedom(m.fs) == 0

        m.fs.state.initialize()

        solver = get_solver()
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
    config = {
        # Specifying components
        "components": {
            "H2O": {
                "type": Solvent,
                "parameter_data": {"mw": (18e-3, pyunits.kg / pyunits.mol)},
            },
            "CO2": {
                "type": Solute,
                "parameter_data": {"mw": (44e-3, pyunits.kg / pyunits.mol)},
            },
            "KHCO3": {
                "type": Apparent,
                "dissociation_species": {"K+": 1, "HCO3-": 1},
                "parameter_data": {"mw": (100.1e-3, pyunits.kg / pyunits.mol)},
            },
            "K+": {
                "type": Cation,
                "charge": +1,
                "parameter_data": {"mw": (39.1e-3, pyunits.kg / pyunits.mol)},
            },
            "HCO3-": {
                "type": Anion,
                "charge": -1,
                "parameter_data": {"mw": (61e-3, pyunits.kg / pyunits.mol)},
            },
            "N2": {
                "type": Component,
                "parameter_data": {"mw": (28e-3, pyunits.kg / pyunits.mol)},
            },
        },
        # Specifying phases
        "phases": {
            "Liq": {
                "type": AqueousPhase,
                "equation_of_state": DummyEoS,
                "equation_of_state_options": {"pH_range": "basic"},
            },
            "Vap": {"type": VaporPhase, "equation_of_state": Ideal},
        },
        # Set base units of measurement
        "base_units": {
            "time": pyunits.s,
            "length": pyunits.m,
            "mass": pyunits.kg,
            "amount": pyunits.mol,
            "temperature": pyunits.K,
        },
        # Specifying state definition
        "state_definition": FcTP,
        "state_bounds": {
            "flow_mol_comp": (0, 100, 1000, pyunits.mol / pyunits.s),
            "temperature": (273.15, 300, 500, pyunits.K),
            "pressure": (5e4, 1e5, 1e6, pyunits.Pa),
        },
        "state_components": StateIndex.true,
        "pressure_ref": (101325, pyunits.Pa),
        "temperature_ref": (298.15, pyunits.K),
    }

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

        assert isinstance(m.fs.state[1].flow_mol, Expression)
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
    config = {
        # Specifying components
        "components": {
            "H2O": {
                "type": Solvent,
                "dens_mol_liq_comp": dens_mol_H2O,
                "parameter_data": {"mw": (18e-3, pyunits.kg / pyunits.mol)},
            },
            "KHCO3": {
                "type": Apparent,
                "dissociation_species": {"K+": 1, "HCO3-": 1},
                "parameter_data": {"mw": (100.1e-3, pyunits.kg / pyunits.mol)},
            },
            "K2CO3": {
                "type": Apparent,
                "dissociation_species": {"K+": 2, "CO3--": 1},
                "parameter_data": {"mw": (138.2e-3, pyunits.kg / pyunits.mol)},
            },
            "KOH": {
                "type": Apparent,
                "dissociation_species": {"K+": 1, "OH-": 1},
                "parameter_data": {"mw": (56.1e-3, pyunits.kg / pyunits.mol)},
            },
            "H+": {
                "type": Cation,
                "charge": +1,
                "dens_mol_liq_comp": dens_mol_H2O,
                "parameter_data": {"mw": (1e-3, pyunits.kg / pyunits.mol)},
            },
            "K+": {
                "type": Cation,
                "charge": +1,
                "dens_mol_liq_comp": dens_mol_H2O,
                "parameter_data": {"mw": (39.1e-3, pyunits.kg / pyunits.mol)},
            },
            "OH-": {
                "type": Anion,
                "charge": -1,
                "dens_mol_liq_comp": dens_mol_H2O,
                "parameter_data": {"mw": (17e-3, pyunits.kg / pyunits.mol)},
            },
            "HCO3-": {
                "type": Anion,
                "charge": -1,
                "dens_mol_liq_comp": dens_mol_H2O,
                "parameter_data": {"mw": (61e-3, pyunits.kg / pyunits.mol)},
            },
            "CO3--": {
                "type": Anion,
                "charge": -2,
                "dens_mol_liq_comp": dens_mol_H2O,
                "parameter_data": {"mw": (60e-3, pyunits.kg / pyunits.mol)},
            },
        },
        # Specifying phases
        "phases": {
            "Liq": {
                "type": AqueousPhase,
                "equation_of_state": DummyEoS,
                "equation_of_state_options": {"pH_range": "basic"},
            }
        },
        # Set base units of measurement
        "base_units": {
            "time": pyunits.s,
            "length": pyunits.m,
            "mass": pyunits.kg,
            "amount": pyunits.mol,
            "temperature": pyunits.K,
        },
        # Specifying state definition
        "state_definition": FcTP,
        "state_bounds": {
            "flow_mol_comp": (0, 100, 1000, pyunits.mol / pyunits.s),
            "temperature": (273.15, 300, 500, pyunits.K),
            "pressure": (5e4, 1e5, 1e6, pyunits.Pa),
        },
        "state_components": StateIndex.true,
        "pressure_ref": (101325, pyunits.Pa),
        "temperature_ref": (298.15, pyunits.K),
        "inherent_reactions": {
            "H2O_si": {
                "stoichiometry": {
                    ("Liq", "H2O"): -1,
                    ("Liq", "H+"): 1,
                    ("Liq", "OH-"): 1,
                },
                "heat_of_reaction": constant_dh_rxn,
                "equilibrium_constant": van_t_hoff,
                "equilibrium_form": power_law_equil,
                "concentration_form": ConcentrationForm.molarity,
                "parameter_data": {
                    "reaction_order": {("Liq", "H+"): 1, ("Liq", "OH-"): 1},
                    "dh_rxn_ref": 1,
                    "k_eq_ref": 1e-14,
                    "T_eq_ref": 350,
                },
            },
            "co3_hco3": {
                "stoichiometry": {
                    ("Liq", "CO3--"): -1,
                    ("Liq", "H2O"): -1,
                    ("Liq", "HCO3-"): 1,
                    ("Liq", "OH-"): 1,
                },
                "heat_of_reaction": constant_dh_rxn,
                "equilibrium_constant": van_t_hoff,
                "equilibrium_form": power_law_equil,
                "concentration_form": ConcentrationForm.molarity,
                "parameter_data": {
                    "reaction_order": {
                        ("Liq", "CO3--"): -1,
                        ("Liq", "HCO3-"): 1,
                        ("Liq", "OH-"): 1,
                    },
                    "dh_rxn_ref": 1,
                    "k_eq_ref": 5e-11,
                    "T_eq_ref": 350,
                },
            },
        },
        "default_scaling_factors": {
            ("mole_frac_comp", "OH-"): 1e8,
            ("mole_frac_comp", "HCO3-"): 1e8,
            ("mole_frac_comp", "H+"): 1e16,
        },
    }

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

        assert isinstance(m.fs.state[1].flow_mol, Expression)
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
            assert i in ["H2O_si", "co3_hco3"]

        assert not hasattr(m.fs.state[1], "apparent_inherent_reaction_extent")

    @pytest.mark.component
    def test_solve_for_true_species(self, frame):
        m = frame

        m.fs.state[1].temperature.fix(350)
        m.fs.state[1].pressure.fix(1e5)

        m.fs.state[1].flow_mol_comp["H2O"].fix(0.5716 * 2)
        m.fs.state[1].flow_mol_comp["K+"].fix(0.2858 * 2)
        m.fs.state[1].flow_mol_comp["H+"].fix(2.9e-16 * 2)
        m.fs.state[1].flow_mol_comp["HCO3-"].fix(1.14e-8 * 2)
        m.fs.state[1].flow_mol_comp["CO3--"].fix(0.1429 * 2)
        m.fs.state[1].flow_mol_comp["OH-"].fix(1.14e-8 * 2)

        assert degrees_of_freedom(m.fs) == 0

        m.fs.state.initialize()

        solver = get_solver()
        res = solver.solve(m.fs)

        # Check for optimal solution
        assert check_optimal_termination(res)

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
