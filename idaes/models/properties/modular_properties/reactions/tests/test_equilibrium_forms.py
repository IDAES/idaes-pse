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
Tests for equilibrium forms
"""

import pytest

from pyomo.environ import Block, ConcreteModel, Param, Var, units as pyunits, value

from idaes.models.properties.modular_properties.base.generic_reaction import (
    GenericReactionParameterBlock,
    ConcentrationForm,
)
from idaes.models.properties.modular_properties.reactions.equilibrium_forms import *

from idaes.core import SolidPhase
from idaes.core.util.testing import PhysicalParameterTestBlock
from idaes.core.util.misc import add_object_reference
from idaes.core.util.math import smooth_max
from idaes.core.util.exceptions import ConfigurationError


@pytest.mark.unit
def test_power_law_equil_no_order():
    m = ConcreteModel()

    # # Add a test thermo package for validation
    m.pparams = PhysicalParameterTestBlock()

    # Add a solid phase for testing
    m.pparams.sol = SolidPhase()

    m.thermo = m.pparams.build_state_block([1])

    # Create a dummy reaction parameter block
    m.rparams = GenericReactionParameterBlock(
        property_package=m.pparams,
        base_units={
            "time": pyunits.s,
            "mass": pyunits.kg,
            "amount": pyunits.mol,
            "length": pyunits.m,
            "temperature": pyunits.K,
        },
        equilibrium_reactions={
            "r1": {
                "stoichiometry": {
                    ("p1", "c1"): -1,
                    ("p1", "c2"): 2,
                    ("sol", "c1"): -3,
                    ("sol", "c2"): 4,
                },
                "equilibrium_form": power_law_equil,
                "concentration_form": ConcentrationForm.moleFraction,
            }
        },
    )

    # Create a dummy state block
    m.rxn = Block([1])
    add_object_reference(
        m.rxn[1], "phase_component_set", m.pparams._phase_component_set
    )
    add_object_reference(m.rxn[1], "params", m.rparams)
    add_object_reference(m.rxn[1], "state_ref", m.thermo[1])

    m.rxn[1].k_eq = Var(["r1"], initialize=1)

    power_law_equil.build_parameters(
        m.rparams.reaction_r1, m.rparams.config.equilibrium_reactions["r1"]
    )

    # Check parameter construction
    assert isinstance(m.rparams.reaction_r1.reaction_order, Var)
    assert len(m.rparams.reaction_r1.reaction_order) == 6

    assert m.rparams.reaction_r1.reaction_order["p1", "c1"].value == -1
    assert m.rparams.reaction_r1.reaction_order["p1", "c2"].value == 2
    assert m.rparams.reaction_r1.reaction_order["p2", "c1"].value == 0
    assert m.rparams.reaction_r1.reaction_order["p2", "c2"].value == 0
    # Solids should have zero order, as they are excluded
    assert m.rparams.reaction_r1.reaction_order["sol", "c1"].value == 0
    assert m.rparams.reaction_r1.reaction_order["sol", "c2"].value == 0

    # Check reaction form
    rform = power_law_equil.return_expression(
        m.rxn[1], m.rparams.reaction_r1, "r1", 300
    )

    assert str(rform) == str(
        m.rxn[1].k_eq["r1"]
        == (
            m.thermo[1].mole_frac_phase_comp["p1", "c1"]
            ** m.rparams.reaction_r1.reaction_order["p1", "c1"]
            * m.thermo[1].mole_frac_phase_comp["p1", "c2"]
            ** m.rparams.reaction_r1.reaction_order["p1", "c2"]
        )
    )


@pytest.mark.unit
def test_power_law_equil_with_order():
    m = ConcreteModel()

    # # Add a test thermo package for validation
    m.pparams = PhysicalParameterTestBlock()

    # Add a solid phase for testing
    m.pparams.sol = SolidPhase()

    m.thermo = m.pparams.build_state_block([1])

    # Create a dummy reaction parameter block
    m.rparams = GenericReactionParameterBlock(
        property_package=m.pparams,
        base_units={
            "time": pyunits.s,
            "mass": pyunits.kg,
            "amount": pyunits.mol,
            "length": pyunits.m,
            "temperature": pyunits.K,
        },
        equilibrium_reactions={
            "r1": {
                "stoichiometry": {
                    ("p1", "c1"): -1,
                    ("p1", "c2"): 2,
                    ("sol", "c1"): -3,
                    ("sol", "c2"): 4,
                },
                "equilibrium_form": power_law_equil,
                "concentration_form": ConcentrationForm.moleFraction,
                "parameter_data": {
                    "reaction_order": {
                        ("p1", "c1"): 1,
                        ("p1", "c2"): 2,
                        ("p2", "c1"): 3,
                        ("p2", "c2"): 4,
                        ("sol", "c1"): 5,
                        ("sol", "c2"): 6,
                    }
                },
            }
        },
    )

    # Create a dummy state block
    m.rxn = Block([1])
    add_object_reference(
        m.rxn[1], "phase_component_set", m.pparams._phase_component_set
    )
    add_object_reference(m.rxn[1], "params", m.rparams)
    add_object_reference(m.rxn[1], "state_ref", m.thermo[1])

    m.rxn[1].k_eq = Var(["r1"], initialize=1)

    power_law_equil.build_parameters(
        m.rparams.reaction_r1, m.rparams.config.equilibrium_reactions["r1"]
    )

    # Check parameter construction
    assert isinstance(m.rparams.reaction_r1.reaction_order, Var)
    assert len(m.rparams.reaction_r1.reaction_order) == 6
    assert m.rparams.reaction_r1.reaction_order["p1", "c1"].value == 1
    assert m.rparams.reaction_r1.reaction_order["p1", "c2"].value == 2
    assert m.rparams.reaction_r1.reaction_order["p2", "c1"].value == 3
    assert m.rparams.reaction_r1.reaction_order["p2", "c2"].value == 4
    assert m.rparams.reaction_r1.reaction_order["sol", "c1"].value == 5
    assert m.rparams.reaction_r1.reaction_order["sol", "c2"].value == 6

    # Check reaction form
    rform = power_law_equil.return_expression(
        m.rxn[1], m.rparams.reaction_r1, "r1", 300
    )

    assert str(rform) == str(
        m.rxn[1].k_eq["r1"]
        == (
            m.thermo[1].mole_frac_phase_comp["p1", "c1"]
            ** m.rparams.reaction_r1.reaction_order["p1", "c1"]
            * m.thermo[1].mole_frac_phase_comp["p1", "c2"]
            ** m.rparams.reaction_r1.reaction_order["p1", "c2"]
            * m.thermo[1].mole_frac_phase_comp["p2", "c1"]
            ** m.rparams.reaction_r1.reaction_order["p2", "c1"]
            * m.thermo[1].mole_frac_phase_comp["p2", "c2"]
            ** m.rparams.reaction_r1.reaction_order["p2", "c2"]
            * m.thermo[1].mole_frac_phase_comp["sol", "c1"]
            ** m.rparams.reaction_r1.reaction_order["sol", "c1"]
            * m.thermo[1].mole_frac_phase_comp["sol", "c2"]
            ** m.rparams.reaction_r1.reaction_order["sol", "c2"]
        )
    )


@pytest.mark.unit
def test_log_power_law_equil_no_order():
    m = ConcreteModel()

    # # Add a test thermo package for validation
    m.pparams = PhysicalParameterTestBlock()

    # Add a solid phase for testing
    m.pparams.sol = SolidPhase()

    m.thermo = m.pparams.build_state_block([1])
    m.thermo[1].log_mole_frac_phase_comp = Var(m.pparams._phase_component_set)

    # Create a dummy reaction parameter block
    m.rparams = GenericReactionParameterBlock(
        property_package=m.pparams,
        base_units={
            "time": pyunits.s,
            "mass": pyunits.kg,
            "amount": pyunits.mol,
            "length": pyunits.m,
            "temperature": pyunits.K,
        },
        equilibrium_reactions={
            "r1": {
                "stoichiometry": {
                    ("p1", "c1"): -1,
                    ("p1", "c2"): 2,
                    ("sol", "c1"): -3,
                    ("sol", "c2"): 4,
                },
                "equilibrium_form": log_power_law_equil,
                "concentration_form": ConcentrationForm.moleFraction,
            }
        },
    )

    # Create a dummy state block
    m.rxn = Block([1])
    add_object_reference(
        m.rxn[1], "phase_component_set", m.pparams._phase_component_set
    )
    add_object_reference(m.rxn[1], "params", m.rparams)
    add_object_reference(m.rxn[1], "state_ref", m.thermo[1])

    m.rxn[1].log_k_eq = Var(["r1"], initialize=1)

    # Check parameter construction
    assert isinstance(m.rparams.reaction_r1.reaction_order, Var)
    assert len(m.rparams.reaction_r1.reaction_order) == 6

    assert m.rparams.reaction_r1.reaction_order["p1", "c1"].value == -1
    assert m.rparams.reaction_r1.reaction_order["p1", "c2"].value == 2
    assert m.rparams.reaction_r1.reaction_order["p2", "c1"].value == 0
    assert m.rparams.reaction_r1.reaction_order["p2", "c2"].value == 0
    # Solids should have zero order, as they are excluded
    assert m.rparams.reaction_r1.reaction_order["sol", "c1"].value == 0
    assert m.rparams.reaction_r1.reaction_order["sol", "c2"].value == 0

    # Check reaction form
    rform = log_power_law_equil.return_expression(
        m.rxn[1], m.rparams.reaction_r1, "r1", 300
    )

    assert str(rform) == str(
        m.rxn[1].log_k_eq["r1"]
        == m.rparams.reaction_r1.reaction_order["p1", "c1"]
        * m.thermo[1].log_mole_frac_phase_comp["p1", "c1"]
        + m.rparams.reaction_r1.reaction_order["p1", "c2"]
        * m.thermo[1].log_mole_frac_phase_comp["p1", "c2"]
    )


@pytest.mark.unit
def test_log_power_law_equil_with_order():
    m = ConcreteModel()

    # # Add a test thermo package for validation
    m.pparams = PhysicalParameterTestBlock()

    # Add a solid phase for testing
    m.pparams.sol = SolidPhase()

    m.thermo = m.pparams.build_state_block([1])
    m.thermo[1].log_mole_frac_phase_comp = Var(m.pparams._phase_component_set)

    # Create a dummy reaction parameter block
    m.rparams = GenericReactionParameterBlock(
        property_package=m.pparams,
        base_units={
            "time": pyunits.s,
            "mass": pyunits.kg,
            "amount": pyunits.mol,
            "length": pyunits.m,
            "temperature": pyunits.K,
        },
        equilibrium_reactions={
            "r1": {
                "stoichiometry": {
                    ("p1", "c1"): -1,
                    ("p1", "c2"): 2,
                    ("sol", "c1"): -3,
                    ("sol", "c2"): 4,
                },
                "equilibrium_form": log_power_law_equil,
                "concentration_form": ConcentrationForm.moleFraction,
                "parameter_data": {
                    "reaction_order": {
                        ("p1", "c1"): 1,
                        ("p1", "c2"): 2,
                        ("p2", "c1"): 3,
                        ("p2", "c2"): 4,
                        ("sol", "c1"): 5,
                        ("sol", "c2"): 6,
                    }
                },
            }
        },
    )

    # Create a dummy state block
    m.rxn = Block([1])
    add_object_reference(
        m.rxn[1], "phase_component_set", m.pparams._phase_component_set
    )
    add_object_reference(m.rxn[1], "params", m.rparams)
    add_object_reference(m.rxn[1], "state_ref", m.thermo[1])

    m.rxn[1].log_k_eq = Var(["r1"], initialize=1)

    # Check parameter construction
    assert isinstance(m.rparams.reaction_r1.reaction_order, Var)
    assert len(m.rparams.reaction_r1.reaction_order) == 6
    assert m.rparams.reaction_r1.reaction_order["p1", "c1"].value == 1
    assert m.rparams.reaction_r1.reaction_order["p1", "c2"].value == 2
    assert m.rparams.reaction_r1.reaction_order["p2", "c1"].value == 3
    assert m.rparams.reaction_r1.reaction_order["p2", "c2"].value == 4
    assert m.rparams.reaction_r1.reaction_order["sol", "c1"].value == 5
    assert m.rparams.reaction_r1.reaction_order["sol", "c2"].value == 6

    # Check reaction form
    rform = log_power_law_equil.return_expression(
        m.rxn[1], m.rparams.reaction_r1, "r1", 300
    )

    assert str(rform) == str(
        m.rxn[1].log_k_eq["r1"]
        == m.rparams.reaction_r1.reaction_order["p1", "c1"]
        * m.thermo[1].log_mole_frac_phase_comp["p1", "c1"]
        + m.rparams.reaction_r1.reaction_order["p1", "c2"]
        * m.thermo[1].log_mole_frac_phase_comp["p1", "c2"]
        + m.rparams.reaction_r1.reaction_order["p2", "c1"]
        * m.thermo[1].log_mole_frac_phase_comp["p2", "c1"]
        + m.rparams.reaction_r1.reaction_order["p2", "c2"]
        * m.thermo[1].log_mole_frac_phase_comp["p2", "c2"]
        + m.rparams.reaction_r1.reaction_order["sol", "c1"]
        * m.thermo[1].log_mole_frac_phase_comp["sol", "c1"]
        + m.rparams.reaction_r1.reaction_order["sol", "c2"]
        * m.thermo[1].log_mole_frac_phase_comp["sol", "c2"]
    )


@pytest.mark.unit
def test_solubility_no_order():
    m = ConcreteModel()

    # # Add a test thermo package for validation
    m.pparams = PhysicalParameterTestBlock()

    # Add a solid phase for testing
    m.pparams.sol = SolidPhase()

    m.thermo = m.pparams.build_state_block([1])

    # Create a dummy reaction parameter block
    m.rparams = GenericReactionParameterBlock(
        property_package=m.pparams,
        base_units={
            "time": pyunits.s,
            "mass": pyunits.kg,
            "amount": pyunits.mol,
            "length": pyunits.m,
            "temperature": pyunits.K,
        },
        equilibrium_reactions={
            "r1": {
                "stoichiometry": {
                    ("p1", "c1"): -1,
                    ("p1", "c2"): 2,
                    ("sol", "c1"): -3,
                    ("sol", "c2"): 4,
                },
                "equilibrium_form": solubility_product,
                "concentration_form": ConcentrationForm.moleFraction,
            }
        },
    )

    # Create a dummy state block
    m.rxn = Block([1])
    add_object_reference(
        m.rxn[1], "phase_component_set", m.pparams._phase_component_set
    )
    add_object_reference(m.rxn[1], "params", m.rparams)
    add_object_reference(m.rxn[1], "state_ref", m.thermo[1])

    m.rxn[1].k_eq = Var(["r1"], initialize=1)

    # Check parameter construction
    assert isinstance(m.rparams.reaction_r1.eps, Param)
    assert value(m.rparams.reaction_r1.eps) == 1e-4
    assert isinstance(m.rparams.reaction_r1.reaction_order, Var)
    assert len(m.rparams.reaction_r1.reaction_order) == 6

    assert m.rparams.reaction_r1.reaction_order["p1", "c1"].value == -1
    assert m.rparams.reaction_r1.reaction_order["p1", "c2"].value == 2
    assert m.rparams.reaction_r1.reaction_order["p2", "c1"].value == 0
    assert m.rparams.reaction_r1.reaction_order["p2", "c2"].value == 0
    # Solids should have zero order, as they are excluded
    assert m.rparams.reaction_r1.reaction_order["sol", "c1"].value == 0
    assert m.rparams.reaction_r1.reaction_order["sol", "c2"].value == 0

    # Check reaction form
    rform = solubility_product.return_expression(
        m.rxn[1], m.rparams.reaction_r1, "r1", 300
    )

    s = (
        m.thermo[1].flow_mol_phase_comp["sol", "c1"]
        + m.thermo[1].flow_mol_phase_comp["sol", "c2"]
    ) / (pyunits.mol / pyunits.s)
    Q = m.rxn[1].k_eq["r1"] - (
        m.thermo[1].mole_frac_phase_comp["p1", "c1"]
        ** m.rparams.reaction_r1.reaction_order["p1", "c1"]
        * m.thermo[1].mole_frac_phase_comp["p1", "c2"]
        ** m.rparams.reaction_r1.reaction_order["p1", "c2"]
    )

    assert str(rform) == str(s - smooth_max(0, s - Q, m.rparams.reaction_r1.eps) == 0)


@pytest.mark.unit
def test_solubility_product_with_order():
    m = ConcreteModel()

    # # Add a test thermo package for validation
    m.pparams = PhysicalParameterTestBlock()

    # Add a solid phase for testing
    m.pparams.sol = SolidPhase()

    m.thermo = m.pparams.build_state_block([1])

    # Create a dummy reaction parameter block
    m.rparams = GenericReactionParameterBlock(
        property_package=m.pparams,
        base_units={
            "time": pyunits.s,
            "mass": pyunits.kg,
            "amount": pyunits.mol,
            "length": pyunits.m,
            "temperature": pyunits.K,
        },
        equilibrium_reactions={
            "r1": {
                "stoichiometry": {
                    ("p1", "c1"): -1,
                    ("p1", "c2"): 2,
                    ("sol", "c1"): -3,
                    ("sol", "c2"): 4,
                },
                "equilibrium_form": solubility_product,
                "concentration_form": ConcentrationForm.moleFraction,
                "parameter_data": {
                    "reaction_order": {
                        ("p1", "c1"): 1,
                        ("p1", "c2"): 2,
                        ("p2", "c1"): 3,
                        ("p2", "c2"): 4,
                        ("sol", "c1"): 5,
                        ("sol", "c2"): 6,
                    }
                },
            }
        },
    )

    # Create a dummy state block
    m.rxn = Block([1])
    add_object_reference(
        m.rxn[1], "phase_component_set", m.pparams._phase_component_set
    )
    add_object_reference(m.rxn[1], "params", m.rparams)
    add_object_reference(m.rxn[1], "state_ref", m.thermo[1])

    m.rxn[1].k_eq = Var(["r1"], initialize=1)

    # Check parameter construction
    assert isinstance(m.rparams.reaction_r1.eps, Param)
    assert value(m.rparams.reaction_r1.eps) == 1e-4
    assert isinstance(m.rparams.reaction_r1.reaction_order, Var)
    assert len(m.rparams.reaction_r1.reaction_order) == 6
    assert m.rparams.reaction_r1.reaction_order["p1", "c1"].value == 1
    assert m.rparams.reaction_r1.reaction_order["p1", "c2"].value == 2
    assert m.rparams.reaction_r1.reaction_order["p2", "c1"].value == 3
    assert m.rparams.reaction_r1.reaction_order["p2", "c2"].value == 4
    assert m.rparams.reaction_r1.reaction_order["sol", "c1"].value == 5
    assert m.rparams.reaction_r1.reaction_order["sol", "c2"].value == 6

    # Check reaction form
    rform = solubility_product.return_expression(
        m.rxn[1], m.rparams.reaction_r1, "r1", 300
    )

    s = (
        m.thermo[1].flow_mol_phase_comp["sol", "c1"]
        + m.thermo[1].flow_mol_phase_comp["sol", "c2"]
    ) / (pyunits.mol / pyunits.s)
    Q = m.rxn[1].k_eq["r1"] - (
        m.thermo[1].mole_frac_phase_comp["p1", "c1"]
        ** m.rparams.reaction_r1.reaction_order["p1", "c1"]
        * m.thermo[1].mole_frac_phase_comp["p1", "c2"]
        ** m.rparams.reaction_r1.reaction_order["p1", "c2"]
        * m.thermo[1].mole_frac_phase_comp["p2", "c1"]
        ** m.rparams.reaction_r1.reaction_order["p2", "c1"]
        * m.thermo[1].mole_frac_phase_comp["p2", "c2"]
        ** m.rparams.reaction_r1.reaction_order["p2", "c2"]
        * m.thermo[1].mole_frac_phase_comp["sol", "c1"]
        ** m.rparams.reaction_r1.reaction_order["sol", "c1"]
        * m.thermo[1].mole_frac_phase_comp["sol", "c2"]
        ** m.rparams.reaction_r1.reaction_order["sol", "c2"]
    )

    assert str(rform) == str(s - smooth_max(0, s - Q, m.rparams.reaction_r1.eps) == 0)


@pytest.mark.unit
def test_solubility_no_solids():
    m = ConcreteModel()

    # # Add a test thermo package for validation
    m.pparams = PhysicalParameterTestBlock()

    m.thermo = m.pparams.build_state_block([1])

    # Create a dummy reaction parameter block
    m.rparams = GenericReactionParameterBlock(
        property_package=m.pparams,
        base_units={
            "time": pyunits.s,
            "mass": pyunits.kg,
            "amount": pyunits.mol,
            "length": pyunits.m,
            "temperature": pyunits.K,
        },
        equilibrium_reactions={
            "r1": {
                "stoichiometry": {("p1", "c1"): -1, ("p1", "c2"): 2},
                "equilibrium_form": solubility_product,
                "concentration_form": ConcentrationForm.moleFraction,
            }
        },
    )

    # Create a dummy state block
    m.rxn = Block([1])
    add_object_reference(
        m.rxn[1], "phase_component_set", m.pparams._phase_component_set
    )
    add_object_reference(m.rxn[1], "params", m.rparams)
    add_object_reference(m.rxn[1], "state_ref", m.thermo[1])

    m.rxn[1].k_eq = Var(["r1"], initialize=1)

    solubility_product.build_parameters(
        m.rparams.reaction_r1, m.rparams.config.equilibrium_reactions["r1"]
    )

    # Check reaction form - should raise exception
    with pytest.raises(
        ConfigurationError,
        match="did not find a solid phase component for "
        "precipitation reaction r1. This is likely due to the "
        "reaction configuration.",
    ):
        solubility_product.return_expression(m.rxn[1], m.rparams.reaction_r1, "r1", 300)


@pytest.mark.unit
def test_log_solubility_no_order():
    m = ConcreteModel()

    # # Add a test thermo package for validation
    m.pparams = PhysicalParameterTestBlock()

    # Add a solid phase for testing
    m.pparams.sol = SolidPhase()

    m.thermo = m.pparams.build_state_block([1])

    # Create a dummy reaction parameter block
    m.rparams = GenericReactionParameterBlock(
        property_package=m.pparams,
        base_units={
            "time": pyunits.s,
            "mass": pyunits.kg,
            "amount": pyunits.mol,
            "length": pyunits.m,
            "temperature": pyunits.K,
        },
        equilibrium_reactions={
            "r1": {
                "stoichiometry": {
                    ("p1", "c1"): -1,
                    ("p1", "c2"): 2,
                    ("sol", "c1"): -3,
                    ("sol", "c2"): 4,
                },
                "equilibrium_form": log_solubility_product,
                "concentration_form": ConcentrationForm.moleFraction,
            }
        },
    )

    # Create a dummy state block
    m.rxn = Block([1])
    add_object_reference(
        m.rxn[1], "phase_component_set", m.pparams._phase_component_set
    )
    add_object_reference(m.rxn[1], "params", m.rparams)
    add_object_reference(m.rxn[1], "state_ref", m.thermo[1])

    m.rxn[1].log_k_eq = Var(["r1"], initialize=1)
    m.thermo[1].log_mole_frac_phase_comp = Var(
        m.pparams._phase_component_set, initialize=1
    )

    # Check parameter construction
    assert isinstance(m.rparams.reaction_r1.eps, Param)
    assert value(m.rparams.reaction_r1.eps) == 1e-4
    assert isinstance(m.rparams.reaction_r1.reaction_order, Var)
    assert len(m.rparams.reaction_r1.reaction_order) == 6

    assert m.rparams.reaction_r1.reaction_order["p1", "c1"].value == -1
    assert m.rparams.reaction_r1.reaction_order["p1", "c2"].value == 2
    assert m.rparams.reaction_r1.reaction_order["p2", "c1"].value == 0
    assert m.rparams.reaction_r1.reaction_order["p2", "c2"].value == 0
    # Solids should have zero order, as they are excluded
    assert m.rparams.reaction_r1.reaction_order["sol", "c1"].value == 0
    assert m.rparams.reaction_r1.reaction_order["sol", "c2"].value == 0

    # Check reaction form
    rform = log_solubility_product.return_expression(
        m.rxn[1], m.rparams.reaction_r1, "r1", 300
    )

    s = (
        m.thermo[1].flow_mol_phase_comp["sol", "c1"]
        + m.thermo[1].flow_mol_phase_comp["sol", "c2"]
    ) / (pyunits.mol / pyunits.s)
    Q = m.rxn[1].log_k_eq["r1"] - (
        m.thermo[1].log_mole_frac_phase_comp["p1", "c1"]
        * m.rparams.reaction_r1.reaction_order["p1", "c1"]
        + m.thermo[1].log_mole_frac_phase_comp["p1", "c2"]
        * m.rparams.reaction_r1.reaction_order["p1", "c2"]
    )

    assert str(rform) == str(s - smooth_max(0, s - Q, m.rparams.reaction_r1.eps) == 0)


@pytest.mark.unit
def test_log_solubility_product_with_order():
    m = ConcreteModel()

    # # Add a test thermo package for validation
    m.pparams = PhysicalParameterTestBlock()

    # Add a solid phase for testing
    m.pparams.sol = SolidPhase()

    m.thermo = m.pparams.build_state_block([1])

    # Create a dummy reaction parameter block
    m.rparams = GenericReactionParameterBlock(
        property_package=m.pparams,
        base_units={
            "time": pyunits.s,
            "mass": pyunits.kg,
            "amount": pyunits.mol,
            "length": pyunits.m,
            "temperature": pyunits.K,
        },
        equilibrium_reactions={
            "r1": {
                "stoichiometry": {
                    ("p1", "c1"): -1,
                    ("p1", "c2"): 2,
                    ("sol", "c1"): -3,
                    ("sol", "c2"): 4,
                },
                "equilibrium_form": log_solubility_product,
                "concentration_form": ConcentrationForm.moleFraction,
                "parameter_data": {
                    "reaction_order": {
                        ("p1", "c1"): 1,
                        ("p1", "c2"): 2,
                        ("p2", "c1"): 3,
                        ("p2", "c2"): 4,
                        ("sol", "c1"): 5,
                        ("sol", "c2"): 6,
                    }
                },
            }
        },
    )

    # Create a dummy state block
    m.rxn = Block([1])
    add_object_reference(
        m.rxn[1], "phase_component_set", m.pparams._phase_component_set
    )
    add_object_reference(m.rxn[1], "params", m.rparams)
    add_object_reference(m.rxn[1], "state_ref", m.thermo[1])

    m.rxn[1].log_k_eq = Var(["r1"], initialize=1)
    m.thermo[1].log_mole_frac_phase_comp = Var(
        m.pparams._phase_component_set, initialize=1
    )

    # Check parameter construction
    assert isinstance(m.rparams.reaction_r1.eps, Param)
    assert value(m.rparams.reaction_r1.eps) == 1e-4
    assert isinstance(m.rparams.reaction_r1.reaction_order, Var)
    assert len(m.rparams.reaction_r1.reaction_order) == 6
    assert m.rparams.reaction_r1.reaction_order["p1", "c1"].value == 1
    assert m.rparams.reaction_r1.reaction_order["p1", "c2"].value == 2
    assert m.rparams.reaction_r1.reaction_order["p2", "c1"].value == 3
    assert m.rparams.reaction_r1.reaction_order["p2", "c2"].value == 4
    assert m.rparams.reaction_r1.reaction_order["sol", "c1"].value == 5
    assert m.rparams.reaction_r1.reaction_order["sol", "c2"].value == 6

    # Check reaction form
    rform = log_solubility_product.return_expression(
        m.rxn[1], m.rparams.reaction_r1, "r1", 300
    )

    s = (
        m.thermo[1].flow_mol_phase_comp["sol", "c1"]
        + m.thermo[1].flow_mol_phase_comp["sol", "c2"]
    ) / (pyunits.mol / pyunits.s)
    Q = m.rxn[1].log_k_eq["r1"] - (
        m.thermo[1].log_mole_frac_phase_comp["p1", "c1"]
        * m.rparams.reaction_r1.reaction_order["p1", "c1"]
        + m.thermo[1].log_mole_frac_phase_comp["p1", "c2"]
        * m.rparams.reaction_r1.reaction_order["p1", "c2"]
        + m.thermo[1].log_mole_frac_phase_comp["p2", "c1"]
        * m.rparams.reaction_r1.reaction_order["p2", "c1"]
        + m.thermo[1].log_mole_frac_phase_comp["p2", "c2"]
        * m.rparams.reaction_r1.reaction_order["p2", "c2"]
        + m.thermo[1].log_mole_frac_phase_comp["sol", "c1"]
        * m.rparams.reaction_r1.reaction_order["sol", "c1"]
        + m.thermo[1].log_mole_frac_phase_comp["sol", "c2"]
        * m.rparams.reaction_r1.reaction_order["sol", "c2"]
    )

    assert str(rform) == str(s - smooth_max(0, s - Q, m.rparams.reaction_r1.eps) == 0)


@pytest.mark.unit
def test_log_solubility_no_solids():
    m = ConcreteModel()

    # # Add a test thermo package for validation
    m.pparams = PhysicalParameterTestBlock()

    m.thermo = m.pparams.build_state_block([1])

    # Create a dummy reaction parameter block
    m.rparams = GenericReactionParameterBlock(
        property_package=m.pparams,
        base_units={
            "time": pyunits.s,
            "mass": pyunits.kg,
            "amount": pyunits.mol,
            "length": pyunits.m,
            "temperature": pyunits.K,
        },
        equilibrium_reactions={
            "r1": {
                "stoichiometry": {("p1", "c1"): -1, ("p1", "c2"): 2},
                "equilibrium_form": log_solubility_product,
                "concentration_form": ConcentrationForm.moleFraction,
            }
        },
    )

    # Create a dummy state block
    m.rxn = Block([1])
    add_object_reference(
        m.rxn[1], "phase_component_set", m.pparams._phase_component_set
    )
    add_object_reference(m.rxn[1], "params", m.rparams)
    add_object_reference(m.rxn[1], "state_ref", m.thermo[1])

    m.rxn[1].log_k_eq = Var(["r1"], initialize=1)
    m.thermo[1].log_mole_frac_phase_comp = Var(
        m.pparams._phase_component_set, initialize=1
    )

    log_solubility_product.build_parameters(
        m.rparams.reaction_r1, m.rparams.config.equilibrium_reactions["r1"]
    )

    # Check reaction form - should raise exception
    with pytest.raises(
        ConfigurationError,
        match="did not find a solid phase component for "
        "precipitation reaction r1. This is likely due to the "
        "reaction configuration.",
    ):
        log_solubility_product.return_expression(
            m.rxn[1], m.rparams.reaction_r1, "r1", 300
        )
