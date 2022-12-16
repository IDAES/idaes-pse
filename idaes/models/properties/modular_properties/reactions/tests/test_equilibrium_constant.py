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
Tests for rate forms
"""

import pytest

from pyomo.environ import ConcreteModel, log, Var, units as pyunits, value
from pyomo.util.check_units import assert_units_equivalent

from idaes.models.properties.modular_properties.base.generic_reaction import (
    GenericReactionParameterBlock,
    ConcentrationForm,
)
from idaes.models.properties.modular_properties.reactions.equilibrium_constant import *
from idaes.models.properties.modular_properties.reactions.dh_rxn import constant_dh_rxn
from idaes.core import MaterialFlowBasis
from idaes.core.util.testing import PhysicalParameterTestBlock
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.constants import Constants as c


@pytest.fixture
def model():
    m = ConcreteModel()

    # # Add a test thermo package for validation
    m.pparams = PhysicalParameterTestBlock()
    m.thermo = m.pparams.build_state_block([1])

    m.rparams = GenericReactionParameterBlock(
        property_package=m.pparams,
        reaction_basis=MaterialFlowBasis.molar,
        equilibrium_reactions={
            "r1": {
                "stoichiometry": {("p1", "c1"): -1, ("p1", "c2"): 2},
                "equilibrium_form": "foo",
                "concentration_form": ConcentrationForm.moleFraction,
            }
        },
        base_units={
            "amount": pyunits.mol,
            "mass": pyunits.kg,
            "time": pyunits.s,
            "length": pyunits.m,
            "temperature": pyunits.K,
        },
    )

    m.rxn = m.rparams.build_reaction_block(
        [1], state_block=m.thermo, has_equilibrium=False
    )

    m.rxn[1].dh_rxn = Var(["r1"], initialize=1, units=pyunits.J / pyunits.mol)

    return m


class TestConstantKeq(object):
    @pytest.mark.unit
    def test_ConstantKeq_mole_frac(self, model):
        model.rparams.config.equilibrium_reactions.r1.parameter_data = {"k_eq_ref": 1}

        ConstantKeq.build_parameters(
            model.rparams.reaction_r1, model.rparams.config.equilibrium_reactions["r1"]
        )

        # Check parameter construction
        assert isinstance(model.rparams.reaction_r1.k_eq_ref, Var)
        assert model.rparams.reaction_r1.k_eq_ref.value == 1

        assert not hasattr(model.rparams.reaction_r1, "T_eq_ref")

        # Check expression
        rform = ConstantKeq.return_expression(
            model.rxn[1], model.rparams.reaction_r1, "r1", 300 * pyunits.K
        )

        assert value(rform) == 1
        assert_units_equivalent(rform, None)

    @pytest.mark.unit
    def test_ConstantKeq_molarity(self, model):
        model.rparams.config.equilibrium_reactions.r1.concentration_form = (
            ConcentrationForm.molarity
        )
        model.rparams.config.equilibrium_reactions.r1.parameter_data = {"k_eq_ref": 1}

        ConstantKeq.build_parameters(
            model.rparams.reaction_r1, model.rparams.config.equilibrium_reactions["r1"]
        )

        # Check parameter construction
        assert isinstance(model.rparams.reaction_r1.k_eq_ref, Var)
        assert model.rparams.reaction_r1.k_eq_ref.value == 1

        assert not hasattr(model.rparams.reaction_r1, "T_eq_ref")

        # Check expression
        rform = ConstantKeq.return_expression(
            model.rxn[1], model.rparams.reaction_r1, "r1", 300 * pyunits.K
        )

        assert value(rform) == 1
        assert_units_equivalent(rform, pyunits.mol / pyunits.m**3)

    @pytest.mark.unit
    def test_ConstantKeq_molarity_convert(self, model):
        model.rparams.config.equilibrium_reactions.r1.concentration_form = (
            ConcentrationForm.molarity
        )
        model.rparams.config.equilibrium_reactions.r1.parameter_data = {
            "k_eq_ref": (1e-3, pyunits.kmol / pyunits.m**3)
        }

        ConstantKeq.build_parameters(
            model.rparams.reaction_r1, model.rparams.config.equilibrium_reactions["r1"]
        )

        # Check parameter construction
        assert isinstance(model.rparams.reaction_r1.k_eq_ref, Var)
        assert model.rparams.reaction_r1.k_eq_ref.value == 1

        assert not hasattr(model.rparams.reaction_r1, "T_eq_ref")

        # Check expression
        rform = ConstantKeq.return_expression(
            model.rxn[1], model.rparams.reaction_r1, "r1", 300 * pyunits.K
        )

        assert value(rform) == 1
        assert_units_equivalent(rform, pyunits.mol / pyunits.m**3)

    @pytest.mark.unit
    def test_ConstantKeq_molality(self, model):
        model.rparams.config.equilibrium_reactions.r1.concentration_form = (
            ConcentrationForm.molality
        )
        model.rparams.config.equilibrium_reactions.r1.parameter_data = {"k_eq_ref": 1}

        ConstantKeq.build_parameters(
            model.rparams.reaction_r1, model.rparams.config.equilibrium_reactions["r1"]
        )

        # Check parameter construction
        assert isinstance(model.rparams.reaction_r1.k_eq_ref, Var)
        assert model.rparams.reaction_r1.k_eq_ref.value == 1

        assert not hasattr(model.rparams.reaction_r1, "T_eq_ref")

        # Check expression
        rform = ConstantKeq.return_expression(
            model.rxn[1], model.rparams.reaction_r1, "r1", 300 * pyunits.K
        )

        assert value(rform) == 1
        assert_units_equivalent(rform, pyunits.mol / pyunits.kg)

    @pytest.mark.unit
    def test_ConstantKeq_molality_convert(self, model):
        model.rparams.config.equilibrium_reactions.r1.concentration_form = (
            ConcentrationForm.molality
        )
        model.rparams.config.equilibrium_reactions.r1.parameter_data = {
            "k_eq_ref": (1e-3, pyunits.kmol / pyunits.kg)
        }

        ConstantKeq.build_parameters(
            model.rparams.reaction_r1, model.rparams.config.equilibrium_reactions["r1"]
        )

        # Check parameter construction
        assert isinstance(model.rparams.reaction_r1.k_eq_ref, Var)
        assert model.rparams.reaction_r1.k_eq_ref.value == 1

        assert not hasattr(model.rparams.reaction_r1, "T_eq_ref")

        # Check expression
        rform = ConstantKeq.return_expression(
            model.rxn[1], model.rparams.reaction_r1, "r1", 300 * pyunits.K
        )

        assert value(rform) == 1
        assert_units_equivalent(rform, pyunits.mol / pyunits.kg)

    @pytest.mark.unit
    def test_ConstantKeq_partial_pressure(self, model):
        model.rparams.config.equilibrium_reactions.r1.concentration_form = (
            ConcentrationForm.partialPressure
        )
        model.rparams.config.equilibrium_reactions.r1.parameter_data = {"k_eq_ref": 1}

        ConstantKeq.build_parameters(
            model.rparams.reaction_r1, model.rparams.config.equilibrium_reactions["r1"]
        )

        # Check parameter construction
        assert isinstance(model.rparams.reaction_r1.k_eq_ref, Var)
        assert model.rparams.reaction_r1.k_eq_ref.value == 1

        assert not hasattr(model.rparams.reaction_r1, "T_eq_ref")

        # Check expression
        rform = ConstantKeq.return_expression(
            model.rxn[1], model.rparams.reaction_r1, "r1", 300 * pyunits.K
        )

        assert value(rform) == 1
        assert_units_equivalent(rform, pyunits.Pa)

    @pytest.mark.unit
    def test_ConstantKeq_partial_pressure_convert(self, model):
        model.rparams.config.equilibrium_reactions.r1.concentration_form = (
            ConcentrationForm.partialPressure
        )
        model.rparams.config.equilibrium_reactions.r1.parameter_data = {
            "k_eq_ref": (1e-3, pyunits.kPa)
        }

        ConstantKeq.build_parameters(
            model.rparams.reaction_r1, model.rparams.config.equilibrium_reactions["r1"]
        )

        # Check parameter construction
        assert isinstance(model.rparams.reaction_r1.k_eq_ref, Var)
        assert model.rparams.reaction_r1.k_eq_ref.value == 1

        assert not hasattr(model.rparams.reaction_r1, "T_eq_ref")

        # Check expression
        rform = ConstantKeq.return_expression(
            model.rxn[1], model.rparams.reaction_r1, "r1", 300 * pyunits.K
        )

        assert value(rform) == 1
        assert_units_equivalent(rform, pyunits.Pa)


class TestVanTHoff(object):
    @pytest.mark.unit
    def test_van_t_hoff_mole_frac(self, model):
        model.rparams.config.equilibrium_reactions.r1.parameter_data = {
            "k_eq_ref": 1,
            "T_eq_ref": 500,
        }

        van_t_hoff.build_parameters(
            model.rparams.reaction_r1, model.rparams.config.equilibrium_reactions["r1"]
        )

        # Check parameter construction
        assert isinstance(model.rparams.reaction_r1.k_eq_ref, Var)
        assert model.rparams.reaction_r1.k_eq_ref.value == 1

        assert isinstance(model.rparams.reaction_r1.T_eq_ref, Var)
        assert model.rparams.reaction_r1.T_eq_ref.value == 500

        # Check expression
        units = model.rparams.get_metadata().derived_units

        model.rxn[1].log_k_eq = Var(["r1"], initialize=0)
        rform1 = van_t_hoff.return_log_expression(
            model.rxn[1], model.rparams.reaction_r1, "r1", 300 * pyunits.K
        )
        assert str(rform1) == str(
            (model.rxn[1].log_k_eq["r1"] - log(model.rparams.reaction_r1.k_eq_ref))
            == (
                -model.rxn[1].dh_rxn["r1"]
                / pyunits.convert(c.gas_constant, to_units=units.GAS_CONSTANT)
                * (1 / (300 * pyunits.K) - 1 / model.rparams.reaction_r1.T_eq_ref)
            )
        )

        rform2 = van_t_hoff.return_expression(
            model.rxn[1], model.rparams.reaction_r1, "r1", 300 * pyunits.K
        )
        assert value(rform2) == pytest.approx(1, rel=1e-3)
        assert_units_equivalent(rform1, None)
        assert_units_equivalent(rform2, None)

    @pytest.mark.unit
    def test_van_t_hoff_mole_frac_convert(self, model):
        model.rparams.config.equilibrium_reactions.r1.parameter_data = {
            "k_eq_ref": (1, None),
            "T_eq_ref": (900, pyunits.degR),
        }

        van_t_hoff.build_parameters(
            model.rparams.reaction_r1, model.rparams.config.equilibrium_reactions["r1"]
        )

        # Check parameter construction
        assert isinstance(model.rparams.reaction_r1.k_eq_ref, Var)
        assert model.rparams.reaction_r1.k_eq_ref.value == 1

        assert isinstance(model.rparams.reaction_r1.T_eq_ref, Var)
        assert model.rparams.reaction_r1.T_eq_ref.value == 500

        # Check expression
        units = model.rparams.get_metadata().derived_units

        model.rxn[1].log_k_eq = Var(["r1"], initialize=0)
        rform1 = van_t_hoff.return_log_expression(
            model.rxn[1], model.rparams.reaction_r1, "r1", 300 * pyunits.K
        )
        assert str(rform1) == str(
            (model.rxn[1].log_k_eq["r1"] - log(model.rparams.reaction_r1.k_eq_ref))
            == (
                -model.rxn[1].dh_rxn["r1"]
                / pyunits.convert(c.gas_constant, to_units=units.GAS_CONSTANT)
                * (1 / (300 * pyunits.K) - 1 / model.rparams.reaction_r1.T_eq_ref)
            )
        )

        rform2 = van_t_hoff.return_expression(
            model.rxn[1], model.rparams.reaction_r1, "r1", 300 * pyunits.K
        )
        assert value(rform2) == pytest.approx(1, rel=1e-3)
        assert_units_equivalent(rform1, None)
        assert_units_equivalent(rform2, None)

    @pytest.mark.unit
    def test_van_t_hoff_molarity(self, model):
        model.rparams.config.equilibrium_reactions.r1.concentration_form = (
            ConcentrationForm.molarity
        )
        model.rparams.config.equilibrium_reactions.r1.parameter_data = {
            "k_eq_ref": 1,
            "T_eq_ref": 500,
        }

        van_t_hoff.build_parameters(
            model.rparams.reaction_r1, model.rparams.config.equilibrium_reactions["r1"]
        )

        # Check parameter construction
        assert isinstance(model.rparams.reaction_r1.k_eq_ref, Var)
        assert model.rparams.reaction_r1.k_eq_ref.value == 1

        assert isinstance(model.rparams.reaction_r1.T_eq_ref, Var)
        assert model.rparams.reaction_r1.T_eq_ref.value == 500

        # Check expression
        units = model.rparams.get_metadata().derived_units

        model.rxn[1].log_k_eq = Var(["r1"], initialize=0)
        rform1 = van_t_hoff.return_log_expression(
            model.rxn[1], model.rparams.reaction_r1, "r1", 300 * pyunits.K
        )
        assert str(rform1) == str(
            (
                model.rxn[1].log_k_eq["r1"]
                - log(
                    1
                    / (pyunits.mol / pyunits.m**3)
                    * model.rparams.reaction_r1.k_eq_ref
                )
            )
            == (
                -model.rxn[1].dh_rxn["r1"]
                / pyunits.convert(c.gas_constant, to_units=units.GAS_CONSTANT)
                * (1 / (300 * pyunits.K) - 1 / model.rparams.reaction_r1.T_eq_ref)
            )
        )

        rform2 = van_t_hoff.return_expression(
            model.rxn[1], model.rparams.reaction_r1, "r1", 300 * pyunits.K
        )
        assert value(rform2) == pytest.approx(1, rel=1e-3)
        assert_units_equivalent(rform1, None)
        assert_units_equivalent(rform2, pyunits.mol / pyunits.m**3)

    @pytest.mark.unit
    def test_van_t_hoff_molarity_convert(self, model):
        model.rparams.config.equilibrium_reactions.r1.concentration_form = (
            ConcentrationForm.molarity
        )
        model.rparams.config.equilibrium_reactions.r1.parameter_data = {
            "k_eq_ref": (1e-3, pyunits.kmol / pyunits.m**3),
            "T_eq_ref": (900, pyunits.degR),
        }

        van_t_hoff.build_parameters(
            model.rparams.reaction_r1, model.rparams.config.equilibrium_reactions["r1"]
        )

        # Check parameter construction
        assert isinstance(model.rparams.reaction_r1.k_eq_ref, Var)
        assert model.rparams.reaction_r1.k_eq_ref.value == 1

        assert isinstance(model.rparams.reaction_r1.T_eq_ref, Var)
        assert model.rparams.reaction_r1.T_eq_ref.value == 500

        # Check expression
        units = model.rparams.get_metadata().derived_units

        model.rxn[1].log_k_eq = Var(["r1"], initialize=0)
        rform1 = van_t_hoff.return_log_expression(
            model.rxn[1], model.rparams.reaction_r1, "r1", 300 * pyunits.K
        )
        assert str(rform1) == str(
            (
                model.rxn[1].log_k_eq["r1"]
                - log(
                    1
                    / (pyunits.mol / pyunits.m**3)
                    * model.rparams.reaction_r1.k_eq_ref
                )
            )
            == (
                -model.rxn[1].dh_rxn["r1"]
                / pyunits.convert(c.gas_constant, to_units=units.GAS_CONSTANT)
                * (1 / (300 * pyunits.K) - 1 / model.rparams.reaction_r1.T_eq_ref)
            )
        )

        rform2 = van_t_hoff.return_expression(
            model.rxn[1], model.rparams.reaction_r1, "r1", 300 * pyunits.K
        )
        assert value(rform2) == pytest.approx(1, rel=1e-3)
        assert_units_equivalent(rform1, None)
        assert_units_equivalent(rform2, pyunits.mol / pyunits.m**3)

    @pytest.mark.unit
    def test_van_t_hoff_molality(self, model):
        model.rparams.config.equilibrium_reactions.r1.concentration_form = (
            ConcentrationForm.molality
        )
        model.rparams.config.equilibrium_reactions.r1.parameter_data = {
            "k_eq_ref": 1,
            "T_eq_ref": 500,
        }

        van_t_hoff.build_parameters(
            model.rparams.reaction_r1, model.rparams.config.equilibrium_reactions["r1"]
        )

        # Check parameter construction
        assert isinstance(model.rparams.reaction_r1.k_eq_ref, Var)
        assert model.rparams.reaction_r1.k_eq_ref.value == 1

        assert isinstance(model.rparams.reaction_r1.T_eq_ref, Var)
        assert model.rparams.reaction_r1.T_eq_ref.value == 500

        # Check expression
        units = model.rparams.get_metadata().derived_units

        model.rxn[1].log_k_eq = Var(["r1"], initialize=0)
        rform1 = van_t_hoff.return_log_expression(
            model.rxn[1], model.rparams.reaction_r1, "r1", 300 * pyunits.K
        )
        assert str(rform1) == str(
            (
                model.rxn[1].log_k_eq["r1"]
                - log(
                    1 / (pyunits.mol / pyunits.kg) * model.rparams.reaction_r1.k_eq_ref
                )
            )
            == (
                -model.rxn[1].dh_rxn["r1"]
                / pyunits.convert(c.gas_constant, to_units=units.GAS_CONSTANT)
                * (1 / (300 * pyunits.K) - 1 / model.rparams.reaction_r1.T_eq_ref)
            )
        )

        rform2 = van_t_hoff.return_expression(
            model.rxn[1], model.rparams.reaction_r1, "r1", 300 * pyunits.K
        )
        assert value(rform2) == pytest.approx(1, rel=1e-3)
        assert_units_equivalent(rform1, None)
        assert_units_equivalent(rform2, pyunits.mol / pyunits.kg)

    @pytest.mark.unit
    def test_van_t_hoff_molality_convert(self, model):
        model.rparams.config.equilibrium_reactions.r1.concentration_form = (
            ConcentrationForm.molality
        )
        model.rparams.config.equilibrium_reactions.r1.parameter_data = {
            "k_eq_ref": (1e-3, pyunits.kmol / pyunits.kg),
            "T_eq_ref": (900, pyunits.degR),
        }

        van_t_hoff.build_parameters(
            model.rparams.reaction_r1, model.rparams.config.equilibrium_reactions["r1"]
        )

        # Check parameter construction
        assert isinstance(model.rparams.reaction_r1.k_eq_ref, Var)
        assert model.rparams.reaction_r1.k_eq_ref.value == 1

        assert isinstance(model.rparams.reaction_r1.T_eq_ref, Var)
        assert model.rparams.reaction_r1.T_eq_ref.value == 500

        # Check expression
        units = model.rparams.get_metadata().derived_units

        model.rxn[1].log_k_eq = Var(["r1"], initialize=0)
        rform1 = van_t_hoff.return_log_expression(
            model.rxn[1], model.rparams.reaction_r1, "r1", 300 * pyunits.K
        )
        assert str(rform1) == str(
            (
                model.rxn[1].log_k_eq["r1"]
                - log(
                    1 / (pyunits.mol / pyunits.kg) * model.rparams.reaction_r1.k_eq_ref
                )
            )
            == (
                -model.rxn[1].dh_rxn["r1"]
                / pyunits.convert(c.gas_constant, to_units=units.GAS_CONSTANT)
                * (1 / (300 * pyunits.K) - 1 / model.rparams.reaction_r1.T_eq_ref)
            )
        )

        rform2 = van_t_hoff.return_expression(
            model.rxn[1], model.rparams.reaction_r1, "r1", 300 * pyunits.K
        )
        assert value(rform2) == pytest.approx(1, rel=1e-3)
        assert_units_equivalent(rform1, None)
        assert_units_equivalent(rform2, pyunits.mol / pyunits.kg)

    @pytest.mark.unit
    def test_van_t_hoff_partial_pressure(self, model):
        model.rparams.config.equilibrium_reactions.r1.concentration_form = (
            ConcentrationForm.partialPressure
        )
        model.rparams.config.equilibrium_reactions.r1.parameter_data = {
            "k_eq_ref": 1,
            "T_eq_ref": 500,
        }

        van_t_hoff.build_parameters(
            model.rparams.reaction_r1, model.rparams.config.equilibrium_reactions["r1"]
        )

        # Check parameter construction
        assert isinstance(model.rparams.reaction_r1.k_eq_ref, Var)
        assert model.rparams.reaction_r1.k_eq_ref.value == 1

        assert isinstance(model.rparams.reaction_r1.T_eq_ref, Var)
        assert model.rparams.reaction_r1.T_eq_ref.value == 500

        # Check expression
        units = model.rparams.get_metadata().derived_units

        model.rxn[1].log_k_eq = Var(["r1"], initialize=0)
        rform1 = van_t_hoff.return_log_expression(
            model.rxn[1], model.rparams.reaction_r1, "r1", 300 * pyunits.K
        )
        assert str(rform1) == str(
            (
                model.rxn[1].log_k_eq["r1"]
                - log(
                    1
                    / (pyunits.kg / pyunits.m / pyunits.s**2)
                    * model.rparams.reaction_r1.k_eq_ref
                )
            )
            == (
                -model.rxn[1].dh_rxn["r1"]
                / pyunits.convert(c.gas_constant, to_units=units.GAS_CONSTANT)
                * (1 / (300 * pyunits.K) - 1 / model.rparams.reaction_r1.T_eq_ref)
            )
        )

        rform2 = van_t_hoff.return_expression(
            model.rxn[1], model.rparams.reaction_r1, "r1", 300 * pyunits.K
        )
        assert value(rform2) == pytest.approx(1, rel=1e-3)
        assert_units_equivalent(rform1, None)
        assert_units_equivalent(rform2, pyunits.Pa)

    @pytest.mark.unit
    def test_van_t_hoff_partial_pressure_convert(self, model):
        model.rparams.config.equilibrium_reactions.r1.concentration_form = (
            ConcentrationForm.partialPressure
        )
        model.rparams.config.equilibrium_reactions.r1.parameter_data = {
            "k_eq_ref": (1e-3, pyunits.kPa),
            "T_eq_ref": (900, pyunits.degR),
        }

        van_t_hoff.build_parameters(
            model.rparams.reaction_r1, model.rparams.config.equilibrium_reactions["r1"]
        )

        # Check parameter construction
        assert isinstance(model.rparams.reaction_r1.k_eq_ref, Var)
        assert model.rparams.reaction_r1.k_eq_ref.value == 1

        assert isinstance(model.rparams.reaction_r1.T_eq_ref, Var)
        assert model.rparams.reaction_r1.T_eq_ref.value == 500

        # Check expression
        units = model.rparams.get_metadata().derived_units

        model.rxn[1].log_k_eq = Var(["r1"], initialize=0)
        rform1 = van_t_hoff.return_log_expression(
            model.rxn[1], model.rparams.reaction_r1, "r1", 300 * pyunits.K
        )
        assert str(rform1) == str(
            (
                model.rxn[1].log_k_eq["r1"]
                - log(
                    1
                    / (pyunits.kg / pyunits.m / pyunits.s**2)
                    * model.rparams.reaction_r1.k_eq_ref
                )
            )
            == (
                -model.rxn[1].dh_rxn["r1"]
                / pyunits.convert(c.gas_constant, to_units=units.GAS_CONSTANT)
                * (1 / (300 * pyunits.K) - 1 / model.rparams.reaction_r1.T_eq_ref)
            )
        )

        rform2 = van_t_hoff.return_expression(
            model.rxn[1], model.rparams.reaction_r1, "r1", 300 * pyunits.K
        )
        assert value(rform2) == pytest.approx(1, rel=1e-3)
        assert_units_equivalent(rform1, None)
        assert_units_equivalent(rform2, pyunits.Pa)


class TestGibbsEnergy(object):
    @pytest.mark.unit
    def test_gibbs_energy_mole_frac(self, model):
        model.rparams.config.equilibrium_reactions.r1.heat_of_reaction = constant_dh_rxn
        model.rparams.config.equilibrium_reactions.r1.parameter_data = {
            "ds_rxn_ref": 1,
            "T_eq_ref": 500,
        }
        model.rparams.reaction_r1.dh_rxn_ref = Var(
            initialize=2, units=pyunits.J / pyunits.mol
        )

        gibbs_energy.build_parameters(
            model.rparams.reaction_r1, model.rparams.config.equilibrium_reactions["r1"]
        )

        # Check parameter construction
        assert isinstance(model.rparams.reaction_r1.ds_rxn_ref, Var)
        assert model.rparams.reaction_r1.ds_rxn_ref.value == 1

        assert isinstance(model.rparams.reaction_r1.T_eq_ref, Var)
        assert model.rparams.reaction_r1.T_eq_ref.value == 500

        # Check expression
        model.rxn[1].log_k_eq = Var(["r1"], initialize=0)

        rform1 = gibbs_energy.return_expression(
            model.rxn[1], model.rparams.reaction_r1, "r1", 300 * pyunits.K
        )
        rform2 = gibbs_energy.return_log_expression(
            model.rxn[1], model.rparams.reaction_r1, "r1", 300 * pyunits.K
        )

        assert str(rform1) == ("exp(rxn[1].log_k_eq[r1])")
        assert str(rform2) == (
            "rxn[1].log_k_eq[r1]  ==  "
            "- rparams.reaction_r1.dh_rxn_ref/("
            + str(
                pyunits.convert(
                    c.gas_constant,
                    to_units=pyunits.kg
                    * pyunits.m**2
                    / pyunits.s**2
                    / pyunits.mol
                    / pyunits.K,
                )
            )
            + "*(300*K)) + 1/("
            + str(
                pyunits.convert(
                    c.gas_constant,
                    to_units=pyunits.kg
                    * pyunits.m**2
                    / pyunits.s**2
                    / pyunits.mol
                    / pyunits.K,
                )
            )
            + ")*rparams.reaction_r1.ds_rxn_ref"
        )

        assert value(rform1) == pytest.approx(1, rel=1e-3)
        assert_units_equivalent(rform1, pyunits.dimensionless)
        assert_units_equivalent(rform2, pyunits.dimensionless)

    @pytest.mark.unit
    def test_gibbs_energy_invalid_dh_rxn(self, model):
        model.rparams.config.equilibrium_reactions.r1.parameter_data = {
            "ds_rxn_ref": 1,
            "T_eq_ref": 500,
        }
        model.rparams.reaction_r1.dh_rxn_ref = Var(
            initialize=2, units=pyunits.J / pyunits.mol
        )

        with pytest.raises(
            ConfigurationError,
            match="rparams.reaction_r1 calculating equilibrium "
            "constants from Gibbs energy assumes constant "
            "heat of reaction. Please ensure you are using the "
            "constant_dh_rxn method for this reaction",
        ):
            gibbs_energy.build_parameters(
                model.rparams.reaction_r1,
                model.rparams.config.equilibrium_reactions["r1"],
            )

    @pytest.mark.unit
    def test_gibbs_energy_molarity(self, model):
        model.rparams.config.equilibrium_reactions.r1.concentration_form = (
            ConcentrationForm.molarity
        )
        model.rparams.config.equilibrium_reactions.r1.heat_of_reaction = constant_dh_rxn
        model.rparams.config.equilibrium_reactions.r1.parameter_data = {
            "ds_rxn_ref": 1,
            "T_eq_ref": 500,
        }
        model.rparams.reaction_r1.dh_rxn_ref = Var(
            initialize=2, units=pyunits.J / pyunits.mol
        )

        with pytest.raises(
            ConfigurationError,
            match="rparams.reaction_r1 calculation of equilibrium constant"
            " based on Gibbs energy is only supported for mole fraction or"
            " activity forms. Currently selected form: "
            "ConcentrationForm.molarity",
        ):
            gibbs_energy.build_parameters(
                model.rparams.reaction_r1,
                model.rparams.config.equilibrium_reactions["r1"],
            )

    @pytest.mark.unit
    def test_gibbs_energy_molality(self, model):
        model.rparams.config.equilibrium_reactions.r1.concentration_form = (
            ConcentrationForm.molality
        )
        model.rparams.config.equilibrium_reactions.r1.heat_of_reaction = constant_dh_rxn
        model.rparams.config.equilibrium_reactions.r1.parameter_data = {
            "ds_rxn_ref": 1,
            "T_eq_ref": 500,
        }
        model.rparams.reaction_r1.dh_rxn_ref = Var(
            initialize=2, units=pyunits.J / pyunits.mol
        )

        with pytest.raises(
            ConfigurationError,
            match="rparams.reaction_r1 calculation of equilibrium constant"
            " based on Gibbs energy is only supported for mole fraction or"
            " activity forms. Currently selected form: "
            "ConcentrationForm.molality",
        ):
            gibbs_energy.build_parameters(
                model.rparams.reaction_r1,
                model.rparams.config.equilibrium_reactions["r1"],
            )

    @pytest.mark.unit
    def test_gibbs_energy_partial_pressure(self, model):
        model.rparams.config.equilibrium_reactions.r1.concentration_form = (
            ConcentrationForm.partialPressure
        )
        model.rparams.config.equilibrium_reactions.r1.heat_of_reaction = constant_dh_rxn
        model.rparams.config.equilibrium_reactions.r1.parameter_data = {
            "ds_rxn_ref": 1,
            "T_eq_ref": 500,
        }
        model.rparams.reaction_r1.dh_rxn_ref = Var(
            initialize=2, units=pyunits.J / pyunits.mol
        )

        with pytest.raises(
            ConfigurationError,
            match="rparams.reaction_r1 calculation of equilibrium constant"
            " based on Gibbs energy is only supported for mole fraction or"
            " activity forms. Currently selected form: "
            "ConcentrationForm.partialPressure",
        ):
            gibbs_energy.build_parameters(
                model.rparams.reaction_r1,
                model.rparams.config.equilibrium_reactions["r1"],
            )
