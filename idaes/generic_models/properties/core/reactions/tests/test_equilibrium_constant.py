##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
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
Tests for rate forms
"""

import pytest

from pyomo.environ import ConcreteModel, Var, units as pyunits, value
from pyomo.util.check_units import assert_units_equivalent

from idaes.generic_models.properties.core.generic.generic_reaction import \
    GenericReactionParameterBlock, ConcentrationForm
from idaes.generic_models.properties.core.reactions.equilibrium_constant import *
from idaes.core import MaterialFlowBasis
from idaes.core.util.testing import PhysicalParameterTestBlock


@pytest.fixture
def model():
    m = ConcreteModel()

    # # Add a test thermo package for validation
    m.pparams = PhysicalParameterTestBlock()
    m.thermo = m.pparams.build_state_block([1])

    m.rparams = GenericReactionParameterBlock(default={
         "property_package": m.pparams,
         "reaction_basis": MaterialFlowBasis.molar,
         "equilibrium_reactions": {
             "r1": {"stoichiometry": {("p1", "c1"): -1,
                                      ("p1", "c2"): 2},
                    "equilibrium_form": "foo",
                    "concentration_form": ConcentrationForm.moleFraction}},
         "base_units": {"amount": pyunits.mol,
                        "mass": pyunits.kg,
                        "time": pyunits.s,
                        "length": pyunits.m,
                        "temperature": pyunits.K}})

    m.rxn = m.rparams.build_reaction_block([1], default={
        "state_block": m.thermo, "has_equilibrium": False})

    m.rxn[1].dh_rxn = Var(["r1"],
                          initialize=1,
                          units=pyunits.J/pyunits.mol)

    return m


class TestVanTHoff(object):
    @pytest.mark.unit
    def test_van_t_hoff_mole_frac(self, model):
        model.rparams.config.equilibrium_reactions.r1.parameter_data = {
            "k_eq_ref": 1,
            "T_eq_ref": 500}

        van_t_hoff.build_parameters(
            model.rparams.reaction_r1,
            model.rparams.config.equilibrium_reactions["r1"])

        # Check parameter construction
        assert isinstance(model.rparams.reaction_r1.k_eq_ref, Var)
        assert model.rparams.reaction_r1.k_eq_ref.value == 1

        assert isinstance(model.rparams.reaction_r1.T_eq_ref, Var)
        assert model.rparams.reaction_r1.T_eq_ref.value == 500

        # Check expression
        rform = van_t_hoff.return_expression(
            model.rxn[1], model.rparams.reaction_r1, "r1", 300*pyunits.K)

        assert value(rform) == pytest.approx(0.99984, rel=1e-3)
        assert_units_equivalent(rform, None)

    @pytest.mark.unit
    def test_van_t_hoff_mole_frac_convert(self, model):
        model.rparams.config.equilibrium_reactions.r1.parameter_data = {
            "k_eq_ref": (1, None),
            "T_eq_ref": (900, pyunits.degR)}

        van_t_hoff.build_parameters(
            model.rparams.reaction_r1,
            model.rparams.config.equilibrium_reactions["r1"])

        # Check parameter construction
        assert isinstance(model.rparams.reaction_r1.k_eq_ref, Var)
        assert model.rparams.reaction_r1.k_eq_ref.value == 1

        assert isinstance(model.rparams.reaction_r1.T_eq_ref, Var)
        assert model.rparams.reaction_r1.T_eq_ref.value == 500

        # Check expression
        rform = van_t_hoff.return_expression(
            model.rxn[1], model.rparams.reaction_r1, "r1", 300*pyunits.K)

        assert value(rform) == pytest.approx(0.99984, rel=1e-3)
        assert_units_equivalent(rform, None)

    @pytest.mark.unit
    def test_van_t_hoff_molarity(self, model):
        model.rparams.config.equilibrium_reactions.r1.concentration_form = \
            ConcentrationForm.molarity
        model.rparams.config.equilibrium_reactions.r1.parameter_data = {
            "k_eq_ref": 1,
            "T_eq_ref": 500}

        van_t_hoff.build_parameters(
            model.rparams.reaction_r1,
            model.rparams.config.equilibrium_reactions["r1"])

        # Check parameter construction
        assert isinstance(model.rparams.reaction_r1.k_eq_ref, Var)
        assert model.rparams.reaction_r1.k_eq_ref.value == 1

        assert isinstance(model.rparams.reaction_r1.T_eq_ref, Var)
        assert model.rparams.reaction_r1.T_eq_ref.value == 500

        # Check expression
        rform = van_t_hoff.return_expression(
            model.rxn[1], model.rparams.reaction_r1, "r1", 300*pyunits.K)

        assert value(rform) == pytest.approx(0.99984, rel=1e-3)
        assert_units_equivalent(rform, pyunits.mol/pyunits.m**3)

    @pytest.mark.unit
    def test_van_t_hoff_molarity_convert(self, model):
        model.rparams.config.equilibrium_reactions.r1.concentration_form = \
            ConcentrationForm.molarity
        model.rparams.config.equilibrium_reactions.r1.parameter_data = {
            "k_eq_ref": (1e-3, pyunits.kmol/pyunits.m**3),
            "T_eq_ref": (900, pyunits.degR)}

        van_t_hoff.build_parameters(
            model.rparams.reaction_r1,
            model.rparams.config.equilibrium_reactions["r1"])

        # Check parameter construction
        assert isinstance(model.rparams.reaction_r1.k_eq_ref, Var)
        assert model.rparams.reaction_r1.k_eq_ref.value == 1

        assert isinstance(model.rparams.reaction_r1.T_eq_ref, Var)
        assert model.rparams.reaction_r1.T_eq_ref.value == 500

        # Check expression
        rform = van_t_hoff.return_expression(
            model.rxn[1], model.rparams.reaction_r1, "r1", 300*pyunits.K)

        assert value(rform) == pytest.approx(0.99984, rel=1e-3)
        assert_units_equivalent(rform, pyunits.mol/pyunits.m**3)

    @pytest.mark.unit
    def test_van_t_hoff_molality(self, model):
        model.rparams.config.equilibrium_reactions.r1.concentration_form = \
            ConcentrationForm.molality
        model.rparams.config.equilibrium_reactions.r1.parameter_data = {
            "k_eq_ref": 1,
            "T_eq_ref": 500}

        van_t_hoff.build_parameters(
            model.rparams.reaction_r1,
            model.rparams.config.equilibrium_reactions["r1"])

        # Check parameter construction
        assert isinstance(model.rparams.reaction_r1.k_eq_ref, Var)
        assert model.rparams.reaction_r1.k_eq_ref.value == 1

        assert isinstance(model.rparams.reaction_r1.T_eq_ref, Var)
        assert model.rparams.reaction_r1.T_eq_ref.value == 500

        # Check expression
        rform = van_t_hoff.return_expression(
            model.rxn[1], model.rparams.reaction_r1, "r1", 300*pyunits.K)

        assert value(rform) == pytest.approx(0.99984, rel=1e-3)
        assert_units_equivalent(rform, pyunits.mol/pyunits.kg)

    @pytest.mark.unit
    def test_van_t_hoff_molality_convert(self, model):
        model.rparams.config.equilibrium_reactions.r1.concentration_form = \
            ConcentrationForm.molality
        model.rparams.config.equilibrium_reactions.r1.parameter_data = {
            "k_eq_ref": (1e-3, pyunits.kmol/pyunits.kg),
            "T_eq_ref": (900, pyunits.degR)}

        van_t_hoff.build_parameters(
            model.rparams.reaction_r1,
            model.rparams.config.equilibrium_reactions["r1"])

        # Check parameter construction
        assert isinstance(model.rparams.reaction_r1.k_eq_ref, Var)
        assert model.rparams.reaction_r1.k_eq_ref.value == 1

        assert isinstance(model.rparams.reaction_r1.T_eq_ref, Var)
        assert model.rparams.reaction_r1.T_eq_ref.value == 500

        # Check expression
        rform = van_t_hoff.return_expression(
            model.rxn[1], model.rparams.reaction_r1, "r1", 300*pyunits.K)

        assert value(rform) == pytest.approx(0.99984, rel=1e-3)
        assert_units_equivalent(rform, pyunits.mol/pyunits.kg)

    @pytest.mark.unit
    def test_van_t_hoff_partial_pressure(self, model):
        model.rparams.config.equilibrium_reactions.r1.concentration_form = \
            ConcentrationForm.partialPressure
        model.rparams.config.equilibrium_reactions.r1.parameter_data = {
            "k_eq_ref": 1,
            "T_eq_ref": 500}

        van_t_hoff.build_parameters(
            model.rparams.reaction_r1,
            model.rparams.config.equilibrium_reactions["r1"])

        # Check parameter construction
        assert isinstance(model.rparams.reaction_r1.k_eq_ref, Var)
        assert model.rparams.reaction_r1.k_eq_ref.value == 1

        assert isinstance(model.rparams.reaction_r1.T_eq_ref, Var)
        assert model.rparams.reaction_r1.T_eq_ref.value == 500

        # Check expression
        rform = van_t_hoff.return_expression(
            model.rxn[1], model.rparams.reaction_r1, "r1", 300*pyunits.K)

        assert value(rform) == pytest.approx(0.99984, rel=1e-3)
        assert_units_equivalent(rform, pyunits.Pa)

    @pytest.mark.unit
    def test_van_t_hoff_partial_pressure_convert(self, model):
        model.rparams.config.equilibrium_reactions.r1.concentration_form = \
            ConcentrationForm.partialPressure
        model.rparams.config.equilibrium_reactions.r1.parameter_data = {
            "k_eq_ref": (1e-3, pyunits.kPa),
            "T_eq_ref": (900, pyunits.degR)}

        van_t_hoff.build_parameters(
            model.rparams.reaction_r1,
            model.rparams.config.equilibrium_reactions["r1"])

        # Check parameter construction
        assert isinstance(model.rparams.reaction_r1.k_eq_ref, Var)
        assert model.rparams.reaction_r1.k_eq_ref.value == 1

        assert isinstance(model.rparams.reaction_r1.T_eq_ref, Var)
        assert model.rparams.reaction_r1.T_eq_ref.value == 500

        # Check expression
        rform = van_t_hoff.return_expression(
            model.rxn[1], model.rparams.reaction_r1, "r1", 300*pyunits.K)

        assert value(rform) == pytest.approx(0.99984, rel=1e-3)
        assert_units_equivalent(rform, pyunits.Pa)


class TestGibbsEnergy(object):
    @pytest.mark.unit
    def test_gibbs_energy_mole_frac(self, model):
        model.rparams.config.equilibrium_reactions.r1.parameter_data = {
            "dh_rxn_ref": 2,
            "ds_rxn_ref": 1,
            "T_eq_ref": 500}

        gibbs_energy.build_parameters(
            model.rparams.reaction_r1,
            model.rparams.config.equilibrium_reactions["r1"])

        # Check parameter construction
        assert isinstance(model.rparams.reaction_r1.dh_rxn_ref, Var)
        assert model.rparams.reaction_r1.dh_rxn_ref.value == 2

        assert isinstance(model.rparams.reaction_r1.ds_rxn_ref, Var)
        assert model.rparams.reaction_r1.ds_rxn_ref.value == 1

        assert isinstance(model.rparams.reaction_r1.T_eq_ref, Var)
        assert model.rparams.reaction_r1.T_eq_ref.value == 500

        assert_units_equivalent(model.rparams.reaction_r1._keq_units,
                                pyunits.dimensionless)

        # Check expression
        rform = gibbs_energy.return_expression(
            model.rxn[1], model.rparams.reaction_r1, "r1", 300*pyunits.K)

        assert str(rform) == (
            'exp(- rparams.reaction_r1.dh_rxn_ref/(kg*m**2/J/s**2*'
            '(8.314462618*(J)/mol/K)*(300*K)) + '
            '1/(kg*m**2/J/s**2*(8.314462618*(J)/mol/K))*'
            'rparams.reaction_r1.ds_rxn_ref)*dimensionless')
        assert value(rform) == pytest.approx(1.1269, rel=1e-3)
        assert_units_equivalent(rform, None)

    @pytest.mark.unit
    def test_gibbs_energy_mole_frac_convert(self, model):
        model.rparams.config.equilibrium_reactions.r1.parameter_data = {
            "dh_rxn_ref": (2000, pyunits.J/pyunits.kmol),
            "ds_rxn_ref": (1000, pyunits.J/pyunits.kmol/pyunits.K),
            "T_eq_ref": (900, pyunits.degR)}

        gibbs_energy.build_parameters(
            model.rparams.reaction_r1,
            model.rparams.config.equilibrium_reactions["r1"])

        # Check parameter construction
        assert isinstance(model.rparams.reaction_r1.dh_rxn_ref, Var)
        assert model.rparams.reaction_r1.dh_rxn_ref.value == 2

        assert isinstance(model.rparams.reaction_r1.ds_rxn_ref, Var)
        assert model.rparams.reaction_r1.ds_rxn_ref.value == 1

        assert isinstance(model.rparams.reaction_r1.T_eq_ref, Var)
        assert model.rparams.reaction_r1.T_eq_ref.value == 500

        assert_units_equivalent(model.rparams.reaction_r1._keq_units,
                                pyunits.dimensionless)

        # Check expression
        rform = gibbs_energy.return_expression(
            model.rxn[1], model.rparams.reaction_r1, "r1", 300*pyunits.K)

        assert str(rform) == (
            'exp(- rparams.reaction_r1.dh_rxn_ref/(kg*m**2/J/s**2*'
            '(8.314462618*(J)/mol/K)*(300*K)) + '
            '1/(kg*m**2/J/s**2*(8.314462618*(J)/mol/K))*'
            'rparams.reaction_r1.ds_rxn_ref)*dimensionless')
        assert value(rform) == pytest.approx(1.1269, rel=1e-3)
        assert_units_equivalent(rform, None)

    @pytest.mark.unit
    def test_gibbs_energy_molarity(self, model):
        model.rparams.config.equilibrium_reactions.r1.concentration_form = \
            ConcentrationForm.molarity
        model.rparams.config.equilibrium_reactions.r1.parameter_data = {
            "dh_rxn_ref": 2,
            "ds_rxn_ref": 1,
            "T_eq_ref": 500}

        gibbs_energy.build_parameters(
            model.rparams.reaction_r1,
            model.rparams.config.equilibrium_reactions["r1"])

        # Check parameter construction
        assert isinstance(model.rparams.reaction_r1.dh_rxn_ref, Var)
        assert model.rparams.reaction_r1.dh_rxn_ref.value == 2

        assert isinstance(model.rparams.reaction_r1.ds_rxn_ref, Var)
        assert model.rparams.reaction_r1.ds_rxn_ref.value == 1

        assert isinstance(model.rparams.reaction_r1.T_eq_ref, Var)
        assert model.rparams.reaction_r1.T_eq_ref.value == 500

        assert_units_equivalent(model.rparams.reaction_r1._keq_units,
                                pyunits.mol*pyunits.m**-3)

        # Check expression
        rform = gibbs_energy.return_expression(
            model.rxn[1], model.rparams.reaction_r1, "r1", 300*pyunits.K)

        assert str(rform) == (
            'exp(- rparams.reaction_r1.dh_rxn_ref/(kg*m**2/J/s**2*'
            '(8.314462618*(J)/mol/K)*(300*K)) + '
            '1/(kg*m**2/J/s**2*(8.314462618*(J)/mol/K))*'
            'rparams.reaction_r1.ds_rxn_ref)*(mol*m**-3)')
        assert value(rform) == pytest.approx(1.1269, rel=1e-3)
        assert_units_equivalent(rform, pyunits.mol*pyunits.m**-3)

    @pytest.mark.unit
    def test_gibbs_energy_molarity_convert(self, model):
        model.rparams.config.equilibrium_reactions.r1.concentration_form = \
            ConcentrationForm.molarity
        model.rparams.config.equilibrium_reactions.r1.parameter_data = {
            "dh_rxn_ref": (2000, pyunits.J/pyunits.kmol),
            "ds_rxn_ref": (1000, pyunits.J/pyunits.kmol/pyunits.K),
            "T_eq_ref": (900, pyunits.degR)}

        gibbs_energy.build_parameters(
            model.rparams.reaction_r1,
            model.rparams.config.equilibrium_reactions["r1"])

        # Check parameter construction
        assert isinstance(model.rparams.reaction_r1.dh_rxn_ref, Var)
        assert model.rparams.reaction_r1.dh_rxn_ref.value == 2

        assert isinstance(model.rparams.reaction_r1.ds_rxn_ref, Var)
        assert model.rparams.reaction_r1.ds_rxn_ref.value == 1

        assert isinstance(model.rparams.reaction_r1.T_eq_ref, Var)
        assert model.rparams.reaction_r1.T_eq_ref.value == 500

        assert_units_equivalent(model.rparams.reaction_r1._keq_units,
                                pyunits.mol/pyunits.m**3)

        # Check expression
        rform = gibbs_energy.return_expression(
            model.rxn[1], model.rparams.reaction_r1, "r1", 300*pyunits.K)

        assert str(rform) == (
            'exp(- rparams.reaction_r1.dh_rxn_ref/(kg*m**2/J/s**2*'
            '(8.314462618*(J)/mol/K)*(300*K)) + '
            '1/(kg*m**2/J/s**2*(8.314462618*(J)/mol/K))*'
            'rparams.reaction_r1.ds_rxn_ref)*(mol*m**-3)')
        assert value(rform) == pytest.approx(1.1269, rel=1e-3)
        assert_units_equivalent(rform, pyunits.mol*pyunits.m**-3)

    @pytest.mark.unit
    def test_gibss_energy_molality(self, model):
        model.rparams.config.equilibrium_reactions.r1.concentration_form = \
            ConcentrationForm.molality
        model.rparams.config.equilibrium_reactions.r1.parameter_data = {
            "dh_rxn_ref": 2,
            "ds_rxn_ref": 1,
            "T_eq_ref": 500}

        gibbs_energy.build_parameters(
            model.rparams.reaction_r1,
            model.rparams.config.equilibrium_reactions["r1"])

        # Check parameter construction
        assert isinstance(model.rparams.reaction_r1.dh_rxn_ref, Var)
        assert model.rparams.reaction_r1.dh_rxn_ref.value == 2

        assert isinstance(model.rparams.reaction_r1.ds_rxn_ref, Var)
        assert model.rparams.reaction_r1.ds_rxn_ref.value == 1

        assert isinstance(model.rparams.reaction_r1.T_eq_ref, Var)
        assert model.rparams.reaction_r1.T_eq_ref.value == 500

        assert_units_equivalent(model.rparams.reaction_r1._keq_units,
                                pyunits.mol/pyunits.kg)

        # Check expression
        rform = gibbs_energy.return_expression(
            model.rxn[1], model.rparams.reaction_r1, "r1", 300*pyunits.K)

        assert str(rform) == (
            'exp(- rparams.reaction_r1.dh_rxn_ref/(kg*m**2/J/s**2*'
            '(8.314462618*(J)/mol/K)*(300*K)) + '
            '1/(kg*m**2/J/s**2*(8.314462618*(J)/mol/K))*'
            'rparams.reaction_r1.ds_rxn_ref)*(mol*kg**-1)')
        assert value(rform) == pytest.approx(1.1269, rel=1e-3)
        assert_units_equivalent(rform, pyunits.mol/pyunits.kg)

    @pytest.mark.unit
    def test_gibbs_energy_molality_convert(self, model):
        model.rparams.config.equilibrium_reactions.r1.concentration_form = \
            ConcentrationForm.molality
        model.rparams.config.equilibrium_reactions.r1.parameter_data = {
            "dh_rxn_ref": (2000, pyunits.J/pyunits.kmol),
            "ds_rxn_ref": (1000, pyunits.J/pyunits.kmol/pyunits.K),
            "T_eq_ref": (900, pyunits.degR)}

        gibbs_energy.build_parameters(
            model.rparams.reaction_r1,
            model.rparams.config.equilibrium_reactions["r1"])

        # Check parameter construction
        assert isinstance(model.rparams.reaction_r1.dh_rxn_ref, Var)
        assert model.rparams.reaction_r1.dh_rxn_ref.value == 2

        assert isinstance(model.rparams.reaction_r1.ds_rxn_ref, Var)
        assert model.rparams.reaction_r1.ds_rxn_ref.value == 1

        assert isinstance(model.rparams.reaction_r1.T_eq_ref, Var)
        assert model.rparams.reaction_r1.T_eq_ref.value == 500

        assert_units_equivalent(model.rparams.reaction_r1._keq_units,
                                pyunits.mol/pyunits.kg)

        # Check expression
        rform = gibbs_energy.return_expression(
            model.rxn[1], model.rparams.reaction_r1, "r1", 300*pyunits.K)

        assert str(rform) == (
            'exp(- rparams.reaction_r1.dh_rxn_ref/(kg*m**2/J/s**2*'
            '(8.314462618*(J)/mol/K)*(300*K)) + '
            '1/(kg*m**2/J/s**2*(8.314462618*(J)/mol/K))*'
            'rparams.reaction_r1.ds_rxn_ref)*(mol*kg**-1)')
        assert value(rform) == pytest.approx(1.1269, rel=1e-3)
        assert_units_equivalent(rform, pyunits.mol*pyunits.kg**-1)

    @pytest.mark.unit
    def test_gibbs_energy_partial_pressure(self, model):
        model.rparams.config.equilibrium_reactions.r1.concentration_form = \
            ConcentrationForm.partialPressure
        model.rparams.config.equilibrium_reactions.r1.parameter_data = {
            "dh_rxn_ref": 2,
            "ds_rxn_ref": 1,
            "T_eq_ref": 500}

        gibbs_energy.build_parameters(
            model.rparams.reaction_r1,
            model.rparams.config.equilibrium_reactions["r1"])

        # Check parameter construction
        assert isinstance(model.rparams.reaction_r1.dh_rxn_ref, Var)
        assert model.rparams.reaction_r1.dh_rxn_ref.value == 2

        assert isinstance(model.rparams.reaction_r1.ds_rxn_ref, Var)
        assert model.rparams.reaction_r1.ds_rxn_ref.value == 1

        assert isinstance(model.rparams.reaction_r1.T_eq_ref, Var)
        assert model.rparams.reaction_r1.T_eq_ref.value == 500

        assert_units_equivalent(model.rparams.reaction_r1._keq_units,
                                pyunits.Pa)

        # Check expression
        rform = gibbs_energy.return_expression(
            model.rxn[1], model.rparams.reaction_r1, "r1", 300*pyunits.K)

        assert str(rform) == (
            'exp(- rparams.reaction_r1.dh_rxn_ref/(kg*m**2/J/s**2*'
            '(8.314462618*(J)/mol/K)*(300*K)) + '
            '1/(kg*m**2/J/s**2*(8.314462618*(J)/mol/K))*'
            'rparams.reaction_r1.ds_rxn_ref)*(kg*m**-1*s**-2)')
        assert value(rform) == pytest.approx(1.1269, rel=1e-3)
        assert_units_equivalent(rform, pyunits.Pa)

    @pytest.mark.unit
    def test_gibbs_energy_partial_pressure_convert(self, model):
        model.rparams.config.equilibrium_reactions.r1.concentration_form = \
            ConcentrationForm.partialPressure
        model.rparams.config.equilibrium_reactions.r1.parameter_data = {
            "dh_rxn_ref": (2000, pyunits.J/pyunits.kmol),
            "ds_rxn_ref": (1000, pyunits.J/pyunits.kmol/pyunits.K),
            "T_eq_ref": (900, pyunits.degR)}

        gibbs_energy.build_parameters(
            model.rparams.reaction_r1,
            model.rparams.config.equilibrium_reactions["r1"])

        # Check parameter construction
        assert isinstance(model.rparams.reaction_r1.dh_rxn_ref, Var)
        assert model.rparams.reaction_r1.dh_rxn_ref.value == 2

        assert isinstance(model.rparams.reaction_r1.ds_rxn_ref, Var)
        assert model.rparams.reaction_r1.ds_rxn_ref.value == 1

        assert isinstance(model.rparams.reaction_r1.T_eq_ref, Var)
        assert model.rparams.reaction_r1.T_eq_ref.value == 500

        assert_units_equivalent(model.rparams.reaction_r1._keq_units,
                                pyunits.Pa)

        # Check expression
        rform = gibbs_energy.return_expression(
            model.rxn[1], model.rparams.reaction_r1, "r1", 300*pyunits.K)

        assert str(rform) == (
            'exp(- rparams.reaction_r1.dh_rxn_ref/(kg*m**2/J/s**2*'
            '(8.314462618*(J)/mol/K)*(300*K)) + '
            '1/(kg*m**2/J/s**2*(8.314462618*(J)/mol/K))*'
            'rparams.reaction_r1.ds_rxn_ref)*(kg*m**-1*s**-2)')
        assert value(rform) == pytest.approx(1.1269, rel=1e-3)
        assert_units_equivalent(rform, pyunits.Pa)
