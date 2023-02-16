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
Tests for methods from CoolProp forms

Authors: Andrew Lee
"""

import pytest

from pyomo.environ import ConcreteModel, Block, value, Var, units as pyunits
from pyomo.util.check_units import assert_units_equivalent

import idaes.models.properties.modular_properties.coolprop.coolprop_forms as cforms

from idaes.core.util.exceptions import ConfigurationError


def between(y, x1, x2):
    return 0 > (y - x1) * (y - x2)


class TestBaseForms:
    # Use parameters for O2 for testing purposes
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.params = Block()
        m.state = Block([0])

        return m

    @pytest.mark.unit
    def test_exponential_parameters_length_mismatch(self, model):
        # Using parameters for O2 saturation pressure as test case
        with pytest.raises(
            ConfigurationError,
            match="params mismatched length between n and t "
            "parameters for CoolProp exponential form for "
            "property pressure_sat. Please ensure the number "
            "of n and t parameters are equal.",
        ):
            cforms.parameters_nt_sum(
                model.params,
                "pressure_sat",
                [
                    -7.645535357,
                    2.42142887,
                    9.642620061,
                    -10.09482287,
                    -1.685668959,
                    1.277640761,
                ],
                [1.019, 1.177, 2.44, 2.493, 5.646],
            )  # Missing one term for t

    @pytest.mark.unit
    def test_exponential_parameters(self, model):
        # Using parameters for O2 saturation pressure as test case
        cforms.parameters_nt_sum(
            model.params,
            "pressure_sat",
            [
                -7.645535357,
                2.42142887,
                9.642620061,
                -10.09482287,
                -1.685668959,
                1.277640761,
            ],
            [1.019, 1.177, 2.44, 2.493, 5.646, 10.887],
        )

        assert isinstance(model.params.pressure_sat_coeff_n1, Var)
        assert isinstance(model.params.pressure_sat_coeff_n2, Var)
        assert isinstance(model.params.pressure_sat_coeff_n3, Var)
        assert isinstance(model.params.pressure_sat_coeff_n4, Var)
        assert isinstance(model.params.pressure_sat_coeff_n5, Var)
        assert isinstance(model.params.pressure_sat_coeff_n6, Var)

        assert isinstance(model.params.pressure_sat_coeff_t1, Var)
        assert isinstance(model.params.pressure_sat_coeff_t2, Var)
        assert isinstance(model.params.pressure_sat_coeff_t3, Var)
        assert isinstance(model.params.pressure_sat_coeff_t4, Var)
        assert isinstance(model.params.pressure_sat_coeff_t5, Var)
        assert isinstance(model.params.pressure_sat_coeff_t6, Var)

        assert model.params.pressure_sat_coeff_n1.value == -7.645535357
        assert model.params.pressure_sat_coeff_n1.fixed
        assert model.params.pressure_sat_coeff_n2.value == 2.42142887
        assert model.params.pressure_sat_coeff_n2.fixed
        assert model.params.pressure_sat_coeff_n3.value == 9.642620061
        assert model.params.pressure_sat_coeff_n3.fixed
        assert model.params.pressure_sat_coeff_n4.value == -10.09482287
        assert model.params.pressure_sat_coeff_n4.fixed
        assert model.params.pressure_sat_coeff_n5.value == -1.685668959
        assert model.params.pressure_sat_coeff_n5.fixed
        assert model.params.pressure_sat_coeff_n6.value == 1.277640761
        assert model.params.pressure_sat_coeff_n6.fixed

        assert model.params.pressure_sat_coeff_t1.value == 1.019
        assert model.params.pressure_sat_coeff_t1.fixed
        assert model.params.pressure_sat_coeff_t2.value == 1.177
        assert model.params.pressure_sat_coeff_t2.fixed
        assert model.params.pressure_sat_coeff_t3.value == 2.44
        assert model.params.pressure_sat_coeff_t3.fixed
        assert model.params.pressure_sat_coeff_t4.value == 2.493
        assert model.params.pressure_sat_coeff_t4.fixed
        assert model.params.pressure_sat_coeff_t5.value == 5.646
        assert model.params.pressure_sat_coeff_t5.fixed
        assert model.params.pressure_sat_coeff_t6.value == 10.887
        assert model.params.pressure_sat_coeff_t6.fixed

    @pytest.mark.unit
    def test_exponential_sum(self, model):
        expr = cforms._nt_sum(model.params, "pressure_sat", 42)

        assert str(expr) == str(
            model.params.pressure_sat_coeff_n1
            * 42**model.params.pressure_sat_coeff_t1
            + model.params.pressure_sat_coeff_n2
            * 42**model.params.pressure_sat_coeff_t2
            + model.params.pressure_sat_coeff_n3
            * 42**model.params.pressure_sat_coeff_t3
            + model.params.pressure_sat_coeff_n4
            * 42**model.params.pressure_sat_coeff_t4
            + model.params.pressure_sat_coeff_n5
            * 42**model.params.pressure_sat_coeff_t5
            + model.params.pressure_sat_coeff_n6
            * 42**model.params.pressure_sat_coeff_t6
        )

        assert_units_equivalent(expr, pyunits.dimensionless)

    @pytest.mark.unit
    def test_expression_exponential_tau(self, model):
        model.params.pressure_crit = Var(initialize=5043000, units=pyunits.Pa)
        model.params.temperature_crit = Var(initialize=154.581, units=pyunits.K)

        data = {
            100: 254007.641,
            110: 543399.5533,
            120: 1022269.512,
            130: 1749056.109,
            140: 2787800.883,
            150: 4218667.236,
        }

        for T, Psat in data.items():
            expr = cforms.expression_exponential(
                model.params,
                "pressure_sat",
                T * pyunits.K,
                model.params.pressure_crit,
                tau=True,
            )

            assert value(expr) == pytest.approx(Psat, rel=1e-8)
            assert_units_equivalent(expr, pyunits.Pa)

    @pytest.mark.unit
    def test_expression_exponential_tau_dT(self, model):

        fdiff = 1e-5
        for T in range(10, 101, 10):
            dT = value(
                cforms.dT_expression_exponential(
                    model.params,
                    "pressure_sat",
                    T * pyunits.K,
                    model.params.pressure_crit,
                    tau=True,
                )
            )

            em = cforms.expression_exponential(
                model.params,
                "pressure_sat",
                (T - fdiff) * pyunits.K,
                model.params.pressure_crit,
                tau=True,
            )
            e = cforms.expression_exponential(
                model.params,
                "pressure_sat",
                T * pyunits.K,
                model.params.pressure_crit,
                tau=True,
            )
            ep = cforms.expression_exponential(
                model.params,
                "pressure_sat",
                (T + fdiff) * pyunits.K,
                model.params.pressure_crit,
                tau=True,
            )

            dT_fe_m = value((e - em) / fdiff)
            dT_fe_p = value((ep - e) / fdiff)

            assert between(dT, dT_fe_m, dT_fe_p)
            assert pytest.approx(dT_fe_m, rel=1e-4) == dT
            assert pytest.approx(dT_fe_p, rel=1e-4) == dT

    @pytest.mark.unit
    def test_expression_exponential(self, model):
        # Use O2 parameters just to check results
        data = {
            100: 729628.576,
            110: 1033154.739,
            120: 1460907.321,
            130: 2069823.39,
            140: 2948110.91,
            150: 4241040.16,
        }

        for T, Psat in data.items():
            expr = cforms.expression_exponential(
                model.params, "pressure_sat", T * pyunits.K, model.params.pressure_crit
            )

            assert value(expr) == pytest.approx(Psat, rel=1e-8)
            assert_units_equivalent(expr, pyunits.Pa)

    @pytest.mark.unit
    def test_expression_nonexponential(self, model):
        # Manufacture results using O2 paramters for Psat
        data = {
            100: -4706232.585007142,
            110: -2952092.3626429834,
            120: -1204991.804869883,
            130: 552017.5795317084,
            140: 2335733.299130539,
            150: 4169589.797731975,
        }

        for T, Psat in data.items():
            expr = cforms.expression_nonexponential(
                model.params, "pressure_sat", T * pyunits.K, model.params.pressure_crit
            )

            assert value(expr) == pytest.approx(Psat, rel=1e-8)
            assert_units_equivalent(expr, pyunits.Pa)

    @pytest.mark.unit
    def test_polynomial_parameters(self, model):
        # Using parameters for O2 saturation pressure as test case
        cforms.parameters_polynomial(
            model.params,
            "enth_mol_liq_comp",
            pyunits.J / pyunits.mol,
            [
                -9025.64288846283,
                -35.0797393389397,
                4.78790950986069,
                -8.95948250127742e-2,
                9.20005110933356e-4,
                -5.57355242293296e-06,
                1.85555931770993e-08,
                -2.63480658094238e-11,
            ],
            [1, -6.39773824359397e-3],
        )

        assert isinstance(model.params.enth_mol_liq_comp_coeff_A0, Var)
        assert isinstance(model.params.enth_mol_liq_comp_coeff_A1, Var)
        assert isinstance(model.params.enth_mol_liq_comp_coeff_A2, Var)
        assert isinstance(model.params.enth_mol_liq_comp_coeff_A3, Var)
        assert isinstance(model.params.enth_mol_liq_comp_coeff_A4, Var)
        assert isinstance(model.params.enth_mol_liq_comp_coeff_A5, Var)
        assert isinstance(model.params.enth_mol_liq_comp_coeff_A6, Var)
        assert isinstance(model.params.enth_mol_liq_comp_coeff_A7, Var)

        assert isinstance(model.params.enth_mol_liq_comp_coeff_B0, Var)
        assert isinstance(model.params.enth_mol_liq_comp_coeff_B1, Var)

        assert model.params.enth_mol_liq_comp_coeff_A0.value == -9025.64288846283
        assert model.params.enth_mol_liq_comp_coeff_A0.fixed
        assert model.params.enth_mol_liq_comp_coeff_A1.value == -35.0797393389397
        assert model.params.enth_mol_liq_comp_coeff_A1.fixed
        assert model.params.enth_mol_liq_comp_coeff_A2.value == 4.78790950986069
        assert model.params.enth_mol_liq_comp_coeff_A2.fixed
        assert model.params.enth_mol_liq_comp_coeff_A3.value == -8.95948250127742e-2
        assert model.params.enth_mol_liq_comp_coeff_A3.fixed
        assert model.params.enth_mol_liq_comp_coeff_A4.value == 9.20005110933356e-4
        assert model.params.enth_mol_liq_comp_coeff_A4.fixed
        assert model.params.enth_mol_liq_comp_coeff_A5.value == -5.57355242293296e-06
        assert model.params.enth_mol_liq_comp_coeff_A5.fixed
        assert model.params.enth_mol_liq_comp_coeff_A6.value == 1.85555931770993e-08
        assert model.params.enth_mol_liq_comp_coeff_A6.fixed
        assert model.params.enth_mol_liq_comp_coeff_A7.value == -2.63480658094238e-11
        assert model.params.enth_mol_liq_comp_coeff_A7.fixed

        assert model.params.enth_mol_liq_comp_coeff_B0.value == 1
        assert model.params.enth_mol_liq_comp_coeff_B0.fixed
        assert model.params.enth_mol_liq_comp_coeff_B1.value == -6.39773824359397e-3
        assert model.params.enth_mol_liq_comp_coeff_B1.fixed

    @pytest.mark.unit
    def test_expression_polynomial(self, model):
        # Use O2 parameters just to check results
        data = {
            100: -3726.193393,
            110: -3156.212331,
            120: -2557.035333,
            130: -1909.402314,
            140: -1173.118775,
            150: -217.8380819,
        }
        h_anchor = 2002.355546

        for T, hL in data.items():
            expr = cforms.expression_polynomial(
                model.params, "enth_mol_liq_comp", T * pyunits.K
            )

            assert value(expr) == pytest.approx(hL - h_anchor, rel=1e-8)
            assert_units_equivalent(expr, pyunits.J / pyunits.mol)
