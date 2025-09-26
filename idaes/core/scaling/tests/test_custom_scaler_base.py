#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Tests for CustomScalerBase.

Author: Andrew Lee
"""
import pytest
import re

from pyomo.environ import (
    Block,
    ComponentMap,
    ConcreteModel,
    Constraint,
    Expression,
    Set,
    Suffix,
    units,
    value,
    Var,
)
from pyomo.common.config import ConfigDict

from idaes.core.scaling.custom_scaler_base import (
    CustomScalerBase,
    ConstraintScalingScheme,
    DefaultScalingRecommendation,
)
from idaes.core.util.constants import Constants
from idaes.core.util.testing import PhysicalParameterTestBlock
import idaes.logger as idaeslog


class DummyScaler(CustomScalerBase):
    """
    Create some dummy methods that will record that they ran in lists
    """

    def variable_scaling_routine(self, model, overwrite, submodel_scalers):
        # Create some lists to show that each method was run
        model._verification = []
        model._subscalers = {}
        model.overwrite = []

        model._verification.append("variables")
        model._subscalers["variables"] = submodel_scalers
        model.overwrite.append(overwrite)

    def constraint_scaling_routine(self, model, overwrite, submodel_scalers):
        model._verification.append("constraints")
        model._subscalers["constraints"] = submodel_scalers
        model.overwrite.append(overwrite)

    def fill_in_1(self, model):
        model._verification.append("fill_in_1")

    def fill_in_2(self, model):
        model._verification.append("fill_in_2")

    def dummy_method(self, model, overwrite, submodel_scalers):
        model._dummy_scaler_test = overwrite


@pytest.fixture
def model():
    m = ConcreteModel()

    m.pressure = Var(units=units.Pa, bounds=(1e4, None))
    m.temperature = Var(units=units.K, bounds=(250, 350))
    m.volume_mol = Var(units=units.m**3 / units.mol, bounds=(None, 10))
    m.enth_mol = Var(units=units.J / units.mol)

    m.ideal_gas = Constraint(
        expr=m.pressure * m.volume_mol == Constants.gas_constant * m.temperature
    )
    m.enthalpy_eq = Constraint(
        expr=4.81 * units.J / units.mol / units.K * m.temperature == m.enth_mol
    )
    m.enthalpy_expr = Expression(
        expr=4.81 * units.J / units.mol / units.K * m.temperature
    )

    m.scaling_factor = Suffix(direction=Suffix.EXPORT)

    return m


class TestCustomScalerBase:
    @pytest.mark.unit
    def test_init(self):
        sb = CustomScalerBase()

        assert isinstance(sb.config, ConfigDict)
        assert sb.config.zero_tolerance == 1e-12
        assert sb.config.max_variable_scaling_factor == 1e10
        assert sb.config.min_variable_scaling_factor == 1e-10
        assert sb.config.max_constraint_scaling_factor == 1e10
        assert sb.config.min_constraint_scaling_factor == 1e-10
        assert not sb.config.overwrite

        assert sb.default_scaling_factors == {}
        assert sb.unit_scaling_factors == {
            "Temperature": (units.K, 1e-2),
            "Pressure": (units.Pa, 1e-5),
        }

    @pytest.mark.unit
    def test_variable_scaling_routine(self):
        sb = CustomScalerBase()

        with pytest.raises(
            NotImplementedError,
            match="Custom Scaler has not implemented a variable_scaling_routine method.",
        ):
            sb.variable_scaling_routine("foo")

    @pytest.mark.unit
    def test_constraint_scaling_routine(self):
        sb = CustomScalerBase()

        with pytest.raises(
            NotImplementedError,
            match="Custom Scaler has not implemented a constraint_scaling_routine method.",
        ):
            sb.constraint_scaling_routine("foo")

    @pytest.mark.unit
    def test_scale_model(self):
        # Dummy object to hold testing data
        model = ConcreteModel()

        sb = DummyScaler()
        sb.scale_model(
            model,
            first_stage_fill_in=[sb.fill_in_1],
            second_stage_fill_in=[sb.fill_in_2],
            submodel_scalers={"foo", "bar"},
        )

        # Check order methods were run in
        assert model._verification == [
            "variables",
            "fill_in_1",
            "constraints",
            "fill_in_2",
        ]

        # Check that overwrite was passed on
        model.overwrite == [False, False]

        # Check that submodel scalers were passed on to methods
        assert model._subscalers == {
            "variables": {"foo", "bar"},
            "constraints": {"foo", "bar"},
        }

    @pytest.mark.unit
    def test_get_default_scaling_factor(self, model, caplog):
        caplog.set_level(idaeslog.DEBUG, logger="idaes")
        sb = CustomScalerBase()

        # No defaults defined yet
        assert sb.get_default_scaling_factor(model.pressure) is None
        assert "No default scaling factor found for pressure" in caplog.text

        # Set a default
        sb.default_scaling_factors["pressure"] = 1e-4
        assert sb.get_default_scaling_factor(model.pressure) == 1e-4

    @pytest.mark.unit
    def test_get_default_scaling_factor_indexed(self, caplog):
        caplog.set_level(idaeslog.DEBUG, logger="idaes")
        m = ConcreteModel()
        m.v = Var([1, 2, 3, 4])

        sb = CustomScalerBase()

        # No defaults defined yet
        assert sb.get_default_scaling_factor(m.v[1]) is None
        assert "No default scaling factor found for v[1]" in caplog.text

        # Set a default for the indexed var
        sb.default_scaling_factors["v"] = 1e-4
        assert sb.get_default_scaling_factor(m.v[1]) == 1e-4

        # Set a default for the specific element
        sb.default_scaling_factors["v[1]"] = 1e-8
        assert sb.get_default_scaling_factor(m.v[1]) == 1e-8

    @pytest.mark.unit
    def test_scale_variable_by_component(self, model, caplog):
        caplog.set_level(idaeslog.DEBUG, logger="idaes")
        sb = CustomScalerBase()

        # No scaling factors set
        sb.scale_variable_by_component(model.pressure, model.temperature)
        assert sb.get_scaling_factor(model.pressure) is None
        assert (
            "Could not set scaling factor for pressure, no scaling factor set for temperature"
            in caplog.text
        )

        # Set a scaling factor for temperature
        model.scaling_factor = Suffix(direction=Suffix.EXPORT)
        model.scaling_factor[model.temperature] = 1e-2

        sb.scale_variable_by_component(model.pressure, model.temperature)
        assert model.scaling_factor[model.pressure] == 1e-2

        # Change temperature scaling and check for overwrite
        model.scaling_factor[model.temperature] = 1e-3

        sb.scale_variable_by_component(
            model.pressure, model.temperature, overwrite=False
        )
        assert model.scaling_factor[model.pressure] == 1e-2

        sb.scale_variable_by_component(
            model.pressure, model.temperature, overwrite=True
        )
        assert model.scaling_factor[model.pressure] == 1e-3

    @pytest.mark.unit
    def test_scale_variable_by_bounds(self, model, caplog):
        caplog.set_level(idaeslog.DEBUG, logger="idaes")
        sb = CustomScalerBase()

        # Set a scaling factor for temperature to test overwrite
        model.scaling_factor[model.temperature] = 1e-2

        sb.scale_variable_by_bounds(model.temperature, overwrite=False)
        assert model.scaling_factor[model.temperature] == 1e-2

        # Both bounds
        sb.scale_variable_by_bounds(model.temperature, overwrite=True)
        assert model.scaling_factor[model.temperature] == 1 / 300

        # Lower bound only
        sb.scale_variable_by_bounds(model.pressure)
        assert model.scaling_factor[model.pressure] == 1e-4

        # Upper bound only
        sb.scale_variable_by_bounds(model.volume_mol)
        assert model.scaling_factor[model.volume_mol] == 0.1

        # No bounds
        sb.scale_variable_by_bounds(model.enth_mol)
        assert (
            "No scaling factor set for enth_mol; variable has no bounds." in caplog.text
        )
        assert model.enth_mol not in model.scaling_factor

        # Test case where bounds give a magnitude of 0
        model.pressure.setlb(-1e4)
        model.pressure.setub(1e4)
        sb.scale_variable_by_bounds(model.pressure, overwrite=True)
        assert model.scaling_factor[model.pressure] == 1

    @pytest.mark.unit
    def test_scale_variable_by_default_no_default(self, model):
        sb = CustomScalerBase()

        # No defaults defined yet
        with pytest.raises(
            ValueError, match=re.escape("No default scaling factor set for pressure.")
        ):
            sb.scale_variable_by_default(model.pressure)
        assert model.pressure not in model.scaling_factor

    @pytest.mark.unit
    def test_scale_variable_by_default_constraint(self, model):
        sb = CustomScalerBase()

        # No defaults defined yet
        with pytest.raises(
            TypeError,
            match=re.escape(
                "ideal_gas is type <class 'pyomo.core.base.constraint.ScalarConstraint'>, "
                "but a variable or expression was expected."
            ),
        ):
            sb.scale_variable_by_default(model.ideal_gas)
        assert model.ideal_gas not in model.scaling_factor

    @pytest.mark.unit
    def test_scale_variable_by_default_indexed(self, model):
        model.mole_frac_comp = Var(["N2", "O2"])

        @model.Constraint(["N2", "O2"])
        def mole_frac_eqn(b, j):
            return b.mole_frac_comp[j] == 0.5

        sb = CustomScalerBase()

        with pytest.raises(
            TypeError,
            match=re.escape(
                "mole_frac_comp is indexed. Call with ComponentData children instead."
            ),
        ):
            sb.scale_variable_by_default(model.mole_frac_comp)
        with pytest.raises(
            TypeError,
            match=re.escape(
                "mole_frac_eqn is indexed. Call with ComponentData children instead."
            ),
        ):
            sb.scale_variable_by_default(model.mole_frac_eqn)
        with pytest.raises(
            ValueError,
            match=re.escape("No default scaling factor set for mole_frac_comp[N2]."),
        ):
            sb.scale_variable_by_default(model.mole_frac_comp["N2"])
        with pytest.raises(
            TypeError,
            match=re.escape(
                "mole_frac_eqn[N2] is type <class 'pyomo.core.base.constraint.ConstraintData'>,"
                " but a variable or expression was expected."
            ),
        ):
            sb.scale_variable_by_default(model.mole_frac_eqn["N2"])
        sb.default_scaling_factors["mole_frac_comp[N2]"] = 7
        sb.scale_variable_by_default(model.mole_frac_comp["N2"])
        assert model.scaling_factor[model.mole_frac_comp["N2"]] == 7
        assert model.mole_frac_comp["O2"] not in model.scaling_factor

    @pytest.mark.unit
    def test_scale_variable_by_default_user_input_required(self, model):
        sb = CustomScalerBase()
        sb.default_scaling_factors["pressure"] = (
            DefaultScalingRecommendation.userInputRequired
        )
        # No defaults defined yet
        with pytest.raises(
            ValueError, match=re.escape("No default scaling factor set for pressure.")
        ):
            sb.scale_variable_by_default(model.pressure)
        assert model.pressure not in model.scaling_factor

        # If a scaling factor is already set, then no exception is raised
        sb.set_component_scaling_factor(model.pressure, 1e-4)
        sb.scale_variable_by_default(model.pressure)
        assert model.scaling_factor[model.pressure] == 1e-4

        # If we tell it to overwrite the scaling factors, the existence of
        # a preexisting scaling factor is no longer sufficient.
        with pytest.raises(
            ValueError, match=re.escape("No default scaling factor set for pressure.")
        ):
            sb.scale_variable_by_default(model.pressure, overwrite=True)
        assert model.scaling_factor[model.pressure] == 1e-4

        # If user certifies that they set the scaling factor manually,
        # then overwrite doesn't raise an exception
        sb.default_scaling_factors["pressure"] = (
            DefaultScalingRecommendation.userSetManually
        )
        sb.scale_variable_by_default(model.pressure, overwrite=True)
        assert model.scaling_factor[model.pressure] == 1e-4

    @pytest.mark.unit
    def test_scale_variable_by_default_user_input_recommended(self, model):
        sb = CustomScalerBase()
        sb.default_scaling_factors["pressure"] = (
            DefaultScalingRecommendation.userInputRecommended
        )
        # The scaling method is going to generate a guess for the scaling
        # factor later, so no guess is necessary now.
        sb.scale_variable_by_default(model.pressure)
        assert model.pressure not in model.scaling_factor

    @pytest.mark.unit
    def test_scale_variable_by_default(self, model, caplog):
        caplog.set_level(idaeslog.DEBUG, logger="idaes")
        sb = CustomScalerBase()

        # Set a default
        sb.default_scaling_factors["pressure"] = 1e-4
        sb.scale_variable_by_default(model.pressure)
        assert model.scaling_factor[model.pressure] == 1e-4

        # Change default to check overwrite
        sb.default_scaling_factors["pressure"] = 1e-5
        sb.scale_variable_by_default(model.pressure, overwrite=False)
        assert model.scaling_factor[model.pressure] == 1e-4
        sb.scale_variable_by_default(model.pressure, overwrite=True)
        assert model.scaling_factor[model.pressure] == 1e-5
        sb.default_scaling_factors["enthalpy_expr"] = 1e-3

        sb.scale_variable_by_default(model.enthalpy_expr)
        assert model.scaling_hint[model.enthalpy_expr] == 1e-3

    @pytest.mark.unit
    def test_scale_variable_by_units(self, model, caplog):
        caplog.set_level(idaeslog.DEBUG, logger="idaes")
        sb = CustomScalerBase()

        # Units in dict, no conversion required
        sb.scale_variable_by_units(model.pressure)
        assert model.scaling_factor[model.pressure] == 1e-5

        # Equivalent units in dict, need conversion
        model.p2 = Var(units=units.bar)
        sb.scale_variable_by_units(model.p2)
        assert model.scaling_factor[model.p2] == pytest.approx(1, rel=1e-8)

        # No match - no sf assigned
        sb.scale_variable_by_units(model.enth_mol)
        assert model.enth_mol not in model.scaling_factor
        assert (
            "No scaling factor set for enth_mol; no match for units J/mol found "
            "in self.unit_scaling_factors" in caplog.text
        )

        # Test for overwrite
        model.scaling_factor[model.temperature] = 42
        sb.scale_variable_by_units(model.temperature, overwrite=False)
        assert model.scaling_factor[model.temperature] == 42
        sb.scale_variable_by_units(model.temperature, overwrite=True)
        assert model.scaling_factor[model.temperature] == 1e-2

    @pytest.mark.unit
    def test_scale_variable_by_definition_constraint(self, model):
        sb = CustomScalerBase()
        model.temperature.value = 7
        sb.set_variable_scaling_factor(model.temperature, 1 / 300)
        sb.scale_variable_by_definition_constraint(
            model.enth_mol,
            model.enthalpy_eq,
        )
        assert model.scaling_factor[model.enth_mol] == pytest.approx(
            1 / (4.81 * 300), rel=1e-15
        )
        assert model.temperature.value == 7

    @pytest.mark.unit
    def test_scale_variable_by_definition_constraint_without_variable(self, model):
        sb = CustomScalerBase()
        with pytest.raises(
            ValueError,
            match=re.escape(
                "Variable pressure does not appear in constraint "
                "enthalpy_eq, cannot calculate scaling factor."
            ),
        ):
            sb.scale_variable_by_definition_constraint(
                model.pressure,
                model.enthalpy_eq,
            )

    @pytest.mark.unit
    def test_scale_variable_by_definition_constraint_not_constraint(self, model):
        sb = CustomScalerBase()
        with pytest.raises(
            TypeError,
            match=re.escape(
                "enthalpy_expr is not a constraint, but instead "
                "<class 'pyomo.core.base.expression.ScalarExpression'>"
            ),
        ):
            sb.scale_variable_by_definition_constraint(
                model.pressure,
                model.enthalpy_expr,
            )

    @pytest.mark.unit
    def test_scale_variable_by_definition_constraint_indexed(self, model):
        sb = CustomScalerBase()

        # The fact that this constraint overdetermines the model is
        # of no consequence. We don't even reach the computation
        # stage by the time an exception is thrown.
        @model.Constraint([1, 2, 3])
        def foo(b, idx):
            return b.pressure == 0

        with pytest.raises(
            TypeError,
            match=re.escape(
                f"Constraint foo is indexed. Call with ConstraintData "
                "children instead."
            ),
        ):
            sb.scale_variable_by_definition_constraint(
                model.pressure,
                model.foo,
            )

    @pytest.mark.unit
    def test_scale_variable_by_definition_constraint_zero_derivative(self, model):
        sb = CustomScalerBase()
        model.foo = Constraint(expr=1 == 0 * model.enth_mol)
        model.enth_mol.value = 42
        with pytest.raises(
            RuntimeError,
            match=re.escape(
                "Could not calculate scaling factor from definition constraint foo. "
                "Does enth_mol appear nonlinearly in it or have a linear coefficient "
                "equal to zero?"
            ),
        ):
            sb.scale_variable_by_definition_constraint(
                model.enth_mol,
                model.foo,
            )
        assert model.enth_mol.value == 42

    @pytest.mark.unit
    def test_scale_variable_by_definition_constraint_zero_derivative(self, model):
        sb = CustomScalerBase()
        model.enth_mol.value = 42
        model.foo = Constraint(expr=0 == model.enth_mol)
        with pytest.raises(
            ValueError,
            match=re.escape(
                "Calculated nominal value of zero from definition constraint."
            ),
        ):
            sb.scale_variable_by_definition_constraint(
                model.enth_mol,
                model.foo,
            )
        assert model.enth_mol.value == 42

    @pytest.mark.unit
    def test_scale_constraint_by_default_no_default(self, model):
        sb = CustomScalerBase()
        with pytest.raises(
            TypeError,
            match=re.escape(
                "pressure is type <class 'pyomo.core.base.var.ScalarVar'>, "
                "but a constraint was expected."
            ),
        ):
            sb.scale_constraint_by_default(model.pressure)
        assert model.pressure not in model.scaling_factor

    @pytest.mark.unit
    def test_scale_constraint_by_default(self, model):
        sb = CustomScalerBase()

        # Set a default
        sb.default_scaling_factors["ideal_gas"] = 1e-3
        sb.scale_constraint_by_default(model.ideal_gas)
        assert model.scaling_factor[model.ideal_gas] == 1e-3

        # Change default to check overwrite
        sb.default_scaling_factors["ideal_gas"] = 1e-5
        sb.scale_constraint_by_default(model.ideal_gas, overwrite=False)
        assert model.scaling_factor[model.ideal_gas] == 1e-3
        sb.scale_constraint_by_default(model.ideal_gas, overwrite=True)
        assert model.scaling_factor[model.ideal_gas] == 1e-5

    @pytest.mark.unit
    def test_scale_constraint_by_default_indexed(self, model):

        model.mole_frac_comp = Var(["N2", "O2"])

        @model.Constraint(["N2", "O2"])
        def mole_frac_eqn(b, j):
            return b.mole_frac_comp[j] == 0.5

        sb = CustomScalerBase()

        with pytest.raises(
            TypeError,
            match=re.escape(
                "mole_frac_comp is indexed. Call with ComponentData children instead."
            ),
        ):
            sb.scale_constraint_by_default(model.mole_frac_comp)
        with pytest.raises(
            TypeError,
            match=re.escape(
                "mole_frac_eqn is indexed. Call with ComponentData children instead."
            ),
        ):
            sb.scale_constraint_by_default(model.mole_frac_eqn)
        with pytest.raises(
            TypeError,
            match=re.escape(
                "mole_frac_comp[N2] is type <class 'pyomo.core.base.var.VarData'>,"
                " but a constraint was expected."
            ),
        ):
            sb.scale_constraint_by_default(model.mole_frac_comp["N2"])
        with pytest.raises(
            ValueError,
            match=re.escape("No default scaling factor set for mole_frac_eqn[N2]."),
        ):
            sb.scale_constraint_by_default(model.mole_frac_eqn["N2"])
        sb.default_scaling_factors["mole_frac_eqn[N2]"] = 7
        sb.scale_constraint_by_default(model.mole_frac_eqn["N2"])
        assert model.scaling_factor[model.mole_frac_eqn["N2"]] == 7
        assert model.mole_frac_eqn["O2"] not in model.scaling_factor

    @pytest.mark.unit
    def test_get_expression_nominal_values(self, model):
        sb = CustomScalerBase()

        # Set variable scaling factors for testing
        model.scaling_factor[model.pressure] = 1e-5
        model.scaling_factor[model.temperature] = 1e-2
        model.scaling_factor[model.volume_mol] = 1e-1

        nominal_value = sb.get_expression_nominal_value(model.ideal_gas.body)

        # Nominal value will be P*V - R*T
        assert nominal_value == pytest.approx(831.446 - 1e6, rel=1e-5)

        # Check redirection for ConstraintData objects
        nominal_value = sb.get_expression_nominal_value(model.ideal_gas)
        assert nominal_value == pytest.approx(831.446 - 1e6, rel=1e-5)

    @pytest.mark.unit
    def test_get_sum_terms_nominal_value(self, model):
        sb = CustomScalerBase()

        # Set variable scaling factors for testing
        model.scaling_factor[model.pressure] = 1e-5
        model.scaling_factor[model.temperature] = 1e-2
        model.scaling_factor[model.volume_mol] = 1e-1

        nominal_values = sb.get_expression_nominal_values(model.ideal_gas.expr)

        # Nominal values will be (R*T, P*V)
        assert nominal_values == [
            pytest.approx(831.446, rel=1e-5),
            pytest.approx(1e6, rel=1e-5),
        ]

        # Check redirection for ConstraintData objects
        nominal_values = sb.get_sum_terms_nominal_values(model.ideal_gas)
        assert nominal_values == [
            pytest.approx(831.446, rel=1e-5),
            pytest.approx(1e6, rel=1e-5),
        ]

    @pytest.mark.unit
    def test_get_sum_terms_nominal_values(self, model):
        sb = CustomScalerBase()

        # Set variable scaling factors for testing
        model.scaling_factor[model.pressure] = 1e-5
        model.scaling_factor[model.temperature] = 1e-2
        model.scaling_factor[model.volume_mol] = 1e-1

        nominal_values = sb.get_sum_terms_nominal_values(model.ideal_gas.expr)

        # Nominal values will be (R*T, P*V)
        assert nominal_values == [
            pytest.approx(831.446, rel=1e-5),
            pytest.approx(1e6, rel=1e-5),
        ]

        # Check redirection for ConstraintData objects
        nominal_values = sb.get_sum_terms_nominal_values(model.ideal_gas)
        assert nominal_values == [
            pytest.approx(831.446, rel=1e-5),
            pytest.approx(1e6, rel=1e-5),
        ]

    @pytest.mark.unit
    def test_scale_constraint_by_nominal_value_harmonic(self, model):
        sb = CustomScalerBase()

        # Set variable scaling factors for testing
        model.scaling_factor[model.pressure] = 1e-5
        model.scaling_factor[model.temperature] = 1e-2
        model.scaling_factor[model.volume_mol] = 1e-1
        model.scaling_factor[model.ideal_gas] = 1

        # overwrite = False, no change
        sb.scale_constraint_by_nominal_value(
            model.ideal_gas,
            scheme=ConstraintScalingScheme.harmonicMean,
            overwrite=False,
        )
        assert model.scaling_factor[model.ideal_gas] == 1
        # overwrite = True
        sb.scale_constraint_by_nominal_value(
            model.ideal_gas, scheme=ConstraintScalingScheme.harmonicMean, overwrite=True
        )
        assert model.scaling_factor[model.ideal_gas] == pytest.approx(
            (1 / 831.446 + 1e-6), rel=1e-5
        )

    @pytest.mark.unit
    def test_scale_constraint_by_nominal_value_max(self, model):
        sb = CustomScalerBase()

        # Set variable scaling factors for testing
        model.scaling_factor[model.pressure] = 1e-5
        model.scaling_factor[model.temperature] = 1e-2
        model.scaling_factor[model.volume_mol] = 1e-1
        model.scaling_factor[model.ideal_gas] = 1

        # overwrite = False, no change
        sb.scale_constraint_by_nominal_value(model.ideal_gas, overwrite=False)
        assert model.scaling_factor[model.ideal_gas] == 1
        # overwrite = True
        sb.scale_constraint_by_nominal_value(model.ideal_gas, overwrite=True)
        assert model.scaling_factor[model.ideal_gas] == pytest.approx(1e-6, rel=1e-5)

    @pytest.mark.unit
    def test_scale_constraint_by_nominal_value_min(self, model):
        sb = CustomScalerBase()

        # Set variable scaling factors for testing
        model.scaling_factor[model.pressure] = 1e-5
        model.scaling_factor[model.temperature] = 1e-2
        model.scaling_factor[model.volume_mol] = 1e-1
        model.scaling_factor[model.ideal_gas] = 1

        # overwrite = False, no change
        sb.scale_constraint_by_nominal_value(
            model.ideal_gas,
            scheme=ConstraintScalingScheme.inverseMinimum,
            overwrite=False,
        )
        assert model.scaling_factor[model.ideal_gas] == 1
        # overwrite = True
        sb.scale_constraint_by_nominal_value(
            model.ideal_gas,
            scheme=ConstraintScalingScheme.inverseMinimum,
            overwrite=True,
        )
        assert model.scaling_factor[model.ideal_gas] == pytest.approx(
            1 / 831.446, rel=1e-5
        )

    @pytest.mark.unit
    def test_scale_constraint_by_nominal_value_inv_sum(self, model):
        sb = CustomScalerBase()

        # Set variable scaling factors for testing
        model.scaling_factor[model.pressure] = 1e-5
        model.scaling_factor[model.temperature] = 1e-2
        model.scaling_factor[model.volume_mol] = 1e-1
        model.scaling_factor[model.ideal_gas] = 1

        # overwrite = False, no change
        sb.scale_constraint_by_nominal_value(
            model.ideal_gas, scheme=ConstraintScalingScheme.inverseSum, overwrite=False
        )
        assert model.scaling_factor[model.ideal_gas] == 1
        # overwrite = True
        sb.scale_constraint_by_nominal_value(
            model.ideal_gas, scheme=ConstraintScalingScheme.inverseSum, overwrite=True
        )
        assert model.scaling_factor[model.ideal_gas] == pytest.approx(
            1 / (831.446 + 1e6), rel=1e-5
        )

    @pytest.mark.unit
    def test_scale_constraint_by_nominal_value_rss(self, model):
        sb = CustomScalerBase()

        # Set variable scaling factors for testing
        model.scaling_factor[model.pressure] = 1e-5
        model.scaling_factor[model.temperature] = 1e-2
        model.scaling_factor[model.volume_mol] = 1e-1
        model.scaling_factor[model.ideal_gas] = 1

        # overwrite = False, no change
        sb.scale_constraint_by_nominal_value(
            model.ideal_gas, scheme=ConstraintScalingScheme.inverseRSS, overwrite=False
        )
        assert model.scaling_factor[model.ideal_gas] == 1
        # overwrite = True
        sb.scale_constraint_by_nominal_value(
            model.ideal_gas, scheme=ConstraintScalingScheme.inverseRSS, overwrite=True
        )
        assert model.scaling_factor[model.ideal_gas] == pytest.approx(
            1 / (831.446**2 + 1e6**2) ** 0.5, rel=1e-5
        )

    @pytest.mark.unit
    def test_scale_constraint_by_nominal_value_invalid_scheme(self, model):
        sb = CustomScalerBase()

        with pytest.raises(
            ValueError,
            match=re.escape(
                "Invalid value for 'scheme' argument (foo) in "
                "scale_constraint_by_nominal_value."
            ),
        ):
            sb.scale_constraint_by_nominal_value(model.ideal_gas, scheme="foo")

    @pytest.mark.unit
    def test_scale_constraint_by_nominal_derivative_1norm(self, model):
        sb = CustomScalerBase()

        # Set variable scaling factors for testing
        model.scaling_factor[model.pressure] = 1e-5
        model.scaling_factor[model.temperature] = 1e-2
        model.scaling_factor[model.volume_mol] = 1e-1
        model.scaling_factor[model.ideal_gas] = 1

        # overwrite = False, no change
        sb.scale_constraint_by_nominal_derivative_norm(
            model.ideal_gas, norm=1, overwrite=False
        )
        assert model.scaling_factor[model.ideal_gas] == 1
        # overwrite = True
        sb.scale_constraint_by_nominal_derivative_norm(
            model.ideal_gas, norm=1, overwrite=True
        )
        assert model.scaling_factor[model.ideal_gas] == pytest.approx(
            4.99792e-7, rel=1e-5  # (1/(8.314+1e6+1e6)
        )

        # Check for clean up
        assert model.pressure.value is None
        assert model.temperature.value is None
        assert model.pressure.value is None

    @pytest.mark.unit
    def test_scale_constraint_by_nominal_derivative_2norm(self, model):
        sb = CustomScalerBase()

        # Set variable scaling factors for testing
        model.scaling_factor[model.pressure] = 1e-5
        model.scaling_factor[model.temperature] = 1e-2
        model.scaling_factor[model.volume_mol] = 1e-1
        model.scaling_factor[model.ideal_gas] = 1

        # overwrite = False, no change
        sb.scale_constraint_by_nominal_derivative_norm(model.ideal_gas, overwrite=False)
        assert model.scaling_factor[model.ideal_gas] == 1
        # overwrite = True
        sb.scale_constraint_by_nominal_derivative_norm(model.ideal_gas, overwrite=True)
        assert model.scaling_factor[model.ideal_gas] == pytest.approx(
            1 / (8.314**2 + 1e6**2 + 1e6**2) ** 0.5, rel=1e-5
        )

        # Check for clean up
        assert model.pressure.value is None
        assert model.temperature.value is None
        assert model.pressure.value is None

    @pytest.mark.unit
    def test_propagate_state_scaling_indexed_to_indexed(self):
        # Dummy up two state blocks
        m = ConcreteModel()

        m.properties = PhysicalParameterTestBlock()

        m.state1 = m.properties.build_state_block([1, 2, 3])
        m.state2 = m.properties.build_state_block([1, 2, 3])

        # Set scaling factors on state1
        for t, sd in m.state1.items():
            sd.scaling_factor = Suffix(direction=Suffix.EXPORT)
            sd.scaling_factor[sd.temperature] = 100 * t
            sd.scaling_factor[sd.pressure] = 1e5 * t

            count = 1
            for j in sd.flow_mol_phase_comp.values():
                sd.scaling_factor[j] = 10 * t * count
                count += 1

        sb = CustomScalerBase()
        sb.propagate_state_scaling(m.state2, m.state1)

        for t, sd in m.state2.items():
            assert sd.scaling_factor[sd.temperature] == 100 * t
            assert sd.scaling_factor[sd.pressure] == 1e5 * t

            count = 1
            for j in sd.flow_mol_phase_comp.values():
                assert sd.scaling_factor[j] == 10 * t * count
                count += 1

        # Test for overwrite
        for t, sd in m.state1.items():
            sd.scaling_factor[sd.temperature] = 200 * t
            sd.scaling_factor[sd.pressure] = 2e5 * t

            count = 1
            for j in sd.flow_mol_phase_comp.values():
                sd.scaling_factor[j] = 20 * t * count
                count += 1

        sb.propagate_state_scaling(m.state2, m.state1, overwrite=False)

        for t, sd in m.state2.items():
            assert sd.scaling_factor[sd.temperature] == 100 * t
            assert sd.scaling_factor[sd.pressure] == 1e5 * t

            count = 1
            for j in sd.flow_mol_phase_comp.values():
                assert sd.scaling_factor[j] == 10 * t * count
                count += 1

    @pytest.mark.unit
    def test_propagate_state_scaling_scalar_to_indexed(self):
        # Dummy up two state blocks
        m = ConcreteModel()

        m.properties = PhysicalParameterTestBlock()

        m.state1 = m.properties.build_state_block([1, 2, 3])
        m.state2 = m.properties.build_state_block([1, 2, 3])

        # Set scaling factors on state1
        for t, sd in m.state1.items():
            sd.scaling_factor = Suffix(direction=Suffix.EXPORT)
            sd.scaling_factor[sd.temperature] = 100 * t
            sd.scaling_factor[sd.pressure] = 1e5 * t

            count = 1
            for j in sd.flow_mol_phase_comp.values():
                sd.scaling_factor[j] = 10 * t * count
                count += 1

        sb = CustomScalerBase()
        sb.propagate_state_scaling(m.state2, m.state1[1])

        for t, sd in m.state2.items():
            assert sd.scaling_factor[sd.temperature] == 100 * 1
            assert sd.scaling_factor[sd.pressure] == 1e5 * 1

            count = 1
            for j in sd.flow_mol_phase_comp.values():
                assert sd.scaling_factor[j] == 10 * 1 * count
                count += 1

        # Test for overwrite=False
        for t, sd in m.state1.items():
            sd.scaling_factor[sd.temperature] = 200 * t
            sd.scaling_factor[sd.pressure] = 2e5 * t

            count = 1
            for j in sd.flow_mol_phase_comp.values():
                sd.scaling_factor[j] = 20 * t * count
                count += 1

        sb.propagate_state_scaling(m.state2, m.state1[1], overwrite=False)

        for t, sd in m.state2.items():
            assert sd.scaling_factor[sd.temperature] == 100 * 1
            assert sd.scaling_factor[sd.pressure] == 1e5 * 1

            count = 1
            for j in sd.flow_mol_phase_comp.values():
                assert sd.scaling_factor[j] == 10 * 1 * count
                count += 1

        # Test for overwrite=True
        sb.propagate_state_scaling(m.state2, m.state1[2], overwrite=True)

        for t, sd in m.state2.items():
            assert sd.scaling_factor[sd.temperature] == 200 * 2
            assert sd.scaling_factor[sd.pressure] == 2e5 * 2

            count = 1
            for j in sd.flow_mol_phase_comp.values():
                assert sd.scaling_factor[j] == 20 * 2 * count
                count += 1

    @pytest.mark.unit
    def test_propagate_state_scaling_scalar_to_indexed_within_block(self):
        # Dummy up two state blocks
        m = ConcreteModel()

        m.properties = PhysicalParameterTestBlock()

        m.state1 = m.properties.build_state_block([1, 2, 3])

        # Set scaling factors on state1
        for t, sd in m.state1.items():
            sd.scaling_factor = Suffix(direction=Suffix.EXPORT)
            sd.scaling_factor[sd.temperature] = 100 * t
            sd.scaling_factor[sd.pressure] = 1e5 * t

            count = 1
            for j in sd.flow_mol_phase_comp.values():
                sd.scaling_factor[j] = 10 * t * count
                count += 1

        sb = CustomScalerBase()
        sb.propagate_state_scaling(m.state1, m.state1[1])

        for t, sd in m.state1.items():
            assert sd.scaling_factor[sd.temperature] == 100 * t
            assert sd.scaling_factor[sd.pressure] == 1e5 * t

            count = 1
            for j in sd.flow_mol_phase_comp.values():
                assert sd.scaling_factor[j] == 10 * t * count
                count += 1

        # Test for overwrite=False
        for t, sd in m.state1.items():
            sd.scaling_factor[sd.temperature] = 200 * t
            sd.scaling_factor[sd.pressure] = 2e5 * t

            count = 1
            for j in sd.flow_mol_phase_comp.values():
                sd.scaling_factor[j] = 20 * t * count
                count += 1

        sb.propagate_state_scaling(m.state1, m.state1[1], overwrite=False)

        for t, sd in m.state1.items():
            assert sd.scaling_factor[sd.temperature] == 200 * t
            assert sd.scaling_factor[sd.pressure] == 2e5 * t

            count = 1
            for j in sd.flow_mol_phase_comp.values():
                assert sd.scaling_factor[j] == 20 * t * count
                count += 1

        # Test for overwrite=True
        sb.propagate_state_scaling(m.state1, m.state1[2], overwrite=True)

        for t, sd in m.state1.items():
            assert sd.scaling_factor[sd.temperature] == 200 * 2
            assert sd.scaling_factor[sd.pressure] == 2e5 * 2

            count = 1
            for j in sd.flow_mol_phase_comp.values():
                assert sd.scaling_factor[j] == 20 * 2 * count
                count += 1

    @pytest.mark.unit
    def test_propagate_state_scaling_scalar_to_scalar(self):
        # Dummy up two state blocks
        m = ConcreteModel()

        m.properties = PhysicalParameterTestBlock()

        m.state1 = m.properties.build_state_block([1, 2, 3])
        m.state2 = m.properties.build_state_block([1, 2, 3])

        # Set scaling factors on state1
        for t, sd in m.state1.items():
            sd.scaling_factor = Suffix(direction=Suffix.EXPORT)
            sd.scaling_factor[sd.temperature] = 100 * t
            sd.scaling_factor[sd.pressure] = 1e5 * t

            count = 1
            for j in sd.flow_mol_phase_comp.values():
                sd.scaling_factor[j] = 10 * t * count
                count += 1

        sb = CustomScalerBase()
        sb.propagate_state_scaling(m.state2[2], m.state1[3])

        for t, sd in m.state2.items():
            if t == 2:
                assert sd.scaling_factor[sd.temperature] == 100 * 3
                assert sd.scaling_factor[sd.pressure] == 1e5 * 3

                count = 1
                for j in sd.flow_mol_phase_comp.values():
                    assert sd.scaling_factor[j] == 10 * 3 * count
                    count += 1
            else:
                assert not hasattr(sd, "scaling_factor")

    @pytest.mark.unit
    def test_propagate_state_scaling_indexed_to_scalar(self):
        # Dummy up two state blocks
        m = ConcreteModel()

        m.properties = PhysicalParameterTestBlock()

        m.state1 = m.properties.build_state_block([1, 2, 3])
        m.state2 = m.properties.build_state_block([1, 2, 3])

        # Set scaling factors on state1
        for t, sd in m.state1.items():
            sd.scaling_factor = Suffix(direction=Suffix.EXPORT)
            sd.scaling_factor[sd.temperature] = 100 * t
            sd.scaling_factor[sd.pressure] = 1e5 * t

            count = 1
            for j in sd.flow_mol_phase_comp.values():
                sd.scaling_factor[j] = 10 * t * count
                count += 1

        sb = CustomScalerBase()
        with pytest.raises(
            ValueError,
            match=re.escape(
                "Source state block is indexed but target state block is not indexed. "
                "It is ambiguous which index should be used."
            ),
        ):
            sb.propagate_state_scaling(m.state2[2], m.state1)

    @pytest.mark.unit
    def test_call_submodel_scaler_method_no_scaler(self, caplog):
        caplog.set_level(idaeslog.DEBUG, logger="idaes")

        # Dummy up a nested model
        m = ConcreteModel()
        m.b = Block([1, 2, 3])

        sb = CustomScalerBase()
        sb.call_submodel_scaler_method(m.b, method="dummy_method", overwrite=True)

        for bd in m.b.values():
            assert not hasattr(bd, "_dummy_scaler_test")

        assert "No default Scaler set for b. Cannot call dummy_method." in caplog.text

    @pytest.mark.unit
    def test_call_submodel_scaler_method_default_scaler(self, caplog):
        caplog.set_level(idaeslog.DEBUG, logger="idaes")

        # Dummy up a nested model
        m = ConcreteModel()
        m.b = Block([1, 2, 3])
        for bd in m.b.values():
            bd.default_scaler = DummyScaler

        sb = CustomScalerBase()
        sb.call_submodel_scaler_method(m.b, method="dummy_method", overwrite=True)

        for bd in m.b.values():
            assert bd._dummy_scaler_test

        assert "Using default Scaler for b." in caplog.text

    @pytest.mark.unit
    def test_call_submodel_scaler_method_user_scaler(self, caplog):
        caplog.set_level(idaeslog.DEBUG, logger="idaes")

        # Dummy up a nested model
        m = ConcreteModel()
        m.b = Block([1, 2, 3])

        scaler_map = ComponentMap()
        scaler_map[m.b] = DummyScaler()

        sb = CustomScalerBase()
        sb.call_submodel_scaler_method(
            m.b,
            method="dummy_method",
            submodel_scalers=scaler_map,
            overwrite=False,
        )

        for bd in m.b.values():
            assert not bd._dummy_scaler_test

        assert "Using user-defined Scaler for b." in caplog.text

    @pytest.mark.unit
    def test_call_submodel_scaler_method_user_scaler_class(self, caplog):
        caplog.set_level(idaeslog.DEBUG, logger="idaes")

        # Dummy up a nested model
        m = ConcreteModel()
        m.b = Block([1, 2, 3])

        scaler_map = ComponentMap()
        scaler_map[m.b] = DummyScaler()

        sb = CustomScalerBase()
        sb.call_submodel_scaler_method(
            m.b,
            method="dummy_method",
            submodel_scalers=scaler_map,
            overwrite=False,
        )

        for bd in m.b.values():
            assert not bd._dummy_scaler_test

        assert "Using user-defined Scaler for b." in caplog.text

    @pytest.mark.unit
    def test_call_submodel_scaler_method_default_scaler_blockdata(self, caplog):
        caplog.set_level(idaeslog.DEBUG, logger="idaes")

        # Dummy up a nested model
        m = ConcreteModel()
        m.b = Block([1, 2, 3])
        for bd in m.b.values():
            bd.default_scaler = DummyScaler

        sb = CustomScalerBase()
        sb.call_submodel_scaler_method(m.b[1], method="dummy_method", overwrite=True)

        for i, bd in m.b.items():
            if i == 1:
                assert bd._dummy_scaler_test
            else:
                assert not hasattr(bd, "_dummy_scaler_test")

        assert "Using default Scaler for b[1]." in caplog.text


# TODO additional tests for nested submodel scalers.
