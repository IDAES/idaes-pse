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

Author: Andrew Lee, Douglas Allan
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
from idaes.core import (
    declare_process_block_class,
    PhysicalParameterBlock,
    Phase,
    Component,
    StateBlock,
    StateBlockData,
)
from idaes.core.util.constants import Constants
from idaes.core.util.testing import PhysicalParameterTestBlock

# Dummy parameter block to test create-on-demand properties
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


# --------------------------------------------------------------------
# Dummy parameter and state blocks to test interactions with
# build-on-demand properties
# --------------------------------------------------------------------


@declare_process_block_class("Parameters")
class _Parameters(PhysicalParameterBlock):
    def build(self):
        super(_Parameters, self).build()

        self.p1 = Phase()
        self.c1 = Component()

    @classmethod
    def define_metadata(cls, obj):
        obj.define_custom_properties(
            {
                "a": {"method": "a_method"},
                "b": {"method": "b_method"},
                "recursion1": {"method": "_recursion1"},
                "recursion2": {"method": "_recursion2"},
                "not_callable": {"method": "test_obj"},
                "raise_exception": {"method": "_raise_exception"},
                "not_supported": {"supported": False},
                "does_not_create_component": {"method": "_does_not_create_component"},
            }
        )
        obj.add_default_units(
            {
                "time": units.s,
                "mass": units.kg,
                "length": units.m,
                "amount": units.mol,
                "temperature": units.K,
            }
        )


@declare_process_block_class("State", block_class=StateBlock)
class _State(StateBlockData):
    def build(self):
        super(StateBlockData, self).build()

        self.test_obj = 1

    def a_method(self):
        self.a = Var(initialize=1)

    def b_method(self):
        self.b = Var(initialize=1)

    def _recursion1(self):
        self.recursive_cons1 = Constraint(expr=self.recursion2 == 1)

    def _recursion2(self):
        self.recursive_cons2 = Constraint(expr=self.recursion1 == 1)

    def _raise_exception(self):
        # PYLINT-TODO
        # pylint: disable-next=broad-exception-raised
        raise Exception()

    def _does_not_create_component(self):
        pass


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
        m.w = Var(["a", "b", "c"], ["x", "y", "z"])

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

        caplog.clear()
        assert sb.get_default_scaling_factor(m.w["c", "y"]) is None
        assert "No default scaling factor found for w[c,y]" in caplog.text

        sb.default_scaling_factors["w[c,y]"] = 17
        assert sb.get_default_scaling_factor(m.w["c", "y"]) == 17

        caplog.clear()
        assert sb.get_default_scaling_factor(m.w["a", "x"]) is None
        assert "No default scaling factor found for w[a,x]" in caplog.text

        # Make sure that entries with spaces between indices are still found
        sb.default_scaling_factors["w[a, x]"] = 23
        assert sb.get_default_scaling_factor(m.w["a", "x"]) == 23

    @pytest.mark.unit
    def test_get_default_scaling_factor_indexed(self):
        m = ConcreteModel()
        m.params = Parameters()
        m.state = State(parameters=m.params)

        # Trigger build-on-demand creation
        m.state.a

        sb = CustomScalerBase()
        sb.default_scaling_factors["a"] = 7
        sb.default_scaling_factors["b"] = 11
        m.state._lock_attribute_creation = True

        assert sb.get_default_scaling_factor(m.state.a) == 7
        # Want to make sure that _lock_attribute_creation doesn't get
        # changed to false upon exiting get_default_scaling_factor
        assert m.state._lock_attribute_creation == True
        assert not m.state.is_property_constructed("b")

        m.state._lock_attribute_creation = False
        assert sb.get_default_scaling_factor(m.state.a) == 7
        assert m.state._lock_attribute_creation == False
        assert not m.state.is_property_constructed("b")

        # Now trigger creation of b
        m.state.b
        assert sb.get_default_scaling_factor(m.state.b) == 11
        assert m.state._lock_attribute_creation == False

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
            ValueError,
            match=re.escape(
                "This scaler requires the user to provide a default "
                "scaling factor for pressure, but no default scaling "
                "factor was set."
            ),
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
            match=re.escape(
                "This scaler requires the user to provide a default "
                f"scaling factor for mole_frac_comp[N2], but no default scaling "
                "factor was set."
            ),
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
            ValueError,
            match=re.escape(
                "This scaler requires the user to provide a default "
                "scaling factor for pressure, but no default scaling "
                "factor was set."
            ),
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
            ValueError,
            match=re.escape(
                "This scaler requires the user to provide a default "
                "scaling factor for pressure, but no default scaling "
                "factor was set."
            ),
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
    def test_scale_variable_by_definition_constraint_not_variable(self, model):
        sb = CustomScalerBase()
        with pytest.raises(
            TypeError,
            match=re.escape(
                f"enthalpy_eq is not a variable, but instead <class 'pyomo.core.base.constraint.ScalarConstraint'>"
            ),
        ):
            sb.scale_variable_by_definition_constraint(
                model.enthalpy_eq,
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
    def test_scale_variable_by_definition_variable_indexed(self, model):
        sb = CustomScalerBase()

        model.foo = Var([1, 2, 3])

        @model.Constraint([1, 2, 3])
        def bar(b, idx):
            return b.foo[idx] == 0

        with pytest.raises(
            TypeError,
            match=re.escape(
                f"Variable foo is indexed. Call with VarData " "children instead."
            ),
        ):
            sb.scale_variable_by_definition_constraint(
                model.foo,
                model.bar,
            )

    @pytest.mark.unit
    def test_scale_variable_by_definition_constraint_indexed(self, model):
        sb = CustomScalerBase()

        model.foo = Var([1, 2, 3])

        @model.Constraint([1, 2, 3])
        def bar(b, idx):
            return b.foo[idx] == idx

        with pytest.raises(
            TypeError,
            match=re.escape(
                f"Constraint bar is indexed. Call with ConstraintData "
                "children instead."
            ),
        ):
            sb.scale_variable_by_definition_constraint(
                model.foo[1],
                model.bar,
            )

    @pytest.mark.unit
    def test_scale_variable_by_definition_data_children(self, model):
        sb = CustomScalerBase()

        model.foo = Var([1, 2, 3])

        @model.Constraint([1, 2, 3])
        def bar(b, idx):
            return b.foo[idx] == idx

        for idx in model.foo:
            sb.scale_variable_by_definition_constraint(
                model.foo[idx],
                model.bar[idx],
            )

        for idx in model.foo:
            assert model.scaling_factor[model.foo[idx]] == 1 / idx
            assert model.foo[idx].value is None

    @pytest.mark.unit
    def test_scale_variable_by_definition_foo_in_expression_scalar(self, model):
        sb = CustomScalerBase()

        model.foo = Var()

        @model.Expression()
        def biz(b):
            return b.foo

        sb.set_component_scaling_factor(model.biz, 42)

        @model.Constraint()
        def bar(b):
            return b.biz == 1

        sb.scale_variable_by_definition_constraint(
            model.foo,
            model.bar,
        )

        assert model.scaling_factor[model.foo] == 1
        assert model.foo.value is None

    @pytest.mark.unit
    def test_scale_variable_by_definition_foo_in_expression_indexed(self, model):
        sb = CustomScalerBase()

        model.foo = Var([1, 2, 3])

        @model.Expression([1, 2, 3])
        def biz(b, idx):
            return b.foo[idx]

        for idx in model.biz:
            sb.set_component_scaling_factor(model.biz[idx], 42)

        @model.Constraint([1, 2, 3])
        def bar(b, idx):
            return b.biz[idx] == idx

        for idx in model.foo:
            sb.scale_variable_by_definition_constraint(
                model.foo[idx],
                model.bar[idx],
            )

        for idx in model.foo:
            assert model.scaling_factor[model.foo[idx]] == 1 / idx
            assert model.foo[idx].value is None

    @pytest.mark.unit
    def test_scale_variable_by_definition_scaling_hint_scalar(self, model):
        sb = CustomScalerBase()

        model.foo = Var()

        @model.Expression()
        def biz(b):
            return 1

        sb.set_component_scaling_factor(model.biz, 42)

        @model.Constraint()
        def bar(b):
            return b.foo == b.biz

        sb.scale_variable_by_definition_constraint(
            model.foo,
            model.bar,
        )

        assert model.scaling_factor[model.foo] == pytest.approx(42)
        assert model.foo.value is None

    @pytest.mark.unit
    def test_scale_variable_by_definition_scaling_hint_indexed(self, model):
        sb = CustomScalerBase()

        model.foo = Var([1, 2, 3])

        @model.Expression([1, 2, 3])
        def biz(b, idx):
            return idx

        scaling_hints = [11, 13, 17]
        for idx in model.biz:
            sb.set_component_scaling_factor(model.biz[idx], scaling_hints[idx - 1])

        @model.Constraint([1, 2, 3])
        def bar(b, idx):
            return b.foo[idx] == b.biz[idx]

        for idx in model.foo:
            sb.scale_variable_by_definition_constraint(
                model.foo[idx],
                model.bar[idx],
            )

        for idx in model.foo:
            assert model.scaling_factor[model.foo[idx]] == pytest.approx(
                scaling_hints[idx - 1]
            )
            assert model.foo[idx].value is None

    @pytest.mark.unit
    def test_scale_variable_by_definition_no_hint_scalar(self, model):
        sb = CustomScalerBase()

        model.foo = Var()

        @model.Expression()
        def biz(b):
            return 1

        @model.Constraint()
        def bar(b):
            return b.foo == b.biz

        sb.scale_variable_by_definition_constraint(
            model.foo,
            model.bar,
        )

        assert model.scaling_factor[model.foo] == pytest.approx(1)
        assert model.foo.value is None

    @pytest.mark.unit
    def test_scale_variable_by_definition_no_hint_indexed(self, model):
        sb = CustomScalerBase()

        model.foo = Var([1, 2, 3])

        @model.Expression([1, 2, 3])
        def biz(b, idx):
            return idx

        @model.Constraint([1, 2, 3])
        def bar(b, idx):
            return b.foo[idx] == b.biz[idx]

        for idx in model.foo:
            sb.scale_variable_by_definition_constraint(
                model.foo[idx],
                model.bar[idx],
            )

        for idx in model.foo:
            assert model.scaling_factor[model.foo[idx]] == pytest.approx(1 / idx)
            assert model.foo[idx].value is None

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
    def test_scale_variable_by_definition_constraint_min_scaling_factors(self, model):
        sb = CustomScalerBase()
        model.foo = Var(initialize=1e8)
        model.bar = Var(initialize=1e16)
        model.biz = Var()

        @model.Constraint()
        def con(b):
            return b.biz == b.foo**2 / b.bar

        sb.scale_variable_by_definition_constraint(
            model.biz,
            model.con,
        )
        assert model.foo.value == 1e8
        assert model.bar.value == 1e16
        assert model.biz.value is None

        assert model.scaling_factor[model.biz] == pytest.approx(1e-6)

    @pytest.mark.unit
    def test_scale_variable_by_definition_constraint_max_scaling_factors(self, model):
        sb = CustomScalerBase()
        model.foo = Var(initialize=1e-8)
        model.bar = Var(initialize=1e-16)
        model.biz = Var()

        @model.Constraint()
        def con(b):
            return b.biz == b.foo**2 / b.bar

        sb.scale_variable_by_definition_constraint(
            model.biz,
            model.con,
        )
        assert model.foo.value == 1e-8
        assert model.bar.value == 1e-16
        assert model.biz.value is None

        assert model.scaling_factor[model.biz] == pytest.approx(1e6)

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
            match=re.escape(
                "This scaler requires the user to provide a default "
                "scaling factor for mole_frac_eqn[N2], but no default scaling "
                "factor was set."
            ),
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
