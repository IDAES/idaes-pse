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
Tests for InitializerBase class
"""
import pytest
import types
import os

from pyomo.environ import Block, ConcreteModel, Constraint, value, Var

from idaes.core.base.process_base import ProcessBaseBlock
from idaes.core.initialization.initializer_base import (
    InitializerBase,
    InitializationStatus,
    ModularInitializerBase,
)
from idaes.core.initialization.block_triangularization import (
    BlockTriangularizationInitializer,
)

from idaes.core.util.exceptions import InitializationError
import idaes.logger as idaeslog

path = os.path.dirname(__file__)
fname = os.path.join(path, "init_example.json")

__author__ = "Andrew Lee"


class TestSubMethods:
    @pytest.mark.unit
    def test_init(self):
        initializer = InitializerBase()

        assert initializer.summary == {}
        assert initializer.initial_state == {}

        assert hasattr(initializer, "config")

        assert "constraint_tolerance" in initializer.config

    @pytest.fixture
    def model(self):
        m = ConcreteModel()

        m.v1 = Var()
        m.v2 = Var()

        m.c1 = Constraint(expr=m.v1 == 2)
        m.c2 = Constraint(expr=m.v1 == m.v2)
        m.c3 = Constraint(expr=m.v2 >= 0)
        m.c4 = Constraint(expr=m.v1 <= 5)

        return m

    @pytest.mark.unit
    def test_get_initial_state(self, model):
        model.v1.fix(42)
        model.v2.set_value(43)

        initializer = InitializerBase()
        state = initializer.get_current_state(model)

        assert state is initializer.initial_state[model]

        expected = {
            "__type__": "<class 'pyomo.core.base.PyomoModel.ConcreteModel'>",
            "__id__": 0,
            "active": True,
            "data": {
                "None": {
                    "__type__": "<class 'pyomo.core.base.PyomoModel.ConcreteModel'>",
                    "__id__": 1,
                    "active": True,
                    "__pyomo_components__": {
                        "v1": {
                            "__type__": "<class 'pyomo.core.base.var.ScalarVar'>",
                            "__id__": 2,
                            "data": {
                                "None": {
                                    "__type__": "<class 'pyomo.core.base.var.ScalarVar'>",
                                    "__id__": 3,
                                    "fixed": True,
                                    "value": 42,
                                }
                            },
                        },
                        "v2": {
                            "__type__": "<class 'pyomo.core.base.var.ScalarVar'>",
                            "__id__": 4,
                            "data": {
                                "None": {
                                    "__type__": "<class 'pyomo.core.base.var.ScalarVar'>",
                                    "__id__": 5,
                                    "fixed": False,
                                    "value": 43,
                                }
                            },
                        },
                        "c1": {
                            "__type__": "<class 'pyomo.core.base.constraint.ScalarConstraint'>",
                            "__id__": 6,
                            "active": True,
                            "data": {
                                "None": {
                                    "__type__": "<class 'pyomo.core.base.constraint.ScalarConstraint'>",
                                    "__id__": 7,
                                    "active": True,
                                }
                            },
                        },
                        "c2": {
                            "__type__": "<class 'pyomo.core.base.constraint.ScalarConstraint'>",
                            "__id__": 8,
                            "active": True,
                            "data": {
                                "None": {
                                    "__type__": "<class 'pyomo.core.base.constraint.ScalarConstraint'>",
                                    "__id__": 9,
                                    "active": True,
                                }
                            },
                        },
                        "c3": {
                            "__type__": "<class 'pyomo.core.base.constraint.ScalarConstraint'>",
                            "__id__": 10,
                            "active": True,
                            "data": {
                                "None": {
                                    "__type__": "<class 'pyomo.core.base.constraint.ScalarConstraint'>",
                                    "__id__": 11,
                                    "active": True,
                                }
                            },
                        },
                        "c4": {
                            "__type__": "<class 'pyomo.core.base.constraint.ScalarConstraint'>",
                            "__id__": 12,
                            "active": True,
                            "data": {
                                "None": {
                                    "__type__": "<class 'pyomo.core.base.constraint.ScalarConstraint'>",
                                    "__id__": 13,
                                    "active": True,
                                }
                            },
                        },
                    },
                }
            },
        }

        assert expected == state["unknown"]

    @pytest.mark.unit
    def test_get_and_restore_initial_state(self, model):
        # Set an initial state
        model.v1.fix(10)
        model.c4.deactivate()

        initializer = InitializerBase()

        # Store the initial state
        initializer.get_current_state(model)

        expected = {
            "__type__": "<class 'pyomo.core.base.PyomoModel.ConcreteModel'>",
            "__id__": 0,
            "active": True,
            "data": {
                "None": {
                    "__type__": "<class 'pyomo.core.base.PyomoModel.ConcreteModel'>",
                    "__id__": 1,
                    "active": True,
                    "__pyomo_components__": {
                        "v1": {
                            "__type__": "<class 'pyomo.core.base.var.ScalarVar'>",
                            "__id__": 2,
                            "data": {
                                "None": {
                                    "__type__": "<class 'pyomo.core.base.var.ScalarVar'>",
                                    "__id__": 3,
                                    "fixed": True,
                                    "value": 10,
                                }
                            },
                        },
                        "v2": {
                            "__type__": "<class 'pyomo.core.base.var.ScalarVar'>",
                            "__id__": 4,
                            "data": {
                                "None": {
                                    "__type__": "<class 'pyomo.core.base.var.ScalarVar'>",
                                    "__id__": 5,
                                    "fixed": False,
                                    "value": None,
                                }
                            },
                        },
                        "c1": {
                            "__type__": "<class 'pyomo.core.base.constraint.ScalarConstraint'>",
                            "__id__": 6,
                            "active": True,
                            "data": {
                                "None": {
                                    "__type__": "<class 'pyomo.core.base.constraint.ScalarConstraint'>",
                                    "__id__": 7,
                                    "active": True,
                                }
                            },
                        },
                        "c2": {
                            "__type__": "<class 'pyomo.core.base.constraint.ScalarConstraint'>",
                            "__id__": 8,
                            "active": True,
                            "data": {
                                "None": {
                                    "__type__": "<class 'pyomo.core.base.constraint.ScalarConstraint'>",
                                    "__id__": 9,
                                    "active": True,
                                }
                            },
                        },
                        "c3": {
                            "__type__": "<class 'pyomo.core.base.constraint.ScalarConstraint'>",
                            "__id__": 10,
                            "active": True,
                            "data": {
                                "None": {
                                    "__type__": "<class 'pyomo.core.base.constraint.ScalarConstraint'>",
                                    "__id__": 11,
                                    "active": True,
                                }
                            },
                        },
                        "c4": {
                            "__type__": "<class 'pyomo.core.base.constraint.ScalarConstraint'>",
                            "__id__": 12,
                            "active": False,
                            "data": {
                                "None": {
                                    "__type__": "<class 'pyomo.core.base.constraint.ScalarConstraint'>",
                                    "__id__": 13,
                                    "active": False,
                                }
                            },
                        },
                    },
                }
            },
        }
        assert expected == initializer.initial_state[model]["unknown"]

        # Make some more changes to the state
        model.v1.set_value(21)
        model.v1.unfix()
        model.v2.fix(22)
        model.c1.deactivate()
        model.c4.activate()

        initializer.restore_model_state(model)

        # Check that state was reverted correctly
        assert model.v1.value == 10  # Value should have changed back to original
        assert model.v1.fixed  # Should be fixed again
        assert model.v2.value == 22  # Value should not have changed
        assert not model.v2.fixed  # Should have been unfixed
        assert model.c1.active
        assert model.c2.active
        assert model.c3.active
        assert not model.c4.active

    @pytest.mark.unit
    def test_load_value_from_file(self):
        m = ConcreteModel()

        m.v1 = Var()
        m.v2 = Var()
        m.v3 = Var()
        m.v4 = Var()

        # Fix value of v2, it should not get overwritten
        m.v2.fix(10)

        initializer = InitializerBase()

        # Load values from example file, all values are 4
        initializer.load_initial_guesses(m, json_file=fname)

        assert m.v1.value == 4
        assert m.v2.value == 10  # Should retain fixed value
        assert m.v3.value == 4
        assert m.v4.value == 4

    @pytest.mark.unit
    def test_load_values_dict(self, model):
        model.v3 = Var(["a", "b"])
        model.v4 = Var(["a", "b"])
        initializer = InitializerBase()

        init_dict = {"v1": 10, "v2": 20, "v3[a]": 30, "v3[b]": 40, "v4": 50}

        initializer.load_initial_guesses(model, initial_guesses=init_dict)

        assert model.v1.value == 10
        assert model.v2.value == 20
        assert model.v3["a"].value == 30
        assert model.v3["b"].value == 40
        assert model.v4["a"].value == 50
        assert model.v4["b"].value == 50

    @pytest.mark.unit
    def test_load_values_dict_and_file(self, model):
        initializer = InitializerBase()

        with pytest.raises(
            ValueError,
            match="Cannot provide both a set of initial guesses and a json file to load.",
        ):
            initializer.load_initial_guesses(
                model, initial_guesses="foo", json_file="bar"
            )
        assert initializer.summary[model]["status"] == InitializationStatus.Error

    @pytest.mark.unit
    def test_load_values_dict_non_var(self, model):
        initializer = InitializerBase()

        init_dict = {"c1": 10}

        with pytest.raises(
            TypeError,
            match="Component c1 is not a Var. Initial guesses should only contain values for variables.",
        ):
            initializer.load_initial_guesses(model, initial_guesses=init_dict)
        assert initializer.summary[model]["status"] == InitializationStatus.Error

    @pytest.mark.unit
    def test_load_value_none(self, caplog):
        caplog.set_level(
            idaeslog.DEBUG,
            logger="idaes",
        )

        m = ConcreteModel()
        m.v1 = Var()

        initializer = InitializerBase()
        initializer.load_initial_guesses(m)

        expected = "No initial guesses provided during initialization of model unknown."
        assert expected in caplog.text

        assert m.v1.value is None

    @pytest.mark.unit
    def test_fix_states(self, model):
        # Create a dummy method to fix states on test model
        def fix_initialization_states(blk):
            blk.v1.fix(12)

        model.fix_initialization_states = types.MethodType(
            fix_initialization_states, model
        )

        initializer = InitializerBase()
        initializer.fix_initialization_states(model)

        assert model.v1.fixed
        assert not model.v2.fixed

    @pytest.mark.unit
    def test_fix_states_no_method(self, model, caplog):
        caplog.set_level(
            idaeslog.DEBUG,
            logger="idaes",
        )

        initializer = InitializerBase()
        initializer.fix_initialization_states(model)

        expected = "Model unknown does not have a fix_initialization_states method - attempting to continue."
        assert expected in caplog.text

        assert not model.v1.fixed
        assert not model.v2.fixed

    @pytest.mark.unit
    def test_initialization_routine(self, model):
        initializer = InitializerBase()

        with pytest.raises(NotImplementedError):
            initializer.initialization_routine(model)
        assert initializer.summary[model]["status"] == InitializationStatus.Error

    @pytest.mark.unit
    def test_precheck_fail(self, model):
        model.v1.fix(10)

        initializer = InitializerBase()
        with pytest.raises(
            InitializationError,
            match="Degrees of freedom for unknown were not equal to zero during "
            "initialization \(DoF = -1\).",
        ):
            initializer.precheck(model)
        assert initializer.summary[model]["status"] == InitializationStatus.DoF
        assert initializer.summary[model]["DoF"] == -1

    @pytest.mark.unit
    def test_precheck_ok(self, model):
        initializer = InitializerBase()
        initializer.precheck(model)

        assert initializer.summary[model]["status"] == InitializationStatus.none
        assert initializer.summary[model]["DoF"] == 0

    @pytest.mark.unit
    def test_postcheck_vars_none(self, model):
        initializer = InitializerBase()

        with pytest.raises(
            InitializationError,
            match=f"unknown failed to initialize successfully: uninitialized variables or "
            "unconverged equality constraints detected. Please check postcheck summary for "
            "more information.",
        ):
            initializer.postcheck(model)

        assert initializer.summary[model]["status"] == InitializationStatus.Failed
        assert len(initializer.summary[model]["uninitialized_vars"]) == 2
        assert len(initializer.summary[model]["unconverged_constraints"]) == 4

    @pytest.mark.unit
    def test_postcheck_vars_wrong(self, model):
        model.v1.set_value(10)
        model.v2.set_value(-10)

        initializer = InitializerBase()

        with pytest.raises(
            InitializationError,
            match=f"unknown failed to initialize successfully: uninitialized variables or "
            "unconverged equality constraints detected. Please check postcheck summary for "
            "more information.",
        ):
            initializer.postcheck(model)

        assert initializer.summary[model]["status"] == InitializationStatus.Failed
        assert len(initializer.summary[model]["uninitialized_vars"]) == 0
        assert len(initializer.summary[model]["unconverged_constraints"]) == 4

    @pytest.mark.unit
    def test_postcheck_partially_satisfied(self, model):
        model.v1.set_value(-10)
        model.v2.set_value(-10)

        initializer = InitializerBase()

        with pytest.raises(
            InitializationError,
            match=f"unknown failed to initialize successfully: uninitialized variables or "
            "unconverged equality constraints detected. Please check postcheck summary for "
            "more information.",
        ):
            initializer.postcheck(model)

        assert initializer.summary[model]["status"] == InitializationStatus.Failed
        assert len(initializer.summary[model]["uninitialized_vars"]) == 0
        assert len(initializer.summary[model]["unconverged_constraints"]) == 2

    @pytest.mark.unit
    def test_postcheck_correct(self, model):
        model.v1.set_value(2)
        model.v2.set_value(2)

        initializer = InitializerBase()
        initializer.postcheck(model)

        assert initializer.summary[model]["status"] == InitializationStatus.Ok
        assert len(initializer.summary[model]["uninitialized_vars"]) == 0
        assert len(initializer.summary[model]["unconverged_constraints"]) == 0

    @pytest.mark.unit
    def test_postcheck_near_tolerance(self, model):
        model.v1.set_value(2)
        model.v2.set_value(2.0001)

        initializer = InitializerBase()
        with pytest.raises(
            InitializationError,
            match=f"unknown failed to initialize successfully: uninitialized variables or "
            "unconverged equality constraints detected. Please check postcheck summary for "
            "more information.",
        ):
            initializer.postcheck(model)

        assert initializer.summary[model]["status"] == InitializationStatus.Failed
        assert len(initializer.summary[model]["uninitialized_vars"]) == 0
        assert len(initializer.summary[model]["unconverged_constraints"]) == 1

        # Change tolerance and recheck
        initializer.config.constraint_tolerance = 1e-3

        initializer.postcheck(model)

        assert initializer.summary[model]["status"] == InitializationStatus.Ok
        assert len(initializer.summary[model]["uninitialized_vars"]) == 0
        assert len(initializer.summary[model]["unconverged_constraints"]) == 0

    @pytest.mark.unit
    def test_postcheck_ignore_unused(self):
        m = ConcreteModel()
        m.v1 = Var(initialize=1)
        m.v2 = Var(initialize=2)
        m.v3 = Var()  # unused
        m.c = Constraint(expr=m.v1 * 2 == m.v2)

        initializer = InitializerBase()

        # By default ,will fail due to m.v3
        with pytest.raises(
            InitializationError,
            match=f"unknown failed to initialize successfully: uninitialized variables or "
            "unconverged equality constraints detected. Please check postcheck summary for "
            "more information.",
        ):
            initializer.postcheck(m)

        assert initializer.summary[m]["status"] == InitializationStatus.Failed
        assert initializer.summary[m]["uninitialized_vars"] == [m.v3]
        assert len(initializer.summary[m]["unconverged_constraints"]) == 0

        # Try again whilst skipping unused vars
        initializer.postcheck(m, exclude_unused_vars=True)

        assert initializer.summary[m]["status"] == InitializationStatus.Ok
        assert len(initializer.summary[m]["uninitialized_vars"]) == 0
        assert len(initializer.summary[m]["unconverged_constraints"]) == 0

    @pytest.mark.unit
    def test_update_summary(self):
        initializer = InitializerBase()

        assert initializer.summary == {}

        initializer._update_summary("foo", "bar", "baz")

        assert initializer.summary["foo"]["bar"] == "baz"
        assert initializer.summary["foo"]["status"] == InitializationStatus.none

    @pytest.mark.unit
    def test_load_values_from_dict_not_var(self):
        m = ConcreteModel()
        m.v = Var()
        m.c = Constraint(expr=m.v == 0)

        initializer = InitializerBase()

        with pytest.raises(
            TypeError,
            match="Component c is not a Var. Initial guesses should only contain values for variables.",
        ):
            initializer._load_values_from_dict(m, {m.c: 10})

        assert initializer.summary[m]["status"] == InitializationStatus.Error

    @pytest.mark.unit
    def test_load_values_from_dict_indexed(self):
        m = ConcreteModel()
        m.v = Var(["a", "b"])

        initializer = InitializerBase()

        initializer._load_values_from_dict(m, {m.v: 10})

        assert m.v["a"].value == 10
        assert m.v["b"].value == 10

    @pytest.mark.unit
    def test_load_values_from_dict_indexed_fixed(self):
        m = ConcreteModel()
        m.v = Var(["a", "b"])
        m.v["b"].fix(20)

        initializer = InitializerBase()

        with pytest.raises(
            InitializationError,
            match="Attempted to change the value of fixed variable v\[b\]. "
            "Initialization from initial guesses does not support changing the value "
            "of fixed variables.",
        ):
            initializer._load_values_from_dict(m, {m.v: 10})

        assert initializer.summary[m]["status"] == InitializationStatus.Error

    @pytest.mark.unit
    def test_load_values_from_dict_indexed_fixed_ignore(self, caplog):
        caplog.set_level(
            idaeslog.DEBUG,
            logger="idaes",
        )

        m = ConcreteModel()
        m.v = Var(["a", "b"])
        m.v["b"].fix(20)

        initializer = InitializerBase()
        initializer._load_values_from_dict(m, {m.v: 10}, exception_on_fixed=False)

        assert m.v["a"].value == 10
        assert m.v["b"].value == 20

        expected = "Found initial guess for fixed Var v[b] - ignoring."
        assert expected in caplog.text

    @pytest.mark.unit
    def test_load_values_from_dict_scalar(self):
        m = ConcreteModel()
        m.v = Var(["a", "b"])

        initializer = InitializerBase()

        initializer._load_values_from_dict(m, {'v["a"]': 10})

        assert m.v["a"].value == 10
        assert m.v["b"].value == None

    @pytest.mark.unit
    def test_load_values_from_dict_scalar_fixed(self):
        m = ConcreteModel()
        m.v = Var(["a", "b"])
        m.v["a"].fix(20)

        initializer = InitializerBase()

        with pytest.raises(
            InitializationError,
            match="Attempted to change the value of fixed variable v\[a\]. "
            "Initialization from initial guesses does not support changing the value "
            "of fixed variables.",
        ):
            initializer._load_values_from_dict(m, {'v["a"]': 10})

        assert initializer.summary[m]["status"] == InitializationStatus.Error

    @pytest.mark.unit
    def test_load_values_from_dict_scalar_fixed_ignore(self, caplog):
        caplog.set_level(
            idaeslog.DEBUG,
            logger="idaes",
        )

        m = ConcreteModel()
        m.v = Var()
        m.v.fix(20)

        initializer = InitializerBase()
        initializer._load_values_from_dict(m, {"v": 10}, exception_on_fixed=False)

        assert m.v.value == 20

        expected = "Found initial guess for fixed Var v - ignoring."
        assert expected in caplog.text


class DummyInit:
    def plugin_prepare(self, plugin, output="default"):
        plugin._test = output

    def plugin_initialize(self, model, **kwargs):
        model._plugin_initialized = True

    def plugin_finalize(self, plugin, output="default"):
        plugin._test2 = output

    def initialize(self, model, **kwargs):
        model._initialized = True


class TestModularInitializerBase:
    @pytest.mark.unit
    def test_base_attributed(self):
        initializer = ModularInitializerBase()

        assert initializer.submodel_initializers == {}
        assert initializer.config.default_submodel_initializer is None

        assert initializer._solver is None
        assert initializer.config.solver is None
        assert initializer.config.solver_options == {}

    @pytest.mark.unit
    def test_get_submodel_initializer_specific_model(self):
        m = ConcreteModel()

        initializer = ModularInitializerBase()
        initializer.submodel_initializers[m] = DummyInit

        assert isinstance(initializer.get_submodel_initializer(m), DummyInit)

    @pytest.mark.unit
    def test_get_submodel_initializer_model_type(self):
        m = ConcreteModel()

        initializer = ModularInitializerBase()
        initializer.submodel_initializers[ConcreteModel] = DummyInit

        assert isinstance(initializer.get_submodel_initializer(m), DummyInit)

    @pytest.mark.unit
    def test_get_submodel_initializer_model_default(self):
        m = ConcreteModel()
        m.default_initializer = DummyInit

        initializer = ModularInitializerBase()

        assert isinstance(initializer.get_submodel_initializer(m), DummyInit)

    @pytest.mark.unit
    def test_get_submodel_initializer_global_default(self):
        m = ConcreteModel()

        initializer = ModularInitializerBase(default_submodel_initializer=DummyInit)
        assert initializer.config.default_submodel_initializer == DummyInit

        assert isinstance(initializer.get_submodel_initializer(m), DummyInit)

    @pytest.mark.unit
    def test_get_submodel_initializer_model_w_params(self):
        class DummyParam:
            pass

        dummy_param = DummyParam()

        m = ConcreteModel()
        m.params = dummy_param

        initializer = ModularInitializerBase()
        initializer.submodel_initializers[dummy_param] = DummyInit

        assert isinstance(initializer.get_submodel_initializer(m), DummyInit)

    @pytest.mark.unit
    def test_get_submodel_initializer_none(self, caplog):
        m = ConcreteModel()

        initializer = ModularInitializerBase()

        assert initializer.get_submodel_initializer(m) is None
        expected = "No Initializer found for submodel unknown - attempting to continue."
        assert expected in caplog.text

    @pytest.mark.unit
    def test_get_submodel_initializer_process_base_default(self):
        m = ConcreteModel()
        m.b = ProcessBaseBlock()

        initializer = ModularInitializerBase()
        assert isinstance(
            initializer.get_submodel_initializer(m.b), BlockTriangularizationInitializer
        )

    @pytest.mark.unit
    def test_get_submodel_initializer_priorit(self):
        # Progressively add higher priority initializers and ensure they are returned
        class DummyParam:
            def __init__(self):
                self.name = "dummy"

        dummy_param = DummyParam()

        m = ConcreteModel()
        m.params = dummy_param

        initializer = ModularInitializerBase()

        # No default set
        assert initializer.get_submodel_initializer(m) is None

        # Add a global default
        initializer.config.default_submodel_initializer = "global"
        assert initializer.get_submodel_initializer(m) == "global"

        # Parameter block
        initializer.submodel_initializers[dummy_param] = "params"
        assert initializer.get_submodel_initializer(m) == "params"

        # Model default
        m.default_initializer = "model_default"
        assert initializer.get_submodel_initializer(m) == "model_default"

        # Model type
        initializer.submodel_initializers[ConcreteModel] = "type"
        assert initializer.get_submodel_initializer(m) == "type"

        # Specific model
        initializer.submodel_initializers[m] = "specific_model"
        assert initializer.get_submodel_initializer(m) == "specific_model"

    @pytest.mark.unit
    def test_add_submodel_initializer(self):
        initializer = ModularInitializerBase()
        assert initializer.submodel_initializers == {}

        initializer.add_submodel_initializer("foo", "bar")
        assert initializer.submodel_initializers == {"foo": "bar"}

    @pytest.mark.unit
    def test_prepare_plugins_none(self):
        m = ConcreteModel()
        m.initialization_order = [m]

        initializer = ModularInitializerBase()

        args = {"foo": "bar"}
        subinit, plugin_args = initializer.prepare_plugins(m, args)

        assert subinit == {}
        assert plugin_args == {"foo": "bar"}
        # Make sure we didn't change the original args
        assert args == {"foo": "bar"}

    @pytest.mark.unit
    def test_prepare_plugins_no_args(self):
        m = ConcreteModel()
        m.b = Block()
        m.initialization_order = [m, m.b]

        initializer = ModularInitializerBase()
        initializer.add_submodel_initializer(m.b, DummyInit)

        args = {"foo": "bar"}
        subinit, plugin_args = initializer.prepare_plugins(m, args)

        assert len(subinit) == 1
        assert isinstance(subinit[m.b], DummyInit)
        assert plugin_args == {"foo": "bar", m.b: {}}
        assert m.b._test == "default"
        # Make sure we didn't change the original args
        assert args == {"foo": "bar"}

    @pytest.mark.unit
    def test_prepare_plugins_w_args(self):
        m = ConcreteModel()
        m.b = Block()
        m.initialization_order = [m, m.b]

        initializer = ModularInitializerBase()
        initializer.add_submodel_initializer(m.b, DummyInit)

        args = {"foo": "bar", m.b: {"output": "checkval"}}
        subinit, plugin_args = initializer.prepare_plugins(m, args)

        assert len(subinit) == 1
        assert isinstance(subinit[m.b], DummyInit)
        assert plugin_args == {"foo": "bar", m.b: {"output": "checkval"}}
        assert m.b._test == "checkval"
        # Make sure we didn't change the original args
        assert args == {"foo": "bar", m.b: {"output": "checkval"}}

    @pytest.mark.unit
    def test_solve_full_model_no_plugins(self):
        m = ConcreteModel()
        m.initialization_order = [m]

        initializer = ModularInitializerBase()

        results = initializer.solve_full_model(m, "test_results")

        assert results == "test_results"

    @pytest.mark.unit
    def test_solve_full_model_w_plugins(self):
        m = ConcreteModel()
        m.b = Block()
        m.b.v = Var()
        m.b.c = Constraint(expr=m.b.v == 12)

        m.initialization_order = [m, m.b]

        initializer = ModularInitializerBase()

        results = initializer.solve_full_model(m, "test_results")

        assert value(m.b.v) == pytest.approx(12, rel=1e-8)
        assert isinstance(results, dict)
        assert "Solution" in results

    @pytest.mark.unit
    def test_cleanup_none(self):
        m = ConcreteModel()
        m.initialization_order = [m]

        initializer = ModularInitializerBase()

        args = {"foo": "bar"}
        subinit = {}
        initializer.cleanup(m, args, subinit)

        # No-op - what to assert?

    @pytest.mark.unit
    def test_cleanup_no_args(self):
        m = ConcreteModel()
        m.b = Block()
        m.initialization_order = [m, m.b]

        initializer = ModularInitializerBase()

        args = {"foo": "bar", m.b: {}}
        subinit = {m.b: DummyInit()}
        initializer.cleanup(m, args, subinit)

        assert m.b._test2 == "default"

    @pytest.mark.unit
    def test_cleanup_w_args(self):
        m = ConcreteModel()
        m.b = Block()
        m.initialization_order = [m, m.b]

        initializer = ModularInitializerBase()

        args = {"foo": "bar", m.b: {"output": "forty-two"}}
        subinit = {m.b: DummyInit()}
        initializer.cleanup(m, args, subinit)

        assert m.b._test2 == "forty-two"
