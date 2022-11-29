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
Tests for InitializerBase class
"""
import pytest
import types
import os

from pyomo.environ import ConcreteModel, Constraint, Var

from idaes.core.initialization.initializer_base import (
    InitializerBase,
    InitializationStatus,
)

from idaes.core.util.exceptions import InitializationError

path = os.path.dirname(__file__)
fname = os.path.join(path, "init_example.json")

__author__ = "Andrew Lee"


class TestSubMethods:
    @pytest.mark.unit
    def test_init(self):
        initializer = InitializerBase()

        assert initializer.postcheck_summary == {}

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

        assert state is initializer.initial_state

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
        assert expected == initializer.initial_state["unknown"]

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
        assert initializer.status == InitializationStatus.Error

    @pytest.mark.unit
    def test_load_values_dict_non_var(self, model):
        initializer = InitializerBase()

        init_dict = {"c1": 10}

        with pytest.raises(
            TypeError,
            match="Component c1 is not a Var. Initial guesses should only contain values for variables.",
        ):
            initializer.load_initial_guesses(model, initial_guesses=init_dict)
        assert initializer.status == InitializationStatus.Error

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
    def test_initialization_routine(self, model):
        initializer = InitializerBase()

        with pytest.raises(NotImplementedError):
            initializer.initialization_routine(model)
        assert initializer.status == InitializationStatus.Error

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
        assert initializer.status == InitializationStatus.DoF

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

        assert initializer.status == InitializationStatus.Failed
        assert len(initializer.postcheck_summary["uninitialized_vars"]) == 2
        assert len(initializer.postcheck_summary["unconverged_constraints"]) == 4

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

        assert initializer.status == InitializationStatus.Failed
        assert len(initializer.postcheck_summary["uninitialized_vars"]) == 0
        assert len(initializer.postcheck_summary["unconverged_constraints"]) == 4

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

        assert initializer.status == InitializationStatus.Failed
        assert len(initializer.postcheck_summary["uninitialized_vars"]) == 0
        assert len(initializer.postcheck_summary["unconverged_constraints"]) == 2

    @pytest.mark.unit
    def test_postcheck_correct(self, model):
        model.v1.set_value(2)
        model.v2.set_value(2)

        initializer = InitializerBase()
        initializer.postcheck(model)

        assert initializer.status == InitializationStatus.Ok
        assert len(initializer.postcheck_summary["uninitialized_vars"]) == 0
        assert len(initializer.postcheck_summary["unconverged_constraints"]) == 0

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

        assert initializer.status == InitializationStatus.Failed
        assert len(initializer.postcheck_summary["uninitialized_vars"]) == 0
        assert len(initializer.postcheck_summary["unconverged_constraints"]) == 1

        # Change tolerance and recheck
        initializer.config.constraint_tolerance = 1e-3

        initializer.postcheck(model)

        assert initializer.status == InitializationStatus.Ok
        assert len(initializer.postcheck_summary["uninitialized_vars"]) == 0
        assert len(initializer.postcheck_summary["unconverged_constraints"]) == 0
