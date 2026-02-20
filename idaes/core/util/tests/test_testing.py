#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Who will test the testers?

This file was created with the assistance of AI
"""

__author__ = "Douglas Allan"


import pytest
from pyomo.environ import Block, ConcreteModel, Var, Expression, Set
from idaes.core.util.testing import assert_solution_equivalent

##############################################################
#### Test basic functionality on a flat Pyomo model
##############################################################


class TestAssertSolutionEquivalent:
    @pytest.fixture
    def model(self):
        """Creates a standard Pyomo model for testing."""
        m = ConcreteModel()
        m.i = Set(initialize=["A", "B", "C"])

        # Indexed Variable
        m.x = Var(m.i, initialize={"A": 1.0, "B": 2.5, "C": 3.0})

        # Unindexed Variable
        m.y = Var(initialize=100.0)

        # For testing absolute tolerance
        m.z = Var(initialize=0.001)

        # Indexed Expression
        m.x_squared = Expression(m.i, rule=lambda model, i: model.x[i] ** 2)

        # A component that is not a Var or Expression
        m.my_set = Set(initialize=[1, 2, 3])
        return m

    # --- Test Cases ---
    @pytest.mark.unit
    def test_all_values_correct(self, model):
        """
        Tests that assert_solution_equivalent passes when all model values
        match the expected results within tolerance.
        """
        expected_results = {
            "x": {
                "A": (1.0, 1e-6, 1e-8),
                "B": (2.5, 1e-6, 1e-8),
                "C": (3.0, 1e-6, 1e-8),
            },
            "y": {None: (100.0, 1e-6, 1e-8)},
            "x_squared": {
                "A": (1.0, 1e-6, 1e-8),
                "B": (6.25, 1e-6, 1e-8),
                "C": (9.0, 1e-6, 1e-8),
            },
        }
        # This should execute without raising an AssertionError
        assert_solution_equivalent(model, expected_results)

    @pytest.mark.unit
    def test_single_failure(self, model):
        """
        Tests that a single incorrect value is correctly identified and
        raises an AssertionError with a proper message.
        """
        expected_results = {
            "x": {
                "A": (1.0, 1e-6, 1e-8),
                "B": (2.5, 1e-6, 1e-8),
                "C": (3.1, 1e-6, 1e-8),  # Incorrect value
            },
        }
        with pytest.raises(AssertionError, match=r"Found 1 mismatch\(es\)"):
            assert_solution_equivalent(model, expected_results)

    @pytest.mark.unit
    def test_multiple_failures(self, model):
        """
        Tests that multiple incorrect values across different components
        are all reported in a single AssertionError.
        """
        expected_results = {
            "x": {
                "A": (1.1, 1e-7, None),  # Incorrect value
                "B": (2.5, 1e-7, None),
                "C": (3.0, 1e-7, None),
            },
            "y": {None: (101.0, 1e-7, None)},  # Incorrect value
            "x_squared": {
                "A": (1.0, 1e-7, None),
                "B": (6.25, 1e-7, None),
                "C": (9.1, 1e-7, None),  # Incorrect value
            },
        }
        with pytest.raises(AssertionError, match=r"Found 3 mismatch\(es\)"):
            assert_solution_equivalent(model, expected_results)

    @pytest.mark.unit
    def test_missing_component(self, model):
        """
        Tests that a component listed in expected_results but not present
        in the model is correctly flagged.
        """
        expected_results = {
            "zed": {  # This variable does not exist on the model
                "A": (1.0, 1e-6, 1e-8),
            },
        }
        with pytest.raises(AssertionError, match="Could not find object zed on model"):
            assert_solution_equivalent(model, expected_results)

    @pytest.mark.unit
    def test_missing_index(self, model):
        """
        Tests that an index listed in expected_results but not present
        in the component is correctly flagged.
        """
        expected_results = {
            "x": {
                "D": (4.0, 1e-6, 1e-8),  # Index 'D' does not exist in model.x
            },
        }
        with pytest.raises(AssertionError, match=r"Index:    D is absent"):
            assert_solution_equivalent(model, expected_results)

    @pytest.mark.unit
    def test_unindexed_component_failure(self, model):
        """
        Tests failure for a single, unindexed component.
        """
        expected_results = {
            "y": {
                None: (99.0, 1e-6, 1e-8),  # Incorrect value, correct is 100.0
            },
        }
        with pytest.raises(AssertionError, match=r"Found 1 mismatch\(es\)"):
            assert_solution_equivalent(model, expected_results)

    @pytest.mark.unit
    def test_wrong_component_type(self, model):
        """
        Tests that providing a component that is not a Var or Expression
        (e.g., a Set) is handled and reported.
        """
        expected_results = {
            "my_set": {  # my_set is a Pyomo Set, not a Var or Expression
                None: (1, None, None)
            }
        }
        with pytest.raises(
            AssertionError, match="Error: object my_set is not a Var or Expression"
        ):
            assert_solution_equivalent(model, expected_results)

    @pytest.mark.unit
    def test_relative_tolerance_pass(self, model):
        """
        Tests that a value passes when it's within the relative tolerance,
        even if it would fail a stricter absolute tolerance.
        - Expected: 100.0, Actual: 100.1
        - Relative tolerance (1e-2) allows a deviation up to 1.0.
        - Absolute tolerance (1e-8) would fail.
        The test should pass.
        """
        model.y.set_value(100.1)
        expected_results = {
            "y": {
                None: (100.0, 1e-2, 1e-8),  # rel=1e-2 (pass), abs=1e-8 (fail)
            },
        }
        # Should pass because it meets the relative tolerance
        assert_solution_equivalent(model, expected_results)

    @pytest.mark.unit
    def test_absolute_tolerance_pass(self, model):
        """
        Tests that a value passes when it's within the absolute tolerance,
        even if it would fail a stricter relative tolerance. This is key for
        comparisons with numbers close to zero.
        - Expected: 0.001, Actual: 0.0011
        - Relative tolerance (1e-2) would fail, as the difference is 10%.
        - Absolute tolerance (1e-3) allows a deviation up to 0.001.
        The test should pass.
        """
        model.z.set_value(0.0011)
        expected_results = {
            "z": {
                None: (0.001, 1e-2, 1e-3),  # rel=1e-2 (fail), abs=1e-3 (pass)
            },
        }
        # Should pass because it meets the absolute tolerance
        assert_solution_equivalent(model, expected_results)

    @pytest.mark.unit
    def test_relative_tolerance_fail(self, model):
        """
        Tests that a value just outside the relative tolerance fails.
        - Expected: 100.0, Actual: 101.1
        - Relative tolerance (1e-2) allows a deviation up to 1.0. The actual
        deviation is 1.1.
        The test should fail.
        """
        model.y.set_value(101.1)
        expected_results = {
            "y": {
                None: (100.0, 1e-2, 1e-8),
            },
        }
        with pytest.raises(AssertionError, match=r"Found 1 mismatch\(es\)"):
            assert_solution_equivalent(model, expected_results)

    @pytest.mark.unit
    def test_absolute_tolerance_fail(self, model):
        """
        Tests that a value just outside the absolute tolerance fails.
        - Expected: 0.001, Actual: 0.0021
        - Absolute tolerance (1e-3) allows a deviation up to 0.001. The actual
        deviation is 0.0011.
        The test should fail.
        """
        model.z.set_value(0.0021)
        expected_results = {
            "z": {
                None: (0.001, 1e-2, 1e-3),
            },
        }
        with pytest.raises(AssertionError, match=r"Found 1 mismatch\(es\)"):
            assert_solution_equivalent(model, expected_results)

    @pytest.mark.unit
    def test_default_tolerances_pass(self, model):
        """
        Tests that pytest.approx defaults are used when rel and abs are None.
        Default rel_tol is 1e-6.
        - Expected: 100.0, Actual: 100.0000001
        The test should pass.
        """
        model.y.set_value(100.0000001)
        expected_results = {
            "y": {
                None: (100.0, None, None),  # Use pytest.approx defaults
            },
        }
        assert_solution_equivalent(model, expected_results)

    @pytest.mark.unit
    def test_default_tolerances_fail(self, model):
        """
        Tests that pytest.approx defaults are used and can fail the test.
        Default rel_tol is 1e-6.
        - Expected: 100.0, Actual: 100.001
        The deviation of 0.001 is greater than the default relative tolerance
        (100.0 * 1e-6 = 0.0001).
        The test should fail.
        """
        model.y.set_value(100.001)
        expected_results = {
            "y": {
                None: (100.0, None, None),  # Use pytest.approx defaults
            },
        }
        with pytest.raises(AssertionError, match=r"Found 1 mismatch\(es\)"):
            assert_solution_equivalent(model, expected_results)


##############################################################
#### Test formatting in reports
##############################################################
class TestAssertSolutionEquivalentFormatting:
    @pytest.fixture
    def reporting_model(self):
        """Creates a Pyomo model with specific values for testing report formats."""
        m = ConcreteModel()
        m.normal_var = Var(initialize=123.456)
        m.large_var = Var(initialize=1.23456789e20)
        m.small_var = Var(initialize=1.23456789e-9)
        return m

    @pytest.mark.unit
    def test_reporting_format_relative_tolerance(self, reporting_model):
        """
        Tests that the 'e' format with correct significant figures is used
        when a relative tolerance is provided.
        - rel = 1e-4 should result in ceil(-log10(1e-4)) + 1 = 5 sig figs.
        - Format spec should be '.5e'.
        """
        expected_results = {
            "normal_var": {
                None: (123.0, 1e-4, None),  # Fails the check
            }
        }
        with pytest.raises(AssertionError) as excinfo:
            assert_solution_equivalent(reporting_model, expected_results)

        # Check that the numbers in the report use the correct 'e' format
        assert "Expected: 1.23000e+02" in str(excinfo.value)
        assert "Actual:   1.23456e+02" in str(excinfo.value)

    @pytest.mark.unit
    def test_reporting_format_absolute_tolerance(self, reporting_model):
        """
        Tests that the 'f' format is used when only absolute tolerance is provided.
        - abs = 1e-2 should result in ceil(-log10(1e-2)) + 1 = 3 decimal places.
        - Format spec should be '.3f'.
        """
        expected_results = {
            "normal_var": {
                None: (123.0, None, 1e-2),  # Fails the check
            }
        }
        with pytest.raises(AssertionError) as excinfo:
            assert_solution_equivalent(reporting_model, expected_results)

        # Check that the numbers in the report use the correct 'f' format
        assert "Expected: 123.000" in str(excinfo.value)
        assert "Actual:   123.456" in str(excinfo.value)

    @pytest.mark.unit
    def test_reporting_format_default_tolerance(self, reporting_model):
        """
        Tests that the default '.7e' format is used when both tolerances are None.
        """
        expected_results = {
            "normal_var": {
                None: (123.0, None, None),  # Fails the check
            }
        }
        with pytest.raises(AssertionError) as excinfo:
            assert_solution_equivalent(reporting_model, expected_results)

        # Check that the numbers in the report use the default '.7e' format
        assert "Expected: 1.2300000e+02" in str(excinfo.value)
        assert "Actual:   1.2345600e+02" in str(excinfo.value)

    @pytest.mark.unit
    def test_reporting_format_large_numbers(self, reporting_model):
        """
        Tests the reporting format for very large numbers using scientific notation.
        - rel = 1e-5 should result in '.6e' format.
        """
        expected_results = {
            "large_var": {
                None: (1.23e20, 1e-5, None),  # Fails the check
            }
        }
        with pytest.raises(AssertionError) as excinfo:
            assert_solution_equivalent(reporting_model, expected_results)

        # Check that the large numbers are formatted correctly
        assert "Expected: 1.230000e+20" in str(excinfo.value)
        assert "Actual:   1.234568e+20" in str(
            excinfo.value
        )  # Note pytest.approx rounding

    @pytest.mark.unit
    def test_reporting_format_small_numbers(self, reporting_model):
        """
        Tests the reporting format for very small numbers using scientific notation.
        - rel = 1e-3 should result in '.4e' format.
        """
        expected_results = {
            "small_var": {
                None: (1.2e-9, 1e-3, None),  # Fails the check
            }
        }
        with pytest.raises(AssertionError) as excinfo:
            assert_solution_equivalent(reporting_model, expected_results)

        # Check that the small numbers are formatted correctly
        assert "Expected: 1.2000e-09" in str(excinfo.value)
        assert "Actual:   1.2346e-09" in str(
            excinfo.value
        )  # Note pytest.approx rounding


##############################################################
#### Test complex block structures
##############################################################


class TestComplexBlockStructures:
    @pytest.fixture
    def complex_model(self):
        """
        Creates a Pyomo model with a variety of nested and indexed block structures
        for comprehensive testing.

        Structure:
        m.simple_b - A non-indexed block
            - m.simple_b.v_simple
            - m.simple_b.indexed_b_nested - An indexed block nested inside
                - m.simple_b.indexed_b_nested[i].v_nested_indexed

        m.indexed_b[i] - An indexed block
            - m.indexed_b[i].v_indexed
            - m.indexed_b[i].simple_b_nested - A simple block nested inside
                - m.indexed_b[i].simple_b_nested.v_nested_simple
        """
        m = ConcreteModel()
        m.i = Set(initialize=["X", "Y"])

        # 1. Simple (non-indexed) Block
        m.simple_b = Block()
        m.simple_b.v_simple = Var(initialize=10)

        # 4. Indexed Block nested within a Simple Block
        def indexed_nested_rule(b, i):
            b.v_nested_indexed = Var(initialize=20 if i == "X" else 21)

        m.simple_b.indexed_b_nested = Block(m.i, rule=indexed_nested_rule)

        # 2. Indexed Block
        def indexed_block_rule(b, i):
            # Variable indexed by the same set as the block
            b.v_indexed = Var(initialize=30 if i == "X" else 31)

            # 3. Simple Block nested within an Indexed Block
            b.simple_b_nested = Block()
            b.simple_b_nested.v_nested_simple = Var(initialize=40 if i == "X" else 41)

        m.indexed_b = Block(m.i, rule=indexed_block_rule)

        return m

    @pytest.mark.unit
    def test_on_simple_block(self, complex_model):
        """
        Tests asserting values on a simple, non-indexed block.
        The function is called on the block 'simple_b' itself.
        """
        target_block = complex_model.simple_b
        expected_results = {
            "v_simple": {None: (10, 1e-6, None)},
        }
        # This should pass
        assert_solution_equivalent(target_block, expected_results)

        # Now, test a failure on the same block
        expected_results_fail = {
            "v_simple": {None: (10.1, 1e-6, None)},
        }
        with pytest.raises(AssertionError, match="Found 1 mismatch"):
            assert_solution_equivalent(target_block, expected_results_fail)

        # Now test with full path
        expected_results = {
            "simple_b.v_simple": {None: (10, 1e-6, None)},
        }
        # This should pass
        assert_solution_equivalent(complex_model, expected_results)

        # Now, test a failure on the same block
        expected_results_fail = {
            "simple_b.v_simple": {None: (10.1, 1e-6, None)},
        }
        with pytest.raises(AssertionError, match="Found 1 mismatch"):
            assert_solution_equivalent(complex_model, expected_results_fail)

    @pytest.mark.unit
    def test_on_indexed_block(self, complex_model):
        """
        Tests asserting values on one instance of an indexed block ('indexed_b[X]').
        """
        target_block = complex_model.indexed_b["X"]
        expected_results = {
            "v_indexed": {None: (30, 1e-6, None)},
        }
        # This should pass
        assert_solution_equivalent(target_block, expected_results)

        # Test failure on a different index
        target_block_Y = complex_model.indexed_b["Y"]
        expected_results_fail = {
            "v_indexed": {None: (30, 1e-6, None)},  # Should be 31
        }
        with pytest.raises(AssertionError, match="Found 1 mismatch"):
            assert_solution_equivalent(target_block_Y, expected_results_fail)

        # Now use full path
        expected_results = {
            "indexed_b[X].v_indexed": {None: (30, 1e-6, None)},
        }
        assert_solution_equivalent(complex_model, expected_results)

        # Test failure on a different index
        expected_results_fail = {
            "indexed_b[Y].v_indexed": {None: (30, 1e-6, None)},  # Should be 31
        }
        with pytest.raises(AssertionError, match="Found 1 mismatch"):
            assert_solution_equivalent(complex_model, expected_results_fail)

    @pytest.mark.unit
    def test_on_nested_simple_in_indexed_block(self, complex_model):
        """
        Tests asserting values on a simple block that is nested inside one
        instance of an indexed block ('indexed_b[Y].simple_b_nested').
        """
        target_block = complex_model.indexed_b["Y"].simple_b_nested
        expected_results = {
            "v_nested_simple": {None: (41, 1e-6, None)},
        }
        # This should pass
        assert_solution_equivalent(target_block, expected_results)

        # Test a failure on the same block
        expected_results_fail = {
            "v_nested_simple": {None: (40, 1e-6, None)},  # Wrong value
        }
        with pytest.raises(AssertionError, match="Found 1 mismatch"):
            assert_solution_equivalent(target_block, expected_results_fail)

        # Now test using full path
        expected_results = {
            "indexed_b[Y].simple_b_nested.v_nested_simple": {None: (41, 1e-6, None)},
        }
        # This should pass
        assert_solution_equivalent(complex_model, expected_results)

        # Now test failure
        expected_results_fail = {
            "indexed_b[Y].simple_b_nested.v_nested_simple": {None: (40, 1e-6, None)},
        }
        with pytest.raises(AssertionError, match="Found 1 mismatch"):
            assert_solution_equivalent(complex_model, expected_results_fail)

    @pytest.mark.unit
    def test_on_nested_indexed_in_simple_block(self, complex_model):
        """
        Tests asserting values on an indexed block nested within a simple block.
        The function must be called on the top-level block containing the component.
        We check a variable within one instance of the nested indexed block.
        """
        target_block = complex_model.simple_b
        # The name must specify the index of the component to find.
        # We are looking for component 'indexed_b_nested[X]' on 'simple_b'.
        # Then we check the variable 'v_nested_indexed' on that component.
        # The utility handles this by first finding 'indexed_b_nested' and then
        # applying the index to it.
        expected_results = {
            "indexed_b_nested": {
                "X": (
                    20,
                    1e-6,
                    None,
                ),  # Special case: value of a VarData on a BlockData
                "Y": (21, 1e-6, None),
            }
        }
        # This test is tricky. `value(m.simple_b.indexed_b_nested['X'])` is not what we want.
        # The component we are checking is `m.simple_b.indexed_b_nested['X'].v_nested_indexed`.
        # Your current utility is designed to find a component by *name* and then apply an index.
        # It finds 'indexed_b_nested' and then accesses index 'X'. This yields the BlockData
        # m.simple_b.indexed_b_nested['X']. Your code then tries to call value() on this,
        # which is not a valid operation for a block.

        # To make this work, we must check the variable itself.
        # Let's see how the utility handles a fully-qualified component name.
        # Pyomo's find_component can handle dot notation.
        expected_results_correct = {
            "indexed_b_nested[X].v_nested_indexed": {None: (20, 1e-6, None)},
            "indexed_b_nested[Y].v_nested_indexed": {None: (21, 1e-6, None)},
        }
        assert_solution_equivalent(target_block, expected_results_correct)

        # And a failure case
        expected_results_fail = {
            "indexed_b_nested[X].v_nested_indexed": {None: (99, 1e-6, None)},
        }
        with pytest.raises(AssertionError, match="Found 1 mismatch"):
            assert_solution_equivalent(target_block, expected_results_fail)
