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
Tests for scaling utility functions.

Author: Andrew Lee, Douglas Allan
"""

import re

import pytest

from pyomo.environ import (
    assert_optimal_termination,
    Constraint,
    ConcreteModel,
    SolverFactory,
    Suffix,
    Var,
)
from pyomo.contrib.pynumero.asl import AmplInterface
from pyomo.contrib.pynumero.interfaces.external_grey_box import ExternalGreyBoxBlock
import pyomo.contrib.pynumero.interfaces.tests.external_grey_box_models as ex_models

from idaes.core.scaling.util import (
    get_jacobian,
    _get_jacobian_greybox_compatible,
    jacobian_cond,
    set_scaling_factor,
)
import idaes.logger as idaeslog


@pytest.mark.component
def test_get_jacobian_greybox_error_with_egb_constraints():
    m = ConcreteModel()
    m.egb = ExternalGreyBoxBlock()
    m.egb.set_external_model(
        ex_models.PressureDropSingleOutput(),
        build_implicit_constraint_objects=True,
    )

    # Add linking variables and constraints
    m.Pin = Var(initialize=1e5)
    m.c = Var(initialize=1e3)
    m.F = Var(initialize=0.5)
    m.Pout = Var(initialize=1e5)

    m.link_P_in = Constraint(expr=m.Pin == m.egb.inputs["Pin"])
    m.link_c = Constraint(expr=m.c == m.egb.inputs["c"])
    m.link_F = Constraint(expr=m.F == m.egb.inputs["F"])
    m.link_P_out = Constraint(expr=m.Pout == m.egb.outputs["Pout"])

    m.Pin.fix()
    m.c.fix()
    m.F.fix()

    # Get Jacobian
    with pytest.raises(
        ValueError,
        match="The model contains components that are not supported by the Pyomo NL writer. "
        "This may be because the model contains a grey-box model. If you want to include "
        "grey-box variables and constraints in the Jacobian, set "
        "include_greybox=True when calling get_jacobian.",
    ):
        get_jacobian(m, include_greybox=False)


# Parity tests to ensure grey box compatible Jacobian matches legacy tools
# when no grey box is present
# Tests based on John Eslick's originals for the old scaling tools
@pytest.mark.skipif(
    not AmplInterface.available(), reason="pynumero_ASL is not available"
)
class TestJacobianMethodsNoGreyBox:
    @pytest.fixture
    def model(self):
        m = ConcreteModel()
        x = m.x = Var(initialize=1e3)
        y = m.y = Var(initialize=1e6)
        z = m.z = Var(initialize=1e4)
        m.c1 = Constraint(expr=0 == -x * y + z)
        m.c2 = Constraint(expr=0 == 3 * x + 4 * y + 2 * z)
        m.c3 = Constraint(expr=0 <= z**3)
        return m

    @pytest.mark.unit
    def test_jacobian(self, model):
        """Make sure the Jacobian from Pynumero matches expectation.  This is
        mostly to ensure we understand the interface and catch if things change.
        """
        m = model
        jac, nlp = _get_jacobian_greybox_compatible(m)
        assert not hasattr(m, "scaling_factor")
        assert not hasattr(m, "scaling_hint")

        c1_row = nlp.constraint_names().index("c1")
        c2_row = nlp.constraint_names().index("c2")
        c3_row = nlp.constraint_names().index("c3")
        x_col = nlp.primals_names().index("x")
        y_col = nlp.primals_names().index("y")
        z_col = nlp.primals_names().index("z")

        assert jac[c1_row, x_col] == pytest.approx(-1e6)
        assert jac[c1_row, y_col] == pytest.approx(-1e3)
        assert jac[c1_row, z_col] == pytest.approx(1)

        assert jac[c2_row, x_col] == pytest.approx(3)
        assert jac[c2_row, y_col] == pytest.approx(4)
        assert jac[c2_row, z_col] == pytest.approx(2)

        assert jac[c3_row, z_col] == pytest.approx(3e8)

        # Make sure scaling factors don't affect the result
        set_scaling_factor(m.c1, 1e-6)
        set_scaling_factor(m.x, 1e-3)
        set_scaling_factor(m.y, 1e-6)
        set_scaling_factor(m.z, 1e-4)
        jac, _ = _get_jacobian_greybox_compatible(m, include_scaling_factors=False)
        assert len(m.scaling_factor) == 4
        assert not hasattr(m, "scaling_hint")
        assert jac[c1_row, x_col] == pytest.approx(-1e6)

        # Check the scaled jacobian calculation
        jac_scaled, _ = _get_jacobian_greybox_compatible(m)
        assert len(m.scaling_factor) == 4
        assert not hasattr(m, "scaling_hint")
        assert jac_scaled[c1_row, x_col] == pytest.approx(-1000)
        assert jac_scaled[c1_row, y_col] == pytest.approx(-1000)
        assert jac_scaled[c1_row, z_col] == pytest.approx(0.01)

    @pytest.mark.unit
    def test_scale_no_var_scale(self, model):
        m = model
        jac_scaled, nlp = _get_jacobian_greybox_compatible(
            m,
            include_scaling_factors=False,
            include_ipopt_autoscaling=True,
            min_scale=1e-6,
        )
        # get_scaling_factor isn't called here so the suffix shouldn't exist
        assert not hasattr(m, "scaling_factor")
        assert not hasattr(m, "scaling_hint")

        c1_row = nlp.constraint_names().index("c1")
        c2_row = nlp.constraint_names().index("c2")
        c3_row = nlp.constraint_names().index("c3")
        x_col = nlp.primals_names().index("x")
        y_col = nlp.primals_names().index("y")
        z_col = nlp.primals_names().index("z")

        assert jac_scaled[c1_row, x_col] == pytest.approx(-100)
        assert jac_scaled[c1_row, y_col] == pytest.approx(-0.1)
        assert jac_scaled[c1_row, z_col] == pytest.approx(1e-4)

        assert jac_scaled[c2_row, x_col] == pytest.approx(3)
        assert jac_scaled[c2_row, y_col] == pytest.approx(4)
        assert jac_scaled[c2_row, z_col] == pytest.approx(2)

        assert jac_scaled[c3_row, z_col] == pytest.approx(3e2)

    @pytest.mark.unit
    def test_scale_with_var_scale(self, model):
        m = model
        m.scaling_factor = Suffix(direction=Suffix.EXPORT)
        m.scaling_factor[m.x] = 1e-3
        m.scaling_factor[m.y] = 1e-6
        m.scaling_factor[m.z] = 1e-4

        # In these tests derived from the old scaling tools, we use
        # min_scale=1e-6 because that's what the old tools use as a
        # default. However, we're using a new default of 1e-8
        # because that's what appears to be IPOPT's actual default.
        jac_scaled, nlp = _get_jacobian_greybox_compatible(
            m,
            include_scaling_factors=True,
            include_ipopt_autoscaling=True,
            min_scale=1e-6,
        )
        assert len(m.scaling_factor) == 3
        assert not hasattr(m, "scaling_hint")

        c1_row = nlp.constraint_names().index("c1")
        c2_row = nlp.constraint_names().index("c2")
        c3_row = nlp.constraint_names().index("c3")
        x_col = nlp.primals_names().index("x")
        y_col = nlp.primals_names().index("y")
        z_col = nlp.primals_names().index("z")

        assert jac_scaled[c1_row, x_col] == pytest.approx(-1000)
        assert jac_scaled[c1_row, y_col] == pytest.approx(-1000)
        assert jac_scaled[c1_row, z_col] == pytest.approx(1e-2)

        assert jac_scaled[c2_row, x_col] == pytest.approx(0.075)
        assert jac_scaled[c2_row, y_col] == pytest.approx(100)
        assert jac_scaled[c2_row, z_col] == pytest.approx(0.5)

        assert jac_scaled[c3_row, z_col] == pytest.approx(3e6)

    @pytest.mark.unit
    def test_exclude_scaling_factors_variables(self, model):
        m = model
        m.scaling_factor = Suffix(direction=Suffix.EXPORT)
        m.scaling_factor[m.x] = 1e-3
        m.scaling_factor[m.y] = 1e-6
        m.scaling_factor[m.z] = 1e-4

        jac_scaled, nlp = _get_jacobian_greybox_compatible(
            m,
            include_scaling_factors=False,
            include_ipopt_autoscaling=True,
            min_scale=1e-6,
        )
        assert len(m.scaling_factor) == 3
        assert not hasattr(m, "scaling_hint")

        c1_row = nlp.constraint_names().index("c1")
        c2_row = nlp.constraint_names().index("c2")
        c3_row = nlp.constraint_names().index("c3")
        x_col = nlp.primals_names().index("x")
        y_col = nlp.primals_names().index("y")
        z_col = nlp.primals_names().index("z")

        assert jac_scaled[c1_row, x_col] == pytest.approx(-100)
        assert jac_scaled[c1_row, y_col] == pytest.approx(-0.1)
        assert jac_scaled[c1_row, z_col] == pytest.approx(1e-4)

        assert jac_scaled[c2_row, x_col] == pytest.approx(3)
        assert jac_scaled[c2_row, y_col] == pytest.approx(4)
        assert jac_scaled[c2_row, z_col] == pytest.approx(2)

        assert jac_scaled[c3_row, z_col] == pytest.approx(3e2)

    @pytest.mark.unit
    def test_condition_number(self, model, caplog):
        """Calculate the condition number of the Jacobian"""
        m = model
        m.scaling_factor = Suffix(direction=Suffix.EXPORT)
        m.scaling_factor[m.x] = 1e-3
        m.scaling_factor[m.y] = 1e-6
        m.scaling_factor[m.z] = 1e-4
        m.scaling_factor[m.c1] = 1e-6
        m.scaling_factor[m.c2] = 1e-6
        m.scaling_factor[m.c3] = 1e-12

        jac, _ = _get_jacobian_greybox_compatible(
            m,
            include_scaling_factors=False,
            include_ipopt_autoscaling=False,
        )
        jac_scaled, _ = _get_jacobian_greybox_compatible(
            m,
            include_scaling_factors=True,
            include_ipopt_autoscaling=False,
            min_scale=1e-6,
        )

        n = jacobian_cond(m, jac=jac_scaled)
        assert n == pytest.approx(687.47, rel=1e-3)
        n = jacobian_cond(m, jac=jac)
        assert n == pytest.approx(7.50567e7, rel=1e-3)

        # Nonsquare condition number
        m.c3.deactivate()
        jac, _ = _get_jacobian_greybox_compatible(
            m,
            include_scaling_factors=False,
            include_ipopt_autoscaling=False,
        )
        jac_scaled, _ = _get_jacobian_greybox_compatible(
            m,
            include_scaling_factors=True,
            include_ipopt_autoscaling=False,
            min_scale=1e-6,
        )

        # Scaled
        with caplog.at_level(idaeslog.INFO):
            n = jacobian_cond(m, jac=jac_scaled)
        assert (
            "Nonsquare Jacobian. Using pseudoinverse to calculate Frobenius norm."
        ) in caplog.text
        assert n == pytest.approx(500.367, rel=1e-3)

        # Unscaled
        with caplog.at_level(idaeslog.INFO):
            n = jacobian_cond(m, jac=jac)
        assert (
            "Nonsquare Jacobian. Using pseudoinverse to calculate Frobenius norm."
        ) in caplog.text
        assert n == pytest.approx(2.23741e5, rel=1e-3)


@pytest.mark.skipif(
    not AmplInterface.available(), reason="pynumero_ASL is not available"
)
class TestJacobianMethodsWithGreyBox:
    @pytest.fixture
    def model(self):
        m = ConcreteModel()
        m.egb = ExternalGreyBoxBlock()
        m.egb.set_external_model(
            ex_models.PressureDropSingleOutput(),
            build_implicit_constraint_objects=True,
        )

        # Add linking variables and constraints
        m.Pin = Var(initialize=1e5)
        m.c = Var(initialize=1e3)
        m.F = Var(initialize=0.5)
        m.Pout = Var(initialize=1e5)

        m.link_P_in = Constraint(expr=m.Pin == m.egb.inputs["Pin"])
        m.link_c = Constraint(expr=m.c == m.egb.inputs["c"])
        m.link_F = Constraint(expr=m.F == m.egb.inputs["F"])
        m.link_P_out = Constraint(expr=m.Pout == m.egb.outputs["Pout"])

        m.Pin.fix()
        m.c.fix()
        m.F.fix()

        solver = SolverFactory("cyipopt")
        results = solver.solve(m, tee=False)
        assert_optimal_termination(results)

        return m

    @pytest.mark.unit
    def test_jacobian(self, model):
        """Make sure the Jacobian from Pynumero matches expectation.  This is
        mostly to ensure we understand the interface and catch if things change.
        """
        m = model
        jac, nlp = _get_jacobian_greybox_compatible(m, include_scaling_factors=False)
        assert not hasattr(m, "scaling_factor")
        assert not hasattr(m, "scaling_hint")

        expected_jac = {
            ("link_P_in", "egb.inputs[Pin]"): -1.0,
            ("link_c", "egb.inputs[c]"): -1.0,
            ("link_F", "egb.inputs[F]"): -1.0,
            ("link_P_out", "Pout"): 1.0,
            ("link_P_out", "egb.outputs[Pout]"): -1.0,
            ("egb.Pout_constraint", "egb.inputs[Pin]"): 1.0,
            ("egb.Pout_constraint", "egb.inputs[c]"): -1.0,  # -4*0.5**2
            ("egb.Pout_constraint", "egb.inputs[F]"): -4000.0,  # -4*1e5*2*0.5
            ("egb.Pout_constraint", "egb.outputs[Pout]"): -1.0,
        }

        jac_coo = jac.tocoo()
        assert len(jac_coo.data) == 9
        for i, j, val in zip(jac_coo.row, jac_coo.col, jac_coo.data):
            if (nlp.constraint_names()[i], nlp.primals_names()[j]) not in expected_jac:
                assert val == 0
            else:
                assert val == pytest.approx(
                    expected_jac[nlp.constraint_names()[i], nlp.primals_names()[j]]
                )

        # Make sure scaling factors don't affect the result
        set_scaling_factor(m.link_P_in, 1e-5)
        set_scaling_factor(m.link_c, 1e-3)
        set_scaling_factor(m.link_F, 10)
        set_scaling_factor(m.link_P_out, 1e-5)
        set_scaling_factor(m.Pin, 1e-5)
        set_scaling_factor(m.c, 1e-3)
        set_scaling_factor(m.F, 10)
        set_scaling_factor(m.Pout, 1e-5)
        jac, _ = _get_jacobian_greybox_compatible(m, include_scaling_factors=False)
        assert len(m.scaling_factor) == 8
        assert not hasattr(m, "scaling_hint")

        # Check the scaled jacobian calculation
        jac_coo = jac.tocoo()
        assert len(jac_coo.data) == 9
        for i, j, val in zip(jac_coo.row, jac_coo.col, jac_coo.data):
            if (nlp.constraint_names()[i], nlp.primals_names()[j]) not in expected_jac:
                assert val == 0
            else:
                assert val == pytest.approx(
                    expected_jac[nlp.constraint_names()[i], nlp.primals_names()[j]]
                )

    @pytest.mark.unit
    def test_jacobian_equality_only(self, model):
        """Make sure the equality_only behaviour works as expected."""
        m = model

        # Add an inequality
        m.ineq = Constraint(expr=m.Pout <= 1e6)

        jac, nlp = _get_jacobian_greybox_compatible(m, equality_constraints_only=False)
        assert not hasattr(m, "scaling_factor")
        assert not hasattr(m, "scaling_hint")

        expected_jac = {
            ("link_P_in", "egb.inputs[Pin]"): -1.0,
            ("link_c", "egb.inputs[c]"): -1.0,
            ("link_F", "egb.inputs[F]"): -1.0,
            ("link_P_out", "Pout"): 1.0,
            ("link_P_out", "egb.outputs[Pout]"): -1.0,
            ("egb.Pout_constraint", "egb.inputs[Pin]"): 1.0,
            ("egb.Pout_constraint", "egb.inputs[c]"): -1.0,  # -4*0.5**2
            ("egb.Pout_constraint", "egb.inputs[F]"): -4000.0,  # -4*1e5*2*0.5
            ("egb.Pout_constraint", "egb.outputs[Pout]"): -1.0,
            ("ineq", "Pout"): 1.0,
        }

        jac_coo = jac.tocoo()
        assert len(jac_coo.data) == 10
        for i, j, val in zip(jac_coo.row, jac_coo.col, jac_coo.data):
            if (nlp.constraint_names()[i], nlp.primals_names()[j]) not in expected_jac:
                assert val == 0
            else:
                assert val == pytest.approx(
                    expected_jac[nlp.constraint_names()[i], nlp.primals_names()[j]]
                )

        # Test with equality_only=True to make sure the inequality is excluded
        jac, nlp = _get_jacobian_greybox_compatible(m, equality_constraints_only=True)
        jac_coo = jac.tocoo()
        assert len(jac_coo.data) == 9

        # Create a new mapping for constraint names that leaves out `ineq`
        eq_constraints = [c for c in nlp.constraint_names() if not c.startswith("ineq")]
        for i, j, val in zip(jac_coo.row, jac_coo.col, jac_coo.data):
            if (eq_constraints[i], nlp.primals_names()[j]) not in expected_jac:
                assert val == 0
            else:
                assert val == pytest.approx(
                    expected_jac[eq_constraints[i], nlp.primals_names()[j]]
                )

    @pytest.mark.unit
    def test_jacobian_w_ipopt_scaling(self, model):
        """Make sure the IPOPT scaling behaviour works as expected."""
        m = model
        jac, nlp = _get_jacobian_greybox_compatible(m, include_ipopt_autoscaling=True)
        assert not hasattr(m, "scaling_factor")
        assert not hasattr(m, "scaling_hint")

        expected_jac = {
            ("link_P_in", "egb.inputs[Pin]"): -1.0,
            ("link_c", "egb.inputs[c]"): -1.0,
            ("link_F", "egb.inputs[F]"): -1.0,
            ("link_P_out", "Pout"): 1.0,
            ("link_P_out", "egb.outputs[Pout]"): -1.0,
            # All EGB constraints get multiplied by 0.025
            ("egb.Pout_constraint", "egb.inputs[Pin]"): 0.025,
            ("egb.Pout_constraint", "egb.inputs[c]"): -0.025,
            ("egb.Pout_constraint", "egb.inputs[F]"): -100.0,
            ("egb.Pout_constraint", "egb.outputs[Pout]"): -0.025,
        }

        jac_coo = jac.tocoo()
        assert len(jac_coo.data) == 9
        for i, j, val in zip(jac_coo.row, jac_coo.col, jac_coo.data):
            if (nlp.constraint_names()[i], nlp.primals_names()[j]) not in expected_jac:
                assert val == 0
            else:
                assert val == pytest.approx(
                    expected_jac[nlp.constraint_names()[i], nlp.primals_names()[j]]
                )

    @pytest.mark.unit
    def test_condition_number(self, model):
        """Calculate the condition number of the Jacobian"""
        m = model
        m.scaling_factor = Suffix(direction=Suffix.EXPORT)
        m.scaling_factor[m.link_P_in] = 1e-5
        m.scaling_factor[m.link_c] = 1e-3
        m.scaling_factor[m.link_F] = 10
        m.scaling_factor[m.link_P_out] = 1e-5
        m.scaling_factor[m.Pin] = 1e-5
        m.scaling_factor[m.c] = 1e-3
        m.scaling_factor[m.F] = 10
        m.scaling_factor[m.Pout] = 1e-5

        jac, _ = _get_jacobian_greybox_compatible(
            m,
            include_scaling_factors=False,
            include_ipopt_autoscaling=False,
        )
        n = jacobian_cond(m, jac=jac)
        assert n == pytest.approx(2.26274262e7, rel=1e-8)

        jac_scaled, _ = _get_jacobian_greybox_compatible(
            m,
            include_scaling_factors=True,
            include_ipopt_autoscaling=False,
            min_scale=1e-6,
        )
        n = jacobian_cond(m, jac=jac_scaled)
        assert n == pytest.approx(5.657178098e8, rel=1e-3)
