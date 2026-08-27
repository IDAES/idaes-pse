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

The tests in this file were generated with the assistance of Google Gemini 3.6 Flash

Author: Douglas Allan
"""

from unittest.mock import patch
import pytest
from pyomo.environ import Block, ConcreteModel, Constraint, TransformationFactory, Var
from pyomo.network import Arc, Port
from idaes.core.scaling.arc_constraint_scaler import ArcConstraintScaler
from idaes.core.scaling.custom_scaler_base import ConstraintScalingScheme
from idaes.core.scaling.util import set_scaling_factor


class TestArcConstraintScaler:
    @pytest.fixture
    def arc_model(self):
        m = ConcreteModel()

        # Block 1 setup
        m.s1 = Block()
        m.s1.x = Var(initialize=10.0)
        m.s1.flow = Var(["a", "b"], initialize={"a": 2.0, "b": 3.0})
        m.s1.port = Port()
        m.s1.port.add(m.s1.x, "x")
        m.s1.port.add(m.s1.flow, "flow", rule=Port.Extensive)

        # Block 2 setup
        m.s2 = Block()
        m.s2.x = Var(initialize=5.0)
        m.s2.flow = Var(["a", "b"], initialize={"a": 1.0, "b": 1.5})
        m.s2.port = Port()
        m.s2.port.add(m.s2.x, "x")
        m.s2.port.add(m.s2.flow, "flow", rule=Port.Extensive)

        # Connect ports via Arc
        m.arc = Arc(src=m.s1.port, dest=m.s2.port)

        # Expand arc to generate expanded_block constraints
        TransformationFactory("network.expand_arcs").apply_to(m)

        # Set variable scaling factors
        set_scaling_factor(m.s1.x, 0.1)  # nominal value = 10.0
        set_scaling_factor(m.s2.x, 0.2)  # nominal value = 5.0

        set_scaling_factor(m.s1.flow["a"], 0.5)  # nominal value = 2.0
        set_scaling_factor(m.s2.flow["a"], 0.4)  # nominal value = 2.5

        set_scaling_factor(m.s1.flow["b"], 0.25)  # nominal value = 4.0
        set_scaling_factor(m.s2.flow["b"], 0.1)  # nominal value = 10.0

        return m

    @pytest.mark.unit
    def test_scale_model_raises_not_implemented(self):
        scaler = ArcConstraintScaler()
        with pytest.raises(
            NotImplementedError,
            match="The scale_model method is not implemented for the ArcConstraintScaler.",
        ):
            scaler.scale_model()

    @pytest.mark.unit
    def test_scale_arc_constraints_by_nominal_value(self, arc_model):
        scaler = ArcConstraintScaler()

        scheme = ConstraintScalingScheme.inverseMaximum
        with patch.object(
            scaler,
            "scale_constraint_by_nominal_value",
            wraps=scaler.scale_constraint_by_nominal_value,
        ) as spy_scale_con:
            scaler.scale_arc_constraints_by_nominal_value(
                arc_model, scheme=scheme, overwrite=True, descend_into=True
            )

            constraints = list(
                arc_model.arc.expanded_block.component_data_objects(ctype=Constraint)
            )
            assert len(constraints) == 3

            # Verify that each constraint was passed to scale_constraint_by_nominal_value with correct kwargs
            for con in constraints:
                spy_scale_con.assert_any_call(con, scheme=scheme, overwrite=True)

        # Expected inverseMaximum scaling factors (1 / max(|nom_src|, |nom_dest|)):
        # x con: max(10, 5) = 10 -> sf = 0.1
        # flow[a] con: max(2, 2.5) = 2.5 -> sf = 0.4
        # flow[b] con: max(4, 10) = 10 -> sf = 0.1
        expected_sfs = {
            "x": 0.1,
            "flow[a]": 0.4,
            "flow[b]": 0.1,
        }

        for con in constraints:
            var_name = con.local_name.replace("_equality", "")
            expected_sf = expected_sfs[var_name]

            sf = scaler.get_scaling_factor(con)
            assert sf == pytest.approx(expected_sf)

            # Verify overwrite behavior
            set_scaling_factor(con, 999.0, overwrite=True)
            scaler.scale_arc_constraints_by_nominal_value(arc_model, overwrite=False)
            assert scaler.get_scaling_factor(con) == 999.0

            scaler.scale_arc_constraints_by_nominal_value(arc_model, overwrite=True)
            assert scaler.get_scaling_factor(con) == pytest.approx(expected_sf)

    @pytest.mark.unit
    def test_scale_arc_constraints_by_nominal_derivative_norm(self, arc_model):
        scaler = ArcConstraintScaler()

        with patch.object(
            scaler,
            "scale_constraint_by_nominal_derivative_norm",
            wraps=scaler.scale_constraint_by_nominal_derivative_norm,
        ) as spy_scale_con:
            scaler.scale_arc_constraints_by_nominal_derivative_norm(
                arc_model, norm=1, overwrite=True, descend_into=True
            )

            constraints = list(
                arc_model.arc.expanded_block.component_data_objects(ctype=Constraint)
            )
            assert len(constraints) == 3

            # Verify that each constraint was passed to scale_constraint_by_nominal_derivative_norm with correct kwargs
            for con in constraints:
                spy_scale_con.assert_any_call(con, norm=1, overwrite=True)

        # Re-run with default norm=2 for value assertion checks
        scaler.scale_arc_constraints_by_nominal_derivative_norm(
            arc_model, norm=2, overwrite=True
        )

        # Expected 2-norm scaling factors (1 / sqrt((1/sf_src)^2 + (-1/sf_dest)^2)):
        # x con: 1 / sqrt(10^2 + 5^2) = 1 / sqrt(125)
        # flow[a] con: 1 / sqrt(2^2 + 2.5^2) = 1 / sqrt(10.25)
        # flow[b] con: 1 / sqrt(4^2 + 10^2) = 1 / sqrt(116)
        expected_sfs = {
            "x": 1.0 / (10.0**2 + 5.0**2) ** 0.5,
            "flow[a]": 1.0 / (2.0**2 + 2.5**2) ** 0.5,
            "flow[b]": 1.0 / (4.0**2 + 10.0**2) ** 0.5,
        }

        for con in constraints:
            var_name = con.local_name.replace("_equality", "")
            expected_sf = expected_sfs[var_name]

            sf = scaler.get_scaling_factor(con)
            assert sf == pytest.approx(expected_sf)

            # Verify overwrite behavior
            set_scaling_factor(con, 999.0, overwrite=True)
            scaler.scale_arc_constraints_by_nominal_derivative_norm(
                arc_model, norm=2, overwrite=False
            )
            assert scaler.get_scaling_factor(con) == 999.0

            scaler.scale_arc_constraints_by_nominal_derivative_norm(
                arc_model, norm=2, overwrite=True
            )
            assert scaler.get_scaling_factor(con) == pytest.approx(expected_sf)

    @pytest.mark.unit
    @pytest.mark.parametrize("descend_into", [False, True])
    def test_descend_into(self, descend_into):
        m = ConcreteModel()

        # Top-level block and arc
        m.b1 = Block()
        m.b1.x = Var(initialize=10.0)
        m.b1.p = Port()
        m.b1.p.add(m.b1.x, "x")

        m.b2 = Block()
        m.b2.x = Var(initialize=5.0)
        m.b2.p = Port()
        m.b2.p.add(m.b2.x, "x")

        m.arc_top = Arc(src=m.b1.p, dest=m.b2.p)

        # Nested subblock and arc
        m.sub = Block()
        m.sub.b1 = Block()
        m.sub.b1.y = Var(initialize=4.0)
        m.sub.b1.p = Port()
        m.sub.b1.p.add(m.sub.b1.y, "y")

        m.sub.b2 = Block()
        m.sub.b2.y = Var(initialize=2.0)
        m.sub.b2.p = Port()
        m.sub.b2.p.add(m.sub.b2.y, "y")

        m.sub.arc_nested = Arc(src=m.sub.b1.p, dest=m.sub.b2.p)

        TransformationFactory("network.expand_arcs").apply_to(m)

        set_scaling_factor(m.b1.x, 0.1)
        set_scaling_factor(m.b2.x, 0.2)
        set_scaling_factor(m.sub.b1.y, 0.25)
        set_scaling_factor(m.sub.b2.y, 0.5)

        scaler = ArcConstraintScaler()

        scaler.scale_arc_constraints_by_nominal_value(m, descend_into=descend_into)

        top_con = next(
            m.arc_top.expanded_block.component_data_objects(ctype=Constraint)
        )
        nested_con = next(
            m.sub.arc_nested.expanded_block.component_data_objects(ctype=Constraint)
        )

        # Top-level arc constraint should always be scaled
        assert scaler.get_scaling_factor(top_con) is not None

        # Nested arc constraint should only be scaled if descend_into is True
        if descend_into:
            assert scaler.get_scaling_factor(nested_con) is not None
        else:
            assert scaler.get_scaling_factor(nested_con) is None
