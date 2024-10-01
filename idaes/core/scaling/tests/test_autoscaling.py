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
Tests for autoscalers.

Author: Andrew Lee
"""
from math import sqrt
import pytest

from pyomo.environ import Block, ConcreteModel, Constraint, Objective, Set, Suffix, Var

from idaes.core.scaling import AutoScaler


@pytest.fixture
def model():
    m = ConcreteModel()
    m.s = Set(initialize=[1, 2, 3, 4])

    m.b = Block(m.s)

    m.v1 = Var(initialize=2)
    m.v2 = Var(m.s, initialize=10)

    for bd in m.b.values():
        bd.v3 = Var(initialize=10)
        bd.v4 = Var(initialize=100)

    @m.Constraint(m.s)
    def c1(blk, i):
        return blk.v1**i == blk.v2[i]

    @m.Constraint(m.s)
    def c2(blk, i):
        return blk.v2[i] == blk.b[i].v3

    for bd in m.b.values():

        @bd.Constraint()
        def c3(blk):
            return blk.v3**2 == blk.v4

    m.o = Objective(expr=m.v1)

    return m


class TestAutoscaleVarMagnitude:
    @pytest.mark.unit
    def test_var_data(self, model):
        scaler = AutoScaler()
        scaler.scale_variables_by_magnitude(model.v1)
        scaler.scale_variables_by_magnitude(model.v2[1])

        assert model.scaling_factor[model.v1] == pytest.approx(1 / 2, rel=1e-8)
        assert model.scaling_factor[model.v2[1]] == pytest.approx(1 / 10, rel=1e-8)
        assert len(model.scaling_factor) == 2

    @pytest.mark.unit
    def test_var_data_no_value(self):
        model = ConcreteModel()
        model.v1 = Var()

        scaler = AutoScaler()
        scaler.scale_variables_by_magnitude(model.v1)

        assert model.scaling_factor[model.v1] == 1

    @pytest.mark.unit
    def test_indexed_var(self, model):
        scaler = AutoScaler()
        scaler.scale_variables_by_magnitude(model.v2)

        for i in model.s:
            assert model.scaling_factor[model.v2[i]] == pytest.approx(1 / 10, rel=1e-8)
        assert len(model.scaling_factor) == 4

    @pytest.mark.unit
    def test_block_data(self, model):
        scaler = AutoScaler()
        scaler.scale_variables_by_magnitude(model.b[1])

        assert model.b[1].scaling_factor[model.b[1].v3] == pytest.approx(
            1 / 10, rel=1e-8
        )
        assert model.b[1].scaling_factor[model.b[1].v4] == pytest.approx(
            1 / 100, rel=1e-8
        )
        assert len(model.b[1].scaling_factor) == 2

    @pytest.mark.unit
    def test_indexed_block(self, model):
        scaler = AutoScaler()
        scaler.scale_variables_by_magnitude(model.b)

        for bd in model.b.values():
            sfx = bd.scaling_factor
            assert sfx[bd.v3] == pytest.approx(1 / 10, rel=1e-8)
            assert sfx[bd.v4] == pytest.approx(1 / 100, rel=1e-8)
            assert len(sfx) == 2

    @pytest.mark.unit
    def test_nested_blocks_descend(self, model):
        scaler = AutoScaler()
        scaler.scale_variables_by_magnitude(model)

        assert model.scaling_factor[model.v1] == pytest.approx(1 / 2, rel=1e-8)
        for i in model.s:
            assert model.scaling_factor[model.v2[i]] == pytest.approx(1 / 10, rel=1e-8)
        assert len(model.scaling_factor) == 5

        for bd in model.b.values():
            sfx = bd.scaling_factor
            assert sfx[bd.v3] == pytest.approx(1 / 10, rel=1e-8)
            assert sfx[bd.v4] == pytest.approx(1 / 100, rel=1e-8)
            assert len(sfx) == 2

    @pytest.mark.unit
    def test_nested_blocks_no_descend(self, model):
        scaler = AutoScaler()
        scaler.scale_variables_by_magnitude(model, descend_into=False)

        assert model.scaling_factor[model.v1] == pytest.approx(1 / 2, rel=1e-8)
        for i in model.s:
            assert model.scaling_factor[model.v2[i]] == pytest.approx(1 / 10, rel=1e-8)
        assert len(model.scaling_factor) == 5

        for bd in model.b.values():
            assert not hasattr(bd, "scaling_factor")

    @pytest.mark.unit
    def test_no_overwrite(self, model):
        # Add some scaling factors to ensure they are not overwritten
        model.scaling_factor = Suffix(direction=Suffix.EXPORT)
        model.scaling_factor[model.v1] = 20

        model.b[1].scaling_factor = Suffix(direction=Suffix.EXPORT)
        model.b[1].scaling_factor[model.b[1].v3] = 20

        scaler = AutoScaler(overwrite=False)
        scaler.scale_variables_by_magnitude(model)

        assert model.scaling_factor[model.v1] == 20
        for i in model.s:
            assert model.scaling_factor[model.v2[i]] == pytest.approx(1 / 10, rel=1e-8)
        assert len(model.scaling_factor) == 5

        for bd in model.b.values():
            sfx = bd.scaling_factor
            if bd is model.b[1]:
                assert sfx[bd.v3] == 20
            else:
                assert sfx[bd.v3] == pytest.approx(1 / 10, rel=1e-8)
            assert sfx[bd.v4] == pytest.approx(1 / 100, rel=1e-8)
            assert len(sfx) == 2

    @pytest.mark.unit
    def test_max_and_min(self, model):
        scaler = AutoScaler(
            max_variable_scaling_factor=1 / 10, min_variable_scaling_factor=1 / 50
        )
        scaler.scale_variables_by_magnitude(model)

        assert model.scaling_factor[model.v1] == pytest.approx(1 / 10, rel=1e-8)
        for i in model.s:
            assert model.scaling_factor[model.v2[i]] == pytest.approx(1 / 10, rel=1e-8)
        assert len(model.scaling_factor) == 5

        for bd in model.b.values():
            sfx = bd.scaling_factor
            assert sfx[bd.v3] == pytest.approx(1 / 10, rel=1e-8)
            assert sfx[bd.v4] == pytest.approx(1 / 50, rel=1e-8)
            assert len(sfx) == 2

    @pytest.mark.unit
    def test_var_fixed(self, model):
        model.v1.fix()
        scaler = AutoScaler()
        scaler.scale_variables_by_magnitude(model.v1)
        scaler.scale_variables_by_magnitude(model.v2[1])

        assert model.scaling_factor[model.v1] == pytest.approx(1 / 2, rel=1e-8)
        assert model.scaling_factor[model.v2[1]] == pytest.approx(1 / 10, rel=1e-8)
        assert len(model.scaling_factor) == 2

    @pytest.mark.unit
    def test_not_block_or_var(self, model):
        scaler = AutoScaler()

        with pytest.raises(TypeError, match="c1 is not a block or variable."):
            scaler.scale_variables_by_magnitude(model.c1)


class TestConstraintsByNorm:
    @pytest.mark.unit
    def test_not_block_or_var(self, model):
        scaler = AutoScaler()

        with pytest.raises(TypeError, match="v1 is not a block or constraint."):
            scaler.scale_constraints_by_jacobian_norm(model.v1)

    @pytest.mark.unit
    def test_block_data_L2(self, model):
        scaler = AutoScaler()

        scaler.scale_constraints_by_jacobian_norm(model.b[1])

        assert model.b[1].scaling_factor[model.b[1].c3] == pytest.approx(
            1 / sqrt(20**2 + 1**2), rel=1e-8
        )
        assert len(model.b[1].scaling_factor) == 1

    @pytest.mark.unit
    def test_block_data_L2_block_data(self, model):
        scaler = AutoScaler()

        scaler.scale_constraints_by_jacobian_norm(model.b)

        for bd in model.b.values():
            assert bd.scaling_factor[bd.c3] == pytest.approx(
                1 / sqrt(20**2 + 1**2), rel=1e-8
            )
            assert len(bd.scaling_factor) == 1

    @pytest.mark.unit
    def test_nested_blocks_L2(self, model):
        scaler = AutoScaler()

        scaler.scale_constraints_by_jacobian_norm(model)

        for i in model.s:
            assert model.scaling_factor[model.c1[i]] == pytest.approx(
                1 / sqrt((i * 2 ** (i - 1)) ** 2 + 1**2), rel=1e-8
            )
            assert model.scaling_factor[model.c2[i]] == pytest.approx(
                1 / sqrt(2), rel=1e-8
            )
        assert len(model.scaling_factor) == 8

        for bd in model.b.values():
            assert bd.scaling_factor[bd.c3] == pytest.approx(
                1 / sqrt(20**2 + 1**2), rel=1e-8
            )
            assert len(bd.scaling_factor) == 1

    @pytest.mark.unit
    def test_nested_blocks_L2_no_descent(self, model):
        scaler = AutoScaler()

        scaler.scale_constraints_by_jacobian_norm(model, descend_into=False)

        for i in model.s:
            assert model.scaling_factor[model.c1[i]] == pytest.approx(
                1 / sqrt((i * 2 ** (i - 1)) ** 2 + 1**2), rel=1e-8
            )
            assert model.scaling_factor[model.c2[i]] == pytest.approx(
                1 / sqrt(2), rel=1e-8
            )
        assert len(model.scaling_factor) == 8

        for bd in model.b.values():
            assert not hasattr(bd, "scaling_factor")

    @pytest.mark.unit
    def test_constraint_data_L2(self, model):
        scaler = AutoScaler()

        scaler.scale_constraints_by_jacobian_norm(model.c1[1])

        assert model.scaling_factor[model.c1[1]] == pytest.approx(
            1 / sqrt((1 * 2 ** (1 - 1)) ** 2 + 1**2), rel=1e-8
        )
        assert len(model.scaling_factor) == 1

    @pytest.mark.unit
    def test_indexed_constraint_L2(self, model):
        scaler = AutoScaler()

        scaler.scale_constraints_by_jacobian_norm(model.c1)

        for i in model.s:
            assert model.scaling_factor[model.c1[i]] == pytest.approx(
                1 / sqrt((i * 2 ** (i - 1)) ** 2 + 1**2), rel=1e-8
            )
        assert len(model.scaling_factor) == 4

    @pytest.mark.unit
    def test_block_data_L1(self, model):
        scaler = AutoScaler()

        scaler.scale_constraints_by_jacobian_norm(model.b[1], norm=1)

        assert model.b[1].scaling_factor[model.b[1].c3] == pytest.approx(
            1 / 21, rel=1e-8
        )
        assert len(model.b[1].scaling_factor) == 1

    @pytest.mark.unit
    def test_nested_blocks_L1(self, model):
        scaler = AutoScaler()

        scaler.scale_constraints_by_jacobian_norm(model, norm=1)

        for i in model.s:
            assert model.scaling_factor[model.c1[i]] == pytest.approx(
                1 / ((i * 2 ** (i - 1)) + 1), rel=1e-8
            )
            assert model.scaling_factor[model.c2[i]] == pytest.approx(1 / 2, rel=1e-8)
        assert len(model.scaling_factor) == 8

        for bd in model.b.values():
            assert bd.scaling_factor[bd.c3] == pytest.approx(1 / 21, rel=1e-8)
            assert len(bd.scaling_factor) == 1

    @pytest.mark.unit
    def test_constraint_data_L1(self, model):
        scaler = AutoScaler()

        scaler.scale_constraints_by_jacobian_norm(model.c1[1], norm=1)

        assert model.scaling_factor[model.c1[1]] == pytest.approx(
            1 / ((1 * 2 ** (1 - 1)) + 1), rel=1e-8
        )
        assert len(model.scaling_factor) == 1

    @pytest.mark.unit
    def test_indexed_constraint_L1(self, model):
        scaler = AutoScaler()

        scaler.scale_constraints_by_jacobian_norm(model.c1, norm=1)

        for i in model.s:
            assert model.scaling_factor[model.c1[i]] == pytest.approx(
                1 / ((i * 2 ** (i - 1)) + 1), rel=1e-8
            )
        assert len(model.scaling_factor) == 4

    @pytest.mark.unit
    def test_nested_blocks_L1_no_overwrite(self, model):
        # Add some scaling factors to ensure they are not overwritten
        model.scaling_factor = Suffix(direction=Suffix.EXPORT)
        model.scaling_factor[model.c1[1]] = 20

        model.b[1].scaling_factor = Suffix(direction=Suffix.EXPORT)
        model.b[1].scaling_factor[model.b[1].c3] = 20

        scaler = AutoScaler()

        scaler.scale_constraints_by_jacobian_norm(model, norm=1)

        for i in model.s:
            if i == 1:
                assert model.scaling_factor[model.c1[i]] == 20
            else:
                assert model.scaling_factor[model.c1[i]] == pytest.approx(
                    1 / ((i * 2 ** (i - 1)) + 1), rel=1e-8
                )
            assert model.scaling_factor[model.c2[i]] == pytest.approx(1 / 2, rel=1e-8)
        assert len(model.scaling_factor) == 8

        for k, bd in model.b.items():
            if k == 1:
                assert bd.scaling_factor[bd.c3] == 20
            else:
                assert bd.scaling_factor[bd.c3] == pytest.approx(1 / 21, rel=1e-8)
            assert len(bd.scaling_factor) == 1


class TestAutoScaleModel:
    @pytest.mark.unit
    def test_scale_model_default(self, model):
        scaler = AutoScaler()

        scaler.scale_model(model)

        c1_sf = {
            "c1[1]": 0.0980580676,
            "c1[2]": 0.0780868809,
            "c1[3]": 0.0384615385,
            "c1[4]": 0.0154376880,
        }

        for k, v in model.scaling_factor.items():
            if str(k).startswith("v2"):
                assert v == pytest.approx(0.1, rel=1e-8)
            elif str(k).startswith("v1"):
                assert v == pytest.approx(0.5, rel=1e-8)
            elif str(k).startswith("c2"):
                assert v == pytest.approx(0.5**0.5 * 0.1, rel=1e-8)
            elif str(k).startswith("c1"):
                assert v == pytest.approx(c1_sf[str(k)], rel=1e-8)
