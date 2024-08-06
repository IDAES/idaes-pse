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
Tests for ScalerBase.

Author: Andrew Lee
"""
import pytest
import re

from pyomo.environ import ConcreteModel, Constraint, Set, Suffix, Var
from pyomo.common.config import ConfigDict

from idaes.core.scaling.scaling_base import ScalerBase
import idaes.logger as idaeslog


@pytest.fixture
def model():
    m = ConcreteModel()

    m.s = Set(initialize=[1, 2, 3, 4])

    m.v = Var(m.s)

    def c_rule(b, i):
        return b.v[i] == i

    m.c = Constraint(m.s, rule=c_rule)

    m.scaling_factor = Suffix(direction=Suffix.EXPORT)

    m.scaling_factor[m.v[1]] = 1
    m.scaling_factor[m.v[2]] = 2
    m.scaling_factor[m.v[3]] = 3
    m.scaling_factor[m.v[4]] = 4
    m.scaling_factor[m.c[1]] = 11
    m.scaling_factor[m.c[2]] = 21
    m.scaling_factor[m.c[3]] = 31
    m.scaling_factor[m.c[4]] = 41

    return m


class TestScalerBase:
    @pytest.mark.unit
    def test_init(self):
        sb = ScalerBase()

        assert isinstance(sb.config, ConfigDict)
        assert sb.config.zero_tolerance == 1e-12
        assert sb.config.max_variable_scaling_factor == 1e10
        assert sb.config.min_variable_scaling_factor == 1e-10
        assert sb.config.max_constraint_scaling_factor == 1e10
        assert sb.config.min_constraint_scaling_factor == 1e-10
        assert not sb.config.overwrite

    @pytest.mark.unit
    def test_get_scaling_factor(self, model):
        sb = ScalerBase()

        assert sb.get_scaling_factor(model.v[1]) == 1

    @pytest.mark.unit
    def test_set_scaling_factor(self, model):
        sb = ScalerBase()

        sb._set_scaling_factor(
            component=model.v[1],
            component_type="variable",
            scaling_factor=42,
            overwrite=False,
        )
        sb._set_scaling_factor(
            component=model.c[1],
            component_type="constraint",
            scaling_factor=42,
            overwrite=False,
        )
        assert model.scaling_factor[model.v[1]] == 1  # Overwrite = False, no change
        assert model.scaling_factor[model.c[1]] == 11  # Overwrite = False, no change

        sb._set_scaling_factor(
            component=model.v[1],
            component_type="variable",
            scaling_factor=42,
            overwrite=True,
        )
        sb._set_scaling_factor(
            component=model.c[1],
            component_type="constraint",
            scaling_factor=42,
            overwrite=True,
        )
        assert model.scaling_factor[model.v[1]] == 42
        assert model.scaling_factor[model.c[1]] == 42

    @pytest.mark.unit
    def test_set_scaling_factor_overwrite_default(self, model):
        sb = ScalerBase(overwrite=False)

        sb._set_scaling_factor(
            component=model.v[1], component_type="variable", scaling_factor=42
        )
        sb._set_scaling_factor(
            component=model.c[1], component_type="constraint", scaling_factor=42
        )
        assert model.scaling_factor[model.v[1]] == 1  # Overwrite = False, no change
        assert model.scaling_factor[model.c[1]] == 11  # Overwrite = False, no change

        # Change default overwrite setting
        sb.config.overwrite = True
        sb._set_scaling_factor(
            component=model.v[1], component_type="variable", scaling_factor=42
        )
        sb._set_scaling_factor(
            component=model.c[1], component_type="constraint", scaling_factor=42
        )
        assert model.scaling_factor[model.v[1]] == 42
        assert model.scaling_factor[model.c[1]] == 42

    @pytest.mark.unit
    def test_set_scaling_factor_invalid_component_type(self, model):
        sb = ScalerBase()

        with pytest.raises(
            ValueError,
            match="Invalid value for component_type.",
        ):
            sb._set_scaling_factor(
                component=model.v[1], component_type="foo", scaling_factor="bar"
            )

    @pytest.mark.unit
    def test_set_scaling_factor_max_limit(self, model, caplog):
        sb = ScalerBase(
            overwrite=True,
            max_variable_scaling_factor=100,
            max_constraint_scaling_factor=200,
        )

        caplog.set_level(
            idaeslog.DEBUG,
            logger="idaes",
        )
        # Below max var limit - value set as given
        sb._set_scaling_factor(
            component=model.v[1], component_type="variable", scaling_factor=42
        )
        assert model.scaling_factor[model.v[1]] == 42

        # Above max var limit - value limited to max
        sb._set_scaling_factor(
            component=model.v[1], component_type="variable", scaling_factor=150
        )
        assert model.scaling_factor[model.v[1]] == 100
        assert (
            "Scaling factor for v[1] limited by maximum value (max_sf: 100.0 < sf: 150)"
            in caplog.text
        )

        # Above max var limit but below max con limit - value set as given
        sb._set_scaling_factor(
            component=model.v[1], component_type="constraint", scaling_factor=150
        )
        assert model.scaling_factor[model.v[1]] == 150

        # Above max con limit - value limited to max
        sb._set_scaling_factor(
            component=model.v[1], component_type="constraint", scaling_factor=250
        )
        assert model.scaling_factor[model.v[1]] == 200
        assert (
            "Scaling factor for v[1] limited by maximum value (max_sf: 200.0 < sf: 250)"
            in caplog.text
        )

    @pytest.mark.unit
    def test_set_scaling_factor_min_limit(self, model, caplog):
        sb = ScalerBase(
            overwrite=True,
            min_variable_scaling_factor=100,
            min_constraint_scaling_factor=200,
        )

        caplog.set_level(idaeslog.DEBUG, logger="idaes")
        # Above both min limits - value set as given
        sb._set_scaling_factor(
            component=model.v[1], component_type="constraint", scaling_factor=400
        )
        assert model.scaling_factor[model.v[1]] == 400

        # Below con min limit - value limited
        sb._set_scaling_factor(
            component=model.v[1], component_type="constraint", scaling_factor=150
        )
        assert model.scaling_factor[model.v[1]] == 200
        assert (
            "Scaling factor for v[1] limited by minimum value (min_sf: 200.0 > sf: 150)"
            in caplog.text
        )

        # Below con min limit but above min var limit - value set as given
        sb._set_scaling_factor(
            component=model.v[1], component_type="variable", scaling_factor=150
        )
        assert model.scaling_factor[model.v[1]] == 150

        # Below both limited
        sb._set_scaling_factor(
            component=model.v[1], component_type="variable", scaling_factor=50
        )
        assert model.scaling_factor[model.v[1]] == 100
        assert (
            "Scaling factor for v[1] limited by minimum value (min_sf: 100.0 > sf: 50)"
            in caplog.text
        )

    @pytest.mark.unit
    def test_set_variable_scaling_factor(self, model, caplog):
        sb = ScalerBase(
            max_variable_scaling_factor=1e3,
            min_variable_scaling_factor=100,
            max_constraint_scaling_factor=2e3,
            min_constraint_scaling_factor=200,
            overwrite=True,
        )

        caplog.set_level(idaeslog.DEBUG, logger="idaes")

        # Scaling factor within limits
        sb.set_variable_scaling_factor(model.v[1], 200)
        assert model.scaling_factor[model.v[1]] == 200

        # Too large
        sb.set_variable_scaling_factor(model.v[1], 2e3)
        assert model.scaling_factor[model.v[1]] == 1e3
        assert (
            "Scaling factor for v[1] limited by maximum value (max_sf: 1000.0 < sf: 2000.0)"
            in caplog.text
        )

        # Too small
        sb.set_variable_scaling_factor(model.v[1], 1)
        assert model.scaling_factor[model.v[1]] == 100
        assert (
            "Scaling factor for v[1] limited by minimum value (min_sf: 100.0 > sf: 1)"
            in caplog.text
        )

    @pytest.mark.unit
    def test_set_variable_scaling_factor_invalid_type(self, model):
        sb = ScalerBase()

        with pytest.raises(
            TypeError, match=re.escape("c[1] is not a variable (or is indexed).")
        ):
            sb.set_variable_scaling_factor(model.c[1], 200)

        with pytest.raises(
            TypeError, match=re.escape("v is not a variable (or is indexed).")
        ):
            sb.set_variable_scaling_factor(model.v, 200)

    @pytest.mark.unit
    def test_set_constraint_scaling_factor(self, model, caplog):
        sb = ScalerBase(
            max_variable_scaling_factor=1e3,
            min_variable_scaling_factor=100,
            max_constraint_scaling_factor=2e3,
            min_constraint_scaling_factor=200,
            overwrite=True,
        )

        caplog.set_level(idaeslog.DEBUG, logger="idaes")

        # Scaling factor within limits
        sb.set_constraint_scaling_factor(model.c[1], 200)
        assert model.scaling_factor[model.c[1]] == 200

        # Too large
        sb.set_constraint_scaling_factor(model.c[1], 3e3)
        assert model.scaling_factor[model.c[1]] == 2e3
        assert (
            "Scaling factor for c[1] limited by maximum value (max_sf: 2000.0 < sf: 3000.0)"
            in caplog.text
        )

        # Too small
        sb.set_constraint_scaling_factor(model.c[1], 1)
        assert model.scaling_factor[model.c[1]] == 200
        assert (
            "Scaling factor for c[1] limited by minimum value (min_sf: 200.0 > sf: 1)"
            in caplog.text
        )

    @pytest.mark.unit
    def test_set_constraint_scaling_factor_invalid_type(self, model):
        sb = ScalerBase()

        with pytest.raises(
            TypeError, match=re.escape("v[1] is not a constraint (or is indexed).")
        ):
            sb.set_constraint_scaling_factor(model.v[1], 200)

        with pytest.raises(
            TypeError, match=re.escape("c is not a constraint (or is indexed).")
        ):
            sb.set_constraint_scaling_factor(model.c, 200)
