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
Tests for scaling utility functions.

Author: Andrew Lee, Douglas Allan
"""
from io import StringIO
import os
import pytest
import re
from copy import deepcopy

from pyomo.environ import Block, Constraint, ConcreteModel, Expression, Set, Suffix, Var
from pyomo.common.fileutils import this_file_dir
from pyomo.common.tempfiles import TempfileManager
from pyomo.contrib.pynumero.asl import AmplInterface

from idaes.core.scaling.util import (
    get_jacobian,
    jacobian_cond,
    get_scaling_factor_suffix,
    get_scaling_hint_suffix,
    get_scaling_factor,
    set_scaling_factor,
    del_scaling_factor,
    _suffix_to_dict,
    _suffix_from_dict,
    _collect_block_suffixes,
    _set_block_suffixes_from_dict,
    list_unscaled_variables,
    list_unscaled_constraints,
    scaling_factors_to_dict,
    scaling_factors_from_dict,
    scaling_factors_to_json_file,
    scaling_factors_from_json_file,
    report_scaling_factors,
    unscaled_variables_generator,
    unscaled_constraints_generator,
)
from idaes.core.util.model_statistics import number_activated_objectives
import idaes.logger as idaeslog

currdir = this_file_dir()


class TestGetScalingFactorSuffix:
    @pytest.mark.unit
    def test_get_scaling_factor_suffix_block_new(self, caplog):
        caplog.set_level(
            idaeslog.DEBUG,
            logger="idaes",
        )

        m = ConcreteModel()
        sfx = get_scaling_factor_suffix(m)

        assert "Created new scaling suffix for model" in caplog.text

        assert isinstance(m.scaling_factor, Suffix)
        assert sfx is m.scaling_factor

    @pytest.mark.unit
    def test_get_scaling_factor_suffix_indexed_component_new(self):
        m = ConcreteModel()
        m.v = Var([1, 2, 3, 4])
        with pytest.raises(
            TypeError,
        ) as einfo:
            _ = get_scaling_factor_suffix(m.v[1])
        assert (
            "Component v[1] was not a BlockData, instead it was a <class 'pyomo.core.base.var.VarData'"
            in str(einfo)
        )

    @pytest.mark.unit
    def test_get_scaling_factor_suffix_indexed_block(self):
        m = ConcreteModel()
        m.b = Block([1, 2, 3, 4])

        with pytest.raises(
            TypeError,
            match="IndexedBlocks cannot have scaling factors attached to them. "
            "Please assign scaling factors to the elements of the IndexedBlock.",
        ):
            get_scaling_factor_suffix(m.b)

    @pytest.mark.unit
    def test_get_scaling_factor_suffix_component_new(self):
        m = ConcreteModel()
        m.v = Var()
        with pytest.raises(TypeError) as e_obj:
            _ = get_scaling_factor_suffix(m.v)
        assert (
            "Component v was not a BlockData, instead it was a <class 'pyomo.core.base.var.ScalarVar'>"
            in str(e_obj)
        )

    @pytest.mark.unit
    def test_get_scaling_factor_suffix_block_existing(self, caplog):
        m = ConcreteModel()
        m.scaling_factor = Suffix(direction=Suffix.EXPORT)
        sfx = get_scaling_factor_suffix(m)

        assert "Created new scaling suffix for model" not in caplog.text

        assert isinstance(m.scaling_factor, Suffix)
        assert sfx is m.scaling_factor

    @pytest.mark.unit
    def test_get_scaling_factor_suffix_component_existing(self):
        m = ConcreteModel()
        m.scaling_factor = Suffix(direction=Suffix.EXPORT)
        m.v = Var()
        with pytest.raises(TypeError) as eobj:
            _ = get_scaling_factor_suffix(m.v)

        assert (
            "Component v was not a BlockData, instead it was a <class 'pyomo.core.base.var.ScalarVar'>"
            in str(eobj)
        )

    @pytest.mark.unit
    def test_get_scaling_factor_suffix_deactivated(self):
        m = ConcreteModel()
        m.scaling_factor = Suffix(direction=Suffix.EXPORT)
        m.scaling_factor.deactivate()
        with pytest.raises(
            RuntimeError,
            match=re.escape(
                f"The scaling suffix on model has been deactivated. "
                "Typically, this means that the user has performed a scaling transformation "
                "on the model."
            ),
        ):
            _ = get_scaling_factor_suffix(m)


class TestGetScalingHintSuffix:
    @pytest.mark.unit
    def test_get_scaling_hint_block_new(self, caplog):
        caplog.set_level(
            idaeslog.DEBUG,
            logger="idaes",
        )

        m = ConcreteModel()
        sfx = get_scaling_hint_suffix(m)

        assert "Created new scaling hint suffix for model" in caplog.text

        assert isinstance(m.scaling_hint, Suffix)
        assert sfx is m.scaling_hint

    @pytest.mark.unit
    def test_get_scaling_hint_suffix_indexed_component_new(self):
        m = ConcreteModel()
        m.v = Var([1, 2, 3, 4])
        with pytest.raises(
            TypeError,
        ) as einfo:
            _ = get_scaling_hint_suffix(m.v[1])
        assert (
            "Component v[1] was not a BlockData, instead it was a <class 'pyomo.core.base.var.VarData'"
            in str(einfo)
        )

    @pytest.mark.unit
    def test_get_scaling_hint_suffix_indexed_block(self):
        m = ConcreteModel()
        m.b = Block([1, 2, 3, 4])

        with pytest.raises(
            TypeError,
            match=re.escape(
                "IndexedBlocks cannot have scaling hints attached to them. "
                "Please assign scaling hints to the elements of the IndexedBlock."
            ),
        ):
            get_scaling_hint_suffix(m.b)

    @pytest.mark.unit
    def test_get_scaling_hint_suffix_component_new(self):
        m = ConcreteModel()
        m.v = Var()
        with pytest.raises(TypeError) as e_obj:
            _ = get_scaling_hint_suffix(m.v)
        assert (
            "Component v was not a BlockData, instead it was a <class 'pyomo.core.base.var.ScalarVar'>"
            in str(e_obj)
        )

    @pytest.mark.unit
    def test_get_scaling_hint_suffix_block_existing(self, caplog):
        m = ConcreteModel()
        m.scaling_hint = Suffix(direction=Suffix.EXPORT)
        sfx = get_scaling_hint_suffix(m)

        assert "Created new scaling suffix for model" not in caplog.text

        assert isinstance(m.scaling_hint, Suffix)
        assert sfx is m.scaling_hint

    @pytest.mark.unit
    def test_get_scaling_hint_suffix_component_existing(self):
        m = ConcreteModel()
        m.scaling_hint = Suffix(direction=Suffix.EXPORT)
        m.v = Var()
        with pytest.raises(TypeError) as eobj:
            _ = get_scaling_hint_suffix(m.v)

        assert (
            "Component v was not a BlockData, instead it was a <class 'pyomo.core.base.var.ScalarVar'>"
            in str(eobj)
        )


def _create_model():
    m = ConcreteModel()
    m.s = Set(initialize=[1, 2, 3, 4])

    m.v = Var(m.s)

    @m.Constraint(m.s)
    def c(b, i):
        return b.v[i] == i

    @m.Expression()
    def e1(b):
        return sum(b.v[k] for k in b.s)

    m.b = Block(m.s)

    def e2_rule(b, k):
        return k * b.v2

    for bd in m.b.values():
        bd.v2 = Var()
        bd.e2 = Expression(m.s, rule=e2_rule)

    return m


class TestSuffixToFromDict:
    @pytest.fixture
    def unscaled_model(self):
        return _create_model()

    @pytest.fixture
    def scaled_model(self):
        m = _create_model()

        for bd in m.b.values():
            bd.scaling_factor = Suffix(direction=Suffix.EXPORT)
            bd.scaling_factor[bd.v2] = 10

            bd.scaling_hint = Suffix(direction=Suffix.EXPORT)
            for k in m.s:
                bd.scaling_hint[bd.e2[k]] = 1 / k

        m.scaling_factor = Suffix(direction=Suffix.EXPORT)
        for i in m.s:
            m.scaling_factor[m.v[i]] = 5 * i
            m.scaling_factor[m.c[i]] = 5**i

        m.scaling_hint = Suffix(direction=Suffix.EXPORT)
        m.scaling_hint[m.e1] = 13
        return m

    @pytest.mark.unit
    def test_suffix_to_dict(self, scaled_model):
        sdict = _suffix_to_dict(scaled_model.scaling_factor)

        assert sdict == {
            "v[1]": 5,
            "v[2]": 10,
            "v[3]": 15,
            "v[4]": 20,
            "c[1]": 5,
            "c[2]": 25,
            "c[3]": 125,
            "c[4]": 625,
        }

        hdict = _suffix_to_dict(scaled_model.scaling_hint)
        assert hdict == {"e1": 13}

    @pytest.mark.unit
    def test_suffix_from_dict(self, scaled_model):
        sdict = {
            "v[1]": 5000,
            "v[2]": 10000,
            "v[3]": 15000,
            "v[4]": 20000,
            "c[1]": 5000,
            "c[2]": 25000,
            "c[3]": 125000,
            "c[4]": 625000,
        }

        _suffix_from_dict(scaled_model.scaling_factor, sdict, overwrite=True)

        assert scaled_model.scaling_factor[scaled_model.v[1]] == 5000
        assert scaled_model.scaling_factor[scaled_model.v[2]] == 10000
        assert scaled_model.scaling_factor[scaled_model.v[3]] == 15000
        assert scaled_model.scaling_factor[scaled_model.v[4]] == 20000

        assert scaled_model.scaling_factor[scaled_model.c[1]] == 5000
        assert scaled_model.scaling_factor[scaled_model.c[2]] == 25000
        assert scaled_model.scaling_factor[scaled_model.c[3]] == 125000
        assert scaled_model.scaling_factor[scaled_model.c[4]] == 625000

        assert len(scaled_model.scaling_factor) == 8

        hdict = {"e1": 130}

        _suffix_from_dict(scaled_model.scaling_hint, hdict, overwrite=True)

        assert scaled_model.scaling_hint[scaled_model.e1] == 130
        assert len(scaled_model.scaling_hint) == 1

    @pytest.mark.unit
    def test_suffix_from_dict_invalid_component_name(self, unscaled_model):
        unscaled_model.scaling_factor = Suffix(direction=Suffix.EXPORT)
        sdict = {
            "v[1]": 5,
            "v[2]": 10,
            "v[3]": 15,
            "v[4]": 20,
            "c[1]": 5,
            "c[2]": 25,
            "c[3]": 125,
            "c[4]": 625,
            "foo": 7,
        }

        with pytest.raises(ValueError) as eobj:
            _suffix_from_dict(unscaled_model.scaling_factor, sdict, overwrite=True)
        assert "Could not find component foo on model." in str(eobj)

        # If we set verify_name=False, it should proceed
        _suffix_from_dict(
            unscaled_model.scaling_factor, sdict, overwrite=True, verify_names=False
        )

        assert unscaled_model.scaling_factor[unscaled_model.v[1]] == 5
        assert unscaled_model.scaling_factor[unscaled_model.v[2]] == 10
        assert unscaled_model.scaling_factor[unscaled_model.v[3]] == 15
        assert unscaled_model.scaling_factor[unscaled_model.v[4]] == 20

        assert unscaled_model.scaling_factor[unscaled_model.c[1]] == 5
        assert unscaled_model.scaling_factor[unscaled_model.c[2]] == 25
        assert unscaled_model.scaling_factor[unscaled_model.c[3]] == 125
        assert unscaled_model.scaling_factor[unscaled_model.c[4]] == 625

        assert len(unscaled_model.scaling_factor) == 8

    @pytest.mark.unit
    def test_collect_block_suffixes_single(self):
        m = ConcreteModel()
        m.s = Set(initialize=[1, 2, 3, 4])

        m.v = Var(m.s)

        @m.Constraint(m.s)
        def c(b, i):
            return b.v[i] == i

        @m.Expression()
        def e1(b):
            return sum(b.v[k] for k in b.s)

        m.scaling_factor = Suffix(direction=Suffix.EXPORT)
        for i in m.s:
            m.scaling_factor[m.v[i]] = 5 * i
            m.scaling_factor[m.c[i]] = 5**i

        m.scaling_hint = Suffix(direction=Suffix.EXPORT)
        m.scaling_hint[m.e1] = 13

        sdict = _collect_block_suffixes(m)

        assert sdict == {
            "scaling_factor_suffix": {
                "v[1]": 5,
                "v[2]": 10,
                "v[3]": 15,
                "v[4]": 20,
                "c[1]": 5,
                "c[2]": 25,
                "c[3]": 125,
                "c[4]": 625,
            },
            "scaling_hint_suffix": {"e1": 13},
            "subblock_suffixes": {},
        }

    @pytest.mark.unit
    def test_collect_block_suffixes_nested(self, scaled_model):
        sdict = _collect_block_suffixes(scaled_model)

        assert sdict == {
            "scaling_factor_suffix": {
                "v[1]": 5,
                "v[2]": 10,
                "v[3]": 15,
                "v[4]": 20,
                "c[1]": 5,
                "c[2]": 25,
                "c[3]": 125,
                "c[4]": 625,
            },
            "scaling_hint_suffix": {"e1": 13},
            "subblock_suffixes": {
                "b[1]": {
                    "scaling_factor_suffix": {
                        "v2": 10,
                    },
                    "scaling_hint_suffix": {
                        "e2[1]": 1,
                        "e2[2]": 1 / 2,
                        "e2[3]": 1 / 3,
                        "e2[4]": 1 / 4,
                    },
                    "subblock_suffixes": {},
                },
                "b[2]": {
                    "scaling_factor_suffix": {
                        "v2": 10,
                    },
                    "scaling_hint_suffix": {
                        "e2[1]": 1,
                        "e2[2]": 1 / 2,
                        "e2[3]": 1 / 3,
                        "e2[4]": 1 / 4,
                    },
                    "subblock_suffixes": {},
                },
                "b[3]": {
                    "scaling_factor_suffix": {
                        "v2": 10,
                    },
                    "scaling_hint_suffix": {
                        "e2[1]": 1,
                        "e2[2]": 1 / 2,
                        "e2[3]": 1 / 3,
                        "e2[4]": 1 / 4,
                    },
                    "subblock_suffixes": {},
                },
                "b[4]": {
                    "scaling_factor_suffix": {
                        "v2": 10,
                    },
                    "scaling_hint_suffix": {
                        "e2[1]": 1,
                        "e2[2]": 1 / 2,
                        "e2[3]": 1 / 3,
                        "e2[4]": 1 / 4,
                    },
                    "subblock_suffixes": {},
                },
            },
        }

    @pytest.mark.unit
    def test_collect_block_suffixes_nested_descend_false(self, scaled_model):
        sdict = _collect_block_suffixes(scaled_model, descend_into=False)

        assert sdict == {
            "scaling_factor_suffix": {
                "v[1]": 5,
                "v[2]": 10,
                "v[3]": 15,
                "v[4]": 20,
                "c[1]": 5,
                "c[2]": 25,
                "c[3]": 125,
                "c[4]": 625,
            },
            "scaling_hint_suffix": {"e1": 13},
        }

    @pytest.mark.unit
    def test_set_block_suffixes_from_dict(self, unscaled_model):
        m = unscaled_model

        # Set suffix values to retrieve
        # Only set values for some subblocks to make sure behaviour is correct
        sdict = {
            "scaling_factor_suffix": {
                "v[3]": 15,
                "v[4]": 20,
                "c[3]": 125,
                "c[4]": 625,
            },
            "scaling_hint_suffix": {"e1": 13},
            "subblock_suffixes": {
                "b[1]": {
                    "scaling_factor_suffix": {
                        "v2": 10,
                    },
                    "scaling_hint_suffix": {
                        "e2[1]": 1,
                        "e2[2]": 1 / 2,
                        "e2[3]": 1 / 3,
                        "e2[4]": 1 / 4,
                    },
                },
                "b[2]": {
                    "scaling_factor_suffix": {
                        "v2": 20,
                    },
                    "scaling_hint_suffix": {
                        "e2[1]": 1,
                        "e2[2]": 1 / 2,
                        "e2[3]": 1 / 3,
                        "e2[4]": 1 / 4,
                    },
                },
            },
        }

        _set_block_suffixes_from_dict(m, sdict)

        assert m.scaling_factor[m.v[3]] == 15
        assert m.scaling_factor[m.v[4]] == 20

        assert m.scaling_factor[m.c[3]] == 125
        assert m.scaling_factor[m.c[4]] == 625

        assert len(m.scaling_factor) == 4

        assert m.b[1].scaling_factor[m.b[1].v2] == 10
        assert len(m.b[1].scaling_factor) == 1

        assert m.b[1].scaling_hint[m.b[1].e2[1]] == 1
        assert m.b[1].scaling_hint[m.b[1].e2[2]] == 1 / 2
        assert m.b[1].scaling_hint[m.b[1].e2[3]] == 1 / 3
        assert m.b[1].scaling_hint[m.b[1].e2[4]] == 1 / 4
        assert len(m.b[1].scaling_hint) == 4

        assert m.b[2].scaling_factor[m.b[2].v2] == 20
        assert len(m.b[2].scaling_factor) == 1

        assert m.b[2].scaling_hint[m.b[2].e2[1]] == 1
        assert m.b[2].scaling_hint[m.b[2].e2[2]] == 1 / 2
        assert m.b[2].scaling_hint[m.b[2].e2[3]] == 1 / 3
        assert m.b[2].scaling_hint[m.b[2].e2[4]] == 1 / 4
        assert len(m.b[2].scaling_hint) == 4

        assert not hasattr(m.b[3], "scaling_factor")
        assert not hasattr(m.b[3], "scaling_hint")
        assert not hasattr(m.b[4], "scaling_factor")
        assert not hasattr(m.b[4], "scaling_hint")

    @pytest.mark.unit
    def test_set_block_suffixes_from_dict_overwrite_false(self, unscaled_model):
        m = unscaled_model

        # Set some existing scaling factors
        m.scaling_factor = Suffix(direction=Suffix.EXPORT)
        m.scaling_hint = Suffix(direction=Suffix.EXPORT)
        m.scaling_factor[m.v[1]] = 7
        m.scaling_factor[m.v[2]] = 17

        m.scaling_factor[m.c[1]] = 23
        m.scaling_factor[m.c[2]] = 29

        m.scaling_hint[m.e1] = 31

        for bd in m.b.values():
            bd.scaling_factor = Suffix(direction=Suffix.EXPORT)
            bd.scaling_factor[bd.v2] = 100

        # Set suffix values to retrieve
        # Only set values for some subblocks to make sure behaviour is correct
        sdict = {
            "scaling_factor_suffix": {
                "v[1]": 5,
                "v[2]": 10,
                "v[3]": 15,
                "v[4]": 20,
                "c[1]": 5,
                "c[2]": 25,
                "c[3]": 125,
                "c[4]": 625,
            },
            "scaling_hint_suffix": {"e1": 13},
            "subblock_suffixes": {
                "b[1]": {
                    "scaling_factor_suffix": {
                        "v2": 10,
                    },
                    "scaling_hint_suffix": {
                        "e2[1]": 1,
                        "e2[2]": 1 / 2,
                        "e2[3]": 1 / 3,
                        "e2[4]": 1 / 4,
                    },
                },
                "b[2]": {
                    "scaling_factor_suffix": {
                        "v2": 20,
                    },
                    "scaling_hint_suffix": {
                        "e2[1]": 1,
                        "e2[2]": 1 / 2,
                        "e2[3]": 1 / 3,
                        "e2[4]": 1 / 4,
                    },
                },
            },
        }
        sdict_copy = deepcopy(sdict)

        _set_block_suffixes_from_dict(m, sdict, overwrite=False)

        assert m.scaling_factor[m.v[1]] == 7
        assert m.scaling_factor[m.v[2]] == 17
        assert m.scaling_factor[m.v[3]] == 15
        assert m.scaling_factor[m.v[4]] == 20

        assert m.scaling_factor[m.c[1]] == 23
        assert m.scaling_factor[m.c[2]] == 29
        assert m.scaling_factor[m.c[3]] == 125
        assert m.scaling_factor[m.c[4]] == 625

        assert len(m.scaling_factor) == 8

        assert m.scaling_hint[m.e1] == 31

        for i in [1, 2]:
            assert m.b[i].scaling_factor[m.b[i].v2] == 100
            assert len(m.b[i].scaling_factor) == 1
            assert len(m.b[i].scaling_hint) == 4
            for k in m.s:
                assert m.b[i].scaling_hint[m.b[i].e2[k]] == 1 / k

        for i in [3, 4]:
            assert m.b[i].scaling_factor[m.b[i].v2] == 100
            assert len(m.b[i].scaling_factor) == 1
            assert not hasattr(m.b[i], "scaling_hint")

        # Check that we did not mutate the original dict
        assert sdict == sdict_copy

    @pytest.mark.unit
    def test_set_block_suffixes_from_dict_verify_names(self, unscaled_model):
        m = unscaled_model

        # Set suffix values to retrieve
        # Only set values for some subblocks to make sure behaviour is correct
        sdict = {
            "scaling_factor_suffix": {
                "v[1]": 5,
                "v[2]": 10,
                "v[3]": 15,
                "v[4]": 20,
                "c[1]": 5,
                "c[2]": 25,
                "c[3]": 125,
                "c[4]": 625,
            },
            "scaling_hint_suffix": {},
            "subblock_suffixes": {
                "b[1]": {
                    "scaling_factor_suffix": {
                        "v2": 10,
                    },
                    "scaling_hint_suffix": {},
                },
                "foo": {
                    "scaling_factor_suffix": {
                        "v2": 20,
                    },
                    "scaling_hint_suffix": {},
                },
            },
        }

        with pytest.raises(
            AttributeError,
            match="Model does not have a subblock named foo.",
        ):
            _set_block_suffixes_from_dict(m, sdict, verify_names=True)

    @pytest.mark.unit
    def test_scaling_factors_to_dict_suffix(self, scaled_model):
        sdict = scaling_factors_to_dict(scaled_model.scaling_factor)

        assert sdict == {
            "suffix": {
                "v[1]": 5,
                "v[2]": 10,
                "v[3]": 15,
                "v[4]": 20,
                "c[1]": 5,
                "c[2]": 25,
                "c[3]": 125,
                "c[4]": 625,
            },
            "block_name": "unknown",
        }

    @pytest.mark.unit
    def test_scaling_hints_to_dict_suffix(self, scaled_model):
        sdict = scaling_factors_to_dict(scaled_model.scaling_hint)

        assert sdict == {
            "suffix": {
                "e1": 13,
            },
            "block_name": "unknown",
        }

    @pytest.mark.unit
    def test_scaling_factors_to_dict_blockdata_descend_false(self, scaled_model):
        sdict = scaling_factors_to_dict(scaled_model, descend_into=False)

        assert sdict == {
            "scaling_factor_suffix": {
                "v[1]": 5,
                "v[2]": 10,
                "v[3]": 15,
                "v[4]": 20,
                "c[1]": 5,
                "c[2]": 25,
                "c[3]": 125,
                "c[4]": 625,
            },
            "scaling_hint_suffix": {"e1": 13},
            "block_name": "unknown",
        }

    @pytest.mark.unit
    def test_scaling_factors_to_dict_blockdata_descend_true(self, scaled_model):
        sdict = scaling_factors_to_dict(scaled_model, descend_into=True)

        assert sdict == {
            "scaling_factor_suffix": {
                "v[1]": 5,
                "v[2]": 10,
                "v[3]": 15,
                "v[4]": 20,
                "c[1]": 5,
                "c[2]": 25,
                "c[3]": 125,
                "c[4]": 625,
            },
            "scaling_hint_suffix": {"e1": 13},
            "subblock_suffixes": {
                "b[1]": {
                    "scaling_factor_suffix": {
                        "v2": 10,
                    },
                    "scaling_hint_suffix": {
                        "e2[1]": 1,
                        "e2[2]": 1 / 2,
                        "e2[3]": 1 / 3,
                        "e2[4]": 1 / 4,
                    },
                    "subblock_suffixes": {},
                },
                "b[2]": {
                    "scaling_factor_suffix": {
                        "v2": 10,
                    },
                    "scaling_hint_suffix": {
                        "e2[1]": 1,
                        "e2[2]": 1 / 2,
                        "e2[3]": 1 / 3,
                        "e2[4]": 1 / 4,
                    },
                    "subblock_suffixes": {},
                },
                "b[3]": {
                    "scaling_factor_suffix": {
                        "v2": 10,
                    },
                    "scaling_hint_suffix": {
                        "e2[1]": 1,
                        "e2[2]": 1 / 2,
                        "e2[3]": 1 / 3,
                        "e2[4]": 1 / 4,
                    },
                    "subblock_suffixes": {},
                },
                "b[4]": {
                    "scaling_factor_suffix": {
                        "v2": 10,
                    },
                    "scaling_hint_suffix": {
                        "e2[1]": 1,
                        "e2[2]": 1 / 2,
                        "e2[3]": 1 / 3,
                        "e2[4]": 1 / 4,
                    },
                    "subblock_suffixes": {},
                },
            },
            "block_name": "unknown",
        }

    @pytest.mark.unit
    def test_scaling_factors_to_dict_indexed_block(self, scaled_model):
        sdict = scaling_factors_to_dict(scaled_model.b, descend_into=True)

        assert sdict == {
            "block_data": {
                "b[1]": {
                    "scaling_factor_suffix": {
                        "v2": 10,
                    },
                    "scaling_hint_suffix": {
                        "e2[1]": 1,
                        "e2[2]": 1 / 2,
                        "e2[3]": 1 / 3,
                        "e2[4]": 1 / 4,
                    },
                    "subblock_suffixes": {},
                },
                "b[2]": {
                    "scaling_factor_suffix": {
                        "v2": 10,
                    },
                    "scaling_hint_suffix": {
                        "e2[1]": 1,
                        "e2[2]": 1 / 2,
                        "e2[3]": 1 / 3,
                        "e2[4]": 1 / 4,
                    },
                    "subblock_suffixes": {},
                },
                "b[3]": {
                    "scaling_factor_suffix": {
                        "v2": 10,
                    },
                    "scaling_hint_suffix": {
                        "e2[1]": 1,
                        "e2[2]": 1 / 2,
                        "e2[3]": 1 / 3,
                        "e2[4]": 1 / 4,
                    },
                    "subblock_suffixes": {},
                },
                "b[4]": {
                    "scaling_factor_suffix": {
                        "v2": 10,
                    },
                    "scaling_hint_suffix": {
                        "e2[1]": 1,
                        "e2[2]": 1 / 2,
                        "e2[3]": 1 / 3,
                        "e2[4]": 1 / 4,
                    },
                    "subblock_suffixes": {},
                },
            },
            "block_name": "b",
        }

    @pytest.mark.unit
    def test_scaling_factors_from_dict_suffix(self, scaled_model):
        # Partial dict of scaling factors to ensure we only touch things in the dict
        sdict = {
            "suffix": {
                "v[1]": 50,
                "v[2]": 100,
                "c[1]": 50,
                "c[2]": 250,
            },
            "block_name": "unknown",
        }
        sdict_copy = deepcopy(sdict)

        m = scaled_model

        scaling_factors_from_dict(
            m.scaling_factor, sdict, overwrite=True, verify_names=True
        )

        assert m.scaling_factor[m.v[1]] == 50
        assert m.scaling_factor[m.v[2]] == 100
        assert m.scaling_factor[m.v[3]] == 15
        assert m.scaling_factor[m.v[4]] == 20
        assert m.scaling_factor[m.c[1]] == 50
        assert m.scaling_factor[m.c[2]] == 250
        assert m.scaling_factor[m.c[3]] == 125
        assert m.scaling_factor[m.c[4]] == 625
        assert len(m.scaling_factor) == 8

        assert m.scaling_hint[m.e1] == 13
        assert len(m.scaling_hint) == 1

        for bd in m.b.values():
            assert bd.scaling_factor[bd.v2] == 10
            assert len(bd.scaling_factor) == 1
            for k in m.s:
                assert bd.scaling_hint[bd.e2[k]] == 1 / k
            assert len(bd.scaling_hint) == 4

        # Ensure we have not mutated original dict
        assert sdict == sdict_copy

    @pytest.mark.unit
    def test_scaling_factors_from_dict_suffix_overwrite_false(self, scaled_model):
        # Partial dict of scaling factors to ensure we only touch things in the dict
        sdict = {
            "suffix": {
                "v[1]": 50,
                "v[2]": 100,
                "c[1]": 50,
                "c[2]": 250,
            },
            "block_name": "unknown",
        }
        sdict_copy = deepcopy(sdict)

        m = scaled_model

        scaling_factors_from_dict(
            m.scaling_factor, sdict, overwrite=False, verify_names=True
        )

        assert m.scaling_factor[m.v[1]] == 5
        assert m.scaling_factor[m.v[2]] == 10
        assert m.scaling_factor[m.v[3]] == 15
        assert m.scaling_factor[m.v[4]] == 20
        assert m.scaling_factor[m.c[1]] == 5
        assert m.scaling_factor[m.c[2]] == 25
        assert m.scaling_factor[m.c[3]] == 125
        assert m.scaling_factor[m.c[4]] == 625
        assert len(m.scaling_factor) == 8

        assert m.scaling_hint[m.e1] == 13
        assert len(m.scaling_hint) == 1

        for bd in m.b.values():
            assert bd.scaling_factor[bd.v2] == 10
            assert len(bd.scaling_factor) == 1
            for k in m.s:
                assert bd.scaling_hint[bd.e2[k]] == 1 / k
            assert len(bd.scaling_hint) == 4

        # Ensure we have not mutated original dict
        assert sdict == sdict_copy

    @pytest.mark.unit
    def test_scaling_factors_from_dict_suffix_verify_fail(self, unscaled_model):
        unscaled_model.scaling_factor = Suffix(direction=Suffix.EXPORT)
        # Partial dict of scaling factors to ensure we only touch things in the dict
        sdict = {
            "v[1]": 50,
            "v[2]": 100,
            "c[1]": 50,
            "c[2]": 250,
            "block_name": "foo",
        }

        with pytest.raises(
            ValueError,
            match=re.escape(
                "Name of parent block (unknown) does not match that recorded in json_dict (foo)"
            ),
        ):
            scaling_factors_from_dict(
                unscaled_model.scaling_factor, sdict, overwrite=True, verify_names=True
            )

        # Ensure we have not mutated original dict
        assert sdict == {
            "v[1]": 50,
            "v[2]": 100,
            "c[1]": 50,
            "c[2]": 250,
            "block_name": "foo",
        }

    @pytest.mark.unit
    def test_missing_scaling_factor_dictionary_model(self, unscaled_model):
        # Partial dict of scaling factors to ensure we only touch things in the dict
        sdict = {
            "block_name": "unknown",
        }
        with pytest.raises(
            KeyError, match=re.escape("Missing scaling factor dictionary for model.")
        ):
            scaling_factors_from_dict(
                unscaled_model, sdict, overwrite=True, verify_names=True
            )

    @pytest.mark.unit
    def test_missing_scaling_hint_dictionary_model(self, unscaled_model):
        # Partial dict of scaling factors to ensure we only touch things in the dict
        sdict = {
            "scaling_factor_suffix": {
                "v[1]": 50,
                "v[2]": 100,
                "c[1]": 50,
                "c[2]": 250,
            },
            "block_name": "unknown",
        }
        with pytest.raises(
            KeyError, match=re.escape("Missing scaling hint dictionary for model.")
        ):
            scaling_factors_from_dict(
                unscaled_model, sdict, overwrite=True, verify_names=True
            )

    @pytest.mark.unit
    def test_scaling_factors_from_dict_block_data(self, scaled_model):
        # Partial dict of scaling factors to ensure we only touch things in the dict
        sdict = {
            "scaling_factor_suffix": {
                "v[1]": 50,
                "v[2]": 100,
                "c[1]": 50,
                "c[2]": 250,
            },
            "scaling_hint_suffix": {"e1": 13},
            "subblock_suffixes": {
                "b[1]": {
                    "scaling_factor_suffix": {
                        "v2": 20,
                    },
                    "scaling_hint_suffix": {
                        "e2[1]": 1,
                        "e2[2]": 1 / 4,
                        "e2[3]": 1 / 9,
                        "e2[4]": 1 / 16,
                    },
                },
                "b[2]": {
                    "scaling_factor_suffix": {
                        "v2": 20,
                    },
                    "scaling_hint_suffix": {
                        "e2[1]": 1,
                        "e2[2]": 1 / 4,
                        "e2[3]": 1 / 9,
                        "e2[4]": 1 / 16,
                    },
                },
            },
            "block_name": "unknown",
        }
        sdict_copy = deepcopy(sdict)

        m = scaled_model

        scaling_factors_from_dict(m, sdict, overwrite=True, verify_names=True)

        assert m.scaling_factor[m.v[1]] == 50
        assert m.scaling_factor[m.v[2]] == 100
        assert m.scaling_factor[m.v[3]] == 15
        assert m.scaling_factor[m.v[4]] == 20
        assert m.scaling_factor[m.c[1]] == 50
        assert m.scaling_factor[m.c[2]] == 250
        assert m.scaling_factor[m.c[3]] == 125
        assert m.scaling_factor[m.c[4]] == 625
        assert len(m.scaling_factor) == 8

        for i in [1, 2]:
            bd = m.b[i]
            assert bd.scaling_factor[bd.v2] == 20
            assert len(bd.scaling_factor) == 1
            assert len(bd.scaling_hint) == 4
            for k in m.s:
                assert bd.scaling_hint[bd.e2[k]] == 1 / k**2
        for i in [3, 4]:
            bd = m.b[i]
            assert bd.scaling_factor[bd.v2] == 10
            assert len(bd.scaling_factor) == 1
            assert len(bd.scaling_hint) == 4
            for k in m.s:
                assert bd.scaling_hint[bd.e2[k]] == 1 / k

        # Ensure we have not mutated original dict
        assert sdict == sdict_copy

    @pytest.mark.unit
    def test_scaling_factors_from_dict_block_data_overwrite_false(self, unscaled_model):
        m = unscaled_model

        # Set some existing scaling factors
        m.scaling_factor = Suffix(direction=Suffix.EXPORT)
        m.scaling_hint = Suffix(direction=Suffix.EXPORT)
        m.scaling_factor[m.v[1]] = 7
        m.scaling_factor[m.v[2]] = 17

        m.scaling_factor[m.c[1]] = 23
        m.scaling_factor[m.c[2]] = 29

        m.scaling_hint[m.e1] = 31

        for i in [1, 2]:
            bd = m.b[i]
            bd.scaling_factor = Suffix(direction=Suffix.EXPORT)
            bd.scaling_factor[bd.v2] = 100
            bd.scaling_hint = Suffix(direction=Suffix.EXPORT)
            for k in m.s:
                bd.scaling_hint[bd.e2[k]] = 1 / (k + 1) ** 2

        # Set suffix values to retrieve
        # Only set values for some subblocks to make sure behaviour is correct
        sdict = {
            "scaling_factor_suffix": {
                "v[1]": 5,
                "v[2]": 10,
                "v[3]": 15,
                "v[4]": 20,
                "c[1]": 5,
                "c[2]": 25,
                "c[3]": 125,
                "c[4]": 625,
            },
            "scaling_hint_suffix": {"e1": 13},
            "subblock_suffixes": {
                "b[1]": {
                    "scaling_factor_suffix": {
                        "v2": 10,
                    },
                    "scaling_hint_suffix": {
                        "e2[1]": 1,
                        "e2[2]": 1 / 2,
                        "e2[3]": 1 / 3,
                        "e2[4]": 1 / 4,
                    },
                },
                "b[2]": {
                    "scaling_factor_suffix": {
                        "v2": 20,
                    },
                    "scaling_hint_suffix": {
                        "e2[1]": 1,
                        "e2[2]": 1 / 2,
                        "e2[3]": 1 / 3,
                        "e2[4]": 1 / 4,
                    },
                },
                "b[3]": {
                    "scaling_factor_suffix": {
                        "v2": 20,
                    },
                    "scaling_hint_suffix": {
                        "e2[1]": 1,
                        "e2[2]": 1 / 2,
                        "e2[3]": 1 / 3,
                        "e2[4]": 1 / 4,
                    },
                },
                "b[4]": {
                    "scaling_factor_suffix": {
                        "v2": 20,
                    },
                    "scaling_hint_suffix": {
                        "e2[1]": 1,
                        "e2[2]": 1 / 2,
                        "e2[3]": 1 / 3,
                        "e2[4]": 1 / 4,
                    },
                },
            },
            "block_name": "unknown",
        }
        sdict_copy = deepcopy(sdict)

        scaling_factors_from_dict(m, sdict, overwrite=False, verify_names=True)

        assert m.scaling_factor[m.v[1]] == 7
        assert m.scaling_factor[m.v[2]] == 17
        assert m.scaling_factor[m.v[3]] == 15
        assert m.scaling_factor[m.v[4]] == 20

        assert m.scaling_factor[m.c[1]] == 23
        assert m.scaling_factor[m.c[2]] == 29
        assert m.scaling_factor[m.c[3]] == 125
        assert m.scaling_factor[m.c[4]] == 625

        assert len(m.scaling_factor) == 8

        assert m.scaling_hint[m.e1] == 31

        for i in [1, 2]:
            assert m.b[i].scaling_factor[m.b[i].v2] == 100
            assert len(m.b[i].scaling_factor) == 1
            assert len(m.b[i].scaling_hint) == 4
            for k in m.s:
                assert m.b[i].scaling_hint[m.b[i].e2[k]] == 1 / (k + 1) ** 2

        for i in [3, 4]:
            assert m.b[i].scaling_factor[m.b[i].v2] == 20
            assert len(m.b[i].scaling_factor) == 1
            assert len(m.b[i].scaling_hint) == 4
            for k in m.s:
                assert m.b[i].scaling_hint[m.b[i].e2[k]] == 1 / k

        # Check that we did not mutate the original dict
        assert sdict == sdict_copy

        # Ensure we have not mutated original dict
        assert sdict == sdict_copy

    @pytest.mark.unit
    def test_scaling_factors_from_dict_block_data_verify_fail(self, scaled_model):
        # Partial dict of scaling factors to ensure we only touch things in the dict
        sdict = {
            "scaling_factor_suffix": {
                "v[1]": 50,
                "v[2]": 100,
                "c[1]": 50,
                "c[2]": 250,
            },
            "scaling_hint_suffix": {},
            "block_name": "foo",
        }
        sdict_copy = deepcopy(sdict)

        with pytest.raises(
            ValueError,
            match=re.escape(
                "Block name (unknown) does not match that recorded in json_dict (foo)"
            ),
        ):
            scaling_factors_from_dict(
                scaled_model, sdict, overwrite=True, verify_names=True
            )

        # Ensure we have not mutated original dict
        assert sdict == sdict_copy

    @pytest.mark.unit
    def test_missing_scaling_factor_dict_block(self, scaled_model):
        # Partial dict of scaling factors to ensure we only touch things in the dict
        sdict = {
            "block_data": {
                "b[1]": {
                    "v2": 42,
                },
                "b[2]": {
                    "v2": 42,
                },
            },
            "block_name": "b",
        }
        with pytest.raises(
            KeyError,
            match=re.escape("Missing scaling factor dictionary for block b[1]."),
        ):
            scaling_factors_from_dict(
                scaled_model.b, sdict, overwrite=True, verify_names=True
            )

    @pytest.mark.unit
    def test_missing_scaling_hint_dict_block(self, scaled_model):
        # Partial dict of scaling factors to ensure we only touch things in the dict
        sdict = {
            "block_data": {
                "b[1]": {
                    "scaling_factor_suffix": {
                        "v2": 42,
                    }
                },
                "b[2]": {
                    "scaling_factor_suffix": {
                        "v2": 42,
                    }
                },
            },
            "block_name": "b",
        }
        with pytest.raises(
            KeyError, match=re.escape("Missing scaling hint dictionary for block b[1].")
        ):
            scaling_factors_from_dict(
                scaled_model.b, sdict, overwrite=True, verify_names=True
            )

    @pytest.mark.unit
    def test_scaling_factors_from_dict_indexed_block(self, scaled_model):
        # Partial dict of scaling factors to ensure we only touch things in the dict
        sdict = {
            "block_data": {
                "b[1]": {
                    "scaling_factor_suffix": {
                        "v2": 42,
                    },
                    "scaling_hint_suffix": {
                        "e2[1]": 1 / 4,
                        "e2[2]": 1 / 9,
                        "e2[3]": 1 / 16,
                        "e2[4]": 1 / 25,
                    },
                },
                "b[2]": {
                    "scaling_factor_suffix": {
                        "v2": 42,
                    },
                    "scaling_hint_suffix": {
                        "e2[1]": 1 / 4,
                        "e2[2]": 1 / 9,
                        "e2[3]": 1 / 16,
                        "e2[4]": 1 / 25,
                    },
                },
            },
            "block_name": "b",
        }
        sdict_copy = deepcopy(sdict)

        m = scaled_model

        scaling_factors_from_dict(m.b, sdict, overwrite=True, verify_names=True)

        assert m.scaling_factor[m.v[1]] == 5
        assert m.scaling_factor[m.v[2]] == 10
        assert m.scaling_factor[m.v[3]] == 15
        assert m.scaling_factor[m.v[4]] == 20
        assert m.scaling_factor[m.c[1]] == 5
        assert m.scaling_factor[m.c[2]] == 25
        assert m.scaling_factor[m.c[3]] == 125
        assert m.scaling_factor[m.c[4]] == 625
        assert len(m.scaling_factor) == 8

        for k in [1, 2]:
            assert m.b[k].scaling_factor[m.b[k].v2] == 42
            assert len(m.b[k].scaling_factor) == 1
            assert len(m.b[k].scaling_hint) == 4
            for i in m.s:
                assert m.b[k].scaling_hint[m.b[k].e2[i]] == 1 / (i + 1) ** 2

        for k in [3, 4]:
            assert m.b[k].scaling_factor[m.b[k].v2] == 10
            assert len(m.b[k].scaling_factor) == 1
            assert len(m.b[k].scaling_hint) == 4
            for i in m.s:
                assert m.b[k].scaling_hint[m.b[k].e2[i]] == 1 / i

        # Ensure we have not mutated original dict
        assert sdict == sdict_copy

    @pytest.mark.unit
    def test_scaling_factors_from_dict_indexed_block_overwrite_false(
        self, scaled_model
    ):
        # Partial dict of scaling factors to ensure we only touch things in the dict
        sdict = {
            "block_data": {
                "b[1]": {
                    "scaling_factor_suffix": {
                        "v2": 42,
                    },
                    "scaling_hint_suffix": {
                        "e2[1]": 1 / 4,
                        "e2[2]": 1 / 9,
                        "e2[3]": 1 / 16,
                        "e2[4]": 1 / 25,
                    },
                },
                "b[2]": {
                    "scaling_factor_suffix": {
                        "v2": 42,
                    },
                    "scaling_hint_suffix": {
                        "e2[1]": 1 / 4,
                        "e2[2]": 1 / 9,
                        "e2[3]": 1 / 16,
                        "e2[4]": 1 / 25,
                    },
                },
            },
            "block_name": "b",
        }
        sdict_copy = deepcopy(sdict)

        m = scaled_model

        scaling_factors_from_dict(m.b, sdict, overwrite=False, verify_names=True)

        assert m.scaling_factor[m.v[1]] == 5
        assert m.scaling_factor[m.v[2]] == 10
        assert m.scaling_factor[m.v[3]] == 15
        assert m.scaling_factor[m.v[4]] == 20
        assert m.scaling_factor[m.c[1]] == 5
        assert m.scaling_factor[m.c[2]] == 25
        assert m.scaling_factor[m.c[3]] == 125
        assert m.scaling_factor[m.c[4]] == 625
        assert len(m.scaling_factor) == 8

        for k in [1, 2, 3, 4]:
            assert m.b[k].scaling_factor[m.b[k].v2] == 10
            assert len(m.b[k].scaling_factor) == 1
            assert len(m.b[k].scaling_hint) == 4
            for i in m.s:
                assert m.b[k].scaling_hint[m.b[k].e2[i]] == 1 / i

        # Ensure we have not mutated original dict
        assert sdict == sdict_copy

    @pytest.mark.unit
    def test_scaling_factors_from_dict_verify_names_failure(self, scaled_model):
        # Partial dict of scaling factors to ensure we only touch things in the dict
        sdict = {
            "block_data": {
                "b[1]": {
                    "scaling_factor_suffix": {
                        "v2": 42,
                    },
                    "scaling_hint_suffix": {},
                },
                "b[2]": {
                    "scaling_factor_suffix": {
                        "v2": 42,
                    },
                    "scaling_hint_suffix": {},
                },
            },
            "block_name": "foo",
        }
        sdict_copy = deepcopy(sdict)

        with pytest.raises(
            ValueError,
            match=re.escape(
                "Block name (b) does not match that recorded in json_dict (foo)"
            ),
        ):
            scaling_factors_from_dict(
                scaled_model.b, sdict, overwrite=True, verify_names=True
            )

        # Ensure we have not mutated original dict
        assert sdict == sdict_copy

    @pytest.mark.unit
    def test_scaling_factors_from_dict_invalid_component(self, scaled_model):
        # Partial dict of scaling factors to ensure we only touch things in the dict
        sdict = {
            "block_data": {
                "b[1]": {
                    "scaling_factor_suffix": {
                        "v2": 42,
                    },
                    "scaling_hint_suffix": {},
                },
                "b[2]": {
                    "scaling_factor_suffix": {
                        "v2": 42,
                    },
                    "scaling_hint_suffix": {},
                },
            },
            "block_name": "foo",
        }

        with pytest.raises(
            TypeError, match=re.escape("v is not an instance of a Block of Suffix.")
        ):
            scaling_factors_from_dict(
                scaled_model.v, sdict, overwrite=True, verify_names=True
            )

    @pytest.mark.unit
    def test_scaling_factors_to_json_file(self, scaled_model):
        temp_context = TempfileManager.new_context()
        tmpfile = temp_context.create_tempfile(suffix=".json")

        scaling_factors_to_json_file(scaled_model, tmpfile)

        with open(tmpfile, "r") as f:
            lines = f.read()
        f.close()

        print(lines)

        expected = """{\n   "scaling_factor_suffix": {\n      "v[1]": 5,\n      "c[1]": 5,\n      "v[2]": 10,\n      "c[2]": 25,\n      "v[3]": 15,\n      "c[3]": 125,\n      "v[4]": 20,\n      "c[4]": 625\n   },\n   "scaling_hint_suffix": {\n      "e1": 13\n   },\n   "subblock_suffixes": {\n      "b[1]": {\n         "scaling_factor_suffix": {\n            "v2": 10\n         },\n         "scaling_hint_suffix": {\n            "e2[1]": 1.0,\n            "e2[2]": 0.5,\n            "e2[3]": 0.3333333333333333,\n            "e2[4]": 0.25\n         },\n         "subblock_suffixes": {}\n      },\n      "b[2]": {\n         "scaling_factor_suffix": {\n            "v2": 10\n         },\n         "scaling_hint_suffix": {\n            "e2[1]": 1.0,\n            "e2[2]": 0.5,\n            "e2[3]": 0.3333333333333333,\n            "e2[4]": 0.25\n         },\n         "subblock_suffixes": {}\n      },\n      "b[3]": {\n         "scaling_factor_suffix": {\n            "v2": 10\n         },\n         "scaling_hint_suffix": {\n            "e2[1]": 1.0,\n            "e2[2]": 0.5,\n            "e2[3]": 0.3333333333333333,\n            "e2[4]": 0.25\n         },\n         "subblock_suffixes": {}\n      },\n      "b[4]": {\n         "scaling_factor_suffix": {\n            "v2": 10\n         },\n         "scaling_hint_suffix": {\n            "e2[1]": 1.0,\n            "e2[2]": 0.5,\n            "e2[3]": 0.3333333333333333,\n            "e2[4]": 0.25\n         },\n         "subblock_suffixes": {}\n      }\n   },\n   "block_name": "unknown"\n}"""

        assert lines == expected

        # Check for clean up
        temp_context.release(remove=True)
        assert not os.path.exists(tmpfile)

    @pytest.mark.unit
    def test_scaling_factors_from_json_file(self, unscaled_model):
        fname = os.path.join(currdir, "load_scaling_factors.json")

        m = unscaled_model

        scaling_factors_from_json_file(m, fname, overwrite=True)

        assert m.scaling_factor[m.v[3]] == 15
        assert m.scaling_factor[m.v[4]] == 20

        assert m.scaling_factor[m.c[3]] == 125
        assert m.scaling_factor[m.c[4]] == 625

        assert len(m.scaling_factor) == 4

        assert m.b[1].scaling_factor[m.b[1].v2] == 10
        assert len(m.b[1].scaling_factor) == 1

        assert m.b[1].scaling_hint[m.b[1].e2[1]] == 1
        assert m.b[1].scaling_hint[m.b[1].e2[2]] == 1 / 2
        assert m.b[1].scaling_hint[m.b[1].e2[3]] == 1 / 3
        assert m.b[1].scaling_hint[m.b[1].e2[4]] == 1 / 4
        assert len(m.b[1].scaling_hint) == 4

        assert m.b[2].scaling_factor[m.b[2].v2] == 20
        assert len(m.b[2].scaling_factor) == 1

        assert m.b[2].scaling_hint[m.b[2].e2[1]] == 1
        assert m.b[2].scaling_hint[m.b[2].e2[2]] == 1 / 2
        assert m.b[2].scaling_hint[m.b[2].e2[3]] == 1 / 3
        assert m.b[2].scaling_hint[m.b[2].e2[4]] == 1 / 4
        assert len(m.b[2].scaling_hint) == 4

        assert not hasattr(m.b[3], "scaling_factor")
        assert not hasattr(m.b[3], "scaling_hint")
        assert not hasattr(m.b[4], "scaling_factor")
        assert not hasattr(m.b[4], "scaling_hint")


class TestGetScalingFactor:
    @pytest.mark.unit
    def test_get_scaling_factor_block(self):
        m = ConcreteModel()

        with pytest.raises(
            TypeError,
            match=re.escape(
                "Can get scaling factors for only VarData, ConstraintData, and (hints from) ExpressionData. Component unknown is instead <class 'pyomo.core.base.PyomoModel.ConcreteModel'>."
            ),
        ):
            get_scaling_factor(m)

    @pytest.mark.unit
    def test_get_scaling_factor(self, caplog):
        m = ConcreteModel()
        m.v = Var()

        m.scaling_factor = Suffix(direction=Suffix.EXPORT)
        m.scaling_factor[m.v] = 10

        m.scaling_hint = Suffix(direction=Suffix.EXPORT)
        m.scaling_hint[m.v] = 13
        with caplog.at_level(idaeslog.WARNING):
            sf = get_scaling_factor(m.v, default=17, warning=True)
        assert len(caplog.text) == 0
        assert sf == 10

    @pytest.mark.unit
    def test_get_scaling_factor_warning_false(self, caplog):
        m = ConcreteModel()
        m.v = Var()

        m.scaling_factor = Suffix(direction=Suffix.EXPORT)
        m.scaling_factor[m.v] = 10

        m.scaling_hint = Suffix(direction=Suffix.EXPORT)
        m.scaling_hint[m.v] = 13

        with caplog.at_level(idaeslog.WARNING):
            sf = get_scaling_factor(m.v, default=17, warning=False)
        assert len(caplog.text) == 0

        assert sf == 10

    @pytest.mark.unit
    def test_get_scaling_factor_none_warning_true(self, caplog):
        m = ConcreteModel()
        m.v = Var()

        m.scaling_hint = Suffix(direction=Suffix.EXPORT)
        m.scaling_hint[m.v] = 13

        m.scaling_factor = Suffix(direction=Suffix.EXPORT)
        with caplog.at_level(idaeslog.WARNING):
            sf = get_scaling_factor(m.v, warning=True)
        assert len(caplog.records) == 1
        assert "Missing scaling factor for v" in caplog.text
        assert sf is None

    @pytest.mark.unit
    def test_get_scaling_factor_none_warning_false(self, caplog):
        m = ConcreteModel()
        m.v = Var()

        m.scaling_hint = Suffix(direction=Suffix.EXPORT)
        m.scaling_hint[m.v] = 13

        m.scaling_factor = Suffix(direction=Suffix.EXPORT)
        with caplog.at_level(idaeslog.WARNING):
            sf = get_scaling_factor(m.v, warning=False)
        assert len(caplog.text) == 0
        assert sf is None

    @pytest.mark.unit
    def test_get_scaling_factor_none_warning_true_default(self, caplog):
        m = ConcreteModel()
        m.v = Var()

        m.scaling_hint = Suffix(direction=Suffix.EXPORT)
        m.scaling_hint[m.v] = 13

        m.scaling_factor = Suffix(direction=Suffix.EXPORT)
        with caplog.at_level(idaeslog.WARNING):
            sf = get_scaling_factor(m.v, default=7, warning=True)
        assert len(caplog.records) == 1
        assert "Missing scaling factor for v" in caplog.text
        assert sf == 7

    @pytest.mark.unit
    def test_get_scaling_factor_none_warning_false_default(self, caplog):
        m = ConcreteModel()
        m.v = Var()

        m.scaling_hint = Suffix(direction=Suffix.EXPORT)
        m.scaling_hint[m.v] = 13

        m.scaling_factor = Suffix(direction=Suffix.EXPORT)
        with caplog.at_level(idaeslog.WARNING):
            sf = get_scaling_factor(m.v, default=7, warning=False)
        assert len(caplog.text) == 0
        assert sf == 7

    @pytest.mark.unit
    def test_get_scaling_factor_no_suffix(self):
        m = ConcreteModel()
        m.v = Var()

        m.scaling_hint = Suffix(direction=Suffix.EXPORT)
        m.scaling_hint[m.v] = 13

        assert get_scaling_factor(m.v) is None

    @pytest.mark.unit
    def test_get_scaling_factor_expression(self, caplog):
        m = ConcreteModel()
        m.e = Expression(expr=4)

        # We don't want expression scaling hints to be
        # stored in the scaling factor suffix, but in
        # the event that one ends up there we want to
        # guarantee good behavior
        m.scaling_factor = Suffix(direction=Suffix.EXPORT)
        m.scaling_factor[m.e] = 13

        m.scaling_hint = Suffix(direction=Suffix.EXPORT)
        m.scaling_hint[m.e] = 10
        with caplog.at_level(idaeslog.WARNING):
            sf = get_scaling_factor(m.e, default=17)
        assert len(caplog.text) == 0
        assert sf == 10

    @pytest.mark.unit
    def test_get_scaling_factor_none(self, caplog):
        m = ConcreteModel()
        m.e = Expression(expr=4)
        m.scaling_factor = Suffix(direction=Suffix.EXPORT)
        m.scaling_factor[m.e] = 13

        m.scaling_hint = Suffix(direction=Suffix.EXPORT)

        with caplog.at_level(idaeslog.WARNING):
            sf = get_scaling_factor(m.e, warning=True)
        assert len(caplog.records) == 1
        assert "Missing scaling factor for e" in caplog.text
        assert sf is None

    @pytest.mark.unit
    def test_get_scaling_factor_expression_warning_false(self, caplog):
        m = ConcreteModel()
        m.e = Expression(expr=4)

        # We don't want expression scaling hints to be
        # stored in the scaling factor suffix, but in
        # the event that one ends up there we want to
        # guarantee good behavior
        m.scaling_factor = Suffix(direction=Suffix.EXPORT)
        m.scaling_factor[m.e] = 13

        m.scaling_hint = Suffix(direction=Suffix.EXPORT)
        m.scaling_hint[m.e] = 10
        with caplog.at_level(idaeslog.WARNING):
            sf = get_scaling_factor(m.e, warning=False)
        assert len(caplog.text) == 0
        assert sf == 10

    @pytest.mark.unit
    def test_get_scaling_factor_none_default(self, caplog):
        m = ConcreteModel()
        m.e = Expression(expr=4)
        m.scaling_factor = Suffix(direction=Suffix.EXPORT)
        m.scaling_factor[m.e] = 13

        m.scaling_hint = Suffix(direction=Suffix.EXPORT)

        with caplog.at_level(idaeslog.WARNING):
            sf = get_scaling_factor(m.e, default=17, warning=True)
        assert len(caplog.records) == 1
        assert "Missing scaling factor for e" in caplog.text
        assert sf == 17

    @pytest.mark.unit
    def test_get_scaling_factor_expression_warning_false_default(self, caplog):
        m = ConcreteModel()
        m.e = Expression(expr=4)

        # We don't want expression scaling hints to be
        # stored in the scaling factor suffix, but in
        # the event that one ends up there we want to
        # guarantee good behavior
        m.scaling_factor = Suffix(direction=Suffix.EXPORT)
        m.scaling_factor[m.e] = 13

        m.scaling_hint = Suffix(direction=Suffix.EXPORT)
        with caplog.at_level(idaeslog.WARNING):
            sf = get_scaling_factor(m.e, default=17, warning=False)
        assert len(caplog.text) == 0
        assert sf == 17

    @pytest.mark.unit
    def test_get_scaling_factor_no_suffix(self):
        m = ConcreteModel()
        m.e = Expression(expr=4)

        m.scaling_factor = Suffix(direction=Suffix.EXPORT)
        m.scaling_factor[m.e] = 13

        assert get_scaling_factor(m.e) is None

    @pytest.mark.unit
    def test_get_scaling_factor_unnamed_expression(self):
        m = ConcreteModel()
        m.v = Var()
        e = 2 * m.v
        with pytest.raises(
            TypeError,
            match=re.escape(
                "Can only get scaling hints for named expressions, but component was an unnamed expression."
            ),
        ):
            get_scaling_factor(e)

    @pytest.mark.unit
    def test_get_scaling_factor_deactivated_suffix(self):
        m = ConcreteModel()
        m.v = Var()
        m.scaling_factor = Suffix(direction=Suffix.EXPORT)
        m.scaling_factor[m.v] = 17
        m.scaling_factor.deactivate()
        assert get_scaling_factor(m.v) is None


class TestSetScalingFactor:
    @pytest.mark.unit
    def test_set_scaling_factor(self):
        m = ConcreteModel()
        m.scaling_factor = Suffix(direction=Suffix.EXPORT)

        m.v = Var()

        set_scaling_factor(m.v, 42)

        assert m.scaling_factor[m.v] == 42

        assert not hasattr(m, "scaling_hint")

    @pytest.mark.unit
    def test_set_scaling_factor_new_suffix(self, caplog):
        caplog.set_level(
            idaeslog.DEBUG,
            logger="idaes",
        )

        m = ConcreteModel()
        m.v = Var()

        set_scaling_factor(m.v, 42)

        assert m.scaling_factor[m.v] == 42.0
        assert not hasattr(m, "scaling_hint")

        assert "Created new scaling suffix for model" in caplog.text

    @pytest.mark.unit
    def test_set_scaling_factor_not_float(self):
        m = ConcreteModel()
        m.v = Var()

        with pytest.raises(
            ValueError, match="could not convert string to float: 'foo'"
        ):
            set_scaling_factor(m.v, "foo")

    @pytest.mark.unit
    def test_set_scaling_factor_negative(self):
        m = ConcreteModel()
        m.v = Var()

        with pytest.raises(
            ValueError,
            match=re.escape(
                "Scaling factor for v is negative (-42.0). "
                "Scaling factors must be strictly positive."
            ),
        ):
            set_scaling_factor(m.v, -42)

    @pytest.mark.unit
    def test_set_scaling_factor_zero(self):
        m = ConcreteModel()
        m.v = Var()

        with pytest.raises(
            ValueError,
            match=re.escape(
                "Scaling factor for v is zero. "
                "Scaling factors must be strictly positive."
            ),
        ):
            set_scaling_factor(m.v, 0)

    @pytest.mark.unit
    def test_set_scaling_factor_infinity(self):
        m = ConcreteModel()
        m.v = Var()

        with pytest.raises(
            ValueError,
            match=re.escape(
                "Scaling factor for v is infinity. " "Scaling factors must be finite."
            ),
        ):
            set_scaling_factor(m.v, float("inf"))

    @pytest.mark.unit
    def test_set_scaling_factor_NaN(self):
        m = ConcreteModel()
        m.v = Var()

        with pytest.raises(
            ValueError,
            match=re.escape("Scaling factor for v is NaN."),
        ):
            set_scaling_factor(m.v, float("NaN"))

    @pytest.mark.unit
    def test_set_scaling_factor_indexed(self, caplog):
        caplog.set_level(
            idaeslog.DEBUG,
            logger="idaes",
        )

        m = ConcreteModel()
        m.v = Var([1, 2, 3])

        with pytest.raises(
            TypeError,
            match=re.escape(
                "Component v is indexed. Set scaling factors for individual indices instead."
            ),
        ):
            set_scaling_factor(m.v, 42)

        set_scaling_factor(m.v[1], 42)

        assert m.scaling_factor[m.v[1]] == 42.0
        assert m.v[2] not in m.scaling_factor
        assert m.v[3] not in m.scaling_factor
        assert not hasattr(m, "scaling_hint")

        assert "Created new scaling suffix for model" in caplog.text

    @pytest.mark.unit
    def test_set_scaling_factor_overwrite_false(self, caplog):
        caplog.set_level(
            idaeslog.DEBUG,
            logger="idaes",
        )

        m = ConcreteModel()
        m.v = Var()

        m.scaling_factor = Suffix(direction=Suffix.EXPORT)
        m.scaling_factor[m.v] = 10

        set_scaling_factor(m.v, 42, overwrite=False)

        assert (
            "Existing scaling factor for v found and overwrite=False. "
            "Scaling factor unchanged." in caplog.text
        )
        assert m.scaling_factor[m.v] == 10

    @pytest.fixture
    def model_expr(self):
        m = ConcreteModel()
        m.e = Expression(expr=4)
        return m

    @pytest.mark.unit
    def test_set_scaling_factor_expr(self, model_expr):
        m = model_expr
        m.scaling_hint = Suffix(direction=Suffix.EXPORT)

        set_scaling_factor(m.e, 42)

        assert m.scaling_hint[m.e] == 42

        assert not hasattr(m, "scaling_factor")

    @pytest.mark.unit
    def test_set_scaling_factor_new_suffix(self, caplog, model_expr):
        caplog.set_level(
            idaeslog.DEBUG,
            logger="idaes",
        )

        m = model_expr

        set_scaling_factor(m.e, 42)

        assert m.scaling_hint[m.e] == 42.0
        assert not hasattr(m, "scaling_factor")

        assert "Created new scaling hint suffix for model" in caplog.text

    @pytest.mark.unit
    def test_set_scaling_factor_not_float(self, model_expr):
        m = model_expr

        with pytest.raises(
            ValueError, match="could not convert string to float: 'foo'"
        ):
            set_scaling_factor(m.e, "foo")

    @pytest.mark.unit
    def test_set_scaling_factor_negative(self, model_expr):
        m = model_expr

        with pytest.raises(
            ValueError,
            match=re.escape(
                "Scaling factor for e is negative (-42.0). "
                "Scaling factors must be strictly positive."
            ),
        ):
            set_scaling_factor(m.e, -42)

    @pytest.mark.unit
    def test_set_scaling_factor_zero(self, model_expr):
        m = model_expr

        with pytest.raises(
            ValueError,
            match=re.escape(
                "Scaling factor for e is zero. "
                "Scaling factors must be strictly positive."
            ),
        ):
            set_scaling_factor(m.e, 0)

    @pytest.mark.unit
    def test_set_scaling_factor_overwrite_false(self, caplog, model_expr):
        caplog.set_level(
            idaeslog.DEBUG,
            logger="idaes",
        )

        m = model_expr

        m.scaling_hint = Suffix(direction=Suffix.EXPORT)
        m.scaling_hint[m.e] = 10

        set_scaling_factor(m.e, 42, overwrite=False)

        assert (
            "Existing scaling factor for e found and overwrite=False. "
            "Scaling factor unchanged." in caplog.text
        )
        assert m.scaling_hint[m.e] == 10

    @pytest.mark.unit
    def test_set_scaling_factor_deactivated_suffix(self):
        m = ConcreteModel()
        m.scaling_factor = Suffix(direction=Suffix.EXPORT)
        m.scaling_factor.deactivate()
        m.v = Var()
        with pytest.raises(
            RuntimeError,
            match=re.escape(
                "Cannot set a scaling factor for v because the scaling_factor suffix has been deactivated."
            ),
        ):
            set_scaling_factor(m.v, 42)

        assert m.v not in m.scaling_factor


class TestDelScalingFactor:
    @pytest.mark.unit
    def test_del_scaling_factor(self):
        m = ConcreteModel()
        m.v = Var()

        m.scaling_factor = Suffix(direction=Suffix.EXPORT)
        m.scaling_hint = Suffix(direction=Suffix.EXPORT)
        m.scaling_factor[m.v] = 10

        del_scaling_factor(m.v)

        assert len(m.scaling_factor) == 0
        assert len(m.scaling_hint) == 0

    @pytest.mark.unit
    def test_del_scaling_factor_not_present(self, caplog):
        caplog.set_level(
            idaeslog.DEBUG,
            logger="idaes",
        )

        m = ConcreteModel()
        m.v = Var()

        m.scaling_factor = Suffix(direction=Suffix.EXPORT)
        m.scaling_hint = Suffix(direction=Suffix.EXPORT)

        del_scaling_factor(m.v)

        assert len(m.scaling_factor) == 0
        assert len(m.scaling_hint) == 0

    @pytest.mark.unit
    def test_del_scaling_factor_delete_empty(self, caplog):
        caplog.set_level(
            idaeslog.DEBUG,
            logger="idaes",
        )

        m = ConcreteModel()
        m.v = Var()

        m.scaling_factor = Suffix(direction=Suffix.EXPORT)
        m.scaling_hint = Suffix(direction=Suffix.EXPORT)
        m.scaling_factor[m.v] = 10

        del_scaling_factor(m.v, delete_empty_suffix=True)

        assert not hasattr(m, "scaling_factor")

        assert "Deleting empty scaling suffix from unknown" in caplog.text

        assert len(m.scaling_hint) == 0

    @pytest.mark.unit
    def test_del_scaling_factor(self):
        m = ConcreteModel()
        m.e = Expression(expr=4)

        m.scaling_factor = Suffix(direction=Suffix.EXPORT)
        m.scaling_hint = Suffix(direction=Suffix.EXPORT)
        m.scaling_hint[m.e] = 10

        del_scaling_factor(m.e)

        assert len(m.scaling_factor) == 0
        assert len(m.scaling_hint) == 0

    @pytest.mark.unit
    def test_del_scaling_factor_not_present(self, caplog):
        caplog.set_level(
            idaeslog.DEBUG,
            logger="idaes",
        )

        m = ConcreteModel()
        m.e = Expression(expr=4)

        m.scaling_factor = Suffix(direction=Suffix.EXPORT)
        m.scaling_hint = Suffix(direction=Suffix.EXPORT)

        del_scaling_factor(m.e)

        assert len(m.scaling_factor) == 0
        assert len(m.scaling_hint) == 0

    @pytest.mark.unit
    def test_del_scaling_factor_delete_empty(self, caplog):
        caplog.set_level(
            idaeslog.DEBUG,
            logger="idaes",
        )

        m = ConcreteModel()
        m.e = Expression(expr=4)

        m.scaling_factor = Suffix(direction=Suffix.EXPORT)
        m.scaling_hint = Suffix(direction=Suffix.EXPORT)
        m.scaling_hint[m.e] = 10

        del_scaling_factor(m.e, delete_empty_suffix=True)

        assert not hasattr(m, "scaling_hint")

        assert "Deleting empty scaling hint suffix from unknown" in caplog.text

        assert len(m.scaling_factor) == 0


class TestReportScalingFactors:
    @pytest.fixture
    def model(self):
        m = _create_model()

        # Need to check all possible behaviours
        # Set values for half the variables (indexes 1 and 3)
        for i in [1, 3]:
            m.v[i].set_value(42)
            m.b[i].v2.set_value(42)
            m.b[i].c2 = Constraint(expr=m.b[i].v2 == m.b[i].e2[i])

        # Set scaling factors for half the components (indexed 1 and 2)
        m.scaling_factor = Suffix(direction=Suffix.EXPORT)
        for i in [1, 2]:
            m.scaling_factor[m.v[i]] = 5 * i
            m.scaling_factor[m.c[i]] = 5**i
            m.b[i].scaling_factor = Suffix(direction=Suffix.EXPORT)
            m.b[i].scaling_factor[m.b[i].v2] = 10

        set_scaling_factor(m.b[1].c2, 89)

        for i in [1, 3]:
            m.b[i].scaling_hint = Suffix(direction=Suffix.EXPORT)
            for k in [2, 4]:
                m.b[i].scaling_hint[m.b[i].e2[k]] = 1 / k

        return m

    @pytest.mark.unit
    def test_report_scaling_factors_all(self, model):
        stream = StringIO()

        report_scaling_factors(model, descend_into=True, stream=stream)

        expected = """Scaling Factors for model

Variable    Scaling Factor    Value        Scaled Value
v[1]        5.000E+00         4.200E+01    2.100E+02
v[2]        1.000E+01         None         None
v[3]        None              4.200E+01    4.200E+01
v[4]        None              None         None
b[1].v2     1.000E+01         4.200E+01    4.200E+02
b[2].v2     1.000E+01         None         None
b[3].v2     None              4.200E+01    4.200E+01
b[4].v2     None              None         None

Constraint    Scaling Factor
c[1]          5.000E+00
c[2]          2.500E+01
c[3]          None
c[4]          None
b[1].c2       8.900E+01
b[3].c2       None

Expression    Scaling Hint
e1            None
b[1].e2[1]    None
b[1].e2[2]    5.000E-01
b[1].e2[3]    None
b[1].e2[4]    2.500E-01
b[2].e2[1]    None
b[2].e2[2]    None
b[2].e2[3]    None
b[2].e2[4]    None
b[3].e2[1]    None
b[3].e2[2]    5.000E-01
b[3].e2[3]    None
b[3].e2[4]    2.500E-01
b[4].e2[1]    None
b[4].e2[2]    None
b[4].e2[3]    None
b[4].e2[4]    None
"""

        print(stream.getvalue())
        assert stream.getvalue() == expected

    @pytest.mark.unit
    def test_report_scaling_factors_descend_false(self, model):
        stream = StringIO()

        report_scaling_factors(model, descend_into=False, stream=stream)

        expected = """Scaling Factors for model

Variable    Scaling Factor    Value        Scaled Value
v[1]        5.000E+00         4.200E+01    2.100E+02
v[2]        1.000E+01         None         None
v[3]        None              4.200E+01    4.200E+01
v[4]        None              None         None

Constraint    Scaling Factor
c[1]          5.000E+00
c[2]          2.500E+01
c[3]          None
c[4]          None

Expression    Scaling Hint
e1            None
"""

        assert stream.getvalue() == expected

    @pytest.mark.unit
    def test_report_scaling_factors_vars_only(self, model):
        stream = StringIO()

        report_scaling_factors(model, descend_into=True, ctype=Var, stream=stream)

        expected = """Scaling Factors for model

Variable    Scaling Factor    Value        Scaled Value
v[1]        5.000E+00         4.200E+01    2.100E+02
v[2]        1.000E+01         None         None
v[3]        None              4.200E+01    4.200E+01
v[4]        None              None         None
b[1].v2     1.000E+01         4.200E+01    4.200E+02
b[2].v2     1.000E+01         None         None
b[3].v2     None              4.200E+01    4.200E+01
b[4].v2     None              None         None
"""

        assert stream.getvalue() == expected

    @pytest.mark.unit
    def test_report_scaling_factors_constraints_only(self, model):
        stream = StringIO()

        report_scaling_factors(
            model, descend_into=True, stream=stream, ctype=Constraint
        )

        expected = """Scaling Factors for model

Constraint    Scaling Factor
c[1]          5.000E+00
c[2]          2.500E+01
c[3]          None
c[4]          None
b[1].c2       8.900E+01
b[3].c2       None
"""

        assert stream.getvalue() == expected

    @pytest.mark.unit
    def test_report_scaling_factors_indexed_block(self, model):
        stream = StringIO()

        report_scaling_factors(model.b, descend_into=True, stream=stream)

        expected = """Scaling Factors for block b

Variable    Scaling Factor    Value        Scaled Value
b[1].v2     1.000E+01         4.200E+01    4.200E+02
b[2].v2     1.000E+01         None         None
b[3].v2     None              4.200E+01    4.200E+01
b[4].v2     None              None         None

Constraint    Scaling Factor
b[1].c2       8.900E+01
b[3].c2       None

Expression    Scaling Hint
b[1].e2[1]    None
b[1].e2[2]    5.000E-01
b[1].e2[3]    None
b[1].e2[4]    2.500E-01
b[2].e2[1]    None
b[2].e2[2]    None
b[2].e2[3]    None
b[2].e2[4]    None
b[3].e2[1]    None
b[3].e2[2]    5.000E-01
b[3].e2[3]    None
b[3].e2[4]    2.500E-01
b[4].e2[1]    None
b[4].e2[2]    None
b[4].e2[3]    None
b[4].e2[4]    None
"""

        assert stream.getvalue() == expected

    @pytest.mark.unit
    def test_report_scaling_factors_not_block(self, model):
        stream = StringIO()

        with pytest.raises(
            TypeError,
            match="report_scaling_factors: blk must be an instance of a Pyomo Block.",
        ):
            report_scaling_factors(
                model.v, descend_into=True, stream=stream, ctype=Constraint
            )

    @pytest.mark.unit
    def test_report_scaling_factors_invalid_ctype(self, model):
        stream = StringIO()

        with pytest.raises(
            ValueError,
            match="report_scaling_factors only supports None, Var, Constraint, or Expression for argument ctype: "
            "received foo.",
        ):
            report_scaling_factors(model, descend_into=True, stream=stream, ctype="foo")


# Adopted from old scaling tools
# originally by John Eslick
@pytest.mark.unit
def test_find_unscaled_vars_and_constraints():
    m = ConcreteModel()
    m.b = Block()
    m.x = Var(initialize=1e6)
    m.y = Var(initialize=1e-8)
    m.z = Var(initialize=1e-20)
    m.c1 = Constraint(expr=m.x == 0)
    m.c2 = Constraint(expr=m.y == 0)
    m.b.w = Var([1, 2, 3], initialize=1e10)
    m.b.c1 = Constraint(expr=m.b.w[1] == 0)
    m.b.c2 = Constraint(expr=m.b.w[2] == 0)
    m.c3 = Constraint(expr=m.z == 0)

    set_scaling_factor(m.x, 1)
    set_scaling_factor(m.b.w[1], 2)
    set_scaling_factor(m.c1, 1)
    set_scaling_factor(m.b.c1, 1)
    set_scaling_factor(m.c3, 1)

    a = [id(v) for v in unscaled_variables_generator(m)]
    # Make sure we pick up the right variables
    assert id(m.x) not in a
    assert id(m.y) in a
    assert id(m.z) in a
    assert id(m.b.w[1]) not in a
    assert id(m.b.w[2]) in a
    assert id(m.b.w[3]) in a
    assert len(a) == 4  # make sure we didn't pick up any other random stuff

    b = [id(v) for v in list_unscaled_variables(m)]
    for foo, bar in zip(a, b):
        assert foo == bar

    a = [id(v) for v in unscaled_constraints_generator(m)]
    assert id(m.c1) not in a
    assert id(m.b.c1) not in a
    assert id(m.c2) in a
    assert id(m.b.c2) in a
    assert id(m.c3) not in a
    assert len(a) == 2  # make sure we didn't pick up any other random stuff

    b = [id(v) for v in list_unscaled_constraints(m)]
    for foo, bar in zip(a, b):
        assert foo == bar


# Tests based on John Eslick's originals for the old scaling tools
@pytest.mark.skipif(
    not AmplInterface.available(), reason="pynumero_ASL is not available"
)
class TestJacobianMethods:
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
        assert number_activated_objectives(m) == 0
        jac, nlp = get_jacobian(m)
        assert not hasattr(m, "scaling_factor")
        assert not hasattr(m, "scaling_hint")
        assert number_activated_objectives(m) == 0

        c1_row = nlp._condata_to_idx[m.c1]
        c2_row = nlp._condata_to_idx[m.c2]
        c3_row = nlp._condata_to_idx[m.c3]
        x_col = nlp._vardata_to_idx[m.x]
        y_col = nlp._vardata_to_idx[m.y]
        z_col = nlp._vardata_to_idx[m.z]

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
        jac, _ = get_jacobian(m, include_scaling_factors=False)
        assert len(m.scaling_factor) == 4
        assert not hasattr(m, "scaling_hint")
        assert jac[c1_row, x_col] == pytest.approx(-1e6)

        # Check the scaled jacobian calculation
        jac_scaled, _ = get_jacobian(m)
        assert len(m.scaling_factor) == 4
        assert not hasattr(m, "scaling_hint")
        assert jac_scaled[c1_row, x_col] == pytest.approx(-1000)
        assert jac_scaled[c1_row, y_col] == pytest.approx(-1000)
        assert jac_scaled[c1_row, z_col] == pytest.approx(0.01)

    @pytest.mark.unit
    def test_scale_no_var_scale(self, model):
        m = model
        jac_scaled, nlp = get_jacobian(
            m,
            include_scaling_factors=False,
            include_ipopt_autoscaling=True,
            min_scale=1e-6,
        )
        # get_scaling_factor isn't called here so the suffix shouldn't exist
        assert not hasattr(m, "scaling_factor")
        assert not hasattr(m, "scaling_hint")

        c1_row = nlp._condata_to_idx[m.c1]
        c2_row = nlp._condata_to_idx[m.c2]
        c3_row = nlp._condata_to_idx[m.c3]
        x_col = nlp._vardata_to_idx[m.x]
        y_col = nlp._vardata_to_idx[m.y]
        z_col = nlp._vardata_to_idx[m.z]

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
        jac_scaled, nlp = get_jacobian(
            m,
            include_scaling_factors=True,
            include_ipopt_autoscaling=True,
            min_scale=1e-6,
        )
        assert len(m.scaling_factor) == 3
        assert not hasattr(m, "scaling_hint")

        c1_row = nlp._condata_to_idx[m.c1]
        c2_row = nlp._condata_to_idx[m.c2]
        c3_row = nlp._condata_to_idx[m.c3]
        x_col = nlp._vardata_to_idx[m.x]
        y_col = nlp._vardata_to_idx[m.y]
        z_col = nlp._vardata_to_idx[m.z]

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

        jac_scaled, nlp = get_jacobian(
            m,
            include_scaling_factors=False,
            include_ipopt_autoscaling=True,
            min_scale=1e-6,
        )
        assert len(m.scaling_factor) == 3
        assert not hasattr(m, "scaling_hint")

        c1_row = nlp._condata_to_idx[m.c1]
        c2_row = nlp._condata_to_idx[m.c2]
        c3_row = nlp._condata_to_idx[m.c3]
        x_col = nlp._vardata_to_idx[m.x]
        y_col = nlp._vardata_to_idx[m.y]
        z_col = nlp._vardata_to_idx[m.z]

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

        n = jacobian_cond(m, scaled=True)
        assert n == pytest.approx(687.47, rel=1e-3)
        n = jacobian_cond(m, scaled=False)
        assert n == pytest.approx(7.50567e7, rel=1e-3)

        # Nonsquare condition number
        m.c3.deactivate()

        # Scaled
        with caplog.at_level(idaeslog.INFO):
            n = jacobian_cond(m, scaled=True)
        assert (
            "Nonsquare Jacobian. Using pseudoinverse to calculate Frobenius norm."
        ) in caplog.text
        assert n == pytest.approx(500.367, rel=1e-3)

        # Unscaled
        with caplog.at_level(idaeslog.INFO):
            n = jacobian_cond(m, scaled=False)
        assert (
            "Nonsquare Jacobian. Using pseudoinverse to calculate Frobenius norm."
        ) in caplog.text
        assert n == pytest.approx(2.23741e5, rel=1e-3)

    @pytest.mark.unit
    def test_condition_number_none(self, model):
        with pytest.raises(
            RuntimeError,
            match=re.escape(
                "User must provide either a Pyomo model or a Jacobian "
                "to calculate the condition number."
            ),
        ):
            _ = jacobian_cond()
