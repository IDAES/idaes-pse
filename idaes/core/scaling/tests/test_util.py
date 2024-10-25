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

Author: Andrew Lee
"""
from io import StringIO
import os
import pytest
import re

from pyomo.environ import Block, Constraint, ConcreteModel, Set, Suffix, Var
from pyomo.common.fileutils import this_file_dir
from pyomo.common.tempfiles import TempfileManager

from idaes.core.scaling.util import (
    get_scaling_suffix,
    get_scaling_factor,
    set_scaling_factor,
    del_scaling_factor,
    _suffix_to_dict,
    _suffix_from_dict,
    _collect_block_suffixes,
    _set_block_suffixes_from_dict,
    scaling_factors_to_dict,
    scaling_factors_from_dict,
    scaling_factors_to_json_file,
    scaling_factors_from_json_file,
    report_scaling_factors,
)
import idaes.logger as idaeslog

currdir = this_file_dir()


class TestGetScalingSuffix:
    @pytest.mark.unit
    def test_get_scaling_suffix_block_new(self, caplog):
        caplog.set_level(
            idaeslog.DEBUG,
            logger="idaes",
        )

        m = ConcreteModel()
        sfx = get_scaling_suffix(m)

        assert "Created new scaling suffix for unknown" in caplog.text

        assert isinstance(m.scaling_factor, Suffix)
        assert sfx is m.scaling_factor

    @pytest.mark.unit
    def test_get_scaling_suffix_indexed_component_new(self, caplog):
        caplog.set_level(
            idaeslog.DEBUG,
            logger="idaes",
        )

        m = ConcreteModel()
        m.v = Var([1, 2, 3, 4])
        sfx = get_scaling_suffix(m.v[1])

        assert "Created new scaling suffix for unknown" in caplog.text

        assert isinstance(m.scaling_factor, Suffix)
        assert sfx is m.scaling_factor

    @pytest.mark.unit
    def test_get_scaling_suffix_indexed_block(self):
        m = ConcreteModel()
        m.b = Block([1, 2, 3, 4])

        with pytest.raises(
            TypeError,
            match="IndexedBlocks cannot have scaling factors attached to them. "
            "Please assign scaling factors to the elements of the IndexedBlock.",
        ):
            get_scaling_suffix(m.b)

    @pytest.mark.unit
    def test_get_scaling_suffix_component_new(self, caplog):
        caplog.set_level(
            idaeslog.DEBUG,
            logger="idaes",
        )

        m = ConcreteModel()
        m.v = Var()
        sfx = get_scaling_suffix(m.v)

        assert "Created new scaling suffix for unknown" in caplog.text

        assert isinstance(m.scaling_factor, Suffix)
        assert sfx is m.scaling_factor

    @pytest.mark.unit
    def test_get_scaling_suffix_block_existing(self, caplog):
        m = ConcreteModel()
        m.scaling_factor = Suffix(direction=Suffix.EXPORT)
        sfx = get_scaling_suffix(m)

        assert "Created new scaling suffix for unknown" not in caplog.text

        assert isinstance(m.scaling_factor, Suffix)
        assert sfx is m.scaling_factor

    @pytest.mark.unit
    def test_get_scaling_suffix_component_existing(self, caplog):
        caplog.set_level(
            idaeslog.DEBUG,
            logger="idaes",
        )

        m = ConcreteModel()
        m.scaling_factor = Suffix(direction=Suffix.EXPORT)
        m.v = Var()
        sfx = get_scaling_suffix(m.v)

        assert "Created new scaling suffix for unknown" not in caplog.text

        assert isinstance(m.scaling_factor, Suffix)
        assert sfx is m.scaling_factor


class TestSuffixToFromDict:
    @pytest.fixture
    def model(self):
        m = ConcreteModel()
        m.s = Set(initialize=[1, 2, 3, 4])

        m.v = Var(m.s)

        @m.Constraint(m.s)
        def c(b, i):
            return b.v[i] == i

        m.b = Block(m.s)

        for bd in m.b.values():
            bd.v2 = Var()

            bd.scaling_factor = Suffix(direction=Suffix.EXPORT)
            bd.scaling_factor[bd.v2] = 10

        m.scaling_factor = Suffix(direction=Suffix.EXPORT)
        for i in m.s:
            m.scaling_factor[m.v[i]] = 5 * i
            m.scaling_factor[m.c[i]] = 5**i

        return m

    @pytest.mark.unit
    def test_suffix_to_dict(self, model):
        sdict = _suffix_to_dict(model.scaling_factor)

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

    @pytest.mark.unit
    def test_suffix_from_dict(self, model):
        sdict = {
            "v[1]": 5,
            "v[2]": 10,
            "v[3]": 15,
            "v[4]": 20,
            "c[1]": 5,
            "c[2]": 25,
            "c[3]": 125,
            "c[4]": 625,
        }

        _suffix_from_dict(model.scaling_factor, sdict, overwrite=True)

        assert model.scaling_factor[model.v[1]] == 5
        assert model.scaling_factor[model.v[2]] == 10
        assert model.scaling_factor[model.v[3]] == 15
        assert model.scaling_factor[model.v[4]] == 20

        assert model.scaling_factor[model.c[1]] == 5
        assert model.scaling_factor[model.c[2]] == 25
        assert model.scaling_factor[model.c[3]] == 125
        assert model.scaling_factor[model.c[4]] == 625

        assert len(model.scaling_factor) == 8

    @pytest.mark.unit
    def test_suffix_from_dict_invalid_component_name(self, model):
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

        with pytest.raises(
            ValueError,
            match=re.escape("Could not find component foo on block unknown."),
        ):
            _suffix_from_dict(model.scaling_factor, sdict, overwrite=True)

        # If we set verify_name=False, it should proceed
        _suffix_from_dict(
            model.scaling_factor, sdict, overwrite=True, verify_names=False
        )

        assert model.scaling_factor[model.v[1]] == 5
        assert model.scaling_factor[model.v[2]] == 10
        assert model.scaling_factor[model.v[3]] == 15
        assert model.scaling_factor[model.v[4]] == 20

        assert model.scaling_factor[model.c[1]] == 5
        assert model.scaling_factor[model.c[2]] == 25
        assert model.scaling_factor[model.c[3]] == 125
        assert model.scaling_factor[model.c[4]] == 625

        assert len(model.scaling_factor) == 8

    @pytest.mark.unit
    def test_collect_block_suffixes_single(self):
        m = ConcreteModel()
        m.s = Set(initialize=[1, 2, 3, 4])

        m.v = Var(m.s)

        @m.Constraint(m.s)
        def c(b, i):
            return b.v[i] == i

        m.scaling_factor = Suffix(direction=Suffix.EXPORT)
        for i in m.s:
            m.scaling_factor[m.v[i]] = 5 * i
            m.scaling_factor[m.c[i]] = 5**i

        sdict = _collect_block_suffixes(m)

        assert sdict == {
            "v[1]": 5,
            "v[2]": 10,
            "v[3]": 15,
            "v[4]": 20,
            "c[1]": 5,
            "c[2]": 25,
            "c[3]": 125,
            "c[4]": 625,
            "subblock_suffixes": {},
        }

    @pytest.mark.unit
    def test_collect_block_suffixes_nested(self, model):
        sdict = _collect_block_suffixes(model)

        assert sdict == {
            "v[1]": 5,
            "v[2]": 10,
            "v[3]": 15,
            "v[4]": 20,
            "c[1]": 5,
            "c[2]": 25,
            "c[3]": 125,
            "c[4]": 625,
            "subblock_suffixes": {
                "b[1]": {
                    "v2": 10,
                    "subblock_suffixes": {},
                },
                "b[2]": {
                    "v2": 10,
                    "subblock_suffixes": {},
                },
                "b[3]": {
                    "v2": 10,
                    "subblock_suffixes": {},
                },
                "b[4]": {
                    "v2": 10,
                    "subblock_suffixes": {},
                },
            },
        }

    @pytest.mark.unit
    def test_collect_block_suffixes_nested_descend_false(self, model):
        sdict = _collect_block_suffixes(model, descend_into=False)

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

    @pytest.mark.unit
    def test_set_block_suffixes_from_dict(self):
        m = ConcreteModel()
        m.s = Set(initialize=[1, 2, 3, 4])

        m.v = Var(m.s)

        @m.Constraint(m.s)
        def c(b, i):
            return b.v[i] == i

        m.b = Block(m.s)

        for bd in m.b.values():
            bd.v2 = Var()

        # Set suffix values to retrieve
        # Only set values for some subblocks to make sure behaviour is correct
        sdict = {
            "v[1]": 5,
            "v[2]": 10,
            "v[3]": 15,
            "v[4]": 20,
            "c[1]": 5,
            "c[2]": 25,
            "c[3]": 125,
            "c[4]": 625,
            "subblock_suffixes": {
                "b[1]": {
                    "v2": 10,
                },
                "b[2]": {
                    "v2": 20,
                },
            },
        }

        _set_block_suffixes_from_dict(m, sdict)

        assert m.scaling_factor[m.v[1]] == 5
        assert m.scaling_factor[m.v[2]] == 10
        assert m.scaling_factor[m.v[3]] == 15
        assert m.scaling_factor[m.v[4]] == 20

        assert m.scaling_factor[m.c[1]] == 5
        assert m.scaling_factor[m.c[2]] == 25
        assert m.scaling_factor[m.c[3]] == 125
        assert m.scaling_factor[m.c[4]] == 625

        assert len(m.scaling_factor) == 8

        assert m.b[1].scaling_factor[m.b[1].v2] == 10
        assert len(m.b[1].scaling_factor) == 1

        assert m.b[2].scaling_factor[m.b[2].v2] == 20
        assert len(m.b[2].scaling_factor) == 1

        assert not hasattr(m.b[3], "scaling_factor")
        assert not hasattr(m.b[4], "scaling_factor")

    @pytest.mark.unit
    def test_set_block_suffixes_from_dict_overwrite_false(self):
        m = ConcreteModel()
        m.s = Set(initialize=[1, 2, 3, 4])

        m.v = Var(m.s)

        @m.Constraint(m.s)
        def c(b, i):
            return b.v[i] == i

        m.b = Block(m.s)

        for bd in m.b.values():
            bd.v2 = Var()

            # Set some existing scaling factors
            bd.scaling_factor = Suffix(direction=Suffix.EXPORT)
            bd.scaling_factor[bd.v2] = 100

        # Set suffix values to retrieve
        # Only set values for some subblocks to make sure behaviour is correct
        sdict = {
            "v[1]": 5,
            "v[2]": 10,
            "v[3]": 15,
            "v[4]": 20,
            "c[1]": 5,
            "c[2]": 25,
            "c[3]": 125,
            "c[4]": 625,
            "subblock_suffixes": {
                "b[1]": {
                    "v2": 10,
                },
                "b[2]": {
                    "v2": 20,
                },
            },
        }

        _set_block_suffixes_from_dict(m, sdict, overwrite=False)

        assert m.scaling_factor[m.v[1]] == 5
        assert m.scaling_factor[m.v[2]] == 10
        assert m.scaling_factor[m.v[3]] == 15
        assert m.scaling_factor[m.v[4]] == 20

        assert m.scaling_factor[m.c[1]] == 5
        assert m.scaling_factor[m.c[2]] == 25
        assert m.scaling_factor[m.c[3]] == 125
        assert m.scaling_factor[m.c[4]] == 625

        assert len(m.scaling_factor) == 8

        for i in [1, 2, 3, 4]:
            assert m.b[i].scaling_factor[m.b[i].v2] == 100
            assert len(m.b[i].scaling_factor) == 1

        # Check that we did not mutate the original dict
        assert sdict == {
            "v[1]": 5,
            "v[2]": 10,
            "v[3]": 15,
            "v[4]": 20,
            "c[1]": 5,
            "c[2]": 25,
            "c[3]": 125,
            "c[4]": 625,
            "subblock_suffixes": {
                "b[1]": {
                    "v2": 10,
                },
                "b[2]": {
                    "v2": 20,
                },
            },
        }

    @pytest.mark.unit
    def test_set_block_suffixes_from_dict_verify_names(self):
        m = ConcreteModel()
        m.s = Set(initialize=[1, 2, 3, 4])

        m.v = Var(m.s)

        @m.Constraint(m.s)
        def c(b, i):
            return b.v[i] == i

        m.b = Block(m.s)

        for bd in m.b.values():
            bd.v2 = Var()

        # Set suffix values to retrieve
        # Only set values for some subblocks to make sure behaviour is correct
        sdict = {
            "v[1]": 5,
            "v[2]": 10,
            "v[3]": 15,
            "v[4]": 20,
            "c[1]": 5,
            "c[2]": 25,
            "c[3]": 125,
            "c[4]": 625,
            "subblock_suffixes": {
                "b[1]": {
                    "v2": 10,
                },
                "foo": {
                    "v2": 20,
                },
            },
        }

        with pytest.raises(
            AttributeError,
            match="Block unknown does not have a subblock named foo.",
        ):
            _set_block_suffixes_from_dict(m, sdict, verify_names=True)

    @pytest.mark.unit
    def test_scaling_factors_to_dict_suffix(self, model):
        sdict = scaling_factors_to_dict(model.scaling_factor)

        assert sdict == {
            "v[1]": 5,
            "v[2]": 10,
            "v[3]": 15,
            "v[4]": 20,
            "c[1]": 5,
            "c[2]": 25,
            "c[3]": 125,
            "c[4]": 625,
            "block_name": "unknown",
        }

    @pytest.mark.unit
    def test_scaling_factors_to_dict_blockdata_descend_false(self, model):
        sdict = scaling_factors_to_dict(model, descend_into=False)

        assert sdict == {
            "v[1]": 5,
            "v[2]": 10,
            "v[3]": 15,
            "v[4]": 20,
            "c[1]": 5,
            "c[2]": 25,
            "c[3]": 125,
            "c[4]": 625,
            "block_name": "unknown",
        }

    @pytest.mark.unit
    def test_scaling_factors_to_dict_blockdata_descend_true(self, model):
        sdict = scaling_factors_to_dict(model, descend_into=True)

        assert sdict == {
            "v[1]": 5,
            "v[2]": 10,
            "v[3]": 15,
            "v[4]": 20,
            "c[1]": 5,
            "c[2]": 25,
            "c[3]": 125,
            "c[4]": 625,
            "subblock_suffixes": {
                "b[1]": {
                    "v2": 10,
                    "subblock_suffixes": {},
                },
                "b[2]": {
                    "v2": 10,
                    "subblock_suffixes": {},
                },
                "b[3]": {
                    "v2": 10,
                    "subblock_suffixes": {},
                },
                "b[4]": {
                    "v2": 10,
                    "subblock_suffixes": {},
                },
            },
            "block_name": "unknown",
        }

    @pytest.mark.unit
    def test_scaling_factors_to_dict_indexed_block(self, model):
        sdict = scaling_factors_to_dict(model.b, descend_into=True)

        assert sdict == {
            "block_datas": {
                "b[1]": {
                    "v2": 10,
                    "subblock_suffixes": {},
                },
                "b[2]": {
                    "v2": 10,
                    "subblock_suffixes": {},
                },
                "b[3]": {
                    "v2": 10,
                    "subblock_suffixes": {},
                },
                "b[4]": {
                    "v2": 10,
                    "subblock_suffixes": {},
                },
            },
            "block_name": "b",
        }

    @pytest.mark.unit
    def test_scaling_factors_from_dict_suffix(self, model):
        # Partial dict of scaling factors to ensure we only touch things in the dict
        sdict = {
            "v[1]": 50,
            "v[2]": 100,
            "c[1]": 50,
            "c[2]": 250,
            "block_name": "unknown",
        }

        scaling_factors_from_dict(
            model.scaling_factor, sdict, overwrite=True, verify_names=True
        )

        assert model.scaling_factor[model.v[1]] == 50
        assert model.scaling_factor[model.v[2]] == 100
        assert model.scaling_factor[model.v[3]] == 15
        assert model.scaling_factor[model.v[4]] == 20
        assert model.scaling_factor[model.c[1]] == 50
        assert model.scaling_factor[model.c[2]] == 250
        assert model.scaling_factor[model.c[3]] == 125
        assert model.scaling_factor[model.c[4]] == 625
        assert len(model.scaling_factor) == 8

        for bd in model.b.values():
            assert bd.scaling_factor[bd.v2] == 10
            assert len(bd.scaling_factor) == 1

        # Ensure we have not mutated original dict
        assert sdict == {
            "v[1]": 50,
            "v[2]": 100,
            "c[1]": 50,
            "c[2]": 250,
            "block_name": "unknown",
        }

    @pytest.mark.unit
    def test_scaling_factors_from_dict_suffix_overwrite_false(self, model):
        # Partial dict of scaling factors to ensure we only touch things in the dict
        sdict = {
            "v[1]": 50,
            "v[2]": 100,
            "c[1]": 50,
            "c[2]": 250,
            "block_name": "unknown",
        }

        scaling_factors_from_dict(
            model.scaling_factor, sdict, overwrite=False, verify_names=True
        )

        assert model.scaling_factor[model.v[1]] == 5
        assert model.scaling_factor[model.v[2]] == 10
        assert model.scaling_factor[model.v[3]] == 15
        assert model.scaling_factor[model.v[4]] == 20
        assert model.scaling_factor[model.c[1]] == 5
        assert model.scaling_factor[model.c[2]] == 25
        assert model.scaling_factor[model.c[3]] == 125
        assert model.scaling_factor[model.c[4]] == 625
        assert len(model.scaling_factor) == 8

        for bd in model.b.values():
            assert bd.scaling_factor[bd.v2] == 10
            assert len(bd.scaling_factor) == 1

        # Ensure we have not mutated original dict
        assert sdict == {
            "v[1]": 50,
            "v[2]": 100,
            "c[1]": 50,
            "c[2]": 250,
            "block_name": "unknown",
        }

    @pytest.mark.unit
    def test_scaling_factors_from_dict_suffix_verify_fail(self, model):
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
                model.scaling_factor, sdict, overwrite=True, verify_names=True
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
    def test_scaling_factors_from_dict_block_data(self, model):
        # Partial dict of scaling factors to ensure we only touch things in the dict
        sdict = {
            "v[1]": 50,
            "v[2]": 100,
            "c[1]": 50,
            "c[2]": 250,
            "block_name": "unknown",
        }

        scaling_factors_from_dict(model, sdict, overwrite=True, verify_names=True)

        assert model.scaling_factor[model.v[1]] == 50
        assert model.scaling_factor[model.v[2]] == 100
        assert model.scaling_factor[model.v[3]] == 15
        assert model.scaling_factor[model.v[4]] == 20
        assert model.scaling_factor[model.c[1]] == 50
        assert model.scaling_factor[model.c[2]] == 250
        assert model.scaling_factor[model.c[3]] == 125
        assert model.scaling_factor[model.c[4]] == 625
        assert len(model.scaling_factor) == 8

        for bd in model.b.values():
            assert bd.scaling_factor[bd.v2] == 10
            assert len(bd.scaling_factor) == 1

        # Ensure we have not mutated original dict
        assert sdict == {
            "v[1]": 50,
            "v[2]": 100,
            "c[1]": 50,
            "c[2]": 250,
            "block_name": "unknown",
        }

    @pytest.mark.unit
    def test_scaling_factors_from_dict_block_data_overwrite_false(self, model):
        # Partial dict of scaling factors to ensure we only touch things in the dict
        sdict = {
            "v[1]": 50,
            "v[2]": 100,
            "c[1]": 50,
            "c[2]": 250,
            "block_name": "unknown",
        }

        scaling_factors_from_dict(model, sdict, overwrite=False, verify_names=True)

        assert model.scaling_factor[model.v[1]] == 5
        assert model.scaling_factor[model.v[2]] == 10
        assert model.scaling_factor[model.v[3]] == 15
        assert model.scaling_factor[model.v[4]] == 20
        assert model.scaling_factor[model.c[1]] == 5
        assert model.scaling_factor[model.c[2]] == 25
        assert model.scaling_factor[model.c[3]] == 125
        assert model.scaling_factor[model.c[4]] == 625
        assert len(model.scaling_factor) == 8

        for bd in model.b.values():
            assert bd.scaling_factor[bd.v2] == 10
            assert len(bd.scaling_factor) == 1

        # Ensure we have not mutated original dict
        assert sdict == {
            "v[1]": 50,
            "v[2]": 100,
            "c[1]": 50,
            "c[2]": 250,
            "block_name": "unknown",
        }

    @pytest.mark.unit
    def test_scaling_factors_from_dict_block_data_verify_fail(self, model):
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
                "Block name (unknown) does not match that recorded in json_dict (foo)"
            ),
        ):
            scaling_factors_from_dict(model, sdict, overwrite=True, verify_names=True)

        # Ensure we have not mutated original dict
        assert sdict == {
            "v[1]": 50,
            "v[2]": 100,
            "c[1]": 50,
            "c[2]": 250,
            "block_name": "foo",
        }

    @pytest.mark.unit
    def test_scaling_factors_from_dict_indexed_block(self, model):
        # Partial dict of scaling factors to ensure we only touch things in the dict
        sdict = {
            "block_datas": {
                "b[1]": {
                    "v2": 42,
                },
                "b[2]": {
                    "v2": 42,
                },
            },
            "block_name": "b",
        }

        scaling_factors_from_dict(model.b, sdict, overwrite=True, verify_names=True)

        assert model.scaling_factor[model.v[1]] == 5
        assert model.scaling_factor[model.v[2]] == 10
        assert model.scaling_factor[model.v[3]] == 15
        assert model.scaling_factor[model.v[4]] == 20
        assert model.scaling_factor[model.c[1]] == 5
        assert model.scaling_factor[model.c[2]] == 25
        assert model.scaling_factor[model.c[3]] == 125
        assert model.scaling_factor[model.c[4]] == 625
        assert len(model.scaling_factor) == 8

        for k in [1, 2]:
            assert model.b[k].scaling_factor[model.b[k].v2] == 42
            assert len(model.b[k].scaling_factor) == 1
        for k in [3, 4]:
            assert model.b[k].scaling_factor[model.b[k].v2] == 10
            assert len(model.b[k].scaling_factor) == 1

        # Ensure we have not mutated original dict
        assert sdict == {
            "block_datas": {
                "b[1]": {
                    "v2": 42,
                },
                "b[2]": {
                    "v2": 42,
                },
            },
            "block_name": "b",
        }

    @pytest.mark.unit
    def test_scaling_factors_from_dict_indexed_block_overwrite_false(self, model):
        # Partial dict of scaling factors to ensure we only touch things in the dict
        sdict = {
            "block_datas": {
                "b[1]": {
                    "v2": 42,
                },
                "b[2]": {
                    "v2": 42,
                },
            },
            "block_name": "b",
        }

        scaling_factors_from_dict(model.b, sdict, overwrite=False, verify_names=True)

        assert model.scaling_factor[model.v[1]] == 5
        assert model.scaling_factor[model.v[2]] == 10
        assert model.scaling_factor[model.v[3]] == 15
        assert model.scaling_factor[model.v[4]] == 20
        assert model.scaling_factor[model.c[1]] == 5
        assert model.scaling_factor[model.c[2]] == 25
        assert model.scaling_factor[model.c[3]] == 125
        assert model.scaling_factor[model.c[4]] == 625
        assert len(model.scaling_factor) == 8

        for k in [1, 2, 3, 4]:
            assert model.b[k].scaling_factor[model.b[k].v2] == 10
            assert len(model.b[k].scaling_factor) == 1

        # Ensure we have not mutated original dict
        assert sdict == {
            "block_datas": {
                "b[1]": {
                    "v2": 42,
                },
                "b[2]": {
                    "v2": 42,
                },
            },
            "block_name": "b",
        }

    @pytest.mark.unit
    def test_scaling_factors_from_dict_verify_names_failure(self, model):
        # Partial dict of scaling factors to ensure we only touch things in the dict
        sdict = {
            "block_datas": {
                "b[1]": {
                    "v2": 42,
                },
                "b[2]": {
                    "v2": 42,
                },
            },
            "block_name": "foo",
        }

        with pytest.raises(
            ValueError,
            match=re.escape(
                "Block name (b) does not match that recorded in json_dict (foo)"
            ),
        ):
            scaling_factors_from_dict(model.b, sdict, overwrite=True, verify_names=True)

        # Ensure we have not mutated original dict
        assert sdict == {
            "block_datas": {
                "b[1]": {
                    "v2": 42,
                },
                "b[2]": {
                    "v2": 42,
                },
            },
            "block_name": "foo",
        }

    @pytest.mark.unit
    def test_scaling_factors_from_dict_invalid_component(self, model):
        # Partial dict of scaling factors to ensure we only touch things in the dict
        sdict = {
            "block_datas": {
                "b[1]": {
                    "v2": 42,
                },
                "b[2]": {
                    "v2": 42,
                },
            },
            "block_name": "foo",
        }

        with pytest.raises(
            TypeError, match=re.escape("v is not an instance of a Block of Suffix.")
        ):
            scaling_factors_from_dict(model.v, sdict, overwrite=True, verify_names=True)

    @pytest.mark.unit
    def test_scaling_factors_to_json_file(self, model):
        temp_context = TempfileManager.new_context()
        tmpfile = temp_context.create_tempfile(suffix=".json")

        scaling_factors_to_json_file(model, tmpfile)

        with open(tmpfile, "r") as f:
            lines = f.read()
        f.close()

        print(lines)

        expected = """{
   "v[1]": 5,
   "c[1]": 5,
   "v[2]": 10,
   "c[2]": 25,
   "v[3]": 15,
   "c[3]": 125,
   "v[4]": 20,
   "c[4]": 625,
   "subblock_suffixes": {
      "b[1]": {
         "v2": 10,
         "subblock_suffixes": {}
      },
      "b[2]": {
         "v2": 10,
         "subblock_suffixes": {}
      },
      "b[3]": {
         "v2": 10,
         "subblock_suffixes": {}
      },
      "b[4]": {
         "v2": 10,
         "subblock_suffixes": {}
      }
   },
   "block_name": "unknown"
}"""

        assert lines == expected

        # Check for clean up
        temp_context.release(remove=True)
        assert not os.path.exists(tmpfile)

    @pytest.mark.unit
    def test_scaling_factors_from_json_file(self, model):
        fname = os.path.join(currdir, "load_scaling_factors.json")

        scaling_factors_from_json_file(model, fname, overwrite=True)

        assert model.scaling_factor[model.v[1]] == 50
        assert model.scaling_factor[model.v[2]] == 100
        assert model.scaling_factor[model.v[3]] == 150
        assert model.scaling_factor[model.v[4]] == 200
        assert model.scaling_factor[model.c[1]] == 50
        assert model.scaling_factor[model.c[2]] == 250
        assert model.scaling_factor[model.c[3]] == 1250
        assert model.scaling_factor[model.c[4]] == 6250
        assert len(model.scaling_factor) == 8

        for k in [1, 2, 3, 4]:
            assert model.b[k].scaling_factor[model.b[k].v2] == 100
            assert len(model.b[k].scaling_factor) == 1


class TestGetScalingFactor:
    @pytest.mark.unit
    def test_get_scaling_factor_block(self):
        m = ConcreteModel()

        with pytest.raises(TypeError, match="Blocks cannot have scaling factors."):
            get_scaling_factor(m)

    @pytest.mark.unit
    def test_get_scaling_factor(self):
        m = ConcreteModel()
        m.v = Var()

        m.scaling_factor = Suffix(direction=Suffix.EXPORT)
        m.scaling_factor[m.v] = 10

        assert get_scaling_factor(m.v) == 10

    @pytest.mark.unit
    def test_get_scaling_factor_none(self):
        m = ConcreteModel()
        m.v = Var()

        m.scaling_factor = Suffix(direction=Suffix.EXPORT)

        assert get_scaling_factor(m.v) is None

    @pytest.mark.unit
    def test_get_scaling_factor_no_suffix(self):
        m = ConcreteModel()
        m.v = Var()

        assert get_scaling_factor(m.v) is None


class TestSetScalingFactor:
    @pytest.mark.unit
    def test_set_scaling_factor(self):
        m = ConcreteModel()
        m.scaling_factor = Suffix(direction=Suffix.EXPORT)

        m.v = Var()

        set_scaling_factor(m.v, 42)

        assert m.scaling_factor[m.v] == 42

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

        assert "Created new scaling suffix for unknown" in caplog.text

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
                "scaling factor for v is negative (-42.0). "
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
                "scaling factor for v is zero. "
                "Scaling factors must be strictly positive."
            ),
        ):
            set_scaling_factor(m.v, 0)

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


class TestDelScalingFactor:
    @pytest.mark.unit
    def test_del_scaling_factor(self):
        m = ConcreteModel()
        m.v = Var()

        m.scaling_factor = Suffix(direction=Suffix.EXPORT)
        m.scaling_factor[m.v] = 10

        del_scaling_factor(m.v)

        assert len(m.scaling_factor) == 0

    @pytest.mark.unit
    def test_del_scaling_factor_not_present(self, caplog):
        caplog.set_level(
            idaeslog.DEBUG,
            logger="idaes",
        )

        m = ConcreteModel()
        m.v = Var()

        m.scaling_factor = Suffix(direction=Suffix.EXPORT)

        del_scaling_factor(m.v)

        assert len(m.scaling_factor) == 0

    @pytest.mark.unit
    def test_del_scaling_factor_delete_empty(self, caplog):
        caplog.set_level(
            idaeslog.DEBUG,
            logger="idaes",
        )

        m = ConcreteModel()
        m.v = Var()

        m.scaling_factor = Suffix(direction=Suffix.EXPORT)
        m.scaling_factor[m.v] = 10

        del_scaling_factor(m.v, delete_empty_suffix=True)

        assert not hasattr(m, "scaling_factor")

        assert "Deleting empty scaling suffix from unknown" in caplog.text


class TestReportScalingFactors:
    @pytest.fixture
    def model(self):
        m = ConcreteModel()
        m.s = Set(initialize=[1, 2, 3, 4])

        m.v = Var(m.s)

        @m.Constraint(m.s)
        def c(b, i):
            return b.v[i] == i

        m.b = Block(m.s)

        for i, bd in m.b.items():
            bd.v2 = Var()

        # Need to check all possible behaviours
        # Set values for half the variables (indexes 1 and 3)
        for i in [1, 3]:
            m.v[i].set_value(42)
            m.b[i].v2.set_value(42)

        # Set scaling factors for half the components (indexed 1 and 2)
        m.scaling_factor = Suffix(direction=Suffix.EXPORT)
        for i in [1, 2]:
            m.scaling_factor[m.v[i]] = 5 * i
            m.scaling_factor[m.c[i]] = 5**i
            m.b[i].scaling_factor = Suffix(direction=Suffix.EXPORT)
            m.b[i].scaling_factor[m.b[i].v2] = 10

        return m

    @pytest.mark.unit
    def test_report_scaling_factors_all(self, model):
        stream = StringIO()

        report_scaling_factors(model, descend_into=True, stream=stream)

        expected = """Scaling Factors for unknown

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
"""

        print(stream.getvalue())
        assert stream.getvalue() == expected

    @pytest.mark.unit
    def test_report_scaling_factors_descend_false(self, model):
        stream = StringIO()

        report_scaling_factors(model, descend_into=False, stream=stream)

        expected = """Scaling Factors for unknown

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
"""

        assert stream.getvalue() == expected

    @pytest.mark.unit
    def test_report_scaling_factors_vars_only(self, model):
        stream = StringIO()

        report_scaling_factors(model, descend_into=True, ctype=Var, stream=stream)

        expected = """Scaling Factors for unknown

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

        expected = """Scaling Factors for unknown

Constraint    Scaling Factor
c[1]          5.000E+00
c[2]          2.500E+01
c[3]          None
c[4]          None
"""

        assert stream.getvalue() == expected

    @pytest.mark.unit
    def test_report_scaling_factors_indexed_block(self, model):
        stream = StringIO()

        report_scaling_factors(model.b, descend_into=True, stream=stream)

        expected = """Scaling Factors for b

Variable    Scaling Factor    Value        Scaled Value
b[1].v2     1.000E+01         4.200E+01    4.200E+02
b[2].v2     1.000E+01         None         None
b[3].v2     None              4.200E+01    4.200E+01
b[4].v2     None              None         None
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
            match="report_scaling_factors only supports None, Var or Constraint for argument ctype: "
            "received foo.",
        ):
            report_scaling_factors(model, descend_into=True, stream=stream, ctype="foo")
