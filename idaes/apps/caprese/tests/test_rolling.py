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
""" Tests for rolling horizon helper classes
"""

import pytest
import sys
from collections import OrderedDict
import pyomo.environ as pyo
from idaes.apps.caprese.rolling import *


class TestTimeList(object):
    @pytest.mark.unit
    def test_constructor(self):
        empty_time_list = TimeList()
        assert hasattr(empty_time_list, "tolerance")
        assert len(empty_time_list) == 0

        tol = 0.49
        good_time_list = TimeList(time=[1, 2, 3], tolerance=tol)
        assert good_time_list.tolerance == tol

        with pytest.raises(ValueError) as err:
            # Error message will be checked later
            bad_time_list = TimeList(time=[3, 1, 2], tolerance=tol)

    @pytest.mark.unit
    def test_validate_time(self):
        increasing_time = [1, 2, 3]
        nonincreasing_time = [3, 1, 2]

        good_tol = 0.4999
        bad_tol = 0.5

        time_list_1 = TimeList(tolerance=good_tol)
        time_list_2 = TimeList(tolerance=bad_tol)

        assert time_list_1.validate_time(increasing_time) == increasing_time
        with pytest.raises(ValueError) as err:
            time = time_list_1.validate_time(nonincreasing_time)
            msg = "Time points must be increasing"
            assert msg in str(err)

        with pytest.raises(ValueError) as err:
            time = time_list_2.validate_time(increasing_time)
            msg = " and separated by more than twice the tolerance."
            assert msg in str(err)

    @pytest.mark.unit
    def test_is_within_bounds(self):
        empty_time_list = TimeList()

        with pytest.raises(ValueError) as err:
            empty_time_list.is_within_bounds(0.0)
            msg = "List is empty"
            assert msg in str(err)

        time = [1, 2, 3]
        t_last = time[-1]
        tol = 0.1
        eps = sys.float_info.epsilon
        root_eps = eps**0.5
        time_list = TimeList(time, tolerance=tol)
        assert time_list.is_within_bounds(2.5)
        assert time_list.is_within_bounds(t_last + tol / 2)
        assert time_list.is_within_bounds(t_last + tol)
        assert not time_list.is_within_bounds(3.5)

        assert not time_list.is_within_bounds(3.1 + root_eps)

        # Of course, we may still violate our tolerance
        # by less than machine precision.
        assert time_list.is_within_bounds(3.1 + eps / 2)

    @pytest.mark.unit
    def test_is_valid_append(self):
        time = [1, 2, 3]
        t_last = time[-1]
        tol = 0.1
        root_eps = sys.float_info.epsilon**0.5
        empty_time_list = TimeList()
        time_list = TimeList(time, tolerance=tol)
        assert empty_time_list.is_valid_append(-1)
        assert not time_list.is_valid_append(0.0)
        assert not time_list.is_valid_append(t_last)
        assert not time_list.is_valid_append(t_last + tol)
        assert time_list.is_valid_append(t_last + 2 * tol + root_eps)

    @pytest.mark.unit
    def test_append(self):
        time = [1, 2, 3]
        t_last = time[-1]
        tol = 0.1
        root_eps = sys.float_info.epsilon**0.5
        time_list = TimeList(time, tolerance=tol)

        t_new = t_last + tol
        with pytest.raises(ValueError) as err:
            time_list.append(t_new)
            msg = "Appended time values must be"
            assert msg in str(err)

        t_new = t_last + 2 * tol + root_eps
        time_list.append(t_new)
        assert time_list == time + [t_new]

    @pytest.mark.unit
    def test_validate_extend(self):
        time = [1, 2, 3]
        t_last = time[-1]
        tol = 0.1
        root_eps = sys.float_info.epsilon**0.5
        time_list = TimeList(time, tolerance=tol)

        assert time_list.validate_extend([]) == []

        new_points = [2.5, 4]
        with pytest.raises(ValueError) as err:
            time_list.validate_extend(new_points)
            msg = "First new point must not be earlier"
            assert msg in str(err)

        new_points = [t_last, 4]
        assert time_list.validate_extend(new_points) == [4]
        new_points = [t_last + tol / 2, 4]
        assert time_list.validate_extend(new_points) == [4]

        new_points = [t_last + 3 * tol / 2, 4]
        with pytest.raises(ValueError) as err:
            time_list.validate_extend(new_points)
            msg = "must be separated by more than twice the tolerance"
            assert msg in str(err)

        new_points = [3.5, 4]
        assert time_list.validate_extend(new_points) == new_points

    @pytest.mark.unit
    def test_extend(self):
        time = [1, 2, 3]
        t_last = time[-1]
        tol = 0.1
        time_list = TimeList(time, tolerance=tol)

        new_points = [4, 5]
        time_list.extend(new_points)
        assert time_list == [1, 2, 3, 4, 5]

        new_points = [5, 6, 7]
        time_list.extend(new_points)
        assert time_list == [1, 2, 3, 4, 5, 6, 7]

    @pytest.mark.unit
    def test_find_nearest_index(self):
        time = [1, 2, 3, 4, 5]
        tol = 0.1
        time_list = TimeList(time, tolerance=tol)

        assert time_list.find_nearest_index(2.5) is None
        assert time_list.find_nearest_index(2.05) == 1
        assert time_list.find_nearest_index(1.0) == 0
        assert time_list.find_nearest_index(5.05) == 4


def make_model():
    m = pyo.ConcreteModel()
    m.time = pyo.Set(initialize=range(11))
    m.space = pyo.Set(initialize=range(3))
    m.comp = pyo.Set(initialize=["a", "b"])

    m.v1 = pyo.Var(m.time, m.space)
    m.v2 = pyo.Var(m.time, m.space, m.comp)

    @m.Block(m.time, m.space)
    def b(b, t, x):
        b.v3 = pyo.Var()
        b.v4 = pyo.Var(m.comp)

    m.v1_refs = [pyo.Reference(m.v1[:, x]) for x in m.space]
    m.v2a_refs = [pyo.Reference(m.v2[:, x, "a"]) for x in m.space]
    m.v2b_refs = [pyo.Reference(m.v2[:, x, "b"]) for x in m.space]
    m.v3_refs = [pyo.Reference(m.b[:, x].v3) for x in m.space]
    m.v4a_refs = [pyo.Reference(m.b[:, x].v4["a"]) for x in m.space]
    m.v4b_refs = [pyo.Reference(m.b[:, x].v4["a"]) for x in m.space]

    for t in m.time:
        for ref in m.v1_refs:
            ref[t] = 1.0 * t
        for ref in m.v2a_refs + m.v2b_refs:
            ref[t] = 2.0 * t
        for ref in m.v3_refs:
            ref[t] = 3.0 * t
        for ref in m.v4a_refs + m.v4b_refs:
            ref[t] = 4.0 * t

    return m


class TestVectorSeries(object):
    @pytest.mark.unit
    def test_constructor(self):
        # VectorSeries is constructed from an OrderedDict with data
        # series (lists) as values. Keys are meant to be identifiers
        # of Pyomo model components, but could technically be anything.
        # To illustrate their intended use, tests use CUIDs as keys.
        m = make_model()
        data = OrderedDict(
            [
                (pyo.ComponentUID(_slice.referent), [_slice[t].value for t in m.time])
                for _slice in m.v1_refs
            ]
        )

        name = "v1"
        tol = 0.1
        history = VectorSeries(data, time=list(m.time), name=name, tolerance=tol)
        assert history.name == name
        assert type(history.time) is TimeList
        assert history.time == m.time
        assert history.time.tolerance == tol

        for _slice in m.v1_refs:
            assert pyo.ComponentUID(_slice.referent) in history

        empty_series = VectorSeries()
        assert empty_series == OrderedDict()

    @pytest.mark.unit
    def test_dim(self):
        # This test also covers consistent_dimension and
        # validate_dimension.
        m = make_model()
        data = OrderedDict(
            [
                (pyo.ComponentUID(_slice.referent), [_slice[t].value for t in m.time])
                for _slice in m.v2a_refs + m.v2b_refs
            ]
        )
        name = "v2"
        tol = 0.1
        history = VectorSeries(data, time=list(m.time), name=name, tolerance=tol)
        assert history.dim() == len(m.space * m.comp)
        assert VectorSeries().dim() == 0

        vals = [1 for _ in m.space * m.comp]
        assert history.consistent_dimension(vals)
        assert not history.consistent_dimension([1])
        assert history.validate_dimension(vals) == vals
        with pytest.raises(ValueError) as err:
            history.validate_dimension([1])
            assert "inconsistent dimension" in str(err)

    @pytest.mark.unit
    def test_len(self):
        m = make_model()
        data = OrderedDict(
            [
                (pyo.ComponentUID(_slice.referent), [_slice[t].value for t in m.time])
                for _slice in m.v3_refs
            ]
        )
        name = "v3"
        tol = 0.1
        history = VectorSeries(data, time=list(m.time), name=name, tolerance=tol)
        assert len(history) == len(m.time)
        assert len(VectorSeries()) == 0

    @pytest.mark.unit
    def test_consistent(self):
        m = make_model()
        data = OrderedDict(
            [
                (pyo.ComponentUID(_slice.referent), [_slice[t].value for t in m.time])
                for _slice in m.v3_refs
            ]
        )
        name = "v3"
        tol = 0.1
        history = VectorSeries(data, time=list(m.time), name=name, tolerance=tol)

        # These vals are consistent because this is how
        # I happened to initialize m.v3
        vals = [3.0 * m.time.last() for _ in m.space]
        bad_vals = [m.time.last() for _ in m.space]
        assert not history.consistent([])
        assert not history.consistent(bad_vals)
        assert history.consistent(vals)

        assert VectorSeries().consistent([])
        assert not VectorSeries().consistent([1])

        empty_data = OrderedDict(
            [(pyo.ComponentUID(_slice.referent), []) for _slice in m.v1_refs]
        )
        empty_series = VectorSeries(empty_data)
        vals = [1 for _ in m.space]
        assert not empty_series.consistent([])
        assert empty_series.consistent(vals)

    @pytest.mark.unit
    def test_append(self):
        m = make_model()
        data = OrderedDict(
            [
                (pyo.ComponentUID(_slice.referent), [_slice[t].value for t in m.time])
                for _slice in m.v3_refs
            ]
        )
        name = "v3"
        tol = 0.1
        history = VectorSeries(data, time=list(m.time), name=name, tolerance=tol)
        new_t = 11
        with pytest.raises(ValueError) as err:
            history.append(new_t, [1])
            assert "inconsistent dimension" in str(err)
        new_vals = [3 * new_t for _ in m.space]
        history.append(new_t, new_vals)
        last = [series[-1] for series in history.values()]
        assert last == new_vals
        assert history.time[-1] == new_t

        with pytest.raises(ValueError) as err:
            history.append(new_t, new_vals)
            msg = "Appended time values must be"
            assert msg in str(err)

    @pytest.mark.unit
    def test_extend(self):
        m = make_model()
        midpoint = len(m.time) // 2
        tlist = list(m.time)
        time_1 = tlist[0:midpoint]
        time_2 = tlist[midpoint:]

        data_1 = OrderedDict(
            [
                (pyo.ComponentUID(_slice.referent), [_slice[t].value for t in time_1])
                for _slice in m.v1_refs
            ]
        )
        data_2 = OrderedDict(
            [
                (pyo.ComponentUID(_slice.referent), [_slice[t].value for t in time_2])
                for _slice in m.v1_refs
            ]
        )
        tol = 0.1
        history_1 = VectorSeries(data_1, time_1, tolerance=tol)
        history_2 = VectorSeries(data_2, time_2, tolerance=tol)
        history_1.extend(history_2.time, history_2)

        vals = [1.0 * t for t in m.time]
        for series in history_1.values():
            assert series == vals
        assert history_1.time == m.time

        new_time = [10, 11, 13]
        new_vals = [1 * t for t in new_time]
        new_data = [new_vals for _ in history_1.values()]
        history_1.extend(new_time, new_data)

        time = list(m.time) + new_time
        time = list(OrderedDict.fromkeys(time))  # No duplicates
        vals = [1 * t for t in time]
        for series in history_1.values():
            assert series == vals
        assert history_1.time == time

        new_time = [13, 14, 15]
        new_vals = [2 * t for t in new_time]
        new_data = [new_vals for _ in history_1.values()]
        with pytest.raises(ValueError) as err:
            history_1.extend(new_time, new_data)
            msg = "data was not consistent"
            assert mst in str(err)

        with pytest.raises(ValueError) as err:
            # Assuming that history_2.dim() != 1 ...
            history_2.extend([], [[]])
            assert "inconsistent dimension" in str(err)
