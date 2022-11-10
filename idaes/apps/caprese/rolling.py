# -*- coding: utf-8 -*-
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
""" Classes used for rolling-horizon simulation.
"""

from pyomo.core.base.componentuid import ComponentUID
from collections import OrderedDict
import bisect


class TimeList(list):
    """
    Maintains a list of time points in sorted order, ensures time points
    are unique, and implements a binary search method.
    """

    # TODO: This class could inherit from ContinuousSet.
    # That would be convenient for a lot of reasons, one of which is that
    # I could use the find_nearest_index method
    def __init__(self, time=None, tolerance=0.0):
        """
        A tolerance is required because I want to ensure time points
        are unique.
        """
        if time is None:
            time = []
        time = list(time)
        self.tolerance = tolerance
        time = self.validate_time(time)
        super().__init__(time)

    def validate_time(self, time):
        """
        Make sure that list of time points used to initialize is
        increasing, and that points are separated by more than twice
        the tolerance, as is required for time points to be "unique-
        within-tolerance."
        """
        # Could sort time before doing this...
        tolerance = self.tolerance
        len_time = len(time)
        if len_time == 0 or len_time == 1:
            return time
        t0 = time[0]
        for t1 in time[1:]:
            if t1 - t0 <= 2 * self.tolerance:
                raise ValueError(
                    "Time points must be increasing and separated by more "
                    "than twice the tolerance."
                )
        return time

    def is_within_bounds(self, t):
        """
        Checks whether point t is in the interval bounded by the first and
        last time points.
        """
        if not self:
            raise ValueError("List is empty. No bounds exist.")
        lower = self[0]
        upper = self[-1]
        tol = self.tolerance
        return lower - tol <= t and t <= upper + tol

    def is_valid_append(self, t):
        """
        Checks whether point t is greater than the last existing time point
        by an amount consistent with the tolerance.
        """
        if not self:
            return True
        t_last = self[-1]
        tol = self.tolerance
        return t - t_last > 2 * tol

    def append(self, t):
        """
        Appends a valid time point.
        """
        if not self.is_valid_append(t):
            raise ValueError(
                "Appended time values must be more than 2*tolerance later "
                "than current values."
            )
        super().append(t)

    def validate_extend(self, tpoints):
        """
        Checks whether new time points overlap with current time points
        at more than a single point. If there is no overlap, the
        separation between current and new points must be more than
        twice the tolerance.
        """
        if not tpoints:
            # May assume input tpoints is a list.
            return tpoints
        tolerance = self.tolerance
        t_last = self[-1]
        t_new = tpoints[0]
        if t_new < t_last - tolerance:
            raise ValueError(
                "First new point must not be earlier than last existing point."
            )
        elif t_new <= t_last + tolerance:
            # Do not want to have duplicates in time list, but want to allow
            # extending with an interval where the boundaries overlap.
            return tpoints[1:]
        elif t_new <= t_last + 2 * tolerance:
            raise ValueError(
                "Distinct time points must be separated by more than twice the "
                "tolerance."
            )
        else:
            # t_new > t_last + 2*tolerance
            return tpoints

    def extend(self, tpoints):
        """
        Extends by a valid iterable of time points.
        """
        tpoints = list(tpoints)
        tpoints = self.validate_time(tpoints)
        tpoints = self.validate_extend(tpoints)
        super().extend(tpoints)

    def find_nearest_index(self, target):
        """
        Returns the index of the nearest point in self.
        """
        # This is a copy of the same method of ContinuousSet.
        #
        tol = self.tolerance
        lo = 0
        hi = len(self)
        i = bisect.bisect_right(self, target, lo=lo, hi=hi)
        # i is index at which target should be inserted if it is
        # to be right of any equal components.

        if i == lo:
            # target is less than every entry of the set
            nearest_index = i
            delta = self[nearest_index] - target
        elif i == hi:
            nearest_index = i - 1
            delta = target - self[nearest_index]
        else:
            delta, nearest_index = min((abs(target - self[_]), _) for _ in [i - 1, i])

        if delta > tol:
            return None
        return nearest_index

    def binary_search(self, t):
        """ """
        # TODO: clean up this implementation
        tolerance = self.tolerance

        i = bisect.bisect_left(self, t)
        t_ge = self[i]
        if t_ge - t < tolerance:
            return True, i

        if i == 0:
            return False, i

        t_l = self[i - 1]
        if t - t_l < tolerance:
            return True, i - 1

        return False, i


class VectorSeries(OrderedDict):
    """
    A time-indexed vector.
    """

    # This class's purpose is suspiciously similar to that of IndexedComponent.
    # Could have an OrderedDict of time-indexed variables/Parameters that I
    # continually update along with the underlying time set.
    def __init__(self, data=None, time=None, name=None, tolerance=0.0):
        """
        Args:
            data: OrderedDict mapping cuids to lists of values
            time: list of time points corresponding to data values
            name: name of this vector
        """
        if data is None:
            data = OrderedDict()
        if time is None:
            time = []

        len_time = len(time)
        for val_list in data.values():
            if len(val_list) != len_time:
                raise ValueError("data lists and time list must have same length")

        self.name = name
        self.time = TimeList(time, tolerance=tolerance)

        super().__init__(data)

    def dim(self):
        """This is the dimension of the vector that is indexed by time."""
        return len(self.keys())

    def __len__(self):
        """This is the length of the time series, each element of which
        is a vector.
        """
        return len(self.time)

    def consistent(self, target):
        """
        target is a vector
        """
        if len(target) != self.dim():
            return False
        if len(self) == 0:
            return True
        last = [series[-1] for series in self.values()]
        return all(new == old for new, old in zip(target, last))

    def consistent_dimension(self, target):
        return len(target) == self.dim()

    def validate_dimension(self, target):
        if not self.consistent_dimension(target):
            raise ValueError(
                "Tried to validate a vector with inconsistent dimension. "
                "Expected %s, got %s." % (self.dim(), len(target))
            )
        return target

    def append(self, t, data):
        """
        data is a vector.
        """
        data = self.validate_dimension(data)
        self.time.append(t)
        for series, val in zip(self.values(), data):
            series.append(val)

    def extend(self, tpoints, data):
        """
        data is a matrix.
        """
        #       How should I infer whether to check for consistency?
        #       Just compare last existing and first new time values
        # TODO: Should I have an option that allows the user to violate
        #       consistency?
        try:
            # This allows data to be an OrderedDict
            # or another VectorSeries.
            data = list(data.values())
        except AttributeError as ae:
            if "values" not in str(ae):
                # If this attribute error is caught, data should
                # behave like a list of lists.
                raise ae
        data = self.validate_dimension(data)
        tolerance = self.time.tolerance
        if len(self) != 0:
            tlast = self.time[-1]
            t0 = tpoints[0]
            if tlast - tolerance <= t0 and t0 <= tlast + tolerance:
                data0 = [series[0] for series in data]
                if not self.consistent(data0):
                    raise ValueError(
                        "Tried to extend with series that overlapped at time "
                        "point %s, but the series data was not consistent "
                        "with pre-existing data." % t0
                    )
                tpoints = tpoints[1:]
                data = [series[1:] for series in data]
        self.time.extend(tpoints)
        for series, new_data in zip(self.values(), data):
            series.extend(new_data)
