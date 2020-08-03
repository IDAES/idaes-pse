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
            if t1 - t0 <= 2*self.tolerance:
                raise ValueError(
                    'Time points must be increasing and separated by at '
                    'least twice the tolerance')
        return time

    def is_within_bounds(self, t):
        """
        Checks whether point t is in the interval bounded by the first and
        last time points.
        """
        if not self:
            raise ValueError('List is empty. No bounds exist.')
        lower = self[0]
        upper = self[-1]
        tolerance = self.tolerance
        return lower - tolerance <= t and t <= upper + tolerance

    def is_valid_append(self, t):
        """
        Checks whether point t is greater than the last existing time point
        by an amount consistent with the tolerance.
        """
        if not self:
            return True
        t_last = self[-1]
        return (t - t_last > 2*tolerance)

    def append(self, t):
        """
        Appends a valid time point.
        """
        if not self.is_valid_append(t):
            raise ValueError(
                'Appended time values must be more than 2*tolerance later '
                'than current values.')
        super().append(t)

    def validate_extend(self, tpoints):
        """
        """
        if not tpoints:
            return tpoints
        tolerance = self.tolerance
        t_last = self[-1]
        t_new = tpoints[0]
        if t_new < t_last - tolerance:
            raise ValueError(
                'First new point must not be earlier than last existing point.')
        elif t_new <= t_last + tolerance:
            # Do not want to have duplicates in time list, but want to allow
            # extending with an interval where the boundaries overlap.
            return tpoints[1:]
        elif t_new <= t_last + 2*tolerance:
            raise ValueError(
                'Distinct time points must be separated by more than twice the '
                'tolerance.')
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

    def binary_search(self, t):
        """
        """
        # TODO: clean up this implementation
        tolerance = self.tolerance

        i = bisect.bisect_left(self, t)
        t_ge = self[i]
        if t_ge - t < tolerance:
            return True, i

        if i == 0:
            return False, i

        t_l = self[i-1]
        if t - t_l < tolerance:
            return True, i-1

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
                raise ValueError(
                    'data lists and time list must have same length')

        self.name = name
        self.time = TimeList(time, tolerance=tolerance)

        super().__init__(data)

    def dim(self):
        return len(self.data)

    def __len__(self):
        return len(self.time)

    # TODO in this class:
    # - check for consistency
    # - append/extend values, including checks in TimeList
    # - interpolate. Then interpolate each component of the vector?

    # This is just my SynchronousHistory class.
    # Same question, again: Should this class contain functions for consistency
    # and interpolation.
    # Question really is: Is this a vector of time-indexed components, or is it
    # a time-indexed vector. I think the former.
    # Know (enforce) that I will do append/extend as a block. 
    # => consistency check can be done as a block as well
    #    I will not mix variables that require a different check for consistency.
    # Interpolation may not even make sense for, e.g. inputs
    # (Interpolation should still be well defined.)
    # Want:
    # - find_nearest_index() function that I can call with a tight tolerance
    #   ^ This should be a method of time_list
    # - interpolate() function that uses find_nearest_index. 
    #   Vector-valued, for now
    # - get_value() function that, by default, calls interpolate()
    #   ^ should be overridden for inputs, objectives
    #   => will be vector-valued
    #      could just have an interpolate function that calls a Series 
    #      interpolate function...
    #      or an interpolate utility function...
            
    def consistent(self, target):
        if len(self) == 0:
            return True
        last = [series[-1] for series in self.values()]
        return all(new == old for new, old in zip(target, last))

    def append(self, t, data):
        self.time.append(t)
        self.append(data)

    def extend(self, tpoints, data):
        """
        data is a matrix.
        """
        #       How should I infer whether to check for consistency?
        #       Just compare last existing and first new time values
        # TODO: Should I have an option that allows the user to violate 
        #       consistency?
        try:
            data = list(data.values())
        except AttributeError as ae:
            if 'values' not in str(ae):
                raise ae
        tolerance = self.time.tolerance
        if len(self) != 0:
            tlast = self.time[-1]
            t0 = tpoints[0]
            if tlast-tolerance <= t0 and t0 <= tlast+tolerance:
                data0 = [series[0] for series in data]
                if not self.consistent(data0):
                    raise ValueError(
                        'Tried to extend with series that overlapped at time '
                        'point %s, but the series data was not consistent '
                        'with pre-existing data.' 
                        % t0)
            tpoints = tpoints[1:]
            data = [series[1:] for series in data]
        self.time.extend(tpoints)
        for series, new_data in zip(self.values(), data):
            series.extend(new_data)


class Series(list):
    """
    Class used to represent data series.
    """
    # Do I intend for this class to be user-facing?
    # I don't think so.
    # Because this class will have no logic for maintaining
    # its TimeList. 
    def __init__(self, data=None, time=None, name=None):
        """
        Args:
            data: List of data points
            time: TimeList object
            name: Name of this data series
        """
        if data is None:
            data = []
        if time is None:
            time = TimeList()

        self.name = name
        self.time = time
        data = list(data)
        if len(time) != len(data):
            raise ValueError(
                'time and data lists must have the same length')
        super().__init__(data)

    def consistent(self, val):
        """
        Used to determine if a new value is consistent with the last
        existing data value. May need to be overridden.
        """
        if not self:
            return True
        last = self[-1]
        return val == last
    
    def interpolate(self, t):
        """
        """
        found, i = self.time.binary_search(t)
        if found:
            return self[i]

        if i == 0 or i == len(self) + 1:
            raise ValueError(
                'Value to interpolate is out of bounds')

        time = self.time
        t_l = time[i-1]
        t_g = time[i]
        val_l = self[i-1]
        val_g = self[i]

        return val_l + (val_g - val_l)/(t_g - t_l)*(t - t_l)


class History(OrderedDict):
    """
    """
    def __init__(self, init=None, name=None, tolerance=0.0):
        """
        Args:
            init: OrderedDict: cuid -> list of time, value tuples
        """
        if init is None:
            init = OrderedDict()
        self.name = name
        for cuid, data in init.items():
            time_list = TimeList([t for t, val in data], tolerance)
            data_list = [val for t, val in data]

            # This should be some subclass of Series that manages its own
            # time_list. TimeSeries or something...
            init[cuid] = Series(
                    data=data_list,
                    time=time_list,
                    )
        super().__init__(init)

    # This class may implement its own append/extend methods, but they will
    # just be wrappers around the "TimeSeries" methods.


class SynchronizedHistory(OrderedDict):
    """
    """
    def __init__(self, init=OrderedDict(), time=TimeList(), name=None):
        """
        Args:
            init: OrderedDict: cuid -> list of data points
            time: TimeList object. Must be consistent with data lists in init.
        """
        # Tolerance should be provided in the TimeList that is used to 
        # construct.
        self.name = name
        self.time = time
        n_time = len(time)
        for cuid, data in init.items():
            if len(data) != n_time:
                raise ValueError(
                    'Time and data lists must have same length')
            init[cuid] = Series(
                    data=data,
                    time=self.time,
                    )
        super().__init__(init)
            

class _History(OrderedDict):
    """
    """
    def __init__(self, init={}, time=None, name=None, tolerance=1e-8):
        """
        init should map: cuid -> list of data points
        """
        # TODO: make a SynchronousHistory subclass of History
        if time is None:
            self.synchronous = False
        else:
            self.synchronous = True
            # Cast to list so instances may be constructed with a ContinuousSet
            time = list(time)
            time_list = time

        for cuid, data in init.items():
            if not self.synchronous:
                # In this case, data must be a list of time-point, value tuples
                # This is necessary as each series must contain its own time
                # list.
                time_list = [t for t, val in data]
                data = [val for t, val in data]

            init[cuid] = Series(
                    data=data,
                    time=time_list,
                    name=str(cuid),
                    tolerance=tolerance)

        self.dim = len(init)
        self.name = name
        self.tolerance = tolerance
        self.time = time

        super().__init__(init)

    def extend(self, data_list, time_list=None):
        if time_list is None and self.synchronous:
            raise ValueError(
                'A single time list must be passed in for synchronous data')
        elif time_list is not None and not self.synchronous:
            raise ValueError(
                'Time list must be omitted for asynchronous data')
        if len(data_list) != self.dim:
            raise ValueError(
                "Must extend with same number of data series as History's dimension")

        if self.synchonous:
            t_last = self.time[-1]
            t0 = time_list[0]

            # For now, time will always be offset from last existing value.
            # Unsure whether I should (add an option to) relax this.
            new_time = [t_last + (t-t0) for t in time_list]
            try:
                # Extending from time-indexed Pyomo variables
                new_data_list = [[container[t].value for t in time_list]
                        for container in data_list]
            except AttributeError as ae:
                if 'value' not in str(ae):
                    raise ae
                # Extending from lists of values
                new_data_list = []
                for data in data_list:
                    if len(data) != len(time_list):
                        raise ValueError(
                            'Data must have same length as time list for '
                            'synchronous data.')
                    new_data_list.append(data)
        
        for data, series in zip(data_list, self.values()):
            series.extend


class _Series(list):
    """
    Old class used to represent data series
    This class is basically a combination of Series and TimeList,
    and is mostly deprecated.
    """
    def __init__(self, data=[], time=None, name=None, tolerance=1e-8):
        """
        Args:
            data: List of data points
            time: List of time points
            name: Name of this data series
            tolerance: Tolerance used to check "equality" between floating
                       point time values. Values with a difference less than
                       or equal to this tolerance are "equal."
        """
        self.name = name
        self.tolerance = tolerance
        if len(time) != len(data):
            raise ValueError(
                'time and data lists must have same length')
        self.time = self.validate_time(time)
        super().__init__(data)


    def is_consistent_with_last_value(self, val):
        """
        Checks whether a value is equal to the last existing
        data value. This may need to be overridden.
        """
        last = self[-1]
        return val == last

    # TODO: Convert to validate_extend function so I can give more specific
    # error messages
    def is_valid_extend(self, tlist, val_0):
        tlist = self.validate_time(tlist)
        if not self:
            return True
        t_last = self.time[-1]
        t0 = tlist[0]
        if t0 - t_last <= -self.tolerance:
            # An earlier time point is never valid
            return False

        elif abs(t0 - t_last) < self.tolerance:
            # A time point at the current end point is valid if 
            # the corresponding value is consistent
            return self.is_consistent_with_last_value(val_0)

        elif t0 - t_last <= 2*self.tolerance:
            # A later time point is valid if tolerance is a small enough radius
            # to separate it from the last existing time point.
            return False

        else:
            # t0 - t_last > 2*self.tolerance
            return True

    def append(self, t, val):
        """
        
        """
        if not self.is_valid_append(t):
            raise ValueError(
                'Appended time values must more than 2*tolerance later than '
                'current values')
        self.time.append(t)
        super().append(val)

    def extend(self, tlist, vals):
        val_0 = vals[0]
        if not self.is_valid_extend(tlist, val_0):
            raise ValueError(
                'Invalid extend')
        self.time.extend(tlist)
        super().extend(vals)

    def binary_search(self, t):
        time = self.time
        tolerance = self.tolerance

        i = bisect.bisect_left(time, t)
        t_ge = time[i]
        if t_ge - t < tolerance:
            return True, i

        if i == 0:
            return False, i

        t_l = time[i-1]
        if t - t_l < tolerance:
            return True, i-1

        return False, i

    def get(self, t):
        found, i = self.binary_search(t)
        if found:
            return self[i]
        else:
            return None

    def interpolate(self, t):
        found, i = self.binary_search(t)
        if found:
            return self[i]

        if i == 0 or i == len(self) + 1:
            raise ValueError(
                'Value to interpolate is out of bounds')

        time = self.time
        t_l = time[i-1]
        t_g = time[i]
        val_l = self[i-1]
        val_g = self[i]

        return val_l + (val_g - val_l)/(t_g - t_l)*(t - t_l)


class InputSeries(Series):
    """
    For consistency with the backwards discretization schemes most commonly
    used, the time corresponding to a value in an input history is the last
    time point at which that input is applied.
    """
    # TODO: revisit this once Series implementation is more stable
    def get(self, t):
        if not self.is_within_bounds(t):
            return None
        found, i = self.binary_search(t)
        return self[i]


