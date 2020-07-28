from collections import OrderedDict
import bisect

class Series(list):
    """
    Class used to represent data series
    """
    def __init__(self, data=[], time=[], name=None, tolerance=1e-8):
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

    def validate_time(self, time):
        """
        Make sure that list of time points used to initialize is 
        increasing, and that points are separated by more than twice
        the tolerance, as is required for time points to be "unique-
        within-tolerance."
        """
        len_time = len(time)
        if len_time == 0 or len_time == 1:
            return time
        t0 = time[0]
        for t1 in time[1:]:
            if t1 - t0 > 2*self.tolerance:
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
        lower = self.time[0]
        upper = self.time[-1]
        tolerance = self.tolerance
        return lower - tolerance <= t and t <= upper + tolerance

    def is_valid_append(self, t):
        """
        Checks whether point t is greater than the last existing time point
        by an amount consistent with the tolerance.
        """
        if not self:
            return True
        t_last = self.time[-1]
        return (t - t_last > 2*tolerance)

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
    def get(self, t):
        if not self.is_within_bounds(t):
            return None
        found, i = self.binary_search(t)
        return self[i]


class History(OrderedDict):
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
