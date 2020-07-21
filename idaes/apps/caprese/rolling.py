from collections import OrderedDict
import bisect

class Series(list):
    """
    Class used to represent data series
    """
    def __init__(self, data=[], time=[], name=None, tolerance=1e-8):
        """
        """
        self.name = name
        self.tolerance = tolerance
        if len(time) != len(data):
            raise ValueError(
                'time and data lists must have same length')
        self.time = self.validate_time(time)
        super().__init__(data)

    def validate_time(self, time):
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

    def is_valid_append(self, t):
        if not self:
            return True
        t_last = self.time[-1]
        return (t - t_last >= 2*tolerance)

    def is_consistent_with_last_value(self, val):
        last = self[-1]
        # Equality may or may not be the right check here. 
        # Should be overridden if necessary.
        return last == val

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
        if not self.is_valid_append(t):
            raise ValueError(
                'Appended time values must at least 2*tolerance later than '
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


class History(OrderedDict):
    """
    """
    def __init__(self, init=OrderedDict(), time=None, name=None, tolerance=1e-8):
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

        self.name = name
        self.tolerance = tolerance
        self.time = time

        super().__init__(init)
