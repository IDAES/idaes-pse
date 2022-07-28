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


def assert_disjoint_intervals(intervals):
    """
    This function takes intervals in the form of tuples and makes sure
    that they are disjoint.

    Arguments
    ---------
    intervals: iterable
        Iterable of tuples, each containing the low and high values of an
        interval.

    """
    intervals = list(sorted(intervals))
    for i, (lo, hi) in enumerate(intervals):
        if not lo <= hi:
            raise RuntimeError(
                "Lower endpoint of interval is higher than upper endpoint"
            )
        if i != 0:
            prev_lo, prev_hi = intervals[i - 1]
            if not prev_hi <= lo:
                raise RuntimeError(
                    "Intervals %s and %s are not disjoint"
                    % ((prev_lo, prev_hi), (lo, hi))
                )


def load_inputs_into_model(model, time, input_data, time_tol=0):
    """
    This function loads piecewise constant values into variables (or
    mutable parameters) of a model.

    Arguments
    ---------
    model: _BlockData
        Pyomo block containing the variables and parameters whose values
        will be set
    time: ContinuousSet
        Pyomo ContinuousSet corresponding to the piecewise constant intervals
    input_data: dict of dicts
        Maps variable names to dictionaries mapping 2-tuples to values.
        Each tuple contains the low and high endpoints of an interval
        on which the variable is to take the specified value.
    time_tol: float
        Optional. Tolerance within which the ContinuousSet will be searched
        for interval endpoints. The default is zero, i.e. the endpoints
        must be within the ContinuousSet exactly.

    """
    for cuid, inputs in input_data.items():
        var = model.find_component(cuid)
        if var is None:
            raise RuntimeError(
                "Could not find a variable on model %s with ComponentUID %s"
                % (model.name, cuid)
            )

        intervals = list(sorted(inputs.keys()))
        assert_disjoint_intervals(intervals)
        for i, interval in enumerate(intervals):
            idx0 = time.find_nearest_index(interval[0], tolerance=time_tol)
            idx1 = time.find_nearest_index(interval[1], tolerance=time_tol)
            if idx0 is None or idx1 is None:
                # One of the interval boundaries is not a valid time index
                # within tolerance. Skip this interval and move on.
                continue
            input_val = inputs[interval]
            idx_iter = range(idx0 + 1, idx1 + 1) if idx0 != idx1 else (idx0,)
            for idx in idx_iter:
                t = time.at(idx)
                var[t].set_value(input_val)


def interval_data_from_time_series(data, use_left_endpoint=False):
    """
    This function converts time series data to piecewise constant
    interval data. A series of N time points and values yields
    N-1 intervals. By default, each interval takes the value of
    its right endpoint.

    In:
        ([t0, ...], {str(cuid): [value0, ...], },)
    Out:
        {str(cuid): {(t0, t1): value0 or value1, ...},}

    Arguments
    ---------
    data: tuple
        First entry is a list of time points, second entry is a dict
        mapping names each to a list of values at the corresponding time
        point
    use_left_endpoint: bool
        Optional. Indicates whether each interval should take the value
        of its left endpoint. Default is False, i.e. each interval takes
        the value of its right endpoint.

    Returns
    -------
    dict of dicts
        Maps names to a dict that maps interval-tuples each to the value
        over that interval

    """
    time, value_dict = data
    n_t = len(time)
    if n_t == 1:
        t0 = time[0]
        return {name: {(t0, t0): values[0]} for name, values in value_dict.items()}
    else:
        # This covers the case of n_t > 1 and n_t == 0
        interval_data = {}
        intervals = [(time[i - 1], time[i]) for i in range(1, n_t)]
        for name, values in value_dict.items():
            interval_values = [
                values[i - 1] if use_left_endpoint else values[i] for i in range(1, n_t)
            ]
            interval_data[name] = dict(zip(intervals, interval_values))
        return interval_data
