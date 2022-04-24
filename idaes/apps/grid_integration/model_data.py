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

from numbers import Real
from abc import ABC, abstractmethod

{
    "bus": "Adams",
    "in_service": True,
    "mbase": 100.0,
    "pg": 76.0,
    "qg": -2.31,
    "p_min": 30.0,
    "p_max": 76.0,
    "q_min": -25.0,
    "q_max": 30.0,
    "ramp_q": 2.0,
    "fuel": "Coal",
    "unit_type": "STEAM",
    "area": "1",
    "zone": "12.0",
    "generator_type": "thermal",
    "p_fuel": {
        "data_type": "fuel_curve",
        "values": [(30.0, 347.73), (45.3, 481.36), (60.7, 633.22), (76.0, 796.18)],
    },
    "startup_fuel": [(4.0, 1111.637669), (10.0, 1599.13457), (12.0, 1738.41)],
    "non_fuel_startup_cost": 0.0,
    "shutdown_cost": 0.0,
    "agc_capable": True,
    "p_min_agc": 30.0,
    "p_max_agc": 76.0,
    "ramp_agc": 2.0,
    "ramp_up_60min": 120.0,
    "ramp_down_60min": 120.0,
    "fuel_cost": 2.11399,
    "startup_capacity": 30.0,
    "shutdown_capacity": 30.0,
    "min_up_time": 8.0,
    "min_down_time": 4.0,
    "initial_status": 9.0,
    "initial_p_output": 30.0,
    "initial_q_output": 0.0,
}


class BaseValidator(ABC):
    def __set_name__(self, cls, prop_name):
        self.prop_name = prop_name

    def __set__(self, instance, value):

        self._validate(instance, value)
        instance.__dict__[self.prop_name] = value

    def __get__(self, instance, cls):

        if instance is None:
            return self
        return instance.__dict__.get(self.prop_name, None)

    @abstractmethod
    def _validate(self, instance, value):
        pass


class RealValueValidator(BaseValidator):
    def __init__(self, min_val=None, max_val=None):
        self.min_val = min_val
        self.max_val = max_val

    def _validate(self, instance, value):

        if not isinstance(value, Real):
            raise TypeError(f"Value for {self.prop_name} shoulde be real numbers.")

        if self.min_val is not None and value < self.min_val:
            raise ValueError(
                f"Value should be greater than or equal to {self.min_val}."
            )

        if self.max_val is not None and value > self.max_val:
            raise ValueError(f"Value should be less than or equal to {self.max_val}.")


class AtLeastPminValidator(BaseValidator):
    def _validate(self, instance, value):

        pmin = getattr(instance, "p_min", None)

        if pmin is None:
            raise RuntimeError(
                f"Value for Pmin does not exist, so {self.prop_name} cannot be validated."
            )

        if not isinstance(value, Real):
            raise TypeError(f"Value for {self.prop_name} shoulde be real numbers.")

        if value < pmin:
            raise ValueError(
                f"Value for {self.prop_name} shoulde be greater or equal to Pmin."
            )


class GeneratorModelData:

    p_min = RealValueValidator(min_val=0)
    min_down_time = RealValueValidator(min_val=0)
    min_up_time = RealValueValidator(min_val=0)
    ramp_up_60min = RealValueValidator(min_val=0)
    ramp_down_60min = RealValueValidator(min_val=0)

    p_max = AtLeastPminValidator()
    shutdown_capacity = AtLeastPminValidator()
    startup_capacity = AtLeastPminValidator()

    def __init__(
        self,
        gen_name,
        gen_type,
        p_min,
        p_max,
        min_down_time,
        min_up_time,
        ramp_up_60min,
        ramp_down_60min,
        shutdown_capacity,
        startup_capacity,
        output_points,
        marginal_costs,
        start_up_lags=None,
        start_up_costs=None,
        fixed_commitment=None,
    ):
        self.gen_name = gen_name
        self.gen_type = gen_type
        self.p_min = p_min
        self.p_max = p_max
        self.min_down_time = min_down_time
        self.min_up_time = min_up_time
        self.ramp_up_60min = ramp_up_60min
        self.ramp_down_60min = ramp_down_60min
        self.shutdown_capacity = shutdown_capacity
        self.startup_capacity = startup_capacity
        self.fixed_commitment = fixed_commitment

        self._assemble_default_cost_bids(output_points, marginal_costs)
        self._assemble_default_start_up_cost_bids(start_up_lags, start_up_costs)

    @property
    def gen_name(self):
        return self._gen_name

    @gen_name.setter
    def gen_name(self, value):
        if not isinstance(value, str):
            raise TypeError(
                f"Value for generator names must be str, but {type(value)} is provided."
            )
        self._gen_name = value

    @property
    def gen_type(self):
        return self._gen_type

    @gen_type.setter
    def gen_type(self, value):
        allowed_types = ["thermal", "renewable"]
        if value not in allowed_types:
            raise ValueError(
                f"Value for generator types must be one of {allowed_types}, but {value} is provided."
            )
        self._gen_type = value

    def _assemble_default_cost_bids(self, output_points, marginal_costs):
        pass

    def _assemble_default_start_up_cost_bids(self, start_up_lags, start_up_costs):
        pass
