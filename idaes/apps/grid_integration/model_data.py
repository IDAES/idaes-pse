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
        generator_type,
        p_min,
        p_max,
        min_down_time,
        min_up_time,
        ramp_up_60min,
        ramp_down_60min,
        shutdown_capacity,
        startup_capacity,
        production_cost_bid_pairs=None,
        startup_cost_pairs=None,
        fixed_commitment=None,
    ):
        self.gen_name = gen_name
        self.generator_type = generator_type
        self.p_min = p_min
        self.p_max = p_max
        self.min_down_time = min_down_time
        self.min_up_time = min_up_time
        self.ramp_up_60min = ramp_up_60min
        self.ramp_down_60min = ramp_down_60min
        self.shutdown_capacity = shutdown_capacity
        self.startup_capacity = startup_capacity

        fixed_commitment_allowed_values = [0, 1, None]
        if fixed_commitment not in fixed_commitment_allowed_values:
            raise ValueError(
                f"Value for generator fixed commitment must be one of {fixed_commitment_allowed_values}, but {fixed_commitment} is provided."
            )
        self.fixed_commitment = fixed_commitment

        self.default_bids = self._assemble_default_cost_bids(production_cost_bid_pairs)
        self.default_startup_bids = self._assemble_default_startup_cost_bids(
            startup_cost_pairs
        )

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
    def generator_type(self):
        return self._generator_type

    @generator_type.setter
    def generator_type(self, value):
        allowed_types = ["thermal", "renewable"]
        if value not in allowed_types:
            raise ValueError(
                f"Value for generator types must be one of {allowed_types}, but {value} is provided."
            )
        self._generator_type = value

    def _check_empty_and_sort_cost_pairs(self, pair_description, pairs):

        if pairs is None or len(pairs) == 0:
            raise ValueError(f"Empty {pair_description} are provided.")

        # sort based on power output
        pairs.sort(key=lambda p: p[0])

        return

    def _assemble_default_cost_bids(self, production_cost_bid_pairs):

        self._check_empty_and_sort_cost_pairs(
            pair_description="production cost pairs", pairs=production_cost_bid_pairs
        )

        if production_cost_bid_pairs[0][0] != self.p_min:
            raise ValueError(
                f"The first power output in the bid should be the Pmin {self.p_min}, but {production_cost_bid_pairs[0][0]} is provided."
            )

        if production_cost_bid_pairs[-1][0] != self.p_max:
            raise ValueError(
                f"The last power output in the bid should be the Pmax {self.p_max}, but {production_cost_bid_pairs[-1][0]} is provided."
            )

        return {
            power: marginal_cost for power, marginal_cost in production_cost_bid_pairs
        }

    def _assemble_default_startup_cost_bids(self, startup_cost_pairs):

        self._check_empty_and_sort_cost_pairs(
            pair_description="startup cost pairs", pairs=startup_cost_pairs
        )

        if startup_cost_pairs[0][0] != self.min_down_time:
            raise ValueError(
                f"The first startup lag should be the same as minimum down time {self.min_down_time}, but {startup_cost_pairs[0][0]} is provided."
            )

        return {
            startup_lag: startup_cost
            for startup_lag, startup_cost in startup_cost_pairs
        }
