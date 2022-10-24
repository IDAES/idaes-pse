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
from math import isclose
from abc import ABC, abstractmethod


class BaseValidator(ABC):

    """
    A base data descriptor that validate values before store them.
    """

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

        """
        Validate the provided value.

        Args:
            instance: an instance that uses this data descriptor.
            value: value that needs to be validated.
        """
        pass


class StrValidator(BaseValidator):
    def _validate(self, instance, value):

        if not isinstance(value, str):
            raise TypeError(
                f"Value for {self.prop_name} must be str, but {type(value)} is provided."
            )


class RealValueValidator(BaseValidator):

    """
    A data descriptor that validate whether a value is a real number and whether
    it is within its lower and upper bounds (if provided).
    """

    def __init__(self, min_val=None, max_val=None):
        self.min_val = min_val
        self.max_val = max_val

    def _validate(self, instance, value):

        """
        Validate whether the provided value is within bounds and is a real number.

        Args:
            instance: an instance that uses this data descriptor.
            value: value that needs to be validated.
        """

        if not isinstance(value, Real):
            raise TypeError(f"Value for {self.prop_name} shoulde be real numbers.")

        if (
            self.min_val is not None
            and value < self.min_val
            and not isclose(value, self.min_val)
        ):
            raise ValueError(
                f"Value should be greater than or equal to {self.min_val}."
            )

        if (
            self.max_val is not None
            and value > self.max_val
            and not isclose(value, self.max_val)
        ):
            raise ValueError(f"Value should be less than or equal to {self.max_val}.")


class AtLeastPminValidator(BaseValidator):

    """
    A data descriptor that validate whether a value is greater or equal to the
    generator's Pmin. These values include Pmax, shut-down capacity, start-up
    capacity, and etc.
    """

    def _validate(self, instance, value):

        """
        Validate whether the provided value is at least Pmin.

        Args:
            instance: an instance that uses this data descriptor.
            value: value that needs to be validated.
        """

        pmin = getattr(instance, "p_min", None)

        if pmin is None:
            raise RuntimeError(
                f"Value for Pmin does not exist, so {self.prop_name} cannot be validated."
            )

        if not isinstance(value, Real):
            raise TypeError(f"Value for {self.prop_name} shoulde be real numbers.")

        if value < pmin and not isclose(value, pmin):
            raise ValueError(
                f"Value for {self.prop_name} shoulde be greater or equal to Pmin."
            )


class GeneratorModelData:

    """
    A class that holds data for generator parameters.
    """

    gen_name = StrValidator()
    bus = StrValidator()

    p_min = RealValueValidator(min_val=0)
    p_max = AtLeastPminValidator()

    def __init__(self, gen_name, bus, p_min, p_max, fixed_commitment=None):

        self.gen_name = gen_name
        self.bus = bus
        self.p_min = p_min
        self.p_max = p_max

        fixed_commitment_allowed_values = [0, 1, None]
        if fixed_commitment not in fixed_commitment_allowed_values:
            raise ValueError(
                f"Value for generator fixed commitment must be one of {fixed_commitment_allowed_values}, but {fixed_commitment} is provided."
            )
        self.fixed_commitment = fixed_commitment

    def __iter__(self):
        """
        Make it iteratble.
        """

        if not hasattr(self, "_collection"):
            self._collection = [
                name
                for name in dir(self)
                if not name.startswith("__")
                and not name.startswith("_")
                and not callable(getattr(self, name))
            ]
            self._index = -1

        return self

    def __next__(self):
        """
        Implement iterator method.
        """

        self._index += 1
        if self._index >= len(self._collection):
            self._index = -1
            raise StopIteration
        else:
            name = self._collection[self._index]
            return name, getattr(self, name)


class ThermalGeneratorModelData(GeneratorModelData):

    """
    Adds run time, ramping, shutdown and startup information for thermal generators.
    Requires as initial status and power production.
    Bids are in ( MWh, $ ) pairs
    """

    min_down_time = RealValueValidator(min_val=0)
    min_up_time = RealValueValidator(min_val=0)
    ramp_up_60min = RealValueValidator(min_val=0)
    ramp_down_60min = RealValueValidator(min_val=0)

    shutdown_capacity = AtLeastPminValidator()
    startup_capacity = AtLeastPminValidator()

    def __init__(
        self,
        gen_name,
        bus,
        p_min,
        p_max,
        min_down_time,
        min_up_time,
        ramp_up_60min,
        ramp_down_60min,
        shutdown_capacity,
        startup_capacity,
        initial_status,
        initial_p_output,
        fixed_commitment=None,
        production_cost_bid_pairs=None,
        startup_cost_pairs=None,
    ):

        super().__init__(
            gen_name=gen_name,
            bus=bus,
            p_min=p_min,
            p_max=p_max,
            fixed_commitment=fixed_commitment,
        )

        self.min_down_time = min_down_time
        self.min_up_time = min_up_time
        self.ramp_up_60min = ramp_up_60min
        self.ramp_down_60min = ramp_down_60min
        self.shutdown_capacity = shutdown_capacity
        self.startup_capacity = startup_capacity
        self.initial_status = initial_status
        self.initial_p_output = initial_p_output

        self.p_cost = self._assemble_default_cost_bids(production_cost_bid_pairs)
        self.startup_cost = self._assemble_default_startup_cost_bids(startup_cost_pairs)

    def _check_empty_and_sort_cost_pairs(self, pair_description, pairs):

        """
        Check whether a list of pairs for production and startup costs is empty,
        and then sort them based on power output and startup lag.

        Args:
            pair_description: a string description for the list of pairs
            pairs: a list of pairs

        Returns:
            None

        """

        if pairs is None or len(pairs) == 0:
            raise ValueError(f"Empty {pair_description} are provided.")

        # sort based on power output
        pairs.sort(key=lambda p: p[0])

        return

    def _assemble_default_cost_bids(self, production_cost_bid_pairs):

        """
        Assemble the default production cost bids.

        Args:
            production_cost_bid_pairs: a list of pairs, which consit of power output
            and corresponding marginal costs

        Returns:
            dict: the assembled production cost bids
        """

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

        return production_cost_bid_pairs

    def _assemble_default_startup_cost_bids(self, startup_cost_pairs):

        """
        Assemble the default startup cost bids.

        Args:
            startup_cost_pairs: a list of pairs, which consit of startup time lag
            and corresponding startup costs

        Returns:
            dict: the assembled startup cost bids
        """

        self._check_empty_and_sort_cost_pairs(
            pair_description="startup cost pairs", pairs=startup_cost_pairs
        )

        if startup_cost_pairs[0][0] != self.min_down_time:
            raise ValueError(
                f"The first startup lag should be the same as minimum down time {self.min_down_time}, but {startup_cost_pairs[0][0]} is provided."
            )

        return startup_cost_pairs

    @property
    def generator_type(self):
        """
        Generator type property: fixed to thermal.
        """
        return "thermal"

    @property
    def initial_status(self):

        """
        Generator initial status proptery. If positive, the number of hours prior
        to (and including) t=0 that the unit has been on. If negative, the number
        of hours prior to (and including) t=0 that the unit has been off. The value
        cannot be 0, by definition.
        """

        return self._initial_status

    @initial_status.setter
    def initial_status(self, value):

        """
        Generator initial status proptery setter. Validate the value before setting.
        """

        if not isinstance(value, Real):
            raise TypeError("Value for initial_status shoulde be real numbers.")

        if isclose(value, 0):
            raise ValueError("Value for initial_status cannot be zero.")

        self._initial_status = value

    @property
    def initial_p_output(self):
        """
        Generator initial power output proptery.
        """
        return self._initial_p_output

    @initial_p_output.setter
    def initial_p_output(self, value):

        """
        Generator initial power output proptery setter.
        """

        if not isinstance(value, Real):
            raise TypeError("Value for initial_p_output shoulde be real numbers.")

        if (
            self.initial_status > 0
            and value < self.p_min
            and not isclose(value, self.p_min)
        ):
            raise ValueError(
                f"The initial status of the generator was on before T0, so the initial power output should at least be p_min {self.p_min}, but {value} is provided."
            )

        if self.initial_status < 0 and (value > 0 and not isclose(value, 0)):
            raise ValueError(
                f"The initial status of the generator was off before T0, so the initial power output should at 0, but {value} is provided."
            )

        self._initial_p_output = value


class RenewableGeneratorModelData(GeneratorModelData):

    """
    Adds a single production cost for renewable generators.
    """

    p_cost = RealValueValidator(min_val=0)

    def __init__(self, gen_name, bus, p_min, p_max, p_cost, fixed_commitment=None):

        super().__init__(
            gen_name=gen_name,
            bus=bus,
            p_min=p_min,
            p_max=p_max,
            fixed_commitment=fixed_commitment,
        )
        self.p_cost = p_cost

    @property
    def generator_type(self):
        """
        Generator type property: fixed to renewable.
        """
        return "renewable"
