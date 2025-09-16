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

from typing import Optional, Union, Callable, List
import pandas as pd
from pyomo.environ import (
    ConcreteModel,
    Block,
    Var,
    Param,
    RangeSet,
    Objective,
    Constraint,
    NonNegativeReals,
    Expression,
    maximize,
)
from pyomo.common.config import (
    ConfigDict,
    ConfigValue,
    NonNegativeInt,
    NonNegativeFloat,
    ListOf,
    PositiveInt,
)

from idaes.apps.grid_integration.pricetaker.design_and_operation_models import (
    DesignModelData,
)
from idaes.apps.grid_integration.pricetaker.clustering import (
    generate_daily_data,
    cluster_lmp_data,
    get_optimal_num_clusters,
)
from idaes.apps.grid_integration.pricetaker.unit_commitment import (
    UnitCommitmentData,
    capacity_limits,
    ramping_limits,
    startup_shutdown_constraints,
)
from idaes.core.util.config import ConfigurationError, is_in_range
import idaes.logger as idaeslog

_logger = idaeslog.getLogger(__name__)


# Defining a ConfigDict to simplify the domain validation of
# arguments needed for all the methods of the PriceTakerModel class
CONFIG = ConfigDict()

# List of arguments needed for the `append_lmp_data` method
CONFIG.declare(
    "seed",
    ConfigValue(
        default=20,
        domain=NonNegativeInt,
        doc="Seed value for clustering techniques",
    ),
)
CONFIG.declare(
    "horizon_length",
    ConfigValue(
        domain=PositiveInt,
        doc="Length of each time horizon",
    ),
)
CONFIG.declare(
    "num_clusters_range",
    ConfigValue(
        default=(5, 30),
        domain=ListOf(int, PositiveInt),
        doc="Range of number of clusters for generating elbow plot",
    ),
)
CONFIG.declare(
    "num_clusters",
    ConfigValue(
        domain=PositiveInt,
        doc="Number of clusters/representaive days/periods",
    ),
)
CONFIG.declare(
    "lmp_data",
    ConfigValue(
        domain=ListOf(float),
        doc="List containing locational marginal prices (LMP)",
    ),
)

# List of arguments needed for adding startup/shutdown constraints
CONFIG.declare(
    "minimum_up_time",
    ConfigValue(
        domain=PositiveInt,
        doc="Minimum uptime [in hours]",
    ),
)
CONFIG.declare(
    "minimum_down_time",
    ConfigValue(
        domain=PositiveInt,
        doc="Minimum downtime [in hours]",
    ),
)

# List of arguments for NPV calculation
CONFIG.declare(
    "lifetime",
    ConfigValue(
        domain=PositiveInt,
        doc="Total lifetime of the system [in years]",
    ),
)
CONFIG.declare(
    "discount_rate",
    ConfigValue(
        domain=is_in_range(0, 1),
        doc="Discount rate for annualization [fraction]",
    ),
)
CONFIG.declare(
    "corporate_tax_rate",
    ConfigValue(
        domain=is_in_range(0, 1),
        doc="Effective corporate tax rate [fraction]",
    ),
)
CONFIG.declare(
    "annualization_factor",
    ConfigValue(
        domain=is_in_range(0, 1),
        doc="Capital cost annualization factor [fraction of investment cost/year]",
    ),
)
CONFIG.declare(
    "cash_inflow_scale_factor",
    ConfigValue(
        domain=NonNegativeFloat,
        doc="Scaling factor for net cash inflow calculations",
    ),
)


# pylint: disable = too-many-ancestors, too-many-instance-attributes
class PriceTakerModel(ConcreteModel):
    """Builds a price-taker model for a given system"""

    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)
        self._config = CONFIG()
        self._has_hourly_cashflows = False
        self._has_overall_cashflows = False
        self._op_blk_uc_data = {}  # Save UC data for op. blocks
        self._op_blk_uptime_downtime = {}
        self._linking_constraint_counter = 1

    @property
    def num_representative_days(self):
        """Returns the number of representative days"""
        return self._config.num_clusters

    @num_representative_days.setter
    def num_representative_days(self, value):
        """Setter for the num_representative_days property"""
        if self._config.num_clusters is not None:
            raise ConfigurationError(
                f"num_representative_days is already defined as "
                f"{self.num_representative_days} and it cannot be overwritten."
                f"\n\tInstantiate a new PriceTakerModel object."
            )

        self._config.num_clusters = value

    @property
    def horizon_length(self):
        """Returns the length of each representative day"""
        if self._config.horizon_length is not None:
            return self._config.horizon_length

        _logger.warning(
            "Attribute horizon_length is not specified. Using 24 hours "
            "as the horizon length."
        )
        return 24

    @horizon_length.setter
    def horizon_length(self, value):
        """Setter for the horizon_length property"""
        if self._config.horizon_length is not None:
            raise ConfigurationError(
                f"horizon_length is already defined as {self.horizon_length} and "
                f"it cannot be overwritten.\n\tInstantiate a new PriceTakerModel object."
            )

        self._config.horizon_length = value

    def _assert_lmp_data_exists(self):
        """Raise an error if LMP data does not exist"""
        if self._config.lmp_data is None:
            raise ConfigurationError(
                "LMP data is missing. Please append the LMP data using the "
                "`append_lmp_data` method."
            )

    def append_lmp_data(
        self,
        lmp_data: Union[list, pd.DataFrame, pd.Series],
        num_representative_days: Optional[int] = None,
        horizon_length: Optional[int] = None,
        seed: Optional[int] = 42,
    ):
        """
        Appends the locational marginal price (LMP) data to the
        PriceTakerModel object. If desired, the method can cluster the data
        into representative periods before appending the data.

        Args:
            lmp_data: Union[list, tuple, pd.DataFrame]
                List of locational marginal prices

            num_representative_days: Optional[int], default=None
                number of clusters or representative periods.

            horizon_length: Optional[int], default=None
                Length of each representative period

            seed: Optional[int], default=42
                Seed value for k-means clustering technique
        """

        # Perform domain validation (ConfigDict performs the validation)
        if self._config.lmp_data is None:
            assert len(lmp_data) >= 2  # Do not remove this check!
            self._config.lmp_data = lmp_data
            self._config.num_clusters = num_representative_days
            self._config.horizon_length = horizon_length
            self._config.seed = seed

        else:
            # We do not allow modification of lmp data after it is defined
            raise ConfigurationError(
                "Attempted to overwrite the LMP data. Instantiate a "
                "new PriceTakerModel object to change the LMP data."
            )

    def get_optimal_representative_days(
        self,
        kmin: int = 4,
        kmax: int = 30,
        method: str = "silhouette",
        generate_elbow_plot: bool = True,
    ):
        """
        Returns the optimal number of representative days for the
        given price signal.
        """
        self._assert_lmp_data_exists()
        # Domain validation of kmin and kmax
        self._config.num_clusters_range = (kmin, kmax)

        daily_data = generate_daily_data(self._config.lmp_data, self.horizon_length)

        return get_optimal_num_clusters(
            daily_data, kmin, kmax, method, generate_elbow_plot, self._config.seed
        )

    def _assert_mp_model_exists(self):
        """Raise an error if the multiperiod model does not exist"""
        if not hasattr(self, "period"):
            raise ConfigurationError(
                "Unable to find the multiperiod model. Please use the "
                "build_multiperiod_model method to construct one."
            )

    def build_multiperiod_model(
        self,
        flowsheet_func: Callable,
        flowsheet_options: Optional[dict] = None,
    ):
        """
        Builds the multiperiod model.

        Args:
            flowsheet_func : Callable
                A function that returns an instance of the flowsheet

            flowsheet_options : dict,
                Optional arguments needed for `flowsheet_func`

        """
        self._assert_lmp_data_exists()
        if hasattr(self, "period"):
            # Object may contain a multiperiod model. so raise an error
            raise ConfigurationError(
                "A multiperiod model might already exist, as the object has "
                "`period` attribute."
            )

        lmp_data = self._config.lmp_data

        # pylint: disable=attribute-defined-outside-init
        if self.num_representative_days is not None:
            # Need to use representative days
            rep_days_lmp, rep_days_weights = cluster_lmp_data(
                raw_data=lmp_data,
                horizon_length=self.horizon_length,
                n_clusters=self.num_representative_days,
                seed=self._config.seed,
            )

        else:
            # Use full year price signal
            self.num_representative_days = 1
            self.horizon_length = len(lmp_data)
            rep_days_lmp = {1: {t + 1: val for t, val in enumerate(lmp_data)}}
            rep_days_weights = {1: 1}

        self.set_days = RangeSet(self.num_representative_days)
        self.set_time = RangeSet(self.horizon_length)
        self.rep_days_lmp = rep_days_lmp
        self.rep_days_weights = rep_days_weights

        flowsheet_blk = ConcreteModel()
        if flowsheet_options is None:
            flowsheet_options = {}
        flowsheet_func(flowsheet_blk, **flowsheet_options)

        # Based on some earlier testing, it is faster to clone the model
        # than to build it from scratch. So, instead of using the `rule`
        # argument, we are cloning the model and then transferring the
        # attributes to the `period` block.
        self.period = Block(self.set_days, self.set_time)
        for d, t in self.period:
            self.period[d, t].transfer_attributes_from(flowsheet_blk.clone())

            # If the flowsheet model has LMP defined, update it.
            if hasattr(self.period[d, t], "LMP"):
                self.period[d, t].LMP = self.rep_days_lmp[d][t]

            # Iterate through model to append LMP data if it's been defined
            for blk in self.period[d, t].component_data_objects(Block):
                if hasattr(blk, "LMP"):
                    blk.LMP = self.rep_days_lmp[d][t]

    def _get_operation_blocks(
        self,
        blk_name: str,
        attribute_list: list,
    ):
        """
        Returns a dictionary of operational blocks named 'blk_name'.
        In addition, it also checks the existence of the operational
        blocks, and the existence of specified attributes.
        """
        # Ensure that the multiperiod model exists
        self._assert_mp_model_exists()

        # pylint: disable=not-an-iterable
        op_blocks = {
            d: {t: self.period[d, t].find_component(blk_name) for t in self.set_time}
            for d in self.set_days
        }

        # NOTE: It is sufficient to perform checks only for one block, because
        # the rest of them are clones.
        blk = op_blocks[1][1]  # This object always exists

        # First, check for the existence of the operational block
        if blk is None:
            raise AttributeError(f"Operational block {blk_name} does not exist.")

        # Next, check for the existence of attributes.
        for attribute in attribute_list:
            if not hasattr(blk, attribute):
                raise AttributeError(
                    f"Required attribute {attribute} is not found in "
                    f"the operational block {blk_name}."
                )

        return op_blocks

    def _get_operation_vars(self, var_name: str):
        """
        Returns a dictionary of pointers to the var_name variable located in each flowsheet
        instance. If the variable is not present, then an error is raised.
        """
        # Ensure that the multiperiod model exists
        self._assert_mp_model_exists()

        # pylint: disable=not-an-iterable
        op_vars = {
            d: {t: self.period[d, t].find_component(var_name) for t in self.set_time}
            for d in self.set_days
        }

        # NOTE: It is sufficient to perform checks only for one variable
        if op_vars[1][1] is None:
            raise AttributeError(
                f"Variable {var_name} does not exist in the multiperiod model."
            )

        return op_vars

    def add_linking_constraints(self, previous_time_var: str, current_time_var: str):
        """
        Adds constraints to relate variables across two consecutive time periods. This method is
        usually needed if the system has storage. Using this method, the holdup at the end of the
        previous time period can be equated to the holdup at the beginning of the current time
        period.

        Args:
            previous_time_var : str,
                Name of the operational variable at the end of the previous time step

        current_time_var : str,
                Name of the operational variable at the beginning of the current time step

        """
        old_time_var = self._get_operation_vars(previous_time_var)
        new_time_var = self._get_operation_vars(current_time_var)

        def _rule_linking_constraints(_, d, t):
            if t == 1:
                return Constraint.Skip
            return old_time_var[d][t - 1] == new_time_var[d][t]

        setattr(
            self,
            "variable_linking_constraints_" + str(self._linking_constraint_counter),
            Constraint(self.set_days, self.set_time, rule=_rule_linking_constraints),
        )
        _logger.info(
            f"Linking constraints are added to the model at "
            f"variable_linking_constraints_{self._linking_constraint_counter}"
        )
        self._linking_constraint_counter += 1  # Increase the counter for the new pair

    def _retrieve_uc_data(self, op_blk: str, commodity: str, **kwargs):
        """
        Retrieves unit commitment data for the given operation
        block if it exists. Otherwise, it creates a new object and
        returns the object containing the unit commitment data
        """
        if (op_blk, commodity) in self._op_blk_uc_data:
            self._op_blk_uc_data[op_blk, commodity].update(**kwargs)

        else:
            self._op_blk_uc_data[op_blk, commodity] = UnitCommitmentData(
                blk_name=op_blk, commodity_name=commodity, **kwargs
            )

        return self._op_blk_uc_data[op_blk, commodity]

    def add_capacity_limits(
        self,
        op_block_name: str,
        commodity: str,
        capacity: Union[float, Param, Var, Expression],
        op_range_lb: float,
    ):
        """
        Adds capacity limit constraints of the form:
        op_range_lb * capacity * op_mode(t) <= commodity(t) <= capacity * op_mode(t)
        ex: 0.3 * P_max * op_mode(t) <= P(t) <= P_max * op_mode(t),
        where P(t) is power at time t and op_mode(t) = 1 if the system is operating
        at time t; and op_mode(t) = 0, otherwise.

        Args:
            op_block_name: str,
                Name of the operation model block, ex: ("fs.ngcc")

            commodity: str,
                Name of the commodity on the model the capacity constraints
                will be applied to, ex: ("total_power")

            capacity: float, or Pyomo Var, Param, or Expression
                Maximum capacity on the commodity, ex: (650.0, or m.min_power_capacity)

            op_range_lb: float,
                Ratio of the capacity at minimum stable operation to the maximum capacity.
                Must be a number in the interval [0, 1]
        """
        # _get_operation_blocks method ensures that the operation block exists, and
        # the commodity variable also exists.
        op_blocks = self._get_operation_blocks(
            blk_name=op_block_name, attribute_list=["op_mode", commodity]
        )
        # This method performs all the necessary data checks
        uc_data = self._retrieve_uc_data(
            op_blk=op_block_name,
            commodity=commodity,
            capacity=capacity,
            op_range_lb=op_range_lb,
        )

        # Create a block for storing capacity limit constraints
        cap_limit_blk_name = op_block_name.split(".")[-1] + f"_{commodity}_limits"
        if hasattr(self, cap_limit_blk_name):
            raise ConfigurationError(
                f"Attempting to overwrite capacity limits for {commodity} in {op_block_name}."
            )

        setattr(self, cap_limit_blk_name, Block(self.set_days))
        cap_limit_blk = getattr(self, cap_limit_blk_name)

        # pylint: disable = not-an-iterable
        for d in self.set_days:
            capacity_limits(
                blk=cap_limit_blk[d],
                op_blocks=op_blocks[d],
                uc_data=uc_data,
                set_time=self.set_time,
            )

        # Logger info for where constraint is located on the model
        _logger.info(
            f"Created capacity limit constraints for commodity {commodity} in "
            f"operation block {op_block_name} at {cap_limit_blk.name}"
        )

    def add_ramping_limits(
        self,
        op_block_name: str,
        commodity: str,
        capacity: Union[float, Param, Var, Expression],
        startup_rate: float,
        shutdown_rate: float,
        rampup_rate: float,
        rampdown_rate: float,
    ):
        """
        Adds ramping constraints of the form:
        ramping_var[t] - ramping_var[t-1] <=
        startup_rate * capacity * startup[t] + rampup_rate * capacity * op_mode[t-1];
        ramping_var[t-1] - ramping_var[t] <=
        shutdown_rate * capacity * shutdown[t] + rampdown_rate * capacity * op_mode[t]

        Args:
            op_block_name: str,
                Name of the operation model block, e.g., ("fs.ngcc")

            commodity: str,
                Name of the variable that the ramping constraints will be applied to,
                e.g., "power"

            capacity: float, or Pyomo Var, or Param, or Expression
                String of the name of the entity on the model the ramping constraints
                will be applied to, ex: ("total_power")

            startup_rate: float,
                Fraction of the maximum capacity that variable ramping_var can
                increase during startup (between 0 and 1)

            shutdown_rate: float,
                Fraction of the maximum capacity that variable ramping_var can
                decrease during shutdown (between 0 and 1)

            rampup_rate: float,
                Fraction of the maximum capacity that variable ramping_var can
                increase during operation (between 0 and 1)

            rampdown_rate: float,
                Fraction of the maximum capacity that variable ramping_var can
                decrease during operation (between 0 and 1)
        """
        # Get operational blocks
        op_blocks = self._get_operation_blocks(
            blk_name=op_block_name,
            attribute_list=["op_mode", "startup", "shutdown", commodity],
        )

        # Perform all the necessary datachecks
        uc_data = self._retrieve_uc_data(
            op_blk=op_block_name,
            commodity=commodity,
            capacity=capacity,
            startup_rate=startup_rate,
            shutdown_rate=shutdown_rate,
            rampup_rate=rampup_rate,
            rampdown_rate=rampdown_rate,
        )
        uc_data.assert_ramping_args_present()

        # Creating the pyomo block
        ramp_blk_name = op_block_name.split(".")[-1] + f"_{commodity}_ramping"
        if hasattr(self, ramp_blk_name):
            raise ConfigurationError(
                f"Attempting to overwrite ramping limits for {commodity} in {op_block_name}."
            )

        setattr(self, ramp_blk_name, Block(self.set_days))
        ramp_blk = getattr(self, ramp_blk_name)

        # pylint: disable=not-an-iterable
        for d in self.set_days:
            ramping_limits(
                blk=ramp_blk[d],
                op_blocks=op_blocks[d],
                uc_data=uc_data,
                set_time=self.set_time,
            )

        # Logger info for where constraint is located on the model
        _logger.info(
            f"Created ramping constraints for variable {commodity} "
            f"on operational block {op_block_name} at {ramp_blk.name}"
        )

    def add_startup_shutdown(
        self,
        op_block_name: str,
        des_block_name: Optional[str] = None,
        minimum_up_time: int = 1,
        minimum_down_time: int = 1,
    ):
        """
        Adds minimum uptime/downtime constraints for a given unit/process

        Args:
            op_block_name: str,
                Name of the operation model block, e.g., "fs.ngcc"

            des_block_name: str, default=None,
                Name of the design model block for the operation block
                op_block_name. This argument is specified if the design is
                being optimized simultaneously, e.g., "ngcc_design"

            minimum_up_time: int, default=1,
                Total uptime (must be >= 1), e.g., 4 time periods
                Uptime must include the minimum uptime and the time required
                for shutdown.

            minimum_down_time: int, default=1,
                Total downtime (must be >= 1), e.g., 4 time periods
                Downtime must include the minimum downtime and the time
                required for startup
        """
        # Check minimum_up_time and minimum_down_time for validity
        self.config.minimum_up_time = minimum_up_time
        self.config.minimum_down_time = minimum_down_time

        op_blocks = self._get_operation_blocks(
            blk_name=op_block_name,
            attribute_list=["op_mode", "startup", "shutdown"],
        )
        if des_block_name is not None:
            install_unit = self.find_component(des_block_name + ".install_unit")

            if install_unit is None:
                raise AttributeError(
                    f"Binary variable associated with unit installation is not found "
                    f"in {des_block_name}. \n\tDo not specify des_block_name argument if "
                    f"installation of the unit is not a decision variable."
                )
        else:
            install_unit = 1

        start_shut_blk_name = op_block_name.split(".")[-1] + "_startup_shutdown"
        if hasattr(self, start_shut_blk_name):
            raise ConfigurationError(
                f"Attempting to overwrite startup/shutdown constraints "
                f"for operation block {op_block_name}."
            )

        setattr(self, start_shut_blk_name, Block(self.set_days))
        start_shut_blk = getattr(self, start_shut_blk_name)

        # pylint: disable=not-an-iterable
        for d in self.set_days:
            startup_shutdown_constraints(
                blk=start_shut_blk[d],
                op_blocks=op_blocks[d],
                install_unit=install_unit,
                minimum_up_time=minimum_up_time,
                minimum_down_time=minimum_down_time,
                set_time=self.set_time,
            )

        # Save the uptime and downtime data for reference
        self._op_blk_uptime_downtime[op_block_name] = {
            "minimum_up_time": minimum_up_time,
            "minimum_down_time": minimum_down_time,
        }

        # Logger info for where constraint is located on the model
        _logger.info(
            f"Created startup/shutdown constraints for operation model "
            f" {op_block_name} at {start_shut_blk.name}."
        )

    def add_hourly_cashflows(
        self,
        revenue_streams: Optional[List] = None,
        operational_costs: Optional[List] = None,
    ):
        """
        Adds an expression for the net cash inflow for each flowsheet
        instance. Operational costs in each flowsheet instance may include
        'non_fuel_vom' (non-fuel variable operating costs), 'fuel_cost'
        (cost of fuel), and 'carbon_price' (cost associated with producing
        carbon; i.e., a carbon tax). Revenue streams in each flowsheet instance
        may include the revenue from the wholesale electricity market, revenue
        from alternative products (e.g., hydrogen), etc.

        The net cash inflow is calculated as:

            Sum(revenue streams) - Sum(costs)

        for each time period.

        Args:
            revenue_streams: List of strings representing the names of the
                revenue streams, default: None

                Example: ::

                    ['elec_revenue', 'H2_revenue', ]

            costs: List of strings representing the names of the
                   costs associated with operating at a time period.
                   default: None
                   Example: ::

                        ['hourly_fixed_cost', 'electricity_cost',]
        """
        # Ensure that multiperiod model exists
        self._assert_mp_model_exists()

        if operational_costs is None:
            _logger.warning(
                "Argument operational_costs is not specified, so the total "
                "operational cost will be set to 0."
            )
            operational_costs = []

        if revenue_streams is None:
            _logger.warning(
                "Argument revenue_streams is not specified, so the total "
                "revenue will be set to 0."
            )
            revenue_streams = []

        for p in self.period:
            # Set net profit contribution expressions to 0
            total_cost_expr = 0
            total_revenue_expr = 0

            for blk in self.period[p].component_data_objects(Block):

                # Add costs for each block. If more than one block, may have
                # costs that exist on one block and not on another. (i.e., coproduction)
                for cost in operational_costs:
                    curr_cost = blk.find_component(cost)
                    total_cost_expr += curr_cost if curr_cost is not None else 0

                # Add revenue streams for each block. If more than one block, may have
                # revenues that exist on one block and not on another. (i.e., coproduction)
                for revenue in revenue_streams:
                    curr_rev = blk.find_component(revenue)
                    total_revenue_expr += curr_rev if curr_rev is not None else 0

            # Account for flowsheet-level operational_costs
            for cost in operational_costs:
                curr_cost = self.period[p].find_component(cost)
                total_cost_expr += curr_cost if curr_cost is not None else 0

            # Account for flowsheet-level revenue_streams
            for revenue in revenue_streams:
                curr_rev = self.period[p].find_component(revenue)
                total_revenue_expr += curr_rev if curr_rev is not None else 0

            # Set total cost expression
            self.period[p].total_hourly_cost = Expression(expr=total_cost_expr)

            # Set total revenue expression
            self.period[p].total_hourly_revenue = Expression(expr=total_revenue_expr)

            # Set net cash inflow expression
            self.period[p].net_hourly_cash_inflow = Expression(
                expr=self.period[p].total_hourly_revenue
                - self.period[p].total_hourly_cost
            )

        # Logger info for where constraint is located on the model
        self._has_hourly_cashflows = True
        _logger.info(
            "Created hourly cashflow expressions at:\n\t period[d, t].total_hourly_cost, "
            "period[d, t].total_hourly_revenue, period[d, t].net_hourly_cash_inflow."
        )

    def add_overall_cashflows(
        self,
        lifetime: int = 30,
        discount_rate: float = 0.08,
        corporate_tax_rate: float = 0.2,
        annualization_factor: Optional[float] = None,
        cash_inflow_scale_factor: Optional[float] = 1.0,
    ):
        """
        Builds overall cashflow expressions.

        Args:
            lifetime: int, default=30,
                Lifetime of the unit/process [in years]

            discount_rate: float, default=0.08,
                Discount rate [fraction] for calculating the current value
                of the cashflow. It is also used to compute the annualization
                factor. Must be between 0 and 1.

            corporate_tax_rate: float, default=0.2,
                Fractional value of corporate tax used in NPV calculations.
                Must be between 0 and 1.

            annualization_factor: float, default=None,
                Annualization factor

            cash_inflow_scale_factor: float, default=1.0,
                Scaling factor for the sum of net hourly cashflows.
                If the specified price signal is for the full year, then the
                value of this argument must be 1.0. However, if the specified
                price signal is for a short period, say 1 month, then an
                appropriate scaling factor can be specified to compute the value
                for a year (12, for the case of 1 month's price signal).

        """
        # Ensure that multiperiod model exists
        self._assert_mp_model_exists()

        # Domain validation of input arguments
        self.config.lifetime = lifetime
        self.config.discount_rate = discount_rate
        self.config.corporate_tax = corporate_tax_rate
        self.config.annualization_factor = annualization_factor
        self.config.cash_inflow_scale_factor = cash_inflow_scale_factor

        if not self._has_hourly_cashflows:
            raise ConfigurationError(
                "Hourly cashflows are not added to the model. Please run "
                "add_hourly_cashflows method before calling the "
                "add_overall_cashflows method."
            )

        capex_expr = 0
        fom_expr = 0

        count_des_blks = 0
        for blk in self.component_data_objects(Block):
            if isinstance(blk, DesignModelData):
                count_des_blks += 1
                capex_expr += blk.capex
                fom_expr += blk.fom

        if count_des_blks == 0:
            _logger.warning(
                "No design blocks were found, so the overall capital cost "
                "(capex) and the fixed O&M cost (fom) are set to 0."
            )

        # pylint: disable=attribute-defined-outside-init
        self.cashflows = Block()
        cf = self.cashflows
        cf.capex = Var(within=NonNegativeReals, doc="Total CAPEX")
        cf.capex_calculation = Constraint(expr=cf.capex == capex_expr)

        cf.fom = Var(within=NonNegativeReals, doc="Yearly Fixed O&M Cost")
        cf.fom_calculation = Constraint(expr=cf.fom == fom_expr)

        cf.depreciation = Var(within=NonNegativeReals, doc="Yearly depreciation")
        cf.depreciation_calculation = Constraint(
            expr=cf.depreciation == cf.capex / lifetime
        )

        cf.net_cash_inflow = Var(doc="Net cash inflow")
        cf.net_cash_inflow_calculation = Constraint(
            expr=cf.net_cash_inflow
            == cash_inflow_scale_factor
            * sum(
                self.rep_days_weights[d] * self.period[d, t].net_hourly_cash_inflow
                for d, t in self.period
            )
        )

        cf.corporate_tax = Var(within=NonNegativeReals, doc="Corporate tax")
        cf.corporate_tax_calculation = Constraint(
            expr=cf.corporate_tax
            >= corporate_tax_rate * (cf.net_cash_inflow - cf.fom - cf.depreciation)
        )

        cf.net_profit = Var(doc="Net profit after taxes")
        cf.net_profit_calculation = Constraint(
            expr=cf.net_profit == cf.net_cash_inflow - cf.fom - cf.corporate_tax
        )

        if annualization_factor is None:
            # If the annualization factor is not specified
            annualization_factor = discount_rate / (
                1 - (1 + discount_rate) ** (-lifetime)
            )

        cf.lifetime_npv = Expression(
            expr=(1 / annualization_factor) * cf.net_profit - cf.capex
        )
        cf.npv = Expression(
            expr=cf.net_profit - annualization_factor * cf.capex,
        )

        self._has_overall_cashflows = True
        _logger.info(f"Overall cashflows are added to the block {cf.name}")

    def add_objective_function(self, objective_type="npv"):
        """
        Appends the objective function to the model.

        Args:
            objective_type: str, default="npv",
                Supported objective functions are
                annualized net present value ("npv"),
                net profit ("net_profit"), and
                lifetime net present value ("lifetime_npv").
        """
        # pylint: disable = attribute-defined-outside-init
        if not self._has_overall_cashflows:
            raise ConfigurationError(
                "Overall cashflows are not appended. Please run the "
                "add_overall_cashflows method."
            )

        try:
            self.obj = Objective(
                expr=getattr(self.cashflows, objective_type),
                sense=maximize,
            )

        except AttributeError as msg:
            raise ConfigurationError(
                f"{objective_type} is not a supported objective function."
                f"Please specify either npv, or lifetime_npv, or net_profit "
                f"as the objective_type."
            ) from msg
