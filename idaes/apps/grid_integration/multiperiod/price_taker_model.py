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
import pandas as pd
import numpy as np
from functools import reduce
from operator import attrgetter
import os.path

from importlib import resources
from pathlib import Path
from pyomo.environ import (
    ConcreteModel,
    Block,
    Var,
    RangeSet,
    Objective,
    Constraint,
    NonNegativeReals,
    Expression,
    maximize,
    Param,
)

from idaes.apps.grid_integration import MultiPeriodModel
from idaes.apps.grid_integration.multiperiod.design_and_operation_models import (
    DesignModelData,
    OperationModelData,
)

from sklearn.cluster import KMeans
from kneed import KneeLocator

import matplotlib.pyplot as plt

import logging

_logger = logging.getLogger(__name__)


def deepgetattr(obj, attr):
    return reduce(getattr, attr.split('.'), obj)


class PriceTakerModel(ConcreteModel):
    def __init__(self, seed=20, horizon_length=24):
        super().__init__()
        self.seed = seed
        self.horizon_length = horizon_length

    @property
    def seed(self):
        return self._seed

    @seed.setter
    def seed(self, value):
        if not isinstance(value, int):
            raise ValueError(f"seed must be an integer, but {value} is not an integer")
        self._seed = value

    @property
    def horizon_length(self):
        return self._horizon_length

    @horizon_length.setter
    def horizon_length(self, value):
        if value <= 0:
            raise ValueError(f"horizon_length must be > 0, but {value} is provided.")
        if not isinstance(value, int):
            raise ValueError(f"horizon_length must be an integer, but {value} is not an integer")
        self._horizon_length = value

    def generate_daily_data(self,raw_data,day_list):
        """
        Function used to generate the daily data in a usable format
        from the raw data provided.

        Args:
            raw_data:   LMP data in raw format
            day_list:   List of days to be used from the raw data
        
        Returns:
            daily_data: Correctly formatted LMP data for later use
        """

        daily_data = pd.DataFrame(columns=day_list)

        # Extracting data to populate empty dataframe
        i = 0
        j = self._horizon_length
        day = 1
        while j <= len(raw_data):
            daily_data[day] = raw_data[i:j].reset_index(drop=True)
            i = j
            j = j + self._horizon_length
            day = day + 1
        
        return daily_data

    def reconfigure_raw_data(self, raw_data):
        """
        Reconfigures the raw price series data into a usable form

        Args:
            raw_data: imported price series data

        Returns:
            daily_data: reconfigured price series data
            scenarios:  ??? Add description ???
        """

        # Get column headings
        column_head = raw_data.columns.tolist()
        # Remove the date/time column
        scenarios = column_head[2:]

        # Creating an empty dataframe to store daily data for clustering
        day_list = list(range(1, (len(raw_data) // self._horizon_length) + 1))

        # Generate daily data
        for i in scenarios[0:1]:
            daily_data = self.generate_daily_data(
                raw_data=raw_data[i], day_list=day_list
            )

        return daily_data, scenarios

    def get_optimal_n_clusters(
        self, 
        daily_data, 
        kmin=None, 
        kmax=None, 
        plot=False,
    ):
        """
        Determines the appropriate number of clusters needed for a
        given price signal.

        Args:
            daily_data: reconfigured price series data
            kmin:       minimum number of clusters
            kmax:       maximum number of clusters
            plot:       flag to determine if an elbow plot should be displayed

        Returns:
            n_clusters:     the optimal number of clusters for the given data
            inertia_values: within-cluster sum-of-squares
        """
        # Check kmin and kmax for validity
        if kmin is None:
            kmin = 4
        if kmax is None:
            kmax = 30
            _logger.warning(f"kmax was not set - using a default value of 30.")
        
        if not isinstance(kmin, int):
            raise ValueError(f"kmin must be an integer, but {kmin} is not an integer")
        if not isinstance(kmax, int):
            raise ValueError(f"kmax must be an integer, but {kmax} is not an integer")
        if kmin < 1:
            raise ValueError(f"kmin must be > 0, but {kmin} is provided.")
        if kmax < 1:
            raise ValueError(f"kmax must be > 0, but {kmax} is provided.")
        if  kmin >= kmax:
            raise ValueError(f"kmin must be less than kmax, but {kmin} >= {kmax}")
        
        k_values = range(kmin, kmax)
        inertia_values = []

        np.random.seed(self._seed)

        # Compute the inertia (SSE) for k clusters

        for k in k_values:
            kmeans = KMeans(n_clusters=k).fit(daily_data.transpose())
            inertia_values.append(kmeans.inertia_)

        # Identify the "elbow point"
        elbow_point = KneeLocator(
            k_values, inertia_values, curve="convex", direction="decreasing"
        )
        n_clusters = elbow_point.knee

        print(f"Optimal # of clusters is: {n_clusters}")

        if n_clusters + 2 >= kmax:
            _logger.warning(
                f"Optimal number of clusters is close to kmax: {kmax}. Consider increasing kmax."
            )

        if plot == True:
            plt.show()
            plt.plot(k_values, inertia_values)
            plt.axvline(x=n_clusters, color="red", linestyle="--", label="Elbow")
            plt.xlabel("Number of clusters")
            plt.ylabel("Inertia")
            plt.title("Elbow Method")
            plt.xlim(kmin, kmax)
            plt.grid()
            plt.show()

        return n_clusters, inertia_values

    def cluster_lmp_data(self, raw_data, n_clusters):
        """
        Clusters the given price signal in n_clusters. This method supports k-means, k-meteiod,...
        techniques for clustering.

        Args:
            raw_data:   imported price series data
            n_clusters: number of clusters desired for the data (representative days)

        Returns:
            lmp_data:   dict of representative day LMP data, indices are indexed
                        by integers starting at 1 (ex: {1: {1: 4, 2: 3, 3: 5}, 
                                                        2: {1: 1, 2: 7, 3: 3}})
                                               Format: {day: {time_period: LMP}}
            weights:    dict of weights for each representative day, indexed the
                        same way as lmp_data      (ex: {1: 45, 2: 56})
                                               Format: {day: weight}
        """
        # reconfiguring raw data 
        daily_data, scenarios  = self.reconfigure_raw_data(raw_data)
        
        # KMeans clustering with the optimal number of clusters
        kmeans = KMeans(n_clusters=n_clusters).fit(daily_data.transpose())
        centroids = kmeans.cluster_centers_
        labels = kmeans.labels_

        # Set any centroid values that are < 1e-4 to 0 to avoid noise
        for d in range(n_clusters):
            for t in range(24):
                if centroids[d][t] < 1e-4:
                    centroids[d][t] = 0

        n_clusters_list = range(0, n_clusters)
        weights_counter = np.zeros(n_clusters)

        # Compute weight for each cluster by counting its occurrences in the dataset
        for j in range(0, len(labels)):
            for k in n_clusters_list:
                if labels[j] == k:
                    weights_counter[k] += 1

        # Create dicts for lmp data and the weight of each cluster
        rep_days_data = {}
        weights_data = {}

        
        rep_days_data = pd.DataFrame(centroids.transpose(), columns = range(1,n_clusters+1))
        lmp_data = rep_days_data.to_dict()
        weights_data = pd.DataFrame(weights_counter)
        weights_data.index = np.arange(1, len(weights_data) + 1)
        weights = weights_data.to_dict()

        return lmp_data, weights

    def append_lmp_data(
        self,
        file_path,
        sheet=None,
        column_name=None,
        n_clusters=None,
        horizon_length=None,
    ):
        """
        This function appends LMP data to the PriceTakerModel using single
        or multiple year signals, as well as using full price data or 
        clustering to use representative days.

        Args:
            file_path:      path to the file containing LMP data
            sheet:          if the file is in axcel file, the sheet name
                            where the LMP data is located - (default: None)
            column_name:    list of integer year names if multiple years of
                            data are to be used, None otherwise - (default: None)
            n_clusters:     number of clusters for the data (representative days)
                            None if representative days not used - (default: None)
            horizon_length: if a value is given, this will be used to set the
                            horizon_length attribute of the PriceTakerModel.
                            (default: None --> use existing horizon_length (default: 24))

        Returns:
            
        """
        if os.path.exists(file_path):
            path_to_file = file_path
        else:
            raise ValueError()

        if '.xls' in path_to_file[-5:]:
            full_data = pd.read_excel( path_to_file, sheet_name=[sheet])[sheet]
        elif '.csv' in path_to_file[-5:]:
            full_data = pd.read_csv(path_to_file, )
        
        if horizon_length is not None:
            self.horizon_length(horizon_length)
        
        # editing the data
        if isinstance(column_name, list) and n_clusters is not None:
            # Multiple years and representative days
            self.set_years = [int(y) for y in column_name]
            self.set_days = RangeSet(1,n_clusters)
            self._n_time_points = self.horizon_length if self.horizon_length is not None else 24
            self.set_time = RangeSet(self._n_time_points)

            self.LMP = {}
            self.WEIGHTS = {}

            for year in column_name:
                price_data = full_data[year]
                lmp_data, weights = self.cluster_lmp_data(
                    price_data, n_clusters
                )
                y = int(year)

                for d in self.set_days:
                    for t in self.set_time:
                        self.LMP[t, d, y] = lmp_data[d][t-1]
                        self.WEIGHTS[d, y] = weights[0][d]

            return

        elif isinstance(column_name, list):
            # Multiple years, use fullyear price signal for each year
            self.set_years = [int(y) for y in column_name]
            self.set_days = None
            self._n_time_points = len(full_data)
            self.set_time = RangeSet(self._n_time_points)

            self.LMP = {}

            for year in column_name:
                price_data = full_data[year]
                y = int(year)

                for t in self.set_time:
                    self.LMP[t, y] = price_data[t - 1]

            return

        elif n_clusters is not None:
            # Single price signal, use reprentative days
            self.set_years = None
            self.set_days = RangeSet(1,n_clusters)
            self._n_time_points = self.horizon_length if self.horizon_length is not None else 24
            self.set_time = RangeSet(self._n_time_points)

            self.LMP = {}
            self.WEIGHTS = {}

            price_data = full_data
            lmp_data, weights = self.cluster_lmp_data(
                price_data, n_clusters
            )
            
            for d in self.set_days:
                for t in self.set_time:
                    self.LMP[t, d] = lmp_data[d][t-1]
                    self.WEIGHTS[d] = weights[0][d]
                    
            return

        else:
            # Single price signal, use full year's price signal
            self.set_years = None
            self.set_days = None
            self._n_time_points = len(full_data)
            self.set_time = RangeSet(self._n_time_points)

            price_data = full_data[column_name].to_list()
            self.LMP = {t: price_data[t - 1] for t in self.set_time}

            return

    def build_multiperiod_model(self, **kwargs):

        # if not self.model_sets_available:
        #     raise Exception(
        #         "Model sets have not been defined. Run get_lmp_data to construct model sets"
        #     )

        self.mp_model = MultiPeriodModel(
            n_time_points=self._n_time_points,
            set_days=self.set_days,
            set_years=self.set_years,
            use_stochastic_build=True,
            **kwargs,
        )

        # If append_lmp_data is automatic, need to append the LMP data.
        # Check if LMP has already been defined with the append_lmp_data
        # function above
        LMP_exists = hasattr(self, 'LMP')

        # Iterate through model to append LMP data if it's been defined
        # and the model says it should be (default)
        period = self.mp_model.period
        for p in period:
            for blk in period[p].component_data_objects(Block):
                if isinstance(blk, OperationModelData):
                    if blk.config.append_lmp_data:
                        if not LMP_exists:
                            raise ValueError(f"OperationModelData has been defined to automatically " +
                                             f"populate LMP data. However, m.LMP does not exist. " +
                                             f"Please run append_lmp_data() first or set the " +
                                             f"append_lmp_data attribute to False when configuring " +
                                             f"your OperationModelData object.")
                        blk.LMP = self.LMP[p]




    #TODO: (1) Need to determine whether the minimum opterating power can be a parameter 
    #      or whether it largly depends on the capacity of the plant. 
    #      (2) Need to determine if the start up and shutdown limits can be modeled as bulk amounts vs
    #      a percentage of capacity 
    #      (3) How to model/code the cosntraints in without knowing what variable
    #      these constraints will be applied to. 
    def add_ramping_constraints(
        self,
        op_blk,
        design_blk,
        capacity_var,
        ramping_var,
        constraint_type,
        linearization=True,
        op_range_lb=0.6, 
        startup_rate=0.7, 
        shutdown_rate=0.7,
        ramp_up_rate=0.7, 
        ramp_down_rate=0.7,
    ):
        """
        Adds ramping constraints of the form:
            -ramp_down_limit <= var(t) - var(t-1) <= ramp_up_limit on var
        

        Args: 
            op_blk:             String of the name of the operation model block, ex: ("fs.op_name")
            design_blk:         String of the name of the design model block, ex: ("m.design_name")
            capacity_var:       String of the name of the entity on the model the ramping constraints 
                                will be applied to, ex: ("total_power") OR integer/float value of the
                                capacity value
            ramping_var:        String of the name of the variable that the ramping constraints will
                                be applied to
            constraint_type:    String to choose between linear and nonlinear constraints. Valid 
                                inputs are in: ["linear", "nonlinear"]
            linearization:      Boolean indicating whether linearization is used when constraint_type is
                                "nonlinear", True to use linearization, False otherwise
            op_range_lb:        The fraction of the maximum capacity that represents the lower operating 
                                bound (between 0 and 1)
            startup_rate:       The fraction of the maximum capacity that variable ramping_var can 
                                increase during startup (between 0 and 1)
            shutdown_rate:      The fraction of the maximum capacity that variable ramping_var can 
                                decrease during shutdown (between 0 and 1)
            ramp_up_rate:       The fraction of the maximum capacity that variable ramping_var can 
                                increase during operation (between 0 and 1)
            ramp_down_rate:     The fraction of the maximum capacity that variable ramping_var can 
                                decrease during operation (between 0 and 1)

    
        Assumptions/relationship:
            total_power_upper_bound >= ramp_up_limit >= startup_limit >= total_power_lower_bound > 0
            total_power_upper_bound  >= ramp_down_limit >= shutdown_limit >= total_power_lower_bound > 0
        """
        # Checking that all ramping rates are between 0 and 1 and that the
        # lower bound for operation is less than the startup/shutdown ramps
        if startup_rate <= 0  or startup_rate > 1:
            raise ValueError(f"startup_rate fraction must be between 0 and 1, but {startup_rate} is not.")
        if shutdown_rate <= 0  or shutdown_rate > 1:
            raise ValueError(f"shutdown_rate fraction must be between 0 and 1, but {shutdown_rate} is not.")
        if ramp_up_rate <= 0  or ramp_up_rate > 1:
            raise ValueError(f"ramp_up_rate fraction must be between 0 and 1, but {ramp_up_rate} is not.")
        if ramp_down_rate <= 0  or ramp_down_rate > 1:
            raise ValueError(f"ramp_down_rate fraction must be between 0 and 1, but {ramp_down_rate} is not.")
        if op_range_lb < 0  or op_range_lb > 1:
            raise ValueError(f"op_range_lb fraction must be between 0 and 1, but {op_range_lb} is not.")
        if op_range_lb >= shutdown_rate:
            raise ValueError(f"op_range_lb fraction must be <= shut_down_rate, otherwise the system cannot reach the off state.")


        start_up = {t: deepgetattr(self.mp_model.period[t], op_blk + ".startup") for t in self.mp_model.period}
        op_mode = {t: deepgetattr(self.mp_model.period[t], op_blk + ".op_mode") for t in self.mp_model.period}
        shut_down = {t: deepgetattr(self.mp_model.period[t], op_blk + ".shutdown") for t in self.mp_model.period}
        var_ramping = {t: deepgetattr(self.mp_model.period[t], op_blk + "." + ramping_var) for t in self.mp_model.period}

        if constraint_type == "linear":
            if isinstance(capacity_var,str):
                var_capacity = deepgetattr(self, design_blk + "." + capacity_var)
                act_startup_rate = {t: var_capacity * start_up[t] for t in self.mp_model.period}
                act_shutdown_rate = {t: var_capacity * shut_down[t] for t in self.mp_model.period}
                act_op_mode_rate = {t: var_capacity * op_mode[t] for t in self.mp_model.period}
            
            if isinstance(capacity_var, (int,float)):
                act_startup_rate = {t: capacity_var * start_up[t] for t in self.mp_model.period}
                act_shutdown_rate = {t: capacity_var * shut_down[t] for t in self.mp_model.period}
                act_op_mode_rate = {t: capacity_var * op_mode[t] for t in self.mp_model.period}
        elif constraint_type == "nonlinear":
            if linearization == True:
                if hasattr(self, "capacity_startup") == False:
                    for t in self.mp_model.period:
                        blk = deepgetattr(self.mp_model.period[t],op_blk)
                        blk._add_capacity_aux_vars()

                act_startup_rate= {t: deepgetattr(self.mp_model.period[t], op_blk + "." + capacity_var + "_startup") for t in self.mp_model.period}
                act_shutdown_rate= {t: deepgetattr(self.mp_model.period[t], op_blk + "." + capacity_var + "_shutdown") for t in self.mp_model.period}
                act_op_mode_rate = {t: deepgetattr(self.mp_model.period[t], op_blk + "." + capacity_var + "_op_mode") for t in self.mp_model.period}
            elif linearization == False:
                var_capacity = deepgetattr(self, design_blk + "." + capacity_var)
                act_startup_rate = {t: var_capacity * start_up[t] for t in self.mp_model.period}
                act_shutdown_rate = {t: var_capacity * shut_down[t] for t in self.mp_model.period}
                act_op_mode_rate = {t: var_capacity * op_mode[t] for t in self.mp_model.period}



       
        # Importing in the necessary variables
        if not hasattr(self, 'range_time_steps'):
            self.range_time_steps = RangeSet(len(self.mp_model.set_period))
        
        #Creating the pyomo block
        blk_name = op_blk.split(".")[-1] + "_rampup_rampdown"
        setattr(self.mp_model, blk_name, Block())
        blk = getattr(self.mp_model,blk_name)
        # The linearized ramping constraints
        @blk.Constraint(self.range_time_steps)
        def ramp_up_con(b,t):
            if t == 1:
                return Constraint.Skip
            else:
                return (var_ramping[self.mp_model.set_period.at(t)] - var_ramping[self.mp_model.set_period.at(t-1)] <= 
                        (startup_rate - ramp_up_rate) * act_startup_rate[self.mp_model.set_period.at(t)]  
                        + ramp_up_rate * act_op_mode_rate[self.mp_model.set_period.at(t)]
                        )

        @blk.Constraint(self.range_time_steps)
        def ramp_down_con(b,t):
            if t == 1:
                return Constraint.Skip
            else:
                return (
                        var_ramping[self.mp_model.set_period.at(t-1)]- var_ramping[self.mp_model.set_period.at(t)] <= 
                       shutdown_rate * act_shutdown_rate[self.mp_model.set_period.at(t)]
                       + ramp_down_rate * act_op_mode_rate[self.mp_model.set_period.at(t)]
                       )
        

    def add_startup_shutdown(
        self, 
        op_blk,
        design_blk, 
        build_binary_var,
        use_min_time=True, 
        up_time=1, 
        down_time=1,
    ):
        """
        Adds startup/shutdown and minimum uptime/downtime constraints on
        a given unit/process

        
        Args:
            op_blk:             String of the name of the operation model block, ex: ("fs.op_name")
            design_blk:         String of the name of the design model block, ex: ("m.design_name")
            build_binary_var:   String of the name of the binary variable which indicates if we 
                                should build (1) or not build (0) the design corresponding to the
                                'design_blk' referenced above
            use_min_time:    Boolean that can bypass constructing minimum up and downtime
                                constraints. (default: True --> constraints are constructed)
            up_time:            Time required for the system to start up fully 
                                    ex: 4 
            down_time:          Time required for the system to shutdown fully 
                                    ex: 4
        
        Returns:

        Assumption:
            up_time >= 1 & down_time >= 1
        """
        # Check up_time and down_time for validity
        if not isinstance(up_time, int):
            raise ValueError(f"up_time must be an integer, but {up_time} is not an integer")
        if up_time < 1:
            raise ValueError(f"up_time must be >= 1, but {up_time} is not")
        if not isinstance(down_time, int):
            raise ValueError(f"down_time must be an integer, but {down_time} is not an integer")
        if down_time < 1:
            raise ValueError(f"down_time must be >= 1, but {down_time} is not")
        
        

        start_up = {t: deepgetattr(self.mp_model.period[t], op_blk + ".startup") for t in self.mp_model.period}
        op_mode = {t: deepgetattr(self.mp_model.period[t], op_blk + ".op_mode") for t in self.mp_model.period}
        shut_down = {t: deepgetattr(self.mp_model.period[t], op_blk + ".shutdown") for t in self.mp_model.period}

        if design_blk is not None:
            build =  deepgetattr(self,design_blk + "." + build_binary_var)

        if not hasattr(self, 'range_time_steps'):
            self.range_time_steps = RangeSet(len(self.mp_model.set_period))
        number_time_steps = len(self.mp_model.set_period)

        blk_name = op_blk.split(".")[-1] + "_startup_shutdown"
        setattr(self.mp_model, blk_name, Block())
        blk = getattr(self.mp_model,blk_name)

        if design_blk is not None:
            @blk.Constraint(self.range_time_steps)
            def design_op_relationship_con(b, t):
                return (op_mode[self.mp_model.set_period.at(t)] <= build)
        
        @blk.Constraint(self.range_time_steps)
        def Binary_relationhsip_con(b, t):
            if t == 1 or t > number_time_steps:
                return Constraint.Skip
            return (op_mode[self.mp_model.set_period.at(t)] - op_mode[self.mp_model.set_period.at(t-1)] == 
                    start_up[self.mp_model.set_period.at(t-1)] - shut_down[self.mp_model.set_period.at(t)]
                   )

        if use_min_time:
            @blk.Constraint(self.range_time_steps)
            def minimum_up_time_con(b,t):
                if t == 1 or t < up_time or t > number_time_steps:
                    return Constraint.Skip
                else:
                    return sum(start_up[self.mp_model.set_period.at(i)] for i in range(t-up_time+1, t+1)) <= op_mode[self.mp_model.set_period.at(t)]
            
            @blk.Constraint(self.range_time_steps)
            def minimum_down_time_con(b, t):
                if t < down_time or t == 1 or t > number_time_steps:
                    return Constraint.Skip
                return (sum(shut_down[self.mp_model.set_period.at(i)] for i in range(t-down_time+1, t+1)) <= 
                    1 - op_mode[self.mp_model.set_period.at(t)])


    def build_hourly_cashflows(self, revenue_streams=['elec_revenue',], additional_costs=None):
        """
        Adds an expression for the net cash inflow for each operational
        block. This is the new cash inflow for each time period of the
        PriceTakerModel. Default costs for each model should include
        'non_fuel_vom' (non-fuel variable operating costs), 'fuel_cost'
        (cost of fuel), and 'carbon_price' (cost associated with producing
        carbon; i.e., a carbon tax). The net cash inflow is caluclated as:

            Sum(revenue streams) - Sum(costs)
        
        for every time period in the PriceTakerModel's MultiPeriodModel
        
        Args:
            revenue_streams:    List of strings representing the names of the 
                                revenue streams coming from the model. 
                                (default: ['elec_revenue',])
                                Coproduction example: ['elec_revenue', 
                                                       'H2_revenue', ]
            additional_costs:   List of strings representing the names of the
                                costs associated with operating at a time period.
                                (default: None)
                                example: ['hourly_fixed_cost',
                                          '',]

        Returns:

        """
        period = self.mp_model.period

        for p in period:
            non_fuel_vom = 0
            fuel_cost = 0
            elec_revenue = 0
            carbon_price = 0

            # ToDo: Add multiple revenue streams to the model
            # ToDo: Add additional costs for time period

            for blk in period[p].component_data_objects(Block):
                if isinstance(blk, OperationModelData):
                    non_fuel_vom += blk.non_fuel_vom
                    fuel_cost += blk.fuel_cost
                    elec_revenue += blk.elec_revenue
                    carbon_price += blk.carbon_price

            period[p].non_fuel_vom = Expression(expr=non_fuel_vom)
            period[p].fuel_cost = Expression(expr=fuel_cost)
            period[p].elec_revenue = Expression(expr=elec_revenue)
            period[p].carbon_price = Expression(expr=carbon_price)

            period[p].net_cash_inflow = Expression(
                expr=period[p].elec_revenue
                - period[p].non_fuel_vom
                - period[p].fuel_cost
                - period[p].carbon_price
            )

    def build_cashflows(
        self,
        lifetime=30,
        discount_rate=0.08,
        corp_tax=0.2,
        other_costs=0,
        other_revenue=0,
        objective="NPV",
    ):
        """
        Builds overall cashflow expressions and appends the objective function
        in terms of cashflows to the PriceTakerModel
        
        Args:
            lifetime:       Number of years (lifetime) to evaluate the equipment
            discount_rate:  Fractional rate of discount used in NPV calculations.
                            Must be between 0 and 1.
            corp_tax:       Fractional value of corporate tax used in NPV calculations.
                            ??? What are the restrictions on tax ???
            other_costs:    Additional costs that are not captured
            other_revenue:  
            objective:      

        Returns:

        """

        capex_expr = 0
        fom_expr = 0
        for blk in self.component_data_objects(Block):
            if isinstance(blk, DesignModelData):
                capex_expr += blk.capex
                fom_expr += blk.fom

        self.CAPEX = Var(within=NonNegativeReals, doc="Total CAPEX")
        self.capex_calculation = Constraint(expr=self.CAPEX == capex_expr)

        self.FOM = Var(within=NonNegativeReals, doc="Yearly Fixed O&M")
        self.fom_calculation = Constraint(expr=self.FOM == fom_expr)

        self.DEPRECIATION = Var(within=NonNegativeReals, doc="Yearly depreciation")
        self.dep_calculation = Constraint(
            expr=self.DEPRECIATION == self.CAPEX / lifetime
        )

        self.NET_CASH_INFLOW = Var(doc="Net cash inflow")
        self.net_cash_inflow_calculation = Constraint(
            expr=self.NET_CASH_INFLOW
            == sum(self.mp_model.period[p].net_cash_inflow for p in self.mp_model.period) # added period block name to the net_cash_inflow callout
        ) # added period block name to mp_model, to kthe len, and made a range list to loop over

        self.CORP_TAX = Var(within=NonNegativeReals, doc="Corporate tax")
        self.corp_tax_calculation = Constraint(
            expr=self.CORP_TAX
            >= corp_tax
            * (
                self.NET_CASH_INFLOW
                + other_revenue
                - other_costs
                - self.FOM
                - self.DEPRECIATION
            )
        )

        self.NET_PROFIT = Var(doc="Net profit after taxes")
        self.net_profit_calculation = Constraint(
            expr=self.NET_PROFIT
            == self.NET_CASH_INFLOW
            + other_revenue
            - other_costs
            - self.FOM
            - self.CORP_TAX
        )

        constant_cf_factor = (1 - (1 + discount_rate) ** (-lifetime)) / discount_rate
        self.NPV = Expression(expr=constant_cf_factor * self.NET_PROFIT - self.CAPEX)
        self.Annualized_NPV = Expression(
            expr=self.NET_PROFIT - (1 / constant_cf_factor) * self.CAPEX,
        )

        if objective == "NPV":
            self.obj = Objective(expr=self.NPV, sense=maximize)

        elif objective == "Annualized NPV":
            self.obj = Objective(expr=self.Annualized_NPV, sense=maximize)

        elif objective == "Net Profit":
            self.obj = Objective(expr=self.NET_PROFIT, sense=maximize)
