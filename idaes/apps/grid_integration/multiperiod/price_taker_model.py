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


class PriceTakerModel(ConcreteModel):
    @staticmethod
    def generate_daily_data(raw_data, day_list):

        daily_data = pd.DataFrame(columns=day_list)

        # Extracting data to populate empty dataframe
        i = 0
        j = 24
        day = 1
        while j <= len(raw_data):
            daily_data[day] = raw_data[i:j].reset_index(drop=True)
            i = j
            j = j + 24
            day = day + 1

        return daily_data

    @staticmethod
    def reconfigure_raw_data(self, raw_data):
        """
        Reconfigures the raw price series data into a usable form

        Args:
            raw_data: imported price series data

        Returns:
            daily_data: reconfigured price series data
        """
        # Get column headings
        column_head = raw_data.columns.tolist()
        # Remove the date/time column
        scenarios = column_head[1:]

        # Creating an empty dataframe to store daily data for clustering
        day_list = list(range(1, (len(raw_data) // 24) + 1))

        # Generate daily data
        for i in scenarios[0:1]:
            daily_data = self.generate_daily_data(
                raw_data=raw_data[i], day_list=day_list
            )

        return daily_data

    @staticmethod
    def get_optimal_n_clusters(daily_data, kmin=None, kmax=None, sample_weight=None):
        """
        Determines the appropriate number of clusters needed for a
        given price signal.

        Args:
            daily_data: reconfigured price series data
            kmin: minimum number of clusters
            kmax: maximum number of clusters
            sample_weight: applies a weight to each entry in daily_data

        Returns:
            n_clusters: the optimal number of clusters for the given data
        """
        if kmin is None:
            kmin = 1
        if kmax is None:
            # setting an arbitrary kmax. Maybe log warning that no kmax was set and default of 14 used?
            kmax = 14
        if kmin > kmax:
            _logger.error(f"kmin:{kmin} needs to be less than kmax:{kmax}.")

        if sample_weight is not None:
            if daily_data.shape != sample_weight.shape:
                _logger.error(
                    f"Ensure that the dimensions of the datasets match or set sample_weight to None."
                )

        k_values = range(kmin, kmax)
        inertia_values = []

        # Compute the inertia (SSE) for k clusters
        for k in k_values:
            kmeans = KMeans(n_clusters=k).fit(
                daily_data.transpose(), sample_weight=sample_weight
            )
            inertia_values.append(kmeans.inertia_)

        # Identify the "elbow point"
        elbow_point = KneeLocator(
            k_values, inertia_values, curve="convex", direction="decreasing"
        )
        n_clusters = elbow_point.knee

        print(f"Optimal # of clusters is: {n_clusters}")

        return n_clusters, inertia_values

    @staticmethod
    def get_elbow_plot(self, daily_data, kmin=None, kmax=None, sample_weight=None):
        """
        Determines the appropriate number of clusters needed for a
        given price signal.

        Args:
            daily_data: reconfigured price series data
            kmin: minimum number of clusters
            kmax: maximum number of clusters
            sample_weight: applies a weight to each entry in daily_data

        Returns:
            k_values: range of clusters from kmin to kmax (x-axis of elbow plot)
            inertia_values: within-cluster sum-of-squares (y-axis of elbow plot)
        """
        if kmin is None:
            kmin = 1
        if kmax is None:
            # setting an arbitrary kmax. Maybe log warning that no kmax was set and default of 14 used?
            kmax = 14

        k_values = range(kmin, kmax)
        n_clusters, inertia_values = self.get_optimal_n_clusters(
            daily_data, kmin, kmax, sample_weight
        )

        plt.show()
        plt.plot(k_values, inertia_values)
        plt.axvline(x=n_clusters, color="red", linestyle="--", label="Elbow")
        plt.xlabel("Number of clusters")
        plt.ylabel("Inertia")
        plt.title("Elbow Method")
        plt.xlim(kmin, kmax)
        plt.grid()
        plt.show()

        return k_values, inertia_values

    @staticmethod
    def cluster_lmp_data(data, n_clusters, horizon_length):
        """
        Clusters the given price signal in n_clusters. This method supports k-means, k-meteiod,...
        techniques for clustering.

        Args:

        Returns:
        """
        lmp_data = {1: {1: 2, 2: 3, 3: 5}, 2: {1: 2, 2: 3, 3: 5}}
        weights = {1: 45, 2: 56}

        return lmp_data, weights

    def append_lmp_data(
        self,
        filename,
        column_name="price",
        n_clusters=None,
        horizon_length=None,
    ):
        full_data = pd.read_csv(filename)

        if isinstance(column_name, list) and n_clusters is not None:
            # Multiple years and representative days
            self.set_years = [int(y) for y in column_name]
            self.set_days = RangeSet(n_clusters)
            self._n_time_points = horizon_length if horizon_length is not None else 24
            self.set_time = RangeSet(self._n_time_points)

            self.LMP = {}
            self.WEIGHTS = {}

            for year in column_name:
                price_data = full_data[year].to_list()
                lmp_data, weights = self.cluster_lmp_data(
                    price_data, n_clusters, horizon_length
                )
                y = int(year)

                for d in self.set_days:
                    for t in self.set_time:
                        self.LMP[t, d, y] = lmp_data[d][t]
                        self.WEIGHTS[d, y] = weights[d]

            return

        elif isinstance(column_name, list):
            # Multiple years, use fullyear price signal for each year
            self.set_years = [int(y) for y in column_name]
            self.set_days = None
            self._n_time_points = len(full_data)
            self.set_time = RangeSet(self._n_time_points)

            self.LMP = {}

            for year in column_name:
                price_data = full_data[year].to_list()
                y = int(year)

                for t in self.set_time:
                    self.LMP[t, y] = price_data[t - 1]

            return

        elif n_clusters is not None:
            # Single price signal, use reprentative days
            self.set_years = None
            self.set_days = RangeSet(n_clusters)
            self._n_time_points = horizon_length if horizon_length is not None else 24
            self.set_time = RangeSet(self._n_time_points)

            self.LMP = {}
            self.WEIGHTS = {}

            price_data = full_data[year].to_list()
            lmp_data, weights = self.cluster_lmp_data(
                price_data, n_clusters, horizon_length
            )

            for d in self.set_days:
                for t in self.set_time:
                    self.LMP[t, d] = lmp_data[d][t]
                    self.WEIGHTS[d] = weights[d]

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

        if not self.model_sets_available:
            raise Exception(
                "Model sets have not been defined. Run get_lmp_data to construct model sets"
            )

        self.mp_model = MultiPeriodModel(
            n_time_points=self._n_time_points,
            set_days=self.set_days,
            set_years=self.set_years,
            use_stochastic_build=True,
            **kwargs,
        )

    def add_ramping_constraints(self, var, ramp_up_limit, ramp_down_limit):
        """
        Adds ramping constraints of the form
        -ramp_down_limit <= var(t) - var(t-1) <= ramp_up_limit on var
        """
        pass

    def add_startup_shutdown(self):
        """
        Adds startup/shutdown and minimum uptime/downtime constraints on
        a given unit/process
        """
        pass

    def build_hourly_cashflows(self):
        period = self.mp_model.period

        for p in period:
            non_fuel_vom = 0
            fuel_cost = 0
            elec_revenue = 0
            carbon_price = 0

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
        Builds overall cashflow expressions and appends objective function
        to the model
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
            == sum(self.mp_model[p].net_cash_inflow for p in self.mp_model)
        )

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
