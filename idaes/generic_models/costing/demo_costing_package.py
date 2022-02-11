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
"""
Demonstration costing package for testing costing proposal
"""
import pyomo.environ as pyo
from pyomo.util.calc_var_value import calculate_variable_from_constraint

from idaes.generic_models.unit_models import Heater
from idaes.core.util.misc import register_units_of_measurement

from idaes.costing.costing_base import CostingPackageBase


class MyCosting(CostingPackageBase):

    currency_units = [
        "USD2010 = [currency]",
        "USD2020 = 0.5 * USD2010"]  # Massive inflation
    # Need to register these units before we continue
    register_units_of_measurement(currency_units)

    # Set the base year for all costs
    base_currency = pyo.units.USD2020
    # Set a base period for all operating costs
    base_period = pyo.units.year

    # Define the list of known flow types this costing package will support
    defined_flows = {"steam": 0.2*pyo.units.USD2010/pyo.units.J}

    @staticmethod
    def build_global_params(blk):
        """
        This is where we can declare any global parameters we need, such as
        Lang factors, or coefficients for costing methods that should be
        shared across the process.

        You can do what you want here, so you could have e.g. sub-Blocks
        for ach costing method to separate the parameters for each method.
        """
        # Add a global capital cost factor
        blk.lang_factor = pyo.Param(initialize=4.2,
                                    units=pyo.units.dimensionless)

        # Add a costing parameter for Heaters
        blk.k1 = pyo.Param(initialize=42,
                           units=pyo.units.USD2020/pyo.units.kW)

    @staticmethod
    def build_process_costs(blk):
        """
        This is where you do all your process wide costing.
        This is completely up to you, but you will have access to the
        following aggregate costs:

            1. blk.aggregate_capital_cost
            2. blk.aggregate_fixed_operating_cost
            3. blk.aggregate_variable_operating_cost
            4. blk.aggregate_flow_costs (indexed by flow type)
        """
        # Add a total capital cost Var & Constraint
        blk.total_capital_cost = pyo.Var(units=MyCosting.base_currency)

        blk.total_capital_cost_constraint = pyo.Constraint(
            expr=blk.total_capital_cost ==
            blk.aggregate_capital_cost*blk.lang_factor)

    @staticmethod
    def initialize(blk):
        """
        Here we can add intializationsteps for the things we built in
        build_process_costs.

        Note that the aggregate costs will be initialized by the framework.
        """
        if hasattr(blk, "total_capital_cost"):
            calculate_variable_from_constraint(
                blk.total_capital_cost,
                blk.total_capital_cost_constraint)

    def cost_heater(blk):
        """
        This is where you do the unit level costing
        You can do pretty much whatever you want here, including
        intermediate Vars and Constraints. However, we do not have tools to
        initialize more complex costing methods (WIP).

        You have access to the follwoing references for getting things from
        the costing pckage parameters and unit model:

            * blk.costing_package is the FlowsheetCostingBlock
            * blk.unit_model is the UnitModel

        The framework will automatically aggregate the following attributes
        for you:

            1. capital_cost (units of [currency])
            2. fixed_operating_cost (units of [currency]/[time])
            2. variable_operating_cost (units of [currency]/[time])

        Note that all of the above are optional - if they are present they
        will be aggregated, but they will be skipped if not present.
        Note 2: currency units do not need to be consistent - the
        aggregator will handle any conversion to the specified
        base_currency units.
        """
        # $$ = k1*Q
        k1 = blk.costing_package.k1
        Q = blk.unit_model.heat_duty[0]
        # TODO : What about time-indexed models?
        # All capital design vars should be un-indexed, but no guarantee
        blk.capital_cost = pyo.Var(bounds=(0, 1e7),
                                   units=MyCosting.base_currency)
        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost ==
            pyo.units.convert(k1*Q, to_units=MyCosting.base_currency))

    # Map costing methods to unit model classes
    # Here we can provide a dict mapping unit model classes to costing methods
    # Even better, this is inheritance aware so e.g. Pump will first look for a
    # method assigned to Pump, and then fall back to PressureChanger
    unit_mapping = {Heater: cost_heater}
