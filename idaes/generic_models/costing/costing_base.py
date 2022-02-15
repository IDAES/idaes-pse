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
Base classes for process costing
"""
from enum import Enum

import pyomo.environ as pyo
from pyomo.common.config import ConfigValue
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from pyomo.contrib.fbbt.fbbt import compute_bounds_on_expr

from idaes.core.process_base import (declare_process_block_class,
                                     ProcessBlockData)
from idaes.core.util.misc import add_object_reference
from idaes.core.util.exceptions import ConfigurationError

import idaes.logger as idaeslog


# Set up logger
_log = idaeslog.getLogger(__name__)


class DefaultCostingComponents(str, Enum):
    capital = "capital_cost"
    fixed = "fixed_operating_cost"
    variable = "variable_operating_cost"

    def __str__(self):
        return self.value


class CostingPackageBase():
    """
    Base Class for defining property packages.
    """

    # Placeholders for expected attributes
    currency_units = []

    # Set the base year for all costs
    base_currency = None
    base_period = pyo.units.year

    # Placeholder dict for defining material flows
    defined_flows = {}

    # Map costing methods to unit model classes
    unit_mapping = {}

    @staticmethod
    def build_global_params(self):
        """
        This is where any global parameters, such as Lang factors or
        coefficients for costing methods that should be shared across the
        process, should be declared. Sub-Blocks may be used ot help organize
        parameters if requried.
        """
        raise NotImplementedError(
            "Derived CostingPackage class has not defined a "
            "build_global_params method")

    @staticmethod
    def build_process_costs(self):
        """
        This is where process wide costing correlations should be declared.
        The following aggregate costs are available for use in calculating
        these process-wide costs:

            1. blk.aggregate_capital_cost
            2. blk.aggregate_fixed_operating_cost
            3. blk.aggregate_variable_operating_cost
            4. blk.aggregate_flow_costs (indexed by flow type)
        """
        raise NotImplementedError(
            "Derived CostingPackage class has not defined a "
            "build_process_costs method")

    @staticmethod
    def initialize(self):
        """
        This method allows users to definean initialization routine for any
        components createby the build_process_Costs method.

        Note that the aggregate costs will be initialized by the framework.
        """
        raise NotImplementedError(
            "Derived CostingPackage class has not defined an "
            "initialize method")


def is_costing_package(val):
    '''Domain validator for property package attributes

    Args:
        val : value to be checked

    Returns:
        ConfigurationError if val is not an instance of PhysicalParameterBlock
        or useDefault
    '''
    if issubclass(val, CostingPackageBase):
        return val
    else:
        _log.error(f"Costing package argument {val} should "
                   "be a sub-class of CostingPackageBase")
        raise ConfigurationError(
            f"Costing package argument {val} should "
            "be a sub-class of CostingPackageBase")


@declare_process_block_class("FlowsheetCostingBlock")
class FlowsheetCostingBlockData(ProcessBlockData):
    """
    IDAES Process Costing class.

    This class is used to create IDAES Flowsheet Costing Blocks, which are
    used to define capital and operating costs for a process. Each Flowsheet
    Costing Block is associated with a user-defined "costing package" module
    which must contain a library of methods for costing capital equipment
    as well as definitions of standard currency conversions and material and
    utility costs.
    """

    # Create Class ConfigBlock
    CONFIG = ProcessBlockData.CONFIG()
    CONFIG.declare("costing_package", ConfigValue(
        default=None,
        domain=is_costing_package,
        description="Costing package to be used"))

    def build(self):
        """
        Base build method for FlowsheetCostingBlocks.

        This method sets up the basic attributes expected for registering
        costing of unit operations and flows and calls the build_global_params
        method from the associated costing package.
        """
        super().build()

        # Check that costing package was assigned
        if self.config.costing_package is None:
            raise ConfigurationError(
                f"{self.name} - no costing_package was assigned - all "
                "FlowsheetCostingBlocks must be assigned a costing package.")

        # Set up attributes for registering units and flows
        self._registered_unit_models = []
        self.flow_types = pyo.Set()
        self._registered_flows = {}

        # Register unit mapping
        self._costing_methods_map = {}
        self._build_costing_methods_map()

        # Verify that costing package has set key attributes
        c_units = self.config.costing_package.base_currency
        if c_units is None:
            raise ValueError(
                f"{self.name} - costing package has not specified the base "
                "currency units to use for costing.")

        # Build global params
        self.config.costing_package.build_global_params(self)

        # Register pre-defined flow types
        for f, c in self.config.costing_package.defined_flows.items():
            self.register_flow_type(f, c)

    def cost_process(self):
        """
        This method constructs the process-level costing components based on
        the registered unit operations and flows.

        This first aggregates the costs from all the registered units and
        flows, and then calls the build_process_costs method from the
        associated costing package.
        """
        self.aggregate_costs()

        self.config.costing_package.build_process_costs(self)

    def cost_flow(self, flow_expr, flow_type):
        """
        This method registers a given flow component (Var or expression) for
        costing. All flows are required to be bounded to be non-negative (i.e.
        a lower bound equal to or greater than 0).

        Args:
            flow_expr - Pyomo Var or expression that represents a material flow
                        that should be included in the process costing. Units
                        are expected to be on a per time basis.
            flow_type - string identifying the material this flow represents.
                        This string must be registered with the
                        FlowsheetCostingBlock as a known flow type.

        Raises:
            ValueError if flow_type is not recognized.
            TypeError if flow_type is an indexed Var.
        """
        if flow_type not in self.flow_types:
            raise ValueError(
                f"{flow_type} is not a recognized flow type. Please check "
                "your spelling and that the flow type has been registered with"
                " the FlowsheetCostingBlock.")

        if isinstance(flow_expr, pyo.Var):
            # First, check that it is not indexed
            if flow_expr.is_indexed():
                raise TypeError(
                    f"{flow_expr.name} is an indexed Var. Flow costing only "
                    "supports unindexed Vars.")
        if hasattr(flow_expr, "lb"):
            if flow_expr.lb is None or flow_expr.lb < 0:
                raise ValueError(
                    f"{flow_expr.name} has a lower bound of less than zero. "
                    "Costing requires that all flows have a lower bound "
                    "equal to or greater than zero to avoid negative costs.")
        else:
            # Get bounds from expression
            ebounds = compute_bounds_on_expr(flow_expr)
            if ebounds[0] is None or ebounds[0] < 0:
                raise ValueError(
                    "flow_expr is an expression with a lower bound of less "
                    "than zero. Costing requires that all flows have a lower "
                    "bound equal to or greater than zero to avoid negative "
                    "costs.")

        self._registered_flows[flow_type].append(flow_expr)

    def cost_unit(self, unit_model, method=None, **kwargs):
        """
        This method registers the given unit model with the
        FlowsheetCostingBlock and constructs the associated Vars and
        Constraints. User may also provide a specific costing method to use for
        the unit, otherwise the framework will attempt to identify the correct
        method from the costing package mappings.

        An indexed entry for the given unit model is made in a
        UnitModelCostingBlock as a peer of the unit model (and a new
        UnitModelCostingBlock is created if required). The identified costing
        method is then called to populate this UnitModelCostingBlock with the
        required Vars and Constraints.

        Args:
            unit_model - a UnitModelBlock for which costing is to be performed.
            method - (optional) specific costing method to use when costing.
                     unit_model (default = None, use costing package mapping).
            kwargs - any additional arguments to be passed to method.

        Raises:
            RuntimeError if unit model already appears in Set of costed units.
            RuntimeError if method not provided and no default can be found in
            costing package mapping.
        """
        if method is None:
            method = self._get_costing_method_for(unit_model)

        parent = unit_model.parent_block()

        # Check to see if parent block has UnitCostingBlock already
        if not hasattr(parent, "unit_costing"):
            parent.unit_costing_set = pyo.Set()
            parent.unit_costing = UnitModelCostingBlock(
                parent.unit_costing_set)
        else:
            # Check if unit model is already in costing Set
            # We can skip this if the Set didn't exist before, as it is
            # obvious that costing hasn't been done yet
            if unit_model.local_name in parent.unit_costing_set:
                raise RuntimeError(
                    f"Unit model {unit_model.name} already appears in the "
                    "Set of costed units. Each unit model can only be costed "
                    "once.")

        # Register unit model with this costing package
        self._registered_unit_models.append(unit_model)

        # Add unit model to costing set
        parent.unit_costing_set.add(unit_model.local_name)

        # Assign obejct references for costing package and unit model to
        # UnitCostingBlockData
        add_object_reference(parent.unit_costing[unit_model.local_name],
                             "costing_package",
                             self)
        add_object_reference(parent.unit_costing[unit_model.local_name],
                             "unit_model",
                             unit_model)

        # Call unit costing method with unit_model
        method(parent.unit_costing[unit_model.local_name], **kwargs)

        # Check that costs are Vars and have lower bound of 0
        cost_vars = DefaultCostingComponents
        for v in cost_vars:
            try:
                cvar = getattr(parent.unit_costing[unit_model.local_name], v)
                if not isinstance(cvar, pyo.Var):
                    raise TypeError(
                        f"{unit_model.name} {v} component must be a Var. "
                        "Please check the costingpackage you are using to "
                        "ensure that all costing components are declared as "
                        "variables.")
                elif cvar.lb is None or cvar.lb < 0:
                    raise ValueError(
                        f"{unit_model.name} {v} component has a lower bound "
                        "less than zero. All costing components are required "
                        "to have lower bounds of 0 or greater to avoid "
                        "negative costs.")
            except AttributeError:
                pass

    def initialize(self):
        """
        This method attempts to initialize all the costing components
        registered with this FlowsheetCostingBlock. First, each registered
        UnitModelCostingBlock is initialized, followed by the aggregate
        costing variables and finally any process level costing components.

        Args:
            None
        """
        for u in self._registered_unit_models:
            # TODO: Make something more general for this to handle intermediate
            # Vars and Constraints
            parent = u.parent_block()
            cblock = parent.unit_costing[u.local_name]

            for c in DefaultCostingComponents:
                if hasattr(cblock, c):
                    var = getattr(cblock, c)
                    cons = getattr(cblock, f"{c}_constraint")
                    calculate_variable_from_constraint(var, cons)

        # Initialize aggregate cost vars
        for c in DefaultCostingComponents:
            var = getattr(self, f"aggregate_{c}")
            cons = getattr(self, f"aggregate_{c}_constraint")
            calculate_variable_from_constraint(var, cons)

        # Initialize aggregate flows and costs
        for f in self.flow_types:
            try:
                fvar = getattr(self, f"aggregate_flow_{f}")
                fconst = getattr(self, f"aggregate_flow_{f}_constraint")
                calculate_variable_from_constraint(fvar, fconst)

                calculate_variable_from_constraint(
                    self.aggregate_flow_costs[f],
                    self.aggregate_flow_costs_constraint[f])
            except AttributeError:
                self.aggregate_flow_costs[f].set_value(0)

        # Call costing package initialization
        try:
            # TODO: More general approach for initializing costing
            self.config.costing_package.initialize(self)
        except AttributeError:
            # Assume the package has no initialize method
            pass

    def register_flow_type(self, flow_type, cost):
        """
        This method allows users to register new material and utility flows
        with the FlowsheetCostingBlock for use when costing flows.

        Args:
            flow_type - string name to represent flow type
            cost - a Pyomo expression with units representing the flow cost
        """
        self.flow_types.add(flow_type)

        # Create a Var to hold the cost
        # Units will be different between flows, so have to use scalar Vars
        fvar = pyo.Var(units=pyo.units.get_units(cost),
                       doc=f"Cost parameter for {flow_type} flow")
        self.add_component(f"{flow_type}_cost", fvar)
        fvar.fix(cost)

        self._registered_flows[flow_type] = []

    def aggregate_costs(self):
        """
        This method aggregates costs from all the unit models and flows
        registered with this FlowsheetCostingBlock and creates aggregate
        variables for these on the FlowsheetCostingBlock that can be used for
        further process-wide costing calculations.

        The following costing variables are aggregated from all the registered
        UnitModelCostingBlocks (if they exist):
            * capital_cost,
            *fixed_operating_cost, and
            *variable_operating_cost

        Additionally, aggregate flow variables are created for all registered
        flow types along with aggregate costs associated with each of these.

        Args:
            None
        """
        c_units = self.config.costing_package.base_currency
        t_units = self.config.costing_package.base_period

        self.aggregate_capital_cost = pyo.Var(units=c_units)

        def agg_cap_cost_rule(blk):
            e = 0
            for u in self._registered_unit_models:
                parent = u.parent_block()
                cblock = parent.unit_costing[u.local_name]

                # Allow for units that might only have a subset of cost Vars
                if hasattr(cblock, "capital_cost"):
                    e += pyo.units.convert(cblock.capital_cost,
                                           to_units=c_units)

            return blk.aggregate_capital_cost == e

        self.aggregate_capital_cost_constraint = pyo.Constraint(
            rule=agg_cap_cost_rule)

        # Aggregate unit operating costs
        self.aggregate_fixed_operating_cost = pyo.Var(
            units=c_units/t_units)

        def agg_fixed_om_cost_rule(blk):
            e = 0
            for u in self._registered_unit_models:
                parent = u.parent_block()
                cblock = parent.unit_costing[u.local_name]

                # Allow for units that might only have a subset of cost Vars
                if hasattr(cblock, "fixed_operating_cost"):
                    e += pyo.units.convert(cblock.fixed_operating_cost,
                                           to_units=c_units/t_units)

            return blk.aggregate_fixed_operating_cost == e

        self.aggregate_fixed_operating_cost_constraint = pyo.Constraint(
            rule=agg_fixed_om_cost_rule)

        self.aggregate_variable_operating_cost = pyo.Var(
            units=c_units/t_units)

        def agg_var_om_cost_rule(blk):
            e = 0
            for u in self._registered_unit_models:
                parent = u.parent_block()
                cblock = parent.unit_costing[u.local_name]

                # Allow for units that might only have a subset of cost Vars
                if hasattr(cblock, "variable_operating_cost"):
                    e += pyo.units.convert(cblock.variable_operating_cost,
                                           to_units=c_units/t_units)

            return blk.aggregate_variable_operating_cost == e

        self.aggregate_variable_operating_cost_constraint = pyo.Constraint(
            rule=agg_var_om_cost_rule)

        # Aggregate flows
        # Units of flows will not be consistent, need separate Vars
        for f in self.flow_types:
            # We will use the first costed flow as representative of the whole
            if len(self._registered_flows[f]) > 0:
                f1 = self._registered_flows[f][0]
                funits = pyo.units.get_units(f1)
                agg_var = pyo.Var(units=funits,
                                  doc=f"Aggregate flow for {f}")
                self.add_component(f"aggregate_flow_{f}", agg_var)

                def agg_flow_rule(blk):
                    e = 0
                    for flow in self._registered_flows[f]:
                        e += pyo.units.convert(flow, to_units=funits)
                    agg_flow = getattr(blk, f"aggregate_flow_{f}")
                    return agg_flow == e

                agg_const = pyo.Constraint(rule=agg_flow_rule)

                self.add_component(f"aggregate_flow_{f}_constraint",
                                   agg_const)

        # TODO : More complex cost functions
        self.aggregate_flow_costs = pyo.Var(self.flow_types,
                                            units=c_units/t_units)

        def agg_flow_cost_rule(blk, ftype):
            try:
                agg_var = getattr(blk, f"aggregate_flow_{ftype}")
                cost_var = getattr(blk, f"{ftype}_cost")
                return blk.aggregate_flow_costs[ftype] == (
                    pyo.units.convert(agg_var * cost_var,
                                      to_units=c_units/t_units))
            except AttributeError:
                return blk.aggregate_flow_costs[ftype] == 0*c_units/t_units

        self.aggregate_flow_costs_constraint = pyo.Constraint(
            self.flow_types, rule=agg_flow_cost_rule)

    def del_unit_costing(self, unit_model):
        """
        This method unregisters the provided unit model from the current
        FlowsheetCostingBlock and delete the associated UnitModelCostingBlock
        element.

        Args:
            unit_model - unit model to be unregistered

        Raises:
            RuntimeError if unit_model is not registered with this
            FlowsheetCostingBlock
        """
        if unit_model not in self._registered_unit_models:
            raise RuntimeError(
                f"{unit_model.name} was not registered with this "
                "FlowsheetCostingBlock. del_unit_costing can only be used "
                "from the block with which the unit model is registered for "
                "costing.")

        parent = unit_model.parent_block()

        del parent.unit_costing[unit_model.local_name]
        parent.unit_costing_set.discard(unit_model.local_name)

        self._registered_unit_models.remove(unit_model)

    def _build_costing_methods_map(self):
        """
        This method takes the mapping of unit model classes to costing methods
        from the costing package and constructs a mapping of the associated
        unit model data classes which is used to look up default costing
        methods.
        """
        package_map = self.config.costing_package.unit_mapping

        for unit_class, meth in package_map.items():
            # Need to get unit data class from unit class
            data_class = unit_class._ComponentDataClass
            self._costing_methods_map[data_class] = meth

    def _get_costing_method_for(self, unit_model):
        """
        This method used the _costing_methods_map mapping to try to identify
        the default costing method to use for a given unit model based on
        the data class of the unit model.
        """
        for base in unit_model.__class__.__mro__:
            if base in self._costing_methods_map:
                return self._costing_methods_map[base]
        raise RuntimeError(
            f"Could not identify default costing method for {unit_model.name}."
            " This implies the unit model's class and parent classes do not "
            "exist in the default mapping provided by the costing package. "
            "Please provide a specific costing method for this unit.")

    # -------------------------------------------------------------------------
    def report(self):
        # Need a clean reporting function
        # This will need input from the costing package
        pass

    def display_registered_unit_models(self):
        # Need a method for cleanly showing registered unit models for costing
        pass

    def display_registered_flow_types(self):
        # Need a method for cleanly showing registered flow types
        pass

    def display_registered_flows(self):
        # Need a method for cleanly showing registered flow for costing
        pass


@declare_process_block_class("UnitModelCostingBlock")
class UnitModelCostingBlockData(ProcessBlockData):
    """
    Class for constructing costing blocks for unit models.

    At this stage, the only purpose of this class it to provide a distinct
    type for type-checking.
    """

    def build(self):
        super().build()
