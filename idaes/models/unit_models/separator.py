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
"""
General purpose separator block for IDAES models
"""

from enum import Enum
from functools import partial
from pandas import DataFrame

from pyomo.environ import (
    Block,
    check_optimal_termination,
    Constraint,
    Param,
    Reals,
    Reference,
    Set,
    Var,
    value,
)
from pyomo.network import Port
from pyomo.common.config import ConfigBlock, ConfigValue, In, ListOf, Bool

from idaes.core import (
    declare_process_block_class,
    UnitModelBlockData,
    useDefault,
    MaterialBalanceType,
    MomentumBalanceType,
    MaterialFlowBasis,
    VarLikeExpression,
)
from idaes.core.util.config import (
    is_physical_parameter_block,
    is_state_block,
)
from idaes.core.util.exceptions import (
    BurntToast,
    ConfigurationError,
    PropertyNotSupportedError,
    InitializationError,
)
from idaes.core.solvers import get_solver
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.logger as idaeslog
import idaes.core.util.scaling as iscale
from idaes.core.util.units_of_measurement import report_quantity
from idaes.core.initialization import ModularInitializerBase

__author__ = "Andrew Lee"


# Set up logger
_log = idaeslog.getLogger(__name__)


# Enumerate options for balances
class SplittingType(Enum):
    """
    Enum of supported material split types.
    """

    totalFlow = 1
    phaseFlow = 2
    componentFlow = 3
    phaseComponentFlow = 4


class EnergySplittingType(Enum):
    """
    Enum of support energy split types.
    """

    none = 0
    equal_temperature = 1
    equal_molar_enthalpy = 2
    enthalpy_split = 3


class SeparatorInitializer(ModularInitializerBase):
    """
    Initializer for Separator blocks.

    """

    def initialization_routine(
        self,
        model: Block,
    ):
        """
        Initialization routine for Separator Blocks.

        This routine starts by initializing the feed stream followed by solving
        for all split fractions. Next, states for each outlet are estimated using the
        following rules:

        1. Intensive states remain unchanged
        2. Extensive states are multiplied by split fractions if index matches, or
            average of split fractions for outlet otherwise

        Each outlet state block is then initialized, and finally the full model is
        solved.

        Args:
            model: model to be initialized

        Returns:
            None

        """
        init_log = idaeslog.getInitLogger(
            model.name, self.get_output_level(), tag="unit"
        )
        solve_log = idaeslog.getSolveLogger(
            model.name, self.get_output_level(), tag="unit"
        )

        # Create solver
        solver = get_solver(self.config.solver, self.config.solver_options)

        # Initialize mixed state block
        if model.config.mixed_state_block is not None:
            mblock = model.config.mixed_state_block
        else:
            mblock = model.mixed_state
        self.get_submodel_initializer(mblock).initialize(mblock)

        # Solve for split fractions only
        component_status = {}
        for c in model.component_objects((Block, Constraint)):
            for i in c:
                if not c[i].local_name == "sum_split_frac":
                    # Record current status of components to restore later
                    component_status[c[i]] = c[i].active
                    c[i].deactivate()

        res = None
        if degrees_of_freedom(model) != 0:
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                res = solver.solve(model, tee=slc.tee)
                init_log.info(
                    "Initialization Step 1 Complete: {}".format(idaeslog.condition(res))
                )

        for c, s in component_status.items():
            if s:
                c.activate()

        if model.config.ideal_separation:
            # If using ideal splitting, initialization should be complete
            return res

        # Initialize outlet StateBlocks
        outlet_list = model.create_outlet_list()

        # Premises for initializing outlet states:
        # 1. Intensive states remain unchanged - this is either a valid premise
        # or the actual state is impossible to calculate without solving the
        # full separator model.
        # 2. Extensive states use split fractions if index matches, or
        # average of split fractions for outlet otherwise
        for o in outlet_list:
            # Get corresponding outlet StateBlock
            o_block = getattr(model, o + "_state")

            # Create dict to store fixed status of state variables
            for t in model.flowsheet().time:
                # Calculate values for state variables
                s_vars = o_block[t].define_state_vars()
                for var_name_port, var_obj in s_vars.items():
                    for k in var_obj:
                        # If fixed, use current value
                        # otherwise calculate guess from mixed state and fix
                        if not var_obj[k].fixed:
                            m_var = getattr(mblock[t], var_obj.local_name)
                            if "flow" in var_name_port:
                                # If a "flow" variable, is extensive
                                # Apply split fraction
                                if model.config.split_basis == SplittingType.totalFlow:
                                    # All flows split by outlet
                                    var_obj[k].set_value(
                                        value(m_var[k] * model.split_fraction[(t, o)])
                                    )
                                elif "_phase_comp" in var_name_port:
                                    # Need to match indices, but use split frac
                                    if (
                                        model.config.split_basis
                                        == SplittingType.phaseComponentFlow
                                    ):
                                        var_obj[k].set_value(
                                            value(
                                                m_var[k]
                                                * model.split_fraction[(t, o) + (k,)]
                                            )
                                        )
                                    elif (
                                        model.config.split_basis
                                        == SplittingType.phaseFlow
                                    ):
                                        var_obj[k].set_value(
                                            value(
                                                m_var[k]
                                                * model.split_fraction[(t, o) + (k[0],)]
                                            )
                                        )
                                    elif (
                                        model.config.split_basis
                                        == SplittingType.componentFlow
                                    ):
                                        var_obj[k].set_value(
                                            value(
                                                m_var[k]
                                                * model.split_fraction[(t, o) + (k[1],)]
                                            )
                                        )
                                    else:
                                        raise BurntToast(
                                            "{} encountered unrecognised "
                                            "SplittingType. This should not "
                                            "occur - please send this bug to "
                                            "the IDAES developers.".format(model.name)
                                        )
                                elif "_phase" in var_name_port:
                                    if (
                                        model.config.split_basis
                                        == SplittingType.phaseComponentFlow
                                    ):
                                        # Need average split fraction
                                        avg_split = value(
                                            sum(
                                                model.split_fraction[t, o, k, j]
                                                for j in mblock.component_list
                                            )
                                            / len(mblock.component_list)
                                        )
                                        var_obj[k].set_value(
                                            value(m_var[k] * avg_split)
                                        )
                                    elif (
                                        model.config.split_basis
                                        == SplittingType.phaseFlow
                                    ):
                                        var_obj[k].set_value(
                                            value(
                                                m_var[k]
                                                * model.split_fraction[(t, o) + (k,)]
                                            )
                                        )
                                    elif (
                                        model.config.split_basis
                                        == SplittingType.componentFlow
                                    ):
                                        # Need average split fraction
                                        avg_split = value(
                                            sum(
                                                model.split_fraction[t, o, j]
                                                for j in mblock.component_list
                                            )
                                            / len(mblock.component_list)
                                        )
                                        var_obj[k].set_value(
                                            value(m_var[k] * avg_split)
                                        )
                                    else:
                                        raise BurntToast(
                                            "{} encountered unrecognised "
                                            "SplittingType. This should not "
                                            "occur - please send this bug to "
                                            "the IDAES developers.".format(model.name)
                                        )
                                elif "_comp" in var_name_port:
                                    if (
                                        model.config.split_basis
                                        == SplittingType.phaseComponentFlow
                                    ):
                                        # Need average split fraction
                                        avg_split = value(
                                            sum(
                                                model.split_fraction[t, o, p, k]
                                                for p in mblock.phase_list
                                            )
                                            / len(mblock.phase_list)
                                        )
                                        var_obj[k].set_value(
                                            value(m_var[k] * avg_split)
                                        )
                                    elif (
                                        model.config.split_basis
                                        == SplittingType.phaseFlow
                                    ):
                                        # Need average split fraction
                                        avg_split = value(
                                            sum(
                                                model.split_fraction[t, o, p]
                                                for p in mblock.phase_list
                                            )
                                            / len(mblock.phase_list)
                                        )
                                        var_obj[k].set_value(
                                            value(m_var[k] * avg_split)
                                        )
                                    elif (
                                        model.config.split_basis
                                        == SplittingType.componentFlow
                                    ):
                                        var_obj[k].set_value(
                                            value(
                                                m_var[k]
                                                * model.split_fraction[(t, o) + (k,)]
                                            )
                                        )
                                    else:
                                        raise BurntToast(
                                            "{} encountered unrecognised "
                                            "SplittingType. This should not "
                                            "occur - please send this bug to "
                                            "the IDAES developers.".format(model.name)
                                        )
                                else:
                                    # Assume unindexed extensive state
                                    # Need average split
                                    if (
                                        model.config.split_basis
                                        == SplittingType.phaseComponentFlow
                                    ):
                                        # Need average split fraction
                                        avg_split = value(
                                            sum(
                                                model.split_fraction[t, o, p, j]
                                                for (p, j) in mblock.phase_component_set
                                            )
                                            / len(mblock.phase_component_set)
                                        )
                                    elif (
                                        model.config.split_basis
                                        == SplittingType.phaseFlow
                                    ):
                                        # Need average split fraction
                                        avg_split = value(
                                            sum(
                                                model.split_fraction[t, o, p]
                                                for p in mblock.phase_list
                                            )
                                            / len(mblock.phase_list)
                                        )
                                    elif (
                                        model.config.split_basis
                                        == SplittingType.componentFlow
                                    ):
                                        # Need average split fraction
                                        avg_split = value(
                                            sum(
                                                model.split_fraction[t, o, j]
                                                for j in mblock.component_list
                                            )
                                            / len(mblock.component_list)
                                        )
                                    else:
                                        raise BurntToast(
                                            "{} encountered unrecognised "
                                            "SplittingType. This should not "
                                            "occur - please send this bug to "
                                            "the IDAES developers.".format(model.name)
                                        )
                                    var_obj[k].set_value(value(m_var[k] * avg_split))
                            else:
                                # Otherwise intensive, equate to mixed stream
                                var_obj[k].set_value(m_var[k].value)

            # Call initialization routine for outlet StateBlock
            self.get_submodel_initializer(o_block).initialize(o_block)

        if model.config.mixed_state_block is None:
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                res = solver.solve(model, tee=slc.tee)

            if not check_optimal_termination(res):
                raise InitializationError(
                    f"{model.name} failed to initialize successfully. Please "
                    f"check the output logs for more information."
                )

            init_log.info(
                "Initialization Step 2 Complete: {}".format(idaeslog.condition(res))
            )
        else:
            init_log.info("Initialization Complete.")

        return res


@declare_process_block_class("Separator")
class SeparatorData(UnitModelBlockData):
    """
    This is a general purpose model for a Separator block with the IDAES
    modeling framework. This block can be used either as a stand-alone
    Separator unit operation, or as a sub-model within another unit operation.

    This model creates a number of StateBlocks to represent the outgoing
    streams, then writes a set of phase-component material balances, an
    overall enthalpy balance (2 options), and a momentum balance (2 options)
    linked to a mixed-state StateBlock. The mixed-state StateBlock can either
    be specified by the user (allowing use as a sub-model), or created by the
    Separator.

    When being used as a sub-model, Separator should only be used when a
    set of new StateBlocks are required for the streams to be separated. It
    should not be used to separate streams to go to multiple ControlVolumes in a
    single unit model - in these cases the unit model developer should write
    their own splitting equations.
    """

    default_initializer = SeparatorInitializer

    CONFIG = ConfigBlock()
    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""Indicates whether this model will be dynamic or not,
**default** = False. Product blocks are always steady-state.""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be False",
            doc="""Product blocks do not contain holdup, thus this must be
False.""",
        ),
    )
    CONFIG.declare(
        "property_package",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for mixer",
            doc="""Property parameter object used to define property
calculations,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PropertyParameterObject** - a PropertyParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "property_package_args",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing property packages",
            doc="""A ConfigBlock with arguments to be passed to a property
block(s) and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}""",
        ),
    )
    CONFIG.declare(
        "outlet_list",
        ConfigValue(
            domain=ListOf(str),
            description="List of outlet names",
            doc="""A list containing names of outlets,
**default** - None.
**Valid values:** {
**None** - use num_outlets argument,
**list** - a list of names to use for outlets.}""",
        ),
    )
    CONFIG.declare(
        "num_outlets",
        ConfigValue(
            domain=int,
            description="Number of outlets to unit",
            doc="""Argument indicating number (int) of outlets to construct,
not used if outlet_list arg is provided,
**default** - None.
**Valid values:** {
**None** - use outlet_list arg instead, or default to 2 if neither argument
provided,
**int** - number of outlets to create (will be named with sequential integers
from 1 to num_outlets).}""",
        ),
    )
    CONFIG.declare(
        "split_basis",
        ConfigValue(
            default=SplittingType.totalFlow,
            domain=SplittingType,
            description="Basis for splitting stream",
            doc="""Argument indicating basis to use for splitting mixed stream,
**default** - SplittingType.totalFlow.
**Valid values:** {
**SplittingType.totalFlow** - split based on total flow (split
fraction indexed only by time and outlet),
**SplittingType.phaseFlow** - split based on phase flows (split fraction
indexed by time, outlet and phase),
**SplittingType.componentFlow** - split based on component flows (split
fraction indexed by time, outlet and components),
**SplittingType.phaseComponentFlow** - split based on phase-component flows (
split fraction indexed by both time, outlet, phase and components).}""",
        ),
    )
    CONFIG.declare(
        "material_balance_type",
        ConfigValue(
            default=MaterialBalanceType.useDefault,
            domain=In(MaterialBalanceType),
            description="Material balance construction flag",
            doc="""Indicates what type of mass balance should be constructed,
**default** - MaterialBalanceType.useDefault.
**Valid values:** {
**MaterialBalanceType.useDefault - refer to property package for default
balance type
**MaterialBalanceType.none** - exclude material balances,
**MaterialBalanceType.componentPhase** - use phase component balances,
**MaterialBalanceType.componentTotal** - use total component balances,
**MaterialBalanceType.elementTotal** - use total element balances,
**MaterialBalanceType.total** - use total material balance.}""",
        ),
    )
    CONFIG.declare(
        "momentum_balance_type",
        ConfigValue(
            default=MomentumBalanceType.pressureTotal,
            domain=In(MomentumBalanceType),
            description="Momentum balance construction flag",
            doc="""Indicates what type of momentum balance should be constructed,
    **default** - MomentumBalanceType.pressureTotal.
    **Valid values:** {
    **MomentumBalanceType.none** - exclude momentum balances,
    **MomentumBalanceType.pressureTotal** - pressure in all outlets is equal,
    **MomentumBalanceType.pressurePhase** - not yet supported,
    **MomentumBalanceType.momentumTotal** - not yet supported,
    **MomentumBalanceType.momentumPhase** - not yet supported.}""",
        ),
    )
    CONFIG.declare(
        "has_phase_equilibrium",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Calculate phase equilibrium in mixed stream",
            doc="""Argument indicating whether phase equilibrium should be
calculated for the resulting mixed stream,
**default** - False.
**Valid values:** {
**True** - calculate phase equilibrium in mixed stream,
**False** - do not calculate equilibrium in mixed stream.}""",
        ),
    )
    CONFIG.declare(
        "energy_split_basis",
        ConfigValue(
            default=EnergySplittingType.equal_temperature,
            domain=EnergySplittingType,
            description="Type of constraint to write for energy splitting",
            doc="""Argument indicating basis to use for splitting energy this
is not used for when ideal_separation == True.
**default** - EnergySplittingType.equal_temperature.
**Valid values:** {
**EnergySplittingType.none** - no energy balance constraints,
**EnergySplittingType.equal_temperature** - outlet temperatures equal inlet,
**EnergySplittingType.equal_molar_enthalpy** - outlet molar enthalpies equal
inlet,
**EnergySplittingType.enthalpy_split** - apply split fractions to enthalpy
flows. Does not work with component or phase-component splitting.}""",
        ),
    )
    CONFIG.declare(
        "ideal_separation",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Ideal splitting flag",
            doc="""Argument indicating whether ideal splitting should be used.
Ideal splitting assumes perfect spearation of material, and attempts to
avoid duplication of StateBlocks by directly partitioning outlet flows to
ports,
**default** - False.
**Valid values:** {
**True** - use ideal splitting methods. Cannot be combined with
has_phase_equilibrium = True,
**False** - use explicit splitting equations with split fractions.}""",
        ),
    )
    CONFIG.declare(
        "ideal_split_map",
        ConfigValue(
            domain=dict,
            description="Ideal splitting partitioning map",
            doc="""Dictionary containing information on how extensive variables
should be partitioned when using ideal splitting (ideal_separation = True).
**default** - None.
**Valid values:** {
**dict** with keys of indexing set members and values indicating which outlet
this combination of keys should be partitioned to.
E.g. {("Vap", "H2"): "outlet_1"}}""",
        ),
    )
    CONFIG.declare(
        "mixed_state_block",
        ConfigValue(
            domain=is_state_block,
            description="Existing StateBlock to use as mixed stream",
            doc="""An existing state block to use as the source stream from the
Separator block,
**default** - None.
**Valid values:** {
**None** - create a new StateBlock for the mixed stream,
**StateBlock** - a StateBock to use as the source for the mixed stream.}""",
        ),
    )
    CONFIG.declare(
        "construct_ports",
        ConfigValue(
            default=True,
            domain=Bool,
            description="Construct inlet and outlet Port objects",
            doc="""Argument indicating whether model should construct Port
objects linked the mixed state and all outlet states,
**default** - True.
**Valid values:** {
**True** - construct Ports for all states,
**False** - do not construct Ports.""",
        ),
    )

    def build(self):
        """
        General build method for SeparatorData. This method calls a number
        of sub-methods which automate the construction of expected attributes
        of unit models.

        Inheriting models should call `super().build`.

        Args:
            None

        Returns:
            None
        """
        # Call super.build()
        super(SeparatorData, self).build()

        self._validate_config_arguments()

        # Call setup methods from ControlVolumeBlockData
        self._get_property_package()
        self._get_indexing_sets()

        # Create list of inlet names
        outlet_list = self.create_outlet_list()

        if self.config.mixed_state_block is None:
            mixed_block = self.add_mixed_state_block()
        else:
            mixed_block = self.get_mixed_state_block()

        # Add inlet port
        self.add_inlet_port_objects(mixed_block)

        # Construct splitter based on ideal_separation argument
        if self.config.ideal_separation:
            # Use ideal partitioning method
            self.partition_outlet_flows(mixed_block, outlet_list)
        else:
            # Otherwise, Build StateBlocks for outlet
            outlet_blocks = self.add_outlet_state_blocks(outlet_list)

            # Add split fractions
            self.add_split_fractions(outlet_list, mixed_block)

            # Construct splitting equations
            self.add_material_splitting_constraints(mixed_block)
            self.add_energy_splitting_constraints(mixed_block)
            self.add_momentum_splitting_constraints(mixed_block)

            # Construct outlet port objects
            self.add_outlet_port_objects(outlet_list, outlet_blocks)

    def _validate_config_arguments(self):
        if self.config.has_phase_equilibrium and self.config.ideal_separation:
            raise ConfigurationError(
                """{} received arguments has_phase_equilibrium = True and
                    ideal_separation = True. These arguments are incompatible
                    with each other, and you should choose one or the other.""".format(
                    self.name
                )
            )

    def create_outlet_list(self):
        """
        Create list of outlet stream names based on config arguments.

        Returns:
            list of strings
        """
        if self.config.outlet_list is not None and self.config.num_outlets is not None:
            # If both arguments provided and not consistent, raise Exception
            if len(self.config.outlet_list) != self.config.num_outlets:
                raise ConfigurationError(
                    "{} Separator provided with both outlet_list and "
                    "num_outlets arguments, which were not consistent ("
                    "length of outlet_list was not equal to num_outlets). "
                    "Please check your arguments for consistency, and "
                    "note that it is only necessry to provide one of "
                    "these arguments.".format(self.name)
                )
        elif self.config.outlet_list is None and self.config.num_outlets is None:
            # If no arguments provided for outlets, default to num_outlets = 2
            self.config.num_outlets = 2

        # Create a list of names for outlet StateBlocks
        if self.config.outlet_list is not None:
            outlet_list = self.config.outlet_list
        else:
            outlet_list = [
                "outlet_" + str(n) for n in range(1, self.config.num_outlets + 1)
            ]

        return outlet_list

    def add_outlet_state_blocks(self, outlet_list):
        """
        Construct StateBlocks for all outlet streams.

        Args:
            list of strings to use as StateBlock names

        Returns:
            list of StateBlocks
        """
        # Setup StateBlock argument dict
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["defined_state"] = False

        # Create empty list to hold StateBlocks for return
        outlet_blocks = []

        # Create an instance of StateBlock for all outlets
        for o in outlet_list:
            o_obj = self.config.property_package.build_state_block(
                self.flowsheet().time, doc="Material properties at outlet", **tmp_dict
            )

            setattr(self, o + "_state", o_obj)

            outlet_blocks.append(getattr(self, o + "_state"))

        return outlet_blocks

    def add_mixed_state_block(self):
        """
        Constructs StateBlock to represent mixed stream.

        Returns:
            New StateBlock object
        """
        # Setup StateBlock argument dict
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["defined_state"] = True

        self.mixed_state = self.config.property_package.build_state_block(
            self.flowsheet().time, doc="Material properties of mixed stream", **tmp_dict
        )

        return self.mixed_state

    def get_mixed_state_block(self):
        """
        Validates StateBlock provided in user arguments for mixed stream.

        Returns:
            The user-provided StateBlock or an Exception
        """
        # Sanity check to make sure method is not called when arg missing
        if self.config.mixed_state_block is None:
            try:
                return self.mixed_state
            except AttributeError:
                raise BurntToast(
                    f"{self.name} get_mixed_state_block method called when the "
                    "mixed_state_block argument is None, and no mixed_state "
                    "block is contained in separator. This should not happen."
                )
        # Check that the user-provided StateBlock uses the same prop pack
        if (
            self.config.mixed_state_block[
                self.flowsheet().time.first()
            ].config.parameters
            != self.config.property_package
        ):
            raise ConfigurationError(
                "{} StateBlock provided in mixed_state_block argument "
                " does not come from the same property package as "
                "provided in the property_package argument. All "
                "StateBlocks within a Separator must use the same "
                "property package.".format(self.name)
            )

        return self.config.mixed_state_block

    def add_inlet_port_objects(self, mixed_block):
        """
        Adds inlet Port object if required.

        Args:
            a mixed state StateBlock object

        Returns:
            None
        """
        if self.config.construct_ports is True:
            self.add_port(name="inlet", block=mixed_block, doc="Inlet Port")

    def add_outlet_port_objects(self, outlet_list, outlet_blocks):
        """
        Adds outlet Port objects if required.

        Args:
            a list of outlet StateBlock objects

        Returns:
            None
        """
        if self.config.construct_ports is True:
            # Add ports
            for p in outlet_list:
                o_state = getattr(self, p + "_state")
                self.add_port(name=p, block=o_state, doc="Outlet Port")

    def add_split_fractions(self, outlet_list, mixed_block):
        """
        Creates outlet Port objects and tries to partition mixed stream flows
        between these

        Args:
            StateBlock representing the mixed flow to be split
            a list of names for outlets

        Returns:
            None
        """
        self.outlet_idx = Set(initialize=outlet_list)
        pc_set = mixed_block.phase_component_set

        if self.config.split_basis == SplittingType.totalFlow:
            sf_idx = [self.flowsheet().time, self.outlet_idx]
            sf_sum_idx = [self.flowsheet().time]
        elif self.config.split_basis == SplittingType.phaseFlow:
            sf_idx = [
                self.flowsheet().time,
                self.outlet_idx,
                mixed_block.phase_list,
            ]
            sf_sum_idx = [
                self.flowsheet().time,
                mixed_block.phase_list,
            ]
        elif self.config.split_basis == SplittingType.componentFlow:
            sf_idx = [
                self.flowsheet().time,
                self.outlet_idx,
                mixed_block.component_list,
            ]
            sf_sum_idx = [
                self.flowsheet().time,
                mixed_block.component_list,
            ]
        elif self.config.split_basis == SplittingType.phaseComponentFlow:
            sf_idx = [
                self.flowsheet().time,
                self.outlet_idx,
                pc_set,
            ]
            sf_sum_idx = [
                self.flowsheet().time,
                pc_set,
            ]
        else:
            raise BurntToast(
                "{} split_basis has unexpected value. This "
                "should not happen.".format(self.name)
            )

        # Create split fraction variable
        self.split_fraction = Var(*sf_idx, initialize=0.5, doc="Outlet split fractions")

        # Add constraint that split fractions sum to 1
        def sum_sf_rule(b, t, *args):
            return 1 == sum(b.split_fraction[t, o, args] for o in self.outlet_idx)

        self.sum_split_frac = Constraint(*sf_sum_idx, rule=sum_sf_rule)

    def add_material_splitting_constraints(self, mixed_block):
        """
        Creates constraints for splitting the material flows
        """
        pc_set = mixed_block.phase_component_set

        def sf(t, o, p, j):
            if self.config.split_basis == SplittingType.totalFlow:
                return self.split_fraction[t, o]
            elif self.config.split_basis == SplittingType.phaseFlow:
                return self.split_fraction[t, o, p]
            elif self.config.split_basis == SplittingType.componentFlow:
                return self.split_fraction[t, o, j]
            elif self.config.split_basis == SplittingType.phaseComponentFlow:
                return self.split_fraction[t, o, p, j]

        mb_type = self.config.material_balance_type
        if mb_type == MaterialBalanceType.useDefault:
            t_ref = self.flowsheet().time.first()
            mb_type = mixed_block[t_ref].default_material_balance_type()

        # Get units from property package
        units_meta = mixed_block.params.get_metadata()
        flow_basis = mixed_block[
            self.flowsheet().time.first()
        ].get_material_flow_basis()
        if flow_basis == MaterialFlowBasis.molar:
            flow_units = units_meta.get_derived_units("flow_mole")
        elif flow_basis == MaterialFlowBasis.mass:
            flow_units = units_meta.get_derived_units("flow_mass")
        else:
            # Let this pass for now with no units
            flow_units = None

        if mixed_block.include_inherent_reactions:
            # Add extents of reaction and stoichiometric constraints
            # TODO: It would be nice if there was a way to recognise that
            # a total flow split does not require calculation of equilibrium
            # (assuming T&P are constant) and turn of inherent reactions in
            # that case, but that feature does not exist in property packages
            self.inherent_reaction_extent = Var(
                self.flowsheet().time,
                self.outlet_idx,
                mixed_block.params.inherent_reaction_idx,
                domain=Reals,
                initialize=0.0,
                doc="Extent of inherent reactions in outlets",
                units=flow_units,
            )

            self.inherent_reaction_generation = Var(
                self.flowsheet().time,
                self.outlet_idx,
                pc_set,
                domain=Reals,
                initialize=0.0,
                doc="Generation due to inherent reactions in outlets",
                units=flow_units,
            )

            @self.Constraint(
                self.flowsheet().time,
                self.outlet_idx,
                pc_set,
            )
            def inherent_reaction_constraint(b, t, o, p, j):
                if (p, j) in pc_set:
                    return b.inherent_reaction_generation[t, o, p, j] == (
                        sum(
                            mixed_block[t].params.inherent_reaction_stoichiometry[
                                r, p, j
                            ]
                            * self.inherent_reaction_extent[t, o, r]
                            for r in mixed_block[t].params.inherent_reaction_idx
                        )
                    )
                return Constraint.Skip

        if mb_type == MaterialBalanceType.componentPhase:
            if self.config.has_phase_equilibrium is True:
                try:
                    self.phase_equilibrium_generation = Var(
                        self.flowsheet().time,
                        self.outlet_idx,
                        self.config.property_package.phase_equilibrium_idx,
                        domain=Reals,
                        doc="Amount of generation in unit by phase equilibria",
                        units=flow_units,
                    )
                except AttributeError:
                    raise PropertyNotSupportedError(
                        "{} Property package does not contain a list of phase "
                        "equilibrium reactions (phase_equilibrium_idx), thus "
                        "does not support phase equilibrium.".format(self.name)
                    )

            @self.Constraint(
                self.flowsheet().time,
                self.outlet_idx,
                pc_set,
                doc="Material splitting equations",
            )
            def material_splitting_eqn(b, t, o, p, j):
                o_block = getattr(self, o + "_state")
                lhs = sf(t, o, p, j) * mixed_block[t].get_material_flow_terms(p, j)

                rhs = o_block[t].get_material_flow_terms(p, j)

                if self.config.has_phase_equilibrium:
                    rhs += -sum(
                        b.phase_equilibrium_generation[t, o, r]
                        for r in b.config.property_package.phase_equilibrium_idx
                        if mixed_block[t].params.phase_equilibrium_list[r][0] == j
                        and mixed_block[t].params.phase_equilibrium_list[r][1][0] == p
                    ) + sum(
                        b.phase_equilibrium_generation[t, o, r]
                        for r in b.config.property_package.phase_equilibrium_idx
                        if mixed_block[t].params.phase_equilibrium_list[r][0] == j
                        and mixed_block[t].params.phase_equilibrium_list[r][1][1] == p
                    )

                if mixed_block.include_inherent_reactions:
                    rhs += b.inherent_reaction_generation[t, o, p, j]

                return lhs == rhs

        elif mb_type == MaterialBalanceType.componentTotal:

            @self.Constraint(
                self.flowsheet().time,
                self.outlet_idx,
                mixed_block.component_list,
                doc="Material splitting equations",
            )
            def material_splitting_eqn(b, t, o, j):
                o_block = getattr(self, o + "_state")

                lhs = sum(
                    sf(t, o, p, j) * mixed_block[t].get_material_flow_terms(p, j)
                    for p in mixed_block.phase_list
                    if (p, j) in pc_set
                )
                rhs = sum(
                    o_block[t].get_material_flow_terms(p, j)
                    for p in o_block.phase_list
                    if (p, j) in pc_set
                )

                if mixed_block.include_inherent_reactions:
                    rhs += sum(
                        b.inherent_reaction_generation[t, o, p, j]
                        for p in o_block.phase_list
                        if (p, j) in pc_set
                    )

                return lhs == rhs

        elif mb_type == MaterialBalanceType.total:

            @self.Constraint(
                self.flowsheet().time,
                self.outlet_idx,
                doc="Material splitting equations",
            )
            def material_splitting_eqn(b, t, o):
                o_block = getattr(self, o + "_state")

                lhs = sum(
                    sum(
                        sf(t, o, p, j) * mixed_block[t].get_material_flow_terms(p, j)
                        for j in mixed_block.component_list
                        if (p, j) in pc_set
                    )
                    for p in mixed_block.phase_list
                )
                rhs = sum(
                    sum(
                        o_block[t].get_material_flow_terms(p, j)
                        for j in mixed_block.component_list
                        if (p, j) in pc_set
                    )
                    for p in o_block.phase_list
                )

                if mixed_block.include_inherent_reactions:
                    rhs += sum(
                        sum(
                            b.inherent_reaction_generation[t, o, p, j]
                            for j in mixed_block.component_list
                            if (p, j) in pc_set
                        )
                        for p in o_block.phase_list
                    )

                return lhs == rhs

        elif mb_type == MaterialBalanceType.elementTotal:
            raise ConfigurationError(
                "{} Separators do not support elemental "
                "material balances.".format(self.name)
            )
        elif mb_type == MaterialBalanceType.none:
            pass
        else:
            raise BurntToast(
                "{} Separator received unrecognised value for "
                "material_balance_type. This should not happen, "
                "please report this bug to the IDAES developers.".format(self.name)
            )

    def add_energy_splitting_constraints(self, mixed_block):
        """
        Creates constraints for splitting the energy flows.
        """
        if self.config.energy_split_basis == EnergySplittingType.none:
            # Do nothing and return
            return
        elif self.config.energy_split_basis == EnergySplittingType.equal_temperature:

            @self.Constraint(
                self.flowsheet().time,
                self.outlet_idx,
                doc="Temperature equality constraint",
            )
            def temperature_equality_eqn(b, t, o):
                o_block = getattr(self, o + "_state")
                return mixed_block[t].temperature == o_block[t].temperature

        elif self.config.energy_split_basis == EnergySplittingType.equal_molar_enthalpy:

            @self.Constraint(
                self.flowsheet().time,
                self.outlet_idx,
                doc="Molar enthalpy equality constraint",
            )
            def molar_enthalpy_equality_eqn(b, t, o):
                o_block = getattr(self, o + "_state")
                return mixed_block[t].enth_mol == o_block[t].enth_mol

        elif self.config.energy_split_basis == EnergySplittingType.enthalpy_split:
            # Validate split fraction type
            if (
                self.config.split_basis == SplittingType.phaseComponentFlow
                or self.config.split_basis == SplittingType.componentFlow
            ):
                raise ConfigurationError(
                    "{} Cannot use energy_split_basis == enthalpy_split "
                    "with split_basis == component or phaseComponent.".format(self.name)
                )

            def sf(t, o, p):
                if self.config.split_basis == SplittingType.totalFlow:
                    return self.split_fraction[t, o]
                elif self.config.split_basis == SplittingType.phaseFlow:
                    return self.split_fraction[t, o, p]

            @self.Constraint(
                self.flowsheet().time,
                self.outlet_idx,
                doc="Molar enthalpy splitting constraint",
            )
            def molar_enthalpy_splitting_eqn(b, t, o):
                o_block = getattr(self, o + "_state")
                return sum(
                    mixed_block[t].get_enthalpy_flow_terms(p) * sf(t, o, p)
                    for p in mixed_block.phase_list
                ) == sum(
                    o_block[t].get_enthalpy_flow_terms(p) for p in o_block.phase_list
                )

        else:
            raise BurntToast(
                "{} received unrecognised value for "
                "energy_split_basis. This should never happen, so"
                " please contact the IDAES developers with this "
                "bug.".format(self.name)
            )

    def add_momentum_splitting_constraints(self, mixed_block):
        """
        Creates constraints for splitting the momentum flows - done by equating
        pressures in outlets.
        """

        if self.config.momentum_balance_type is MomentumBalanceType.none:
            pass
        elif self.config.momentum_balance_type is MomentumBalanceType.pressureTotal:

            @self.Constraint(
                self.flowsheet().time,
                self.outlet_idx,
                doc="Pressure equality constraint",
            )
            def pressure_equality_eqn(b, t, o):
                o_block = getattr(self, o + "_state")
                return mixed_block[t].pressure == o_block[t].pressure

        else:
            raise NotImplementedError(
                f"Separators do not yet support momentum balances of type "
                f"{self.config.momentum_balance_type}"
            )

    def partition_outlet_flows(self, mb, outlet_list):
        """
        Creates outlet Port objects and tries to partition mixed stream flows
        between these

        Args:
            StateBlock representing the mixed flow to be split
            a list of names for outlets

        Returns:
            None
        """
        # Check arguments
        if self.config.construct_ports is False:
            raise ConfigurationError(
                "{} cannot have and ideal separator "
                "(ideal_separation = True) with "
                "construct_ports = False.".format(self.name)
            )
        if self.config.split_basis == SplittingType.totalFlow:
            raise ConfigurationError(
                "{} cannot do an ideal separation based "
                "on total flow. Either use ideal_separation = False or a "
                "different separation basis.".format(self.name)
            )
        if self.config.ideal_split_map is None:
            raise ConfigurationError(
                "{} was not provided with an "
                "ideal_split_map argument which is "
                "necessary for doing an ideal_separation.".format(self.name)
            )

        # Validate split map
        split_map = self.config.ideal_split_map
        idx_list = []
        if self.config.split_basis == SplittingType.phaseFlow:
            for p in mb.phase_list:
                idx_list.append((p))

            if len(idx_list) != len(split_map):
                raise ConfigurationError(
                    "{} ideal_split_map does not match with "
                    "split_basis chosen. ideal_split_map must"
                    " have a key for each combination of indices.".format(self.name)
                )
            for k in idx_list:
                if k not in split_map:
                    raise ConfigurationError(
                        "{} ideal_split_map does not match with "
                        "split_basis chosen. ideal_split_map must"
                        " have a key for each combination of indices.".format(self.name)
                    )

        elif self.config.split_basis == SplittingType.componentFlow:
            for j in mb.component_list:
                idx_list.append((j))

            if len(idx_list) != len(split_map):
                raise ConfigurationError(
                    "{} ideal_split_map does not match with "
                    "split_basis chosen. ideal_split_map must"
                    " have a key for each component.".format(self.name)
                )
        elif self.config.split_basis == SplittingType.phaseComponentFlow:
            for p in mb.phase_list:
                for j in mb.component_list:
                    idx_list.append((p, j))

            if len(idx_list) != len(split_map):
                raise ConfigurationError(
                    "{} ideal_split_map does not match with "
                    "split_basis chosen. ideal_split_map must"
                    " have a key for each phase-component pair.".format(self.name)
                )

        # Check that no. outlets matches split_basis
        if len(outlet_list) != len(idx_list):
            raise ConfigurationError(
                "{} Cannot perform ideal separation. Must have one "
                "outlet for each possible combination of the "
                "chosen split_basis.".format(self.name)
            )

        # Get units metadata
        units_meta = self.config.property_package.get_metadata()

        flow_basis = mb[self.flowsheet().time.first()].get_material_flow_basis()
        if flow_basis == MaterialFlowBasis.molar:
            flow_units = units_meta.get_derived_units("flow_mole")
        elif flow_basis == MaterialFlowBasis.mass:
            flow_units = units_meta.get_derived_units("flow_mass")
        else:
            # Let this pass for now with no units
            flow_units = None

        # Create tolerance Parameter for 0 flow outlets
        self.eps = Param(default=1e-8, mutable=True, units=flow_units)

        # Get list of port members
        s_vars = mb[self.flowsheet().time.first()].define_port_members()

        # Get phase component list(s)
        pc_set = mb.phase_component_set

        # Add empty Port objects
        for o in outlet_list:
            p_obj = Port(noruleinit=True, doc="Outlet Port")
            setattr(self, o, p_obj)

            # Iterate over members to create References or Expressions
            for var_name_port, var_obj in s_vars.items():
                # Get local variable name of component
                var_name_local = var_obj.local_name
                if var_name_local == "pressure" or var_name_local == "temperature":
                    # Assume outlets same as mixed flow - make Reference
                    e_obj = Reference(mb[:].component(var_name_local))

                elif var_name_local.startswith(
                    "mole_frac"
                ) or var_name_local.startswith("mass_frac"):
                    # Mole and mass frac need special handling
                    if "_phase" in var_name_local:
                        e_obj = VarLikeExpression(
                            self.flowsheet().time,
                            pc_set,
                            rule=partial(
                                _e_rule_mole_mass_frac_phase_comp,
                                mixed_block=mb,
                                var_obj=var_obj,
                                outlet=o,
                            ),
                        )

                    else:
                        e_obj = VarLikeExpression(
                            self.flowsheet().time,
                            mb.component_list,
                            rule=partial(
                                _e_rule_mole_mass_frac_comp,
                                mixed_block=mb,
                                var_obj=var_obj,
                                outlet=o,
                            ),
                        )

                elif var_name_local.endswith("_phase_comp"):
                    e_obj = VarLikeExpression(
                        self.flowsheet().time,
                        pc_set,
                        rule=partial(
                            _e_rule_phase_comp,
                            mixed_block=mb,
                            var_obj=var_obj,
                            outlet=o,
                        ),
                    )

                elif var_name_local.endswith("_phase"):
                    e_obj = VarLikeExpression(
                        self.flowsheet().time,
                        mb.phase_list,
                        rule=partial(
                            _e_rule_phase, mixed_block=mb, var_obj=var_obj, outlet=o
                        ),
                    )

                elif var_name_local.endswith("_comp"):
                    e_obj = VarLikeExpression(
                        self.flowsheet().time,
                        mb.component_list,
                        rule=partial(
                            _e_rule_comp, mixed_block=mb, var_obj=var_obj, outlet=o
                        ),
                    )

                else:
                    # Variable unindexed by phase or component
                    e_obj = VarLikeExpression(
                        self.flowsheet().time,
                        rule=partial(
                            _e_rule_other, mixed_block=mb, var_obj=var_obj, outlet=o
                        ),
                    )

                # Add Reference/Expression object to Separator model object
                setattr(self, "_" + o + "_" + var_name_local + "_ref", e_obj)

                # Add member to Port object
                p_obj.add(e_obj, var_name_port)

    def model_check(blk):
        """
        This method executes the model_check methods on the associated state
        blocks (if they exist). This method is generally called by a unit model
        as part of the unit's model_check method.

        Args:
            None

        Returns:
            None
        """
        # Try property block model check
        for t in blk.flowsheet().time:
            try:
                if blk.config.mixed_state_block is None:
                    blk.mixed_state[t].model_check()
                else:
                    blk.config.mixed_state_block.model_check()
            except AttributeError:
                _log.warning(
                    "{} Separator inlet state block has no "
                    "model check. To correct this, add a "
                    "model_check method to the associated "
                    "StateBlock class.".format(blk.name)
                )

            try:
                outlet_list = blk.create_outlet_list()
                for o in outlet_list:
                    o_block = getattr(blk, o + "_state")
                    o_block[t].model_check()
            except AttributeError:
                _log.warning(
                    "{} Separator outlet state block has no "
                    "model checks. To correct this, add a model_check"
                    " method to the associated StateBlock class.".format(blk.name)
                )

    def initialize_build(
        blk, outlvl=idaeslog.NOTSET, optarg=None, solver=None, hold_state=False
    ):
        """
        Initialization routine for separator

        Keyword Arguments:
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None, use
                     default solver options)
            solver : str indicating which solver to use during
                     initialization (default = None, use default solver)
            hold_state : flag indicating whether the initialization routine
                     should unfix any state variables fixed during
                     initialization, **default** - False. **Valid values:**
                     **True** - states variables are not unfixed, and a dict of
                     returned containing flags for which states were fixed
                     during initialization, **False** - state variables are
                     unfixed after initialization by calling the release_state
                     method.

        Returns:
            If hold_states is True, returns a dict containing flags for which
            states were fixed during initialization.
        """
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")

        # Create solver
        opt = get_solver(solver, optarg)

        # Initialize mixed state block
        if blk.config.mixed_state_block is not None:
            mblock = blk.config.mixed_state_block
        else:
            mblock = blk.mixed_state
        flags = mblock.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            hold_state=True,
        )

        # Solve for split fractions only
        component_status = {}
        for c in blk.component_objects((Block, Constraint)):
            for i in c:
                if not c[i].local_name == "sum_split_frac":
                    # Record current status of components to restore later
                    component_status[c[i]] = c[i].active
                    c[i].deactivate()

        if degrees_of_freedom(blk) != 0:
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                res = opt.solve(blk, tee=slc.tee)
                init_log.info(
                    "Initialization Step 1 Complete: {}".format(idaeslog.condition(res))
                )

        for c, s in component_status.items():
            if s:
                c.activate()

        if blk.config.ideal_separation:
            # If using ideal splitting, initialization should be complete
            return flags

        # Initialize outlet StateBlocks
        outlet_list = blk.create_outlet_list()

        # Premises for initializing outlet states:
        # 1. Intensive states remain unchanged - this is either a valid premise
        # or the actual state is impossible to calculate without solving the
        # full separator model.
        # 2. Extensive states are use split fractions if index matches, or
        # average of split fractions for outlet otherwise
        for o in outlet_list:
            # Get corresponding outlet StateBlock
            o_block = getattr(blk, o + "_state")

            # Create dict to store fixed status of state variables
            o_flags = {}
            for t in blk.flowsheet().time:

                # Calculate values for state variables
                s_vars = o_block[t].define_state_vars()
                for v in s_vars:
                    for k in s_vars[v]:
                        # Record whether variable was fixed or not
                        o_flags[t, v, k] = s_vars[v][k].fixed

                        # If fixed, use current value
                        # otherwise calculate guess from mixed state and fix
                        if not s_vars[v][k].fixed:
                            m_var = getattr(mblock[t], s_vars[v].local_name)
                            if "flow" in v:
                                # If a "flow" variable, is extensive
                                # Apply split fraction
                                if blk.config.split_basis == SplittingType.totalFlow:
                                    # All flows split by outlet
                                    s_vars[v][k].fix(
                                        value(m_var[k] * blk.split_fraction[(t, o)])
                                    )
                                elif "_phase_comp" in v:
                                    # Need to match indices, but use split frac
                                    if (
                                        blk.config.split_basis
                                        == SplittingType.phaseComponentFlow
                                    ):
                                        s_vars[v][k].fix(
                                            value(
                                                m_var[k]
                                                * blk.split_fraction[(t, o) + (k,)]
                                            )
                                        )
                                    elif (
                                        blk.config.split_basis
                                        == SplittingType.phaseFlow
                                    ):
                                        s_vars[v][k].fix(
                                            value(
                                                m_var[k]
                                                * blk.split_fraction[(t, o) + (k[0],)]
                                            )
                                        )
                                    elif (
                                        blk.config.split_basis
                                        == SplittingType.componentFlow
                                    ):
                                        s_vars[v][k].fix(
                                            value(
                                                m_var[k]
                                                * blk.split_fraction[(t, o) + (k[1],)]
                                            )
                                        )
                                    else:
                                        raise BurntToast(
                                            "{} encountered unrecognised "
                                            "SplittingType. This should not "
                                            "occur - please send this bug to "
                                            "the IDAES developers.".format(blk.name)
                                        )
                                elif "_phase" in v:
                                    if (
                                        blk.config.split_basis
                                        == SplittingType.phaseComponentFlow
                                    ):
                                        # Need average split fraction
                                        avg_split = value(
                                            sum(
                                                blk.split_fraction[t, o, k, j]
                                                for j in mblock.component_list
                                            )
                                            / len(mblock.component_list)
                                        )
                                        s_vars[v][k].fix(value(m_var[k] * avg_split))
                                    elif (
                                        blk.config.split_basis
                                        == SplittingType.phaseFlow
                                    ):
                                        s_vars[v][k].fix(
                                            value(
                                                m_var[k]
                                                * blk.split_fraction[(t, o) + (k,)]
                                            )
                                        )
                                    elif (
                                        blk.config.split_basis
                                        == SplittingType.componentFlow
                                    ):
                                        # Need average split fraction
                                        avg_split = value(
                                            sum(
                                                blk.split_fraction[t, o, j]
                                                for j in mblock.component_list
                                            )
                                            / len(mblock.component_list)
                                        )
                                        s_vars[v][k].fix(value(m_var[k] * avg_split))
                                    else:
                                        raise BurntToast(
                                            "{} encountered unrecognised "
                                            "SplittingType. This should not "
                                            "occur - please send this bug to "
                                            "the IDAES developers.".format(blk.name)
                                        )
                                elif "_comp" in v:
                                    if (
                                        blk.config.split_basis
                                        == SplittingType.phaseComponentFlow
                                    ):
                                        # Need average split fraction
                                        avg_split = value(
                                            sum(
                                                blk.split_fraction[t, o, p, k]
                                                for p in mblock.phase_list
                                            )
                                            / len(mblock.phase_list)
                                        )
                                        s_vars[v][k].fix(value(m_var[k] * avg_split))
                                    elif (
                                        blk.config.split_basis
                                        == SplittingType.phaseFlow
                                    ):
                                        # Need average split fraction
                                        avg_split = value(
                                            sum(
                                                blk.split_fraction[t, o, p]
                                                for p in mblock.phase_list
                                            )
                                            / len(mblock.phase_list)
                                        )
                                        s_vars[v][k].fix(value(m_var[k] * avg_split))
                                    elif (
                                        blk.config.split_basis
                                        == SplittingType.componentFlow
                                    ):
                                        s_vars[v][k].fix(
                                            value(
                                                m_var[k]
                                                * blk.split_fraction[(t, o) + (k,)]
                                            )
                                        )
                                    else:
                                        raise BurntToast(
                                            "{} encountered unrecognised "
                                            "SplittingType. This should not "
                                            "occur - please send this bug to "
                                            "the IDAES developers.".format(blk.name)
                                        )
                                else:
                                    # Assume unindexed extensive state
                                    # Need average split
                                    if (
                                        blk.config.split_basis
                                        == SplittingType.phaseComponentFlow
                                    ):
                                        # Need average split fraction
                                        avg_split = value(
                                            sum(
                                                blk.split_fraction[t, o, p, j]
                                                for (p, j) in mblock.phase_component_set
                                            )
                                            / len(mblock.phase_component_set)
                                        )
                                    elif (
                                        blk.config.split_basis
                                        == SplittingType.phaseFlow
                                    ):
                                        # Need average split fraction
                                        avg_split = value(
                                            sum(
                                                blk.split_fraction[t, o, p]
                                                for p in mblock.phase_list
                                            )
                                            / len(mblock.phase_list)
                                        )
                                    elif (
                                        blk.config.split_basis
                                        == SplittingType.componentFlow
                                    ):
                                        # Need average split fraction
                                        avg_split = value(
                                            sum(
                                                blk.split_fraction[t, o, j]
                                                for j in mblock.component_list
                                            )
                                            / len(mblock.component_list)
                                        )
                                    else:
                                        raise BurntToast(
                                            "{} encountered unrecognised "
                                            "SplittingType. This should not "
                                            "occur - please send this bug to "
                                            "the IDAES developers.".format(blk.name)
                                        )
                                    s_vars[v][k].fix(value(m_var[k] * avg_split))
                            else:
                                # Otherwise intensive, equate to mixed stream
                                s_vars[v][k].fix(m_var[k].value)

            # Call initialization routine for outlet StateBlock
            o_block.initialize(
                outlvl=outlvl,
                optarg=optarg,
                solver=solver,
                hold_state=False,
            )

            # Revert fixed status of variables to what they were before
            for t in blk.flowsheet().time:
                s_vars = o_block[t].define_state_vars()
                for v in s_vars:
                    for k in s_vars[v]:
                        s_vars[v][k].fixed = o_flags[t, v, k]

        if blk.config.mixed_state_block is None:
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                res = opt.solve(blk, tee=slc.tee)

            if not check_optimal_termination(res):
                raise InitializationError(
                    f"{blk.name} failed to initialize successfully. Please "
                    f"check the output logs for more information."
                )

            init_log.info(
                "Initialization Step 2 Complete: {}".format(idaeslog.condition(res))
            )
        else:
            init_log.info("Initialization Complete.")

        if hold_state is True:
            return flags
        else:
            blk.release_state(flags, outlvl=outlvl)

    def release_state(blk, flags, outlvl=idaeslog.NOTSET):
        """
        Method to release state variables fixed during initialization.

        Keyword Arguments:
            flags : dict containing information of which state variables
                    were fixed during initialization, and should now be
                    unfixed. This dict is returned by initialize if
                    hold_state = True.
            outlvl : sets output level of logging

        Returns:
            None
        """
        if blk.config.mixed_state_block is None:
            mblock = blk.mixed_state
        else:
            mblock = blk.config.mixed_state_block

        mblock.release_state(flags, outlvl=outlvl)

    def calculate_scaling_factors(self):
        mb_type = self.config.material_balance_type
        mixed_state = self.get_mixed_state_block()
        if mb_type == MaterialBalanceType.useDefault:
            t_ref = self.flowsheet().time.first()
            mb_type = mixed_state[t_ref].default_material_balance_type()
        super().calculate_scaling_factors()

        if hasattr(self, "temperature_equality_eqn"):
            for (t, i), c in self.temperature_equality_eqn.items():
                s = iscale.get_scaling_factor(
                    mixed_state[t].temperature, default=1, warning=True
                )
                iscale.constraint_scaling_transform(c, s)

        if hasattr(self, "pressure_equality_eqn"):
            for (t, i), c in self.pressure_equality_eqn.items():
                s = iscale.get_scaling_factor(
                    mixed_state[t].pressure, default=1, warning=True
                )
                iscale.constraint_scaling_transform(c, s)

        if hasattr(self, "material_splitting_eqn"):
            if mb_type == MaterialBalanceType.componentPhase:
                for (t, _, p, j), c in self.material_splitting_eqn.items():
                    flow_term = mixed_state[t].get_material_flow_terms(p, j)
                    s = iscale.get_scaling_factor(flow_term, default=1)
                    iscale.constraint_scaling_transform(c, s)
            elif mb_type == MaterialBalanceType.componentTotal:
                for (t, _, j), c in self.material_splitting_eqn.items():
                    s = None
                    for p in mixed_state.phase_list:
                        try:
                            ft = mixed_state[t].get_material_flow_terms(p, j)
                        except KeyError:
                            # This component does not exist in this phase
                            continue
                        if s is None:
                            s = iscale.get_scaling_factor(ft, default=1)
                        else:
                            _s = iscale.get_scaling_factor(ft, default=1)
                            s = _s if _s < s else s
                    iscale.constraint_scaling_transform(c, s)
            elif mb_type == MaterialBalanceType.total:
                pc_set = mixed_state.phase_component_set
                for (t, _), c in self.material_splitting_eqn.items():
                    for i, (p, j) in enumerate(pc_set):
                        ft = mixed_state[t].get_material_flow_terms(p, j)
                        if i == 0:
                            s = iscale.get_scaling_factor(ft, default=1)
                        else:
                            _s = iscale.get_scaling_factor(ft, default=1)
                            s = _s if _s < s else s
                    iscale.constraint_scaling_transform(c, s)

    def _get_performance_contents(self, time_point=0):
        if hasattr(self, "split_fraction"):
            var_dict = {}
            for k in self.split_fraction.keys():
                if k[0] == time_point:
                    var_dict[f"Split Fraction [{str(k[1:])}]"] = self.split_fraction[k]
            return {"vars": var_dict}
        else:
            return None

    def _get_stream_table_contents(self, time_point=0):
        outlet_list = self.create_outlet_list()

        if not self.config.ideal_separation:
            io_dict = {}
            if self.config.mixed_state_block is None:
                io_dict["Inlet"] = self.mixed_state
            else:
                io_dict["Inlet"] = self.config.mixed_state_block

            for o in outlet_list:
                io_dict[o] = getattr(self, o + "_state")

            return create_stream_table_dataframe(io_dict, time_point=time_point)

        else:
            stream_attributes = {}
            stream_attributes["Units"] = {}

            for n in ["inlet"] + outlet_list:
                port_obj = getattr(self, n)

                stream_attributes[n] = {}

                for k in port_obj.vars:
                    for i in port_obj.vars[k]:
                        if isinstance(i, float):
                            quant = report_quantity(port_obj.vars[k][time_point])
                            stream_attributes[n][k] = quant.m
                            stream_attributes["Units"][k] = quant.u
                        else:
                            if len(i) == 2:
                                kname = str(i[1])
                            else:
                                kname = str(i[1:])
                            quant = report_quantity(port_obj.vars[k][time_point, i[1:]])
                            stream_attributes[n][k + " " + kname] = quant.m
                            stream_attributes["Units"][k + " " + kname] = quant.u

            return DataFrame.from_dict(stream_attributes, orient="columns")


def _raise_split_indexed_fail_err(name, s):
    raise ConfigurationError(
        "{} Cannot use ideal splitting with "
        "this property package. Package uses "
        "indexed port member {} which cannot "
        "be partitioned. Please set "
        "configuration argument "
        "ideal_separation = False for this "
        "property package.".format(name, s)
    )


def _raise_split_unindexed_fail_err(name, s):
    raise ConfigurationError(
        "{} Cannot use ideal splitting with "
        "this property package. Package uses "
        "unindexed port member {} which does "
        "not have an equivalent indexed form. ".format(name, s)
    )


# Mole and mass frac need special handling
def _e_rule_mole_mass_frac_phase_comp(b, t, p, j, mixed_block, var_obj, outlet):
    pc_set = mixed_block.phase_component_set
    split_map = b.config.ideal_split_map

    if (p, j) in pc_set:
        if b.config.split_basis == SplittingType.phaseFlow:
            return var_obj[p, j]
        elif b.config.split_basis == SplittingType.componentFlow:
            if split_map[j] == outlet:
                return 1
            else:
                return b.eps
        elif b.config.split_basis == SplittingType.phaseComponentFlow:
            for ps in mixed_block.phase_list:
                if split_map[ps, j] == outlet:
                    return 1
            else:  # pylint: disable=W0120
                return b.eps
    else:
        raise BurntToast(
            "{} This should not happen. Please "
            "report this bug to the IDAES "
            "developers.".format(b.name)
        )


def _e_rule_mole_mass_frac_comp(b, t, j, mixed_block, var_obj, outlet):
    pc_set = mixed_block.phase_component_set
    l_name = var_obj.local_name
    split_map = b.config.ideal_split_map

    if b.config.split_basis == SplittingType.componentFlow:
        if split_map[j] == outlet:
            return 1
        # else:
        return b.eps

    elif b.config.split_basis == SplittingType.phaseComponentFlow:
        if any(split_map[phase, j] == outlet for phase in mixed_block.phase_list):
            return 1
        # else:
        return b.eps

    else:
        mfp = mixed_block[t].component(l_name.replace("_comp", "_phase_comp"))

        if mfp is None:
            raise ConfigurationError(
                "{} Cannot use ideal splitting with "
                "this property package. Package uses "
                "indexed port member {} which cannot "
                "be partitioned. Please set "
                "configuration argument "
                "ideal_separation = False for this "
                "property package.".format(b.name, l_name)
            )

        for phase in mixed_block.phase_list:
            if b.config.split_basis == SplittingType.phaseFlow:
                s_check = split_map[phase]
            else:
                raise BurntToast(
                    "{} This should not happen. Please"
                    " report this bug to the IDAES "
                    "developers.".format(b.name)
                )

            if s_check == outlet and (phase, j) in pc_set:
                return mfp[phase, j]
        # else:
        return b.eps


def _e_rule_phase_comp(b, t, p, j, mixed_block, var_obj, outlet):
    l_name = var_obj.local_name
    split_map = b.config.ideal_split_map

    if b.config.split_basis == SplittingType.phaseFlow:
        s_check = split_map[p]
    elif b.config.split_basis == SplittingType.componentFlow:
        s_check = split_map[j]
    elif b.config.split_basis == SplittingType.phaseComponentFlow:
        s_check = split_map[p, j]
    else:
        raise BurntToast(
            "{} This should not happen. Please"
            " report this bug to the IDAES "
            "developers.".format(b.name)
        )

    if s_check == outlet:
        return mixed_block[t].component(l_name)[p, j]
    else:
        return b.eps


def _e_rule_phase(b, t, p, mixed_block, var_obj, outlet):
    l_name = var_obj.local_name
    split_map = b.config.ideal_split_map

    if b.config.split_basis == SplittingType.phaseFlow:
        if split_map[p] == outlet:
            return mixed_block[t].component(l_name)[p]
        else:
            return b.eps

    else:
        mfp = mixed_block[t].component(l_name + "_comp")

        if mfp is None:
            _raise_split_indexed_fail_err(b.name, l_name)

        for j in mixed_block.component_list:
            if b.config.split_basis == SplittingType.componentFlow:
                s_check = split_map[j]
            elif b.config.split_basis == SplittingType.phaseComponentFlow:
                s_check = split_map[p, j]
            else:
                raise BurntToast(
                    "{} This should not happen. Please"
                    " report this bug to the IDAES "
                    "developers.".format(b.name)
                )

            if s_check == outlet:
                return mfp[p, j]
        # else:
        return b.eps


def _e_rule_comp(b, t, j, mixed_block, var_obj, outlet):
    pc_set = mixed_block.phase_component_set
    l_name = var_obj.local_name
    split_map = b.config.ideal_split_map

    if b.config.split_basis == SplittingType.componentFlow:
        if split_map[j] == outlet:
            return mixed_block[t].component(l_name)[j]
        else:
            return b.eps

    else:
        mfp = mixed_block[t].component("{0}_phase{1}".format(l_name[:-5], l_name[-5:]))

        if mfp is None:
            _raise_split_indexed_fail_err(b.name, l_name)

        for p in mixed_block.phase_list:
            if (p, j) in pc_set:
                if b.config.split_basis == SplittingType.phaseFlow:
                    s_check = split_map[p]
                elif b.config.split_basis == SplittingType.phaseComponentFlow:
                    s_check = split_map[p, j]
                else:
                    raise BurntToast(
                        "{} This should not happen. Please"
                        " report this bug to the IDAES "
                        "developers.".format(b.name)
                    )

                if s_check == outlet:
                    return mfp[p, j]
        return b.eps


def _e_rule_other(b, t, mixed_block, var_obj, outlet):
    pc_set = mixed_block.phase_component_set
    l_name = var_obj.local_name
    split_map = b.config.ideal_split_map

    if b.config.split_basis == SplittingType.phaseFlow:
        ivar = mixed_block[t].component(l_name + "_phase")
        if ivar is not None:
            for p in mixed_block.phase_list:
                if split_map[p] == outlet:
                    return ivar[p]
                else:
                    continue
        else:
            ivar = mixed_block[t].component(l_name + "_phase_comp")
            if ivar is not None:
                for p in mixed_block.phase_list:
                    if split_map[p] == outlet:
                        return sum(
                            ivar[p, j]
                            for j in mixed_block.component_list
                            if (p, j) in pc_set
                        )
                    else:
                        continue
            else:
                _raise_split_unindexed_fail_err(b.name, l_name)

    elif b.config.split_basis == SplittingType.componentFlow:
        ivar = mixed_block[t].component(l_name + "_comp")
        if ivar is not None:
            for j in mixed_block.component_list:
                if split_map[j] == outlet:
                    return ivar[j]
                else:
                    continue
        else:
            ivar = mixed_block[t].component(l_name + "_phase_comp")
            if ivar is not None:
                for j in mixed_block.component_list:
                    if split_map[j] == outlet:
                        return sum(
                            ivar[p, j]
                            for p in mixed_block.phase_list
                            if (p, j) in pc_set
                        )
                    else:
                        continue
            else:
                _raise_split_unindexed_fail_err(b.name, l_name)
    elif b.config.split_basis == SplittingType.phaseComponentFlow:
        ivar = mixed_block[t].component(l_name + "_phase_comp")
        if ivar is not None:
            for p in mixed_block.phase_list:
                for j in mixed_block.component_list:
                    if split_map[p, j] == outlet and (p, j) in pc_set:
                        return ivar[p, j]
                    else:
                        continue
        else:
            _raise_split_unindexed_fail_err(b.name, l_name)
    else:
        # Unrecognised split tupe
        raise BurntToast(
            "{} received unrecognised value for "
            "split_basis argument. This should never "
            "happen, so please contact the IDAES "
            "developers with this bug.".format(b.name)
        )
