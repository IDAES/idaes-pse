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
Base class for unit models
"""

from pyomo.environ import check_optimal_termination
from pyomo.common.config import ConfigValue

from idaes.core.base.process_base import (
    declare_process_block_class,
    ProcessBlockData,
    useDefault,
)
from idaes.core.base.property_base import StateBlock
from idaes.core.base.control_volume_base import (
    ControlVolumeBlockData,
    FlowDirection,
    MaterialBalanceType,
)
from idaes.core.util.exceptions import (
    BurntToast,
    ConfigurationError,
    BalanceTypeNotSupportedError,
    InitializationError,
)
from idaes.core.util.tables import create_stream_table_dataframe
import idaes.logger as idaeslog
from idaes.core.solvers import get_solver
from idaes.core.util.config import DefaultBool


__author__ = "John Eslick, Qi Chen, Andrew Lee"


__all__ = ["UnitModelBlockData", "UnitModelBlock"]

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("UnitModelBlock")
class UnitModelBlockData(ProcessBlockData):
    """
    This is the class for process unit operations models. These are models that
    would generally appear in a process flowsheet or superstructure.
    """

    # Create Class ConfigBlock
    CONFIG = ProcessBlockData.CONFIG()
    CONFIG.declare(
        "dynamic",
        ConfigValue(
            default=useDefault,
            domain=DefaultBool,
            description="Dynamic model flag",
            doc="""Indicates whether this model will be dynamic or not,
**default** = useDefault.
**Valid values:** {
**useDefault** - get flag from parent (default = False),
**True** - set as a dynamic model,
**False** - set as a steady-state model.}""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=useDefault,
            domain=DefaultBool,
            description="Holdup construction flag",
            doc="""Indicates whether holdup terms should be constructed or not.
Must be True if dynamic = True,
**default** - False.
**Valid values:** {
**useDefault** - get flag from parent (default = False),
**True** - construct holdup terms,
**False** - do not construct holdup terms}""",
        ),
    )

    def build(self):
        """
        General build method for UnitModelBlockData. This method calls a number
        of sub-methods which automate the construction of expected attributes
        of unit models.

        Inheriting models should call `super().build`.

        Args:
            None

        Returns:
            None
        """
        super(UnitModelBlockData, self).build()

        # Add a placeholder for initialization order
        self._initialization_order = []

        # Set up dynamic flag and time domain
        self._setup_dynamics()

    def model_check(blk):
        """
        This is a general purpose initialization routine for simple unit
        models. This method assumes a single ControlVolume block called
        controlVolume and tries to call the model_check method of the
        controlVolume block. If an AttributeError is raised, the check is
        passed.

        More complex models should overload this method with a model_check
        suited to the particular application, especially if there are multiple
        ControlVolume blocks present.

        Args:
            None

        Returns:
            None
        """
        # Run control volume block model checks
        try:
            blk.controlVolume.model_check()
        except AttributeError:
            pass

    def add_port(self, name, block, doc=None):
        """
        This is a method to build Port objects in a unit model and
        connect these to a specified StateBlock.

        Keyword Args:
            name : name to use for Port object.
            block : an instance of a StateBlock to use as the source to
                    populate the Port object
            doc : doc string for Port object, default = None.

        Returns:
            A Pyomo Port object and associated components.
        """
        # Create Port
        try:
            port, ref_name_list = block.build_port(doc)
        except AttributeError:
            raise ConfigurationError(
                f"{self.name} block object provided to add_port method is not an "
                f"instance of a StateBlock object (does not have a build_port method)."
            )

        # Add Port and References to unit moodel
        self.add_component(name, port)
        for ref, cname in ref_name_list:
            ref_name = block.get_port_reference_name(cname, name)
            self.add_component(ref_name, ref)

        return port

    def add_inlet_port(self, name=None, block=None, doc=None):
        """
        This is a method to build inlet Port objects in a unit model and
        connect these to a specified control volume or state block.

        The name and block arguments are optional, but must be used together.
        i.e. either both arguments are provided or neither.

        Keyword Args:
            name : name to use for Port object (default = "inlet").
            block : an instance of a ControlVolume or StateBlock to use as the
                    source to populate the Port object. If a ControlVolume is
                    provided, the method will use the inlet state block as
                    defined by the ControlVolume. If not provided, method will
                    attempt to default to an object named control_volume.
            doc : doc string for Port object (default = "Inlet Port")

        Returns:
            A Pyomo Port object and associated components.
        """
        if block is None:
            # Check that name is None
            if name is not None:
                raise ConfigurationError(
                    "{} add_inlet_port was called without a block argument"
                    " but a name argument was provided. Either both "
                    "a name and a block must be provided or neither.".format(self.name)
                )
            else:
                name = "inlet"
            # Try for default ControlVolume name
            try:
                block = self.control_volume
            except AttributeError:
                raise ConfigurationError(
                    "{} add_inlet_port was called without a block argument"
                    " but no default ControlVolume exists "
                    "(control_volume). Please provide block to which the "
                    "Port should be associated.".format(self.name)
                )
        else:
            # Check that name is not None
            if name is None:
                raise ConfigurationError(
                    "{} add_inlet_port was called with a block argument, "
                    "but a name argument was not provided. Either both "
                    "a name and a block must be provided or neither.".format(self.name)
                )

        if doc is None:
            doc = "Inlet Port"

        # Determine type of sourcce block
        if isinstance(block, ControlVolumeBlockData):
            # Work out if this is a 0D or 1D block
            try:
                # Try 0D first
                sblock = block.properties_in
                port, ref_name_list = sblock.build_port(doc)
            except AttributeError:
                # Otherwise a 1D control volume
                try:
                    sblock = block.properties

                    # Need to determine correct subset of indices to add to Port
                    if block._flow_direction == FlowDirection.forward:
                        _idx = block.length_domain.first()
                    elif block._flow_direction == FlowDirection.backward:
                        _idx = block.length_domain.last()

                    port, ref_name_list = sblock.build_port(
                        doc, slice_index=(slice(None), _idx)
                    )

                except AttributeError:
                    raise ConfigurationError(
                        f"{self.name} - control volume does not have expected "
                        f"names for StateBlocks. Please check that the control "
                        f"volume was constructed correctly."
                    )
        else:
            # Assume a StateBlock indexed only by time
            sblock = block
            port, ref_name_list = sblock.build_port(doc)

        # Add Port and References to unit moodel
        self.add_component(name, port)
        for ref, cname in ref_name_list:
            ref_name = sblock.get_port_reference_name(cname, name)
            self.add_component(ref_name, ref)

        return port

    def add_outlet_port(self, name=None, block=None, doc=None):
        """
        This is a method to build outlet Port objects in a unit model and
        connect these to a specified control volume or state block.

        The name and block arguments are optional, but must be used together.
        i.e. either both arguments are provided or neither.

        Keyword Args:
            name : name to use for Port object (default = "outlet").
            block : an instance of a ControlVolume or StateBlock to use as the
                    source to populate the Port object. If a ControlVolume is
                    provided, the method will use the outlet state block as
                    defined by the ControlVolume. If not provided, method will
                    attempt to default to an object named control_volume.
            doc : doc string for Port object (default = "Outlet Port")

        Returns:
            A Pyomo Port object and associated components.
        """
        if block is None:
            # Check that name is None
            if name is not None:
                raise ConfigurationError(
                    "{} add_outlet_port was called without a block "
                    "argument but a name argument was provided. Either "
                    "both a name and a block must be provided or neither.".format(
                        self.name
                    )
                )
            else:
                name = "outlet"
            # Try for default ControlVolume name
            try:
                block = self.control_volume
            except AttributeError:
                raise ConfigurationError(
                    "{} add_outlet_port was called without a block "
                    "argument but no default ControlVolume exists "
                    "(control_volume). Please provide block to which the "
                    "Port should be associated.".format(self.name)
                )
        else:
            # Check that name is not None
            if name is None:
                raise ConfigurationError(
                    "{} add_outlet_port was called with a block argument, "
                    "but a name argument was not provided. Either both "
                    "a name and a block must be provided or neither.".format(self.name)
                )

        if doc is None:
            doc = "Outlet Port"

        # Determine type of sourcce block
        if isinstance(block, ControlVolumeBlockData):
            # Work out if this is a 0D or 1D block
            try:
                # Try 0D first
                sblock = block.properties_out
                port, ref_name_list = sblock.build_port(doc)
            except AttributeError:
                # Otherwise a 1D control volume
                try:
                    sblock = block.properties

                    # Need to determine correct subset of indices to add to Port
                    if block._flow_direction == FlowDirection.backward:
                        _idx = block.length_domain.first()
                    elif block._flow_direction == FlowDirection.forward:
                        _idx = block.length_domain.last()

                    port, ref_name_list = sblock.build_port(
                        doc, slice_index=(slice(None), _idx)
                    )

                except AttributeError:
                    raise ConfigurationError(
                        f"{self.name} - control volume does not have expected "
                        f"names for StateBlocks. Please check that the control "
                        f"volume was constructed correctly."
                    )
        else:
            # Assume a StateBlock indexed only by time
            sblock = block
            port, ref_name_list = sblock.build_port(doc)

        # Add Port and References to unit moodel
        self.add_component(name, port)
        for ref, cname in ref_name_list:
            ref_name = sblock.get_port_reference_name(cname, name)
            self.add_component(ref_name, ref)

        return port

    def add_state_material_balances(self, balance_type, state_1, state_2):
        """
        Method to add material balances linking two State Blocks in a Unit
        Model. This method is not intended to replace Control Volumes, but
        to automate writing material balances linking isolated State Blocks
        in those models where this is required.

        Args:
            balance_type - a MaterialBalanceType Enum indicating the type
                            of material balances to write
            state_1 - first State Block to be linked by balances
            state_2 - second State Block to be linked by balances

        Returns:
            None
        """
        # Confirm that both state blocks stem from the same parameter block
        if not isinstance(state_1, StateBlock):
            raise ConfigurationError(
                "{} state_1 argument to add_state_material_balances "
                "was not an instance of a State Block.".format(self.name)
            )

        if not isinstance(state_2, StateBlock):
            raise ConfigurationError(
                "{} state_2 argument to add_state_material_balances "
                "was not an instance of a State Block.".format(self.name)
            )

        # Check that no constraint with the same name exists
        # We will only support using this method once per Block
        if hasattr(self, "state_material_balances"):
            raise AttributeError(
                "{} a set of constraints named state_material_balances "
                "already exists in the current UnitModel. To avoid "
                "confusion, add_state_material_balances is only supported "
                "once per UnitModel.".format(self.name)
            )

        # Get a representative time point for testing
        rep_time = self.flowsheet().time.first()
        if state_1[rep_time].params is not state_2[rep_time].params:
            raise ConfigurationError(
                "{} add_state_material_balances method was provided with "
                "State Blocks are not linked to the same "
                "instance of a Physical Parameter Block. This method "
                "only supports linking State Blocks from the same "
                "Physical Parameter Block.".format(self.name)
            )

        if balance_type == MaterialBalanceType.useDefault:
            balance_type = state_1[rep_time].default_material_balance_type()

        phase_list = state_1.phase_list
        component_list = state_1.component_list
        pc_set = state_1.phase_component_set

        if balance_type == MaterialBalanceType.componentPhase:
            # TODO : Should we include an optional phase equilibrium term here
            # to allow for systems where a phase-transition may occur?

            @self.Constraint(
                self.flowsheet().time,
                pc_set,
                doc="State material balances",
            )
            def state_material_balances(b, t, p, j):
                return state_1[t].get_material_flow_terms(p, j) == state_2[
                    t
                ].get_material_flow_terms(p, j)

        elif balance_type == MaterialBalanceType.componentTotal:

            @self.Constraint(
                self.flowsheet().time,
                component_list,
                doc="State material balances",
            )
            def state_material_balances(b, t, j):
                return sum(
                    state_1[t].get_material_flow_terms(p, j)
                    for p in phase_list
                    if (p, j) in pc_set
                ) == sum(
                    state_2[t].get_material_flow_terms(p, j)
                    for p in phase_list
                    if (p, j) in pc_set
                )

        elif balance_type == MaterialBalanceType.total:

            @self.Constraint(
                self.flowsheet().time,
                doc="State material balances",
            )
            def state_material_balances(b, t):
                return sum(
                    state_1[t].get_material_flow_terms(p, j) for p, j in pc_set
                ) == sum(state_2[t].get_material_flow_terms(p, j) for p, j in pc_set)

        elif balance_type == MaterialBalanceType.elementTotal:
            raise BalanceTypeNotSupportedError(
                "{} add_state_material_balances does not support "
                "MaterialBalanceType.elementTotal.".format(self.name)
            )
        elif balance_type == MaterialBalanceType.none:
            raise BalanceTypeNotSupportedError(
                "{} add_state_material_balances does not support "
                "MaterialBalanceType.None.".format(self.name)
            )
        else:
            raise BurntToast(
                "{} add_state_material_balances received an unexpected "
                "argument for balance_type. This should never happen. Please "
                "contact the IDAES developers with this bug.".format(self.name)
            )

    def _get_stream_table_contents(self, time_point=0):
        """
        Assume unit has standard configuration of 1 inlet and 1 outlet.

        Developers should overload this as appropriate.
        """
        try:
            return create_stream_table_dataframe(
                {"Inlet": self.inlet, "Outlet": self.outlet}, time_point=time_point
            )
        except AttributeError:
            raise ConfigurationError(
                f"Unit model {self.name} does not have the standard Port "
                f"names (inet and outlet). Please contact the unit model "
                f"developer to develop a unit specific stream table."
            )

    def initialize(blk, *args, **kwargs):
        """
        Initialization routine for Unit Model objects and associated
        components.

        This method is intended for initializing IDAES unit models and any
        modeling components attached to them, such as costing blocks. This
        method iterates through all objects in blk._initialization_order and
        deactivates them, followed by calling blk.initialize_build. Finally,
        it iterates through all objects in blk._initialization_order in reverse
        and re-activates these whilst calling the associated initialize method.

        Currently, parsing of arguments to the initialize method of attached
        blocks is hard coded - this will be addressed in a future PR.
        Currently, the method supports the following attached components:

        * UnitModelCostingBlocks

        Args:
            costing_args - dict arguments to be passed to costing block
                           initialize method

        For other arguments, see the initilize_unit method.
        """
        # Get any arguments for costing if provided
        cost_args = kwargs.pop("costing_args", {})

        # Get the costing block if present
        # TODO: Clean up in IDAES v2.0
        init_order = blk._initialization_order
        if hasattr(blk, "costing") and blk.costing not in blk._initialization_order:
            # Fallback for older style costing
            init_order.append(blk.costing)

        # If costing block exists, deactivate
        for c in init_order:
            c.deactivate()

        # Remember to collect flags for fixed vars
        flags = blk.initialize_build(*args, **kwargs)

        # If costing block exists, activate and initialize
        for c in init_order:
            c.activate()

            if hasattr(c, "initialize"):
                c.initialize(**cost_args)

        # Return any flags returned by initialize_build
        return flags

    def initialize_build(
        blk, state_args=None, outlvl=idaeslog.NOTSET, solver=None, optarg=None
    ):
        """
        This is a general purpose initialization routine for simple unit
        models. This method assumes a single ControlVolume block called
        controlVolume, and first initializes this and then attempts to solve
        the entire unit.

        More complex models should overload this method with their own
        initialization routines,

        Keyword Arguments:
            state_args : a dict of arguments to be passed to the property
                           package(s) to provide an initial state for
                           initialization (see documentation of the specific
                           property package) (default = {}).
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None, use
                     default solver options)
            solver : str indicating which solver to use during
                     initialization (default = None, use default IDAES solver)

        Returns:
            None
        """
        if optarg is None:
            optarg = {}

        # Set solver options
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")

        opt = get_solver(solver, optarg)

        # ---------------------------------------------------------------------
        # Initialize control volume block
        flags = blk.control_volume.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
        )

        init_log.info_high("Initialization Step 1 Complete.")

        # ---------------------------------------------------------------------
        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            results = opt.solve(blk, tee=slc.tee)

        init_log.info_high(
            "Initialization Step 2 {}.".format(idaeslog.condition(results))
        )

        # ---------------------------------------------------------------------
        # Release Inlet state
        blk.control_volume.release_state(flags, outlvl)

        if not check_optimal_termination(results):
            raise InitializationError(
                f"{blk.name} failed to initialize successfully. Please check "
                f"the output logs for more information."
            )

        init_log.info("Initialization Complete: {}".format(idaeslog.condition(results)))

        return None

    def del_component(self, name_or_object):
        """
        Delete a component from this block. Need to introduce code to handle
        un-registering costing blocks
        """
        obj = self.component(name_or_object)

        # TODO: See if Pyomo can give us a call-back to do this
        try:
            # If this is a costing block, need to unregister it
            obj.del_costing()
        except AttributeError:
            pass

        super().del_component(obj)
