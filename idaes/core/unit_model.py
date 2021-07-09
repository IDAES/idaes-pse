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

from pyomo.environ import Reference
from pyomo.network import Port
from pyomo.common.config import ConfigValue, In

from .process_base import (declare_process_block_class,
                           ProcessBlockData,
                           useDefault)
from .property_base import StateBlock
from .control_volume_base import (ControlVolumeBlockData,
                                  FlowDirection,
                                  MaterialBalanceType)
from idaes.core.util.exceptions import (BurntToast,
                                        ConfigurationError,
                                        PropertyPackageError,
                                        BalanceTypeNotSupportedError)
from idaes.core.util.tables import create_stream_table_dataframe
import idaes.core.util.unit_costing
import idaes.logger as idaeslog
from idaes.core.util import get_solver

__author__ = "John Eslick, Qi Chen, Andrew Lee"


__all__ = ['UnitModelBlockData', 'UnitModelBlock']

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
    CONFIG.declare("dynamic", ConfigValue(
        default=useDefault,
        domain=In([useDefault, True, False]),
        description="Dynamic model flag",
        doc="""Indicates whether this model will be dynamic or not,
**default** = useDefault.
**Valid values:** {
**useDefault** - get flag from parent (default = False),
**True** - set as a dynamic model,
**False** - set as a steady-state model.}"""))
    CONFIG.declare("has_holdup", ConfigValue(
        default=useDefault,
        domain=In([useDefault, True, False]),
        description="Holdup construction flag",
        doc="""Indicates whether holdup terms should be constructed or not.
Must be True if dynamic = True,
**default** - False.
**Valid values:** {
**useDefault** - get flag from parent (default = False),
**True** - construct holdup terms,
**False** - do not construct holdup terms}"""))

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

    def add_port(blk, name=None, block=None, doc=None):
        """
        This is a method to build Port objects in a unit model and
        connect these to a specified StateBlock.

        Keyword Args:
            name : name to use for Port object.
            block : an instance of a StateBlock to use as the source to
                    populate the Port object
            doc : doc string for Port object

        Returns:
            A Pyomo Port object and associated components.
        """
        # Validate block object
        if not isinstance(block, StateBlock):
            raise ConfigurationError("{} block object provided to add_port "
                                     "method is not an instance of a "
                                     "StateBlock object. IDAES port objects "
                                     "should only be associated with "
                                     "StateBlocks.".format(blk.name))

        # Create empty Port
        p = Port(noruleinit=True, doc=doc)
        setattr(blk, name, p)

        p._state_block = (block, )

        # Get dict of Port members and names
        member_list = block[
                blk.flowsheet().time.first()].define_port_members()

        # Create References for port members
        for s in member_list:
            if not member_list[s].is_indexed():
                slicer = block[:].component(member_list[s].local_name)
            else:
                slicer = block[:].component(member_list[s].local_name)[...]

            r = Reference(slicer)
            setattr(blk, "_"+s+"_"+name+"_ref", r)

            # Add Reference to Port
            p.add(r, s)

        return p

    def add_inlet_port(blk, name=None, block=None, doc=None):
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
                        "a name and a block must be provided or neither."
                        .format(blk.name))
            else:
                name = "inlet"
            # Try for default ControlVolume name
            try:
                block = blk.control_volume
            except AttributeError:
                raise ConfigurationError(
                        "{} add_inlet_port was called without a block argument"
                        " but no default ControlVolume exists "
                        "(control_volume). Please provide block to which the "
                        "Port should be associated.".format(blk.name))
        else:
            # Check that name is not None
            if name is None:
                raise ConfigurationError(
                        "{} add_inlet_port was called with a block argument, "
                        "but a name argument was not provided. Either both "
                        "a name and a block must be provided or neither."
                        .format(blk.name))

        if doc is None:
            doc = "Inlet Port"

        # Create empty Port
        p = Port(noruleinit=True, doc=doc)
        setattr(blk, name, p)

        # Get dict of Port members and names
        if isinstance(block, ControlVolumeBlockData):
            try:
                member_list = (block.properties_in[
                                    block.flowsheet().time.first()]
                               .define_port_members())
                p._state_block = (block.properties_in, )
            except AttributeError:
                try:
                    member_list = (block.properties[
                                    block.flowsheet().time.first(), 0]
                                   .define_port_members())
                    if block._flow_direction == FlowDirection.forward:
                        p._state_block = (block.properties,
                                          block.length_domain.first())
                    elif block._flow_direction == FlowDirection.backward:
                        p._state_block = (block.properties,
                                          block.length_domain.last())
                except AttributeError:
                    raise PropertyPackageError(
                            "{} property package does not appear to have "
                            "implemented a define_port_memebers method. "
                            "Please contact the developer of the property "
                            "package.".format(blk.name))
        elif isinstance(block, StateBlock):
            member_list = block[
                    blk.flowsheet().time.first()].define_port_members()
            p._state_block = (block, )
        else:
            raise ConfigurationError(
                    "{} block provided to add_inlet_port "
                    "method was not an instance of a "
                    "ControlVolume or a StateBlock."
                    .format(blk.name))

        # Create References for port members
        for s in member_list:
            if not member_list[s].is_indexed():
                if isinstance(block, ControlVolumeBlockData):
                    try:
                        slicer = block.properties_in[:].component(
                                        member_list[s].local_name)
                    except AttributeError:
                        if block._flow_direction == FlowDirection.forward:
                            _idx = block.length_domain.first()
                        elif block._flow_direction == FlowDirection.backward:
                            _idx = block.length_domain.last()
                        else:
                            raise BurntToast(
                                    "{} flow_direction argument received "
                                    "invalid value. This should never "
                                    "happen, so please contact the IDAES "
                                    "developers with this bug."
                                    .format(blk.name))
                        slicer = (block.properties[:, _idx]
                                  .component(member_list[s].local_name))
                elif isinstance(block, StateBlock):
                    slicer = block[:].component(member_list[s].local_name)
                else:
                    raise ConfigurationError(
                            "{} block provided to add_inlet_port "
                            "method was not an instance of a "
                            "ControlVolume or a StateBlock."
                            .format(blk.name))
            else:
                if isinstance(block, ControlVolumeBlockData):
                    try:
                        slicer = block.properties_in[:].component(
                                        member_list[s].local_name)[...]
                    except AttributeError:
                        if block._flow_direction == FlowDirection.forward:
                            _idx = block.length_domain.first()
                        elif block._flow_direction == FlowDirection.backward:
                            _idx = block.length_domain.last()
                        else:
                            raise BurntToast(
                                    "{} flow_direction argument received "
                                    "invalid value. This should never "
                                    "happen, so please contact the IDAES "
                                    "developers with this bug."
                                    .format(blk.name))
                        slicer = (block.properties[:, _idx].component(
                                    member_list[s].local_name))[...]
                elif isinstance(block, StateBlock):
                    slicer = block[:].component(member_list[s].local_name)[...]
                else:
                    raise ConfigurationError(
                            "{} block provided to add_inlet_port "
                            "method was not an instance of a "
                            "ControlVolume or a StateBlock."
                            .format(blk.name))

            r = Reference(slicer)
            setattr(blk, "_"+s+"_"+name+"_ref", r)

            # Add Reference to Port
            p.add(r, s)

        return p

    def add_outlet_port(blk, name=None, block=None, doc=None):
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
                        "argument  but a name argument was provided. Either "
                        "both a name and a block must be provided or neither."
                        .format(blk.name))
            else:
                name = "outlet"
            # Try for default ControlVolume name
            try:
                block = blk.control_volume
            except AttributeError:
                raise ConfigurationError(
                        "{} add_outlet_port was called without a block "
                        "argument but no default ControlVolume exists "
                        "(control_volume). Please provide block to which the "
                        "Port should be associated.".format(blk.name))
        else:
            # Check that name is not None
            if name is None:
                raise ConfigurationError(
                        "{} add_outlet_port was called with a block argument, "
                        "but a name argument was not provided. Either both "
                        "a name and a block must be provided or neither."
                        .format(blk.name))

        if doc is None:
            doc = "Outlet Port"

        # Create empty Port
        p = Port(noruleinit=True, doc=doc)
        setattr(blk, name, p)

        # Get dict of Port members and names
        if isinstance(block, ControlVolumeBlockData):
            try:
                member_list = (block.properties_out[
                                    block.flowsheet().time.first()]
                               .define_port_members())
                p._state_block = (block.properties_out, )
            except AttributeError:
                try:
                    member_list = (block.properties[
                                    block.flowsheet().time.first(), 0]
                                   .define_port_members())
                    if block._flow_direction == FlowDirection.forward:
                        p._state_block = (block.properties,
                                          block.length_domain.last())
                    elif block._flow_direction == FlowDirection.backward:
                        p._state_block = (block.properties,
                                          block.length_domain.first())
                except AttributeError:
                    raise PropertyPackageError(
                            "{} property package does not appear to have "
                            "implemented a define_port_members method. "
                            "Please contact the developer of the property "
                            "package.".format(blk.name))
        elif isinstance(block, StateBlock):
            member_list = block[
                    blk.flowsheet().time.first()].define_port_members()
            p._state_block = (block, )
        else:
            raise ConfigurationError(
                    "{} block provided to add_inlet_port "
                    "method was not an instance of a "
                    "ControlVolume or a StateBlock."
                    .format(blk.name))

        # Create References for port members
        for s in member_list:
            if not member_list[s].is_indexed():
                if isinstance(block, ControlVolumeBlockData):
                    try:
                        slicer = block.properties_out[:].component(
                                        member_list[s].local_name)
                    except AttributeError:
                        if block._flow_direction == FlowDirection.forward:
                            _idx = block.length_domain.last()
                        elif block._flow_direction == FlowDirection.backward:
                            _idx = block.length_domain.first()
                        else:
                            raise BurntToast(
                                    "{} flow_direction argument received "
                                    "invalid value. This should never "
                                    "happen, so please contact the IDAES "
                                    "developers with this bug."
                                    .format(blk.name))
                        slicer = (block.properties[:, _idx]
                                  .component(member_list[s].local_name))
                elif isinstance(block, StateBlock):
                    slicer = block[:].component(member_list[s].local_name)
                else:
                    raise ConfigurationError(
                            "{} block provided to add_inlet_port "
                            "method was not an instance of a "
                            "ControlVolume or a StateBlock."
                            .format(blk.name))
            else:
                # Need to use slice notation on indexed comenent as well
                if isinstance(block, ControlVolumeBlockData):
                    try:
                        slicer = block.properties_out[:].component(
                                        member_list[s].local_name)[...]
                    except AttributeError:
                        if block._flow_direction == FlowDirection.forward:
                            _idx = block.length_domain.last()
                        elif block._flow_direction == FlowDirection.backward:
                            _idx = block.length_domain.first()
                        else:
                            raise BurntToast(
                                    "{} flow_direction argument received "
                                    "invalid value. This should never "
                                    "happen, so please contact the IDAES "
                                    "developers with this bug."
                                    .format(blk.name))
                        slicer = (block.properties[:, _idx].component(
                                    member_list[s].local_name))[...]
                elif isinstance(block, StateBlock):
                    slicer = block[:].component(member_list[s].local_name)[...]
                else:
                    raise ConfigurationError(
                            "{} block provided to add_inlet_port "
                            "method was not an instance of a "
                            "ControlVolume or a StateBlock."
                            .format(blk.name))

            r = Reference(slicer)
            setattr(blk, "_"+s+"_"+name+"_ref", r)

            # Add Reference to Port
            p.add(r, s)

        return p

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
                    "was not an instance of a State Block.".format(self.name))

        if not isinstance(state_2, StateBlock):
            raise ConfigurationError(
                    "{} state_2 argument to add_state_material_balances "
                    "was not an instance of a State Block.".format(self.name))

        # Check that no constraint with the same name exists
        # We will only support using this method once per Block
        if hasattr(self, "state_material_balances"):
            raise AttributeError(
                    "{} a set of constraints named state_material_balances "
                    "already exists in the current UnitModel. To avoid "
                    "confusion, add_state_material_balances is only supported "
                    "once per UnitModel.".format(self.name))

        # Get a representative time point for testing
        rep_time = self.flowsheet().time.first()
        if state_1[rep_time].params is not state_2[rep_time].params:
            raise ConfigurationError(
                    "{} add_state_material_balances method was provided with "
                    "State Blocks are not linked to the same "
                    "instance of a Physical Parameter Block. This method "
                    "only supports linking State Blocks from the same "
                    "Physical Parameter Block.".format(self.name))

        if balance_type == MaterialBalanceType.useDefault:
            balance_type = (
                state_1[rep_time].default_material_balance_type()
            )

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
                return state_1[t].get_material_flow_terms(
                    p, j
                ) == state_2[t].get_material_flow_terms(p, j)

        elif balance_type == MaterialBalanceType.componentTotal:

            @self.Constraint(
                self.flowsheet().time,
                component_list,
                doc="State material balances",
            )
            def state_material_balances(b, t, j):
                return sum(
                    state_1[t].get_material_flow_terms(p, j)
                    for p in phase_list if (p, j) in pc_set
                ) == sum(
                    state_2[t].get_material_flow_terms(p, j)
                    for p in phase_list if (p, j) in pc_set
                )

        elif balance_type == MaterialBalanceType.total:

            @self.Constraint(
                self.flowsheet().time,
                doc="State material balances",
            )
            def state_material_balances(b, t):
                return sum(
                    state_1[t].get_material_flow_terms(p, j) for p, j in pc_set
                ) == sum(
                    state_2[t].get_material_flow_terms(p, j) for p, j in pc_set
                )

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
            return create_stream_table_dataframe({"Inlet": self.inlet,
                                                  "Outlet": self.outlet},
                                                 time_point=time_point)
        except AttributeError:
            raise ConfigurationError(
                    f"Unit model {self.name} does not have the standard Port "
                    f"names (inet and outlet). Please contact the unit model "
                    f"developer to develop a unit specific stream table.")

    def initialize(blk, state_args=None, outlvl=idaeslog.NOTSET,
                   solver=None, optarg=None):
        '''
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
        '''
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

        init_log.info_high('Initialization Step 1 Complete.')

        # ---------------------------------------------------------------------
        # Solve unit

        # if costing block exists, deactivate
        if hasattr(blk, "costing"):
            blk.costing.deactivate()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            results = opt.solve(blk, tee=slc.tee)

        init_log.info_high(
            "Initialization Step 2 {}.".format(idaeslog.condition(results))
        )

        # if costing block exists, activate and initialize
        if hasattr(blk, "costing"):
            blk.costing.activate()
            idaes.core.util.unit_costing.initialize(blk.costing)
        # ---------------------------------------------------------------------
        # Release Inlet state
        blk.control_volume.release_state(flags, outlvl)

        init_log.info('Initialization Complete: {}'
                      .format(idaeslog.condition(results)))
