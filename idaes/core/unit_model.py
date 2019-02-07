##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
"""
Base class for unit models
"""
from __future__ import absolute_import  # disable implicit relative imports
from __future__ import division, print_function

import logging

from pyomo.environ import Reference, SolverFactory
from pyomo.network import Port
from pyomo.opt import TerminationCondition
from pyomo.common.config import ConfigValue, In

from .process_base import (declare_process_block_class,
                           ProcessBlockData,
                           useDefault)
from .property_base import StateBlockBase
from .control_volume_base import ControlVolumeBase, FlowDirection
from idaes.core.util.exceptions import (BurntToast,
                                        ConfigurationError,
                                        DynamicError,
                                        PropertyPackageError)
from idaes.core.util.misc import add_object_reference

__author__ = "John Eslick, Qi Chen, Andrew Lee"


__all__ = ['UnitBlockData', 'UnitBlock']

# Set up logger
_log = logging.getLogger(__name__)


@declare_process_block_class("UnitBlock")
class UnitBlockData(ProcessBlockData):
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

    def build(self):
        """
        General build method for UnitBlockData. This method calls a number
        of sub-methods which automate the construction of expected attributes
        of unit models.

        Inheriting models should call `super().build`.

        Args:
            None

        Returns:
            None
        """
        super(UnitBlockData, self).build()

        # Set up dynamic flag and time domain
        self._setup_dynamics()

    def _setup_dynamics(self):
        """
        This method automates the setting of the dynamic flag and time domain
        for unit models.

        Performs the following:
         1) Determines if this is a top level flowsheet
         2) Gets dynamic flag from parent if not top level, or checks validity
            of argument provided
         3) Gets time domain from parent, or creates domain if top level model
         4) Checks has_holdup flag if present and dynamic = True

        Args:
            None

        Returns:
            None
        """
        # Check the dynamic flag, and retrieve if necessary
        if self.config.dynamic == useDefault:
            # Get dynamic flag from parent
            try:
                self.config.dynamic = self.parent_block().config.dynamic
            except AttributeError:
                # If parent does not have dynamic flag, raise Exception
                raise DynamicError('{} has a parent model '
                                   'with no dynamic attribute.'
                                   .format(self.name))

        # Check for case when dynamic=True, but parent dynamic=False
        if (self.config.dynamic and not self.parent_block().config.dynamic):
            raise DynamicError('{} trying to declare a dynamic model within '
                               'a steady-state flowsheet. This is not '
                               'supported by the IDAES framework. Try '
                               'creating a dynamic flowsheet instead, and '
                               'declaring some models as steady-state.'
                               .format(self.name))

        # Try to get reference to time object from parent
        try:
            # Guess parent is top level flowsheet, has time domain
            add_object_reference(self, "time_ref", self.parent_block().time)
        except AttributeError:
            # Try to look for a reference to time domain (time_ref)
            try:
                add_object_reference(self,
                                     "time_ref",
                                     self.parent_block().time_ref)
            except AttributeError:
                # Can't find time domain
                raise DynamicError('{} has a parent model '
                                   'with no time domain'.format(self.name))

        # Check has_holdup, if present
        if self.config.dynamic:
            if hasattr(self.config, "has_holdup"):
                if not self.config.has_holdup:
                    # Dynamic model must have has_holdup = True
                    raise ConfigurationError(
                            "{} invalid arguments for dynamic and has_holdup. "
                            "If dynamic = True, has_holdup must also be True "
                            "(was False)".format(self.name))

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
            name = name to use for Port object.
            block = an instance of a StateBlock to use as the source to
                    populate the Port object
            doc = doc string for Port object
        Returns:
            A Pyomo Port object and associated components.
        """
        # Validate block object
        if not isinstance(block, StateBlockBase):
            raise ConfigurationError("{} block object provided to add_port "
                                     "method is not an instance of a "
                                     "StateBlock object. IDAES port objects "
                                     "should only be associated with "
                                     "StateBlocks.".format(blk.name))

        # Create empty Port
        p = Port(noruleinit=True, doc=doc)
        setattr(blk, name, p)

        # Get dict of Port members and names
        member_list = block[blk.time_ref.first()].define_port_members()

        # Create References for port members
        for s in member_list:
            if not member_list[s].is_indexed():
                slicer = block[:].component(member_list[s].local_name)
            else:
                slicer = block[:].component(member_list[s].local_name)[...]

            r = Reference(slicer)
            setattr(blk, "_"+name+"_"+s+"_ref", r)

            # Add Reference to Port
            p.add(r, s)

        return p

    def add_inlet_port(blk, name=None, block=None, doc=None):
        """
        This is a method to build inlet Port objects in a unit model and
        connect these to a specified control volume or state block.

        Keyword Args:
            name = name to use for Port object (default = "inlet").
            block = an instance of a ControlVolume or StateBlock to use as the
                    source to populate the Port object. If a ControlVolume is
                    provided, the method will use the inlet state block as
                    defined by the ControlVolume. If not provided, method will
                    attempt to default to an object named control_volume.
            doc = doc string for Port object (default = "Inlet Port")

        Returns:
            A Pyomo Port object and associated components.
        """
        if block is None:
            # Try for default ControlVolume name
            try:
                block = blk.control_volume
            except AttributeError:
                raise ConfigurationError(
                        "{} add_inlet_port was called without a block argument"
                        " but no default ControlVolume exists "
                        "(control_volume). Please provide block to which the "
                        "Port should be associated.".format(blk.name))

        # Set default name and doc string if not provided
        # TODO: This does not work when the block has multiple control volumes.
        # Better to raise exception when user does not provide a name than us
        # setting a name automatically.
        if name is None:
            name = "inlet"
        if doc is None:
            doc = "Inlet Port"

        # Create empty Port
        p = Port(noruleinit=True, doc=doc)
        setattr(blk, name, p)

        # Get dict of Port members and names
        if isinstance(block, ControlVolumeBase):
            try:
                member_list = (block.properties_in[block.time_ref.first()]
                               .define_port_members())
            except AttributeError:
                try:
                    member_list = (block.properties[block.time_ref.first(), 0]
                                   .define_port_members())
                except AttributeError:
                    raise PropertyPackageError(
                            "{} property package does not appear to have "
                            "implemented a define_port_memebers method. "
                            "Please contact the developer of the property "
                            "package.".format(blk.name))
        elif isinstance(block, StateBlockBase):
            member_list = block[blk.time_ref.first()].define_port_members()
        else:
            raise ConfigurationError(
                    "{} block provided to add_inlet_port "
                    "method was not an instance of a "
                    "ControlVolume or a StateBlock."
                    .format(blk.name))

        # Create References for port members
        for s in member_list:
            if not member_list[s].is_indexed():
                if isinstance(block, ControlVolumeBase):
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
                elif isinstance(block, StateBlockBase):
                    slicer = block[:].component(member_list[s].local_name)
                else:
                    raise ConfigurationError(
                            "{} block provided to add_inlet_port "
                            "method was not an instance of a "
                            "ControlVolume or a StateBlock."
                            .format(blk.name))
            else:
                if isinstance(block, ControlVolumeBase):
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
                elif isinstance(block, StateBlockBase):
                    slicer = block[:].component(member_list[s].local_name)[...]
                else:
                    raise ConfigurationError(
                            "{} block provided to add_inlet_port "
                            "method was not an instance of a "
                            "ControlVolume or a StateBlock."
                            .format(blk.name))

            r = Reference(slicer)
            setattr(blk, "_"+name+"_"+s+"_ref", r)

            # Add Reference to Port
            p.add(r, s)

        return p

    def add_outlet_port(blk, name=None, block=None, doc=None):
        """
        This is a method to build outlet Port objects in a unit model and
        connect these to a specified control volume or state block.

        Keyword Args:
            name = name to use for Port object (default = "outlet").
            block = an instance of a ControlVolume or StateBlock to use as the
                    source to populate the Port object. If a ControlVolume is
                    provided, the method will use the outlet state block as
                    defined by the ControlVolume. If not provided, method will
                    attempt to default to an object named control_volume.
            doc = doc string for Port object (default = "Outlet Port")

        Returns:
            A Pyomo Port object and associated components.
        """
        if block is None:
            # Try for default ControlVolume name
            try:
                block = blk.control_volume
            except AttributeError:
                raise ConfigurationError(
                        "{} add_outlet_port was called without a block "
                        "argument but no default ControlVolume exists "
                        "(control_volume). Please provide block to which the "
                        "Port should be associated.".format(blk.name))

        # Set default name and doc string if not provided
        if name is None:
            name = "outlet"
        if doc is None:
            doc = "Outlet Port"

        # Create empty Port
        p = Port(noruleinit=True, doc=doc)
        setattr(blk, name, p)

        # Get dict of Port members and names
        if isinstance(block, ControlVolumeBase):
            try:
                member_list = (block.properties_out[block.time_ref.first()]
                               .define_port_members())
            except AttributeError:
                try:
                    member_list = (block.properties[block.time_ref.first(), 0]
                                   .define_port_members())
                except AttributeError:
                    raise PropertyPackageError(
                            "{} property package does not appear to have "
                            "implemented a define_port_memebers method. "
                            "Please contact the developer of the property "
                            "package.".format(blk.name))
        elif isinstance(block, StateBlockBase):
            member_list = block[blk.time_ref.first()].define_port_members()
        else:
            raise ConfigurationError(
                    "{} block provided to add_inlet_port "
                    "method was not an instance of a "
                    "ControlVolume or a StateBlock."
                    .format(blk.name))

        # Create References for port members
        for s in member_list:
            if not member_list[s].is_indexed():
                if isinstance(block, ControlVolumeBase):
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
                elif isinstance(block, StateBlockBase):
                    slicer = block[:].component(member_list[s].local_name)
                else:
                    raise ConfigurationError(
                            "{} block provided to add_inlet_port "
                            "method was not an instance of a "
                            "ControlVolume or a StateBlock."
                            .format(blk.name))
            else:
                # Need to use slice notation on indexed comenent as well
                if isinstance(block, ControlVolumeBase):
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
                elif isinstance(block, StateBlockBase):
                    slicer = block[:].component(member_list[s].local_name)[...]
                else:
                    raise ConfigurationError(
                            "{} block provided to add_inlet_port "
                            "method was not an instance of a "
                            "ControlVolume or a StateBlock."
                            .format(blk.name))

            r = Reference(slicer)
            setattr(blk, "_"+name+"_"+s+"_ref", r)

            # Add Reference to Port
            p.add(r, s)

        return p

    def initialize(blk, state_args=None, outlvl=0,
                   solver='ipopt', optarg={'tol': 1e-6}):
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
            outlvl : sets output level of initialisation routine

                     * 0 = no output (default)
                     * 1 = return solver state for each step in routine
                     * 2 = return solver state for each step in subroutines
                     * 3 = include solver output infomation (tee=True)

            optarg : solver options dictionary object (default={'tol': 1e-6})
            solver : str indicating which solver to use during
                     initialization (default = 'ipopt')

        Returns:
            None
        '''
        # Set solver options
        if outlvl > 3:
            stee = True
        else:
            stee = False

        opt = SolverFactory(solver)
        opt.options = optarg

        # ---------------------------------------------------------------------
        # Initialize control volume block
        flags = blk.control_volume.initialize(outlvl=outlvl-1,
                                              optarg=optarg,
                                              solver=solver,
                                              state_args=state_args)

        if outlvl > 0:
            _log.info('{} Initialisation Step 1 Complete.'.format(blk.name))

        # ---------------------------------------------------------------------
        # Solve unit
        results = opt.solve(blk, tee=stee)

        if outlvl > 0:
            if results.solver.termination_condition == \
                    TerminationCondition.optimal:
                _log.info('{} Initialisation Step 2 Complete.'
                          .format(blk.name))
            else:
                _log.warning('{} Initialisation Step 2 Failed.'
                             .format(blk.name))

        # ---------------------------------------------------------------------
        # Release Inlet state
        blk.control_volume.release_state(flags, outlvl-1)

        if outlvl > 0:
            _log.info('{} Initialisation Complete.'.format(blk.name))
