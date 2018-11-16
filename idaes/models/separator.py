##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes".
##############################################################################
"""
General purpose separator block for IDAES models
"""
from __future__ import absolute_import  # disable implicit relative imports
from __future__ import division, print_function

import logging

from pyomo.environ import Constraint, Param, PositiveReals, Reals, Set, Var
from pyomo.common.config import ConfigBlock, ConfigValue, In
from pyutilib.enum import Enum

from idaes.core import (declare_process_block_class,
                        UnitBlockData,
                        useDefault)
from idaes.core.util.config import (is_physical_parameter_block,
                                    is_state_block,
                                    list_of_strings)
from idaes.core.util.exceptions import (BurntToast,
                                        ConfigurationError,
                                        PropertyNotSupportedError)
from idaes.core.util.math import smooth_min

__author__ = "Andrew Lee"


# Set up logger
_log = logging.getLogger(__name__)


@declare_process_block_class("SeparatorBlock")
class SeparatorBlockData(UnitBlockData):
    """
    This is a general purpose model for a Separator block with the IDAES
    modeling framework. This block can be used either as a stand-alone
    Separator unit operation, or as a sub-model within another unit operation.

    This model creates a number of StateBlocks to represent the outgoing
    streams, then writes a set of phase-component material balances, an
    overall enthalpy balance (2 options), and a momentum balance (2 options)
    linked to a mixed-state StateBlock. The mixed-state StateBlock can either
    be specified by the user (allowing use as a sub-model), or created by the
    SeparatorBlock.

    When being used as a sub-model, SeparatorBlock should only be used when a
    set of new StateBlocks are required for the streams to be separated. It
    should not be used to separate streams to go to mutiple ControlVolumes in a
    single unit model - in these cases the unit model developer should write
    their own splitting equations.
    """
    CONFIG = UnitBlockData.CONFIG()
    CONFIG.declare("property_package", ConfigValue(
        default=useDefault,
        domain=is_physical_parameter_block,
        description="Property package to use for mixer",
        doc="""Property parameter object used to define property calculations,
**default** - useDefault. **Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PropertyParameterObject** - a PropertyParameterBlock object.}"""))
    CONFIG.declare("property_package_args", ConfigBlock(
        implicit=True,
        description="Arguments to use for constructing property packages",
        doc="""A ConfigBlock with arguments to be passed to a property block(s)
and used when constructing these, **default** - None. **Valid values:** {
see property package for documentation.}"""))
    CONFIG.declare("outlet_list", ConfigValue(
        domain=list_of_strings,
        description="List of outlet names",
        doc="""A list containing names of outlets, **default** - None.
**Valid values:** {**None** - use num_outlets argument, **list** - a list of
names to use for outlets.}"""))
    CONFIG.declare("num_outlets", ConfigValue(
        domain=int,
        description="Number of outlets to unit",
        doc="""Argument indicating number (int) of outlets to construct, not
used if outlet_list arg is provided, **default** - None. **Valid values:** {
**None** - use outlet_list arg instead, or default to 2 if neither argument
provided, **int** - number of outlets to create (will be named with sequential
integers from 1 to num_outlets).}"""))
    CONFIG.declare("calculate_phase_equilibrium", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="Calculate phase equilibrium in outlet streams",
        doc="""Argument indicating whether phase equilibrium should be
calculated for the resulting outlet streams, **default** - False. **Valid
values: ** {**True** - calculate phase equilibrium in outlet streams,
**False** - do not calculate equilibrium in outlet streams.}"""))
#    CONFIG.declare("momentum_mixing_type", ConfigValue(
#        default=MomentumMixingType.minimize,
#        domain=MomentumMixingType,
#        description="Method to use when mxing momentum/pressure",
#        doc="""Argument indicating what method to use when mixing momentum/
#pressure of incoming streams, **default** - MomentumMixingType.minimize.
#**Valid values:** {**MomentumMixingType.minimize** - mixed stream has
#pressure equal to the minimimum pressure of the incoming streams (uses
#smoothMin operator), **MomentumMixingType.equality** - enforces equality of
#pressure in mixed and all incoming streams.}"""))
    CONFIG.declare("mixed_state_block", ConfigValue(
        domain=is_state_block,
        description="Existing StateBlock to use as mixed stream",
        doc="""An existing state block to use as the source stream from the
Separator block, **default** - None. **Valid values:** {**None** - create a new
StateBlock for the mixed stream, **StateBlock** - a StateBock to use as the
source for the mixed stream.}"""))
    CONFIG.declare("construct_ports", ConfigValue(
        default=True,
        domain=In([True, False]),
        description="Construct inlet and outlet Port objects",
        doc="""Argument indicating whether model should construct Port objects
linked the mixed state and all outlet states, **default** - True.
**Valid values:** {**True** - construct Ports for all states, **False** - do
not construct Ports."""))

    def build(self):
        """
        General build method for SeparatorBlockData. This method calls a number
        of sub-methods which automate the construction of expected attributes
        of unit models.

        Inheriting models should call `super().build`.

        Args:
            None

        Returns:
            None
        """
        # Call super.build()
        super(SeparatorBlockData, self).build()

        # Call setup methods from ControlVolumeBase
        self._get_property_package()
        self._get_indexing_sets()

        # Create list of inlet names
        outlet_list = self.create_outlet_list()

        # Build StateBlocks
        outlet_blocks = self.add_outlet_state_blocks(outlet_list)

        if self.config.mixed_state_block is None:
            mixed_block = self.add_mixed_state_block()
        else:
            mixed_block = self.get_mixed_state_block()

        self.add_port_objects(outlet_list, outlet_blocks, mixed_block)

    def create_outlet_list(self):
        """
        Create list of outlet stream names based on config arguments.

        Returns:
            list of strings
        """
        if (self.config.outlet_list is not None and
                self.config.num_outlets is not None):
            # If both arguments provided and not consistent, raise Exception
            if len(self.config.outlet_list) != self.config.num_outlets:
                raise ConfigurationError(
                        "{} SeparatorBlock provided with both outlet_list and "
                        "num_outlets arguments, which were not consistent ("
                        "length of outlet_list was not equal to num_outlets). "
                        "PLease check your arguments for consistency, and "
                        "note that it is only necessry to provide one of "
                        "these arguments.".format(self.name))
        elif (self.config.outlet_list is None and
              self.config.num_outlets is None):
            # If no arguments provided for outlets, default to num_outlets = 2
            self.config.num_outlets = 2

        # Create a list of names for outlet StateBlocks
        if self.config.outlet_list is not None:
            outlet_list = self.config.outlet_list
        else:
            outlet_list = ['outlet_' + str(n)
                           for n in range(1, self.config.num_outlets+1)]

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
        tmp_dict = self.config.property_package_args
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = False

        # Create empty list to hold StateBlocks for return
        outlet_blocks = []

        # Create an instance of StateBlock for all outlets
        for o in outlet_list:
            o_obj = self._property_module.StateBlock(
                        self.time,
                        doc="Material properties at outlet",
                        default=tmp_dict)

            setattr(self, o+"_state", o_obj)

            outlet_blocks.append(getattr(self, o+"_state"))

        return outlet_blocks

    def add_mixed_state_block(self):
        """
        Constructs StateBlock to represent mixed stream.

        Returns:
            New StateBlock object
        """
        # Setup StateBlock argument dict
        tmp_dict = self.config.property_package_args
        tmp_dict["has_phase_equilibrium"] = \
            self.config.calculate_phase_equilibrium
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = True

        self.mixed_state = self._property_module.StateBlock(
                                self.time,
                                doc="Material properties of mixed stream",
                                default=tmp_dict)

        return self.mixed_state

    def get_mixed_state_block(self):
        """
        Validates StateBlock provided in user arguments for mixed stream.

        Returns:
            The user-provided StateBlock or an Exception
        """
        # Sanity check to make sure method is not called when arg missing
        if self.config.mixed_state_block is None:
            raise BurntToast("{} get_mixed_state_block method called when "
                             "mixed_state_block argument is None. This should "
                             "not happen.".format(self.name))

        # Check that the user-provided StateBlock uses the same prop pack
        if (self.config.mixed_state_block[self.time.first()].config.parameters
                != self.config.property_package):
            raise ConfigurationError(
                    "{} StateBlock provided in mixed_state_block argument "
                    " does not come from the same property package as "
                    "provided in the property_package argument. All "
                    "StateBlocks within a SeparatorBlock must use the same "
                    "property package.".format(self.name))

        return self.config.mixed_state_block

    def add_port_objects(self, outlet_list, outlet_blocks, mixed_block):
        """
        Adds Port objects if required.

        Args:
            a list of outlet StateBlock objects
            a mixed state StateBlock object

        Returns:
            None
        """
        if self.config.construct_ports is True:
            # Add ports
            for p in outlet_list:
                o_state = getattr(self, p+"_state")
                self.add_port(name=p, block=o_state, doc="Outlet Port")
            self.add_port(name="inlet", block=mixed_block, doc="Inlet Port")

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
        for t in blk.time:
            try:
                if blk.config.mixed_state_block is None:
                    blk.mixed_state[t].model_check()
                else:
                    blk.config.mixed_state_block.model_check()
            except AttributeError:
                _log.warning('{} SeparatorBlock inlet state block has no '
                             'model check. To correct this, add a '
                             'model_check method to the associated '
                             'StateBlock class.'.format(blk.name))

            try:
                outlet_list = blk.create_outlet_list()
                for o in outlet_list:
                    o_block = getattr(blk, o+"_state")
                    o_block[t].model_check()
            except AttributeError:
                _log.warning('{} SeparatorBlock outlet state block has no '
                             'model checks. To correct this, add a model_check'
                             ' method to the associated StateBlock class.'
                             .format(blk.name))

    def initialize(blk, outlvl=0, optarg=None,
                   solver='ipopt', hold_state=True):
        '''
        Initialisation routine for holdup (default solver ipopt)

        Keyword Arguments:
            outlvl : sets output level of initialisation routine. **Valid
                     values:** **0** - no output (default), **1** - return
                     solver state for each step in routine, **2** - include
                     solver output infomation (tee=True)
            optarg : solver options dictionary object (default=None)
            solver : str indicating whcih solver to use during
                     initialization (default = 'ipopt')
            hold_state : flag indicating whether the initialization routine
                     should unfix any state variables fixed during
                     initialization, **default** - True. **Valid values:**
                     **True** - states variables are not unfixed, and a dict of
                     returned containing flags for which states were fixed
                     during initialization, **False** - state variables are
                     unfixed after initialization by calling the release_state
                     method.

        Returns:
            If hold_states is True, returns a dict containing flags for which
            states were fixed during initialization.
        '''
        # Initialize inlet state blocks
        flags = {}
        inlet_list = blk.create_inlet_list()
        for i in inlet_list:
            i_block = getattr(blk, i+"_state")
            flags[i] = {}
            flags[i] = i_block.initialize(outlvl=outlvl-1,
                                          optarg=optarg,
                                          solver=solver,
                                          hold_state=hold_state)

        if blk.config.mixed_state_block is None:
            mblock = blk.mixed_state
        else:
            mblock = blk.config.mixed_state_block

        mblock.initialize(outlvl=outlvl-1,
                          optarg=optarg,
                          solver=solver,
                          hold_state=False)

        if outlvl > 0:
            _log.info('{} Initialisation Complete'.format(blk.name))

        return flags

    def release_state(blk, flags, outlvl=0):
        '''
        Method to release state variables fixed during initialisation.

        Keyword Arguments:
            flags : dict containing information of which state variables
                    were fixed during initialization, and should now be
                    unfixed. This dict is returned by initialize if
                    hold_state = True.
            outlvl : sets output level of logging

        Returns:
            None
        '''
        inlet_list = blk.create_inlet_list()
        for i in inlet_list:
            i_block = getattr(blk, i+"_state")
            i_block.release_state(flags[i], outlvl=outlvl-1)
