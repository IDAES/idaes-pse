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
Base class for unit models
"""
from __future__ import absolute_import  # disable implicit relative imports
from __future__ import division, print_function

import logging

from pyomo.environ import SolverFactory
from pyomo.network import Port
from pyomo.opt import TerminationCondition
from pyomo.common.config import ConfigValue, In

from .process_base import declare_process_block_class, ProcessBlockData
from idaes.core.util.exceptions import DynamicError

__author__ = "John Eslick, Qi Chen, Andrew Lee"


__all__ = ['UnitBlockData', 'UnitBlock']

# Set up logger
logger = logging.getLogger('idaes.core')
unit_logger = logging.getLogger('idaes.unit_model')


@declare_process_block_class("UnitBlock")
class UnitBlockData(ProcessBlockData):
    """
    This is the class for process unit operations models. These are models that
    would generally appear in a process flowsheet or superstructure.
    """
    # Create Class ConfigBlock
    CONFIG = ProcessBlockData.CONFIG()
    CONFIG.declare("dynamic", ConfigValue(
        domain=In([True, False]),
        description="Dynamic model flag",
        doc="""Indicates whether this model will be dynamic or not
(default = None).
None - get flag from parent (default = False)
True - set as a dynamic model
False - set as a steady-state model"""))

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
         4) Checks include_holdup flag if present and dynamic = True

        Args:
            None

        Returns:
            None
        """
        # Check the dynamic flag, and retrieve if necessary
        if self.config.dynamic is None:
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
            object.__setattr__(self, "time", self.parent_block().time)
        except AttributeError:
            raise DynamicError('{} has a parent model '
                               'with no time domain'.format(self.name))

        # Check include_holdup, if present
        if self.config.dynamic:
            if hasattr(self.config, "include_holdup"):
                if not self.config.include_holdup:
                    # Dynamic model must have include_holdup = True
                    logger.warning('{} Dynamic models must have '
                                   'include_holdup = True. '
                                   'Overwritting argument.'
                                   .format(self.name))
                    self.config.include_holdup = True

#    def build_inlets(self, holdup=None, inlets=None, num_inlets=None):
#        """
#        This is a method to build inlet Port objects in a unit model and
#        connect these to holdup blocks as needed. This method supports an
#        arbitary number of inlets and holdup blocks, and works for both
#        simple (0D) and 1D IDAES holdup blocks.
#
#        Keyword Args:
#            holdup = holdup block to which inlets are associated. If left None,
#                    assumes a default holdup (default = None).
#            inlets = argument defining inlet names (default: None). inlets may
#                    be None or list.
#                    - None - assumes a single inlet.
#                    - list - use names provided in list for inlets (can be
#                            other iterables, but not a string or dict)
#            num_inlets = argument indication number (int) of inlets to
#                    construct (default = None). Not used if inlets arg is
#                    provided.
#                    - None - use inlets arg instead
#                    - int - Inlets will be named with sequential numbers from 1
#                            to num_inlets.
#
#        Returns:
#            A Pyomo Port object and assoicated components.
#        """
#        # Check holdup argument, and get holdup block
#        if holdup is None:
#            # If None, assume default names
#            try:
#                hblock = self.holdup
#                inlet_name = 'inlet'
#            except AttributeError:
#                raise AttributeError('{} Invalid holdup argument. Unit model '
#                                     'contains no attribute holdup.'
#                                     .format(self.name))
#        else:
#            # Otherwise, use named holdup and inlet
#            try:
#                hblock = getattr(self, holdup)
#                inlet_name = holdup+'_inlet'
#            except AttributeError:
#                raise AttributeError('{} Invalid holdup argument. Unit model '
#                                     'contains no attribute {}.'
#                                     .format(self.name, holdup))
#
#        # Validate inlets argument and create inlet_list if needed
#        inlet_list = []
#        if inlets is None:
#            if num_inlets is not None:
#                # Create list of integers as strings and put in inlet_list
#                inlet_list = list(map(str, range(1, num_inlets+1)))
#            else:
#                # No arguments provided, assume single default inlet
#                inlet_list = [None]
#        else:
#            # List of named inlets, should be iterable but not string or dict
#            try:
#                # Test iterable
#                iter(inlets)
#                # Raise exception if string or dict
#                if isinstance(inlets, (str, dict)):
#                    raise TypeError('{} Invalid inlets argument. Must be '
#                                    'iterable but not a string or a dict.'
#                                    .format(self.name))
#                # Ensure that all elements are strings
#                for i in inlets:
#                    inlet_list += [str(i)]
#            except TypeError:
#                # Unrecognised type for inlets
#                raise TypeError('{} Unrecognised type for inlets argument.'
#                                .format(self.name))
#
#            # Check if both inlets and num_inlets args were provided
#            if num_inlets is not None:
#                # Log a warning
#                if num_inlets != len(inlets):
#                    # If args don't agree, raise Exception
#                    raise ValueError("{} provided with both inlets and "
#                                     "num_inlets arguments which do "
#                                     "not agree. Check the values assigned to "
#                                     "these. It is also only necessary to "
#                                     "specify one of these arguments."
#                                     .format(self.name))
#                else:
#                    # Otherwise just log an debug message
#                    logger.debug("{} provided with both inlets and "
#                                 "num_inlets arguments - num_inlets will be "
#                                 "ignored.".format(self.name))
#
#        # Check for duplicate inlet names
#        if inlets is not None:
#            if any(inlet_list.count(x) > 1 for x in inlet_list):
#                raise ValueError('{} Duplicate inlet name found.'
#                                 .format(self.name))
#
#        if len(inlet_list) > 1:
#            # Multiple inlets to holdup, so call InletMixer
#            hblock.inlet_mixer = InletMixer(doc="Mixer for multiple inlets.",
#                                            inlets=inlet_list)
#
#        # Build inlet port object and populate
#        if len(inlet_list) == 1:
#            # Only one inlet, index only by time
#            def inlet_rule(b, t):
#                try:
#                    return hblock.properties_in[t].declare_port_members()
#                except AttributeError:
#                    if hasattr(hblock, "ldomain"):
#                        if hblock.config.flow_direction is "forward":
#                            prop_key = (t, hblock.ldomain.first())
#                        else:
#                            prop_key = (t, hblock.ldomain.last())
#                    else:
#                        prop_key = t
#                    return hblock.properties[prop_key].declare_port_members()
#            i = Port(self.time,
#                     rule=inlet_rule,
#                     doc="Inlet port object")
#            setattr(self, inlet_name, i)
#        else:
#            # Multiple inlets, need to index conenctor
#            # Create inlet port and populate later
#            i = Port(self.time,
#                     hblock.inlet_mixer.inlet_idx,
#                     noruleinit=True,
#                     doc="Inlet port object")
#            setattr(self, inlet_name, i)
#            inlet_obj = getattr(self, inlet_name)
#
#            # Get the assoicated property block from inlet_mixer
#            pblock = hblock.inlet_mixer.properties
#
#            # Iterate over time and inlets
#            for i in hblock.inlet_mixer.inlet_idx:
#                for t in self.time:
#                    # Get the member of the port to add
#                    members = pblock[t, i].declare_port_members()
#                    for obj in members:
#                        inlet_obj[t, i].add(members[obj], obj)
#
#    def build_outlets(self, holdup=None, outlets=None,
#                      num_outlets=None, material_split_type='flow',
#                      energy_split_type='temperature'):
#        """
#        This is a method to build outlet Port objects in a unit model and
#        connect these to holdup blocks as needed. This method supports an
#        arbitary number of outlets and holdup blocks, and works for both
#        simple (0D) and 1D IDAES holdup blocks.
#
#        Keyword Args:
#            holdup = holdup block to which inlets are associated. If left None,
#                    assumes a default holdup (default = None).
#            outlets = argument defining outlet names (default: None). outlets
#                    may be None or list.
#                    - None - assumes a single outlet.
#                    - list - use names provided in list for outlets (can be
#                            other iterables, but not a string or dict)
#            num_outlets = argument indication number (int) of outlets to
#                    construct (default = None). Not used if outlets arg is
#                    provided.
#                    - None - use outlets arg instead
#                    - int - Outlets will be named with sequential numbers from
#                            1 to num_outlets.
#            material_split_type = argument defining method to use to split
#                    outlet material flow in case of multiple outlets
#                    (default = 'flow').
#                        - 'flow' - outlets are split by total flow
#                        - 'phase' - outlets are split by phase
#                        - 'component' - outlets are split by component
#                        - 'total' - outlets are split by both phase and
#                                    component
#                        - 'duplicate' - all outlets are duplicates of the
#                                        total outlet stream.
#            energy_split_type = argument defining method to use to split
#                    outlet energy flow in case of multiple outlets
#                    (default = 'temperature').
#                        - 'temperature' - equate temperatures in outlets
#                        - 'enth_mol' - equate molar enthalpies in outlets
#                        - 'enth_mass' - equate mass enthalpies in outlets
#                        - 'energy_balance' - outlets energy split using split
#                                fractions
#
#        Returns:
#            A Pyomo Port object and assoicated components.
#        """
#        # Check holdup argument, and get holdup block
#        if holdup is None:
#            # If None, assume default names
#            try:
#                hblock = self.holdup
#                outlet_name = 'outlet'
#            except AttributeError:
#                raise AttributeError('{} Invalid holdup argument. Unit model '
#                                     'contains no attribute holdup.'
#                                     .format(self.name))
#        else:
#            # Otherwise, use named holdup and outlet
#            try:
#                hblock = getattr(self, holdup)
#                outlet_name = holdup+'_outlet'
#            except AttributeError:
#                raise AttributeError('{} Invalid holdup argument. Unit model '
#                                     'contains no attribute {}.'
#                                     .format(self.name, holdup))
#
#        # Validate inlets argument and create outlet_list
#        outlet_list = []
#        if outlets is None:
#            if num_outlets is not None:
#                # Create list of integers as strings and put in outlet_list
#                outlet_list = list(map(str, range(1, num_outlets+1)))
#            else:
#                # No arguments provided, assume single default outlet
#                outlet_list = [None]
#        else:
#            # List of named outlets, should be iterable but not string or dict
#            try:
#                # Test iterable
#                iter(outlets)
#                # Raise exception if string or dict
#                if isinstance(outlets, (str, dict)):
#                    raise TypeError('{} Invalid outlets argument. Must be '
#                                    'iterable but not a string or a dict.'
#                                    .format(self.name))
#                # Ensure that all elements are strings
#                for o in outlets:
#                    outlet_list += [str(o)]
#            except TypeError:
#                # Unrecognised type for outlets
#                raise TypeError('{} Unrecognised type for outlets argument.'
#                                .format(self.name))
#
#            # Check if both outlets and num_outlets args were provided
#            if num_outlets is not None:
#                # Log a warning
#                if num_outlets != len(outlets):
#                    # If args don't agree, raise Excpetion
#                    raise ValueError("{} provided with both outlets and "
#                                     "num_outlets arguments which do "
#                                     "not agree. Check the values assigned to "
#                                     "these. It is also only necessary to "
#                                     "specify one of these arguments."
#                                     .format(self.name))
#                else:
#                    # Otherwise just log an info message
#                    logger.debug("{} provided with both outlets and "
#                                 "num_outlets arguments - num_outlets will be "
#                                 "ignored.".format(self.name))
#
#        # Check for duplicate outlet names
#        if outlets is not None:
#            if any(outlet_list.count(x) > 1 for x in outlet_list):
#                raise ValueError('{} Duplicate outlet name found.'
#                                 .format(self.name))
#
#        # Check that material_split_type argument is valid
#        valid_args = ['flow', 'phase', 'component', 'total', 'duplicate']
#        if material_split_type not in valid_args:
#            raise ValueError('{} Unrecognised value for material_split_type '
#                             'argument.'.format(self.name))
#
#        if len(outlet_list) > 1:
#            # Multiple outlets to holdup, so check if splitter is required
#            # If material_split_type is duplicate, no splitter is required
#            if material_split_type != 'duplicate':
#                # Call OutletSplitter
#                hblock.outlet_splitter = OutletSplitter(
#                        doc="Splitter for multiple outlets.",
#                        outlets=outlet_list,
#                        material_split_type=material_split_type,
#                        energy_split_type=energy_split_type)
#
#        # Build outlet port object and populate
#        if len(outlet_list) == 1:
#            # Only one outlet, index only by time
#            def outlet_rule(b, t):
#                try:
#                    return hblock.properties_out[t].declare_port_members()
#                except AttributeError:
#                    if hasattr(hblock, "ldomain"):
#                        if hblock.config.flow_direction is "forward":
#                            prop_key = (t, hblock.ldomain.last())
#                        else:
#                            prop_key = (t, hblock.ldomain.first())
#                    else:
#                        prop_key = t
#                    return hblock.properties[prop_key].declare_port_members()
#            o = Port(self.time,
#                     rule=outlet_rule,
#                     doc="Outlet port object")
#            setattr(self, outlet_name, o)
#        else:
#            # Multiple inlets, need to index conenctor
#            # Create outlet port and populate later
#            o = Port(self.time,
#                     hblock.outlet_splitter.outlet_idx,
#                     noruleinit=True,
#                     doc="Outlet port object")
#            setattr(self, outlet_name, o)
#            outlet_obj = getattr(self, outlet_name)
#
#            # Check material_split_type to populate outlet
#            if material_split_type != 'duplicate':
#                # Get the assoicated property block from outlet_mixer
#                pblock = hblock.outlet_splitter.properties
#            else:
#                # Outlets are duplicates, so no splitter
#                pblock = hblock.properties_out
#
#            # Iterate over time and inlets
#            for i in hblock.outlet_splitter.outlet_idx:
#                for t in self.time:
#                    # Get the member of the port to add
#                    members = pblock[t, i].declare_port_members()
#                    for obj in members:
#                        outlet_obj[t, i].add(members[obj], obj)

    def model_check(blk):
        """
        This is a general purpose initialization routine for simple unit
        models. This method assumes a single Holdup block called holdup and
        tries to call the model_check method of the holdup block. If an
        AttributeError is raised, the check is passed.

        More complex models should overload this method with a model_check
        suited to the particular application, especially if there are multiple
        Holdup blocks present.

        Args:
            None

        Returns:
            None
        """
        # Run holdup block model checks
        try:
            blk.holdup.model_check()
        except AttributeError:
            pass

    def initialize(blk, state_args=None, outlvl=0,
                   solver='ipopt', optarg={'tol': 1e-6}):
        '''
        This is a general purpose initialization routine for simple unit
        models. This method assumes a single Holdup block called holdup, and
        first initializes this and then attempts to solve the entire unit.

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
            solver : str indicating whcih solver to use during
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
        # Initialize holdup block
        flags = blk.holdup.initialize(outlvl=outlvl-1,
                                      optarg=optarg,
                                      solver=solver,
                                      state_args=state_args)

        if outlvl > 0:
            unit_logger.info('{} Initialisation Step 1 Complete.'
                             .format(blk.name))

        # ---------------------------------------------------------------------
        # Solve unit
        results = opt.solve(blk, tee=stee)

        if outlvl > 0:
            if results.solver.termination_condition == \
                    TerminationCondition.optimal:
                unit_logger.info('{} Initialisation Step 2 Complete.'
                                 .format(blk.name))
            else:
                unit_logger.warning('{} Initialisation Step 2 Failed.'
                                    .format(blk.name))

        # ---------------------------------------------------------------------
        # Release Inlet state
        blk.holdup.release_state(flags, outlvl-1)

        if outlvl > 0:
            unit_logger.info('{} Initialisation Complete.'.format(blk.name))
