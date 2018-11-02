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

from .process_base import (declare_process_block_class,
                           ProcessBlockData,
                           useDefault)
from idaes.core.util.exceptions import ConfigurationError, DynamicError

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
        doc="""Indicates whether this model will be dynamic or not
(default = useDefault).
useDefault - get flag from parent (default = False)
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
            # TODO : Replace with Reference
            object.__setattr__(self, "time", self.parent_block().time)
        except AttributeError:
            raise DynamicError('{} has a parent model '
                               'with no time domain'.format(self.name))

        # Check include_holdup, if present
        if self.config.dynamic:
            if hasattr(self.config, "include_holdup"):
                if not self.config.include_holdup:
                    # Dynamic model must have include_holdup = True
                    raise ConfigurationError(
                            "{} invalid arguments for dynamic and has_holdup. "
                            "If dynamic = True, has_holdup must also be True "
                            "(was False)".format(self.name))

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
        blk.holdup.release_state(flags, outlvl-1)

        if outlvl > 0:
            _log.info('{} Initialisation Complete.'.format(blk.name))
