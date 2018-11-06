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
This module contains the base class for constructing flowsheet models in the
IDAES modeling framework.
"""

from __future__ import division, print_function

import logging

import pyomo.environ as pe
from pyomo.dae import ContinuousSet
from pyomo.common.config import ConfigValue, In

from idaes.core import (ProcessBlockData, declare_process_block_class,
                        UnitBlockData, useDefault)
from idaes.core.util.config import is_property_parameter_block, list_of_floats
from idaes.core.util.exceptions import ConfigurationError, DynamicError

# Some more information about this module
__author__ = "John Eslick, Qi Chen, Andrew Lee"


__all__ = ['FlowsheetBlock, FlowsheetBlockData']

# Set up logger
logger = logging.getLogger(__name__)


@declare_process_block_class("FlowsheetBlock", doc="""
    FlowsheetBlock is a specialized Pyomo block for IDAES flowsheet models, and
    contains instances of FlowsheetBlockData.""")
class FlowsheetBlockData(ProcessBlockData):
    """
    The FlowsheetBlockData Class forms the base class for all IDAES process
    flowsheet models. The main purpose of this class is to automate the tasks
    common to all flowsheet models and ensure that the necessary attributes of
    a flowsheet model are present.

    The most signfiicant role of the FlowsheetBlockData class is to
    automatically create the time domain for the flowsheet.
    """

    # Create Class ConfigBlock
    CONFIG = ProcessBlockData.CONFIG()
    CONFIG.declare("dynamic", ConfigValue(
        default=useDefault,
        domain=In([useDefault, True, False]),
        description="Dynamic model flag",
        doc="""Indicates whether this model will be dynamic,
**default** - useDefault.  **Valid values:** {
**useDefault** - get flag from parent or False,
**True** - set as a dynamic model,
**False** - set as a steady-state model}"""))
    CONFIG.declare("time_set", ConfigValue(
        default=[0],
        domain=list_of_floats,
        description="Set of points for initializing time domain",
        doc="""Set of points for initializing time domain. This should be a
        list of floating point numbers, **default** - [0]."""))
    CONFIG.declare("default_property_package", ConfigValue(
        default=None,
        domain=is_property_parameter_block,
        description="Default property package to use in flowsheet",
        doc="""Indicates the default property package to be used by models
within this flowsheet if not otherwise specified, **default** - None.
**Valid values:** {**None** - no default property package,
**a ParameterBlock object**.}"""))

    def build(self):
        """
        General build method for FlowsheetBlockData. This method calls a number
        of sub-methods which automate the construction of expected attributes
        of flowsheets.

        Inheriting models should call `super().build`.

        Args:
            None

        Returns:
            None
        """
        super(FlowsheetBlockData, self).build()

        # Set up dynamic flag and time domain
        self._setup_dynamics()

    # TODO [Qi]: this should be implemented as a transformation
    def model_check(self):
        """
        This method runs model checks on all unit models in a flowsheet.

        This method searches for objects which inherit from UnitBlockData and
        executes the model_check method if it exists.

        Args:
            None

        Returns:
            None
        """
        logger.info("Executing model checks.")
        for o in self.component_objects(descend_into=False):
            if isinstance(o, UnitBlockData):
                try:
                    o.model_check()
                except AttributeError:
                    logger.warning('{} Model/block has no model check. To '
                                   'correct this, add a model_check method to '
                                   'the associated unit model class'
                                   .format(o.name))

    def _setup_dynamics(self):
        """
        This method automates the setting of the dynamic flag and time domain
        for flowsheet models.

        Performs the following:
         1) Determines if this is a top level flowsheet
         2) Gets dynamic flag from parent if not top level, or checks validity
            of argument provided
         3) Gets time domain from parent, or creates domain if top level model

        Args:
            None

        Returns:
            None
        """
        # Test to ensure model is constructed
        if not self._constructed:
            raise ConfigurationError('{} flowsheet has no parent object but '
                                     'has not yet been constructed. Either '
                                     'attach the flowsheet to a ConcreteModel '
                                     'or use the argument concrete = True.'
                                     .format(self.name))

        # Determine if this is top level flowsheet
        if self.parent_block() is None:
            # Flowsheet has no parent, so top level model
            # Check that flowsheet has been set as concrete
            if not self._constructed:
                raise DynamicError('{} flowsheet has no parent object but '
                                   'has not yet been constructed. Either '
                                   'attach the flowsheet to a ConcreteModel '
                                   'or use the argument concrete = True.'
                                   .format(self.name))
            top_level = True
        elif isinstance(self.parent_block(), pe.ConcreteModel):
            # If flowsheet is attached to a ConcreteModel, this is a top level
            # flowsheet. However, check if the ConcreteModel has time
            try:
                # If parent has time, it must be a ContinuousSet or Set
                if isinstance(self.parent_block().time,
                              (ContinuousSet, pe.Set)):
                    logger.warning('{} parent ConcreteModel has an attribute '
                                   'time. Flowsheet will use this as its time '
                                   'domain, however this may be unexpected '
                                   'behaviour'.format(self.name))
                    top_level = False
                else:
                    raise DynamicError('{} has an attribute time which is not '
                                       'a Set or ContinuousSet.'
                                       .format(self.name))
            except AttributeError:
                # Set top level flag
                top_level = True
        else:
            top_level = False

        # Check the dynamic flag, and retrieve if necessary
        if self.config.dynamic == useDefault:
            # Check to see if this is a top level model
            if top_level:
                # If there is no parent, set dynamic to False by default and
                # warn the user
                logger.warning('{} is a top level flowhseet, but dynamic flag '
                               'set to "use_parent"value". Dynamic '
                               'flag set to False by default'
                               .format(self.name))
                self.config.dynamic = False
            else:
                # Get dynamic flag from parent
                try:
                    self.config.dynamic = self.parent_block().config.dynamic
                except AttributeError:
                    # If parent does not have dynamic flag, raise Exception
                    raise DynamicError('{} has a parent model '
                                       'with no dynamic attribute.'
                                       .format(self.name))

        # Check for case when dynamic=True, but parent dynamic=False
        if (not top_level and self.config.dynamic is True and
                not self.parent_block().config.dynamic):
            raise DynamicError('{} trying to declare a dynamic model within '
                               'a steady-state flowsheet. This is not '
                               'supported by the IDAES framework. Try '
                               'creating a dynamic flowsheet instead, and '
                               'declaring some models as  steady-state.'
                               .format(self.name))

        # Set up time domain
        if not top_level:
            # Try to get reference to time object from parent
            try:
                object.__setattr__(self, "time", self.parent_block().time)
            except AttributeError:
                raise DynamicError('{} has a parent model '
                                   'with no time domain'.format(self.name))
        else:
            # Create time domain
            if self.config.dynamic:
                # Check if time_set has at least two points
                if len(self.config.time_set) < 2:
                    # Check if time_set is at default value
                    if self.config.time_set == [0.0]:
                        # If default, set default end point to be 1.0
                        self.config.time_set = [0.0, 1.0]
                    else:
                        # Invalid user input, raise Excpetion
                        raise DynamicError("Flowsheet provided with invalid "
                                           "time_set attribute - must have at "
                                           "least two values (start and end).")
                # For dynamics, need a ContinuousSet
                self.time = ContinuousSet(initialize=self.config.time_set)
            else:
                # For steady-state, use an ordered Set
                self.time = pe.Set(initialize=self.config.time_set,
                                   ordered=True)
