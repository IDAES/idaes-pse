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
This module contains the base class for constructing flowsheet models in the
IDAES modeling framework.
"""

import pyomo.environ as pe
from pyomo.dae import ContinuousSet
from pyomo.network import Arc
from pyomo.common.config import ConfigValue, ListOf
from pyomo.core.base.units_container import _PyomoUnit

from idaes.core import (
    ProcessBlockData,
    declare_process_block_class,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.util.config import (
    is_physical_parameter_block,
    is_time_domain,
    DefaultBool,
)
from idaes.core.util.misc import add_object_reference
from idaes.core.util.exceptions import DynamicError, ConfigurationError
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.ui.fsvis.fsvis import visualize

import idaes.logger as idaeslog

# Some more information about this module
__author__ = "John Eslick, Qi Chen, Andrew Lee"


__all__ = ["FlowsheetBlock", "FlowsheetBlockData"]

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class(
    "FlowsheetBlock",
    doc="""
    FlowsheetBlock is a specialized Pyomo block for IDAES flowsheet models, and
    contains instances of FlowsheetBlockData.""",
)
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
    CONFIG.declare(
        "dynamic",
        ConfigValue(
            default=useDefault,
            domain=DefaultBool,
            description="Dynamic model flag",
            doc="""Indicates whether this model will be dynamic,
**default** - useDefault.
**Valid values:** {
**useDefault** - get flag from parent or False,
**True** - set as a dynamic model,
**False** - set as a steady-state model.}""",
        ),
    )
    CONFIG.declare(
        "time",
        ConfigValue(
            default=None,
            domain=is_time_domain,
            description="Flowsheet time domain",
            doc="""Pointer to the time domain for the flowsheet. Users may provide
an existing time domain from another flowsheet, otherwise the flowsheet will
search for a parent with a time domain or create a new time domain and
reference it here.""",
        ),
    )
    CONFIG.declare(
        "time_set",
        ConfigValue(
            default=[0],
            domain=ListOf(float),
            description="Set of points for initializing time domain",
            doc="""Set of points for initializing time domain. This should be a
list of floating point numbers,
**default** - [0].""",
        ),
    )
    CONFIG.declare(
        "time_units",
        ConfigValue(
            description="Units for time domain",
            doc="""Pyomo Units object describing the units of the time domain.
This must be defined for dynamic simulations, default = None.""",
        ),
    )
    CONFIG.declare(
        "default_property_package",
        ConfigValue(
            default=None,
            domain=is_physical_parameter_block,
            description="Default property package to use in flowsheet",
            doc="""Indicates the default property package to be used by models
within this flowsheet if not otherwise specified,
**default** - None.
**Valid values:** {
**None** - no default property package,
**a ParameterBlock object**.}""",
        ),
    )

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

        self._time_units = None

        # Set up dynamic flag and time domain
        self._setup_dynamics()

    @property
    def time(self):
        # _time will be created by the _setup_dynamics method
        return self._time

    @property
    def time_units(self):
        return self._time_units

    def is_flowsheet(self):
        """
        Method which returns True to indicate that this component is a
        flowsheet.

        Args:
            None

        Returns:
            True
        """
        return True

    def model_check(self):
        """
        This method runs model checks on all unit models in a flowsheet.

        This method searches for objects which inherit from UnitModelBlockData
        and executes the model_check method if it exists.

        Args:
            None

        Returns:
            None
        """
        _log.info("Executing model checks.")
        for o in self.component_objects(descend_into=False):
            if isinstance(o, UnitModelBlockData):
                try:
                    o.model_check()
                except AttributeError:
                    # This should never happen, but just in case
                    _log.warning(
                        "{} Model/block has no model_check method.".format(o.name)
                    )

    def stream_table(self, true_state=False, time_point=0, orient="columns"):
        """
        Method to generate a stream table by iterating over all Arcs in the
        flowsheet.

        Args:
            true_state : whether the state variables (True) or display
                         variables (False, default) from the StateBlocks should
                         be used in the stream table.
            time_point : point in the time domain at which to create stream
                         table (default = 0)
            orient : whether stream should be shown by columns ("columns") or
                     rows ("index")

        Returns:
            A pandas dataframe containing stream table information
        """
        dict_arcs = {}

        for a in self.component_objects(ctype=Arc, descend_into=False):
            dict_arcs[a.local_name] = a

        return create_stream_table_dataframe(
            dict_arcs, time_point=time_point, orient=orient, true_state=true_state
        )

    def visualize(self, model_name, **kwargs):
        """
        Starts up a flask server that serializes the model and pops up a
        webpage with the visualization

        Args:
            model_name : The name of the model that flask will use as an argument
                         for the webpage
        Keyword Args:
            **kwargs: Additional keywords for :func:`idaes.core.ui.fsvis.visualize()`

        Returns:
            None
        """
        visualize(self, model_name, **kwargs)

    def _get_stream_table_contents(self, time_point=0):
        """
        Calls stream_table method and returns result
        """
        return self.stream_table(time_point)

    def _setup_dynamics(self):
        # Look for parent flowsheet
        fs = self.flowsheet()

        # Check the dynamic flag, and retrieve if necessary
        if self.config.dynamic == useDefault:
            if fs is None:
                # No parent, so default to steady-state and warn user
                _log.warning(
                    "{} is a top level flowsheet, but dynamic flag "
                    "set to useDefault. Dynamic "
                    "flag set to False by default".format(self.name)
                )
                self.config.dynamic = False

            else:
                # Get dynamic flag from parent flowsheet
                self.config.dynamic = fs.config.dynamic

        # Check for case when dynamic=True, but parent dynamic=False
        elif self.config.dynamic is True:
            if fs is not None and fs.config.dynamic is False:
                raise DynamicError(
                    "{} trying to declare a dynamic model within "
                    "a steady-state flowsheet. This is not "
                    "supported by the IDAES framework. Try "
                    "creating a dynamic flowsheet instead, and "
                    "declaring some models as steady-state.".format(self.name)
                )

        # Validate units for time domain
        if self.config.time is None and fs is not None:
            # We will get units from parent
            pass
        elif self.config.time_units is None and self.config.dynamic:
            raise ConfigurationError(
                f"{self.name} - no units were specified for the time domain. "
                f"Units must be be specified for dynamic models."
            )
        elif self.config.time_units is None and not self.config.dynamic:
            _log.debug("No units specified for stady-state time domain.")
        elif not isinstance(self.config.time_units, _PyomoUnit):
            raise ConfigurationError(
                "{} unrecognised value for time_units argument. This must be "
                "a Pyomo Unit object (not a compound unit).".format(self.name)
            )

        if self.config.time is not None:
            # Validate user provided time domain
            if self.config.dynamic is True and not isinstance(
                self.config.time, ContinuousSet
            ):
                raise DynamicError(
                    "{} was set as a dynamic flowsheet, but time domain "
                    "provided was not a ContinuousSet.".format(self.name)
                )
            add_object_reference(self, "_time", self.config.time)
            self._time_units = self.config.time_units
        else:
            # If no parent flowsheet, set up time domain
            if fs is None:
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
                            raise DynamicError(
                                "Flowsheet provided with invalid "
                                "time_set attribute - must have at "
                                "least two values (start and end)."
                            )
                    # For dynamics, need a ContinuousSet
                    self._time = ContinuousSet(initialize=self.config.time_set)
                else:
                    # For steady-state, use an ordered Set
                    self._time = pe.Set(initialize=self.config.time_set, ordered=True)
                self._time_units = self.config.time_units

                # Set time config argument as reference to time domain
                self.config.time = self._time
            else:
                # Set time config argument to parent time
                self.config.time = fs.time
                add_object_reference(self, "_time", fs.time)
                self._time_units = fs._time_units
