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
Standard IDAES Feed block.
"""
# Import Pyomo libraries
from pyomo.environ import Reference
from pyomo.common.config import ConfigBlock, ConfigValue, In

# Import IDAES cores
from idaes.core import declare_process_block_class, UnitModelBlockData, useDefault
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.tables import create_stream_table_dataframe
import idaes.logger as idaeslog


__author__ = "Andrew Lee"


# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("Feed")
class FeedData(UnitModelBlockData):
    """
    Standard Feed Block Class
    """

    CONFIG = ConfigBlock()
    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""Indicates whether this model will be dynamic or not,
**default** = False. Feed blocks are always steady-state.""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be False",
            doc="""Feed blocks do not contain holdup, thus this must be False.""",
        ),
    )
    CONFIG.declare(
        "property_package",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for control volume",
            doc="""Property parameter object used to define property calculations,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PhysicalParameterObject** - a PhysicalParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "property_package_args",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing property packages",
            doc="""A ConfigBlock with arguments to be passed to a property block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}""",
        ),
    )

    def build(self):
        """
        Begin building model.

        Args:
            None

        Returns:
            None
        """
        # Call UnitModel.build to setup dynamics
        super(FeedData, self).build()

        # Add State Block
        self.properties = self.config.property_package.build_state_block(
            self.flowsheet().time,
            doc="Material properties in feed",
            defined_state=True,
            has_phase_equilibrium=False,
            **self.config.property_package_args
        )

        # Add references to all state vars
        s_vars = self.properties[self.flowsheet().time.first()].define_state_vars()
        for s in s_vars:
            l_name = s_vars[s].local_name
            if s_vars[s].is_indexed():
                slicer = self.properties[:].component(l_name)[...]
            else:
                slicer = self.properties[:].component(l_name)

            r = Reference(slicer)
            setattr(self, s, r)

        # Add outlet port
        self.add_port(name="outlet", block=self.properties, doc="Outlet Port")

    def initialize_build(
        blk, state_args=None, outlvl=idaeslog.NOTSET, solver=None, optarg=None
    ):
        """
        This method calls the initialization method of the state block.

        Keyword Arguments:
            state_args : a dict of arguments to be passed to the property
                           package(s) to provide an initial state for
                           initialization (see documentation of the specific
                           property package) (default = None).
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None, use
                     default solver options)
            solver : str indicating which solver to use during
                     initialization (default = None, use default solver)

        Returns:
            None
        """
        # ---------------------------------------------------------------------
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")

        if state_args is None:
            state_args = {}

        # Initialize state block
        blk.properties.initialize(
            outlvl=outlvl, optarg=optarg, solver=solver, state_args=state_args
        )

        init_log.info("Initialization Complete.")

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {"Outlet": self.outlet}, time_point=time_point
        )
