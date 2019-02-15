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
Generic template for a translator block.
"""
from __future__ import division

import logging

# Import Pyomo libraries
from pyomo.common.config import ConfigBlock, ConfigValue, In

# Import IDAES cores
from idaes.core import declare_process_block_class, UnitModelBlockData
from idaes.core.util.config import is_physical_parameter_block

__author__ = "Andrew Lee"


# Set up logger
_log = logging.getLogger(__name__)


@declare_process_block_class("Translator")
class TranslatorData(UnitModelBlockData):
    """
    Standard Translator Block Class
    """
    CONFIG = ConfigBlock()
    CONFIG.declare("dynamic", ConfigValue(
        domain=In([False]),
        default=False,
        description="Dynamic model flag - must be False",
        doc="""Translator blocks are always steady-state."""))
    CONFIG.declare("has_holdup", ConfigValue(
        default=False,
        domain=In([False]),
        description="Holdup construction flag - must be False",
        doc="""Translator blocks do not contain holdup."""))
    CONFIG.declare("inlet_property_package", ConfigValue(
        default=None,
        domain=is_physical_parameter_block,
        description="Property package to use for incoming stream",
        doc="""Property parameter object used to define property calculations
for the incoming stream,
**default** - None.
**Valid values:** {
**PhysicalParameterObject** - a PhysicalParameterBlock object.}"""))
    CONFIG.declare("inlet_property_package_args", ConfigBlock(
        implicit=True,
        description="Arguments to use for constructing property package of "
                    "the incoming stream",
        doc="""A ConfigBlock with arguments to be passed to the property block
associated with the incoming stream,
**default** - None.
**Valid values:** {
see property package for documentation.}"""))
    CONFIG.declare("outlet_property_package", ConfigValue(
        default=None,
        domain=is_physical_parameter_block,
        description="Property package to use for outgoing stream",
        doc="""Property parameter object used to define property calculations
for the outgoing stream,
**default** - None.
**Valid values:** {
**PhysicalParameterObject** - a PhysicalParameterBlock object.}"""))
    CONFIG.declare("outlet_property_package_args", ConfigBlock(
        implicit=True,
        description="Arguments to use for constructing property package of "
                    "the outgoing stream",
        doc="""A ConfigBlock with arguments to be passed to the property block
associated with the outgoing stream,
**default** - None.
**Valid values:** {
see property package for documentation.}"""))

    def build(self):
        """
        Begin building model.

        Args:
            None

        Returns:
            None
        """
        # Call UnitModel.build to setup dynamics
        super(TranslatorData, self).build()

        # Add State Blocks
        self.properties_in = (
                    self.config.inlet_property_package.state_block_class(
                        self.time_ref,
                        doc="Material properties in incoming stream",
                        default={
                            "defined_state": True,
                            "parameters": self.config.inlet_property_package,
                            "has_phase_equilibrium": False,
                            **self.config.inlet_property_package_args}))

        self.properties_out = (
                    self.config.outlet_property_package.state_block_class(
                        self.time_ref,
                        doc="Material properties in outgoing stream",
                        default={
                            "defined_state": True,
                            "parameters": self.config.outlet_property_package,
                            "has_phase_equilibrium": False,
                            **self.config.outlet_property_package_args}))

        # Add outlet port
        self.add_port(name="inlet",
                      block=self.properties_in,
                      doc="Inlet Port")
        self.add_port(name="outlet",
                      block=self.properties_out,
                      doc="Outlet Port")

    def initialize(blk, state_args_in={}, state_args_out={}, outlvl=0,
                   solver='ipopt', optarg={'tol': 1e-6}):
        '''
        This method calls the initialization method of the state blocks.

        Keyword Arguments:
            state_args_in : a dict of arguments to be passed to the inlet
                            property package (to provide an initial state for
                            initialization (see documentation of the specific
                            property package) (default = {}).
            state_args_out : a dict of arguments to be passed to the outlet
                             property package (to provide an initial state for
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
        # ---------------------------------------------------------------------
        # Initialize state block
        flags = blk.properties_in.initialize(outlvl=outlvl-1,
                                             optarg=optarg,
                                             solver=solver,
                                             **state_args_in,
                                             hold_state=True)

        blk.properties_out.initialize(outlvl=outlvl-1,
                                      optarg=optarg,
                                      solver=solver,
                                      **state_args_out)

        blk.properties_in.release_state(flags)

        if outlvl > 0:
            _log.info('{} Initialisation Complete.'.format(blk.name))
