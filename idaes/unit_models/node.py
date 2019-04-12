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
Standard IDAES Node model.
"""
from __future__ import division
import logging

# Import Pyomo libraries
from pyomo.environ import SolverFactory
from pyomo.common.config import ConfigBlock, ConfigValue, In

# Import IDAES cores
from idaes.core import (declare_process_block_class,
                        UnitModelBlockData,
                        useDefault)
from idaes.core.util.config import is_physical_parameter_block

__author__ = "Andrew Lee"

# Set up logger
_log = logging.getLogger(__name__)


@declare_process_block_class("Node")
class NodeData(UnitModelBlockData):
    """
    Standard Node Unit Model Class
    """
    CONFIG = ConfigBlock()
    CONFIG.declare("dynamic", ConfigValue(
        domain=In([False]),
        default=False,
        description="Dynamic model flag - must be False",
        doc="""Indicates whether this model will be dynamic or not,
**default** = False."""))
    CONFIG.declare("has_holdup", ConfigValue(
        default=False,
        domain=In([False]),
        description="Holdup construction flag - must be False",
        doc="""Indicates whether holdup terms should be constructed or not.
**default** - False. Nodes do not have defined volume, thus
this must be False."""))
    CONFIG.declare("property_package", ConfigValue(
        default=useDefault,
        domain=is_physical_parameter_block,
        description="Property package to use for control volume",
        doc="""Property parameter object used to define property calculations,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PhysicalParameterObject** - a PhysicalParameterBlock object.}"""))
    CONFIG.declare("property_package_args", ConfigBlock(
        implicit=True,
        description="Arguments to use for constructing property packages",
        doc="""A ConfigBlock with arguments to be passed to a property block(s)
and used when constructing these,
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
        super(NodeData, self).build()

        self.properties = \
                self.config.property_package.state_block_class(
                        self.time_ref,
                        doc="Material properties",
                        default={"has_phase_equilibrium": False,
                                 "parameters": self.config.property_package,
                                 "defined_state": True,
                                 **self.config.property_package_args})

        # Add Ports
        self.add_inlet_port(name="inlet",
                            block=self.properties,
                            doc="Inlet block")
        self.add_outlet_port(name="outlet",
                             block=self.properties,
                             doc="Outlet block")

    def initialize(blk, state_args={}, outlvl=0,
                   solver='ipopt', optarg={'tol': 1e-6}):
        '''
        This method initializes the Node block by calling the initialize method
        on the properties block.

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
        opt = SolverFactory(solver)
        opt.options = optarg

        # ---------------------------------------------------------------------
        # Initialize control volume block
        blk.properties.initialize(outlvl=outlvl-1,
                                  optarg=optarg,
                                  solver=solver,
                                  hold_state=False,
                                  **state_args)

        if outlvl > 0:
            _log.info('{} Initialisation Step Complete.'.format(blk.name))
