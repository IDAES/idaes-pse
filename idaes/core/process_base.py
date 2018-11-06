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
Base for IDAES process model objects.
"""
from __future__ import absolute_import  # disable implicit relative imports
from __future__ import division  # No integer division
from __future__ import print_function  # Python 3 style print

from idaes.core.process_block import declare_process_block_class
from pyomo.core.base.block import _BlockData
from pyomo.environ import Block
from pyomo.common.config import ConfigBlock


# Some more inforation about this module
__author__ = "John Eslick, Qi Chen, Andrew Lee"


__all__ = ['ProcessBlockData']


useDefault = object()


@declare_process_block_class("ProcessBaseBlock")
class ProcessBlockData(_BlockData):
    """
    Base class for most IDAES process models and classes.

    The primary purpose of this class is to create the local config block to
    handle arguments provided by the user when constructing an object and to
    ensure that these arguments are stored in the config block.

    Additionally, this class contains a number of methods common to all IDAES
    classes.
    """
    CONFIG = ConfigBlock("ProcessBlockData", implicit=False)

    def __init__(self, component):
        """
        Initialize the object. Anything inheriting from process base should
        impliment a build() function that creates the Pyomo comonents. This
        build function is called by the blocks default rule.

        Args:
            component: container Block instance to which this _BlockData
                       belongs.

        Returns:
            None
        """
        super(ProcessBlockData, self).__init__(component=component)

    def build(self):
        """
        Default build method for all Classes inheriting from ProcessBlockData.
        Currently empty, but left in place to allow calls to super().build and
        for future compatability.

        Args:
            None

        Returns:
            None
        """
        try:
            idx = self.index()
        except:
            idx = None
        kwargs = self.parent_component()._block_data_config_initialize.get(
            idx, self.parent_component()._block_data_config_default)
        self.config = self.CONFIG(kwargs)

    def fix_initial_conditions(self, state="steady-state"):
        """This method fixes the initial conditions for dynamic models.

        Args:
            state : initial state to use for simulation (default =
                    'steady-state')

        Returns :
            None
        """
        if state == 'steady-state':
            for obj in self.component_objects((Block, Disjunct),
                                              descend_into=True):
                # Try to fix material_accumulation @ first time point
                try:
                    obj.material_accumulation[obj.time.first(), ...].fix(0.0)
                except AttributeError:
                    pass

                # Try to fix energy_accumulation @ first time point
                try:
                    obj.energy_accumulation[obj.time.first(), ...].fix(0.0)
                except AttributeError:
                    pass

        else:
            raise ValueError("Unrecognised value for argument 'state'. "
                             "Valid values are 'steady-state'.")

    def unfix_initial_conditions(self):
        """This method unfixed the initial conditions for dynamic models.

        Args:
            None

        Returns :
            None
        """
        for obj in self.component_objects(Block, descend_into=True):
            # Try to unfix material_accumulation @ first time point
            try:
                obj.material_accumulation[obj.time.first(), ...].unfix()
            except AttributeError:
                pass

            # Try to unfix energy_accumulation @ first time point
            try:
                obj.energy_accumulation[obj.time.first(), ...].unfix()
            except AttributeError:
                pass
