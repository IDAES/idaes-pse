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
General purpose mixer block for IDAES models
"""
from __future__ import absolute_import  # disable implicit relative imports
from __future__ import division, print_function

import logging

from pyomo.common.config import ConfigBlock, ConfigValue, In

from idaes.core import declare_process_block_class, UnitBlockData, useDefault
from idaes.core.util.config import is_physical_parameter_block, list_of_strings

__author__ = "Andrew Lee"


# Set up logger
_log = logging.getLogger(__name__)


@declare_process_block_class("MixerBlock")
class MixerBlockData(UnitBlockData):
    """
    This is a general purpose model for a Mixer block with the IDAES modeling
    framework. This block can be used either as a stand-alone Mixer unit
    operation, or as a sub-model within another unit operation.

    This model creates a number of StateBlocks to represent the incoming
    streams, then writes a set of phase-component material balances, an
    overall enthalpy balance and a momentum balance (2 options) linked to a
    mixed-state StateBlock. The mixed-state StateBlock can either be specified
    by the user (allowing use as a sub-model), or created by the MixerBlock.

    When being used as a sub-model, MixerBlock should only be used when a set
    of new StateBlocks are required for the streams to be mixed. It should not
    be used to mix streams from mutiple ControlVolumes in a single unit model -
    in these cases the unit model developer should write their own mixing
    equations.
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
    CONFIG.declare("inlet_list", ConfigValue(
        domain=list_of_strings,
        description="List of inlet names",
        doc="""A list containing names of inlets, **default** - None.
**Valid values:** {**None** - use num_inlets argument, **list** - a list of
names to use for inlets.}"""))
    CONFIG.declare("num_inlets", ConfigValue(
        domain=int,
        description="Number of inlets to unit",
        doc="""Argument indication number (int) of inlets to construct, not
used if inlet_list arg is provided, **default** - None. **Valid values:** {
**None** - use inlet_list arg instead, or default to 2 if neither argument
provided, **int** - number of inlets to creat (will be named with sequential
integers from 1 to num_inlets).}"""))
    CONFIG.declare("calculate_phase_equilibrium", ConfigValue(
        default=False,
        domain=In(True, False),
        description="Calculate phase equilibrium in mixed stream",
        doc="""Argument indicating whether phase equilibrium should be
calculated for the resulting mixed stream, **default** - False. **Valid values:
** {**True** - calculate phase equilibrium in mixed stream, **False** -
do not calculate equilibrium in mixed stream.}"""))

    def build(self):
        """
        General build method for MixerBlockData. This method calls a number
        of sub-methods which automate the construction of expected attributes
        of unit models.

        Inheriting models should call `super().build`.

        Args:
            None

        Returns:
            None
        """
        super(MixerBlockData, self).build()
