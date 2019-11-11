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
Tray Mixer model for distillation.
"""

__author__ = "Jaffer Ghouse"

import logging

# Import Pyomo libraries
from pyomo.common.config import ConfigBlock, ConfigValue, In
from pyomo.network import Port
from pyomo.environ import Reference, Expression, Var, Constraint, \
    TerminationCondition

# Import IDAES cores
from idaes.core import (ControlVolume0DBlock,
                        declare_process_block_class,
                        EnergyBalanceType,
                        MomentumBalanceType,
                        MaterialBalanceType,
                        UnitModelBlockData,
                        useDefault)
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.misc import add_object_reference
from idaes.core.util.exceptions import PropertyPackageError

_log = logging.getLogger(__name__)


@declare_process_block_class("TrayMixer")
class TrayMixerData(UnitModelBlockData):
    """
    Tray Mixer unit for distillation model.
    """
    CONFIG = UnitModelBlockData.CONFIG()
    CONFIG.declare("has_feed", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="Feed inlet construction flag.",
        doc="""Indicates if there is a feed inlet to the tray,
**default** - False.
**Valid values:** {
**True** - include a feed inlet to the tray,
**False** - exclude a feed inlet to the tray.}"""))
    CONFIG.declare("property_package", ConfigValue(
        default=useDefault,
        domain=is_physical_parameter_block,
        description="Property package to use for control volume",
        doc="""Property parameter object used to define property calculations,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PropertyParameterObject** - a PropertyParameterBlock object.}"""))
    CONFIG.declare("property_package_args", ConfigBlock(
        implicit=True,
        description="Arguments to use for constructing property packages",
        doc="""A ConfigBlock with arguments to be passed to a property block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}"""))

    def build(self):
        """Build the model.

        Args:
            None
        Returns:
            None
        """
        # Call UnitModel.build to setup dynamics
        super(TrayMixerData, self).build()

        # Create the inlets for the tray mixer
        inlet_list = ["liq", "vap"]

        # add inlet state blocks

        # Setup StateBlock argument dict
        state_block_args = dict(**self.config.property_package_args)
        state_block_args["has_phase_equilibrium"] = False
        state_block_args["parameters"] = self.config.property_package
        state_block_args["defined_state"] = True

        for i in inlet_list:
            state_obj = self.config.property_package.state_block_class(
                self.flowsheet().config.time,
                doc="State block for inlet to tray",
                default=state_block_args)

            setattr(self, i + "_in_properties", state_obj)

        # add mixed outlet state blocks which is the feed to the tray

        # Setup StateBlock argument dict
        mixed_block_args = dict(**self.config.property_package_args)
        mixed_block_args["has_phase_equilibrium"] = True
        mixed_block_args["parameters"] = self.config.property_package
        mixed_block_args["defined_state"] = False

        self.mixed_feed_properties = self.config.property_package.\
            state_block_class(self.flowsheet().config.time,
                              doc="State block for inlet to tray",
                              default=mixed_block_args)

        self._add_material_balance()
        self._add_energy_balance()

    def _add_material_balance(self):

        @self.Constraint(self.flowsheet().config.time,
                         self.config.property_package.component_list,
                         doc="Material mixing equations")
        def material_mixing_equations(b, t, j):
            return 0 == sum(
                self.liq_in_properties[t].get_material_flow_terms(p, j) +
                self.vap_in_properties[t].get_material_flow_terms(p, j) -
                self.mixed_feed_properties[t].get_material_flow_terms(p, j)
                for p in b.config.property_package.phase_list)

    def _add_energy_balance(self):

        @self.Constraint(self.flowsheet().config.time, doc="Energy balances")
        def enthalpy_mixing_equations(b, t):
            return 0 == (
                sum(self.liq_in_properties[t].get_enthalpy_flow_terms(p)
                    for p in b.config.property_package.phase_list) +
                sum(self.vap_in_properties[t].get_enthalpy_flow_terms(p)
                    for p in b.config.property_package.phase_list) -
                sum(self.mixed_feed_properties[t].get_enthalpy_flow_terms(p)
                    for p in b.config.property_package.phase_list))

    def initialize(self, solver=None, outlvl=None):

        # Initialize the inlet state blocks
        self.liq_in_properties.initialize(outlvl=outlvl)
        self.vap_in_properties.initialize(outlvl=outlvl)

        if solver is not None:
            if outlvl > 2:
                tee = True
            else:
                tee = False

            solver_output = solver.solve(self, tee=tee)
            if solver_output.solver.termination_condition == \
                    TerminationCondition.optimal:
                _log.info('{} Tray Mixer Initialisation Complete.'
                          .format(self.name))
