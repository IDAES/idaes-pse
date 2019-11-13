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
    TerminationCondition, value

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
#     CONFIG.declare("has_feed", ConfigValue(
#         default=False,
#         domain=In([True, False]),
#         description="Feed inlet construction flag.",
#         doc="""Indicates if there is a feed inlet to the tray,
# **default** - False.
# **Valid values:** {
# **True** - include a feed inlet to the tray,
# **False** - exclude a feed inlet to the tray.}"""))
    CONFIG.declare("has_side_liquid_draw", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="Side liquid draw construction flag.",
        doc="""Indicates if there is a side liquid draw from the tray,
**default** - False.
**Valid values:** {
**True** - include a side liquid draw from the tray,
**False** - exclude a side liquid draw from the tray.}"""))
    CONFIG.declare("has_side_vapor_draw", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="Side vapor draw construction flag.",
        doc="""Indicates if there is a side vapor draw from the tray,
**default** - False.
**Valid values:** {
**True** - include a side vapor draw from the tray,
**False** - exclude a side vapor draw from the tray.}"""))
    CONFIG.declare("has_heat_transfer", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="Heat transfer to/from tray construction flag.",
        doc="""Indicates if there is heat transfer to/from the tray,
**default** - False.
**Valid values:** {
**True** - include a heat transfer term,
**False** - exclude a heat transfer term.}"""))
    CONFIG.declare("has_pressure_change", ConfigValue(
        default=True,
        domain=In([True, False]),
        description="Pressure change term construction flag",
        doc="""Indicates whether terms for pressure change should be
    constructed,
    **default** - False.
    **Valid values:** {
    **True** - include pressure change terms,
    **False** - exclude pressure change terms.}"""))
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

        self._add_pressure_balance()
        self._add_ports()

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

        if self.config.has_heat_transfer:
            self.heat_duty = Var(initialize=0,
                                 doc="heat duty for the tray")

        @self.Constraint(self.flowsheet().config.time, doc="Energy balances")
        def enthalpy_mixing_equations(b, t):
            return 0 == (
                sum(self.liq_in_properties[t].get_enthalpy_flow_terms(p)
                    for p in b.config.property_package.phase_list) +
                sum(self.vap_in_properties[t].get_enthalpy_flow_terms(p)
                    for p in b.config.property_package.phase_list) -
                sum(self.mixed_feed_properties[t].get_enthalpy_flow_terms(p)
                    for p in b.config.property_package.phase_list))

    def _add_pressure_balance(self):
        if self.config.has_pressure_change:
            self.deltaP = Var(initialize=0,
                              doc="pressure drop across tray")
        @self.Constraint(self.flowsheet().config.time,
                         doc="Pressure drop constraint for tray")
        def pressure_drop_equation(self, t):
            if self.config.has_pressure_change:
                return self.mixed_feed_properties[t].pressure == \
                    self.vap_in_properties[t].pressure - self.deltaP
            else:
                return self.mixed_feed_properties[t].pressure == \
                    self.vap_in_properties[t].pressure

    def _add_ports(self):

        if self.config.has_side_liquid_draw:
            self.side_liq = Port(noruleinit=True, doc="side liquid draw.")
        if self.config.has_side_vapor_draw:
            self.side_vap = Port(noruleinit=True, doc="side vapor draw.")

    def initialize(self, solver=None, outlvl=None):

        # Initialize the inlet state blocks
        self.liq_in_properties.initialize(outlvl=outlvl)
        self.vap_in_properties.initialize(outlvl=outlvl)

        # Initialize the mixed outlet state block
        self.mixed_feed_properties.initialize(outlvl=outlvl)

        # Deactivate energy balance
        self.enthalpy_mixing_equations.deactivate()
        average_temperature = \
            0.5 * (self.liq_in_properties[0].temperature.value
                   + self.vap_in_properties[0].temperature.value)

        self.mixed_feed_properties[:].temperature.fix(average_temperature)

        # Deactivate pressure balance
        self.pressure_drop_equation.deactivate()
        self.mixed_feed_properties[:].pressure.\
            fix(self.vap_in_properties[0].pressure.value)

        if solver is not None:
            if outlvl > 2:
                tee = True
            else:
                tee = False

        solver_output = solver.solve(self, tee=tee)

        if solver_output.solver.termination_condition == \
                TerminationCondition.optimal:
            _log.info('{} Mass balance solve successful.'
                      .format(self.name))
        else:
            _log.info('{} Mass balance solve failed.'
                      .format(self.name))

        # Activate energy balance
        self.enthalpy_mixing_equations.activate()
        self.mixed_feed_properties[:].temperature.unfix()

        solver_output = solver.solve(self, tee=tee)

        if solver_output.solver.termination_condition == \
                TerminationCondition.optimal:
            _log.info('{} Mass/Energy balance solve successful.'
                      .format(self.name))
        else:
            _log.info('{} Mass/Energy balance solve failed.'
                      .format(self.name))

        # Activate pressure balance
        self.pressure_drop_equation.activate()
        self.mixed_feed_properties[:].pressure.unfix()

        if solver_output.solver.termination_condition == \
                TerminationCondition.optimal:
            _log.info('{} Tray initialisation complete.'
                      .format(self.name))
        else:
            _log.info('{} Mass/Energy/Pressure balance solve failed.'
                      .format(self.name))
