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
Tray column model for distillation.
"""

__author__ = "Jaffer Ghouse"

import idaes.logger as idaeslog

# Import Pyomo libraries
from pyomo.common.config import ConfigBlock, ConfigValue, In
from pyomo.network import Port
from pyomo.environ import Reference, Expression, Var, Constraint, \
    TerminationCondition, value, Integers, RangeSet

# Import IDAES cores
from idaes.generic_models.unit_models.distillation import Tray, Condenser, \
    Reboiler
from idaes.generic_models.unit_models.distillation.condenser \
    import CondenserType, TemperatureSpec
from idaes.core import (ControlVolume0DBlock,
                        declare_process_block_class,
                        EnergyBalanceType,
                        MomentumBalanceType,
                        MaterialBalanceType,
                        UnitModelBlockData,
                        useDefault)
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import ConfigurationError, \
    PropertyPackageError, PropertyNotSupportedError
from idaes.core.util.testing import get_default_solver

_log = idaeslog.getLogger(__name__)


@declare_process_block_class("TrayColumn")
class TrayColumnData(UnitModelBlockData):
    """
    Tray Column model.
    """
    CONFIG = ConfigBlock()
    CONFIG.declare("dynamic", ConfigValue(
        domain=In([False]),
        default=False,
        description="Dynamic model flag - must be False",
        doc="""Indicates whether this model will be dynamic or not,
**default** = False. Tray column units do not support dynamic behavior."""))
    CONFIG.declare("has_holdup", ConfigValue(
        default=False,
        domain=In([False]),
        description="Holdup construction flag - must be False",
        doc="""Indicates whether holdup terms should be constructed or not.
**default** - False. Tray column units do not have defined volume, thus
this must be False."""))
    CONFIG.declare("number_of_trays", ConfigValue(
        default=None,
        domain=In(Integers),
        description="Number of trays in the column",
        doc="""Indicates the number of trays to be constructed.
**default** - None.
**Valid values:**{
Must be integer number.}"""))
    CONFIG.declare("feed_tray_location", ConfigValue(
        default=None,
        domain=In(Integers),
        description="Feed tray location in the column",
        doc="""Indicates the number of trays to be constructed.
**default** - None.
**Valid values:**{
Must be integer number.}"""))
    CONFIG.declare("condenser_type", ConfigValue(
        default=CondenserType.totalCondenser,
        domain=In(CondenserType),
        description="Type of condenser flag",
        doc="""Indicates what type of condenser should be constructed,
**default** - CondenserType.totalCondenser.
**Valid values:** {
**CondenserType.totalCondenser** - Incoming vapor from top tray is condensed
to all liquid,
**CondenserType.partialCondenser** - Incoming vapor from top tray is
partially condensed to a vapor and liquid stream.}"""))
    CONFIG.declare("has_heat_transfer", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="heat duty to/from tray construction flag.",
        doc="""indicates if there is heat duty to/from the tray,
**default** - False.
**Valid values:** {
**True** - include a heat duty term,
**False** - exclude a heat duty term.}"""))
    CONFIG.declare("has_pressure_change", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="pressure change term construction flag",
        doc="""indicates whether terms for pressure change should be
constructed,
**default** - False.
**Valid values:** {
**True** - include pressure change terms,
**False** - exclude pressure change terms.}"""))
    CONFIG.declare("property_package", ConfigValue(
        default=useDefault,
        domain=is_physical_parameter_block,
        description="property package to use for control volume",
        doc="""property parameter object used to define property calculations,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PropertyParameterObject** - a PropertyParameterBlock object.}"""))
    CONFIG.declare("property_package_args", ConfigBlock(
        implicit=True,
        description="arguments to use for constructing property packages",
        doc="""a ConfigBlock with arguments to be passed to a property block(s)
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
        super(TrayColumnData, self).build()

        # Create set for constructing indexed trays
        tray_index = RangeSet(1, self.config.number_of_trays)

        # Add trays
        self.tray = Tray(tray_index, default={})
        
        # Add condenser
        self.condenser = Condenser(
            default={"property_package": self.config.property_package,
                     "property_package_args":
                     self.config.property_package_args,
                     "condenser_type": CondenserType.totalCondenser,
                     "temperature_spec": TemperatureSpec.atBubblePoint,
                     "has_pressure_change": self.config.has_pressure_change})

        # Add reboiler
        self.reboiler = Reboiler(
            default={"property_package": self.config.property_package,
                     "property_package_args":
                     self.config.property_package_args,
                     "has_boilup_ratio": True,
                     "has_pressure_change":
                     self.config.has_pressure_change})

        raise Exception(tray_index)


    def initialize(self, state_args_feed=None, state_args_liq=None,
                   state_args_vap=None, solver=None, outlvl=idaeslog.NOTSET):

        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        init_log.info("Begin initialization.")

        if solver is None:
            init_log.warning("Solver not provided. Default solver(ipopt) "
                             " being used for initialization.")
            solver = get_default_solver()

        if self.config.has_liquid_side_draw:
            if not self.liq_side_sf.fixed:
                raise ConfigurationError(
                    "Liquid side draw split fraction not fixed but "
                    "has_liquid_side_draw set to True.")

        if self.config.has_vapor_side_draw:
            if not self.vap_side_sf.fixed:
                raise ConfigurationError(
                    "Vapor side draw split fraction not fixed but "
                    "has_vapor_side_draw set to True.")

        # Initialize the inlet state blocks
        if self.config.is_feed_tray:
            self.properties_in_feed.initialize(state_args=state_args_feed,
                                               solver=solver,
                                               outlvl=outlvl)
        self.properties_in_liq.initialize(state_args=state_args_liq,
                                          solver=solver,
                                          outlvl=outlvl)
        self.properties_in_vap.initialize(state_args=state_args_vap,
                                          solver=solver,
                                          outlvl=outlvl)

        # Deactivate energy balance
        self.enthalpy_mixing_equations.deactivate()

        # state args to initialize the mixed outlet state block
        if self.config.is_feed_tray and state_args_feed is not None:
            # if initial guess provided for the feed stream, initialize the
            # mixed state block at the same condition.
            self.properties_out.initialize(state_args=state_args_feed,
                                           solver=solver,
                                           outlvl=outlvl)
        else:
            state_args_mixed = {}
            state_dict = \
                self.properties_out[self.flowsheet().config.time.first()].\
                define_port_members()
            if self.config.is_feed_tray:
                for k in state_dict.keys():
                    if state_dict[k].is_indexed():
                        state_args_mixed[k] = {}
                        for m in state_dict[k].keys():
                            state_args_mixed[k][m] = \
                                self.properties_in_feed[0].\
                                component(state_dict[k].local_name)[m].value
                    else:
                        state_args_mixed[k] = \
                            self.properties_in_feed[0].\
                            component(state_dict[k].local_name).value

            else:
                # if not feed tray, initialize mixed state block at average of
                # vap/liq inlets except pressure. While this is crude, it
                # will work for most combination of state vars.
                for k in state_dict.keys():
                    if "pressure" in k:
                        # Take the lowest pressure and this is the liq inlet
                        state_args_mixed[k] = self.properties_in_liq[0].\
                            component(state_dict[k].local_name).value
                    elif state_dict[k].is_indexed():
                        state_args_mixed[k] = {}
                        for m in state_dict[k].keys():
                            state_args_mixed[k][m] = \
                                0.5 * (self.properties_in_liq[0].
                                       component(state_dict[k].local_name)[m].
                                       value + self.properties_in_vap[0].
                                       component(state_dict[k].local_name)[m].
                                       value)
                    else:
                        state_args_mixed[k] = \
                            0.5 * (self.properties_in_liq[0].
                                   component(state_dict[k].local_name).value +
                                   self.properties_in_vap[0].
                                   component(state_dict[k].local_name).value)

        # Deactivate pressure balance
        self.pressure_drop_equation.deactivate()

        self.properties_out.initialize(state_args=state_args_mixed,
                                       solver=solver,
                                       outlvl=outlvl)

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = solver.solve(self, tee=slc.tee)
        init_log.info_high(
            "Mass balance solve {}.".format(idaeslog.condition(res))
        )

        # Activate energy balance
        self.enthalpy_mixing_equations.activate()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = solver.solve(self, tee=slc.tee)
        init_log.info_high(
            "Mass and energy balance solve {}.".format(idaeslog.condition(res))
        )

        # Activate pressure balance
        self.pressure_drop_equation.activate()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = solver.solve(self, tee=slc.tee)
        init_log.info_high(
            "Mass, energy and pressure balance solve {}.".
            format(idaeslog.condition(res)))
        init_log.info(
            "Initialization complete, status {}.".
            format(idaeslog.condition(res)))
