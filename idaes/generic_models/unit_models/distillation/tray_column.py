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
from pyomo.network import Arc, Port
from pyomo.environ import value, Integers, RangeSet, TransformationFactory

# Import IDAES cores
from idaes.generic_models.unit_models.distillation import Tray, Condenser, \
    Reboiler
from idaes.generic_models.unit_models.distillation.condenser \
    import CondenserType, TemperatureSpec
from idaes.core import (declare_process_block_class,
                        UnitModelBlockData,
                        useDefault)
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.config import is_physical_parameter_block
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
**Valid values:** {
Must be integer number.}"""))
    CONFIG.declare("feed_tray_location", ConfigValue(
        default=None,
        domain=In(Integers),
        description="Feed tray location in the column",
        doc="""Indicates the number of trays to be constructed.
**default** - None.
**Valid values:** {
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
    CONFIG.declare("condenser_temperature_spec", ConfigValue(
        default=None,
        domain=In(TemperatureSpec),
        description="Temperature spec for the condenser",
        doc="""Temperature specification for the condenser,
**default** - TemperatureSpec.none
**Valid values:** {
**TemperatureSpec.none** - No spec is selected,
**TemperatureSpec.atBubblePoint** - Condenser temperature set at
bubble point i.e. total condenser,
**TemperatureSpec.customTemperature** - Condenser temperature at
user specified temperature.}"""))
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
    CONFIG.declare("has_liquid_side_draw", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="liquid side draw construction flag.",
        doc="""indicates if there is a liquid side draw from all trays,
**default** - False.
**Valid values:** {
**True** - include a liquid side draw for all trays,
**False** - exclude a liquid side draw for all trays.}"""))
    CONFIG.declare("has_vapor_side_draw", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="vapor side draw construction flag.",
        doc="""indicates if there is a vapor side draw from all trays,
**default** - False.
**Valid values:** {
**True** - include a vapor side draw for all trays,
**False** - exclude a vapor side draw for all trays.}"""))
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
        if self.config.number_of_trays is not None:
            self.tray_index = RangeSet(1, self.config.number_of_trays)
        else:
            raise ConfigurationError("The config argument number_of_trays "
                                     "needs to specified and cannot be None.")

        # Add trays

        # TODO:
        # 1. Add support for multiple feed tray locations
        # 2. Add support for specifying which trays have side draws
        self.tray = Tray(self.tray_index,
                         default={"has_liquid_side_draw":
                                  self.config.has_liquid_side_draw,
                                  "has_vapor_side_draw":
                                  self.config.has_vapor_side_draw,
                                  "has_heat_transfer":
                                  self.config.has_heat_transfer,
                                  "has_pressure_change":
                                  self.config.has_pressure_change,
                                  "property_package":
                                  self.config.property_package,
                                  "property_package_args":
                                  self.config.property_package_args},
                         initialize={self.config.feed_tray_location:
                                     {"is_feed_tray": True,
                                      "has_liquid_side_draw":
                                      self.config.has_liquid_side_draw,
                                      "has_vapor_side_draw":
                                      self.config.has_vapor_side_draw,
                                      "has_heat_transfer":
                                      self.config.has_heat_transfer,
                                      "has_pressure_change":
                                      self.config.has_pressure_change,
                                      "property_package":
                                      self.config.property_package,
                                      "property_package_args":
                                      self.config.property_package_args}})

        # Add condenser
        self.condenser = Condenser(
            default={"condenser_type": self.config.condenser_type,
                     "temperature_spec":
                         self.config.condenser_temperature_spec,
                     "property_package": self.config.property_package,
                     "property_package_args":
                     self.config.property_package_args})

        # Add reboiler
        self.reboiler = Reboiler(
            default={"has_boilup_ratio": True,
                     "has_pressure_change": self.config.has_pressure_change,
                     "property_package": self.config.property_package,
                     "property_package_args":
                     self.config.property_package_args})

        # Add extension to the feed port
        self.feed = Port(extends=self.tray[self.config.feed_tray_location].
                         feed)

        # Construct arcs between trays, condenser, and reboiler
        self._make_arcs()
        TransformationFactory("network.expand_arcs").apply_to(self)

    def _make_arcs(self):
        # make arcs
        self.liq_stream_index = RangeSet(0, self.config.number_of_trays)
        self.vap_stream_index = RangeSet(1, self.config.number_of_trays + 1)

        def rule_liq_stream(self, i):
            if i == 0:
                return {"source": self.condenser.reflux,
                        "destination": self.tray[i + 1].liq_in}
            elif i == self.config.number_of_trays:
                return {"source": self.tray[i].liq_out,
                        "destination": self.reboiler.inlet}
            else:
                return {"source": self.tray[i].liq_out,
                        "destination": self.tray[i + 1].liq_in}

        def rule_vap_stream(self, i):
            if i == 1:
                return {"source": self.tray[i].vap_out,
                        "destination": self.condenser.inlet}
            elif i == self.config.number_of_trays + 1:
                return {"source": self.reboiler.vapor_reboil,
                        "destination": self.tray[i - 1].vap_in}
            else:
                return {"source": self.tray[i].vap_out,
                        "destination": self.tray[i - 1].vap_in}

        self.liq_stream = Arc(self.liq_stream_index, rule=rule_liq_stream)
        self.vap_stream = Arc(self.vap_stream_index, rule=rule_vap_stream)

    def propagate_stream_state(self, source=None,
                               destination=None):
        """
        This method is used during initialization to propage state values
        between any two ports.

        Args:
            source : source port
            destination : destination port
        """
        for v in source.vars:
            for i in destination.vars[v]:
                if not destination.vars[v][i].fixed:
                    destination.vars[v][i].value = value(source.vars[v][i])

    def initialize(self, state_args_feed=None, state_args_liq=None,
                   state_args_vap=None, solver=None, outlvl=idaeslog.NOTSET):

        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        init_log.info("Begin initialization.")

        if solver is None:
            init_log.warning("Solver not provided. Default solver(ipopt) "
                             "being used for initialization.")
            solver = get_default_solver()

        feed_flags = self.tray[self.config.feed_tray_location].initialize()

        self.propagate_stream_state(
            source=self.tray[self.config.feed_tray_location].vap_out,
            destination=self.condenser.inlet)

        self.condenser.initialize()

        self.propagate_stream_state(
            source=self.tray[self.config.feed_tray_location].liq_out,
            destination=self.reboiler.inlet)

        self.reboiler.initialize()

        for i in self.tray_index:
            self.propagate_stream_state(
                source=self.condenser.reflux,
                destination=self.tray[i].liq_in)
            self.propagate_stream_state(
                source=self.reboiler.vapor_reboil,
                destination=self.tray[i].vap_in)
            self.tray[i].initialize()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = solver.solve(self, tee=slc.tee)
        init_log.info(
            "Column initialization status {}.".format(idaeslog.condition(res))
        )

        # release feed tray state once initialization is complete
        self.tray[self.config.feed_tray_location].properties_in_feed.\
            release_state(flags=feed_flags, outlvl=outlvl)
