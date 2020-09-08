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
from pyomo.network import Arc
from pyomo.environ import Reference, Expression, Var, Constraint, \
    TerminationCondition, value, Integers, RangeSet, TransformationFactory

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
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import propagate_state, \
    solve_indexed_blocks

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
        self.tray_index = RangeSet(1, self.config.number_of_trays)

        # Add trays
        self.tray = Tray(self.tray_index,
                         default={"property_package":
                                  self.config.property_package,
                                  "has_heat_transfer":
                                      self.config.has_heat_transfer,
                                  "has_pressure_change":
                                      self.config.has_pressure_change},
                         initialize={self.config.feed_tray_location:
                                     {"property_package":
                                         self.config.property_package,
                                      "is_feed_tray": True,
                                      "has_heat_transfer":
                                          self.config.has_heat_transfer,
                                      "has_pressure_change":
                                          self.config.has_pressure_change}})

        # Add condenser
        self.condenser = Condenser(
            default={"property_package": self.config.property_package,
                     "property_package_args":
                     self.config.property_package_args,
                     "condenser_type": self.config.condenser_type,
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

        self._make_pressure_balance()

        self._make_arcs()
        TransformationFactory("network.expand_arcs").apply_to(self)

    def _make_pressure_balance(self):
        if self.config.has_pressure_change:
            self.deltaP = Var(self.flowsheet().config.time,
                              self.tray_index,
                              initialize=0,
                              doc="pressure drop across tray")

        @self.Constraint(self.flowsheet().config.time,
                         self.tray_index,
                         doc="pressure balance for tray")
        def pressure_drop_equation(self, t, i):
            if self.config.has_pressure_change:
                if i == 1:
                    return self.tray[i].properties_out[t].pressure == \
                        self.condenser.control_volume.properties_out[t].\
                        pressure + self.deltaP[t, i]
                else:
                    return self.tray[i].properties_out[t].pressure == \
                        self.tray[i - 1].properties_out[t].\
                        pressure + self.deltaP[t, i]
            else:
                if i == 1:
                    return self.tray[i].properties_out[t].pressure == \
                        self.condenser.control_volume.properties_out[t].\
                        pressure
                else:
                    return self.tray[i].properties_out[t].pressure == \
                        self.tray[i - 1].properties_out[t].\
                        pressure            

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

    def my_propagate_state(self, source=None,
                           destination=None, direction="forward"):
        """
        This method propagates values between Ports along Arcs. Values can be
        propagated in either direction using the direction argument.

        Args:
            stream : Arc object along which to propagate values
            direction: direction in which to propagate values. Default = 'forward'
                    Valid value: 'forward', 'backward'.
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
                             " being used for initialization.")
            solver = get_default_solver()

        self.tray[self.config.feed_tray_location].initialize()

        self.my_propagate_state(
            source=self.tray[self.config.feed_tray_location].vap_out,
            destination=self.condenser.inlet)

        self.condenser.initialize()

        self.my_propagate_state(
            source=self.tray[self.config.feed_tray_location].liq_out,
            destination=self.reboiler.inlet)

        self.reboiler.initialize()

        for i in self.tray_index:
            if i < self.config.feed_tray_location:
                self.my_propagate_state(
                    source=self.condenser.reflux,
                    destination=self.tray[i].liq_in)
                self.my_propagate_state(
                    source=self.tray[self.config.feed_tray_location].vap_out,
                    destination=self.tray[i].vap_in)
                self.tray[i].initialize()

            elif i > self.config.feed_tray_location:
                self.my_propagate_state(
                    source=self.tray[self.config.feed_tray_location].liq_out,
                    destination=self.tray[i].liq_in)
                self.my_propagate_state(
                    source=self.reboiler.vapor_reboil,
                    destination=self.tray[i].vap_in)
                self.tray[i].initialize()

        # self.reboiler.control_volume.properties_in[0].flow_mol.fix(100)
        # self.reboiler.control_volume.properties_in[0].temperature.fix(368)
        # self.reboiler.control_volume.properties_in[0].pressure.fix(101325)
        # self.reboiler.control_volume.properties_in[0].mole_frac_comp["benzene"].fix(0.5)
        # self.reboiler.control_volume.properties_in[0].mole_frac_comp["toluene"].fix(0.5)

        # self.vap_stream.display()
        print(degrees_of_freedom(self))

        self.vap_stream[self.config.number_of_trays + 1].deactivate()
        # self.vap_stream[1].deactivate()

        # self.liq_stream[0].deactivate()
        self.liq_stream[3].deactivate()

        self.reboiler.deactivate()
        # self.condenser.deactivate()

        # self.tray[1].properties_in_liq[0].flow_mol.fix()
        # self.tray[1].properties_in_liq[0].temperature.fix()
        # self.tray[1].properties_in_liq[0].pressure.fix()
        # self.tray[1].properties_in_liq[0].mole_frac_comp["benzene"].fix()
        # self.tray[1].properties_in_liq[0].mole_frac_comp["toluene"].fix()

        self.tray[3].properties_in_vap[0].flow_mol.fix()
        self.tray[3].properties_in_vap[0].temperature.fix()
        self.tray[3].properties_in_vap[0].pressure.fix()
        self.tray[3].properties_in_vap[0].mole_frac_comp["benzene"].fix()
        self.tray[3].properties_in_vap[0].mole_frac_comp["toluene"].fix()


        # raise Exception(degrees_of_freedom(self))

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = solver.solve(self, tee=slc.tee)
        init_log.info_high(
            "Column initialization status {}.".format(idaeslog.condition(res))
        )

        print(res)

        # self.tray[3].properties_in_vap[0].pressure.unfix()
        # self.eq_pressure = Constraint(
        #     expr=self.reboiler.control_volume.properties_out[0].pressure ==
        #     self.tray[3].properties_in_vap[0].pressure)

        # self.tray[3].properties_in_vap[0].temperature.unfix()
        # self.eq_temperature = Constraint(
        #     expr=self.reboiler.control_volume.properties_out[0].temperature ==
        #     self.tray[3].properties_in_vap[0].temperature)


        # with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
        #     res = solver.solve(self, tee=slc.tee)
        # init_log.info_high(
        #     "Column initialization status {}.".format(idaeslog.condition(res))
        # )

        for i in self.tray_index:
            self.tray[i].liq_in.display()
            self.tray[i].vap_in.display()

            self.tray[i].liq_out.display()
            self.tray[i].vap_out.display()

        self.reboiler.report()
        self.condenser.report()
        # self.tray[3].properties_in_vap[0].display()
        # self.reboiler.control_volume.properties_out[0].display()
        raise Exception(res)
