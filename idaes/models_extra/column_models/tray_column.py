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
Tray column model for distillation.
"""

__author__ = "Jaffer Ghouse"

import idaes.logger as idaeslog

# Import Pyomo libraries
from pyomo.common.config import ConfigBlock, ConfigValue, In, Bool
from pyomo.network import Arc, Port
from pyomo.environ import (
    value,
    Integers,
    RangeSet,
    TransformationFactory,
    Block,
    Reference,
)

# Import IDAES cores
from idaes.models_extra.column_models import Tray, Condenser, Reboiler
from idaes.models_extra.column_models.condenser import CondenserType, TemperatureSpec
from idaes.core import declare_process_block_class, UnitModelBlockData, useDefault
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.solvers import get_solver

_log = idaeslog.getLogger(__name__)


@declare_process_block_class("TrayColumn")
class TrayColumnData(UnitModelBlockData):
    """
    Tray Column model.
    """

    CONFIG = ConfigBlock()
    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""Indicates whether this model will be dynamic or not,
**default** = False. Tray column units do not support dynamic behavior.""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be False",
            doc="""Indicates whether holdup terms should be constructed or not.
**default** - False. Tray column units do not have defined volume, thus
this must be False.""",
        ),
    )
    CONFIG.declare(
        "number_of_trays",
        ConfigValue(
            default=None,
            domain=In(Integers),
            description="Number of trays in the column",
            doc="""Indicates the number of trays to be constructed.
**default** - None.
**Valid values:** {
Must be integer number.}""",
        ),
    )
    CONFIG.declare(
        "feed_tray_location",
        ConfigValue(
            default=None,
            domain=In(Integers),
            description="Feed tray location in the column",
            doc="""Indicates the number of trays to be constructed.
**default** - None.
**Valid values:** {
Must be integer number.}""",
        ),
    )
    CONFIG.declare(
        "condenser_type",
        ConfigValue(
            default=CondenserType.totalCondenser,
            domain=In(CondenserType),
            description="Type of condenser flag",
            doc="""Indicates what type of condenser should be constructed,
**default** - CondenserType.totalCondenser.
**Valid values:** {
**CondenserType.totalCondenser** - Incoming vapor from top tray is condensed
to all liquid,
**CondenserType.partialCondenser** - Incoming vapor from top tray is
partially condensed to a vapor and liquid stream.}""",
        ),
    )
    CONFIG.declare(
        "condenser_temperature_spec",
        ConfigValue(
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
user specified temperature.}""",
        ),
    )
    CONFIG.declare(
        "has_heat_transfer",
        ConfigValue(
            default=False,
            domain=Bool,
            description="heat duty to/from tray construction flag.",
            doc="""indicates if there is heat duty to/from the tray,
**default** - False.
**Valid values:** {
**True** - include a heat duty term,
**False** - exclude a heat duty term.}""",
        ),
    )
    CONFIG.declare(
        "has_pressure_change",
        ConfigValue(
            default=False,
            domain=Bool,
            description="pressure change term construction flag",
            doc="""indicates whether terms for pressure change should be
constructed,
**default** - False.
**Valid values:** {
**True** - include pressure change terms,
**False** - exclude pressure change terms.}""",
        ),
    )
    CONFIG.declare(
        "has_liquid_side_draw",
        ConfigValue(
            default=False,
            domain=Bool,
            description="liquid side draw construction flag.",
            doc="""indicates if there is a liquid side draw from all trays,
**default** - False.
**Valid values:** {
**True** - include a liquid side draw for all trays,
**False** - exclude a liquid side draw for all trays.}""",
        ),
    )
    CONFIG.declare(
        "has_vapor_side_draw",
        ConfigValue(
            default=False,
            domain=Bool,
            description="vapor side draw construction flag.",
            doc="""indicates if there is a vapor side draw from all trays,
**default** - False.
**Valid values:** {
**True** - include a vapor side draw for all trays,
**False** - exclude a vapor side draw for all trays.}""",
        ),
    )
    CONFIG.declare(
        "property_package",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="property package to use for control volume",
            doc="""property parameter object used to define property calculations,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PropertyParameterObject** - a PropertyParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "property_package_args",
        ConfigBlock(
            implicit=True,
            description="arguments to use for constructing property packages",
            doc="""a ConfigBlock with arguments to be passed to a property block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}""",
        ),
    )

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
            self._rectification_index = RangeSet(1, self.config.feed_tray_location - 1)
            self._stripping_index = RangeSet(
                self.config.feed_tray_location + 1, self.config.number_of_trays
            )

        else:
            raise ConfigurationError(
                "The config argument number_of_trays "
                "needs to specified and cannot be None."
            )

        # Add trays
        # TODO:
        # 1. Add support for multiple feed tray locations
        # 2. Add support for specifying which trays have side draws
        self.rectification_section = Tray(
            self._rectification_index,
            has_liquid_side_draw=self.config.has_liquid_side_draw,
            has_vapor_side_draw=self.config.has_vapor_side_draw,
            has_heat_transfer=self.config.has_heat_transfer,
            has_pressure_change=self.config.has_pressure_change,
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
        )

        self.feed_tray = Tray(
            is_feed_tray=True,
            has_liquid_side_draw=self.config.has_liquid_side_draw,
            has_vapor_side_draw=self.config.has_vapor_side_draw,
            has_heat_transfer=self.config.has_heat_transfer,
            has_pressure_change=self.config.has_pressure_change,
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
        )

        self.stripping_section = Tray(
            self._stripping_index,
            has_liquid_side_draw=self.config.has_liquid_side_draw,
            has_vapor_side_draw=self.config.has_vapor_side_draw,
            has_heat_transfer=self.config.has_heat_transfer,
            has_pressure_change=self.config.has_pressure_change,
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
        )

        # Add condenser
        self.condenser = Condenser(
            condenser_type=self.config.condenser_type,
            temperature_spec=self.config.condenser_temperature_spec,
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
        )

        # Add reboiler
        self.reboiler = Reboiler(
            has_boilup_ratio=True,
            has_pressure_change=self.config.has_pressure_change,
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
        )

        # Add extension to the feed port
        self.feed = Port(extends=self.feed_tray.feed)

        # Add extensions to rectification section

        # Add extensions to stripping section

        # Construct arcs between trays, condenser, and reboiler
        self._make_rectification_arcs()
        self._make_stripping_arcs()
        self._make_feed_arcs()
        self._make_condenser_arcs()
        self._make_reboiler_arcs()

        TransformationFactory("network.expand_arcs").apply_to(self)

    def _make_rectification_arcs(self):
        self._rectification_stream_index = RangeSet(
            1, self.config.feed_tray_location - 2
        )

        def rule_liq_stream(self, i):
            return {
                "source": self.rectification_section[i].liq_out,
                "destination": self.rectification_section[i + 1].liq_in,
            }

        def rule_vap_stream(self, i):
            return {
                "source": self.rectification_section[i + 1].vap_out,
                "destination": self.rectification_section[i].vap_in,
            }

        self.rectification_liq_stream = Arc(
            self._rectification_stream_index, rule=rule_liq_stream
        )
        self.rectification_vap_stream = Arc(
            self._rectification_stream_index, rule=rule_vap_stream
        )

    def _make_stripping_arcs(self):

        self._stripping_stream_index = RangeSet(
            self.config.feed_tray_location + 1, self.config.number_of_trays - 1
        )

        def rule_liq_stream(self, i):
            return {
                "source": self.stripping_section[i].liq_out,
                "destination": self.stripping_section[i + 1].liq_in,
            }

        def rule_vap_stream(self, i):
            return {
                "source": self.stripping_section[i + 1].vap_out,
                "destination": self.stripping_section[i].vap_in,
            }

        self.stripping_liq_stream = Arc(
            self._stripping_stream_index, rule=rule_liq_stream
        )
        self.stripping_vap_stream = Arc(
            self._stripping_stream_index, rule=rule_vap_stream
        )

    def _make_feed_arcs(self):

        self.feed_liq_in = Arc(
            source=self.rectification_section[
                self.config.feed_tray_location - 1
            ].liq_out,
            destination=self.feed_tray.liq_in,
        )

        self.feed_liq_out = Arc(
            source=self.feed_tray.liq_out,
            destination=self.stripping_section[
                self.config.feed_tray_location + 1
            ].liq_in,
        )

        self.feed_vap_in = Arc(
            source=self.stripping_section[self.config.feed_tray_location + 1].vap_out,
            destination=self.feed_tray.vap_in,
        )

        self.feed_vap_out = Arc(
            source=self.feed_tray.vap_out,
            destination=self.rectification_section[
                self.config.feed_tray_location - 1
            ].vap_in,
        )

    def _make_condenser_arcs(self):

        self.condenser_vap_in = Arc(
            source=self.rectification_section[1].vap_out,
            destination=self.condenser.inlet,
        )

        self.condenser_reflux_out = Arc(
            source=self.condenser.reflux,
            destination=self.rectification_section[1].liq_in,
        )

    def _make_reboiler_arcs(self):

        self.reboiler_liq_in = Arc(
            source=self.stripping_section[self.config.number_of_trays].liq_out,
            destination=self.reboiler.inlet,
        )

        self.reboiler_vap_out = Arc(
            source=self.reboiler.vapor_reboil,
            destination=self.stripping_section[self.config.number_of_trays].vap_in,
        )

    def propagate_stream_state(self, source=None, destination=None):
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

    def initialize(
        self,
        state_args_feed=None,
        state_args_liq=None,
        state_args_vap=None,
        solver=None,
        optarg=None,
        outlvl=idaeslog.NOTSET,
    ):

        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        init_log.info("Begin initialization.")

        solverobj = get_solver(solver, optarg)

        feed_flags = self.feed_tray.initialize(
            solver=solver, optarg=optarg, outlvl=outlvl
        )

        self.propagate_stream_state(
            source=self.feed_tray.vap_out, destination=self.condenser.inlet
        )

        self.condenser.initialize(solver=solver, optarg=optarg, outlvl=outlvl)

        self.propagate_stream_state(
            source=self.feed_tray.liq_out, destination=self.reboiler.inlet
        )

        self.reboiler.initialize(solver=solver, optarg=optarg, outlvl=outlvl)

        # initialize the rectification section
        for i in self._rectification_index:
            self.propagate_stream_state(
                source=self.condenser.reflux,
                destination=self.rectification_section[i].liq_in,
            )
            self.propagate_stream_state(
                source=self.feed_tray.vap_out,
                destination=self.rectification_section[i].vap_in,
            )
            if i == 1:
                rect_liq_flags = self.rectification_section[i].initialize(
                    hold_state_liq=True,
                    hold_state_vap=False,
                    solver=solver,
                    optarg=optarg,
                    outlvl=outlvl,
                )
            elif i == len(self._rectification_index):
                rect_vap_flags = self.rectification_section[i].initialize(
                    hold_state_liq=False,
                    hold_state_vap=True,
                    solver=solver,
                    optarg=optarg,
                    outlvl=outlvl,
                )
            else:
                self.rectification_section[i].initialize(
                    solver=solver, optarg=optarg, outlvl=outlvl
                )

        # initialize the stripping section
        for i in self._stripping_index:
            self.propagate_stream_state(
                source=self.feed_tray.liq_out,
                destination=self.stripping_section[i].liq_in,
            )
            self.propagate_stream_state(
                source=self.reboiler.vapor_reboil,
                destination=self.stripping_section[i].vap_in,
            )
            if i == self.config.feed_tray_location + 1:
                strip_liq_flags = self.stripping_section[i].initialize(
                    hold_state_liq=True,
                    hold_state_vap=False,
                    solver=solver,
                    optarg=optarg,
                    outlvl=outlvl,
                )
            elif i == self.config.number_of_trays:
                strip_vap_flags = self.stripping_section[i].initialize(
                    hold_state_liq=False,
                    hold_state_vap=True,
                    solver=solver,
                    optarg=optarg,
                    outlvl=outlvl,
                )
            else:
                self.stripping_section[i].initialize(
                    solver=None, optarg=optarg, outlvl=outlvl
                )

        # For initialization purposes and to enable solving individual sections
        # creating a temp block. Note that this temp block is a reference to
        # the rectification, stripping, and feed sections. Also, expanded arcs
        # are added to the temp block as the initialization solve proceeds.
        self._temp_block = Block()

        self._temp_block.rectification = Block()

        # adding reference to the rectification section and the expanded
        # vapor and liquid arcs
        self._temp_block.rectification.trays = Reference(self.rectification_section)
        self._temp_block.rectification.expanded_liq_stream = Reference(
            self.rectification_liq_stream[:].expanded_block
        )
        self._temp_block.rectification.expanded_vap_stream = Reference(
            self.rectification_vap_stream[:].expanded_block
        )

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = solverobj.solve(self._temp_block.rectification, tee=slc.tee)
        init_log.info(
            "Rectification section initialization status {}.".format(
                idaeslog.condition(res)
            )
        )

        self._temp_block.stripping = Block()

        # adding reference to the stripping section and the expanded
        # vapor and liquid arcs
        self._temp_block.stripping.trays = Reference(self.stripping_section)
        self._temp_block.stripping.expanded_liq_stream = Reference(
            self.stripping_liq_stream[:].expanded_block
        )
        self._temp_block.stripping.expanded_vap_stream = Reference(
            self.stripping_vap_stream[:].expanded_block
        )

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = solverobj.solve(self._temp_block.stripping, tee=slc.tee)
        init_log.info(
            "Stripping section initialization status {}.".format(
                idaeslog.condition(res)
            )
        )

        # releasing the fixed inlets for the vap in to the rectification
        # to enable connection with the feed tray vap out
        self.rectification_section[
            len(self._rectification_index)
        ].properties_in_vap.release_state(flags=rect_vap_flags, outlvl=outlvl)

        # releasing the fixed inlets for the liq in to the stripping
        # to enable connection with the feed tray liq out
        self.stripping_section[
            self.config.feed_tray_location + 1
        ].properties_in_liq.release_state(flags=strip_liq_flags, outlvl=outlvl)

        # Adding the feed tray to temp block solve
        self._temp_block.feed_tray = Reference(self.feed_tray)
        self._temp_block.expanded_feed_liq_stream_in = Reference(
            self.feed_liq_in.expanded_block
        )
        self._temp_block.expanded_feed_liq_stream_out = Reference(
            self.feed_liq_out.expanded_block
        )
        self._temp_block.expanded_feed_vap_stream_in = Reference(
            self.feed_vap_in.expanded_block
        )
        self._temp_block.expanded_feed_vap_stream_out = Reference(
            self.feed_vap_out.expanded_block
        )

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = solverobj.solve(self._temp_block, tee=slc.tee)
        init_log.info(
            "Column section initialization status {}.".format(idaeslog.condition(res))
        )

        self.rectification_section[1].properties_in_liq.release_state(
            flags=rect_liq_flags, outlvl=outlvl
        )

        # Adding the condenser to the temp block solve
        self._temp_block.condenser = Reference(self.condenser)
        self._temp_block.expanded_condenser_vap_in = Reference(
            self.condenser_vap_in.expanded_block
        )
        self._temp_block.expanded_condenser_reflux_out = Reference(
            self.condenser_reflux_out.expanded_block
        )

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = solverobj.solve(self._temp_block, tee=slc.tee)
        init_log.info(
            "Column section + Condenser initialization status {}.".format(
                idaeslog.condition(res)
            )
        )

        self.stripping_section[
            self.config.number_of_trays
        ].properties_in_vap.release_state(flags=strip_vap_flags, outlvl=outlvl)

        # Delete the _temp_block as next solve is solving the entire column.
        # If we add the reboiler to the temp block, it will be similar to
        # solving the original block. This ensures that if
        # initialize is triggered twice, there is no implicit replacing
        # component error.
        self.del_component(self._temp_block)

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = solverobj.solve(self, tee=slc.tee)
        init_log.info(
            "Column section + Condenser + Reboiler initialization status {}.".format(
                idaeslog.condition(res)
            )
        )

        # release feed tray state once initialization is complete
        self.feed_tray.properties_in_feed.release_state(flags=feed_flags, outlvl=outlvl)

        # TODO : This fails in the current model
        # if not check_optimal_termination(res):
        #     raise InitializationError(
        #         f"{self.name} failed to initialize successfully. Please check "
        #         f"the output logs for more information.")
