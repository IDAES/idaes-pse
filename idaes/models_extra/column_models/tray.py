#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Tray model for distillation. Can build the following:
1. Conventional tray with liq/vap inlet and liq/vap outlet
2. Feed tray with single feed inlet, liq/vap inlet, and liq/vap outlet
3. Conventional tray with side liq/vap draws (with or without feed inlet)

NOTE:
1. Does not use the IDAES control volume blocks.
2. Side vapor draw is unconventional but the model allows the user to
have a side vapor draw if required.
"""

__author__ = "Jaffer Ghouse"

# Import Pyomo libraries
from pyomo.common.config import ConfigBlock, ConfigValue, In, Bool
from pyomo.network import Port
from pyomo.environ import (
    Var,
    Set,
    value,
    check_optimal_termination,
)

# Import IDAES cores
from idaes.core import declare_process_block_class, UnitModelBlockData, useDefault
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import (
    ConfigurationError,
    InitializationError,
)
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.logger as idaeslog
from idaes.models_extra.column_models.util import make_phase_split

_log = idaeslog.getLogger(__name__)


@declare_process_block_class("Tray")
class TrayData(UnitModelBlockData):
    """
    Tray unit for distillation model.
    """

    CONFIG = ConfigBlock()
    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""Indicates whether this model will be dynamic or not,
**default** = False. Flash units do not support dynamic behavior.""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be False",
            doc="""Indicates whether holdup terms should be constructed or not.
**default** - False. Flash units do not have defined volume, thus
this must be False.""",
        ),
    )
    CONFIG.declare(
        "is_feed_tray",
        ConfigValue(
            default=False,
            domain=Bool,
            description="flag to indicate feed tray.",
            doc="""indicates if this is a feed tray and constructs
corresponding ports,
**default** - False.
**Valid values:** {
**True** - feed tray,
**False** - conventional tray with no feed inlet}""",
        ),
    )
    CONFIG.declare(
        "has_liquid_side_draw",
        ConfigValue(
            default=False,
            domain=Bool,
            description="liquid side draw construction flag.",
            doc="""indicates if there is a liquid side draw from the tray,
**default** - False.
**Valid values:** {
**True** - include a liquid side draw from the tray,
**False** - exclude a liquid side draw from the tray.}""",
        ),
    )
    CONFIG.declare(
        "has_vapor_side_draw",
        ConfigValue(
            default=False,
            domain=Bool,
            description="vapor side draw construction flag.",
            doc="""indicates if there is a vapor side draw from the tray,
**default** - False.
**Valid values:** {
**True** - include a vapor side draw from the tray,
**False** - exclude a vapor side draw from the tray.}""",
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
        # Setup model build logger
        model_log = idaeslog.getModelLogger(self.name, tag="unit")

        # Call UnitModel.build to setup dynamics
        super(TrayData, self).build()

        # Create the inlet list to build inlet state blocks
        if self.config.is_feed_tray:
            inlet_list = ["feed", "liq", "vap"]
        else:
            inlet_list = ["liq", "vap"]

        # Create a dict to set up the inlet state blocks
        state_block_args = dict(**self.config.property_package_args)
        state_block_args["has_phase_equilibrium"] = True
        state_block_args["defined_state"] = True

        # Check for default property package
        self._get_property_package()

        for i in inlet_list:
            state_obj = self.config.property_package.build_state_block(
                self.flowsheet().time,
                doc="State block for " + i + "_inlet to tray",
                **state_block_args,
            )

            setattr(self, "properties_in_" + i, state_obj)

        # Create a dict to set up the mixed outlet state blocks
        mixed_block_args = dict(**self.config.property_package_args)
        mixed_block_args["has_phase_equilibrium"] = True
        mixed_block_args["defined_state"] = False

        self.properties_out = self.config.property_package.build_state_block(
            self.flowsheet().time,
            doc="State block for mixed outlet from tray",
            **mixed_block_args,
        )

        self._add_material_balance()
        self._add_energy_balance()
        self._add_pressure_balance()

        # Get liquid and vapor phase objects from the property package
        # to be used below. Avoids repetition.
        _liquid_list = []
        _vapor_list = []
        for p in self.config.property_package.phase_list:
            pobj = self.config.property_package.get_phase(p)
            if pobj.is_vapor_phase():
                _vapor_list.append(p)
            elif pobj.is_liquid_phase():
                _liquid_list.append(p)
            else:
                _liquid_list.append(p)
                model_log.warning(
                    "A non-liquid/non-vapor phase was detected but will "
                    "be treated as a liquid."
                )

        # Create a pyomo set for indexing purposes. This set is appended to
        # model otherwise results in an abstract set.
        self._liquid_set = Set(initialize=_liquid_list)
        self._vapor_set = Set(initialize=_vapor_list)

        self._add_ports()

    def _add_material_balance(self):
        """Method to construct the mass balance equation."""

        @self.Constraint(
            self.flowsheet().time,
            self.config.property_package.component_list,
            doc="material balance",
        )
        def material_mixing_equations(b, t, j):
            if self.config.is_feed_tray:
                return 0 == sum(
                    self.properties_in_feed[t].get_material_flow_terms(p, j)
                    + self.properties_in_liq[t].get_material_flow_terms(p, j)
                    + self.properties_in_vap[t].get_material_flow_terms(p, j)
                    - self.properties_out[t].get_material_flow_terms(p, j)
                    for p in b.config.property_package.phase_list
                )
            else:
                return 0 == sum(
                    self.properties_in_liq[t].get_material_flow_terms(p, j)
                    + self.properties_in_vap[t].get_material_flow_terms(p, j)
                    - self.properties_out[t].get_material_flow_terms(p, j)
                    for p in b.config.property_package.phase_list
                )

    def _add_energy_balance(self):
        """Method to construct the energy balance equation."""

        if self.config.has_heat_transfer:
            units_meta = self.config.property_package.get_metadata()
            self.heat_duty = Var(
                self.flowsheet().time,
                initialize=0,
                doc="Heat duty for the tray",
                units=units_meta.get_derived_units("power"),
            )

        @self.Constraint(self.flowsheet().time, doc="energy balance")
        def enthalpy_mixing_equations(b, t):
            if self.config.is_feed_tray:
                if self.config.has_heat_transfer:
                    return (
                        0
                        == (
                            sum(
                                self.properties_in_feed[t].get_enthalpy_flow_terms(p)
                                for p in b.config.property_package.phase_list
                            )
                            + sum(
                                self.properties_in_liq[t].get_enthalpy_flow_terms(p)
                                for p in b.config.property_package.phase_list
                            )
                            + sum(
                                self.properties_in_vap[t].get_enthalpy_flow_terms(p)
                                for p in b.config.property_package.phase_list
                            )
                            - sum(
                                self.properties_out[t].get_enthalpy_flow_terms(p)
                                for p in b.config.property_package.phase_list
                            )
                        )
                        + self.heat_duty[t]
                    )
                else:
                    return 0 == (
                        sum(
                            self.properties_in_feed[t].get_enthalpy_flow_terms(p)
                            for p in b.config.property_package.phase_list
                        )
                        + sum(
                            self.properties_in_liq[t].get_enthalpy_flow_terms(p)
                            for p in b.config.property_package.phase_list
                        )
                        + sum(
                            self.properties_in_vap[t].get_enthalpy_flow_terms(p)
                            for p in b.config.property_package.phase_list
                        )
                        - sum(
                            self.properties_out[t].get_enthalpy_flow_terms(p)
                            for p in b.config.property_package.phase_list
                        )
                    )
            else:
                if self.config.has_heat_transfer:
                    return (
                        0
                        == (
                            sum(
                                self.properties_in_liq[t].get_enthalpy_flow_terms(p)
                                for p in b.config.property_package.phase_list
                            )
                            + sum(
                                self.properties_in_vap[t].get_enthalpy_flow_terms(p)
                                for p in b.config.property_package.phase_list
                            )
                            - sum(
                                self.properties_out[t].get_enthalpy_flow_terms(p)
                                for p in b.config.property_package.phase_list
                            )
                        )
                        + self.heat_duty[t]
                    )
                else:
                    return 0 == (
                        sum(
                            self.properties_in_liq[t].get_enthalpy_flow_terms(p)
                            for p in b.config.property_package.phase_list
                        )
                        + sum(
                            self.properties_in_vap[t].get_enthalpy_flow_terms(p)
                            for p in b.config.property_package.phase_list
                        )
                        - sum(
                            self.properties_out[t].get_enthalpy_flow_terms(p)
                            for p in b.config.property_package.phase_list
                        )
                    )

    def _add_pressure_balance(self):
        """Method to construct the pressure balance."""
        if self.config.has_pressure_change:
            units_meta = self.config.property_package.get_metadata()
            self.deltaP = Var(
                self.flowsheet().time,
                initialize=0,
                doc="Pressure drop across tray",
                units=units_meta.get_derived_units("pressure"),
            )

        @self.Constraint(self.flowsheet().time, doc="pressure balance for tray")
        def pressure_drop_equation(self, t):
            if self.config.has_pressure_change:
                return (
                    self.properties_out[t].pressure
                    == self.properties_in_liq[t].pressure - self.deltaP[t]
                )
            else:
                return (
                    self.properties_out[t].pressure
                    == self.properties_in_liq[t].pressure
                )

    def _add_ports(self):
        """Method to construct the ports for the tray."""

        # Add feed inlet port
        if self.config.is_feed_tray:
            self.add_inlet_port(name="feed", block=self.properties_in_feed)

        # Add liquid and vapor inlet ports
        self.add_inlet_port(name="liq_in", block=self.properties_in_liq)
        self.add_inlet_port(name="vap_in", block=self.properties_in_vap)

        # Add liquid outlet port
        self.liq_out = Port(noruleinit=True, doc="liquid outlet from tray")

        # Add liquid side draw port if selected
        if self.config.has_liquid_side_draw:
            self.liq_side_sf = Var(
                initialize=0.01, doc="split fraction for the liquid side draw"
            )
            self.liq_side_draw = Port(noruleinit=True, doc="liquid side draw.")
            make_phase_split(
                self,
                port=self.liq_side_draw,
                phase=self._liquid_set,
                side_sf=self.liq_side_sf,
            )

            # Populate the liquid outlet port with the remaining liquid
            # after the side draw
            make_phase_split(
                self,
                port=self.liq_out,
                phase=self._liquid_set,
                side_sf=1 - self.liq_side_sf,
            )
        else:
            # Populate the liquid outlet port when no liquid side draw
            make_phase_split(self, port=self.liq_out, phase=self._liquid_set, side_sf=1)

        # Add the vapor outlet port
        self.vap_out = Port(noruleinit=True, doc="vapor outlet from tray")

        # Add vapor side draw port if selected
        if self.config.has_vapor_side_draw:
            self.vap_side_sf = Var(
                initialize=0.01, doc="split fraction for the vapor side draw"
            )
            self.vap_side_draw = Port(noruleinit=True, doc="vapor side draw.")
            make_phase_split(
                self,
                port=self.vap_side_draw,
                phase=self._vapor_set,
                side_sf=self.vap_side_sf,
            )
            # Populate the vapor outlet port with the remaining vapor
            # after the vapor side draw
            make_phase_split(
                self,
                port=self.vap_out,
                phase=self._vapor_set,
                side_sf=1 - self.vap_side_sf,
            )
        else:
            # Populate the vapor outlet port when no vapor side draw
            make_phase_split(self, port=self.vap_out, phase=self._vapor_set, side_sf=1)

    def initialize(
        self,
        state_args_feed=None,
        state_args_liq=None,
        state_args_vap=None,
        hold_state_liq=False,
        hold_state_vap=False,
        solver=None,
        optarg=None,
        outlvl=idaeslog.NOTSET,
    ):

        # TODO:
        # 1. Initialization for dynamic mode. Currently not supported.
        # 2. Handle unfixed side split fraction vars
        # 3. Better logic to handle and fix state vars.

        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        init_log.info("Begin initialization.")

        solverobj = get_solver(solver, optarg)

        if self.config.has_liquid_side_draw:
            if not self.liq_side_sf.fixed:
                raise ConfigurationError(
                    "Liquid side draw split fraction not fixed but "
                    "has_liquid_side_draw set to True."
                )

        if self.config.has_vapor_side_draw:
            if not self.vap_side_sf.fixed:
                raise ConfigurationError(
                    "Vapor side draw split fraction not fixed but "
                    "has_vapor_side_draw set to True."
                )

        # Create initial guess if not provided by using current values
        if self.config.is_feed_tray and state_args_feed is None:
            state_args_feed = {}
            state_args_liq = {}
            state_args_vap = {}
            state_dict = self.properties_in_feed[
                self.flowsheet().time.first()
            ].define_port_members()

            for k in state_dict.keys():
                if "flow" in k:
                    if state_dict[k].is_indexed():
                        state_args_feed[k] = {}
                        state_args_liq[k] = {}
                        state_args_vap[k] = {}
                        for m in state_dict[k].keys():
                            state_args_feed[k][m] = value(state_dict[k][m])
                            state_args_liq[k][m] = value(0.1 * state_dict[k][m])
                            state_args_vap[k][m] = value(0.1 * state_dict[k][m])

                    else:
                        state_args_feed[k] = value(state_dict[k])
                        state_args_liq[k] = 0.1 * value(state_dict[k])
                        state_args_vap[k] = 0.1 * value(state_dict[k])
                else:
                    if state_dict[k].is_indexed():
                        state_args_feed[k] = {}
                        state_args_liq[k] = {}
                        state_args_vap[k] = {}
                        for m in state_dict[k].keys():
                            state_args_feed[k][m] = value(state_dict[k][m])
                            state_args_liq[k][m] = value(state_dict[k][m])
                            state_args_vap[k][m] = value(state_dict[k][m])

                    else:
                        state_args_feed[k] = value(state_dict[k])
                        state_args_liq[k] = value(state_dict[k])
                        state_args_vap[k] = value(state_dict[k])

        # Create initial guess if not provided by using current values
        if not self.config.is_feed_tray and state_args_liq is None:
            state_args_liq = {}
            state_dict = self.properties_in_liq[
                self.flowsheet().time.first()
            ].define_port_members()

            for k in state_dict.keys():
                if state_dict[k].is_indexed():
                    state_args_liq[k] = {}
                    for m in state_dict[k].keys():
                        state_args_liq[k][m] = value(state_dict[k][m])
                else:
                    state_args_liq[k] = value(state_dict[k])

        # Create initial guess if not provided by using current values
        if not self.config.is_feed_tray and state_args_vap is None:
            state_args_vap = {}
            state_dict = self.properties_in_vap[
                self.flowsheet().time.first()
            ].define_port_members()

            for k in state_dict.keys():
                if state_dict[k].is_indexed():
                    state_args_vap[k] = {}
                    for m in state_dict[k].keys():
                        state_args_vap[k][m] = value(state_dict[k][m])
                else:
                    state_args_vap[k] = value(state_dict[k])

        if self.config.is_feed_tray:
            feed_flags = self.properties_in_feed.initialize(
                outlvl=outlvl,
                solver=solver,
                optarg=optarg,
                hold_state=True,
                state_args=state_args_feed,
                state_vars_fixed=False,
            )

        liq_in_flags = self.properties_in_liq.initialize(
            outlvl=outlvl,
            solver=solver,
            optarg=optarg,
            hold_state=True,
            state_args=state_args_liq,
            state_vars_fixed=False,
        )

        vap_in_flags = self.properties_in_vap.initialize(
            outlvl=outlvl,
            solver=solver,
            optarg=optarg,
            hold_state=True,
            state_args=state_args_vap,
            state_vars_fixed=False,
        )

        # state args to initialize the mixed outlet state block
        state_args_mixed = {}

        if self.config.is_feed_tray:

            # if feed tray, initialize the mixed state block at
            # the same condition.
            state_args_mixed = state_args_feed
        else:
            # if not feed tray, initialize mixed state block at average of
            # vap/liq inlets except pressure. While this is crude, it
            # will work for most combination of state vars.
            state_dict = self.properties_in_liq[
                self.flowsheet().time.first()
            ].define_port_members()
            for k in state_dict.keys():
                if k == "pressure":
                    # Take the lowest pressure and this is the liq inlet
                    state_args_mixed[k] = value(
                        self.properties_in_liq[0].component(state_dict[k].local_name)
                    )
                elif state_dict[k].is_indexed():
                    state_args_mixed[k] = {}
                    for m in state_dict[k].keys():
                        if "flow" in k:
                            state_args_mixed[k][m] = value(
                                self.properties_in_liq[0].component(
                                    state_dict[k].local_name
                                )[m]
                            ) + value(
                                self.properties_in_vap[0].component(
                                    state_dict[k].local_name
                                )[m]
                            )

                        else:
                            state_args_mixed[k][m] = 0.5 * (
                                value(
                                    self.properties_in_liq[0].component(
                                        state_dict[k].local_name
                                    )[m]
                                )
                                + value(
                                    self.properties_in_vap[0].component(
                                        state_dict[k].local_name
                                    )[m]
                                )
                            )

                else:
                    if "flow" in k:
                        state_args_mixed[k] = value(
                            self.properties_in_liq[0].component(
                                state_dict[k].local_name
                            )
                        ) + value(
                            self.properties_in_vap[0].component(
                                state_dict[k].local_name
                            )
                        )
                    else:
                        state_args_mixed[k] = 0.5 * (
                            value(
                                self.properties_in_liq[0].component(
                                    state_dict[k].local_name
                                )
                            )
                            + value(
                                self.properties_in_vap[0].component(
                                    state_dict[k].local_name
                                )
                            )
                        )

        # Initialize the mixed outlet state block
        self.properties_out.initialize(
            outlvl=outlvl,
            solver=solver,
            optarg=optarg,
            hold_state=False,
            state_args=state_args_mixed,
            state_vars_fixed=False,
        )

        # Deactivate energy balance
        self.enthalpy_mixing_equations.deactivate()

        # Try fixing the outlet temperature if else pass
        # NOTE: if passed then there would probably be a degree of freedom
        try:
            self.properties_out[:].temperature.fix(state_args_mixed["temperature"])
        except AttributeError:
            init_log.warning(
                "Trying to fix outlet temperature "
                "during initialization but temperature attribute "
                "unavailable in the state block. Initialization "
                "proceeding with a potential degree of freedom."
            )

        # Deactivate pressure balance
        self.pressure_drop_equation.deactivate()

        # Try fixing the outlet temperature if else pass
        # NOTE: if passed then there would probably be a degree of freedom
        try:
            self.properties_out[:].pressure.fix(state_args_mixed["pressure"])
        except AttributeError:
            init_log.warning(
                "Trying to fix outlet pressure "
                "during initialization but pressure attribute "
                "unavailable in the state block. Initialization "
                "proceeding with a potential degree of freedom."
            )

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = solverobj.solve(self, tee=slc.tee)
        init_log.info("Mass balance solve {}.".format(idaeslog.condition(res)))

        # Activate energy balance
        self.enthalpy_mixing_equations.activate()
        try:
            self.properties_out[:].temperature.unfix()
        except AttributeError:
            pass

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = solverobj.solve(self, tee=slc.tee)
        init_log.info(
            "Mass and energy balance solve {}.".format(idaeslog.condition(res))
        )

        # Activate pressure balance
        self.pressure_drop_equation.activate()
        try:
            self.properties_out[:].pressure.unfix()
        except AttributeError:
            pass

        if degrees_of_freedom(self) == 0:
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                res = solverobj.solve(self, tee=slc.tee)
            init_log.info(
                "Mass, energy and pressure balance solve {}.".format(
                    idaeslog.condition(res)
                )
            )
        else:
            raise Exception(
                "State vars fixed but degrees of freedom "
                "for tray block is not zero during "
                "initialization."
            )

        if not check_optimal_termination(res):
            raise InitializationError(
                f"{self.name} failed to initialize successfully. Please check "
                f"the output logs for more information."
            )

        init_log.info(
            "Initialization complete, status {}.".format(idaeslog.condition(res))
        )

        if not self.config.is_feed_tray:
            if not hold_state_vap:
                self.properties_in_vap.release_state(flags=vap_in_flags, outlvl=outlvl)
            if not hold_state_liq:
                self.properties_in_liq.release_state(flags=liq_in_flags, outlvl=outlvl)
            if hold_state_liq and hold_state_vap:
                return liq_in_flags, vap_in_flags
            elif hold_state_vap:
                return vap_in_flags
            elif hold_state_liq:
                return liq_in_flags
        else:
            self.properties_in_liq.release_state(flags=liq_in_flags, outlvl=outlvl)
            self.properties_in_vap.release_state(flags=vap_in_flags, outlvl=outlvl)
            return feed_flags
