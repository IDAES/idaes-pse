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
Heat Exchanger Models.
"""

__author__ = "John Eslick"

from enum import Enum

# Import Pyomo libraries
from pyomo.environ import (
    Var,
    Param,
    log,
    Reference,
    PositiveReals,
    ExternalFunction,
    units as pyunits,
    check_optimal_termination,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    UnitModelBlockData,
)

import idaes.logger as idaeslog
from idaes.core.util.functions import functions_lib
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.models.unit_models.heater import (
    _make_heater_config_block,
    _make_heater_control_volume,
)

from idaes.core.util.misc import add_object_reference
from idaes.core.util import scaling as iscale
from idaes.core.solvers import get_solver
from idaes.core.util.exceptions import ConfigurationError, InitializationError


_log = idaeslog.getLogger(__name__)


class HeatExchangerFlowPattern(Enum):
    countercurrent = 1
    cocurrent = 2
    crossflow = 3


def _make_heat_exchanger_config(config):
    """
    Declare configuration options for HeatExchangerData block.
    """
    config.declare(
        "hot_side_name",
        ConfigValue(
            default=None,
            domain=str,
            doc="Hot side name, sets control volume and inlet and outlet names",
        ),
    )
    config.declare(
        "cold_side_name",
        ConfigValue(
            default=None,
            domain=str,
            doc="Cold side name, sets control volume and inlet and outlet names",
        ),
    )
    config.declare(
        "hot_side",
        ConfigBlock(
            description="Config block for hot side",
            doc="""A config block used to construct the hot side control volume.
This config can be given by the hot side name instead of hot_side.""",
        ),
    )
    config.declare(
        "cold_side",
        ConfigBlock(
            description="Config block for cold side",
            doc="""A config block used to construct the cold side control volume.
This config can be given by the cold side name instead of cold_side.""",
        ),
    )
    _make_heater_config_block(config.hot_side)
    _make_heater_config_block(config.cold_side)
    config.declare(
        "delta_temperature_callback",
        ConfigValue(
            default=delta_temperature_lmtd_callback,
            description="Callback for for temperature difference calculations",
        ),
    )
    config.declare(
        "flow_pattern",
        ConfigValue(
            default=HeatExchangerFlowPattern.countercurrent,
            domain=In(HeatExchangerFlowPattern),
            description="Heat exchanger flow pattern",
            doc="""Heat exchanger flow pattern,
**default** - HeatExchangerFlowPattern.countercurrent.
**Valid values:** {
**HeatExchangerFlowPattern.countercurrent** - countercurrent flow,
**HeatExchangerFlowPattern.cocurrent** - cocurrent flow,
**HeatExchangerFlowPattern.crossflow** - cross flow, factor times
countercurrent temperature difference.}""",
        ),
    )


def delta_temperature_lmtd_smooth_callback(b):
    r"""
    This is a callback for a temperature difference expression to calculate
    :math:`\Delta T` in the heat exchanger model using log-mean temperature
    difference (LMTD) approximation from Kazi, S., M. Short, A.J. Isafiade,
    L.T. Biegler. "Heat exchanger network synthesus with detailed exchanger
    designs - 2. Hybrid optimization strategy for synthesis of heat exchanger
    networks". AIChE Journal, 2020.

    .. math::

        \alpha = \frac{\Delta T_2}{\Delta T_1}

    .. math::

        \Delta T = \frac{\Delta T_1 \sqrt{(\alpha - 1)^2 + \varepsilon}}{\sqrt{(\log_e{\alpha})^2 + \varepsilon}}

    where :math:`\Delta T_1` is the temperature difference at the hot inlet end,
    :math:`\Delta T_2` is the temperature difference at the hot outlet end, and
    :math:`\varepsilon` is the smoothing parameter.

    The smoothing parameter (``eps_lmtd_smoothing``) is mutable and the
    default is 1e-10.
    """
    b.eps_lmtd_smoothing = Param(mutable=True, initialize=1e-10)
    eps = b.eps_lmtd_smoothing
    dT1 = b.delta_temperature_in
    dT2 = b.delta_temperature_out

    @b.Expression(b.flowsheet().time)
    def delta_temperature(b, t):
        alpha = dT2[t] / dT1[t]
        return (
            dT1[t]
            * ((alpha - 1) ** 2 + eps) ** (1 / 2)
            / (log(alpha) ** 2 + eps) ** (1 / 2)
        )


def delta_temperature_lmtd_callback(b):
    r"""
    This is a callback for a temperature difference expression to calculate
    :math:`\Delta T` in the heat exchanger model using log-mean temperature
    difference (LMTD).  It can be supplied to "delta_temperature_callback"
    HeatExchanger configuration option.  This form is

    .. math::

        \Delta T = \frac{\Delta T_1 - \Delta T_2}{
            \log_e\left(\frac{\Delta T_1}{\Delta T_2}\right)}

    where :math:`\Delta T_1` is the temperature difference at the hot inlet end
    and :math:`\Delta T_2` is the temperature difference at the hot outlet end.
    """
    dT1 = b.delta_temperature_in
    dT2 = b.delta_temperature_out

    @b.Expression(b.flowsheet().time)
    def delta_temperature(b, t):
        return (dT1[t] - dT2[t]) / log(dT1[t] / dT2[t])


def delta_temperature_lmtd2_callback(b):
    r"""
    This is a callback for a temperature difference expression to calculate
    :math:`\Delta T` in the heat exchanger model using log-mean temperature
    difference (LMTD).  It can be supplied to "delta_temperature_callback"
    HeatExchanger configuration option. This form is

    .. math::

        \Delta T = \frac{\Delta T_2 - \Delta T_1}{
            \log_e\left(\frac{\Delta T_2}{\Delta T_1}\right)}

    where :math:`\Delta T_1` is the temperature difference at the hot inlet end
    and :math:`\Delta T_2` is the temperature difference at the hot outlet end.
    """
    dT1 = b.delta_temperature_in
    dT2 = b.delta_temperature_out

    @b.Expression(b.flowsheet().time)
    def delta_temperature(b, t):
        return (dT2[t] - dT1[t]) / log(dT2[t] / dT1[t])


def delta_temperature_lmtd3_callback(b):
    r"""
    This is a callback for a temperature difference expression to calculate
    :math:`\Delta T` in the heat exchanger model using log-mean temperature
    difference (LMTD).  It can be supplied to "delta_temperature_callback"
    HeatExchanger configuration option. This form is

    .. math::

        \Delta T = \frac{\Delta T_1 - \Delta T_2}{
            \log_e\left(\Delta T_2\right) - \log_e\left(\Delta T_1\right)}

    where :math:`\Delta T_1` is the temperature difference at the hot inlet end
    and :math:`\Delta T_2` is the temperature difference at the hot outlet end.
    """
    dT1 = b.delta_temperature_in
    dT2 = b.delta_temperature_out

    @b.Expression(b.flowsheet().time)
    def delta_temperature(b, t):
        return (dT2[t] - dT1[t]) / (log(dT2[t]) - log(dT1[t]))


def delta_temperature_amtd_callback(b):
    r"""
    This is a callback for a temperature difference expression to calculate
    :math:`\Delta T` in the heat exchanger model using arithmetic-mean
    temperature difference (AMTD).  It can be supplied to
    "delta_temperature_callback" HeatExchanger configuration option. This form is

    .. math::

        \Delta T = \frac{\Delta T_1 + \Delta T_2}{2}

    where :math:`\Delta T_1` is the temperature difference at the hot inlet end
    and :math:`\Delta T_2` is the temperature difference at the hot outlet end.
    """
    dT1 = b.delta_temperature_in
    dT2 = b.delta_temperature_out

    @b.Expression(b.flowsheet().time)
    def delta_temperature(b, t):
        return (dT1[t] + dT2[t]) * 0.5


def delta_temperature_underwood_callback(b):
    r"""
    This is a callback for a temperature difference expression to calculate
    :math:`\Delta T` in the heat exchanger model using log-mean temperature
    difference (LMTD) approximation given by Underwood (1970).  It can be
    supplied to "delta_temperature_callback" HeatExchanger configuration option.
    This uses a cube root function that works with negative numbers returning
    the real negative root. This should always evaluate successfully. This form
    is

    .. math::

        \Delta T = \left(\frac{
            \Delta T_1^\frac{1}{3} + \Delta T_2^\frac{1}{3}}{2}\right)^3

    where :math:`\Delta T_1` is the temperature difference at the hot inlet end
    and :math:`\Delta T_2` is the temperature difference at the hot outlet end.
    """
    dT1 = b.delta_temperature_in
    dT2 = b.delta_temperature_out
    temp_units = pyunits.get_units(dT1[dT1.index_set().first()])

    # external function that ruturns the real root, for the cuberoot of negitive
    # numbers, so it will return without error for positive and negitive dT.
    b.cbrt = ExternalFunction(
        library=functions_lib(), function="cbrt", arg_units=[temp_units]
    )

    @b.Expression(b.flowsheet().time)
    def delta_temperature(b, t):
        return ((b.cbrt(dT1[t]) + b.cbrt(dT2[t])) / 2.0) ** 3 * temp_units


def hx_process_config(self):
    """Check for configuration errors and alternate config option names."""
    config = self.config

    if config.cold_side_name in ["hot_side", "cold_side"]:
        raise ConfigurationError(f"cold_side_name cannot be '{config.cold_side_name}'.")
    if config.hot_side_name in ["hot_side", "cold_side"]:
        raise ConfigurationError(f"hot_side_name cannot be '{config.hot_side_name}'.")

    if (
        config.hot_side_name is not None
        and config.cold_side_name is not None
        and config.hot_side_name == config.cold_side_name
    ):
        raise NameError(
            f"HeatExchanger hot and cold side cannot have the same name "
            f"'{config.hot_side_name}'."
        )

    for o in config:
        if not (o in self.CONFIG or o in [config.hot_side_name, config.cold_side_name]):
            raise KeyError("HeatExchanger config option {} not defined".format(o))

    if config.hot_side_name is not None and config.hot_side_name in config:
        config.hot_side.set_value(config[config.hot_side_name])
        # Allow access to hot_side under the hot_side_name, backward
        # compatible with the tube and shell notation
        setattr(config, config.hot_side_name, config.hot_side)
    if config.cold_side_name is not None and config.cold_side_name in config:
        config.cold_side.set_value(config[config.cold_side_name])
        # Allow access to hot_side under the cold_side_name, backward
        # compatible with the tube and shell notation
        setattr(config, config.cold_side_name, config.cold_side)


def add_hx_references(self):
    """
    Method to add common references for hot and cold sides in heat exchangers.

    Args:
        None

    Returns:
        None
    """
    # Add references to the user provided aliases (if applicable).
    # Using add_object_reference keeps these from showing up when you
    # iterate through pyomo components in a model
    if self.config.hot_side_name is not None:
        if not hasattr(self, self.config.hot_side_name):
            add_object_reference(self, self.config.hot_side_name, self.hot_side)
        else:
            raise ValueError(
                f"{self.name} could not assign hot side alias {self.config.hot_side_name} "
                f"as an attribute of that name already exists."
            )
        if not hasattr(self, self.config.hot_side_name + "_inlet"):
            add_object_reference(
                self, self.config.hot_side_name + "_inlet", self.hot_side_inlet
            )
        else:
            raise ValueError(
                f"{self.name} could not assign hot side inlet alias {self.config.hot_side_name}_inlet "
                f"as an attribute of that name already exists."
            )
        if not hasattr(self, self.config.hot_side_name + "_outlet"):
            add_object_reference(
                self, self.config.hot_side_name + "_outlet", self.hot_side_outlet
            )
        else:
            raise ValueError(
                f"{self.name} could not assign hot side outlet alias {self.config.hot_side_name}_outlet "
                f"as an attribute of that name already exists."
            )
    if self.config.cold_side_name is not None:
        if not hasattr(self, self.config.cold_side_name):
            add_object_reference(self, self.config.cold_side_name, self.cold_side)
        else:
            raise ValueError(
                f"{self.name} could not assign cold side alias {self.config.cold_side_name} "
                f"as an attribute of that name already exists."
            )
        if self.config.cold_side_name is not None:
            if not hasattr(self, self.config.cold_side_name + "_inlet"):
                add_object_reference(
                    self, self.config.cold_side_name + "_inlet", self.cold_side_inlet
                )
            else:
                raise ValueError(
                    f"{self.name} could not assign cold side inlet alias {self.config.cold_side_name}_inlet "
                    f"as an attribute of that name already exists."
                )
            if not hasattr(self, self.config.cold_side_name + "_outlet"):
                add_object_reference(
                    self, self.config.cold_side_name + "_outlet", self.cold_side_outlet
                )
            else:
                raise ValueError(
                    f"{self.name} could not assign cold side outlet alias {self.config.cold_side_name}_outlet "
                    f"as an attribute of that name already exists."
                )


@declare_process_block_class("HeatExchanger", doc="Simple 0D heat exchanger model.")
class HeatExchangerData(UnitModelBlockData):
    """
    Simple 0D heat exchange unit.
    Unit model to transfer heat from one material to another.
    """

    CONFIG = UnitModelBlockData.CONFIG(implicit=True)
    _make_heat_exchanger_config(CONFIG)

    def build(self):
        """
        Building model

        Args:
            None
        Returns:
            None
        """
        ########################################################################
        #  Call UnitModel.build to setup dynamics and configure                #
        ########################################################################
        super().build()
        hx_process_config(self)
        config = self.config

        ########################################################################
        # Add control volumes                                                  #
        ########################################################################
        hot_side = _make_heater_control_volume(
            self,
            "hot_side",
            config.hot_side,
            dynamic=config.dynamic,
            has_holdup=config.has_holdup,
        )
        cold_side = _make_heater_control_volume(
            self,
            "cold_side",
            config.cold_side,
            dynamic=config.dynamic,
            has_holdup=config.has_holdup,
        )

        ########################################################################
        # Add variables                                                        #
        ########################################################################
        # Use hot side units as basis
        s1_metadata = config.hot_side.property_package.get_metadata()

        q_units = s1_metadata.get_derived_units("power")
        u_units = s1_metadata.get_derived_units("heat_transfer_coefficient")
        a_units = s1_metadata.get_derived_units("area")
        temp_units = s1_metadata.get_derived_units("temperature")

        u = self.overall_heat_transfer_coefficient = Var(
            self.flowsheet().time,
            domain=PositiveReals,
            initialize=100.0,
            doc="Overall heat transfer coefficient",
            units=u_units,
        )
        a = self.area = Var(
            domain=PositiveReals,
            initialize=1000.0,
            doc="Heat exchange area",
            units=a_units,
        )
        self.delta_temperature_in = Var(
            self.flowsheet().time,
            initialize=10.0,
            doc="Temperature difference at the hot inlet end",
            units=temp_units,
        )
        self.delta_temperature_out = Var(
            self.flowsheet().time,
            initialize=10.1,
            doc="Temperature difference at the hot outlet end",
            units=temp_units,
        )
        if self.config.flow_pattern == HeatExchangerFlowPattern.crossflow:
            self.crossflow_factor = Var(
                self.flowsheet().time,
                initialize=1.0,
                doc="Factor to adjust countercurrent flow heat "
                "transfer calculation for cross flow.",
            )
            f = self.crossflow_factor

        self.heat_duty = Reference(cold_side.heat)
        ########################################################################
        # Add ports                                                            #
        ########################################################################
        self.add_inlet_port(name="hot_side_inlet", block=hot_side, doc="Hot side inlet")
        self.add_inlet_port(
            name="cold_side_inlet",
            block=cold_side,
            doc="Cold side inlet",
        )
        self.add_outlet_port(
            name="hot_side_outlet", block=hot_side, doc="Hot side outlet"
        )
        self.add_outlet_port(
            name="cold_side_outlet",
            block=cold_side,
            doc="Cold side outlet",
        )
        ########################################################################
        # Add aliases                                                          #
        ########################################################################
        add_hx_references(self)

        ########################################################################
        # Add end temperature difference constraints                           #
        ########################################################################

        @self.Constraint(self.flowsheet().time)
        def delta_temperature_in_equation(b, t):
            if b.config.flow_pattern == HeatExchangerFlowPattern.cocurrent:
                return b.delta_temperature_in[t] == hot_side.properties_in[
                    t
                ].temperature - pyunits.convert(
                    cold_side.properties_in[t].temperature, to_units=temp_units
                )
            else:
                return b.delta_temperature_in[t] == hot_side.properties_in[
                    t
                ].temperature - pyunits.convert(
                    cold_side.properties_out[t].temperature, to_units=temp_units
                )

        @self.Constraint(self.flowsheet().time)
        def delta_temperature_out_equation(b, t):
            if b.config.flow_pattern == HeatExchangerFlowPattern.cocurrent:
                return b.delta_temperature_out[t] == hot_side.properties_out[
                    t
                ].temperature - pyunits.convert(
                    cold_side.properties_out[t].temperature, to_units=temp_units
                )
            else:
                return b.delta_temperature_out[t] == hot_side.properties_out[
                    t
                ].temperature - pyunits.convert(
                    cold_side.properties_in[t].temperature, to_units=temp_units
                )

        ########################################################################
        # Add a unit level energy balance                                      #
        ########################################################################
        @self.Constraint(self.flowsheet().time)
        def unit_heat_balance(b, t):
            return 0 == (
                hot_side.heat[t] + pyunits.convert(cold_side.heat[t], to_units=q_units)
            )

        ########################################################################
        # Add delta T calculations using callack function, lots of options,    #
        #   and users can provide their own if needed                          #
        ########################################################################
        config.delta_temperature_callback(self)
        ########################################################################
        # Add Heat transfer equation                                           #
        ########################################################################
        deltaT = self.delta_temperature

        @self.Constraint(self.flowsheet().time)
        def heat_transfer_equation(b, t):
            if self.config.flow_pattern == HeatExchangerFlowPattern.crossflow:
                return pyunits.convert(self.heat_duty[t], to_units=q_units) == (
                    f[t] * u[t] * a * deltaT[t]
                )
            else:
                return pyunits.convert(self.heat_duty[t], to_units=q_units) == (
                    u[t] * a * deltaT[t]
                )

        ########################################################################
        # Add symbols for LaTeX equation rendering                             #
        ########################################################################
        self.overall_heat_transfer_coefficient.latex_symbol = "U"
        self.area.latex_symbol = "A"
        hot_side.heat.latex_symbol = "Q_1"
        cold_side.heat.latex_symbol = "Q_2"
        self.delta_temperature.latex_symbol = "\\Delta T"

    def initialize_build(
        self,
        state_args_1=None,
        state_args_2=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
        duty=None,
    ):
        """
        Heat exchanger initialization method.

        Args:
            state_args_1 : a dict of arguments to be passed to the property
                initialization for the hot side (see documentation of the specific
                property package) (default = {}).
            state_args_2 : a dict of arguments to be passed to the property
                initialization for the cold side (see documentation of the specific
                property package) (default = {}).
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None, use
                     default solver options)
            solver : str indicating which solver to use during
                     initialization (default = None, use default solver)
            duty : an initial guess for the amount of heat transfered. This
                should be a tuple in the form (value, units), (default
                = (1000 J/s))

        Returns:
            None

        """
        # Set solver options
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        # Create solver
        opt = get_solver(solver, optarg)

        flags1 = self.hot_side.initialize(
            outlvl=outlvl, optarg=optarg, solver=solver, state_args=state_args_1
        )

        init_log.info_high("Initialization Step 1a (hot side) Complete.")

        flags2 = self.cold_side.initialize(
            outlvl=outlvl, optarg=optarg, solver=solver, state_args=state_args_2
        )

        init_log.info_high("Initialization Step 1b (cold side) Complete.")
        # ---------------------------------------------------------------------
        # Solve unit without heat transfer equation
        self.heat_transfer_equation.deactivate()

        # Get side 1 and side 2 heat units, and convert duty as needed
        s1_units = self.hot_side.heat.get_units()
        s2_units = self.cold_side.heat.get_units()

        # Check to see if heat duty is fixed
        # WE will assume that if the first point is fixed, it is fixed at all points
        if not self.cold_side.heat[self.flowsheet().time.first()].fixed:
            cs_fixed = False
            if duty is None:
                # Assume 1000 J/s and check for unitless properties
                if s1_units is None and s2_units is None:
                    # Backwards compatibility for unitless properties
                    s1_duty = -1000
                    s2_duty = 1000
                else:
                    s1_duty = pyunits.convert_value(
                        -1000, from_units=pyunits.W, to_units=s1_units
                    )
                    s2_duty = pyunits.convert_value(
                        1000, from_units=pyunits.W, to_units=s2_units
                    )
            else:
                # Duty provided with explicit units
                s1_duty = -pyunits.convert_value(
                    duty[0], from_units=duty[1], to_units=s1_units
                )
                s2_duty = pyunits.convert_value(
                    duty[0], from_units=duty[1], to_units=s2_units
                )

            self.cold_side.heat.fix(s2_duty)
            for i in self.hot_side.heat:
                self.hot_side.heat[i].value = s1_duty
        else:
            cd_fixed = True
            for i in self.hot_side.heat:
                self.hot_side.heat[i].set_value(self.cold_side.heat[i])

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        init_log.info_high("Initialization Step 2 {}.".format(idaeslog.condition(res)))
        if not cs_fixed:
            self.cold_side.heat.unfix()
        self.heat_transfer_equation.activate()
        # ---------------------------------------------------------------------
        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        init_log.info_high("Initialization Step 3 {}.".format(idaeslog.condition(res)))
        # ---------------------------------------------------------------------
        # Release Inlet state
        self.hot_side.release_state(flags1, outlvl=outlvl)
        self.cold_side.release_state(flags2, outlvl=outlvl)

        init_log.info("Initialization Completed, {}".format(idaeslog.condition(res)))

        if not check_optimal_termination(res):
            raise InitializationError(
                f"{self.name} failed to initialize successfully. Please check "
                f"the output logs for more information."
            )

    def _get_performance_contents(self, time_point=0):
        var_dict = {
            "HX Coefficient": self.overall_heat_transfer_coefficient[time_point]
        }
        var_dict["HX Area"] = self.area
        var_dict["Heat Duty"] = self.heat_duty[time_point]
        if self.config.flow_pattern == HeatExchangerFlowPattern.crossflow:
            var_dict["Crossflow Factor"] = self.crossflow_factor[time_point]

        expr_dict = {}
        expr_dict["Delta T Driving"] = self.delta_temperature[time_point]
        expr_dict["Delta T In"] = self.delta_temperature_in[time_point]
        expr_dict["Delta T Out"] = self.delta_temperature_out[time_point]

        return {"vars": var_dict, "exprs": expr_dict}

    def _get_stream_table_contents(self, time_point=0):
        # Get names for hot and cold sides
        hot_name = self.config.hot_side_name
        if hot_name is None:
            hot_name = "Hot Side"
        cold_name = self.config.cold_side_name
        if cold_name is None:
            cold_name = "Cold Side"
        return create_stream_table_dataframe(
            {
                f"{hot_name} Inlet": self.hot_side_inlet,
                f"{hot_name} Outlet": self.hot_side_outlet,
                f"{cold_name} Inlet": self.cold_side_inlet,
                f"{cold_name} Outlet": self.cold_side_outlet,
            },
            time_point=time_point,
        )

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        # We have a pretty good idea that the delta Ts will be between about
        # 1 and 100 regardless of process of temperature units, so a default
        # should be fine, so don't warn.  Guessing a typical delta t around 10
        # the default scaling factor is set to 0.1
        sf_dT1 = dict(
            zip(
                self.delta_temperature_in.keys(),
                [
                    iscale.get_scaling_factor(v, default=0.1)
                    for v in self.delta_temperature_in.values()
                ],
            )
        )
        sf_dT2 = dict(
            zip(
                self.delta_temperature_out.keys(),
                [
                    iscale.get_scaling_factor(v, default=0.1)
                    for v in self.delta_temperature_out.values()
                ],
            )
        )

        # U depends a lot on the process and units of measure so user should set
        # this one.
        sf_u = dict(
            zip(
                self.overall_heat_transfer_coefficient.keys(),
                [
                    iscale.get_scaling_factor(v, default=1, warning=True)
                    for v in self.overall_heat_transfer_coefficient.values()
                ],
            )
        )

        # Since this depends on the process size this is another scaling factor
        # the user should always set.
        sf_a = iscale.get_scaling_factor(self.area, default=1, warning=True)

        for t, c in self.heat_transfer_equation.items():
            iscale.constraint_scaling_transform(
                c, sf_dT1[t] * sf_u[t] * sf_a, overwrite=False
            )

        for t, c in self.unit_heat_balance.items():
            iscale.constraint_scaling_transform(
                c, sf_dT1[t] * sf_u[t] * sf_a, overwrite=False
            )

        for t, c in self.delta_temperature_in_equation.items():
            iscale.constraint_scaling_transform(c, sf_dT1[t], overwrite=False)

        for t, c in self.delta_temperature_out_equation.items():
            iscale.constraint_scaling_transform(c, sf_dT2[t], overwrite=False)
