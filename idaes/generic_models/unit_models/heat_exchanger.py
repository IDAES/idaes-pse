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
    log,
    Reference,
    PositiveReals,
    SolverFactory,
    ExternalFunction,
    Block,
    units as pyunits
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
from idaes.generic_models.unit_models.heater import (
    _make_heater_config_block,
    _make_heater_control_volume,
)

import idaes.core.util.unit_costing as costing
from idaes.core.util.misc import add_object_reference
from idaes.core.util import get_solver, scaling as iscale
from idaes.core.util.exceptions import ConfigurationError

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
            default="shell",
            domain=str,
            doc="Hot side name, sets control volume and inlet and outlet names",
        ),
    )
    config.declare(
        "cold_side_name",
        ConfigValue(
            default="tube",
            domain=str,
            doc="Cold side name, sets control volume and inlet and outlet names",
        ),
    )
    config.declare(
        "hot_side_config",
        ConfigBlock(
            implicit=True,
            description="Config block for hot side",
            doc="""A config block used to construct the hot side control volume.
This config can be given by the hot side name instead of hot_side_config.""",
        ),
    )
    config.declare(
        "cold_side_config",
        ConfigBlock(
            implicit=True,
            description="Config block for cold side",
            doc="""A config block used to construct the cold side control volume.
This config can be given by the cold side name instead of cold_side_config.""",
        ),
    )
    _make_heater_config_block(config.hot_side_config)
    _make_heater_config_block(config.cold_side_config)
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


def delta_temperature_lmtd_callback(b):
    """
    This is a callback for a temperature difference expression to calculate
    :math:`\Delta T` in the heat exchanger model using log-mean temperature
    difference (LMTD).  It can be supplied to "delta_temperature_callback"
    HeatExchanger configuration option.
    """
    dT1 = b.delta_temperature_in
    dT2 = b.delta_temperature_out

    @b.Expression(b.flowsheet().time)
    def delta_temperature(b, t):
        return (dT1[t] - dT2[t]) / log(dT1[t] / dT2[t])


def delta_temperature_amtd_callback(b):
    """
    This is a callback for a temperature difference expression to calculate
    :math:`\Delta T` in the heat exchanger model using arithmetic-mean
    temperature difference (AMTD).  It can be supplied to
    "delta_temperature_callback" HeatExchanger configuration option.
    """
    dT1 = b.delta_temperature_in
    dT2 = b.delta_temperature_out

    @b.Expression(b.flowsheet().time)
    def delta_temperature(b, t):
        return (dT1[t] + dT2[t]) * 0.5


def delta_temperature_underwood_callback(b):
    """
    This is a callback for a temperature difference expression to calculate
    :math:`\Delta T` in the heat exchanger model using log-mean temperature
    difference (LMTD) approximation given by Underwood (1970).  It can be
    supplied to "delta_temperature_callback" HeatExchanger configuration option.
    This uses a cube root function that works with negative numbers returning
    the real negative root. This should always evaluate successfully.
    """
    dT1 = b.delta_temperature_in
    dT2 = b.delta_temperature_out
    temp_units = pyunits.get_units(dT1)

    # external function that ruturns the real root, for the cuberoot of negitive
    # numbers, so it will return without error for positive and negitive dT.
    b.cbrt = ExternalFunction(
        library=functions_lib(),
        function="cbrt",
        arg_units=[temp_units])

    @b.Expression(b.flowsheet().time)
    def delta_temperature(b, t):
        return ((b.cbrt(dT1[t]) + b.cbrt(dT2[t])) / 2.0) ** 3 * temp_units


@declare_process_block_class("HeatExchanger", doc="Simple 0D heat exchanger model.")
class HeatExchangerData(UnitModelBlockData):
    """
    Simple 0D heat exchange unit.
    Unit model to transfer heat from one material to another.
    """

    CONFIG = UnitModelBlockData.CONFIG(implicit=True)
    _make_heat_exchanger_config(CONFIG)

    def _process_config(self):
        """Check for configuration errors and alternate config option names.
        """
        config = self.config

        if config.hot_side_name == config.cold_side_name:
            raise NameError(
                "Heatexchanger hot and cold side cannot have the same name '{}'."
                " Be sure to set both the hot_side_name and cold_side_name.".format(
                    config.hot_side_name
                )
            )

        for o in config:
            if not (
                o in self.CONFIG or o in [config.hot_side_name, config.cold_side_name]
            ):
                raise KeyError("Heatexchanger config option {} not defined".format(o))

        if config.hot_side_name in config:
            config.hot_side_config.set_value(config[config.hot_side_name])
            # Allow access to hot_side_config under the hot_side_name, backward
            # compatible with the tube and shell notation
            setattr(config, config.hot_side_name, config.hot_side_config)
        if config.cold_side_name in config:
            config.cold_side_config.set_value(config[config.cold_side_name])
            # Allow access to hot_side_config under the cold_side_name, backward
            # compatible with the tube and shell notation
            setattr(config, config.cold_side_name, config.cold_side_config)

        if config.cold_side_name in ["hot_side", "side_1"]:
            raise ConfigurationError("Cold side name cannot be in ['hot_side', 'side_1'].")
        if config.hot_side_name in ["cold_side", "side_2"]:
            raise ConfigurationError("Hot side name cannot be in ['cold_side', 'side_2'].")

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
        self._process_config()
        config = self.config

        ########################################################################
        # Add control volumes                                                  #
        ########################################################################
        hot_side = _make_heater_control_volume(
            self,
            config.hot_side_name,
            config.hot_side_config,
            dynamic=config.dynamic,
            has_holdup=config.has_holdup,
        )
        cold_side = _make_heater_control_volume(
            self,
            config.cold_side_name,
            config.cold_side_config,
            dynamic=config.dynamic,
            has_holdup=config.has_holdup,
        )
        # Add references to the hot side and cold side, so that we have solid
        # names to refer to internally.  side_1 and side_2 also maintain
        # compatability with older models.  Using add_object_reference keeps
        # these from showing up when you iterate through pyomo compoents in a
        # model, so only the user specified control volume names are "seen"
        if not hasattr(self, "side_1"):
            add_object_reference(self, "side_1", hot_side)
        if not hasattr(self, "side_2"):
            add_object_reference(self, "side_2", cold_side)
        if not hasattr(self, "hot_side"):
            add_object_reference(self, "hot_side", hot_side)
        if not hasattr(self, "cold_side"):
            add_object_reference(self, "cold_side", cold_side)

        ########################################################################
        # Add variables                                                        #
        ########################################################################
        # Use hot side units as basis
        s1_metadata = config.hot_side_config.property_package.get_metadata()

        q_units = s1_metadata.get_derived_units("power")
        u_units = s1_metadata.get_derived_units("heat_transfer_coefficient")
        a_units = s1_metadata.get_derived_units("area")
        temp_units = s1_metadata.get_derived_units("temperature")

        u = self.overall_heat_transfer_coefficient = Var(
            self.flowsheet().time,
            domain=PositiveReals,
            initialize=100.0,
            doc="Overall heat transfer coefficient",
            units=u_units
        )
        a = self.area = Var(
            domain=PositiveReals,
            initialize=1000.0,
            doc="Heat exchange area",
            units=a_units
        )
        self.delta_temperature_in = Var(
            self.flowsheet().time,
            initialize=10.0,
            doc="Temperature difference at the hot inlet end",
            units=temp_units
        )
        self.delta_temperature_out = Var(
            self.flowsheet().time,
            initialize=10.1,
            doc="Temperature difference at the hot outlet end",
            units=temp_units
        )
        if self.config.flow_pattern == HeatExchangerFlowPattern.crossflow:
            self.crossflow_factor = Var(
                self.flowsheet().time,
                initialize=1.0,
                doc="Factor to adjust coutercurrent flow heat "
                "transfer calculation for cross flow.",
            )
            f = self.crossflow_factor

        self.heat_duty = Reference(cold_side.heat)
        ########################################################################
        # Add ports                                                            #
        ########################################################################
        i1 = self.add_inlet_port(
            name=f"{config.hot_side_name}_inlet",
            block=hot_side,
            doc="Hot side inlet")
        i2 = self.add_inlet_port(
            name=f"{config.cold_side_name}_inlet",
            block=cold_side,
            doc="Cold side inlet")
        o1 = self.add_outlet_port(
            name=f"{config.hot_side_name}_outlet",
            block=hot_side,
            doc="Hot side outlet")
        o2 = self.add_outlet_port(
            name=f"{config.cold_side_name}_outlet",
            block=cold_side,
            doc="Cold side outlet")

        # Using Andrew's function for now.  I want these port names for backward
        # compatablity, but I don't want them to appear if you iterate throught
        # components and add_object_reference hides them from Pyomo.
        if not hasattr(self, "inlet_1"):
            add_object_reference(self, "inlet_1", i1)
        if not hasattr(self, "inlet_2"):
            add_object_reference(self, "inlet_2", i2)
        if not hasattr(self, "outlet_1"):
            add_object_reference(self, "outlet_1", o1)
        if not hasattr(self, "outlet_2"):
            add_object_reference(self, "outlet_2", o2)

        if not hasattr(self, "hot_inlet"):
            add_object_reference(self, "hot_inlet", i1)
        if not hasattr(self, "cold_inlet"):
            add_object_reference(self, "cold_inlet", i2)
        if not hasattr(self, "hot_outlet"):
            add_object_reference(self, "hot_outlet", o1)
        if not hasattr(self, "cold_outlet"):
            add_object_reference(self, "cold_outlet", o2)
        ########################################################################
        # Add end temperature differnece constraints                           #
        ########################################################################

        @self.Constraint(self.flowsheet().time)
        def delta_temperature_in_equation(b, t):
            if b.config.flow_pattern == HeatExchangerFlowPattern.cocurrent:
                return (
                    b.delta_temperature_in[t]
                    == hot_side.properties_in[t].temperature
                    - pyunits.convert(cold_side.properties_in[t].temperature,
                                      to_units=temp_units)
                )
            else:
                return (
                    b.delta_temperature_in[t]
                    == hot_side.properties_in[t].temperature
                    - pyunits.convert(cold_side.properties_out[t].temperature,
                                      to_units=temp_units)
                )

        @self.Constraint(self.flowsheet().time)
        def delta_temperature_out_equation(b, t):
            if b.config.flow_pattern == HeatExchangerFlowPattern.cocurrent:
                return (
                    b.delta_temperature_out[t]
                    == hot_side.properties_out[t].temperature
                    - pyunits.convert(cold_side.properties_out[t].temperature,
                                      to_units=temp_units)
                )
            else:
                return (
                    b.delta_temperature_out[t]
                    == hot_side.properties_out[t].temperature
                    - pyunits.convert(cold_side.properties_in[t].temperature,
                                      to_units=temp_units)
                )

        ########################################################################
        # Add a unit level energy balance                                      #
        ########################################################################
        @self.Constraint(self.flowsheet().time)
        def unit_heat_balance(b, t):
            return 0 == (hot_side.heat[t] +
                         pyunits.convert(cold_side.heat[t],
                                         to_units=q_units))

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
                    f[t] * u[t] * a * deltaT[t])
            else:
                return pyunits.convert(self.heat_duty[t], to_units=q_units) == (
                    u[t] * a * deltaT[t])

        ########################################################################
        # Add symbols for LaTeX equation rendering                             #
        ########################################################################
        self.overall_heat_transfer_coefficient.latex_symbol = "U"
        self.area.latex_symbol = "A"
        hot_side.heat.latex_symbol = "Q_1"
        cold_side.heat.latex_symbol = "Q_2"
        self.delta_temperature.latex_symbol = "\\Delta T"

    def initialize(
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

        hot_side = getattr(self, self.config.hot_side_name)
        cold_side = getattr(self, self.config.cold_side_name)

        # Create solver
        opt = get_solver(solver, optarg)

        flags1 = hot_side.initialize(
            outlvl=outlvl, optarg=optarg, solver=solver, state_args=state_args_1
        )

        init_log.info_high("Initialization Step 1a (hot side) Complete.")

        flags2 = cold_side.initialize(
            outlvl=outlvl, optarg=optarg, solver=solver, state_args=state_args_2
        )

        init_log.info_high("Initialization Step 1b (cold side) Complete.")
        # ---------------------------------------------------------------------
        # Solve unit without heat transfer equation
        # if costing block exists, deactivate
        if hasattr(self, "costing"):
            self.costing.deactivate()

        self.heat_transfer_equation.deactivate()

        # Get side 1 and side 2 heat units, and convert duty as needed
        s1_units = hot_side.heat.get_units()
        s2_units = cold_side.heat.get_units()

        if duty is None:
            # Assume 1000 J/s and check for unitless properties
            if s1_units is None and s2_units is None:
                # Backwards compatability for unitless properties
                s1_duty = - 1000
                s2_duty = 1000
            else:
                s1_duty = pyunits.convert_value(-1000,
                                                from_units=pyunits.W,
                                                to_units=s1_units)
                s2_duty = pyunits.convert_value(1000,
                                                from_units=pyunits.W,
                                                to_units=s2_units)
        else:
            # Duty provided with explicit units
            s1_duty = -pyunits.convert_value(duty[0],
                                             from_units=duty[1],
                                             to_units=s1_units)
            s2_duty = pyunits.convert_value(duty[0],
                                            from_units=duty[1],
                                            to_units=s2_units)

        cold_side.heat.fix(s2_duty)
        for i in hot_side.heat:
            hot_side.heat[i].value = s1_duty

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        init_log.info_high("Initialization Step 2 {}.".format(idaeslog.condition(res)))
        cold_side.heat.unfix()
        self.heat_transfer_equation.activate()
        # ---------------------------------------------------------------------
        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        init_log.info_high("Initialization Step 3 {}.".format(idaeslog.condition(res)))
        # ---------------------------------------------------------------------
        # Release Inlet state
        hot_side.release_state(flags1, outlvl=outlvl)
        cold_side.release_state(flags2, outlvl=outlvl)

        init_log.info("Initialization Completed, {}".format(idaeslog.condition(res)))
        # if costing block exists, activate and initialize
        if hasattr(self, "costing"):
            self.costing.activate()
            costing.initialize(self.costing)

    def _get_performance_contents(self, time_point=0):
        var_dict = {
            "HX Coefficient": self.overall_heat_transfer_coefficient[time_point]
        }
        var_dict["HX Area"] = self.area
        var_dict["Heat Duty"] = self.heat_duty[time_point]
        if self.config.flow_pattern == HeatExchangerFlowPattern.crossflow:
            var_dict = {"Crossflow Factor": self.crossflow_factor[time_point]}

        expr_dict = {}
        expr_dict["Delta T Driving"] = self.delta_temperature[time_point]
        expr_dict["Delta T In"] = self.delta_temperature_in[time_point]
        expr_dict["Delta T Out"] = self.delta_temperature_out[time_point]

        return {"vars": var_dict, "exprs": expr_dict}

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {
                "Hot Inlet": self.inlet_1,
                "Hot Outlet": self.outlet_1,
                "Cold Inlet": self.inlet_2,
                "Cold Outlet": self.outlet_2,
            },
            time_point=time_point,
        )

    def get_costing(self, module=costing, year=None, **kwargs):
        if not hasattr(self.flowsheet(), "costing"):
            self.flowsheet().get_costing(year=year)

        self.costing = Block()
        module.hx_costing(self.costing, **kwargs)

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        # We have a pretty good idea that the delta Ts will be between about
        # 1 and 100 regardless of process of temperature units, so a default
        # should be fine, so don't warn.  Guessing a typical delta t around 10
        # the default scaling factor is set to 0.1
        sf_dT1 = dict(zip(
            self.delta_temperature_in.keys(),
            [iscale.get_scaling_factor(v, default=0.1)
                for v in self.delta_temperature_in.values()]))
        sf_dT2 = dict(zip(
            self.delta_temperature_out.keys(),
            [iscale.get_scaling_factor(v, default=0.1)
                for v in self.delta_temperature_out.values()]))

        # U depends a lot on the process and units of measure so user should set
        # this one.
        sf_u = dict(zip(
            self.overall_heat_transfer_coefficient.keys(),
            [iscale.get_scaling_factor(v, default=1, warning=True)
                for v in self.overall_heat_transfer_coefficient.values()]))

        # Since this depends on the process size this is another scaling factor
        # the user should always set.
        sf_a = iscale.get_scaling_factor(self.area, default=1, warning=True)

        for t, c in self.heat_transfer_equation.items():
            iscale.constraint_scaling_transform(
                c, sf_dT1[t]*sf_u[t]*sf_a, overwrite=False)

        for t, c in self.unit_heat_balance.items():
            iscale.constraint_scaling_transform(
                c, sf_dT1[t]*sf_u[t]*sf_a, overwrite=False)

        for t, c in self.delta_temperature_in_equation.items():
            iscale.constraint_scaling_transform(c, sf_dT1[t], overwrite=False)

        for t, c in self.delta_temperature_out_equation.items():
            iscale.constraint_scaling_transform(c, sf_dT2[t], overwrite=False)

        if hasattr(self, "costing"):
            # import costing scaling factors
            costing.calculate_scaling_factors(self.costing)
