##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
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
from idaes.core.util import scaling as iscale

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

    @b.Expression(b.flowsheet().config.time)
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

    @b.Expression(b.flowsheet().config.time)
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
    # external function that ruturns the real root, for the cuberoot of negitive
    # numbers, so it will return without error for positive and negitive dT.
    b.cbrt = ExternalFunction(library=functions_lib(), function="cbrt")
    dT1 = b.delta_temperature_in
    dT2 = b.delta_temperature_out

    @b.Expression(b.flowsheet().config.time)
    def delta_temperature(b, t):
        return ((b.cbrt(dT1[t]) + b.cbrt(dT2[t])) / 2.0) ** 3


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
        # Add variables                                                        #
        ########################################################################
        # Use side 1 units as basis
        s1_metadata = config.hot_side_config.property_package.get_metadata()

        q_units = s1_metadata.get_derived_units("power")
        u_units = s1_metadata.get_derived_units("heat_transfer_coefficient")
        a_units = s1_metadata.get_derived_units("area")
        t_units = s1_metadata.get_derived_units("temperature")

        u = self.overall_heat_transfer_coefficient = Var(
            self.flowsheet().config.time,
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
            self.flowsheet().config.time,
            initialize=10.0,
            doc="Temperature difference at the hot inlet end",
            units=t_units
        )
        self.delta_temperature_out = Var(
            self.flowsheet().config.time,
            initialize=10.1,
            doc="Temperature difference at the hot outlet end",
            units=t_units
        )
        if self.config.flow_pattern == HeatExchangerFlowPattern.crossflow:
            self.crossflow_factor = Var(
                self.flowsheet().config.time,
                initialize=1.0,
                doc="Factor to adjust coutercurrent flow heat "
                "transfer calculation for cross flow.",
            )
            f = self.crossflow_factor
        ########################################################################
        # Add control volumes                                                  #
        ########################################################################
        _make_heater_control_volume(
            self,
            "side_1",
            config.hot_side_config,
            dynamic=config.dynamic,
            has_holdup=config.has_holdup,
        )
        _make_heater_control_volume(
            self,
            "side_2",
            config.cold_side_config,
            dynamic=config.dynamic,
            has_holdup=config.has_holdup,
        )
        # Add named references to side_1 and side_2, side 1 and 2 maintain
        # backward compatability and are names the user doesn't need to worry
        # about. The sign convention for duty is heat from side 1 to side 2 is
        # positive
        add_object_reference(self, config.hot_side_name, self.side_1)
        add_object_reference(self, config.cold_side_name, self.side_2)

        # Add convienient references to heat duty.
        self.heat_duty = Reference(self.side_2.heat)
        # Need to do this as Reference does not have units (Pyomo fixing this)
        q = self.side_2.heat
        ########################################################################
        # Add ports                                                            #
        ########################################################################
        # Keep old port names, just for backward compatability
        self.add_inlet_port(name="inlet_1", block=self.side_1, doc="Hot side inlet")
        self.add_inlet_port(name="inlet_2", block=self.side_2, doc="Cold side inlet")
        self.add_outlet_port(name="outlet_1", block=self.side_1, doc="Hot side outlet")
        self.add_outlet_port(name="outlet_2", block=self.side_2, doc="Cold side outlet")

        # Using Andrew's function for now, I think Pyomo's refrence has trouble
        # with scalar (pyomo) components.
        add_object_reference(self, config.hot_side_name + "_inlet", self.inlet_1)
        add_object_reference(self, config.cold_side_name + "_inlet", self.inlet_2)
        add_object_reference(self, config.hot_side_name + "_outlet", self.outlet_1)
        add_object_reference(self, config.cold_side_name + "_outlet", self.outlet_2)
        ########################################################################
        # Add end temperature differnece constraints                            #
        ########################################################################
        @self.Constraint(self.flowsheet().config.time)
        def delta_temperature_in_equation(b, t):
            if b.config.flow_pattern == HeatExchangerFlowPattern.cocurrent:
                return (
                    b.delta_temperature_in[t]
                    == b.side_1.properties_in[t].temperature
                    - pyunits.convert(b.side_2.properties_in[t].temperature,
                                      to_units=t_units)
                )
            else:
                return (
                    b.delta_temperature_in[t]
                    == b.side_1.properties_in[t].temperature
                    - pyunits.convert(b.side_2.properties_out[t].temperature,
                                      to_units=t_units)
                )

        @self.Constraint(self.flowsheet().config.time)
        def delta_temperature_out_equation(b, t):
            if b.config.flow_pattern == HeatExchangerFlowPattern.cocurrent:
                return (
                    b.delta_temperature_out[t]
                    == b.side_1.properties_out[t].temperature
                    - pyunits.convert(b.side_2.properties_out[t].temperature,
                                      to_units=t_units)
                )
            else:
                return (
                    b.delta_temperature_out[t]
                    == b.side_1.properties_out[t].temperature
                    - pyunits.convert(b.side_2.properties_in[t].temperature,
                                      to_units=t_units)
                )

        ########################################################################
        # Add a unit level energy balance                                      #
        ########################################################################
        @self.Constraint(self.flowsheet().config.time)
        def unit_heat_balance(b, t):
            return 0 == (self.side_1.heat[t] +
                         pyunits.convert(self.side_2.heat[t],
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

        @self.Constraint(self.flowsheet().config.time)
        def heat_transfer_equation(b, t):
            if self.config.flow_pattern == HeatExchangerFlowPattern.crossflow:
                return pyunits.convert(q[t], to_units=q_units) == (
                    f[t] * u[t] * a * deltaT[t])
            else:
                return pyunits.convert(q[t], to_units=q_units) == (
                    u[t] * a * deltaT[t])

        ########################################################################
        # Add symbols for LaTeX equation rendering                             #
        ########################################################################
        self.overall_heat_transfer_coefficient.latex_symbol = "U"
        self.area.latex_symbol = "A"
        self.side_1.heat.latex_symbol = "Q_1"
        self.side_2.heat.latex_symbol = "Q_2"
        self.delta_temperature.latex_symbol = "\\Delta T"

    def initialize(
        self,
        state_args_1=None,
        state_args_2=None,
        outlvl=idaeslog.NOTSET,
        solver="ipopt",
        optarg={"tol": 1e-6},
        duty=None,
    ):
        """
        Heat exchanger initialization method.

        Args:
            state_args_1 : a dict of arguments to be passed to the property
                initialization for side_1 (see documentation of the specific
                property package) (default = {}).
            state_args_2 : a dict of arguments to be passed to the property
                initialization for side_2 (see documentation of the specific
                property package) (default = {}).
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default={'tol': 1e-6})
            solver : str indicating which solver to use during
                     initialization (default = 'ipopt')
            duty : an initial guess for the amount of heat transfered. This
                should be a tuple in the form (value, units), (default
                = (1000 J/s))

        Returns:
            None

        """
        # Set solver options
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        opt = SolverFactory(solver)
        opt.options = optarg
        flags1 = self.side_1.initialize(
            outlvl=outlvl, optarg=optarg, solver=solver, state_args=state_args_1
        )

        init_log.info_high("Initialization Step 1a (side_1) Complete.")

        flags2 = self.side_2.initialize(
            outlvl=outlvl, optarg=optarg, solver=solver, state_args=state_args_2
        )

        init_log.info_high("Initialization Step 1b (side_2) Complete.")
        # ---------------------------------------------------------------------
        # Solve unit without heat transfer equation
        # if costing block exists, deactivate
        try:
            self.costing.deactivate()
        except AttributeError:
            pass

        self.heat_transfer_equation.deactivate()

        # Get side 1 and side 2 heat units, and convert duty as needed
        s1_units = self.side_1.heat.get_units()
        s2_units = self.side_2.heat.get_units()

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

        self.side_2.heat.fix(s2_duty)
        for i in self.side_1.heat:
            self.side_1.heat[i].value = s1_duty

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        init_log.info_high("Initialization Step 2 {}.".format(idaeslog.condition(res)))
        self.side_2.heat.unfix()
        self.heat_transfer_equation.activate()
        # ---------------------------------------------------------------------
        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        init_log.info_high("Initialization Step 3 {}.".format(idaeslog.condition(res)))
        # ---------------------------------------------------------------------
        # Release Inlet state
        self.side_1.release_state(flags1, outlvl=outlvl)
        self.side_2.release_state(flags2, outlvl=outlvl)

        init_log.info("Initialization Completed, {}".format(idaeslog.condition(res)))
        # if costing block exists, activate
        try:
            self.costing.activate()
        except AttributeError:
            pass

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

    def get_costing(self, module=costing, hx_type='U-tube',
                    Mat_factor='stainless steel/stainless steel',
                    length_factor='12ft', year=None):
        if not hasattr(self.flowsheet(), "costing"):
            self.flowsheet().get_costing(year=year)

        self.costing = Block()
        module.hx_costing(self.costing, hx_type=hx_type, Mat_factor=Mat_factor,
                          length_factor=length_factor)

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
            iscale.constraint_scaling_transform(c, sf_dT1[t]*sf_u[t]*sf_a)

        for t, c in self.unit_heat_balance.items():
            iscale.constraint_scaling_transform(c, sf_dT1[t]*sf_u[t]*sf_a)

        for t,c in self.delta_temperature_in_equation.items():
            iscale.constraint_scaling_transform(c, sf_dT1[t])

        for t,c in self.delta_temperature_out_equation.items():
            iscale.constraint_scaling_transform(c, sf_dT2[t])
