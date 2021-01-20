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
This provides valve models for steam and liquid water.  These are for
steam cycle control valves and the turbine throttle valves.
"""

__Author__ = "John Eslick"

from enum import Enum

from pyomo.common.config import ConfigValue, In
from pyomo.environ import Var, Expression, value, Constraint, sqrt, Param

from idaes.core import declare_process_block_class
from idaes.generic_models.unit_models.pressure_changer import (
    PressureChangerData,
    ThermodynamicAssumption,
    MaterialBalanceType,
)
from idaes.core.util import from_json, to_json, StoreSpec
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)


class ValveFunctionType(Enum):
    linear = 1
    quick_opening = 2
    equal_percentage = 3
    custom = 4


def _define_config(config):
    config.compressor = False
    config.get("compressor")._default = False
    config.get("compressor")._domain = In([False])
    config.material_balance_type = MaterialBalanceType.componentTotal
    config.get("material_balance_type")._default = \
        MaterialBalanceType.componentTotal
    config.thermodynamic_assumption = ThermodynamicAssumption.adiabatic
    config.get("thermodynamic_assumption")._default = \
        ThermodynamicAssumption.adiabatic
    config.get("thermodynamic_assumption")._domain = In(
        [ThermodynamicAssumption.adiabatic]
    )
    config.declare(
        "valve_function",
        ConfigValue(
            default=ValveFunctionType.linear,
            domain=In(ValveFunctionType),
            description="Valve function type, if custom provide an expression rule",
            doc="""The type of valve function, if custom provide an expression
rule with the valve_function_rule argument.
**default** - ValveFunctionType.linear
**Valid values** - {
ValveFunctionType.linear,
ValveFunctionType.quick_opening,
ValveFunctionType.equal_percentage,
ValveFunctionType.custom}""",
        ),
    )
    config.declare(
        "valve_function_rule",
        ConfigValue(
            default=None,
            description="This is a rule that returns a time indexed valve function expression.",
            doc="""This is a rule that returns a time indexed valve function
expression. This is required only if valve_function==ValveFunctionType.custom""",
        ),
    )
    config.declare(
        "phase",
        ConfigValue(
            default="Vap",
            domain=In(("Vap", "Liq")),
            description='Expected phase of fluid in valve in {"Liq", "Vap"}',
        ),
    )


def _linear_rule(b, t):
    return b.valve_opening[t]


def _quick_open_rule(b, t):
    return sqrt(b.valve_opening[t])


def _equal_percentage_rule(b, t):
    return b.alpha ** (b.valve_opening[t] - 1)


def _liquid_pressure_flow_rule(b, t):
    """
    For liquid F = Cv*sqrt(Pi**2 - Po**2)*f(x)
    """
    Po = b.control_volume.properties_out[t].pressure
    Pi = b.control_volume.properties_in[t].pressure
    F = b.control_volume.properties_in[t].flow_mol
    Cv = b.Cv
    fun = b.valve_function[t]
    return ((1 / b.flow_scale ** 2) * F ** 2 ==
            (1 / b.flow_scale ** 2) * Cv ** 2 *
            (Pi - Po) * fun ** 2)


def _vapor_pressure_flow_rule(b, t):
    """
    For vapor F = Cv*sqrt(Pi**2 - Po**2)*f(x)
    """
    Po = b.control_volume.properties_out[t].pressure
    Pi = b.control_volume.properties_in[t].pressure
    F = b.control_volume.properties_in[t].flow_mol
    Cv = b.Cv
    fun = b.valve_function[t]
    return ((1 / b.flow_scale ** 2) * F ** 2 ==
            (1 / b.flow_scale ** 2) * Cv ** 2 *
            (Pi ** 2 - Po ** 2) * fun ** 2)


@declare_process_block_class("SteamValve", doc="Basic steam valve models")
class SteamValveData(PressureChangerData):
    # Same settings as the default pressure changer, but force to expander with
    # isentropic efficiency
    CONFIG = PressureChangerData.CONFIG()
    _define_config(CONFIG)

    def build(self):
        super().build()

        self.valve_opening = Var(
            self.flowsheet().config.time,
            initialize=1,
            doc="Fraction open for valve from 0 to 1",
        )

        umeta = self.config.property_package.get_metadata().get_derived_units
        if self.config.phase == "Liq":
            cv_units = umeta("amount")/umeta("time")/umeta("pressure")**0.5
        else:
            cv_units = umeta("amount")/umeta("time")/umeta("pressure")

        self.Cv = Var(
            initialize=0.1,
            doc="Valve flow coefficent",
            units=cv_units
        )
        self.flow_scale = Param(
            mutable=True,
            default=1e3,
            doc="Scaling factor for pressure flow relation should be "
            "approximatly the same order of magnitude as the expected flow.",
        )
        self.Cv.fix()
        self.valve_opening.fix()

        # set up the valve function rule.  I'm not sure these matter too much
        # for us, but the options are easy enough to provide.
        if self.config.valve_function == ValveFunctionType.linear:
            rule = _linear_rule
        elif self.config.valve_function == ValveFunctionType.quick_opening:
            rule = _quick_open_rule
        elif self.config.valve_function == ValveFunctionType.equal_percentage:
            self.alpha = Var(initialize=1, doc="Valve function parameter")
            self.alpha.fix()
            rule = _equal_percentage_rule
        else:
            rule = self.config.valve_function_rule

        self.valve_function = Expression(
            self.flowsheet().config.time,
            rule=rule,
            doc="Valve function expression"
        )

        if self.config.phase == "Liq":
            rule = _liquid_pressure_flow_rule
        else:
            rule = _vapor_pressure_flow_rule

        self.pressure_flow_equation = Constraint(
            self.flowsheet().config.time, rule=rule
        )

    def initialize(
        self,
        state_args={},
        outlvl=idaeslog.NOTSET,
        solver="ipopt",
        optarg={"tol": 1e-6, "max_iter": 30},
    ):
        """
        Initialize the turbine stage model.  This deactivates the
        specialized constraints, then does the isentropic turbine
        initialization, then reactivates the constraints and solves.

        Args:
            state_args (dict): Initial state for property initialization
            outlvl : sets output level of initialization routine
            solver (str): Solver to use for initialization
            optarg (dict): Solver arguments dictionary
        """
        # sp is what to save to make sure state after init is same as the start
        #   saves value, fixed, and active state, doesn't load originally free
        #   values, this makes sure original problem spec is same but
        #   initializes the values of free vars
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")

        sp = StoreSpec.value_isfixed_isactive(only_fixed=True)
        istate = to_json(self, return_dict=True, wts=sp)

        self.deltaP[:].unfix()
        self.ratioP[:].unfix()

        # fix inlet and free outlet
        for t in self.flowsheet().config.time:
            for k, v in self.inlet.vars.items():
                v[t].fix()
            for k, v in self.outlet.vars.items():
                v[t].unfix()
            # to calculate outlet pressure
            Pout = self.outlet.pressure[t]
            Pin = self.inlet.pressure[t]
            if self.deltaP[t].value is not None:
                prdp = value((self.deltaP[t] - Pin) / Pin)
            else:
                prdp = -100  # crazy number to say don't use deltaP as guess
            if value(Pout / Pin) > 1 or value(Pout / Pin) < 0.0:
                if value(self.ratioP[t]) <= 1 and value(self.ratioP[t]) >= 0:
                    Pout.value = value(Pin * self.ratioP[t])
                elif prdp <= 1 and prdp >= 0:
                    Pout.value = value(prdp * Pin)
                else:
                    Pout.value = value(Pin * 0.95)
            self.deltaP[t] = value(Pout - Pin)
            self.ratioP[t] = value(Pout / Pin)

        # Make sure the initialization problem has no degrees of freedom
        # This shouldn't happen here unless there is a bug in this
        dof = degrees_of_freedom(self)
        try:
            assert dof == 0
        except:
            init_log.exception("degrees_of_freedom = {}".format(dof))
            raise

        # one bad thing about reusing this is that the log messages aren't
        # really compatible with being nested inside another initialization
        super().initialize(
            state_args=state_args, outlvl=outlvl, solver=solver, optarg=optarg
        )

        # reload original spec
        from_json(self, sd=istate, wts=sp)

    def _get_performance_contents(self, time_point=0):
        pc = super()._get_performance_contents(time_point=time_point)

        pc["vars"]["Opening"] = self.valve_opening[time_point]
        pc["vars"]["Valve Coefficient"] = self.Cv
        if self.config.valve_function == ValveFunctionType.equal_percentage:
            pc["vars"]["alpha"] = self.alpha

        pc["params"] = {}
        pc["params"]["Flow Scaling"] = self.flow_scale

        return pc
