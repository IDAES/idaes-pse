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
import pyomo.environ as pyo
from pyomo.common.config import ConfigValue, In
from idaes.core import declare_process_block_class
from idaes.models_extra.power_generation.unit_models.balance import BalanceBlockData
from idaes.core.util import from_json, to_json, StoreSpec
from idaes.core.solvers import get_solver
import idaes.models.properties.helmholtz.helmholtz as hltz
from idaes.models.properties.helmholtz.helmholtz import (
    HelmholtzThermoExpressions as ThermoExpr,
)
import idaes.logger as idaeslog
import idaes.core.util.scaling as iscale
from idaes.core.util.exceptions import ConfigurationError
from enum import Enum


_log = idaeslog.getLogger(__name__)


class ValveFunctionType(Enum):
    linear = 1
    quick_opening = 2
    equal_percentage = 3
    custom = 4


def _assert_properties(pb):
    """Assert that the properies parameter block conforms to the requirements"""
    try:
        assert isinstance(pb, hltz.HelmholtzParameterBlockData)
        assert pb.config.phase_presentation in {
            hltz.PhaseType.MIX,
            hltz.PhaseType.L,
            hltz.PhaseType.G,
        }
        assert pb.config.state_vars == hltz.StateVars.PH
    except AssertionError:
        _log.error(
            "helm.HelmValve requires a Helmholtz EOS with "
            "a single or mixed phase and pressure-enthalpy state vars."
        )
        raise


def _linear_callback(blk):
    @blk.Expression(blk.flowsheet().time)
    def valve_function(b, t):
        return b.valve_opening[t]


def _quick_open_callback(blk):
    @blk.Expression(blk.flowsheet().time)
    def valve_function(b, t):
        return pyo.sqrt(b.valve_opening[t])


def _equal_percentage_callback(blk):
    blk.alpha = pyo.Var(initialize=1, doc="Valve function parameter")
    blk.alpha.fix()

    @blk.Expression(blk.flowsheet().time)
    def valve_function(b, t):
        return b.alpha ** (b.valve_opening[t] - 1)


def _liquid_pressure_flow_rule(b, t):
    """
    For liquid F = Cv*sqrt(Pi**2 - Po**2)*f(x)
    """
    Po = b.control_volume.properties_out[t].pressure
    Pi = b.control_volume.properties_in[t].pressure
    F = b.control_volume.properties_in[t].flow_mol
    fun = b.valve_function[t]
    return F**2 == b.Cv**2 * (Pi - Po) * fun**2


def _vapor_pressure_flow_rule(b, t):
    """
    For vapor F = Cv*sqrt(Pi**2 - Po**2)*f(x)
    """
    Po = b.control_volume.properties_out[t].pressure
    Pi = b.control_volume.properties_in[t].pressure
    F = b.control_volume.properties_in[t].flow_mol
    fun = b.valve_function[t]
    return F**2 == b.Cv**2 * (Pi**2 - Po**2) * fun**2


@declare_process_block_class("HelmValve")
class HelmValveData(BalanceBlockData):
    """
    Basic adiabatic 0D valve model.  This inherits the balance block to get
    a lot of unit model boilerplate and the mass balance, enegy balance and
    pressure equations.  This model is intended to be used only with Helmholtz
    EOS property pacakges in mixed or single phase mode with P-H state vars.

    Since this inherits BalanceBlockData, and only operates in steady-state or
    pseudo-steady-state (for dynamic models) the following mass, energy and
    pressure equations are implicitly writen.

    1) Mass Balance:
        0 = flow_mol_in[t] - flow_mol_out[t]
    2) Energy Balance:
        0 = (flow_mol[t]*h_mol[t])_in - (flow_mol[t]*h_mol[t])_out
    3) Pressure:
        0 = P_in[t] + deltaP[t] - P_out[t]
    """

    CONFIG = BalanceBlockData.CONFIG()
    # For dynamics assume pseudo-steady-state
    CONFIG.dynamic = False
    CONFIG.get("dynamic")._default = False
    CONFIG.get("dynamic")._domain = In([False])
    CONFIG.has_holdup = False
    CONFIG.get("has_holdup")._default = False
    CONFIG.get("has_holdup")._domain = In([False])
    # Rest of config to make this function like a turbine
    CONFIG.has_pressure_change = True
    CONFIG.get("has_pressure_change")._default = True
    CONFIG.get("has_pressure_change")._domain = In([True])
    CONFIG.has_work_transfer = False
    CONFIG.get("has_work_transfer")._default = False
    CONFIG.get("has_work_transfer")._domain = In([False])
    CONFIG.has_heat_transfer = False
    CONFIG.get("has_heat_transfer")._default = False
    CONFIG.get("has_heat_transfer")._domain = In([False])
    CONFIG.declare(
        "valve_function",
        ConfigValue(
            default=ValveFunctionType.linear,
            domain=In(ValveFunctionType),
            description="Valve function type, if custom provide an expression rule",
            doc="""The type of valve function, if custom provide an expression rule
with the valve_function_rule argument.
**default** - ValveFunctionType.linear
**Valid values** - {
ValveFunctionType.linear,
ValveFunctionType.quick_opening,
ValveFunctionType.equal_percentage,
ValveFunctionType.custom}""",
        ),
    )
    CONFIG.declare(
        "valve_function_callback",
        ConfigValue(
            default=None,
            description="This is a callback that adds a valve function.  The "
            "callback function takes the valve bock data argument.",
        ),
    )
    CONFIG.declare(
        "phase",
        ConfigValue(
            default="Vap",
            domain=In(("Vap", "Liq")),
            description='Expected phase of fluid in valve in {"Liq", "Vap"}',
        ),
    )

    def build(self):
        """
        Add model equations to the unit model.  This is called by a default block
        construnction rule when the unit model is created.
        """
        super().build()  # Basic unit model build/read config
        config = self.config  # shorter config pointer

        # The thermodynamic expression writer object, te, writes expressions
        # including external function calls to calculate thermodynamic quantities
        # from a set of state variables.
        _assert_properties(config.property_package)
        te = ThermoExpr(blk=self, parameters=config.property_package)

        self.valve_opening = pyo.Var(
            self.flowsheet().time,
            initialize=1,
            doc="Fraction open for valve from 0 to 1",
        )
        self.Cv = pyo.Var(
            initialize=0.1,
            doc="Valve flow coefficent, for vapor " "[mol/s/Pa] for liquid [mol/s/Pa]",
            units=pyo.units.mol / pyo.units.s / pyo.units.Pa,
        )
        # self.Cv.fix()

        # set up the valve function rule.  I'm not sure these matter too much
        # for us, but the options are easy enough to provide.
        vfcb = self.config.valve_function_callback
        vfselect = self.config.valve_function
        if vfselect is not ValveFunctionType.custom and vfcb is not None:
            _log.warning(
                f"A valve function callback was provided but the valve "
                "function type is not custom."
            )

        if vfselect == ValveFunctionType.linear:
            _linear_callback(self)
        elif vfselect == ValveFunctionType.quick_opening:
            _quick_open_callback(self)
        elif vfselect == ValveFunctionType.equal_percentage:
            _equal_percentage_callback(self)
        else:
            if vfcb is None:
                raise ConfigurationError("No custom valve function callback provided")
            vfcb(self)

        if self.config.phase == "Liq":
            rule = _liquid_pressure_flow_rule
        else:
            rule = _vapor_pressure_flow_rule

        self.pressure_flow_equation = pyo.Constraint(self.flowsheet().time, rule=rule)

    def _get_performance_contents(self, time_point=0):
        """This returns a dictionary of quntities to be used in IDAES unit model
        report generation routines.
        """
        pc = super()._get_performance_contents(time_point=time_point)
        return pc

    def initialize_build(
        self,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
        calculate_cv=False,
        calculate_opening=False,
    ):
        """
        For simplicity this initialization requires you to set values for the
        efficency, inlet, and one of pressure ratio, pressure change or outlet
        pressure.
        """
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")
        init_log.info("Steam valve intialization started")

        # Create solver
        opt = get_solver(solver, optarg)

        # Store original specification so initialization doesn't change the model
        # This will only resore the values of varaibles that were originally fixed
        sp = StoreSpec.value_isfixed_isactive(only_fixed=True)
        istate = to_json(self, return_dict=True, wts=sp)
        # Check for alternate pressure specs
        for t in self.flowsheet().time:
            if self.outlet.pressure[t].fixed:
                self.deltaP[t].fix(
                    pyo.value(self.outlet.pressure[t] - self.inlet.pressure[t])
                )
                self.outlet.pressure[t].unfix()
            elif self.deltaP[t].fixed:
                # No outlet pressure specified guess a small pressure drop
                self.outlet.pressure[t] = pyo.value(
                    self.inlet.pressure[t] + self.deltaP[t]
                )

        self.inlet.fix()
        self.outlet.unfix()
        for t, v in self.deltaP.items():
            if calculate_cv:
                self.Cv.unfix()
            elif calculate_opening:
                self.valve_opening.unfix()
            elif v.fixed and self.pressure_flow_equation.active:
                self.inlet.flow_mol[t].unfix()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)

        init_log.info("Steam valve intialization complete")

        from_json(self, sd=istate, wts=sp)

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        for t, c in self.pressure_flow_equation.items():
            s = iscale.get_scaling_factor(self.control_volume.properties_in[t].flow_mol)
            s = s**2
            iscale.constraint_scaling_transform(c, s, overwrite=False)
