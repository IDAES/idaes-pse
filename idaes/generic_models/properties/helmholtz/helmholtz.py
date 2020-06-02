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
"""Generic Helmholtz EOS Base Class
"""
__author__ = "John Eslick"

# Import Python libraries
import os
import enum

# Import Pyomo libraries
from pyomo.environ import (
    Constraint,
    Expression,
    Param,
    PositiveReals,
    Set,
    value,
    Var,
    NonNegativeReals,
    ConcreteModel,
    Suffix,
)
from pyomo.environ import ExternalFunction as EF
from pyomo.core.kernel.component_set import ComponentSet
from pyomo.common.config import ConfigValue, In

# Import IDAES
from idaes.core import (
    StateBlock,
    StateBlockData,
    PhysicalParameterBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    LiquidPhase,
    VaporPhase,
    Phase,
    Component
)
from idaes.core.util.math import smooth_max
from idaes.core.util.exceptions import ConfigurationError
import idaes
import idaes.logger as idaeslog

# Logger
_log = idaeslog.getLogger(__name__)


def _available(shared_lib):
    """Make sure the compiled library functions are available. Yes, in Windows
    the .so extension is still used.
    """
    return os.path.isfile(shared_lib)

#python external function object name to ASL user function name
_external_function_map = {
    "func_p": "p",
    "func_u": "u",
    "func_s": "s",
    "func_h": "h",
    "func_hvpt": "hvpt",
    "func_hlpt": "hlpt",
    "func_svpt": "svpt",
    "func_slpt": "slpt",
    "func_uvpt": "uvpt",
    "func_ulpt": "ulpt",
    "func_tau": "tau",
    "memo_test_tau": "memo_test_tau",
    "func_tau_sp": "tau_sp",
    "func_tau_up": "tau_up",
    "func_p_stau": "p_stau",
    "func_vf": "vf",
    "func_vfs": "vfs",
    "func_vfu": "vfu",
    "func_g": "g",
    "func_f": "f",
    "func_cv": "cv",
    "func_cp": "cp",
    "func_w": "w",
    "func_delta_liq": "delta_liq",
    "func_delta_vap": "delta_vap",
    "func_delta_sat_l": "delta_sat_l",
    "func_delta_sat_v": "delta_sat_v",
    "func_p_sat": "p_sat",
    "func_tau_sat": "tau_sat",
    "func_phi0": "phi0",
    "func_phi0_delta": "phi0_delta",
    "func_phi0_delta2": "phi0_delta2",
    "func_phi0_tau": "phi0_tau",
    "func_phi0_tau2": "phi0_tau2",
    "func_phir": "phir",
    "func_phir_delta": "phir_delta",
    "func_phir_delta2": "phir_delta2",
    "func_phir_tau": "phir_tau",
    "func_phir_tau2": "phir_tau2",
    "func_phir_delta_tau": "phir_delta_tau",
}


def _available(shared_lib):
    """Make sure the compiled library functions are available. Yes, in Windows
    the .so extension is still used.
    """
    return os.path.isfile(shared_lib)


def _add_external_functions(blk, eos_tag, shared_lib, names=None):
    """Create ExternalFunction components."""
    def _fnc(x):
        x = "_".join([x, eos_tag])
        return x
    if names is None:
        names = _external_function_map.keys()
    for name in names:
        if hasattr(blk, name):
            continue
        setattr(
            blk,
            name,
            EF(library=shared_lib, function=_fnc(_external_function_map[name]))
        )


class StateVars(enum.Enum):
    """
    State variable set options
    """

    PH = 1  # Pressure-Enthalpy
    TPX = 2  # Temperature-Pressure-Quality


class PhaseType(enum.Enum):
    """
    Ways to present phases to the framework
    """

    MIX = 1  # Looks like a single phase called mixed with a vapor fraction
    LG = 2  # Looks like two phases vapor and liquid
    L = 3  # Assume only liquid is present
    G = 4  # Assume only vapor is pressent


def _htpx(T, prop=None, P=None, x=None, Tmin=200, Tmax=1200, Pmin=1, Pmax=1e9):
    """
    Convenience function to calculate enthalpy from temperature and either
    pressure or vapor fraction. This function can be used for inlet streams and
    initialization where temperature is known instead of enthalpy.

    User must provide values for one (and only one) of arguments P and x.

    Args:
        T: Temperature [K] (between Tmin and Tmax)
        Tmin: Lower bound on allowed temperatures
        Tmax: Upper bound on allowed temperatures
        Pmin: Lower bound on allowed pressures
        PmaxL Upper bound on allowed pressures
        prop: Property block to use for the enthalpy calcuations
        P: Pressure [Pa] (between Pmin and Pmax), None if saturated
        x: Vapor fraction [mol vapor/mol total] (between 0 and 1), None if
        superheated or subcooled

    Returns:
        Total molar enthalpy [J/mol].
    """
    if not (P is None) ^ (x is None):
        raise ConfigurationError(
            "htpx must be provided with one (and only one) of arguments P and x."
        )
    if not Tmin <= T <= Tmax:
        raise ConfigurationError("T = {}, ({} <= T <= {})".format(T, Tmin, Tmax))
    if P is not None and not Pmin <= P <= Pmax:
        raise ConfigurationError("P = {}, ({} <= P <= {})".format(P, Pmin, Pmax))
    if x is not None and not 0 <= x <= 1:
        raise ConfigurationError("x = {}, (0 <= x <= 1)".format(x))

    model = ConcreteModel()
    model.param = prop.config.parameters
    model.prop = prop
    Tc = model.param.temperature_crit
    Pc = model.param.pressure_crit

    if x is None:
        Tsat = Tc / value(prop.func_tau_sat(P / 1000))
        if value(T) < Tsat or value(P) > Pc:  # liquid
            return value(prop.func_hlpt(P / 1000, Tc / T) * prop.mw * 1000.0)
        else:  # vapor
            return value(prop.func_hvpt(P / 1000, Tc / T) * prop.mw * 1000.0)
    if P is None:
        Psat = value(prop.func_p_sat(Tc / T))  # kPa
        return (
            value(prop.func_hlpt(Psat, Tc / T) * prop.mw * 1000.0) * (1 - x)
            + value(prop.func_hvpt(Psat, Tc / T) * prop.mw * 1000.0) * x
        )


class HelmholtzThermoExpressions(object):
    """Class to write thermodynamic property expressions.  Take one of these
    possible sets of state variables: {h, p}, {u, p}, {s, p}, {s, T}, {T, x},
    {P, x}, or {T, P, x}, and return an expression for thermo property.
    This works by converting the given state varaibles and writing expressions
    for liquid and vapor density, vapor fraction, and temerature.  Then those
    can be used to calculate any other property.  You can specify an integer
    value 1 or 0 for x with any state varaibles to get liquid or vapor properies.
    """
    def __init__(self, blk, parameters):
        self.param = parameters
        self.blk=blk

    @staticmethod
    def _sv_str(**kwargs):
        a = [x for x in kwargs if kwargs[x] is not None]
        return ", ".join(a)

    def add_funcs(self, names=None):
        _add_external_functions(
            self.blk,
            eos_tag=self.param.eos_tag,
            shared_lib=self.param.plib,
            names=names
        )

    def basic_calculations(self, h=None, s=None, p=None, T=None, u=None, x=None):
        """This function is called as this basis for most thermo expression writer
        functions.  It takes the given state variables and returns expressions for
        liqid density, vapor density, vapor fraction and temperature, which can be
        used to write an expression for any thermo quantity.
        """
        mw = self.param.mw
        # 1.) convert units to those expected by external functions
        if h is not None:
            h = h/mw/1000 # J/mol -> kJ/kg
        if u is not None:
            u = u/mw/1000 # J/mol -> kJ/kg
        if s is not None:
            s = s/mw/1000 # J/mol/K -> kJ/kg/K
        if p is not None:
            p = p/1000 # Pa -> kPa
        if T is not None:
            tau = self.param.temperature_crit/T

        # 2.) find the block with the external functions
        blk = self.blk

        # 3.) Take given state varaibles and convert to density, T, and x

        if h is not None and p is not None:
            # h, p
            self.add_funcs(names=["func_tau", "func_vf"])
            tau = blk.func_tau(h, p)
            if x is None:
                x = blk.func_vf(h, p)
        elif s is not None and p is not None:
            # s, p
            self.add_funcs(names=["func_tau_sp", "func_vfs"])
            tau = blk.func_tau_sp(s, p)
            if x is None:
                x = blk.func_vfs(s, p)
        elif u is not None and p is not None:
            # u, p
            self.add_funcs(names=["func_tau_up", "func_vfu"])
            tau = blk.func_tau_up(u, p)
            if x is None:
                x = blk.func_vfu(u, p)
        elif s is not None and T is not None:
            # s, p
            self.add_funcs(names=["func_p_stau", "func_vfs"])
            p = blk.func_p_stau(s, tau)
            if x is None:
                x = blk.func_vfs(s, p)
        elif x is not None and T is not None and p is not None:
            # T, P, x (okay, but I hope you know what you're doing)
            pass
        elif x is not None and p is not None:
            # x, p
            self.add_funcs(names=["func_tau_sat"])
            tau = blk.func_tau_sat(p)
        elif x is not None and T is not None:
            # x, T
            self.add_funcs(names=["func_p_sat"])
            p = blk.func_p_sat(tau)
        else:
            m = "This choice of state variables ({}) are not yet supported.".format(
                self._sv_str(h=h, s=s, p=p, T=T, u=u, x=x, y=y)
            )
            _log.error(m)
            raise NotImplementedError(m)

        # 4.) Calculate density
        self.add_funcs(names=["func_delta_liq", "func_delta_vap"])
        delta_liq = blk.func_delta_liq(p, tau)
        delta_vap = blk.func_delta_vap(p, tau)

        # 5.) From here its straight forward to calculate any property
        return blk, delta_liq, delta_vap, tau, x

    def s(self, **kwargs):
        blk, delta_liq, delta_vap, tau, x = self.basic_calculations(**kwargs)
        self.add_funcs(names=["func_s"])
        s = blk.func_s(delta_liq, tau)*(1-x) + blk.func_s(delta_vap, tau)*x
        s = s*self.param.mw*1000 # kJ/kg/K -> J/mol/K
        return s

    def h(self, **kwargs):
        blk, delta_liq, delta_vap, tau, x = self.basic_calculations(**kwargs)
        self.add_funcs(names=["func_h"])
        h = blk.func_h(delta_liq, tau)*(1-x) + blk.func_h(delta_vap, tau)*x
        h = h*self.param.mw*1000 # kJ/kg -> J/mol
        return h

    def u(self, **kwargs):
        blk, delta_liq, delta_vap, tau, x = self.basic_calculations(**kwargs)
        self.add_funcs(names=["func_u"])
        u = blk.func_u(delta_liq, tau)*(1-x) + blk.func_u(delta_vap, tau)*x
        u = u*self.param.mw*1000 # kJ/kg -> J/mol
        return u

    def g(self, **kwargs):
        blk, delta_liq, delta_vap, tau, x = self.basic_calculations(**kwargs)
        self.add_funcs(names=["func_g"])
        g = blk.func_g(delta_liq, tau)*(1-x) + blk.func_g(delta_vap, tau)*x
        g = g*self.param.mw*1000 # kJ/kg -> J/mol
        return g

    def f(self, **kwargs):
        blk, delta_liq, delta_vap, tau, x = self.basic_calculations(**kwargs)
        self.add_funcs(names=["func_f"])
        f = blk.func_f(delta_liq, tau)*(1-x) + blk.func_f(delta_vap, tau)*x
        f = f*self.param.mw*1000 # kJ/kg -> J/mol
        return f

    def p(self, **kwargs):
        blk, delta_liq, delta_vap, tau, x = self.basic_calculations(**kwargs)
        self.add_funcs(names=["func_p"])
        # The following line looks a bit weird, but it is okay.  When in the
        # two-phase region the pressure for both phases is the same
        p = blk.func_p(delta_liq, tau)*(1-x) + blk.func_p(delta_vap, tau)*x
        p = p*1000 # kPa -> Pa
        return p

    def v(self, **kwargs):
        blk, delta_liq, delta_vap, tau, x = self.basic_calculations(**kwargs)
        self.add_funcs(names=["func_p"])
        v = ((1-x)/delta_liq + x/delta_vap)/self.param.dens_mass_crit*self.mw
        return v

    def x(self, **kwargs):
        blk, delta_liq, delta_vap, tau, x = self.basic_calculations(**kwargs)
        return x

    def T(self, **kwargs):
        blk, delta_liq, delta_vap, tau, x = self.basic_calculations(**kwargs)
        return self.param.temperature_crit/tau

    def tau(self, **kwargs):
        blk, delta_liq, delta_vap, tau, x = self.basic_calculations(**kwargs)
        return tau

    def delta_liq(self, **kwargs):
        blk, delta_liq, delta_vap, tau, x = self.basic_calculations(**kwargs)
        return delta_liq

    def rho_liq(self, **kwargs):
        blk, delta_liq, delta_vap, tau, x = self.basic_calculations(**kwargs)
        return delta_liq*self.param.dens_mass_crit

    def rho_mol_liq(self, **kwargs):
        blk, delta_liq, delta_vap, tau, x = self.basic_calculations(**kwargs)
        return delta_liq*self.param.dens_mass_crit/self.param.mw

    def delta_vap(self, **kwargs):
        blk, delta_liq, delta_vap, tau, x = self.basic_calculations(**kwargs)
        return delta_vap

    def rho_vap(self, **kwargs):
        blk, delta_liq, delta_vap, tau, x = self.basic_calculations(**kwargs)
        return delta_vap*self.param.dens_mass_crit

    def rho_mol_vap(self, **kwargs):
        blk, delta_liq, delta_vap, tau, x = self.basic_calculations(**kwargs)
        return delta_vap*self.param.dens_mass_crit/self.param.mw


class HelmholtzParameterBlockData(PhysicalParameterBlock):
    CONFIG = PhysicalParameterBlock.CONFIG()

    CONFIG.declare(
        "phase_presentation",
        ConfigValue(
            default=PhaseType.MIX,
            domain=In(PhaseType),
            description="Set the way phases are presented to models",
            doc="""Set the way phases are presented to models. The MIX option
appears to the framework to be a mixed phase containing liquid and/or vapor.
The mixed option can simplify calculations at the unit model level since it can
be treated as a single phase, but unit models such as flash vessels will not be
able to treat the phases independently. The LG option presents as two separate
phases to the framework. The L or G options can be used if it is known for sure
that only one phase is present.
**default** - PhaseType.MIX
**Valid values:** {
**PhaseType.MIX** - Present a mixed phase with liquid and/or vapor,
**PhaseType.LG** - Present a liquid and vapor phase,
**PhaseType.L** - Assume only liquid can be present,
**PhaseType.G** - Assume only vapor can be present}""",
        ),
    )

    CONFIG.declare(
        "state_vars",
        ConfigValue(
            default=StateVars.PH,
            domain=In(StateVars),
            description="State variable set",
            doc="""The set of state variables to use. Depending on the use, one
state variable set or another may be better computationally. Usually pressure
and enthalpy are the best choice because they are well behaved during a phase
change.
**default** - StateVars.PH
**Valid values:** {
**StateVars.PH** - Pressure-Enthalpy,
**StateVars.TPX** - Temperature-Pressure-Quality}""",
        ),
    )
    def _set_parameters(
        self,
        library,
        state_block_class,
        component_list,
        phase_equilibrium_idx,
        phase_equilibrium_list,
        mw,
        temperature_crit,
        pressure_crit,
        dens_mass_crit,
        specific_gas_constant,
        pressure_bounds,
        temperature_bounds,
        enthalpy_bounds,
        eos_tag,
        pressure_value=1e5,
        temperature_value=300,
        enthalpy_value=1000,
    ):
        """This function sets the parameters that are required for a Helmholtz
        equation of state parameter block, and ensures that all required parameters
        are set.

        """
        # Location of the *.so or *.dll file for external functions
        self.plib = library
        self.eos_tag = eos_tag
        self._state_block_class = state_block_class
        self.component_list = component_list
        self.phase_equilibrium_idx = phase_equilibrium_idx
        self.phase_equilibrium_list = phase_equilibrium_list
        # Parameters, these should match what's in the C code
        self.temperature_crit = temperature_crit
        self.pressure_crit = pressure_crit
        self.dens_mass_crit = dens_mass_crit
        self.specific_gas_constant = specific_gas_constant
        self.mw = mw
        self.default_pressure_bounds = pressure_bounds
        self.default_enthalpy_bounds = enthalpy_bounds
        self.default_temperature_bounds = temperature_bounds
        self.default_pressure_value = pressure_value
        self.default_temperature_value = temperature_value
        self.default_enthalpy_value = enthalpy_value

    def build(self):
        super().build()
        # Location of the *.so or *.dll file for external functions
        # Phase list
        self.available = _available(self.plib)

        # Create Component objects
        for c in self.component_list:
            setattr(self, str(c), Component(default={"_component_list_exists": True}))

        # Create Phase objects
        self.private_phase_list = Set(initialize=["Vap", "Liq"])
        if self.config.phase_presentation == PhaseType.MIX:
            self.Mix = Phase()

        if self.config.phase_presentation == PhaseType.LG or \
                self.config.phase_presentation == PhaseType.L:
            self.Liq = LiquidPhase()

        if self.config.phase_presentation == PhaseType.LG or \
                self.config.phase_presentation == PhaseType.G:
            self.Vap = VaporPhase()

        # State var set
        self.state_vars = self.config.state_vars

        self.smoothing_pressure_over = Param(
            mutable=True, initialize=1e-4, doc="Smooth max parameter (pressure over)"
        )
        self.smoothing_pressure_under = Param(
            mutable=True, initialize=1e-4, doc="Smooth max parameter (pressure under)"
        )

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties(
            {
                "temperature_crit": {"method": None, "units": "K"},
                "pressure_crit": {"method": None, "units": "Pa"},
                "dens_mass_crit": {"method": None, "units": "kg/m^3"},
                "specific_gas_constant": {"method": None, "units": "J/kg.K"},
                "mw": {"method": None, "units": "kg/mol"},
                "temperature_sat": {"method": "None", "units": "K"},
                "flow_mol": {"method": None, "units": "mol/s"},
                "flow_mass": {"method": None, "units": "kg/s"},
                "flow_vol": {"method": None, "units": "m^3/s"},
                "temperature": {"method": None, "units": "K"},
                "pressure": {"method": None, "units": "Pa"},
                "vapor_frac": {"method": None, "units": None},
                "dens_mass_phase": {"method": None, "units": "kg/m^3"},
                "temperature_red": {"method": None, "units": None},
                "pressure_sat": {"method": None, "units": "kPa"},
                "energy_internal_mol_phase": {"method": None, "units": "J/mol"},
                "enth_mol_phase": {"method": None, "units": "J/mol"},
                "entr_mol_phase": {"method": None, "units": "J/mol.K"},
                "cp_mol_phase": {"method": None, "units": "J/mol.K"},
                "cv_mol_phase": {"method": None, "units": "J/mol.K"},
                "speed_sound_phase": {"method": None, "units": "m/s"},
                "dens_mol_phase": {"method": None, "units": "mol/m^3"},
                "therm_cond_phase": {"method": None, "units": "W/m.K"},
                "visc_d_phase": {"method": None, "units": "Pa.s"},
                "visc_k_phase": {"method": None, "units": "m^2/s"},
                "phase_frac": {"method": None, "units": None},
                "flow_mol_comp": {"method": None, "units": "mol/s"},
                "energy_internal_mol": {"method": None, "units": "J/mol"},
                "enth_mol": {"method": None, "units": "J/mol"},
                "entr_mol": {"method": None, "units": "J/mol.K"},
                "cp_mol": {"method": None, "units": "J/mol.K"},
                "cv_mol": {"method": None, "units": "J/mol.K"},
                "heat_capacity_ratio": {"method": None, "units": None},
                "dens_mass": {"method": None, "units": "kg/m^3"},
                "dens_mol": {"method": None, "units": "mol/m^3"},
                "dh_vap_mol": {"method": None, "units": "J/mol"},
            }
        )

        obj.add_default_units(
            {
                "time": "s",
                "length": "m",
                "mass": "kg",
                "amount": "mol",
                "temperature": "K",
                "energy": "J",
                "holdup": "mol",
            }
        )


class _StateBlock(StateBlock):
    """
    This class contains methods which should be applied to Property Blocks as a
    whole, rather than individual elements of indexed Property Blocks.
    """

    @staticmethod
    def _set_fixed(v, f):
        if f:
            v.fix()
        else:
            v.unfix()

    def initialize(self, *args, **kwargs):
        flags = {}
        hold_state = kwargs.pop("hold_state", False)
        state_args = kwargs.pop("state_args", None)

        for i, v in self.items():
            pp = self[i].config.parameters.config.phase_presentation
            if self[i].state_vars == StateVars.PH:
                # hold the P-H vars
                flags[i] = (v.flow_mol.fixed, v.enth_mol.fixed, v.pressure.fixed)

                if state_args is not None:
                    if not v.flow_mol.fixed:
                        try:
                            v.flow_mol.value = state_args["flow_mol"]
                        except KeyError:
                            pass
                    if not v.enth_mol.fixed:
                        try:
                            v.enth_mol.value = state_args["enth_mol"]
                        except KeyError:
                            pass
                    if not v.pressure.fixed:
                        try:
                            v.pressure.value = state_args["pressure"]
                        except KeyError:
                            pass

                if hold_state:
                    v.flow_mol.fix()
                    v.enth_mol.fix()
                    v.pressure.fix()

            elif self[i].state_vars == StateVars.TPX:
                # Hold the T-P-x vars
                if pp in (PhaseType.MIX, PhaseType.LG):
                    flags[i] = (
                        v.flow_mol.fixed,
                        v.temperature.fixed,
                        v.pressure.fixed,
                        v.vapor_frac.fixed,
                    )

                    if state_args is not None:
                        if not v.flow_mol.fixed:
                            try:
                                v.flow_mol.value = state_args["flow_mol"]
                            except KeyError:
                                pass
                        if not v.temperature.fixed:
                            try:
                                v.temperature.value = state_args["temperature"]
                            except KeyError:
                                pass
                        if not v.pressure.fixed:
                            try:
                                v.pressure.value = state_args["pressure"]
                            except KeyError:
                                pass
                        if not v.vapor_frac.fixed:
                            try:
                                v.vapor_frac.value = state_args["vapor_frac"]
                            except KeyError:
                                pass

                    if hold_state:
                        v.flow_mol.fix()
                        v.temperature.fix()
                        v.pressure.fix()
                        v.vapor_frac.fix()
                else:
                    flags[i] = (v.flow_mol.fixed, v.temperature.fixed, v.pressure.fixed)

                    if state_args is not None:
                        if not v.flow_mol.fixed:
                            try:
                                v.flow_mol.value = state_args["flow_mol"]
                            except KeyError:
                                pass
                        if not v.temperature.fixed:
                            try:
                                v.temperature.value = state_args["temperature"]
                            except KeyError:
                                pass
                        if not v.pressure.fixed:
                            try:
                                v.pressure.value = state_args["pressure"]
                            except KeyError:
                                pass

                    if hold_state:
                        v.flow_mol.fix()
                        v.temperature.fix()
                        v.pressure.fix()

        # Call initialize on each data element
        for i in self:
            self[i].initialize(*args, **kwargs)
        return flags

    def release_state(self, flags, **kwargs):
        for i, f in flags.items():
            pp = self[i].config.parameters.config.phase_presentation
            if self[i].state_vars == StateVars.PH:
                self._set_fixed(self[i].flow_mol, f[0])
                self._set_fixed(self[i].enth_mol, f[1])
                self._set_fixed(self[i].pressure, f[2])
            elif self[i].state_vars == StateVars.TPX:
                self._set_fixed(self[i].flow_mol, f[0])
                self._set_fixed(self[i].temperature, f[1])
                self._set_fixed(self[i].pressure, f[2])
                if pp in (PhaseType.MIX, PhaseType.LG):
                    self._set_fixed(self[i].vapor_frac, f[3])


class HelmholtzStateBlockData(StateBlockData):
    """
    This is a base clase for Helmholtz equations of state using IDAES standard
    Helmholtz EOS external functions written in C++.
    """
    def initialize(self, *args, **kwargs):
        # With this particualr property pacakage there is not need for
        # initialization
        pass

    def _external_functions(self):
        """Create ExternalFunction components.  This includes some external
        functions that are not usually used for testing purposes."""
        _add_external_functions(
            blk=self,
            eos_tag=self.config.parameters.eos_tag,
            shared_lib=self.config.parameters.plib
        )

    def _state_vars(self):
        """ Create the state variables
        """
        params = self.config.parameters

        self.flow_mol = Var(
            initialize=1,
            doc="Total flow [mol/s]"
        )
        self.scaling_factor[self.flow_mol] = 1e-3

        if self.state_vars == StateVars.PH:
            self.pressure = Var(
                domain=PositiveReals,
                initialize=params.default_pressure_value,
                doc="Pressure [Pa]",
                bounds=params.default_pressure_bounds,
            )
            self.enth_mol = Var(
                initialize=params.default_enthalpy_value,
                doc="Total molar enthalpy (J/mol)",
                bounds=params.default_enthalpy_bounds,
            )
            self.scaling_factor[self.enth_mol] = 1e-3

            P = self.pressure / 1000.0  # Pressure expr [kPA] (for external func)
            h_mass = self.enth_mol / self.mw / 1000  # enthalpy expr [kJ/kg]
            phase_set = params.config.phase_presentation

            self.temperature = Expression(
                expr=self.temperature_crit / self.func_tau(h_mass, P),
                doc="Temperature (K)",
            )
            if phase_set == PhaseType.MIX or phase_set == PhaseType.LG:
                self.vapor_frac = Expression(
                    expr=self.func_vf(h_mass, P),
                    doc="Vapor mole fraction (mol vapor/mol total)",
                )
            elif phase_set == PhaseType.L:
                self.vapor_frac = Expression(
                    expr=0.0,
                    doc="Vapor mole fraction (mol vapor/mol total)",
                )
            elif phase_set == PhaseType.G:
                self.vapor_frac = Expression(
                    expr=1.0,
                    doc="Vapor mole fraction (mol vapor/mol total)",
                )

            # For variables that show up in ports specify extensive/intensive
            self.extensive_set = ComponentSet((self.flow_mol,))
            self.intensive_set = ComponentSet((self.enth_mol, self.pressure))

        elif self.state_vars == StateVars.TPX:
            self.temperature = Var(
                domain=PositiveReals,
                initialize=params.default_temperature_value,
                doc="Temperature [K]",
                bounds=params.default_temperature_bounds
            )
            self.pressure = Var(
                domain=PositiveReals,
                initialize=params.default_pressure_value,
                doc="Pressure [Pa]",
                bounds=params.default_pressure_bounds,
            )
            self.vapor_frac = Var(
                initialize=0.0,
                doc="Vapor fraction [none]"
                # No bounds here, since it is often (usually) on it's bound
                # and that's not the best for IPOPT
            )

            # enth_mol is defined later, since in this case it needs
            # enth_mol_phase to be defined first

            # For variables that show up in ports specify extensive/intensive
            self.extensive_set = ComponentSet((self.flow_mol,))
            self.intensive_set = ComponentSet(
                (self.temperature, self.pressure, self.vapor_frac)
            )
        self.scaling_factor[self.temperature] = 1e-1
        self.scaling_factor[self.pressure] = 1e-6
        self.scaling_factor[self.vapor_frac] = 1e1

    def _tpx_phase_eq(self):
        # Saturation pressure
        eps_pu = self.config.parameters.smoothing_pressure_under
        eps_po = self.config.parameters.smoothing_pressure_over
        priv_plist = self.config.parameters.private_phase_list
        plist = self.config.parameters.phase_list
        rhoc = self.config.parameters.dens_mass_crit

        P = self.pressure / 1000  # expression for pressure in kPa
        Psat = self.pressure_sat / 1000.0  # expression for Psat in kPA
        vf = self.vapor_frac
        tau = self.tau

        # Terms for determining if you are above, below, or at the Psat
        self.P_under_sat = Expression(
            expr=smooth_max(0, Psat - P, eps_pu),
            doc="pressure above Psat, 0 if liqid exists [kPa]",
        )
        self.P_over_sat = Expression(
            expr=smooth_max(0, P - Psat, eps_po),
            doc="pressure below Psat, 0 if vapor exists [kPa]",
        )

        # Calculate liquid and vapor density.  If the phase doesn't exist,
        # density will be calculated at the saturation or critical pressure
        def rule_dens_mass(b, p):
            if p == "Liq":
                self.scaling_factor[self.dens_mass_phase[p]] = 1e-2
                return rhoc * self.func_delta_liq(P + self.P_under_sat, tau)
            else:
                self.scaling_factor[self.dens_mass_phase[p]] = 1e1
                return rhoc * self.func_delta_vap(P - self.P_over_sat, tau)

        self.dens_mass_phase = Expression(priv_plist, rule=rule_dens_mass)

        # Reduced Density (no _mass_ identifier because mass or mol is same)
        def rule_dens_red(b, p):
            return self.dens_mass_phase[p] / rhoc

        self.dens_phase_red = Expression(
            priv_plist, rule=rule_dens_red, doc="reduced density [unitless]"
        )

        # If there is only one phase fix the vapor fraction appropriately
        if len(plist) == 1:
            if "Vap" in plist:
                self.vapor_frac.fix(1.0)
            else:
                self.vapor_frac.fix(0.0)
        elif not self.config.defined_state:
            self.eq_complementarity = Constraint(
                expr=0 == (vf * self.P_over_sat - (1 - vf) * self.P_under_sat)
            )
            self.scaling_expression[self.eq_complementarity] = 10 / self.pressure

        # eq_sat can activated to force the pressure to be the saturation
        # pressure, if you use this constraint deactivate eq_complementarity
        self.eq_sat = Constraint(expr=P / 1000.0 == Psat / 1000.0)
        self.scaling_expression[self.eq_sat] = 1000 / self.pressure
        self.eq_sat.deactivate()


    def build(self, *args):
        """
        Callable method for Block construction
        """
        super().build(*args)

        # Create the scaling suffixes for the state block
        self.scaling_factor = Suffix(direction=Suffix.EXPORT)
        self.scaling_expression = Suffix()

        # Check if the library is available.
        self.available = self.config.parameters.available
        if not self.available:
            _log.error("Library file '{}' not found. Was it installed?".format(
                    self.config.parameters.plib
                )
            )

        # Add external functions
        self._external_functions()

        # Thermo expression writer
        te = HelmholtzThermoExpressions(
            blk=self,
            parameters=self.config.parameters
        )

        # Which state vars to use
        self.state_vars = self.config.parameters.state_vars
        # The private phase list contains phases that may be present and is
        # used internally.  If using the single mixed phase option the phase
        # list would be mixed while the private phase list would be ["Liq, "Vap"]
        phlist = self.config.parameters.private_phase_list
        pub_phlist = self.config.parameters.phase_list
        component_list = self.config.parameters.component_list
        phase_set = self.config.parameters.config.phase_presentation
        self.phase_equilibrium_list = self.config.parameters.phase_equilibrium_list

        # Expressions that link to some parameters in the param block, which
        # are commonly needed, this lets you get the parameters with scale
        # factors directly from the state block
        self.temperature_crit = Expression(expr=self.config.parameters.temperature_crit)
        self.scaling_factor[self.temperature_crit] = 1e-2
        self.pressure_crit = Expression(expr=self.config.parameters.pressure_crit)
        self.scaling_factor[self.pressure_crit] = 1e-6
        self.dens_mass_crit = Expression(expr=self.config.parameters.dens_mass_crit)
        self.scaling_factor[self.dens_mass_crit] = 1e-2
        self.mw = Expression(
            expr=self.config.parameters.mw, doc="molecular weight [kg/mol]"
        )
        self.scaling_factor[self.mw] = 1e3

        # create the appropriate state variables
        self._state_vars()

        # Some parameters/variables show up in several expressions, so to enhance
        # readability and compactness, give them short aliases
        Tc = self.config.parameters.temperature_crit
        rhoc = self.config.parameters.dens_mass_crit
        mw = self.mw
        P = self.pressure / 1000.0  # Pressure expr [kPA] (for external func)
        T = self.temperature
        vf = self.vapor_frac

        # Saturation temperature expression
        self.temperature_sat = Expression(
            expr=Tc / self.func_tau_sat(P), doc="Stauration temperature (K)"
        )
        self.scaling_factor[self.temperature_sat] = 1e-2

        # Saturation tau (tau = Tc/T)
        self.tau_sat = Expression(expr=self.func_tau_sat(P))

        # Reduced temperature
        self.temperature_red = Expression(
            expr=T / Tc, doc="reduced temperature T/Tc (unitless)"
        )
        self.scaling_factor[self.temperature_red] = 1

        self.tau = Expression(expr=Tc / T, doc="Tc/T (unitless)")
        tau = self.tau

        # Saturation pressure
        self.pressure_sat = Expression(
            expr=1000 * self.func_p_sat(tau), doc="Saturation pressure (Pa)"
        )
        self.scaling_factor[self.pressure_sat] = 1e-5

        if self.state_vars == StateVars.PH:
            # If TPx state vars the expressions are given in _tpx_phase_eq
            # Calculate liquid and vapor density.  If the phase doesn't exist,
            # density will be calculated at the saturation or critical pressure
            # depending on whether the temperature is above the critical
            # temperature supercritical fluid is considered to be the liquid
            # phase
            def rule_dens_mass(b, p):
                if p == "Liq":
                    self.scaling_factor[self.dens_mass_phase[p]] = 1e-2
                    return rhoc * self.func_delta_liq(P, tau)
                else:
                    self.scaling_factor[self.dens_mass_phase[p]] = 1e1
                    return rhoc * self.func_delta_vap(P, tau)

            self.dens_mass_phase = Expression(
                phlist, rule=rule_dens_mass, doc="Mass density by phase (kg/m3)"
            )

            # Reduced Density (no _mass_ identifier as mass or mol is same)
            def rule_dens_red(b, p):
                self.scaling_factor[self.dens_phase_red[p]] = 1
                return self.dens_mass_phase[p] / rhoc

            self.dens_phase_red = Expression(
                phlist, rule=rule_dens_red, doc="reduced density (unitless)"
            )

        elif self.state_vars == StateVars.TPX:
            self._tpx_phase_eq()
        delta = self.dens_phase_red

        # Phase property expressions all converted to SI

        # Saturated Enthalpy
        def rule_enth_mol_sat_phase(b, p):
            if p == "Liq":
                self.scaling_factor[self.enth_mol_sat_phase[p]] = 1e-2
                return 1000 * mw * self.func_hlpt(P, self.tau_sat)
            elif p == "Vap":
                self.scaling_factor[self.enth_mol_sat_phase[p]] = 1e-4
                return 1000 * mw * self.func_hvpt(P, self.tau_sat)

        self.enth_mol_sat_phase = Expression(
            phlist,
            rule=rule_enth_mol_sat_phase,
            doc="Saturated enthalpy of the phases at pressure (J/mol)",
        )

        self.dh_vap_mol = Expression(
            expr=self.enth_mol_sat_phase["Vap"] - self.enth_mol_sat_phase["Liq"],
            doc="Enthaply of vaporization at pressure and saturation (J/mol)",
        )
        self.scaling_factor[self.dh_vap_mol] = 1e-4

        # Phase Internal Energy
        def rule_energy_internal_mol_phase(b, p):
            if p == "Liq":
                self.scaling_factor[self.energy_internal_mol_phase[p]] = 1e-2
            else:
                self.scaling_factor[self.energy_internal_mol_phase[p]] = 1e-4
            return 1000 * mw * self.func_u(delta[p], tau)

        self.energy_internal_mol_phase = Expression(
            phlist,
            rule=rule_energy_internal_mol_phase,
            doc="Phase internal energy or saturated if phase doesn't exist [J/mol]",
        )

        # Phase Enthalpy
        def rule_enth_mol_phase(b, p):
            if p == "Liq":
                self.scaling_factor[self.enth_mol_phase[p]] = 1e-2
            elif p == "Vap":
                self.scaling_factor[self.enth_mol_phase[p]] = 1e-4
            return 1000 * mw * self.func_h(delta[p], tau)

        self.enth_mol_phase = Expression(
            phlist,
            rule=rule_enth_mol_phase,
            doc="Phase enthalpy or saturated if phase doesn't exist [J/mol]",
        )

        # Phase Entropy
        def rule_entr_mol_phase(b, p):
            if p == "Liq":
                self.scaling_factor[self.entr_mol_phase[p]] = 1e-1
            elif p == "Vap":
                self.scaling_factor[self.entr_mol_phase[p]] = 1e-1
            return 1000 * mw * self.func_s(delta[p], tau)

        self.entr_mol_phase = Expression(
            phlist,
            rule=rule_entr_mol_phase,
            doc="Phase entropy or saturated if phase doesn't exist [J/mol/K]",
        )

        # Phase constant pressure heat capacity, cp
        def rule_cp_mol_phase(b, p):
            if p == "Liq":
                self.scaling_factor[self.cp_mol_phase[p]] = 1e-2
            elif p == "Vap":
                self.scaling_factor[self.cp_mol_phase[p]] = 1e-2
            return 1000 * mw * self.func_cp(delta[p], tau)

        self.cp_mol_phase = Expression(
            phlist,
            rule=rule_cp_mol_phase,
            doc="Phase cp or saturated if phase doesn't exist [J/mol/K]",
        )

        # Phase constant pressure heat capacity, cv
        def rule_cv_mol_phase(b, p):
            if p == "Liq":
                self.scaling_factor[self.cv_mol_phase[p]] = 1e-2
            elif p == "Vap":
                self.scaling_factor[self.cv_mol_phase[p]] = 1e-2
            return 1000 * mw * self.func_cv(delta[p], tau)

        self.cv_mol_phase = Expression(
            phlist,
            rule=rule_cv_mol_phase,
            doc="Phase cv or saturated if phase doesn't exist [J/mol/K]",
        )

        # Phase speed of sound
        def rule_speed_sound_phase(b, p):
            if p == "Liq":
                self.scaling_factor[self.speed_sound_phase[p]] = 1e-2
            elif p == "Vap":
                self.scaling_factor[self.speed_sound_phase[p]] = 1e-2
            return self.func_w(delta[p], tau)

        self.speed_sound_phase = Expression(
            phlist,
            rule=rule_speed_sound_phase,
            doc="Phase speed of sound or saturated if phase doesn't exist [m/s]",
        )

        # Phase Mole density
        def rule_dens_mol_phase(b, p):
            if p == "Liq":
                self.scaling_factor[self.dens_mol_phase[p]] = 1e-2
            elif p == "Vap":
                self.scaling_factor[self.dens_mol_phase[p]] = 1e-4
            return self.dens_mass_phase[p] / mw

        self.dens_mol_phase = Expression(
            phlist,
            rule=rule_dens_mol_phase,
            doc="Phase mole density or saturated if phase doesn't exist [mol/m3]",
        )

        # Phase fraction
        def rule_phase_frac(b, p):
            self.scaling_factor[self.phase_frac[p]] = 10
            if p == "Vap":
                return vf
            elif p == "Liq":
                return 1.0 - vf

        self.phase_frac = Expression(
            phlist, rule=rule_phase_frac, doc="Phase fraction [unitless]"
        )

        # Component flow (for units that need it)
        def component_flow(b, i):
            self.scaling_factor[self.flow_mol_comp[i]] = 1e-3
            return self.flow_mol

        self.flow_mol_comp = Expression(
            component_list,
            rule=component_flow,
            doc="Total flow (both phases) of component [mol/s]",
        )

        # Total (mixed phase) properties

        # Enthalpy
        if self.state_vars == StateVars.TPX:
            self.enth_mol = Expression(
                expr=sum(self.phase_frac[p] * self.enth_mol_phase[p] for p in phlist)
            )
            self.scaling_factor[self.enth_mol] = 1e-3
        # Internal Energy
        self.energy_internal_mol = Expression(
            expr=sum(
                self.phase_frac[p] * self.energy_internal_mol_phase[p] for p in phlist
            )
        )
        self.scaling_factor[self.energy_internal_mol] = 1e-3
        # Entropy
        self.entr_mol = Expression(expr=te.s(h=self.enth_mol, p=self.pressure))
        self.scaling_factor[self.entr_mol] = 1e-1
        # cp
        self.cp_mol = Expression(
            expr=sum(self.phase_frac[p] * self.cp_mol_phase[p] for p in phlist)
        )
        self.scaling_factor[self.cp_mol] = 1e-2
        # cv
        self.cv_mol = Expression(
            expr=sum(self.phase_frac[p] * self.cv_mol_phase[p] for p in phlist)
        )
        self.scaling_factor[self.cv_mol] = 1e-2
        # mass density
        self.dens_mass = Expression(
            expr=1.0
            / sum(self.phase_frac[p] * 1.0 / self.dens_mass_phase[p] for p in phlist)
        )
        self.scaling_factor[self.dens_mass] = 1e0
        # mole density
        self.dens_mol = Expression(
            expr=1.0
            / sum(self.phase_frac[p] * 1.0 / self.dens_mol_phase[p] for p in phlist)
        )
        self.scaling_factor[self.dens_mol] = 1e-3
        # heat capacity ratio
        self.heat_capacity_ratio = Expression(expr=self.cp_mol / self.cv_mol)
        self.scaling_factor[self.heat_capacity_ratio] = 1e1
        # Flows
        self.flow_vol = Expression(
            expr=self.flow_mol / self.dens_mol,
            doc="Total liquid + vapor volumetric flow (m3/s)",
        )
        self.scaling_factor[self.flow_vol] = 100

        self.flow_mass = Expression(
            expr=self.mw * self.flow_mol, doc="mass flow rate [kg/s]"
        )
        self.scaling_factor[self.flow_mass] = 1

        self.enth_mass = Expression(expr=self.enth_mol / mw, doc="Mass enthalpy (J/kg)")
        self.scaling_factor[self.enth_mass] = 1

        # Set the state vars dictionary
        if self.state_vars == StateVars.PH:
            self._state_vars_dict = {
                "flow_mol": self.flow_mol,
                "enth_mol": self.enth_mol,
                "pressure": self.pressure,
            }
        elif self.state_vars == StateVars.TPX and phase_set in (
            PhaseType.MIX,
            PhaseType.LG,
        ):
            self._state_vars_dict = {
                "flow_mol": self.flow_mol,
                "temperature": self.temperature,
                "pressure": self.pressure,
                "vapor_frac": self.vapor_frac,
            }
        elif self.state_vars == StateVars.TPX and phase_set in (
            PhaseType.G,
            PhaseType.L,
        ):
            self._state_vars_dict = {
                "flow_mol": self.flow_mol,
                "temperature": self.temperature,
                "pressure": self.pressure,
            }

        # Define some expressions for the balance terms returned by functions
        # This is just to allow assigning scale factors to the expressions
        # returned
        #
        # Marterial flow term exprsssions
        def rule_material_flow_terms(b, p):
            self.scaling_expression[b.material_flow_terms[p]] = 1 / self.flow_mol
            if p == "Mix":
                return self.flow_mol
            else:
                return self.flow_mol * self.phase_frac[p]

        self.material_flow_terms = Expression(pub_phlist, rule=rule_material_flow_terms)

        # Enthaply flow term expressions
        def rule_enthalpy_flow_terms(b, p):
            if p == "Mix":
                self.scaling_expression[b.enthalpy_flow_terms[p]] = 1 / (
                    self.enth_mol * self.flow_mol
                )
                return self.enth_mol * self.flow_mol
            else:
                self.scaling_expression[b.enthalpy_flow_terms[p]] = 1 / (
                    self.enth_mol_phase[p] * self.phase_frac[p] * self.flow_mol
                )
                return self.enth_mol_phase[p] * self.phase_frac[p] * self.flow_mol

        self.enthalpy_flow_terms = Expression(pub_phlist, rule=rule_enthalpy_flow_terms)

        # Energy density term expressions
        def rule_energy_density_terms(b, p):
            if p == "Mix":
                self.scaling_expression[b.energy_density_terms[p]] = 1 / (
                    self.energy_internal_mol * self.flow_mol
                )
                return self.dens_mol * self.energy_internal_mol
            else:
                self.scaling_expression[b.energy_density_terms[p]] = 1 / (
                    self.dens_mol_phase[p] * self.energy_internal_mol_phase[p]
                )
                return self.dens_mol_phase[p] * self.energy_internal_mol_phase[p]

        self.energy_density_terms = Expression(
            pub_phlist, rule=rule_energy_density_terms
        )

    def get_material_flow_terms(self, p, j):
        return self.material_flow_terms[p]

    def get_enthalpy_flow_terms(self, p):
        return self.enthalpy_flow_terms[p]

    def get_material_density_terms(self, p, j):
        if p == "Mix":
            return self.dens_mol
        else:
            return self.dens_mol_phase[p]

    def get_energy_density_terms(self, p):
        return self.energy_density_terms[p]

    def default_material_balance_type(self):
        return MaterialBalanceType.componentTotal

    def default_energy_balance_type(self):
        return EnergyBalanceType.enthalpyTotal

    def define_state_vars(self):
        return self._state_vars_dict

    def define_display_vars(self):
        return {
            "Molar Flow (mol/s)": self.flow_mol,
            "Mass Flow (kg/s)": self.flow_mass,
            "T (K)": self.temperature,
            "P (Pa)": self.pressure,
            "Vapor Fraction": self.vapor_frac,
            "Molar Enthalpy (J/mol)": self.enth_mol_phase,
        }

    def extensive_state_vars(self):
        return self.extensive_set

    def intensive_state_vars(self):
        return self.intensive_set

    def model_check(self):
        pass
