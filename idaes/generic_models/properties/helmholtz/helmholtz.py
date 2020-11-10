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
    ConcreteModel,
    units as pyunits
)
from pyomo.environ import ExternalFunction as EF
from pyomo.common.collections import ComponentSet
from pyomo.common.config import ConfigValue, In

# Import IDAES
from idaes.core import (
    StateBlock,
    StateBlockData,
    PhysicalParameterBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MaterialFlowBasis,
    LiquidPhase,
    VaporPhase,
    Phase,
    Component
)
from idaes.core.util.math import smooth_max
from idaes.core.util.exceptions import ConfigurationError
import idaes
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

# Logger
_log = idaeslog.getLogger(__name__)


def _available(shared_lib):
    """Make sure the compiled library functions are available. Yes, in Windows
    the .so extension is still used.
    """
    return os.path.isfile(shared_lib)


# python external function object name to ASL user function name
_external_function_map = {
    "func_p": {"fname": "p",
               "units": pyunits.kPa,
               "arg_units": [pyunits.dimensionless, pyunits.dimensionless]},
    "func_u": {"fname": "u",
               "units": pyunits.kJ/pyunits.kg,
               "arg_units": [pyunits.dimensionless, pyunits.dimensionless]},
    "func_s": {"fname": "s",
               "units": pyunits.kJ/pyunits.kg/pyunits.K,
               "arg_units": [pyunits.dimensionless, pyunits.dimensionless]},
    "func_h": {"fname": "h",
               "units": pyunits.kJ/pyunits.kg,
               "arg_units": [pyunits.dimensionless, pyunits.dimensionless]},
    "func_hvpt": {"fname": "hvpt",
                  "units": pyunits.kJ/pyunits.kg,
                  "arg_units": [pyunits.kPa, pyunits.dimensionless]},
    "func_hlpt": {"fname": "hlpt",
                  "units": pyunits.kJ/pyunits.kg,
                  "arg_units": [pyunits.kPa, pyunits.dimensionless]},
    "func_svpt": {"fname": "svpt", "units": None, "arg_units": [None, None]},
    "func_slpt": {"fname": "slpt", "units": None, "arg_units": [None, None]},
    "func_uvpt": {"fname": "uvpt", "units": None, "arg_units": [None, None]},
    "func_ulpt": {"fname": "ulpt", "units": None, "arg_units": [None, None]},
    "func_tau": {"fname": "tau",
                 "units": pyunits.dimensionless,
                 "arg_units": [pyunits.kJ/pyunits.kg, pyunits.kPa]},
    "memo_test_tau": {"fname": "memo_test_tau", "units": None, "arg_units": [None, None]},
    "func_tau_sp": {"fname": "tau_sp",
                    "units": pyunits.dimensionless,
                    "arg_units": [pyunits.kJ/pyunits.kg/pyunits.K,
                                  pyunits.kPa]},
    "func_tau_up": {"fname": "tau_up",
                    "units": pyunits.dimensionless,
                    "arg_units": [pyunits.kJ/pyunits.kg, pyunits.kPa]},
    "func_p_stau": {"fname": "p_stau",
                    "units": pyunits.kPa,
                    "arg_units": [pyunits.kJ/pyunits.kg/pyunits.K,
                                  pyunits.dimensionless]},
    "func_vf": {"fname": "vf",
                "units": pyunits.dimensionless,
                "arg_units": [pyunits.kJ/pyunits.kg, pyunits.kPa]},
    "func_vfs": {"fname": "vfs",
                 "units": pyunits.dimensionless,
                 "arg_units": [pyunits.kJ/pyunits.kg/pyunits.K, pyunits.kPa]},
    "func_vfu": {"fname": "vfu",
                 "units": pyunits.dimensionless,
                 "arg_units": [pyunits.kJ/pyunits.kg, pyunits.kPa]},
    "func_g": {"fname": "g",
               "units": pyunits.kJ/pyunits.kg,
               "arg_units": [pyunits.dimensionless, pyunits.dimensionless]},
    "func_f": {"fname": "f",
               "units": pyunits.kJ/pyunits.kg,
               "arg_units": [pyunits.dimensionless, pyunits.dimensionless]},
    "func_cv": {"fname": "cv",
                "units": pyunits.kJ/pyunits.kg/pyunits.K,
                "arg_units": [pyunits.dimensionless, pyunits.dimensionless]},
    "func_cp": {"fname": "cp",
                "units": pyunits.kJ/pyunits.kg/pyunits.K,
                "arg_units": [pyunits.dimensionless, pyunits.dimensionless]},
    "func_w": {"fname": "w",
               "units": pyunits.m/pyunits.s,
               "arg_units": [pyunits.dimensionless, pyunits.dimensionless]},
    "func_delta_liq": {"fname": "delta_liq",
                       "units": pyunits.dimensionless,
                       "arg_units": [pyunits.kPa, None]},
    "func_delta_vap": {"fname": "delta_vap",
                       "units": pyunits.dimensionless,
                       "arg_units": [pyunits.kPa, None]},
    "func_delta_sat_l": {"fname": "delta_sat_l", "units": None, "arg_units": [None, None]},
    "func_delta_sat_v": {"fname": "delta_sat_v", "units": None, "arg_units": [None, None]},
    "func_p_sat": {"fname": "p_sat",
                   "units": pyunits.kPa,
                   "arg_units": [pyunits.dimensionless]},
    "func_tau_sat": {"fname": "tau_sat",
                     "units": pyunits.dimensionless,
                     "arg_units": [pyunits.kPa]},
    "func_phi0": {"fname": "phi0", "units": None, "arg_units": [None, None]},
    "func_phi0_delta": {"fname": "phi0_delta", "units": None, "arg_units": [None, None]},
    "func_phi0_delta2": {"fname": "phi0_delta2", "units": None, "arg_units": [None, None]},
    "func_phi0_tau": {"fname": "phi0_tau", "units": None, "arg_units": [None, None]},
    "func_phi0_tau2": {"fname": "phi0_tau2", "units": None, "arg_units": [None, None]},
    "func_phir": {"fname": "phir", "units": None, "arg_units": [None, None]},
    "func_phir_delta": {"fname": "phir_delta", "units": None, "arg_units": [None, None]},
    "func_phir_delta2": {"fname": "phir_delta2", "units": None, "arg_units": [None, None]},
    "func_phir_tau": {"fname": "phir_tau", "units": None, "arg_units": [None, None]},
    "func_phir_tau2": {"fname": "phir_tau2", "units": None, "arg_units": [None, None]},
    "func_phir_delta_tau": {"fname": "phir_delta_tau", "units": None, "arg_units": [None, None]},
}


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
            EF(library=shared_lib,
               function=_fnc(_external_function_map[name]["fname"]),
               units=_external_function_map[name]["units"],
               arg_units=_external_function_map[name]["arg_units"])
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


def _htpx(T, prop=None, P=None, x=None,
          Tmin=200, Tmax=1200, Pmin=1e-3, Pmax=1e6):
    """
    Convenience function to calculate enthalpy from temperature and either
    pressure or vapor fraction. This function can be used for inlet streams and
    initialization where temperature is known instead of enthalpy.

    User must provide values for one (and only one) of arguments P and x.

    Args:
        T: Temperature (between Tmin and Tmax)
        Tmin: Lower bound on allowed temperatures [K]
        Tmax: Upper bound on allowed temperatures [K]
        Pmin: Lower bound on allowed pressures [kPa]
        PmaxL Upper bound on allowed pressures [kPa]
        prop: Property block to use for the enthalpy calcuations
        P: Pressure (between Pmin and Pmax), None if saturated
        x: Vapor fraction [mol vapor/mol total] (between 0 and 1), None if
        superheated or subcooled

    Returns:
        Total molar enthalpy [J/mol].
    """
    if not (P is None) ^ (x is None):
        raise ConfigurationError("htpx must be provided with one (and only "
                                 "one) of arguments P and x.")

    T = pyunits.convert(T, to_units=pyunits.K)
    if not Tmin <= value(T) <= Tmax:
        raise ConfigurationError("T = {}, ({} <= T <= {} [K])"
                                 .format(value(T), Tmin, Tmax))
    if P is not None:
        P = pyunits.convert(P, to_units=pyunits.kPa)
        if not (Pmin <= value(P) <= Pmax):
            raise ConfigurationError("P = {}, ({} <= P <= {} [kPa])"
                                     .format(value(P), Pmin, Pmax))
    if x is not None and not 0 <= x <= 1:
        raise ConfigurationError("x = {}, (0 <= x <= 1)".format(x))

    model = ConcreteModel()
    model.param = prop.config.parameters
    model.prop = prop
    Tc = model.param.temperature_crit
    Pc = model.param.pressure_crit

    if x is None:
        Tsat = Tc / prop.func_tau_sat(P)
        if value(T) < value(Tsat) or value(P) > value(Pc):  # liquid
            return value(pyunits.convert(prop.func_hlpt(P, Tc/T) * prop.mw,
                                         to_units=pyunits.J/pyunits.mol))
        else:  # vapor
            return value(pyunits.convert(prop.func_hvpt(P, Tc/T) * prop.mw,
                                         to_units=pyunits.J/pyunits.mol))
    if P is None:
        Psat = prop.func_p_sat(Tc/T)  # kPa
        return (
            value(pyunits.convert(
                prop.func_hlpt(Psat, Tc/T) * prop.mw * (1 - x) +
                prop.func_hvpt(Psat, Tc/T) * prop.mw * x,
                to_units=pyunits.J/pyunits.mol)))


class HelmholtzThermoExpressions(object):
    """Class to write thermodynamic property expressions.  Take one of these
    possible sets of state variables: {h, p}, {u, p}, {s, p}, {s, T}, {T, x},
    {P, x}, or {T, P, x}, and return an expression for thermo property.
    This works by converting the given state varaibles and writing expressions
    for liquid and vapor density, vapor fraction, and temerature.  Then those
    can be used to calculate any other property.  You can specify an integer
    value 1 or 0 for x with any state varaibles to get liquid or vapor
    properies.
    """

    def __init__(self, blk, parameters):
        self.param = parameters
        self.blk = blk

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

    def basic_calculations(
            self, h=None, s=None, p=None, T=None, u=None, x=None):
        """This function is called as this basis for most thermo expression
        writer functions.  It takes the given state variables and returns
        expressions for liqid density, vapor density, vapor fraction and
        temperature, which can be used to write an expression for any thermo
        quantity.
        """
        mw = self.param.mw
        # 1.) convert units to those expected by external functions
        if h is not None:
            h = pyunits.convert(h/mw, to_units=pyunits.kJ/pyunits.kg)
        if u is not None:
            u = pyunits.convert(u/mw, to_units=pyunits.kJ/pyunits.kg)
        if s is not None:
            s = pyunits.convert(s/mw, to_units=pyunits.kJ/pyunits.kg/pyunits.K)
        if p is not None:
            p = pyunits.convert(p, to_units=pyunits.kPa)
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
            m = ("This choice of state variables ({}) is not yet supported."
                 .format(self._sv_str(h=h, s=s, p=p, T=T, u=u, x=x)))
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
        s = pyunits.convert(s*self.param.mw,
                            to_units=pyunits.J/pyunits.mol/pyunits.K)
        return s

    def h(self, **kwargs):
        blk, delta_liq, delta_vap, tau, x = self.basic_calculations(**kwargs)
        self.add_funcs(names=["func_h"])
        h = blk.func_h(delta_liq, tau)*(1-x) + blk.func_h(delta_vap, tau)*x
        h = pyunits.convert(h*self.param.mw,
                            to_units=pyunits.J/pyunits.mol)
        return h

    def u(self, **kwargs):
        blk, delta_liq, delta_vap, tau, x = self.basic_calculations(**kwargs)
        self.add_funcs(names=["func_u"])
        u = blk.func_u(delta_liq, tau)*(1-x) + blk.func_u(delta_vap, tau)*x
        u = pyunits.convert(u*self.param.mw,
                            to_units=pyunits.J/pyunits.mol)
        return u

    def g(self, **kwargs):
        blk, delta_liq, delta_vap, tau, x = self.basic_calculations(**kwargs)
        self.add_funcs(names=["func_g"])
        g = blk.func_g(delta_liq, tau)*(1-x) + blk.func_g(delta_vap, tau)*x
        g = pyunits.convert(g*self.param.mw,
                            to_units=pyunits.J/pyunits.mol)
        return g

    def f(self, **kwargs):
        blk, delta_liq, delta_vap, tau, x = self.basic_calculations(**kwargs)
        self.add_funcs(names=["func_f"])
        f = blk.func_f(delta_liq, tau)*(1-x) + blk.func_f(delta_vap, tau)*x
        f = pyunits.convert(f*self.param.mw,
                            to_units=pyunits.J/pyunits.mol)
        return f

    def p(self, **kwargs):
        blk, delta_liq, delta_vap, tau, x = self.basic_calculations(**kwargs)
        self.add_funcs(names=["func_p"])
        # The following line looks a bit weird, but it is okay.  When in the
        # two-phase region the pressure for both phases is the same
        p = blk.func_p(delta_liq, tau)*(1-x) + blk.func_p(delta_vap, tau)*x
        p = pyunits.convert(p, to_units=pyunits.Pa)
        return p

    def v(self, **kwargs):
        blk, delta_liq, delta_vap, tau, x = self.basic_calculations(**kwargs)
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
        equation of state parameter block, and ensures that all required
        parameters are set.
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
            setattr(self, str(c), Component(
                default={"_component_list_exists": True}))

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
            mutable=True,
            initialize=1e-4,
            doc="Smooth max parameter (pressure over)",
            units=pyunits.kPa)

        self.smoothing_pressure_under = Param(
            mutable=True,
            initialize=1e-4,
            doc="Smooth max parameter (pressure under)",
            units=pyunits.kPa)

        # TODO<jce> leaving flow to not break things, but plan to remove
        self.set_default_scaling("flow_mol", 1e-4)
        self.set_default_scaling("flow_mol_comp", 1e-4)
        self.set_default_scaling("flow_vol", 100)
        self.set_default_scaling("flow_mass", 1)

        # Set some scalings with reasonable a priori values
        self.set_default_scaling("temperature_crit", 1e-2)
        self.set_default_scaling("enth_mol", 1e-3)
        self.set_default_scaling("temperature", 1e-1)
        self.set_default_scaling("pressure", 1e-6)
        self.set_default_scaling("vapor_frac", 1e1)
        self.set_default_scaling("dens_mass_phase", 1e-2, index="Liq")
        self.set_default_scaling("dens_mass_phase", 1e1)  # Not liq
        self.set_default_scaling("pressure_crit", 1e-6)
        self.set_default_scaling("dens_mass_crit", 1e-2)
        self.set_default_scaling("mw", 1e3)
        self.set_default_scaling("temperature_sat", 1e-2)
        self.set_default_scaling("temperature_red", 1)
        self.set_default_scaling("pressure_sat", 1e-5)
        self.set_default_scaling("dens_mass_phase", 1e-2, index="Liq")
        self.set_default_scaling("dens_mass_phase", 1e1)
        self.set_default_scaling("dens_phase_red", 1)
        self.set_default_scaling("enth_mol_sat_phase", 1e-2, index="Liq")
        self.set_default_scaling("enth_mol_sat_phase", 1e-4, index="Vap")
        self.set_default_scaling("dh_vap_mol", 1e-4)
        self.set_default_scaling(
            "energy_internal_mol_phase", 1e-2, index="Liq")
        self.set_default_scaling("energy_internal_mol_phase", 1e-4)
        self.set_default_scaling("enth_mol_phase", 1e-2, index="Liq")
        self.set_default_scaling("enth_mol_phase", 1e-4, index="Vap")
        self.set_default_scaling("entr_mol_phase", 1e-1, index="Liq")
        self.set_default_scaling("entr_mol_phase", 1e-1, index="Vap")
        self.set_default_scaling("cp_mol_phase", 1e-2, index="Liq")
        self.set_default_scaling("cp_mol_phase", 1e-2, index="Vap")
        self.set_default_scaling("cv_mol_phase", 1e-2, index="Liq")
        self.set_default_scaling("cv_mol_phase", 1e-2, index="Vap")
        self.set_default_scaling("dens_mol_phase", 1e-2, index="Liq")
        self.set_default_scaling("dens_mol_phase", 1e-4, index="Vap")
        self.set_default_scaling("phase_frac", 10)
        self.set_default_scaling("enth_mol", 1e-3)
        self.set_default_scaling("energy_internal_mol", 1e-3)
        self.set_default_scaling("entr_mol", 1e-1)
        self.set_default_scaling("cp_mol", 1e-2)
        self.set_default_scaling("cv_mol", 1e-2)
        self.set_default_scaling("dens_mass", 1)
        self.set_default_scaling("dens_mol", 1e-3)
        self.set_default_scaling("heat_capacity_ratio", 1e1)
        self.set_default_scaling("enth_mass", 1)

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
                "energy_internal_mol_phase": {
                    "method": None, "units": "J/mol"},
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
                "time": pyunits.s,
                "length": pyunits.m,
                "mass": pyunits.kg,
                "amount": pyunits.mol,
                "temperature": pyunits.K,
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
                flags[i] = (v.flow_mol.fixed,
                            v.enth_mol.fixed,
                            v.pressure.fixed)

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
                    flags[i] = (v.flow_mol.fixed,
                                v.temperature.fixed,
                                v.pressure.fixed)

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
            doc="Total flow [mol/s]",
            units=pyunits.mol/pyunits.s
        )

        if self.state_vars == StateVars.PH:
            self.pressure = Var(
                domain=PositiveReals,
                initialize=params.default_pressure_value,
                doc="Pressure [Pa]",
                bounds=params.default_pressure_bounds,
                units=pyunits.Pa)
            self.enth_mol = Var(
                initialize=params.default_enthalpy_value,
                doc="Total molar enthalpy (J/mol)",
                bounds=params.default_enthalpy_bounds,
                units=pyunits.J/pyunits.mol)

            P = pyunits.convert(self.pressure, to_units=pyunits.kPa)
            h_mass = pyunits.convert(self.enth_mol / self.mw,
                                     to_units=pyunits.kJ/pyunits.kg)
            phase_set = params.config.phase_presentation

            self.temperature = Expression(
                expr=self.temperature_crit / self.func_tau(h_mass, P),
                doc="Temperature (K)")

            if phase_set == PhaseType.MIX or phase_set == PhaseType.LG:
                self.vapor_frac = Expression(
                    expr=self.func_vf(h_mass, P),
                    doc="Vapor mole fraction (mol vapor/mol total)")
            elif phase_set == PhaseType.L:
                self.vapor_frac = Expression(
                    expr=0.0,
                    doc="Vapor mole fraction (mol vapor/mol total)")
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
                bounds=params.default_temperature_bounds,
                units=pyunits.K)

            self.pressure = Var(
                domain=PositiveReals,
                initialize=params.default_pressure_value,
                doc="Pressure [Pa]",
                bounds=params.default_pressure_bounds,
                units=pyunits.Pa)

            self.vapor_frac = Var(
                initialize=0.0,
                units=pyunits.dimensionless,
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

    def _tpx_phase_eq(self):
        # Saturation pressure
        eps_pu = self.config.parameters.smoothing_pressure_under
        eps_po = self.config.parameters.smoothing_pressure_over
        priv_plist = self.config.parameters.private_phase_list
        plist = self.config.parameters.phase_list
        rhoc = self.config.parameters.dens_mass_crit

        P = pyunits.convert(self.pressure, to_units=pyunits.kPa)
        Psat = pyunits.convert(self.pressure_sat, to_units=pyunits.kPa)
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
                return rhoc * self.func_delta_liq(P + self.P_under_sat, tau)
            else:
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

        # eq_sat can activated to force the pressure to be the saturation
        # pressure, if you use this constraint deactivate eq_complementarity
        self.eq_sat = Constraint(expr=P / 1000.0 == Psat / 1000.0)
        self.eq_sat.deactivate()

    def build(self, *args):
        """
        Callable method for Block construction
        """
        super().build(*args)

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
        # list would be mixed while the private phase list would be
        # ["Liq, "Vap"]
        phlist = self.config.parameters.private_phase_list
        pub_phlist = self.config.parameters.phase_list
        component_list = self.config.parameters.component_list
        phase_set = self.config.parameters.config.phase_presentation
        self.phase_equilibrium_list = \
            self.config.parameters.phase_equilibrium_list

        # Expressions that link to some parameters in the param block, which
        # are commonly needed, this lets you get the parameters with scale
        # factors directly from the state block
        self.temperature_crit = Expression(
            expr=self.config.parameters.temperature_crit)
        self.pressure_crit = Expression(
            expr=self.config.parameters.pressure_crit)
        self.dens_mass_crit = Expression(
            expr=self.config.parameters.dens_mass_crit)
        self.mw = Expression(
            expr=self.config.parameters.mw, doc="molecular weight [kg/mol]"
        )

        # create the appropriate state variables
        self._state_vars()

        # Some parameters/variables show up in several expressions, so to
        # enhance readability and compactness, give them short aliases
        Tc = self.config.parameters.temperature_crit
        rhoc = self.config.parameters.dens_mass_crit
        mw = self.mw
        P = pyunits.convert(self.pressure, to_units=pyunits.kPa)
        T = self.temperature
        vf = self.vapor_frac

        # Saturation temperature expression
        self.temperature_sat = Expression(
            expr=Tc / self.func_tau_sat(P), doc="Stauration temperature (K)"
        )

        # Saturation tau (tau = Tc/T)
        self.tau_sat = Expression(expr=self.func_tau_sat(P))

        # Reduced temperature
        self.temperature_red = Expression(
            expr=T / Tc, doc="reduced temperature T/Tc (unitless)"
        )

        self.tau = Expression(expr=Tc / T, doc="Tc/T (unitless)")
        tau = self.tau

        # Saturation pressure
        self.pressure_sat = Expression(
            expr=pyunits.convert(self.func_p_sat(tau), to_units=pyunits.Pa),
            doc="Saturation pressure (Pa)"
        )

        if self.state_vars == StateVars.PH:
            # If TPx state vars the expressions are given in _tpx_phase_eq
            # Calculate liquid and vapor density.  If the phase doesn't exist,
            # density will be calculated at the saturation or critical pressure
            # depending on whether the temperature is above the critical
            # temperature supercritical fluid is considered to be the liquid
            # phase
            def rule_dens_mass(b, p):
                if p == "Liq":
                    return rhoc * self.func_delta_liq(P, tau)
                else:
                    return rhoc * self.func_delta_vap(P, tau)

            self.dens_mass_phase = Expression(
                phlist,
                rule=rule_dens_mass,
                doc="Mass density by phase (kg/m3)"
            )

            # Reduced Density (no _mass_ identifier as mass or mol is same)
            def rule_dens_red(b, p):
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
                return pyunits.convert(mw * self.func_hlpt(P, self.tau_sat),
                                       to_units=pyunits.J/pyunits.mol)
            elif p == "Vap":
                return pyunits.convert(mw * self.func_hvpt(P, self.tau_sat),
                                       to_units=pyunits.J/pyunits.mol)

        self.enth_mol_sat_phase = Expression(
            phlist,
            rule=rule_enth_mol_sat_phase,
            doc="Saturated enthalpy of the phases at pressure (J/mol)",
        )

        self.dh_vap_mol = Expression(
            expr=(self.enth_mol_sat_phase["Vap"] -
                  self.enth_mol_sat_phase["Liq"]),
            doc="Enthaply of vaporization at pressure and saturation (J/mol)",
        )

        # Phase Internal Energy
        def rule_energy_internal_mol_phase(b, p):
            return pyunits.convert(mw * self.func_u(delta[p], tau),
                                   to_units=pyunits.J/pyunits.mol)

        self.energy_internal_mol_phase = Expression(
            phlist,
            rule=rule_energy_internal_mol_phase,
            doc="Phase internal energy or saturated if phase doesn't exist",
        )

        # Phase Enthalpy
        def rule_enth_mol_phase(b, p):
            return pyunits.convert(mw * self.func_h(delta[p], tau),
                                   to_units=pyunits.J/pyunits.mol)

        self.enth_mol_phase = Expression(
            phlist,
            rule=rule_enth_mol_phase,
            doc="Phase enthalpy or saturated if phase doesn't exist [J/mol]",
        )

        # Phase Entropy
        def rule_entr_mol_phase(b, p):
            return pyunits.convert(mw * self.func_s(delta[p], tau),
                                   to_units=pyunits.J/pyunits.mol/pyunits.K)

        self.entr_mol_phase = Expression(
            phlist,
            rule=rule_entr_mol_phase,
            doc="Phase entropy or saturated if phase doesn't exist [J/mol/K]",
        )

        # Phase constant pressure heat capacity, cp
        def rule_cp_mol_phase(b, p):
            return pyunits.convert(mw * self.func_cp(delta[p], tau),
                                   to_units=pyunits.J/pyunits.mol/pyunits.K)

        self.cp_mol_phase = Expression(
            phlist,
            rule=rule_cp_mol_phase,
            doc="Phase cp or saturated if phase doesn't exist [J/mol/K]",
        )

        # Phase constant pressure heat capacity, cv
        def rule_cv_mol_phase(b, p):
            return pyunits.convert(mw * self.func_cv(delta[p], tau),
                                   to_units=pyunits.J/pyunits.mol/pyunits.K)

        self.cv_mol_phase = Expression(
            phlist,
            rule=rule_cv_mol_phase,
            doc="Phase cv or saturated if phase doesn't exist [J/mol/K]",
        )

        # Phase speed of sound
        def rule_speed_sound_phase(b, p):
            return self.func_w(delta[p], tau)

        self.speed_sound_phase = Expression(
            phlist,
            rule=rule_speed_sound_phase,
            doc="Phase speed of sound or saturated if phase doesn't exist",
        )

        # Phase Mole density
        def rule_dens_mol_phase(b, p):
            return self.dens_mass_phase[p] / mw

        self.dens_mol_phase = Expression(
            phlist,
            rule=rule_dens_mol_phase,
            doc="Phase mole density or saturated if phase doesn't exist")

        # Phase fraction
        def rule_phase_frac(b, p):
            if p == "Vap":
                return vf
            elif p == "Liq":
                return 1.0 - vf

        self.phase_frac = Expression(
            phlist, rule=rule_phase_frac, doc="Phase fraction [unitless]"
        )

        # Component flow (for units that need it)
        def component_flow(b, i):
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
                expr=sum(self.phase_frac[p] * self.enth_mol_phase[p]
                         for p in phlist)
            )
        # Internal Energy
        self.energy_internal_mol = Expression(
            expr=sum(
                self.phase_frac[p] * self.energy_internal_mol_phase[p]
                for p in phlist
            )
        )
        # Entropy
        self.entr_mol = Expression(expr=te.s(h=self.enth_mol, p=self.pressure))
        # cp
        self.cp_mol = Expression(
            expr=sum(self.phase_frac[p] * self.cp_mol_phase[p]
                     for p in phlist)
        )
        # cv
        self.cv_mol = Expression(
            expr=sum(self.phase_frac[p] * self.cv_mol_phase[p] for p in phlist)
        )

        # mass density
        self.dens_mass = Expression(
            expr=1.0
            / sum(self.phase_frac[p] * 1.0 / self.dens_mass_phase[p]
                  for p in phlist)
        )
        # mole density
        self.dens_mol = Expression(
            expr=1.0
            / sum(self.phase_frac[p] * 1.0 / self.dens_mol_phase[p]
                  for p in phlist)
        )
        # heat capacity ratio
        self.heat_capacity_ratio = Expression(expr=self.cp_mol / self.cv_mol)
        # Flows
        self.flow_vol = Expression(
            expr=self.flow_mol / self.dens_mol,
            doc="Total liquid + vapor volumetric flow (m3/s)",
        )

        self.flow_mass = Expression(
            expr=self.mw * self.flow_mol, doc="mass flow rate [kg/s]"
        )

        self.enth_mass = Expression(expr=self.enth_mol / mw,
                                    doc="Mass enthalpy (J/kg)")

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

        def rule_mole_frac_phase_comp(b, p, j):
            if p == "Mix":
                return 1.0
            else:
                return self.phase_frac[p]
        self.mole_frac_phase_comp = Expression(
            pub_phlist,
            component_list,
            rule=rule_mole_frac_phase_comp
        )

        # Define some expressions for the balance terms returned by functions
        # This is just to allow assigning scale factors to the expressions
        # returned
        #
        # Marterial flow term exprsssions
        def rule_material_flow_terms(b, p):
            if p == "Mix":
                return self.flow_mol
            else:
                return self.flow_mol * self.phase_frac[p]

        self.material_flow_terms = Expression(pub_phlist,
                                              rule=rule_material_flow_terms)

        # Enthaply flow term expressions
        def rule_enthalpy_flow_terms(b, p):
            if p == "Mix":
                return self.enth_mol * self.flow_mol
            else:
                return (self.enth_mol_phase[p] * self.phase_frac[p] *
                        self.flow_mol)

        self.enthalpy_flow_terms = Expression(pub_phlist,
                                              rule=rule_enthalpy_flow_terms)

        # Energy density term expressions
        def rule_energy_density_terms(b, p):
            if p == "Mix":
                return self.dens_mol * self.energy_internal_mol
            else:
                return (self.dens_mol_phase[p] *
                        self.energy_internal_mol_phase[p])

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

    def get_material_flow_basis(b):
        return MaterialFlowBasis.molar

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

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        sf_flow = iscale.get_scaling_factor(self.flow_mol, default=1)
        sf_enth = iscale.get_scaling_factor(self.enth_mol, default=1)
        sf_inte = iscale.get_scaling_factor(self.energy_internal_mol,
                                            default=1)
        sf_dens = iscale.get_scaling_factor(self.dens_mol, default=1)
        sf_pres = iscale.get_scaling_factor(self.pressure, default=1)
        for v in self.material_flow_terms.values():
            iscale.set_scaling_factor(v, sf_flow)
        for v in self.enthalpy_flow_terms.values():
            iscale.set_scaling_factor(v, sf_enth*sf_flow)
        for k, v in self.energy_density_terms.items():
            if k == "Mix":
                iscale.set_scaling_factor(v, sf_inte*sf_dens)
            else:
                sf_inte_p = iscale.get_scaling_factor(
                    self.energy_internal_mol_phase[k], default=1)
                sf_dens_p = iscale.get_scaling_factor(
                    self.dens_mol_phase[k], default=1)
                iscale.set_scaling_factor(v, sf_inte_p*sf_dens_p)
        try:
            iscale.set_scaling_factor(self.eq_sat, sf_pres/1000.0)
        except AttributeError:
            pass  # may not have eq_sat, and that's ok
        try:
            iscale.set_scaling_factor(self.eq_complementarity, sf_pres/10)
        except AttributeError:
            pass  # may not have eq_complementarity which is fine
