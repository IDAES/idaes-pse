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
"""IDAES IAPWS-95 Steam properties

Dropped all critical enhancments and non-analytic terms ment to improve
accruacy near the critical point. These tend to cause singularities in the
equations, and it is assumend that we will try to avoid operating very close to
the critical point.

 References: (some of this is only used in the C++ part)
   International Association for the Properties of Water and Steam (2016).
       IAPWS R6-95 (2016), "Revised Release on the IAPWS Formulation 1995 for
       the Properties of Ordinary Water Substance for General Scientific Use,"
       URL: http://iapws.org/relguide/IAPWS95-2016.pdf
   Wagner, W.,  A. Pruss (2002). "The IAPWS Formulation 1995 for the
       Thermodynamic Properties of Ordinary Water Substance for General and
       Scientific Use." J. Phys. Chem. Ref. Data, 31, 387-535.
   Wagner, W. et al. (2000). "The IAPWS Industrial Formulation 1997 for the
       Thermodynamic Properties of Water and Steam," ASME J. Eng. Gas Turbines
       and Power, 122, 150-182.
   Akasaka, R. (2008). "A Reliable and Useful Method to Determine the
       Saturation State from Helmholtz Energy Equations of State." Journal of
       Thermal Science and Technology, 3(3), 442-451.
   International Association for the Properties of Water and Steam (2011).
       IAPWS R15-11, "Release on the IAPWS Formulation 2011 for the
       Thermal Conductivity of Ordinary Water Substance,"
       URL: http://iapws.org/relguide/ThCond.pdf
   International Association for the Properties of Water and Steam (2008).
       IAPWS R12-08, "Release on the IAPWS Formulation 2008 for the Viscosity
       of Ordinary Water Substance,"
       URL: http://iapws.org/relguide/visc.pdf
"""
__author__ = "John Eslick"

# Import Python libraries
import os
import enum

# Import Pyomo libraries
from pyomo.environ import Constraint, Expression, Param, PositiveReals,\
                          RangeSet, Set, value, Var, NonNegativeReals,\
                          exp, sqrt, ConcreteModel
from pyomo.environ import ExternalFunction as EF
from pyomo.core.kernel.component_set import ComponentSet
from pyomo.common.config import ConfigValue, In

# Import IDAES
from idaes.core import declare_process_block_class, \
                       StateBlock, StateBlockData, PhysicalParameterBlock, \
                       MaterialBalanceType, EnergyBalanceType
from idaes.core.util.math import smooth_max
from idaes.core.util.exceptions import ConfigurationError
import idaes
from idaes.logger import getIdaesLogger, getInitLogger, init_tee, condition

# Logger
_log = getIdaesLogger(__name__)
_so = os.path.join(idaes.lib_directory, "iapws95_external.so")


def iapws95_available():
    """Make sure the compiled IAPWS-95 functions are available. Yes, in Windows
    the .so extention is still used.
    """
    return os.path.isfile(_so)


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


def htpx(T, P=None, x=None):
    """
    Convenience function to calculate steam enthalpy from temperature and
    either pressure or vapor fraction. This function can be used for inlet
    streams and initialization where temperature is known instead of enthalpy.

    User must provided values for one (and only one) of arguments P and x.

    Args:
        T: Temperature [K] (between 200 and 3000)
        P: Pressure [Pa] (between 1 and 1e9), None if saturated steam
        x: Vapor fraction [mol vapor/mol total] (between 0 and 1), None if
        superheated or subcooled

    Returns:
        Total molar enthalpy [J/mol].
    """
    if not (P is None) ^ (x is None):
        raise ConfigurationError(
                "htpx must be provided with one (and only one) of "
                "arguments P and x.")
    if not 200 <= T <= 3e3:
        raise ConfigurationError("T out of range. Must be between 2e2 and 3e3")
    if P is not None and not 1 <= P <= 1e9:
        raise ConfigurationError("P out of range. Must be between 1 and 1e9")
    if x is not None and not 0 <= x <= 1:
        raise ConfigurationError("x must be between 0 and 1")

    model = ConcreteModel()
    model.prop_param = Iapws95ParameterBlock()
    prop = model.prop = Iapws95StateBlock(
            default={"parameters": model.prop_param})

    if x is None:
        Tsat = 647.096/value(prop.func_tau_sat(P/1000))
        if value(T) < Tsat or value(P/1000) > 22064:  # liquid
            return value(prop.func_hlpt(P/1000, 647.096/T)*prop.mw*1000.0)
        else:  # vapor
            return value(prop.func_hvpt(P/1000, 647.096/T)*prop.mw*1000.0)
    if P is None:
        Psat = value(prop.func_p_sat(647.096/T))  # kPa
        return value(prop.func_hlpt(Psat, 647.096/T)*prop.mw*1000.0)*(1-x) +\
            value(prop.func_hvpt(Psat, 647.096/T)*prop.mw*1000.0)*x


@declare_process_block_class("Iapws95ParameterBlock")
class Iapws95ParameterBlockData(PhysicalParameterBlock):
    CONFIG = PhysicalParameterBlock.CONFIG()

    CONFIG.declare("phase_presentation", ConfigValue(
        default=PhaseType.MIX,
        domain=In(PhaseType),
        description="Set the way phases are presented to models",
        doc="""Set the way phases are presented to models. The MIX option
appears to the framework to be a mixed phase containing liquid and/or vapor.
The mixed option can simplify calculations at the unit model level since it can
be treated as a single phase, but unit models such as flash vessels will not be
able to treate the phases indepedently. The LG option presents as two sperate
phases to the framework. The L or G options can be used if it is known for sure
that only one phase is present.
**default** - PhaseType.MIX
**Valid values:** {
**PhaseType.MIX** - Present a mixed phase with liquid and/or vapor,
**PhaseType.LG** - Present a liquid and vapor phase,
**PhaseType.L** - Assume only liquid can be present,
**PhaseType.G** - Assume only vapor can be present}"""))

    CONFIG.declare("state_vars", ConfigValue(
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
**StateVars.TPX** - Temperature-Pressure-Quality}"""))

    def build(self):
        super(Iapws95ParameterBlockData, self).build()
        self.state_block_class = Iapws95StateBlock
        # Location of the *.so or *.dll file for external functions
        self.plib = _so
        self.available = os.path.isfile(self.plib)
        # Phase list
        self.private_phase_list = Set(initialize=["Vap", "Liq"])
        if self.config.phase_presentation == PhaseType.MIX:
            self.phase_list = Set(initialize=["Mix"])
        elif self.config.phase_presentation == PhaseType.LG:
            self.phase_list = Set(initialize=["Vap", "Liq"])
        elif self.config.phase_presentation == PhaseType.L:
            self.phase_list = Set(initialize=["Liq"])
        elif self.config.phase_presentation == PhaseType.G:
            self.phase_list = Set(initialize=["Vap"])

        # State var set
        self.state_vars = self.config.state_vars

        # Component list - a list of component identifiers
        self.component_list = Set(initialize=['H2O'])

        # List of phase equilibrium
        self.phase_equilibrium_idx = Set(initialize=[1])
        self.phase_equilibrium_list = {1: ["H2O", ("Vap", "Liq")]}

        # Parameters, these should match what's in the C code
        self.temperature_crit = Param(initialize=647.096,
                                      doc='Critical temperature [K]')
        self.pressure_crit = Param(initialize=2.2064e7,
                                   doc='Critical pressure [Pa]')
        self.dens_mass_crit = Param(initialize=322,
                                    doc='Critical density [kg/m3]')
        self.gas_const = Param(initialize=8.3144598,
                               doc='Gas Constant [J/mol/K]')
        self.mw = Param(initialize=0.01801528,
                        doc='Molecular weight [kg/mol]')
        # Thermal conductivity parameters.
        # "Release on the IAPWS Formulation 2011 for the Thermal Conductivity
        # of Ordinary Water Substance"
        self.tc_L0 = Param(RangeSet(0, 5), initialize={
            0: 2.443221e-3,
            1: 1.323095e-2,
            2: 6.770357e-3,
            3: -3.454586e-3,
            4: 4.096266e-4},
            doc="0th order thermal conductivity paramters")

        self.tc_L1 = Param(RangeSet(0, 5), RangeSet(0, 6), initialize={
            (0, 0): 1.60397357,
            (1, 0): 2.33771842,
            (2, 0): 2.19650529,
            (3, 0): -1.21051378,
            (4, 0): -2.7203370,
            (0, 1): -0.646013523,
            (1, 1): -2.78843778,
            (2, 1): -4.54580785,
            (3, 1): 1.60812989,
            (4, 1): 4.57586331,
            (0, 2): 0.111443906,
            (1, 2): 1.53616167,
            (2, 2): 3.55777244,
            (3, 2): -0.621178141,
            (4, 2): -3.18369245,
            (0, 3): 0.102997357,
            (1, 3): -0.463045512,
            (2, 3): -1.40944978,
            (3, 3): 0.0716373224,
            (4, 3): 1.1168348,
            (0, 4): -0.0504123634,
            (1, 4): 0.0832827019,
            (2, 4): 0.275418278,
            (3, 4): 0.0,
            (4, 4): -0.19268305,
            (0, 5): 0.00609859258,
            (1, 5): -0.00719201245,
            (2, 5): -0.0205938816,
            (3, 5): 0.0,
            (4, 5): 0.012913842},
            doc="1st order thermal conductivity paramters")
        # Viscosity paramters
        # "Release on the IAPWS Formulation 2008 for the Viscosity of
        # Ordinary Water Substance "
        self.visc_H0 = Param(RangeSet(0, 4), initialize={
            0: 1.67752,
            1: 2.20462,
            2: 0.6366564,
            3: -0.241605},
            doc="0th order viscosity parameters")

        self.visc_H1 = Param(RangeSet(0, 6), RangeSet(0, 7), initialize={
            (0, 0): 5.20094e-1,
            (1, 0): 8.50895e-2,
            (2, 0): -1.08374,
            (3, 0): -2.89555e-1,
            (4, 0): 0.0,
            (5, 0): 0.0,
            (0, 1): 2.22531e-1,
            (1, 1): 9.99115e-1,
            (2, 1): 1.88797,
            (3, 1): 1.26613,
            (4, 1): 0.0,
            (5, 1): 1.20573e-1,
            (0, 2): -2.81378e-1,
            (1, 2): -9.06851e-1,
            (2, 2): -7.72479e-1,
            (3, 2): -4.89837e-1,
            (4, 2): -2.57040e-1,
            (5, 2):  0.0,
            (0, 3): 1.61913e-1,
            (1, 3): 2.57399e-1,
            (2, 3): 0.0,
            (3, 3): 0.0,
            (4, 3): 0.0,
            (5, 3): 0.0,
            (0, 4): -3.25372e-2,
            (1, 4): 0.0,
            (2, 4): 0.0,
            (3, 4): 6.98452e-2,
            (4, 4): 0.0,
            (5, 4): 0.0,
            (0, 5): 0.0,
            (1, 5): 0.0,
            (2, 5): 0.0,
            (3, 5): 0.0,
            (4, 5): 8.72102e-3,
            (5, 5): 0.0,
            (0, 6): 0.0,
            (1, 6): 0.0,
            (2, 6): 0.0,
            (3, 6): -4.35673e-3,
            (4, 6): 0.0,
            (5, 6): -5.93264e-4},
            doc="1st order viscosity parameters")

        self.smoothing_pressure_over = Param(
                mutable=True,
                initialize=1e-4,
                doc='Smooth max parameter (pressure over)')
        self.smoothing_pressure_under = Param(
                mutable=True,
                initialize=1e-4,
                doc='Smooth max parameter (pressure under)')

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties({
            'temperature_crit': {'method': None, 'units': 'K'},
            'pressure_crit': {'method': None, 'units': 'Pa'},
            'dens_mass_crit': {'method': None, 'units': 'kg/m^3'},
            'gas_const': {'method': None, 'units': 'J/mol.K'},
            'mw': {'method': None, 'units': 'kg/mol'},
            'temperature_sat': {'method': 'None', 'units': 'K'},
            'flow_mol': {'method': None, 'units': 'mol/s'},
            'flow_mass': {'method': None, 'units': 'kg/s'},
            'temperature': {'method': None, 'units': 'K'},
            'pressure': {'method': None, 'units': 'Pa'},
            'vapor_frac': {'method': None, 'units': None},
            'dens_mass_phase': {'method': None, 'units': 'kg/m^3'},
            'temperature_red': {'method': None, 'units': None},
            'pressure_sat': {'method': None, 'units': 'kPa'},
            'energy_internal_mol_phase':  {'method': None, 'units': 'J/mol'},
            'enth_mol_phase': {'method': None, 'units': 'J/mol'},
            'entr_mol_phase': {'method': None, 'units': 'J/mol.K'},
            'cp_mol_phase': {'method': None, 'units': 'J/mol.K'},
            'cv_mol_phase': {'method': None, 'units': 'J/mol.K'},
            'speed_sound_phase': {'method': None, 'units': 'm/s'},
            'dens_mol_phase': {'method': None, 'units': 'mol/m^3'},
            'therm_cond_phase': {'method': None, 'units': 'W/m.K'},
            'visc_d_phase': {'method': None, 'units': 'Pa.s'},
            'visc_k_phase': {'method': None, 'units': 'm^2/s'},
            'phase_frac': {'method': None, 'units': None},
            'flow_mol_comp': {'method': None, 'units': 'mol/s'},
            'energy_internal_mol':  {'method': None, 'units': 'J/mol'},
            'enth_mol': {'method': None, 'units': 'J/mol'},
            'entr_mol': {'method': None, 'units': 'J/mol.K'},
            'cp_mol': {'method': None, 'units': 'J/mol.K'},
            'cv_mol': {'method': None, 'units': 'J/mol.K'},
            'heat_capacity_ratio': {'method': None, 'units': None},
            'dens_mass': {'method': None, 'units': 'kg/m^3'},
            'dens_mol': {'method': None, 'units': 'mol/m^3'},
            'dh_vap_mol': {'method': None, 'units': 'J/mol'}})

        obj.add_default_units({
            'time': 's',
            'length': 'm',
            'mass': 'kg',
            'amount': 'mol',
            'temperature': 'K',
            'energy': 'J',
            'holdup': 'mol'})


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
                    flags[i] = (v.flow_mol.fixed,
                                v.temperature.fixed,
                                v.pressure.fixed,
                                v.vapor_frac.fixed)

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


@declare_process_block_class("Iapws95StateBlock", block_class=_StateBlock,
                             doc="""This is some placeholder doc.""")
class Iapws95StateBlockData(StateBlockData):
    """
    This is a property package for calculating thermophysical properties of
    water
    """
    def initialize(self, *args, **kwargs):
        # With this particualr property pacakage there is not need for
        # initialization
        pass

    def _state_vars(self):
        """ Create the state variables
        """
        self.flow_mol = Var(initialize=1, domain=NonNegativeReals,
                            doc="Total flow [mol/s]")
        self.flow_mol.latex_symbol = "F"

        if self.state_vars == StateVars.PH:
            self.pressure = Var(domain=PositiveReals,
                                initialize=1e5,
                                doc="Pressure [Pa]",
                                bounds=(1, 1e9))
            self.pressure.latex_symbol = "P"

            self.enth_mol = Var(initialize=1000,
                                doc="Total molar enthalpy (J/mol)",
                                bounds=(1, 1e5))
            self.enth_mol.latex_symbol = "h"

            # For variables that show up in ports specify extensive/intensive
            self.extensive_set = ComponentSet((self.flow_mol,))
            self.intensive_set = ComponentSet((self.enth_mol, self.pressure))

        elif self.state_vars == StateVars.TPX:
            self.temperature = Var(initialize=300,
                                   doc="Temperature [K]",
                                   bounds=(200, 3e3))
            self.temperature.latex_symbol = "T"

            self.pressure = Var(domain=PositiveReals,
                                initialize=1e5,
                                doc="Pressure [Pa]",
                                bounds=(1, 1e9))
            self.pressure.latex_symbol = "P"

            self.vapor_frac = Var(initialize=0.0,
                                  doc="Vapor fraction [none]")
            self.vapor_frac.latex_symbol = "x"

            # For variables that show up in ports specify extensive/intensive
            self.extensive_set = ComponentSet((self.flow_mol,))
            self.intensive_set = ComponentSet((self.temperature,
                                               self.pressure,
                                               self.vapor_frac))

    def _tpx_phase_eq(self):
        # Saturation pressure
        eps_pu = self.config.parameters.smoothing_pressure_under
        eps_po = self.config.parameters.smoothing_pressure_over
        priv_plist = self.config.parameters.private_phase_list
        plist = self.config.parameters.phase_list
        rhoc = self.config.parameters.dens_mass_crit

        P = self.pressure/1000  # expression for pressure in kPa
        Psat = self.pressure_sat/1000.0  # expression for Psat in kPA
        vf = self.vapor_frac
        tau = self.tau

        # Terms for determining if you are above, below, or at the Psat
        self.P_under_sat = Expression(
                expr=smooth_max(0, Psat - P, eps_pu),
                doc="pressure above Psat, 0 if liqid exists [kPa]")
        self.P_over_sat = Expression(
                expr=smooth_max(0, P - Psat, eps_po),
                doc="pressure below Psat, 0 if vapor exists [kPa]")

        # Calculate liquid and vapor density.  If the phase doesn't exist,
        # density will be calculated at the saturation or critical pressure
        def rule_dens_mass(b, i):
            if i == "Liq":
                return rhoc*self.func_delta_liq(P + self.P_under_sat, tau)
            else:
                return rhoc*self.func_delta_vap(P - self.P_over_sat, tau)
        self.dens_mass_phase = Expression(priv_plist, rule=rule_dens_mass)

        # Reduced Density (no _mass_ identifier because mass or mol is same)
        def rule_dens_red(b, p):
            return self.dens_mass_phase[p]/rhoc
        self.dens_phase_red = Expression(priv_plist,
                                         rule=rule_dens_red,
                                         doc="reduced density [unitless]")

        # If there is only one phase fix the vapor fraction appropriatly
        if len(plist) == 1:
            if "Vap" in plist:
                self.vapor_frac.fix(1.0)
            else:
                self.vapor_frac.fix(0.0)
        elif not self.config.defined_state:
            self.eq_complementarity = Constraint(
                expr=0 == (vf*self.P_over_sat - (1 - vf)*self.P_under_sat))

        # eq_sat can activated to force the pressure to be the saturation
        # pressure, if you use this constraint deactivate eq_complementarity
        self.eq_sat = Constraint(expr=P/1000.0 == Psat/1000.0)
        self.eq_sat.deactivate()

    def build(self, *args):
        """
        Callable method for Block construction
        """
        super(Iapws95StateBlockData, self).build(*args)

        self.state_vars = self.config.parameters.state_vars
        # parameter aliases for convienient use later
        component_list = self.config.parameters.component_list
        phlist = self.config.parameters.private_phase_list
        Tc = self.config.parameters.temperature_crit
        rhoc = self.config.parameters.dens_mass_crit
        phase_set = self.config.parameters.config.phase_presentation

        self.phase_equilibrium_list = {1: ["H2O", ("Vap", "Liq")]}

        self.config.parameters.config.phase_presentation == PhaseType.MIX

        # Set if the IAPWS library is available.
        self.available = self.config.parameters.available
        if not self.available:
            _log.error("IAPWS library file not found. Was it compiled?")

        self._state_vars()  # create the appropriate state variables

        # External Functions (some of these are included only for testing)
        plib = self.config.parameters.plib
        self.func_p = EF(library=plib, function="p")
        self.func_u = EF(library=plib, function="u")
        self.func_s = EF(library=plib, function="s")
        self.func_h = EF(library=plib, function="h")
        self.func_hvpt = EF(library=plib, function="hvpt")
        self.func_hlpt = EF(library=plib, function="hlpt")
        self.func_tau = EF(library=plib, function="tau")
        self.func_vf = EF(library=plib, function="vf")
        self.func_g = EF(library=plib, function="g")
        self.func_f = EF(library=plib, function="f")
        self.func_cv = EF(library=plib, function="cv")
        self.func_cp = EF(library=plib, function="cp")
        self.func_w = EF(library=plib, function="w")
        self.func_delta_liq = EF(library=plib, function="delta_liq")
        self.func_delta_vap = EF(library=plib, function="delta_vap")
        self.func_delta_sat_l = EF(library=plib, function="delta_sat_l")
        self.func_delta_sat_v = EF(library=plib, function="delta_sat_v")
        self.func_p_sat = EF(library=plib, function="p_sat")
        self.func_tau_sat = EF(library=plib, function="tau_sat")
        self.func_phi0 = EF(library=plib, function="phi0")
        self.func_phi0_delta = EF(library=plib, function="phi0_delta")
        self.func_phi0_delta2 = EF(library=plib, function="phi0_delta2")
        self.func_phi0_tau = EF(library=plib, function="phi0_tau")
        self.func_phi0_tau2 = EF(library=plib, function="phi0_tau2")
        self.func_phir = EF(library=plib, function="phir")
        self.func_phir_delta = EF(library=plib, function="phir_delta")
        self.func_phir_delta2 = EF(library=plib, function="phir_delta2")
        self.func_phir_tau = EF(library=plib, function="phir_tau")
        self.func_phir_tau2 = EF(library=plib, function="phir_tau2")
        self.func_phir_delta_tau = EF(library=plib, function="phir_delta_tau")

        # Calculations

        # molecular weight
        self.mw = Expression(expr=self.config.parameters.mw,
                             doc="molecular weight [kg/mol]")
        mw = self.mw
        mw.latex_symbol = "M"
        P = self.pressure/1000.0  # Pressure expr [kPA] (for external func)

        if self.state_vars == StateVars.PH:
            # if the state vars are P and H create expressions for T and x
            h_mass = self.enth_mol/mw/1000  # enthalpy expr [kJ/kg]
            self.temperature = Expression(
                expr=Tc/self.func_tau(h_mass, P),
                doc="Temperature (K)")
            self.temperature.latex_symbol = "T"
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
                    doc="Vapor mole fraction (mol vapor/mol total)")
            self.vapor_frac.latex_symbol = "x"
        elif self.state_vars == StateVars.TPX:
            # Need to get enthalpy expressions in here
            pass

        # Convenient shorter names and expressions
        T = self.temperature
        vf = self.vapor_frac

        # Saturation temperature expression
        self.temperature_sat = Expression(expr=Tc/self.func_tau_sat(P),
                                          doc="Stauration temperature (K)")
        self.temperature_sat.latex_symbol = "T_{sat}"

        # Saturation tau (tau = Tc/T)
        self.tau_sat = Expression(expr=self.func_tau_sat(P))

        # Reduced temperature
        self.temperature_red = Expression(
                expr=T/Tc,
                doc="reduced temperature T/Tc (unitless)")
        self.temperature_red.latex_symbol = "T_r"

        self.tau = Expression(expr=Tc/T, doc="Tc/T (unitless)")
        tau = self.tau
        self.tau.latex_symbol = "\\tau"

        # Saturation pressure
        self.pressure_sat = Expression(expr=1000*self.func_p_sat(tau),
                                       doc="Saturation pressure (Pa)")
        self.pressure_sat.latex_symbol = "P_{sat}"
        self.pressure_sat/1000.0  # expression for Psat in kPA

        if self.state_vars == StateVars.PH:
            # If TPx state vars the expressions are given in _tpx_phase_eq
            # Calculate liquid and vapor density.  If the phase doesn't exist,
            # density will be calculated at the saturation or critical pressure
            # depending on whether the temperature is above the critical
            # temperature supercritical fluid is considered to be the liquid
            # phase
            def rule_dens_mass(b, i):
                if i == "Liq":
                    return rhoc*self.func_delta_liq(P, tau)
                else:
                    return rhoc*self.func_delta_vap(P, tau)
            self.dens_mass_phase = Expression(
                    phlist,
                    rule=rule_dens_mass,
                    doc="Mass density by phase (kg/m3)")
            self.dens_mass_phase.latex_symbol = "\\rho"

            # Reduced Density (no _mass_ identifier as mass or mol is same)
            def rule_dens_red(b, p):
                return self.dens_mass_phase[p]/rhoc
            self.dens_phase_red = Expression(phlist,
                                             rule=rule_dens_red,
                                             doc="reduced density (unitless)")
        elif self.state_vars == StateVars.TPX:
            self._tpx_phase_eq()
        delta = self.dens_phase_red  # this shorter name is from IAPWS
        self.dens_phase_red.latex_symbol = "\\delta"

        # Phase property expressions all converted to SI

        # Saturated Enthalpy
        def rule_enth_mol_sat_phase(b, p):
            if p == "Liq":
                return 1000*mw*self.func_hlpt(P, self.tau_sat)
            elif p == "Vap":
                return 1000*mw*self.func_hvpt(P, self.tau_sat)
        self.enth_mol_sat_phase = Expression(
                phlist,
                rule=rule_enth_mol_sat_phase,
                doc="Saturated enthalpy of the phases at pressure (J/mol)")

        self.dh_vap_mol = Expression(
            expr=self.enth_mol_sat_phase["Vap"] -
            self.enth_mol_sat_phase["Liq"],
            doc="Enthaply of vaporization at pressure and saturation (J/mol)")

        # Phase Internal Energy
        def rule_energy_internal_mol_phase(b, p):
            return 1000*mw*self.func_u(delta[p], tau)
        self.energy_internal_mol_phase = Expression(
                phlist,
                rule=rule_energy_internal_mol_phase,
                doc="Phase internal energy or saturated if phase "
                "doesn't exist [J/mol]")

        # Phase Enthalpy
        def rule_enth_mol_phase(b, p):
            return 1000*mw*self.func_h(delta[p], tau)
        self.enth_mol_phase = Expression(
            phlist,
            rule=rule_enth_mol_phase,
            doc="Phase enthalpy or saturated if phase doesn't exist [J/mol]")

        # Phase Entropy
        def rule_entr_mol_phase(b, p):
            return 1000*mw*self.func_s(delta[p], tau)
        self.entr_mol_phase = Expression(
            phlist,
            rule=rule_entr_mol_phase,
            doc="Phase entropy or saturated if phase doesn't exist [J/mol/K]")

        # Phase constant pressure heat capacity, cp
        def rule_cp_mol_phase(b, p):
            return 1000*mw*self.func_cp(delta[p], tau)
        self.cp_mol_phase = Expression(
            phlist,
            rule=rule_cp_mol_phase,
            doc="Phase cp or saturated if phase doesn't exist [J/mol/K]")

        # Phase constant pressure heat capacity, cv
        def rule_cv_mol_phase(b, p):
            return 1000*mw*self.func_cv(delta[p], tau)
        self.cv_mol_phase = Expression(
            phlist,
            rule=rule_cv_mol_phase,
            doc="Phase cv or saturated if phase doesn't exist [J/mol/K]")

        # Phase speed of sound
        def rule_speed_sound_phase(b, p):
            return self.func_w(delta[p], tau)
        self.speed_sound_phase = Expression(
            phlist,
            rule=rule_speed_sound_phase,
            doc="Phase speed of sound or saturated if phase "
            "doesn't exist [m/s]")

        # Phase Mole density
        def rule_dens_mol_phase(b, p):
            return self.dens_mass_phase[p]/mw
        self.dens_mol_phase = Expression(
            phlist,
            rule=rule_dens_mol_phase,
            doc="Phase mole density or saturated if phase "
            "doesn't exist [mol/m3]")

        # Phase Thermal conductiviy
        def rule_tc(b, p):
            L0 = self.config.parameters.tc_L0
            L1 = self.config.parameters.tc_L1
            return (1e-3*sqrt(1.0/tau)/sum(L0[i]*tau**i for i in L0) *
                    exp(delta[p]*sum((tau - 1)**i *
                        sum(L1[i, j]*(delta[p] - 1)**j
                            for j in range(0, 6)) for i in range(0, 5))))
        self.therm_cond_phase = Expression(
                phlist,
                rule=rule_tc,
                doc="Thermal conductivity [W/K/m]")

        # Phase dynamic viscosity
        def rule_mu(b, p):
            H0 = self.config.parameters.visc_H0
            H1 = self.config.parameters.visc_H1
            return (1e-4*sqrt(1.0/tau)/sum(H0[i]*tau**i for i in H0) *
                    exp(delta[p]*sum((tau - 1)**i *
                        sum(H1[i, j]*(delta[p] - 1)**j
                            for j in range(0, 7)) for i in range(0, 6))))
        self.visc_d_phase = Expression(
                phlist,
                rule=rule_mu,
                doc="Viscosity (dynamic) [Pa*s]")

        # Phase kinimatic viscosity
        def rule_nu(b, p):
            return self.visc_d_phase[p]/self.dens_mass_phase[p]
        self.visc_k_phase = Expression(phlist,
                                       rule=rule_nu,
                                       doc="Kinematic viscosity [m^2/s]")

        # Phase fraction
        def rule_phase_frac(b, p):
            if p == "Vap":
                return vf
            elif p == "Liq":
                return 1.0 - vf
        self.phase_frac = Expression(phlist,
                                     rule=rule_phase_frac,
                                     doc="Phase fraction [unitless]")

        # Component flow (for units that need it)
        def component_flow(b, i):
            return self.flow_mol
        self.flow_mol_comp = Expression(
                component_list,
                rule=component_flow,
                doc="Total flow (both phases) of component [mol/s]")

        # Total (mixed phase) properties

        # Enthalpy
        if self.state_vars == StateVars.TPX:
            self.enth_mol = Expression(
                    expr=sum(self.phase_frac[p]*self.enth_mol_phase[p]
                             for p in phlist))
            self.enth_mol.latex_symbol = "h"
        # Internal Energy
        self.energy_internal_mol = Expression(
                expr=sum(self.phase_frac[p]*self.energy_internal_mol_phase[p]
                         for p in phlist))
        self.energy_internal_mol.latex_symbol = "u"
        # Entropy
        self.entr_mol = Expression(
                expr=sum(self.phase_frac[p]*self.entr_mol_phase[p]
                         for p in phlist))
        self.entr_mol.latex_symbol = "s"
        # cp
        self.cp_mol = Expression(
                expr=sum(self.phase_frac[p]*self.cp_mol_phase[p]
                         for p in phlist))
        self.cp_mol.latex_symbol = "c_p"
        # cv
        self.cv_mol = Expression(
                expr=sum(self.phase_frac[p]*self.cv_mol_phase[p]
                         for p in phlist))
        self.cv_mol.latex_symbol = "c_v"
        # mass density
        self.dens_mass = Expression(
                expr=1.0/sum(self.phase_frac[p]*1.0/self.dens_mass_phase[p]
                             for p in phlist))
        # mole density
        self.dens_mol = Expression(
                expr=1.0/sum(self.phase_frac[p]*1.0/self.dens_mol_phase[p]
                             for p in phlist))
        # heat capacity ratio
        self.heat_capacity_ratio = Expression(expr=self.cp_mol/self.cv_mol)
        # Flows
        self.flow_vol = Expression(
                expr=self.flow_mol/self.dens_mol,
                doc="Total liquid + vapor volumetric flow (m3/s)")

        self.flow_mass = Expression(expr=self.mw*self.flow_mol,
                                    doc="mass flow rate [kg/s]")

        self.enth_mass = Expression(expr=self.enth_mol/mw,
                                    doc="Mass enthalpy (J/kg)")

        # Set the state vars dictionary
        if self.state_vars == StateVars.PH:
            self._state_vars_dict = {
                "flow_mol": self.flow_mol,
                "enth_mol": self.enth_mol,
                "pressure": self.pressure}
        elif self.state_vars == StateVars.TPX and \
                phase_set in (PhaseType.MIX, PhaseType.LG):
            self._state_vars_dict = {
                "flow_mol": self.flow_mol,
                "temperature": self.temperature,
                "pressure": self.pressure,
                "vapor_frac": self.vapor_frac}
        elif self.state_vars == StateVars.TPX and \
                phase_set in (PhaseType.G, PhaseType.L):
            self._state_vars_dict = {
                "flow_mol": self.flow_mol,
                "temperature": self.temperature,
                "pressure": self.pressure}

    def get_material_flow_terms(self, p, j):
        if p == "Mix":
            return self.flow_mol
        else:
            return self.flow_mol*self.phase_frac[p]

    def get_enthalpy_flow_terms(self, p):
        if p == "Mix":
            return self.enth_mol*self.flow_mol
        else:
            return self.enth_mol_phase[p]*self.phase_frac[p]*self.flow_mol

    def get_material_density_terms(self, p, j):
        if p == "Mix":
            return self.dens_mol
        else:
            return self.dens_mol_phase[p]

    def get_energy_density_terms(self, p):
        if p == "Mix":
            return self.dens_mol*self.energy_internal_mol
        else:
            return self.dens_mol_phase[p]*self.energy_internal_mol_phase[p]

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
            "Molar Enthalpy (J/mol)": self.enth_mol_phase}

    def extensive_state_vars(self):
        return self.extensive_set

    def intensive_state_vars(self):
        return self.intensive_set

    def model_check(self):
        pass
