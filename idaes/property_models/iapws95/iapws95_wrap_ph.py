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

Dropped all critical enhancments and non-analytic terms ment to imporve accruacy
near the critical point. These tend to cause singularities in the equations, and
it is assumend that we will try to avoid operating very close to the critical
point.

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
   Akasaka, R. (2008). "A Reliable and Useful Method to Determine the Saturation
       State from Helmholtz Energy Equations of State." Journal of Thermal
       Science and Technology, 3(3), 442-451.
   International Association for the Properties of Water and Steam (2011).
       IAPWS R15-11, "Release on the IAPWS Formulation 2011 for the
       Thermal Conductivity of Ordinary Water Substance,"
       URL: http://iapws.org/relguide/ThCond.pdf
   International Association for the Properties of Water and Steam (2008).
       IAPWS R12-08, "Release on the IAPWS Formulation 2008 for the Viscosity of
       Ordinary Water Substance,"
       URL: http://iapws.org/relguide/visc.pdf
"""
from __future__ import division

__author__ = "John Eslick"

# Import Python libraries
import logging
import os

# Import Pyomo libraries
from pyomo.environ import Constraint, Expression, Param, PositiveReals,\
                          RangeSet, Reals, Set, value, Var, ExternalFunction,\
                          NonNegativeReals, exp, sqrt, log, tanh, ConcreteModel
from pyomo.opt import SolverFactory, TerminationCondition
from pyutilib.misc.config import ConfigValue

# Import IDAES
from idaes.core import declare_process_block_class, ProcessBlock, \
                       StateBlock, StateBlockDataBase, PhysicalParameterBlock

# Logger
_log = logging.getLogger(__name__)

def htpx(T, P=None, x=None):
    """
    Conveneince function to calculate steam enthalpy from temperature and
    either pressure or vapor fraction. This function can be used for inlet
    streams and initialization where temperature is known instread of enthalpy.
    Args:
        T: Temperature [K]
        P: Pressure [Pa], None if saturated steam
        x: Vapor fraction [mol vapor/mol total], None if superheated or subcooled

    Returns:
        Molar enthalpy [J/mol].
    """
    model = ConcreteModel()
    model.prop_param = Iapws95ParameterBlock()
    prop = model.prop = Iapws95StateBlock(default={"parameters":model.prop_param})

    if x is None:
        Tsat = 647.096/value(prop.func_tau_sat(P/1000))
        if value(T) < Tsat or value(P) > 22064: #liquid
            return value(prop.func_hlpt(P/1000, 647.096/T)*prop.mw*1000.0)
        else: #vapor
            return value(prop.func_hvpt(P/1000, 647.096/T)*prop.mw*1000.0)
    if P is None:
        Psat = value(prop.func_p_sat(647.096/T))
        return value(prop.func_hlpt(Psat, 647.096/T)*prop.mw*1000.0)*(x-1) +\
            value(prop.func_hlpt(Psat, 647.096/T)*prop.mw*1000.0)*x

@declare_process_block_class("Iapws95ParameterBlock")
class Iapws95ParameterBlockData(PhysicalParameterBlock):

    def build(self):
        super(Iapws95ParameterBlockData, self).build()
        self.state_block_class = Iapws95StateBlock
        # Location of the *.so or *.dll file for external functions
        self.plib = os.path.dirname(__file__)
        self.plib = os.path.join(self.plib, "iapws95.so")
        # Phase list
        self.private_phase_list = Set(initialize=["Vap", "Liq"])
        self.phase_list = Set(initialize=["Mix"])
        # Component list - a list of component identifiers
        self.component_list = Set(initialize=['H2O'])
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
        #Thermal conductivity parameters.
        # "Release on the IAPWS Formulation 2011 for the Thermal Conductivity of
        # Ordinary Water Substance"
        self.tc_L0 = Param(RangeSet(0,5), initialize={
            0:2.443221e-3,
            1:1.323095e-2,
            2:6.770357e-3,
            3:-3.454586e-3,
            4:4.096266e-4},
            doc="0th order themalcondutivity paramters")

        self.tc_L1 = Param(RangeSet(0,5), RangeSet(0,6), initialize={
            (0,0):1.60397357,
            (1,0):2.33771842,
            (2,0):2.19650529,
            (3,0):-1.21051378,
            (4,0):-2.7203370,
            (0,1):-0.646013523,
            (1,1):-2.78843778,
            (2,1):-4.54580785,
            (3,1):1.60812989,
            (4,1):4.57586331,
            (0,2):0.111443906,
            (1,2):1.53616167,
            (2,2):3.55777244,
            (3,2):-0.621178141,
            (4,2):-3.18369245,
            (0,3):0.102997357,
            (1,3):-0.463045512,
            (2,3):-1.40944978,
            (3,3):0.0716373224,
            (4,3):1.1168348,
            (0,4):-0.0504123634,
            (1,4):0.0832827019,
            (2,4):0.275418278,
            (3,4):0.0,
            (4,4):-0.19268305,
            (0,5):0.00609859258,
            (1,5):-0.00719201245,
            (2,5):-0.0205938816,
            (3,5):0.0,
            (4,5):0.012913842},
            doc="1st order themalcondutivity paramters")
        #Viscosity paramters
        #"Release on the IAPWS Formulation 2008 for the Viscosity of
        # Ordinary Water Substance "
        self.visc_H0 = Param(RangeSet(0,4), initialize={
            0:1.67752,
            1:2.20462,
            2:0.6366564,
            3:-0.241605},
            doc="0th order viscosity parameters")

        self.visc_H1 = Param(RangeSet(0,6), RangeSet(0,7), initialize={
            (0,0):5.20094e-1,
            (1,0):8.50895e-2,
            (2,0):-1.08374,
            (3,0):-2.89555e-1,
            (4,0):0.0,
            (5,0):0.0,
            (0,1):2.22531e-1,
            (1,1):9.99115e-1,
            (2,1):1.88797,
            (3,1):1.26613,
            (4,1):0.0,
            (5,1):1.20573e-1,
            (0,2):-2.81378e-1,
            (1,2):-9.06851e-1,
            (2,2):-7.72479e-1,
            (3,2):-4.89837e-1,
            (4,2):-2.57040e-1,
            (5,2):0.0,
            (0,3):1.61913e-1,
            (1,3):2.57399e-1,
            (2,3):0.0,
            (3,3):0.0,
            (4,3):0.0,
            (5,3):0.0,
            (0,4):-3.25372e-2,
            (1,4):0.0,
            (2,4):0.0,
            (3,4):6.98452e-2,
            (4,4):0.0,
            (5,4):0.0,
            (0,5):0.0,
            (1,5):0.0,
            (2,5):0.0,
            (3,5):0.0,
            (4,5):8.72102e-3,
            (5,5):0.0,
            (0,6):0.0,
            (1,6):0.0,
            (2,6):0.0,
            (3,6):-4.35673e-3,
            (4,6):0.0,
            (5,6):-5.93264e-4},
            doc="1st order viscosity parameters")

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties({
            'temperature_crit': {'method': None, 'units': 'K'},
            'pressure_crit': {'method': None, 'units': 'Pa'},
            'dens_mass_crit': {'method': None, 'units': 'kg/m^3'},
            'gas_const': {'method': None, 'units': 'J/mol.K'},
            'mw': {'method': None, 'units': 'kg/mol'},
            'temperature_sat': {'method': '_temperature_sat', 'units': 'K'},
            'flow_mol': {'method': None, 'units': 'mol/s'},
            'temperature': {'method': None, 'units': 'K'},
            'pressure': {'method': None, 'units': 'Pa'},
            'vapor_frac': {'method': None, 'units': None},
            'dens_mass_phase': {'method': None, 'units': 'kg/m^3'},
            'temperature_red': {'method': None, 'units': None},
            'pressure_sat': {'method': None, 'units': 'kPa'},
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
            'enth_mol': {'method': None, 'units': 'J/mol'},
            'entr_mol': {'method': None, 'units': 'J/mol.K'},
            'cp_mol': {'method': None, 'units': 'J/mol.K'},
            'cv_mol': {'method': None, 'units': 'J/mol.K'},
            'heat_capacity_ratio': {'method': None, 'units': None},
            'dens_mass': {'method': None, 'units': 'kg/m^3'},
            'dens_mol': {'method': None, 'units': 'mol/m^3'}})

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
    def initialize(self, *args, **kwargs):
        for i in self.keys():
            self[i].initialize(*args, **kwargs)

    def release_state(*args, **kwargs):
        pass


@declare_process_block_class("Iapws95StateBlock", block_class=_StateBlock, doc="""
    This is some placeholder doc.
    """)
class Iapws95StateBlockData(StateBlockDataBase):
    """
    This is a property package for calcuating thermophysical properties of water
    """

    def initialize(self, *args, **kwargs):
        Fold = self.flow_mol.value
        Hold = self.enth_mol.value
        Pold = self.pressure.value
        self.flow_mol.value = kwargs.get("flow_mol", Fold)
        self.enth_mol.value = kwargs.get("enth_mol", Hold)
        self.pressure.value = kwargs.get("presssure", Pold)
        #
        # This is where stuff would get initialized if there was anything to do
        #
        # reset the values of fixed variables so this doesn't propblem
        if self.flow_mol.is_fixed():
            self.flow_mol.value = Fold
        if self.pressure.is_fixed():
            self.pressure.value = Pold
        if self.enth_mol.is_fixed():
            self.enth_mol.value = Hold

    def build(self, *args):
        """
        Callable method for Block construction
        """
        super(Iapws95StateBlockData, self).build(*args)

        # parameters
        component_list = self.config.parameters.component_list
        phase_list = self.config.parameters.phase_list
        private_phase_list = self.config.parameters.private_phase_list
        temperature_crit = self.config.parameters.temperature_crit
        pressure_crit = self.config.parameters.pressure_crit
        dens_mass_crit = self.config.parameters.dens_mass_crit
        gas_const = self.config.parameters.gas_const

        # Variables
        self.flow_mol = Var(initialize=1, domain=NonNegativeReals,
            doc="Total flow [mol/s]")
        self.flow_mol.latex_symbol = "F"

        self.pressure = Var(domain=PositiveReals, initialize=1e5,
            doc="Pressure [Pa]")
        self.pressure.latex_symbol = "P"

        self.enth_mol = Var(initialize=1000,
            doc="Total molar enthalpy (J/mol)")
        self.enth_mol.latex_symbol = "h"

        # External Functions (some of these are included only for testing)
        plib = self.config.parameters.plib
        self.func_p = ExternalFunction(library=plib, function="p")
        self.func_u = ExternalFunction(library=plib, function="u")
        self.func_s = ExternalFunction(library=plib, function="s")
        self.func_h = ExternalFunction(library=plib, function="h")
        self.func_hvpt = ExternalFunction(library=plib, function="hvpt")
        self.func_hlpt = ExternalFunction(library=plib, function="hlpt")
        self.func_tau = ExternalFunction(library=plib, function="tau")
        self.func_vf = ExternalFunction(library=plib, function="vf")
        self.func_g = ExternalFunction(library=plib, function="g")
        self.func_f = ExternalFunction(library=plib, function="f")
        self.func_cv = ExternalFunction(library=plib, function="cv")
        self.func_cp = ExternalFunction(library=plib, function="cp")
        self.func_w = ExternalFunction(library=plib, function="w")
        self.func_delta_liq = ExternalFunction(library=plib,
            function="delta_liq")
        self.func_delta_vap = ExternalFunction(library=plib,
            function="delta_vap")
        self.func_delta_sat_l = ExternalFunction(library=plib,
            function="delta_sat_l")
        self.func_delta_sat_v = ExternalFunction(library=plib,
            function="delta_sat_v")
        self.func_p_sat = ExternalFunction(library=plib,
            function="p_sat")
        self.func_tau_sat = ExternalFunction(library=plib,
            function="tau_sat")
        self.func_phi0 = ExternalFunction(library=plib,
            function="phi0")
        self.func_phi0_delta = ExternalFunction(library=plib,
            function="phi0_delta")
        self.func_phi0_delta2 = ExternalFunction(library=plib,
            function="phi0_delta2")
        self.func_phi0_tau = ExternalFunction(library=plib,
            function="phi0_tau")
        self.func_phi0_tau2 = ExternalFunction(library=plib,
            function="phi0_tau2")
        self.func_phir = ExternalFunction(library=plib,
            function="phir")
        self.func_phir_delta = ExternalFunction(library=plib,
            function="phir_delta")
        self.func_phir_delta2 = ExternalFunction(library=plib,
            function="phir_delta2")
        self.func_phir_tau = ExternalFunction(library=plib,
            function="phir_tau")
        self.func_phir_tau2 = ExternalFunction(library=plib,
            function="phir_tau2")
        self.func_phir_delta_tau = ExternalFunction(library=plib,
            function="phir_delta_tau")

        # Calcuations (In this case all expressions no constraints)

        # molecular weight
        self.mw = Expression(expr=self.config.parameters.mw,
            doc="molecular weight [kg/mol]")
        mw = self.mw
        mw.latex_symbol = "M"

        self.enth_mass = Expression(expr = self.enth_mol/mw,
            doc="Mass enthalpy (J/kg)")

        # Conveneint shorter names and expressions
        P = self.pressure/1000.0 # Pressure expr [kPA] (for external func)
        h_mass = self.enth_mass/1000 #enthalpy expr [kJ/kg] (for external func)
        Tc = temperature_crit # critical temperature

        # Temperature expression, external function takes
        self.temperature = Expression(expr=Tc/self.func_tau(h_mass, P),
            doc="Temperature (K)")
        self.temperature.latex_symbol = "T"

        self.vapor_frac = Expression(expr=self.func_vf(h_mass, P),
            doc="Vapor mole fraction (mol vapor/mol total)")
        self.vapor_frac.latex_symbol = "y"

        # More conveneint shorter names and expressions
        rhoc = dens_mass_crit
        T = self.temperature
        vf = self.vapor_frac
        phlist = private_phase_list

        # Saturation temperature expression
        self.temperature_sat = Expression(expr=Tc/self.func_tau_sat(P),
            doc="Stauration temperature (K)")
        self.temperature_sat.latex_symbol = "T_\{sat\}"

        # Saturation tau (tau = Tc/T)
        self.tau_sat = Expression(expr=self.func_tau_sat(P))

        # Reduced temperature
        self.temperature_red = Expression(expr=T/Tc,
            doc="reduced temperature T/Tc (unitless)")
        self.temperature_red.latex_symbol = "T_r"

        self.tau = Expression(expr=Tc/T, doc="Tc/T (unitless)")
        tau = self.tau
        self.tau.latex_symbol = "\\tau"

        # Saturation pressure
        self.pressure_sat = Expression(expr=1000*self.func_p_sat(tau),
            doc="Saturation pressure (Pa)")
        self.pressure_sat.latex_symbol = "P_\{sat\}"
        Psat = self.pressure_sat/1000.0 # expression for Psat in kPA

        # Calculate liquid and vapor density.  If the phase doesn't exist,
        # density will be calculated at the saturation or critical pressure
        def rule_dens_mass(b, i):
            if i=="Liq":
                return rhoc*self.func_delta_liq(P, tau)
            else:
                return rhoc*self.func_delta_vap(P, tau)
        self.dens_mass_phase = Expression(phlist, rule=rule_dens_mass,
            doc="Mass density by phase (kg/m3)")
        self.dens_mass_phase.latex_symbol = "\\rho"

        # Reduced Density (no _mass_ identifier because mass or mol is same)
        def rule_dens_red(b, p):
            return self.dens_mass_phase[p]/rhoc
        self.dens_phase_red = Expression(phlist, rule=rule_dens_red,
            doc="reduced density (unitless)")
        self.dens_phase_red.latex_symbol = "\\delta"
        delta = self.dens_phase_red # this shorter name is from IAPWS

        # Phase property expressions all converted to SI

        # Saturated Enthalpy
        def rule_enth_mol_sat_phase(b, p):
            if p == "Liq":
                return 1000*mw*self.func_hlpt(P, self.tau_sat)
            elif p == "Vap":
                return 1000*mw*self.func_hvpt(P, self.tau_sat)
        self.enth_mol_sat_phase = Expression(phlist,
            rule=rule_enth_mol_sat_phase,
            doc="Saturated enthalpy of the phases at pressure (J/mol)")

        # Phase Enthalpy
        def rule_enth_mol_phase(b, p):
            return 1000*mw*self.func_h(delta[p], tau)
        self.enth_mol_phase = Expression(phlist,
            rule=rule_enth_mol_phase,
            doc="Phase enthalpy or saturated if phase doesn't exist [J/mol]")

        # Phase Entropy
        def rule_entr_mol_phase(b, p):
            return 1000*mw*self.func_s(delta[p], tau)
        self.entr_mol_phase = Expression(phlist,
            rule=rule_entr_mol_phase,
            doc="Phase entropy or saturated if phase doesn't exist [J/mol/K]")

        # Phase constant pressure heat capacity, cp
        def rule_cp_mol_phase(b, p):
            return 1000*mw*self.func_cp(delta[p], tau)
        self.cp_mol_phase = Expression(phlist,
            rule=rule_cp_mol_phase,
            doc="Phase cp or saturated if phase doesn't exist [J/mol/K]")

        # Phase constant pressure heat capacity, cv
        def rule_cv_mol_phase(b, p):
            return 1000*mw*self.func_cv(delta[p], tau)
        self.cv_mol_phase = Expression(phlist,
            rule=rule_cv_mol_phase,
            doc="Phase cv or saturated if phase doesn't exist [J/mol/K]")

        # Phase speed of sound
        def rule_speed_sound_phase(b, p):
            return self.func_w(delta[p], tau)
        self.speed_sound_phase = Expression(phlist,
            rule=rule_speed_sound_phase,
            doc="Phase speed of sound or saturated if phase doesn't exist [m/s]")

        # Phase Mole density
        def rule_dens_mol_phase(b, p):
            return self.dens_mass_phase[p]/mw
        self.dens_mol_phase = Expression(phlist,
            rule=rule_dens_mol_phase,
            doc="Phase mole density or saturated if phase doesn't exist [mol/m3]")

        # Phase Thermal conductiviy
        def rule_tc(b, p):
            L0 = self.config.parameters.tc_L0
            L1 = self.config.parameters.tc_L1
            return 1e-3*sqrt(1.0/tau)/sum(L0[i]*tau**i for i in L0)*\
                exp(delta[p]*sum((tau - 1)**i*sum(L1[i,j]*(delta[p] - 1)**j\
                    for j in range(0,6)) for i in range(0,5)))
        self.therm_cond_phase = Expression(phlist, rule=rule_tc,
            doc="Thermal conductivity [W/K/m]")

        # Phase dynamic viscosity
        def rule_mu(b, p):
            H0 = self.config.parameters.visc_H0
            H1 = self.config.parameters.visc_H1
            return 1e-4*sqrt(1.0/tau)/sum(H0[i]*tau**i for i in H0)*\
                exp(delta[p]*sum((tau - 1)**i*sum(H1[i,j]*(delta[p] - 1)**j\
                    for j in range(0,7)) for i in range(0,6)))
        self.visc_d_phase = Expression(phlist, rule=rule_mu,
            doc="Viscosity (dynamic) [Pa*s]")

        # Phase kinimatic viscosity
        def rule_nu(b, p):
            return self.visc_d_phase[p]/self.dens_mass_phase[p]
        self.visc_k_phase = Expression(phlist, rule=rule_nu,
            doc="Kinematic viscosity [m^2/s]")

        #Phase fraction
        def rule_phase_frac(b, p):
            if p == "Vap":
                return vf
            elif p == "Liq":
                return 1.0 - vf
        self.phase_frac = Expression(phlist,
            rule=rule_phase_frac, doc="Phase fraction [unitless]")

        # Component flow (for units that need it)
        def component_flow(b, i):
            return self.flow_mol
        self.flow_mol_comp = Expression(component_list,
            rule=component_flow,
            doc="Total flow (both phases) of component [mol/s]")

        # Total (mixed phase) properties

        #Entropy
        self.entr_mol = Expression(expr=
            sum(self.phase_frac[p]*self.entr_mol_phase[p] for p in phlist))
        self.entr_mol.latex_symbol = "s"
        #cp
        self.cp_mol = Expression(expr=
            sum(self.phase_frac[p]*self.cp_mol_phase[p] for p in phlist))
        self.cp_mol.latex_symbol = "c_p"
        #cv
        self.cv_mol = Expression(expr=
            sum(self.phase_frac[p]*self.cv_mol_phase[p] for p in phlist))
        self.cv_mol.latex_symbol = "c_v"
        #mass density
        self.dens_mass = Expression(expr=
            1.0/sum(self.phase_frac[p]*1.0/self.dens_mass_phase[p]
                    for p in phlist))
        #mole density
        self.dens_mol = Expression(expr=
            1.0/sum(self.phase_frac[p]*1.0/self.dens_mol_phase[p]
                    for p in phlist))
        #heat capacity ratio
        self.heat_capacity_ratio = Expression(expr=self.cp_mol/self.cv_mol)
        #Flows
        self.flow_vol = Expression(expr=self.flow_mol/self.dens_mol,
            doc="Total liquid + vapor volumetric flow (m3/s)")

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

    def get_enthalpy_density_terms(self, p):
        if p == "Mix":
            return self.dens_mol*self.enth_mol
        else:
            return self.dens_mol_phase[p]*self.enth_mol_phase[p]

    def define_state_vars(self):
        return {"flow_mol": self.flow_mol,
                "enth_mol": self.enth_mol,
                "pressure": self.pressure}

    def model_check(self):
        pass
