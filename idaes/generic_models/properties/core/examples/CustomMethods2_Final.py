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
Custom Methods for calculating pure component properties for CO2
its equilibrium equation, and molar liquid phase density of solution:
"""
from pyomo.environ import exp, log, Var, units as pyunits
import idaes.generic_models.properties.core.pure.NIST as NIST
import idaes.generic_models.properties.core.pure.Perrys as Perrys
from idaes.core.util.exceptions import PropertyNotSupportedError
from idaes.generic_models.properties.core.generic.utility import (
    get_method, get_component_object as cobj)
from idaes.generic_models.properties.core.examples.CustomMethods_Final import \
    vol_mol_liq_comp

from idaes.core.util.misc import set_param_from_config

# Molar Enthalpy of CO2 dissolved in liquid phase (J/mol)
class enth_mol_liq_comp():

    @staticmethod
    def build_parameters(cobj):
        NIST.cp_mol_ig_comp.build_parameters(cobj)

    @staticmethod
    def return_expression(b, cobj, T):
        
        # Specific enthalpy via the Shomate equation
        t = pyunits.convert(T, to_units=pyunits.kiloK)
        
        h_form = b.params.config.include_enthalpy_of_formation
        H = (cobj.cp_mol_ig_comp_coeff_H if not h_form
              else 0*pyunits.kJ*pyunits.mol**-1)

        h = (cobj.cp_mol_ig_comp_coeff_A*(t) +
              (cobj.cp_mol_ig_comp_coeff_B/2)*(t**2) +
              (cobj.cp_mol_ig_comp_coeff_C/3)*(t**3) +
              (cobj.cp_mol_ig_comp_coeff_D/4)*(t**4) -
              cobj.cp_mol_ig_comp_coeff_E*(1/t) +
              cobj.cp_mol_ig_comp_coeff_F -
              H)

        units = b.params.get_metadata().derived_units
        h1 = pyunits.convert(h, units["energy_mole"])
        
        # Hliq(T) = Hig(T) + R*T*ln(H_CO2_H2O*xco2/Pref)
        # Data source [4], Eqn 13-33
        htot = h1 + 8.314*(pyunits.J/(pyunits.mol*pyunits.K))*T*\
            log(henrys_constant_co2_h2o.return_expression(b,cobj,T)*\
                b.mole_frac_phase_comp['Liq', 'CO2']/\
                (b.params.pressure_ref))
        
        return htot
    
class pressure_sat_comp():
    
    @staticmethod
    def build_parameters(cobj):
        pass
    
    @staticmethod
    def return_expression(b, cobj, T, dT=False):
        return 0

# Henrys constant of CO2 dissolved in liquid phase (Pa.m3/mol)
# Data source [3], N2O Analogy 
class henrys_constant():
    
    @staticmethod
    def build_parameters(cobj):
        pass
    
    @staticmethod
    def return_expression(b, cobj, T, xmea, xh2o, xco2):
        T = pyunits.convert(T, to_units=pyunits.K)
        t      = T - (273.15*pyunits.K)
        den = (xmea*61.0831E-3*(pyunits.kg/pyunits.mol)) + \
            (xh2o*18.0153E-3*(pyunits.kg/pyunits.mol))
        wt_MEA = xmea*61.0831E-3*(pyunits.kg/pyunits.mol)/den
        wt_H2O = xh2o*18.0153E-3*(pyunits.kg/pyunits.mol)/den
        H_N2O_MEA = 2.448e5*exp(-1348*pyunits.K/T)
        H_CO2_H2O = 3.52e6*exp(-2113*pyunits.K/T)
        H_N2O_H2O = 8.449e6*exp(-2283*pyunits.K/T)
        H_CO2_MEA = H_N2O_MEA*(H_CO2_H2O/H_N2O_H2O)
        lwm = 1.70981 + 0.03972*(1/pyunits.K)*t - \
            4.3e-4*(1/pyunits.K**2)*t**2 - 2.20377*wt_H2O
        
        return (exp(wt_MEA*log(H_CO2_MEA) +wt_H2O*log(H_CO2_H2O)
                    + wt_MEA*wt_H2O*lwm))*pyunits.Pa*pyunits.m**3/pyunits.mol

# Henrys constant of CO2 in H2O (Pa/mole fraction)
# Data source [4], table 6.2-7
class henrys_constant_co2_h2o():
    
    @staticmethod
    def build_parameters(cobj):
        pass
    
    @staticmethod
    def return_expression(b, cobj, T):
        T = pyunits.convert(T, to_units=pyunits.K)
        H_CO2_H2O = exp(170.7126 + (-8477.711*pyunits.K/T) + \
                        (-21.95743*log(T))
            + (0.005781*(1/pyunits.K)*T))
        
        return H_CO2_H2O*pyunits.Pa

# Build fugacity based equilibrium equation for CO2    
# y.P = x.rho.He_CO2
# rho - molar density of liquid phase solution
# He_CO2 - Henrys constant of CO2 in liquid phase (N2O Analogy)
class fugacity():
    
    @staticmethod
    def build_parameters(cobj):
        pass
    
    def return_expression(b, phase1, phase2, j):
        pp = (phase1, phase2)
        return (b.params.get_phase('Vap')
                .config.equation_of_state.fug_phase_comp_eq(b, 'Vap', j, pp) ==
                fug_phase_comp_eq(b, 'Liq', j, pp))
    
def fug_phase_comp_eq(b, p, j, pp):
    return _fug_phase_comp(b, p, j, b._teq[pp], \
                           b.mole_frac_phase_comp['Liq', 'MEA'], \
                               b.mole_frac_phase_comp['Liq', 'H2O'], \
                                   b.mole_frac_phase_comp['Liq', 'CO2'])

def _fug_phase_comp(b, p, j, T, xmea, xh2o, xco2):
    pobj = b.params.get_phase(p)
    if pobj.is_liquid_phase():
        return (b.mole_frac_phase_comp[p, j] *\
                dens_mol_phase.return_expression(b,p)*
                    henrys_constant.return_expression(
                        b, cobj(b, j), T, xmea, xh2o, xco2))
    else:
        raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))

# Molar density of liquid phase (m3/mol) 
# Data source [5], based on Eqn 24, table 5      
class dens_mol_phase():
    
    @staticmethod
    def build_parameters(cobj):
        pass
    
    @staticmethod
    def return_expression(b, p):
        pobj = b.params.get_phase(p)
        if pobj.is_liquid_phase():
            v_soln = ((b.mole_frac_phase_comp['Liq','MEA']*\
                      vol_mol_liq_comp.return_expression(
                          b, cobj(b,'MEA'), b.temperature)) + \
                          (b.mole_frac_phase_comp['Liq','H2O']*\
                              vol_mol_liq_comp.return_expression(
                                   b, cobj(b,'H2O'), b.temperature))+
                      (10.2074*(pyunits.cm**3/pyunits.mol)* \
                          b.mole_frac_phase_comp['Liq','CO2']) + \
                              ((-2.2642*(pyunits.cm**3/pyunits.mol) + \
                               3.0059*(pyunits.cm**3/pyunits.mol)* \
                                   b.mole_frac_phase_comp['Liq','MEA'])
                      * b.mole_frac_phase_comp['Liq','MEA']*\
                          b.mole_frac_phase_comp['Liq','H2O']) +
                      ((207*(pyunits.cm**3/pyunits.mol) + \
                       (-563.3701)*(pyunits.cm**3/pyunits.mol) * \
                           b.mole_frac_phase_comp['Liq','MEA']) *
                      b.mole_frac_phase_comp['Liq','MEA']*\
                          b.mole_frac_phase_comp['Liq','CO2']))*\
                      1e-6*pyunits.m**3/pyunits.cm**3
        return 1/v_soln  