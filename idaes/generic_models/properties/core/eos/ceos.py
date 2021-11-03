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
Methods for cubic equations of state.

Currently only supports liquid and vapor phases
"""
import os
from enum import Enum

from pyomo.environ import (exp,
                           Expression,
                           ExternalFunction,
                           log,
                           Param,
                           Reals,
                           sqrt,
                           Var)
from pyomo.common.config import ConfigBlock, ConfigValue, In

from idaes.core.util.exceptions import PropertyNotSupportedError
from idaes.generic_models.properties.core.generic.utility import (
    get_method, get_component_object as cobj)
from idaes.core.util.math import safe_log
from .eos_base import EoSBase
from idaes import bin_directory
import idaes.logger as idaeslog
from idaes.core.util.exceptions import \
    BurntToast, ConfigurationError, PropertyNotSupportedError


# Set up logger
_log = idaeslog.getLogger(__name__)


# Set path to root finder .so file
_so = os.path.join(bin_directory, "cubic_roots.so")

"""
References:
[1]. Poling, B.E., Prausnitz, J.M. and O’connell, J.P., 2001. 
     Properties of gases and liquids. McGraw-Hill Education.

[2]. Trujillo, M.F., O'Rourke, P. and Torres, D., 2002.
     Generalizing the Thermodynamics State Relationships in KIVA-3V
     https://www.osti.gov/servlets/purl/809947 (Last accessed: 08/13/2021)
"""


def cubic_roots_available():
    """Make sure the compiled cubic root functions are available. Yes, in
    Windows the .so extention is still used.
    """
    return os.path.isfile(_so)


class CubicType(Enum):
    PR = 0
    SRK = 1


class MixingRuleA(Enum):
    # Rule to calculate am for cubic equations of state
    default = 0


class MixingRuleB(Enum):
    # Rule to calculate bm for cubic equations of state
    default = 0


EoS_param = {
        CubicType.PR: {'u': 2, 'w': -1, 'omegaA': 0.45724, 'coeff_b': 0.07780},
        CubicType.SRK: {'u': 1, 'w': 0, 'omegaA': 0.42748, 'coeff_b': 0.08664}
        }


CubicConfig = ConfigBlock()
CubicConfig.declare("type", ConfigValue(
    domain=In(CubicType),
    description="Equation of state to use",
    doc="Enum indicating type of cubic equation of state to use."))


class Cubic(EoSBase):

    @staticmethod
    def common(b, pobj):
        # TODO: determine if Henry's Law applies to Cubic EoS systems
        # For now, raise an exception if found
        # Follow on questions:
        # If Henry's law is used for a component, how does that effect
        # calculating A, B and phi?
        for j in b.component_list:
            cobj = b.params.get_component(j)
            if (cobj.config.henry_component is not None and
                    pobj.local_name in cobj.config.henry_component):
                raise PropertyNotSupportedError(
                    "{} Cubic equations of state do not support Henry's "
                    "components [{}, {}].".format(b.name, pobj.local_name, j))

        ctype = pobj._cubic_type
        cname = pobj.config.equation_of_state_options["type"].name
        mixing_rule_a = pobj._mixing_rule_a
        mixing_rule_b = pobj._mixing_rule_b

        if hasattr(b, cname+"_fw"):
            # Common components already constructed by previous phase
            return

        # Create expressions for coefficients
        def func_fw(m, j):
            cobj = m.params.get_component(j)
            if ctype == CubicType.PR:
                return 0.37464 + 1.54226*cobj.omega - \
                       0.26992*cobj.omega**2
            elif ctype == CubicType.SRK:
                return 0.48 + 1.574*cobj.omega - \
                       0.176*cobj.omega**2
            else:
                raise BurntToast(
                        "{} received unrecognized cubic type. This should "
                        "never happen, so please contact the IDAES developers "
                        "with this bug.".format(b.name))

        b.add_component(cname+'_fw',
                        Expression(b.component_list,
                                   rule=func_fw,
                                   doc='EoS S factor'))
        
        def func_a(m, j):
            cobj = m.params.get_component(j)
            fw = getattr(m, cname+"_fw")
            return (EoS_param[ctype]['omegaA']*(
                       (Cubic.gas_constant(b) *
                        cobj.temperature_crit)**2/cobj.pressure_crit) *
                    ((1+fw[j]*(1-sqrt(m.temperature /
                                      cobj.temperature_crit)))**2))
        b.add_component(cname+'_a',
                        Expression(b.component_list,
                                   rule=func_a,
                                   doc='Component a coefficient'))
        
        def func_da_dT(m, j):
            cobj = m.params.get_component(j)
            fw = getattr(m, cname+"_fw")
            
            ac = (EoS_param[ctype]['omegaA']*(
                       (Cubic.gas_constant(b) *
                        cobj.temperature_crit)**2/cobj.pressure_crit))
            Tr = m.temperature/cobj.temperature_crit
            
            return -ac*fw[j]/(cobj.temperature_crit*sqrt(Tr))*(1+fw[j]*(1-sqrt(Tr)))
        b.add_component(cname+'_da_dT',
                        Expression(b.component_list,
                                   rule=func_da_dT,
                                   doc='Temperature derivative of component a'))
        
        def func_d2a_dT2(m, j):
            cobj = m.params.get_component(j)
            fw = getattr(m, cname+"_fw")
            
            ac = (EoS_param[ctype]['omegaA']*(
                       (Cubic.gas_constant(b) *
                        cobj.temperature_crit)**2/cobj.pressure_crit))
            Tr = m.temperature/cobj.temperature_crit
            
            return 1/cobj.temperature_crit**2*ac*fw[j]*(1+fw[j])/(2*Tr**(3/2))
                    
        b.add_component(cname+'_d2a_dT2',
                        Expression(b.component_list,
                                   rule=func_d2a_dT2,
                                   doc='Second temperature derivative'
                                       'of component a'))

        def func_b(m, j):
            cobj = m.params.get_component(j)
            return (EoS_param[ctype]['coeff_b'] * Cubic.gas_constant(b) *
                    cobj.temperature_crit/cobj.pressure_crit)
        b.add_component(cname+'_b',
                        Expression(b.component_list,
                                   rule=func_b,
                                   doc='Component b coefficient'))

        if mixing_rule_a == MixingRuleA.default:
            def rule_am(m, p):
                a = getattr(m, cname+"_a")
                return rule_am_default(m, cname, a, p)
            b.add_component(cname+'_am',
                            Expression(b.phase_list, rule=rule_am))
            
            def rule_daij_dT(m, i, j):
                a = getattr(m, cname+"_a")
                da_dT = getattr(m, cname+"_da_dT")
                k = getattr(m.params, cname+"_kappa")
                
                # Include temperature derivative of k for future extension
                dk_ij_dT = 0
                
                return sqrt(a[i]*a[j])*(-dk_ij_dT
                                        + (1-k[i,j])/2
                                        *(da_dT[i]/a[i]+da_dT[j]/a[j])
                                        )
            b.add_component(cname+'_daij_dT',
                            Expression(b.component_list,
                                       b.component_list,
                                       rule=rule_daij_dT))

            def rule_dam_dT(m, p):
                daij_dT = getattr(m, cname+"_daij_dT")
                return sum(sum(m.mole_frac_phase_comp[p, i]
                                *  m.mole_frac_phase_comp[p, j]
                                 * daij_dT[i,j]
                                 for j in m.components_in_phase(p))
                             for i in m.components_in_phase(p))
            
            b.add_component(cname+"_dam_dT",
                            Expression(b.phase_list,
                                       rule=rule_dam_dT))
            
            
            def rule_d2am_dT2(m,p):
                 k = getattr(m.params, cname+"_kappa")
                 a = getattr(m, cname+"_a")        
                 da_dT = getattr(m, cname+"_da_dT")
                 d2a_dT2 = getattr(m, cname+"_d2a_dT2")
                # Placeholders for if temperature dependent k is needed
                 dk_dT = 0
                 d2k_dT2 = 0
                 
                 # Initialize loop variable
                 d2am_dT2 = 0
                 
                 for i in m.components_in_phase(p):
                     for j in m.components_in_phase(p):
                         d2aij_dT2 = (
                             sqrt(a[i]*a[j])
                             * (-d2k_dT2- dk_dT*(da_dT[i]/a[i] + da_dT[j]/a[j])
                                + (1-k[i,j])/2
                                * (d2a_dT2[i]/a[i] + d2a_dT2[j]/a[j]
                                   - 1/2*(da_dT[i]/a[i]-da_dT[j]/a[j])**2
                                   )
                                )
                             )
                         d2am_dT2 += (m.mole_frac_phase_comp[p, i]
                                 * m.mole_frac_phase_comp[p, j]
                                 * d2aij_dT2)
                 return d2am_dT2
            b.add_component(cname+"_d2am_dT2",
                            Expression(b.phase_list,
                                       rule=rule_d2am_dT2))                                
                 
            
            def rule_delta(m, p, i):
                # See pg. 145 in Properties of Gases and Liquids
                a = getattr(m, cname+"_a")
                am = getattr(m, cname+"_am")
                kappa = getattr(m.params, cname+"_kappa")
                return (2*sqrt(a[i])/am[p] *
                        sum(m.mole_frac_phase_comp[p, j]*sqrt(a[j]) *
                            (1-kappa[i, j])
                            for j in b.components_in_phase(p)))
            b.add_component(cname+"_delta",
                            Expression(b.phase_component_set,
                                       rule=rule_delta))
            
        else:
            raise ConfigurationError(
                "{} Unrecognized option for Equation of State "
                "mixing_rule_a: {}. Must be an instance of MixingRuleA "
                "Enum.".format(b.name, mixing_rule_a))        
        
        
        if mixing_rule_b == MixingRuleB.default:
            def rule_bm(m, p):
                b = getattr(m, cname+"_b")
                return rule_bm_default(m, b, p)
            b.add_component(cname+'_bm',
                            Expression(b.phase_list, rule=rule_bm))
        else:
            raise ConfigurationError(
                "{} Unrecognized option for Equation of State "
                "mixing_rule_a: {}. Must be an instance of MixingRuleB "
                "Enum.".format(b.name, mixing_rule_b))

        def rule_A(m, p):
            am = getattr(m, cname+"_am")
            return (am[p]*m.pressure /
                    (Cubic.gas_constant(b)*m.temperature)**2)
        b.add_component(cname+'_A',
                        Expression(b.phase_list, rule=rule_A))

        def rule_B(m, p):
            bm = getattr(m, cname+"_bm")
            return (bm[p]*m.pressure /
                    (Cubic.gas_constant(b)*m.temperature))
        b.add_component(cname+'_B',
                        Expression(b.phase_list, rule=rule_B))
        
        

        # Add components at equilibrium state if required
        if (b.params.config.phases_in_equilibrium is not None and
                (not b.config.defined_state or b.always_flash)):
            def func_a_eq(m, p1, p2, j):
                cobj = m.params.get_component(j)
                fw = getattr(m, cname+"_fw")
                return (EoS_param[ctype]['omegaA']*(
                            (Cubic.gas_constant(b) *
                             cobj.temperature_crit)**2/cobj.pressure_crit) *
                        ((1+fw[j]*(1-sqrt(m._teq[p1, p2] /
                                          cobj.temperature_crit)))**2))
            b.add_component('_'+cname+'_a_eq',
                            Expression(b.params._pe_pairs,
                                       b.component_list,
                                       rule=func_a_eq,
                                       doc='Component a coefficient at Teq'))

            def rule_am_eq(m, p1, p2, p3):
                try:
                    rule = m.params.get_phase(p3).config.equation_of_state_options[
                        "mixing_rule_a"]
                except (KeyError, TypeError):
                    rule = MixingRuleA.default

                a = getattr(m, "_"+cname+"_a_eq")
                if rule == MixingRuleA.default:
                    return rule_am_default(m, cname, a, p3, (p1, p2))
                else:
                    raise ConfigurationError(
                        "{} Unrecognized option for Equation of State "
                        "mixing_rule_a: {}. Must be an instance of MixingRuleA "
                        "Enum.".format(m.name, rule))
            b.add_component('_'+cname+'_am_eq',
                            Expression(b.params._pe_pairs,
                                       b.phase_list,
                                       rule=rule_am_eq))

            def rule_A_eq(m, p1, p2, p3):
                am_eq = getattr(m, "_"+cname+"_am_eq")
                return (am_eq[p1, p2, p3]*m.pressure /
                        (Cubic.gas_constant(b)*m._teq[p1, p2])**2)
            b.add_component('_'+cname+'_A_eq',
                            Expression(b.params._pe_pairs,
                                       b.phase_list,
                                       rule=rule_A_eq))

            def rule_B_eq(m, p1, p2, p3):
                bm = getattr(m, cname+"_bm")
                return (bm[p3]*m.pressure /
                        (Cubic.gas_constant(b)*m._teq[p1, p2]))
            b.add_component('_'+cname+'_B_eq',
                            Expression(b.params._pe_pairs,
                                       b.phase_list,
                                       rule=rule_B_eq))

            def rule_delta_eq(m, p1, p2, p3, i):
                # See pg. 145 in Properties of Gases and Liquids
                a = getattr(m, "_"+cname+"_a_eq")
                am = getattr(m, "_"+cname+"_am_eq")
                kappa = getattr(m.params, cname+"_kappa")
                return (2*sqrt(a[p1, p2, i])/am[p1, p2, p3] *
                        sum(m.mole_frac_phase_comp[p3, j]*sqrt(a[p1, p2, j]) *
                            (1-kappa[i, j])
                            for j in m.components_in_phase(p3)))
            b.add_component("_"+cname+"_delta_eq",
                            Expression(b.params._pe_pairs,
                                       b.phase_component_set,
                                       rule=rule_delta_eq))

        # Set up external function calls
        b.add_component("_"+cname+"_ext_func_param",
                        Param(default=ctype.value))
        b.add_component("_"+cname+"_proc_Z_liq",
                        ExternalFunction(library=_so,
                                         function="ceos_z_liq"))
        b.add_component("_"+cname+"_proc_Z_vap",
                        ExternalFunction(library=_so,
                                         function="ceos_z_vap"))

    @staticmethod
    def calculate_scaling_factors(b, pobj):
        pass

    @staticmethod
    def build_parameters(b):
        b._cubic_type = b.config.equation_of_state_options["type"]
  
        cname = b._cubic_type.name
        try:
            mixing_rule_a = b.config.equation_of_state_options["mixing_rule_a"]
        except (KeyError, TypeError):
            mixing_rule_a = MixingRuleA.default
        
        if mixing_rule_a == MixingRuleA.default:
            b._mixing_rule_a = mixing_rule_a
        else:
            raise ConfigurationError(
                "{} Unrecognized option for Equation of State "
                "mixing_rule_a: {}. Must be an instance of MixingRuleA "
                "Enum.".format(b.name, mixing_rule_a))
        
        try:
            mixing_rule_b = b.config.equation_of_state_options["mixing_rule_b"]
        except (KeyError, TypeError):
            mixing_rule_b = MixingRuleB.default

        if mixing_rule_b == MixingRuleB.default:
            b._mixing_rule_b = mixing_rule_b
        else:
            raise ConfigurationError(
                "{} Unrecognized option for Equation of State "
                "mixing_rule_a: {}. Must be an instance of MixingRuleB "
                "Enum.".format(b.name, mixing_rule_b))
        
        param_block = b.parent_block()

        if hasattr(param_block, cname+"_kappa"):
            # Common components already constructed by previous phase
            return

        kappa_data = param_block.config.parameter_data[cname+"_kappa"]
        param_block.add_component(
            cname+'_kappa',
            Var(param_block.component_list,
                param_block.component_list,
                within=Reals,
                initialize=kappa_data,
                doc=cname+' binary interaction parameters',
                units=None))

    @staticmethod
    def compress_fact_phase(b, p):
        pobj = b.params.get_phase(p)
        cname = pobj._cubic_type.name
        A = getattr(b, cname+"_A")
        B = getattr(b, cname+"_B")
        f = getattr(b, "_"+cname+"_ext_func_param")
        if pobj.is_vapor_phase():
            proc = getattr(b, "_"+cname+"_proc_Z_vap")
        elif pobj.is_liquid_phase():
            proc = getattr(b, "_"+cname+"_proc_Z_liq")
        else:
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))
        return proc(f, A[p], B[p])   
    
    @staticmethod
    def cp_mol_phase(blk, p):
        pobj = blk.params.get_phase(p)
        cname = pobj._cubic_type.name

        am = getattr(blk, cname+"_am")[p]
        bm = getattr(blk, cname+"_bm")[p]
        B = getattr(blk, cname+"_B")[p]

        dam_dT = getattr(blk, cname+"_dam_dT")[p]
        d2am_dT2 = getattr(blk, cname+"_d2am_dT2")[p]
        
        T = blk.temperature
        R = Cubic.gas_constant(blk)
        Z = blk.compress_fact_phase[p]
        dZdT = _dZ_dT(blk,p)

        EoS_u = EoS_param[pobj._cubic_type]['u']
        EoS_w = EoS_param[pobj._cubic_type]['w']
        EoS_p = sqrt(EoS_u**2 - 4*EoS_w)

        # d2adT2 = - (0.5/T) * dam_dT + \
        #          ((R**2 * omegaA) / (2*T)) * (
        #              sum(sum(
        #                  blk.mole_frac_phase_comp[p, i] * blk.mole_frac_phase_comp[p, j]
        #                  * (1 - kappa[i, j]) * fw[i] * fw[j] * sqrt(
        #                      (blk.params.get_component(i).temperature_crit * blk.params.get_component(j).temperature_crit) /
        #                      (blk.params.get_component(i).pressure_crit * blk.params.get_component(j).pressure_crit)
        #                  )
        #                  for j in blk.components_in_phase(p))
        #                  for i in blk.components_in_phase(p))
        #          )

        expression1 = 2*Z + (EoS_u + EoS_p)*B
        expression2 = 2*Z + (EoS_u - EoS_p)*B
        expression3 = B*(dZdT + Z/T)/(Z**2 + Z*EoS_u*B + EoS_w*B**2)

        cp_ideal_gas =  sum(blk.mole_frac_phase_comp[p, j] 
             * get_method(blk, "cp_mol_ig_comp", j)(blk, cobj(blk, j), T)
             for j in blk.components_in_phase(p))
        
        # Derived from the relations in Chapter 6 of [1]
        cp_departure =  (R*(T*dZdT + Z - 1) 
                      + (T*d2am_dT2/(EoS_p*bm))*safe_log(expression1/expression2,
                                                         eps = 1e-6)
                      + ((am - T*dam_dT)*expression3/bm))

        
        return cp_ideal_gas + cp_departure

    @staticmethod
    def cv_mol_phase(blk, p):
        pobj = blk.params.get_phase(p)
        cname = pobj._cubic_type.name
        am = getattr(blk, cname+"_am")[p]
        bm = getattr(blk, cname+"_bm")[p]
        cp = blk.cp_mol_phase[p]
        V = 1 / blk.dens_mol_phase[p]
        dam_dT = getattr(blk, cname+"_dam_dT")[p]

        EoS_u = EoS_param[pobj._cubic_type]['u']
        EoS_w = EoS_param[pobj._cubic_type]['w']

        dPdV = - ((Cubic.gas_constant(blk) * blk.temperature) / (V - bm)**2 ) + \
                (am * (2 * V + EoS_u * bm) / (V**2 + EoS_u * bm * V + EoS_w * bm**2)**2)

        dPdT = (Cubic.gas_constant(blk) / (V - bm)) - \
               (1 / (V**2 + EoS_u * bm * V + EoS_w * bm**2)) * dam_dT

        # See Chapter 6 in [1]
        return (
            cp + blk.temperature * dPdT**2 / dPdV
        )


    @staticmethod
    def dens_mass_phase(b, p):
        return b.dens_mol_phase[p]*b.mw_phase[p]
    

    @staticmethod
    def dens_mol_phase(b, p):
        pobj = b.params.get_phase(p)
        if pobj.is_vapor_phase() or pobj.is_liquid_phase():
            return b.pressure/(
                Cubic.gas_constant(b)*b.temperature*b.compress_fact_phase[p])
        else:
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))
    

    @staticmethod
    def energy_internal_mol_phase(blk, p):
        pobj = blk.params.get_phase(p)
        if not (pobj.is_vapor_phase() or pobj.is_liquid_phase()):
            raise PropertyNotSupportedError(_invalid_phase_msg(blk.name, p))

        cname = pobj._cubic_type.name
        am = getattr(blk, cname+"_am")[p]
        bm = getattr(blk, cname+"_bm")[p]
        B = getattr(blk, cname+"_B")[p]
        dam_dT = getattr(blk, cname+"_dam_dT")[p]
        Z = blk.compress_fact_phase[p]

        EoS_u = EoS_param[pobj._cubic_type]['u']
        EoS_w = EoS_param[pobj._cubic_type]['w']
        EoS_p = sqrt(EoS_u**2 - 4*EoS_w)

        # Derived from equation on pg. 120 in Properties of Gases and Liquids
        # Departure function for U is similar to H minus the RT(Z-1) term
        return (((blk.temperature*dam_dT - am) *
                 safe_log((2*Z + B*(EoS_u+EoS_p)) / (2*Z + B*(EoS_u-EoS_p)),
                          eps=1e-6)) / (bm*EoS_p) +
                sum(blk.mole_frac_phase_comp[p, j] *
                    EoSBase.energy_internal_mol_ig_comp_pure(blk, j)
                    for j in blk.components_in_phase(p)))

    @staticmethod
    def energy_internal_mol_phase_comp(blk, p, j):
        pobj = blk.params.get_phase(p)
        if not (pobj.is_vapor_phase() or pobj.is_liquid_phase()):
            raise PropertyNotSupportedError(_invalid_phase_msg(blk.name, p))

        cname = pobj._cubic_type.name
        am = getattr(blk, cname+"_am")[p]
        bm = getattr(blk, cname+"_bm")[p]
        B = getattr(blk, cname+"_B")[p]
        dam_dT = getattr(blk, cname+"_dam_dT")[p]
        Z = blk.compress_fact_phase[p]

        EoS_u = EoS_param[pobj._cubic_type]['u']
        EoS_w = EoS_param[pobj._cubic_type]['w']
        EoS_p = sqrt(EoS_u**2 - 4*EoS_w)

        # Derived from equation on pg. 120 in Properties of Gases and Liquids
        # Departure function for U is similar to H minus the RT(Z-1) term
        return (((blk.temperature*dam_dT - am) *
                 safe_log((2*Z + B*(EoS_u+EoS_p)) / (2*Z + B*(EoS_u-EoS_p)),
                          eps=1e-6)) / (bm*EoS_p) +
                EoSBase.energy_internal_mol_ig_comp_pure(blk, j))

    @staticmethod
    def enth_mol_phase(blk, p):
        pobj = blk.params.get_phase(p)
        if not (pobj.is_vapor_phase() or pobj.is_liquid_phase()):
            raise PropertyNotSupportedError(_invalid_phase_msg(blk.name, p))

        cname = pobj._cubic_type.name
        am = getattr(blk, cname+"_am")[p]
        bm = getattr(blk, cname+"_bm")[p]
        B = getattr(blk, cname+"_B")[p]
        dam_dT = getattr(blk, cname+"_dam_dT")[p]
        Z = blk.compress_fact_phase[p]

        EoS_u = EoS_param[pobj._cubic_type]['u']
        EoS_w = EoS_param[pobj._cubic_type]['w']
        EoS_p = sqrt(EoS_u**2 - 4*EoS_w)

        # Derived from equation on pg. 120 in Properties of Gases and Liquids
        return (((blk.temperature*dam_dT - am) *
                 safe_log((2*Z + B*(EoS_u+EoS_p)) / (2*Z + B*(EoS_u-EoS_p)),
                          eps=1e-6) +
                 Cubic.gas_constant(blk)*blk.temperature*(Z-1)*bm*EoS_p) /
                (bm*EoS_p) + sum(blk.mole_frac_phase_comp[p, j] *
                                 get_method(blk, "enth_mol_ig_comp", j)(
                                            blk, cobj(blk, j), blk.temperature)
                                 for j in blk.components_in_phase(p)))

    @staticmethod
    def enth_mol_phase_comp(blk, p, j):
        pobj = blk.params.get_phase(p)
        if not (pobj.is_vapor_phase() or pobj.is_liquid_phase()):
            raise PropertyNotSupportedError(_invalid_phase_msg(blk.name, p))

        cname = pobj._cubic_type.name
        am = getattr(blk, cname+"_am")[p]
        bm = getattr(blk, cname+"_bm")[p]
        B = getattr(blk, cname+"_B")[p]
        dam_dT = getattr(blk, cname+"_dam_dT")[p]
        Z = blk.compress_fact_phase[p]

        EoS_u = EoS_param[pobj._cubic_type]['u']
        EoS_w = EoS_param[pobj._cubic_type]['w']
        EoS_p = sqrt(EoS_u**2 - 4*EoS_w)

        # Derived from equation on pg. 120 in Properties of Gases and Liquids
        return (((blk.temperature*dam_dT - am) *
                 safe_log((2*Z + B*(EoS_u+EoS_p)) / (2*Z + B*(EoS_u-EoS_p)),
                          eps=1e-6) +
                 Cubic.gas_constant(blk)*blk.temperature*(Z-1)*bm*EoS_p) /
                (bm*EoS_p) + get_method(blk, "enth_mol_ig_comp", j)(
                                        blk, cobj(blk, j), blk.temperature))

    @staticmethod
    def entr_mol_phase(blk, p):
        pobj = blk.params.get_phase(p)
        if not (pobj.is_vapor_phase() or pobj.is_liquid_phase()):
            raise PropertyNotSupportedError(_invalid_phase_msg(blk.name, p))

        cname = pobj._cubic_type.name
        bm = getattr(blk, cname+"_bm")[p]
        B = getattr(blk, cname+"_B")[p]
        dam_dT = getattr(blk, cname+"_dam_dT")[p]
        Z = blk.compress_fact_phase[p]

        EoS_u = EoS_param[pobj._cubic_type]['u']
        EoS_w = EoS_param[pobj._cubic_type]['w']
        EoS_p = sqrt(EoS_u**2 - 4*EoS_w)
        
        R = Cubic.gas_constant(blk)
        
        entr_ideal_gas = -R*safe_log(blk.pressure
                                 /blk.params.pressure_ref, eps=1e-6)
        for j in blk.components_in_phase(p):
            entr_j = get_method(blk, "entr_mol_ig_comp", j)(
                blk, cobj(blk, j), blk.temperature)
            xj = blk.mole_frac_phase_comp[p, j]
            
            entr_ideal_gas += xj*(entr_j - R*safe_log(xj,eps=1e-6))
        
        # See pg. 102 in Properties of Gases and Liquids
        # or pg. 208 of Sandler, 4th Ed.
        entr_departure = (
            R*safe_log((Z-B), eps=1e-6)
            +  dam_dT/(bm*EoS_p)*safe_log((2*Z + B*(EoS_u + EoS_p)) /
                                         (2*Z + B*(EoS_u - EoS_p)),eps=1e-6)
            )
        
        return entr_ideal_gas + entr_departure

    @staticmethod
    def entr_mol_phase_comp(blk, p, j):
        pobj = blk.params.get_phase(p)
        if not (pobj.is_vapor_phase() or pobj.is_liquid_phase()):
            raise PropertyNotSupportedError(_invalid_phase_msg(blk.name, p))

        cname = pobj._cubic_type.name
        
        bm = getattr(blk, cname+"_bm")[p]
        B = getattr(blk, cname+"_B")[p]
        dam_dT = getattr(blk, cname+"_dam_dT")[p]
        
        b = getattr(blk, cname+"_b")
        bm = getattr(blk, cname+"_bm")
        A = getattr(blk, cname+"_A")
        delta = getattr(blk, cname+"_delta")
        Z = blk.compress_fact_phase[p]
        
        logphi_j = _log_fug_coeff_phase_comp(blk,p,j)
        dlogphi_j_dT = _d_log_fug_coeff_dT_phase_comp(blk,p,j)
        Z = blk.compress_fact_phase[p]

        EoS_u = EoS_param[pobj._cubic_type]['u']
        EoS_w = EoS_param[pobj._cubic_type]['w']
        EoS_p = sqrt(EoS_u**2 - 4*EoS_w)
        
        R = Cubic.gas_constant(blk)
        
        entr_ideal_gas = (get_method(blk, "entr_mol_ig_comp", j)(
                              blk, cobj(blk, j), blk.temperature)
                            - R*(safe_log(blk.pressure
                                 /blk.params.pressure_ref, eps=1e-6)
                                 + safe_log(blk.mole_frac_phase_comp[p, j],
                                            eps=1e-6)
                                 )
                            )

        # Departure function for partial molar entropy can be obtained by
        # writing the departure function for the Gibbs energy in terms of the
        # log of the fugacity coefficient, and differentiating that expression.
        # This is mixing-rule dependent
        
        entr_departure = -R*logphi_j - R*blk.temperature*dlogphi_j_dT
        
        return entr_ideal_gas + entr_departure

    @staticmethod
    def fug_phase_comp(b, p, j):
        pobj = b.params.get_phase(p)
        if pobj.is_vapor_phase() or pobj.is_liquid_phase():
            return (b.mole_frac_phase_comp[p, j] *
                    b.pressure *
                    b.fug_coeff_phase_comp[p, j])
        else:
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))
    
    @staticmethod
    def fug_phase_comp_eq(b, p, j, pp):
        pobj = b.params.get_phase(p)
        if pobj.is_vapor_phase() or pobj.is_liquid_phase():
            return (b.mole_frac_phase_comp[p, j] *
                    b.pressure *
                    exp(_log_fug_coeff_phase_comp_eq(b, p, j, pp)))
        else:
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))


    @staticmethod
    def log_fug_phase_comp_eq(b, p, j, pp):
        pobj = b.params.get_phase(p)
        if pobj.is_vapor_phase() or pobj.is_liquid_phase():
            return (log(b.mole_frac_phase_comp[p, j]) +
                    log(b.pressure/b.params.pressure_ref) +
                    _log_fug_coeff_phase_comp_eq(b, p, j, pp))
        else:
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))
    
    @staticmethod
    def fug_coeff_phase_comp(blk, p, j):
        pobj = blk.params.get_phase(p)
        ctype = pobj._cubic_type

        if not (pobj.is_vapor_phase() or pobj.is_liquid_phase()):
            raise PropertyNotSupportedError(_invalid_phase_msg(blk.name, p))

        cname = pobj._cubic_type.name
        b = getattr(blk, cname+"_b")[j]
        bm = getattr(blk, cname+"_bm")[p]
        A = getattr(blk, cname+"_A")[p]
        B = getattr(blk, cname+"_B")[p]
        delta = getattr(blk, cname+"_delta")[p, j]
        Z = blk.compress_fact_phase[p]

        return exp(_log_fug_coeff_method(A, b, bm, B, delta, Z, ctype))

    @staticmethod
    def fug_coeff_phase_comp_eq(blk, p, j, pp):
        return exp(_log_fug_coeff_phase_comp_eq(blk, p, j, pp))

    @staticmethod
    def log_fug_phase_comp_Tbub(blk, p, j, pp):
        pobj = blk.params.get_phase(p)
        ctype = pobj._cubic_type
        cname = pobj.config.equation_of_state_options["type"].name

        if pobj.is_liquid_phase():
            x = blk.mole_frac_comp
            xidx = ()
        elif pobj.is_vapor_phase():
            x = blk._mole_frac_tbub
            xidx = pp
        else:
            raise BurntToast("{} non-vapor or liquid phase called for bubble "
                             "temperature calculation. This should never "
                             "happen, so please contact the IDAES developers "
                             "with this bug.".format(blk.name))

        def a(k):
            cobj = blk.params.get_component(k)
            fw = getattr(blk, cname+"_fw")[k]
            return (EoS_param[ctype]['omegaA'] *
                    ((Cubic.gas_constant(blk) * cobj.temperature_crit)**2 /
                     cobj.pressure_crit) *
                    ((1+fw*(1-sqrt(blk.temperature_bubble[pp] /
                                   cobj.temperature_crit)))**2))

        kappa = getattr(blk.params, cname+"_kappa")
        am = sum(sum(x[xidx, i]*x[xidx, j]*sqrt(a(i)*a(j))*(1-kappa[i, j])
                     for j in blk.component_list)
                 for i in blk.component_list)

        b = getattr(blk, cname+"_b")
        bm = sum(x[xidx, i]*b[i] for i in blk.component_list)

        A = am*blk.pressure/(Cubic.gas_constant(blk) *
                             blk.temperature_bubble[pp])**2
        B = bm*blk.pressure/(Cubic.gas_constant(blk) *
                             blk.temperature_bubble[pp])

        delta = (2*sqrt(a(j))/am * sum(x[xidx, i]*sqrt(a(i))*(1-kappa[j, i])
                                       for i in blk.component_list))

        f = getattr(blk, "_"+cname+"_ext_func_param")
        if pobj.is_vapor_phase():
            proc = getattr(blk, "_"+cname+"_proc_Z_vap")
        elif pobj.is_liquid_phase():
            proc = getattr(blk, "_"+cname+"_proc_Z_liq")

        Z = proc(f, A, B)

        if pobj.is_vapor_phase():
            mole_frac = blk._mole_frac_tbub[pp[0], pp[1], j]
        else:
            mole_frac = blk.mole_frac_comp[j]

        return (_log_fug_coeff_method(A, b[j], bm, B, delta, Z, ctype) +
                log(mole_frac) + log(blk.pressure/blk.pressure._units))

    @staticmethod
    def log_fug_phase_comp_Tdew(blk, p, j, pp):
        pobj = blk.params.get_phase(p)
        ctype = pobj._cubic_type
        cname = pobj.config.equation_of_state_options["type"].name

        if pobj.is_liquid_phase():
            x = blk._mole_frac_tdew
            xidx = pp
        elif pobj.is_vapor_phase():
            x = blk.mole_frac_comp
            xidx = ()
        else:
            raise BurntToast("{} non-vapor or liquid phase called for bubble "
                             "temperature calculation. This should never "
                             "happen, so please contact the IDAES developers "
                             "with this bug.".format(blk.name))

        def a(k):
            cobj = blk.params.get_component(k)
            fw = getattr(blk, cname+"_fw")[k]
            return (EoS_param[ctype]['omegaA'] *
                    ((Cubic.gas_constant(blk) * cobj.temperature_crit)**2 /
                     cobj.pressure_crit) *
                    ((1+fw*(1-sqrt(blk.temperature_dew[pp] /
                                   cobj.temperature_crit)))**2))

        kappa = getattr(blk.params, cname+"_kappa")
        am = sum(sum(x[xidx, i]*x[xidx, j]*sqrt(a(i)*a(j))*(1-kappa[i, j])
                     for j in blk.component_list)
                 for i in blk.component_list)

        b = getattr(blk, cname+"_b")
        bm = sum(x[xidx, i]*b[i] for i in blk.component_list)

        A = am*blk.pressure/(Cubic.gas_constant(blk) *
                             blk.temperature_dew[pp])**2
        B = bm*blk.pressure/(Cubic.gas_constant(blk) *
                             blk.temperature_dew[pp])

        delta = (2*sqrt(a(j))/am * sum(x[xidx, i]*sqrt(a(i))*(1-kappa[j, i])
                                       for i in blk.component_list))

        f = getattr(blk, "_"+cname+"_ext_func_param")
        if pobj.is_vapor_phase():
            proc = getattr(blk, "_"+cname+"_proc_Z_vap")
        elif pobj.is_liquid_phase():
            proc = getattr(blk, "_"+cname+"_proc_Z_liq")

        Z = proc(f, A, B)

        if pobj.is_vapor_phase():
            mole_frac = blk.mole_frac_comp[j]
        else:
            mole_frac = blk._mole_frac_tdew[pp[0], pp[1], j]

        return (_log_fug_coeff_method(A, b[j], bm, B, delta, Z, ctype) +
                log(mole_frac) + log(blk.pressure/blk.pressure._units))

    @staticmethod
    def log_fug_phase_comp_Pbub(blk, p, j, pp):
        pobj = blk.params.get_phase(p)
        ctype = pobj._cubic_type
        cname = pobj.config.equation_of_state_options["type"].name

        if pobj.is_liquid_phase():
            x = blk.mole_frac_comp
            xidx = ()
        elif pobj.is_vapor_phase():
            x = blk._mole_frac_pbub
            xidx = pp
        else:
            raise BurntToast("{} non-vapor or liquid phase called for bubble "
                             "temperature calculation. This should never "
                             "happen, so please contact the IDAES developers "
                             "with this bug.".format(blk.name))

        a = getattr(blk, cname+"_a")
        kappa = getattr(blk.params, cname+"_kappa")
        am = sum(sum(x[xidx, i]*x[xidx, j] *
                     sqrt(a[i]*a[j])*(1-kappa[i, j])
                     for j in blk.component_list)
                 for i in blk.component_list)

        b = getattr(blk, cname+"_b")
        bm = sum(x[xidx, i]*b[i] for i in blk.component_list)

        A = am*blk.pressure_bubble[pp]/(Cubic.gas_constant(blk) *
                                        blk.temperature)**2
        B = bm*blk.pressure_bubble[pp]/(Cubic.gas_constant(blk) *
                                        blk.temperature)

        delta = (2*sqrt(a[j])/am * sum(x[xidx, i]*sqrt(a[i])*(1-kappa[j, i])
                                       for i in blk.component_list))

        f = getattr(blk, "_"+cname+"_ext_func_param")
        if pobj.is_vapor_phase():
            proc = getattr(blk, "_"+cname+"_proc_Z_vap")
        elif pobj.is_liquid_phase():
            proc = getattr(blk, "_"+cname+"_proc_Z_liq")

        Z = proc(f, A, B)

        if pobj.is_vapor_phase():
            mole_frac = blk._mole_frac_pbub[pp[0], pp[1], j]
        else:
            mole_frac = blk.mole_frac_comp[j]

        return (_log_fug_coeff_method(A, b[j], bm, B, delta, Z, ctype) +
                log(mole_frac) + log(blk.pressure_bubble[pp] /
                                     blk.pressure_bubble._units))

    @staticmethod
    def log_fug_phase_comp_Pdew(blk, p, j, pp):
        pobj = blk.params.get_phase(p)
        ctype = pobj._cubic_type
        cname = pobj.config.equation_of_state_options["type"].name

        if pobj.is_liquid_phase():
            x = blk._mole_frac_pdew
            xidx = pp
        elif pobj.is_vapor_phase():
            x = blk.mole_frac_comp
            xidx = ()
        else:
            raise BurntToast("{} non-vapor or liquid phase called for bubble "
                             "temperature calculation. This should never "
                             "happen, so please contact the IDAES developers "
                             "with this bug.".format(blk.name))

        a = getattr(blk, cname+"_a")
        kappa = getattr(blk.params, cname+"_kappa")
        am = sum(sum(x[xidx, i]*x[xidx, j] *
                     sqrt(a[i]*a[j])*(1-kappa[i, j])
                     for j in blk.component_list)
                 for i in blk.component_list)

        b = getattr(blk, cname+"_b")
        bm = sum(x[xidx, i]*b[i] for i in blk.component_list)

        A = am*blk.pressure_dew[pp]/(Cubic.gas_constant(blk) *
                                     blk.temperature)**2
        B = bm*blk.pressure_dew[pp]/(Cubic.gas_constant(blk)*blk.temperature)

        delta = (2*sqrt(a[j])/am * sum(x[xidx, i]*sqrt(a[i])*(1-kappa[j, i])
                                       for i in blk.component_list))

        f = getattr(blk, "_"+cname+"_ext_func_param")
        if pobj.is_vapor_phase():
            proc = getattr(blk, "_"+cname+"_proc_Z_vap")
        elif pobj.is_liquid_phase():
            proc = getattr(blk, "_"+cname+"_proc_Z_liq")

        Z = proc(f, A, B)

        if pobj.is_vapor_phase():
            mole_frac = blk.mole_frac_comp[j]
        else:
            mole_frac = blk._mole_frac_pdew[pp[0], pp[1], j]

        return (_log_fug_coeff_method(A, b[j], bm, B, delta, Z, ctype) +
                log(mole_frac) + log(blk.pressure_dew[pp] /
                                     blk.pressure_dew._units))

    @staticmethod
    def gibbs_mol_phase(b, p):
        return (b.enth_mol_phase[p] - b.entr_mol_phase[p]*b.temperature)

    @staticmethod
    def gibbs_mol_phase_comp(b, p, j):
        return (b.enth_mol_phase_comp[p, j] -
                b.entr_mol_phase_comp[p, j] *
                b.temperature)


    @staticmethod
    def isentropic_speed_sound_phase(blk, p):
        # See Reference [2]
        return sqrt(blk.heat_capacity_ratio_phase[p]) * blk.isothermal_speed_sound_phase[p]

    
    @staticmethod
    def isothermal_speed_sound_phase(blk, p):
        pobj = blk.params.get_phase(p)
        cname = pobj._cubic_type.name
        am = getattr(blk, cname+"_am")[p]
        bm = getattr(blk, cname+"_bm")[p]
        V = 1 / blk.dens_mol_phase[p]
        mw = blk.mw
        rho = blk.dens_mass_phase[p]

        EoS_u = EoS_param[pobj._cubic_type]['u']
        EoS_w = EoS_param[pobj._cubic_type]['w']

        dPdV = - ((Cubic.gas_constant(blk) * blk.temperature) / (V - bm)**2 ) + \
                (am * (2 * V + EoS_u * bm) / (V**2 + EoS_u * bm * V + EoS_w * bm**2)**2)

        # see reference [2]
        return sqrt(- dPdV * mw / rho**2)

    @staticmethod    
    def vol_mol_phase(b, p):
        pobj = b.params.get_phase(p)
        if pobj.is_vapor_phase() or pobj.is_liquid_phase():
            return (Cubic.gas_constant(b)*b.temperature *
                    b.compress_fact_phase[p] /
                    b.pressure)
        else:
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))
    
    @staticmethod
    def vol_mol_phase_comp(b,p,j):
        pobj = b.params.get_phase(p)
        if pobj.is_vapor_phase() or pobj.is_liquid_phase():
            return (Cubic.gas_constant(b)*b.temperature *
                    _dZ_dxj(b,p,j) /
                    b.pressure)
        else:
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))


def _invalid_phase_msg(name, phase):
    return ("{} received unrecognized phase name {}. Ideal property "
            "libray only supports Vap and Liq phases."
            .format(name, phase))



def _dZ_dT(blk,p):
    pobj = blk.params.get_phase(p)
    cname = pobj._cubic_type.name
    am = getattr(blk, cname+"_am")[p]
    A = getattr(blk, cname+"_A")[p]
    B = getattr(blk, cname+"_B")[p]
    dam_dT = getattr(blk, cname+"_dam_dT")[p]
    Z = blk.compress_fact_phase[p]
    T = blk.temperature

    EoS_u = EoS_param[pobj._cubic_type]['u']
    EoS_w = EoS_param[pobj._cubic_type]['w']
    
    dBdT = -B/T
    dAdT = (A/am) * dam_dT - (2*A/T)
    
    K2 = (EoS_u - 1) * B - 1
    K3 = A - EoS_u*B - (EoS_u - EoS_w)*B**2
    K4 = - (A*B + EoS_w*B**2 + EoS_w*B**3)
    
    dK2dT = (EoS_u - 1) * dBdT
    dK3dT = dAdT - EoS_u*dBdT - 2*(EoS_u - EoS_w)*B*dBdT
    dK4dT = -(A*dBdT + B*dAdT + 2*EoS_w*B*dBdT + 3*EoS_w*B**2*dBdT) 
    
    return -(Z**2*dK2dT + Z*dK3dT + dK4dT)/(3*Z**2 + 2*K2*Z + K3)

def _dZ_dxj(blk,p,j):
    pobj = blk.params.get_phase(p)
    cname = pobj._cubic_type.name
    
    if not (pobj._mixing_rule_a == MixingRuleA.default
            and pobj._mixing_rule_b == MixingRuleB.default):
        # Any user adding more mixing rules will need to explicitly add
        # support for these functions
        raise NotImplementedError("Block {} called for a property "
                                  "that is not supported by this choice "
                                  "of mixing rules.".format(blk.name))
    

    a = getattr(blk, cname+"_a")
    b = getattr(blk, cname+"_b")
    k = getattr(blk, cname+"_kappa")
    am = getattr(blk, cname+"_am")[p]
    bm = getattr(blk, cname+"_bm")[p]
    A = getattr(blk, cname+"_A")[p]
    B = getattr(blk, cname+"_B")[p]
    
    # This derivative takes a somewhat unexpected form because it's a
    # partial molar quantity
    if pobj._mixing_rule_a == MixingRuleA.default:
        dam_dxj = (2/(1-blk.mole_frac_phase_comp[p, j])
                   *(-am + sum(blk.mole_frac_phase_comp[p, i]
                       * (1 - k[i,j])*sqrt(a[i]*a[j])
                               for i in blk.components_in_phase(p))
                     )
                   )
    # Same thing here. 
    if pobj._mixing_rule_b == MixingRuleB.default:
        dbm_dxj = (b[j]-bm)/(1-blk.mole_frac_phase_comp[p, j])
    
    Z = blk.compress_fact_phase[p]
    R = Cubic.gas_constant(blk)
    T = blk.temperature
    P = blk.pressure

    EoS_u = EoS_param[pobj._cubic_type]['u']
    EoS_w = EoS_param[pobj._cubic_type]['w']
    
    dA_dxj = P/(R*T)**2*dam_dxj
    dB_dxj = P/(R*T)*dbm_dxj
    
    
    K2 = (EoS_u - 1) * B - 1
    K3 = A - EoS_u*B - (EoS_u - EoS_w)*B**2
    K4 = - (A*B + EoS_w*B**2 + EoS_w*B**3)
    
    dK2_dxj = (EoS_u - 1) * dB_dxj
    dK3_dxj = dA_dxj - EoS_u*dB_dxj - 2*(EoS_u - EoS_w)*B*dB_dxj
    dK4_dxj = -(A*dB_dxj + B*dA_dxj + 2*EoS_w*B*dB_dxj + 3*EoS_w*B**2*dB_dxj) 
    
    return -(Z**2*dK2_dxj + Z*dK3_dxj + dK4_dxj)/(3*Z**2 + 2*K2*Z + K3)

def _log_fug_coeff_phase_comp_eq(blk, p, j, pp):
    pobj = blk.params.get_phase(p)
    if not (pobj.is_vapor_phase() or pobj.is_liquid_phase()):
        raise PropertyNotSupportedError(_invalid_phase_msg(blk.name, p))

    cname = pobj._cubic_type.name
    b = getattr(blk, cname+"_b")
    bm = getattr(blk, cname+"_bm")
    Aeq = getattr(blk, "_"+cname+"_A_eq")
    Beq = getattr(blk, "_"+cname+"_B_eq")
    delta_eq = getattr(blk, "_"+cname+"_delta_eq")

    f = getattr(blk, "_"+cname+"_ext_func_param")
    if pobj.is_vapor_phase():
        proc = getattr(blk, "_"+cname+"_proc_Z_vap")
    elif pobj.is_liquid_phase():
        proc = getattr(blk, "_"+cname+"_proc_Z_liq")

    def Zeq(p):
        return proc(f, Aeq[pp, p], Beq[pp, p])

    return _log_fug_coeff_method(Aeq[pp, p], b[j], bm[p], Beq[pp, p],
                                 delta_eq[pp, p, j], Zeq(p), pobj._cubic_type)


def _log_fug_coeff_phase_comp(blk,p,j):
    pobj = blk.params.get_phase(p)
    if not (pobj.is_vapor_phase() or pobj.is_liquid_phase()):
        raise PropertyNotSupportedError(_invalid_phase_msg(blk.name, p))
    
    cname = pobj._cubic_type.name
    b = getattr(blk, cname+"_b")
    bm = getattr(blk, cname+"_bm")
    A = getattr(blk, cname+"_A")
    B = getattr(blk, cname+"_B")
    delta = getattr(blk, cname+"_delta")

    f = getattr(blk, "_"+cname+"_ext_func_param")
    if pobj.is_vapor_phase():
        proc = getattr(blk, "_"+cname+"_proc_Z_vap")
    elif pobj.is_liquid_phase():
        proc = getattr(blk, "_"+cname+"_proc_Z_liq")

    def Z(p):
        return proc(f, A[p], B[p])

    return _log_fug_coeff_method(A[p], b[j], bm[p], B[p],
                                 delta[p, j], Z(p), pobj._cubic_type)

def _log_fug_coeff_method(A, b, bm, B, delta, Z, cubic_type):
    u = EoS_param[cubic_type]['u']
    w = EoS_param[cubic_type]['w']
    p = sqrt(u**2 - 4*w)

    return ((b/bm*(Z-1)*(B*p) - safe_log(Z-B, eps=1e-6)*(B*p) +
             A*(b/bm - delta)*safe_log((2*Z + B*(u + p))/(2*Z + B*(u - p)),
                                       eps=1e-6)) /
            (B*p))

def _d_log_fug_coeff_dT_phase_comp(blk,p,j):
    pobj = blk.params.get_phase(p)
    # Is this check necessary?
    if not (pobj.is_vapor_phase() or pobj.is_liquid_phase()):
        raise PropertyNotSupportedError(_invalid_phase_msg(blk.name, p))
    
    if not (pobj._mixing_rule_a == MixingRuleA.default
            and pobj._mixing_rule_b == MixingRuleB.default):
        # Any 
        raise NotImplementedError("Block {} called for a property "
                                  "that is not supported by this choice "
                                  "of mixing rules.".format(blk.name))
    
    cname = pobj._cubic_type.name
    am = getattr(blk, cname+"_am")[p]
    daij_dT = getattr(blk, cname+"_daij_dT")
    dam_dT = getattr(blk, cname+"_dam_dT")[p]
    b = getattr(blk, cname+"_b")[j]
    bm = getattr(blk, cname+"_bm")[p]
    A = getattr(blk, cname+"_A")[p]
    B = getattr(blk, cname+"_B")[p]
    delta = getattr(blk, cname+"_delta")[p,j]
    
    Z = blk.compress_fact_phase[p]
    dZ_dT = _dZ_dT(blk,p)
    T = blk.temperature
    
    u = EoS_param[pobj._cubic_type]['u']
    w = EoS_param[pobj._cubic_type]['w']
    EoS_p = sqrt(u**2 - 4*w)
    
    expr = A/(EoS_p*B)*(b/bm*(1/am*dam_dT-1/T) + delta/T
                    - 2/am*sum(blk.mole_frac_phase_comp[p,i]*daij_dT[i,j]
                               for i in blk.component_list))
    
    return (b/bm*dZ_dT - (dZ_dT+B/T)/(Z-B)
             - A*(b/bm-delta)*(Z/T+dZ_dT)/(Z**2+u*B*Z+w*B**2)
             + log((2*Z+B*(u+EoS_p))/(2*Z+B*(u-EoS_p)))*expr)

# -----------------------------------------------------------------------------
# Mixing rules
def rule_am_default(m, cname, a, p, pp=()):
    k = getattr(m.params, cname+"_kappa")
    return sum(sum(
        m.mole_frac_phase_comp[p, i]*m.mole_frac_phase_comp[p, j] *
        sqrt(a[pp, i]*a[pp, j])*(1-k[i, j])
        for j in m.components_in_phase(p))
        for i in m.components_in_phase(p))

def rule_bm_default(m, b, p):
    return sum(m.mole_frac_phase_comp[p, i]*b[i]
               for i in m.components_in_phase(p))