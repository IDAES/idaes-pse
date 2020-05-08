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
from idaes.generic_models.properties.core.generic.generic_property import (
    get_method, get_component_object as cobj)
from .eos_base import EoSBase
from idaes import bin_directory
from idaes.core.util.constants import Constants as const
import idaes.logger as idaeslog
from idaes.core.util.exceptions import BurntToast


# Set up logger
_log = idaeslog.getLogger(__name__)


# Set path to root finder .so file
_so = os.path.join(bin_directory, "cubic_roots.so")


def cubic_roots_available():
    """Make sure the compiled cubic root functions are available. Yes, in
    Windows the .so extention is still used.
    """
    return os.path.isfile(_so)


class CubicType(Enum):
    PR = 0
    SRK = 1


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
    def common(b, pobj):
        ctype = pobj._cubic_type
        cname = pobj.config.equation_of_state_options["type"].name

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
                        "{} received unrecognised cubic type. This should "
                        "never happen, so please contact the IDAES developers "
                        "with this bug.".format(b.name))

        b.add_component(cname+'_fw',
                        Expression(b.params.component_list,
                                   rule=func_fw,
                                   doc='EoS S factor'))

        def func_a(m, j):
            cobj = m.params.get_component(j)
            fw = getattr(m, cname+"_fw")
            return (EoS_param[ctype]['omegaA']*(
                       (const.gas_constant *
                        cobj.temperature_crit)**2/cobj.pressure_crit) *
                    ((1+fw[j]*(1-sqrt(m.temperature /
                                      cobj.temperature_crit)))**2))
        b.add_component(cname+'_a',
                        Expression(b.params.component_list,
                                   rule=func_a,
                                   doc='Component a coefficient'))

        def func_b(m, j):
            cobj = m.params.get_component(j)
            return (EoS_param[ctype]['coeff_b'] * const.gas_constant *
                    cobj.temperature_crit/cobj.pressure_crit)
        b.add_component(cname+'_b',
                        Expression(b.params.component_list,
                                   rule=func_b,
                                   doc='Component b coefficient'))

        def rule_am(m, p):
            a = getattr(m, cname+"_a")
            k = getattr(m.params, cname+"_kappa")
            return sum(sum(
                m.mole_frac_phase_comp[p, i]*m.mole_frac_phase_comp[p, j] *
                sqrt(a[i]*a[j])*(1-k[i, j])
                for j in m.params.component_list)
                for i in m.params.component_list)
        b.add_component(cname+'_am',
                        Expression(b.params.phase_list, rule=rule_am))

        def rule_bm(m, p):
            b = getattr(m, cname+"_b")
            return sum(m.mole_frac_phase_comp[p, i]*b[i]
                       for i in m.params.component_list)
        b.add_component(cname+'_bm',
                        Expression(b.params.phase_list, rule=rule_bm))

        def rule_A(m, p):
            am = getattr(m, cname+"_am")
            return (am[p]*m.pressure /
                    (const.gas_constant*m.temperature)**2)
        b.add_component(cname+'_A',
                        Expression(b.params.phase_list, rule=rule_A))

        def rule_B(m, p):
            bm = getattr(m, cname+"_bm")
            return (bm[p]*m.pressure /
                    (const.gas_constant*m.temperature))
        b.add_component(cname+'_B',
                        Expression(b.params.phase_list, rule=rule_B))

        def rule_delta(m, p, i):
            # See pg. 145 in Properties of Gases and Liquids
            a = getattr(m, cname+"_a")
            am = getattr(m, cname+"_am")
            kappa = getattr(m.params, cname+"_kappa")
            return (2*sqrt(a[i])/am[p] *
                    sum(m.mole_frac_phase_comp[p, j]*sqrt(a[j]) *
                        (1-kappa[i, j])
                        for j in m.params.component_list))
        b.add_component(cname+"_delta",
                        Expression(b.params.phase_list,
                                   b.params.component_list,
                                   rule=rule_delta))

        def rule_dadT(m, p):
            # See pg. 102 in Properties of Gases and Liquids
            a = getattr(m, cname+"_a")
            fw = getattr(m, cname+"_fw")
            kappa = getattr(m.params, cname+"_kappa")
            return -((const.gas_constant/2)*sqrt(EoS_param[ctype]['omegaA']) *
                     sum(sum(m.mole_frac_phase_comp[p, i] *
                             m.mole_frac_phase_comp[p, j] *
                             (1-kappa[i, j]) *
                             (fw[j]*sqrt(a[i] *
                              m.params.get_component(j).temperature_crit /
                              m.params.get_component(j).pressure_crit) +
                              fw[i]*sqrt(a[j] *
                              m.params.get_component(i).temperature_crit /
                              m.params.get_component(i).pressure_crit))
                             for j in m.params.component_list)
                         for i in m.params.component_list) /
                     sqrt(m.temperature))
        b.add_component(cname+"_dadT",
                        Expression(b.params.phase_list,
                                   rule=rule_dadT))

        # Add components at equilibrium state if required
        if (b.params.config.phases_in_equilibrium is not None and
                (not b.config.defined_state or b.always_flash)):
            def func_a_eq(m, p1, p2, j):
                cobj = m.params.get_component(j)
                fw = getattr(m, cname+"_fw")
                return (EoS_param[ctype]['omegaA']*(
                            (const.gas_constant *
                             cobj.temperature_crit)**2/cobj.pressure_crit) *
                        ((1+fw[j]*(1-sqrt(m._teq[p1, p2] /
                                          cobj.temperature_crit)))**2))
            b.add_component('_'+cname+'_a_eq',
                            Expression(b.params._pe_pairs,
                                       b.params.component_list,
                                       rule=func_a_eq,
                                       doc='Component a coefficient at Teq'))

            def rule_am_eq(m, p1, p2, p3):
                a = getattr(m, '_'+cname+"_a_eq")
                k = getattr(m.params, cname+"_kappa")
                return sum(sum(
                    m.mole_frac_phase_comp[p3, i] *
                    m.mole_frac_phase_comp[p3, j] *
                    sqrt(a[p1, p2, i]*a[p1, p2, j])*(1-k[i, j])
                    for j in m.params.component_list)
                    for i in m.params.component_list)
            b.add_component('_'+cname+'_am_eq',
                            Expression(b.params._pe_pairs,
                                       b.params.phase_list,
                                       rule=rule_am_eq))

            def rule_A_eq(m, p1, p2, p3):
                am_eq = getattr(m, "_"+cname+"_am_eq")
                return (am_eq[p1, p2, p3]*m.pressure /
                        (const.gas_constant*m._teq[p1, p2])**2)
            b.add_component('_'+cname+'_A_eq',
                            Expression(b.params._pe_pairs,
                                       b.params.phase_list,
                                       rule=rule_A_eq))

            def rule_B_eq(m, p1, p2, p3):
                bm = getattr(m, cname+"_bm")
                return (bm[p3]*m.pressure /
                        (const.gas_constant*m._teq[p1, p2]))
            b.add_component('_'+cname+'_B_eq',
                            Expression(b.params._pe_pairs,
                                       b.params.phase_list,
                                       rule=rule_B_eq))

            def rule_delta_eq(m, p1, p2, p3, i):
                # See pg. 145 in Properties of Gases and Liquids
                a = getattr(m, "_"+cname+"_a_eq")
                am = getattr(m, "_"+cname+"_am_eq")
                kappa = getattr(m.params, cname+"_kappa")
                return (2*sqrt(a[p1, p2, i])/am[p1, p2, p3] *
                        sum(m.mole_frac_phase_comp[p3, j]*sqrt(a[p1, p2, j]) *
                            (1-kappa[i, j])
                            for j in m.params.component_list))
            b.add_component("_"+cname+"_delta_eq",
                            Expression(b.params._pe_pairs,
                                       b.params.phase_list,
                                       b.params.component_list,
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

    def build_parameters(b):
        b._cubic_type = b.config.equation_of_state_options["type"]
        cname = b._cubic_type.name
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
                doc=cname+' binary interaction parameters'))

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

    def dens_mass_phase(b, p):
        return b.dens_mol_phase[p]*b.mw_phase[p]

    def dens_mol_phase(b, p):
        pobj = b.params.get_phase(p)
        if pobj.is_vapor_phase():
            return b.pressure/(const.gas_constant*b.temperature)
        elif pobj.is_liquid_phase():
            return sum(b.mole_frac_phase_comp[p, j] *
                       get_method(b, "dens_mol_liq_comp", j)(
                           b, cobj(b, j), b.temperature)
                       for j in b.components_in_phase(p))
        else:
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))

    def enth_mol_phase(b, p):
        return sum(b.mole_frac_phase_comp[p, j]*b.enth_mol_phase_comp[p, j]
                   for j in b.components_in_phase(p))

    def enth_mol_phase_comp(b, p, j):
        pobj = b.params.get_phase(p)
        if pobj.is_vapor_phase():
            return get_method(b, "enth_mol_ig_comp", j)(
                b, cobj(b, j), b.temperature)
        elif pobj.is_liquid_phase():
            return get_method(b, "enth_mol_liq_comp", j)(
                b, cobj(b, j), b.temperature)
        else:
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))

    def entr_mol_phase(b, p):
        return sum(b.mole_frac_phase_comp[p, j]*b.entr_mol_phase_comp[p, j]
                   for j in b.components_in_phase(p))

    def entr_mol_phase_comp(b, p, j):
        pobj = b.params.get_phase(p)
        if pobj.is_vapor_phase():
            return get_method(b, "entr_mol_ig_comp", j)(
                b, cobj(b, j), b.temperature)
        elif pobj.is_liquid_phase():
            return get_method(b, "entr_mol_liq_comp", j)(
                b, cobj(b, j), b.temperature)
        else:
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))

    def fug_phase_comp(b, p, j, pp):
        pobj = b.params.get_phase(p)
        if pobj.is_vapor_phase():
            return b.mole_frac_phase_comp[p, j]*b.pressure
        elif pobj.is_liquid_phase():
            return b.mole_frac_phase_comp[p, j] * \
                   get_method(b, "pressure_sat_comp", j)(
                       b, cobj(b, j), b._teq[pp])
        else:
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))

    def fug_coeff_phase_comp(blk, p, j):
        pobj = blk.params.get_phase(p)
        if not (pobj.is_vapor_phase() or pobj.is_liquid_phase()):
            raise PropertyNotSupportedError(_invalid_phase_msg(blk.name, p))

        cname = pobj._cubic_type.name
        b = getattr(blk, cname+"_b")[j]
        bm = getattr(blk, cname+"_bm")[p]
        A = getattr(blk, cname+"_A")[p]
        B = getattr(blk, cname+"_B")[p]
        delta = getattr(blk, cname+"_delta")[p, j]
        Z = blk.compress_fact_phase[p]

        EoS_u = EoS_param[pobj._cubic_type]['u']
        EoS_w = EoS_param[pobj._cubic_type]['w']
        EoS_p = sqrt(EoS_u**2 - 4*EoS_w)

        return exp((b/bm*(Z-1)*(B*EoS_p) - log(Z-B)*(B*EoS_p) +
                    A*(b/bm - delta) *
                    log((2*Z + B*(EoS_u + EoS_p))/(2*Z + B*(EoS_u - EoS_p)))) /
                   (B*EoS_p))

    # def fug_coeff_phase_comp(blk, p, j):
    #     pobj = blk.params.get_phase(p)
    #     if not (pobj.is_vapor_phase() or pobj.is_liquid_phase()):
    #         raise PropertyNotSupportedError(_invalid_phase_msg(blk.name, p))

    #     cname = pobj._cubic_type.name
    #     b = getattr(blk, cname+"_b")
    #     bm = getattr(blk, cname+"_bm")
    #     Aeq = getattr(blk, "_"+cname+"_A_eq")
    #     Beq = getattr(blk, "_"+cname+"_B_eq")
    #     delta_eq = getattr(blk, "_"+cname+"_delta_eq")

    #     EoS_u = EoS_param[pobj._cubic_type]['u']
    #     EoS_w = EoS_param[pobj._cubic_type]['w']
    #     EoS_p = sqrt(EoS_u**2 - 4*EoS_w)

    #     f = getattr(blk, "_"+cname+"_ext_func_param")
    #     if pobj.is_vapor_phase():
    #         proc = getattr(blk, "_"+cname+"_proc_Z_vap")
    #     elif pobj.is_liquid_phase():
    #         proc = getattr(blk, "_"+cname+"_proc_Z_liq")

    #     def Zeq(p):
    #         return proc(f, Aeq[p], Beq[p])

    #     return exp((b[j]/bm[p]*(Zeq(p)-1)*(Beq[p]*EoS_p) -
    #                 log(Zeq(p)-Beq[p])*(Beq[p]*EoS_p) +
    #                 Aeq[p]*(b[j]/bm[p] - delta_eq[p, j]) *
    #                 log((2*Zeq(p) + Beq[p]*(EoS_u + EoS_p)) /
    #                     (2*Zeq(p) + Beq[p]*(EoS_u - EoS_p)))) /
    #                (Beq[p]*EoS_p))

    def gibbs_mol_phase(b, p):
        return sum(b.mole_frac_phase_comp[p, j]*b.gibbs_mol_phase_comp[p, j]
                   for j in b.components_in_phase(p))

    def gibbs_mol_phase_comp(b, p, j):
        return (b.enth_mol_phase_comp[p, j] -
                b.entr_mol_phase_comp[p, j] *
                b.temperature)


def _invalid_phase_msg(name, phase):
    return ("{} received unrecognised phase name {}. Ideal property "
            "libray only supports Vap and Liq phases."
            .format(name, phase))
