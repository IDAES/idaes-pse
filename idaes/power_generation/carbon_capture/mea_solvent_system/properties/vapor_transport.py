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
Methods for Vapor Phase Transport Properties

References:

[1]. Perry and Green Handbook; McGraw Hill,8th edition 2008
[2].Wilke, C. R. "A viscosity equation for gas mixtures.
    " The journal of chemical physics 18.4 (1950): 517-519.
[3]. Fuller, Edward N., Paul D. Schettler, and J. Calvin Giddings.
   "New method for prediction of binary gas-phase diffusion coefficients."
   IEC 58.5 (1966): 18-27.

"""
from pyomo.environ import exp, log, Var, Expression, units as pyunits

from idaes.core.util.misc import set_param_from_config
from idaes.generic_models.properties.core.generic.utility import get_method


class ViscSutherland():
    # Sutherland Method for computing vapor component viscosity
    @staticmethod
    def build_parameters(cobj):
        cobj.visc_d_comp_coeff_a = Var(
            doc="Parameter a for vapor phase viscosity model",
            units=pyunits.Pa * pyunits.s)
        set_param_from_config(cobj, param="visc_d_comp_coeff", index="a")

        cobj.visc_d_comp_coeff_b = Var(
            doc="Parameter b for vapor phase viscosity model",
            units=pyunits.K)
        set_param_from_config(cobj, param="visc_d_comp_coeff", index="b")

        cobj.visc_d_comp_coeff_c = Var(
            doc="Parameter c for vapor phase viscosity model",
            units=pyunits.K)
        set_param_from_config(cobj, param="visc_d_comp_coeff", index="c")

    @staticmethod
    def return_expression(blk, cobj, T):
        a = cobj.visc_d_comp_coeff_a
        b = cobj.visc_d_comp_coeff_b
        c = cobj.visc_d_comp_coeff_c

        visc = a * (b + c) / (T + c) * (T / b)**1.5
        return visc


class ViscDIPPR():
    # DIPPR method for computing Dynamic vicosity of vapor components Ref: #1
    @staticmethod
    def build_parameters(cobj):
        cobj.visc_d_comp_coeff_a = Var(
            doc="Parameter a for vapor phase viscosity model",
            units=pyunits.Pa * pyunits.s)
        set_param_from_config(cobj, param="visc_d_comp_coeff", index="a")

        cobj.visc_d_comp_coeff_b = Var(
            doc="Parameter b for vapor phase viscosity model",
            units=pyunits.dimensionless)
        set_param_from_config(cobj, param="visc_d_comp_coeff", index="b")

        cobj.visc_d_comp_coeff_c = Var(
            doc="Parameter c for vapor phase viscosity model",
            units=pyunits.K)
        set_param_from_config(cobj, param="visc_d_comp_coeff", index="c")

    @staticmethod
    def return_expression(blk, cobj, T):
        a = cobj.visc_d_comp_coeff_a
        b = cobj.visc_d_comp_coeff_b
        c = cobj.visc_d_comp_coeff_c

        visc = a * T**b / (1 + c / T)
        return visc / pyunits.K**b


class ViscVapor():
    # Vapor phase dynamic viscosity (Wilke,1950), Ref: #2
    @staticmethod
    def return_expression(blk, p):
        pobj = blk.params.get_phase(p)

        bin_set = []
        for i in pobj.parent_block().component_list:
            for j in pobj.parent_block().component_list:
                bin_set.append((i, j))

        pobj.thetha_ij = Expression(bin_set,
                                    doc='viscosity interaction parameter')

        comp = [i for i in pobj.parent_block().component_list]
        o = dict()
        for (i, j) in enumerate(comp, 1):
            o[i] = j

        def thetha_ij_expression(mw, visc):
            for i in range(1, len(comp)):
                for j in range(i + 1, len(comp) + 1):
                    pobj.thetha_ij[o[i], o[j]] =\
                        (1 + 2 * (visc[p, o[i]] / visc[p, o[j]])**0.5 *
                         (mw[o[j]] / mw[o[i]])**0.25 + visc[p, o[i]] / visc[p, o[j]] *
                         (mw[o[j]] / mw[o[i]])**0.5) /\
                        (8 + 8 * mw[o[i]] / mw[o[j]])**0.5

                    pobj.thetha_ij[o[j], o[i]] = (visc[p, o[j]] / visc[p, o[i]] *
                                                  mw[o[i]] / mw[o[j]] *
                                                  pobj.thetha_ij[o[i], o[j]])
            for i in pobj.parent_block().component_list:
                for j in pobj.parent_block().component_list:
                    if i == j:
                        pobj.thetha_ij[i, j] = 1

        thetha_ij_expression(blk.mw_comp, blk.visc_d_comp)

        mu_vap = sum(blk.mole_frac_comp[i] * blk.visc_d_comp[p, i] /
                     sum(blk.mole_frac_comp[j] * pobj.thetha_ij[i, j]
                         for j in pobj.parent_block().component_list)
                     for i in pobj.parent_block().component_list)
        return mu_vap


class ThermCondComp():
    # Method for computing vapor component thermal conductivity Ref: #1
    @staticmethod
    def build_parameters(cobj):
        cobj.therm_cond_comp_coeff_a = Var(
            doc="Parameter a for vapor phase thermal conductivity model",
            units=pyunits.W / pyunits.m / pyunits.K)
        set_param_from_config(cobj, param="therm_cond_comp_coeff", index="a")

        cobj.therm_cond_comp_coeff_b = Var(
            doc="Parameter b for vapor phase thermal conductivity model",
            units=pyunits.dimensionless)
        set_param_from_config(cobj, param="therm_cond_comp_coeff", index="b")

        cobj.therm_cond_comp_coeff_c = Var(
            doc="Parameter c for vapor phase thermal conductivity model",
            units=pyunits.K)
        set_param_from_config(cobj, param="therm_cond_comp_coeff", index="c")

        cobj.therm_cond_comp_coeff_d = Var(
            doc="Parameter d for vapor phase thermal conductivity model",
            units=pyunits.K**2)
        set_param_from_config(cobj, param="therm_cond_comp_coeff", index="d")

    @staticmethod
    def return_expression(blk, cobj, T):
        a = cobj.therm_cond_comp_coeff_a
        b = cobj.therm_cond_comp_coeff_b
        c = cobj.therm_cond_comp_coeff_c
        d = cobj.therm_cond_comp_coeff_d

        therm_cond = a * (T**b) / ((1 + c / T) + d / T**2)
        return therm_cond / pyunits.K**b


class ThermCond():
    # Vapor phase Thermal Conductivty
    # Wassiljewa-Mason-Saxena mixing rule(low pressure)
    @staticmethod
    def return_expression(blk, p):
        pobj = blk.params.get_phase(p)

        def Wassiljewa_Mason_Saxena_mixing_rule(x, mw, visc, kcomp):
            k_vap = 0
            for i in pobj.parent_block().component_list:
                sumij = 0
                for j in pobj.parent_block().component_list:
                    Aij = ((1 + (visc[p, i] / visc[p, j])**0.5 *
                            (mw[j] / mw[i])**0.25)**2 *
                           (8 * (1 + mw[i] / mw[j]))**-0.5)
                    sumij += x[p, j] * Aij
                k_vap += x[p, i] * kcomp[p, i] / sumij
            return k_vap

        return Wassiljewa_Mason_Saxena_mixing_rule(blk.mole_frac_phase_comp,
                                                   blk.mw_comp,
                                                   blk.visc_d_comp,
                                                   blk.therm_cond_comp)


class Diffus():
    @staticmethod
    def build_parameters(cobj, phase):
        pass

    @staticmethod
    def return_expression(blk, p, i, T):
        # Diffusion Coefficient(binary) constants - Ref: #3
        # Diffusion volumes in Fuller-Schettler-Giddings correlation
        # for estimating binary diffusivities of components in vapor phase
        # Reference: Table 3.1 pp 71 Seader Henley (2006)

        pobj = blk.params.get_phase(p)

        binary_set = []
        for ii in pobj.parent_block().component_list:
            for jj in pobj.parent_block().component_list:
                if ii != jj and (jj, ii) not in binary_set:
                    binary_set.append((ii, jj))

        def diffus_binary(i, j):
            return (1.013e-2 * pyunits.Pa / pyunits.K**1.75 / pyunits.s *
                    T**1.75 / blk.pressure *
                    (1e-3 * (1 / blk.mw_comp[i] + 1 / blk.mw_comp[j]))**0.5 /
                    (blk.params.get_component(i).diffus_volume**(0.333333) +
                     blk.params.get_component(j).diffus_volume**(0.333333))**2)

        return ((1 - blk.mole_frac_phase_comp[p, i]) /
                (sum(blk.mole_frac_phase_comp[p, j] / diffus_binary(i, j)
                     for j in blk.component_list if (i, j) in binary_set) +
                 sum(blk.mole_frac_phase_comp[p, j] / diffus_binary(j, i)
                     for j in blk.component_list if (j, i) in binary_set)))

class TransportMethod():
    ViscDIPPR = ViscDIPPR
    ViscSutherland = ViscSutherland
    ViscVapor = ViscVapor
    ThermCondComp = ThermCondComp
    ThermCond = ThermCond
    Diffus = Diffus
