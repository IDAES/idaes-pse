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
Methods for eNRTL activity coefficient method.

Only applicable to liquid/electrolyte phases

Reference:

Song, Y. and Chen, C.-C., Symmetric Electrolyte Nonrandom Two-Liquid Activity
Coefficient Model, Ind. Eng. Chem. Res., 2009, Vol. 48, pgs. 7788–7797
"""
from pyomo.environ import Expression, exp, log, Set, units as pyunits

from .eos_base import EoSBase
from .enrtl_parameters import ConstantAlpha, ConstantTau
from idaes.generic_models.properties.core.generic.utility import (
    get_method, get_component_object as cobj)
from idaes.core.util.constants import Constants
from idaes.core.util.exceptions import BurntToast
import idaes.logger as idaeslog


# Set up logger
_log = idaeslog.getLogger(__name__)


DefaultAlphaRule = ConstantAlpha
DefaultTauRule = ConstantTau

# Closest appraoch parameter - implemented as a global constant for now
# This is not something the user should be changing in most cases
ClosestApproach = 14.9


class ENRTL(EoSBase):
    # Add attribute indicating support for electrolyte systems
    electrolyte_support = True

    def build_parameters(b):
        # Build additional indexing sets
        pblock = b.parent_block()
        ion_pair = []
        for i in pblock.cation_set:
            for j in pblock.anion_set:
                ion_pair.append(i+", "+j)
        b.ion_pair_set = Set(initialize=ion_pair)

        comps = pblock.solvent_set | pblock.solute_set | b.ion_pair_set
        comp_pairs = []
        comp_pairs_sym = []
        for i in comps:
            for j in comps:
                if i in pblock.solvent_set | pblock.solute_set or i != j:
                    comp_pairs.append((i, j))
                    if (j, i) not in comp_pairs_sym:
                        comp_pairs_sym.append((i, j))
        b.component_pair_set = Set(initialize=comp_pairs)
        b.component_pair_set_symmetric = Set(initialize=comp_pairs_sym)

        # Check options for alpha rule
        if (b.config.equation_of_state_options is not None and
                "alpha_rule" in b.config.equation_of_state_options):
            b.config.equation_of_state_options[
                "alpha_rule"].build_parameters(b)
        else:
            DefaultAlphaRule.build_parameters(b)

        # Check options for tau rule
        if (b.config.equation_of_state_options is not None and
                "tau_rule" in b.config.equation_of_state_options):
            b.config.equation_of_state_options[
                "tau_rule"].build_parameters(b)
        else:
            DefaultTauRule.build_parameters(b)

    def common(b, pobj):
        pname = pobj.local_name

        molecular_set = b.params.solvent_set | b.params.solute_set

        # Check options for alpha rule
        if (pobj.config.equation_of_state_options is not None and
                "alpha_rule" in pobj.config.equation_of_state_options):
            alpha_rule = pobj.config.equation_of_state_options[
                "alpha_rule"].return_expression
        else:
            alpha_rule = DefaultAlphaRule.return_expression

        # Check options for tau rule
        if (pobj.config.equation_of_state_options is not None and
                "tau_rule" in pobj.config.equation_of_state_options):
            tau_rule = pobj.config.equation_of_state_options[
                "tau_rule"].return_expression
        else:
            tau_rule = DefaultTauRule.return_expression

        # Ionic Strength
        def rule_I(b):
            return (0.5*sum(b.mole_frac_phase_comp_true[pname, c] *
                            b.params.get_component(c).config.charge**2
                            for c in b.params.cation_set) +
                    0.5*sum(b.mole_frac_phase_comp_true[pname, a] *
                            b.params.get_component(a).config.charge**2
                            for a in b.params.anion_set))
        b.add_component(pname+"_ionic_strength",
                        Expression(rule=rule_I,
                                   doc="Ionic strength"))

        # Average molar volume of solvent
        def rule_vol_mol_solvent(b):
            return (sum(b.mole_frac_phase_comp_true[pname, s] /
                        get_method(b, "dens_mol_liq_comp", s)(
                            b, cobj(b, s), b.temperature)
                        for s in molecular_set) /
                    sum(b.mole_frac_phase_comp_true[pname, s]
                        for s in molecular_set))
        b.add_component(pname+"_vol_mol_solvent",
                        Expression(rule=rule_vol_mol_solvent,
                                   doc="Mean molar volume of solvent"))

        # Mean relative permitivity of solvent
        def rule_eps_solvent(b):
            return (sum(b.mole_frac_phase_comp_true[pname, s] *
                        get_method(b, "relative_permittivity_liq_comp", s)(
                            b, cobj(b, s), b.temperature) *
                        b.params.get_component(s).mw
                        for s in molecular_set) /
                    sum(b.mole_frac_phase_comp_true[pname, s] *
                        b.params.get_component(s).mw
                        for s in molecular_set))
        b.add_component(pname+"_relative_permittivity_solvent",
                        Expression(
                            rule=rule_eps_solvent,
                            doc="Mean relative permittivity  of solvent"))

        # Debye-Huckel parameter
        def rule_A_DH(b):
            # Note: Where the paper referes to the dielectric constant, it
            # actually means the electric permittivity of the solver
            # eps = eps_r*eps_0 (units F/m)
            v = pyunits.convert(getattr(b, pname+"_vol_mol_solvent"),
                                pyunits.m**3/pyunits.mol)
            eps = getattr(b, pname+"_relative_permittivity_solvent")
            eps0 = Constants.vacuum_electric_permittivity
            return ((1/3)*(2*Constants.pi*Constants.avogadro_number/v)**0.5 *
                    (Constants.elemental_charge**2/(
                        eps*eps0*Constants.boltzmann_constant *
                        b.temperature))**(3/2)
                    )
        b.add_component(pname+"_A_DH",
                        Expression(
                            rule=rule_A_DH,
                            doc="Debye-Huckel parameter"))

        # Long-range (PDH) contribution to activity coefficient
        def rule_log_gamma_pdh(b, j):
            A = getattr(b, pname+"_A_DH")
            I = getattr(b, pname+"_ionic_strength")
            rho = ClosestApproach
            if j in molecular_set:
                # Note typo in original paper. Correct power for I is (3/2)
                return (2*A*I**(3/2)/(1+rho*I**(1/2)))
            elif j in b.params.cation_set:
                return 200
            elif j in b.params.anion_set:
                return 300
            else:
                raise BurntToast(
                    "{} eNRTL model encountered unexpected component."
                    .format(b.name))
        b.add_component(
            pname+"_log_gamma_pdh",
            Expression(b.params.true_species_set,
                       rule=rule_log_gamma_pdh,
                       doc="Long-range contribution to activity coefficient"))

        # Calculate mixing factors
        def rule_X(b, j):
            if (pname, j) not in b.params.true_phase_component_set:
                return Expression.Skip
            elif j in b.params.cation_set or j in b.params.anion_set:
                cobj = b.params.get_component(j)
                return (b.mole_frac_phase_comp_true[pname, j] *
                        cobj.config.charge)
            else:
                return b.mole_frac_phase_comp_true[pname, j]

        b.add_component(pname+"_X",
                        Expression(b.params.true_species_set,
                                   rule=rule_X,
                                   doc="Charge x mole fraction term"))

        def rule_Y(b, j):
            cobj = b.params.get_component(j)
            if cobj.config.charge < 0:
                # Anion
                dom = b.params.anion_set
            else:
                dom = b.params.cation_set

            X = getattr(b, pname+"_X")
            return X[j]/sum(X[i] for i in dom)

        b.add_component(pname+"_Y",
                        Expression(b.params.ion_set,
                                   rule=rule_Y,
                                   doc="Charge composition"))

        # Calculate alphas for all true species pairings
        def rule_alpha_expr(b, i, j):
            Y = getattr(b, pname+"_Y")
            if ((pname, i) not in b.params.true_phase_component_set or
                    (pname, j) not in b.params.true_phase_component_set):
                return Expression.Skip
            elif ((i in molecular_set) and
                    (j in molecular_set)):
                return alpha_rule(b, pobj, i, j, b.temperature)
            elif (i in b.params.cation_set and j in molecular_set):
                return sum(
                    Y[k] *
                    alpha_rule(
                        b, pobj, (i+", "+k), j, b.temperature)
                    for k in b.params.anion_set)
            elif (j in b.params.cation_set and i in molecular_set):
                return sum(
                    Y[k] *
                    alpha_rule(
                        b, pobj, (j+", "+k), i, b.temperature)
                    for k in b.params.anion_set)
            elif (i in b.params.anion_set and j in molecular_set):
                return sum(
                    Y[k] *
                    alpha_rule(
                        b, pobj, (k+", "+i), j, b.temperature)
                    for k in b.params.cation_set)
            elif (j in b.params.anion_set and i in molecular_set):
                return sum(
                    Y[k] *
                    alpha_rule(
                        b, pobj, (k+", "+j), i, b.temperature)
                    for k in b.params.cation_set)
            elif (i in b.params.cation_set and j in b.params.anion_set):
                return sum(Y[k]*alpha_rule(
                               b, pobj, (i+", "+j), (k+", "+j), b.temperature)
                           for k in b.params.cation_set
                           if i != k)
            elif (i in b.params.anion_set and j in b.params.cation_set):
                return sum(Y[k]*alpha_rule(
                               b, pobj, (j+", "+i), (j+", "+k), b.temperature)
                           for k in b.params.anion_set
                           if i != k)
            elif ((i in b.params.cation_set and j in b.params.cation_set) or
                  (i in b.params.anion_set and j in b.params.anion_set)):
                # No like-ion interactions
                return Expression.Skip
            else:
                raise BurntToast(
                    "{} eNRTL model encountered unexpected component pair {}."
                    .format(b.name, (i, j)))

        b.add_component(pname+"_alpha",
                        Expression(b.params.true_species_set,
                                   b.params.true_species_set,
                                   rule=rule_alpha_expr,
                                   doc="Non-randomness parameters"))

        # Calculate G terms
        def rule_G_expr(b, i, j):
            Y = getattr(b, pname+"_Y")

            def _G_appr(b, pobj, i, j, T):
                return exp(-alpha_rule(b, pobj, i, j, T) *
                           tau_rule(b, pobj, i, j, T))

            if ((pname, i) not in b.params.true_phase_component_set or
                    (pname, j) not in b.params.true_phase_component_set):
                return Expression.Skip
            elif ((i in molecular_set) and
                    (j in molecular_set)):
                return _G_appr(b, pobj, i, j, b.temperature)
            elif (i in b.params.cation_set and j in molecular_set):
                return sum(
                    Y[k] * _G_appr(b, pobj, (i+", "+k), j,  b.temperature)
                    for k in b.params.anion_set)
            elif (j in b.params.cation_set and i in molecular_set):
                return sum(
                    Y[k] * _G_appr(b, pobj, (j+", "+k), i, b.temperature)
                    for k in b.params.anion_set)
            elif (i in b.params.anion_set and j in molecular_set):
                return sum(
                    Y[k] * _G_appr(b, pobj, (k+", "+i), j, b.temperature)
                    for k in b.params.cation_set)
            elif (j in b.params.anion_set and i in molecular_set):
                return sum(
                    Y[k] * _G_appr(b, pobj, (k+", "+j), i, b.temperature)
                    for k in b.params.cation_set)
            elif (i in b.params.cation_set and j in b.params.anion_set):
                return sum(Y[k] * _G_appr(
                    b, pobj, (i+", "+j), (k+", "+j), b.temperature)
                    for k in b.params.cation_set
                    if i != k)
            elif (i in b.params.anion_set and j in b.params.cation_set):
                return sum(Y[k] * _G_appr(
                    b, pobj, (j+", "+i), (j+", "+k), b.temperature)
                    for k in b.params.anion_set
                    if i != k)
            elif ((i in b.params.cation_set and j in b.params.cation_set) or
                  (i in b.params.anion_set and j in b.params.anion_set)):
                # No like-ion interactions
                return Expression.Skip
            else:
                raise BurntToast(
                    "{} eNRTL model encountered unexpected component pair {}."
                    .format(b.name, (i, j)))
        b.add_component(pname+"_G",
                        Expression(b.params.true_species_set,
                                   b.params.true_species_set,
                                   rule=rule_G_expr,
                                   doc="Local interaction G term"))

        # Calculate tau terms
        def rule_tau_expr(b, i, j):
            if ((pname, i) not in b.params.true_phase_component_set or
                    (pname, j) not in b.params.true_phase_component_set):
                return Expression.Skip
            elif ((i in molecular_set) and
                    (j in molecular_set)):
                return tau_rule(b, pobj, i, j, b.temperature)
            elif ((i in b.params.cation_set and j in b.params.cation_set) or
                  (i in b.params.anion_set and j in b.params.anion_set)):
                # No like-ion interactions
                return Expression.Skip
            else:
                alpha = getattr(b, pname+"_alpha")
                G = getattr(b, pname+"_G")
                return log(G[i, j])/alpha[i, j]
        b.add_component(pname+"_tau",
                        Expression(b.params.true_species_set,
                                   b.params.true_species_set,
                                   rule=rule_tau_expr,
                                   doc="Binary interaction energy parameters"))

        # Local contribution to activity coefficient
        # Indicies in expressions use same names as source paper
        # mp = m'
        aqu_species = b.params.true_species_set - b.params._non_aqueous_set

        def rule_log_gamma_lc(b, s):
            X = getattr(b, pname+"_X")
            G = getattr(b, pname+"_G")
            tau = getattr(b, pname+"_tau")
            if (pname, s) not in b.params.true_phase_component_set:
                # Non-aqueous component
                return Expression.Skip
            if s in b.params.cation_set:
                c = s
                Z = b.params.get_component(c).config.charge
                return Z*(
                    sum((X[m]*G[c, m] /
                         sum(X[i]*G[i, m] for i in aqu_species)) *
                        (tau[c, m] -
                         (sum(X[i]*G[i, m]*tau[i, m] for i in aqu_species) /
                          sum(X[i]*G[i, m] for i in aqu_species)))
                        for m in molecular_set) +
                    sum(X[i]*G[i, c]*tau[i, c]
                        for i in (aqu_species-b.params.cation_set)) /
                    sum(X[i]*G[i, c]
                        for i in (aqu_species-b.params.cation_set)) +
                    sum((X[a]*G[c, a] /
                         sum(X[i]*G[i, a]
                             for i in (aqu_species-b.params.anion_set))) *
                        (tau[c, a] -
                         sum(X[i]*G[i, a]*tau[i, a]
                             for i in (aqu_species-b.params.anion_set)) /
                         sum(X[i]*G[i, a]
                             for i in (aqu_species-b.params.anion_set)))
                        for a in b.params.anion_set))
            elif s in b.params.anion_set:
                a = s
                Z = b.params.get_component(a).config.charge
                return Z*(
                    sum((X[m]*G[a, m] /
                         sum(X[i]*G[i, m] for i in aqu_species)) *
                        (tau[a, m] -
                         (sum(X[i]*G[i, m]*tau[i, m] for i in aqu_species) /
                          sum(X[i]*G[i, m] for i in aqu_species)))
                        for m in molecular_set) +
                    sum(X[i]*G[i, a]*tau[i, a]
                        for i in (aqu_species-b.params.anion_set)) /
                    sum(X[i]*G[i, a]
                        for i in (aqu_species-b.params.anion_set)) +
                    sum((X[c]*G[a, c] /
                         sum(X[i]*G[i, c]
                             for i in (aqu_species-b.params.cation_set))) *
                        (tau[a, c] -
                         sum(X[i]*G[i, c]*tau[i, c]
                             for i in (aqu_species-b.params.cation_set)) /
                         sum(X[i]*G[i, c]
                             for i in (aqu_species-b.params.cation_set)))
                        for c in b.params.cation_set))
            else:
                m = s
                return (sum(X[m]*G[i, m]*tau[i, m] for i in aqu_species) /
                        sum(X[m]*G[i, m] for i in aqu_species) +
                        sum((X[mp]*G[m, mp] /
                             sum(X[i]*G[i, mp] for i in aqu_species)) *
                            (tau[m, mp] -
                             (sum(X[i]*G[i, mp]*tau[i, mp]
                                  for i in aqu_species) /
                              sum(X[i]*G[i, mp] for i in aqu_species)))
                            for mp in molecular_set) +
                        sum((X[c]*G[m, c] /
                             sum(X[i]*G[i, c]
                                 for i in (aqu_species-b.params.cation_set))) *
                            (tau[m, c] -
                             (sum(X[i]*G[i, c]*tau[i, c]
                                  for i in (aqu_species-b.params.cation_set)) /
                              sum(X[i]*G[i, c]
                                  for i in (aqu_species-b.params.cation_set))))
                            for c in b.params.cation_set) +
                        sum((X[a]*G[m, a] /
                             sum(X[i]*G[i, a]
                                 for i in (aqu_species-b.params.anion_set))) *
                            (tau[m, a] -
                             (sum(X[i]*G[i, a]*tau[i, a]
                                  for i in (aqu_species-b.params.anion_set)) /
                              sum(X[i]*G[i, a]
                                  for i in (aqu_species-b.params.anion_set))))
                            for a in b.params.anion_set)
                        )
        b.add_component(pname+"_log_gamma_lc",
                        Expression(
                            b.params.true_species_set,
                            rule=rule_log_gamma_lc,
                            doc="Local contribution to activity coefficient"))

    @staticmethod
    def calculate_scaling_factors(b, pobj):
        pass

    @staticmethod
    def dens_mol_phase(b, p):
        return 55e3

    @staticmethod
    def enth_mol_phase(b, p):
        return 1e2*b.temperature

    @staticmethod
    def enth_mol_phase_comp(b, p, j):
        return 1e2*b.temperature
