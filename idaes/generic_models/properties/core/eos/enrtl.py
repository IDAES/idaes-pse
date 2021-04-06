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
"""
from pyomo.environ import Expression, exp, log

from .eos_base import EoSBase
from .enrtl_submethods import ConstantAlpha, ConstantTau
from idaes.core.util.exceptions import BurntToast
import idaes.logger as idaeslog


# Set up logger
_log = idaeslog.getLogger(__name__)


DefaultAlphaRule = ConstantAlpha
DefaultTauRule = ConstantTau


class ENRTL(EoSBase):
    # Add attribute indicating support for electrolyte systems
    electrolyte_support = True

    def build_parameters(b):

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

        # Calculate mixing factors
        def rule_X(b, j):
            if j in b.params.cation_set or j in b.params.anion_set:
                cobj = b.params.get_component(j)
                return (b.mole_frac_phase_comp_true[pobj.local_name, j] *
                        cobj.config.charge)
            else:
                return b.mole_frac_phase_comp_true[pobj.local_name, j]

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
            if ((i in molecular_set) and
                    (j in molecular_set)):
                return alpha_rule(b, pobj, i, j, b.temperature)
            elif (i in b.params.cation_set and j in molecular_set):
                return sum(
                    Y[k] *
                    alpha_rule(
                        b, pobj, _get_salt(b, pname, i, k), j, b.temperature)
                    for k in b.params.anion_set)
            elif (j in b.params.cation_set and i in molecular_set):
                return sum(
                    Y[k] *
                    alpha_rule(
                        b, pobj, _get_salt(b, pname, j, k), i, b.temperature)
                    for k in b.params.anion_set)
            elif (i in b.params.anion_set and j in molecular_set):
                return sum(
                    Y[k] *
                    alpha_rule(
                        b, pobj, _get_salt(b, pname, k, i), j, b.temperature)
                    for k in b.params.cation_set)
            elif (j in b.params.anion_set and i in molecular_set):
                return sum(
                    Y[k] *
                    alpha_rule(
                        b, pobj, _get_salt(b, pname, k, j), i, b.temperature)
                    for k in b.params.cation_set)
            elif (i in b.params.cation_set and j in b.params.anion_set):
                ipair = _get_salt(b, pname, i, j)
                return sum(Y[k]*alpha_rule(b,
                                           pobj,
                                           ipair,
                                           _get_salt(b, pname, k, j),
                                           b.temperature)
                           for k in b.params.cation_set)
            elif (i in b.params.anion_set and j in b.params.cation_set):
                ipair = _get_salt(b, pname, j, i)
                return sum(Y[k]*alpha_rule(b,
                                           pobj,
                                           ipair,
                                           _get_salt(b, pname, j, k),
                                           b.temperature)
                           for k in b.params.anion_set)
            elif ((i in b.params.cation_set and j in b.params.cation_set) or
                  (i in b.params.anion_set and j in b.params.anion_set)):
                # No like ion interactions
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

            if ((i in molecular_set) and
                    (j in molecular_set)):
                return _G_appr(b, pobj, i, j, b.temperature)
            elif (i in b.params.cation_set and j in molecular_set):
                return sum(
                    Y[k] *
                    _G_appr(
                        b, pobj, _get_salt(b, pname, i, k), j,  b.temperature)
                    for k in b.params.anion_set)
            elif (j in b.params.cation_set and i in molecular_set):
                return sum(
                    Y[k] *
                    _G_appr(
                        b, pobj, _get_salt(b, pname, j, k), i, b.temperature)
                    for k in b.params.anion_set)
            elif (i in b.params.anion_set and j in molecular_set):
                return sum(
                    Y[k] *
                    _G_appr(
                        b, pobj, _get_salt(b, pname, k, i), j, b.temperature)
                    for k in b.params.cation_set)
            elif (j in b.params.anion_set and i in molecular_set):
                return sum(
                    Y[k] *
                    _G_appr(
                        b, pobj, _get_salt(b, pname, k, j), i, b.temperature)
                    for k in b.params.cation_set)
            elif (i in b.params.cation_set and j in b.params.anion_set):
                ipair = _get_salt(b, pname, i, j)
                return sum(Y[k]*_G_appr(b,
                                        pobj,
                                        ipair,
                                        _get_salt(b, pname, k, j),
                                        b.temperature)
                           for k in b.params.cation_set)
            elif (i in b.params.anion_set and j in b.params.cation_set):
                ipair = _get_salt(b, pname, j, i)
                return sum(Y[k]*_G_appr(b,
                                        pobj,
                                        ipair,
                                        _get_salt(b, pname, j, k),
                                        b.temperature)
                           for k in b.params.anion_set)
            elif ((i in b.params.cation_set and j in b.params.cation_set) or
                  (i in b.params.anion_set and j in b.params.anion_set)):
                # No like ion interactions
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
            if ((i in molecular_set) and
                    (j in molecular_set)):
                return tau_rule(b, pobj, i, j, b.temperature)
            elif ((i in b.params.cation_set and j in b.params.cation_set) or
                  (i in b.params.anion_set and j in b.params.anion_set)):
                # No like ion interactions
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


def _get_salt(b, p, c, a):
    # First, check apparent species
    for app in b.params._apparent_set:
        app_obj = b.params.get_component(app)
        if (c in app_obj.config.dissociation_species and
                a in app_obj.config.dissociation_species):
            return app
    # Next, check for weak acids and bases
    # TODO: Implement something properly for weak acids and bases
    if c == "H+" or a == "OH-":
        for r in b.params.config.inherent_reactions:
            stoic = b.params.config.inherent_reactions[r].stoichiometry
            # First, length of stoich must be == 3: 2 ions and ion-pair
            # Second, stoic mist contain both (p, c) and (p, a)
            if (len(stoic) == 3 and
                    (p, c) in stoic and
                    (p, a) in stoic):
                # Next, stoic coeff. for c and a must have same sign
                if stoic[p, c]*stoic[p, a] > 0:
                    # Need to find remaining component (the ion-pair)
                    for q, s in stoic:
                        # Ion-pair must not be a or c
                        # Stoic coeff. must have opposite sign to c and a
                        if (stoic[p, s] != stoic[p, c] and
                                stoic[p, s] != stoic[p, a] and
                                stoic[p, s]*stoic[p, c] < 0):
                            return s

    raise BurntToast("{} eNRTL error. Could not find ion pair for {}."
                     .format(b.name, (c, a)))
