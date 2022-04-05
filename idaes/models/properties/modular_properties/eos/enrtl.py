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
Methods for eNRTL activity coefficient method.

Only applicable to liquid/electrolyte phases

Many thanks to C.-C. Chen for his assistance and suggestions on testing and
verifying the model.

Reference:

Song, Y. and Chen, C.-C., Symmetric Electrolyte Nonrandom Two-Liquid Activity
Coefficient Model, Ind. Eng. Chem. Res., 2009, Vol. 48, pgs. 7788–7797

Note that "charge number" in the paper referes to the absolute value of the
ionic charge.
"""
from pyomo.environ import Expression, exp, log, Set, units as pyunits

from .ideal import Ideal
from .enrtl_reference_states import Symmetric
from .enrtl_parameters import ConstantAlpha, ConstantTau
from idaes.models.properties.modular_properties.base.utility import (
    get_method,
    get_component_object as cobj,
)
from idaes.models.properties.modular_properties.base.generic_property import StateIndex
from idaes.core.util.constants import Constants
from idaes.core.util.exceptions import BurntToast
import idaes.logger as idaeslog


# Set up logger
_log = idaeslog.getLogger(__name__)


DefaultAlphaRule = ConstantAlpha
DefaultTauRule = ConstantTau
DefaultRefState = Symmetric

# Closest appraoch parameter - implemented as a global constant for now
# This is not something the user should be changing in most cases
ClosestApproach = 14.9


class ENRTL(Ideal):
    # Add attribute indicating support for electrolyte systems
    electrolyte_support = True

    @staticmethod
    def build_parameters(b):
        # Build additional indexing sets
        pblock = b.parent_block()
        ion_pair = []
        for i in pblock.cation_set:
            for j in pblock.anion_set:
                ion_pair.append(i + ", " + j)
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
        if (
            b.config.equation_of_state_options is not None
            and "alpha_rule" in b.config.equation_of_state_options
        ):
            b.config.equation_of_state_options["alpha_rule"].build_parameters(b)
        else:
            DefaultAlphaRule.build_parameters(b)

        # Check options for tau rule
        if (
            b.config.equation_of_state_options is not None
            and "tau_rule" in b.config.equation_of_state_options
        ):
            b.config.equation_of_state_options["tau_rule"].build_parameters(b)
        else:
            DefaultTauRule.build_parameters(b)

    @staticmethod
    def common(b, pobj):
        pname = pobj.local_name

        molecular_set = b.params.solvent_set | b.params.solute_set

        # Check options for alpha rule
        if (
            pobj.config.equation_of_state_options is not None
            and "alpha_rule" in pobj.config.equation_of_state_options
        ):
            alpha_rule = pobj.config.equation_of_state_options[
                "alpha_rule"
            ].return_expression
        else:
            alpha_rule = DefaultAlphaRule.return_expression

        # Check options for tau rule
        if (
            pobj.config.equation_of_state_options is not None
            and "tau_rule" in pobj.config.equation_of_state_options
        ):
            tau_rule = pobj.config.equation_of_state_options[
                "tau_rule"
            ].return_expression
        else:
            tau_rule = DefaultTauRule.return_expression

        # Check options for reference state
        if (
            pobj.config.equation_of_state_options is not None
            and "reference_state" in pobj.config.equation_of_state_options
        ):
            ref_state = pobj.config.equation_of_state_options["reference_state"]
        else:
            ref_state = DefaultRefState

        # ---------------------------------------------------------------------
        # Calculate composition terms at both actual and reference state
        # Calculate reference state mole fractions
        ref_state.ref_state(b, pname)

        # Ionic Strength
        def rule_I(b):  # Eqn 62
            return 0.5 * sum(
                b.mole_frac_phase_comp_true[pname, c]
                * b.params.get_component(c).config.charge ** 2
                for c in b.params.ion_set
            )

        b.add_component(
            pname + "_ionic_strength", Expression(rule=rule_I, doc="Ionic strength")
        )

        def rule_I_ref(b):  # Eqn 62 evaluated at reference state
            x = getattr(b, pname + "_x_ref")
            return 0.5 * sum(
                x[c] * b.params.get_component(c).config.charge ** 2
                for c in b.params.ion_set
            )

        b.add_component(
            pname + "_ionic_strength_ref",
            Expression(rule=rule_I_ref, doc="Ionic strength at reference state"),
        )

        # Calculate mixing factors
        def rule_X(b, j):  # Eqn 21
            if (pname, j) not in b.params.true_phase_component_set:
                return Expression.Skip
            elif j in b.params.cation_set or j in b.params.anion_set:
                return b.mole_frac_phase_comp_true[pname, j] * abs(
                    cobj(b, j).config.charge
                )
            else:
                return b.mole_frac_phase_comp_true[pname, j]

        b.add_component(
            pname + "_X",
            Expression(
                b.params.true_species_set,
                rule=rule_X,
                doc="Charge x mole fraction term",
            ),
        )

        def rule_X_ref(b, j):  # Eqn 21 evaluated at reference state
            x = getattr(b, pname + "_x_ref")
            if (pname, j) not in b.params.true_phase_component_set:
                return Expression.Skip
            elif j in b.params.cation_set or j in b.params.anion_set:
                return x[j] * abs(cobj(b, j).config.charge)
            else:
                return x[j]

        b.add_component(
            pname + "_X_ref",
            Expression(
                b.params.true_species_set,
                rule=rule_X_ref,
                doc="Charge x mole fraction term at reference state",
            ),
        )

        def rule_Y(b, j):
            if cobj(b, j).config.charge < 0:
                # Anion
                dom = b.params.anion_set
            else:
                dom = b.params.cation_set

            X = getattr(b, pname + "_X")
            return X[j] / sum(X[i] for i in dom)  # Eqns 36 and 37

        # Y is a charge ratio, and thus independent of x for symmetric state
        # TODO: This may need to change for the unsymmetric state
        b.add_component(
            pname + "_Y",
            Expression(b.params.ion_set, rule=rule_Y, doc="Charge composition"),
        )

        # ---------------------------------------------------------------------
        # Long-range terms
        # Average molar volume of solvent
        def rule_vol_mol_solvent(b):  # Eqn 77
            if len(b.params.solvent_set) == 1:
                s = b.params.solvent_set.first()
                return ENRTL.get_vol_mol_pure(b, "liq", s, b.temperature)
            else:
                return sum(
                    b.mole_frac_phase_comp_true[pname, s]
                    * ENRTL.get_vol_mol_pure(b, "liq", s, b.temperature)
                    for s in b.params.solvent_set
                ) / sum(
                    b.mole_frac_phase_comp_true[pname, s] for s in b.params.solvent_set
                )

        b.add_component(
            pname + "_vol_mol_solvent",
            Expression(rule=rule_vol_mol_solvent, doc="Mean molar volume of solvent"),
        )

        # Mean relative permitivity of solvent
        def rule_eps_solvent(b):  # Eqn 78
            if len(b.params.solvent_set) == 1:
                s = b.params.solvent_set.first()
                return get_method(b, "relative_permittivity_liq_comp", s)(
                    b, cobj(b, s), b.temperature
                )
            else:
                return sum(
                    b.mole_frac_phase_comp_true[pname, s]
                    * get_method(b, "relative_permittivity_liq_comp", s)(
                        b, cobj(b, s), b.temperature
                    )
                    * b.params.get_component(s).mw
                    for s in b.params.solvent_set
                ) / sum(
                    b.mole_frac_phase_comp_true[pname, s] * b.params.get_component(s).mw
                    for s in b.params.solvent_set
                )

        b.add_component(
            pname + "_relative_permittivity_solvent",
            Expression(
                rule=rule_eps_solvent, doc="Mean relative permittivity  of solvent"
            ),
        )

        # Debye-Huckel parameter
        def rule_A_DH(b):  # Eqn 61
            # Note: Where the paper refers to the dielectric constant, it
            # actually means the electric permittivity of the solvent
            # eps = eps_r*eps_0 (units F/m)
            # Note that paper is missing a required 4*pi term
            v = pyunits.convert(
                getattr(b, pname + "_vol_mol_solvent"), pyunits.m**3 / pyunits.mol
            )
            eps = getattr(b, pname + "_relative_permittivity_solvent")
            eps0 = Constants.vacuum_electric_permittivity
            return (
                (1 / 3)
                * (2 * Constants.pi * Constants.avogadro_number / v) ** 0.5
                * (
                    Constants.elemental_charge**2
                    / (
                        4
                        * Constants.pi
                        * eps
                        * eps0
                        * Constants.boltzmann_constant
                        * b.temperature
                    )
                )
                ** (3 / 2)
            )

        b.add_component(
            pname + "_A_DH", Expression(rule=rule_A_DH, doc="Debye-Huckel parameter")
        )

        # Long-range (PDH) contribution to activity coefficient
        def rule_log_gamma_pdh(b, j):
            A = getattr(b, pname + "_A_DH")
            Ix = getattr(b, pname + "_ionic_strength")
            I0 = getattr(b, pname + "_ionic_strength_ref")
            rho = ClosestApproach
            if j in molecular_set:
                # Eqn 69
                # Note typo in original paper. Correct power for I is (3/2)
                return 2 * A * Ix ** (3 / 2) / (1 + rho * Ix ** (1 / 2))
            elif j in b.params.ion_set:
                # Eqn 70
                z = abs(cobj(b, j).config.charge)
                return -A * (
                    (2 * z**2 / rho)
                    * log((1 + rho * Ix**0.5) / (1 + rho * I0**0.5))
                    + (z**2 * Ix**0.5 - 2 * Ix ** (3 / 2)) / (1 + rho * Ix**0.5)
                    - (2 * Ix * I0**-0.5)
                    / (1 + rho * I0**0.5)
                    * ref_state.ndIdn(b, pname, j)
                )
            else:
                raise BurntToast(
                    "{} eNRTL model encountered unexpected component.".format(b.name)
                )

        b.add_component(
            pname + "_log_gamma_pdh",
            Expression(
                b.params.true_species_set,
                rule=rule_log_gamma_pdh,
                doc="Long-range contribution to activity coefficient",
            ),
        )

        # ---------------------------------------------------------------------
        # Local Contribution Terms
        # For the symmetric state, all of these are independent of composition
        # TODO: For the unsymmetric state, it may be necessary to recalculate
        # Calculate alphas for all true species pairings
        def rule_alpha_expr(b, i, j):
            Y = getattr(b, pname + "_Y")
            if (pname, i) not in b.params.true_phase_component_set or (
                pname,
                j,
            ) not in b.params.true_phase_component_set:
                return Expression.Skip
            elif (i in molecular_set) and (j in molecular_set):
                # alpha equal user provided parameters
                return alpha_rule(b, pobj, i, j, b.temperature)
            elif i in b.params.cation_set and j in molecular_set:
                # Eqn 32
                return sum(
                    Y[k] * alpha_rule(b, pobj, (i + ", " + k), j, b.temperature)
                    for k in b.params.anion_set
                )
            elif j in b.params.cation_set and i in molecular_set:
                # Eqn 32
                return sum(
                    Y[k] * alpha_rule(b, pobj, (j + ", " + k), i, b.temperature)
                    for k in b.params.anion_set
                )
            elif i in b.params.anion_set and j in molecular_set:
                # Eqn 33
                return sum(
                    Y[k] * alpha_rule(b, pobj, (k + ", " + i), j, b.temperature)
                    for k in b.params.cation_set
                )
            elif j in b.params.anion_set and i in molecular_set:
                # Eqn 33
                return sum(
                    Y[k] * alpha_rule(b, pobj, (k + ", " + j), i, b.temperature)
                    for k in b.params.cation_set
                )
            elif i in b.params.cation_set and j in b.params.anion_set:
                # Eqn 34
                if len(b.params.cation_set) > 1:
                    return sum(
                        Y[k]
                        * alpha_rule(
                            b, pobj, (i + ", " + j), (k + ", " + j), b.temperature
                        )
                        for k in b.params.cation_set
                    )
                else:
                    return 0.2
            elif i in b.params.anion_set and j in b.params.cation_set:
                # Eqn 35
                if len(b.params.anion_set) > 1:
                    return sum(
                        Y[k]
                        * alpha_rule(
                            b, pobj, (j + ", " + i), (j + ", " + k), b.temperature
                        )
                        for k in b.params.anion_set
                    )
                else:
                    return 0.2
            elif (i in b.params.cation_set and j in b.params.cation_set) or (
                i in b.params.anion_set and j in b.params.anion_set
            ):
                # No like-ion interactions
                return Expression.Skip
            else:
                raise BurntToast(
                    "{} eNRTL model encountered unexpected component pair {}.".format(
                        b.name, (i, j)
                    )
                )

        b.add_component(
            pname + "_alpha",
            Expression(
                b.params.true_species_set,
                b.params.true_species_set,
                rule=rule_alpha_expr,
                doc="Non-randomness parameters",
            ),
        )

        # Calculate G terms
        def rule_G_expr(b, i, j):
            Y = getattr(b, pname + "_Y")

            def _G_appr(b, pobj, i, j, T):  # Eqn 23
                if i != j:
                    return exp(
                        -alpha_rule(b, pobj, i, j, T) * tau_rule(b, pobj, i, j, T)
                    )
                else:
                    return 1

            if (pname, i) not in b.params.true_phase_component_set or (
                pname,
                j,
            ) not in b.params.true_phase_component_set:
                return Expression.Skip
            elif (i in molecular_set) and (j in molecular_set):
                # G comes directly from parameters
                return _G_appr(b, pobj, i, j, b.temperature)
            elif i in b.params.cation_set and j in molecular_set:
                # Eqn 38
                return sum(
                    Y[k] * _G_appr(b, pobj, (i + ", " + k), j, b.temperature)
                    for k in b.params.anion_set
                )
            elif i in molecular_set and j in b.params.cation_set:
                # Eqn 40
                return sum(
                    Y[k] * _G_appr(b, pobj, i, (j + ", " + k), b.temperature)
                    for k in b.params.anion_set
                )
            elif i in b.params.anion_set and j in molecular_set:
                # Eqn 39
                return sum(
                    Y[k] * _G_appr(b, pobj, (k + ", " + i), j, b.temperature)
                    for k in b.params.cation_set
                )
            elif i in molecular_set and j in b.params.anion_set:
                # Eqn 41
                return sum(
                    Y[k] * _G_appr(b, pobj, i, (k + ", " + j), b.temperature)
                    for k in b.params.cation_set
                )
            elif i in b.params.cation_set and j in b.params.anion_set:
                # Eqn 42
                if len(b.params.cation_set) > 1:
                    return sum(
                        Y[k]
                        * _G_appr(
                            b, pobj, (i + ", " + j), (k + ", " + j), b.temperature
                        )
                        for k in b.params.cation_set
                    )
                else:
                    # This term does not exist for single cation systems
                    # However, need a valid result to calculate tau
                    return 1
            elif i in b.params.anion_set and j in b.params.cation_set:
                # Eqn 43
                if len(b.params.anion_set) > 1:
                    return sum(
                        Y[k]
                        * _G_appr(
                            b, pobj, (j + ", " + i), (j + ", " + k), b.temperature
                        )
                        for k in b.params.anion_set
                    )
                else:
                    # This term does not exist for single anion systems
                    # However, need a valid result to calculate tau
                    return 1
            elif (i in b.params.cation_set and j in b.params.cation_set) or (
                i in b.params.anion_set and j in b.params.anion_set
            ):
                # No like-ion interactions
                return Expression.Skip
            else:
                raise BurntToast(
                    "{} eNRTL model encountered unexpected component pair {}.".format(
                        b.name, (i, j)
                    )
                )

        b.add_component(
            pname + "_G",
            Expression(
                b.params.true_species_set,
                b.params.true_species_set,
                rule=rule_G_expr,
                doc="Local interaction G term",
            ),
        )

        # Calculate tau terms
        def rule_tau_expr(b, i, j):
            if (pname, i) not in b.params.true_phase_component_set or (
                pname,
                j,
            ) not in b.params.true_phase_component_set:
                return Expression.Skip
            elif (i in molecular_set) and (j in molecular_set):
                # tau equal to parameter
                return tau_rule(b, pobj, i, j, b.temperature)
            elif (i in b.params.cation_set and j in b.params.cation_set) or (
                i in b.params.anion_set and j in b.params.anion_set
            ):
                # No like-ion interactions
                return Expression.Skip
            else:
                alpha = getattr(b, pname + "_alpha")
                G = getattr(b, pname + "_G")
                # Eqn 44
                return -log(G[i, j]) / alpha[i, j]

        b.add_component(
            pname + "_tau",
            Expression(
                b.params.true_species_set,
                b.params.true_species_set,
                rule=rule_tau_expr,
                doc="Binary interaction energy parameters",
            ),
        )

        # Local contribution to activity coefficient
        def rule_log_gamma_lc_I(b, s):
            X = getattr(b, pname + "_X")
            G = getattr(b, pname + "_G")
            tau = getattr(b, pname + "_tau")

            return log_gamma_lc(b, pname, s, X, G, tau)

        b.add_component(
            pname + "_log_gamma_lc_I",
            Expression(
                b.params.true_species_set,
                rule=rule_log_gamma_lc_I,
                doc="Local contribution at actual state",
            ),
        )

        def rule_log_gamma_lc_I0(b, s):
            X = getattr(b, pname + "_X_ref")
            G = getattr(b, pname + "_G")
            tau = getattr(b, pname + "_tau")

            return log_gamma_lc(b, pname, s, X, G, tau)

        b.add_component(
            pname + "_log_gamma_lc_I0",
            Expression(
                b.params.ion_set,
                rule=rule_log_gamma_lc_I0,
                doc="Local contribution at reference state",
            ),
        )

        def rule_log_gamma_lc(b, s):
            log_gamma_lc_I = getattr(b, pname + "_log_gamma_lc_I")
            if s in molecular_set:
                return log_gamma_lc_I[s]
            else:
                log_gamma_lc_I0 = getattr(b, pname + "_log_gamma_lc_I0")
                return log_gamma_lc_I[s] - log_gamma_lc_I0[s]

        b.add_component(
            pname + "_log_gamma_lc",
            Expression(
                b.params.true_species_set,
                rule=rule_log_gamma_lc,
                doc="Local contribution contribution to activity coefficient",
            ),
        )

        # Overall log gamma
        def rule_log_gamma(b, j):
            pdh = getattr(b, pname + "_log_gamma_pdh")
            lc = getattr(b, pname + "_log_gamma_lc")
            return pdh[j] + lc[j]

        b.add_component(
            pname + "_log_gamma",
            Expression(
                b.params.true_species_set,
                rule=rule_log_gamma,
                doc="Log of activity coefficient",
            ),
        )

        # Activity coefficient of apparent species
        def rule_log_gamma_pm(b, j):
            cobj = b.params.get_component(j)

            if "dissociation_species" in cobj.config:
                dspec = cobj.config.dissociation_species

                n = 0
                d = 0
                for s in dspec:
                    dobj = b.params.get_component(s)
                    ln_g = getattr(b, pname + "_log_gamma")[s]
                    n += abs(dobj.config.charge) * ln_g
                    d += abs(dobj.config.charge)

                return n / d
            else:
                return getattr(b, pname + "_log_gamma")[j]

        b.add_component(
            pname + "_log_gamma_appr",
            Expression(
                b.params.apparent_species_set,
                rule=rule_log_gamma_pm,
                doc="Log of mean activity coefficient",
            ),
        )

    @staticmethod
    def calculate_scaling_factors(b, pobj):
        pass

    @staticmethod
    def act_phase_comp(b, p, j):
        return b.mole_frac_phase_comp[p, j] * b.act_coeff_phase_comp[p, j]

    @staticmethod
    def act_phase_comp_true(b, p, j):
        ln_gamma = getattr(b, p + "_log_gamma")
        return b.mole_frac_phase_comp_true[p, j] * exp(ln_gamma[j])

    @staticmethod
    def act_phase_comp_appr(b, p, j):
        ln_gamma = getattr(b, p + "_log_gamma_appr")
        return b.mole_frac_phase_comp_apparent[p, j] * exp(ln_gamma[j])

    @staticmethod
    def act_coeff_phase_comp(b, p, j):
        if b.params.config.state_components == StateIndex.true:
            ln_gamma = getattr(b, p + "_log_gamma")
        else:
            ln_gamma = getattr(b, p + "_log_gamma_appr")
        return exp(ln_gamma[j])

    @staticmethod
    def act_coeff_phase_comp_true(b, p, j):
        ln_gamma = getattr(b, p + "_log_gamma")
        return exp(ln_gamma[j])

    @staticmethod
    def act_coeff_phase_comp_appr(b, p, j):
        ln_gamma = getattr(b, p + "_log_gamma_appr")
        return exp(ln_gamma[j])

    @staticmethod
    def pressure_osm_phase(b, p):
        return (
            -ENRTL.gas_constant(b)
            * b.temperature
            * b.log_act_phase_solvents[p]
            / b.vol_mol_phase[p]
        )

    @staticmethod
    def vol_mol_phase(b, p):
        # eNRTL model uses apparent species for calculating molar volume
        # TODO : Need something more rigorus to handle concentrated solutions
        v_expr = 0
        for j in b.params.apparent_species_set:
            v_comp = ENRTL.get_vol_mol_pure(b, "liq", j, b.temperature)
            v_expr += b.mole_frac_phase_comp_apparent[p, j] * v_comp

        return v_expr


def log_gamma_lc(b, pname, s, X, G, tau):
    # General function for calculating local contributions
    # The same method can be used for both actual state and reference state
    # by providing different X, G and tau expressions.

    # Indicies in expressions use same names as source paper
    # mp = m'
    molecular_set = b.params.solvent_set | b.params.solute_set
    aqu_species = b.params.true_species_set - b.params._non_aqueous_set

    if (pname, s) not in b.params.true_phase_component_set:
        # Non-aqueous component
        return Expression.Skip
    if s in b.params.cation_set:
        c = s
        Z = b.params.get_component(c).config.charge

        # Eqn 26
        return Z * (
            sum(
                (X[m] * G[c, m] / sum(X[i] * G[i, m] for i in aqu_species))
                * (
                    tau[c, m]
                    - (
                        sum(X[i] * G[i, m] * tau[i, m] for i in aqu_species)
                        / sum(X[i] * G[i, m] for i in aqu_species)
                    )
                )
                for m in molecular_set
            )
            + sum(
                X[i] * G[i, c] * tau[i, c] for i in (aqu_species - b.params.cation_set)
            )
            / sum(X[i] * G[i, c] for i in (aqu_species - b.params.cation_set))
            + sum(
                (
                    X[a]
                    * G[c, a]
                    / sum(X[i] * G[i, a] for i in (aqu_species - b.params.anion_set))
                )
                * (
                    tau[c, a]
                    - sum(
                        X[i] * G[i, a] * tau[i, a]
                        for i in (aqu_species - b.params.anion_set)
                    )
                    / sum(X[i] * G[i, a] for i in (aqu_species - b.params.anion_set))
                )
                for a in b.params.anion_set
            )
        )
    elif s in b.params.anion_set:
        a = s
        Z = abs(b.params.get_component(a).config.charge)

        # Eqn 27
        return Z * (
            sum(
                (X[m] * G[a, m] / sum(X[i] * G[i, m] for i in aqu_species))
                * (
                    tau[a, m]
                    - (
                        sum(X[i] * G[i, m] * tau[i, m] for i in aqu_species)
                        / sum(X[i] * G[i, m] for i in aqu_species)
                    )
                )
                for m in molecular_set
            )
            + sum(
                X[i] * G[i, a] * tau[i, a] for i in (aqu_species - b.params.anion_set)
            )
            / sum(X[i] * G[i, a] for i in (aqu_species - b.params.anion_set))
            + sum(
                (
                    X[c]
                    * G[a, c]
                    / sum(X[i] * G[i, c] for i in (aqu_species - b.params.cation_set))
                )
                * (
                    tau[a, c]
                    - sum(
                        X[i] * G[i, c] * tau[i, c]
                        for i in (aqu_species - b.params.cation_set)
                    )
                    / sum(X[i] * G[i, c] for i in (aqu_species - b.params.cation_set))
                )
                for c in b.params.cation_set
            )
        )
    else:
        m = s
        # Eqn 25
        return (
            sum(X[i] * G[i, m] * tau[i, m] for i in aqu_species)
            / sum(X[i] * G[i, m] for i in aqu_species)
            + sum(
                (X[mp] * G[m, mp] / sum(X[i] * G[i, mp] for i in aqu_species))
                * (
                    tau[m, mp]
                    - (
                        sum(X[i] * G[i, mp] * tau[i, mp] for i in aqu_species)
                        / sum(X[i] * G[i, mp] for i in aqu_species)
                    )
                )
                for mp in molecular_set
            )
            + sum(
                (
                    X[c]
                    * G[m, c]
                    / sum(X[i] * G[i, c] for i in (aqu_species - b.params.cation_set))
                )
                * (
                    tau[m, c]
                    - (
                        sum(
                            X[i] * G[i, c] * tau[i, c]
                            for i in (aqu_species - b.params.cation_set)
                        )
                        / sum(
                            X[i] * G[i, c] for i in (aqu_species - b.params.cation_set)
                        )
                    )
                )
                for c in b.params.cation_set
            )
            + sum(
                (
                    X[a]
                    * G[m, a]
                    / sum(X[i] * G[i, a] for i in (aqu_species - b.params.anion_set))
                )
                * (
                    tau[m, a]
                    - (
                        sum(
                            X[i] * G[i, a] * tau[i, a]
                            for i in (aqu_species - b.params.anion_set)
                        )
                        / sum(
                            X[i] * G[i, a] for i in (aqu_species - b.params.anion_set)
                        )
                    )
                )
                for a in b.params.anion_set
            )
        )
