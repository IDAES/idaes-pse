#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Class implementing a variation of the enhancement factor model presented in Gaspar and Fosbol
(2015, https://doi.org/10.1016/j.ces.2015.08.023). An expression giving the enhancement factor
at equilibrium is used instead of solving the full system of equations. A paper containing a
full explanation of this model has been submitted to CAChe
"""

from pyomo.environ import (
    value,
    Var,
    Constraint,
    Expression,
    exp,
    log,
    units as pyunits,
)
from pyomo.util.calc_var_value import calculate_variable_from_constraint

import idaes.logger as idaeslog

__author__ = "Douglas Allan"


class PseudoSecondOrderExplicit(object):
    """
    Class implementing a variation of the enhancement factor model presented in Gaspar and Fosbol
    (2015, https://doi.org/10.1016/j.ces.2015.08.023). An expression giving the enhancement factor
    at equilibrium is used instead of solving the full system of equations. A paper containing a
    full explanation of this model has been submitted to CAChe
    """

    @staticmethod
    def make_model(blk, kinetics="Putta"):
        """
        Enhancement factor based liquid phase mass transfer model.

        Arguments:
            blk: MEASolventColumn object whose enhancement factor model needs to be created
            kinetics : str giving which set of kinetics to use. Valid options are "Putta"
                and "Luo"
        """
        lunits = blk.liquid_phase.properties.params.get_metadata().get_derived_units

        blk.log_rate_constant_MEA = Var(
            blk.flowsheet().time,
            blk.liquid_phase.length_domain,
            doc="Logarithm of rate constant for MEA mechanism",
            initialize=0,
        )

        @blk.Constraint(blk.flowsheet().time, blk.liquid_phase.length_domain)
        def log_rate_constant_MEA_eqn(b, t, x):
            if x == b.liquid_phase.length_domain.last():
                return Constraint.Skip
            else:
                if kinetics == "Putta":
                    # Putta, Svendsen, Knuutila 2017 Eqn. 42
                    reduced_activation_energy_MEA = pyunits.convert(
                        4936.6 * pyunits.K, to_units=lunits("temperature")
                    )
                    preexponential_factor_MEA = pyunits.convert(
                        3.1732e3 * ((pyunits.m) ** 6 / (pyunits.mol**2 * pyunits.s)),
                        to_units=1 / (lunits("time") * lunits("density_mole") ** 2),
                    )
                elif kinetics == "Luo":
                    reduced_activation_energy_MEA = pyunits.convert(
                        4742.0 * pyunits.K, to_units=lunits("temperature")
                    )
                    preexponential_factor_MEA = pyunits.convert(
                        2.003e4 * ((pyunits.m) ** 6 / (pyunits.mol**2 * pyunits.s)),
                        to_units=1 / (lunits("time") * lunits("density_mole") ** 2),
                    )
                else:
                    return ValueError(
                        "The kinetics option can take on values of 'Luo' and 'Putta', but "
                        f"an unknown option {kinetics} was passed instead."
                    )

                log_preexponential_factor_MEA = log(value(preexponential_factor_MEA))

                return b.log_rate_constant_MEA[t, x] == (
                    log_preexponential_factor_MEA
                    + b.liquid_phase.properties[t, x].log_conc_mol_phase_comp_true[
                        "Liq", "MEA"
                    ]
                    - reduced_activation_energy_MEA
                    / b.liquid_phase.properties[t, x].temperature
                )

        blk.log_rate_constant_H2O = Var(
            blk.flowsheet().time,
            blk.liquid_phase.length_domain,
            doc="Logarithm of rate constant for H2O mechanism",
            initialize=0,
        )

        @blk.Constraint(blk.flowsheet().time, blk.liquid_phase.length_domain)
        def log_rate_constant_H2O_eqn(b, t, x):
            if x == b.liquid_phase.length_domain.last():
                return Constraint.Skip
            else:
                if kinetics == "Putta":
                    # Putta, Svendsen, Knuutila 2017 Eqn. 42
                    reduced_activation_energy_H2O = pyunits.convert(
                        3900 * pyunits.K, to_units=lunits("temperature")
                    )
                    preexponential_factor_H2O = pyunits.convert(
                        1.0882e2 * ((pyunits.m) ** 6 / (pyunits.mol**2 * pyunits.s)),
                        to_units=1 / (lunits("time") * lunits("density_mole") ** 2),
                    )

                elif kinetics == "Luo":
                    reduced_activation_energy_H2O = pyunits.convert(
                        3110 * pyunits.K, to_units=lunits("temperature")
                    )
                    preexponential_factor_H2O = pyunits.convert(
                        4.147 * ((pyunits.m) ** 6 / (pyunits.mol**2 * pyunits.s)),
                        to_units=1 / (lunits("time") * lunits("density_mole") ** 2),
                    )
                else:
                    return ValueError(
                        "The kinetics option can take on values of 'Luo' and 'Putta', but "
                        f"an unknown option {kinetics} was passed instead."
                    )

                log_preexponential_factor_H2O = log(value(preexponential_factor_H2O))
                return b.log_rate_constant_H2O[t, x] == (
                    log_preexponential_factor_H2O
                    + b.liquid_phase.properties[t, x].log_conc_mol_phase_comp_true[
                        "Liq", "H2O"
                    ]
                    - reduced_activation_energy_H2O
                    / b.liquid_phase.properties[t, x].temperature
                )

        @blk.Expression(
            blk.flowsheet().time,
            blk.liquid_phase.length_domain,
            doc="Second order rate constant [m3/(mol.s)]",
        )
        def log_rate_constant(b, t, x):
            if x == b.liquid_phase.length_domain.last():
                return Expression.Skip
            else:
                return log(
                    exp(b.log_rate_constant_MEA[t, x])
                    + exp(b.log_rate_constant_H2O[t, x])
                )

        blk.log_hatta_number = Var(
            blk.flowsheet().time,
            blk.liquid_phase.length_domain,
            initialize=3,
            doc="Hatta number",
        )

        @blk.Constraint(
            blk.flowsheet().time,
            blk.liquid_phase.length_domain,
            doc="Hatta number constraint",
        )
        def hatta_number_eqn(b, t, x):
            if x == b.liquid_phase.length_domain.last():
                return Constraint.Skip
            else:
                return (
                    b.log_hatta_number[t, x]
                    == 0.5
                    * (
                        b.log_rate_constant[t, x]
                        + b.liquid_phase.properties[t, x].log_conc_mol_phase_comp_true[
                            "Liq", "MEA"
                        ]
                        + b.log_diffus_liq_comp[t, x, "CO2"]
                    )
                    - b.log_mass_transfer_coeff_liq[t, x, "CO2"]
                )

        @blk.Expression(
            blk.flowsheet().time,
            blk.liquid_phase.length_domain,
            ["MEA_+", "MEACOO_-", "CO2"],
            doc="Quotient of diffusion rates",
        )
        def diffus_ratio(b, t, x, j):
            # When evaluated for CO2, we get instant_E_hat_minus_one
            if x == b.liquid_phase.length_domain.last():
                return Expression.Skip
            props = b.liquid_phase.properties[t, x]
            return (
                props.diffus_phase_comp_true["Liq", "MEA"]
                * props.conc_mol_phase_comp_true["Liq", "MEA"]
                / (
                    2
                    * props.diffus_phase_comp_true["Liq", j]
                    * props.conc_mol_phase_comp_true["Liq", j]
                )
            )

        @blk.Expression(
            blk.flowsheet().time,
            blk.liquid_phase.length_domain,
        )
        def E_hat(b, t, x):
            if x == b.liquid_phase.length_domain.last():
                return Expression.Skip
            Ha = exp(b.log_hatta_number[t, x])
            R_plus = b.diffus_ratio[t, x, "MEA_+"]
            R_minus = b.diffus_ratio[t, x, "MEACOO_-"]
            Estar_hat_minus_one = b.diffus_ratio[t, x, "CO2"]
            return 1 + (Ha - 1) / (
                Ha * (R_plus + R_minus + 2) / Estar_hat_minus_one + 1
            )

        @blk.Constraint(
            blk.flowsheet().time,
            blk.liquid_phase.length_domain,
        )
        def enhancement_factor_eqn(b, t, x):
            if x == b.liquid_phase.length_domain.last():
                return Constraint.Skip
            return b.log_enhancement_factor[t, x] == log(b.E_hat[t, x])

        enhancement_factor_vars = [
            blk.log_enhancement_factor,
            blk.log_hatta_number,
            blk.log_rate_constant_MEA,
            blk.log_rate_constant_H2O,
        ]

        enhancement_factor_constraints = [
            blk.enhancement_factor_eqn,
            blk.hatta_number_eqn,
            blk.log_rate_constant_MEA_eqn,
            blk.log_rate_constant_H2O_eqn,
        ]

        return enhancement_factor_vars, enhancement_factor_constraints

    @staticmethod
    def initialize_model(
        blk,
        outlvl=idaeslog.NOTSET,
        optarg=None,
        solver=None,
    ):
        """
        Initialization method for enhancement factor submodel for MEASolventColumn

        Arguments:
            blk: MEASolventColumn object whose enhancement factor model needs to be initialized
            optarg : solver options dictionary object (default=None, use
                        default solver options)
            solver : str indicating which solver to use during initialization
                    (default = None, use IDAES default solver)
            outlvl : output level for logging
        """
        # Set solver options
        if optarg is None:
            optarg = {}

        for t in blk.flowsheet().time:
            for x in blk.liquid_phase.length_domain:
                if x == blk.liquid_phase.length_domain.last():
                    continue

                calculate_variable_from_constraint(
                    blk.log_rate_constant_MEA[t, x], blk.log_rate_constant_MEA_eqn[t, x]
                )
                calculate_variable_from_constraint(
                    blk.log_rate_constant_H2O[t, x], blk.log_rate_constant_H2O_eqn[t, x]
                )
                calculate_variable_from_constraint(
                    blk.log_hatta_number[t, x], blk.hatta_number_eqn[t, x]
                )
                calculate_variable_from_constraint(
                    blk.log_enhancement_factor[t, x], blk.enhancement_factor_eqn[t, x]
                )
