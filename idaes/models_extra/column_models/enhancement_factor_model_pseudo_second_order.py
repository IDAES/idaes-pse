# Import Python libraries
import logging
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Import Pyomo libraries
import pyomo.opt
from pyomo.environ import (
    Block,
    ConcreteModel,
    value,
    Var,
    Reals,
    NonNegativeReals,
    Param,
    TransformationFactory,
    Constraint,
    Expression,
    Objective,
    SolverStatus,
    TerminationCondition,
    check_optimal_termination,
    assert_optimal_termination,
    exp,
    log,
    sqrt,
    units as pyunits,
    Set,
    Reference
)
from pyomo.common.collections import ComponentSet, ComponentMap

from pyomo.util.calc_var_value import calculate_variable_from_constraint
from pyomo.common.config import ConfigValue, Bool

# Import IDAES Libraries
from idaes.core.util.constants import Constants as CONST
from idaes.models_extra.column_models.solvent_column import PackedColumnData

from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import _fix_vars, _restore_fixedness
from idaes.core import declare_process_block_class, FlowsheetBlock, StateBlock
from idaes.core.util.exceptions import InitializationError
from idaes.core.solvers.get_solver import get_solver
import idaes.logger as idaeslog

from idaes.core.solvers import use_idaes_solver_configuration_defaults
import idaes.core.util.scaling as iscale
from pyomo.util.subsystems import (
    create_subsystem_block,
)
from idaes.core.solvers.petsc import (
    _sub_problem_scaling_suffix,
)

def make_enhancement_factor_model(blk, lunits, kinetics="Putta"):
    """
    Enhancement factor based liquid phase mass transfer model.
    """
    assert kinetics in {"Luo", "Putta"}
    blk.log_rate_constant_MEA = Var(
        blk.flowsheet().time,
        blk.liquid_phase.length_domain,
        doc="Logarithm of rate constant for MEA mechanism",
        initialize = 0,
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
                    3.1732e3 * ((pyunits.m) ** 6 / (pyunits.mol ** 2 * pyunits.s)),
                    to_units=1 / (lunits("time") * lunits("density_mole")**2),
                )
            elif kinetics == "Luo":
                reduced_activation_energy_MEA = pyunits.convert(
                    4742.0 * pyunits.K, to_units=lunits("temperature")
                )
                preexponential_factor_MEA = pyunits.convert(
                    2.003e4 * ((pyunits.m) ** 6 / (pyunits.mol ** 2 * pyunits.s)),
                    to_units=1 / (lunits("time") * lunits("density_mole")**2),
                )
            else:
                return AssertionError
                
            log_preexponential_factor_MEA = log(value(preexponential_factor_MEA))

            return b.log_rate_constant_MEA[t, x] == (
                log_preexponential_factor_MEA 
                + b.liquid_phase.properties[t, x].log_conc_mol_phase_comp_true["Liq", "MEA"]
                - reduced_activation_energy_MEA / b.liquid_phase.properties[t, x].temperature
            )
    

    blk.log_rate_constant_H2O = Var(
        blk.flowsheet().time,
        blk.liquid_phase.length_domain,
        doc="Logarithm of rate constant for H2O mechanism",
        initialize = 0,
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
                    1.0882e2 * ((pyunits.m) ** 6 / (pyunits.mol ** 2 * pyunits.s)),
                    to_units=1 / (lunits("time") * lunits("density_mole")**2),
                )
                
            elif kinetics == "Luo":
                reduced_activation_energy_H2O = pyunits.convert(
                    3110 * pyunits.K, to_units=lunits("temperature")
                )
                preexponential_factor_H2O = pyunits.convert(
                    4.147 * ((pyunits.m) ** 6 / (pyunits.mol ** 2 * pyunits.s)),
                    to_units=1 / (lunits("time") * lunits("density_mole")**2),
                )
            else:
                return AssertionError

            log_preexponential_factor_H2O = log(value(preexponential_factor_H2O))
            return b.log_rate_constant_H2O[t, x] == (
                log_preexponential_factor_H2O 
                + b.liquid_phase.properties[t, x].log_conc_mol_phase_comp_true["Liq", "H2O"]
                - reduced_activation_energy_H2O / b.liquid_phase.properties[t, x].temperature
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
            return b.log_hatta_number[t, x] == 0.5 * (
                    b.log_rate_constant[t, x]
                    + b.liquid_phase.properties[t, x].log_conc_mol_phase_comp_true[
                        "Liq", "MEA"
                    ]
                    + b.log_diffus_liq_comp[t, x, "CO2"]
            ) - b.log_mass_transfer_coeff_liq[t, x, "CO2"]

    blk.conc_CO2_bulk = Var(
        blk.flowsheet().time,
        blk.liquid_phase.length_domain,
        initialize=0.95,
        units=pyunits.dimensionless,
        bounds=(0, None),
        doc="""Dimensionless concentration of CO2,
                                     Driving force term where
                                     Absorption implies conc_CO2_bulk < 1 and 
                                     Desorption implies conc_CO2_bulk > 1 """,
    )
    blk.log_conc_CO2_bulk = Var(
        blk.flowsheet().time,
        blk.liquid_phase.length_domain,
        initialize=-0.25,
        units=pyunits.dimensionless,
        doc="Natural logarithm of conc_CO2_bulk",
    )

    @blk.Constraint(
        blk.flowsheet().time,
        blk.liquid_phase.length_domain,
        doc="Equation defining log_conc_CO2_bulk",
    )
    def log_conc_CO2_bulk_eqn(b, t, x):
        if x == b.liquid_phase.length_domain.last():
            return Constraint.Skip
        else:
            return exp(b.log_conc_CO2_bulk[t, x]) == b.conc_CO2_bulk[t, x]

    @blk.Constraint(
        blk.flowsheet().time,
        blk.liquid_phase.length_domain,
        doc="""Dimensionless concentration of CO2""",
    )
    def conc_CO2_bulk_eqn(b, t, x):
        if x == b.liquid_phase.length_domain.last():
            return Constraint.Skip
        else:
            zf = b.liquid_phase.length_domain.next(x)
            Pressure = pyunits.convert(
                b.vapor_phase.properties[t, zf].pressure,
                to_units=lunits("pressure"),
            )
            return b.log_conc_CO2_bulk[t, x] + log(
                (
                    b.vapor_phase.properties[t, zf].mole_frac_comp["CO2"]
                    * Pressure
                    / b.psi[t, zf]
                    + exp(b.liquid_phase.properties[t, x].log_conc_mol_phase_comp_true[
                        "Liq", "CO2"
                    ])
                ) / lunits("density_mole")
            ) == b.liquid_phase.properties[t, x].log_conc_mol_phase_comp_true[
                       "Liq", "CO2"
                   ] + log(
                b.liquid_phase.properties[t, x].henry["Liq", "CO2"]
                / b.psi[t, zf]
                + 1
            )
            # return (
            #     b.conc_CO2_bulk[t, x] * (
            #         b.vapor_phase.properties[t, zf].mole_frac_comp["CO2"]
            #         * Pressure
            #         / b.psi[t, zf]
            #         + b.liquid_phase.properties[t, x].conc_mol_phase_comp_true["Liq", "CO2"]
            #     )
            #     == b.liquid_phase.properties[t, x].conc_mol_phase_comp_true["Liq", "CO2"]
            #     * (
            #         b.liquid_phase.properties[t, x].henry["Liq", "CO2"]
            #         / b.psi[t, zf]
            #         + 1
            #     )
            # )

    @blk.Expression(
        blk.flowsheet().time,
        blk.liquid_phase.length_domain,
        doc="""Instantaneous Enhancement factor""",
    )
    def log_instant_E_minus_one(b, t, x):
        if x == b.liquid_phase.length_domain.last():
            return Expression.Skip
        else:
            zf = b.liquid_phase.length_domain.next(x)
            return (
                    b.log_diffus_liq_comp[t, x, "MEA"]
                    + b.liquid_phase.properties[t, x].log_conc_mol_phase_comp_true[
                        "Liq", "MEA"
                    ]
                    + b.log_conc_CO2_bulk[t, x]
                    - log(2)
                    - b.log_diffus_liq_comp[t, x, "CO2"]
                    - b.liquid_phase.properties[t, x].log_conc_mol_phase_comp_true[
                        "Liq", "CO2"
                    ]
            )
    # ======================================================================
    # Enhancement factor model
    # Reference: Jozsef Gaspar,Philip Loldrup Fosbol, (2015)

    blk.conc_interface_MEA = Var(
        blk.flowsheet().time,
        blk.liquid_phase.length_domain,
        bounds=(0, None),
        initialize=0.95,
        units=pyunits.dimensionless,
        doc="""Dimensionless concentration of MEA
                                    at interface """,
    )

    blk.log_conc_interface_MEA = Var(
        blk.flowsheet().time,
        blk.liquid_phase.length_domain,
        bounds=(None, 100),
        initialize=-0.05,
        units=pyunits.dimensionless,
        doc="""Logarithm of conc_interface_MEA""",
    )
    # =============================================================================
    # ------------------------ ORIGINAL -----------------------------------
    @blk.Expression(
        blk.flowsheet().time,
        blk.liquid_phase.length_domain,
        doc="Dimensionless concentration of MEACOO-",
    )
    def conc_interface_MEACOO(b, t, x):
        if x == b.liquid_phase.length_domain.last():
            return Expression.Skip
        else:
            return 1 + (
                    b.liquid_phase.properties[t, x].diffus_phase_comp_true[
                        "Liq", "MEA"
                    ]
                    * b.liquid_phase.properties[t, x].conc_mol_phase_comp_true[
                        "Liq", "MEA"
                    ]
            ) * (1 - b.conc_interface_MEA[t, x]) / (
                           2
                           * b.liquid_phase.properties[t, x].diffus_phase_comp_true[
                               "Liq", "MEACOO_-"
                           ]
                           * b.liquid_phase.properties[t, x].conc_mol_phase_comp_true[
                               "Liq", "MEACOO_-"
                           ]
                   )

    @blk.Expression(
        blk.flowsheet().time,
        blk.liquid_phase.length_domain,
        doc="Dimensionless concentration of MEA+",
    )
    def conc_interface_MEAH(b, t, x):
        if x == b.liquid_phase.length_domain.last():
            return Expression.Skip
        else:
            return 1 + (
                    b.liquid_phase.properties[t, x].diffus_phase_comp_true[
                        "Liq", "MEA"
                    ]
                    * b.liquid_phase.properties[t, x].conc_mol_phase_comp_true[
                        "Liq", "MEA"
                    ]
            ) * (1 - b.conc_interface_MEA[t, x]) / (
                           2
                           * b.liquid_phase.properties[t, x].diffus_phase_comp_true[
                               "Liq", "MEA_+"
                           ]
                           * b.liquid_phase.properties[t, x].conc_mol_phase_comp_true[
                               "Liq", "MEA_+"
                           ]
                   )

    # =============================================================================
    blk.conc_CO2_equil_bulk = Var(
        blk.flowsheet().time,
        blk.liquid_phase.length_domain,
        initialize=1,
        bounds=(0, None),
        doc="""Dimensionless concentration of CO2
                                           at equilibrium with the bulk """,
    )

    @blk.Constraint(
        blk.flowsheet().time,
        blk.liquid_phase.length_domain,
        doc="""Constraint for dimensionless concentration of CO2
                              at equilibrium with the bulk """,
    )
    def conc_CO2_equil_bulk_eqn(b, t, x):
        if x == b.liquid_phase.length_domain.last():
            return Constraint.Skip
        else:
            return (
                    b.conc_CO2_equil_bulk[t, x] * b.conc_interface_MEA[t, x] ** 2
                    == b.conc_CO2_bulk[t, x]
                    * b.conc_interface_MEAH[t, x]
                    * b.conc_interface_MEACOO[t, x]
            )

    @blk.Constraint(
        blk.flowsheet().time,
        blk.liquid_phase.length_domain,
        doc="Defining log_conc_interface_MEA",
    )
    def log_conc_interface_MEA_eqn(b, t, x):
        if x == b.liquid_phase.length_domain.last():
            return Constraint.Skip
        else:
            return (
                    exp(b.log_conc_interface_MEA[t, x])
                    == b.conc_interface_MEA[t, x]
            )

    blk.log_singular_CO2_CO2_ratio = Var(
        blk.flowsheet().time,
        blk.liquid_phase.length_domain,
        initialize=-3,
        bounds=(-12, 12),
        doc="""Logarithm of the ratio of one minus dimensionless concentration CO2 at equilibrium
                       with the bulk to one minus dimensionless concentration of CO2 actually in the bulk.
                       Individually, numerator and denominator go to zero at equilibrium, but the ratio
                       remains (hopefully) well defined mathematically, if not numerically.""",
    )

    @blk.Constraint(
        blk.flowsheet().time,
        blk.liquid_phase.length_domain,
        doc="Equation defining log_singular_CO2_CO2_ratio",
    )
    def log_singular_CO2_CO2_ratio_eqn(b, t, x):
        if x == b.liquid_phase.length_domain.last():
            return Constraint.Skip
        else:
            # Constraint written presently to lose meaning when dimensionless concentrations approach 1
            # If this form doesn't work, we can try one with division in it
            return exp(b.log_singular_CO2_CO2_ratio[t, x]) * (
                    1 - b.conc_CO2_bulk[t, x]
            ) == (1 - b.conc_CO2_equil_bulk[t, x])

    @blk.Constraint(
        blk.flowsheet().time,
        blk.liquid_phase.length_domain,
        doc="Enhancement factor - function of Hatta number",
    )
    def enhancement_factor_eqn1(b, t, x):
        if x == b.liquid_phase.length_domain.last():
            return Constraint.Skip
        else:
            return (
                    b.log_enhancement_factor[t, x]
                    == b.log_hatta_number[t, x]
                    + 0.5 * b.log_conc_interface_MEA[t, x]
                    + b.log_singular_CO2_CO2_ratio[t, x]
            )
            # return (b.enhancement_factor[t, x] * (
            #     1 - b.conc_CO2_bulk[t, x]
            # ) == b.Hatta[t, x] * b.sqrt_conc_interface_MEA[t, x] * (
            #     1 - b.conc_CO2_equil_bulk[t, x]
            # )

    # blk.log_singular_MEA_CO2_ratio = Var(
    #     blk.flowsheet().time,
    #     blk.liquid_phase.length_domain,
    #     initialize=-5,
    #     bounds=(
    #         -12,
    #         3,
    #     ),  # Should be straight up nonpositive, but give it a bit of slack
    #     doc="""Logarithm of the ratio of one minus dimensionless concentration MEA at interface
    #                    to one minus dimensionless concentration of CO2 at equilibrium with the bulk.
    #                    Individually, numerator and denominator go to zero at equilibrium, but the ratio
    #                    remains (hopefully) well defined mathematically, if not numerically.""",
    # )

    # @blk.Constraint(
    #     blk.flowsheet().time,
    #     blk.liquid_phase.length_domain,
    #     doc="Equation defining log_singular_MEA_CO2_ratio",
    # )
    # def log_singular_MEA_CO2_ratio_eqn(b, t, x):
    #     if x == b.liquid_phase.length_domain.last():
    #         return Constraint.Skip
    #     else:
    #         # Constraint written presently to lose meaning when dimensionless concentrations approach 1
    #         # If this form doesn't work, we can try one with division in it
    #         return exp(b.log_singular_MEA_CO2_ratio[t, x]) * (
    #                 1 - b.conc_CO2_bulk[t, x]
    #         ) == (1 - b.conc_interface_MEA[t, x])

    @blk.Constraint(
        blk.flowsheet().time,
        blk.liquid_phase.length_domain,
        doc="Enhancement factor - function of instantaneous enhancement factor",
    )
    def enhancement_factor_eqn2(b, t, x):
        if x == b.liquid_phase.length_domain.last():
            return Constraint.Skip
        else:
            # return b.log_enhancement_factor_minus_one[t, x] == (
            #         b.log_instant_E_minus_one[t, x]
            #         + b.log_singular_MEA_CO2_ratio[t, x]
            # )
            return (b.enhancement_factor[t, x] - 1) * (
                1 - b.conc_CO2_bulk[t, x]
            ) == exp(b.log_instant_E_minus_one[t, x]) * (1 - b.conc_interface_MEA[t, x])

    enhancement_factor_vars = [
        blk.log_enhancement_factor,
        blk.log_hatta_number,
        blk.conc_CO2_bulk,
        blk.log_conc_CO2_bulk,
        blk.conc_interface_MEA,
        blk.log_conc_interface_MEA,
        blk.conc_CO2_equil_bulk,
        blk.log_singular_CO2_CO2_ratio,
        blk.log_rate_constant_MEA,
        blk.log_rate_constant_H2O,
    ]

    enhancement_factor_constraints = [
        blk.hatta_number_eqn,
        blk.log_conc_CO2_bulk_eqn,
        blk.conc_CO2_bulk_eqn,
        blk.conc_CO2_equil_bulk_eqn,
        blk.log_conc_interface_MEA_eqn,
        blk.log_singular_CO2_CO2_ratio_eqn,
        blk.enhancement_factor_eqn1,
        blk.enhancement_factor_eqn2,
        blk.log_rate_constant_MEA_eqn,
        blk.log_rate_constant_H2O_eqn,
    ]
    
    return enhancement_factor_vars, enhancement_factor_constraints

def initialize_enhancement_factor_model(
        blk,
        state_args=None,
        outlvl=idaeslog.NOTSET,
        optarg=None,
        solver=None,
    ):
    # Set up logger for initialization and solve
    init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
    solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")
    # Set solver options
    if optarg is None:
        optarg = {}

    solver_obj = get_solver(solver, optarg)

    long_var_list = []
    long_eqn_list = []

    for t in blk.flowsheet().time:
        for x in blk.liquid_phase.length_domain:
            if x == blk.liquid_phase.length_domain.last():
                continue
            zf = blk.liquid_phase.length_domain.next(x)
            
            calculate_variable_from_constraint(
                blk.log_rate_constant_MEA[t, x], blk.log_rate_constant_MEA_eqn[t, x]
            )
            calculate_variable_from_constraint(
                blk.log_rate_constant_H2O[t, x], blk.log_rate_constant_H2O_eqn[t, x]
            )
            calculate_variable_from_constraint(
                blk.log_conc_CO2_bulk[t, x], blk.conc_CO2_bulk_eqn[t, x]
            )
            calculate_variable_from_constraint(
                blk.conc_CO2_bulk[t, x], blk.log_conc_CO2_bulk_eqn[t, x]
            )
            calculate_variable_from_constraint(
                blk.log_hatta_number[t, x], blk.hatta_number_eqn[t, x]
            )

            if value(blk.conc_CO2_bulk[t, x]) < 1:
                blk.conc_interface_MEA[t, x].value = 0.95
                blk.log_conc_interface_MEA[t, x].value = log(0.95)
                blk.log_singular_CO2_CO2_ratio[t, x].value = log(1.1)
            else:
                blk.conc_interface_MEA[t, x].value = 1.05
                blk.log_conc_interface_MEA[t, x].value = log(1.05)
                blk.log_singular_CO2_CO2_ratio[t, x].value = log(0.9)
            calculate_variable_from_constraint(
                blk.conc_CO2_equil_bulk[t, x],
                blk.conc_CO2_equil_bulk_eqn[t, x],
            )

            entangled_vars = [
                blk.conc_interface_MEA[t, x],
                blk.log_conc_interface_MEA[t, x],
                blk.log_enhancement_factor[t, x],
                blk.log_singular_CO2_CO2_ratio[t, x],
                blk.conc_CO2_equil_bulk[t, x],
                blk.conc_CO2_bulk[t, x],
                blk.log_conc_CO2_bulk[t, x],
            ]
            entangled_eqns = [
                blk.log_conc_interface_MEA_eqn[t, x],
                blk.enhancement_factor_eqn1[t, x],
                blk.log_singular_CO2_CO2_ratio_eqn[t, x],
                blk.enhancement_factor_eqn2[t, x],
                blk.conc_CO2_equil_bulk_eqn[t, x],
                blk.conc_CO2_bulk_eqn[t, x],
                blk.log_conc_CO2_bulk_eqn[t, x],
            ]
            long_var_list += entangled_vars
            long_eqn_list += entangled_eqns

    tmp_blk = create_subsystem_block(long_eqn_list, long_var_list)
    _sub_problem_scaling_suffix(blk, tmp_blk)
    flags = _fix_vars([var for var in tmp_blk.input_vars.values()])
    assert degrees_of_freedom(tmp_blk) == 0
    with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
        res = solver_obj.solve(tmp_blk, tee=slc.tee, symbolic_solver_labels=True)
    assert_optimal_termination(res)
    # if not check_optimal_termination(res):
    #     # IPOPT's solve failed, so let's try a sketchy iterative solution instead
    #     for t in blk.flowsheet().time:
    #         for x in blk.liquid_phase.length_domain:
    #             if x == blk.liquid_phase.length_domain.last():
    #                 continue
    #             zf = blk.liquid_phase.length_domain.next(x)
    #             for j in blk.log_diffus_liq_comp_list:
    #                 blk.log_diffus_liq_comp[t, x, j].value = log(value(blk.liquid_phase.properties[t, x].diffus_phase_comp["Liq", j]))
    #             # CO2 is the only index
    #             blk.log_mass_transfer_coeff_liq[t, x, "CO2"].value = log(value(blk.mass_transfer_coeff_liq[t, x, "CO2"]))

    #             calculate_variable_from_constraint(
    #                 blk.log_hatta_number[t, x], blk.hatta_number_eqn[t, x]
    #             )
    #             calculate_variable_from_constraint(
    #                 blk.log_conc_CO2_bulk[t, x], blk.conc_CO2_bulk_eqn[t, x]
    #             )
    #             calculate_variable_from_constraint(
    #                 blk.conc_CO2_bulk[t, x], blk.log_conc_CO2_bulk_eqn[t, x]
    #             )


    #             if value(blk.conc_CO2_bulk[t, x]) < 1:
    #                 blk.conc_interface_MEA[t, x].value = 0.95
    #                 blk.log_conc_interface_MEA[t, x].value = log(0.95)
    #                 blk.log_singular_CO2_CO2_ratio[t, x].value = log(1.1)
    #             else:
    #                 blk.conc_interface_MEA[t, x].value = 2
    #                 blk.log_conc_interface_MEA[t, x].value = log(2)
    #                 blk.log_singular_CO2_CO2_ratio[t, x].value = log(0.9)
    #             blk.log_singular_MEA_CO2_ratio[t, x].value = value(
    #                 log(
    #                     (1 - blk.conc_interface_MEA[t, x]) / (1 - blk.conc_CO2_bulk[t, x])
    #                 )
    #             )
    #             calculate_variable_from_constraint(
    #                 blk.conc_CO2_equil_bulk[t, x],
    #                 blk.conc_CO2_equil_bulk_eqn[t, x],
    #             )
    #             for k in range(8):
    #             # Solving Eq. 24 from Gaspar et. al.
    #                 instant_E = value(exp(blk.log_instant_E_minus_one[t, x]) + 1)
    #                 a = value(1 - instant_E)
    #                 b = value(exp(blk.log_hatta_number[t, x]) * (blk.conc_CO2_equil_bulk[t, x] - 1))
    #                 c = value(instant_E - blk.conc_CO2_bulk[t, x])
    #                 Yplus = (-b + sqrt(b**2 - 4*a*c))/(2*a)
    #                 Yminus = (-b - sqrt(b**2 - 4*a*c))/(2*a)
    #                 if Yplus > 0:
    #                     assert Yminus < 0
    #                     blk.conc_interface_MEA[t, x].value = Yplus**2
    #                 elif Yminus > 0:
    #                     blk.conc_interface_MEA[t, x].value = Yminus**2
    #                 else:
    #                     raise AssertionError
                    
    #                 blk.log_conc_interface_MEA[t, x].value = log(blk.conc_interface_MEA[t, x].value)

    #                 # Use new value for conc_interface_MEA to calculate enhancement factor
    #                 calculate_variable_from_constraint(
    #                     blk.conc_CO2_equil_bulk[t, x],
    #                     blk.conc_CO2_equil_bulk_eqn[t, x],
    #                 )

    #                 blk.log_singular_CO2_CO2_ratio[t, x].value = value(
    #                     log(
    #                         (1 - blk.conc_CO2_equil_bulk[t, x]) / (1 - blk.conc_CO2_bulk[t, x])
    #                     )
    #                 )

    #                 enhancement_factor = value(
    #                     exp(
    #                         blk.log_hatta_number[t, x]
    #                         + 0.5 * blk.log_conc_interface_MEA[t, x]
    #                         + blk.log_singular_CO2_CO2_ratio[t, x]
    #                     )
    #                 )
    #                 if enhancement_factor > 1:
    #                     blk.log_enhancement_factor_minus_one[t, x].value = log(enhancement_factor - 1)
    #                 else:
    #                     blk.log_enhancement_factor_minus_one[t, x].value = -3
    #                 # Now update values with new value for enhancement factor
    #                 calculate_variable_from_constraint(
    #                     blk.log_conc_CO2_bulk[t, x], blk.conc_CO2_bulk_eqn[t, x]
    #                 )
    #                 calculate_variable_from_constraint(
    #                     blk.conc_CO2_bulk[t, x], blk.log_conc_CO2_bulk_eqn[t, x]
    #                 )
    #             blk.log_singular_MEA_CO2_ratio[t, x].value = log(
    #                 (1 - blk.conc_interface_MEA[t, x]) / (1 - blk.conc_CO2_bulk[t, x])
    #             )

    #     # Resolve after (hopefully) better initialization
    #     assert degrees_of_freedom(tmp_blk) == 0
    #     with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
    #         res = solver_obj.solve(tmp_blk, tee=slc.tee, symbolic_solver_labels=True)
    #     assert_optimal_termination(res)


    _restore_fixedness(flags)