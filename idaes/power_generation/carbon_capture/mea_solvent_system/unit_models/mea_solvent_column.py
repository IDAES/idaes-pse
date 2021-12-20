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
Packed Solvent Column Model for MEA systems
"""

from pyomo.environ import Constraint, exp, SolverStatus, TerminationCondition, Var

from idaes.generic_models.unit_models.column_models.solvent_column import \
    PackedColumnData

from idaes.core import declare_process_block_class
from idaes.core.util import get_solver
import idaes.logger as idaeslog
from idaes.core.util.exceptions import InitializationError


__author__ = "Paul Akula, John Eslick, Anuja Deshpande, Andrew Lee"


@declare_process_block_class("MEAColumn")
class MEAColumnData(PackedColumnData):
    """
    MEA Packed Column Model Class.
    """

    CONFIG = PackedColumnData.CONFIG()

    def build(self):

        super().build()

        # ---------------------------------------------------------------------
        # Unit level sets
        vap_comp = self.config.vapor_side.property_package.component_list
        liq_comp = self.config.liquid_side.property_package.component_list
        equilibrium_comp = vap_comp & liq_comp
        solvent_comp_list = \
            self.config.liquid_side.property_package.solvent_set
        solute_comp_list = self.config.liquid_side.property_package.solute_set

    #     # ---------------------------------------------------------------------
    #     # Vapor-liquid heat transfer coeff modified by Ackmann factor
    #     self.heat_transfer_coeff_base = Var(
    #         self.flowsheet().time,
    #         self.vapor_phase.length_domain,
    #         initialize=100,
    #         doc='Uncorrected vapor-liquid heat transfer coefficient')

    #     def rule_heat_transfer_coeff_Ack(blk, t, x):
    #         if x == self.vapor_phase.length_domain.first():
    #             return Constraint.Skip
    #         else:
    #             Ackmann_factor =\
    #                 sum(blk.vapor_phase.properties[
    #                         t, x].cp_mol_phase_comp['Vap', j] *
    #                     blk.interphase_mass_transfer[t, x, j]
    #                     for j in equilibrium_comp)
    #             return blk.heat_transfer_coeff_base[t, x] == (
    #                 Ackmann_factor /
    #                 (1 - exp(-Ackmann_factor /
    #                          (blk.heat_transfer_coeff_base[t, x] *
    #                           blk.area_interfacial[t, x] *
    #                           blk.area_column))))
    #     self.heat_transfer_Ackmann_correction = Constraint(
    #         self.flowsheet().time,
    #         self.vapor_phase.length_domain,
    #         rule=rule_heat_transfer_coeff_Ack,
    #         doc='Vap-Liq heat transfer correction by Ackmann factor')

    # # =========================================================================
    # # Model initialization routine
    # def initialize(blk,
    #                vapor_phase_state_args=None,
    #                liquid_phase_state_args=None,
    #                state_vars_fixed=False,
    #                outlvl=idaeslog.NOTSET,
    #                solver=None,
    #                optarg=None):
    #     """
    #     Standard Packed Column initialization.

    #     Arguments:
    #         state_args : a dict of arguments to be passed to the property
    #                      package(s) to provide an initial state for
    #                      initialization (see documentation of the specific
    #                      property package) (default = None).
    #         optarg : solver options dictionary object (default=None, use
    #                  default solver options)
    #         solver : str indicating which solver to use during initialization
    #                 (default = None, use IDAES default solver)
    #     """

    #     # Set up logger for initialization and solve
    #     init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
    #     solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")

    #     # Set solver options
    #     opt = get_solver(solver, optarg)

    #     unit_constraints = [
    #         "pressure_at_interface",
    #         "interphase_mass_transfer_eqn",
    #         "liquid_mass_transfer_eqn",
    #         "vapor_mass_transfer_eqn",
    #         "heat_transfer_eqn1",
    #         "heat_transfer_eqn2",
    #         "enthalpy_transfer_eqn1",
    #         "enthalpy_transfer_eqn2"]

    #     # ---------------------------------------------------------------------
    #     # Deactivate unit model level constraints
    #     for c in blk.component_objects(Constraint, descend_into=True):
    #         if c.local_name in unit_constraints:
    #             c.deactivate()

    #     # Fix variables

    #     # Interface pressure
    #     blk.pressure_equil.fix()

    #     # Molar flux
    #     blk.interphase_mass_transfer.fix(0.0)
    #     blk.vapor_phase.mass_transfer_term.fix(0.0)
    #     blk.liquid_phase.mass_transfer_term.fix(0.0)

    #     # Heat transfer rate
    #     blk.vapor_phase.heat.fix(0.0)
    #     blk.liquid_phase.heat.fix(0.0)
    #     blk.vapor_phase.enthalpy_transfer.fix(0.0)
    #     blk.liquid_phase.enthalpy_transfer.fix(0.0)

    #     # ---------------------------------------------------------------------
    #     # Provide state arguments for property package initialization

    #     init_log.info("Step 1: Property Package initialization")
    #     # Initialize vapor_phase properties block
    #     vflag = blk.vapor_phase.initialize(
    #         state_args=vapor_phase_state_args,
    #         outlvl=outlvl,
    #         optarg=optarg,
    #         solver=solver,
    #         hold_state=True)

    #     # Initialize liquid_phase properties block
    #     lflag = blk.liquid_phase.initialize(
    #         state_args=liquid_phase_state_args,
    #         outlvl=outlvl,
    #         optarg=optarg,
    #         solver=solver,
    #         hold_state=True)

    #     init_log.info("Step 2: Steady-State isothermal mass balance")

    #     with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
    #         res = opt.solve(blk, tee=slc.tee)
    #     init_log.info_high("Step 2: {}.".format(idaeslog.condition(res)))

    #     # ---------------------------------------------------------------------
    #     init_log.info('Step 3: Interface equilibrium')

    #     # Activate interface pressure constraint
    #     blk.pressure_equil.unfix()
    #     blk.pressure_at_interface.activate()

    #     with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
    #         res = opt.solve(blk, tee=slc.tee)
    #     init_log.info_high(
    #         "Step 3 complete: {}.".format(idaeslog.condition(res)))

    #     # ---------------------------------------------------------------------
    #     init_log.info('Step 4a: Isothermal absoption')
    #     init_log.info_high("Calculating mass flux")

    #     # Unfix mass transfer terms
    #     blk.interphase_mass_transfer.unfix()

    #     # Activate mass transfer equation in vapor phase
    #     blk.interphase_mass_transfer_eqn.activate()

    #     with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
    #         res = opt.solve(blk, tee=slc.tee)

    #     init_log.info('Step 4b: Isothermal chemical absoption')
    #     init_log.info_high("Adding mass transfer to material balances")

    #     blk.vapor_phase.mass_transfer_term.unfix()
    #     blk.liquid_phase.mass_transfer_term.unfix()
    #     blk.vapor_mass_transfer_eqn.activate()
    #     blk.liquid_mass_transfer_eqn.activate()

    #     # Fix this
    #     with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
    #         res = opt.solve(blk, tee=slc.tee)
    #     init_log.info_high(
    #         "Step 4 complete: {}.".format(idaeslog.condition(res)))

    #     # ---------------------------------------------------------------------
    #     init_log.info('Step 5: Adiabatic chemical absoption')
    #     init_log.info_high("Isothermal to Adiabatic ")

    #     # Unfix heat transfer terms
    #     blk.vapor_phase.heat.unfix()
    #     blk.liquid_phase.heat.unfix()
    #     blk.vapor_phase.enthalpy_transfer.unfix()
    #     blk.liquid_phase.enthalpy_transfer.unfix()

    #     # Activate heat transfer equations
    #     for c in ["heat_transfer_eqn1",
    #               "heat_transfer_eqn2",
    #               "enthalpy_transfer_eqn1",
    #               "enthalpy_transfer_eqn2"]:
    #         getattr(blk, c).activate()

    #     with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
    #         res = opt.solve(blk, tee=slc.tee)

    #     init_log.info_high(
    #         "Step 5 complete: {}.".format(idaeslog.condition(res)))

    #     # ---------------------------------------------------------------------
    #     blk.vapor_phase.release_state(flags=vflag)
    #     blk.liquid_phase.release_state(flags=lflag)

    #     if (res.solver.termination_condition != TerminationCondition.optimal or
    #             res.solver.status != SolverStatus.ok):
    #         raise InitializationError(
    #             f"{blk.name} failed to initialize successfully. Please check "
    #             f"the output logs for more information.")
