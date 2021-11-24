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
Gas-solvent contactor column model for carbon capture using MEA.

"""

# Import Pyomo libraries
from pyomo.environ import (Constraint,
                           exp,
                           Expression,
                           NonNegativeReals,
                           Param,
                           Reals,
                           units as pyunits,
                           Var)
from pyomo.util.calc_var_value import calculate_variable_from_constraint

# Import IDAES Libraries
from idaes.core import declare_process_block_class
from idaes.core.util.constants import Constants

from idaes.generic_models.unit_models.column_models.solvent_column import \
    PackedColumnData
from idaes.core.util import get_solver
import idaes.logger as idaeslog


__author__ = "Andrew Lee, Paul Akula, John Eslick, Anuja Deshpande"


# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("MEAColumn")
class MEAColumnData(PackedColumnData):
    """
    Model for carbon capture column using MEA solvent

    """
    CONFIG = PackedColumnData.CONFIG()

    # -------------------------------------------------------------------------
    def build(self):
        """
        Begin building model (pre-DAE transformation).

        Args:
            None

        Returns:
            None
        """
        # Call UnitModel.build to build default attributes
        super().build()

        # ---------------------------------------------------------------------
        # Define component lists for all species and those in equilibrium
        vap_comp = self.config.vapor_side.property_package.component_list
        liq_comp = self.config.liquid_side.property_package.component_list
        equilibrium_comp = vap_comp & liq_comp

        # ---------------------------------------------------------------------
        # Hydrodynamics
        self.velocity_liq = Var(self.flowsheet().time,
                                self.liquid_phase.length_domain,
                                units=pyunits.m/pyunits.s,
                                domain=NonNegativeReals,
                                initialize=0.01,
                                doc='Liquid superficial velocity')

        @self.Constraint(self.flowsheet().time,
                         self.liquid_phase.length_domain,
                         doc="Liquid superficial velocity")
        def eq_velocity_liq(blk, t, x):
            return blk.liquid_phase.properties[t, x].flow_mol == (
                blk.velocity_liq[t, x] * blk.area_column *
                blk.liquid_phase.properties[t, x].dens_mol)

        self.velocity_vap = Var(self.flowsheet().time,
                                self.vapor_phase.length_domain,
                                units=pyunits.m/pyunits.s,
                                domain=NonNegativeReals,
                                initialize=2,
                                doc='Vapor superficial velocity')

        @self.Constraint(self.flowsheet().time,
                         self.vapor_phase.length_domain,
                         doc="Vapor superficial velocity")
        def eq_velocity_vap(blk, t, x):
            return blk.vapor_phase.properties[t, x].flow_mol == (
                blk.velocity_vap[t, x] * blk.area_column *
                blk.vapor_phase.properties[t, x].dens_mol)

        # ---------------------------------------------------------------------
        # Interfacial Area model ([m2/m3]):
        # Reference: Tsai correlation,regressed by Chinen et al. 2018
        self.area_interfacial_parameter_A = Var(
            initialize=0.6486,
            units=pyunits.dimensionless,
            doc='Interfacial area parameter A')

        self.area_interfacial_parameter_B = Var(
            initialize=0.12,
            units=pyunits.dimensionless,
            doc='Interfacial area parameter B')

        self.area_interfacial_parameter_A.fix()
        self.area_interfacial_parameter_B.fix()

        def rule_interfacial_area(blk, t, x):
            # Assume this is purely empirical and drop units
            return blk.area_interfacial[t, x] == (
                blk.packing_specific_surface_area *
                blk.area_interfacial_parameter_A *
                (blk.liquid_phase.properties[t, x].dens_mass /
                 pyunits.get_units(
                     blk.liquid_phase.properties[t, x].dens_mass) /
                 blk.liquid_phase.properties[t, x].surf_tens_phase["Liq"] *
                 pyunits.get_units(
                     blk.liquid_phase.properties[t, x].surf_tens_phase["Liq"]) *
                 (blk.velocity_liq[t, x] /
                  pyunits.get_units(blk.velocity_liq[t, x]))**(4/3)
                 )**blk.area_interfacial_parameter_B)

        self.area_interfacial_model = Constraint(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            rule=rule_interfacial_area,
            doc='Specific interfacial area model')

        # ---------------------------------------------------------------------
        # Liquid holdup model
        # Reference: Tsai correlation, regressed by Chinen et al. 2018
        self.holdup_parameter_A = Var(
            initialize=24.2355,
            units=pyunits.dimensionless,
            doc='Holdup parameter A')
        self.holdup_parameter_A.fix()

        self.holdup_parameter_B = Var(
            initialize=0.6471,
            units=pyunits.dimensionless,
            doc='Holdup parameter B')
        self.holdup_parameter_B.fix()

        def rule_holdup_liq(blk, t, x):
            # Assume this is purely empirical and drop units
            return blk.liquid_holdup_fraction[t, x] == (
                blk.holdup_parameter_A *
                (blk.velocity_liq[t, x] /
                 pyunits.get_units(blk.velocity_liq[t, x]) *
                 (blk.liquid_phase.properties[t, x].visc_d_phase["Liq"] /
                  pyunits.get_units(
                      blk.liquid_phase.properties[t, x].visc_d_phase["Liq"]) /
                  blk.liquid_phase.properties[t, x].dens_mass *
                  pyunits.get_units(
                      blk.liquid_phase.properties[t, x].dens_mass))**(1/3)
                 )**blk.holdup_parameter_B)

        self.liquid_holdup_model = Constraint(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            rule=rule_holdup_liq,
            doc='Volumetric liquid holdup fraction model')

        # ---------------------------------------------------------------------
        # Heat transfer coefficients, Chilton-Colburn analogy
        # Vapor-liquid heat transfer coefficient [J/m2.s.K]

        def rule_heat_transfer_coeff(blk, t, x):
            if x == self.vapor_phase.length_domain.first():
                return Expression.Skip
            else:
                return (
                    blk.mass_transfer_coeff_vap_comp[t, x, 'CO2'] *
                    blk.vapor_phase.properties[t, x].pressure *
                    blk.vapor_phase.properties[t, x].cp_mol_phase["Vap"] *
                    (blk.vapor_phase.properties[t, x].therm_cond_phase["Vap"] /
                     (blk.vapor_phase.properties[t, x].dens_mol *
                      blk.vapor_phase.properties[t, x].cp_mol_phase["Vap"] *
                      blk.vapor_phase.properties[t, x].diffus_phase_comp[
                          "Vap", 'CO2'])
                     )**(2/3))

        self.h_v = Expression(self.flowsheet().time,
                              self.vapor_phase.length_domain,
                              rule=rule_heat_transfer_coeff,
                              doc='Vapor-liquid heat transfer coefficient')

        # Overall heat transfer coefficient modified by Ackmann factor
        def rule_heat_transfer_coeff_Ack(blk, t, x):
            if x == self.vapor_phase.length_domain.first():
                return Constraint.Skip
            else:
                xl = self.liquid_phase.length_domain.prev(x)
                Ackmann_factor =\
                    (blk.vapor_phase.properties[t, x].cp_mol_phase_comp[
                        'Vap', 'CO2'] *
                     blk.liquid_phase.mass_transfer_term[t, xl, 'Liq', 'CO2'] +
                     blk.vapor_phase.properties[t, x].cp_mol_phase_comp[
                         'Vap', 'H2O'] *
                     blk.liquid_phase.mass_transfer_term[t, xl, 'Liq', 'H2O'])
                return blk.heat_transfer_coeff[t, x] == (
                    Ackmann_factor /
                    (1 -
                     exp(-Ackmann_factor /
                         (blk.h_v[t, x] * blk.area_interfacial[t, x] *
                          blk.area_column))))
        self.heat_transfer_model = Constraint(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            rule=rule_heat_transfer_coeff_Ack,
            doc='Overall heat transfer coefficient corrected by Ackmann factor'
            )

    # -------------------------------------------------------------------------
    # Model initialization routine
    def initialize(blk,
                   vapor_state_args=None,
                   liquid_state_args=None,
                   state_vars_fixed=False,
                   outlvl=idaeslog.NOTSET,
                   solver=None,
                   optarg=None):
        """
        Column initialization.

        Arguments:
            vapor_state_args : a dict of arguments to be passed to the vapor
                property package to provide an initial state for initialization
                (see documentation of the specific property package)
                (default = None).
            liquid_state_args : a dict of arguments to be passed to the liquid
                property package to provide an initial state for initialization
                (see documentation of the specific property package)
                (default = None).
            optarg : solver options dictionary object (default=None, use
                      default solver options)
            solver : str indicating which solver to use during initialization
                    (default = None, use IDAES default solver)
        """

        # Set up logger for initialization and solve
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")

        # Set solver options
        opt = get_solver(solver, optarg)

        # ---------------------------------------------------------------------
        # Initialize vapor phase control volume block
        vflags = blk.vapor_phase.initialize(
                outlvl=outlvl,
                optarg=optarg,
                solver=solver,
                state_args=vapor_state_args,
        )

        init_log.info_high('Initialization Step 1 Complete.')

        # ---------------------------------------------------------------------
        # Initialize liquid phase control volume block
        lflags = blk.liquid_phase.initialize(
                outlvl=outlvl,
                optarg=optarg,
                solver=solver,
                state_args=liquid_state_args,
        )

        init_log.info_high('Initialization Step 2 Complete.')

        # ---------------------------------------------------------------------
        # Deactivate heat & mass transfer and solve for intermediate variables
        blk.mass_transfer_liq.deactivate()
        blk.mass_transfer_vap.deactivate()
        blk.vapor_phase_heat_transfer.deactivate()
        blk.column_energy_balance.deactivate()
        blk.vapor_enthalpy_transfer_equation.deactivate()
        blk.enthalpy_transfer_balance.deactivate()

        blk.vapor_phase.mass_transfer_term.fix(0)
        blk.liquid_phase.mass_transfer_term.fix(0)
        blk.vapor_phase.heat.fix(0)
        blk.liquid_phase.heat.fix(0)
        blk.vapor_phase.enthalpy_transfer.fix(0)
        blk.liquid_phase.enthalpy_transfer.fix(0)

        blk.heat_transfer_coeff.fix(6000)
        blk.heat_transfer_model.deactivate()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            results = opt.solve(blk, tee=slc.tee)

        init_log.info_high(
            "Initialization Step 3 {}.".format(idaeslog.condition(results))
        )

        # ---------------------------------------------------------------------
        # Reactivate mass transfer and solve
        blk.mass_transfer_liq.activate()
        blk.mass_transfer_vap.activate()

        blk.vapor_phase.mass_transfer_term.unfix()
        blk.liquid_phase.mass_transfer_term.unfix()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            results = opt.solve(blk, tee=slc.tee)

        init_log.info_high(
            "Initialization Step 4 {}.".format(idaeslog.condition(results))
            )

        # ---------------------------------------------------------------------
        # Reactivate enthalpy transfer and solve
        blk.vapor_enthalpy_transfer_equation.activate()
        blk.enthalpy_transfer_balance.activate()

        blk.vapor_phase.enthalpy_transfer.unfix()
        blk.liquid_phase.enthalpy_transfer.unfix()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            results = opt.solve(blk, tee=slc.tee)

        init_log.info_high(
            "Initialization Step 5 {}.".format(idaeslog.condition(results))
        )

        # ---------------------------------------------------------------------
        # Reactivate heat transfer and solve
        blk.vapor_phase_heat_transfer.activate()
        blk.column_energy_balance.activate()

        blk.vapor_phase.heat.unfix()
        blk.liquid_phase.heat.unfix()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            results = opt.solve(blk, tee=slc.tee)

        init_log.info_high(
            "Initialization Step 6 {}.".format(idaeslog.condition(results))
        )

        # ---------------------------------------------------------------------
        # Add heat transfer model
        blk.heat_transfer_model.activate()
        blk.heat_transfer_coeff.unfix()
        for k in blk.heat_transfer_model:
            calculate_variable_from_constraint(
                blk.heat_transfer_coeff[k],
                blk.heat_transfer_model[k])

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            results = opt.solve(blk, tee=slc.tee)

        init_log.info_high(
            "Initialization Step 7 {}.".format(idaeslog.condition(results))
        )

        # ---------------------------------------------------------------------
        # Release states
        blk.vapor_phase.release_state(vflags, outlvl)
        blk.liquid_phase.release_state(lflags, outlvl)

        init_log.info('Initialization Complete: {}'
                      .format(idaeslog.condition(results)))
