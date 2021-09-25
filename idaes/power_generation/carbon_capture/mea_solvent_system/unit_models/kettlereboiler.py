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
IDAES Kettle Reboiler

Author : Paul Akula
"""

# Import Pyomo libraries
from pyomo.environ import Expression, Constraint, value, Var, units
from pyomo.common.config import ConfigBlock, ConfigValue, Bool

# Import IDAES cores
from idaes.core import (ControlVolume0DBlock,
                        declare_process_block_class,
                        MaterialBalanceType,
                        MomentumBalanceType,
                        UnitModelBlockData,
                        useDefault)
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util import get_solver
import idaes.logger as idaeslog

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("KettleReboiler")
class KettleReboilerData(UnitModelBlockData):
    """Kettle Reboiler Unit Model."""

    CONFIG = UnitModelBlockData.CONFIG()

    # Configuration template for Phase specific  arguments
    _PhaseCONFIG = ConfigBlock()

    CONFIG.declare("specify_heat_duty", ConfigValue(
        default=True,
        domain=Bool,
        description="Indicates if reboiler duty is specified",
        doc="""Indicates if the user specifies the reboiler duty or not
       **default** - True.
        **Valid values:** {
        **True** - heat duty is specified,
        **False** - heat duty is not specified}"""))

    CONFIG.declare("heat_duty", ConfigValue(
        default=430.61,
        domain=float,
        description="Reboiler heat duty in kW",
        doc="Reboiler heat duty in kW"))

    _PhaseCONFIG.declare("property_package", ConfigValue(
        default=useDefault,
        domain=is_physical_parameter_block,
        description="Property package to use ",
        doc="""Property parameter object used to define property calculations
        **default** - useDefault.
        **Valid values:** {
        **useDefault** - use default package from parent model or flowsheet,
        **PhysicalParameterObject** - a PhysicalParameterBlock object.}"""))

    _PhaseCONFIG.declare("property_package_args", ConfigBlock(
        implicit=True,
        description="Arguments to use for constructing property package",
        doc="""A ConfigBlock with arguments to be passed to
        property block(s) and used when constructing these,
        **default** - None.
        **Valid values:** {
        see property package for documentation.}"""))

    # Create individual config blocks for vapor and liquid phases
    CONFIG.declare("vapor_side",
                   _PhaseCONFIG(doc="Vapor phase config arguments"))
    CONFIG.declare("liquid_side",
                   _PhaseCONFIG(doc="liquid phase config arguments"))

    def build(self):
        # Call UnitModel.build to setup model
        super(KettleReboilerData, self).build()

        # Build Control Volume for liquid Phase
        self.liquid_phase = ControlVolume0DBlock(default={
            "dynamic": self.config.dynamic,
            "property_package": self.config.liquid_side.property_package,
            "property_package_args": self.config.liquid_side.property_package_args})

        self.liquid_phase.add_state_blocks(has_phase_equilibrium=False)

        # Add only momentum balance
        # Combined  Mass and Energy balances for  Liquid and Vapor phases
        # are added as performance equations

        self.liquid_phase.add_momentum_balances(
            balance_type=MomentumBalanceType.pressureTotal,
            has_pressure_change=False)

        # ======================================================================
        # Build Control Volume for vapor Phase
        self.vapor_phase = ControlVolume0DBlock(default={
            "dynamic": self.config.dynamic,
            "property_package": self.config.vapor_side.property_package,
            "property_package_args": self.config.vapor_side.property_package_args})

        self.vapor_phase.add_state_blocks(has_phase_equilibrium=False)

        # Energy balance not added: Exit temperature is the equilibruim Temp.

        self.vapor_phase.add_material_balances(
            balance_type=MaterialBalanceType.componentTotal,
            has_mass_transfer=False,
            has_phase_equilibrium=False,
            has_rate_reactions=False)

        self.vapor_phase.add_momentum_balances(
            balance_type=MomentumBalanceType.pressureTotal,
            has_pressure_change=False)

        # ======================================================================
        # Add Ports to control volumes
        # Liquid Phase
        self.add_inlet_port(name="liquid_inlet",
                            block=self.liquid_phase, doc='inlet Port')
        self.add_outlet_port(name="liquid_outlet",
                             block=self.liquid_phase, doc='outlet Port')

        # Vapor Phase : No inlet port
        self.add_outlet_port(name="vapor_outlet",
                             block=self.vapor_phase, doc='outlet Port')
        # ======================================================================
        # Add performace equation method
        self._make_params()
        self._make_performance_method()

    def _make_params(self):

        # Internal variables for custom mass and energy balance

        self.y = Var(self.flowsheet().time,
                     self.config.vapor_side.property_package.component_list,
                     initialize=0.5,
                     units=None,
                     doc="Vapor component composition")

        self.vf = Var(self.flowsheet().time,
                      initialize=0.5,
                      units=None,
                      doc="Vapor phase fraction")

        self.FV = Var(self.flowsheet().time,
                      initialize=1,
                      units=units.kmol,
                      doc="Vapor exit flow")

        self.heat_duty = Var(self.flowsheet().time,
                             initialize=1,
                             units=units.kW,
                             doc="Reboiler heat duty")

    def _make_performance_method(self):

        # =====================================================================
        # Custom Material Balance
        @self.Constraint(self.flowsheet().time,
                         self.config.liquid_side.property_package.component_list,
                         doc="Reboiler Material Balance")
        def reb_material_balance(blk, t, i):
            if i == 'MEA':
                return blk.liquid_phase.properties_in[t].mole_frac_comp[i] ==\
                    (1 - blk.vf[t]) * blk.liquid_phase.properties_out[t].mole_frac_comp[i]
            else:
                return blk.liquid_phase.properties_in[t].mole_frac_comp[i] ==\
                    (1 - blk.vf[t]) * blk.liquid_phase.properties_out[t].mole_frac_comp[i] + blk.vf[t] * blk.y[t, i]

        @self.Constraint(self.flowsheet().time,
                         doc="Vapor exit flow")
        def vap_exit_flow(blk, t):
            return blk.FV[t] == blk.vf[t] *\
                blk.liquid_phase.properties_in[t].flow_mol*1e-3

        @self.Constraint(self.flowsheet().time,
                         doc="Liquid exit flow")
        def liq_exit_flow(blk, t):
            return blk.liquid_phase.properties_out[t].flow_mol ==\
                   blk.liquid_phase.properties_in[t].flow_mol - blk.FV[t]*1e3

        # =====================================================================
        # Custom Energy Balance

        def rule_heat_vaporization(blk, t):
            return 40650
        self.heat_vaporization = Expression(self.flowsheet().time,
                                            rule=rule_heat_vaporization,
                                            doc='heat of water vaporization in kJ/mol')

        def rule_heat_desorption(blk, t):
            return 95000
        self.heat_desorption = Expression(self.flowsheet().time,
                                          rule=rule_heat_desorption,
                                          doc='heat of CO2 desorption in kJ/mol')

        def rule_heat_latent(blk, t):
            return blk.vf[t] * (blk.y[t, 'H2O'] * blk.heat_vaporization[t] +
                                blk.y[t, 'CO2'] * blk.heat_desorption[t])
        self.heat_latent = Expression(self.flowsheet().time,
                                      rule=rule_heat_latent,
                                      doc='latent heat ')

        def rule_heat_sens(blk, t):
            return (1 - blk.vf[t]) * blk.liquid_phase.properties_out[t].cp_mol *\
                   (blk.liquid_phase.properties_out[t].temperature -
                    blk.liquid_phase.properties_in[t].temperature)
        self.heat_sens = Expression(self.flowsheet().time,
                                    rule=rule_heat_sens,
                                    doc='sensible heat ')

        @self.Constraint(self.flowsheet().time,
                         doc="Energy balance")
        def reb_energy_balance(blk, t):
            return blk.heat_duty[t] == \
                (blk.heat_latent[t] + blk.heat_sens[t]) *\
                blk.liquid_phase.properties_in[t].flow_mol*1e-3

        # =====================================================================
        # Vapor -Liquid Equilibruim
        @self.Constraint(self.flowsheet().time,
                         self.config.vapor_side.property_package.component_list,
                         doc="Vapor-Liquid Equilibruim ")
        def VLE(blk, t, j):
            if j == 'H2O':
                return blk.y[t, j] * blk.liquid_phase.properties_out[t].pressure ==\
                    (blk.liquid_phase.properties_out[t].vol_mol *
                     blk.liquid_phase.properties_out[t].conc_mol_comp_true[j] *
                     blk.liquid_phase.properties_out[t].pressure_sat[j])
            elif j == 'CO2':
                return blk.y[t, j] * blk.liquid_phase.properties_out[t].pressure ==\
                    (blk.liquid_phase.properties_out[t].conc_mol_comp_true[j] *
                     blk.liquid_phase.properties_out[t].henry_N2O_analogy)

        # ----------------------------------------------------------------------
        # Vapor Phase Inlet Connections and Exit Temperature Relation

        @self.Constraint(self.flowsheet().time,
                         self.config.vapor_side.property_package.component_list,
                         doc="Vapor phase inlet composition relation")
        def vapor_inlet_comp(blk, t, i):
            return blk.vapor_phase.properties_in[t].mole_frac_comp[i] == blk.y[t, i]

        @self.Constraint(self.flowsheet().time,
                         doc="Vapor phase inlet flow relation")
        def vapor_inlet_flow(blk, t):
            return blk.vapor_phase.properties_in[t].flow_mol == blk.FV[t]*1e3

        @self.Constraint(self.flowsheet().time,
                         doc="Vapor phase inlet temperature relation")
        def vapor_inlet_temp(blk, t):
            return blk.vapor_phase.properties_in[t].temperature ==\
                blk.liquid_phase.properties_out[t].temperature

        @self.Constraint(self.flowsheet().time,
                         doc="Vapor phase exit temperature relation")
        def vapor_exit_temp(blk, t):
            return blk.vapor_phase.properties_in[t].temperature ==\
                blk.vapor_phase.properties_out[t].temperature

        @self.Constraint(self.flowsheet().time,
                         doc="Vapor phase inlet pressure relation")
        def vapor_inlet_pressure(blk, t):
            return blk.vapor_phase.properties_in[t].pressure ==\
                  blk.liquid_phase.properties_in[t].pressure

        # fix heat duty if provided
        if self.config.specify_heat_duty:
            for t in self.flowsheet().time:
                self.heat_duty[t].fix(self.config.heat_duty)

    def initialize(blk, liquid_state_args=None,
                   outlvl=idaeslog.NOTSET, solver=None, optarg=None):
        '''
        Initialisation routine for reboiler unit (default solver ipopt)

        Keyword Arguments:
            state_args : a dict of arguments to be passed to the property
                           package(s) to provide an initial state for
                           initialization (see documentation of the specific
                           property package) (default = {}).
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None, use
                     default solver options)
            solver : str indicating which solver to use during
                     initialization (default = None)
        Returns:
            None
        '''
        # Set solver options
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag='unit')
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")

        # Create solver
        opt = get_solver(solver, optarg)

        liquid_state_args = {
            'flow_mol': value(blk.liquid_inlet.flow_mol[0]),
            'temperature': value(blk.liquid_inlet.temperature[0]),
            'pressure': value(blk.liquid_inlet.pressure[0]),
            'mole_frac_comp':
            {'H2O': value(blk.liquid_inlet.mole_frac_comp[0, 'H2O']),
             'CO2': value(blk.liquid_inlet.mole_frac_comp[0, 'CO2']),
             'MEA': value(blk.liquid_inlet.mole_frac_comp[0, 'MEA'])}}

        # ---------------------------------------------------------------------
        # Initialize the liquid INLET properties
        init_log.info('STEP 1: PROPERTY INITIALIZATION')
        init_log.info_high("INLET Properties initialization")
        blk.liquid_phase.properties_in.initialize(
            state_args=liquid_state_args,
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            hold_state=True)

        blk.vf[0].value = 0.5
        flow_out = value((1 - blk.vf[0])
                         *blk.liquid_phase.properties_in[0].flow_mol)

        liquid_state_args_out = {
            'flow_mol': flow_out,
            'temperature': value(blk.liquid_inlet.temperature[0]),
            'pressure': value(blk.liquid_inlet.pressure[0]),
            'mole_frac_comp':
            {'H2O': value(blk.liquid_inlet.mole_frac_comp[0, 'H2O']),
             'CO2': value(blk.liquid_inlet.mole_frac_comp[0, 'CO2']),
             'MEA': value(blk.liquid_inlet.mole_frac_comp[0, 'MEA'])}}

        # Initialize the Liquid OUTLET properties
        init_log.info_high("Liquid OUTLET Properties initialization")
        blk.liquid_phase.properties_out.initialize(
            state_args=liquid_state_args_out,
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            hold_state=False)

        blk.heat_duty[0].fix(blk.config.heat_duty)

        blk.FV[0].value = value(blk.vf[0] *
                                blk.liquid_phase.properties_in[0].flow_mol*1e-3)

        # ----------------------------------------------------------------------
        init_log.info('STEP 2: Reboiler INITIALIZATION')
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high(
            "STEP 2 Complete: {}.".format(idaeslog.condition(res)))
        if not blk.config.specify_heat_duty:
            blk.heat_duty.unfix()
        init_log.info('INITIALIZATION COMPLETED')
