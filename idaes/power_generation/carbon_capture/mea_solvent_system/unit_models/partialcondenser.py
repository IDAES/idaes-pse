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
IDAES Partial Condenser

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


@declare_process_block_class("PartialCondenser")
class PartialCondenserData(UnitModelBlockData):
    """Partial Condenser Unit Model."""

    CONFIG = UnitModelBlockData.CONFIG()

    # Configuration template for Phase specific  arguments
    _PhaseCONFIG = ConfigBlock()

    CONFIG.declare("specify_condenser_temperature", ConfigValue(
        default=True,
        domain=Bool,
        description="Indicates if the outlet temperature of the condenser is specified",
        doc="""Indicates if the user specifies the condenser temperature or not
       **default** - True.
        **Valid values:** {
        **True** - Condenser outlet temperature is specified,
        **False** - Condenser outlet temperature is not specified}"""))

    CONFIG.declare("condenser_temperature_spec", ConfigValue(
        default=303.15,
        domain=float,
        description="Condenser outlet temperature in Kelvin",
        doc="Condenser temperature in Kelvin"))

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
        super(PartialCondenserData, self).build()

        # Build Control Volume for Vapor Phase
        self.vapor_phase = ControlVolume0DBlock(default={
            "dynamic": False,
            "property_package": self.config.vapor_side.property_package,
            "property_package_args": self.config.vapor_side.property_package_args})

        self.vapor_phase.add_state_blocks(has_phase_equilibrium=False)

        # Add only momentum balance
        # Combined  Mass and Energy balances for  Liquid and Vapor phases
        # are added as performance equations

        self.vapor_phase.add_momentum_balances(
            balance_type=MomentumBalanceType.pressureTotal,
            has_pressure_change=False)

        # ======================================================================
        # Build Control Volume for liquid Phase
        self.liquid_phase = ControlVolume0DBlock(default={
            "dynamic": False,
            "property_package": self.config.liquid_side.property_package,
            "property_package_args": self.config.liquid_side.property_package_args})

        self.liquid_phase.add_state_blocks(has_phase_equilibrium=False)

        # Energy balance not added: Exit temperature is the equilibruim Temp.

        self.liquid_phase.add_material_balances(
            balance_type=MaterialBalanceType.componentTotal,
            has_mass_transfer=False,
            has_phase_equilibrium=False,
            has_rate_reactions=False)

        self.liquid_phase.add_momentum_balances(
            balance_type=MomentumBalanceType.pressureTotal,
            has_pressure_change=False)

        # ======================================================================
        # Add Ports to control volumes
        # Vapor Phase
        self.add_inlet_port(name="vapor_inlet",
                            block=self.vapor_phase, doc='inlet Port')
        self.add_outlet_port(name="vapor_outlet",
                             block=self.vapor_phase, doc='outlet Port')

        # Liquid Phase : No inlet port
        self.add_outlet_port(name="liquid_outlet",
                             block=self.liquid_phase, doc='outlet Port')
        # ======================================================================
        # Add performace equation method
        self._make_params()
        self._make_performance_method()

    def _make_params(self):

        # Internal variables for custom mass and energy balance
        self.vf = Var(self.flowsheet().time,
                      initialize=0.5,
                      units=None,
                      doc="Vapor phase fraction")

        self.FL = Var(self.flowsheet().time,
                      initialize=1,
                      units=units.kmol,
                      doc="Liquid exit flow")

        self.condenser_duty = Var(self.flowsheet().time,
                             initialize=1,
                             units=units.kW,
                             doc="Condenser heat duty")

    def _make_performance_method(self):

        # =====================================================================
        # Custom Material Balance
        @self.Constraint(self.flowsheet().time,
                         doc="Condenser Material Balance")
        def cond_material_balance(blk, t):
            return blk.vapor_phase.properties_in[t].mole_frac_comp['CO2'] ==\
                   blk.vf[t] * blk.vapor_phase.properties_out[t].mole_frac_comp['CO2']

        @self.Constraint(self.flowsheet().time,
                         doc="Vapor exit flow")
        def vap_exit_flow(blk, t):
            return blk.vapor_phase.properties_out[t].flow_mol == blk.vf[t] *\
                   blk.vapor_phase.properties_in[t].flow_mol

        @self.Constraint(self.flowsheet().time,
                         doc="Liquid exit flow")
        def liq_exit_flow(blk, t):
            return blk.FL[t]*1e3 ==blk.vapor_phase.properties_in[t].flow_mol-\
                   blk.vapor_phase.properties_out[t].flow_mol
        # =====================================================================
        # Custom Energy Balance

        def rule_cp_mol_comp_avg(blk, t, i):
            return 0.5*(blk.vapor_phase.properties_in[t].cp_mol_comp[i] +
                        blk.vapor_phase.properties_out[t].cp_mol_comp[i])
        self.cp_mol_comp_avg = Expression(
                           self.flowsheet().time,
                           self.config.vapor_side.property_package.component_list,
                           rule=rule_cp_mol_comp_avg,
                           doc='''Average molar heat capacity of vapor components
                                  btw Inlet and Outlet temperature''')

        def rule_cp_mol_avg(blk, t):
            return sum(blk.vapor_phase.properties_out[t].mole_frac_comp[i]*
                       blk.cp_mol_comp_avg[t,i] for i in
                       blk.config.vapor_side.property_package.component_list)
        self.cp_mol_avg = Expression(
                           self.flowsheet().time,
                           rule=rule_cp_mol_avg,
                           doc='''Average molar heat capacity of vapor phase
                                  btw Inlet and Outlet temperature''')

        def rule_heat_latent(blk, t):
            return (1-blk.vf[t]) *blk.liquid_phase.properties_out[t].hvap
        self.heat_latent = Expression(self.flowsheet().time,
                                      rule=rule_heat_latent,
                                      doc='latent heat ')

        def rule_heat_sens(blk, t):
            return ((1 - blk.vf[t])*blk.cp_mol_comp_avg[t,'H2O'] +
                    blk.vf[t]*blk.cp_mol_avg[t])*\
                   (blk.vapor_phase.properties_out[t].temperature -
                    blk.vapor_phase.properties_in[t].temperature)
        self.heat_sens = Expression(self.flowsheet().time,
                                    rule=rule_heat_sens,
                                    doc='sensible heat ')

        @self.Constraint(self.flowsheet().time,
                         doc="Energy balance")
        def cond_energy_balance(blk, t):
            return blk.condenser_duty[t] == (blk.heat_sens[t] -
                                             blk.heat_latent[t]) *\
                   blk.vapor_phase.properties_in[t].flow_mol*1e-3

        # =====================================================================
        # Vapor -Liquid Equilibruim- Only water condensation
        @self.Constraint(self.flowsheet().time,
                         doc="Vapor-Liquid Equilibruim ")
        def VLE(blk, t):
            return blk.liquid_phase.properties_out[t].pressure_sat['H2O']==\
                   blk.vapor_phase.properties_out[t].pressure *\
                   blk.vapor_phase.properties_out[t].mole_frac_comp['H2O']

        # ----------------------------------------------------------------------
        # Liquid Phase Inlet Connections and Exit Temperature Relation
        @self.Constraint(self.flowsheet().time,
                         self.config.liquid_side.property_package.component_list,
                         doc="Liquid inlet composition relation")
        def liquid_inlet_comp_mol_frac(blk, t, i):
            if i == 'H2O':
                return blk.liquid_phase.properties_in[t].mole_frac_comp[i] == 1
            else:
                return blk.liquid_phase.properties_in[t].mole_frac_comp[i] == 0

        @self.Constraint(self.flowsheet().time,
                         doc="Liquid phase inlet flow relation")
        def liquid_inlet_flow_mol(blk, t):
            return blk.liquid_phase.properties_in[t].flow_mol == blk.FL[t]*1e3

        @self.Constraint(self.flowsheet().time,
                         doc="Liquid phase inlet temperature relation")
        def liquid_inlet_temp(blk, t):
            return blk.liquid_phase.properties_in[t].temperature ==\
                blk.vapor_phase.properties_out[t].temperature

        @self.Constraint(self.flowsheet().time,
                         doc="Liquid phase exit temperature relation")
        def liquid_exit_temp(blk, t):
            return blk.liquid_phase.properties_in[t].temperature ==\
                   blk.liquid_phase.properties_out[t].temperature

        @self.Constraint(self.flowsheet().time,
                         doc="Liquid phase inlet pressure relation")
        def liquid_inlet_pressure(blk, t):
            return blk.vapor_phase.properties_in[t].pressure ==\
                   blk.liquid_phase.properties_in[t].pressure

        # fix condenser temperature if provided
        if self.config.specify_condenser_temperature:
            for t in self.flowsheet().time:
                self.vapor_phase.properties_out[t].temperature.fix(
                    self.config.condenser_temperature_spec)

        # This is a rare case where the speciation model is not required
        for t in self.flowsheet().time:
            self.liquid_phase.properties_out[t].speciation_model.deactivate()
            self.liquid_phase.properties_out[t].conc_mol_comp_true.fix()
            self.liquid_phase.properties_in[t].speciation_model.deactivate()
            self.liquid_phase.properties_in[t].conc_mol_comp_true.fix()

    def initialize(blk, vapor_state_args=None,
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

        vapor_state_args = {
            'flow_mol': value(blk.vapor_inlet.flow_mol[0]),
            'temperature': value(blk.vapor_inlet.temperature[0]),
            'pressure': value(blk.vapor_inlet.pressure[0]),
            'mole_frac_comp':
            {'H2O': value(blk.vapor_inlet.mole_frac_comp[0, 'H2O']),
             'CO2': value(blk.vapor_inlet.mole_frac_comp[0, 'CO2'])}}

        # ---------------------------------------------------------------------
        # Initialize the vapor INLET properties
        init_log.info('STEP 1: PROPERTY INITIALIZATION')
        init_log.info_high("INLET Properties initialization")
        blk.vapor_phase.properties_in.initialize(
            state_args=vapor_state_args,
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            hold_state=True)

        blk.vf[0].value = 0.5
        flow_out = value(blk.vf[0]*blk.vapor_phase.properties_in[0].flow_mol)

        vapor_state_args_out = {
            'flow_mol': flow_out,
            'temperature': value(blk.vapor_inlet.temperature[0]),
            'pressure': value(blk.vapor_inlet.pressure[0]),
            'mole_frac_comp':
            {'H2O': value(blk.vapor_inlet.mole_frac_comp[0, 'H2O']),
             'CO2': value(blk.vapor_inlet.mole_frac_comp[0, 'CO2'])}}

        # Initialize the Vapor OUTLET properties
        init_log.info_high("Vapor OUTLET Properties initialization")
        blk.vapor_phase.properties_out.initialize(
            state_args=vapor_state_args_out,
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            hold_state=False)

        blk.vapor_phase.properties_out[0].temperature.fix(
                                         blk.config.condenser_temperature_spec)

        blk.FL[0].value = value((1- blk.vf[0]) *
                                blk.vapor_phase.properties_in[0].flow_mol*1e-3)

        # ----------------------------------------------------------------------
        init_log.info('STEP 2: Condenser INITIALIZATION')
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high(
            "STEP 2 Complete: {}.".format(idaeslog.condition(res)))
        if not blk.config.specify_condenser_temperature:
            for t in blk.flowsheet().time:
                blk.vapor_phase.properties_out[t].temperature.unfix()
        init_log.info('INITIALIZATION COMPLETED')
