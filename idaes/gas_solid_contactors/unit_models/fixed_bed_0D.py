##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
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
IDAES 0D Fixed Bed Reactor model.
"""

# Import Pyomo libraries
from pyomo.environ import (Var, Param, Reals, value, SolverFactory,
                           TransformationFactory, Constraint,
                           TerminationCondition)
from pyomo.dae import DerivativeVar, ContinuousSet # added ContinuousSet
from pyomo.common.config import ConfigBlock, ConfigValue, In
from pyomo.util.calc_var_value import calculate_variable_from_constraint # added

# Import IDAES cores
from idaes.core import (ControlVolume0DBlock,
                        declare_process_block_class,
                        MaterialBalanceType,
                        EnergyBalanceType,
                        MomentumBalanceType,
                        UnitModelBlockData,
                        useDefault)
from idaes.core.util.config import (is_physical_parameter_block,
                                    is_reaction_parameter_block)
from idaes.core.util.exceptions import (ConfigurationError,
                                        BurntToast) # added
from idaes.core.util.tables import create_stream_table_dataframe # added
from idaes.core.util.misc import add_object_reference
from idaes.core.util.constants import Constants as constants # added
from idaes.core.util.math import smooth_abs # added
from idaes.core.util import get_solver # added
from idaes.logger import (getInitLogger, condition, solver_log,
                          DEBUG, getSolveLogger)

__author__ = "Andrew Lee" # edited by Brandon Paul


# Assumptions:
# Only Vap and Sol phases, which are explicitly named as such
# Perfect mixing in Vap phase
# Static solid phase

# Need to have build-on-demand properties in ractions for this to work

@declare_process_block_class("FixedBed0D")
class FixedBed0DData(UnitModelBlockData):
    """
    0D Fixed Bed Reactor Model Class
    """
    CONFIG = ConfigBlock()
    CONFIG.declare("dynamic", ConfigValue(
        domain=In([True]),
        default=True,
        description="Dynamic model flag - must be True",
        doc="""Indicates whether this model will be dynamic or not,
**default** = True. Fixed beds must be dynamicr."""))
    CONFIG.declare("has_holdup", ConfigValue(
        default=True,
        domain=In([True]),
        description="Holdup construction flag - must be True",
        doc="""Indicates whether holdup terms should be constructed or not.
**default** - True. Fixed bed reactors must be dynamic, thus this must be
True."""))
    CONFIG.declare("material_balance_type", ConfigValue(
        default=MaterialBalanceType.useDefault,
        domain=In(MaterialBalanceType),
        description="Material balance construction flag",
        doc="""Indicates what type of mass balance should be constructed,
**default** - MaterialBalanceType.useDefault.
**Valid values:** {
**MaterialBalanceType.useDefault - refer to property package for default
balance type
**MaterialBalanceType.none** - exclude material balances,
**MaterialBalanceType.componentPhase** - use phase component balances,
**MaterialBalanceType.componentTotal** - use total component balances,
**MaterialBalanceType.elementTotal** - use total element balances,
**MaterialBalanceType.total** - use total material balance.}"""))
    CONFIG.declare("energy_balance_type", ConfigValue(
        default=EnergyBalanceType.useDefault,
        domain=In(EnergyBalanceType),
        description="Energy balance construction flag",
        doc="""Indicates what type of energy balance should be constructed,
**default** - EnergyBalanceType.useDefault.
**Valid values:** {
**EnergyBalanceType.useDefault - refer to property package for default
balance type
**EnergyBalanceType.none** - exclude energy balances,
**EnergyBalanceType.enthalpyTotal** - single enthalpy balance for material,
**EnergyBalanceType.enthalpyPhase** - enthalpy balances for each phase,
**EnergyBalanceType.energyTotal** - single energy balance for material,
**EnergyBalanceType.energyPhase** - energy balances for each phase.}"""))
    CONFIG.declare("momentum_balance_type", ConfigValue(
        default=MomentumBalanceType.pressureTotal,
        domain=In(MomentumBalanceType),
        description="Momentum balance construction flag",
        doc="""Indicates what type of momentum balance should be constructed,
**default** - MomentumBalanceType.pressureTotal.
**Valid values:** {
**MomentumBalanceType.none** - exclude momentum balances,
**MomentumBalanceType.pressureTotal** - single pressure balance for material,
**MomentumBalanceType.pressurePhase** - pressure balances for each phase,
**MomentumBalanceType.momentumTotal** - single momentum balance for material,
**MomentumBalanceType.momentumPhase** - momentum balances for each phase.}"""))
#    CONFIG.declare("has_heat_transfer", ConfigValue(
#        default=False,
#        domain=In([True, False]),
#        description="External Heat transfer term construction flag",
#        doc="""Indicates whether terms for external heat transfer should be
#constructed, **default** - False.
#**Valid values:** {
#**True** - include external heat transfer terms,
#**False** - exclude external heat transfer terms.}"""))
#    CONFIG.declare("has_equilibrium_reactions", ConfigValue(
#        default=False,
#        domain=In([True, False]),
#        description="Equilibrium reaction construction flag",
#        doc="""Indicates whether terms for equilibrium controlled reactions
#should be constructed,
#**default** - True.
#**Valid values:** {
#**True** - include equilibrium reaction terms,
#**False** - exclude equilibrium reaction terms.}"""))
#    CONFIG.declare("has_phase_equilibrium", ConfigValue(
#        default=False,
#        domain=In([True, False]),
#        description="Phase equilibrium construction flag",
#        doc="""Indicates whether terms for phase equilibrium should be
#constructed,
#**default** = False.
#**Valid values:** {
#**True** - include phase equilibrium terms
#**False** - exclude phase equilibrium terms.}"""))
#    CONFIG.declare("has_heat_of_reaction", ConfigValue(
#        default=False,
#        domain=In([True, False]),
#        description="Heat of reaction term construction flag",
#        doc="""Indicates whether terms for heat of reaction terms should be
#constructed,
#**default** - False.
#**Valid values:** {
#**True** - include heat of reaction terms,
#**False** - exclude heat of reaction terms.}"""))
    CONFIG.declare("has_pressure_change", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="Pressure change term construction flag",
        doc="""Indicates whether terms for pressure change should be
constructed,
**default** - False.
**Valid values:** {
**True** - include pressure change terms,
**False** - exclude pressure change terms.}"""))
    CONFIG.declare("gas_property_package", ConfigValue(
        default=useDefault,
        domain=is_physical_parameter_block,
        description="Property package to use for gas phase",
        doc="""Property parameter object used to define property calculations
for the gas phase, **default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PhysicalParameterObject** - a PhysicalParameterBlock object.}"""))
    CONFIG.declare("gas_property_package_args", ConfigBlock(
        implicit=True,
        description="Arguments to use for constructing gas phase "
        "property packages",
        doc="""A ConfigBlock with arguments to be passed to a gas phase
property block(s) and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}"""))
    CONFIG.declare("solid_property_package", ConfigValue(
        default=useDefault,
        domain=is_physical_parameter_block,
        description="Property package to use for solid phase",
        doc="""Property parameter object used to define property calculations
for the solid phase, **default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PhysicalParameterObject** - a PhysicalParameterBlock object.}"""))
    CONFIG.declare("solid_property_package_args", ConfigBlock(
        implicit=True,
        description="Arguments to use for constructing solid phase "
        "property packages",
        doc="""A ConfigBlock with arguments to be passed to a solid phase
property block(s) and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}"""))
    CONFIG.declare("reaction_package", ConfigValue(
        default=None,
        # TODO: Had to remove domain due to limitations in ReactionBase for heterogeneous systems
        description="Reaction package to use for unit",
        doc="""Reaction parameter object used to define reaction calculations,
**default** - None.
**Valid values:** {
**None** - no reaction package,
**ReactionParameterBlock** - a ReactionParameterBlock object.}"""))
    CONFIG.declare("reaction_package_args", ConfigBlock(
        implicit=True,
        description="Arguments to use for constructing reaction packages",
        doc="""A ConfigBlock with arguments to be passed to a reaction block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see reaction package for documentation.}"""))

    def build(self):
        # Call UnitModel.build to setup dynamics
        super(FixedBed0DData, self).build()

        # Build Gas Phase Control Volume
        self.gas_phase = ControlVolume0DBlock(default={
                "dynamic": self.config.dynamic,
                "has_holdup": self.config.has_holdup,
                "property_package": self.config.gas_property_package,
                "property_package_args": self.config.gas_property_package_args})#,
#                "reaction_package": self.config.reaction_package,
#                "reaction_package_args": self.config.reaction_package_args})

        self.gas_phase.add_geometry()

        self.gas_phase.add_state_blocks(
                has_phase_equilibrium=False)

#        self.gas_phase.add_reaction_blocks(
#                has_equilibrium=self.config.has_equilibrium_reactions)

        self.gas_phase.add_material_balances(
            balance_type=self.config.material_balance_type,
            has_mass_transfer=True)

        self.gas_phase.add_energy_balances(
            balance_type=self.config.energy_balance_type,
            has_heat_of_reaction=False,
            has_heat_transfer=True)

        self.gas_phase.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=self.config.has_pressure_change)

        # Build Solid Phase StateBlock
        # As there is no solid flow, only need a single state block
        self.solids = self.config.solid_property_package.state_block_class(
                self.flowsheet().time,
                default={"parameters": self.config.solid_property_package})

        # Create Reaction Block
#        self.reactions = self.config.reaction_package.reaction_block_class(
#                self.flowsheet().time,
#                default={"parameters": self.config.reaction_package,
#                         # Had to change these keys to get 'working' with reduction packages
#                         #"solid_state_block": self.solids,
#                         # ^ old
#                         "state_block": self.solids,
#                         # ^ new
#                         "gas_state_block": self.gas_phase.properties_out})

        tmp_dict = dict(**self.config.reaction_package_args)
        tmp_dict["gas_state_block"] = self.gas_phase.properties_out
        tmp_dict["solid_state_block"] = self.solids
        tmp_dict["parameters"] = (self.config.reaction_package)
        self.reactions = (
                self.config.reaction_package.
                reaction_block_class(
                    self.flowsheet().time,
                    doc="Reaction properties in control volume",
                    default=tmp_dict))

        # Solid phase material balance
        # Volume of solid
        self.volume_solid = Var(self.flowsheet().config.time,
                                initialize=1.0,
                                doc="Solids phase volume")

        # Accumulation equal to mass transfer/reaction
        self.solids_material_holdup = Var(
                self.flowsheet().time,
                self.config.solid_property_package.component_list,
                initialize=1,
                doc="Solid phase component holdups")

        @self.Constraint(self.flowsheet().config.time,
                         self.config.solid_property_package.component_list,
                         doc="Solid phase material holdup constraints")
        def solids_material_holdup_constraints(b, t, j):
            return b.solids_material_holdup[t, j] == (
                  self.volume_solid[t] *
                  b.solids[t].get_material_density_terms("Sol", j))

        self.solids_material_accumulation = DerivativeVar(
                    self.solids_material_holdup,
                    wrt=self.flowsheet().config.time,
                    doc="Solids material accumulation")

        @self.Constraint(self.flowsheet().config.time,
                         self.config.solid_property_package.component_list,
                         doc="Solid phase material accumulation constraints")
        def solids_material_accumulation_constraints(b, t, j):
            if t == self.flowsheet().config.time.first():
                return Constraint.Skip
            else:
                return (b.solids_material_accumulation[t, j]*1e5 ==
                        1e5*b.volume_solid[t] * b.solids[t]._params.mw[j] *
                        sum(b.reactions[t].reaction_rate[r] *
                            b.config.reaction_package
                            .rate_reaction_stoichiometry[r, "Sol", j] for r in
                            b.config.reaction_package.rate_reaction_idx))

        @self.Constraint(self.flowsheet().config.time,
                         self.config.gas_property_package.component_list,
                         doc="Gas phase mass transfer constraints")
        def gas_phase_mass_transfer_constraints(b, t, j):
            return (b.gas_phase.mass_transfer_term[t, "Vap", j]*1e5 ==
                    1e5*b.volume_solid[t] *
                    sum(b.reactions[t].reaction_rate[r] *
                        b.config.reaction_package
                        .rate_reaction_stoichiometry[r, "Vap", j]
                        for r in b.config.reaction_package.rate_reaction_idx))

        # Add solid mass variable and constraint for TGA tracking
        self.mass_solids = Var(self.flowsheet().config.time,
                               doc="Total mass of solids")

        @self.Constraint(self.flowsheet().config.time,
                         doc="Calculating total mass of solids")
        def mass_solids_eq(b, t):
            return b.mass_solids[t] == b.volume_solid[t]*sum(
                    b.solids[t].get_material_density_terms("Sol", j)
                    for j in b.config.solid_property_package.component_list)

        # Sum of mass fractions at all time equals 1
        @self.Constraint(self.flowsheet().config.time,
                         doc="Sum of mass fractions at all time")
        def sum_component_eqn(b, t):
            if t == b.flowsheet().config.time.first():
                return Constraint.Skip
            else:
                return 1e2 == 1e2*sum(b.solids[t].mass_frac[j]
                                      for j in b.solids[t].
                                      _params.component_list)

        # Create total volume variable
        self.volume_reactor = Var(initialize=0.1,
                                  doc="Total reactor volume")

        # Add object reference to gas volume
        add_object_reference(self,
                             "volume_gas",
                             self.gas_phase.volume)

        # Calculate volume of gas
        @self.Constraint(self.flowsheet().config.time,
                         doc="Volume of gas constraint")
        def total_volume_constraint(b, t):
            return b.volume_reactor == b.volume_gas[t] + b.volume_solid[t]

        # Solid phase energy balance
        # Accumulation equal to heat transfer
        self.solids_energy_holdup = Var(
                self.flowsheet().time,
                initialize=1,
                doc="Solid phase energy holdup")

        @self.Constraint(self.flowsheet().config.time,
                         doc="Solid phase energy holdup constraints")
        def solids_energy_holdup_constraints(b, t):
            return b.solids_energy_holdup[t] == (
                  self.volume_solid[t] *
                  b.solids[t].get_energy_density_terms("Sol"))

        self.solids_energy_accumulation = DerivativeVar(
                    self.solids_energy_holdup,
                    doc="Solids energy accumulation")

        @self.Constraint(self.flowsheet().config.time,
                         doc="Solid phase energy accumulation constraints")
        def solids_energy_accumulation_constraints(b, t):
            if t == self.flowsheet().config.time.first():
                return Constraint.Skip
            else:
                return b.solids_energy_accumulation[t] == \
                    - self.gas_phase.heat[t]

#        # Heat Transfer Model
#        # TODO: For now assume thermal equilibrium
#        @self.Constraint(self.flowsheet().config.time,
#                         doc="Gas-solids temperature equilibrium constraint")
#        def solids_temperature_constraint(b, t):
#            return b.solids[t].temperature == \
#                b.gas_phase.properties_out[t].temperature

        # Solids pressure equal to gas pressure
        @self.Constraint(self.flowsheet().config.time,
                         doc="Gas-solids pressure equilibrium constraint")
        def solids_pressure_constraint(b, t):
            return b.solids[t].pressure == \
                b.gas_phase.properties_out[t].pressure

        # Add Ports for gas phase
        self.add_inlet_port("inlet", self.gas_phase)
        self.add_outlet_port("outlet", self.gas_phase)

    def initialize(blk, state_args=None, outlvl=6,
                   solver='ipopt', optarg={'tol': 1e-6}):
        '''
        This is a general purpose initialization routine for simple unit
        models. This method assumes a single ControlVolume block called
        controlVolume, and first initializes this and then attempts to solve
        the entire unit.

        More complex models should overload this method with their own
        initialization routines,

        Keyword Arguments:
            state_args : a dict of arguments to be passed to the property
                           package(s) to provide an initial state for
                           initialization (see documentation of the specific
                           property package) (default = {}).
            outlvl : sets output level of initialization routine
                 * 0 = Use default idaes.init logger setting
                 * 1 = Maximum output
                 * 2 = Include solver output
                 * 3 = Return solver state for each step in subroutines
                 * 4 = Return solver state for each step in routine
                 * 5 = Final initialization status and exceptions
                 * 6 = No output
            optarg : solver options dictionary object (default={'tol': 1e-6})
            solver : str indicating which solver to use during
                     initialization (default = 'ipopt')

        Returns:
            None
        '''
        # Set solver options
        init_log = getInitLogger(blk.name, outlvl)
        solve_log = getSolveLogger(blk.name, outlvl, tag="unit")
        opt = get_solver(solver, optarg) # create solver
        opt.options = optarg

        # ---------------------------------------------------------------------
        # Initialize control volume block
        flags = blk.gas_phase.initialize(outlvl=outlvl+1,
                                         optarg=optarg,
                                         solver=solver,
                                         state_args=state_args)

        init_log.log(4, 'Initialization Step 1 Complete.')

        # ---------------------------------------------------------------------
        # Solve unit
        with solver_log(solve_log, DEBUG) as slc:
            results = opt.solve(blk, tee=slc.tee)
#        results = opt.solve(blk, tee=init_tee(init_log))

        init_log.log(4, "Initialization Step 2 {}."
                     .format(condition(results)))

        # ---------------------------------------------------------------------
        # Release Inlet state
        blk.gas_phase.release_state(flags, outlvl+1)

        init_log.log(5, 'Initialization Complete: {}'
                     .format(condition(results)))
