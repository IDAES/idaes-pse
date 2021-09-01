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
from pyomo.environ import Constraint, Var
from pyomo.dae import DerivativeVar
from pyomo.common.config import ConfigBlock, ConfigValue, In

# Import IDAES cores
from idaes.core import (declare_process_block_class,
                        EnergyBalanceType,
                        UnitModelBlockData,
                        useDefault)
from idaes.core.util.config import (is_physical_parameter_block)
import idaes.logger as idaeslog
from idaes.core.util import get_solver

__author__ = "Chinedu Okoli, Andrew Lee"


# Assumptions:
# Only Vap and Sol phases, which are explicitly named as such
# Perfect mixing in Vap phase
# Static solid phase

# Need to have build-on-demand properties in reactions for this to work

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
**default** = True. Fixed beds must be dynamic."""))
    CONFIG.declare("has_holdup", ConfigValue(
        default=True,
        domain=In([True]),
        description="Holdup construction flag - must be True",
        doc="""Indicates whether holdup terms should be constructed or not.
**default** - True. Fixed bed reactors must be dynamic, thus this must be
True."""))
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
        # TODO: Had to remove domain due to limitations in ReactionBase for
        # heterogeneous systems
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

        # Build Gas Phase StateBlock
        # This block only needed so gas conc. to rxn block is calculated
        # Defined state is set to True so that the "sum(mole_frac)=1" eqn in
        # the gas state block is deactivated.
        self.gas = self.config.gas_property_package.state_block_class(
                self.flowsheet().time,
                default={"parameters": self.config.gas_property_package,
                         "defined_state": True})

        # Build Solid Phase StateBlock
        # As there is no solid flow, only need a single state block
        # Defined state is set to True so that the "sum(mass_frac)=1" eqn in
        # the solid state block is deactivated. This is done here as there is
        # currently no way to deactivate the constraint at the initial time
        # for batch systems (i.e. no inlet or outlet ports).
        # The "sum(mass_frac)=1 for all t neq 0" eqn is written in the
        # unit model instead
        self.solids = self.config.solid_property_package.state_block_class(
                self.flowsheet().time,
                default={"parameters": self.config.solid_property_package,
                         "defined_state": True})

        tmp_dict = dict(**self.config.reaction_package_args)
        tmp_dict["gas_state_block"] = self.gas
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
        self.volume_solid = Var(initialize=1.0,
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
                b.volume_solid *
                b.solids[t].get_material_density_terms("Sol", j))

        self.solids_material_accumulation = DerivativeVar(
                    self.solids_material_holdup,
                    wrt=self.flowsheet().config.time,
                    doc="Solids material accumulation")

        @self.Constraint(self.flowsheet().config.time,
                         self.config.solid_property_package.component_list,
                         doc="Solid phase material accumulation constraints")
        def solids_material_accumulation_constraints(b, t, j):
            # if t == self.flowsheet().config.time.first():
            #     return Constraint.Skip
            # else:
            return (b.solids_material_accumulation[t, j] ==
                    b.volume_solid *
                    b.solids[t]._params.mw_comp[j] *
                    sum(b.reactions[t].reaction_rate[r] *
                        b.config.reaction_package
                        .rate_reaction_stoichiometry[r, "Sol", j] for r in
                        b.config.reaction_package.rate_reaction_idx))

        # Add solid mass variable and constraint for TGA tracking
        self.mass_solids = Var(self.flowsheet().config.time,
                               doc="Total mass of solids")

        @self.Constraint(self.flowsheet().config.time,
                         doc="Calculating total mass of solids")
        def mass_solids_eq(b, t):
            return 1e2*b.mass_solids[t] == 1e2*b.volume_solid * sum(
                    b.solids[t].get_material_density_terms("Sol", j)
                    for j in b.config.solid_property_package.component_list)

        # Sum of mass fractions at all time equals 1
        @self.Constraint(self.flowsheet().config.time,
                         doc="Sum of mass fractions at all time")
        def sum_component_eqn(b, t):
            if t == b.flowsheet().config.time.first():
                return Constraint.Skip
            else:
                return 1e2 == 1e2*sum(b.solids[t].mass_frac_comp[j]
                                      for j in b.solids[t].
                                      _params.component_list)

        if self.config.energy_balance_type != EnergyBalanceType.none:
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
                      self.volume_solid *
                      b.solids[t].get_energy_density_terms("Sol"))

            self.solids_energy_accumulation = DerivativeVar(
                        self.solids_energy_holdup,
                        doc="Solids energy accumulation")

            @self.Constraint(self.flowsheet().config.time,
                             doc="Solid phase energy accumulation constraints")
            def solids_energy_accumulation_constraints(b, t):
                # if t == self.flowsheet().config.time.first():
                #     return Constraint.Skip
                # else:
                return b.solids_energy_accumulation[t] == \
                    - sum(b.reactions[t].reaction_rate[r] *
                          b.volume_solid *
                          b.reactions[t].dh_rxn[r]
                          for r in b.config.reaction_package.
                          rate_reaction_idx)

    def initialize(blk, state_args=None, outlvl=idaeslog.NOTSET,
                   solver=None, optarg=None):
        """
        This is a general purpose initialization routine for simple unit
        models. This method assumes a single ControlVolume block called
        controlVolume, and first initializes this and then attempts to solve
        the entire unit.
        More complex models should overload this method with their own
        initialization routines.

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
        """

        # Set solver options
        init_log = idaeslog.getInitLogger(blk.name, outlvl)
        opt = get_solver(solver, optarg)  # create solver
        opt.options = optarg

        # ---------------------------------------------------------------------
        # Initialize control volume block
        init_log.info('Initialize Thermophysical Properties')

        # Initialize gas_phase block
        gas_phase_flags = blk.gas.initialize(
                                outlvl=outlvl,
                                optarg=optarg,
                                solver=solver)

        # Initialize solid_phase properties block
        solid_phase_flags = blk.solids.initialize(
                                outlvl=outlvl,
                                optarg=optarg,
                                solver=solver)

        print()
        print('Initialize reaction properties')
        # Initialize reactions
        blk.reactions.initialize(outlvl=outlvl,
                                 optarg=optarg,
                                 solver=solver)

        # TODO - maybe set this up in a nicer way
        # Fix the solids temperature to initial value if isothermal conditions
        if blk.config.energy_balance_type == EnergyBalanceType.none:
            for t in blk.flowsheet().config.time:
                blk.solids[t].temperature.fix(blk.solids[0].temperature.value)

        init_log.info_high("Initialization Step 1 Complete.")

        # ---------------------------------------------------------------------
        # Release Inlet state
        blk.gas.release_state(gas_phase_flags, outlvl+1)
        blk.solids.release_state(solid_phase_flags, outlvl+1)
