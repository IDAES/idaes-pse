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
Reboiler model for solvent columns.

This is a simple model for a rboiler in the case where liquid and vapor phases
have separate proeprty packages, suchas the case of solvent columns.

Assumptions:
     * Liquid phase property pakcage has a single phase named Liq
     * Vapor phase property pakcagehas a single phase named Vap
     * Liquid and vapor phase proeprtes need not have the same component lists
"""

__author__ = "andrew Lee, Paul Akula"

# Import Pyomo libraries
from pyomo.environ import Constraint, Param
from pyomo.common.config import ConfigBlock, ConfigValue, In, Bool

# Import IDAES cores
import idaes.logger as idaeslog
from idaes.core import (ControlVolume0DBlock,
                        declare_process_block_class,
                        EnergyBalanceType,
                        MomentumBalanceType,
                        MaterialBalanceType,
                        UnitModelBlockData,
                        useDefault)
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import PropertyPackageError, \
    PropertyNotSupportedError, ConfigurationError
from idaes.core.util import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom


_log = idaeslog.getIdaesLogger(__name__)


@declare_process_block_class("SolventReboiler")
class SolventReboilerData(UnitModelBlockData):
    """
    Reboiler unit for solvent column models using separate property packages
    for liquid and vpor phases.

    Unit model to reboil the liquid from the bottom tray of
    the distillation column.
    """
    CONFIG = ConfigBlock()
    # TOOO: Add dynamics in future
    CONFIG.declare("dynamic", ConfigValue(
        domain=In([False]),
        default=False,
        description="Dynamic model flag - must be False",
        doc="""Indicates whether this model will be dynamic or not,
**default** = False. Equilibrium Reactors do not support dynamic behavior."""))
    CONFIG.declare("has_holdup", ConfigValue(
        default=False,
        domain=In([False]),
        description="Holdup construction flag - must be False",
        doc="""Indicates whether holdup terms should be constructed or not.
**default** - False. Equilibrium reactors do not have defined volume, thus
this must be False."""))
    # TODO :Add boilup ratio back later if needed
#     CONFIG.declare("has_boilup_ratio", ConfigValue(
#         default=False,
#         domain=Bool,
#         description="Boilup ratio term construction flag",
#         doc="""Indicates whether terms for boilup ratio should be
# constructed,
# **default** - False.
# **Valid values:** {
# **True** - include construction of boilup ratio constraint,
# **False** - exclude construction of boilup ratio constraint}"""))
    CONFIG.declare("material_balance_type", ConfigValue(
        default=MaterialBalanceType.useDefault,
        domain=In(MaterialBalanceType),
        description="Material balance construction flag",
        doc="""Indicates what type of mass balance should be constructed,
**default** - MaterialBalanceType.componentPhase.
**Valid values:** {
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
**default** - EnergyBalanceType.enthalpyTotal.
**Valid values:** {
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
    CONFIG.declare("has_pressure_change", ConfigValue(
        default=False,
        domain=Bool,
        description="Pressure change term construction flag",
        doc="""Indicates whether terms for pressure change should be
constructed,
**default** - False.
**Valid values:** {
**True** - include pressure change terms,
**False** - exclude pressure change terms.}"""))
    CONFIG.declare("liquid_property_package", ConfigValue(
        default=useDefault,
        domain=is_physical_parameter_block,
        description="Property package to use for liquid phase",
        doc="""Property parameter object used to define property calculations
for the liquid phase,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PropertyParameterObject** - a PropertyParameterBlock object.}"""))
    CONFIG.declare("liquid_property_package_args", ConfigBlock(
        implicit=True,
        description="Arguments to use for constructing liquid phase properties",
        doc="""A ConfigBlock with arguments to be passed to liquid phase
property block(s) and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}"""))
    CONFIG.declare("vapor_property_package", ConfigValue(
        default=useDefault,
        domain=is_physical_parameter_block,
        description="Property package to use for vapor phase",
        doc="""Property parameter object used to define property calculations
for the vapor phase,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PropertyParameterObject** - a PropertyParameterBlock object.}"""))
    CONFIG.declare("vapor_property_package_args", ConfigBlock(
        implicit=True,
        description="Arguments to use for constructing vapor phase properties",
        doc="""A ConfigBlock with arguments to be passed to vapor phase
property block(s) and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}"""))

    def build(self):
        """Build the model.

        Args:
            None
        Returns:
            None
        """
        # Call UnitModel.build to setup dynamics
        super().build()

        # TODO: Add checks for compatability of property packages
        # Flow basis, phase naming assumptions

        # ---------------------------------------------------------------------
        # Add Control Volume for the Liquid Phase
        self.liquid_phase = ControlVolume0DBlock(default={
            "dynamic": self.config.dynamic,
            "has_holdup": self.config.has_holdup,
            "property_package": self.config.liquid_property_package,
            "property_package_args": self.config.liquid_property_package_args})

        self.liquid_phase.add_state_blocks(
            has_phase_equilibrium=True)

        # Separate liquid and vapor phases measn that phase equilibrium will
        # be handled at the unit model level, thus has_phase_equilibrium is
        # Flase, but has_mass_transfer is True.
        self.liquid_phase.add_material_balances(
            balance_type=self.config.material_balance_type,
            has_mass_transfer=True,
            has_phase_equilibrium=False)

        # Need to include enthalpy transfer term for the mass transfer
        self.liquid_phase.add_energy_balances(
            balance_type=self.config.energy_balance_type,
            has_heat_transfer=True,
            has_enthalpy_transfer=True)

        self.liquid_phase.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=self.config.has_pressure_change)

        # ---------------------------------------------------------------------
        # Add single state block for vapor phase
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["defined_state"] = False
        self.vapor_phase = self.config.vapor_property_package.buid_state_block(
                self.flowsheet().time,
                doc="Vapor phase properties",
                default=tmp_dict)

        # ---------------------------------------------------------------------
        # Add Ports for the reboiler
        self.add_inlet_port()
        self.add_outlet_port(name="bottoms", doc="Bottoms stream")
        self.add_outlet_port(name="vapor_reboil",
                             block=self.vapor_phase,
                             doc="Vapor stream from reboiler")

        # ---------------------------------------------------------------------
        # Add unit level constraints
        # First, need the union and intersection of component lists
        all_comps = (self.vapor_phase.component_list |
                     self.liquid_phase_properties_out.component_list)
        common_comps = (self.vapor_phase.component_list &
                        self.liquid_phase_properties_out.component_list)

        if any(j not in common_comps for j in self.vapor_phase.component_list):
            # We have non-condensable components present, need zero-flow param
            units = self.config.vapor_property_package.get_metadata(
                ).get_derived_units
            self.zero_flow_param = Param(
                mutable=True,
                default=1e-8,
                units=units["flow_mole"])

        # Material balances
        # TODO : Need to make sure flow bases are the same for both prop packs
        # TODO: Need unit conversions
        def rule_material_balance(blk, t, j):
            if j in common_comps:
                # Component is in equilibrium
                # Mass transfer equals vapor flowrate
                return (blk.liquid_phase.mass_transfer_term[t, "Liq", j] ==
                        blk.vapor_phase[t].get_material_flow_term("Vap", j))
            elif j in self.vapor_phase.component_list:
                # Non-condensable component
                # No mass transfer term
                # Set vapor flowrate to an arbitary small value
                return (blk.vapor_phase[t].get_material_flow_term("Vap", j) ==
                        blk.zero_flow_param)
            else:
                # Non-vaporisable comonent
                # Mass transfer term is zero, no vapor flowrate
                return blk.liquid_phase.mass_transfer_term[t, "Liq", j] == 0
        self.unit_mass_balance = Constraint(
            self.flowsheet().time,
            all_comps,
            rule=rule_material_balance,
            doc="Unit level material balances")

        # Phase equilibrium constraints
        # For all common components, equate fugacity in vapor and liquid
        def rule_phase_equilibrium(blk, t, j):
            return (
                blk.liquid_phase.properties_out[t].log_fugacity_phase_comp[
                    "Liq", j] ==
                blk.vapor_phase[t].log_fugacity_phase_comp["Vap", j])
        self.unit_phase_equilibrium = Constraint(
            self.flowsheet().time,
            common_comps,
            rule=rule_phase_equilibrium,
            doc="Unit levelphase equilibrium constraints")

        # Temperature equality constraint
        def rule_temperature_balance(blk, t):
            return (blk.liquid_phase.properties_out[t].temperature ==
                    blk.vapor_phase[t].temperature)
        self.unit_temperature_equality = Constraint(
            self.flowsheet().time,
            rule=rule_temperature_balance,
            doc="Unit level temperature equality")

        # Unit level energy balance
        # Energy leaving in vapor phase must be equal and opposite to enthalpy
        # transfer from liquid phase
        # TODO: How does this need to account for dynamics?
        def rule_energy_balance(blk, t):
            return (blk.liquid_phase.enthalpy_transfer_term[t] ==
                    blk.vapor_phase[t].get_enthalpy_flow_term("Vap"))
        self.unit_enthalpy_balance = Constraint(
            self.flowsheet().time,
            rule=rule_energy_balance,
            doc="Unit level enthalpy_balance")

        # Pressure balance constraint
        def rule_pressure_balance(blk, t):
            return (blk.liquid_phase.properties_out[t].pressure ==
                    blk.vapor_phase[t].pressure)
        self.unit_pressure_balance = Constraint(
            self.flowsheet().time,
            rule=rule_pressure_balance,
            doc="Unit level pressure balance")

        # TODO: Unit level references (heat duty, deltaP)

    # def initialize(self, state_args=None, solver=None, optarg=None,
    #                outlvl=idaeslog.NOTSET):

    #     init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
    #     solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

    #     solverobj = get_solver(solver, optarg)

    #     # Initialize the inlet and outlet state blocks. Calling the state
    #     # blocks initialize methods directly so that custom set of state args
    #     # can be passed to the inlet and outlet state blocks as control_volume
    #     # initialize method initializes the state blocks with the same
    #     # state conditions.
    #     flags = self.control_volume.properties_in. \
    #         initialize(state_args=state_args,
    #                    solver=solver,
    #                    optarg=optarg,
    #                    outlvl=outlvl,
    #                    hold_state=True)

    #     # Initialize outlet state block at same conditions of inlet except
    #     # the temperature. Set the temperature to a temperature guess based
    #     # on the desired boilup_ratio.

    #     # Get index for bubble point temperature and and assume it
    #     # will have only a single phase equilibrium pair. This is to
    #     # support the generic property framework where the T_bubble
    #     # is indexed by the phases_in_equilibrium. In distillation,
    #     # the assumption is that there will only be a single pair
    #     # i.e. vap-liq. 
    #     idx = next(iter(self.control_volume.properties_in[0].
    #                     temperature_bubble))
    #     temp_guess = 0.5 * (
    #         value(self.control_volume.properties_in[0].temperature_dew[idx]) -
    #         value(self.control_volume.properties_in[0].
    #               temperature_bubble[idx])) + \
    #         value(self.control_volume.properties_in[0].temperature_bubble[idx])

    #     state_args_outlet = {}
    #     state_dict_outlet = (
    #         self.control_volume.properties_in[
    #             self.flowsheet().time.first()]
    #         .define_port_members())

    #     for k in state_dict_outlet.keys():
    #         if state_dict_outlet[k].is_indexed():
    #             state_args_outlet[k] = {}
    #             for m in state_dict_outlet[k].keys():
    #                 state_args_outlet[k][m] = value(state_dict_outlet[k][m])
    #         else:
    #             if k != "temperature":
    #                 state_args_outlet[k] = value(state_dict_outlet[k])
    #             else:
    #                 state_args_outlet[k] = temp_guess

    #     self.control_volume.properties_out.initialize(
    #         state_args=state_args_outlet,
    #         solver=solver,
    #         optarg=optarg,
    #         outlvl=outlvl,
    #         hold_state=False)

    #     if degrees_of_freedom(self) == 0:
    #         with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
    #             res = solverobj.solve(self, tee=slc.tee)
    #         init_log.info(
    #             "Initialization Complete, {}.".format(idaeslog.condition(res))
    #         )
    #     else:
    #         raise ConfigurationError(
    #             "State vars fixed but degrees of freedom "
    #             "for reboiler is not zero during "
    #             "initialization. Please ensure that the boilup_ratio "
    #             "or the outlet temperature is fixed.")

    #     self.control_volume.properties_in.\
    #         release_state(flags=flags, outlvl=outlvl)

    # def _get_performance_contents(self, time_point=0):
    #     var_dict = {}
    #     if hasattr(self, "heat_duty"):
    #         var_dict["Heat Duty"] = self.heat_duty[time_point]

    #     return {"vars": var_dict}

    # def _get_stream_table_contents(self, time_point=0):
    #     stream_attributes = {}

    #     stream_dict = {"Inlet": "inlet",
    #                    "Vapor Reboil": "vapor_reboil",
    #                    "Bottoms": "bottoms"}

    #     for n, v in stream_dict.items():
    #         port_obj = getattr(self, v)

    #         stream_attributes[n] = {}

    #         for k in port_obj.vars:
    #             for i in port_obj.vars[k].keys():
    #                 if isinstance(i, float):
    #                     stream_attributes[n][k] = value(
    #                         port_obj.vars[k][time_point])
    #                 else:
    #                     if len(i) == 2:
    #                         kname = str(i[1])
    #                     else:
    #                         kname = str(i[1:])
    #                     stream_attributes[n][k + " " + kname] = \
    #                         value(port_obj.vars[k][time_point, i[1:]])

    #     return DataFrame.from_dict(stream_attributes, orient="columns")
