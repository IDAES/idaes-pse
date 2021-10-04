##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
##############################################################################
# DISPATCHES was produced under the DOE Design Integration and Synthesis
# Platform to Advance Tightly Coupled Hybrid Energy Systems program (DISPATCHES),
# and is copyright (c) 2021 by the software owners: The Regents of the University
# of California, through Lawrence Berkeley National Laboratory, National
# Technology & Engineering Solutions of Sandia, LLC, Alliance for Sustainable
# Energy, LLC, Battelle Energy Alliance, LLC, University of Notre Dame du Lac, et
# al. All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. Both files are also available online at the URL:
# "https://github.com/gmlc-dispatches/dispatches".
#
##############################################################################

"""
Basic IDAES unit model for compressed gas tank.
Uses custom material and energy balance equations.
Suitable for use with steady state flowsheets indexed with time.

A previous state block (previous_state) is defined
to pass the previous time point or starting condition of the tank,
 i.e., Pressure and Temperature.
For dynamic operations, multiple instances of this model indexed with time
can be used. For this, equality constraints must be written at
flowsheet level, i.e.,
the final state of tank (Pres & Temp) at time (t-1) ==
    intial state of tank (Pres & Temp) at time (t)

Adiabatic operations are assumed.
TODO: add constraints to enable isothermal operations
TODO: add capability to allow liquid phase storage
"""

# Import Pyomo libraries
from pyomo.environ import (Var,
                           Reals,
                           NonNegativeReals,
                           Constraint,
                           )
# from pyomo.network import Port
from pyomo.common.config import ConfigBlock, ConfigValue, In

# Import IDAES modules
from idaes.core import (ControlVolume0DBlock,
                        declare_process_block_class,
                        MaterialFlowBasis,
                        MomentumBalanceType,
                        UnitModelBlockData,
                        useDefault)
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.constants import Constants as const
from idaes.core.util import get_solver
import idaes.logger as idaeslog
import idaes.core.util.scaling as iscale

__author__ = "Naresh Susarla"

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("CompressedGasTank",
                             doc="Simple compressed gas tank model")
class CompressedGasTankData(UnitModelBlockData):
    """
    Simple gas tank model.
    Unit model to store or supply compressed gas.

    """
    CONFIG = ConfigBlock()

    # This model is based on steady state material & energy balances.
    # The accumulation term is computed based on the tank state at
    # previous time step. Thus, dynamic option is turned off.
    # However, a dynamic analysis can be performed by creating
    # an instance of this model for every time step.
    CONFIG.declare("dynamic", ConfigValue(
        domain=In([False]),
        default=False,
        description="Dynamic model flag - must be False",
        doc="""Indicats if Gas tank model is dynamic,
**default** = False. Equilibrium Reactors do not support dynamic behavior."""))
    CONFIG.declare("has_holdup", ConfigValue(
        default=False,
        domain=In([False]),
        description="Holdup construction flag - must be False",
        doc="""Indicates whether holdup terms should be constructed or not.
**default** - False. Gas tank model uses custom equations for holdup."""))
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
    CONFIG.declare("property_package", ConfigValue(
        default=useDefault,
        domain=is_physical_parameter_block,
        description="Property package to use for control volume",
        doc="""Property parameter object used to define property calculations,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PhysicalParameterObject** - a PhysicalParameterBlock object.}"""))
    CONFIG.declare("property_package_args", ConfigBlock(
        implicit=True,
        description="Arguments to use for constructing property packages",
        doc="""A ConfigBlock with arguments to be passed to a property block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}"""))

    def build(self):
        """Building model
        Args:
            None
        Returns:
            None
        """
        super().build()

        # Build Control Volume
        self.control_volume = ControlVolume0DBlock(default={
            "dynamic": self.config.dynamic,
            "property_package": self.config.property_package,
            "property_package_args": self.config.property_package_args})

        # add inlet and outlet states
        self.control_volume.add_state_blocks(has_phase_equilibrium=False)

        # add tank volume
        self.control_volume.add_geometry()

        # add phase fractions
        self.control_volume._add_phase_fractions()

        # add pressure balance
        self.control_volume.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=True)

        # add a state block 'previous_state' for the storage tank
        # this state block is needed to compute material and energy holdup
        # at previous or the starting time step using the property package
        # for a given P_prev and T_prev
        # NOTE: there is no flow in the previous state so,
        # flow_mol state variable is fixed to 0
        self.previous_state = (
            self.config.property_package.build_state_block(
                self.flowsheet().config.time,
                doc="Tank state at previous time"))

        # previous state should not have any flow
        self.previous_state[:].flow_mol.fix(0)

        # add local lists for easy use
        phase_list = self.control_volume.properties_in.phase_list
        pc_set = self.control_volume.properties_in.phase_component_set
        component_list = self.control_volume.properties_in.component_list

        # Get units from property package
        units = self.config.property_package.\
            get_metadata().get_derived_units
        if (self.control_volume.properties_in[
                self.flowsheet().config.time.first()]
                .get_material_flow_basis() == MaterialFlowBasis.molar):
            flow_units = units("flow_mole")
            material_units = units("amount")
        elif (self.control_volume.properties_in[
                self.flowsheet().config.time.first()]
              .get_material_flow_basis() == MaterialFlowBasis.mass):
            flow_units = units("flow_mass")
            material_units = units("mass")

        # Add Inlet and Outlet Ports
        self.add_inlet_port()
        self.add_outlet_port()

        # Define Vars for Tank volume calculations
        self.tank_diameter = Var(
            self.flowsheet().config.time,
            within=NonNegativeReals,
            initialize=1.0,
            doc="Diameter of storage tank in m",
            units=units("length"))

        self.tank_length = Var(
            self.flowsheet().config.time,
            within=NonNegativeReals,
            initialize=1.0,
            doc="Length of storage tank in m",
            units=units("length"))

        # Tank volume calculation
        @self.Constraint(self.flowsheet().config.time)
        def volume_cons(b, t):
            return (b.control_volume.volume[t] ==
                    const.pi * b.tank_length[t] *
                    ((b.tank_diameter[t]/2)**2))

        # define Vars for the model
        self.dt = Var(
            self.flowsheet().config.time,
            domain=NonNegativeReals,
            initialize=100,
            doc="Time step for holdup calculation",
            units=units("time"))

        self.heat_duty = Var(
            self.flowsheet().config.time,
            domain=Reals,
            initialize=0.0,
            doc="Heat transferred from surroundings, 0 for adiabatic",
            units=units("power"))

        self.material_accumulation = Var(
            self.flowsheet().config.time,
            pc_set,
            within=Reals,
            initialize=1.0,
            doc="Accumulation of material in tank",
            units=flow_units)

        self.energy_accumulation = Var(
            self.flowsheet().config.time,
            phase_list,
            within=Reals,
            initialize=1.0,
            doc="Energy accumulation",
            units=units("power"))

        self.material_holdup = Var(
            self.flowsheet().config.time,
            pc_set,
            within=Reals,
            initialize=1.0,
            doc="Material holdup in tank",
            units=material_units)

        self.energy_holdup = Var(
            self.flowsheet().config.time,
            phase_list,
            within=Reals,
            initialize=1.0,
            doc="Energy holdup in tank",
            units=units("energy"))

        self.previous_material_holdup = Var(
            self.flowsheet().config.time,
            pc_set,
            within=Reals,
            initialize=1.0,
            doc="Tank material holdup at previous time",
            units=material_units)

        self.previous_energy_holdup = Var(
            self.flowsheet().config.time,
            phase_list,
            within=Reals,
            initialize=1.0,
            doc="Tank energy holdup at previous time",
            units=units("energy"))

        # Adiabatic operations are assumed
        # Fixing the heat_duty to 0 here to avoid any misakes at use
        # TODO: remove this once the isothermal constraints are added
        self.heat_duty.fix(0)

        # Computing material and energy holdup in the tank at previous time
        # using previous state Pressure and Temperature of the tank
        @self.Constraint(self.flowsheet().config.time,
                         pc_set,
                         doc="Material holdup at previous time")
        def previous_material_holdup_rule(b, t, p, j):
            return (
                b.previous_material_holdup[t, p, j]
                == b.control_volume.volume[t] *
                b.control_volume.phase_fraction[t, p] *
                b.previous_state[t].get_material_density_terms(p, j)
                )

        @self.Constraint(self.flowsheet().config.time,
                         phase_list,
                         doc="Energy holdup at previous time")
        def previous_energy_holdup_rule(b, t, p):
            if (self.control_volume.properties_in[t]
                    .get_material_flow_basis() == MaterialFlowBasis.molar):
                return (
                    b.previous_energy_holdup[t, p] == (
                        sum(b.previous_material_holdup[t, p, j]
                            for j in component_list) *
                        b.previous_state[t].energy_internal_mol_phase[p])
                    )
            if (self.control_volume.properties_in[t]
                    .get_material_flow_basis() == MaterialFlowBasis.mass):
                return (
                    b.previous_energy_holdup[t, p] == (
                        sum(b.previous_material_holdup[t, p, j]
                            for j in component_list) *
                        (b.previous_state[t].energy_internal_mol_phase[p] /
                         b.previous_state[t].mw))
                    )

        # component material balances
        @self.Constraint(self.flowsheet().config.time,
                         pc_set,
                         doc="Material balances")
        def material_balances(b, t, p, j):
            if (p, j) in pc_set:
                return (
                    b.material_accumulation[t, p, j] == (
                        b.control_volume.properties_in[t].
                        get_material_flow_terms(p, j) -
                        b.control_volume.properties_out[t].
                        get_material_flow_terms(p, j))
                    )
            else:
                return Constraint.Skip

        # integration of material accumulation
        @self.Constraint(self.flowsheet().config.time,
                         pc_set,
                         doc="Material holdup integration")
        def material_holdup_integration(b, t, p, j):
            if (p, j) in pc_set:
                return b.material_holdup[t, p, j] == (
                      b.dt[t]
                      * b.material_accumulation[t, p, j]
                      + b.previous_material_holdup[t, p, j])

        # material holdup calculation
        @self.Constraint(self.flowsheet().config.time,
                         pc_set,
                         doc="Material holdup calculations")
        def material_holdup_calculation(b, t, p, j):
            if (p, j) in pc_set:
                return (
                    b.material_holdup[t, p, j] == (
                        b.control_volume.volume[t] *
                        b.control_volume.phase_fraction[t, p] *
                        b.control_volume.properties_out[t].
                        get_material_density_terms(p, j)))

        # energy accumulation
        @self.Constraint(self.flowsheet().config.time,
                         doc="Energy accumulation")
        def energy_accumulation_equation(b, t):
            return (
                sum(b.energy_accumulation[t, p] for p in phase_list)
                * b.dt[t] == sum(b.energy_holdup[t, p] for p in phase_list) -
                sum(b.previous_energy_holdup[t, p] for p in phase_list)
                )

        # energy holdup calculation
        @self.Constraint(self.flowsheet().config.time,
                         phase_list,
                         doc="Energy holdup calculation")
        def energy_holdup_calculation(b, t, p):
            if (self.control_volume.properties_in[t].
                    get_material_flow_basis() == MaterialFlowBasis.molar):
                return (
                    b.energy_holdup[t, p] == (
                        sum(b.material_holdup[t, p, j] for j in component_list)
                        * b.control_volume.properties_out[t].
                        energy_internal_mol_phase[p])
                    )
            if (self.control_volume.properties_in[t].
                    get_material_flow_basis() == MaterialFlowBasis.mass):
                return (
                    b.energy_holdup[t, p] == (
                        sum(b.material_holdup[t, p, j] for j in component_list)
                        * (b.control_volume.properties_out[t].
                           energy_internal_mol_phase[p] /
                           b.control_volume.properties_out[t].mw))
                    )

        # Energy balance based on internal energy, as follows:
        #     n_final * U_final =
        #                        n_previous * U_previous +
        #                        n_inlet * H_inlet - n_outlet * H_outlet
        #     where, n is number of moles, U is internal energy, H is enthalpy
        @self.Constraint(self.flowsheet().config.time, doc="Energy balance")
        def energy_balances(b, t):
            return (
                    sum(b.energy_holdup[t, p] for p in phase_list) ==
                    sum(b.previous_energy_holdup[t, p] for p in phase_list) +
                    b.dt[t] *
                    (sum(b.control_volume.properties_in[t].
                         get_enthalpy_flow_terms(p) for p in phase_list) -
                     sum(b.control_volume.properties_out[t].
                         get_enthalpy_flow_terms(p) for p in phase_list))
                    )

    def initialize(blk, state_args=None, outlvl=idaeslog.NOTSET,
                   solver=None, optarg=None):
        '''
        Gas tank model initialization routine.

        Keyword Arguments:
            state_args : a dict of arguments to be passed to the property
                           package(s) for the control_volume of the model to
                           provide an initial state for initialization
                           (see documentation of the specific property package)
                           (default = None).
            outlvl : sets output level of initialisation routine

                     * 0 = no output (default)
                     * 1 = return solver state for each step in routine
                     * 2 = return solver state for each step in subroutines
                     * 3 = include solver output infomation (tee=True)

            optarg : solver options dictionary object (default={'tol': 1e-6})
            solver : str indicating whcih solver to use during
                     initialization (default = 'ipopt')

        Returns:
            None
        '''

        if state_args is None:
            state_args = dict()

        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")

        opt = get_solver(solver, optarg)

        init_log.info_low("Starting initialization...")

        flags = blk.control_volume.initialize(state_args=state_args,
                                              outlvl=outlvl,
                                              optarg=optarg,
                                              solver=solver)
        flag_previous_state = blk.previous_state.initialize(
                outlvl=outlvl,
                optarg=optarg,
                solver=solver,
                hold_state=True,
                state_args=state_args,
        )
        blk.previous_state[0].sum_mole_frac_out.deactivate()

        init_log.info_high("Initialization Step 1 Complete.")

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high(
                "Initialization Step 2 {}.".format(idaeslog.condition(res))
            )
        blk.previous_state[0].sum_mole_frac_out.activate()
        blk.control_volume.release_state(flags, outlvl)
        blk.previous_state.release_state(flag_previous_state, outlvl)
        init_log.info("Initialization Complete.")

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        if hasattr(self, "previous_state"):
            for t, v in self.previous_state.items():
                iscale.set_scaling_factor(v.flow_mol, 1e-3)
                iscale.set_scaling_factor(v.pressure, 1e-5)
                iscale.set_scaling_factor(v.temperature, 1e-1)

        if hasattr(self, "tank_diameter"):
            for t, v in self.tank_diameter.items():
                iscale.set_scaling_factor(v, 1)

        if hasattr(self, "tank_length"):
            for t, v in self.tank_length.items():
                iscale.set_scaling_factor(v, 1)

        if hasattr(self, "heat_duty"):
            for t, v in self.heat_duty.items():
                iscale.set_scaling_factor(v, 1e-5)

        if hasattr(self, "material_accumulation"):
            for (t, p, j), v in self.material_accumulation.items():
                iscale.set_scaling_factor(v, 1e-3)

        if hasattr(self, "energy_accumulation"):
            for (t, p), v in self.energy_accumulation.items():
                iscale.set_scaling_factor(v, 1e-3)

        if hasattr(self, "material_holdup"):
            for (t, p, j), v in self.material_holdup.items():
                iscale.set_scaling_factor(v, 1e-5)

        if hasattr(self, "energy_holdup"):
            for (t, p), v in self.energy_holdup.items():
                iscale.set_scaling_factor(v, 1e-5)

        if hasattr(self, "previous_material_holdup"):
            for (t, p, j), v in self.previous_material_holdup.items():
                iscale.set_scaling_factor(v, 1e-5)

        if hasattr(self, "previous_energy_holdup"):
            for (t, p), v in self.previous_energy_holdup.items():
                iscale.set_scaling_factor(v, 1e-5)

        # Volume constraint
        if hasattr(self, "volume_cons"):
            for t, c in self.volume_cons.items():
                iscale.constraint_scaling_transform(
                    c, iscale.get_scaling_factor(
                        self.tank_length[t],
                        default=1, warning=True))

        # Previous time Material Holdup Rule
        if hasattr(self, "previous_material_holdup_rule"):
            for (t, i, j), c in self.previous_material_holdup_rule.items():
                iscale.constraint_scaling_transform(
                    c, iscale.get_scaling_factor(
                        self.material_holdup[t, i, j],
                        default=1, warning=True))

        # Previous time Energy Holdup Rule
        if hasattr(self, "previous_energy_holdup_rule"):
            for (t, i), c in self.previous_energy_holdup_rule.items():
                iscale.constraint_scaling_transform(
                    c, iscale.get_scaling_factor(
                        self.energy_holdup[t, i],
                        default=1, warning=True))

        # Material Balances
        if hasattr(self, "material_balances"):
            for (t, i, j), c in self.material_balances.items():
                iscale.constraint_scaling_transform(
                    c, iscale.get_scaling_factor(
                        self.material_accumulation[t, i, j],
                        default=1, warning=True))

        # Material Holdup Integration
        if hasattr(self, "material_holdup_integration"):
            for (t, i, j), c in self.material_holdup_integration.items():
                iscale.constraint_scaling_transform(
                    c, iscale.get_scaling_factor(
                        self.material_holdup[t, i, j],
                        default=1, warning=True))

        # Material Holdup Constraints
        if hasattr(self, "material_holdup_calculation"):
            for (t, i, j), c in self.material_holdup_calculation.items():
                iscale.constraint_scaling_transform(
                    c, iscale.get_scaling_factor(
                        self.material_holdup[t, i, j],
                        default=1, warning=True))

        # Enthalpy Balances
        if hasattr(self, "energy_accumulation_equation"):
            for t, c in self.energy_accumulation_equation.items():
                iscale.constraint_scaling_transform(
                    c, iscale.get_scaling_factor(
                        self.energy_accumulation[t, p],
                        default=1, warning=True))

        # Energy Holdup Integration
        if hasattr(self, "energy_holdup_calculation"):
            for (t, i), c in self.energy_holdup_calculation.items():
                iscale.constraint_scaling_transform(
                    c, iscale.get_scaling_factor(
                        self.energy_holdup[t, i],
                        default=1, warning=True))

        # Energy Balance Equation
        if hasattr(self, "energy_balances"):
            for t, c in self.energy_balances.items():
                iscale.constraint_scaling_transform(
                    c, iscale.get_scaling_factor(
                        self.energy_holdup[t, i],
                        default=1, warning=True))
