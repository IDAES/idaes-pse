##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes".
##############################################################################

"""
Standard IDAES pressure changer model.
"""
from __future__ import division

# Import Python libraries
import logging

# Import Pyomo libraries
from pyomo.environ import SolverFactory, value, Var
from pyomo.opt import TerminationCondition
from pyomo.common.config import ConfigBlock, ConfigValue, In

# Import IDAES cores
from idaes.core import (ControlVolume0D,
                        declare_process_block_class,
                        EnergyBalanceType,
                        MomentumBalanceType,
                        UnitBlockData,
                        useDefault)
from idaes.core.util.config import is_physical_parameter_block, list_of_strings
from idaes.core.util.misc import add_object_ref

__author__ = "Emmanuel Ogbe, Andrew Lee"


@declare_process_block_class("GibbsReactor")
class PressureChangerData(UnitBlockData):
    """
    Standard Compressor/Expander Unit Model Class
    """
    CONFIG = ConfigBlock()
    # Set default values of inherited attributes
    CONFIG.declare("dynamic", ConfigValue(
        domain=In([False]),
        default=False,
        description="Dynamic model flag - must be False",
        doc="""Indicates whether or not dynamics should be considered. 
        **default** is False."""))
    CONFIG.declare("energy_balance_type", ConfigValue(
        default=EnergyBalanceType.enthalpyTotal,
        domain=In(EnergyBalanceType),
        description="Energy balance construction flag",
        doc="""Indicates what type of energy balance should be constructed,
c - EnergyBalanceType.enthalpyTotal.
**Valid values:** {
**EnergyBalanceType.none** - exclude energy balances,
**EnergyBalanceType.enthalpyTotal** - single ethalpy balance for material,
**EnergyBalanceType.enthalpyPhase** - ethalpy balances for each phase,
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
    CONFIG.declare("has_heat_transfer", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="Heat transfer term construction flag",
        doc="""Indicates whether terms for heat transfer should be constructed,
**default** - False.
**Valid values:** {
**True** - include heat transfer terms,
**False** - exclude heat transfer terms.}"""))
    CONFIG.declare("has_phase_equilibrium", ConfigValue(
     default=True,
     domain=In([True, False]),
     description="Phase equilibrium construction flag",
     doc="""Indicates whether terms for phase equilibrium should be
 constructed,
 **default** = False.
 **Valid values:** {
 **True** - include phase equilibrium terms
**False** - exclude phase equilibrium terms.}"""))
    CONFIG.declare("has_work_transfer", ConfigValue(
     default=True,
     domain=In([True, False]),
     description="Work transfer term construction flag",
     doc="""Indicates whether terms for work transfer should be constructed,
 **default** - False.
 **Valid values** {
 **True** - include work transfer terms,
 **False** - exclude work transfer terms.}"""))
 CONFIG.declare("has_pressure_change", ConfigValue(
     default=True,
     domain=In([True, False]),
     description="Pressure change term construction flag",
     doc="""Indicates whether terms for pressure change should be
 constructed,
 **default** - False.
 **Valid values:** {
 **True** - include pressure change terms,
**False** - exclude pressure change terms.}"""))
    #CONFIG.get("has_phase_equilibrium")._default = True	
    #CONFIG.get("has_work_transfer")._default = True
    #CONFIG.get("has_pressure_change")._default = True
    # Add unit model attributes
    CONFIG.declare("dynamic", ConfigValue(
        domain=In([False]),
        default=False,
        description="Dynamic model flag - must be False",
        doc="""Pressure changer is dynamic or not; Default is
False."""))
    CONFIG.declare("inlet_list", ConfigValue(
        domain=list_of_strings,
        description="List of inlet names",
        doc="""A list containing names of inlets (default = None)
                - None - default single inlet
                - list - a list of names for inlets"""))
    CONFIG.declare("num_inlets", ConfigValue(
        domain=int,
        description="Number of inlets to unit",
        doc="""Argument indication number (int) of inlets to construct
            (default = None). Not used if inlet_list arg is provided.
                - None - use inlet_list arg instead
                - int - Inlets will be named with sequential numbers from 1
                        to num_inlets"""))
    CONFIG.declare("outlet_list", ConfigValue(
        domain=list_of_strings,
        description="List of outlet names",
        doc="""A list containing names of outlets (default = None)
                - None - default single outlet
                - list - a list of names for outlets"""))
    CONFIG.declare("num_outlets", ConfigValue(
        domain=int,
        description="Number of outlets to unit",
        doc="""Argument indication number (int) of outlets to construct
            (default = None). Not used if outlet_list arg is provided.
                - None - use outlet_list arg instead
                - int - Outlets will be named with sequential numbers from 1
                        to num_outlets"""))
    CONFIG.declare("compressor", ConfigValue(
        default=True,
        domain=In([True, False]),
        description="Compressor flag",
        doc="""Indicates whether this unit should be considered a
            compressor (True (default), pressure increase) or an expander
            (False, pressure decrease)."""))
    CONFIG.declare("thermodynamic_assumption", ConfigValue(
        default='isothermal',
        domain=In(['isothermal', 'isentropic', 'pump', 'adiabatic']),
        description="Thermodynamic assumption to use",
        doc="""Flag to set the thermodynamic assumption to use for the unit.
                - 'isothermal' (default)
                - 'isentropic'
                - 'pump'
                - 'adiabatic'"""))
    CONFIG.declare("property_package", ConfigValue(
        default=useDefault,
        domain=is_physical_parameter_block,
        description="Property package to use for control volume",
        doc="""Property parameter object used to define property calculations,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PropertyParameterObject** - a PropertyParameterBlock object.}"""))
    CONFIG.declare("property_package_args", ConfigBlock(
        implicit=True,
        description="Arguments to use for constructing property packages",
        doc="""A ConfigBlock with arguments to be passed to a property block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}"""))

    def build(self):
        """

        Args:
            None

        Returns:
            None
        """
        # Call UnitModel.build to setup dynamics
        super(PressureChangerData, self).build()

        # Build Holdup Block
        self.control_volume = controlvolume0D()

        # Set Unit Geometry and holdup Volume
        self.add_set_geometry()

        # Construct performance equations
        self.add_make_performance()

        # Construct equations for thermodynamic assumption
        if self.config.thermodynamic_assumption == "isothermal":
            self.add_make_isothermal()
        elif self.config.thermodynamic_assumption == "isentropic":
            self.add_make_isentropic()
        elif self.config.thermodynamic_assumption == "pump":
            self.add_make_pump()
        elif self.config.thermodynamic_assumption == "adiabatic":
            self.add_make_adiabatic()

    def add_set_geometry(self):
        """
        Define the geometry of the unit as necessary, and link to control volume

        Args:
            None

        Returns:
            None
        """
        # For this case, just create a reference to control volume
        if self.config.has_holdup is True:
            add_object_ref(self, "volume", self.control_volume.volume)

    def add_make_performance(self):
        """
        Define constraints which describe the behaviour of the unit model.

        Args:
            None

        Returns:
            None
        """
        # Set references to balance terms at unit level
        add_object_ref(self, "work_mechanical", self.control_volume.work)
        add_object_ref(self, "deltaP", self.control_volume.deltaP)

        # Performance Variables
        self.ratioP = Var(self.time, initialize=1.0,
                          doc="Pressure Ratio")

        # Pressure Ratio
        @self.Constraint(self.time, doc="Pressure ratio constraint")
        def ratioP_calculation(b, t):
            sf = b.control_volume.scaling_factor_pressure
            return sf*b.ratioP[t]*b.control_volume.properties_in[t].pressure == \
                sf*b.control_volume.properties_out[t].pressure

    def add_make_pump(self):
        """
        Add constraints for the incompresisble fluid assumption

        Args:
            None

        Returns:
            None
        """
        self.work_fluid = Var(
                self.time,
                initialize=1.0,
                doc="Work required to increase the pressure of the liquid")
        self.efficiency_pump = Var(
                self.time,
                initialize=1.0,
                doc="Pump efficiency")

        @self.Constraint(self.time, doc="Pump fluid work constraint")
        def fluid_work_calculation(b, t):
            return b.work_fluid[t] == (
                    (b.control_volume.properties_out[t].pressure -
                     b.control_volume.properties_in[t].pressure) *
                    b.control_volume.properties_out[t].flow_vol)

        # Actual work
        @self.Constraint(self.time, doc="Actual mechanical work calculation")
        def actual_work(b, t):
            sf = b.control_volume.scaling_factor_energy
            if b.config.compressor:
                return sf*b.work_fluid[t] == sf*(
                            b.work_mechanical[t]*b.efficiency_pump[t])
            else:
                return sf*b.work_mechanical[t] == sf*(
                            b.work_fluid[t]*b.efficiency_pump[t])

    def add_make_isothermal(self):
        """
        Add constraints for isothermal assumption.

        Args:
            None

        Returns:
            None
        """
        # Isothermal constraint
        @self.Constraint(self.time, doc="Equate inlet and outlet temperature")
        def isothermal(b, t):
            return b.control_volume.properties_in[t].temperature == \
                       b.control_volume.properties_out[t].temperature

    def add_make_adiabatic(self):
        """
        Add constraints for adiabatic assumption.

        Args:
            None

        Returns:
            None
        """
        # Isothermal constraint
        @self.Constraint(self.time, doc="Equate inlet and outlet enthalpy")
        def adiabatic(b, t):
            return b.control_volume.properties_in[t].enth_mol == \
                       b.control_volume.properties_out[t].enth_mol

    def add_make_isentropic(self):
        """
        Add constraints for isentropic assumption.

        Args:
            None

        Returns:
            None
        """
        # Get indexing sets from holdup block
        add_object_ref(self, "phase_list", self.control_volume.phase_list)
        add_object_ref(self, "component_list", self.control_volume.component_list)

        # Add isentropic variables
        self.efficiency_isentropic = Var(self.time,
                                         initialize=0.8,
                                         doc="Efficiency with respect to an "
                                         "isentropic process [-]")
        self.work_isentropic = Var(self.time,
                                   initialize=0.0,
                                   doc="Work input to unit if isentropic "
                                   "process [-]")

        # Build Isentropic Property block
        self.properties_isentropic = \
            self.control_volume.property_module.PropertyBlock(
                self.time,
                parameters=self.control_volume.config.property_package,
                has_sum_fractions=True,
                calculate_equilibrium_reactions=
                    self.config.has_equilibrium_reactions,
                calculate_phase_equilibrium=self.config.has_phase_equilibrium,
                **self.config.property_package_args)

        # Connect isentropic property states
        @self.Constraint(self.time, doc="Pressure for isentropic calculations")
        def isentropic_pressure(b, t):
            sf = b.control_volume.scaling_factor_pressure
            return sf*b.properties_isentropic[t].pressure == \
                sf*b.ratioP[t]*b.control_volume.properties_in[t].pressure

        # TODO: This assumes isentropic composition is the same as outlet
        @self.Constraint(self.time,
                         self.component_list,
                         doc="Material flows for isentropic properties")
        def isentropic_material(b, t, j):
            return b.properties_isentropic[t].flow_mol_comp[j] == \
                        b.control_volume.properties_out[t].flow_mol_comp[j]

        @self.Constraint(self.time, doc="Isentropic assumption")
        def isentropic(b, t):
            return b.properties_isentropic[t].entr_mol == \
                       b.control_volume.properties_in[t].entr_mol

        # Isentropic work
        @self.Constraint(self.time, doc="Calculate work of isentropic process")
        def isentropic_energy_balance(b, t):
            sf = b.control_volume.scaling_factor_energy
            return sf*b.work_isentropic[t] == sf*(
                    sum(b.properties_isentropic[t].energy_balance_term[p]
                        for p in b.phase_list) -
                    sum(b.control_volume.properties_in[t].energy_balance_term[p]
                        for p in b.phase_list))

        # Actual work
        @self.Constraint(self.time, doc="Actual mechanical work calculation")
        def actual_work(b, t):
            sf = b.control_volume.scaling_factor_energy
            if b.config.compressor:
                return sf*b.work_isentropic[t] == sf*(
                            b.work_mechanical[t]*b.efficiency_isentropic[t])
            else:
                return sf*b.work_mechanical[t] == sf*(
                        b.work_isentropic[t]*b.efficiency_isentropic[t])

    def model_check(blk):
        """
        Check that pressure change matches with compressor argument (i.e. if
        compressor = True, pressure should increase or work should be positive)

        Args:
            None

        Returns:
            None
        """
        if blk.config.compressor:
            # Compressor
            # Check that pressure does not decrease
            if any(blk.deltaP[t].fixed and
                    (value(blk.deltaP[t]) < 0.0) for t in blk.time):
                logger.warning('{} Compressor set with negative deltaP.'
                               .format(blk.name))
            if any(blk.ratioP[t].fixed and
                    (value(blk.ratioP[t]) < 1.0) for t in blk.time):
                logger.warning('{} Compressor set with ratioP less than 1.'
                               .format(blk.name))
            if any(blk.control_volume.properties_out[t].pressure.fixed and
                    (value(blk.control_volume.properties_in[t].pressure) >
                     value(blk.control_volume.properties_out[t].pressure))
                    for t in blk.time):
                logger.warning('{} Compressor set with pressure decrease.'
                               .format(blk.name))
            # Check that work is not negative
            if any(blk.work_mechanical[t].fixed and
                   (value(blk.work_mechanical[t]) < 0.0) for t in blk.time):
                logger.warning('{} Compressor maybe set with negative work.'
                               .format(blk.name))
        else:
            # Expander
            # Check that pressure does not increase
            if any(blk.deltaP[t].fixed and
                    (value(blk.deltaP[t]) > 0.0) for t in blk.time):
                logger.warning('{} Expander/turbine set with positive deltaP.'
                               .format(blk.name))
            if any(blk.ratioP[t].fixed and
                    (value(blk.ratioP[t]) > 1.0) for t in blk.time):
                logger.warning('{} Expander/turbine set with ratioP greater '
                               'than 1.'.format(blk.name))
            if any(blk.control_volume.properties_out[t].pressure.fixed and
                    (value(blk.control_volume.properties_in[t].pressure) <
                     value(blk.control_volume.properties_out[t].pressure))
                    for t in blk.time):
                logger.warning('{} Expander/turbine maybe set with pressure ',
                               'increase.'.format(blk.name))
                # Check that work is not positive
            if any(blk.work_mechanical[t].fixed and
                   (value(blk.work_mechanical[t]) > 0.0) for t in blk.time):
                logger.warning('{} Expander/turbine set with positive work.'
                               .format(blk.name))

        # Run holdup block model checks
        blk.control_volume.model_check()

        # Run model checks on isentropic property block
        try:
            for t in blk.time:
                blk.properties_isentropic[t].model_check()
        except AttributeError:
            pass

    def initialize(blk, state_args={}, routine=None, outlvl=0,
                   solver='ipopt', optarg={'tol': 1e-6}):
        '''
        General wrapper for pressure changer initialisation routines

        Keyword Arguments:
            routine : str stating which initialization routine to execute
                        * None - use routine matching thermodynamic_assumption
                        * 'isentropic' - use isentropic initialization routine
                        * 'isothermal' - use isothermal initialization routine
            state_args : a dict of arguments to be passed to the property
                         package(s) to provide an initial state for
                         initialization (see documentation of the specific
                         property package) (default = {}).
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
        if routine is None:
            # Use routine for specific type of unit
            routine = blk.config.thermodynamic_assumption

        # Call initialisation routine
        if routine is "isentropic":
            blk.init_isentropic(state_args=state_args,
                                outlvl=outlvl,
                                solver=solver,
                                optarg=optarg)
        else:
            # Call the general initialization routine in UnitBlockData
            super(PressureChangerData, blk).initialize(state_args=state_args,
                                                       outlvl=outlvl,
                                                       solver=solver,
                                                       optarg=optarg)

    def init_isentropic(blk, state_args, outlvl, solver, optarg):
        '''
        Initialisation routine for unit (default solver ipopt)

        Keyword Arguments:
            state_args : a dict of arguments to be passed to the property
                         package(s) to provide an initial state for
                         initialization (see documentation of the specific
                         property package) (default = {}).
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
        # Set solver options
        if outlvl > 3:
            stee = True
        else:
            stee = False

        opt = SolverFactory(solver)
        opt.options = optarg

        # ---------------------------------------------------------------------
        # Initialize Isentropic block
        blk.properties_isentropic.initialize(outlvl=outlvl-1,
                                             optarg=optarg,
                                             solver=solver,
                                             **state_args)

        if outlvl > 0:
            logger.info('{} Initialisation Step 1 Complete.'.format(blk.name))

        # ---------------------------------------------------------------------
        # Initialize holdup block
        flags = blk.control_volume.initialize(outlvl=outlvl-1,
                                      optarg=optarg,
                                      solver=solver,
                                      state_args=state_args)

        if outlvl > 0:
            logger.info('{} Initialisation Step 2 Complete.'.format(blk.name))

        # ---------------------------------------------------------------------
        # Solve for isothermal conditions
        if isinstance(blk.properties_isentropic[blk.time[1]].temperature, Var):
            for t in blk.time:
                blk.properties_isentropic[t].temperature.fix()
            blk.isentropic.deactivate()
            results = opt.solve(blk, tee=stee)
            if outlvl > 0:
                if results.solver.termination_condition == \
                        TerminationCondition.optimal:
                    logger.info('{} Initialisation Step 3 Complete.'
                                .format(blk.name))
                else:
                    logger.warning('{} Initialisation Step 3 Failed.'
                                   .format(blk.name))
            for t in blk.time:
                blk.properties_isentropic[t].temperature.unfix()
                blk.isentropic.activate()
        elif outlvl > 0:
            logger.info('{} Initialisation Step 3 Skipped.'.format(blk.name))

        # ---------------------------------------------------------------------
        # Solve unit
        results = opt.solve(blk, tee=stee)

        if outlvl > 0:
            if results.solver.termination_condition == \
                    TerminationCondition.optimal:
                logger.info('{} Initialisation Step 4 Complete.'
                            .format(blk.name))
            else:
                logger.warning('{} Initialisation Step 4 Failed.'
                               .format(blk.name))

        # ---------------------------------------------------------------------
        # Release Inlet state
        blk.control_volume.release_state(flags, outlvl-1)

        if outlvl > 0:
            logger.info('{} Initialisation Complete.'.format(blk.name))
