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

"""
Standard IDAES pressure changer model.
"""

# Import Python libraries
from enum import Enum

# Import Pyomo libraries
from pyomo.environ import SolverFactory, value, Var, Block, Expression,\
    Constraint, Reference
from pyomo.common.config import ConfigBlock, ConfigValue, In

# Import IDAES cores
from idaes.core import (
    ControlVolume0DBlock,
    declare_process_block_class,
    EnergyBalanceType,
    MomentumBalanceType,
    MaterialBalanceType,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.util.config import is_physical_parameter_block
import idaes.logger as idaeslog
import idaes.core.util.unit_costing as costing
from idaes.core.util import scaling as iscale


__author__ = "Emmanuel Ogbe, Andrew Lee"
_log = idaeslog.getLogger(__name__)


class ThermodynamicAssumption(Enum):
    isothermal = 1
    isentropic = 2
    pump = 3
    adiabatic = 4


@declare_process_block_class("PressureChanger")
class PressureChangerData(UnitModelBlockData):
    """
    Standard Compressor/Expander Unit Model Class
    """

    CONFIG = UnitModelBlockData.CONFIG()

    CONFIG.declare(
        "material_balance_type",
        ConfigValue(
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
**MaterialBalanceType.total** - use total material balance.}""",
        ),
    )
    CONFIG.declare(
        "energy_balance_type",
        ConfigValue(
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
**EnergyBalanceType.energyPhase** - energy balances for each phase.}""",
        ),
    )
    CONFIG.declare(
        "momentum_balance_type",
        ConfigValue(
            default=MomentumBalanceType.pressureTotal,
            domain=In(MomentumBalanceType),
            description="Momentum balance construction flag",
            doc="""Indicates what type of momentum balance should be
constructed, **default** - MomentumBalanceType.pressureTotal.
**Valid values:** {
**MomentumBalanceType.none** - exclude momentum balances,
**MomentumBalanceType.pressureTotal** - single pressure balance for material,
**MomentumBalanceType.pressurePhase** - pressure balances for each phase,
**MomentumBalanceType.momentumTotal** - single momentum balance for material,
**MomentumBalanceType.momentumPhase** - momentum balances for each phase.}""",
        ),
    )
    CONFIG.declare(
        "has_phase_equilibrium",
        ConfigValue(
            default=False,
            domain=In([True, False]),
            description="Phase equilibrium construction flag",
            doc="""Indicates whether terms for phase equilibrium should be
constructed, **default** = False.
**Valid values:** {
**True** - include phase equilibrium terms
**False** - exclude phase equilibrium terms.}""",
        ),
    )
    CONFIG.declare(
        "compressor",
        ConfigValue(
            default=True,
            domain=In([True, False]),
            description="Compressor flag",
            doc="""Indicates whether this unit should be considered a
            compressor (True (default), pressure increase) or an expander
            (False, pressure decrease).""",
        ),
    )
    CONFIG.declare(
        "thermodynamic_assumption",
        ConfigValue(
            default=ThermodynamicAssumption.isothermal,
            domain=In(ThermodynamicAssumption),
            description="Thermodynamic assumption to use",
            doc="""Flag to set the thermodynamic assumption to use for the unit.
                - ThermodynamicAssumption.isothermal (default)
                - ThermodynamicAssumption.isentropic
                - ThermodynamicAssumption.pump
                - ThermodynamicAssumption.adiabatic""",
        ),
    )
    CONFIG.declare(
        "property_package",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for control volume",
            doc="""Property parameter object used to define property
calculations, **default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PropertyParameterObject** - a PropertyParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "property_package_args",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing property packages",
            doc="""A ConfigBlock with arguments to be passed to a property
block(s) and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}""",
        ),
    )

    def build(self):
        """

        Args:
            None

        Returns:
            None
        """
        # Call UnitModel.build
        super(PressureChangerData, self).build()

        # Add a control volume to the unit including setting up dynamics.
        self.control_volume = ControlVolume0DBlock(
            default={
                "dynamic": self.config.dynamic,
                "has_holdup": self.config.has_holdup,
                "property_package": self.config.property_package,
                "property_package_args": self.config.property_package_args,
            }
        )

        # Add geomerty variables to control volume
        if self.config.has_holdup:
            self.control_volume.add_geometry()

        # Add inlet and outlet state blocks to control volume
        self.control_volume.add_state_blocks(
            has_phase_equilibrium=self.config.has_phase_equilibrium
        )

        # Add mass balance
        # Set has_equilibrium is False for now
        # TO DO; set has_equilibrium to True
        self.control_volume.add_material_balances(
            balance_type=self.config.material_balance_type,
            has_phase_equilibrium=self.config.has_phase_equilibrium,
        )

        # Add energy balance
        self.control_volume.add_energy_balances(
            balance_type=self.config.energy_balance_type,
            has_work_transfer=True
        )

        # add momentum balance
        self.control_volume.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=True
        )

        # Add Ports
        self.add_inlet_port()
        self.add_outlet_port()

        # Set Unit Geometry and holdup Volume
        if self.config.has_holdup is True:
            self.volume = Reference(self.control_volume.volume[:])

        # Construct performance equations
        # Set references to balance terms at unit level
        # Add Work transfer variable 'work'
        self.work_mechanical = Reference(self.control_volume.work[:])

        # Add Momentum balance variable 'deltaP'
        self.deltaP = Reference(self.control_volume.deltaP[:])

        # Performance Variables
        self.ratioP = Var(
            self.flowsheet().config.time, initialize=1.0, doc="Pressure Ratio"
        )

        # Pressure Ratio
        @self.Constraint(self.flowsheet().config.time,
                         doc="Pressure ratio constraint")
        def ratioP_calculation(b, t):
            return (
                b.ratioP[t] * b.control_volume.properties_in[t].pressure
                == b.control_volume.properties_out[t].pressure
            )

        # Construct equations for thermodynamic assumption
        if (self.config.thermodynamic_assumption ==
                ThermodynamicAssumption.isothermal):
            self.add_isothermal()
        elif (self.config.thermodynamic_assumption ==
              ThermodynamicAssumption.isentropic):
            self.add_isentropic()
        elif (self.config.thermodynamic_assumption ==
              ThermodynamicAssumption.pump):
            self.add_pump()
        elif (self.config.thermodynamic_assumption ==
              ThermodynamicAssumption.adiabatic):
            self.add_adiabatic()

    def add_pump(self):
        """
        Add constraints for the incompressible fluid assumption

        Args:
            None

        Returns:
            None
        """
        units_meta = self.config.property_package.get_metadata()

        self.work_fluid = Var(
            self.flowsheet().config.time,
            initialize=1.0,
            doc="Work required to increase the pressure of the liquid",
            units=units_meta.get_derived_units("power")
        )
        self.efficiency_pump = Var(
            self.flowsheet().config.time, initialize=1.0, doc="Pump efficiency"
        )

        @self.Constraint(self.flowsheet().config.time,
                         doc="Pump fluid work constraint")
        def fluid_work_calculation(b, t):
            return b.work_fluid[t] == (
                (
                    b.control_volume.properties_out[t].pressure
                    - b.control_volume.properties_in[t].pressure
                )
                * b.control_volume.properties_out[t].flow_vol
            )

        # Actual work
        @self.Constraint(
            self.flowsheet().config.time,
            doc="Actual mechanical work calculation"
        )
        def actual_work(b, t):
            if b.config.compressor:
                return b.work_fluid[t] == (
                    b.work_mechanical[t] * b.efficiency_pump[t]
                )
            else:
                return b.work_mechanical[t] == (
                    b.work_fluid[t] * b.efficiency_pump[t]
                )

    def add_isothermal(self):
        """
        Add constraints for isothermal assumption.

        Args:
            None

        Returns:
            None
        """
        # Isothermal constraint
        @self.Constraint(
            self.flowsheet().config.time,
            doc="For isothermal condition: Equate inlet and "
            "outlet temperature",
        )
        def isothermal(b, t):
            return (
                b.control_volume.properties_in[t].temperature
                == b.control_volume.properties_out[t].temperature
            )

    def add_adiabatic(self):
        """
        Add constraints for adiabatic assumption.

        Args:
            None

        Returns:
            None
        """
        # Isothermal constraint
        @self.Constraint(
            self.flowsheet().config.time,
            doc="For isothermal condition: Equate inlet and outlet enthalpy",
        )
        def adiabatic(b, t):
            return (
                b.control_volume.properties_in[t].enth_mol
                == b.control_volume.properties_out[t].enth_mol
            )

    def add_isentropic(self):
        """
        Add constraints for isentropic assumption.

        Args:
            None

        Returns:
            None
        """
        units_meta = self.config.property_package.get_metadata()

        # Get indexing sets from control volume
        # Add isentropic variables
        self.efficiency_isentropic = Var(
            self.flowsheet().config.time,
            initialize=0.8,
            doc="Efficiency with respect to an isentropic process [-]",
        )
        self.work_isentropic = Var(
            self.flowsheet().config.time,
            initialize=0.0,
            doc="Work input to unit if isentropic process",
            units=units_meta.get_derived_units("power")
        )

        # Build isentropic state block
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = self.config.has_phase_equilibrium
        tmp_dict["defined_state"] = False

        self.properties_isentropic = (
            self.config.property_package.build_state_block(
                self.flowsheet().config.time,
                doc="isentropic properties at outlet",
                default=tmp_dict)
        )

        # Connect isentropic state block properties
        @self.Constraint(
            self.flowsheet().config.time,
            doc="Pressure for isentropic calculations"
        )
        def isentropic_pressure(b, t):
            return (
                b.properties_isentropic[t].pressure
                == b.control_volume.properties_out[t].pressure
            )

        # This assumes isentropic composition is the same as outlet
        self.add_state_material_balances(self.config.material_balance_type,
                                         self.properties_isentropic,
                                         self.control_volume.properties_out)

        # This assumes isentropic entropy is the same as inlet
        @self.Constraint(self.flowsheet().config.time,
                         doc="Isentropic assumption")
        def isentropic(b, t):
            return (
                b.properties_isentropic[t].entr_mol
                == b.control_volume.properties_in[t].entr_mol
            )

        # Isentropic work
        @self.Constraint(
            self.flowsheet().config.time,
            doc="Calculate work of isentropic process"
        )
        def isentropic_energy_balance(b, t):
            return b.work_isentropic[t] == (
                sum(
                    b.properties_isentropic[t].get_enthalpy_flow_terms(p)
                    for p in b.properties_isentropic.phase_list
                )
                - sum(
                    b.control_volume.properties_in[
                        t].get_enthalpy_flow_terms(p)
                    for p in b.control_volume.properties_in.phase_list
                )
            )

        # Actual work
        @self.Constraint(
            self.flowsheet().config.time,
            doc="Actual mechanical work calculation"
        )
        def actual_work(b, t):
            if b.config.compressor:
                return b.work_isentropic[t] == (
                    b.work_mechanical[t] * b.efficiency_isentropic[t]
                )
            else:
                return b.work_mechanical[t] == (
                    b.work_isentropic[t] * b.efficiency_isentropic[t]
                )

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
            if any(
                blk.deltaP[t].fixed and (value(blk.deltaP[t]) < 0.0)
                for t in blk.flowsheet().config.time
            ):
                _log.warning("{} Compressor set with negative deltaP."
                             .format(blk.name))
            if any(
                blk.ratioP[t].fixed and (value(blk.ratioP[t]) < 1.0)
                for t in blk.flowsheet().config.time
            ):
                _log.warning(
                    "{} Compressor set with ratioP less than 1."
                    .format(blk.name)
                )
            if any(
                blk.control_volume.properties_out[t].pressure.fixed
                and (
                    value(blk.control_volume.properties_in[t].pressure)
                    > value(blk.control_volume.properties_out[t].pressure)
                )
                for t in blk.flowsheet().config.time
            ):
                _log.warning(
                    "{} Compressor set with pressure decrease."
                    .format(blk.name)
                )
            # Check that work is not negative
            if any(
                blk.work_mechanical[t].fixed and (
                    value(blk.work_mechanical[t]) < 0.0)
                for t in blk.flowsheet().config.time
            ):
                _log.warning(
                    "{} Compressor maybe set with negative work."
                    .format(blk.name)
                )
        else:
            # Expander
            # Check that pressure does not increase
            if any(
                blk.deltaP[t].fixed and (value(blk.deltaP[t]) > 0.0)
                for t in blk.flowsheet().config.time
            ):
                _log.warning(
                    "{} Expander/turbine set with positive deltaP."
                    .format(blk.name)
                )
            if any(
                blk.ratioP[t].fixed and (value(blk.ratioP[t]) > 1.0)
                for t in blk.flowsheet().config.time
            ):
                _log.warning(
                    "{} Expander/turbine set with ratioP greater "
                    "than 1.".format(blk.name)
                )
            if any(
                blk.control_volume.properties_out[t].pressure.fixed
                and (
                    value(blk.control_volume.properties_in[t].pressure)
                    < value(blk.control_volume.properties_out[t].pressure)
                )
                for t in blk.flowsheet().config.time
            ):
                _log.warning(
                    "{} Expander/turbine maybe set with pressure ",
                    "increase.".format(blk.name),
                )
            # Check that work is not positive
            if any(
                blk.work_mechanical[t].fixed and (
                    value(blk.work_mechanical[t]) > 0.0)
                for t in blk.flowsheet().config.time
            ):
                _log.warning(
                    "{} Expander/turbine set with positive work."
                    .format(blk.name)
                )

        # Run holdup block model checks
        blk.control_volume.model_check()

        # Run model checks on isentropic property block
        try:
            for t in blk.flowsheet().config.time:
                blk.properties_in[t].model_check()
        except AttributeError:
            pass

    def initialize(
        blk,
        state_args=None,
        routine=None,
        outlvl=idaeslog.NOTSET,
        solver="ipopt",
        optarg={"tol": 1e-6},
    ):
        """
        General wrapper for pressure changer initialization routines

        Keyword Arguments:
            routine : str stating which initialization routine to execute
                        * None - use routine matching thermodynamic_assumption
                        * 'isentropic' - use isentropic initialization routine
                        * 'isothermal' - use isothermal initialization routine
            state_args : a dict of arguments to be passed to the property
                         package(s) to provide an initial state for
                         initialization (see documentation of the specific
                         property package) (default = {}).
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default={'tol': 1e-6})
            solver : str indicating whcih solver to use during
                     initialization (default = 'ipopt')

        Returns:
            None
        """
        # if costing block exists, deactivate
        try:
            blk.costing.deactivate()
        except AttributeError:
            pass

        if routine is None:
            # Use routine for specific type of unit
            routine = blk.config.thermodynamic_assumption

        # Call initialization routine
        if routine is ThermodynamicAssumption.isentropic:
            blk.init_isentropic(
                state_args=state_args,
                outlvl=outlvl,
                solver=solver,
                optarg=optarg
            )
        else:
            # Call the general initialization routine in UnitModelBlockData
            super(PressureChangerData, blk).initialize(
                state_args=state_args,
                outlvl=outlvl,
                solver=solver,
                optarg=optarg
            )
        # if costing block exists, activate
        try:
            blk.costing.activate()
            costing.initialize(blk.costing)
        except AttributeError:
            pass

    def init_isentropic(blk, state_args, outlvl, solver, optarg):
        """
        Initialization routine for unit (default solver ipopt)

        Keyword Arguments:
            state_args : a dict of arguments to be passed to the property
                         package(s) to provide an initial state for
                         initialization (see documentation of the specific
                         property package) (default = {}).
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default={'tol': 1e-6})
            solver : str indicating whcih solver to use during
                     initialization (default = 'ipopt')

        Returns:
            None
        """
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")
        # Set solver options
        opt = SolverFactory(solver)
        opt.options = optarg

        cv = blk.control_volume
        t0 = blk.flowsheet().config.time.first()
        state_args_out = {}
        if state_args is None:
            state_args = {}
            state_dict = (
                cv.properties_in[t0].define_port_members())

            for k in state_dict.keys():
                if state_dict[k].is_indexed():
                    state_args[k] = {}
                    for m in state_dict[k].keys():
                        state_args[k][m] = state_dict[k][m].value
                else:
                    state_args[k] = state_dict[k].value

            # Get initialisation guesses for outlet and isentropic states
            for k in state_args:
                if k == "pressure":
                    # Work out how to estimate outlet pressure
                    if cv.properties_out[t0].pressure.fixed:
                        # Fixed outlet pressure, use this value
                        state_args_out[k] = value(
                            cv.properties_out[t0].pressure)
                    elif blk.deltaP[t0].fixed:
                        state_args_out[k] = value(
                            state_args[k] + blk.deltaP[t0])
                    elif blk.ratioP[t0].fixed:
                        state_args_out[k] = value(
                            state_args[k] * blk.ratioP[t0])
                    else:
                        # Not obvious what to do, use inlet state
                        state_args_out[k] = state_args[k]
                else:
                    state_args_out[k] = state_args[k]

        # Initialize state blocks
        flags = cv.properties_in.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            hold_state=True,
            state_args=state_args,
        )
        cv.properties_out.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            hold_state=False,
            state_args=state_args_out,
        )

        init_log.info_high("Initialization Step 1 Complete.")
        # ---------------------------------------------------------------------
        # Initialize Isentropic block

        blk.properties_isentropic.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_out,
        )

        init_log.info_high("Initialization Step 2 Complete.")

        # ---------------------------------------------------------------------
        # Solve for isothermal conditions
        if isinstance(
            blk.properties_isentropic[
                blk.flowsheet().config.time.first()].temperature,
            Var,
        ):
            blk.properties_isentropic[:].temperature.fix()
        elif isinstance(
            blk.properties_isentropic[
                blk.flowsheet().config.time.first()].enth_mol,
            Var,
        ):
            blk.properties_isentropic[:].enth_mol.fix()
        elif isinstance(
            blk.properties_isentropic[
                blk.flowsheet().config.time.first()].temperature,
            Expression,
        ):
            def tmp_rule(b, t):
                return blk.properties_isentropic[t].temperature == \
                    blk.control_volume.properties_in[t].temperature
            blk.tmp_init_constraint = Constraint(
                blk.flowsheet().config.time, rule=tmp_rule)

        blk.isentropic.deactivate()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high("Initialization Step 3 {}."
                           .format(idaeslog.condition(res)))

        if isinstance(
            blk.properties_isentropic[
                blk.flowsheet().config.time.first()].temperature,
            Var,
        ):
            blk.properties_isentropic[:].temperature.unfix()
        elif isinstance(
            blk.properties_isentropic[
                blk.flowsheet().config.time.first()].enth_mol,
            Var,
        ):
            blk.properties_isentropic[:].enth_mol.unfix()
        elif isinstance(
            blk.properties_isentropic[
                blk.flowsheet().config.time.first()].temperature,
            Expression,
        ):
            blk.del_component(blk.tmp_init_constraint)

        blk.isentropic.activate()

        # ---------------------------------------------------------------------
        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high("Initialization Step 4 {}."
                           .format(idaeslog.condition(res)))

        # ---------------------------------------------------------------------
        # Release Inlet state
        blk.control_volume.release_state(flags, outlvl)
        init_log.info(
            "Initialization Complete: {}"
            .format(idaeslog.condition(res))
        )

    def _get_performance_contents(self, time_point=0):
        var_dict = {}
        if hasattr(self, "deltaP"):
            var_dict["Mechanical Work"] = self.work_mechanical[time_point]
        if hasattr(self, "deltaP"):
            var_dict["Pressure Change"] = self.deltaP[time_point]
        if hasattr(self, "ratioP"):
            var_dict["Pressure Ratio"] = self.ratioP[time_point]
        if hasattr(self, "efficiency_pump"):
            var_dict["Efficiency"] = self.efficiency_pump[time_point]
        if hasattr(self, "efficiency_isentropic"):
            var_dict["Isentropic Efficiency"] = \
                self.efficiency_isentropic[time_point]

        return {"vars": var_dict}

    def get_costing(self, module=costing, year=None, **kwargs):
        if not hasattr(self.flowsheet(), "costing"):
            self.flowsheet().get_costing(year=year)

        self.costing = Block()
        module.pressure_changer_costing(
            self.costing,
            **kwargs)

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        if hasattr(self, "work_fluid"):
            for t, v in self.work_fluid.items():
                iscale.set_scaling_factor(
                    v,
                    iscale.get_scaling_factor(
                        self.control_volume.work[t],
                        default=1,
                        warning=True))

        if hasattr(self, "work_mechanical"):
            for t, v in self.work_mechanical.items():
                iscale.set_scaling_factor(
                    v,
                    iscale.get_scaling_factor(
                        self.control_volume.work[t],
                        default=1,
                        warning=True))

        if hasattr(self, "work_isentropic"):
            for t, v in self.work_isentropic.items():
                iscale.set_scaling_factor(
                    v,
                    iscale.get_scaling_factor(
                        self.control_volume.work[t],
                        default=1,
                        warning=True))

        if hasattr(self, "ratioP_calculation"):
            for t, c in self.ratioP_calculation.items():
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(
                        self.control_volume.properties_in[t].pressure,
                        default=1,
                        warning=True))

        if hasattr(self, "fluid_work_calculation"):
            for t, c in self.fluid_work_calculation.items():
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(
                        self.control_volume.deltaP[t],
                        default=1,
                        warning=True))

        if hasattr(self, "actual_work"):
            for t, c in self.actual_work.items():
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(
                        self.control_volume.work[t],
                        default=1,
                        warning=True))

        if hasattr(self, "adiabatic"):
            for t, c in self.adiabatic.items():
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(
                        self.control_volume.properties_in[t].enth_mol,
                        default=1,
                        warning=True))

        if hasattr(self, "isentropic_pressure"):
            for t, c in self.isentropic_pressure.items():
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(
                        self.control_volume.properties_in[t].pressure,
                        default=1,
                        warning=True))

        if hasattr(self, "isentropic"):
            for t, c in self.isentropic.items():
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(
                        self.control_volume.properties_in[t].entr_mol,
                        default=1,
                        warning=True))

        if hasattr(self, "isentropic_energy_balance"):
            for t, c in self.isentropic_energy_balance.items():
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(
                        self.control_volume.work[t],
                        default=1,
                        warning=True))

        if hasattr(self, "costing"):
            # import costing scaling factors
            costing.calculate_scaling_factors(self.costing)


@declare_process_block_class("Turbine", doc="Isentropic turbine model")
class TurbineData(PressureChangerData):
    # Pressure changer with isentropic turbine options
    CONFIG = PressureChangerData.CONFIG()
    CONFIG.compressor = False
    CONFIG.get("compressor")._default = False
    CONFIG.get("compressor")._domain = In([False])
    CONFIG.thermodynamic_assumption = ThermodynamicAssumption.isentropic
    CONFIG.get("thermodynamic_assumption")._default = \
        ThermodynamicAssumption.isentropic


@declare_process_block_class("Compressor", doc="Isentropic compressor model")
class CompressorData(PressureChangerData):
    # Pressure changer with isentropic turbine options
    CONFIG = PressureChangerData.CONFIG()
    CONFIG.compressor = True
    CONFIG.get("compressor")._default = True
    CONFIG.get("compressor")._domain = In([True])
    CONFIG.thermodynamic_assumption = ThermodynamicAssumption.isentropic
    CONFIG.get("thermodynamic_assumption")._default = \
        ThermodynamicAssumption.isentropic


@declare_process_block_class("Pump", doc="Pump model")
class PumpData(PressureChangerData):
    # Pressure changer with isentropic turbine options
    CONFIG = PressureChangerData.CONFIG()
    CONFIG.compressor = True
    CONFIG.get("compressor")._default = True
    CONFIG.get("compressor")._domain = In([True])
    CONFIG.thermodynamic_assumption = ThermodynamicAssumption.pump
    CONFIG.get("thermodynamic_assumption")._default = \
        ThermodynamicAssumption.pump
