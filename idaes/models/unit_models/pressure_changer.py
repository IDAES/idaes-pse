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
Standard IDAES pressure changer model.
"""

# Import Python libraries
from enum import Enum

# Import Pyomo libraries
from pyomo.environ import (
    value,
    Var,
    Expression,
    Constraint,
    Reference,
    check_optimal_termination,
    Reals,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In, Bool

# Import IDAES cores
from idaes.core import (
    ControlVolume0DBlock,
    declare_process_block_class,
    EnergyBalanceType,
    MomentumBalanceType,
    MaterialBalanceType,
    ProcessBlockData,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.util.exceptions import PropertyNotSupportedError, InitializationError
from idaes.core.util.config import is_physical_parameter_block
import idaes.logger as idaeslog
from idaes.core.util import scaling as iscale
from idaes.core.solvers import get_solver


__author__ = "Emmanuel Ogbe, Andrew Lee"
_log = idaeslog.getLogger(__name__)


class ThermodynamicAssumption(Enum):
    isothermal = 1
    isentropic = 2
    pump = 3
    adiabatic = 4


@declare_process_block_class("IsentropicPerformanceCurve")
class IsentropicPerformanceCurveData(ProcessBlockData):
    """Block that holds performance curves. Typically, these are in the form of
    constraints that relate head, efficiency, or pressure ratio to volumetric
    or mass flow.  Additional variables can be included if needed, such as
    speed. For convenience an option is provided to add head expressions to the
    block. performance curves, and any additional variables, constraints, or
    expressions can be added to this block either via callback provided to the
    configuration, or after the model is constructed."""

    CONFIG = ProcessBlockData.CONFIG(
        doc="Configuration dictionary for the performance curve block."
    )
    CONFIG.declare(
        "build_callback",
        ConfigValue(
            default=None, doc="Optional callback to add performance curve constraints"
        ),
    )
    CONFIG.declare(
        "build_head_expressions",
        ConfigValue(
            default=True,
            domain=bool,
            doc="If true add expressions for 'head' and 'head_isentropic'. "
            "These expressions can be used in performance curve constraints.",
        ),
    )

    def has_constraints(self):
        for o in self.component_data_objects(Constraint):
            return True
        return False

    def build(self):
        super().build()
        if self.config.build_head_expressions:
            try:

                @self.Expression(self.flowsheet().time)
                def head_isentropic(self, t):  # units are energy/mass
                    b = self.parent_block()
                    if hasattr(b.control_volume.properties_in[t], "flow_mass"):
                        return (
                            b.work_isentropic[t]
                            / b.control_volume.properties_in[t].flow_mass
                        )
                    else:
                        return (
                            b.work_isentropic[t]
                            / b.control_volume.properties_in[t].flow_mol
                            / b.control_volume.properties_in[t].mw
                        )

                @self.Expression(self.flowsheet().time)
                def head(self, t):  # units are energy/mass
                    b = self.parent_block()
                    if hasattr(b.control_volume.properties_in[t], "flow_mass"):
                        return (
                            b.work_mechanical[t]
                            / b.control_volume.properties_in[t].flow_mass
                        )
                    else:
                        return (
                            b.work_mechanical[t]
                            / b.control_volume.properties_in[t].flow_mol
                            / b.control_volume.properties_in[t].mw
                        )

            except PropertyNotSupportedError:
                _log.exception(
                    "flow_mass or flow_mol and mw are not supported by the "
                    "property package but are required for isentropic pressure"
                    " changer head calculation"
                )
                raise

        if self.config.build_callback is not None:
            self.config.build_callback(self)


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
            domain=Bool,
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
            domain=Bool,
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
    CONFIG.declare(
        "support_isentropic_performance_curves",
        ConfigValue(
            default=False,
            domain=Bool,
            doc="Include a block for performance curves, configure via"
            " isentropic_performance_curves.",
        ),
    )
    CONFIG.declare(
        "isentropic_performance_curves",
        IsentropicPerformanceCurveData.CONFIG(),
        # doc included in IsentropicPerformanceCurveData
    )

    def build(self):
        """

        Args:
            None

        Returns:
            None
        """
        # Call UnitModel.build
        super().build()

        # Add a control volume to the unit including setting up dynamics.
        self.control_volume = ControlVolume0DBlock(
            dynamic=self.config.dynamic,
            has_holdup=self.config.has_holdup,
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
        )

        # Add geometry variables to control volume
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
        eb = self.control_volume.add_energy_balances(
            balance_type=self.config.energy_balance_type, has_work_transfer=True
        )

        # add momentum balance
        self.control_volume.add_momentum_balances(
            balance_type=self.config.momentum_balance_type, has_pressure_change=True
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
        # If the 'work' variable wasn't already built on the control volume but is needed, create it now.
        if (
            not hasattr(self.control_volume, "work")
            and self.config.thermodynamic_assumption == ThermodynamicAssumption.pump
            and eb is None
        ):
            units = self.config.property_package.get_metadata().get_derived_units
            self.control_volume.work = Var(
                self.flowsheet().time,
                domain=Reals,
                initialize=0.0,
                doc="Work transferred into control volume",
                units=units("power"),
            )
        self.work_mechanical = Reference(self.control_volume.work[:])

        # Add Momentum balance variable 'deltaP'
        self.deltaP = Reference(self.control_volume.deltaP[:])

        # Performance Variables
        self.ratioP = Var(self.flowsheet().time, initialize=1.0, doc="Pressure Ratio")

        # Pressure Ratio
        @self.Constraint(self.flowsheet().time, doc="Pressure ratio constraint")
        def ratioP_calculation(self, t):
            return (
                self.ratioP[t] * self.control_volume.properties_in[t].pressure
                == self.control_volume.properties_out[t].pressure
            )

        # Construct equations for thermodynamic assumption
        if self.config.thermodynamic_assumption == ThermodynamicAssumption.isothermal:
            self.add_isothermal()
        elif self.config.thermodynamic_assumption == ThermodynamicAssumption.isentropic:
            self.add_isentropic()
        elif self.config.thermodynamic_assumption == ThermodynamicAssumption.pump:
            self.add_pump()
        elif self.config.thermodynamic_assumption == ThermodynamicAssumption.adiabatic:
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
            self.flowsheet().time,
            initialize=1.0,
            doc="Work required to increase the pressure of the liquid",
            units=units_meta.get_derived_units("power"),
        )
        self.efficiency_pump = Var(
            self.flowsheet().time, initialize=1.0, doc="Pump efficiency"
        )

        @self.Constraint(self.flowsheet().time, doc="Pump fluid work constraint")
        def fluid_work_calculation(self, t):
            return self.work_fluid[t] == (
                (
                    self.control_volume.properties_out[t].pressure
                    - self.control_volume.properties_in[t].pressure
                )
                * self.control_volume.properties_out[t].flow_vol
            )

        # Actual work
        @self.Constraint(
            self.flowsheet().time, doc="Actual mechanical work calculation"
        )
        def actual_work(self, t):
            if self.config.compressor:
                return self.work_fluid[t] == (
                    self.work_mechanical[t] * self.efficiency_pump[t]
                )
            else:
                return self.work_mechanical[t] == (
                    self.work_fluid[t] * self.efficiency_pump[t]
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
            self.flowsheet().time,
            doc="For isothermal condition: Equate inlet and " "outlet temperature",
        )
        def isothermal(self, t):
            return (
                self.control_volume.properties_in[t].temperature
                == self.control_volume.properties_out[t].temperature
            )

    def add_adiabatic(self):
        """
        Add constraints for adiabatic assumption.

        Args:
            None

        Returns:
            None
        """

        @self.Constraint(self.flowsheet().time)
        def zero_work_equation(self, t):
            return self.control_volume.work[t] == 0

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
            self.flowsheet().time,
            initialize=0.8,
            doc="Efficiency with respect to an isentropic process [-]",
        )
        self.work_isentropic = Var(
            self.flowsheet().time,
            initialize=0.0,
            doc="Work input to unit if isentropic process",
            units=units_meta.get_derived_units("power"),
        )

        # Build isentropic state block
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = self.config.has_phase_equilibrium
        tmp_dict["defined_state"] = False

        self.properties_isentropic = self.config.property_package.build_state_block(
            self.flowsheet().time, doc="isentropic properties at outlet", **tmp_dict
        )

        # Connect isentropic state block properties
        @self.Constraint(
            self.flowsheet().time, doc="Pressure for isentropic calculations"
        )
        def isentropic_pressure(self, t):
            return (
                self.properties_isentropic[t].pressure
                == self.control_volume.properties_out[t].pressure
            )

        # This assumes isentropic composition is the same as outlet
        self.add_state_material_balances(
            self.config.material_balance_type,
            self.properties_isentropic,
            self.control_volume.properties_out,
        )

        # This assumes isentropic entropy is the same as inlet
        @self.Constraint(self.flowsheet().time, doc="Isentropic assumption")
        def isentropic(self, t):
            return (
                self.properties_isentropic[t].entr_mol
                == self.control_volume.properties_in[t].entr_mol
            )

        # Isentropic work
        @self.Constraint(
            self.flowsheet().time, doc="Calculate work of isentropic process"
        )
        def isentropic_energy_balance(self, t):
            return self.work_isentropic[t] == (
                sum(
                    self.properties_isentropic[t].get_enthalpy_flow_terms(p)
                    for p in self.properties_isentropic.phase_list
                )
                - sum(
                    self.control_volume.properties_in[t].get_enthalpy_flow_terms(p)
                    for p in self.control_volume.properties_in.phase_list
                )
            )

        # Actual work
        @self.Constraint(
            self.flowsheet().time, doc="Actual mechanical work calculation"
        )
        def actual_work(self, t):
            if self.config.compressor:
                return self.work_isentropic[t] == (
                    self.work_mechanical[t] * self.efficiency_isentropic[t]
                )
            else:
                return self.work_mechanical[t] == (
                    self.work_isentropic[t] * self.efficiency_isentropic[t]
                )

        if self.config.support_isentropic_performance_curves:
            self.performance_curve = IsentropicPerformanceCurve(
                **self.config.isentropic_performance_curves
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
                for t in blk.flowsheet().time
            ):
                _log.warning("{} Compressor set with negative deltaP.".format(blk.name))
            if any(
                blk.ratioP[t].fixed and (value(blk.ratioP[t]) < 1.0)
                for t in blk.flowsheet().time
            ):
                _log.warning(
                    "{} Compressor set with ratioP less than 1.".format(blk.name)
                )
            if any(
                blk.control_volume.properties_out[t].pressure.fixed
                and (
                    value(blk.control_volume.properties_in[t].pressure)
                    > value(blk.control_volume.properties_out[t].pressure)
                )
                for t in blk.flowsheet().time
            ):
                _log.warning(
                    "{} Compressor set with pressure decrease.".format(blk.name)
                )
            # Check that work is not negative
            if any(
                blk.work_mechanical[t].fixed and (value(blk.work_mechanical[t]) < 0.0)
                for t in blk.flowsheet().time
            ):
                _log.warning(
                    "{} Compressor maybe set with negative work.".format(blk.name)
                )
        else:
            # Expander
            # Check that pressure does not increase
            if any(
                blk.deltaP[t].fixed and (value(blk.deltaP[t]) > 0.0)
                for t in blk.flowsheet().time
            ):
                _log.warning(
                    "{} Expander/turbine set with positive deltaP.".format(blk.name)
                )
            if any(
                blk.ratioP[t].fixed and (value(blk.ratioP[t]) > 1.0)
                for t in blk.flowsheet().time
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
                for t in blk.flowsheet().time
            ):
                _log.warning(
                    "{} Expander/turbine maybe set with pressure "
                    "increase.".format(blk.name),
                )
            # Check that work is not positive
            if any(
                blk.work_mechanical[t].fixed and (value(blk.work_mechanical[t]) > 0.0)
                for t in blk.flowsheet().time
            ):
                _log.warning(
                    "{} Expander/turbine set with positive work.".format(blk.name)
                )

        # Run holdup block model checks
        blk.control_volume.model_check()

        # Run model checks on isentropic property block
        try:
            for t in blk.flowsheet().time:
                blk.properties_in[t].model_check()
        except AttributeError:
            pass

    def initialize_build(
        blk,
        state_args=None,
        routine=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
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
            optarg : solver options dictionary object (default=None, use
                     default solver options)
            solver : str indicating which solver to use during
                     initialization (default = None, use default solver)

        Returns:
            None
        """
        if routine is None:
            # Use routine for specific type of unit
            routine = blk.config.thermodynamic_assumption

        # Call initialization routine
        if routine is ThermodynamicAssumption.isentropic:
            blk.init_isentropic(
                state_args=state_args, outlvl=outlvl, solver=solver, optarg=optarg
            )
        elif routine is ThermodynamicAssumption.adiabatic:
            blk.init_adiabatic(
                state_args=state_args, outlvl=outlvl, solver=solver, optarg=optarg
            )
        else:
            # Call the general initialization routine in UnitModelBlockData
            super().initialize_build(
                state_args=state_args, outlvl=outlvl, solver=solver, optarg=optarg
            )

    def init_adiabatic(blk, state_args, outlvl, solver, optarg):
        """
        Initialization routine for adiabatic pressure changers.

        Keyword Arguments:
            state_args : a dict of arguments to be passed to the property
                         package(s) to provide an initial state for
                         initialization (see documentation of the specific
                         property package) (default = {}).
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default={})
            solver : str indicating which solver to use during
                     initialization (default = None)

        Returns:
            None
        """
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")

        # Create solver
        opt = get_solver(solver, optarg)

        cv = blk.control_volume
        t0 = blk.flowsheet().time.first()
        state_args_out = {}

        if state_args is None:
            state_args = {}
            state_dict = cv.properties_in[t0].define_port_members()

            for k in state_dict.keys():
                if state_dict[k].is_indexed():
                    state_args[k] = {}
                    for m in state_dict[k].keys():
                        state_args[k][m] = state_dict[k][m].value
                else:
                    state_args[k] = state_dict[k].value

        # Get initialisation guesses for outlet and isentropic states
        for k in state_args:
            if k == "pressure" and k not in state_args_out:
                # Work out how to estimate outlet pressure
                if cv.properties_out[t0].pressure.fixed:
                    # Fixed outlet pressure, use this value
                    state_args_out[k] = value(cv.properties_out[t0].pressure)
                elif blk.deltaP[t0].fixed:
                    state_args_out[k] = value(state_args[k] + blk.deltaP[t0])
                elif blk.ratioP[t0].fixed:
                    state_args_out[k] = value(state_args[k] * blk.ratioP[t0])
                else:
                    # Not obvious what to do, use inlet state
                    state_args_out[k] = state_args[k]
            elif k not in state_args_out:
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
        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high("Initialization Step 2 {}.".format(idaeslog.condition(res)))

        # ---------------------------------------------------------------------
        # Release Inlet state
        blk.control_volume.release_state(flags, outlvl)

        if not check_optimal_termination(res):
            raise InitializationError(
                f"{blk.name} failed to initialize successfully. Please check "
                f"the output logs for more information."
            )

        init_log.info(f"Initialization Complete: {idaeslog.condition(res)}")

    def init_isentropic(blk, state_args, outlvl, solver, optarg):
        """
        Initialization routine for isentropic pressure changers.

        Keyword Arguments:
            state_args : a dict of arguments to be passed to the property
                         package(s) to provide an initial state for
                         initialization (see documentation of the specific
                         property package) (default = {}).
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default={})
            solver : str indicating which solver to use during
                     initialization (default = None)

        Returns:
            None
        """
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")

        # Create solver
        opt = get_solver(solver, optarg)

        cv = blk.control_volume
        t0 = blk.flowsheet().time.first()
        state_args_out = {}

        # performance curves exist and are active so initialize with them
        activate_performance_curves = (
            hasattr(blk, "performance_curve")
            and blk.performance_curve.has_constraints()
            and blk.performance_curve.active
        )
        if activate_performance_curves:
            blk.performance_curve.deactivate()
            # The performance curves will provide (maybe indirectly) efficiency
            # and/or pressure ratio. To get through the standard isentropic
            # pressure changer init, we'll see if the user provided a guess for
            # pressure ratio or isentropic efficiency and fix them if needed. If
            # not fixed and no guess provided, fill in something reasonable
            # until the performance curves are turned on.
            unfix_eff = {}
            unfix_ratioP = {}
            for t in blk.flowsheet().time:
                if not (
                    blk.ratioP[t].fixed
                    or blk.deltaP[t].fixed
                    or cv.properties_out[t].pressure.fixed
                ):
                    if blk.config.compressor:
                        if not (
                            value(blk.ratioP[t]) >= 1.01 and value(blk.ratioP[t]) <= 50
                        ):
                            blk.ratioP[t] = 1.8
                    else:
                        if not (
                            value(blk.ratioP[t]) >= 0.01
                            and value(blk.ratioP[t]) <= 0.999
                        ):
                            blk.ratioP[t] = 0.7
                    blk.ratioP[t].fix()
                    unfix_ratioP[t] = True
                if not blk.efficiency_isentropic[t].fixed:
                    if not (
                        value(blk.efficiency_isentropic[t]) >= 0.05
                        and value(blk.efficiency_isentropic[t]) <= 1.0
                    ):
                        blk.efficiency_isentropic[t] = 0.8
                    blk.efficiency_isentropic[t].fix()
                    unfix_eff[t] = True

        if state_args is None:
            state_args = {}
            state_dict = cv.properties_in[t0].define_port_members()

            for k in state_dict.keys():
                if state_dict[k].is_indexed():
                    state_args[k] = {}
                    for m in state_dict[k].keys():
                        state_args[k][m] = state_dict[k][m].value
                else:
                    state_args[k] = state_dict[k].value

        # Get initialisation guesses for outlet and isentropic states
        for k in state_args:
            if k == "pressure" and k not in state_args_out:
                # Work out how to estimate outlet pressure
                if cv.properties_out[t0].pressure.fixed:
                    # Fixed outlet pressure, use this value
                    state_args_out[k] = value(cv.properties_out[t0].pressure)
                elif blk.deltaP[t0].fixed:
                    state_args_out[k] = value(state_args[k] + blk.deltaP[t0])
                elif blk.ratioP[t0].fixed:
                    state_args_out[k] = value(state_args[k] * blk.ratioP[t0])
                else:
                    # Not obvious what to do, use inlet state
                    state_args_out[k] = state_args[k]
            elif k not in state_args_out:
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
            blk.properties_isentropic[blk.flowsheet().time.first()].temperature,
            Var,
        ):
            blk.properties_isentropic[:].temperature.fix()
        elif isinstance(
            blk.properties_isentropic[blk.flowsheet().time.first()].enth_mol,
            Var,
        ):
            blk.properties_isentropic[:].enth_mol.fix()
        elif isinstance(
            blk.properties_isentropic[blk.flowsheet().time.first()].temperature,
            Expression,
        ):

            def tmp_rule(self, t):
                return (
                    blk.properties_isentropic[t].temperature
                    == blk.control_volume.properties_in[t].temperature
                )

            blk.tmp_init_constraint = Constraint(blk.flowsheet().time, rule=tmp_rule)

        blk.isentropic.deactivate()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high("Initialization Step 3 {}.".format(idaeslog.condition(res)))

        if isinstance(
            blk.properties_isentropic[blk.flowsheet().time.first()].temperature,
            Var,
        ):
            blk.properties_isentropic[:].temperature.unfix()
        elif isinstance(
            blk.properties_isentropic[blk.flowsheet().time.first()].enth_mol,
            Var,
        ):
            blk.properties_isentropic[:].enth_mol.unfix()
        elif isinstance(
            blk.properties_isentropic[blk.flowsheet().time.first()].temperature,
            Expression,
        ):
            blk.del_component(blk.tmp_init_constraint)

        blk.isentropic.activate()

        # ---------------------------------------------------------------------
        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high("Initialization Step 4 {}.".format(idaeslog.condition(res)))

        if activate_performance_curves:
            blk.performance_curve.activate()
            for t, v in unfix_eff.items():
                if v:
                    blk.efficiency_isentropic[t].unfix()
            for t, v in unfix_ratioP.items():
                if v:
                    blk.ratioP[t].unfix()
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                res = opt.solve(blk, tee=slc.tee)
            init_log.info_high(f"Initialization Step 5 {idaeslog.condition(res)}.")

        # ---------------------------------------------------------------------
        # Release Inlet state
        blk.control_volume.release_state(flags, outlvl)

        if not check_optimal_termination(res):
            raise InitializationError(
                f"{blk.name} failed to initialize successfully. Please check "
                f"the output logs for more information."
            )

        init_log.info(f"Initialization Complete: {idaeslog.condition(res)}")

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
            var_dict["Isentropic Efficiency"] = self.efficiency_isentropic[time_point]

        return {"vars": var_dict}

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        if hasattr(self, "work_fluid"):
            for t, v in self.work_fluid.items():
                iscale.set_scaling_factor(
                    v,
                    iscale.get_scaling_factor(
                        self.control_volume.work[t], default=1, warning=True
                    ),
                )

        if hasattr(self, "work_mechanical"):
            for t, v in self.work_mechanical.items():
                iscale.set_scaling_factor(
                    v,
                    iscale.get_scaling_factor(
                        self.control_volume.work[t], default=1, warning=True
                    ),
                )

        if hasattr(self, "work_isentropic"):
            for t, v in self.work_isentropic.items():
                iscale.set_scaling_factor(
                    v,
                    iscale.get_scaling_factor(
                        self.control_volume.work[t], default=1, warning=True
                    ),
                )

        if hasattr(self, "ratioP_calculation"):
            for t, c in self.ratioP_calculation.items():
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(
                        self.control_volume.properties_in[t].pressure,
                        default=1,
                        warning=True,
                    ),
                    overwrite=False,
                )

        if hasattr(self, "fluid_work_calculation"):
            for t, c in self.fluid_work_calculation.items():
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(
                        self.control_volume.deltaP[t], default=1, warning=True
                    ),
                    overwrite=False,
                )

        if hasattr(self, "actual_work"):
            for t, c in self.actual_work.items():
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(
                        self.control_volume.work[t], default=1, warning=True
                    ),
                    overwrite=False,
                )

        if hasattr(self, "isentropic_pressure"):
            for t, c in self.isentropic_pressure.items():
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(
                        self.control_volume.properties_in[t].pressure,
                        default=1,
                        warning=True,
                    ),
                    overwrite=False,
                )

        if hasattr(self, "isentropic"):
            for t, c in self.isentropic.items():
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(
                        self.control_volume.properties_in[t].entr_mol,
                        default=1,
                        warning=True,
                    ),
                    overwrite=False,
                )

        if hasattr(self, "isentropic_energy_balance"):
            for t, c in self.isentropic_energy_balance.items():
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(
                        self.control_volume.work[t], default=1, warning=True
                    ),
                    overwrite=False,
                )

        if hasattr(self, "zero_work_equation"):
            for t, c in self.zero_work_equation.items():
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(
                        self.control_volume.work[t], default=1, warning=True
                    ),
                )

        if hasattr(self, "state_material_balances"):
            cvol = self.control_volume
            phase_list = cvol.properties_in.phase_list
            phase_component_set = cvol.properties_in.phase_component_set
            mb_type = cvol._constructed_material_balance_type
            if mb_type == MaterialBalanceType.componentPhase:
                for (t, p, j), c in self.state_material_balances.items():
                    sf = iscale.get_scaling_factor(
                        cvol.properties_in[t].get_material_flow_terms(p, j),
                        default=1,
                        warning=True,
                    )
                    iscale.constraint_scaling_transform(c, sf)
            elif mb_type == MaterialBalanceType.componentTotal:
                for (t, j), c in self.state_material_balances.items():
                    sf = iscale.min_scaling_factor(
                        [
                            cvol.properties_in[t].get_material_flow_terms(p, j)
                            for p in phase_list
                            if (p, j) in phase_component_set
                        ]
                    )
                    iscale.constraint_scaling_transform(c, sf)
            else:
                # There are some other material balance types but they create
                # constraints with different names.
                _log.warning(f"Unknown material balance type {mb_type}")


@declare_process_block_class("Turbine", doc="Isentropic turbine model")
class TurbineData(PressureChangerData):
    # Pressure changer with isentropic turbine options
    CONFIG = PressureChangerData.CONFIG()
    CONFIG.compressor = False
    CONFIG.get("compressor")._default = False
    CONFIG.get("compressor")._domain = In([False])
    CONFIG.thermodynamic_assumption = ThermodynamicAssumption.isentropic
    CONFIG.get("thermodynamic_assumption")._default = ThermodynamicAssumption.isentropic


@declare_process_block_class("Compressor", doc="Isentropic compressor model")
class CompressorData(PressureChangerData):
    # Pressure changer with isentropic turbine options
    CONFIG = PressureChangerData.CONFIG()
    CONFIG.compressor = True
    CONFIG.get("compressor")._default = True
    CONFIG.get("compressor")._domain = In([True])
    CONFIG.thermodynamic_assumption = ThermodynamicAssumption.isentropic
    CONFIG.get("thermodynamic_assumption")._default = ThermodynamicAssumption.isentropic


@declare_process_block_class("Pump", doc="Pump model")
class PumpData(PressureChangerData):
    # Pressure changer with isentropic turbine options
    CONFIG = PressureChangerData.CONFIG()
    CONFIG.compressor = True
    CONFIG.get("compressor")._default = True
    CONFIG.get("compressor")._domain = In([True])
    CONFIG.thermodynamic_assumption = ThermodynamicAssumption.pump
    CONFIG.get("thermodynamic_assumption")._default = ThermodynamicAssumption.pump
