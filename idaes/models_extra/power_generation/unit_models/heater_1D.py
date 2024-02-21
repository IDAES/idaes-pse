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
1-D Cross Flow Heat Exchanger Model With Wall Temperatures

Discretization based on tube rows
"""
from __future__ import division

# Import Python libraries
import math

# Import Pyomo libraries
from pyomo.environ import (
    SolverFactory,
    Var,
    Param,
    Constraint,
    value,
    TerminationCondition,
    exp,
    sqrt,
    log,
    sin,
    cos,
    SolverStatus,
    units as pyunits,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In

# Import IDAES cores
from idaes.core import (
    ControlVolume1DBlock,
    UnitModelBlockData,
    declare_process_block_class,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    FlowDirection,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.util.constants import Constants as const
import idaes.core.util.scaling as iscale
from pyomo.dae import DerivativeVar
from idaes.core.solvers import get_solver
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.misc import add_object_reference
import idaes.logger as idaeslog
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.model_statistics import degrees_of_freedom

from heat_exchanger_common import _make_geometry_common, _make_performance_common, _scale_common

__author__ = "Jinliang Ma, Douglas Allan"

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("Heater1D")
class Heater1DData(UnitModelBlockData):
    """Standard Heat Exchanger Cross Flow Unit Model Class."""

    # Template for config arguments for shell and tube side
    _SideTemplate = ConfigBlock()
    _SideTemplate.declare(
        "dynamic",
        ConfigValue(
            default=useDefault,
            domain=In([useDefault, True, False]),
            description="Dynamic model flag",
            doc="""Indicates whether this model will be dynamic or not,
**default** = useDefault.
**Valid values:** {
**useDefault** - get flag from parent (default = False),
**True** - set as a dynamic model,
**False** - set as a steady-state model.}""",
        ),
    )
    _SideTemplate.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([True, False]),
            description="Holdup construction flag",
            doc="""Indicates whether holdup terms should be constructed or not.
Must be True if dynamic = True,
**default** - False.
**Valid values:** {
**True** - construct holdup terms,
**False** - do not construct holdup terms}""",
        ),
    )
    _SideTemplate.declare(
        "has_fluid_holdup",
        ConfigValue(
            default=False,
            domain=In([True, False]),
            description="Holdup construction flag",
            doc="""Indicates whether holdup terms for the fluid should be constructed or not.
            **default** - False.
            **Valid values:** {
            **True** - construct holdup terms,
            **False** - do not construct holdup terms}""",
        ),
    )
    _SideTemplate.declare(
        "material_balance_type",
        ConfigValue(
            default=MaterialBalanceType.componentTotal,
            domain=In(MaterialBalanceType),
            description="Material balance construction flag",
            doc="""Indicates what type of mass balance should be constructed,
**default** - MaterialBalanceType.componentTotal.
**Valid values:** {
**MaterialBalanceType.none** - exclude material balances,
**MaterialBalanceType.componentPhase** - use phase component balances,
**MaterialBalanceType.componentTotal** - use total component balances,
**MaterialBalanceType.elementTotal** - use total element balances,
**MaterialBalanceType.total** - use total material balance.}""",
        ),
    )
    _SideTemplate.declare(
        "energy_balance_type",
        ConfigValue(
            default=EnergyBalanceType.enthalpyTotal,
            domain=In(EnergyBalanceType),
            description="Energy balance construction flag",
            doc="""Indicates what type of energy balance should be constructed,
**default** - EnergyBalanceType.enthalpyTotal.
**Valid values:** {
**EnergyBalanceType.none** - exclude energy balances,
**EnergyBalanceType.enthalpyTotal** - single enthalpy balance for material,
**EnergyBalanceType.enthalpyPhase** - enthalpy balances for each phase,
**EnergyBalanceType.energyTotal** - single energy balance for material,
**EnergyBalanceType.energyPhase** - energy balances for each phase.}""",
        ),
    )
    _SideTemplate.declare(
        "momentum_balance_type",
        ConfigValue(
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
**MomentumBalanceType.momentumPhase** - momentum balances for each phase.}""",
        ),
    )
    _SideTemplate.declare(
        "has_pressure_change",
        ConfigValue(
            default=False,
            domain=In([True, False]),
            description="Pressure change term construction flag",
            doc="""Indicates whether terms for pressure change should be
constructed,
**default** - False.
**Valid values:** {
**True** - include pressure change terms,
**False** - exclude pressure change terms.}""",
        ),
    )
    _SideTemplate.declare(
        "has_phase_equilibrium",
        ConfigValue(
            default=False,
            domain=In([True, False]),
            description="Phase equilibrium term construction flag",
            doc="""Argument to enable phase equilibrium on the shell side.
- True - include phase equilibrium term
- False - do not include phase equilibrium term""",
        ),
    )
    _SideTemplate.declare(
        "property_package",
        ConfigValue(
            domain=is_physical_parameter_block,
            description="Property package to use for control volume",
            doc="""Property parameter object used to define property calculations
(default = 'use_parent_value')
- 'use_parent_value' - get package from parent (default = None)
- a ParameterBlock object""",
        ),
    )
    _SideTemplate.declare(
        "property_package_args",
        ConfigValue(
            default=None,
            description="Arguments for constructing shell property package",
            doc="""A dict of arguments to be passed to the PropertyBlockData
and used when constructing these
(default = 'use_parent_value')
- 'use_parent_value' - get package from parent (default = None)
- a dict (see property package for documentation)""",
        ),
    )
    # TODO : We should probably think about adding a consistency check for the
    # TODO : discretization methods as well.
    _SideTemplate.declare(
        "transformation_method",
        ConfigValue(
            default=useDefault,
            description="Discretization method to use for DAE transformation",
            doc="""Discretization method to use for DAE transformation. See Pyomo
documentation for supported transformations.""",
        ),
    )
    _SideTemplate.declare(
        "transformation_scheme",
        ConfigValue(
            default=useDefault,
            description="Discretization scheme to use for DAE transformation",
            doc="""Discretization scheme to use when transforming domain. See Pyomo
documentation for supported schemes.""",
        ),
    )

    CONFIG = _SideTemplate

    # Common config args for both sides
    CONFIG.declare(
        "finite_elements",
        ConfigValue(
            default=5,
            domain=int,
            description="Number of finite elements length domain",
            doc="""Number of finite elements to use when discretizing length
domain (default=5). Should set to the number of tube rows""",
        ),
    )
    CONFIG.declare(
        "collocation_points",
        ConfigValue(
            default=3,
            domain=int,
            description="Number of collocation points per finite element",
            doc="""Number of collocation points to use per finite element when
discretizing length domain (default=3)""",
        ),
    )
    CONFIG.declare(
        "tube_arrangement",
        ConfigValue(
            default="in-line",
            domain=In(["in-line", "staggered"]),
            description="tube configuration",
            doc="tube arrangement could be in-line or staggered",
        ),
    )
    CONFIG.declare(
        "has_radiation",
        ConfigValue(
            default=False,
            domain=In([False, True]),
            description="Has side 2 gas radiation",
            doc="define if shell side gas radiation is to be considered",
        ),
    )

    def build(self):
        """
        Begin building model (pre-DAE transformation).

        Args:
            None

        Returns:
            None
        """
        # Call UnitModel.build to setup dynamics
        super(Heater1DData, self).build()

        # Set flow directions for the control volume blocks and specify
        # dicretization if not specified.
        set_direction_shell = FlowDirection.forward
        if self.config.transformation_method is useDefault:
            self.config.transformation_method = "dae.finite_difference"
        if self.config.transformation_scheme is useDefault:
            self.config.transformation_scheme = "FORWARD"

        if self.config.property_package_args is None:
            self.config.property_package_args = {}

        # Control volume 1D for shell, set to steady-state for fluid
        self.control_volume = ControlVolume1DBlock(
            dynamic=self.config.dynamic and self.config.has_fluid_holdup,
            has_holdup=self.config.has_fluid_holdup,
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
            transformation_method=self.config.transformation_method,
            transformation_scheme=self.config.transformation_scheme,
            finite_elements=self.config.finite_elements,
            collocation_points=self.config.collocation_points,
        )

        self.control_volume.add_geometry(flow_direction=set_direction_shell)

        self.control_volume.add_state_blocks(
            information_flow=set_direction_shell,
            has_phase_equilibrium=self.config.has_phase_equilibrium,
        )

        # Populate shell
        self.control_volume.add_material_balances(
            balance_type=self.config.material_balance_type,
            has_phase_equilibrium=self.config.has_phase_equilibrium,
        )

        self.control_volume.add_energy_balances(
            balance_type=self.config.energy_balance_type, has_heat_transfer=True
        )

        self.control_volume.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=self.config.has_pressure_change,
        )

        self.control_volume.apply_transformation()

        # Populate tube

        # Add Ports for shell side
        self.add_inlet_port(name="inlet", block=self.control_volume)
        self.add_outlet_port(name="outlet", block=self.control_volume)

        self._make_geometry()

        self._make_performance()

    def _make_geometry(self):
        """
        Constraints for unit model.

        Args:
            None

        Returns:
            None
        """
        units = self.config.property_package.get_metadata().derived_units
        # Add reference to control volume geometry
        add_object_reference(self, "area_flow_shell", self.control_volume.area)
        add_object_reference(self, "length_flow_shell", self.control_volume.length)
        _make_geometry_common(self, shell_units=units)
        @self.Expression(
            doc="Common performance equations expect this expression to be here"
        )
        def length_flow_tube(b):
            return b.nseg_tube * b.length_tube_seg


    def _make_performance(self):
        """
        Constraints for unit model.

        Args:
            None

        Returns:
            None
        """
        self.electric_heat_duty = Var(
            self.flowsheet().config.time,
            initialize=1e6,
            units=pyunits.W,
            doc="Heat duty provided to heater " "through resistive heating",
        )
        units = self.config.property_package.get_metadata().derived_units
        _make_performance_common(
            self,
            shell=self.control_volume,
            shell_units=units,
            shell_has_pressure_change=self.config.has_pressure_change,
            make_reynolds=True,
            make_nusselt=True,
        )

        def heat_accumulation_term(b, t, x):
            return b.heat_accumulation[t, x] if b.config.dynamic else 0

        # Nusselt number, currently assume Re>300
        @self.Constraint(
            self.flowsheet().config.time,
            self.control_volume.length_domain,
            doc="Nusselts number equation",
        )
        def N_Nu_shell_eqn(b, t, x):
            return (
                b.N_Nu_shell[t, x]
                == b.f_arrangement
                * 0.33
                * b.N_Re_shell[t, x] ** 0.6
                * b.control_volume.properties[t, x].prandtl_number_phase["Vap"]
                ** 0.333333
            )

        # Energy balance with tube wall
        # ------------------------------------
        # Heat to wall per length
        @self.Constraint(
            self.flowsheet().config.time,
            self.control_volume.length_domain,
            doc="heat per length",
        )
        def heat_shell_eqn(b, t, x):
            return b.control_volume.heat[t, x] * b.length_flow_shell == (
                b.hconv_shell_total[t, x]
                * b.total_heat_transfer_area
                * (b.temp_wall_shell[t, x] - b.control_volume.properties[t, x].temperature)
            )

        # Shell side wall temperature
        @self.Constraint(
            self.flowsheet().config.time,
            self.control_volume.length_domain,
            doc="shell side wall temperature",
        )
        def temp_wall_shell_eqn(b, t, x):
            return (
                b.hconv_shell_total[t, x]
                * (b.control_volume.properties[t, x].temperature - b.temp_wall_shell[t, x])
                * (b.thickness_tube / b.therm_cond_wall + b.rfouling_shell)
                == b.temp_wall_shell[t, x] - b.temp_wall_center[t, x]
            )

        @self.Constraint(
            self.flowsheet().config.time,
            self.control_volume.length_domain,
            doc="wall temperature",
        )
        def temp_wall_center_eqn(b, t, x):
            return heat_accumulation_term(b, t, x) == (
                -b.control_volume.heat[t, x] + b.electric_heat_duty[t] / b.length_flow_shell
            )

    def set_initial_condition(self):
        if self.config.dynamic is True:
            self.heat_accumulation[:, :].value = 0
            self.heat_accumulation[0, :].fix(0)
            # no accumulation term for fluid side models to avoid pressure waves

    def initialize_build(blk, state_args=None, outlvl=0, solver="ipopt", optarg=None):
        """
        HeatExchangerCrossFlow1D initialization routine

        Keyword Arguments:
            state_args : a dict of arguments to be passed to the property
                         package(s) to provide an initial state for
                         initialization (see documentation of the specific
                         property package) (default = None).
            outlvl : sets output level of initialization routine

                     * 0 = no output (default)
                     * 1 = return solver state for each step in routine
                     * 2 = return solver state for each step in subroutines
                     * 3 = include solver output information (tee=True)

            optarg : solver options dictionary object (default={'tol': 1e-6})
            solver : str indicating which solver to use during
                     initialization (default = 'ipopt')

        Returns:
            None
        """
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")

        if optarg is None:
            optarg = {}
        opt = get_solver(solver, optarg)

        # ---------------------------------------------------------------------
        # Initialize shell block

        flags = blk.control_volume.initialize(
            outlvl=0, optarg=optarg, solver=solver, state_args=state_args
        )

        init_log.info_high("Initialization Step 1 Complete.")

        # mcp = value(blk.control_volume.properties[0,0].flow_mol * blk.control_volume.properties[0,0].cp_mol)
        # tout_guess = value(blk.tube.properties[0,0].temperature) + value(blk.electric_heat_duty[0]/blk.length_flow)
        calc_var = calculate_variable_from_constraint

        calc_var(blk.length_flow_shell, blk.length_flow_shell_eqn)
        calc_var(blk.area_flow_shell, blk.area_flow_shell_eqn)
        calc_var(blk.area_flow_shell_min, blk.area_flow_shell_min_eqn)

        for t in blk.flowsheet().config.time:
            for x in blk.control_volume.length_domain:
                blk.control_volume.heat[t, x].fix(
                    value(blk.electric_heat_duty[t] / blk.length_flow_shell)
                )
        blk.control_volume.length.fix()
        assert degrees_of_freedom(blk.control_volume) == 0
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk.control_volume, tee=slc.tee)

        assert res.solver.termination_condition == TerminationCondition.optimal
        assert res.solver.status == SolverStatus.ok

        init_log.info_high("Initialization Step 2 Complete.")
        blk.control_volume.length.unfix()
        blk.control_volume.heat.unfix()

        for t in blk.flowsheet().config.time:
            for x in blk.control_volume.length_domain:
                blk.temp_wall_center[t, x].fix(
                    value(blk.control_volume.properties[t, x].temperature) + 10
                )
                calc_var(blk.heat_holdup[t, x], blk.heat_holdup_eqn[t, x])
                blk.temp_wall_center[t, x].unfix()

                # fixed = blk.control_volume.temperature[t,x].fixed
                # blk.control_volume.temperature[t,x].fix()
                # calc_var()

        assert degrees_of_freedom(blk) == 0

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)

        assert res.solver.termination_condition == TerminationCondition.optimal
        assert res.solver.status == SolverStatus.ok

        init_log.info_high("Initialization Step 3 Complete.")

        blk.control_volume.release_state(flags)

    def calculate_scaling_factors(self):
        def gsf(obj):
            return iscale.get_scaling_factor(obj, default=1, warning=True)

        def ssf(obj, sf):
            iscale.set_scaling_factor(obj, sf, overwrite=False)

        def cst(con, sf):
            iscale.constraint_scaling_transform(con, sf, overwrite=False)

        sgsf = iscale.set_and_get_scaling_factor

        _scale_common(
            self,
            self.control_volume,
            self.config.has_pressure_change,
            make_reynolds=True,
            make_nusselt=True
        )

        sf_d_tube = iscale.get_scaling_factor(
            self.do_tube, default=1 / value(self.do_tube)
        )

        for t in self.flowsheet().time:
            for z in self.control_volume.length_domain:
                sf_hconv_conv = gsf(self.hconv_shell_conv[t, z])
                cst(self.hconv_shell_conv_eqn[t, z], sf_hconv_conv * sf_d_tube)

                if self.config.has_radiation:
                    sf_hconv_rad = 1  # FIXME Placeholder
                    sf_hconv_total = 1 / (1 / sf_hconv_conv + 1 / sf_hconv_rad)
                else:
                    sf_hconv_total = sf_hconv_conv

                # FIXME try to do this rigorously later on
                sf_T = gsf(self.control_volume.properties[t, z].temperature)
                ssf(self.temp_wall_shell[t, z], sf_T)
                ssf(self.temp_wall_center[t, z], sf_T)

                sf_area_per_length = value(
                    self.length_flow_shell / self.total_heat_transfer_area
                )
                s_Q = sgsf(
                    self.control_volume.heat[t, z],
                    sf_hconv_total * sf_area_per_length * sf_T,
                )
                ssf(self.electric_heat_duty[t], s_Q / value(self.length_flow_shell))
                cst(self.heat_shell_eqn[t, z], s_Q * value(self.length_flow_shell))
                ssf(self.temp_wall_center[t, z], sf_T)
                cst(self.temp_wall_shell_eqn[t, z], sf_T)
                cst(self.temp_wall_center_eqn[t, z], s_Q)

    def _get_performance_contents(self, time_point=0):
        var_dict = {}
        # var_dict = {
        #     "HX Coefficient": self.overall_heat_transfer_coefficient[time_point]
        # }
        # var_dict["HX Area"] = self.area
        # var_dict["Heat Duty"] = self.heat_duty[time_point]
        # if self.config.flow_pattern == HeatExchangerFlowPattern.crossflow:
        #     var_dict = {"Crossflow Factor": self.crossflow_factor[time_point]}

        expr_dict = {}
        expr_dict["HX Area"] = self.total_heat_transfer_area
        expr_dict["Electric Heat Duty"] = self.electric_heat_duty[time_point]

        return {"vars": var_dict, "exprs": expr_dict}

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {
                "Inlet": self.inlet,
                "Outlet": self.outlet,
            },
            time_point=time_point,
        )
