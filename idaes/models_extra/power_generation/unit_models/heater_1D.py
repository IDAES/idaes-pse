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
1-D Electric Trim Heater Model With Wall Temperatures

Discretization based on tube rows
"""
# Import Pyomo libraries
from pyomo.environ import (
    assert_optimal_termination,
    Block,
    units as pyunits,
    Var,
    value,
)
from pyomo.common.config import ConfigValue, In
from pyomo.util.calc_var_value import calculate_variable_from_constraint

# Import IDAES cores
from idaes.core import (
    ControlVolume1DBlock,
    declare_process_block_class,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    FlowDirection,
    UnitModelBlockData,
    useDefault,
)
import idaes.core.util.scaling as iscale
from idaes.core.solvers import get_solver
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.misc import add_object_reference
import idaes.logger as idaeslog
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models_extra.power_generation.unit_models.heat_exchanger_common import (
    make_geometry_common,
    make_performance_common,
    scale_common,
)
from idaes.core.initialization import SingleControlVolumeUnitInitializer

__author__ = "Jinliang Ma, Douglas Allan"


class Heater1DInitializer(SingleControlVolumeUnitInitializer):
    """
    Initializer for Heater 1D units.
    """

    def initialize_main_model(
        self,
        model: Block,
        copy_inlet_state: bool = False,
    ):
        """
        Initialization routine for the main Heater 1D model
        (as opposed to submodels like costing, which presently do not exist).

        Args:
            model: Pyomo Block to be initialized.
            copy_inlet_state: bool (default=False). Whether to copy inlet state to other states or not
                (0-D control volumes only). Copying will generally be faster, but inlet states may not contain
                all properties required elsewhere.
            duty: initial guess for heat duty to assist with initialization. Can be a Pyomo expression with units.

        Returns:
            Pyomo solver results object.

        """
        # Set solver options
        init_log = idaeslog.getInitLogger(
            model.name, self.get_output_level(), tag="unit"
        )
        solve_log = idaeslog.getSolveLogger(
            model.name, self.get_output_level(), tag="unit"
        )

        # Create solver
        solver_obj = get_solver(self.config.solver, self.config.solver_options)

        # ---------------------------------------------------------------------
        # Initialize shell block

        self.initialize_control_volume(model.control_volume, copy_inlet_state)

        init_log.info_high("Initialization Step 1 Complete.")

        calc_var = calculate_variable_from_constraint

        calc_var(model.length_flow_shell, model.length_flow_shell_eqn)
        calc_var(model.area_flow_shell, model.area_flow_shell_eqn)
        calc_var(model.area_flow_shell_min, model.area_flow_shell_min_eqn)

        for t in model.flowsheet().config.time:
            for x in model.control_volume.length_domain:
                model.control_volume.heat[t, x].fix(
                    value(model.electric_heat_duty[t] / model.length_flow_shell)
                )

        if model.config.has_pressure_change:
            model.control_volume.pressure.fix()

        model.control_volume.length.fix()
        assert degrees_of_freedom(model.control_volume) == 0
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = solver_obj.solve(model.control_volume, tee=slc.tee)

        assert_optimal_termination(res)

        init_log.info_high("Initialization Step 2 Complete.")
        model.control_volume.length.unfix()
        model.control_volume.heat.unfix()

        for t in model.flowsheet().config.time:
            for x in model.control_volume.length_domain:
                model.temp_wall_center[t, x].fix(
                    value(model.control_volume.properties[t, x].temperature) + 10
                )
                calc_var(model.heat_holdup[t, x], model.heat_holdup_eqn[t, x])
                model.temp_wall_center[t, x].unfix()

        if model.config.has_pressure_change:
            model.control_volume.pressure.unfix()
            model.control_volume.pressure[:, 0].fix()

        assert degrees_of_freedom(model) == 0

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = solver_obj.solve(model, tee=slc.tee)

        assert_optimal_termination(res)

        init_log.info_high("Initialization Step 3 Complete.")

        return res


@declare_process_block_class("Heater1D")
class Heater1DData(UnitModelBlockData):
    """Standard Trim Heater Model Class Class."""

    default_initializer = Heater1DInitializer

    CONFIG = UnitModelBlockData.CONFIG()
    CONFIG.declare(
        "has_fluid_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag",
            doc="""Indicates whether holdup terms for the fluid should be constructed or not.
            **default** - False.
            **Valid values:** {
            **False** - do not construct holdup terms}""",
        ),
    )
    CONFIG.declare(
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
    CONFIG.declare(
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
    CONFIG.declare(
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
    CONFIG.declare(
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
    CONFIG.declare(
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
    CONFIG.declare(
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
    CONFIG.declare(
        "transformation_method",
        ConfigValue(
            default=useDefault,
            description="Discretization method to use for DAE transformation",
            doc="""Discretization method to use for DAE transformation. See Pyomo
documentation for supported transformations.""",
        ),
    )
    CONFIG.declare(
        "transformation_scheme",
        ConfigValue(
            default=useDefault,
            description="Discretization scheme to use for DAE transformation",
            doc="""Discretization scheme to use when transforming domain. See Pyomo
documentation for supported schemes.""",
        ),
    )

    CONFIG = CONFIG

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
            doc="""If using collocation, number of collocation points to use
            per finite element when discretizing length domain (default=3)""",
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

    def build(self):
        """
        Begin building model (pre-DAE transformation).

        Args:
            None

        Returns:
            None
        """
        # Call UnitModel.build to setup dynamics
        super().build()

        # Set flow directions for the control volume blocks and specify
        # dicretization if not specified.
        set_direction_shell = FlowDirection.forward
        if self.config.transformation_method is useDefault:
            self.config.transformation_method = "dae.finite_difference"
        if self.config.transformation_scheme is useDefault:
            self.config.transformation_scheme = "BACKWARD"

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
            has_phase_equilibrium=False,
        )

        # Populate shell
        self.control_volume.add_material_balances(
            balance_type=self.config.material_balance_type,
            has_phase_equilibrium=False,
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
        make_geometry_common(self, shell_units=units)

        @self.Expression(
            doc="Common performance equations expect this expression to be here"
        )
        def length_flow_tube(b):
            return b.number_passes * b.length_tube_seg

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
            doc="Heat duty provided to heater through resistive heating",
        )
        units = self.config.property_package.get_metadata().derived_units
        make_performance_common(
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
                b.total_heat_transfer_coeff_shell[t, x]
                * b.total_heat_transfer_area
                * (
                    b.temp_wall_shell[t, x]
                    - b.control_volume.properties[t, x].temperature
                )
            )

        # Shell side wall temperature
        @self.Constraint(
            self.flowsheet().config.time,
            self.control_volume.length_domain,
            doc="shell side wall temperature",
        )
        def temp_wall_shell_eqn(b, t, x):
            return (
                b.total_heat_transfer_coeff_shell[t, x]
                * (
                    b.control_volume.properties[t, x].temperature
                    - b.temp_wall_shell[t, x]
                )
                # Divide thickness by 2 in order to represent center of hollow tube instead of
                # interior edge of hollow tube
                * (b.thickness_tube / (2 * b.therm_cond_wall) + b.rfouling_shell)
                == b.temp_wall_shell[t, x] - b.temp_wall_center[t, x]
            )

        @self.Constraint(
            self.flowsheet().config.time,
            self.control_volume.length_domain,
            doc="wall temperature",
        )
        def temp_wall_center_eqn(b, t, x):
            return heat_accumulation_term(b, t, x) == (
                -b.control_volume.heat[t, x]
                + b.electric_heat_duty[t] / b.length_flow_shell
            )

    def calculate_scaling_factors(self):
        def gsf(obj):
            return iscale.get_scaling_factor(obj, default=1, warning=True)

        def ssf(obj, sf):
            iscale.set_scaling_factor(obj, sf, overwrite=False)

        def cst(con, sf):
            iscale.constraint_scaling_transform(con, sf, overwrite=False)

        scale_common(
            self,
            self.control_volume,
            self.config.has_pressure_change,
            make_reynolds=True,
            make_nusselt=True,
        )

        sf_d_tube = iscale.get_scaling_factor(
            self.do_tube, default=1 / value(self.do_tube)
        )

        for t in self.flowsheet().time:
            for z in self.control_volume.length_domain:
                sf_hconv_conv = gsf(self.conv_heat_transfer_coeff_shell[t, z])
                cst(
                    self.conv_heat_transfer_coeff_shell_eqn[t, z],
                    sf_hconv_conv * sf_d_tube,
                )

                sf_T = gsf(self.control_volume.properties[t, z].temperature)
                ssf(self.temp_wall_shell[t, z], sf_T)
                ssf(self.temp_wall_center[t, z], sf_T)

                s_Q = gsf(self.control_volume.heat[t, z])
                ssf(self.electric_heat_duty[t], s_Q / value(self.length_flow_shell))
                cst(self.heat_shell_eqn[t, z], s_Q * value(self.length_flow_shell))
                ssf(self.temp_wall_center[t, z], sf_T)
                cst(self.temp_wall_shell_eqn[t, z], sf_T)
                cst(self.temp_wall_center_eqn[t, z], s_Q)

    def _get_performance_contents(self, time_point=0):
        var_dict = {}
        var_dict["Electric Heat Duty"] = self.electric_heat_duty[time_point]

        expr_dict = {}
        expr_dict["HX Area"] = self.total_heat_transfer_area

        return {"vars": var_dict, "exprs": expr_dict}

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {
                "Inlet": self.inlet,
                "Outlet": self.outlet,
            },
            time_point=time_point,
        )
