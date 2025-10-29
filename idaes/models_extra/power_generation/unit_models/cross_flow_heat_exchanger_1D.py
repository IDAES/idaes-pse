#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
1-D Cross Flow Heat Exchanger Model With Wall Temperatures

Discretization based on tube rows
"""


# Import Pyomo libraries
from pyomo.environ import (
    assert_optimal_termination,
    Block,
    value,
    Var,
    log,
    Param,
    Reference,
    units as pyunits,
)
from pyomo.common.config import ConfigValue, In, Bool
from pyomo.network import Port
from pyomo.dae.flatten import slice_component_along_sets

# Import IDAES cores
from idaes.core import declare_process_block_class
from idaes.core.util.constants import Constants as const
import idaes.core.util.scaling as iscale
from idaes.core.solvers import get_solver
from idaes.core.util.exceptions import InitializationError
from idaes.core.util.misc import add_object_reference
from idaes.core.util.exceptions import ConfigurationError, BurntToast
import idaes.logger as idaeslog
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.models.unit_models.heat_exchanger import (
    HeatExchangerFlowPattern,
)
from idaes.models.unit_models.heat_exchanger_1D import HeatExchanger1DData
from idaes.models_extra.power_generation.unit_models import heat_exchanger_common
from idaes.core.initialization import SingleControlVolumeUnitInitializer


__author__ = "Jinliang Ma, Douglas Allan"


class CrossFlowHeatExchanger1DInitializer(SingleControlVolumeUnitInitializer):
    """
    Initializer for Cross Flow Heat Exchanger 1D units.

    First, the shell and tube control volumes are initialized without heat transfer. Next
    the total possible heat transfer between streams is estimated based on heat capacity,
    flow rate, and inlet/outlet temperatures. The actual temperature change is set to be
    half the theoretical maximum, and the shell and tube are initialized with linear
    temperature profiles. Finally, temperatures besides the inlets are unfixed and
    the performance equations are activated before a full solve of the system model.
    """

    def initialize_main_model(
        self,
        model: Block,
        copy_inlet_state: bool = False,
    ):
        """
        Initialization routine for the main Cross Flow Heat Exchanger 1D model
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

        solver_obj = get_solver(self.config.solver, self.config.solver_options)

        hot_side = model.hot_side
        cold_side = model.cold_side
        t0 = model.flowsheet().time.first()
        if not (
            "temperature" in hot_side.properties[t0, 0].define_state_vars().keys()
            and "temperature" in cold_side.properties[t0, 0].define_state_vars().keys()
        ):
            raise NotImplementedError(
                "Presently, initialization of the CrossFlowHeatExchanger1D requires "
                "temperature to be a state variable of both hot side and cold side "
                "property packages. Extension to enth_mol or enth_mass as state variables "
                "is straightforward---feel free to open a pull request implementing it."
            )

        hot_units = model.config.hot_side.property_package.get_metadata().derived_units
        cold_units = (
            model.config.cold_side.property_package.get_metadata().derived_units
        )

        if model.config.shell_is_hot:
            shell = model.hot_side
            tube = model.cold_side
            shell_has_pressure_change = model.config.hot_side.has_pressure_change
            tube_has_pressure_change = model.config.cold_side.has_pressure_change
            shell_units = (
                model.config.hot_side.property_package.get_metadata().derived_units
            )
            tube_units = (
                model.config.cold_side.property_package.get_metadata().derived_units
            )
        else:
            shell = model.cold_side
            tube = model.hot_side
            shell_has_pressure_change = model.config.cold_side.has_pressure_change
            tube_has_pressure_change = model.config.hot_side.has_pressure_change
            shell_units = (
                model.config.cold_side.property_package.get_metadata().derived_units
            )
            tube_units = (
                model.config.hot_side.property_package.get_metadata().derived_units
            )

        # Trigger creation of cp for use in future initialization
        # Important to do before initializing property packages in
        # case it is implemented as Var-Constraint pair instead of
        # an Expression
        hot_side.properties[t0, 0].cp_mol  # pylint: disable=pointless-statement
        cold_side.properties[t0, 0].cp_mol  # pylint: disable=pointless-statement

        # ---------------------------------------------------------------------
        # Initialize shell block
        self.initialize_control_volume(tube, copy_inlet_state)
        self.initialize_control_volume(shell, copy_inlet_state)

        init_log.info_high("Initialization Step 1 Complete.")

        # Set tube thermal conductivity to a small value to avoid IPOPT unable to solve initially
        therm_cond_wall_save = model.therm_cond_wall.value
        model.therm_cond_wall.set_value(0.05)
        # In Step 2, fix tube metal temperatures fix fluid state variables (enthalpy/temperature and pressure)
        # calculate maximum heat duty assuming infinite area and use half of the maximum duty as initial guess to calculate outlet temperature

        for t in model.flowsheet().config.time:
            mcp_hot_side = value(
                pyunits.convert(
                    hot_side.properties[t, 0].flow_mol
                    * hot_side.properties[t, 0].cp_mol,
                    to_units=shell_units["power"] / shell_units["temperature"],
                )
            )
            T_in_hot_side = value(
                pyunits.convert(
                    hot_side.properties[t, 0].temperature,
                    to_units=shell_units["temperature"],
                )
            )
            P_in_hot_side = value(hot_side.properties[t, 0].pressure)
            if model.config.flow_type == HeatExchangerFlowPattern.cocurrent:
                mcp_cold_side = value(
                    pyunits.convert(
                        cold_side.properties[t, 0].flow_mol
                        * cold_side.properties[t, 0].cp_mol,
                        to_units=shell_units["power"] / shell_units["temperature"],
                    )
                )
                T_in_cold_side = value(
                    pyunits.convert(
                        cold_side.properties[t, 0].temperature,
                        to_units=shell_units["temperature"],
                    )
                )
                P_in_cold_side = value(cold_side.properties[t, 0].pressure)

                T_out_max = (
                    mcp_cold_side * T_in_cold_side + mcp_hot_side * T_in_hot_side
                ) / (mcp_cold_side + mcp_hot_side)

                q_guess = mcp_cold_side * (T_out_max - T_in_cold_side) / 2

                temp_out_cold_side_guess = T_in_cold_side + q_guess / mcp_cold_side

                cold_side.properties[t, 1].temperature.fix(
                    pyunits.convert_value(
                        temp_out_cold_side_guess,
                        from_units=shell_units["temperature"],
                        to_units=cold_units["temperature"],
                    )
                )

                temp_out_hot_side_guess = T_in_cold_side - q_guess / mcp_hot_side
                hot_side.properties[t, 1].temperature.fix(
                    pyunits.convert_value(
                        temp_out_hot_side_guess,
                        from_units=shell_units["temperature"],
                        to_units=hot_units["temperature"],
                    )
                )

            elif model.config.flow_type == HeatExchangerFlowPattern.countercurrent:
                mcp_cold_side = value(
                    pyunits.convert(
                        cold_side.properties[t, 1].flow_mol
                        * cold_side.properties[t, 1].cp_mol,
                        to_units=shell_units["power"] / shell_units["temperature"],
                    )
                )
                T_in_cold_side = value(
                    pyunits.convert(
                        cold_side.properties[t, 1].temperature,
                        to_units=shell_units["temperature"],
                    )
                )
                P_in_cold_side = value(cold_side.properties[t, 1].pressure)

                if mcp_cold_side < mcp_hot_side:
                    q_guess = mcp_cold_side * (T_in_hot_side - T_in_cold_side) / 2
                else:
                    q_guess = mcp_hot_side * (T_in_hot_side - T_in_cold_side) / 2

                temp_out_cold_side_guess = T_in_cold_side + q_guess / mcp_cold_side
                cold_side.properties[t, 0].temperature.fix(
                    pyunits.convert_value(
                        temp_out_cold_side_guess,
                        from_units=shell_units["temperature"],
                        to_units=cold_units["temperature"],
                    )
                )

                temp_out_hot_side_guess = T_in_hot_side - q_guess / mcp_hot_side
                hot_side.properties[t, 1].temperature.fix(
                    pyunits.convert_value(
                        temp_out_hot_side_guess,
                        from_units=shell_units["temperature"],
                        to_units=hot_units["temperature"],
                    )
                )

            else:
                raise BurntToast(
                    "HeatExchangerFlowPattern should be limited to cocurrent "
                    "or countercurrent flow by parent model. Please open an "
                    "issue on the IDAES Github so this error can be fixed."
                )

            for z in cold_side.length_domain:
                hot_side.properties[t, z].temperature.fix(
                    value(
                        (1 - z) * hot_side.properties[t, 0].temperature
                        + z * hot_side.properties[t, 1].temperature
                    )
                )
                cold_side.properties[t, z].temperature.fix(
                    value(
                        (1 - z) * cold_side.properties[t, 0].temperature
                        + z * cold_side.properties[t, 1].temperature
                    )
                )
                model.temp_wall_center[t, z].fix(
                    value(
                        pyunits.convert(
                            hot_side.properties[t, z].temperature,
                            to_units=shell_units["temperature"],
                        )
                        + pyunits.convert(
                            cold_side.properties[t, z].temperature,
                            to_units=shell_units["temperature"],
                        )
                    )
                    / 2
                )

                model.temp_wall_shell[t, z].set_value(
                    model.temp_wall_center[t, z].value
                )
                model.temp_wall_tube[t, z].set_value(
                    pyunits.convert_value(
                        model.temp_wall_center[t, z].value,
                        from_units=shell_units["temperature"],
                        to_units=tube_units["temperature"],
                    )
                )

                if model.config.cold_side.has_pressure_change:
                    cold_side.properties[t, z].pressure.fix(P_in_cold_side)
                if model.config.hot_side.has_pressure_change:
                    hot_side.properties[t, z].pressure.fix(P_in_hot_side)

        model.temp_wall_center_eqn.deactivate()
        if tube_has_pressure_change:
            model.deltaP_tube_eqn.deactivate()
        if shell_has_pressure_change:
            model.deltaP_shell_eqn.deactivate()
        model.heat_tube_eqn.deactivate()
        model.heat_shell_eqn.deactivate()

        if (
            str.upper(model.config.hot_side.transformation_scheme)
            == "LAGRANGE-LEGENDRE"
        ):
            model.lagrange_legendre_deactivation()
            for t in model.flowsheet().time:
                for x in model.hot_side.length_domain.get_finite_elements():
                    hot_side.properties[t, x].temperature.unfix()
                    cold_side.properties[t, x].temperature.unfix()
            model.fix_initialization_states()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = solver_obj.solve(model, tee=slc.tee)
        try:
            assert_optimal_termination(res)
        except AssertionError:
            raise InitializationError("Initialization solve failed.")
        init_log.info_high("Initialization Step 2 {}.".format(idaeslog.condition(res)))

        # In Step 3, unfix fluid state variables (enthalpy/temperature and pressure)
        # keep the inlet state variables fixed, otherwise, the degree of freedom > 0
        hot_side.properties[:, :].temperature.unfix()
        hot_side.properties[:, :].pressure.unfix()
        hot_side.properties[:, 0].temperature.fix()
        hot_side.properties[:, 0].pressure.fix()

        cold_side.properties[:, :].temperature.unfix()
        cold_side.properties[:, :].pressure.unfix()
        if model.config.flow_type == HeatExchangerFlowPattern.cocurrent:
            cold_side.properties[:, 0].temperature.fix()
            cold_side.properties[:, 0].pressure.fix()
        elif model.config.flow_type == HeatExchangerFlowPattern.countercurrent:
            cold_side.properties[:, 1].temperature.fix()
            cold_side.properties[:, 1].pressure.fix()
        else:
            raise BurntToast(
                "HeatExchangerFlowPattern should be limited to cocurrent "
                "or countercurrent flow by parent model. Please open an "
                "issue on the IDAES Github so this error can be fixed."
            )

        if tube_has_pressure_change:
            model.deltaP_tube_eqn.activate()
        if shell_has_pressure_change:
            model.deltaP_shell_eqn.activate()

        model.hot_side.pressure_dx_disc_eq.activate()
        model.cold_side.pressure_dx_disc_eq.activate()
        model.heat_tube_eqn.activate()
        model.heat_shell_eqn.activate()

        if (
            str.upper(model.config.hot_side.transformation_scheme)
            == "LAGRANGE-LEGENDRE"
        ):
            model.lagrange_legendre_deactivation()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = solver_obj.solve(model, tee=slc.tee)
        try:
            assert_optimal_termination(res)
        except AssertionError:
            raise InitializationError("Initialization solve failed.")

        init_log.info_high("Initialization Step 3 {}.".format(idaeslog.condition(res)))

        model.temp_wall_center.unfix()
        model.temp_wall_center_eqn.activate()

        if (
            str.upper(model.config.hot_side.transformation_scheme)
            == "LAGRANGE-LEGENDRE"
        ):
            model.lagrange_legendre_deactivation()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = solver_obj.solve(model, tee=slc.tee)
        try:
            assert_optimal_termination(res)
        except AssertionError:
            raise InitializationError("Initialization solve failed.")

        init_log.info_high("Initialization Step 4 {}.".format(idaeslog.condition(res)))

        # set the wall thermal conductivity back to the user specified value
        model.therm_cond_wall.set_value(therm_cond_wall_save)

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = solver_obj.solve(model, tee=slc.tee)
        init_log.info_high("Initialization Step 5 {}.".format(idaeslog.condition(res)))
        init_log.info("Initialization Complete.")
        return res


@declare_process_block_class("CrossFlowHeatExchanger1D")
class CrossFlowHeatExchanger1DData(HeatExchanger1DData):
    """Standard Cross Flow Heat Exchanger Unit Model Class."""

    default_initializer = CrossFlowHeatExchanger1DInitializer

    CONFIG = HeatExchanger1DData.CONFIG()
    CONFIG.declare(
        "shell_is_hot",
        ConfigValue(
            default=True,
            domain=Bool,
            description="Shell side contains hot fluid",
            doc="""Boolean flag indicating whether shell side contains hot fluid (default=True).
        If True, shell side will be the hot_side, if False shell side will be cold_side.""",
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

    def _process_config(self):
        # Copy and pasted from ShellAndTube1D
        super()._process_config()

        # Check for custom names, and if not present assign defaults
        if self.config.hot_side_name is None:
            if self.config.shell_is_hot:
                self.config.hot_side_name = "Shell"
            else:
                self.config.hot_side_name = "Tube"

        if self.config.cold_side_name is None:
            if self.config.shell_is_hot:
                self.config.cold_side_name = "Tube"
            else:
                self.config.cold_side_name = "Shell"

    def build(self):
        """
        Begin building model (pre-DAE transformation).

        Args:
            None

        Returns:
            None
        """

        # Call HeatExchanger1DData build to make common components
        super().build()

        # The HeatExchanger1DData equates the heat lost by the hot side and heat gained by the cold side.
        # That equation is deleted here because heat can accumulate in the wall.
        self.del_component(self.heat_conservation)

        if (
            str.upper(self.config.hot_side.transformation_scheme) == "LAGRANGE-LEGENDRE"
            or str.upper(self.config.cold_side.transformation_scheme)
            == "LAGRANGE-LEGENDRE"
        ) and (
            self.config.hot_side.has_pressure_change
            or self.config.cold_side.has_pressure_change
        ):
            raise NotImplementedError(
                "Pressure change is not implemented for Lagrange-Legendre collocation."
            )

        # Create aliases for ports
        if self.config.shell_is_hot:
            self.shell_inlet = Port(extends=self.hot_side_inlet)
            self.shell_outlet = Port(extends=self.hot_side_outlet)
            self.tube_inlet = Port(extends=self.cold_side_inlet)
            self.tube_outlet = Port(extends=self.cold_side_outlet)
        else:
            self.shell_inlet = Port(extends=self.cold_side_inlet)
            self.shell_outlet = Port(extends=self.cold_side_outlet)
            self.tube_inlet = Port(extends=self.hot_side_inlet)
            self.tube_outlet = Port(extends=self.hot_side_outlet)

    def _make_geometry(self):
        """
        Constraints for unit model.

        Args:
            None

        Returns:
            None
        """
        if self.config.shell_is_hot:
            shell = self.hot_side
            tube = self.cold_side
            shell_units = (
                self.config.hot_side.property_package.get_metadata().derived_units
            )
        else:
            shell = self.cold_side
            tube = self.hot_side
            shell_units = (
                self.config.cold_side.property_package.get_metadata().derived_units
            )
        # Use add_object_reference for scalar vars
        add_object_reference(self, "area_flow_shell", shell.area)
        add_object_reference(self, "area_flow_tube", tube.area)
        add_object_reference(self, "length_flow_shell", shell.length)
        add_object_reference(self, "length_flow_tube", tube.length)

        heat_exchanger_common.make_geometry_common(self, shell_units=shell_units)

        # Important that these values about tube geometry are in shell units!
        @self.Constraint(doc="Length of tube side flow")
        def length_flow_tube_eqn(b):
            return (
                pyunits.convert(b.length_flow_tube, to_units=shell_units["length"])
                == b.number_passes * b.length_tube_seg
            )

        @self.Constraint(doc="Total area of tube flow")
        def area_flow_tube_eqn(b):
            return (
                b.area_flow_tube
                == 0.25
                * const.pi
                * b.di_tube**2.0
                * b.number_columns_per_pass
                * b.number_rows_per_pass
            )

    def _make_performance(self):
        """
        Constraints for unit model.

        Args:
            None

        Returns:
            None
        """
        if self.config.shell_is_hot:
            shell = self.hot_side
            tube = self.cold_side
            shell_units = (
                self.config.hot_side.property_package.get_metadata().derived_units
            )
            tube_units = (
                self.config.cold_side.property_package.get_metadata().derived_units
            )
        else:
            shell = self.cold_side
            tube = self.hot_side
            shell_units = (
                self.config.cold_side.property_package.get_metadata().derived_units
            )
            tube_units = (
                self.config.hot_side.property_package.get_metadata().derived_units
            )
        if (
            len(self.config.hot_side.property_package.phase_list) != 1
            or len(self.config.cold_side.property_package.phase_list) != 1
        ):
            raise ConfigurationError(
                "The CrossFlowHeatExchanger1D model is valid only for property packages "
                f"with a single phase. Found {len(self.config.hot_side.property_package.phase_list)} "
                f"phases on the hot side and {len(self.config.cold_side.property_package.phase_list)} "
                "phases on the cold side."
            )

        p_hot = self.config.hot_side.property_package.phase_list.at(1)
        pobj_hot = self.config.hot_side.property_package.get_phase(p_hot)
        p_cold = self.config.cold_side.property_package.phase_list.at(1)
        pobj_cold = self.config.cold_side.property_package.get_phase(p_cold)
        if not pobj_hot.is_vapor_phase():
            raise ConfigurationError(
                "The CrossFlowHeatExchanger1D model is valid only for property packages "
                "whose single phase is a vapor phase. The hot side phase is not a vapor phase."
            )
        if not pobj_cold.is_vapor_phase():
            raise ConfigurationError(
                "The CrossFlowHeatExchanger1D model is valid only for property packages "
                "whose single phase is a vapor phase. The cold side phase is not a vapor phase."
            )

        self.heat_tube = Reference(tube.heat)
        self.heat_shell = Reference(shell.heat)

        shell_has_pressure_change = False
        tube_has_pressure_change = False

        if self.config.cold_side.has_pressure_change:
            if self.config.shell_is_hot:
                self.deltaP_tube = Reference(tube.deltaP)
                tube_has_pressure_change = True
            else:
                self.deltaP_shell = Reference(shell.deltaP)
                shell_has_pressure_change = True
        if self.config.hot_side.has_pressure_change:
            if self.config.shell_is_hot:
                self.deltaP_shell = Reference(shell.deltaP)
                shell_has_pressure_change = True
            else:
                self.deltaP_tube = Reference(tube.deltaP)
                tube_has_pressure_change = True

        heat_exchanger_common.make_performance_common(
            self,
            shell=shell,
            shell_units=shell_units,
            shell_has_pressure_change=shell_has_pressure_change,
            make_reynolds=True,
            make_nusselt=True,
        )

        self.heat_transfer_coeff_tube = Var(
            self.flowsheet().config.time,
            tube.length_domain,
            initialize=100.0,
            units=tube_units["heat_transfer_coefficient"],
            doc="tube side convective heat transfer coefficient",
        )

        # Heat transfer resistance due to the fouling on tube side
        self.rfouling_tube = Param(
            initialize=0.0,
            mutable=True,
            units=1 / tube_units["heat_transfer_coefficient"],
            doc="fouling resistance on tube side",
        )
        # Correction factor for convective heat transfer coefficient on tube side
        self.fcorrection_htc_tube = Var(
            initialize=1.0, doc="correction factor for convective HTC on tube side"
        )
        # Correction factor for tube side pressure drop due to friction
        if tube_has_pressure_change:
            # Loss coefficient for a 180 degree bend (u-turn), usually related to radius to inside diameter ratio
            self.kloss_uturn = Param(
                initialize=0.5, mutable=True, doc="loss coefficient of a tube u-turn"
            )
            self.fcorrection_dp_tube = Var(
                initialize=1.0, doc="correction factor for tube side pressure drop"
            )

        # Boundary wall temperature on tube side
        self.temp_wall_tube = Var(
            self.flowsheet().config.time,
            tube.length_domain,
            initialize=500,
            units=tube_units[
                "temperature"
            ],  # Want to be in shell units for consistency in equations
            doc="boundary wall temperature on tube side",
        )
        # Tube side heat transfer coefficient and pressure drop
        # -----------------------------------------------------
        # Velocity on tube side
        self.v_tube = Var(
            self.flowsheet().config.time,
            tube.length_domain,
            initialize=1.0,
            units=tube_units["velocity"],
            doc="velocity on tube side",
        )

        # Reynalds number on tube side
        self.N_Re_tube = Var(
            self.flowsheet().config.time,
            tube.length_domain,
            initialize=10000.0,
            units=pyunits.dimensionless,
            doc="Reynolds number on tube side",
            bounds=(1e-7, None),
        )

        if tube_has_pressure_change:
            # Friction factor on tube side
            self.friction_factor_tube = Var(
                self.flowsheet().config.time,
                tube.length_domain,
                initialize=1.0,
                doc="friction factor on tube side",
            )

            # Pressure drop due to friction on tube side
            self.deltaP_tube_friction = Var(
                self.flowsheet().config.time,
                tube.length_domain,
                initialize=-10.0,
                units=tube_units["pressure"] / tube_units["length"],
                doc="pressure drop due to friction on tube side",
            )

            # Pressure drop due to 180 degree turn on tube side
            self.deltaP_tube_uturn = Var(
                self.flowsheet().config.time,
                tube.length_domain,
                initialize=-10.0,
                units=tube_units["pressure"] / tube_units["length"],
                doc="pressure drop due to u-turn on tube side",
            )
        # Nusselt number on tube side
        self.N_Nu_tube = Var(
            self.flowsheet().config.time,
            tube.length_domain,
            initialize=1,
            units=pyunits.dimensionless,
            doc="Nusselts number on tube side",
            bounds=(1e-7, None),
        )

        # Velocity equation
        @self.Constraint(
            self.flowsheet().config.time,
            tube.length_domain,
            doc="tube side velocity equation",
        )
        def v_tube_eqn(b, t, x):
            return (
                b.v_tube[t, x]
                * pyunits.convert(b.area_flow_tube, to_units=tube_units["area"])
                * tube.properties[t, x].dens_mol_phase["Vap"]
                == tube.properties[t, x].flow_mol
            )

        # Reynolds number
        @self.Constraint(
            self.flowsheet().config.time,
            tube.length_domain,
            doc="Reynolds number equation on tube side",
        )
        def N_Re_tube_eqn(b, t, x):
            return (
                b.N_Re_tube[t, x] * tube.properties[t, x].visc_d_phase["Vap"]
                == pyunits.convert(b.di_tube, to_units=tube_units["length"])
                * b.v_tube[t, x]
                * tube.properties[t, x].dens_mol_phase["Vap"]
                * tube.properties[t, x].mw
            )

        if tube_has_pressure_change:
            # Friction factor
            @self.Constraint(
                self.flowsheet().config.time,
                tube.length_domain,
                doc="Darcy friction factor on tube side",
            )
            def friction_factor_tube_eqn(b, t, x):
                return (
                    b.friction_factor_tube[t, x] * b.N_Re_tube[t, x] ** 0.25
                    == 0.3164 * b.fcorrection_dp_tube
                )

            # Pressure drop due to friction per tube length
            @self.Constraint(
                self.flowsheet().config.time,
                tube.length_domain,
                doc="pressure drop due to friction per tube length",
            )
            def deltaP_tube_friction_eqn(b, t, x):
                return (
                    b.deltaP_tube_friction[t, x]
                    * pyunits.convert(b.di_tube, to_units=tube_units["length"])
                    == -0.5
                    * tube.properties[t, x].dens_mass_phase["Vap"]
                    * b.v_tube[t, x] ** 2
                    * b.friction_factor_tube[t, x]
                )

            # Pressure drop due to u-turn
            @self.Constraint(
                self.flowsheet().config.time,
                tube.length_domain,
                doc="pressure drop due to u-turn on tube side",
            )
            def deltaP_tube_uturn_eqn(b, t, x):
                return (
                    b.deltaP_tube_uturn[t, x]
                    * pyunits.convert(b.length_tube_seg, to_units=tube_units["length"])
                    == -0.5
                    * tube.properties[t, x].dens_mass_phase["Vap"]
                    * b.v_tube[t, x] ** 2
                    * b.kloss_uturn
                )

            # Total pressure drop on tube side
            @self.Constraint(
                self.flowsheet().config.time,
                tube.length_domain,
                doc="total pressure drop on tube side",
            )
            def deltaP_tube_eqn(b, t, x):
                return b.deltaP_tube[t, x] == (
                    b.deltaP_tube_friction[t, x] + b.deltaP_tube_uturn[t, x]
                )

        @self.Constraint(
            self.flowsheet().config.time,
            tube.length_domain,
            doc="convective heat transfer coefficient equation on tube side",
        )
        def heat_transfer_coeff_tube_eqn(b, t, x):
            return (
                b.heat_transfer_coeff_tube[t, x] * b.di_tube
                == b.N_Nu_tube[t, x]
                * tube.properties[t, x].therm_cond_phase["Vap"]
                * b.fcorrection_htc_tube
            )

        # Nusselts number
        @self.Constraint(
            self.flowsheet().config.time,
            tube.length_domain,
            doc="Nusselts number equation on tube side",
        )
        def N_Nu_tube_eqn(b, t, x):
            return (
                b.N_Nu_tube[t, x]
                == 0.023
                * b.N_Re_tube[t, x] ** 0.8
                * tube.properties[t, x].prandtl_number_phase["Vap"] ** 0.4
            )

        @self.Constraint(
            self.flowsheet().config.time,
            shell.length_domain,
            doc="Nusselts number equation on shell side",
        )
        def N_Nu_shell_eqn(b, t, x):
            return (
                b.N_Nu_shell[t, x]
                == b.f_arrangement
                * 0.33
                * b.N_Re_shell[t, x] ** 0.6
                * shell.properties[t, x].prandtl_number_phase["Vap"] ** 0.333333
            )

        # Energy balance with tube wall
        # ------------------------------------
        # Heat to wall per length on tube side
        @self.Constraint(
            self.flowsheet().config.time,
            tube.length_domain,
            doc="heat per length on tube side",
        )
        def heat_tube_eqn(b, t, x):
            return b.heat_tube[t, x] == (
                b.heat_transfer_coeff_tube[t, x]
                * const.pi
                * pyunits.convert(b.di_tube, to_units=tube_units["length"])
                * b.number_rows_per_pass
                * b.number_columns_per_pass
                * (b.temp_wall_tube[t, x] - tube.properties[t, x].temperature)
            )

        # Heat to wall per length on shell side
        @self.Constraint(
            self.flowsheet().config.time,
            shell.length_domain,
            doc="heat per length on shell side",
        )
        def heat_shell_eqn(b, t, x):
            return b.heat_shell[t, x] * b.length_flow_shell == pyunits.convert(
                b.length_flow_tube, to_units=shell_units["length"]
            ) * b.total_heat_transfer_coeff_shell[
                t, x
            ] * const.pi * b.do_tube * b.number_rows_per_pass * b.number_columns_per_pass * (
                b.temp_wall_shell[t, x] - shell.properties[t, x].temperature
            )

        # Tube side wall temperature
        @self.Constraint(
            self.flowsheet().config.time,
            tube.length_domain,
            doc="tube side wall temperature",
        )
        def temp_wall_tube_eqn(b, t, x):
            return b.heat_transfer_coeff_tube[t, x] * (
                tube.properties[t, x].temperature - b.temp_wall_tube[t, x]
            ) * (
                pyunits.convert(
                    b.thickness_tube / 2 / b.therm_cond_wall,
                    to_units=1 / tube_units["heat_transfer_coefficient"],
                )
                + b.rfouling_tube
            ) == b.temp_wall_tube[
                t, x
            ] - pyunits.convert(
                b.temp_wall_center[t, x], to_units=tube_units["temperature"]
            )

        # Shell side wall temperature
        @self.Constraint(
            self.flowsheet().config.time,
            shell.length_domain,
            doc="shell side wall temperature",
        )
        def temp_wall_shell_eqn(b, t, x):
            return (
                b.total_heat_transfer_coeff_shell[t, x]
                * (shell.properties[t, x].temperature - b.temp_wall_shell[t, x])
                * (b.thickness_tube / 2 / b.therm_cond_wall + b.rfouling_shell)
                == b.temp_wall_shell[t, x] - b.temp_wall_center[t, x]
            )

        def heat_accumulation_term(b, t, x):
            return b.heat_accumulation[t, x] if b.config.dynamic else 0

        # Center point wall temperature based on energy balance for tube wall heat holdup
        @self.Constraint(
            self.flowsheet().config.time,
            tube.length_domain,
            doc="center point wall temperature",
        )
        def temp_wall_center_eqn(b, t, x):
            # heat_shell and heat_tube are positive when heat flows into those
            # control volumes (and out of the wall), hence the negative sign
            # on heat_accumulation_term
            return -heat_accumulation_term(b, t, x) == (
                b.heat_shell[t, x]
                + pyunits.convert(
                    b.heat_tube[t, x],
                    to_units=shell_units["power"] / shell_units["length"],
                )
                * pyunits.convert(b.length_flow_tube, to_units=shell_units["length"])
                / b.length_flow_shell
            )

        if not self.config.dynamic:
            z0 = shell.length_domain.first()
            z1 = shell.length_domain.last()

            @self.Expression(self.flowsheet().config.time)
            def total_heat_duty(b, t):
                enth_in = b.hot_side.properties[t, z0].get_enthalpy_flow_terms(p_hot)
                enth_out = b.hot_side.properties[t, z1].get_enthalpy_flow_terms(p_hot)

                return pyunits.convert(
                    enth_in - enth_out,  # Hot side loses enthalpy
                    to_units=shell_units["power"],  # Hot side isn't always the shell
                )

            @self.Expression(self.flowsheet().config.time)
            def log_mean_delta_temperature(b, t):
                dT0 = b.hot_side.properties[t, z0].temperature - pyunits.convert(
                    b.cold_side.properties[t, z0].temperature,
                    to_units=shell_units["temperature"],
                )
                dT1 = b.hot_side.properties[t, z1].temperature - pyunits.convert(
                    b.cold_side.properties[t, z1].temperature,
                    to_units=shell_units["temperature"],
                )
                return (dT0 - dT1) / log(dT0 / dT1)

            @self.Expression(self.flowsheet().config.time)
            def overall_heat_transfer_coefficient(b, t):
                return b.total_heat_duty[t] / (
                    b.total_heat_transfer_area * b.log_mean_delta_temperature[t]
                )

    def lagrange_legendre_deactivation(self):
        """
        Deactivates variables and fixes constraints at finite element
        boundaries. This process is necessary so that the system of
        equations is not structurally singular if Lagrange-Legendre
        collocation is used.
        """

        element_bounds = self.cold_side.length_domain.get_finite_elements()

        # This section should be expanded into a method on the
        # ControlVolume1D, but that requires including many more
        # variables and constraints
        for side in [self.cold_side, self.hot_side]:
            material_balances_dict = dict(
                (key, Reference(slc))
                for key, slc in slice_component_along_sets(
                    side.material_balances, (side.length_domain,)
                )
            )
            for x in element_bounds:
                for con in material_balances_dict.values():
                    try:
                        con[x].deactivate()
                    except KeyError:
                        # Constraint.skip
                        pass
                if (
                    x != self.cold_side.length_domain.first()
                    and x != self.cold_side.length_domain.last()
                ):
                    pass

                for t in self.flowsheet().time:
                    if (
                        x != self.cold_side.length_domain.first()
                        and x != self.cold_side.length_domain.last()
                    ):
                        # Need to calculate properties at the first
                        # and last elements in order to communicate
                        # state variables to ports
                        pass
                    try:
                        side.enthalpy_balances[t, x].deactivate()
                    except KeyError:
                        pass

    def calculate_scaling_factors(self):
        def gsf(obj):
            return iscale.get_scaling_factor(obj, default=1, warning=True)

        def ssf(obj, sf):
            iscale.set_scaling_factor(obj, sf, overwrite=False)

        def cst(con, sf):
            iscale.constraint_scaling_transform(con, sf, overwrite=False)

        sgsf = iscale.set_and_get_scaling_factor
        if self.config.shell_is_hot:
            shell = self.hot_side
            tube = self.cold_side
            shell_units = (
                self.config.hot_side.property_package.get_metadata().derived_units
            )
        else:
            shell = self.cold_side
            tube = self.hot_side
            shell_units = (
                self.config.cold_side.property_package.get_metadata().derived_units
            )

        shell_has_pressure_change = hasattr(self, "deltaP_shell")

        heat_exchanger_common.scale_common(
            self,
            shell,
            shell_has_pressure_change,
            make_reynolds=True,
            make_nusselt=True,
        )

        sf_di_tube = iscale.get_scaling_factor(
            self.do_tube, default=1 / value(self.di_tube)
        )

        sf_area_per_length_shell = value(
            self.length_flow_shell
            / (
                pyunits.convert(self.length_flow_tube, to_units=shell_units["length"])
                * const.pi
                * self.do_tube
                * self.number_rows_per_pass
                * self.number_columns_per_pass
            )
        )

        for t in self.flowsheet().time:
            for z in shell.length_domain:
                # FIXME get better scaling later
                ssf(self.v_tube[t, z], 1 / 20)
                sf_flow_mol_tube = gsf(tube.properties[t, z].flow_mol)

                cst(self.v_tube_eqn[t, z], sf_flow_mol_tube)

                # FIXME should get scaling of N_Re from defining eqn
                sf_N_Re_tube = sgsf(self.N_Re_tube[t, z], 1e-4)

                sf_visc_d_tube = gsf(tube.properties[t, z].visc_d_phase["Vap"])
                cst(self.N_Re_tube_eqn[t, z], sf_N_Re_tube * sf_visc_d_tube)

                sf_k_tube = gsf(tube.properties[t, z].therm_cond_phase["Vap"])

                sf_N_Nu_tube = sgsf(self.N_Nu_tube[t, z], 1 / 0.023 * sf_N_Re_tube**0.8)
                cst(self.N_Nu_tube_eqn[t, z], sf_N_Nu_tube)

                sf_heat_transfer_coeff_tube = sgsf(
                    self.heat_transfer_coeff_tube[t, z],
                    sf_N_Nu_tube * sf_k_tube / sf_di_tube,
                )
                cst(
                    self.heat_transfer_coeff_tube_eqn[t, z],
                    sf_heat_transfer_coeff_tube * sf_di_tube,
                )

                sf_T_tube = gsf(tube.properties[t, z].temperature)
                ssf(self.temp_wall_tube[t, z], sf_T_tube)
                cst(self.temp_wall_tube_eqn[t, z], sf_T_tube)

                sf_Q_tube = gsf(tube.heat[t, z])
                cst(self.heat_tube_eqn[t, z], sf_Q_tube)

                sf_T_shell = gsf(shell.properties[t, z].temperature)
                ssf(self.temp_wall_shell[t, z], sf_T_shell)
                cst(self.temp_wall_shell_eqn[t, z], sf_T_shell)

                sf_conv_heat_transfer_coeff_shell = gsf(
                    self.conv_heat_transfer_coeff_shell[t, z]
                )
                s_Q_shell = sgsf(
                    shell.heat[t, z],
                    sf_conv_heat_transfer_coeff_shell
                    * sf_area_per_length_shell
                    * sf_T_shell,
                )
                cst(
                    self.heat_shell_eqn[t, z], s_Q_shell * value(self.length_flow_shell)
                )
                # Geometric mean is overkill for most reasonable cases, but it mitigates
                # damage done when one stream has an unset scaling factor
                ssf(self.temp_wall_center[t, z], (sf_T_shell * sf_T_tube) ** 0.5)
                cst(self.temp_wall_center_eqn[t, z], (sf_Q_tube * s_Q_shell) ** 0.5)

    def _get_performance_contents(self, time_point=0):
        var_dict = {}

        expr_dict = {}
        expr_dict["HX Area"] = self.total_heat_transfer_area
        expr_dict["Delta T Driving"] = self.log_mean_delta_temperature[time_point]
        expr_dict["Total Heat Duty"] = self.total_heat_duty[time_point]
        expr_dict["Average HX Coefficient"] = self.overall_heat_transfer_coefficient[
            time_point
        ]

        return {"vars": var_dict, "exprs": expr_dict}

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {
                "Hot Inlet": self.shell_inlet,
                "Hot Outlet": self.shell_outlet,
                "Cold Inlet": self.tube_inlet,
                "Cold Outlet": self.tube_outlet,
            },
            time_point=time_point,
        )
