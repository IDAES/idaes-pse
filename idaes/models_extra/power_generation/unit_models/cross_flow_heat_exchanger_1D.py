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

import pyomo.common.config
import pyomo.opt

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
)
from pyomo.common.config import ConfigBlock, ConfigValue, In, Bool

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
from pyomo.network import Port
from idaes.core.solvers import get_solver
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.misc import add_object_reference
from idaes.core.util.exceptions import ConfigurationError, BurntToast
import idaes.logger as idaeslog
from idaes.core.util.tables import create_stream_table_dataframe

from idaes.models.unit_models.heater import _make_heater_config_block
from idaes.models.unit_models.heat_exchanger import (
    HeatExchangerFlowPattern,
    hx_process_config,
    add_hx_references,
)
from idaes.models.unit_models.heat_exchanger_1D import HeatExchanger1DData
from idaes.models.unit_models.shell_and_tube_1d import ShellAndTube1DData
import heat_exchanger_common

__author__ = "Jinliang Ma, Douglas Allan"

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("HeatExchangerCrossFlow1D")
class CrossFlowHeatExchanger1DData(HeatExchanger1DData):
    """Standard Heat Exchanger Cross Flow Unit Model Class."""

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
        # Add reference to control volume geometry
        add_object_reference(self, "area_flow_shell", shell.area)
        add_object_reference(self, "length_flow_shell", shell.length)
        add_object_reference(self, "area_flow_tube", tube.area)
        # total tube length of flow path
        add_object_reference(self, "length_flow_tube", tube.length)
        heat_exchanger_common._make_geometry_common(self, shell_units=shell_units)
        heat_exchanger_common._make_geometry_tube(self, shell_units=shell_units)

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
                "The HeatExchangerCrossFlow1D model is valid only for property packages "
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
                "The HeatExchangerCrossFlow1D model is valid only for property packages "
                "whose single phase is a vapor phase. The hot side phase is not a vapor phase."
            )
        if not pobj_cold.is_vapor_phase():
            raise ConfigurationError(
                "The HeatExchangerCrossFlow1D model is valid only for property packages "
                "whose single phase is a vapor phase. The cold side phase is not a vapor phase."
            )

        # Reference
        add_object_reference(self, "heat_tube", tube.heat)
        add_object_reference(self, "heat_shell", shell.heat)

        shell_has_pressure_change = False
        tube_has_pressure_change = False

        if self.config.cold_side.has_pressure_change:
            if self.config.shell_is_hot:
                add_object_reference(self, "deltaP_tube", tube.deltaP)
                tube_has_pressure_change = True
            else:
                add_object_reference(self, "deltaP_shell", shell.deltaP)
                shell_has_pressure_change = True
        if self.config.hot_side.has_pressure_change:
            if self.config.shell_is_hot:
                add_object_reference(self, "deltaP_shell", shell.deltaP)
                shell_has_pressure_change = True
            else:
                add_object_reference(self, "deltaP_tube", tube.deltaP)
                tube_has_pressure_change = True

        heat_exchanger_common._make_performance_common(
            self,
            shell=shell,
            shell_units=shell_units,
            shell_has_pressure_change=shell_has_pressure_change,
            make_reynolds=True,
            make_nusselt=True,
        )
        heat_exchanger_common._make_performance_tube(
            self,
            tube=tube,
            tube_units=tube_units,
            tube_has_pressure_change=tube_has_pressure_change,
            make_reynolds=True,
            make_nusselt=True,
        )

        def heat_accumulation_term(b, t, x):
            return b.heat_accumulation[t, x] if b.config.dynamic else 0

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
            return b.heat_tube[t, x] == b.hconv_tube[
                t, x
            ] * const.pi * b.di_tube * b.nrow_inlet * b.ncol_tube * (
                b.temp_wall_tube[t, x] - tube.properties[t, x].temperature
            )

        # Heat to wall per length on shell side
        @self.Constraint(
            self.flowsheet().config.time,
            shell.length_domain,
            doc="heat per length on shell side",
        )
        def heat_shell_eqn(b, t, x):
            return b.heat_shell[
                t, x
            ] * b.length_flow_shell == b.length_flow_tube * b.hconv_shell_total[
                t, x
            ] * const.pi * b.do_tube * b.nrow_inlet * b.ncol_tube * (
                b.temp_wall_shell[t, x] - shell.properties[t, x].temperature
            )

        # Tube side wall temperature
        # FIXME replace with deviation variables
        @self.Constraint(
            self.flowsheet().config.time,
            tube.length_domain,
            doc="tube side wall temperature",
        )
        def temp_wall_tube_eqn(b, t, x):
            return (
                b.hconv_tube[t, x]
                * (tube.properties[t, x].temperature - b.temp_wall_tube[t, x])
                * (b.thickness_tube / 2 / b.therm_cond_wall + b.rfouling_tube)
                == b.temp_wall_tube[t, x] - b.temp_wall_center[t, x]
            )

        # Shell side wall temperature
        @self.Constraint(
            self.flowsheet().config.time,
            shell.length_domain,
            doc="shell side wall temperature",
        )
        def temp_wall_shell_eqn(b, t, x):
            return (
                b.hconv_shell_total[t, x]
                * (shell.properties[t, x].temperature - b.temp_wall_shell[t, x])
                * (b.thickness_tube / 2 / b.therm_cond_wall + b.rfouling_shell)
                == b.temp_wall_shell[t, x] - b.temp_wall_center[t, x]
            )

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
                b.heat_shell[t, x] * b.length_flow_shell / b.length_flow_tube
                + b.heat_tube[t, x]
            )

        if not self.config.dynamic:
            z0 = shell.length_domain.first()
            z1 = shell.length_domain.last()

            @self.Expression(self.flowsheet().config.time)
            def total_heat_duty(b, t):
                enth_in = b.hot_side.properties[t, z0].get_enthalpy_flow_terms(p_hot)
                enth_out = b.hot_side.properties[t, z1].get_enthalpy_flow_terms(p_hot)

                return enth_out - enth_in

            @self.Expression(self.flowsheet().config.time)
            def log_mean_delta_temperature(b, t):
                dT0 = (
                    b.hot_side.properties[t, z0].temperature
                    - b.cold_side.properties[t, z0].temperature
                )
                dT1 = (
                    b.hot_side.properties[t, z1].temperature
                    - b.cold_side.properties[t, z1].temperature
                )
                return (dT0 - dT1) / log(dT0 / dT1)

            @self.Expression(self.flowsheet().config.time)
            def overall_heat_transfer_coefficient(b, t):
                return b.total_heat_duty[t] / (
                    b.total_heat_transfer_area * b.log_mean_delta_temperature[t]
                )

    def set_initial_condition(self):
        if self.config.dynamic is True:
            self.heat_accumulation[:, :].value = 0
            self.heat_accumulation[0, :].fix(0)
            # no accumulation term for fluid side models to avoid pressure waves

    def initialize_build(
        blk,
        shell_state_args=None,
        tube_state_args=None,
        outlvl=idaeslog.NOTSET,
        solver="ipopt",
        optarg=None,
    ):
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

        hot_side = blk.hot_side
        cold_side = blk.cold_side
        t0 = blk.flowsheet().time.first()
        if not (
            "temperature" in hot_side.properties[t0, 0].define_state_vars().keys()
            and "temperature" in cold_side.properties[t0, 0].define_state_vars().keys()
        ):
            raise NotImplementedError(
                "Presently, initialization of the HeatExchangerCrossFlow1D requires "
                "temperature to be a state variable of both hot side and cold side "
                "property packages. Extension to enth_mol or enth_mass as state variables "
                "is straightforward---feel free to open a pull request implementing it."
            )

        if blk.config.shell_is_hot:
            shell = blk.hot_side
            tube = blk.cold_side
            shell_has_pressure_change = blk.config.hot_side.has_pressure_change
            tube_has_pressure_change = blk.config.cold_side.has_pressure_change
        else:
            shell = blk.cold_side
            tube = blk.hot_side
            shell_has_pressure_change = blk.config.cold_side.has_pressure_change
            tube_has_pressure_change = blk.config.hot_side.has_pressure_change

        # ---------------------------------------------------------------------
        # Initialize shell block

        flags_tube = tube.initialize(
            outlvl=outlvl, optarg=optarg, solver=solver, state_args=tube_state_args
        )

        flags_shell = shell.initialize(
            outlvl=outlvl, optarg=optarg, solver=solver, state_args=shell_state_args
        )

        init_log.info_high("Initialization Step 1 Complete.")

        # Set tube thermal conductivity to a small value to avoid IPOPT unable to solve initially
        therm_cond_wall_save = blk.therm_cond_wall.value
        blk.therm_cond_wall = 0.05
        # In Step 2, fix tube metal temperatures fix fluid state variables (enthalpy/temperature and pressure)
        # calculate maximum heat duty assuming infinite area and use half of the maximum duty as initial guess to calculate outlet temperature

        for t in blk.flowsheet().config.time:
            # TODO we first access cp during initialization. That could pose a problem if it is 
            # converted to a Var-Constraint pair instead of being a giant Expression like it is
            # presently.
            mcp_hot_side = value(
                hot_side.properties[t, 0].flow_mol * hot_side.properties[t, 0].cp_mol
            )
            T_in_hot_side = value(hot_side.properties[t, 0].temperature)
            P_in_hot_side = value(hot_side.properties[t, 0].pressure)
            if blk.config.flow_type == HeatExchangerFlowPattern.cocurrent:
                mcp_cold_side = value(
                    cold_side.properties[t, 0].flow_mol * cold_side.properties[t, 0].cp_mol
                )
                T_in_cold_side = value(cold_side.properties[t, 0].temperature)
                P_in_cold_side = value(cold_side.properties[t, 0].pressure)

                T_out_max = (
                    mcp_cold_side * T_in_cold_side
                    + mcp_hot_side * T_in_hot_side
                ) / (mcp_cold_side + mcp_hot_side)

                q_guess = mcp_cold_side * (T_out_max - T_in_cold_side) / 2

                temp_out_cold_side_guess = (
                    T_in_cold_side + q_guess / mcp_cold_side
                )
                cold_side.properties[t, 1].temperature.fix(temp_out_cold_side_guess)

                temp_out_hot_side_guess = (
                    T_in_cold_side - q_guess / mcp_hot_side
                )
                hot_side.properties[t, 1].temperature.fix(temp_out_hot_side_guess)

            elif blk.config.flow_type == HeatExchangerFlowPattern.countercurrent:
                mcp_cold_side = value(
                    cold_side.properties[t, 1].flow_mol * cold_side.properties[t, 1].cp_mol
                )
                T_in_cold_side = value(cold_side.properties[t, 1].temperature)
                P_in_cold_side = value(cold_side.properties[t, 1].pressure)

                if mcp_cold_side < mcp_hot_side:
                    q_guess = mcp_cold_side * (T_in_hot_side - T_in_cold_side) / 2
                else:
                    q_guess = mcp_hot_side * (T_in_hot_side - T_in_cold_side) / 2
                
                temp_out_cold_side_guess = (
                    T_in_cold_side + q_guess / mcp_cold_side
                )
                cold_side.properties[t, 0].temperature.fix(temp_out_cold_side_guess)

                temp_out_hot_side_guess = (
                    T_in_hot_side - q_guess / mcp_hot_side
                )
                hot_side.properties[t, 1].temperature.fix(temp_out_hot_side_guess)

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
                blk.temp_wall_center[t, z].fix(
                    value(hot_side.properties[t, z].temperature + cold_side.properties[t, z].temperature) / 2
                )

                blk.temp_wall_shell[t, z].set_value(blk.temp_wall_center[t, z].value)
                blk.temp_wall_tube[t, z].set_value(blk.temp_wall_center[t, z].value)

                if blk.config.cold_side.has_pressure_change:
                    cold_side.properties[t, z].pressure.fix(P_in_cold_side)
                if blk.config.hot_side.has_pressure_change:
                    hot_side.properties[t, z].pressure.fix(P_in_hot_side)

        blk.temp_wall_center_eqn.deactivate()
        if tube_has_pressure_change == True:
            blk.deltaP_tube_eqn.deactivate()
        if shell_has_pressure_change == True:
            blk.deltaP_shell_eqn.deactivate()
        blk.heat_tube_eqn.deactivate()
        blk.heat_shell_eqn.deactivate()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        pyomo.opt.assert_optimal_termination(res)
        init_log.info_high("Initialization Step 2 {}.".format(idaeslog.condition(res)))

        # In Step 3, unfix fluid state variables (enthalpy/temperature and pressure)
        # keep the inlet state variables fixed, otherwise, the degree of freedom > 0
        hot_side.properties[:, :].temperature.unfix()
        hot_side.properties[:, :].pressure.unfix()
        hot_side.properties[:, 0].temperature.fix()
        hot_side.properties[:, 0].pressure.fix()

        cold_side.properties[:, :].temperature.unfix()
        cold_side.properties[:, :].pressure.unfix()
        if blk.config.flow_type == HeatExchangerFlowPattern.cocurrent:
            cold_side.properties[:, 0].temperature.fix()
            cold_side.properties[:, 0].pressure.fix()
        elif blk.config.flow_type == HeatExchangerFlowPattern.countercurrent:
            cold_side.properties[:, 1].temperature.fix()
            cold_side.properties[:, 1].pressure.fix()
        else:
            raise BurntToast(
                "HeatExchangerFlowPattern should be limited to cocurrent "
                "or countercurrent flow by parent model. Please open an "
                "issue on the IDAES Github so this error can be fixed."
            )

        if tube_has_pressure_change == True:
            blk.deltaP_tube_eqn.activate()
        if shell_has_pressure_change == True:
            blk.deltaP_shell_eqn.activate()
        blk.heat_tube_eqn.activate()
        blk.heat_shell_eqn.activate()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        pyomo.opt.assert_optimal_termination(res)

        init_log.info_high("Initialization Step 3 {}.".format(idaeslog.condition(res)))

        blk.temp_wall_center.unfix()
        blk.temp_wall_center_eqn.activate()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        pyomo.opt.assert_optimal_termination(res)

        init_log.info_high("Initialization Step 4 {}.".format(idaeslog.condition(res)))

        # set the wall thermal conductivity back to the user specified value
        blk.therm_cond_wall = therm_cond_wall_save

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high("Initialization Step 5 {}.".format(idaeslog.condition(res)))
        tube.release_state(flags_tube)
        shell.release_state(flags_shell)
        init_log.info("Initialization Complete.")

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
        else:
            shell = self.cold_side
            tube = self.hot_side

        tube_has_pressure_change = hasattr(self, "deltaP_tube")
        shell_has_pressure_change = hasattr(self, "deltaP_shell")

        heat_exchanger_common._scale_common(
            self,
            shell,
            shell_has_pressure_change,
            make_reynolds=True,
            make_nusselt=True
        )
        heat_exchanger_common._scale_tube(
            self,
            tube,
            tube_has_pressure_change,
            make_reynolds=True,
            make_nusselt=True
        )

        for t in self.flowsheet().time:
            for z in shell.length_domain:
                # FIXME try to do this rigorously later on
                sf_area_per_length_tube = 1 / value(
                    const.pi * self.di_tube * self.nrow_inlet * self.ncol_tube
                )
                sf_T_tube = gsf(tube.properties[t, z].temperature)
                ssf(self.temp_wall_tube[t, z], sf_T_tube)
                cst(self.temp_wall_tube_eqn[t, z], sf_T_tube)

                sf_Q_tube = gsf(tube.heat[t, z])
                cst(self.heat_tube_eqn[t, z], sf_Q_tube)

                sf_T_shell = gsf(shell.properties[t, z].temperature)
                ssf(self.temp_wall_shell[t, z], sf_T_shell)
                cst(self.temp_wall_shell_eqn[t, z], sf_T_shell)

                sf_area_per_length_shell = value(
                    self.length_flow_shell
                    / (
                        self.length_flow_tube
                        * const.pi
                        * self.do_tube
                        * self.nrow_inlet
                        * self.ncol_tube
                    )
                )
                sf_hconv_shell_conv = gsf(self.hconv_shell_conv[t, z])
                s_Q_shell = sgsf(
                    shell.heat[t, z],
                    sf_hconv_shell_conv * sf_area_per_length_shell * sf_T_shell,
                )
                cst(self.heat_shell_eqn[t, z], s_Q_shell * value(self.length_flow_shell))
                # Geometric mean is overkill for most reasonable cases, but it mitigates
                # damage done when one stream has an unset scaling factor
                ssf(self.temp_wall_center[t, z], (sf_T_shell * sf_T_tube)**0.5)
                cst(self.temp_wall_center_eqn[t, z], (sf_Q_tube * s_Q_shell) ** 0.5)

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
        expr_dict["Delta T Driving"] = self.log_mean_delta_temperature[time_point]
        expr_dict["Total Heat Duty"] = self.total_heat_duty[time_point]
        expr_dict["HX Coefficient"] = self.overall_heat_transfer_coefficient[time_point]

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
