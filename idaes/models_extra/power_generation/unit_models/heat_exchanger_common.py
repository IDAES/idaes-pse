#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Methods shared by the CrossFlowHeatExchanger1D and Heater1D models.
"""
# Presently the tube methods are not shared. I'm not sure why I chose to extract
# them here, but I don't want to change things this late into development. --Doug

# Import Pyomo libraries
from pyomo.environ import (
    Var,
    Param,
    value,
    units as pyunits,
)
from pyomo.util.calc_var_value import calculate_variable_from_constraint

from pyomo.dae import DerivativeVar

# Import IDAES cores
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.constants import Constants as const
import idaes.core.util.scaling as iscale
from idaes.core.util.misc import add_object_reference

__author__ = "Jinliang Ma, Douglas Allan"


def make_geometry_common(blk, shell_units):
    """Function to create variables, constraints, and expressions regarding
    unit geometry shared between the CrossFlowHeatExchanger1D and Heater1D models.

    Args:
        blk : unit model on which components are being generated
        shell_units : derived units for property package of shell control volume

    Returns:
        None
    """
    # Number of tube columns in the cross section plane perpendicular to shell side fluid flow (y direction)
    blk.number_columns_per_pass = Var(
        initialize=10.0, doc="Number of tube columns", units=pyunits.dimensionless
    )

    # Number of segments of tube bundles
    blk.number_passes = Var(
        initialize=10.0,
        doc="Number of segments of tube bundles",
        units=pyunits.dimensionless,
    )

    # Number of inlet tube rows
    blk.number_rows_per_pass = Var(
        initialize=1, doc="Number of inlet tube rows", units=pyunits.dimensionless
    )

    # Inner diameter of tubes
    blk.di_tube = Var(
        initialize=0.05, doc="Inner diameter of tube", units=shell_units["length"]
    )

    blk.thickness_tube = Var(
        initialize=0.005, doc="Tube thickness", units=shell_units["length"]
    )

    blk.pitch_y = Var(
        initialize=0.1,
        doc="Pitch between tubes perpendicular to direction of flow",
        units=shell_units["length"],
    )

    blk.pitch_x = Var(
        initialize=0.1,
        doc="Pitch between tubes in line with direction of flow",
        units=shell_units["length"],
    )

    # Length of tube per segment in z direction
    blk.length_tube_seg = Var(
        initialize=1.0, doc="Length of tube per segment", units=shell_units["length"]
    )

    # Minimum cross-sectional area on shell side
    blk.area_flow_shell_min = Var(
        initialize=1.0, doc="Minimum flow area on shell side", units=shell_units["area"]
    )

    # total number of tube rows
    @blk.Expression(doc="Total number of tube rows")
    def nrow_tube(b):
        return b.number_passes * b.number_rows_per_pass

    # Tube outside diameter
    @blk.Expression(doc="Outside diameter of tube")
    def do_tube(b):
        return b.di_tube + b.thickness_tube * 2.0

    # Ratio of pitch_x/do_tube
    @blk.Expression(doc="Ratio of pitch in x direction to tube outside diameter")
    def pitch_x_to_do(b):
        return b.pitch_x / b.do_tube

    # Ratio of pitch_y/do_tube
    @blk.Expression(doc="Ratio of pitch in y direction to tube outside diameter")
    def pitch_y_to_do(b):
        return b.pitch_y / b.do_tube

    # Total cross-sectional area of tube metal per segment
    @blk.Expression(doc="Total cross section area of tube metal per segment")
    def area_wall_seg(b):
        return (
            0.25
            * const.pi
            * (b.do_tube**2 - b.di_tube**2)
            * b.number_columns_per_pass
            * b.number_rows_per_pass
        )

    @blk.Expression(doc="Total heat transfer area on outer surface of tubes")
    def total_heat_transfer_area(b):
        return (
            const.pi
            * b.do_tube
            * b.number_rows_per_pass
            * b.number_columns_per_pass
            * b.number_passes
            * b.length_tube_seg
        )

    # Length of shell side flow
    @blk.Constraint(doc="Length of shell side flow")
    def length_flow_shell_eqn(b):
        return b.length_flow_shell == b.nrow_tube * b.pitch_x

    # Average flow area on shell side
    @blk.Constraint(doc="Average cross section area of shell side flow")
    def area_flow_shell_eqn(b):
        return (
            b.length_flow_shell * b.area_flow_shell
            == b.length_tube_seg
            * b.length_flow_shell
            * b.pitch_y
            * b.number_columns_per_pass
            - b.number_columns_per_pass
            * b.nrow_tube
            * 0.25
            * const.pi
            * b.do_tube**2
            * b.length_tube_seg
        )

    # Minimum flow area on shell side
    @blk.Constraint(doc="Minimum flow area on shell side")
    def area_flow_shell_min_eqn(b):
        return (
            b.area_flow_shell_min
            == b.length_tube_seg * (b.pitch_y - b.do_tube) * b.number_columns_per_pass
        )


def make_performance_common(
    blk,
    shell,
    shell_units,
    shell_has_pressure_change: bool,
    make_reynolds: bool,
    make_nusselt: bool,
):
    """Function to create variables, constraints, and expressions regarding
    performance constraints shared between the CrossFlowHeatExchanger1D and
    Heater1D models.

    Args:
        blk : unit model on which components are being generated
        shell: shell control volume
        shell_units : derived units for property package of shell control volume
        shell_has_pressure_change: bool about whether to make pressure change components
        make_reynolds: bool about whether to create the Reynolds number
        make_nusselt: bool about whether to create Nusselt number

    Returns:
        None
    """
    # We need the Reynolds number for pressure change, even if we don't need it for heat transfer
    if shell_has_pressure_change:
        make_reynolds = True

    add_object_reference(blk, "heat_shell", shell.heat)

    if shell_has_pressure_change:
        add_object_reference(blk, "deltaP_shell", shell.deltaP)

    # Parameters
    # Wall thermal conductivity
    blk.therm_cond_wall = Param(
        initialize=1.0,
        mutable=True,
        units=shell_units["thermal_conductivity"],
        doc="Thermal conductivity of tube wall",
    )

    # Wall heat capacity
    blk.cp_wall = Param(
        initialize=502.4,
        mutable=True,
        units=shell_units["heat_capacity_mass"],
        doc="Tube wall heat capacity",
    )

    # Wall density
    blk.density_wall = Param(
        initialize=7800.0,
        mutable=True,
        units=shell_units["density_mass"],
        doc="Tube wall density",
    )

    # Heat transfer resistance due to the fouling on shell side
    blk.rfouling_shell = Param(
        units=1 / shell_units["heat_transfer_coefficient"],
        initialize=0.0001,
        mutable=True,
        doc="Fouling resistance on tube side",
    )

    # Correction factor for convective heat transfer coefficient on shell side
    blk.fcorrection_htc_shell = Var(
        initialize=1.0, doc="Correction factor for convective HTC on shell"
    )

    # Correction factor for shell side pressure drop due to friction
    if shell_has_pressure_change:
        blk.fcorrection_dp_shell = Var(
            initialize=1.0, doc="Correction factor for shell side pressure drop"
        )

    # Performance variables
    # Shell side convective heat transfer coefficient due to convection only
    blk.conv_heat_transfer_coeff_shell = Var(
        blk.flowsheet().config.time,
        shell.length_domain,
        initialize=100.0,
        bounds=(0, None),
        units=shell_units["heat_transfer_coefficient"],
        doc="Shell side convective heat transfer coefficient due to convection",
    )

    # Boundary wall temperature on shell side
    blk.temp_wall_shell = Var(
        blk.flowsheet().config.time,
        shell.length_domain,
        initialize=500,
        units=shell_units["temperature"],
        doc="Boundary wall temperature on shell side",
    )

    # Central wall temperature of tube metal, used to calculate energy contained by tube metal
    blk.temp_wall_center = Var(
        blk.flowsheet().config.time,
        shell.length_domain,
        initialize=500,
        units=shell_units["temperature"],
        doc="Tube wall temperature at center",
    )

    # Tube wall heat holdup per length of shell
    if blk.config.has_holdup:
        blk.heat_holdup = Var(
            blk.flowsheet().config.time,
            shell.length_domain,
            initialize=1e6,
            units=shell_units["energy"] / shell_units["length"],
            doc="Tube wall heat holdup per length of shell",
        )

        @blk.Constraint(
            blk.flowsheet().config.time,
            shell.length_domain,
            doc="Heat holdup of tube metal",
        )
        def heat_holdup_eqn(b, t, x):
            return (
                b.heat_holdup[t, x]
                == b.cp_wall
                * b.density_wall
                * b.area_wall_seg
                * pyunits.convert(b.length_flow_tube, to_units=shell_units["length"])
                / b.length_flow_shell
                * b.temp_wall_center[t, x]
            )

    # Tube wall heat accumulation term
    if blk.config.dynamic:
        blk.heat_accumulation = DerivativeVar(
            blk.heat_holdup,
            initialize=0,
            wrt=blk.flowsheet().config.time,
            units=shell_units["energy"] / shell_units["length"] / shell_units["time"],
            doc="Tube wall heat accumulation per unit length of shell",
        )

    # Pressure drop and heat transfer coefficient on shell side
    # ----------------------------------------------------------
    # Tube arrangement factor
    if blk.config.tube_arrangement == "in-line":
        blk.f_arrangement = Param(
            initialize=0.788, doc="in-line tube arrangement factor"
        )
    elif blk.config.tube_arrangement == "staggered":
        blk.f_arrangement = Param(
            initialize=1.0, doc="staggered tube arrangement factor"
        )
    else:
        raise ConfigurationError()

    if make_reynolds:
        # Velocity on shell side
        blk.v_shell = Var(
            blk.flowsheet().config.time,
            shell.length_domain,
            initialize=1.0,
            units=shell_units["velocity"],
            doc="velocity on shell side",
        )

        # Reynalds number on shell side
        blk.N_Re_shell = Var(
            blk.flowsheet().config.time,
            shell.length_domain,
            bounds=(1e-7, None),
            initialize=10000.0,
            units=pyunits.dimensionless,
            doc="Reynolds number on shell side",
        )

    if shell_has_pressure_change:
        # Friction factor on shell side
        blk.friction_factor_shell = Var(
            blk.flowsheet().config.time,
            shell.length_domain,
            initialize=1.0,
            doc="friction factor on shell side",
        )

    if make_nusselt:
        # Nusselt number on shell side
        blk.N_Nu_shell = Var(
            blk.flowsheet().config.time,
            shell.length_domain,
            initialize=1,
            units=pyunits.dimensionless,
            doc="Nusselts number on shell side",
            bounds=(1e-7, None),
        )

    if make_reynolds:
        # Velocity equation on shell side
        @blk.Constraint(
            blk.flowsheet().config.time,
            shell.length_domain,
            doc="velocity on shell side",
        )
        def v_shell_eqn(b, t, x):
            return (
                b.v_shell[t, x]
                * shell.properties[t, x].dens_mol_phase["Vap"]
                * b.area_flow_shell_min
                == shell.properties[t, x].flow_mol
            )

        # Reynolds number
        @blk.Constraint(
            blk.flowsheet().config.time,
            shell.length_domain,
            doc="Reynolds number equation on shell side",
        )
        def N_Re_shell_eqn(b, t, x):
            return (
                b.N_Re_shell[t, x] * shell.properties[t, x].visc_d_phase["Vap"]
                == b.do_tube
                * b.v_shell[t, x]
                * shell.properties[t, x].dens_mass_phase["Vap"]
            )

    if shell_has_pressure_change:
        # Friction factor on shell side
        if blk.config.tube_arrangement == "in-line":

            @blk.Constraint(
                blk.flowsheet().config.time,
                shell.length_domain,
                doc="in-line friction factor on shell side",
            )
            def friction_factor_shell_eqn(b, t, x):
                return (
                    b.friction_factor_shell[t, x] * b.N_Re_shell[t, x] ** 0.15
                    == (
                        0.044
                        + 0.08
                        * b.pitch_x_to_do
                        / (b.pitch_y_to_do - 1.0) ** (0.43 + 1.13 / b.pitch_x_to_do)
                    )
                    * b.fcorrection_dp_shell
                )

        elif blk.config.tube_arrangement == "staggered":

            @blk.Constraint(
                blk.flowsheet().config.time,
                shell.length_domain,
                doc="staggered friction factor on shell side",
            )
            def friction_factor_shell_eqn(b, t, x):
                return (
                    b.friction_factor_shell[t, x] * b.N_Re_shell[t, x] ** 0.16
                    == (0.25 + 0.118 / (b.pitch_y_to_do - 1.0) ** 1.08)
                    * b.fcorrection_dp_shell
                )

        else:
            raise ConfigurationError()

        # Pressure drop on shell side
        @blk.Constraint(
            blk.flowsheet().config.time,
            shell.length_domain,
            doc="pressure change on shell side",
        )
        def deltaP_shell_eqn(b, t, x):
            return (
                b.deltaP_shell[t, x] * b.pitch_x
                == -1.4
                * b.friction_factor_shell[t, x]
                * shell.properties[t, x].dens_mass_phase["Vap"]
                * b.v_shell[t, x] ** 2
            )

    if make_nusselt:
        # The actual Nusselt number correlation needs to be made by the particular heat exchanger
        @blk.Constraint(
            blk.flowsheet().config.time,
            shell.length_domain,
            doc="Convective heat transfer coefficient equation on shell side due to convection",
        )
        def conv_heat_transfer_coeff_shell_eqn(b, t, x):
            return (
                b.conv_heat_transfer_coeff_shell[t, x] * b.do_tube
                == b.N_Nu_shell[t, x]
                * shell.properties[t, x].therm_cond_phase["Vap"]
                * b.fcorrection_htc_shell
            )

    # Total convective heat transfer coefficient on shell side
    @blk.Expression(
        blk.flowsheet().config.time,
        shell.length_domain,
        doc="Total convective heat transfer coefficient on shell side",
    )
    def total_heat_transfer_coeff_shell(b, t, x):
        # Retain in case we add back radiation
        return b.conv_heat_transfer_coeff_shell[t, x]


def scale_common(blk, shell, shell_has_pressure_change, make_reynolds, make_nusselt):
    """Function to scale variables and constraints shared between
    the CrossFlowHeatExchanger1D and Heater1D models.

    Args:
        blk : unit model on which components are being generated
        shell: shell control volume
        shell_units : derived units for property package of shell control volume
        shell_has_pressure_change: bool about whether to make pressure change components
        make_reynolds: bool about whether to create the Reynolds number
        make_nusselt: bool about whether to create Nusselt number

    Returns:
        None
    """

    def gsf(obj):
        return iscale.get_scaling_factor(obj, default=1, warning=True)

    def ssf(obj, sf):
        iscale.set_scaling_factor(obj, sf, overwrite=False)

    def cst(con, sf):
        iscale.constraint_scaling_transform(con, sf, overwrite=False)

    sgsf = iscale.set_and_get_scaling_factor

    sf_do_tube = iscale.get_scaling_factor(blk.do_tube, default=1 / value(blk.do_tube))
    calculate_variable_from_constraint(
        blk.area_flow_shell_min, blk.area_flow_shell_min_eqn
    )
    for t in blk.flowsheet().time:
        for z in shell.length_domain:
            sf_flow_mol_shell = gsf(shell.properties[t, z].flow_mol)

            if make_reynolds:
                # FIXME get better scaling later
                ssf(blk.v_shell[t, z], 1 / 10)
                cst(blk.v_shell_eqn[t, z], sf_flow_mol_shell)

                # FIXME should get scaling of N_Re from defining eqn
                sf_N_Re_shell = sgsf(blk.N_Re_shell[t, z], 1e-4)

                sf_visc_d_shell = gsf(shell.properties[t, z].visc_d_phase["Vap"])
                cst(blk.N_Re_shell_eqn[t, z], sf_N_Re_shell * sf_visc_d_shell)
            if make_nusselt:
                sf_k_shell = gsf(shell.properties[t, z].therm_cond_phase["Vap"])

                sf_N_Nu_shell = sgsf(
                    blk.N_Nu_shell[t, z], 1 / 0.33 * sf_N_Re_shell**0.6
                )
                cst(blk.N_Nu_shell_eqn[t, z], sf_N_Nu_shell)

                sf_conv_heat_transfer_coeff_shell = sgsf(
                    blk.conv_heat_transfer_coeff_shell[t, z],
                    sf_N_Nu_shell * sf_k_shell / sf_do_tube,
                )
                cst(
                    blk.conv_heat_transfer_coeff_shell_eqn[t, z],
                    sf_conv_heat_transfer_coeff_shell * sf_do_tube,
                )

            # FIXME estimate from parameters
            if blk.config.has_holdup:
                s_U_holdup = gsf(blk.heat_holdup[t, z])
                cst(blk.heat_holdup_eqn[t, z], s_U_holdup)
