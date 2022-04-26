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
IDAES plate heat exchanger model (PHE) using effectiveness-NTU method, derived
from the core HeatExchangerNTU model.

Notes:
Let P = number of passes(same for both fluids)
if P is even PHE is in counter-current mode
if P is odd PHE  is in co-current mode
Hot and cold fluids flow alternatively in the channels.
Divider plates may be used to partition the passes in separate sections.
The heat transfer area  of a single plate depends on the total heat transfer
area as specified by the manufacturer and the total number of active plates.

Reference
Detailed model equations can be found in the paper :
[1] Akula, P., Eslick, J., Bhattacharyya, D. and Miller, D.C., 2019.
Modelling and Parameter Estimation of a Plate Heat Exchanger
as Part of a Solvent-Based Post-Combustion CO2 Capture System.
In Computer Aided Chemical Engineering (Vol. 47, pp. 47-52). Elsevier.

"""

# Import Pyomo libraries
from pyomo.environ import (
    Constraint,
    exp,
    Expression,
    NonNegativeIntegers,
    Param,
    PositiveIntegers,
    PositiveReals,
    units as pyunits,
    Var,
)
from pyomo.common.config import ConfigValue, In, Integer
from pyomo.util.calc_var_value import calculate_variable_from_constraint

# Import IDAES cores
from idaes.core import declare_process_block_class
from idaes.models.unit_models.heat_exchanger_ntu import HeatExchangerNTUData
from idaes.core.util.constants import Constants
from idaes.core.solvers import get_solver
import idaes.logger as idaeslog

__author__ = "Paul Akula, Andrew Lee"


# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("PlateHeatExchanger")
class PlateHeatExchangerData(HeatExchangerNTUData):
    """Plate Heat Exchanger(PHE) Unit Model."""

    CONFIG = HeatExchangerNTUData.CONFIG()

    CONFIG.declare(
        "passes",
        ConfigValue(
            default=4,
            domain=Integer,
            description="Number of passes",
            doc="""Number of passes of the fluids through the heat exchanger""",
        ),
    )
    CONFIG.declare(
        "channels_per_pass",
        ConfigValue(
            default=12,
            domain=Integer,
            description="Number of channels for each pass",
            doc="""Number of channels to be used in each pass where a channel
               is the space between two plates with a flowing fluid""",
        ),
    )
    CONFIG.declare(
        "number_of_divider_plates",
        ConfigValue(
            default=0,
            domain=Integer,
            description="Number of divider plates in heat exchanger",
            doc="""Divider plates are used to create separate partitions in the
        unit. Each pass can be separated by a divider plate""",
        ),
    )

    # Update config block setting for pressure change to always be true
    CONFIG.hot_side.has_pressure_change = True
    CONFIG.hot_side.get("has_pressure_change").set_domain(In([True]))
    CONFIG.hot_side.get(
        "has_pressure_change"
    )._description = "Pressure change term construction flag - must be True"
    CONFIG.hot_side.get("has_pressure_change")._doc = (
        "Plate Heat Exchanger model includes correlations for pressure drop "
        "thus has_pressure_change must be True"
    )

    CONFIG.cold_side.has_pressure_change = True
    CONFIG.cold_side.get("has_pressure_change").set_domain(In([True]))
    CONFIG.cold_side.get(
        "has_pressure_change"
    )._description = "Pressure change term construction flag - must be True"
    CONFIG.cold_side.get("has_pressure_change")._doc = (
        "Plate Heat Exchanger model includes correlations for pressure drop "
        "thus has_pressure_change must be True"
    )

    def build(self):
        # Call super.build to setup model
        # This will create the control volumes, ports and basic equations
        super().build()

        # Units will be based on hot side properties
        units_meta = (
            self.config.hot_side.property_package.get_metadata().get_derived_units
        )

        # ---------------------------------------------------------------------
        # Plate design variables and parameter
        self.number_of_passes = Param(
            initialize=self.config.passes,
            units=pyunits.dimensionless,
            domain=PositiveIntegers,
            doc="Number of hot/cold fluid passes",
            mutable=True,
        )

        # Assuming number of channels is equal in all plates
        self.channels_per_pass = Param(
            initialize=self.config.channels_per_pass,
            units=pyunits.dimensionless,
            domain=PositiveIntegers,
            doc="Number of channels in each pass",
            mutable=True,
        )

        self.number_of_divider_plates = Param(
            initialize=self.config.number_of_divider_plates,
            units=pyunits.dimensionless,
            domain=NonNegativeIntegers,
            doc="Number of divider plates in heat exchanger",
            mutable=True,
        )

        self.plate_length = Var(
            initialize=1.6925,
            units=units_meta("length"),
            domain=PositiveReals,
            doc="Length of heat exchanger plate",
        )
        self.plate_width = Var(
            initialize=0.6135,
            units=units_meta("length"),
            domain=PositiveReals,
            doc="Width of heat exchanger plate",
        )
        self.plate_thickness = Var(
            initialize=0.0006,
            units=units_meta("length"),
            domain=PositiveReals,
            doc="Thickness of heat exchanger plate",
        )
        self.plate_pact_length = Var(
            initialize=0.381,
            units=units_meta("length"),
            domain=PositiveReals,
            doc="Compressed plate pact length",
        )
        self.port_diameter = Var(
            initialize=0.2045,
            units=units_meta("length"),
            domain=PositiveReals,
            doc="Port diamter",
        )

        self.plate_therm_cond = Var(
            initialize=16.2,
            units=units_meta("thermal_conductivity"),
            domain=PositiveReals,
            doc="Thermal conductivity heat exchanger plates",
        )

        # Set default value of total heat transfer area
        self.area.set_value(114.3)

        # ---------------------------------------------------------------------
        # Derived geometric quantities
        total_plates = (
            2 * self.channels_per_pass * self.number_of_passes
            + 1
            + self.number_of_divider_plates
        )
        total_active_plates = 2 * self.channels_per_pass * self.number_of_passes - (
            1 + self.number_of_divider_plates
        )

        self.plate_gap = Expression(
            expr=self.plate_pact_length / total_plates - self.plate_thickness
        )

        self.plate_area = Expression(
            expr=self.area / total_active_plates,
            doc="Heat transfer area of single plate",
        )

        self.surface_enlargement_factor = Expression(
            expr=self.plate_area / (self.plate_length * self.plate_width)
        )

        # Channel equivalent diameter
        self.channel_diameter = Expression(
            expr=2 * self.plate_gap / self.surface_enlargement_factor,
            doc="Channel equivalent diameter",
        )

        # ---------------------------------------------------------------------
        # Fluid velocities
        def rule_port_vel_hot(blk, t):
            return (
                4
                * blk.hot_side.properties_in[t].flow_vol
                / (Constants.pi * blk.port_diameter**2)
            )

        self.hot_port_velocity = Expression(
            self.flowsheet().time, rule=rule_port_vel_hot, doc="Hot side port velocity"
        )

        def rule_port_vel_cold(blk, t):
            return (
                4
                * pyunits.convert(
                    blk.cold_side.properties_in[t].flow_vol,
                    to_units=units_meta("flow_vol"),
                )
                / (Constants.pi * blk.port_diameter**2)
            )

        self.cold_port_velocity = Expression(
            self.flowsheet().time,
            rule=rule_port_vel_cold,
            doc="Cold side port velocity",
        )

        def rule_channel_vel_hot(blk, t):
            return blk.hot_side.properties_in[t].flow_vol / (
                blk.channels_per_pass * blk.plate_width * blk.plate_gap
            )

        self.hot_channel_velocity = Expression(
            self.flowsheet().time,
            rule=rule_channel_vel_hot,
            doc="Hot side channel velocity",
        )

        def rule_channel_vel_cold(blk, t):
            return pyunits.convert(
                blk.cold_side.properties_in[t].flow_vol, to_units=units_meta("flow_vol")
            ) / (blk.channels_per_pass * blk.plate_width * blk.plate_gap)

        self.cold_channel_velocity = Expression(
            self.flowsheet().time,
            rule=rule_channel_vel_cold,
            doc="Cold side channel velocity",
        )

        # ---------------------------------------------------------------------
        # Reynolds & Prandtl numbers
        # Density cancels out of Reynolds number if mass flow rate is used
        def rule_Re_h(blk, t):
            return (
                blk.hot_side.properties_in[t].flow_mass
                * blk.channel_diameter
                / (
                    blk.channels_per_pass
                    * blk.plate_width
                    * blk.plate_gap
                    * blk.hot_side.properties_in[t].visc_d_phase["Liq"]
                )
            )

        self.Re_hot = Expression(
            self.flowsheet().time, rule=rule_Re_h, doc="Hot side Reynolds number"
        )

        def rule_Re_c(blk, t):
            return (
                pyunits.convert(
                    blk.cold_side.properties_in[t].flow_mass
                    / blk.cold_side.properties_in[t].visc_d_phase["Liq"],
                    to_units=units_meta("length"),
                )
                * blk.channel_diameter
                / (blk.channels_per_pass * blk.plate_width * blk.plate_gap)
            )

        self.Re_cold = Expression(
            self.flowsheet().time, rule=rule_Re_c, doc="Cold side Reynolds number"
        )

        def rule_Pr_h(blk, t):
            return (
                blk.hot_side.properties_in[t].cp_mol
                / blk.hot_side.properties_in[t].mw
                * blk.hot_side.properties_in[t].visc_d_phase["Liq"]
                / blk.hot_side.properties_in[t].therm_cond_phase["Liq"]
            )

        self.Pr_hot = Expression(
            self.flowsheet().time, rule=rule_Pr_h, doc="Hot side Prandtl number"
        )

        def rule_Pr_c(blk, t):
            return (
                blk.cold_side.properties_in[t].cp_mol
                / blk.cold_side.properties_in[t].mw
                * blk.cold_side.properties_in[t].visc_d_phase["Liq"]
                / blk.cold_side.properties_in[t].therm_cond_phase["Liq"]
            )

        self.Pr_cold = Expression(
            self.flowsheet().time, rule=rule_Pr_c, doc="Cold side Prandtl number"
        )

        # ---------------------------------------------------------------------
        # Heat transfer coefficients
        # Parameters for Nusselt number correlation
        self.Nusselt_param_a = Param(
            initialize=0.4,
            domain=PositiveReals,
            units=pyunits.dimensionless,
            mutable=True,
            doc="Nusselt parameter A",
        )
        self.Nusselt_param_b = Param(
            initialize=0.663,
            domain=PositiveReals,
            units=pyunits.dimensionless,
            mutable=True,
            doc="Nusselt parameter B",
        )
        self.Nusselt_param_c = Param(
            initialize=0.333,
            domain=PositiveReals,
            units=pyunits.dimensionless,
            mutable=True,
            doc="Nusselt parameter C",
        )

        # Film heat transfer coefficients
        def rule_hotside_transfer_coeff(blk, t):
            return (
                blk.hot_side.properties_in[t].therm_cond_phase["Liq"]
                / blk.channel_diameter
                * blk.Nusselt_param_a
                * blk.Re_hot[t] ** blk.Nusselt_param_b
                * blk.Pr_hot[t] ** blk.Nusselt_param_c
            )

        self.heat_transfer_coefficient_hot_side = Expression(
            self.flowsheet().time,
            rule=rule_hotside_transfer_coeff,
            doc="Hot side heat transfer coefficient",
        )

        def rule_coldside_transfer_coeff(blk, t):
            return (
                pyunits.convert(
                    blk.cold_side.properties_in[t].therm_cond_phase["Liq"],
                    to_units=units_meta("thermal_conductivity"),
                )
                / blk.channel_diameter
                * blk.Nusselt_param_a
                * blk.Re_cold[t] ** blk.Nusselt_param_b
                * blk.Pr_cold[t] ** blk.Nusselt_param_c
            )

        self.heat_transfer_coefficient_cold_side = Expression(
            self.flowsheet().time,
            rule=rule_coldside_transfer_coeff,
            doc="Cold side heat transfer coefficient",
        )

        # Overall heat transfer coefficient
        def rule_U(blk, t):
            return blk.heat_transfer_coefficient[t] == (
                1.0
                / (
                    1.0 / blk.heat_transfer_coefficient_hot_side[t]
                    + blk.plate_gap / blk.plate_therm_cond
                    + 1.0 / blk.heat_transfer_coefficient_cold_side[t]
                )
            )

        self.overall_heat_transfer_eq = Constraint(
            self.flowsheet().time,
            rule=rule_U,
            doc="Calculations of overall heat transfer coefficient",
        )

        # Effectiveness based on sub-heat exchangers
        # Divide NTU by number of channels per pass
        def rule_Ecf(blk, t):
            if blk.number_of_passes.value % 2 == 0:
                return blk.effectiveness[t] == (
                    1 - exp(-blk.NTU[t] / blk.channels_per_pass * (1 - blk.Cratio[t]))
                ) / (
                    1
                    - blk.Cratio[t]
                    * exp(-blk.NTU[t] / blk.channels_per_pass * (1 - blk.Cratio[t]))
                )
            elif blk.pass_num.value % 2 == 1:
                return blk.effectiveness[t] == (
                    1 - exp(-blk.NTU[t] / blk.channels_per_pass * (1 + blk.Cratio[t]))
                ) / (1 + blk.Cratio[t])

        self.effectiveness_correlation = Constraint(
            self.flowsheet().time,
            rule=rule_Ecf,
            doc="Correlation for effectiveness factor",
        )

        # ---------------------------------------------------------------------
        # Pressure drop correlations
        # Friction factor calculation
        self.friction_factor_param_a = Param(
            initialize=0.0,
            units=pyunits.dimensionless,
            doc="Friction factor parameter A",
            mutable=True,
        )
        self.friction_factor_param_b = Param(
            initialize=18.29,
            units=pyunits.dimensionless,
            doc="Friction factor parameter B",
            mutable=True,
        )
        self.friction_factor_param_c = Param(
            initialize=-0.652,
            units=pyunits.dimensionless,
            doc="Friction factor parameter C",
            mutable=True,
        )

        def rule_fric_h(blk, t):
            return (
                blk.friction_factor_param_a
                + blk.friction_factor_param_b
                * blk.Re_hot[t] ** (blk.friction_factor_param_c)
            )

        self.friction_factor_hot = Expression(
            self.flowsheet().time, rule=rule_fric_h, doc="Hot side friction factor"
        )

        def rule_fric_c(blk, t):
            return (
                blk.friction_factor_param_a
                + blk.friction_factor_param_b
                * blk.Re_cold[t] ** (blk.friction_factor_param_c)
            )

        self.friction_factor_cold = Expression(
            self.flowsheet().time, rule=rule_fric_c, doc="Cold side friction factor"
        )

        def rule_hotside_dP(blk, t):
            return blk.hot_side.deltaP[t] == -(
                (
                    2
                    * blk.friction_factor_hot[t]
                    * (blk.plate_length + blk.port_diameter)
                    * blk.number_of_passes
                    * blk.hot_channel_velocity[t] ** 2
                    * blk.hot_side.properties_in[t].dens_mass
                    / blk.channel_diameter
                )
                + (
                    0.7
                    * blk.number_of_passes
                    * blk.hot_port_velocity[t] ** 2
                    * blk.hot_side.properties_in[t].dens_mass
                )
                + (
                    blk.hot_side.properties_in[t].dens_mass
                    * pyunits.convert(
                        Constants.acceleration_gravity,
                        to_units=units_meta("acceleration"),
                    )
                    * (blk.plate_length + blk.port_diameter)
                )
            )

        self.hot_side_deltaP_eq = Constraint(
            self.flowsheet().time, rule=rule_hotside_dP
        )

        def rule_coldside_dP(blk, t):
            return blk.cold_side.deltaP[t] == -(
                (
                    2
                    * blk.friction_factor_cold[t]
                    * (blk.plate_length + blk.port_diameter)
                    * blk.number_of_passes
                    * blk.cold_channel_velocity[t] ** 2
                    * pyunits.convert(
                        blk.cold_side.properties_in[t].dens_mass,
                        to_units=units_meta("density_mass"),
                    )
                    / blk.channel_diameter
                )
                + (
                    0.7
                    * blk.number_of_passes
                    * blk.cold_port_velocity[t] ** 2
                    * pyunits.convert(
                        blk.cold_side.properties_in[t].dens_mass,
                        to_units=units_meta("density_mass"),
                    )
                )
                + (
                    pyunits.convert(
                        blk.cold_side.properties_in[t].dens_mass,
                        to_units=units_meta("density_mass"),
                    )
                    * pyunits.convert(
                        Constants.acceleration_gravity,
                        to_units=units_meta("acceleration"),
                    )
                    * (blk.plate_length + blk.port_diameter)
                )
            )

        self.cold_side_deltaP_eq = Constraint(
            self.flowsheet().time, rule=rule_coldside_dP
        )

    def initialize(
        self,
        hot_side_state_args=None,
        cold_side_state_args=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
        duty=None,
    ):
        """
        Heat exchanger initialization method.

        Args:
            hot_side_state_args : a dict of arguments to be passed to the
                property initialization for the hot side (see documentation of
                the specific property package) (default = None).
            cold_side_state_args : a dict of arguments to be passed to the
                property initialization for the cold side (see documentation of
                the specific property package) (default = None).
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None, use
                     default solver options)
            solver : str indicating which solver to use during
                     initialization (default = None, use default solver)
            duty : an initial guess for the amount of heat transfered. This
                should be a tuple in the form (value, units), (default
                = (1000 J/s))

        Returns:
            None

        """
        # Set solver options
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        hot_side = self.hot_side
        cold_side = self.cold_side

        # Create solver
        opt = get_solver(solver, optarg)

        flags1 = hot_side.initialize(
            outlvl=outlvl, optarg=optarg, solver=solver, state_args=hot_side_state_args
        )

        init_log.info_high("Initialization Step 1a (hot side) Complete.")

        flags2 = cold_side.initialize(
            outlvl=outlvl, optarg=optarg, solver=solver, state_args=cold_side_state_args
        )

        init_log.info_high("Initialization Step 1b (cold side) Complete.")

        # ---------------------------------------------------------------------
        # Solve unit without heat transfer equation
        self.energy_balance_constraint.deactivate()
        self.effectiveness_correlation.deactivate()
        self.effectiveness.fix(0.68)

        # Get side 1 and side 2 heat units, and convert duty as needed
        s1_units = hot_side.heat.get_units()
        s2_units = cold_side.heat.get_units()

        if duty is None:
            # Assume 1000 J/s and check for unitless properties
            if s1_units is None and s2_units is None:
                # Backwards compatability for unitless properties
                s1_duty = -1000
                s2_duty = 1000
            else:
                s1_duty = pyunits.convert_value(
                    -1000, from_units=pyunits.W, to_units=s1_units
                )
                s2_duty = pyunits.convert_value(
                    1000, from_units=pyunits.W, to_units=s2_units
                )
        else:
            # Duty provided with explicit units
            s1_duty = -pyunits.convert_value(
                duty[0], from_units=duty[1], to_units=s1_units
            )
            s2_duty = pyunits.convert_value(
                duty[0], from_units=duty[1], to_units=s2_units
            )

        cold_side.heat.fix(s2_duty)
        for i in hot_side.heat:
            hot_side.heat[i].value = s1_duty

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)

        init_log.info_high("Initialization Step 2 {}.".format(idaeslog.condition(res)))

        cold_side.heat.unfix()
        self.energy_balance_constraint.activate()

        for t in self.effectiveness:
            calculate_variable_from_constraint(
                self.effectiveness[t], self.effectiveness_correlation[t]
            )

        # ---------------------------------------------------------------------
        # Solve unit with new effectiveness factor
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        init_log.info_high("Initialization Step 3 {}.".format(idaeslog.condition(res)))

        self.effectiveness_correlation.activate()
        self.effectiveness.unfix()

        # ---------------------------------------------------------------------
        # Final solve of full modelr
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        init_log.info_high("Initialization Step 4 {}.".format(idaeslog.condition(res)))

        # ---------------------------------------------------------------------
        # Release Inlet state
        hot_side.release_state(flags1, outlvl=outlvl)
        cold_side.release_state(flags2, outlvl=outlvl)

        init_log.info("Initialization Completed, {}".format(idaeslog.condition(res)))
