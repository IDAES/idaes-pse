################################################################################
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
Packed Solvent Column Model for MEA systems
"""
# TODO: Missing docstrings
# pylint: disable=missing-function-docstring

# TODO: look into protected access issues
# pylint: disable=protected-access

# Import Python libraries
import numpy as np

# Import Pyomo libraries
from pyomo.environ import (
    value,
    Var,
    Param,
    Constraint,
    Expression,
    check_optimal_termination,
    assert_optimal_termination,
    exp,
    log,
    units as pyunits,
    Set,
)
from pyomo.common.collections import ComponentMap

from pyomo.util.calc_var_value import calculate_variable_from_constraint
from pyomo.common.config import ConfigValue

# Import IDAES Libraries
from idaes.core.util.constants import Constants as CONST
from idaes.models_extra.column_models.solvent_column import PackedColumnData

from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core import declare_process_block_class
from idaes.core.util.exceptions import InitializationError
import idaes.core.util.scaling as iscale
from idaes.core.solvers.get_solver import get_solver
import idaes.logger as idaeslog

from idaes.models_extra.column_models.enhancement_factor_model_pseudo_second_order_explicit import PseudoSecondOrderExplicit
__author__ = "Paul Akula, John Eslick, Anuja Deshpande, Andrew Lee, Douglas Allan"


@declare_process_block_class("MEAColumn")
class MEAColumnData(PackedColumnData):
    """
    MEA Packed Column Model Class.
    """

    CONFIG = PackedColumnData.CONFIG()

    CONFIG.declare(
        "enhancement_factor_model",
        ConfigValue(
            default=PseudoSecondOrderExplicit,
            description="Enhancement factor model to use in column",
        ),
    )
    CONFIG.declare(
        "enhancement_factor_kwargs",
        ConfigValue(
            default=None,
            description="Keyword arguments to pass when making enhancement factor model",
        ),
    )

    def liquid_phase_mass_transfer_model(self):
        """
        Enhancement factor based liquid phase mass transfer model.
        """
        time = self.flowsheet().time
        vap_comp = self.config.vapor_phase.property_package.component_list
        liq_comp = self.config.liquid_phase.property_package.component_list
        equilibrium_comp = vap_comp & liq_comp
        solute_comp_list = ["CO2"]
        log_diffus_liq_comp_list = self.log_diffus_liq_comp_list = [
            "CO2",
            "MEA",
        ]  # Can add ions if we want them

        lunits = (
            self.config.liquid_phase.property_package.get_metadata().get_derived_units
        )

        # Liquid phase equilibrium pressure via Enhancement factor
        self.mass_transfer_coeff_liq = Var(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            solute_comp_list,
            bounds=(0, None),
            units=lunits("length") / lunits("time"),
            doc="Liquid phase mass transfer coefficient",
        )

        self.log_enhancement_factor = Var(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            units=pyunits.dimensionless,
            bounds=(None, 100), #  100 is about where we start getting AMPL evaluation errors due to overflow
            initialize=5,
            doc="Natural logarithm of the enhancement factor",
        )

        @self.Expression(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            doc="Enhancement factor",
        )
        def enhancement_factor(blk, t, x):
            if x == blk.liquid_phase.length_domain.last():
                return Expression.Skip
            else:
                return exp(blk.log_enhancement_factor[t, x])

        @self.Expression(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            doc="Intermediate for calculating CO2 equilibrium pressure",
        )
        def psi(blk, t, x):
            if x == blk.vapor_phase.length_domain.first():
                return Expression.Skip
            else:
                zb = self.liquid_phase.length_domain.prev(x)
                return (
                    blk.enhancement_factor[t, zb]
                    * blk.mass_transfer_coeff_liq[t, zb, "CO2"]
                    / blk.mass_transfer_coeff_vap[t, x, "CO2"]
                )

        @self.Constraint(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            equilibrium_comp,
            doc="Equilibrium partial pressure at interface",
        )
        def pressure_at_interface(blk, t, x, j):
            if x == blk.vapor_phase.length_domain.first():
                return Constraint.Skip
            else:
                zb = self.liquid_phase.length_domain.prev(x)
                lprops = blk.liquid_phase.properties[t, zb]
                Pressure = pyunits.convert(
                    blk.vapor_phase.properties[t, x].pressure,
                    to_units=lunits("pressure"),
                )
                if j == "CO2":
                    return blk.pressure_equil[t, x, j] == (
                        (
                            blk.vapor_phase.properties[t, x].mole_frac_comp[j]
                            * Pressure
                            + blk.psi[t, x]
                            * lprops.conc_mol_phase_comp_true["Liq", j]
                        )
                        / (
                            1
                            + blk.psi[t, x]
                            / blk.liquid_phase.properties[t, zb].henry["Liq", j]
                        )
                    )
                else:
                    return blk.pressure_equil[t, x, j] == (
                        lprops.mole_frac_phase_comp_true["Liq", j]
                        * lprops.pressure_sat_comp[j]
                    )
    def build(self):
        super().build()
        if self.config.enhancement_factor_kwargs is None:
            enhancement_factor_kwargs = {}
        else:
            enhancement_factor_kwargs = self.config.enhancement_factor_kwargs

        # ---------------------------------------------------------------------
        # Unit level sets
        vap_comp = self.config.vapor_phase.property_package.component_list
        liq_comp = self.config.liquid_phase.property_package.component_list
        self.equilibirum_comp = equilibrium_comp = vap_comp & liq_comp
        solute_comp_list = self.solute_comp_list = Set(
            initialize=["CO2"],
            ordered=True,
            doc="Set of solutes",
        )

        lunits = (
            self.config.liquid_phase.property_package.get_metadata().get_derived_units
        )

        # Map from log variable of property to equation that defines it.
        self.log_property_var_eqn_map = ComponentMap()
        # Map from log variable of parameter (fixed var) to equation that defines it.
        self.log_parameter_var_eqn_map = ComponentMap()

        # Log variables from property package
        self.log_dens_mass_liq = Var(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            initialize=7,
            units=pyunits.dimensionless,
            doc="Logarithm of liquid mass density",
        )

        @self.Constraint(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            doc="Constraint defining log density variable",
        )
        def log_dens_mass_liq_eqn(blk, t, x):
            if x == blk.liquid_phase.length_domain.last():
                return Constraint.Skip
            return (
                exp(blk.log_dens_mass_liq[t, x]) * lunits("density_mass")
                == blk.liquid_phase.properties[t, x].dens_mass_phase["Liq"]
            )

        self.log_property_var_eqn_map[
            self.log_dens_mass_liq
        ] = self.log_dens_mass_liq_eqn

        self.log_surf_tens_liq = Var(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            initialize=-3,
            units=pyunits.dimensionless,
            doc="Logarithm of liquid surface tension",
        )

        @self.Constraint(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            doc="Constraint defining log surface tension variable",
        )
        def log_surf_tens_liq_eqn(blk, t, x):
            if x == blk.liquid_phase.length_domain.last():
                return Constraint.Skip
            return (
                exp(blk.log_surf_tens_liq[t, x]) * lunits("surface_tension")
                == blk.liquid_phase.properties[t, x].surf_tens_phase["Liq"]
            )

        self.log_property_var_eqn_map[
            self.log_surf_tens_liq
        ] = self.log_surf_tens_liq_eqn

        self.log_visc_d_liq = Var(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            bounds=(None, 100),
            initialize=1,
            units=pyunits.dimensionless,
            doc="Logarithm of liquid viscosity",
        )

        @self.Constraint(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            doc="Defines log variable for liquid viscosity",
        )
        def log_visc_d_liq_eqn(blk, t, x):
            if x == blk.liquid_phase.length_domain.last():
                return Constraint.Skip
            else:
                return (
                    exp(blk.log_visc_d_liq[t, x]) * lunits("dynamic_viscosity")
                    == blk.liquid_phase.properties[t, x].visc_d_phase["Liq"]
                )

        self.log_property_var_eqn_map[self.log_visc_d_liq] = self.log_visc_d_liq_eqn

        self.log_diffus_liq_comp = Var(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            self.log_diffus_liq_comp_list,
            bounds=(None, 100),
            initialize=1,
            units=pyunits.dimensionless,
            doc="""Logarithm of the liquid phase diffusivity of a species""",
        )

        @self.Constraint(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            self.log_diffus_liq_comp_list,
            doc="Defines log variable for liquid phase diffusivity",
        )
        def log_diffus_liq_comp_eqn(blk, t, x, j):
            if x == blk.liquid_phase.length_domain.last():
                return Constraint.Skip
            else:
                return exp(blk.log_diffus_liq_comp[t, x, j]) * lunits(
                    "diffusivity"
                ) == (blk.liquid_phase.properties[t, x].diffus_phase_comp["Liq", j])

        self.log_property_var_eqn_map[
            self.log_diffus_liq_comp
        ] = self.log_diffus_liq_comp_eqn

        self.log_dens_mass_vap = Var(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            bounds=(None, 100),
            initialize=1,
            units=pyunits.dimensionless,
            doc="Logarithm of the vapor phase mass density",
        )

        @self.Constraint(self.flowsheet().time, self.vapor_phase.length_domain)
        def log_dens_mass_vap_eqn(blk, t, x):
            if x == blk.vapor_phase.length_domain.first():
                return Constraint.Skip
            return exp(blk.log_dens_mass_vap[t, x]) * lunits(
                "density_mass"
            ) == pyunits.convert(
                blk.vapor_phase.properties[t, x].dens_mass_phase["Vap"],
                to_units=lunits("density_mass"),
            )

        self.log_property_var_eqn_map[
            self.log_dens_mass_vap
        ] = self.log_dens_mass_vap_eqn

        self.log_visc_d_vap = Var(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            initialize=-10,
            units=pyunits.dimensionless,
            doc="Logarithm of vapor viscosity",
        )

        @self.Constraint(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            doc="Defines log variable for vapor viscosity",
        )
        def log_visc_d_vap_eqn(blk, t, x):
            if x == blk.vapor_phase.length_domain.first():
                return Constraint.Skip
            else:
                return exp(blk.log_visc_d_vap[t, x]) * lunits(
                    "dynamic_viscosity"
                ) == pyunits.convert(
                    blk.vapor_phase.properties[t, x].visc_d_phase["Vap"],
                    to_units=lunits("dynamic_viscosity"),
                )

        self.log_property_var_eqn_map[self.log_visc_d_vap] = self.log_visc_d_vap_eqn

        self.log_diffus_vap_comp = Var(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            equilibrium_comp,
            initialize=1,
            units=pyunits.dimensionless,
            doc="""Logarithm of the vapor mass transfer coeff""",
        )

        @self.Constraint(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            equilibrium_comp,
            doc="""Defines log variable for the vapor species diffusivity""",
        )
        def log_diffus_vap_comp_eqn(blk, t, x, j):
            if x == blk.vapor_phase.length_domain.first():
                return Constraint.Skip
            else:
                diffus_vap_comp = pyunits.convert(
                    blk.vapor_phase.properties[t, x].diffus_phase_comp["Vap", j],
                    to_units=lunits("diffusivity"),
                )
                return (
                    exp(blk.log_diffus_vap_comp[t, x, j]) * lunits("diffusivity")
                    == diffus_vap_comp
                )

        self.log_property_var_eqn_map[
            self.log_diffus_vap_comp
        ] = self.log_diffus_vap_comp_eqn

        self.log_pressure_vap = Var(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            initialize=11.5,
            units=pyunits.dimensionless,
        )

        @self.Constraint(self.flowsheet().time, self.vapor_phase.length_domain)
        def log_pressure_vap_eqn(blk, t, x):
            if x == blk.vapor_phase.length_domain.first():
                return Constraint.Skip
            else:
                pressure = pyunits.convert(
                    blk.vapor_phase.properties[t, x].pressure,
                    to_units=lunits("pressure"),
                )
                return exp(blk.log_pressure_vap[t, x]) * lunits("pressure") == pressure

        self.log_property_var_eqn_map[self.log_pressure_vap] = self.log_pressure_vap_eqn

        self.log_dens_mol_vap = Var(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            initialize=0,
            units=pyunits.dimensionless,
        )

        @self.Constraint(self.flowsheet().time, self.vapor_phase.length_domain)
        def log_dens_mol_vap_eqn(blk, t, x):
            if x == blk.vapor_phase.length_domain.first():
                return Constraint.Skip
            else:
                dens_mol_vap = pyunits.convert(
                    blk.vapor_phase.properties[t, x].dens_mol_phase["Vap"],
                    to_units=lunits("density_mole"),
                )
                return (
                    exp(blk.log_dens_mol_vap[t, x]) * lunits("density_mole")
                    == dens_mol_vap
                )

        self.log_property_var_eqn_map[self.log_dens_mol_vap] = self.log_dens_mol_vap_eqn

        self.log_therm_cond_vap = Var(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            initialize=-3.5,
            units=pyunits.dimensionless,
        )

        @self.Constraint(self.flowsheet().time, self.vapor_phase.length_domain)
        def log_therm_cond_vap_eqn(blk, t, x):
            if x == blk.vapor_phase.length_domain.first():
                return Constraint.Skip
            else:
                therm_cond_vap = pyunits.convert(
                    blk.vapor_phase.properties[t, x].therm_cond_phase["Vap"],
                    to_units=lunits("thermal_conductivity"),
                )
                return (
                    exp(blk.log_therm_cond_vap[t, x]) * lunits("thermal_conductivity")
                    == therm_cond_vap
                )

        self.log_property_var_eqn_map[
            self.log_therm_cond_vap
        ] = self.log_therm_cond_vap_eqn

        self.log_cp_mol_vap = Var(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            initialize=3.5,
            units=pyunits.dimensionless,
        )

        @self.Constraint(self.flowsheet().time, self.vapor_phase.length_domain)
        def log_cp_mol_vap_eqn(blk, t, x):
            if x == blk.vapor_phase.length_domain.first():
                return Constraint.Skip
            else:
                cp_mol_vap = pyunits.convert(
                    blk.vapor_phase.properties[t, x].cp_mol_phase["Vap"],
                    to_units=lunits("heat_capacity_mole"),
                )
                return (
                    exp(blk.log_cp_mol_vap[t, x]) * lunits("heat_capacity_mole")
                    == cp_mol_vap
                )

        self.log_property_var_eqn_map[self.log_cp_mol_vap] = self.log_cp_mol_vap_eqn

        # Velocity calculations
        self.velocity_vap = Var(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            bounds=(0, None),
            initialize=2,
            units=lunits("velocity"),
            doc="Vapor superficial velocity",
        )

        self.velocity_liq = Var(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            bounds=(0, None),
            initialize=0.01,
            units=pyunits.m / pyunits.s,
            doc="Liquid superficial velocity",
        )

        @self.Constraint(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            doc="Vapor superficial velocity",
        )
        def velocity_vap_eqn(blk, t, x):
            if x == blk.vapor_phase.length_domain.first():
                return Constraint.Skip
            return (
                blk.velocity_vap[t, x]
                * blk.area_column
                * blk.vapor_phase.properties[t, x].dens_mol
                == blk.vapor_phase.properties[t, x].flow_mol
            )

        @self.Constraint(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            doc="Liquid superficial velocity",
        )
        def velocity_liq_eqn(blk, t, x):
            if x == blk.liquid_phase.length_domain.last():
                return Constraint.Skip
            return (
                blk.velocity_liq[t, x]
                * blk.area_column
                * blk.liquid_phase.properties[t, x].dens_mol
                == blk.liquid_phase.properties[t, x].flow_mol
            )

        # Log-form variables and constraints for vapor and liquid velocities
        # used in reformulations of correlations
        self.log_velocity_vap = Var(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            bounds=(None, 100),
            initialize=0.69,
            doc="Logarithm of vapor velocity",
        )

        @self.Constraint(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            doc="Defines log_velocity_vap",
        )
        def log_velocity_vap_eqn(blk, t, x):
            if x == blk.vapor_phase.length_domain.first():
                return Constraint.Skip
            else:
                return (
                    exp(blk.log_velocity_vap[t, x]) * lunits("velocity")
                    == blk.velocity_vap[t, x]
                )

        self.log_velocity_liq = Var(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            bounds=(None, 100),
            initialize=-4.6,
            doc="""Logarithm of liquid velocity""",
        )

        @self.Constraint(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            doc="Defines log_velocity_liq",
        )
        def log_velocity_liq_eqn(blk, t, x):
            if x == blk.liquid_phase.length_domain.last():
                return Constraint.Skip
            else:
                return (
                    exp(blk.log_velocity_liq[t, x]) * lunits("velocity")
                    == blk.velocity_liq[t, x]
                )

        self.log_property_var_eqn_map[self.log_velocity_liq] = self.log_velocity_liq_eqn
        self.log_property_var_eqn_map[self.log_velocity_vap] = self.log_velocity_vap_eqn

        @self.Expression(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            doc="Partial pressure of CO2 Henry's law [Pa]",
        )
        def PpCO2_equil_He(blk, t, x):
            if x == blk.vapor_phase.length_domain.first():
                return 1000
            else:
                zb = self.liquid_phase.length_domain.prev(x)
                return (
                    blk.liquid_phase.properties[t, zb].conc_mol_phase_comp_true[
                        "Liq", "CO2"
                    ]
                    * blk.liquid_phase.properties[t, zb].henry["Liq", "CO2"]
                )

        # #####################################################################
        # Interfacial Area model ([m2/m3]):
        # Reference: Tsai correlation,regressed by Chinen et al. 2018

        @self.Expression()
        def wetted_perimeter(blk):
            return blk.area_column * blk.packing_specific_area / blk.eps_ref

        # Original value with weird fractional power units was 0.6486
        # Corrected value was calculated by unlumping parameters to
        # permit for unit handling
        self.area_interfacial_parA = Var(
            initialize=1.43914,
            units=pyunits.dimensionless,
            doc="Interfacial area parameter A",
        )

        self.area_interfacial_parB = Var(
            initialize=0.12,
            units=pyunits.dimensionless,
            doc="Interfacial area parameter B",
        )

        self.log_area_interfacial = Var(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            bounds=(None, 100),
            initialize=5,
            doc="Logarithm of interfacial area",
        )

        @self.Constraint(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            doc="Defines log variable for interfacial area",
        )
        def log_area_interfacial_eqn(blk, t, x):
            if x == blk.vapor_phase.length_domain.first():
                return Constraint.Skip
            else:
                return (
                    exp(blk.log_area_interfacial[t, x])
                    * lunits("length") ** 2
                    / lunits("length") ** 3
                    == blk.area_interfacial[t, x]
                )

        self.log_area_interfacial_parA = Var(
            bounds=(None, 100),
            initialize=0.3640457058674088,
            doc="Logarithm of interfacial area model param A",
        )

        @self.Constraint(doc="Defines log variable for interfacial area param A")
        def log_area_interfacial_parA_eqn(blk):
            return exp(blk.log_area_interfacial_parA) == blk.area_interfacial_parA

        self.log_parameter_var_eqn_map[
            self.log_area_interfacial_parA
        ] = self.log_area_interfacial_parA_eqn

        @self.Constraint(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            doc="Defines specific interfacial area",
        )
        def area_interfacial_eqn(blk, t, x):
            if x == blk.vapor_phase.length_domain.first():
                return Constraint.Skip
            else:
                x_liq = blk.liquid_phase.length_domain.prev(x)
                log_packing_specific_area = log(
                    blk.packing_specific_area
                    / (lunits("length") ** 2 / lunits("length") ** 3)
                )
                log_g = log(
                    pyunits.convert(
                        CONST.acceleration_gravity, to_units=lunits("acceleration")
                    )
                    / lunits("acceleration")
                )
                log_cross_sectional_area = log(blk.area_column / lunits("area"))
                log_wetted_perimeter = log(blk.wetted_perimeter / lunits("length"))
                return blk.log_area_interfacial[t, x] == (
                    log_packing_specific_area
                    + blk.log_area_interfacial_parA
                    + blk.area_interfacial_parB
                    * (
                        blk.log_dens_mass_liq[t, x_liq]
                        - blk.log_surf_tens_liq[t, x_liq]
                        + (1 / 3) * log_g
                        + (4 / 3)
                        * (
                            blk.log_velocity_liq[t, x_liq]
                            + log_cross_sectional_area
                            - log_wetted_perimeter
                        )
                    )
                )

        # #####################################################################
        # Liquid holdup model
        # Reference: Tsai correlation,regressed by Chinen et al. 2018

        # Original, dimensional value is 24.23. Corrected by subtracting the contribution of alpha parameter
        # to make this parameter dimensionless
        self.holdup_parA = Var(
            initialize=11.4474, units=pyunits.dimensionless, doc="Holdup parameter A"
        )

        self.holdup_parB = Var(
            initialize=0.6471, units=pyunits.dimensionless, doc="Holdup parameter B"
        )

        # Value taken from CCSI model. Most of the factors lumped together into this factor can be calculated
        # from information we have, but the channel size appears only in this equation for the set of correlations
        # we are using. Technically that parameter exists in the SolventColumn model, but we have nothing on which
        # to base its value.
        self.log_holdup_parAlpha = Param(
            initialize=log(3.185966), # Original units were 1 / m * (m / s**2) ** (2/3)
            units=pyunits.dimensionless,
            doc="Natural logarithm of holdup parameter alpha",
            mutable=True,
        )

        self.log_holdup_liq = Var(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            bounds=(None, 100),
            initialize=-1,
            units=pyunits.dimensionless,
            doc="Logarithm of liquid holdup",
        )

        @self.Constraint(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            doc="Defines log variable for liquid holdup",
        )
        def log_holdup_liq_eqn(blk, t, x):
            if x == blk.liquid_phase.length_domain.last():
                return Constraint.Skip
            else:
                return exp(blk.log_holdup_liq[t, x]) == blk.holdup_liq[t, x]

        self.log_holdup_parA = Var(
            bounds=(None, 100),
            initialize=3.1875915348284343,
            doc="Logarithm of liquid holdup model param A",
        )

        @self.Constraint(doc="Defines log variable for liquid holdup param A")
        def log_holdup_parA_eqn(blk):
            return exp(blk.log_holdup_parA) == blk.holdup_parA

        self.log_parameter_var_eqn_map[self.log_holdup_parA] = self.log_holdup_parA_eqn

        @self.Constraint(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            doc="Volumetric liquid holdup [-]",
        )
        def holdup_liq_eqn(blk, t, x):
            if x == blk.liquid_phase.length_domain.last():
                return Constraint.Skip
            else:
                return blk.log_holdup_liq[t, x] == (
                    blk.log_holdup_parA
                    + blk.holdup_parB
                    * (
                        blk.log_holdup_parAlpha
                        + blk.log_velocity_liq[t, x]
                        + (1 / 3)
                        * (blk.log_visc_d_liq[t, x] - blk.log_dens_mass_liq[t, x])
                    )
                )

        # #####################################################################

        # Mass transfer coefficients of diffusing components in vapor phase [mol/m2.s.Pa]
        # Mass transfer coefficients, Billet and Schultes (1999) correlation,
        # where parameters are regressed by Chinen et al. (2018).

        self.Cv_ref = Var(
            initialize=0.357,
            units=pyunits.dimensionless,
            doc="""Vapor packing specific constant in Billet and Schultes
                      volumetric mass transfer coefficient correlation""",
        )

        self.log_mass_transfer_coeff_vap = Var(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            equilibrium_comp,
            initialize=-10,
            units=pyunits.dimensionless,
            doc="Logarithm of the vapor mass transfer coeff",
        )

        @self.Constraint(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            equilibrium_comp,
            doc="Defines log variable for the vapor mass transfer coeff",
        )
        def log_mass_transfer_coeff_vap_eqn(blk, t, x, j):
            if x == blk.vapor_phase.length_domain.first():
                return Constraint.Skip
            else:
                return (
                    exp(blk.log_mass_transfer_coeff_vap[t, x, j])
                    * lunits("amount")
                    / lunits("pressure")
                    / lunits("length") ** 2
                    / lunits("time")
                    == blk.mass_transfer_coeff_vap[t, x, j]
                )

        self.log_holdup_vap = Var(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            bounds=(None, 100),
            initialize=-1,
            units=pyunits.dimensionless,
            doc="""Logarithm of vapor holdup""",
        )

        @self.Constraint(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            doc="""Defines log variable for vapor holdup""",
        )
        def log_holdup_vap_eqn(blk, t, x):
            if x == blk.vapor_phase.length_domain.first():
                return Constraint.Skip
            else:
                return exp(blk.log_holdup_vap[t, x]) == blk.holdup_vap[t, x]

        self.log_Cv_ref = Var(
            bounds=(None, 100),
            initialize=-1.030019497202498,
            doc="Logarithm of the vapor mass transfer coeff param",
            units=pyunits.dimensionless,
        )

        @self.Constraint(
            doc="Defines log variable for the vapor mass transfer coeff param"
        )
        def log_Cv_ref_eqn(blk):
            return exp(blk.log_Cv_ref) == blk.Cv_ref

        self.log_parameter_var_eqn_map[self.log_Cv_ref] = self.log_Cv_ref_eqn

        @self.Constraint(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            equilibrium_comp,
            doc="Vapor phase mass transfer coefficient",
        )
        def mass_transfer_coeff_vap_eqn(blk, t, x, j):
            if x == blk.vapor_phase.length_domain.first():
                return Constraint.Skip
            else:
                log_R_gas = log(
                    pyunits.convert(CONST.gas_constant, lunits("gas_constant"))
                    / lunits("gas_constant")
                )
                log_vapor_temp = log(
                    pyunits.convert(
                        blk.vapor_phase.properties[t, x].temperature,
                        to_units=lunits("temperature"),
                    )
                    / lunits("temperature")
                )
                log_packing_specific_area = log(
                    blk.packing_specific_area
                    / (lunits("length") ** 2 / lunits("length") ** 3)
                )
                log_hydraulic_diameter = log(blk.hydraulic_diameter / lunits("length"))
                return blk.log_mass_transfer_coeff_vap[t, x, j] == (
                    blk.log_Cv_ref
                    - log_R_gas
                    - log_vapor_temp
                    - 0.5 * blk.log_holdup_vap[t, x]
                    + 0.5 * (log_packing_specific_area - log_hydraulic_diameter)
                    + (2 / 3) * blk.log_diffus_vap_comp[t, x, j]
                    + (1 / 3) * (blk.log_visc_d_vap[t, x] - blk.log_dens_mass_vap[t, x])
                    + (3 / 4)
                    * (
                        blk.log_velocity_vap[t, x]
                        + blk.log_dens_mass_vap[t, x]
                        - log_packing_specific_area
                        - blk.log_visc_d_vap[t, x]
                    )
                )

        # #####################################################################

        # Mass transfer coefficients of diffusing components in liquid phase  [m/s]
        # Mass transfer coefficients, Billet and Schultes (1999) correlation,
        # where parameters are regressed by Chinen et al. (2018).

        self.Cl_ref = Var(
            initialize=0.5,
            units=pyunits.dimensionless,
            doc="""Liquid packing specific constant in Billet and Schultes 
                      volumetric mass transfer coefficient correlation""",
        )

        self.log_mass_transfer_coeff_liq = Var(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            solute_comp_list,
            initialize=-9,
            units=pyunits.dimensionless,
            doc="""Logarithm of the liquid mass transfer coeff""",
        )

        @self.Constraint(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            solute_comp_list,
            doc="""Defines log variable for the liquid mass transfer coeff""",
        )
        def log_mass_transfer_coeff_liq_eqn(blk, t, x, j):
            if x == blk.liquid_phase.length_domain.last():
                return Constraint.Skip
            else:
                return (
                    exp(blk.log_mass_transfer_coeff_liq[t, x, j]) * lunits("velocity")
                    == blk.mass_transfer_coeff_liq[t, x, j]
                )

        self.log_Cl_ref = Var(
            bounds=(None, 100),
            initialize=-0.30102999566,
            units=pyunits.dimensionless,
            doc="""Logarithm of the liquid mass transfer coeff param""",
        )

        @self.Constraint(
            doc="""Defines log variable for the liquid mass transfer coeff param"""
        )
        def log_Cl_ref_eqn(blk):
            return exp(blk.log_Cl_ref) == blk.Cl_ref

        self.log_parameter_var_eqn_map[self.log_Cl_ref] = self.log_Cl_ref_eqn

        @self.Constraint(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            solute_comp_list,
            doc="Liquid phase mass transfer coefficient",
        )
        def mass_transfer_coeff_liq_eqn(blk, t, x, j):
            if x == blk.liquid_phase.length_domain.last():
                return Constraint.Skip
            else:
                log_hydraulic_diameter = log(blk.hydraulic_diameter / lunits("length"))
                return blk.log_mass_transfer_coeff_liq[t, x, j] == (
                    blk.log_Cl_ref
                    + (1 / 6) * log(12)
                    + 0.5
                    * (
                        blk.log_velocity_liq[t, x]
                        + blk.log_diffus_liq_comp[t, x, j]
                        - blk.log_holdup_liq[t, x]
                        - log_hydraulic_diameter
                    )
                )

        # #####################################################################

        # TODO We shouldn't be fixing variables in a construction method
        # Fix mass transfer parameters
        self.Cv_ref.fix()
        self.Cl_ref.fix()

        # Fix interfacial area parameters
        self.area_interfacial_parA.fix()
        self.area_interfacial_parB.fix()

        # Fix liquid holdup parameters
        self.holdup_parA.fix()
        self.holdup_parB.fix()

        # ---------------------------------------------------------------------
        # Vapor-liquid heat transfer coeff modified by Ackmann factor
        self.heat_transfer_coeff_base = Var(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            initialize=100,
            units=lunits("power") / lunits("temperature") / lunits("length") ** 2,
            doc="Uncorrected vapor-liquid heat transfer coefficient",
        )

        self.log_heat_transfer_coeff_base = Var(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            initialize=5,
            units=pyunits.dimensionless,
        )

        @self.Constraint(self.flowsheet().time, self.vapor_phase.length_domain)
        def log_heat_transfer_coeff_base_eqn(blk, t, x):
            if x == blk.vapor_phase.length_domain.first():
                return Constraint.Skip
            else:
                return (
                    exp(blk.log_heat_transfer_coeff_base[t, x])
                    * lunits("heat_transfer_coefficient")
                    == blk.heat_transfer_coeff_base[t, x]
                )

        # Heat transfer coefficients, Chilton Colburn  analogy
        # Vapor-liquid heat transfer coefficient [J/m2.s.K]
        @self.Constraint(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            doc="""Define vap-liq heat transfer coefficient""",
        )
        def heat_transfer_coeff_base_eqn(blk, t, x):
            if x == blk.vapor_phase.length_domain.first():
                return Constraint.Skip
            else:
                return 3 * blk.log_heat_transfer_coeff_base[t, x] == 3 * (
                    blk.log_mass_transfer_coeff_vap[t, x, "CO2"] + blk.log_pressure_vap[t, x]
                ) + 2 * blk.log_therm_cond_vap[t, x] + blk.log_cp_mol_vap[t, x] - 2 * (
                    blk.log_dens_mol_vap[t, x] + blk.log_diffus_vap_comp[t, x, "CO2"]
                )

        @self.Constraint(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            doc="""Vap-Liq heat transfer correction""",
        )        
        def heat_transfer_coeff_eqn(blk, t, x):
            if x == blk.vapor_phase.length_domain.first():
                return Constraint.Skip
            else:
                return (
                    blk.heat_transfer_coeff[t, x]
                    == blk.heat_transfer_coeff_base[t, x]
                    * blk.area_interfacial[t, x]
                    * blk.area_column
                )


        self.enhancement_factor_vars, self.enhancement_factor_constraints = self.config.enhancement_factor_model.make_model(
            self,
            **enhancement_factor_kwargs
        )
        # Flood point calculations

        self.log_flow_mass_Liq_Vap = Var(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            bounds=(None, 100),
            initialize=1,
        )

        @self.Constraint(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            doc="Constraint for log variable",
        )
        def log_flow_mass_Liq_Vap_eqn(blk, t, x):
            if x == blk.vapor_phase.length_domain.first():
                return Constraint.Skip
            else:
                x_liq = blk.liquid_phase.length_domain.prev(x)
                return (
                    exp(blk.log_flow_mass_Liq_Vap[t, x])
                    * blk.vapor_phase.properties[t, x].flow_mass_phase["Vap"]
                    == blk.liquid_phase.properties[t, x_liq].flow_mass_phase["Liq"]
                )

        self.log_flood_H = Var(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            bounds=(None, 100),
            initialize=-1,
            units=pyunits.dimensionless,
        )

        @self.Constraint(self.flowsheet().time, self.vapor_phase.length_domain)
        def flood_H_eqn(blk, t, x):
            if x == blk.vapor_phase.length_domain.first():
                return Constraint.Skip
            else:
                x_liq = blk.liquid_phase.length_domain.prev(x)
                return blk.log_flood_H[t, x] == blk.log_flow_mass_Liq_Vap[
                    t, x
                ] + 0.5 * (
                    blk.log_dens_mass_vap[t, x] - blk.log_dens_mass_liq[t, x_liq]
                )

        self.fourth_root_flood_H = Var(
            self.flowsheet().time, self.vapor_phase.length_domain, initialize=1
        )

        @self.Constraint(self.flowsheet().time, self.vapor_phase.length_domain)
        def fourth_root_flood_H_eqn(blk, t, x):
            if x == blk.vapor_phase.length_domain.first():
                return Constraint.Skip
            return blk.fourth_root_flood_H[t, x] == exp(0.25 * blk.log_flood_H[t, x])

        self.gas_velocity_fld = Var(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            initialize=1,
            units=lunits("velocity"),
            doc="Flooding velocity",
        )

        self.log_gas_velocity_fld = Var(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            bounds=(None, 100),
            initialize=0.01,
        )

        @self.Constraint(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            doc="Constraint for log variable",
        )
        def log_gas_velocity_fld_eqn(blk, t, x):
            if x == blk.vapor_phase.length_domain.first():
                return Constraint.Skip
            else:
                return (
                    exp(blk.log_gas_velocity_fld[t, x]) * lunits("velocity")
                    == blk.gas_velocity_fld[t, x]
                )

        @self.Constraint(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
        )
        def gas_velocity_fld_eqn(blk, t, x):
            if x == blk.vapor_phase.length_domain.first():
                return Expression.Skip
            else:
                x_liq = blk.liquid_phase.length_domain.prev(x)
                log_packing_specific_area = log(
                    blk.packing_specific_area
                    / (lunits("length") ** 2 / lunits("length") ** 3)
                )
                log_g = log(
                    pyunits.convert(
                        CONST.acceleration_gravity, to_units=lunits("acceleration")
                    )
                    / lunits("acceleration")
                )
                # Value of 0.001 Pa*s for viscosity of water was here when I found it: Doug
                log_visc_d_liq_H2O = log(
                    pyunits.convert(
                        0.001 * pyunits.Pa * pyunits.s,
                        to_units=lunits("dynamic_viscosity"),
                    )
                    / lunits("dynamic_viscosity")
                )
                return blk.log_gas_velocity_fld[t, x] == 0.5 * (
                    log_g
                    + 3 * log(blk.eps_ref)
                    - log_packing_specific_area
                    + blk.log_dens_mass_liq[t, x_liq]
                    - blk.log_dens_mass_vap[t, x]
                    - 0.2 * (blk.log_visc_d_liq[t, x_liq] - log_visc_d_liq_H2O)
                    - 4 * blk.fourth_root_flood_H[t, x]
                )

        self.flood_fraction = Var(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            initialize=0.7,
            units=pyunits.dimensionless,
            doc="Dimensionless flooding fraction",
        )

        @self.Constraint(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            doc="Flooding fraction (expected to be less than 0.8)",
        )
        def flood_fraction_eqn(blk, t, x):
            if x == blk.vapor_phase.length_domain.first():
                return Constraint.Skip

            else:
                return (
                    blk.flood_fraction[t, x] * blk.gas_velocity_fld[t, x]
                    == blk.velocity_vap[t, x]
                )

    # Scaling Routine
    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        vunits = (
            self.config.vapor_phase.property_package.get_metadata().get_derived_units
        )
        lunits = (
            self.config.liquid_phase.property_package.get_metadata().get_derived_units
        )

        def gsf(var):
            return iscale.get_scaling_factor(var, default=1, warning=True)

        def ssf(var, s):
            iscale.set_scaling_factor(var, s, overwrite=False)

        def cst(con, s):
            iscale.constraint_scaling_transform(con, s, overwrite=False)

        for t in self.flowsheet().time:
            for x_liq in self.liquid_phase.length_domain:
                if x_liq == self.liquid_phase.length_domain.last():
                    continue
                x_vap = self.vapor_phase.length_domain.next(x_liq)
                sf_flow_mol = gsf(self.liquid_phase.properties[t, x_liq].flow_mol)

                cst(self.velocity_liq_eqn[t, x_liq], sf_flow_mol)
                # Decent initial guess for liquid velocity
                sf_v = iscale.get_scaling_factor(
                    self.velocity_liq[t, x_liq], default=0.016, warning=False
                )
                ssf(self.log_velocity_liq[t, x_liq], 1)
                cst(self.log_velocity_liq_eqn[t, x_liq], sf_v)

                sf_A_interfacial = iscale.get_scaling_factor(
                    self.area_interfacial[t, x_vap],
                    default=1 / value(self.packing_specific_area),
                    warning=False,
                )
                # In the (common) case where it didn't have anything set, set the default value
                ssf(self.area_interfacial[t, x_vap], sf_A_interfacial)
                cst(self.log_area_interfacial_eqn[t, x_vap], sf_A_interfacial)
                ssf(self.log_area_interfacial[t, x_vap], 1)

                sf_rho_l = iscale.get_scaling_factor(
                    self.liquid_phase.properties[t, x_liq].dens_mass_phase["Liq"],
                    default=1e-3,
                    warning=False,
                )
                cst(self.log_dens_mass_liq_eqn[t, x_liq], sf_rho_l)

                sf_rho_v = iscale.get_scaling_factor(
                    self.vapor_phase.properties[t, x_vap].dens_mass_phase["Vap"],
                    default=1,
                    warning=False,
                )
                cst(self.log_dens_mass_vap_eqn[t, x_vap], sf_rho_v)

                sf = gsf(self.vapor_phase.properties[t, x_vap].pressure)
                cst(self.log_pressure_vap_eqn[t, x_vap], sf)

                sf = iscale.get_scaling_factor(
                    self.vapor_phase.properties[t, x_vap].dens_mol_phase["Vap"],
                    default=1 / 50,
                    warning=False,
                )
                cst(self.log_dens_mol_vap_eqn[t, x_vap], sf)

                sf = iscale.get_scaling_factor(
                    self.vapor_phase.properties[t, x_vap].therm_cond_phase["Vap"],
                    default=100,
                    warning=False,
                )
                cst(self.log_therm_cond_vap_eqn[t, x_vap], sf)

                sf = iscale.get_scaling_factor(
                    self.vapor_phase.properties[t, x_vap].cp_mol_phase["Vap"],
                    default=1 / 50,
                    warning=False,
                )
                cst(self.log_cp_mol_vap_eqn[t, x_vap], sf)

                sf = iscale.get_scaling_factor(
                    self.heat_transfer_coeff_base[t, x_vap],
                    default=1 / 200,
                    warning=False,
                )
                ssf(self.heat_transfer_coeff_base[t, x_vap], sf)
                cst(self.log_heat_transfer_coeff_base_eqn[t, x_vap], sf)

                for j in self.log_diffus_liq_comp_list:
                    sf_diffus_liq_comp = iscale.get_scaling_factor(
                        self.liquid_phase.properties[t, x_liq].diffus_phase_comp[
                            "Liq", j
                        ],
                        default=2e8,
                        warning=False,
                    )
                    cst(self.log_diffus_liq_comp_eqn[t, x_liq, j], sf_diffus_liq_comp)

                for j in self.solute_comp_list:
                    sf = iscale.get_scaling_factor(
                        self.mass_transfer_coeff_liq[t, x_liq, j],
                        default=1e4,
                        warning=False,
                    )
                    ssf(self.mass_transfer_coeff_liq[t, x_liq, j], sf)
                    cst(self.log_mass_transfer_coeff_liq_eqn[t, x_liq, j], sf)

                sf = gsf(self.liquid_phase.properties[t, x_liq].visc_d_phase["Liq"])
                cst(self.log_visc_d_liq_eqn[t, x_liq], sf)

                sf = gsf(self.vapor_phase.properties[t, x_vap].visc_d_phase["Vap"])
                cst(self.log_visc_d_vap_eqn[t, x_vap], sf)

                sf = gsf(self.vapor_phase.properties[t, x_vap].flow_mol)
                cst(self.velocity_vap_eqn[t, x_vap], sf)

                sf = gsf(self.liquid_phase.properties[t, x_liq].flow_mass_phase["Liq"])
                cst(self.log_flow_mass_Liq_Vap_eqn[t, x_vap], sf)

                for j in self.equilibirum_comp:
                    sf = iscale.get_scaling_factor(
                        self.mass_transfer_coeff_vap[t, x_vap, j],
                        default=25000,
                        warning=False,
                    )
                    ssf(self.mass_transfer_coeff_vap[t, x_vap, j], sf)
                    cst(self.log_mass_transfer_coeff_vap_eqn[t, x_vap, j], sf)

                    sf = iscale.get_scaling_factor(
                        self.vapor_phase.properties[t, x_vap].diffus_phase_comp[
                            "Vap", j
                        ],
                        default=5e4,
                        warning=False,
                    )
                    cst(self.log_diffus_vap_comp_eqn[t, x_vap, j], sf)

        # TODO bring this into new form later
        for (t, x), con in self.heat_transfer_coeff_eqn.items():
            # TODO Figure out how to calculate this
            ssf(self.heat_transfer_coeff[t, x], 3e-6)
            # Formulated to have the same scale as the heat transfer coefficient
            cst(con, gsf(self.heat_transfer_coeff[t, x]))

        for (t, x), v in self.heat_transfer_eqn1.items():
            iscale.constraint_scaling_transform(
                v,
                iscale.get_scaling_factor(
                    self.vapor_phase.heat[t, x], default=1, warning=True
                ),
            )

        for (t, x), v in self.heat_transfer_eqn2.items():
            iscale.constraint_scaling_transform(
                v,
                iscale.get_scaling_factor(
                    self.liquid_phase.heat[t, x], default=1, warning=True
                ),
            )

        for (t, x), v in self.enthalpy_transfer_eqn1.items():
            iscale.constraint_scaling_transform(
                v,
                iscale.get_scaling_factor(
                    self.vapor_phase.enthalpy_transfer[t, x], default=1, warning=True
                ),
            )

        for (t, x), v in self.enthalpy_transfer_eqn2.items():
            iscale.constraint_scaling_transform(
                v,
                iscale.get_scaling_factor(
                    self.liquid_phase.enthalpy_transfer[t, x], default=1, warning=True
                ),
            )


    def set_init_values_correlation_vars(blk, nfe, mode):
        """
        This method fixes the initial values of mass and heat transfer coefficients
        interfacial area, enhancement factor, and liquid holdup, based on
        the number of finite elements required by the user.

        A known set of values is used as the basis for fixing the initial values.
        If the number of finite elements required by the user differs from the
        length of the base set of values, linear interpolation is used to
        obtain the initial values.
        """

        # Define vapor phase mass transfer coefficient known values
        k_v_co2_values = [
            0,
            2.81735e-05,
            2.82192e-05,
            2.82675e-05,
            2.83159e-05,
            2.83646e-05,
            2.84139e-05,
            2.84643e-05,
            2.85162e-05,
            2.85701e-05,
            2.86267e-05,
            2.86865e-05,
            2.87505e-05,
            2.88196e-05,
            2.88949e-05,
            2.89779e-05,
            2.90706e-05,
            2.91755e-05,
            2.92957e-05,
            2.94356e-05,
            2.96012e-05,
            2.98008e-05,
            3.00461e-05,
            3.03537e-05,
            3.07478e-05,
            3.12618e-05,
            3.19377e-05,
            3.28062e-05,
            3.37808e-05,
            3.42191e-05,
            3.1977e-05,
        ]

        k_v_h2o_values = [
            0,
            3.70671e-05,
            3.7119e-05,
            3.71728e-05,
            3.72266e-05,
            3.72807e-05,
            3.73355e-05,
            3.73914e-05,
            3.7449e-05,
            3.75089e-05,
            3.75716e-05,
            3.76381e-05,
            3.7709e-05,
            3.77855e-05,
            3.7869e-05,
            3.7961e-05,
            3.80637e-05,
            3.81798e-05,
            3.83129e-05,
            3.84679e-05,
            3.86512e-05,
            3.88721e-05,
            3.91436e-05,
            3.94843e-05,
            3.99212e-05,
            4.0492e-05,
            4.12453e-05,
            4.22199e-05,
            4.33301e-05,
            4.38732e-05,
            4.14247e-05,
        ]

        # Define liquid phase mass transfer coefficient known values
        k_l_co2_values = [
            9.40406e-05,
            9.45346e-05,
            9.50364e-05,
            9.55322e-05,
            9.60246e-05,
            9.65175e-05,
            9.70151e-05,
            9.75218e-05,
            9.80421e-05,
            9.85809e-05,
            9.91439e-05,
            9.97374e-05,
            0.000100369,
            0.000101047,
            0.000101783,
            0.000102591,
            0.000103486,
            0.000104491,
            0.000105632,
            0.000106947,
            0.000108484,
            0.000110305,
            0.000112497,
            0.000115169,
            0.000118462,
            0.000122521,
            0.000127414,
            0.000132803,
            0.000136779,
            0.000131874,
            0.001,
        ]

        # Define heat transfer coefficient values
        h_v_values = [
            100,
            101.6962933,
            101.8538124,
            102.001834,
            102.148048,
            102.2944317,
            102.4423311,
            102.5930583,
            102.7479984,
            102.9086651,
            103.0767604,
            103.2542451,
            103.4434301,
            103.6470939,
            103.8686406,
            104.1123158,
            104.3835061,
            104.6891628,
            105.0384114,
            105.4434399,
            105.9208118,
            106.4934206,
            107.1933964,
            108.066322,
            109.176855,
            110.614142,
            112.487677,
            114.87132,
            117.5153935,
            118.6672687,
            112.4694233,
        ]

        # Define interfacial area known values
        interfacial_area_values = [
            0,
            197.8185136,
            197.8671837,
            197.9188717,
            197.9708021,
            198.0230645,
            198.076072,
            198.130292,
            198.1862321,
            198.2444542,
            198.3055974,
            198.3704047,
            198.4397597,
            198.5147328,
            198.5966454,
            198.6871572,
            198.7883889,
            198.9030972,
            199.0349284,
            199.1887917,
            199.3714147,
            199.5921766,
            199.8643474,
            200.2068651,
            200.6465777,
            201.219827,
            201.9679499,
            202.9052757,
            203.8836586,
            204.1277654,
            201.2884216,
        ]

        # Define enhancement factor known values
        if mode == "absorber":
            enhancement_factor_values = [
                10.59208019,
                10.83945557,
                11.08225811,
                11.3273241,
                11.5773751,
                11.83479677,
                12.10215397,
                12.38234194,
                12.6787182,
                12.99526603,
                13.33681033,
                13.70931137,
                14.12027367,
                14.57932831,
                15.09908116,
                15.69637908,
                16.39424972,
                17.2249598,
                18.23499171,
                19.49343044,
                21.10665852,
                23.24523241,
                26.19540321,
                30.46301154,
                36.99441863,
                47.67140404,
                66.4626304,
                101.9815504,
                169.7954706,
                258.2290563,
                1,
            ]
        elif mode == "stripper":
            enhancement_factor_values = [
                23.68,
                22.39,
                21.25,
                20.24,
                19.78,
                18.93,
                18.15,
                17.45,
                17.12,
                16.80,
                16.21,
                15.93,
                15.41,
                15.16,
                14.92,
                14.47,
                14.25,
                14.05,
                13.85,
                13.66,
                13.48,
                13.29,
                13.12,
                12.8,
                12.64,
                12.5,
                12.36,
                12.25,
                12.21,
                12.38,
                1,
            ]
        else:
            raise RuntimeError

        nfe_base = len(interfacial_area_values) - 1
        x_nfe_list_base = [i / nfe_base for i in range(nfe_base + 1)]

        # Linear Interpolation
        def interpolate_init_values(numdiscretepts, xdata, ydata):

            x_interpolate = [i / numdiscretepts for i in range(numdiscretepts)]
            y_interpolate = np.interp(x_interpolate, xdata, ydata)

            return y_interpolate

        # Interpolate values based on nfe required by user
        if nfe != nfe_base:
            kvco2_interp = interpolate_init_values(
                nfe + 1, x_nfe_list_base, k_v_co2_values
            )
            k_v_co2_values = kvco2_interp

            kvh2o_interp = interpolate_init_values(
                nfe + 1, x_nfe_list_base, k_v_h2o_values
            )
            k_v_h2o_values = kvh2o_interp

            klco2_interp = interpolate_init_values(
                nfe + 1, x_nfe_list_base, k_l_co2_values
            )
            k_l_co2_values = klco2_interp

            hv_interp = interpolate_init_values(nfe + 1, x_nfe_list_base, h_v_values)
            h_v_values = hv_interp

            interfarea_interp = interpolate_init_values(
                nfe + 1, x_nfe_list_base, interfacial_area_values
            )
            interfacial_area_values = interfarea_interp

            enhfactor_interp = interpolate_init_values(
                nfe + 1, x_nfe_list_base, enhancement_factor_values
            )
            enhancement_factor_values = enhfactor_interp

        # Fix values for all correlation variables
        for i, x in enumerate(blk.vapor_phase.length_domain):
            if x == blk.vapor_phase.length_domain.first():
                blk.mass_transfer_coeff_vap[0, x, "CO2"].fix(0.001)
                blk.mass_transfer_coeff_vap[0, x, "H2O"].fix(0.001)
                blk.heat_transfer_coeff_base[0, x].fix(100)
                blk.area_interfacial[0, x].fix(0.001)
            else:
                blk.mass_transfer_coeff_vap[0, x, "CO2"].fix(k_v_co2_values[i])
                blk.mass_transfer_coeff_vap[0, x, "H2O"].fix(k_v_h2o_values[i])
                blk.heat_transfer_coeff_base[0, x].fix(h_v_values[i])
                blk.area_interfacial[0, x].fix(interfacial_area_values[i])

        for i, x in enumerate(blk.liquid_phase.length_domain):
            if x == blk.liquid_phase.length_domain.last():
                blk.mass_transfer_coeff_liq[0, x, "CO2"].fix(0.001)
                blk.holdup_liq[0, x].fix(0.001)
                blk.log_enhancement_factor[0, x].fix(0)
            else:
                blk.mass_transfer_coeff_liq[0, x, "CO2"].fix(k_l_co2_values[i])
                blk.holdup_liq[0, x].fix(0.01)
                blk.log_enhancement_factor[0, x].fix(
                    log(enhancement_factor_values[i])
                )

    # =========================================================================
    # Model initialization routine
    def initialize(
        blk,
        vapor_phase_state_args=None,
        liquid_phase_state_args=None,
        state_vars_fixed=False,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
        mode="absorber",
        monolithic_enhancement_factor_system=True,
    ):
        """
        MEA solvent Packed Column initialization.

        Arguments:
            state_args : a dict of arguments to be passed to the property
                          package(s) to provide an initial state for
                          initialization (see documentation of the specific
                          property package) (default = None).
            optarg : solver options dictionary object (default=None, use
                      default solver options)
            solver : str indicating which solver to use during initialization
                    (default = None, use IDAES default solver)
        """

        # Set up logger for initialization and solve
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")

        # Set solver options
        if optarg is None:
            optarg = {}

        opt = get_solver(solver, optarg)

        initial_dof = degrees_of_freedom(blk)

        lunits = (
            blk.config.liquid_phase.property_package.get_metadata().get_derived_units
        )

        unit_constraints = [
            "pressure_at_interface",
            "interphase_mass_transfer_eqn",
            "liquid_mass_transfer_eqn",
            "vapor_mass_transfer_eqn",
            "Ack_intermediate_eqn",
            "heat_transfer_coeff_eqn",
            "heat_transfer_eqn1",
            "heat_transfer_eqn2",
            "enthalpy_transfer_eqn1",
            "enthalpy_transfer_eqn2",
        ]

        interfacial_area_constraints = [
            "area_interfacial_eqn",
            "log_area_interfacial_eqn",
        ]

        liquid_holdup_constraints = [
            "holdup_liq_eqn",
            "log_holdup_liq_eqn",
        ]

        mass_transfer_coeff_vap_constraints = [
            "log_mass_transfer_coeff_vap_eqn",
            "log_holdup_vap_eqn",
            "log_Cv_ref_eqn",
            "mass_transfer_coeff_vap_eqn",
        ]

        mass_transfer_coeff_liq_constraints = [
            "log_mass_transfer_coeff_liq_eqn",
            "log_Cl_ref_eqn",
            "mass_transfer_coeff_liq_eqn",
        ]

        heat_transfer_coefficient_constraint = [
            "heat_transfer_coeff_base_eqn",
            "log_heat_transfer_coeff_base_eqn",
        ]

        flooding_constraints = [
            "log_flood_H_eqn",
            "log_flow_mass_Liq_Vap_eqn",
            "flood_H_eqn",
            "fourth_root_flood_H_eqn",
            "log_gas_velocity_fld_eqn",
            "gas_velocity_fld_eqn",
            "flood_fraction_eqn",
        ]

        # ---------------------------------------------------------------------
        # Deactivate unit model level constraints
        for c in blk.component_objects(Constraint, descend_into=True):
            if c.local_name in unit_constraints:
                c.deactivate()
            if c.local_name in interfacial_area_constraints:
                c.deactivate()
            if c.local_name in liquid_holdup_constraints:
                c.deactivate()
            if c.local_name in mass_transfer_coeff_vap_constraints:
                c.deactivate()
            if c.local_name in mass_transfer_coeff_liq_constraints:
                c.deactivate()
            if c.local_name in heat_transfer_coefficient_constraint:
                c.deactivate()
            if c.local_name in flooding_constraints:
                c.deactivate()

        for constraint in blk.enhancement_factor_constraints:
            constraint.deactivate()

        # Fix variables

        nfe = blk.config.finite_elements
        blk.set_init_values_correlation_vars(nfe, mode)

        # Interface pressure
        blk.pressure_equil.fix()

        # Molar flux
        blk.interphase_mass_transfer.fix(0.0)
        blk.vapor_phase.mass_transfer_term.fix(0.0)
        blk.liquid_phase.mass_transfer_term.fix(0.0)

        # Heat transfer rate
        blk.heat_transfer_coeff.fix()
        blk.vapor_phase.heat.fix(0.0)
        blk.liquid_phase.heat.fix(0.0)
        blk.vapor_phase.enthalpy_transfer.fix(0.0)
        blk.liquid_phase.enthalpy_transfer.fix(0.0)

        # Velocity
        blk.velocity_vap.fix()
        blk.velocity_liq.fix()

        blk.area_interfacial_parA.fix()
        blk.area_interfacial.fix()
        blk.log_area_interfacial.fix()

        blk.holdup_parA.fix()
        blk.holdup_liq.fix()
        blk.log_holdup_liq.fix()

        # Enhancement factor
        for var in blk.enhancement_factor_vars:
            var.fix()

        # Flood fraction
        blk.log_flood_H.fix()
        blk.log_flow_mass_Liq_Vap.fix()
        blk.fourth_root_flood_H.fix()
        blk.gas_velocity_fld.fix()
        blk.log_gas_velocity_fld.fix()
        blk.flood_fraction.fix()

        # Mass transfer coefficients
        blk.Cl_ref.fix()
        blk.log_mass_transfer_coeff_liq.fix()

        blk.log_mass_transfer_coeff_vap.fix()
        blk.log_holdup_vap.fix()
        blk.Cv_ref.fix()

        init_log.info("Step 1: Property Package initialization")

        # Initialize vapor_phase properties block
        vflag = blk.vapor_phase.initialize(
            state_args=vapor_phase_state_args,
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            hold_state=True,
        )

        # Initialize liquid_phase properties block
        lflag = blk.liquid_phase.initialize(
            state_args=liquid_phase_state_args,
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            hold_state=True,
        )

        # Calculate log variables of (supposed to be fixed) parameters
        for var, eqn in blk.log_parameter_var_eqn_map.items():
            for k in eqn:
                calculate_variable_from_constraint(var[k], eqn[k])

        # Calculate these velocities first because their log counterparts are in log_properties_var_eqn_map
        blk.velocity_liq.unfix()
        blk.velocity_liq_eqn.activate()
        for k in blk.velocity_liq_eqn:
            calculate_variable_from_constraint(
                blk.velocity_liq[k], blk.velocity_liq_eqn[k]
            )

        blk.velocity_vap.unfix()
        blk.velocity_vap_eqn.activate()
        for k in blk.velocity_vap_eqn:
            calculate_variable_from_constraint(
                blk.velocity_vap[k], blk.velocity_vap_eqn[k]
            )

        # Calculate log variables that can be directly calculated from properties from the property packages
        # (and log variables for superficial velocity, which itself can be calculated from the initialized properties)
        for var, eqn in blk.log_property_var_eqn_map.items():
            var.unfix()
            eqn.activate()
            for k in eqn:
                calculate_variable_from_constraint(var[k], eqn[k])

        calculate_variable_from_constraint(
            blk.area_column, blk.column_cross_section_area_eqn
        )

        # ---------------------------------------------------------------------
        init_log.info("Step 2: Steady-State isothermal mass balance")

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        assert_optimal_termination(res)
        init_log.info_high("Step 2: {}.".format(idaeslog.condition(res)))

        # ---------------------------------------------------------------------
        init_log.info("Step 3: Interface equilibrium")

        blk.pressure_equil.unfix()
        blk.pressure_at_interface.activate()
        for idx in blk.pressure_at_interface:
            calculate_variable_from_constraint(blk.pressure_equil[idx], blk.pressure_at_interface[idx])

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        assert_optimal_termination(res)
        init_log.info_high("Step 3 complete: {}.".format(idaeslog.condition(res)))

        # ---------------------------------------------------------------------
        init_log.info("Step 4a: Isothermal absorption")
        init_log.info_high("Calculating mass flux")

        # Unfix mass transfer terms
        blk.interphase_mass_transfer.unfix()

        # Activate mass transfer equation in vapor phase
        blk.interphase_mass_transfer_eqn.activate()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        assert_optimal_termination(res)
        init_log.info("Step 4b: Isothermal chemical absorption")
        init_log.info_high("Adding mass transfer to material balances")

        blk.liquid_phase.mass_transfer_term.unfix()
        blk.liquid_mass_transfer_eqn.activate()

        blk.vapor_phase.mass_transfer_term.unfix()
        blk.vapor_mass_transfer_eqn.activate()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        assert_optimal_termination(res)
        init_log.info_high("Step 4 complete: {}.".format(idaeslog.condition(res)))

        # ---------------------------------------------------------------------
        init_log.info("Step 5: Adiabatic chemical absorption")
        init_log.info_high("Isothermal to Adiabatic ")

        # Unfix heat transfer terms
        blk.heat_transfer_coeff.unfix()
        blk.heat_transfer_coeff_eqn.activate()

        for k in blk.heat_transfer_coeff_eqn:
            calculate_variable_from_constraint(
                blk.heat_transfer_coeff[k], blk.heat_transfer_coeff_eqn[k]
            )

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        assert_optimal_termination(res)
        blk.vapor_phase.heat.unfix()
        blk.liquid_phase.heat.unfix()
        blk.vapor_phase.enthalpy_transfer.unfix()
        blk.liquid_phase.enthalpy_transfer.unfix()

        # Activate heat transfer equations
        for c in [
            "heat_transfer_eqn1",
            "heat_transfer_eqn2",
            "enthalpy_transfer_eqn1",
            "enthalpy_transfer_eqn2",
        ]:
            getattr(blk, c).activate()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        assert_optimal_termination(res)
        init_log.info_high("Step 5 complete: {}.".format(idaeslog.condition(res)))

        init_log.info("Step 6: Flooding velocity constraints")

        blk.log_flood_H.unfix()
        blk.log_flow_mass_Liq_Vap.unfix()
        blk.fourth_root_flood_H.unfix()
        blk.gas_velocity_fld.unfix()
        blk.log_gas_velocity_fld.unfix()
        blk.flood_fraction.unfix()

        blk.flood_H_eqn.activate()
        blk.log_flow_mass_Liq_Vap_eqn.activate()
        blk.flood_H_eqn.activate()
        blk.fourth_root_flood_H_eqn.activate()
        blk.log_gas_velocity_fld_eqn.activate()
        blk.gas_velocity_fld_eqn.activate()
        blk.flood_fraction_eqn.activate()

        for t in blk.flowsheet().time:
            for x in blk.vapor_phase.length_domain:
                if x == blk.vapor_phase.length_domain.first():
                    pass
                else:
                    calculate_variable_from_constraint(
                        blk.log_flood_H[t, x], blk.flood_H_eqn[t, x]
                    )
                    calculate_variable_from_constraint(
                        blk.log_flow_mass_Liq_Vap[t, x],
                        blk.log_flow_mass_Liq_Vap_eqn[t, x],
                    )
                    calculate_variable_from_constraint(
                        blk.log_gas_velocity_fld[t, x],
                        blk.log_gas_velocity_fld_eqn[t, x],
                    )

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee, symbolic_solver_labels=True)
        assert_optimal_termination(res)
        init_log.info_high("Step 6 complete: {}.".format(idaeslog.condition(res)))

        # ---------------------------------------------------------------------
        init_log.info("Step 7: Interfacial area constraint")

        blk.area_interfacial.unfix()
        blk.log_area_interfacial.unfix()

        blk.area_interfacial_eqn.activate()
        blk.log_area_interfacial_eqn.activate()

        init_log.info(
            "Initializing interfacial area - degrees_of_freedom = {}".format(
                degrees_of_freedom(blk)
            )
        )           
        
        # Confusing naming convention: log_var_eqn are always of the form exp(log_var) == var.
        # log_area_interfacial is defined from the performance equation area_interfacial_eqn
        # and then area_interfacial is back-calculated from the exponential relationship
        for k in blk.log_area_interfacial_eqn:
            calculate_variable_from_constraint(
                blk.log_area_interfacial[k], blk.area_interfacial_eqn[k]
            )
            calculate_variable_from_constraint(
                blk.area_interfacial[k], blk.log_area_interfacial_eqn[k]
            )

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee, symbolic_solver_labels=True)
        assert_optimal_termination(res)
        init_log.info_high("Step 7 complete: {}.".format(idaeslog.condition(res)))

        # ---------------------------------------------------------------------
        init_log.info("Step 8: Liquid holdup constraint")
        init_log.info(
            "Initializing liquid holdup- degrees_of_freedom = {}".format(
                degrees_of_freedom(blk)
            )
        )

        for c in liquid_holdup_constraints:
            getattr(blk, c).activate()

        blk.log_holdup_liq.unfix()
        blk.holdup_liq.unfix()

        # Same thing with holdup_liq
        for k in blk.log_holdup_liq_eqn:
            calculate_variable_from_constraint(
                blk.log_holdup_liq[k], blk.holdup_liq_eqn[k]
            )
            calculate_variable_from_constraint(
                blk.holdup_liq[k], blk.log_holdup_liq_eqn[k]
            )

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee, symbolic_solver_labels=True)
        assert_optimal_termination(res)
        init_log.info_high("Step 8 complete: {}.".format(idaeslog.condition(res)))

        # ---------------------------------------------------------------------
        init_log.info("Step 9: Vapor phase mass transfer coefficient constraint")
        init_log.info(
            "Initializing mass transfer coefficient (vapor phase) - degrees_of_freedom = {}".format(
                degrees_of_freedom(blk)
            )
        )

        blk.log_holdup_vap.unfix()
        blk.log_mass_transfer_coeff_vap.unfix()
        blk.mass_transfer_coeff_vap.unfix()

        for c in mass_transfer_coeff_vap_constraints:
            getattr(blk, c).activate()

        for t in  blk.flowsheet().time:
            for x in blk.vapor_phase.length_domain:
                if x == blk.liquid_phase.length_domain.first():
                    pass
                else:
                    calculate_variable_from_constraint(
                        blk.log_holdup_vap[t, x], blk.log_holdup_vap_eqn[t, x]
                    )
                    calculate_variable_from_constraint(
                        blk.log_mass_transfer_coeff_vap[t, x, "CO2"],
                        blk.mass_transfer_coeff_vap_eqn[t, x, "CO2"],
                    )
                    calculate_variable_from_constraint(
                        blk.mass_transfer_coeff_vap[t, x, "H2O"],
                        blk.log_mass_transfer_coeff_vap_eqn[t, x, "H2O"],
                    )

        calculate_variable_from_constraint(blk.Cv_ref, blk.log_Cv_ref_eqn)

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee, symbolic_solver_labels=True)
        assert_optimal_termination(res)
        init_log.info_high("Step 9 complete: {}.".format(idaeslog.condition(res)))

        # ---------------------------------------------------------------------
        init_log.info("Step 10: Liquid phase mass transfer coefficient constraint")
        init_log.info(
            "Initializing mass transfer coefficient (liquid phase) - degrees_of_freedom = {}".format(
                degrees_of_freedom(blk)
            )
        )
        blk.mass_transfer_coeff_liq.unfix()

        # blk.log_diffus_liq_comp.unfix()
        blk.log_mass_transfer_coeff_liq.unfix()
        # blk.Cl_ref.unfix()
        blk.mass_transfer_coeff_liq.unfix()

        for c in mass_transfer_coeff_liq_constraints:
            getattr(blk, c).activate()

        for t in blk.flowsheet().time:
            for x in blk.liquid_phase.length_domain:
                if x == blk.liquid_phase.length_domain.last():
                    pass
                else:
                    calculate_variable_from_constraint(
                        blk.log_mass_transfer_coeff_liq[t, x, "CO2"],
                        blk.mass_transfer_coeff_liq_eqn[t, x, "CO2"],
                    )
                    calculate_variable_from_constraint(
                        blk.mass_transfer_coeff_liq[t, x, "CO2"],
                        blk.log_mass_transfer_coeff_liq_eqn[t, x, "CO2"],
                    )

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee, symbolic_solver_labels=True)
        assert_optimal_termination(res)
        init_log.info_high("Step 10 complete: {}.".format(idaeslog.condition(res)))

        # ---------------------------------------------------------------------
        init_log.info(
            "Step 11a Initializing enhancement factor model - degrees_of_freedom = {}".format(
                degrees_of_freedom(blk)
            )
        )

        for constraint in blk.enhancement_factor_constraints:
            constraint.activate()
        for var in blk.enhancement_factor_vars:
            var.unfix()

        blk.config.enhancement_factor_model.initialize_model(
            blk, 
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
        )
                        
        

        init_log.info("Step 11b: Solve model with initialized enhancement factor")

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee, symbolic_solver_labels=True)
        assert_optimal_termination(res)
        init_log.info_high("Step 11 complete: {}.".format(idaeslog.condition(res)))

        # ---------------------------------------------------------------------
        init_log.info("Step 12: Heat transfer coefficient constraint")
        init_log.info(
            "Initializing heat transfer coefficent - degrees_of_freedom = {}".format(
                degrees_of_freedom(blk)
            )
        )
        blk.heat_transfer_coeff_base.unfix()
        blk.log_heat_transfer_coeff_base.unfix()
        blk.heat_transfer_coeff_base_eqn.activate()
        blk.log_heat_transfer_coeff_base_eqn.activate()

        for x in blk.vapor_phase.length_domain:
            if x == blk.vapor_phase.length_domain.first():
                pass
            else:
                calculate_variable_from_constraint(
                    blk.log_heat_transfer_coeff_base[0, x],
                    blk.heat_transfer_coeff_base_eqn[0, x],
                )
                calculate_variable_from_constraint(
                    blk.heat_transfer_coeff_base[0, x],
                    blk.log_heat_transfer_coeff_base_eqn[0, x],
                )

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee, symbolic_solver_labels=True)
        assert_optimal_termination(res)
        init_log.info_high("Step 12 complete: {}.".format(idaeslog.condition(res)))

        # ---------------------------------------------------------------------
        # Release state
        blk.vapor_phase.release_state(flags=vflag)
        blk.liquid_phase.release_state(flags=lflag)

        init_log.info_high(
            "Initialization completed - degrees_of_freedom = {}".format(
                degrees_of_freedom(blk)
            )
        )

        # Check DOF
        if degrees_of_freedom(blk) != initial_dof:
            raise InitializationError(
                f"Degrees of freedom were not returned to their initial value of {initial_dof} at the end "
                f"of initialization. DoF = {degrees_of_freedom(blk)}"
            )

        # Check solver status
        if not check_optimal_termination(res):
            raise InitializationError(
                f"Model failed to initialize successfully. Please check "
                f"the output logs for more information."
            )
