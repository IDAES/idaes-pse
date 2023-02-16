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
Packed Solvent Column Model for MEA systems
"""

# Import Python libraries
import numpy as np

# Import Pyomo libraries
from pyomo.environ import (
    value,
    Var,
    NonNegativeReals,
    Constraint,
    Expression,
    check_optimal_termination,
    exp,
    units as pyunits,
)
from pyomo.util.calc_var_value import calculate_variable_from_constraint

# Import IDAES Libraries
from idaes.core.util.constants import Constants as CONST
from idaes.models_extra.column_models.solvent_column import PackedColumnData

from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core import declare_process_block_class
from idaes.core.util.exceptions import InitializationError
from idaes.core.solvers.get_solver import get_solver
import idaes.logger as idaeslog

import idaes.core.util.scaling as iscale


__author__ = "Paul Akula, John Eslick, Anuja Deshpande, Andrew Lee"


@declare_process_block_class("MEAColumn")
class MEAColumnData(PackedColumnData):
    """
    MEA Packed Column Model Class.
    """

    CONFIG = PackedColumnData.CONFIG()

    def liquid_phase_mass_transfer_model(self):
        """
        Enhancement factor based liquid phase mass transfer model.
        """
        vap_comp = self.config.vapor_phase.property_package.component_list
        liq_comp = self.config.liquid_phase.property_package.component_list
        equilibrium_comp = vap_comp & liq_comp
        solute_comp_list = ["CO2"]

        lunits = (
            self.config.liquid_phase.property_package.get_metadata().get_derived_units
        )

        # Liquid phase equilibrium pressure via Enhancement factor
        self.mass_transfer_coeff_liq = Var(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            solute_comp_list,
            units=lunits("length") / lunits("time"),
            doc="Liquid phase mass transfer coefficient",
        )

        self.enhancement_factor = Var(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            units=pyunits.dimensionless,
            initialize=160,
            doc="Enhancement factor",
        )

        # Intermediate term
        def rule_phi(blk, t, x, j):
            if x == self.vapor_phase.length_domain.first():
                return Expression.Skip
            else:
                zb = self.liquid_phase.length_domain.prev(x)
                return (
                    blk.enhancement_factor[t, zb]
                    * blk.mass_transfer_coeff_liq[t, zb, j]
                    / blk.mass_transfer_coeff_vap[t, x, j]
                )

        self.phi = Expression(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            solute_comp_list,
            rule=rule_phi,
            doc="Equilibrium partial pressure intermediate term",
        )

        @self.Constraint(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            equilibrium_comp,
            doc="Equilibruim partial pressure at interface",
        )
        def pressure_at_interface(blk, t, x, j):
            if x == self.vapor_phase.length_domain.first():
                return Constraint.Skip
            else:
                zb = self.liquid_phase.length_domain.prev(x)
                lprops = blk.liquid_phase.properties[t, zb]
                henrycomp = lprops.params.get_component(j).config.henry_component
                if henrycomp is not None and "Liq" in henrycomp:
                    return blk.pressure_equil[t, x, j] == (
                        (
                            blk.vapor_phase.properties[t, x].mole_frac_comp[j]
                            * pyunits.convert(
                                blk.vapor_phase.properties[t, x].pressure,
                                to_units=lunits("pressure"),
                            )
                            + blk.phi[t, x, j]
                            * lprops.conc_mol_phase_comp_true["Liq", j]
                        )
                        / (
                            1
                            + blk.phi[t, x, j]
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

        # ---------------------------------------------------------------------
        # Unit level sets
        vap_comp = self.config.vapor_phase.property_package.component_list
        liq_comp = self.config.liquid_phase.property_package.component_list
        equilibrium_comp = vap_comp & liq_comp
        solute_comp_list = ["CO2"]

        lunits = (
            self.config.liquid_phase.property_package.get_metadata().get_derived_units
        )

        # Velocity calculations
        self.velocity_vap = Var(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            domain=NonNegativeReals,
            initialize=2,
            units=pyunits.m / pyunits.s,
            doc="Vapor superficial velocity",
        )

        self.velocity_liq = Var(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            domain=NonNegativeReals,
            initialize=0.01,
            units=pyunits.m / pyunits.s,
            doc="Liquid superficial velocity",
        )

        def eq_velocity_vap(blk, t, x):
            return (
                blk.velocity_vap[t, x]
                * blk.area_column
                * blk.vapor_phase.properties[t, x].dens_mol
                == blk.vapor_phase.properties[t, x].flow_mol
            )

        self.velocity_vap_eqn = Constraint(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            rule=eq_velocity_vap,
            doc="Vapor superficial velocity",
        )

        def eq_velocity_liq(blk, t, x):
            return (
                blk.velocity_liq[t, x]
                * blk.area_column
                * blk.liquid_phase.properties[t, x].dens_mol
                == blk.liquid_phase.properties[t, x].flow_mol
            )

        self.velocity_liq_eqn = Constraint(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            rule=eq_velocity_liq,
            doc="Liquid superficial velocity",
        )

        # ---------------------------------------------------------------------
        # Vapor-liquid heat transfer coeff modified by Ackmann factor
        self.heat_transfer_coeff_base = Var(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            initialize=100,
            units=lunits("power") / lunits("temperature") / lunits("length") ** 2,
            doc="Uncorrected vapor-liquid heat transfer coefficient",
        )

        def rule_heat_transfer_coeff_Ack(blk, t, x):
            if x == self.vapor_phase.length_domain.first():
                return Constraint.Skip
            else:
                Ackmann_factor = sum(
                    blk.vapor_phase.properties[t, x].cp_mol_phase_comp["Vap", j]
                    * blk.interphase_mass_transfer[t, x, j]
                    for j in equilibrium_comp
                )
                return blk.heat_transfer_coeff[t, x] == Ackmann_factor / (
                    1
                    - exp(
                        -Ackmann_factor
                        / (
                            blk.heat_transfer_coeff_base[t, x]
                            * blk.area_interfacial[t, x]
                            * blk.area_column
                        )
                    )
                )

        self.heat_transfer_coeff_corr = Constraint(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            rule=rule_heat_transfer_coeff_Ack,
            doc="Vap-Liq heat transfer correction by Ackmann factor",
        )

        # Interfacial Area model ([m2/m3]):
        # Reference: Tsai correlation,regressed by Chinen et al. 2018

        self.area_interfacial_parA = Var(
            initialize=0.6486, units=None, doc="Interfacial area parameter A"
        )

        self.area_interfacial_parB = Var(
            initialize=0.12,
            units=pyunits.dimensionless,
            doc="Interfacial area parameter B",
        )

        def rule_interfacial_area(blk, t, x):
            if x == self.vapor_phase.length_domain.first():
                return Constraint.Skip
            else:
                return (
                    blk.area_interfacial[t, x]
                    == blk.packing_specific_area
                    * blk.area_interfacial_parA
                    * (
                        (
                            pyunits.newton
                            * (pyunits.m ** (2 / 3))
                            * (pyunits.s ** (4 / 3))
                            / pyunits.kilogram
                        )
                        ** blk.area_interfacial_parB
                    )
                    * (
                        (
                            blk.liquid_phase.properties[t, x].mw
                            / blk.liquid_phase.properties[t, x].vol_mol_phase["Liq"]
                            / blk.liquid_phase.properties[t, x].surf_tens_phase["Liq"]
                        )
                        * (blk.velocity_liq[t, x]) ** (4.0 / 3.0)
                    )
                    ** blk.area_interfacial_parB
                )

        self.area_interfacial_constraint = Constraint(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            rule=rule_interfacial_area,
            doc="Specific interfacial area",
        )

        # Liquid holdup model
        # Reference: Tsai correlation,regressed by Chinen et al. 2018

        self.holdup_parA = Var(initialize=24.23, units=None, doc="Holdup parameter A")

        self.holdup_parB = Var(
            initialize=0.6471, units=pyunits.dimensionless, doc="Holdup parameter B"
        )

        def rule_holdup_liq(blk, t, x):
            if x == self.liquid_phase.length_domain.last():
                return Constraint.Skip
            else:
                return (
                    blk.holdup_liq[t, x]
                    == blk.holdup_parA
                    * (
                        (
                            (pyunits.pascal * (pyunits.s ** (10 / 3)))
                            / (pyunits.kilogram * (pyunits.m ** (2 / 3)))
                        )
                        ** blk.holdup_parB
                    )
                    * (
                        blk.velocity_liq[t, x]
                        * (
                            blk.liquid_phase.properties[t, x].visc_d_phase["Liq"]
                            / (
                                blk.liquid_phase.properties[t, x].mw
                                / blk.liquid_phase.properties[t, x].vol_mol_phase["Liq"]
                            )
                        )
                        ** (1 / 3)
                    )
                    ** blk.holdup_parB
                )

        self.holdup_liq_constraint = Constraint(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            rule=rule_holdup_liq,
            doc="Volumetric liquid holdup [-]",
        )

        # Mass transfer coefficients of diffusing components in vapor phase [mol/m2.s.Pa]
        # Mass transfer coefficients, Billet and Schultes (1999) correlation,
        # where parameters are regressed by Chinen et al. (2018).

        self.Cv_ref = Var(
            initialize=0.357,
            doc="""Vapor packing specific constant in Billet and Schultes 
                      volumetric mass transfer coefficient correlation""",
        )

        def rule_mass_transfer_coeff_vap(blk, t, x, j):
            if x == self.vapor_phase.length_domain.first():
                return Constraint.Skip
            else:
                return blk.mass_transfer_coeff_vap[t, x, j] == (
                    1
                    / (
                        pyunits.convert(
                            CONST.gas_constant,
                            pyunits.kilogram
                            * pyunits.m**2
                            / (pyunits.s**2 * pyunits.mol * pyunits.K),
                        )
                        * blk.vapor_phase.properties[t, x].temperature
                    )
                ) * (blk.Cv_ref / (blk.holdup_vap[t, x]) ** 0.5) * (
                    (blk.packing_specific_area / blk.hydraulic_diameter) ** 0.5
                ) * (
                    (blk.vapor_phase.properties[t, x].diffus_phase_comp["Vap", j])
                    ** (2 / 3)
                ) * (
                    (
                        blk.vapor_phase.properties[t, x].visc_d_phase["Vap"]
                        / (
                            blk.vapor_phase.properties[t, x].mw
                            / blk.vapor_phase.properties[t, x].vol_mol_phase["Vap"]
                        )
                    )
                    ** (1 / 3)
                ) * (
                    (
                        (
                            blk.velocity_vap[t, x]
                            * (
                                blk.vapor_phase.properties[t, x].mw
                                / blk.vapor_phase.properties[t, x].vol_mol_phase["Vap"]
                            )
                        )
                        / (
                            blk.packing_specific_area
                            * blk.vapor_phase.properties[t, x].visc_d_phase["Vap"]
                        )
                    )
                    ** (3 / 4)
                )

        self.mass_transfer_coeff_vap_constraint = Constraint(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            equilibrium_comp,
            rule=rule_mass_transfer_coeff_vap,
            doc=" Vapor phase mass transfer coefficient",
        )

        # Mass transfer coefficients of diffusing components in liquid phase  [m/s]
        # Mass transfer coefficients, Billet and Schultes (1999) correlation,
        # where parameters are regressed by Chinen et al. (2018).

        self.Cl_ref = Var(
            initialize=0.5,
            doc="""Liquid packing specific constant in Billet and Schultes 
                      volumetric mass transfer coefficient correlation""",
        )

        def rule_mass_transfer_coeff_liq(blk, t, x, j):
            if x == self.liquid_phase.length_domain.last():
                return Constraint.Skip
            else:
                return (
                    blk.mass_transfer_coeff_liq[t, x, j]
                    == blk.Cl_ref
                    * (12 ** (1 / 6))
                    * (
                        blk.velocity_liq[t, x]
                        * blk.liquid_phase.properties[t, x].diffus_phase_comp["Liq", j]
                        / (blk.hydraulic_diameter * blk.holdup_liq[t, x])
                    )
                    ** 0.5
                )

        self.mass_transfer_coeff_liq_constraint = Constraint(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            solute_comp_list,
            rule=rule_mass_transfer_coeff_liq,
            doc="Liquid phase mass transfer coefficient",
        )

        # Fix mass transfer parameters
        self.Cv_ref.fix()
        self.Cl_ref.fix()

        # Fix interfacial area parameters
        self.area_interfacial_parA.fix()
        self.area_interfacial_parB.fix()

        # Fix liquid holdup parameters
        self.holdup_parA.fix()
        self.holdup_parB.fix()

        # CO2 molar concentration at interface
        def rule_co2_conc_interface(blk, t, x):

            if x == blk.liquid_phase.length_domain.last():
                return Expression.Skip

            else:
                zf = blk.liquid_phase.length_domain.next(x)

                return (
                    blk.pressure_equil[t, zf, "CO2"]
                    / blk.liquid_phase.properties[t, x].henry["Liq", "CO2"]
                )

        self.co2_conc_interface = Expression(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            rule=rule_co2_conc_interface,
            doc="""CO2 concentration at interface""",
        )

        # Enhancement factor model
        # Reference: Jozsef Gaspar,Philip Loldrup Fosbol, (2015)

        self.conc_interface_MEA = Var(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            bounds=(0.5, 1),
            initialize=1,
            units=pyunits.dimensionless,
            doc="""Dimensionless concentration of MEA
                                      at interface """,
        )

        self.sqrt_conc_interface_MEA = Var(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            bounds=(0.5, 1),
            initialize=0.97,
            units=pyunits.dimensionless,
            doc="""Substitute for conc_interface_MEA""",
        )

        def rule_rate_constant(blk, t, x):
            if x == blk.liquid_phase.length_domain.last():
                return Expression.Skip
            else:
                T = blk.liquid_phase.properties[t, x].temperature
                C_MEA = blk.liquid_phase.properties[t, x].conc_mol_phase_comp_true[
                    "Liq", "MEA"
                ]
                C_H2O = blk.liquid_phase.properties[t, x].conc_mol_phase_comp_true[
                    "Liq", "H2O"
                ]

                # Reference: X.Luo et al., Chem. Eng. Sci. (2015)

                return (
                    (
                        2.003e10 * exp(-4742 * pyunits.K / T) * C_MEA
                        + 4.147e6 * exp(-3110 * pyunits.K / T) * C_H2O
                    )
                    * 1e-6
                    * ((pyunits.m) ** 6 / (pyunits.mol**2 * pyunits.s))
                )

        self.rate_constant = Expression(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            rule=rule_rate_constant,
            doc="Second order rate contant [m3/(mol.s)]",
        )

        def rule_hatta_number(blk, t, x):

            if x == blk.liquid_phase.length_domain.last():
                return Expression.Skip
            else:
                return (
                    blk.rate_constant[t, x]
                    * (
                        blk.liquid_phase.properties[t, x].conc_mol_phase_comp_true[
                            "Liq", "MEA"
                        ]
                    )
                    * blk.liquid_phase.properties[t, x].diffus_phase_comp_true[
                        "Liq", "CO2"
                    ]
                ) ** 0.5 / blk.mass_transfer_coeff_liq[t, x, "CO2"]

        self.Hatta = Expression(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            rule=rule_hatta_number,
            doc="""Hatta number""",
        )

        def rule_conc_CO2_bulk(blk, t, x):
            if x == blk.liquid_phase.length_domain.last():
                return Expression.Skip
            else:
                return (
                    blk.liquid_phase.properties[t, x].conc_mol_phase_comp_true[
                        "Liq", "CO2"
                    ]
                    / blk.co2_conc_interface[t, x]
                )

        self.conc_CO2_bulk = Expression(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            rule=rule_conc_CO2_bulk,
            doc="""Dimensionless concentration of CO2,
                                          Driving force term where
                                          Absorption implies conc_CO2_bulk < 1 and
                                          Desorption impies conc_CO2_bulk > 1 """,
        )

        def rule_instantaneous_E(blk, t, x):
            if x == blk.liquid_phase.length_domain.last():
                return Expression.Skip
            else:
                return 1 + (
                    blk.liquid_phase.properties[t, x].diffus_phase_comp_true[
                        "Liq", "MEA"
                    ]
                    * blk.liquid_phase.properties[t, x].conc_mol_phase_comp_true[
                        "Liq", "MEA"
                    ]
                ) / (
                    2
                    * blk.liquid_phase.properties[t, x].diffus_phase_comp_true[
                        "Liq", "CO2"
                    ]
                    * blk.co2_conc_interface[t, x]
                )

        self.instant_E = Expression(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            rule=rule_instantaneous_E,
            doc="Instantaneous Enhancement factor",
        )

        def rule_conc_interface_MEACOO(blk, t, x):
            if x == blk.liquid_phase.length_domain.last():
                return Expression.Skip
            else:
                return 1 + (
                    blk.liquid_phase.properties[t, x].diffus_phase_comp_true[
                        "Liq", "MEA"
                    ]
                    * blk.liquid_phase.properties[t, x].conc_mol_phase_comp_true[
                        "Liq", "MEA"
                    ]
                ) * (
                    1
                    - blk.sqrt_conc_interface_MEA[t, x]
                    * blk.sqrt_conc_interface_MEA[t, x]
                ) / (
                    2
                    * blk.liquid_phase.properties[t, x].diffus_phase_comp_true[
                        "Liq", "MEACOO_-"
                    ]
                    * blk.liquid_phase.properties[t, x].conc_mol_phase_comp_true[
                        "Liq", "MEACOO_-"
                    ]
                )

        self.conc_interface_MEACOO = Expression(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            rule=rule_conc_interface_MEACOO,
            doc="Dimensionless concentration of MEACOO-",
        )

        def rule_conc_interface_MEAH(blk, t, x):
            if x == blk.liquid_phase.length_domain.last():
                return Expression.Skip
            else:
                return 1 + (
                    blk.liquid_phase.properties[t, x].diffus_phase_comp_true[
                        "Liq", "MEA"
                    ]
                    * blk.liquid_phase.properties[t, x].conc_mol_phase_comp_true[
                        "Liq", "MEA"
                    ]
                ) * (
                    1
                    - blk.sqrt_conc_interface_MEA[t, x]
                    * blk.sqrt_conc_interface_MEA[t, x]
                ) / (
                    2
                    * blk.liquid_phase.properties[t, x].diffus_phase_comp_true[
                        "Liq", "MEA_+"
                    ]
                    * blk.liquid_phase.properties[t, x].conc_mol_phase_comp_true[
                        "Liq", "MEA_+"
                    ]
                )

        self.conc_interface_MEAH = Expression(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            rule=rule_conc_interface_MEAH,
            doc="Dimensionless concentration of MEA+",
        )

        @self.Expression(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            doc="""Dimensionless concentration of CO2
                                at equilibruim with the bulk """,
        )
        def conc_CO2_equil_bulk(blk, t, x):
            if x == blk.liquid_phase.length_domain.last():
                return Expression.Skip
            else:
                return (
                    blk.conc_CO2_bulk[t, x]
                    * blk.conc_interface_MEAH[t, x]
                    * blk.conc_interface_MEACOO[t, x]
                    / (blk.sqrt_conc_interface_MEA[t, x] ** 4)
                )

        @self.Constraint(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            doc="""Defining a substitute for conc_interface_MEA""",
        )
        def sqrt_conc_interface_MEA_eqn(blk, t, x):
            if x == blk.liquid_phase.length_domain.last():
                return Constraint.Skip
            else:
                return (
                    blk.sqrt_conc_interface_MEA[t, x]
                    == blk.conc_interface_MEA[t, x] ** 0.5
                )

        @self.Constraint(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            doc="""Enhancement factor - function of Hatta number""",
        )
        def enhancement_factor_eqn1(blk, t, x):
            if x == blk.liquid_phase.length_domain.last():
                return Constraint.Skip
            else:
                return blk.enhancement_factor[t, x] * (
                    1 - blk.conc_CO2_bulk[t, x]
                ) == blk.Hatta[t, x] * (blk.sqrt_conc_interface_MEA[t, x]) * (
                    1 - blk.conc_CO2_equil_bulk[t, x]
                )

        @self.Constraint(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            doc="""Enhancement factor - function of instantaneous enhancement factor""",
        )
        def enhancement_factor_eqn2(blk, t, x):
            if x == blk.liquid_phase.length_domain.last():
                return Constraint.Skip
            else:
                return (blk.enhancement_factor[t, x] - 1) * (
                    1 - blk.conc_CO2_bulk[t, x]
                ) == (blk.instant_E[t, x] - 1) * (
                    1
                    - blk.sqrt_conc_interface_MEA[t, x]
                    * blk.sqrt_conc_interface_MEA[t, x]
                )

        @self.Objective()
        def enhancement_factor_obj(blk):
            time_set = self.flowsheet().time
            x_set = blk.liquid_phase.length_domain
            return sum(
                sum(
                    (
                        (blk.enhancement_factor[t, x] - 1)
                        * (1 - blk.conc_CO2_bulk[t, x])
                        - (blk.instant_E[t, x] - 1)
                        * (
                            1
                            - blk.sqrt_conc_interface_MEA[t, x]
                            * blk.sqrt_conc_interface_MEA[t, x]
                        )
                    )
                    ** 2
                    for t in time_set
                )
                for x in x_set
                if x != blk.liquid_phase.length_domain.last()
            )

        self.enhancement_factor_obj.deactivate()
        # Note: the objective function is only activated in the
        # initialization routine in an intermediate step, to improve
        # model convergence. It is not included in the unit model.

        # Heat transfer coefficients, Chilton Colburn  analogy
        # Vapor-liquid heat transfer coefficient [J/m2.s.K]

        def rule_heat_transfer_coeff(blk, t, x):
            if x == blk.vapor_phase.length_domain.first():
                return Constraint.Skip
            else:
                return blk.heat_transfer_coeff_base[
                    t, x
                ] == blk.mass_transfer_coeff_vap[
                    t, x, "CO2"
                ] * blk.vapor_phase.properties[
                    t, x
                ].pressure * blk.vapor_phase.properties[
                    t, x
                ].cp_mol_phase[
                    "Vap"
                ] * (
                    blk.vapor_phase.properties[t, x].therm_cond_phase["Vap"]
                    / (
                        blk.vapor_phase.properties[t, x].dens_mol_phase["Vap"]
                        * blk.vapor_phase.properties[t, x].cp_mol_phase["Vap"]
                        * blk.vapor_phase.properties[t, x].diffus_phase_comp[
                            "Vap", "CO2"
                        ]
                    )
                ) ** (
                    2 / 3
                )

        self.heat_transfer_coeff_base_constraint = Constraint(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            rule=rule_heat_transfer_coeff,
            doc="""vap-liq heat transfer coefficient""",
        )

        # Flood point calculations
        def rule_flood_velocity(blk, t, x):

            if x == blk.vapor_phase.length_domain.first():
                return Expression.Skip

            else:
                x_liq = blk.vapor_phase.length_domain.prev(x)

                dens_liq = (
                    blk.liquid_phase.properties[t, x_liq].mw
                    / blk.liquid_phase.properties[t, x_liq].vol_mol_phase["Liq"]
                )

                dens_vap = (
                    blk.vapor_phase.properties[t, x].mw
                    / blk.vapor_phase.properties[t, x].vol_mol_phase["Vap"]
                )

                H = (
                    blk.liquid_phase.properties[t, x_liq].flow_mass_phase["Liq"]
                    / blk.vapor_phase.properties[t, x].flow_mass_phase["Vap"]
                ) * (dens_vap / dens_liq) ** 0.5

                mu_water = 0.001 * pyunits.Pa * pyunits.s

                return (
                    (
                        CONST.acceleration_gravity
                        * (blk.eps_ref) ** 3
                        / blk.packing_specific_area
                    )
                    * (dens_liq / dens_vap)
                    * (
                        blk.liquid_phase.properties[t, x_liq].visc_d_phase["Liq"]
                        / mu_water
                    )
                    ** (-0.2)
                    * exp(-4 * (H) ** 0.25)
                ) ** 0.5

        self.gas_velocity_fld = Expression(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            rule=rule_flood_velocity,
            doc="Gas velocity at flooding point",
        )

        self.flood_fraction = Var(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            initialize=0.7,
            units=pyunits.dimensionless,
            doc="""Dimensionless flooding fraction""",
        )

        def rule_flood_fraction(blk, t, x):
            if x == blk.vapor_phase.length_domain.first():
                return Constraint.Skip

            else:
                return (
                    blk.flood_fraction[t, x] * blk.gas_velocity_fld[t, x]
                    == blk.velocity_vap[t, x]
                )

        self.flood_fraction_constr = Constraint(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            rule=rule_flood_fraction,
            doc="Flooding fraction (expected to be less than 0.8)",
        )

    # Scaling Routine
    def calculate_scaling_factors_props(blk):
        for x in blk.vapor_phase.length_domain:
            iscale.set_scaling_factor(blk.vapor_phase.properties[0, x].pressure, 1e-5)

            iscale.set_scaling_factor(blk.vapor_phase.properties[0, x].flow_mol, 1e-4)

        for x in blk.liquid_phase.length_domain:
            iscale.set_scaling_factor(blk.liquid_phase.properties[0, x].pressure, 1e-5)

            iscale.set_scaling_factor(blk.liquid_phase.properties[0, x].flow_mol, 1e-4)

            iscale.set_scaling_factor(
                blk.liquid_phase.properties[0, x].mole_frac_phase_comp_true[
                    "Liq", "CO2"
                ],
                1e4,
            )

            iscale.set_scaling_factor(
                blk.liquid_phase.properties[0, x].mole_frac_phase_comp_true[
                    "Liq", "HCO3_-"
                ],
                100,
            )

            iscale.set_scaling_factor(
                blk.liquid_phase.properties[0, x].mole_frac_phase_comp_true[
                    "Liq", "MEACOO_-"
                ],
                10,
            )

            iscale.set_scaling_factor(
                blk.liquid_phase.properties[0, x].mole_frac_phase_comp_true[
                    "Liq", "MEA_+"
                ],
                10,
            )

            iscale.set_scaling_factor(
                blk.liquid_phase.properties[0, x].mole_frac_phase_comp_true[
                    "Liq", "MEA"
                ],
                1,
            )

            iscale.set_scaling_factor(
                blk.liquid_phase.properties[0, x].mole_frac_phase_comp_true[
                    "Liq", "H2O"
                ],
                1,
            )

            iscale.set_scaling_factor(
                blk.liquid_phase.properties[0, x].flow_mol_phase_comp_true[
                    "Liq", "CO2"
                ],
                100,
            )

            iscale.set_scaling_factor(
                blk.liquid_phase.properties[0, x].flow_mol_phase_comp_true[
                    "Liq", "H2O"
                ],
                1e-3,
            )

            iscale.set_scaling_factor(
                blk.liquid_phase.properties[0, x].flow_mol_phase_comp_true[
                    "Liq", "MEA"
                ],
                0.1,
            )

            iscale.set_scaling_factor(
                blk.liquid_phase.properties[0, x].log_k_eq["carbamate"], 1
            )

            iscale.set_scaling_factor(
                blk.liquid_phase.properties[0, x].log_k_eq["bicarbonate"], 1
            )

        for (t, x), v in blk.vapor_phase.properties.items():
            iscale.constraint_scaling_transform(
                v.total_flow_balance,
                iscale.get_scaling_factor(v.flow_mol, default=1, warning=True),
            )

        for (t, x), v in blk.liquid_phase.properties.items():
            for (p, j), c in v.appr_to_true_species.items():
                iscale.constraint_scaling_transform(
                    v.appr_to_true_species[p, j],
                    iscale.get_scaling_factor(
                        v.flow_mol_phase_comp_true[p, j], default=1, warning=True
                    ),
                )
            for (p, j), c in v.true_mole_frac_constraint.items():
                iscale.constraint_scaling_transform(
                    v.true_mole_frac_constraint[p, j],
                    iscale.get_scaling_factor(
                        v.flow_mol_phase_comp_true[p, j], default=1, warning=True
                    ),
                )

            for j, c in v.component_flow_balances.items():
                iscale.constraint_scaling_transform(
                    v.component_flow_balances["CO2"], 100
                )
                iscale.constraint_scaling_transform(v.component_flow_balances["H2O"], 1)
                iscale.constraint_scaling_transform(v.component_flow_balances["MEA"], 1)

        for (t, x), v in blk.liquid_phase.properties.items():
            iscale.constraint_scaling_transform(
                v.total_flow_balance,
                iscale.get_scaling_factor(v.flow_mol, default=1, warning=True),
            )

        for (t, x), v in blk.liquid_phase.properties.items():
            for rxn, c in v.log_k_eq_constraint.items():
                iscale.constraint_scaling_transform(
                    v.log_k_eq_constraint[rxn],
                    iscale.get_scaling_factor(v.log_k_eq[rxn], default=1, warning=True),
                )

                iscale.constraint_scaling_transform(
                    v.inherent_equilibrium_constraint[rxn],
                    iscale.get_scaling_factor(v.log_k_eq[rxn], default=1, warning=True),
                )

    def calculate_scaling_factors_control_vol(blk):

        # Scale control volume level variables
        for x in blk.vapor_phase.length_domain:
            iscale.set_scaling_factor(blk.vapor_phase.heat[0, x], 1e-6)

            iscale.set_scaling_factor(blk.vapor_phase.enthalpy_transfer[0, x], 1e-6)

            iscale.set_scaling_factor(blk.vapor_phase._enthalpy_flow[0, x, "Vap"], 1e-6)

            iscale.set_scaling_factor(
                blk.vapor_phase.enthalpy_flow_dx[0, x, "Vap"], 1e-6
            )

        for x in blk.liquid_phase.length_domain:
            iscale.set_scaling_factor(
                blk.liquid_phase._enthalpy_flow[0, x, "Liq"], 1e-8
            )

            iscale.set_scaling_factor(
                blk.liquid_phase.enthalpy_flow_dx[0, x, "Liq"], 1e-6
            )

            iscale.set_scaling_factor(blk.liquid_phase.enthalpy_transfer[0, x], 1e-6)

            iscale.set_scaling_factor(blk.liquid_phase.heat[0, x], 1e-6)

        for (t, x, p, j), v in blk.vapor_phase._flow_terms.items():
            iscale.set_scaling_factor(blk.vapor_phase._flow_terms[t, x, p, j], 0.01)

        for (t, x, p, j), v in blk.vapor_phase.material_flow_dx.items():
            if x != 0:
                iscale.set_scaling_factor(
                    blk.vapor_phase.material_flow_dx[t, x, p, j], 0.01
                )
            else:
                pass

        for (t, x, p, j), v in blk.liquid_phase._flow_terms.items():
            iscale.set_scaling_factor(blk.liquid_phase._flow_terms[t, x, p, j], 0.01)

        for (t, x, p, j), v in blk.liquid_phase.material_flow_dx.items():
            if x != 1:
                if j != "MEA":
                    iscale.set_scaling_factor(
                        blk.liquid_phase.material_flow_dx[t, x, p, j], 0.01
                    )
            else:
                pass

        for (
            t,
            x,
            p,
            j,
        ), v in blk.vapor_phase.material_flow_linking_constraints.items():
            iscale.constraint_scaling_transform(
                v,
                iscale.get_scaling_factor(
                    blk.vapor_phase._flow_terms[t, x, p, j], default=1, warning=True
                ),
            )

        for (t, x, j), v in blk.vapor_phase.material_balances.items():
            iscale.constraint_scaling_transform(
                v,
                iscale.get_scaling_factor(
                    blk.vapor_phase.material_flow_dx[t, x, "Vap", j],
                    default=1,
                    warning=True,
                ),
            )

        for (t, x, p), v in blk.vapor_phase.enthalpy_flow_linking_constraint.items():
            iscale.constraint_scaling_transform(
                v,
                iscale.get_scaling_factor(
                    blk.vapor_phase._enthalpy_flow[t, x, p], default=1, warning=True
                ),
            )

        for (t, x), v in blk.vapor_phase.enthalpy_balances.items():
            iscale.constraint_scaling_transform(
                v,
                iscale.get_scaling_factor(
                    blk.vapor_phase.enthalpy_flow_dx[t, x, "Vap"],
                    default=1,
                    warning=True,
                ),
            )

        for (t, x), v in blk.vapor_phase.pressure_balance.items():
            iscale.constraint_scaling_transform(v, 1)

        for (t, x, p, j), v in blk.vapor_phase.material_flow_dx_disc_eq.items():
            iscale.constraint_scaling_transform(
                v,
                iscale.get_scaling_factor(
                    blk.vapor_phase.material_flow_dx[t, x, p, j],
                    default=1,
                    warning=True,
                ),
            )

        for (t, x, p), v in blk.vapor_phase.enthalpy_flow_dx_disc_eq.items():
            iscale.constraint_scaling_transform(
                v,
                iscale.get_scaling_factor(
                    blk.vapor_phase.enthalpy_flow_dx[t, x, p], default=1, warning=True
                ),
            )

        for (t, x), v in blk.vapor_phase.pressure_dx_disc_eq.items():
            iscale.constraint_scaling_transform(v, 1)

        for (
            t,
            x,
            p,
            j,
        ), v in blk.liquid_phase.material_flow_linking_constraints.items():
            iscale.constraint_scaling_transform(
                v,
                iscale.get_scaling_factor(
                    blk.liquid_phase._flow_terms[t, x, p, j], default=1, warning=True
                ),
            )

        for (t, x, j), v in blk.liquid_phase.material_balances.items():
            iscale.constraint_scaling_transform(
                v,
                iscale.get_scaling_factor(
                    blk.liquid_phase.material_flow_dx[t, x, "Liq", j],
                    default=1,
                    warning=True,
                ),
            )

        for (t, x, p), v in blk.liquid_phase.enthalpy_flow_linking_constraint.items():
            iscale.constraint_scaling_transform(
                v,
                iscale.get_scaling_factor(
                    blk.liquid_phase._enthalpy_flow[t, x, p], default=1, warning=True
                ),
            )

        for (t, x), v in blk.liquid_phase.enthalpy_balances.items():
            iscale.constraint_scaling_transform(
                v,
                iscale.get_scaling_factor(
                    blk.liquid_phase.enthalpy_flow_dx[t, x, "Liq"],
                    default=1,
                    warning=True,
                ),
            )

        for (t, x, p, j), v in blk.liquid_phase.material_flow_dx_disc_eq.items():
            iscale.constraint_scaling_transform(
                v,
                iscale.get_scaling_factor(
                    blk.liquid_phase.material_flow_dx[t, x, p, j],
                    default=1,
                    warning=True,
                ),
            )

        for (t, x, p), v in blk.liquid_phase.enthalpy_flow_dx_disc_eq.items():
            iscale.constraint_scaling_transform(
                v,
                iscale.get_scaling_factor(
                    blk.liquid_phase.enthalpy_flow_dx[t, x, p], default=1, warning=True
                ),
            )

    def calculate_scaling_factors_unit_model(blk):

        # Scale unit model level variables
        for (t, x, j), v in blk.pressure_equil.items():
            if x != 0:
                iscale.set_scaling_factor(
                    v,
                    1
                    / value(blk.liquid_phase.properties[t, x].fug_phase_comp["Liq", j]),
                )
            else:
                iscale.set_scaling_factor(v, 1)

        for (t, x, j), v in blk.interphase_mass_transfer.items():
            if x != 0:
                iscale.set_scaling_factor(
                    v, 1 / value(blk.interphase_mass_transfer[t, x, j])
                )
            else:
                iscale.set_scaling_factor(v, 1)

        # ---------------------------------------------------------------------
        # Scale constraints

        for (t, x), v in blk.mechanical_equilibrium.items():
            iscale.constraint_scaling_transform(
                v,
                iscale.get_scaling_factor(
                    blk.liquid_phase.properties[t, x].pressure, default=1, warning=False
                ),
            )

        for (t, x, j), v in blk.pressure_at_interface.items():
            iscale.constraint_scaling_transform(
                v,
                iscale.get_scaling_factor(
                    blk.pressure_equil[t, x, j], default=1, warning=False
                ),
            )

        for (t, x), v in blk.heat_transfer_eqn1.items():
            iscale.constraint_scaling_transform(
                v,
                iscale.get_scaling_factor(
                    blk.vapor_phase.heat[t, x], default=1, warning=True
                ),
            )

        for (t, x), v in blk.heat_transfer_eqn2.items():
            iscale.constraint_scaling_transform(
                v,
                iscale.get_scaling_factor(
                    blk.liquid_phase.heat[t, x], default=1, warning=True
                ),
            )

        for (t, x), v in blk.enthalpy_transfer_eqn1.items():
            iscale.constraint_scaling_transform(
                v,
                iscale.get_scaling_factor(
                    blk.vapor_phase.enthalpy_transfer[t, x], default=1, warning=True
                ),
            )

        for (t, x), v in blk.enthalpy_transfer_eqn2.items():
            iscale.constraint_scaling_transform(
                v,
                iscale.get_scaling_factor(
                    blk.liquid_phase.enthalpy_transfer[t, x], default=1, warning=True
                ),
            )

    def set_init_values_correlation_vars(blk, nfe):
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
                blk.enhancement_factor[0, x].fix(1)
            else:
                blk.mass_transfer_coeff_liq[0, x, "CO2"].fix(k_l_co2_values[i])
                blk.holdup_liq[0, x].fix(0.01)
                blk.enhancement_factor[0, x].fix(enhancement_factor_values[i])

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
        if "bound_push" not in optarg:
            optarg["bound_push"] = 1e-10
        if "max_iter" not in optarg:
            optarg["max_iter"] = 100

        opt = get_solver(solver, optarg)

        blk.calculate_scaling_factors_props()
        blk.calculate_scaling_factors_control_vol()

        unit_constraints = [
            "pressure_at_interface",
            "interphase_mass_transfer_eqn",
            "liquid_mass_transfer_eqn",
            "vapor_mass_transfer_eqn",
            "heat_transfer_coeff_corr",
            "heat_transfer_eqn1",
            "heat_transfer_eqn2",
            "enthalpy_transfer_eqn1",
            "enthalpy_transfer_eqn2",
        ]

        velocity_constraints = ["velocity_vap_eqn", "velocity_liq_eqn"]

        interfacial_area_constraint = ["area_interfacial_constraint"]

        liquid_holdup_constraint = ["holdup_liq_constraint"]

        mass_transfer_coeff_vap_constraints = ["mass_transfer_coeff_vap_constraint"]

        mass_transfer_coeff_liq_constraints = ["mass_transfer_coeff_liq_constraint"]

        enhancement_factor_constraints = [
            "sqrt_conc_interface_MEA_eqn",
            "enhancement_factor_eqn1",
            "enhancement_factor_eqn2",
        ]

        enhancement_factor_obj = ["enhancement_factor_obj"]

        heat_transfer_coefficient_constraint = ["heat_transfer_coeff_base_constraint"]

        flooding_constraint = ["flood_fraction_constr"]

        # ---------------------------------------------------------------------
        # Deactivate unit model level constraints
        for c in blk.component_objects(Constraint, descend_into=True):
            if c.local_name in unit_constraints:
                c.deactivate()
            if c.local_name in velocity_constraints:
                c.deactivate()
            if c.local_name in interfacial_area_constraint:
                c.deactivate()
            if c.local_name in liquid_holdup_constraint:
                c.deactivate()
            if c.local_name in mass_transfer_coeff_vap_constraints:
                c.deactivate()
            if c.local_name in mass_transfer_coeff_liq_constraints:
                c.deactivate()
            if c.local_name in enhancement_factor_constraints:
                c.deactivate()
            if c.local_name in heat_transfer_coefficient_constraint:
                c.deactivate()
            if c.local_name in flooding_constraint:
                c.deactivate()

        # Fix variables

        nfe = blk.config.finite_elements
        blk.set_init_values_correlation_vars(nfe)

        # Interface pressure
        blk.pressure_equil.fix()

        # Molar flux
        blk.interphase_mass_transfer.fix(0.0)
        blk.vapor_phase.mass_transfer_term.fix(0.0)
        blk.liquid_phase.mass_transfer_term.fix(0.0)

        # Heat transfer rate
        blk.heat_transfer_coeff.fix(0.0)
        blk.vapor_phase.heat.fix(0.0)
        blk.liquid_phase.heat.fix(0.0)
        blk.vapor_phase.enthalpy_transfer.fix(0.0)
        blk.liquid_phase.enthalpy_transfer.fix(0.0)

        # Velocity
        blk.velocity_vap.fix()
        blk.velocity_liq.fix()

        # Enhancement factor
        blk.conc_interface_MEA.fix()
        blk.sqrt_conc_interface_MEA.fix()

        # Flood fraction
        blk.flood_fraction.fix()

        # ---------------------------------------------------------------------
        # Provide state arguments for property package initialization
        vap_comp = blk.config.vapor_phase.property_package.component_list
        liq_apparent_comp = [
            c[1] for c in blk.liquid_phase.properties.phase_component_set
        ]

        if vapor_phase_state_args is None:
            vapor_phase_state_args = {
                "flow_mol": blk.vapor_inlet.flow_mol[0].value,
                "temperature": blk.vapor_inlet.temperature[0].value,
                "pressure": blk.vapor_inlet.pressure[0].value,
                "mole_frac_comp": {
                    j: blk.vapor_inlet.mole_frac_comp[0, j].value for j in vap_comp
                },
            }

        if liquid_phase_state_args is None:
            liquid_phase_state_args = {
                "flow_mol": blk.liquid_inlet.flow_mol[0].value,
                "temperature": blk.liquid_inlet.temperature[0].value,
                "pressure": blk.liquid_inlet.pressure[0].value,
                "mole_frac_comp": {
                    j: blk.liquid_inlet.mole_frac_comp[0, j].value
                    for j in liq_apparent_comp
                },
            }

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

        init_log.info("Step 2: Steady-State isothermal mass balance")

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high("Step 2: {}.".format(idaeslog.condition(res)))

        # ---------------------------------------------------------------------
        init_log.info("Step 3: Interface equilibrium")

        # Activate interface pressure constraint
        blk.pressure_equil.unfix()
        blk.pressure_at_interface.activate()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high("Step 3 complete: {}.".format(idaeslog.condition(res)))

        # ---------------------------------------------------------------------
        init_log.info("Step 4a: Isothermal absoption")
        init_log.info_high("Calculating mass flux")

        # Unfix mass transfer terms
        blk.interphase_mass_transfer.unfix()

        # Activate mass transfer equation in vapor phase
        blk.interphase_mass_transfer_eqn.activate()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)

        init_log.info("Step 4b: Isothermal chemical absoption")
        init_log.info_high("Adding mass transfer to material balances")

        blk.vapor_phase.mass_transfer_term.unfix()
        blk.liquid_phase.mass_transfer_term.unfix()
        blk.vapor_mass_transfer_eqn.activate()
        blk.liquid_mass_transfer_eqn.activate()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high("Step 4 complete: {}.".format(idaeslog.condition(res)))

        # ---------------------------------------------------------------------
        init_log.info("Step 5: Adiabatic chemical absoption")
        init_log.info_high("Isothermal to Adiabatic ")

        # Unfix heat transfer terms
        blk.heat_transfer_coeff.unfix()
        for c in ["heat_transfer_coeff_corr"]:
            getattr(blk, c).activate()

        for k in blk.heat_transfer_coeff_corr:
            calculate_variable_from_constraint(
                blk.heat_transfer_coeff[k], blk.heat_transfer_coeff_corr[k]
            )

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)

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

        init_log.info_high("Step 5 complete: {}.".format(idaeslog.condition(res)))

        # ---------------------------------------------------------------------

        blk.calculate_scaling_factors_unit_model()
        iscale.calculate_scaling_factors(blk)

        # ---------------------------------------------------------------------
        init_log.info("Step 6: Velocity constraints")
        blk.velocity_vap.unfix()
        blk.velocity_liq.unfix()

        for c in velocity_constraints:
            getattr(blk, c).activate()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)

        # ---------------------------------------------------------------------
        init_log.info("Step 7: Interfacial area constraint")
        init_log.info(
            "Initializing interfacial area - degrees_of_freedom = {}".format(
                degrees_of_freedom(blk)
            )
        )
        blk.area_interfacial.unfix()

        for c in interfacial_area_constraint:
            getattr(blk, c).activate()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)

        # Scale variable
        for (t, x), v in blk.area_interfacial.items():
            iscale.set_scaling_factor(v, 1 / value(blk.area_interfacial[t, x]))

        # Scale constraint
        for (t, x), v in blk.area_interfacial_constraint.items():
            iscale.constraint_scaling_transform(
                v, iscale.get_scaling_factor(blk.area_interfacial[t, x])
            )

        # ---------------------------------------------------------------------
        init_log.info("Step 8: Liquid holdup constraint")
        init_log.info(
            "Initializing liquid holdup- degrees_of_freedom = {}".format(
                degrees_of_freedom(blk)
            )
        )
        blk.holdup_liq.unfix()

        for c in liquid_holdup_constraint:
            getattr(blk, c).activate()

        for x in blk.liquid_phase.length_domain:
            if x == blk.liquid_phase.length_domain.last():
                pass
            else:
                calculate_variable_from_constraint(
                    blk.holdup_liq[0, x], blk.holdup_liq_constraint[0, x]
                )

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)

        # Scale variable
        for (t, x), v in blk.holdup_liq.items():
            iscale.set_scaling_factor(v, 1 / value(blk.holdup_liq[t, x]))

        # Scale constraint
        for (t, x), v in blk.holdup_liq_constraint.items():
            iscale.constraint_scaling_transform(
                v, iscale.get_scaling_factor(blk.holdup_liq[t, x])
            )

        # ---------------------------------------------------------------------
        init_log.info("Step 9: Vapor phase mass transfer coefficient constraint")
        init_log.info(
            "Initializing mass transfer coefficient (vapor phase) - degrees_of_freedom = {}".format(
                degrees_of_freedom(blk)
            )
        )
        blk.mass_transfer_coeff_vap.unfix()

        for c in mass_transfer_coeff_vap_constraints:
            getattr(blk, c).activate()

        for x in blk.vapor_phase.length_domain:
            if x == blk.vapor_phase.length_domain.first():
                pass
            else:
                calculate_variable_from_constraint(
                    blk.mass_transfer_coeff_vap[0, x, "CO2"],
                    blk.mass_transfer_coeff_vap_constraint[0, x, "CO2"],
                )
                calculate_variable_from_constraint(
                    blk.mass_transfer_coeff_vap[0, x, "H2O"],
                    blk.mass_transfer_coeff_vap_constraint[0, x, "H2O"],
                )

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)

        # Scale variable
        for (t, x, j), v in blk.mass_transfer_coeff_vap.items():
            iscale.set_scaling_factor(
                v, 1 / value(blk.mass_transfer_coeff_vap[t, x, j])
            )

        # Scale constraint
        for (t, x, j), v in blk.mass_transfer_coeff_vap_constraint.items():
            iscale.constraint_scaling_transform(
                v, iscale.get_scaling_factor(blk.mass_transfer_coeff_vap[t, x, j])
            )

        # ---------------------------------------------------------------------
        init_log.info("Step 10: Liquid phase mass transfer coefficient constraint")
        init_log.info(
            "Initializing mass transfer coefficent (liquid phase) - degrees_of_freedom = {}".format(
                degrees_of_freedom(blk)
            )
        )
        blk.mass_transfer_coeff_liq.unfix()

        for c in mass_transfer_coeff_liq_constraints:
            getattr(blk, c).activate()

        for x in blk.liquid_phase.length_domain:
            if x == blk.liquid_phase.length_domain.last():
                pass
            else:
                calculate_variable_from_constraint(
                    blk.mass_transfer_coeff_liq[0, x, "CO2"],
                    blk.mass_transfer_coeff_liq_constraint[0, x, "CO2"],
                )

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)

        # Scale variable
        for (t, x, j), v in blk.mass_transfer_coeff_liq.items():
            iscale.set_scaling_factor(
                v, 1 / value(blk.mass_transfer_coeff_liq[t, x, "CO2"])
            )

        # Scale constraint
        for (t, x, j), v in blk.mass_transfer_coeff_liq_constraint.items():
            iscale.constraint_scaling_transform(
                v, iscale.get_scaling_factor(blk.mass_transfer_coeff_liq[t, x, "CO2"])
            )

        # ---------------------------------------------------------------------
        init_log.info("Step 11: Enhancement factor model constraints")
        init_log.info(
            "Initializing enhancement factor (liquid phase) - degrees_of_freedom = {}".format(
                degrees_of_freedom(blk)
            )
        )
        blk.conc_interface_MEA.unfix()
        blk.sqrt_conc_interface_MEA.unfix()
        blk.enhancement_factor.unfix()

        for c in enhancement_factor_constraints:
            getattr(blk, c).activate()

        blk.enhancement_factor_eqn2.deactivate()
        blk.enhancement_factor_obj.activate()

        for (t, x), v in blk.enhancement_factor_eqn1.items():
            if x != 1:
                iscale.constraint_scaling_transform(v, 100)
            else:
                pass

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)

        blk.enhancement_factor_eqn2.activate()

        for (t, x), v in blk.enhancement_factor_eqn2.items():
            if x != 1:
                iscale.constraint_scaling_transform(v, 100)
            else:
                pass

        blk.enhancement_factor_obj.deactivate()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)

        # ---------------------------------------------------------------------
        init_log.info("Step 12: Heat transfer coefficient constraint")
        init_log.info(
            "Initializing heat transfer coefficent - degrees_of_freedom = {}".format(
                degrees_of_freedom(blk)
            )
        )
        blk.heat_transfer_coeff_base.unfix()
        blk.heat_transfer_coeff_base_constraint.activate()

        for x in blk.vapor_phase.length_domain:
            if x == blk.vapor_phase.length_domain.first():
                pass
            else:
                calculate_variable_from_constraint(
                    blk.heat_transfer_coeff_base[0, x],
                    blk.heat_transfer_coeff_base_constraint[0, x],
                )

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)

        # Scale variable
        for (t, x), v in blk.heat_transfer_coeff_base.items():
            iscale.set_scaling_factor(v, 1 / value(blk.heat_transfer_coeff_base[t, x]))

        # Scale constraint
        for (t, x), v in blk.heat_transfer_coeff_base_constraint.items():
            iscale.constraint_scaling_transform(
                v, iscale.get_scaling_factor(blk.heat_transfer_coeff_base[t, x])
            )

        # ---------------------------------------------------------------------
        init_log.info("Step 13: Flooding constraint")
        init_log.info(
            "Initializing flooding calculations - degrees_of_freedom = {}".format(
                degrees_of_freedom(blk)
            )
        )

        blk.flood_fraction.unfix()
        blk.flood_fraction_constr.activate()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)

        # ---------------------------------------------------------------------
        # Release state
        blk.vapor_phase.release_state(flags=vflag)
        blk.liquid_phase.release_state(flags=lflag)

        init_log.info(
            "Step 14: Initialization completed - degrees_of_freedom = {}".format(
                degrees_of_freedom(blk)
            )
        )

        # Check DOF
        if degrees_of_freedom(blk) != 0:
            raise InitializationError(
                f"Degrees of freedom were not 0 at the end "
                f"of initialization. DoF = {degrees_of_freedom(blk)}"
            )

        # Check solver status
        if not check_optimal_termination(res):
            raise InitializationError(
                f"Model failed to initialize successfully. Please check "
                f"the output logs for more information."
            )
        else:
            init_log.info("Initialization successful")
