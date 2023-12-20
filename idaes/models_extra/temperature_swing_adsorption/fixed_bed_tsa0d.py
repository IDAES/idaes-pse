#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################

"""
0D Fixed Bed TSA unit model.

The Fixed Bed TSA0d model is implemented as a Temperature Swing Adsorption
(TSA) cycle. The model is an 0D equilibrium-based shortcut model. The model
assumes local adsorption equilibrium and takes into account heat transfer
mechanisms but neglects mass transfer resistances.

TSA cycle is composed of the following four steps:

1. Heating step: Column is in a homogeneous state in saturation with
   the feed. No axial gradient in composition, temperature, or
   pressure, and no pressure drop. The column is heated with one
   open end.
2. Cooling step: Column is cooled down with no in/out flows, modeled
   as a batch with varying pressure and temperature.
3. Pressurization step: Due to a pressure decrease in the cooing
   step, the column is repressurized to atmospheric pressure using
   the feed.
4. Adsorption step: The regenerated column is loaded with CO2 in
   adsorption step producing a N2 stream with high purity. The
   adsorption step is over when the CO2 front reaches the end of
   the column. Therefore, the adsorption time is determined by the
   time needed for a shock wave to travel from the column inlet to
   the column outlet.

Equations written in this model were derived from:

L. Joss, M. Gazzani, M. Hefti, D. Marx, and M. Mazzotti. Temperature
swing adsorption for the recovery of the heavy component: an equilibrium-based
shortcut model. Industrial & Engineering Chemistry Research, 2015, 54(11),
3027-3038
"""
from enum import Enum
from pandas import DataFrame

# Import Pyomo libraries
from pyomo.network import Port
from pyomo.common.config import ConfigValue, In, Bool
from pyomo.environ import (
    Constraint,
    Var,
    Param,
    value,
    Set,
    exp,
    log,
    PositiveReals,
    NonPositiveReals,
    TransformationFactory,
    units,
    Block,
)
from pyomo.dae import ContinuousSet, DerivativeVar, Integral
from pyomo.util.calc_var_value import calculate_variable_from_constraint

# import IDAES core libraries
from idaes.core.util.constants import Constants as const
import idaes.core.util.scaling as iscale
from idaes.core.util.config import is_transformation_method, is_physical_parameter_block
from idaes.core import declare_process_block_class, UnitModelBlockData, useDefault
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.math import smooth_max

from idaes.models.unit_models import SkeletonUnitModel, Heater
from idaes.models.unit_models.pressure_changer import (
    PressureChanger,
    ThermodynamicAssumption,
)

import idaes.logger as idaeslog

__author__ = "Daison Yancy Caballero, Alex Noring"

# Set up logger
_log = idaeslog.getLogger(__name__)


class Adsorbent(Enum):
    """
    Enum for supported adsorbents.
    """

    zeolite_13x = 1
    mmen_mg_mof_74 = 2
    polystyrene_amine = 3


class SteamCalculationType(Enum):
    """
    Enum for steam calculation types.
    """

    none = 0
    simplified = 1
    rigorous = 2


class TransformationScheme(Enum):
    """
    Enum for transformation scheme.
    """

    useDefault = 0
    backward = 1
    forward = 2
    lagrangeRadau = 3


@declare_process_block_class("FixedBedTSA0D")
class FixedBedTSA0DData(UnitModelBlockData):
    """
    Standard FixedBedTSA0D Unit Model Class
    """

    # Create Class ConfigBlock
    CONFIG = UnitModelBlockData.CONFIG()
    CONFIG.declare(
        "adsorbent",
        ConfigValue(
            default=Adsorbent.zeolite_13x,
            domain=In(Adsorbent),
            description="Adsorbent flag",
            doc="""Flag to set adsorbent related properties and isotherms.
- Adsorbent.zeolite_13x (default)
- Adsorbent.mmen_mg_mof_74
- Adsorbent.polystyrene_amine""",
        ),
    )
    CONFIG.declare(
        "number_of_beds",
        ConfigValue(
            default=120,
            domain=int,
            description="Number of beds in fixed bed TSA system",
            doc="""Number of beds in fixed bed TSA system used to split the
mole flow rate at feed,
**default** - 120.
**Valid values:** {
**Int** - number of beds is fixed and defined by user,
**None** - number of beds is calculated for a specific pressure drop
in the column}""",
        ),
    )
    CONFIG.declare(
        "transformation_method",
        ConfigValue(
            default="dae.collocation",
            domain=is_transformation_method,
            description="Method to use for DAE transformation",
            doc="""Method to use to transform domain. Must be a method recognised
by the Pyomo TransformationFactory,
**default** - "dae.collocation".
**Valid values:** {
**"dae.finite_difference"** - Use a finite difference method,
**"dae.collocation"** - Use orthogonal collocation method on finite
elements}""",
        ),
    )
    CONFIG.declare(
        "transformation_scheme",
        ConfigValue(
            default=TransformationScheme.lagrangeRadau,
            domain=In(TransformationScheme),
            description="Scheme to use for DAE transformation",
            doc="""Scheme to use when transforming domain. See Pyomo
documentation for supported schemes,
**default** - TransformationScheme.useDefault.
**Valid values:** {
**TransformationScheme.useDefault** - defaults to
TransformationScheme.backward for finite difference transformation
method, and to TransformationScheme.lagrangeRadau for collocation
transformation method,
**TransformationScheme.backward** - Use a finite difference
transformation method,
**TransformationScheme.forward** - use a finite difference
transformation method,
**TransformationScheme.lagrangeRadau** - use a collocation
transformation method}""",
        ),
    )
    CONFIG.declare(
        "finite_elements",
        ConfigValue(
            default=20,
            domain=int,
            description="Number of finite elements for time domain",
            doc="""Number of finite elements to use when discretizing time
domain (default=20)""",
        ),
    )
    CONFIG.declare(
        "collocation_points",
        ConfigValue(
            default=6,
            domain=int,
            description="Number of collocation points per finite element",
            doc="""Number of collocation points to use per finite element when
discretizing domain with "dae.collocation" (default=6)""",
        ),
    )
    CONFIG.declare(
        "compressor",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Compressor flag",
            doc="""Indicates whether a compressor unit should be added to the
fixed bed TSA system to calculate the required energy needed to provide
the pressure drop in the system.
- False (default): compressor not included
- True           : compressor included""",
        ),
    )
    CONFIG.declare(
        "compressor_properties",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package for compressor",
            doc="""The property package to be used by the compressor unit.
The property package must have N2 and CO2 as the only components.
(default=None)""",
        ),
    )
    CONFIG.declare(
        "steam_calculation",
        ConfigValue(
            default=SteamCalculationType.none,
            domain=In(SteamCalculationType),
            description="Steam calculation flag",
            doc="""Indicates whether a method to estimate the steam flow rate
required in the desorption step should be included.
- SteamCalculationType.none (default): steam calculation method not included
- SteamCalculationType.simplified: a surrogate model is used to
estimate the mass flow rate of steam.
- SteamCalculationType.rigorous: a heater unit model is included in the
TSA system assuming total saturation""",
        ),
    )
    CONFIG.declare(
        "steam_properties",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package for heater unit model",
            doc="""The property package to be used by the heater unit.
The property package must be iapws95.
(default=None)""",
        ),
    )

    def build(self):
        """
        General build method for FixedBedTSA0DData. This method calls a
        number of sub-methods which automate the construction of expected
        attributes of unit models.

        Inheriting models should call `super().build`.

        Args:
            None

        Returns:
            None

        """
        # call UnitModel.build to build default attributes
        super().build()

        # consistency check for transformation method and transformation scheme
        if (
            self.config.transformation_method == "dae.finite_difference"
            and self.config.transformation_scheme == TransformationScheme.useDefault
        ):
            self.config.transformation_scheme = TransformationScheme.backward
        elif (
            self.config.transformation_method == "dae.collocation"
            and self.config.transformation_scheme == TransformationScheme.useDefault
        ):
            self.config.transformation_scheme = TransformationScheme.lagrangeRadau
        elif (
            self.config.transformation_method == "dae.finite_difference"
            and self.config.transformation_scheme != TransformationScheme.backward
            and self.config.transformation_scheme != TransformationScheme.forward
        ):
            raise ConfigurationError(
                "{} invalid value for transformation_scheme argument. "
                "Must be TransformationScheme.backward or "
                "TransformationScheme.forward if transformation_method is "
                "dae.finite_difference.".format(self.name)
            )
        elif (
            self.config.transformation_method == "dae.collocation"
            and self.config.transformation_scheme != TransformationScheme.lagrangeRadau
        ):
            raise ConfigurationError(
                "{} invalid value for transformation_scheme argument."
                "Must be TransformationScheme.lagrangeRadau if "
                "transformation_method is dae.collocation.".format(self.name)
            )

        # consistency check for properties
        if self.config.compressor and self.config.compressor_properties == useDefault:
            raise ConfigurationError(
                "Compressor flag was set to true "
                "but no property package was provided"
            )

        if (
            self.config.steam_calculation == SteamCalculationType.rigorous
            and self.config.steam_properties == useDefault
        ):
            raise ConfigurationError(
                "Steam calculations are set to rigorous "
                "but no property package was provided"
            )

        # determine if number of beds is fixed or calculated
        if self.config.number_of_beds is None:
            self.calculate_beds = True
        else:
            self.calculate_beds = False

        # add required parameters
        self._add_general_parameters()

        # add adsorbent parameters
        if self.config.adsorbent == Adsorbent.zeolite_13x:
            self._add_parameters_zeolite_13x()
        elif self.config.adsorbent == Adsorbent.mmen_mg_mof_74:
            self._add_parameters_mmen_Mg_MOF_74()
        elif self.config.adsorbent == Adsorbent.polystyrene_amine:
            self._add_parameters_polystyrene_amine()

        # add design and operating variables
        self.flow_mol_in_total = Var(
            initialize=0.014,
            domain=PositiveReals,
            units=units.mol / units.s,
            doc="Total mole flow rate at inlet",
        )
        self.mole_frac_in = Var(
            self.isotherm_components,
            initialize=0.05,
            domain=PositiveReals,
            units=units.dimensionless,
            doc="Mole fraction at inlet",
        )
        self.pressure_adsorption = Var(
            initialize=1.0 * 1e5,
            domain=PositiveReals,
            units=units.Pa,
            doc="Adsorption pressure",
        )
        self.temperature_adsorption = Var(
            initialize=310,
            domain=PositiveReals,
            units=units.K,
            doc="Adsorption temperature",
        )
        self.temperature_desorption = Var(
            initialize=430,
            domain=PositiveReals,
            units=units.K,
            doc="Desorption temperature",
        )
        self.temperature_heating = Var(
            initialize=440,
            domain=PositiveReals,
            units=units.K,
            doc="Temperature of heating fluid",
        )
        self.temperature_cooling = Var(
            initialize=300,
            domain=PositiveReals,
            units=units.K,
            doc="Temperature of cooling fluid",
        )
        self.bed_diameter = Var(
            initialize=0.03,
            domain=PositiveReals,
            units=units.m,
            doc="Inner column diameter",
        )
        self.bed_height = Var(
            initialize=1.2,
            domain=PositiveReals,
            units=units.m,
            doc="Column length",
        )
        self.pressure_drop = Var(
            initialize=-3500,
            domain=NonPositiveReals,
            units=units.Pa,
            doc="Pressure drop in the column",
        )
        if self.calculate_beds:
            self.velocity_in = Var(
                initialize=0.5,
                bounds=(0.01, 2.0),
                units=units.m / units.s,
                doc="Interstitial velocity of gas at feed",
            )
        self.velocity_mf = Var(
            initialize=1.0,
            domain=PositiveReals,
            units=units.m / units.s,
            doc="Minimum fluidization velocity",
        )

        # calculate required parameters as expressions
        @self.Expression(doc="Total voidage - inter and intraparticle porosity [-]")
        def total_voidage(b):
            return b.bed_voidage + b.particle_voidage * (1 - b.bed_voidage)

        @self.Expression(doc="Bulk density of bed [kg/m3]")
        def bed_bulk_dens_mass(b):
            return b.dens_mass_sol * (1 - b.total_voidage)

        @self.Expression(doc="Outer column diameter [m]")
        def bed_diameter_outer(b):
            return b.bed_diameter + 0.18 / 100 * units.m

        @self.Expression(doc="Column cross-sectional area [m2]")
        def bed_area(b):
            return const.pi * b.bed_diameter**2 / 4

        @self.Expression(doc="Column volume [m3]")
        def bed_volume(b):
            return b.bed_area * b.bed_height

        @self.Expression(doc="Specific heat exchange surface area [m2/m3]")
        def area_heat_transfer(b):
            return const.pi * b.bed_diameter * b.bed_height / b.bed_volume

        @self.Expression(doc="Volume of the column wall [m3]")
        def wall_volume(b):
            return (
                const.pi
                * (b.bed_diameter_outer**2 - b.bed_diameter**2)
                * b.bed_height
                / 4
            )

        @self.Expression(doc="Mass of adsorbent [tonne of adsorbent]")
        def mass_adsorbent(b):
            return units.convert(
                b.bed_bulk_dens_mass * b.bed_volume, to_units=units.tonne
            )

        @self.Expression(doc="Number of adsorption beds")
        def number_beds_ads(b):
            if b.calculate_beds:
                return (
                    b.flow_mol_in_total
                    * const.gas_constant
                    * b.temperature_adsorption
                    / b.bed_area
                    / b.pressure_adsorption
                    / b.velocity_in
                    / b.total_voidage
                )
            else:
                return b.config.number_of_beds

        if not self.calculate_beds:

            @self.Expression(doc="Interstitial velocity of gas at inlet [m/s]")
            def velocity_in(b):
                return (
                    (b.flow_mol_in_total / b.number_beds_ads)
                    * const.gas_constant
                    * b.temperature_adsorption
                    / b.bed_area
                    / b.pressure_adsorption
                    / b.total_voidage
                )

        @self.Expression(doc="Total mole flow rate per bed at inlet [mol/s]")
        def flow_mol_in_total_bed(b):
            return b.flow_mol_in_total / b.number_beds_ads

        @self.Expression(doc="Mole density of gas at inlet [mol/m3]")
        def dens_mol_gas(b):
            return b.pressure_adsorption / const.gas_constant / b.temperature_adsorption

        @self.Expression(doc="Molar mass of mixture at inlet [kg/mol]")
        def mw_mixture_in(b):
            return sum(b.mw[j] * b.mole_frac_in[j] for j in b.isotherm_components)

        @self.Expression(doc="Total mass flow rate at inlet [kg/s]")
        def flow_mass_in_total(b):
            return b.flow_mol_in_total * b.mw_mixture_in

        @self.Expression(doc="Total mass flow rate per bed at inlet [kg/s]")
        def flow_mass_in_total_bed(b):
            return b.flow_mass_in_total / b.number_beds_ads

        @self.Expression(doc="Mass density of gas at inlet [kg/m3]")
        def dens_mass_gas(b):
            return b.dens_mol_gas * b.mw_mixture_in

        @self.Expression(doc="Reynolds at minimum fluidization [-]")
        def Reynolds_mf(b):
            return b.particle_dia * b.velocity_mf * b.dens_mass_gas / b.visc_d

        # constraint for minimum fluidization velocity
        @self.Constraint(doc="Kunii and Levenspiel equation")
        def velocity_mf_eq(b):
            return (
                1.75
                / b.particle_sphericity
                / b.bed_voidage_mf**3
                * b.Reynolds_mf**2
                + 150
                * (1 - b.bed_voidage_mf)
                / b.particle_sphericity**2
                / b.bed_voidage_mf**3
                * b.Reynolds_mf
                == b.particle_dia**3
                * b.dens_mass_gas
                * (b.dens_mass_sol - b.dens_mass_gas)
                * const.acceleration_gravity
                / b.visc_d**2
            )

        # constraint for pressure drop
        @self.Constraint(doc="Ergun equation for pressure drop")
        def pressure_drop_eq(b):
            return b.pressure_drop == (
                (
                    -150
                    * b.visc_d
                    * (1 - b.bed_voidage) ** 2
                    * b.velocity_in
                    / b.bed_voidage**3
                    / b.particle_dia**2
                    / b.particle_sphericity**2
                    - 1.75
                    * (1 - b.bed_voidage)
                    * b.velocity_in**2
                    * b.dens_mass_gas
                    / b.bed_voidage**3
                    / b.particle_dia
                    / b.particle_sphericity
                )
                * b.bed_height
            )

        # add inlet port for fixed bed TSA
        self._add_inlet_port()

        # build up variables and constraints for heating step
        self._add_heating_step()

        # build up variables and constraints for cooling step
        self._add_cooling_step()

        # build up variables and constraints for pressurization step
        self._add_pressurization_step()

        # build up variables and constraints for adsorption step
        self._add_adsorption_step()

        # discretization of time domain
        self._apply_transformation()

        # add performance equations
        self._make_performance()

        # add outlet ports for fixed bed TSA
        self._add_outlet_port()

        # add emissions equations
        self._emissions()

        # add compressor
        if self.config.compressor:
            self._add_compressor()

        # add steam calculation
        if self.config.steam_calculation != SteamCalculationType.none:
            self._add_steam_calc()

    def _add_general_parameters(self):
        """
        Method to add general required parameters to run fixed bed TSA model.
        This include fluid and wall column properties, component set,
        etc.

        """
        self.component_list = Set(
            initialize=["N2", "CO2", "H2O", "O2"], doc="List of components"
        )
        self.isotherm_components = Set(
            initialize=["CO2", "N2"],
            within=self.component_list,
            doc="List of components: CO2 and N2",
        )
        self.cp_wall = Param(
            initialize=4.0 * 1e6,
            units=units.J / units.m**3 / units.K,
            doc="Heat capacity of wall",
        )
        self.mw = Param(
            self.isotherm_components,
            initialize={"CO2": 0.04401, "N2": 0.02801},
            units=units.kg / units.mol,
            doc="Molar mass",
        )
        self.visc_d = Param(
            initialize=1.9e-5,
            units=units.Pa * units.s,
            doc="Gas phase viscosity",
        )
        self.particle_sphericity = Param(
            initialize=1, units=units.dimensionless, doc="Particle sphericity"
        )
        self.bed_voidage_mf = Param(
            initialize=0.5,
            units=units.dimensionless,
            doc="Minimum fluidization bed voidage",
        )

    def _add_parameters_zeolite_13x(self):
        """
        Method to add adsorbent related parameters to run fixed bed TSA model.
        This method is to add parameters for Zeolite 13X.

        Reference: Hefti, M.; Marx, D.; Joss, L.; Mazzotti, M. Adsorption
        Equilibrium of Binary Mixtures of Carbon Dioxide and Nitrogen on
        Zeolites ZSM-5 and 13X. Microporous Mesoporous Materials, 215, 2014.

        """
        # adsorbent parameters
        self.bed_voidage = Param(
            initialize=0.4,
            units=units.dimensionless,
            doc="Bed voidage - external or interparticle porosity",
        )
        self.particle_voidage = Param(
            initialize=0.5,
            units=units.dimensionless,
            doc="Particle voidage - internal or intraparticle porosity",
        )
        self.heat_transfer_coeff = Param(
            initialize=16.8,
            units=units.J / units.m**2 / units.s / units.K,
            doc="Global heat transfer coefficient bed-wall",
        )
        self.cp_mass_sol = Param(
            initialize=920,
            units=units.J / units.kg / units.K,
            doc="Heat capacity of adsorbent",
        )
        self.dens_mass_sol = Param(
            initialize=2360,
            units=units.kg / units.m**3,
            doc="Density of adsorbent",
        )
        self.particle_dia = Param(
            initialize=2e-3, units=units.m, doc="Particle diameter"
        )
        # isotherm parameters
        self.dh_ads = Param(
            self.isotherm_components,
            initialize={"CO2": -37000, "N2": -18510},
            units=units.J / units.mol,
            doc="Heat of adsorption",
        )
        self.temperature_ref = Param(
            initialize=298.15,
            units=units.K,
            doc="Reference temperature",
        )
        self.saturation_capacity_ref = Param(
            self.isotherm_components,
            initialize={"CO2": 7.268, "N2": 4.051},
            units=units.mol / units.kg,
            doc="Saturation capacity at reference temperature",
        )
        self.saturation_capacity_exponential_factor = Param(
            self.isotherm_components,
            initialize={"CO2": -0.61684, "N2": 0.0},
            units=units.dimensionless,
            doc="Isotherm fitting parameter",
        )
        self.affinity_parameter_preexponential_factor = Param(
            self.isotherm_components,
            initialize={"CO2": 1.129e-4, "N2": 5.8470e-5},
            units=units.bar**-1,
            doc="Affinity parameter pre-exponential factor",
        )
        self.affinity_parameter_characteristic_energy = Param(
            self.isotherm_components,
            initialize={"CO2": 28.389, "N2": 18.4740},
            units=units.kJ / units.mol,
            doc="Characteristic energy for the affinity parameter",
        )
        self.heterogeneity_parameter_ref = Param(
            self.isotherm_components,
            initialize={"CO2": 0.42456, "N2": 0.98624},
            units=units.dimensionless,
            doc="Heterogeneity parameter at reference temperature",
        )
        self.heterogeneity_parameter_alpha = Param(
            self.isotherm_components,
            initialize={"CO2": 0.72378, "N2": 0.0},
            units=units.dimensionless,
            doc="Isotherm fitting parameter",
        )

    def _add_parameters_mmen_Mg_MOF_74(self):
        """
        Method to add adsorbent related parameters to run fixed bed TSA model.
        This method is to add parameters for mmen-Mg-MOF-74.

        Reference: Joss, L.; Hefti, M.; Bjelobrk, Z.; Mazzotti, M.
        Investigating the potential of phase-change materials for CO2
        capture. Faraday Disc, 192, 2016

        """
        # adsorbent parameters
        self.bed_voidage = Param(
            initialize=0.4,
            units=units.dimensionless,
            doc="Bed voidage - external or interparticle porosity",
        )
        self.particle_voidage = Param(
            initialize=0.7,
            units=units.dimensionless,
            doc="Particle voidage - internal or intraparticle porosity",
        )
        self.heat_transfer_coeff = Param(
            initialize=16.8,
            units=units.J / units.m**2 / units.s / units.K,
            doc="Global heat transfer coefficient bed-wall",
        )
        self.cp_mass_sol = Param(
            initialize=1500,
            units=units.J / units.kg / units.K,
            doc="Heat capacity of adsorbent",
        )
        self.dens_mass_sol = Param(
            initialize=3200,
            units=units.kg / units.m**3,
            doc="Density of adsorbent",
        )
        self.particle_dia = Param(
            initialize=2e-3, units=units.m, doc="Particle diameter"
        )
        # isotherm parameters
        self.dh_ads = Param(
            self.isotherm_components,
            initialize={"CO2": -62900, "N2": 0.0},
            units=units.J / units.mol,
            doc="Heat of adsorption - estimated with" "Clausius-Clapeyron relation",
        )
        self.temperature_ref = Param(
            initialize=313.15,
            units=units.K,
            doc="Reference temperature",
        )
        self.lower_saturtion_capacity = Param(
            self.isotherm_components,
            initialize={"CO2": 0.146, "N2": 0.0},
            units=units.mol / units.kg,
            doc="Lower isotherm saturation capacity",
        )
        self.upper_saturtion_capacity = Param(
            self.isotherm_components,
            initialize={"CO2": 3.478, "N2": 0.0},
            units=units.mol / units.kg,
            doc="Upper isotherm saturation capacity",
        )
        self.lower_affinity_preexponential_factor = Param(
            self.isotherm_components,
            initialize={"CO2": 0.009, "N2": 0.0},
            units=units.bar**-1,
            doc="Pre-exponential factor for the lower isotherm affinity parameter",
        )
        self.upper_affinity_preexponential_factor_1 = Param(
            self.isotherm_components,
            initialize={"CO2": 9.00e-07, "N2": 0.0},
            units=units.bar**-1,
            doc="Pre-exponential factor for the upper isotherm site 1 affinity parameter",
        )
        self.upper_affinity_preexponential_factor_2 = Param(
            self.isotherm_components,
            initialize={"CO2": 5.00e-04, "N2": 0.0},
            units=units.mol / units.kg / units.bar,
            doc="Pre-exponential factor for the upper isotherm site 2 affinity parameter",
        )
        self.lower_affinity_characteristic_energy = Param(
            self.isotherm_components,
            initialize={"CO2": 31.0, "N2": 0.0},
            units=units.kJ / units.mol,
            doc="Characteristic energy for the lower isotherm affinity parameter",
        )
        self.upper_affinity_characteristic_energy_1 = Param(
            self.isotherm_components,
            initialize={"CO2": 59.0, "N2": 0.0},
            units=units.kJ / units.mol,
            doc="Characteristic energy for the upper isotherm site 1 affinity parameter",
        )
        self.upper_affinity_characteristic_energy_2 = Param(
            self.isotherm_components,
            initialize={"CO2": 18.0, "N2": 0.0},
            units=units.kJ / units.mol,
            doc="Characteristic energy for the upper isotherm site 2 affinity parameter",
        )
        self.step_width_preexponential_factor = Param(
            self.isotherm_components,
            initialize={"CO2": 1.24e-01, "N2": 0.0},
            units=units.dimensionless,
            doc="Pre-exponential factor for step width",
        )
        self.step_width_exponential_factor = Param(
            self.isotherm_components,
            initialize={"CO2": 0.0, "N2": 0.0},
            units=units.dimensionless,
            doc="Exponential factor for step width",
        )
        self.weighting_function_exponent = Param(
            self.isotherm_components,
            initialize={"CO2": 4.00, "N2": 0.0},
            units=units.dimensionless,
            doc="Isotherm weighting function exponent",
        )
        self.step_partial_pressure_ref = Param(
            self.isotherm_components,
            initialize={"CO2": 0.5 * 1e-3, "N2": 0.0},
            units=units.bar,
            doc="Step partial pressure at reference temperature",
        )
        self.step_enthalpy = Param(
            self.isotherm_components,
            initialize={"CO2": -74.1, "N2": 0.0},
            units=units.kJ / units.mol,
            doc="Enthalpy of phase transition",
        )

    def _add_parameters_polystyrene_amine(self):
        """
        Method to add adsorbent related parameters to run fixed bed TSA model.
        This method is to add parameters for polystyrene functionalized
        with primary amine.

        Elfvinga, J.; Bajamundia, C.; Kauppinena, J.; Sainiob, T. Modelling
        of equilibrium working capacity of PSA, TSA and TVSA processes for
        CO2 adsorption under direct air capture conditions. Journal of CO2
        Utilization, 22, 2017.

        """
        # adsorbent parameters
        self.bed_voidage = Param(
            initialize=0.4,
            units=units.dimensionless,
            doc="Bed voidage - external or interparticle porosity",
        )
        self.particle_voidage = Param(
            initialize=0.5,
            units=units.dimensionless,
            doc="Particle voidage - internal or intraparticle porosity",
        )
        self.heat_transfer_coeff = Param(
            initialize=16.8,
            units=units.J / units.m**2 / units.s / units.K,
            doc="Global heat transfer coefficient bed-wall",
        )
        self.cp_mass_sol = Param(
            initialize=920,
            units=units.J / units.kg / units.K,
            doc="Heat capacity of adsorbent",
        )
        self.dens_mass_sol = Param(
            initialize=2360,
            units=units.kg / units.m**3,
            doc="Density of adsorbent",
        )
        self.particle_dia = Param(
            initialize=2e-3, units=units.m, doc="Particle diameter"
        )
        # isotherm parameters
        self.dh_ads = Param(
            self.isotherm_components,
            initialize={"CO2": -114123, "N2": 0.0},
            units=units.J / units.mol,
            doc="Heat of adsorption",
        )
        self.temperature_ref = Param(
            initialize=298.15,
            units=units.K,
            doc="Reference temperature",
        )
        self.saturation_capacity_ref = Param(
            self.isotherm_components,
            initialize={"CO2": 1.71, "N2": 0.0},
            units=units.mol / units.kg,
            doc="Saturation capacity at reference temperature",
        )
        self.affinity_preexponential_factor = Param(
            self.isotherm_components,
            initialize={"CO2": 1.13e5, "N2": 0.0},
            units=units.bar**-1,
            doc="Pre-exponential factor for the affinity parameter",
        )
        self.toth_constant_ref = Param(
            self.isotherm_components,
            initialize={"CO2": 0.265, "N2": 0.0},
            units=units.dimensionless,
            doc="Toth constant at reference temperature",
        )
        self.toth_constant_alpha = Param(
            self.isotherm_components,
            initialize={"CO2": 0.601, "N2": 0.0},
            units=units.dimensionless,
            doc="Isotherm parameter",
        )
        self.saturation_capacity_exponential_factor = Param(
            self.isotherm_components,
            initialize={"CO2": 4.53, "N2": 0.0},
            units=units.dimensionless,
            doc="Exponential factor for saturation capacity",
        )
        self.affinity_characteristic_energy = Param(
            self.isotherm_components,
            initialize={"CO2": 62.2, "N2": 0.0},
            units=units.kJ / units.mol,
            doc="Characteristic energy for the affinity parameter",
        )

    def _add_inlet_port(self):
        """
        Method to build inlet port object of fixed bed TSA. This method allow
        the unit model to be connected with flue gas streams (stream composed
        by N2, CO2, H2O, O2) contained in other IDAES unit models. This method
        add a constraint to force the fixed bed TSA model to use mole flow
        rates of N2 and CO2 only:
        flow_mol_in_total = flow_mol_in["CO2"] + flow_mol_in["N2"]

        """
        # set up inlet port for fixed bed TSA (exhaust gas stream to TSA inlet)
        self.flow_mol_in = Var(
            self.flowsheet().time,
            self.component_list,
            initialize=10,
            units=units.mol / units.s,
            doc="Component mole flow rate at inlet",
        )
        self.temperature_in = Var(
            self.flowsheet().time,
            initialize=300,
            units=units.K,
            doc="Temperature at inlet",
        )
        self.pressure_in = Var(
            self.flowsheet().time,
            initialize=101325,
            units=units.Pa,
            doc="Pressure at inlet",
        )

        # dictionary containing state of inlet stream
        inlet_dict = {
            "flow_mol_comp": self.flow_mol_in,
            "temperature": self.temperature_in,
            "pressure": self.pressure_in,
        }
        # create port
        self.inlet = Port(
            noruleinit=True, initialize=inlet_dict, doc="Inlet port for fixed bed TSA"
        )
        # add constraints to link variables of CCS and inlet of fixed bed TSA
        tf = self.flowsheet().time.last()

        @self.Constraint(doc="Constraint for mole flow rate at inlet")
        def flow_mol_in_total_eq(b):
            return b.flow_mol_in_total == sum(
                b.flow_mol_in[tf, i] for i in b.isotherm_components
            )

        @self.Constraint(
            self.isotherm_components, doc="Constraint for mole fraction at inlet"
        )
        def mole_frac_in_eq(b, i):
            return b.mole_frac_in[i] == b.flow_mol_in[tf, i] / sum(
                b.flow_mol_in[tf, j] for j in b.isotherm_components
            )

        @self.Constraint(doc="Constraint for pressure at inlet")
        def pressure_in_eq(b):
            return b.pressure_adsorption == b.pressure_in[tf]

    def _add_outlet_port(self):
        """
        Method to build outlet port objects of fixed bed TSA. These outlet
        ports are intended to connect the fixed bed TSA unit model with
        other IDAES units. The ports are populated from mass balance through
        the unit. There are three outlets in the fixed bed TSA system:
        1) a stream that is rich in co2 - to compression (co2 captured)
        2) a stream that is rich in n2 - emissions to atmosphere
        3) a stream that contains H2O and O2 that are not accounted
            for in the fixed bed TSA model

        """

        # set up oulet port #1 for fixed bed TSA (CO2 rich stream)
        self.flow_mol_co2_rich_stream = Var(
            self.flowsheet().time,
            self.isotherm_components,
            initialize=10,
            units=units.mol / units.s,
            doc="Component mole flow rate at CO2 rich stream",
        )
        self.temperature_co2_rich_stream = Var(
            self.flowsheet().time,
            initialize=300,
            units=units.K,
            doc="Temperature at CO2 rich stream",
        )
        self.pressure_co2_rich_stream = Var(
            self.flowsheet().time,
            initialize=101325,
            units=units.Pa,
            doc="Pressure at CO2 rich stream",
        )

        # dictionary containing state of outlet co2_rich_stream
        co2_rich_stream_dict = {
            "flow_mol_comp": self.flow_mol_co2_rich_stream,
            "temperature": self.temperature_co2_rich_stream,
            "pressure": self.pressure_co2_rich_stream,
        }
        # create port
        self.co2_rich_stream = Port(
            noruleinit=True,
            initialize=co2_rich_stream_dict,
            doc="Outlet port for CO2 rich stream",
        )

        # add constraints to populate co2_rich_stream
        @self.Constraint(
            self.flowsheet().time,
            self.isotherm_components,
            doc="Constraint for mole flow rate at co2 rich stream",
        )
        def flow_mol_co2_rich_stream_eq(b, t, i):
            if i == "CO2":
                return (
                    b.flow_mol_co2_rich_stream[t, i]
                    == (
                        b.flue_gas_processed_year_target
                        * 1e9
                        * units.mol
                        / units.Gmol
                        / 8760
                        * units.year
                        / units.hr
                        / 3600
                        * units.hr
                        / units.s
                    )
                    * b.recovery
                    * b.mole_frac_in[i]
                )
            elif i == "N2":
                return b.flow_mol_co2_rich_stream[t, i] == (
                    b.flue_gas_processed_year_target
                    * 1e9
                    * units.mol
                    / units.Gmol
                    / 8760
                    * units.year
                    / units.hr
                    / 3600
                    * units.hr
                    / units.s
                ) * b.recovery * b.mole_frac_in["CO2"] * (1 / b.purity - 1)

        @self.Constraint(
            self.flowsheet().time, doc="Constraint for temperature at co2 rich stream"
        )
        def temperature_co2_rich_stream_eq(b, t):
            return b.temperature_co2_rich_stream[t] == b.temperature_desorption

        @self.Constraint(
            self.flowsheet().time, doc="Constraint for pressure at co2 rich stream"
        )
        def pressure_co2_rich_stream_eq(b, t):
            return b.pressure_co2_rich_stream[t] == b.pressure_adsorption

        # set up oulet port #2 for fixed bed TSA (N2 rich stream)
        self.flow_mol_n2_rich_stream = Var(
            self.flowsheet().time,
            self.isotherm_components,
            initialize=10,
            units=units.mol / units.s,
            doc="Component mole flow rate at N2 rich stream",
        )
        self.temperature_n2_rich_stream = Var(
            self.flowsheet().time,
            initialize=300,
            units=units.K,
            doc="Temperature at N2 rich stream",
        )
        self.pressure_n2_rich_stream = Var(
            self.flowsheet().time,
            initialize=101325,
            units=units.Pa,
            doc="Pressure at N2 rich stream",
        )

        # dictionary containing state of outlet n2_rich_stream
        n2_rich_stream_dict = {
            "flow_mol_comp": self.flow_mol_n2_rich_stream,
            "temperature": self.temperature_n2_rich_stream,
            "pressure": self.pressure_n2_rich_stream,
        }
        # create port
        self.n2_rich_stream = Port(
            noruleinit=True,
            initialize=n2_rich_stream_dict,
            doc="Outlet port for n2 rich stream",
        )

        # add constraints to populate n2_rich_stream
        @self.Constraint(
            self.flowsheet().time,
            self.isotherm_components,
            doc="Constraint for mole flow rate at n2 rich stream",
        )
        def flow_mol_n2_rich_stream_eq(b, t, i):
            return (
                b.flow_mol_n2_rich_stream[t, i]
                == b.flow_mol_in[t, i] - b.flow_mol_co2_rich_stream[t, i]
            )

        @self.Constraint(
            self.flowsheet().time, doc="Constraint for temperature at n2 rich stream"
        )
        def temperature_n2_rich_stream_eq(b, t):
            return b.temperature_n2_rich_stream[t] == b.temperature_adsorption

        @self.Constraint(
            self.flowsheet().time, doc="Constraint for pressure at n2 rich stream"
        )
        def pressure_n2_rich_stream_eq(b, t):
            return b.pressure_n2_rich_stream[t] == b.pressure_adsorption

        # set up oulet port #3 for fixed bed TSA (stream containing H2O and O2)
        # The fixed bed TSA assumes an internal magic dryer and remover of o2
        self.flow_mol_h2o_o2_stream = Var(
            self.flowsheet().time,
            self.component_list - self.isotherm_components,
            initialize=1,
            units=units.mol / units.s,
            doc="Component mole flow rate at H2O+O2 stream",
        )
        self.temperature_h2o_o2_stream = Var(
            self.flowsheet().time,
            initialize=300,
            units=units.K,
            doc="Temperature at H2O+O2 stream",
        )
        self.pressure_h2o_o2_stream = Var(
            self.flowsheet().time,
            initialize=101325,
            units=units.Pa,
            doc="Pressure at H2O+O2 stream",
        )

        # dictionary containing state of outlet h2o_o2_stream
        h2o_o2_stream_dict = {
            "flow_mol_comp": self.flow_mol_h2o_o2_stream,
            "temperature": self.temperature_h2o_o2_stream,
            "pressure": self.pressure_h2o_o2_stream,
        }
        # create port
        self.h2o_o2_stream = Port(
            noruleinit=True,
            initialize=h2o_o2_stream_dict,
            doc="Outlet port for H2O+O2 stream",
        )

        # add constraints to populate h2o_o2_stream
        @self.Constraint(
            self.flowsheet().time,
            self.component_list - self.isotherm_components,
            doc="Constraint for mole flow rate at H2O+O2 rich stream",
        )
        def flow_mol_h2o_o2_stream_eq(b, t, i):
            return b.flow_mol_h2o_o2_stream[t, i] == b.flow_mol_in[t, i]

        @self.Constraint(
            self.flowsheet().time,
            doc="Constraint for temperature at H2O+O2 rich stream",
        )
        def temperature_h2o_o2_stream_eq(b, t):
            return b.temperature_h2o_o2_stream[t] == b.temperature_in[t]

        @self.Constraint(
            self.flowsheet().time, doc="Constraint for pressure at H2O+O2 rich stream"
        )
        def pressure_h2o_o2_stream_eq(b, t):
            return b.pressure_h2o_o2_stream[t] == b.pressure_in[t]

    def _add_heating_step(self):
        """
        Method to add variables and constraints for heating step.

        """
        self.heating = SkeletonUnitModel()

        self.heating.pressure = Param(
            initialize=value(self.pressure_adsorption),
            units=units.Pa,
            doc="Pressure in heating step",
        )

        # time domain and total heating time
        self.heating.time = Var(
            initialize=1e3, units=units.s, doc="Time of heating step"
        )
        self.heating.time_domain = ContinuousSet(
            bounds=(0, 1), doc="Normalized time domain in heating step"
        )

        # variables
        def mole_frac_init(b, t, i):
            return self.mole_frac_in[i]

        self.heating.mole_frac = Var(
            self.heating.time_domain,
            self.isotherm_components,
            bounds=(0, 1),
            initialize=mole_frac_init,
            units=units.dimensionless,
            doc="Mole fraction in heating step",
        )
        self.heating.temperature = Var(
            self.heating.time_domain,
            initialize=self.temperature_adsorption,
            units=units.K,
            doc="Temperature in heating step",
        )
        self.heating.velocity_out = Var(
            self.heating.time_domain,
            initialize=self.velocity_in,
            units=units.m / units.s,
            doc="Outlet velocity in heating step",
        )
        self.heating.loading = Var(
            self.heating.time_domain,
            self.isotherm_components,
            units=units.mol / units.kg,
            doc="Equilibrium loading in heating step",
        )

        # time derivatives
        self.heating.mole_frac_dt = DerivativeVar(
            self.heating.mole_frac,
            wrt=self.heating.time_domain,
            units=units.dimensionless,
            doc="Time derivative of mole fraction in heating step",
        )
        self.heating.temperature_dt = DerivativeVar(
            self.heating.temperature,
            wrt=self.heating.time_domain,
            units=units.K,
            doc="Time derivative of temperature in heating step",
        )
        self.heating.loading_dt = DerivativeVar(
            self.heating.loading,
            wrt=self.heating.time_domain,
            units=units.mol / units.kg,
            doc="Time derivative of equilibrium loading in heating step",
        )

        # partial pressure
        @self.heating.Expression(
            self.heating.time_domain,
            self.isotherm_components,
            doc="Partial pressure in heating step [Pa]",
        )
        def partial_pressure(b, t, i):
            return self._partial_pressure(b.pressure, b.mole_frac[t, i])

        # sum of mole fractions
        @self.heating.Constraint(
            self.heating.time_domain, doc="Sum of mole fraction constraint"
        )
        def sum_mole_frac(b, t):
            return sum(b.mole_frac[t, j] for j in self.isotherm_components) == 1.0

        # component mass balance ODE
        @self.heating.Constraint(
            self.heating.time_domain,
            self.isotherm_components,
            doc="Component mass balance - only for CO2",
        )
        def component_mass_balance_ode(b, t, i):
            if t == b.time_domain.first() or i == "N2":
                return Constraint.Skip
            return (
                self.total_voidage
                * b.pressure
                / const.gas_constant
                / b.temperature[t]
                * b.mole_frac_dt[t, i]
                - self.total_voidage
                * b.pressure
                * b.mole_frac[t, i]
                / const.gas_constant
                / b.temperature[t] ** 2
                * b.temperature_dt[t]
                + self.bed_bulk_dens_mass * b.loading_dt[t, i]
                + b.velocity_out[t]
                * b.mole_frac[t, i]
                * b.pressure
                * b.time
                / self.bed_height
                / const.gas_constant
                / b.temperature[t]
                == 0
            )

        # overall mass balance ODE
        @self.heating.Constraint(self.heating.time_domain, doc="Overall mass balance")
        def overall_mass_balance_ode(b, t):
            if t == b.time_domain.first():
                return Constraint.Skip
            return (
                -self.total_voidage
                * b.pressure
                / const.gas_constant
                / b.temperature[t] ** 2
                * b.temperature_dt[t]
                + self.bed_bulk_dens_mass
                * sum(b.loading_dt[t, j] for j in self.isotherm_components)
                + b.velocity_out[t]
                * b.pressure
                * b.time
                / self.bed_height
                / const.gas_constant
                / b.temperature[t]
                == 0
            )

        # energy balance ODE
        @self.heating.Constraint(self.heating.time_domain, doc="Energy balance")
        def energy_balance_ode(b, t):
            if t == b.time_domain.first():
                return Constraint.Skip
            return (
                self.cp_mass_sol * self.bed_bulk_dens_mass * b.temperature_dt[t]
                - self.bed_bulk_dens_mass
                * sum(
                    -self.dh_ads[j] * b.loading_dt[t, j]
                    for j in self.isotherm_components
                )
                == self.heat_transfer_coeff
                * self.area_heat_transfer
                * (self.temperature_heating - b.temperature[t])
                * b.time
            )

        # component loadings
        @self.heating.Constraint(
            self.heating.time_domain,
            self.isotherm_components,
            doc="Isotherm equation",
        )
        def equil_loading_eq(b, t, i):
            p = {}
            for j in self.isotherm_components:
                p[j] = b.partial_pressure[t, j]
            return b.loading[t, i] == self._equil_loading(i, p, b.temperature[t])

        # initial condition constraints
        t0 = self.heating.time_domain.first()

        @self.heating.Constraint(doc="Initial condition for mole fraction")
        def ic_mole_frac_eq(b):
            return b.mole_frac[t0, "CO2"] == self.mole_frac_in["CO2"]

        @self.heating.Constraint(doc="Initial condition for temperature")
        def ic_temperature_eq(b):
            return b.temperature[t0] == self.temperature_adsorption

        @self.heating.Constraint(doc="Initial condition for outlet velocity")
        def ic_velocity_out_eq(b):
            return b.velocity_out[t0] == 0

        # final condition constraint
        tf = self.heating.time_domain.last()

        @self.heating.Constraint(doc="Final condition for temperature")
        def fc_temperature_eq(b):
            return b.temperature[tf] == self.temperature_desorption

    def _add_cooling_step(self):
        """
        Method to add variables and constraints for cooling step.

        """
        self.cooling = SkeletonUnitModel()

        # time domain and total cooling time
        self.cooling.time = Var(
            initialize=1e3, units=units.s, doc="Time of cooling step"
        )
        self.cooling.time_domain = ContinuousSet(
            bounds=(0, 1), doc="Normalized time domain in cooling step"
        )

        # variables
        def mole_frac_init(b, t, i):
            if i == "CO2":
                return 1.0
            else:
                return 0.0

        self.cooling.mole_frac = Var(
            self.cooling.time_domain,
            self.isotherm_components,
            bounds=(0, 1),
            initialize=mole_frac_init,
            units=units.dimensionless,
            doc="Mole fraction in cooling step",
        )
        self.cooling.temperature = Var(
            self.cooling.time_domain,
            initialize=self.temperature_desorption,
            units=units.K,
            doc="Temperature in cooling step",
        )
        self.cooling.pressure = Var(
            self.cooling.time_domain,
            bounds=(1e1, 1e6),
            initialize=self.pressure_adsorption,
            units=units.Pa,
            doc="Pressure in cooling step",
        )
        self.cooling.loading = Var(
            self.cooling.time_domain,
            self.isotherm_components,
            units=units.mol / units.kg,
            doc="Equilibrium loading in cooling step",
        )

        # auxiliary variables to connect initial state of cooling step
        # with final state of heating step
        self.cooling.mole_frac_heating_end = Var(
            initialize=1.0,
            units=units.dimensionless,
            doc="CO2 mole fraction at end of heating step",
        )

        # time derivatives
        self.cooling.mole_frac_dt = DerivativeVar(
            self.cooling.mole_frac,
            wrt=self.cooling.time_domain,
            units=units.dimensionless,
            doc="Time derivative of CO2 mole fraction in cooling step",
        )
        self.cooling.temperature_dt = DerivativeVar(
            self.cooling.temperature,
            wrt=self.cooling.time_domain,
            units=units.K,
            doc="Time derivative of temperature in cooling step",
        )
        self.cooling.pressure_dt = DerivativeVar(
            self.cooling.pressure,
            wrt=self.cooling.time_domain,
            units=units.Pa,
            doc="Time derivative of pressure in cooling step",
        )
        self.cooling.loading_dt = DerivativeVar(
            self.cooling.loading,
            wrt=self.cooling.time_domain,
            units=units.mol / units.kg,
            doc="Time derivative of equilibrium loading in cooling step",
        )

        # partial pressure
        @self.cooling.Expression(
            self.cooling.time_domain,
            self.isotherm_components,
            doc="Partial pressure in cooling step [Pa]",
        )
        def partial_pressure(b, t, i):
            return self._partial_pressure(b.pressure[t], b.mole_frac[t, i])

        # sum of mole fractions
        @self.cooling.Constraint(
            self.cooling.time_domain, doc="Sum of mole fraction constraint"
        )
        def sum_mole_frac(b, t):
            return sum(b.mole_frac[t, j] for j in self.isotherm_components) == 1.0

        # component mass balance ODE
        @self.cooling.Constraint(
            self.cooling.time_domain,
            self.isotherm_components,
            doc="Component mass balance for CO2",
        )
        def component_mass_balance_ode(b, t, i):
            if t == b.time_domain.first() or i == "N2":
                return Constraint.Skip
            return (
                self.total_voidage
                * b.pressure[t]
                / const.gas_constant
                / b.temperature[t]
                * b.mole_frac_dt[t, i]
                + self.total_voidage
                * b.mole_frac[t, i]
                / const.gas_constant
                / b.temperature[t]
                * b.pressure_dt[t]
                - self.total_voidage
                * b.pressure[t]
                * b.mole_frac[t, i]
                / const.gas_constant
                / b.temperature[t] ** 2
                * b.temperature_dt[t]
                + self.bed_bulk_dens_mass * b.loading_dt[t, i]
                == 0
            )

        # overall mass balance ODE
        @self.cooling.Constraint(self.cooling.time_domain, doc="Overall mass balance")
        def overall_mass_balance_ode(b, t):
            if t == b.time_domain.first():
                return Constraint.Skip
            return (
                self.total_voidage
                / const.gas_constant
                / b.temperature[t]
                * b.pressure_dt[t]
                - self.total_voidage
                * b.pressure[t]
                / const.gas_constant
                / b.temperature[t] ** 2
                * b.temperature_dt[t]
                + self.bed_bulk_dens_mass
                * sum(b.loading_dt[t, j] for j in self.isotherm_components)
                == 0
            )

        # energy balance ODE
        @self.cooling.Constraint(self.cooling.time_domain, doc="Energy balance")
        def energy_balance_ode(b, t):
            if t == b.time_domain.first():
                return Constraint.Skip
            return (
                self.cp_mass_sol * self.bed_bulk_dens_mass * b.temperature_dt[t]
                - self.bed_bulk_dens_mass
                * sum(
                    -self.dh_ads[j] * b.loading_dt[t, j]
                    for j in self.isotherm_components
                )
                == self.heat_transfer_coeff
                * self.area_heat_transfer
                * (self.temperature_cooling - b.temperature[t])
                * b.time
            )

        # component loadings
        @self.cooling.Constraint(
            self.cooling.time_domain,
            self.isotherm_components,
            doc="Isotherm equation",
        )
        def equil_loading_eq(b, t, i):
            p = {}
            for j in self.isotherm_components:
                p[j] = b.partial_pressure[t, j]
            return b.loading[t, i] == self._equil_loading(i, p, b.temperature[t])

        # auxiliary constraints to connect initial state of cooling step
        # with final state of heating step
        @self.cooling.Constraint(
            doc="Constraint to calculate mole fraction at end of heating step"
        )
        def mole_frac_heating_end_eq(b):
            tf = self.heating.time_domain.last()
            return b.mole_frac_heating_end == self.heating.mole_frac[tf, "CO2"]

        # initial condition constraints
        t0 = self.cooling.time_domain.first()

        @self.cooling.Constraint(doc="Initial condition for mole fraction")
        def ic_mole_frac_eq(b):
            return b.mole_frac[t0, "CO2"] == b.mole_frac_heating_end

        @self.cooling.Constraint(doc="Initial condition for temperature")
        def ic_temperature_eq(b):
            return b.temperature[t0] == self.temperature_desorption

        @self.cooling.Constraint(doc="Initial condition for pressure")
        def ic_pressure_eq(b):
            return b.pressure[t0] == self.pressure_adsorption

        # final condition constraint
        tf = self.cooling.time_domain.last()

        @self.cooling.Constraint(doc="Final condition for temperature")
        def fc_temperature_eq(b):
            return b.temperature[tf] == self.temperature_adsorption

    def _add_pressurization_step(self):
        """
        Method to add variables and constraints for pressurization step.

        """
        self.pressurization = SkeletonUnitModel()

        # variables
        self.pressurization.mole_frac = Var(
            self.isotherm_components,
            initialize={"CO2": 0.01, "N2": 0.99},
            units=units.dimensionless,
            doc="Mole fraction in pressurization step",
        )
        self.pressurization.time = Var(
            initialize=1e2, units=units.s, doc="Time of pressurization step"
        )
        self.pressurization.loading = Var(
            self.isotherm_components,
            units=units.mol / units.kg,
            doc="Equilibrium loading in pressurization step",
        )

        # auxiliary variables to connect initial state of pressurization step
        # with final state of cooling step
        self.pressurization.mole_frac_cooling_end = Var(
            initialize=1.0,
            units=units.dimensionless,
            doc="Mole fraction at end of cooling step",
        )
        self.pressurization.pressure_cooling_end = Var(
            initialize=1e3, units=units.Pa, doc="Pressure at end of cooling step"
        )
        self.pressurization.loading_cooling_end = Var(
            self.isotherm_components,
            initialize=1.0,
            units=units.mol / units.kg,
            doc="Equilibrium loading at end of cooling step",
        )

        # partial pressure
        @self.pressurization.Expression(
            self.isotherm_components,
            doc="Partial pressure in pressurization step [Pa]",
        )
        def partial_pressure(b, i):
            return self._partial_pressure(self.pressure_adsorption, b.mole_frac[i])

        # sum of mole fractions
        @self.pressurization.Constraint(doc="Sum of mole fraction constraint")
        def sum_mole_frac(b):
            return sum(b.mole_frac[j] for j in self.isotherm_components) == 1.0

        # mass balance
        @self.pressurization.Constraint(doc="Mass balance in pressurization step")
        def mass_balance_eq(b):
            return (
                self.bed_bulk_dens_mass
                * const.gas_constant
                * self.temperature_adsorption
                / self.total_voidage
                * (
                    (b.loading["CO2"] - b.loading_cooling_end["CO2"])
                    - self.mole_frac_in["CO2"]
                    * sum(
                        (b.loading[j] - b.loading_cooling_end[j])
                        for j in self.isotherm_components
                    )
                )
                + (
                    b.mole_frac["CO2"] * self.pressure_adsorption
                    - b.mole_frac_cooling_end * b.pressure_cooling_end
                )
                - self.mole_frac_in["CO2"]
                * (self.pressure_adsorption - b.pressure_cooling_end)
                == 0
            )

        # component loadings
        @self.pressurization.Constraint(
            self.isotherm_components, doc="Isotherm equation"
        )
        def equil_loading_eq(b, i):
            p = {}
            for j in self.isotherm_components:
                p[j] = b.partial_pressure[j]
            return b.loading[i] == self._equil_loading(
                i, p, self.temperature_adsorption
            )

        # pressurization time
        @self.pressurization.Constraint(doc="Equation for pressurization time")
        def pressurization_time_eq(b):
            return (
                self.bed_height
                / self.velocity_in
                * (
                    const.gas_constant
                    * self.temperature_adsorption
                    * self.bed_bulk_dens_mass
                    / self.total_voidage
                    / self.pressure_adsorption
                    * sum(
                        (b.loading[j] - b.loading_cooling_end[j])
                        for j in self.isotherm_components
                    )
                    + (1 - b.pressure_cooling_end / self.pressure_adsorption)
                )
                == b.time
            )

        # auxiliary constraints to connect initial state of pressurization step
        # with final state of cooling step
        tf = self.cooling.time_domain.last()

        @self.pressurization.Constraint(
            doc="Constraint to calculate mole fraction at end of cooling step"
        )
        def mole_frac_cooling_end_eq(b):
            return b.mole_frac_cooling_end == self.cooling.mole_frac[tf, "CO2"]

        @self.pressurization.Constraint(
            doc="Constraint to calculate pressure at end of cooling step"
        )
        def pressure_cooling_end_eq(b):
            return b.pressure_cooling_end == self.cooling.pressure[tf]

        @self.pressurization.Constraint(
            self.isotherm_components,
            doc="Constraint to calculate equilibrium loading at end of " "cooling step",
        )
        def loading_cooling_end_eq(b, i):
            return b.loading_cooling_end[i] == self.cooling.loading[tf, i]

    def _add_adsorption_step(self):
        """
        Method to add variables and constraints for adsorption step.

        """
        self.adsorption = SkeletonUnitModel()

        # variables
        self.adsorption.time = Var(
            initialize=1e2, units=units.s, doc="Time of adsorption step"
        )
        self.adsorption.loading = Var(
            self.isotherm_components,
            units=units.mol / units.kg,
            doc="Equilibrium loading in adsorption step",
        )

        # auxiliary variables to connect initial state of adsorption step
        # with final state of pressurization step
        self.adsorption.mole_frac_pressurization_end = Var(
            initialize=0.01,
            units=units.dimensionless,
            doc="Mole fraction at end of pressurization step",
        )
        self.adsorption.loading_pressurization_end = Var(
            self.isotherm_components,
            initialize=1.0,
            units=units.mol / units.kg,
            doc="Equilibrium loading at end of pressurization step",
        )

        # partial pressure
        @self.adsorption.Expression(
            self.isotherm_components, doc="Partial pressure in adsorption step [Pa]"
        )
        def partial_pressure(b, i):
            return self._partial_pressure(
                self.pressure_adsorption, self.mole_frac_in[i]
            )

        # component loadings
        @self.adsorption.Constraint(self.isotherm_components, doc="Isotherm equation")
        def equil_loading_eq(b, i):
            p = {}
            for j in self.isotherm_components:
                p[j] = b.partial_pressure[j]
            return b.loading[i] == self._equil_loading(
                i, p, self.temperature_adsorption
            )

        # propagation velocity (shock eq. for two adsorbed components)
        @self.adsorption.Expression(doc="Equation for propagation velocity")
        def ads_prop_vel(b):
            return self.velocity_in / (
                1
                + const.gas_constant
                * self.temperature_adsorption
                * self.bed_bulk_dens_mass
                / self.total_voidage
                / self.pressure_adsorption
                * (
                    (b.loading["CO2"] - b.loading_pressurization_end["CO2"])
                    - b.mole_frac_pressurization_end
                    * sum(
                        (b.loading[j] - b.loading_pressurization_end[j])
                        for j in self.isotherm_components
                    )
                )
                / (self.mole_frac_in["CO2"] * 1.0 - b.mole_frac_pressurization_end)
            )

        # propagation velocity (shock eq. for one adsorbed components)
        @self.adsorption.Expression(doc="Equation for propagation velocity")
        def ads_prop_vel_1(b):
            return self.velocity_in / (
                1
                + const.gas_constant
                * self.temperature_adsorption
                * self.bed_bulk_dens_mass
                / self.total_voidage
                / self.pressure_adsorption
                * ((b.loading["CO2"] - b.loading_pressurization_end["CO2"]))
                / (self.mole_frac_in["CO2"] * 1.0 - b.mole_frac_pressurization_end)
            )

        # adsorption time
        @self.adsorption.Constraint(doc="Equation for adsorption time")
        def adsorption_time_eq(b):
            return self.bed_height / b.ads_prop_vel == b.time

        # auxiliary constraints to connect initial state of adsorption step
        # with final state of pressurization step
        @self.adsorption.Constraint(
            doc="Constraint to calculate mole fraction at end of " "pressurization step"
        )
        def mole_frac_pressurization_end_eq(b):
            return (
                b.mole_frac_pressurization_end == self.pressurization.mole_frac["CO2"]
            )

        @self.adsorption.Constraint(
            self.isotherm_components,
            doc="Constraint to calculate equilibrium loading at end of "
            "pressurization step",
        )
        def loading_pressurization_end_eq(b, i):
            return b.loading_pressurization_end[i] == self.pressurization.loading[i]

    def _apply_transformation(self):
        """
        Method to apply DAE transformation to the time domain.

        """
        schemes = {
            TransformationScheme.backward: "BACKWARD",
            TransformationScheme.forward: "FORWARD",
            TransformationScheme.lagrangeRadau: "LAGRANGE-RADAU",
        }

        if self.config.transformation_method == "dae.finite_difference":
            self.discretizer = TransformationFactory(self.config.transformation_method)
            for time_domain in [
                self.heating.time_domain,
                self.cooling.time_domain,
            ]:
                self.discretizer.apply_to(
                    self,
                    wrt=time_domain,
                    nfe=self.config.finite_elements,
                    scheme=schemes[self.config.transformation_scheme],
                )

        elif self.config.transformation_method == "dae.collocation":
            self.discretizer = TransformationFactory(self.config.transformation_method)
            for time_domain in [
                self.heating.time_domain,
                self.cooling.time_domain,
            ]:
                self.discretizer.apply_to(
                    self,
                    wrt=time_domain,
                    nfe=self.config.finite_elements,
                    ncp=self.config.collocation_points,
                    scheme=schemes[self.config.transformation_scheme],
                )

        if self.config.finite_elements is None:
            raise ConfigurationError(
                "{} was not provided a value for the finite_elements"
                " configuration argument. Please provide a valid value.".format(
                    self.name
                )
            )

        if self.config.collocation_points is None:
            raise ConfigurationError(
                "{} was not provided a value for the collocation_points"
                " configuration argument. Please provide a valid value.".format(
                    self.name
                )
            )

    def _equil_loading(self, i, pressure, temperature):
        """
        Method to add equilibrium loading depending on the adsorbent.

        Keyword Arguments:
            i : component
            pressure : partial pressure of components
            temperature : temperature

        """
        if self.config.adsorbent == Adsorbent.zeolite_13x:
            return self._isotherm_zeolite_13x(i, pressure, temperature)

        elif self.config.adsorbent == Adsorbent.mmen_mg_mof_74:
            return self._isotherm_mmen_Mg_MOF_74(i, pressure, temperature)

        elif self.config.adsorbent == Adsorbent.polystyrene_amine:
            return self._isotherm_polystyrene_amine(i, pressure, temperature)

    # TODO: develop a property package framework for adsorbents
    def _isotherm_zeolite_13x(self, i, pressure, temperature):
        """
        Method to add isotherm for components.
        Isotherm equation: Extended Sips

        Keyword Arguments:
            i : component
            pressure : partial pressure of components
            temperature : temperature

        """
        T = temperature

        saturation_capacity = {}
        affinity_parameter = {}
        heterogeneity_parameter = {}
        p = {}

        for j in self.isotherm_components:

            p[j] = units.convert(pressure[j], to_units=units.bar)

            saturation_capacity[j] = self.saturation_capacity_ref[j] * exp(
                self.saturation_capacity_exponential_factor[j]
                * (T / self.temperature_ref - 1)
            )
            affinity_parameter[j] = self.affinity_parameter_preexponential_factor[
                j
            ] * exp(
                units.convert(
                    self.affinity_parameter_characteristic_energy[j],
                    to_units=units.J / units.mol,
                )
                / const.gas_constant
                / T
            )
            heterogeneity_parameter[j] = self.heterogeneity_parameter_ref[
                j
            ] + self.heterogeneity_parameter_alpha[j] * (T / self.temperature_ref - 1)

        loading = {}

        for j in self.isotherm_components:

            loading[j] = (
                saturation_capacity[j]
                * (affinity_parameter[j] * p[j]) ** heterogeneity_parameter[j]
                / (
                    1
                    + sum(
                        (affinity_parameter[k] * p[k]) ** heterogeneity_parameter[k]
                        for k in self.isotherm_components
                    )
                )
            )

        return loading[i]

    def _isotherm_mmen_Mg_MOF_74(self, i, pressure, temperature):
        """
        Method to add isotherm for components.
        Isotherm equation: Weighted dual site Langmuir (w-DSL) isotherm
        Adsorbent: mmen-Mg(dobpdc): mmen-Mg-MOF-74
                    mmen = N,N"-dimethylethylenediamine
                    dobpdc4^- = 4,4"-dioxido-3,3"-biphenyldicarboxylate

        NOTE: the affinity of the material towards CO2 is not reduced by
              neither N2 nor H2O. N2 adsorption is negligible in
              mmen-Mg(dobpdc). Therefore, CO2 is considered as the only
              adsorbing component

        Keyword Arguments:
            i : component
            pressure : partial pressure of components
            temperature : temperature

        """
        T = temperature
        p = {}
        loading = {}

        for j in self.isotherm_components:
            p[j] = units.convert(pressure[j], to_units=units.bar)

        if i == "CO2":

            lower_affinity_parameter = self.lower_affinity_preexponential_factor[
                i
            ] * exp(
                units.convert(
                    self.lower_affinity_characteristic_energy[i],
                    to_units=units.J / units.mol,
                )
                / const.gas_constant
                / T
            )
            upper_affinity_parameter_1 = self.upper_affinity_preexponential_factor_1[
                i
            ] * exp(
                units.convert(
                    self.upper_affinity_characteristic_energy_1[i],
                    to_units=units.J / units.mol,
                )
                / const.gas_constant
                / T
            )
            upper_affinity_parameter_2 = self.upper_affinity_preexponential_factor_2[
                i
            ] * exp(
                units.convert(
                    self.upper_affinity_characteristic_energy_2[i],
                    to_units=units.J / units.mol,
                )
                / const.gas_constant
                / T
            )

            # lower portion of the isotherm
            lower_loading = (
                self.lower_saturtion_capacity[i]
                * lower_affinity_parameter
                * p[i]
                / (1 + lower_affinity_parameter * p[i])
            )

            # upper portion of the isotherm
            upper_loading = (
                self.upper_saturtion_capacity[i]
                * upper_affinity_parameter_1
                * p[i]
                / (1 + upper_affinity_parameter_1 * p[i])
            ) + upper_affinity_parameter_2 * p[i]

            #  weighting function
            sigma = self.step_width_preexponential_factor[i] * exp(
                self.step_width_exponential_factor[i]
                * (1 / self.temperature_ref - 1 / T)
                * units.K
            )

            pstep = (
                self.step_partial_pressure_ref[i]
                / units.bar
                * exp(
                    (
                        -self.step_enthalpy[i]
                        / const.gas_constant
                        * 1e3
                        * units.J
                        / units.kJ
                    )
                    * (1 / self.temperature_ref - 1 / T)
                )
            )

            w = (
                exp((log(p[i] / units.bar) - log(pstep)) / sigma)
                / (1 + exp((log(p[i] / units.bar) - log(pstep)) / sigma))
            ) ** self.weighting_function_exponent[i]

            loading[i] = lower_loading * (1 - w) + upper_loading * w  # [mol/kg]

        elif i == "N2":
            # no adsorption is assumed of N2 in mmen-Mg(dobpdc)
            loading[i] = 1e-10 * units.mol / units.kg

        return loading[i]

    def _isotherm_polystyrene_amine(self, i, pressure, temperature):
        """
        Method to add isotherm for components.
        Isotherm equation: Toth isotherm
        Adsorbent: polystyrene functionalized with primary amine

        NOTE: the affinity of the material towards CO2 is not reduced by
              neither N2 nor H2O. N2 adsorption is negligible in
              this polystyrene functionalized with primary amine. Therefore,
              CO2 is considered as the only adsorbing component

        Keyword Arguments:
            i : component
            pressure : partial pressure of components
            temperature : temperature

        """

        T = temperature
        p = {}
        loading = {}

        for j in self.isotherm_components:
            p[j] = units.convert(pressure[j], to_units=units.bar)

        if i == "CO2":

            saturation_capacity = self.saturation_capacity_ref[i] * exp(
                self.saturation_capacity_exponential_factor[i]
                * (1 - T / self.temperature_ref)
            )
            affinity = self.affinity_preexponential_factor[i] * exp(
                units.convert(
                    self.affinity_characteristic_energy[i], to_units=units.J / units.mol
                )
                / const.gas_constant
                / self.temperature_ref
                * (self.temperature_ref / T - 1)
            )
            toth_constant = self.toth_constant_ref[i] + self.toth_constant_alpha[i] * (
                1 - self.temperature_ref / T
            )

            loading[i] = (
                saturation_capacity
                * affinity
                * p[i]
                / (1 + (affinity * p[i]) ** toth_constant) ** (1 / toth_constant)
            )

        elif i == "N2":
            # no adsorption is assumed of N2 in this adsorbent
            loading[i] = 1e-10 * units.mol / units.kg

        return loading[i]

    def _partial_pressure(self, P, y):
        """
        Method to calculate partial pressure. Partial pressure
        approximated by the smoothed max operator:
        max(0, x) = 0.5*(x + (x^2 + eps)^0.5)

        """
        return P * smooth_max(0, y)

    def _make_performance(self):
        """
        Method to add all performance equations (performance indicators)
        to unit model.

        Args:
            None

        Returns:
            None

        """
        # performance variables
        self.mole_co2_in = Var(
            initialize=10,
            units=units.mol,
            doc="Amount CO2 moles that goes into the column during "
            "pressurization and adsorption steps",
        )
        self.purity = Var(
            initialize=0.9,
            units=units.dimensionless,
            doc="CO2 purity of stream that goes out the column in heating step",
        )
        self.recovery = Var(
            initialize=0.9,
            units=units.dimensionless,
            doc="Percentage of CO2 recovered",
        )
        self.productivity = Var(
            initialize=10.0,
            units=units.kg / units.tonne / units.hr,
            doc="Productivity",
        )
        self.cycle_time = Var(initialize=1.0, units=units.hr, doc="Total cycle time")
        self.thermal_energy = Var(
            initialize=1.0, units=units.MJ, doc="Thermal energy consumption"
        )
        self.specific_energy = Var(
            initialize=1.0,
            units=units.MJ / units.kg,
            doc="Specific thermal energy consumption",
        )

        # expressions for integrals
        def _mole_co2_out(b, t):
            return (
                self.bed_area
                * self.pressure_adsorption
                / const.gas_constant
                * (b.mole_frac[t, "CO2"] * b.velocity_out[t] / b.temperature[t])
                * b.time
            )

        self.heating.mole_co2_out = Integral(
            self.heating.time_domain,
            rule=_mole_co2_out,
            doc="Moles of CO2 that goes out the column in heating step",
        )

        def _mole_out(b, t):
            return (
                self.bed_area
                * self.pressure_adsorption
                / const.gas_constant
                * (b.velocity_out[t] / b.temperature[t])
                * b.time
            )

        self.heating.mole_out = Integral(
            self.heating.time_domain,
            rule=_mole_out,
            doc="Total moles that goes out the column in heating step",
        )

        # constraints for performance indicators
        @self.Constraint(
            doc="Equation to calculate the amount CO2 moles that goes "
            "into the column during pressurization and adsorption steps [mol]"
        )
        def mole_co2_in_eq(b):
            return (
                b.mole_co2_in
                == (b.pressurization.time + b.adsorption.time)
                * b.bed_area
                * b.velocity_in
                * b.total_voidage
                * b.pressure_adsorption
                * b.mole_frac_in["CO2"]
                / const.gas_constant
                / b.temperature_adsorption
            )

        @self.Constraint(doc="Equation to calculate purity")
        def purity_eq(b):
            return b.purity * b.heating.mole_out == b.heating.mole_co2_out

        @self.Constraint(doc="Equation to calculate recovery")
        def recovery_eq(b):
            return b.recovery * b.mole_co2_in == b.heating.mole_co2_out

        @self.Constraint(doc="Equation to calculate total cycle time [h]")
        def cycle_time_eq(b):
            return b.cycle_time == units.convert(
                (
                    b.heating.time
                    + b.cooling.time
                    + b.pressurization.time
                    + b.adsorption.time
                ),
                to_units=units.hr,
            )

        @self.Constraint(doc="Equation to calculate productivity")
        def productivity_eq(b):
            return (
                b.productivity * b.cycle_time * b.mass_adsorbent
                == b.heating.mole_co2_out * b.mw["CO2"]
            )

        @self.Constraint(doc="Equation to calculate thermal energy consumption")
        def thermal_energy_eq(b):
            tf = b.heating.time_domain.last()
            return (
                units.convert(
                    b.bed_volume
                    * (
                        b.cp_mass_sol
                        * b.bed_bulk_dens_mass
                        * (b.temperature_desorption - b.temperature_adsorption)
                        - b.bed_bulk_dens_mass
                        * sum(
                            b.dh_ads[j]
                            * (b.adsorption.loading[j] - b.heating.loading[tf, j])
                            for j in b.isotherm_components
                        )
                    )
                    + b.wall_volume
                    * b.cp_wall
                    * (b.temperature_desorption - b.temperature_adsorption),
                    to_units=units.MJ,
                )
                == b.thermal_energy
            )

        @self.Constraint(
            doc="Equation to calculate specific thermal energy consumption"
        )
        def specific_thermal_energy_eq(b):
            return (
                b.specific_energy * b.heating.mole_co2_out * b.mw["CO2"]
                == b.thermal_energy
            )

        # expressions for other quantities
        @self.Expression(doc="Number of desorption beds")
        def number_beds_des(b):
            return (
                b.number_beds_ads
                * (b.heating.time + b.cooling.time)
                / (b.adsorption.time + b.pressurization.time)
            )

        @self.Expression(doc="Number of beds")
        def number_beds(b):
            return b.number_beds_ads + b.number_beds_des

        @self.Expression(doc="Heat duty per bed required during the heating step [MW]")
        def heat_duty_bed_heating_step(b):
            return b.thermal_energy / b.heating.time

        @self.Expression(doc="Heat duty per bed required during the cycle [MW]")
        def heat_duty_bed(b):
            return b.thermal_energy / units.convert(b.cycle_time, to_units=units.s)

        @self.Expression(
            doc="Total heat duty of fixed bed TSA system required during "
            "the heating step [MW]"
        )
        def heat_duty_total_heating_step(b):
            return b.heat_duty_bed_heating_step * b.number_beds

        @self.Expression(
            doc="Total heat duty of fixed bed TSA system required during "
            "the cycle [MW]"
        )
        def heat_duty_total(b):
            return b.heat_duty_bed * b.number_beds

        @self.Expression(doc="CO2 captured in one cycle per bed [kg/cycle]")
        def CO2_captured_bed_cycle(b):
            return b.heating.mole_co2_out * b.mw["CO2"]

        @self.Expression(doc="Number of cycles per year")
        def cycles_year(b):
            return 8760 * units.hr / units.year / b.cycle_time

        @self.Expression(doc="Total CO2 captured per year [tonne/year]")
        def total_CO2_captured_year(b):
            return (
                units.convert(b.CO2_captured_bed_cycle, to_units=units.tonne)
                * b.cycles_year
                * b.number_beds
            )

        @self.Expression(doc="Amount of flue gas processed per year [Gmol/year]")
        def flue_gas_processed_year(b):
            return (
                units.convert(
                    b.flow_mol_in_total_bed
                    * (b.pressurization.time + b.adsorption.time)
                    * b.number_beds,
                    to_units=units.Gmol,
                )
                * b.cycles_year
            )

        @self.Expression(
            doc="Amount of flue gas needed to be processed year to keep "
            "continuous process [Gmol/year]"
        )
        def flue_gas_processed_year_target(b):
            return (
                b.flow_mol_in_total
                * 3600
                * units.s
                / units.hr
                * 8760
                * units.hr
                / units.year
                * 1e-9
                * units.Gmol
                / units.mol
            )

    def _emissions(self):
        """
        Method to calculate emissions in fixed bed TSA system.
        to unit model.

        Args:
            None

        Returns:
            None

        """

        @self.Expression(
            doc="Emissions of CO2 vented to atmosphere per year [Gmol/year]"
        )
        def emissions_co2_year(b):
            return b.flue_gas_processed_year * (1 - b.recovery) * b.mole_frac_in["CO2"]

        @self.Expression(doc="Emissions of CO2 vented to atmosphere [mol/s]")
        def emissions_co2(b):
            return (
                b.emissions_co2_year
                * 1e9
                * units.mol
                / units.Gmol
                / 8760
                * units.year
                / units.hr
                / 3600
                * units.hr
                / units.s
            )

        @self.Expression(
            self.isotherm_components,
            doc="Mole fraction of N2 rich stream: stream vented to " "atmosphere [-]",
        )
        def mole_frac_n2_rich_stream(b, i):
            return b.flow_mol_n2_rich_stream[0, i] / sum(
                b.flow_mol_n2_rich_stream[0, j] for j in self.isotherm_components
            )

        @self.Expression(
            doc="Concentration of CO2 in stream emitted to atmosphere [ppm]"
        )
        def emissions_co2_ppm(b):
            return b.mole_frac_n2_rich_stream["CO2"] * 100 * 1e4

    def _calculate_and_fix_variable_from_constraint(
        self,
        variable_list=None,
        constraint_list=None,
        obj=None,
        obj_var=None,
        obj_con=None,
    ):
        """
        Method to calculate variables from a constraints, then fix the
        variables and deactivate the constraints.

        Keyword Arguments:
            variable_list: a list of str with names of variables to be
                calculated.
            constraint_list: a list of str with names of constraints
                used to calculate the variables in variable_list.

        """
        if variable_list is None:
            variable_list = []
        if constraint_list is None:
            constraint_list = []

        if obj is None and obj_var is None and obj_con is None:
            obj = self

        if obj is not None and (obj_var is not None or obj_con is not None):
            raise ConfigurationError(
                "{} only provide obj argument or both obj_var "
                "and obj_con. When only obj is provided, this "
                "method looks for variable_list and constraint_list "
                "in the same object passed. When both obj_var and"
                "obj_con are provided, this method looks for "
                "variable_list in obj_var and constraint_list "
                "in obj_con.".format(self.name)
            )

        if (obj_var is not None and obj_con is None) or (
            obj_con is not None and obj_var is None
        ):
            raise ConfigurationError(
                "{} both obj_var and obj_con need to be provided together.".format(
                    self.name
                )
            )

        v_list = []
        c_list = []

        if obj_var is None and obj_var is None:
            for v in obj.component_objects(Var, descend_into=True):
                if v.local_name in variable_list:
                    v_list.append(v)
            for c in obj.component_objects(Constraint, descend_into=True):
                if c.local_name in constraint_list:
                    c_list.append(c)
        else:
            for v in obj_var.component_objects(Var, descend_into=True):
                if v.local_name in variable_list:
                    v_list.append(v)
            for c in obj_con.component_objects(Constraint, descend_into=True):
                if c.local_name in constraint_list:
                    c_list.append(c)

        v_c_tuple = list(zip(v_list, c_list))

        for k in v_c_tuple:
            if k[1].dim() == 0:
                calculate_variable_from_constraint(k[0], k[1])
            elif k[1].dim() == 1:
                for var_index in k[1].index_set():
                    calculate_variable_from_constraint(k[0][var_index], k[1][var_index])

        for v in v_list:
            v.fix()

        for c in c_list:
            c.deactivate()

    def _add_compressor(self):

        self.compressor = Block()

        self.compressor.unit = PressureChanger(
            property_package=self.config.compressor_properties,
            thermodynamic_assumption=ThermodynamicAssumption.isentropic,
            compressor=True,
        )

        # add constraints to link variables of CCS and compressor systems
        @self.compressor.Constraint(
            self.flowsheet().time,
            self.config.compressor_properties.component_list,
            doc="Constraint for component flow rate at inlet of compressor",
        )
        def flow_mol_in_compressor_eq(b, t, i):
            if i in self.component_list:
                return (
                    b.unit.inlet.flow_mol_comp[t, i] == self.inlet.flow_mol_comp[t, i]
                )
            else:
                return b.unit.inlet.flow_mol_comp[t, i] == 1e-8

        @self.compressor.Constraint(
            self.flowsheet().time,
            doc="Constraint for temperature at inlet of compressor",
        )
        def temperature_in_compressor_eq(b, t):
            return b.unit.inlet.temperature[t] == self.temperature_in[t]

        @self.compressor.Constraint(
            self.flowsheet().time, doc="Constraint for pressure at inlet of compressor"
        )
        def pressure_in_compressor_eq(b, t):
            return b.unit.inlet.pressure[t] == self.pressure_in[t]

        @self.compressor.Constraint(
            self.flowsheet().time,
            doc="Constraint to match pressure drop in fixed bed TSA model and "
            "compressor",
        )
        def pressure_drop_tsa_compressor_eqn(b, t):
            return self.compressor.unit.deltaP[t] * 1e-3 == -self.pressure_drop * 1e-3

        # scaling factors for compression work
        for t in self.flowsheet().time:
            iscale.set_scaling_factor(self.compressor.unit.work_mechanical[t], 1e-5)
            iscale.set_scaling_factor(self.compressor.unit.work_isentropic[t], 1e-5)
            iscale.set_scaling_factor(self.compressor.unit.control_volume.work[t], 1e-5)

    def _add_steam_calc(self):

        self.flow_mass_steam = Var(
            initialize=1.0,
            domain=PositiveReals,
            units=units.kg / units.s,
            doc="Mass flow rate of steam",
        )

        if self.config.steam_calculation == SteamCalculationType.rigorous:

            # add empty block for heater
            self.steam_heater = Block()

            # add heater unit model
            self.steam_heater.unit = Heater(
                property_package=self.config.steam_properties, has_pressure_change=True
            )

            # add constraint for total saturation condition
            @self.steam_heater.unit.Constraint(
                doc="Constraint outlet stream to be liquid (total saturation)"
            )
            def vapor_frac_out_eq(b):
                return b.control_volume.properties_out[0].vapor_frac == 1e-9

            # add constraints to link heat duty of TSA and heater systems
            @self.steam_heater.unit.Constraint(doc="Constraint for heat duty of heater")
            def heat_duty_heater_eq(b):
                return b.heat_duty[0] == -units.convert(
                    self.heat_duty_total, to_units=units.J / units.s
                )

            # fix inlet pressure, enthalpy and pressure drop
            Pin = 510e3 * units.Pa
            hin = 52465.8 * units.J / units.mol  # 500 K

            self.steam_heater.unit.inlet.enth_mol[0].fix(hin)
            self.steam_heater.unit.inlet.pressure[0].fix(Pin)
            self.steam_heater.unit.deltaP.fix(0)

            # scaling factor for heat duty in heater
            for t in self.flowsheet().time:
                iscale.set_scaling_factor(
                    self.steam_heater.unit.control_volume.heat[t], 1e-5
                )

        @self.Constraint(doc="Constraint for mass flow rate of steam")
        def flow_mass_steam_eq(b):
            if self.config.steam_calculation == SteamCalculationType.simplified:
                return (
                    self.flow_mass_steam
                    == (21.594 * self.heat_duty_total / units.MW + 200.42)
                    * units.mol
                    / units.s
                    * 0.018015
                    * units.kg
                    / units.mol
                )

            if self.config.steam_calculation == SteamCalculationType.rigorous:
                return (
                    self.flow_mass_steam
                    == self.steam_heater.unit.control_volume.properties_out[0].flow_mass
                )

    def initialize_build(blk, outlvl=idaeslog.NOTSET, solver=None, optarg=None):
        """
        This method is deprecated.
        """
        raise DeprecationWarning(
            "The initialize_build method has been deprecated, use the "
            "FixedBedTSA0DInitializer class instead."
        )

    def fix_initialization_states(self):
        self.inlet.flow_mol_comp.fix()
        self.inlet.temperature.fix()
        self.inlet.pressure.fix()

    def calculate_scaling_factors(self):

        super().calculate_scaling_factors()

        # scaling factors for heating step
        if iscale.get_scaling_factor(self.heating.mole_frac) is None:
            iscale.set_scaling_factor(self.heating.mole_frac, 1e1)

        if iscale.get_scaling_factor(self.heating.velocity_out) is None:
            iscale.set_scaling_factor(self.heating.velocity_out, 1e2)

        if iscale.get_scaling_factor(self.heating.temperature) is None:
            iscale.set_scaling_factor(self.heating.temperature, 1e-2)

        if iscale.get_scaling_factor(self.heating.temperature_dt) is None:
            iscale.set_scaling_factor(self.heating.temperature_dt, 1e-2)

        for t in self.heating.time_domain:
            if t != self.heating.time_domain.first():
                if (
                    iscale.get_scaling_factor(
                        self.heating.component_mass_balance_ode[t, "CO2"]
                    )
                    is None
                ):
                    iscale.set_scaling_factor(
                        self.heating.component_mass_balance_ode[t, "CO2"], 1e-5
                    )
                if (
                    iscale.get_scaling_factor(self.heating.overall_mass_balance_ode[t])
                    is None
                ):
                    iscale.set_scaling_factor(
                        self.heating.overall_mass_balance_ode[t], 1e-5
                    )
                if (
                    iscale.get_scaling_factor(self.heating.energy_balance_ode[t])
                    is None
                ):
                    iscale.set_scaling_factor(self.heating.energy_balance_ode[t], 1e-6)

        # scaling factors for cooling step
        if iscale.get_scaling_factor(self.cooling.mole_frac) is None:
            iscale.set_scaling_factor(self.cooling.mole_frac, 1e1)

        if iscale.get_scaling_factor(self.cooling.temperature) is None:
            iscale.set_scaling_factor(self.cooling.temperature, 1e-2)

        if iscale.get_scaling_factor(self.cooling.pressure) is None:
            iscale.set_scaling_factor(self.cooling.pressure, 1e-4)

        if iscale.get_scaling_factor(self.cooling.temperature_dt) is None:
            iscale.set_scaling_factor(self.cooling.temperature_dt, 1e-2)

        if iscale.get_scaling_factor(self.cooling.pressure_dt) is None:
            iscale.set_scaling_factor(self.cooling.pressure_dt, 1e-4)

        for t in self.cooling.time_domain:
            if t != self.cooling.time_domain.first():
                if (
                    iscale.get_scaling_factor(
                        self.cooling.component_mass_balance_ode[t, "CO2"]
                    )
                    is None
                ):
                    iscale.set_scaling_factor(
                        self.cooling.component_mass_balance_ode[t, "CO2"], 1e-5
                    )
                if (
                    iscale.get_scaling_factor(self.cooling.overall_mass_balance_ode[t])
                    is None
                ):
                    iscale.set_scaling_factor(
                        self.cooling.overall_mass_balance_ode[t], 1e-5
                    )
                if (
                    iscale.get_scaling_factor(self.cooling.energy_balance_ode[t])
                    is None
                ):
                    iscale.set_scaling_factor(self.cooling.energy_balance_ode[t], 1e-6)

        # scaling for inlet states
        iscale.set_scaling_factor(self.flow_mol_in, 1e-3)
        iscale.set_scaling_factor(self.temperature_in, 1e-2)
        iscale.set_scaling_factor(self.pressure_in, 1e-5)

        # constraint scaling factors
        iscale.constraint_scaling_transform(self.pressure_drop_eq, 1e-5)
        iscale.constraint_scaling_transform(self.velocity_mf_eq, 10)

    def get_var_dict(self):
        """
        Returns a dictionary of important variables

        """

        var_dict = {}

        var_dict["Adsorption temperature [K]"] = self.temperature_adsorption
        var_dict["Desorption temperature [K]"] = self.temperature_desorption
        var_dict["Heating temperature [K]"] = self.temperature_heating
        var_dict["Cooling temperature [K]"] = self.temperature_cooling
        var_dict["Column diameter [m]"] = self.bed_diameter
        var_dict["Column length [m]"] = self.bed_height
        var_dict["Column volume [m3]"] = self.bed_volume
        var_dict["CO2 mole fraction at feed [%]"] = self.mole_frac_in["CO2"] * 100
        var_dict["Feed flow rate [mol/s]"] = self.flow_mol_in_total
        var_dict["Feed velocity [m/s]"] = self.velocity_in
        var_dict["Minimum fluidization velocity [m/s]"] = self.velocity_mf

        var_dict["Time of heating step [h]"] = units.convert(
            self.heating.time, to_units=units.hr
        )
        var_dict["Time of cooling step [h]"] = units.convert(
            self.cooling.time, to_units=units.hr
        )
        var_dict["Time of pressurization step [h]"] = units.convert(
            self.pressurization.time, to_units=units.hr
        )
        var_dict["Time of adsorption step [h]"] = units.convert(
            self.adsorption.time, to_units=units.hr
        )
        var_dict["Cycle time [h]"] = self.cycle_time

        var_dict["Purity [-]"] = self.purity
        var_dict["Recovery [-]"] = self.recovery
        var_dict["Productivity [kg CO2/ton/h]"] = self.productivity
        var_dict["Specific energy [MJ/kg CO2]"] = self.specific_energy
        var_dict["Heat duty per bed [MW]"] = self.heat_duty_bed
        var_dict["Heat duty total [MW]"] = self.heat_duty_total

        var_dict["Pressure drop [Pa]"] = -self.pressure_drop
        var_dict["Number of beds"] = self.number_beds

        var_dict[
            "CO2 captured in one cycle per bed [kg/cycle]"
        ] = self.CO2_captured_bed_cycle

        var_dict["Cycles per year"] = self.cycles_year

        var_dict[
            "Total CO2 captured per year [tonne/year]"
        ] = self.total_CO2_captured_year

        var_dict[
            "Amount of flue gas processed per year [Gmol/year]"
        ] = self.flue_gas_processed_year

        var_dict[
            "Amount of flue gas processed per year (target) [Gmol/year]"
        ] = self.flue_gas_processed_year_target

        var_dict["Amount of CO2 to atmosphere [mol/s]"] = self.emissions_co2

        var_dict[
            "Concentration of CO2 emitted to atmosphere [ppm]"
        ] = self.emissions_co2_ppm

        return var_dict

    def _get_performance_contents(self, time_point=0):

        var_dict = self.get_var_dict()

        rem_from_var_dict = []
        for key, val in var_dict.items():
            if isinstance(val, Var):
                pass
            else:
                rem_from_var_dict.append(key)

        for r in rem_from_var_dict:
            del var_dict[r]

        return {"vars": var_dict}

    def _get_stream_table_contents(self, time_point=0):

        stream_attributes = {}
        stream_attributes["Units"] = {}

        for n, v in {
            "Inlet": "inlet",
            "CO2 Rich Stream": "co2_rich_stream",
            "N2 Rich Stream": "n2_rich_stream",
            "H2O Stream": "h2o_o2_stream",
        }.items():
            port_obj = getattr(self, v)

            stream_attributes[n] = {}

            for k in port_obj.vars:
                for i in port_obj.vars[k].keys():
                    if isinstance(i, float):
                        var = port_obj.vars[k][time_point]
                        stream_attributes[n][k] = value(var)
                        stream_attributes["Units"][k] = units.get_units(var)
                    else:
                        if len(i) == 2:
                            kname = str(i[1])
                        else:
                            kname = str(i[1:])
                        var = port_obj.vars[k][time_point, i[1:]]
                        stream_attributes[n][k + " " + kname] = value(var)
                        stream_attributes["Units"][k + " " + kname] = units.get_units(
                            var
                        )

        return DataFrame.from_dict(stream_attributes, orient="columns")
