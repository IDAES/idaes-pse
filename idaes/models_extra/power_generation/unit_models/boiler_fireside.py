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
Main Assumptions:

- Fireside boiler surrogate model with mass and energy balance enforced
- Flue gas compositions of major species O2, N2, CO2, H2O and SO2 are
  calculated based on elemental mass balance
- Flue gas composition of NO is calculated based on algebraic surrogate model
- Unburned organic elements are calculated based on algebraic surrogate model
- Heat duties in waterwall zones and platen SH and roof are calculated based
  on surrogate model functions
- Flue Gas Exit Temperature (FGET) is calculated based on energy balance

The model has two inlets:

* primary air inlet: PA, flow rate and composition

* secondary air inlet: SA, flow rate and composition

* primary_air_moist: primary air minus vaporized moisture (calculated)


The model has two outlets:

* Flue gas (gas only): flow rate and composition calculated by mass balance
* flyash: Fly ash with unburned organic elements (solid phase only)


Notes:

1. SR versus raw coal flow rate can be specified by a flowsheet constraint
   and should be consistent with the surrogate model.
2. PA to raw coal flow rate ratio can be specified by a flowsheet constraint
3. Fraction of coal moisture vaporized in mill should be consistent with
   the surrogate model (typically 0.6)
4. Moisture mass fraction in raw coal is a user input and should be consistent
   with the values used to train the surrogate model


Surrogate models need to be provided by the user:

* heat duty, one for each water wall zone
* 1 for platen super heater
* 1 for roof and backpass
* 1 for NOx prediction
* 1 for flyash

Note that the user can provide either surrogate models or fixed numeric values

The surrogate models are typically a function of:

* raw coal mass flow rate,
* moisture mass fraction of raw coal,
* overall stoichiometric ratio,
* lower furnace Stoichiometric ratio,
* Secondary air temperature

"""
# Import Pyomo libraries
from pyomo.common.config import ConfigBlock, ConfigValue, Bool

# Import IDAES cores
from idaes.core import declare_process_block_class, UnitModelBlockData, useDefault

from idaes.core.util.model_statistics import degrees_of_freedom

from idaes.core.util.config import is_physical_parameter_block, DefaultBool
from idaes.core.util.exceptions import ConfigurationError
import idaes.logger as idaeslog


# Additional import for the unit operation
from pyomo.environ import Var, Param, exp, RangeSet, Constraint, log
import idaes.core.util.scaling as iscale
from idaes.core.util.constants import Constants as const
from idaes.core.solvers import get_solver

__author__ = "Boiler Team (J. Ma, M. Zamarripa)"
__version__ = "1.0.0"


# ----------------------------------------------------------------------------
@declare_process_block_class("BoilerFireside")
class BoilerFiresideData(UnitModelBlockData):
    """
    Boiler fire-side surrogate model with enforced mass and energy balance
    """

    CONFIG = ConfigBlock()
    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=DefaultBool,
            default=useDefault,
            description="Dynamic model flag",
            doc="""Indicates whether this model will be dynamic or not,
**default** = useDefault.
**Valid values:** {
**useDefault** - get flag from parent (default = False),
**True** - set as dynamic model,
**False** - set as a steady-state model.}""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Holdup construction flag",
            doc="""Indicates whether holdup terms should be constructed or not.
Must be True if dynamic = True,
**default** - False.
**Valid values:** {
**True** - construct holdup terms,
**False** - do not construct holdup terms}""",
        ),
    )
    CONFIG.declare(
        "property_package",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for control volume",
            doc="""Property parameter object used to define property calculations,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PhysicalParameterObject** - a PhysicalParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "property_package_args",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing property packages",
            doc="""A ConfigBlock with arguments to be passed to a property block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}""",
        ),
    )
    CONFIG.declare(
        "calculate_PA_SA_flows",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Pressure change term construction flag",
            doc="""Indicates whether the primary air and secondary air calculations
have to be constructed or not,
**default** - False.
**Valid values:** {
**True** - calculate primary air and secondary air based on stoichiometric
ratio, PA to coal ratio, and lower stoichiometric ratio,
**False** - primary air and secondary air are inputs of the model.}""",
        ),
    )
    CONFIG.declare(
        "number_of_zones",
        ConfigValue(
            default=16,
            description="number of boiler zones",
            doc="number of boiler zones",
        ),
    )
    CONFIG.declare(
        "has_platen_superheater",
        ConfigValue(
            default=True,
            domain=Bool,
            description="True if boiler includes a platen superheater",
            doc="""Indicates if boiler includes a platen superheater,
**default** - True.
**Valid values:** {
**True** - include heat duty to platen superheater,
**False** - exclude heat duty to platen superheater.}""",
        ),
    )
    CONFIG.declare(
        "has_roof_superheater",
        ConfigValue(
            default=True,
            domain=Bool,
            description="True if roof is treated as a superheater",
            doc="""Indicates if roof is a section of superheater,
**default** - True.
**Valid values:** {
**True** - include heat duty to roof superheater,
**False** - exclude heat duty to roof superheater.}""",
        ),
    )
    CONFIG.declare(
        "surrogate_dictionary",
        ConfigValue(
            default=None,
            description="surrogate model dictionary",
            doc="""user must provide surrogate models or values for heat duty,
**default** - False.
**Valid values:** """,
        ),
    )

    def build(self):
        super(BoilerFiresideData, self).build()
        # Insert user provided custom surrogate model function
        # self.config.surrogate_function(self)
        if self.config.surrogate_dictionary is None:
            raise ConfigurationError(
                "{} - User needs to provide a dictionary of surrogates".format(
                    self.name
                )
            )

        self.primary_air = self.config.property_package.build_state_block(
            self.flowsheet().time, **self.config.property_package_args
        )
        self.primary_air_moist = self.config.property_package.build_state_block(
            self.flowsheet().time, **self.config.property_package_args
        )
        self.secondary_air = self.config.property_package.build_state_block(
            self.flowsheet().time, **self.config.property_package_args
        )
        self.flue_gas = self.config.property_package.build_state_block(
            self.flowsheet().time, **self.config.property_package_args
        )

        self.add_port("primary_air_inlet", self.primary_air)
        self.add_port("secondary_air_inlet", self.secondary_air)
        self.add_port("flue_gas_outlet", self.flue_gas)
        # Construct performance equations
        self._make_params()
        self._make_vars()
        self._make_mass_balance()
        self._make_energy_balance()
        self._make_momentum_balance()
        self._import_surrogate_models()

    def _import_surrogate_models(self):

        data_dict = self.config.surrogate_dictionary
        if len(data_dict) != (
            len(self.zones)
            + self.config.has_platen_superheater
            + self.config.has_roof_superheater
            + 2  # flyash surrogate and NOx surrogate
        ):
            raise ConfigurationError(
                "User needs to provide the same number "
                "of surrogate models and water wall zones"
                "and/or platen sh and/or boiler roof"
            )

        # Surrogate model predictions
        # Constraints for heat duty in boiler water wall zones
        @self.Constraint(
            self.flowsheet().time,
            self.zones,
            doc="Surrogate model for heat loss" " to water wall zones",
        )
        def eq_surr_waterwall_heat(b, t, z):
            return b.waterwall_heat[t, z] * b.fcorrection_heat_ww[t] == eval(
                data_dict[z]
            )

        if self.config.has_platen_superheater is True:

            @self.Constraint(
                self.flowsheet().time,
                doc="Surrogate model for heat loss" " to platen superheater",
            )
            def eq_surr_platen_heat(b, t):
                return b.platen_heat[t] * b.fcorrection_heat_platen[t] == eval(
                    data_dict["pl"]
                )

        if self.config.has_roof_superheater is True:

            @self.Constraint(
                self.flowsheet().time,
                doc="Surrogate model for heat loss in " " the roof and backpass heater",
            )
            def eq_surr_roof_heat(b, t):
                return b.roof_heat[t] * b.fcorrection_heat_ww[t] == eval(
                    data_dict["roof"]
                )

        # Constraints for unburned carbon
        @self.Constraint(
            self.flowsheet().time,
            doc="Surrogate model for" " mass fraction of unburned carbon",
        )
        def eq_surr_ln_ubc(b, t):
            return b.ubc_in_flyash[t] == eval(data_dict["flyash"])

        # Constraints for NOx in mol fraction, surrogate model in PPM,
        # converted to mass fraction
        @self.Constraint(
            self.flowsheet().time,
            doc="NOx in mol fraction" "surrogate model must be in PPM",
        )
        def eq_surr_nox(b, t):
            return b.frac_mol_NOx_fluegas[t] * 1e6 == eval(data_dict["NOx"])

        #                    # 1e6 conversion factor from PPM to mol fract

    def _make_params(self):
        """This section is for parameters within this model."""
        # atomic mass of elements involved in kg/mol
        self.atomic_mass_C = Param(initialize=0.0120107, doc="Atomic mass of C")
        self.atomic_mass_H = Param(initialize=0.00100795, doc="Atomic mass of H")
        self.atomic_mass_O = Param(initialize=0.0159994, doc="Atomic mass of O")
        self.atomic_mass_N = Param(initialize=0.0140067, doc="Atomic mass of N")
        self.atomic_mass_S = Param(initialize=0.0320652, doc="Atomic mass of S")
        # mole fractions of O2, H2O, CO2, SO2, and NO in air
        # at 25 C and relative humidity of 25% with SO2 and NO at 0.1 ppm
        # change air compostion accordingly if relative humidity changes
        # N2 mole fraction not needed since summation of mole fractions of
        # all six species is 1. It is 0.7839958
        if self.config.calculate_PA_SA_flows is True:
            self.mole_frac_air = Param(
                self.config.property_package.component_list,
                mutable=True,
                initialize={
                    "O2": 0.20784,
                    "N2": 0.783994,
                    "NO": 0.0000001,
                    "CO2": 0.0003373,
                    "H2O": 0.0078267,
                    "SO2": 0.0000001,
                },
                doc="Mole fraction of air species",
            )

    def _make_vars(self):
        """This section is for variables within this model."""
        # Number of waterwall zones is given by config,
        # Optionally the model could contain a platen and a roof superheater

        self.deltaP = Var(
            self.flowsheet().time,
            initialize=1000,
            doc="Pressure drop of secondary air " "through windbox and burner [Pa]",
        )

        self.zones = RangeSet(self.config.number_of_zones)

        self.fcorrection_heat_ww = Var(
            self.flowsheet().time,
            initialize=1,
            doc="Correction factor " "for waterwall heat duty",
        )

        if self.config.has_platen_superheater is True:
            self.fcorrection_heat_platen = Var(
                self.flowsheet().time,
                initialize=1,
                doc="Correction factor for " "platen SH heat duty",
            )

        # wall temperatures of water wall zones
        self.wall_temperature_waterwall = Var(
            self.flowsheet().time,
            self.zones,
            initialize=700.0,
            doc="Wall temperature [K] in " "Waterwall zones",
        )

        # heat duties for water wall zones
        self.waterwall_heat = Var(
            self.flowsheet().time,
            self.zones,
            initialize=2.0e7,
            doc="Heat duty [W] or heat loss " "to waterwall zones",
        )

        if self.config.has_platen_superheater is True:
            # heat duty to platen super heater
            self.platen_heat = Var(
                self.flowsheet().time,
                initialize=6.0e7,
                doc="Platen superheater heat duty [W]",
            )
            # wall temperature of platen superheater
            self.wall_temperature_platen = Var(
                self.flowsheet().time,
                initialize=800.0,
                doc="Platen superheater" " wall temperature [K]",
            )

        if self.config.has_roof_superheater is True:
            # heat duty to roof superheater
            self.roof_heat = Var(
                self.flowsheet().time,
                initialize=5e6,
                doc="Roof superheater heat duty [W]",
            )
            # wall temperature of roof superheater
            self.wall_temperature_roof = Var(
                self.flowsheet().time,
                initialize=800.0,
                doc="Roof superheater " "wall temperature [K]",
            )

        # PA/coal temperature, usually fixed around 150 F
        self.temperature_coal = Var(
            self.flowsheet().time, initialize=338.7, doc="Coal temperature in K"
        )

        # overall Stoichiometric ratio, used to calculate total
        # combustion air flowrate
        # If SR is a function of load or raw coal flow rate,
        # specify constraint in flowsheet model
        self.SR = Var(
            self.flowsheet().time,
            initialize=1.15,
            doc="Overall furnace Stoichiometric ratio - SR",
        )

        # lower furnace Stoichiometric ratio
        self.SR_lf = Var(
            self.flowsheet().time,
            initialize=1.15,
            doc="Lower furnace Stoichiometric ratio" " - SR excluding overfire air",
        )

        # PA to coal mass flow ratio,
        # typically 2.0 depending on load or mill curve
        # usually this is an input of surrogate model
        self.ratio_PA2coal = Var(
            self.flowsheet().time, initialize=2.0, doc="Primary air to coal ratio"
        )

        # mass flow rate of raw coal fed to mill
        # (raw coal flow without moisture vaporization in mill)
        self.flowrate_coal_raw = Var(
            self.flowsheet().time, initialize=25.0, doc="Raw coal mass flowrate [kg/s]"
        )

        # mass flow rate of coal to burners after moisture vaporization in mill
        self.flowrate_coal_burner = Var(
            self.flowsheet().time,
            initialize=20.0,
            doc="Mass flowrate coal to burners "
            "with moisture"
            " partially vaporized in mill [kg/s]",
        )

        # raw coal moisture mass fraction (on as received basis),
        # can change with time to represent HHV change on as received basis
        self.mf_H2O_coal_raw = Var(
            self.flowsheet().time,
            initialize=0.15,
            doc="Raw coal mass fraction of" " moisture on as received basis",
        )

        # moisture mass fraction of coal to burners after mill,
        # calculated based on fraction of moistures vaporized in mill
        self.mf_H2O_coal_burner = Var(
            self.flowsheet().time,
            initialize=0.15,
            doc="Mass fraction of moisture" " on as received basis",
        )

        # Fraction of moisture vaporized in mill,
        # set in flowsheet as a function of coal flow rate, default is 0.6
        self.frac_moisture_vaporized = Var(
            self.flowsheet().time,
            initialize=0.6,
            doc="Fraction of coal" " moisture vaporized in mill",
        )

        # Vaporized moisture mass flow rate
        @self.Expression(
            self.flowsheet().time, doc="Vaporized moisture mass flowrate [kg/s]"
        )
        def flowrate_moist_vaporized(b, t):
            return (
                b.flowrate_coal_raw[t]
                * b.mf_H2O_coal_raw[t]
                * b.frac_moisture_vaporized[t]
            )

        # Elemental composition of dry coal, assuming fixed over time
        # They are usually fixed as user inputs
        self.mf_C_coal_dry = Var(initialize=0.5, doc="Mass fraction of C on dry basis")

        self.mf_H_coal_dry = Var(initialize=0.05, doc="Mass fraction of H on dry basis")

        self.mf_O_coal_dry = Var(initialize=0.1, doc="Mass fraction of O on dry basis")

        self.mf_N_coal_dry = Var(initialize=0.1, doc="Mass fraction of N on dry basis")

        self.mf_S_coal_dry = Var(initialize=0.05, doc="Mass fraction of S on dry basis")

        self.mf_Ash_coal_dry = Var(
            initialize=0.2, doc="Mass fraction of Ash on dry basis"
        )

        # High heating value of dry coal, usally as a fixed user input
        self.hhv_coal_dry = Var(
            initialize=1e7, doc="High heating value (HHV)" "of coal on dry basis J/kg"
        )

        # Fraction of unburned carbon (actually all organic elements) in
        # flyash, predicted by surrogate model.
        # When generating the surrogate, sum up all elements in
        # flyash from fireside boiler model outputs
        self.ubc_in_flyash = Var(
            self.flowsheet().time,
            initialize=0.01,
            doc="Unburned carbon and" " other organic elements in fly ash",
        )

        # mole fraction of NO in flue gas, predicted by surrogate model
        # usually mole fraction is very close to mass fraction of NO
        self.frac_mol_NOx_fluegas = Var(
            self.flowsheet().time,
            initialize=1e-4,
            doc="Mole fraction of NOx in flue gas",
        )

        # NOx in lb/MMBTU, NO is converted to NO2 as NOx
        @self.Expression(self.flowsheet().time, doc="NOx in lb/MMBTU")
        def nox_lb_mmbtu(b, t):
            return (
                b.flue_gas_outlet.flow_mol_comp[t, "NO"]
                * 0.046
                * 2.20462
                / (
                    b.flowrate_coal_raw[t]
                    * (1 - b.mf_H2O_coal_raw[t])
                    * b.hhv_coal_dry
                    / 1.054e9
                )
            )

        # coal flow rate after mill
        @self.Constraint(
            self.flowsheet().time, doc="Coal flow rate to burners after mill"
        )
        def flowrate_coal_burner_eqn(b, t):
            return (
                b.flowrate_coal_burner[t]
                == b.flowrate_coal_raw[t] - b.flowrate_moist_vaporized[t]
            )

        # moisture mass fraction in coal after mill
        @self.Constraint(
            self.flowsheet().time, doc="Moisture mass fraction for coal after mill"
        )
        def mf_H2O_coal_burner_eqn(b, t):
            return (
                b.flowrate_coal_raw[t] * b.mf_H2O_coal_raw[t]
                == b.flowrate_coal_burner[t] * b.mf_H2O_coal_burner[t]
                + b.flowrate_moist_vaporized[t]
            )

        # fraction of daf elements on dry basis
        @self.Expression(doc="Mass fraction of dry ash free (daf)")
        def mf_daf_dry(b):
            return 1 - b.mf_Ash_coal_dry

        @self.Expression(self.flowsheet().time, doc="Ash mass flow rate kg/s")
        def flowrate_ash(b, t):
            return (
                b.flowrate_coal_raw[t] * (1 - b.mf_H2O_coal_raw[t]) * b.mf_Ash_coal_dry
            )

        @self.Expression(
            self.flowsheet().time,
            doc="Dry ash free - daf_coal flow rate " "in fuel fed to the boiler kg/s",
        )
        def flowrate_daf_fuel(b, t):
            return b.flowrate_coal_raw[t] * (1 - b.mf_H2O_coal_raw[t]) * b.mf_daf_dry

        @self.Expression(
            self.flowsheet().time, doc="Dry ash free (daf) coal flow rate in fly ash"
        )
        def flowrate_daf_flyash(b, t):
            return b.flowrate_ash[t] / (1.0 - b.ubc_in_flyash[t]) * b.ubc_in_flyash[t]

        @self.Expression(self.flowsheet().time, doc="Burned daf coal flow rate")
        def flowrate_daf_burned(b, t):
            return b.flowrate_daf_fuel[t] - b.flowrate_daf_flyash[t]

        @self.Expression(self.flowsheet().time, doc="Coal burnout")
        def coal_burnout(b, t):
            return b.flowrate_daf_burned[t] / b.flowrate_daf_fuel[t]

        @self.Expression(doc="Dry ash free - daf C mass fraction")
        def mf_C_daf(b):
            return b.mf_C_coal_dry / b.mf_daf_dry

        @self.Expression(doc="Dry ash free - daf H mass fraction")
        def mf_H_daf(b):
            return b.mf_H_coal_dry / b.mf_daf_dry

        @self.Expression(doc="Dry ash free - daf O mass fraction")
        def mf_O_daf(b):
            return b.mf_O_coal_dry / b.mf_daf_dry

        @self.Expression(doc="Dry ash free - daf N mass fraction")
        def mf_N_daf(b):
            return b.mf_N_coal_dry / b.mf_daf_dry

        @self.Expression(doc="Dry ash free - daf S mass fraction")
        def mf_S_daf(b):
            return b.mf_S_coal_dry / b.mf_daf_dry

        @self.Expression(doc="Dry ash free high heating value in J/kg")
        def hhv_daf(b):
            return b.hhv_coal_dry / b.mf_daf_dry

        @self.Expression(
            doc="Heat of combustion at constant " "pressure on daf basis in J/kg"
        )
        # Note that measure HHV is the heat of combustion at constant volume
        def dhcoal(b):
            return -b.hhv_daf + const.gas_constant * 298.15 / 2 * (
                -b.mf_H_daf / 2 / b.atomic_mass_H
                + b.mf_O_daf / b.atomic_mass_O
                + b.mf_N_daf / b.atomic_mass_N
            )

        @self.Expression(doc="Carbon heat of reaction - dhc J/kg")
        def dhc(b):
            return -94052 * 4.184 / b.atomic_mass_C * b.mf_C_daf

        @self.Expression(doc="Hydrogen heat of reaction - dhh in J/kg")
        def dhh(b):
            return -68317.4 * 4.184 / b.atomic_mass_H / 2 * b.mf_H_daf

        @self.Expression(doc="Sulfur heat of reaction - dhs in J/kg")
        def dhs(b):
            return -70940 * 4.184 / b.atomic_mass_S * b.mf_S_daf

        @self.Expression(doc="Heat of formation for daf coal in J/kg")
        def hf_daf(b):
            return b.dhc + b.dhh + b.dhs - b.dhcoal

        @self.Expression(
            self.flowsheet().time,
            doc="Heat of formation of " "moisture-containing coal to burners",
        )
        def hf_coal(b, t):
            return b.hf_daf * (
                1 - b.mf_H2O_coal_burner[t]
            ) * b.mf_daf_dry + b.mf_H2O_coal_burner[t] * (-68317.4) * 4.184 / (
                b.atomic_mass_H * 2 + b.atomic_mass_O
            )

        @self.Expression(doc="Average atomic mass of daf coal")
        def am_daf(b):
            return 1 / (
                b.mf_C_daf / b.atomic_mass_C
                + b.mf_H_daf / b.atomic_mass_H
                + b.mf_O_daf / b.atomic_mass_O
                + b.mf_N_daf / b.atomic_mass_N
                + b.mf_S_daf / b.atomic_mass_S
            )

        @self.Expression(
            self.flowsheet().time,
            doc="Gt1 or Einstein quantum theory function for daf "
            "coal sensible heat calculation",
        )
        def gt1(b, t):
            return 1 / (exp(380 / b.temperature_coal[t]) - 1)

        # since flyash and flue gas outlet temperature is not fixed
        self.gt1_flyash = Var(
            self.flowsheet().time,
            initialize=1,
            doc="Gt1 or Einstein quantum theory function" " for daf part of flyash",
        )

        @self.Constraint(
            self.flowsheet().time,
            doc="Gt1 or Einstein quantum theory " "for daf part of flyash",
        )
        def gt1_flyash_eqn(b, t):
            return b.gt1_flyash[t] * (exp(380 / b.flue_gas[t].temperature) - 1) == 1

        @self.Expression(
            self.flowsheet().time,
            doc="Gt2 or Einstein quantum theory "
            "function for daf coal "
            "sensible heat calculation for high temperature",
        )
        def gt2(b, t):
            return 1 / (exp(1800 / b.temperature_coal[t]) - 1)

        # declare variable rather than use expression
        # since flyash and flue gas outlet temperature is not fixed
        self.gt2_flyash = Var(
            self.flowsheet().time,
            initialize=1,
            doc="Gt2 or Einstein quantum theory function " "for daf part of flyash",
        )

        @self.Constraint(
            self.flowsheet().time,
            doc="Gt2 or Einstein quantum theory function " "for daf part of flyash",
        )
        def gt2_flyash_eqn(b, t):
            return b.gt2_flyash[t] * (exp(1800 / b.flue_gas[t].temperature) - 1) == 1

        @self.Expression(self.flowsheet().time, doc="Sensible heat of daf coal")
        def hs_daf(b, t):
            return (
                const.gas_constant
                / b.am_daf
                * (380 * (b.gt1[t] - 0.3880471566) + 3600 * (b.gt2[t] - 0.002393883))
            )

        @self.Expression(
            self.flowsheet().time, doc="Sensible heat for daf part of flyash"
        )
        def hs_daf_flyash(b, t):
            return (
                const.gas_constant
                / b.am_daf
                * (
                    380 * (b.gt1_flyash[t] - 0.3880471566)
                    + 3600 * (b.gt2_flyash[t] - 0.002393883)
                )
            )

        @self.Expression(self.flowsheet().time, doc="Sensible heat of coal to burners")
        def hs_coal(b, t):
            return (
                (1 - b.mf_H2O_coal_burner[t]) * b.mf_daf_dry * b.hs_daf[t]
                + b.mf_H2O_coal_burner[t] * 4184 * (b.temperature_coal[t] - 298.15)
                + (1 - b.mf_H2O_coal_burner[t])
                * b.mf_Ash_coal_dry
                * (
                    593 * (b.temperature_coal[t] - 298.15)
                    + 0.293 * (b.temperature_coal[t] ** 2 - 298.15**2)
                )
            )

        @self.Expression(
            self.flowsheet().time,
            doc="Total enthalpy of " "moisture-containing coal to burners",
        )
        def h_coal(b, t):
            return b.hs_coal[t] + b.hf_coal[t]

        # variable for O2 mole fraction in flue gas on dry basis
        # it can be used for data reconciliation and parameter estimation
        self.fluegas_o2_pct_dry = Var(
            self.flowsheet().time, initialize=3, doc="Mol percent of O2 on dry basis"
        )

    def _make_mass_balance(self):
        # PA flow rate calculated based on PA/coal ratio
        @self.Constraint(
            self.flowsheet().time, doc="Ratio of PA mass flow to raw coal mass flow"
        )
        def ratio_PA2coal_eqn(b, t):
            return (
                b.primary_air[t].flow_mass
                == b.flowrate_coal_raw[t] * b.ratio_PA2coal[t]
            )

        @self.Expression(self.flowsheet().time, doc="C molar flow from coal")
        def molflow_C_fuel(b, t):
            return b.flowrate_daf_fuel[t] * b.mf_C_daf / b.atomic_mass_C

        @self.Expression(self.flowsheet().time, doc="H molar flow from coal")
        def molflow_H_fuel(b, t):
            return b.flowrate_daf_fuel[t] * b.mf_H_daf / b.atomic_mass_H

        @self.Expression(self.flowsheet().time, doc="O molar flow from coal")
        def molflow_O_fuel(b, t):
            return b.flowrate_daf_fuel[t] * b.mf_O_daf / b.atomic_mass_O

        @self.Expression(self.flowsheet().time, doc="N molar flow from coal")
        def molflow_N_fuel(b, t):
            return b.flowrate_daf_fuel[t] * b.mf_N_daf / b.atomic_mass_N

        @self.Expression(self.flowsheet().time, doc="S molar flow from coal")
        def molflow_S_fuel(b, t):
            return b.flowrate_daf_fuel[t] * b.mf_S_daf / b.atomic_mass_S

        if self.config.calculate_PA_SA_flows is True:
            # Calculate PA component molar flow based on air mol fractions
            # specified as model parameter
            # Overwrite the composition from inlets
            # Let N2 molar flow calculated by balance
            @self.Constraint(
                self.flowsheet().time,
                self.config.property_package.component_list,
                doc="PA component molar flow",
            )
            def molar_flow_PA_eqn(b, t, j):
                if j == "N2":
                    return Constraint.Skip
                else:
                    return (
                        b.primary_air[t].flow_mol_comp[j]
                        == b.primary_air[t].flow_mol * b.mole_frac_air[j]
                    )

            # Calculate SA component molar flow based on air mol fractions
            # specified as model parameter
            # Overwrite the composition from inlets
            # Let N2 molar flow calculated by balance
            @self.Constraint(
                self.flowsheet().time,
                self.config.property_package.component_list,
                doc="SA component molar flow",
            )
            def molar_flow_SA_eqn(b, t, j):
                if j == "N2":
                    return Constraint.Skip
                else:
                    return (
                        b.secondary_air[t].flow_mol_comp[j]
                        == b.secondary_air[t].flow_mol * b.mole_frac_air[j]
                    )

        # calculate molar flow of primary_air_moist stream
        @self.Constraint(
            self.flowsheet().time, self.config.property_package.component_list
        )
        def primary_air_moist_comp_flow_eqn(b, t, j):
            if j == "H2O":
                return (
                    b.primary_air[t].flow_mol_comp["H2O"]
                    + b.flowrate_moist_vaporized[t]
                    / (b.atomic_mass_H * 2 + b.atomic_mass_O)
                ) == b.primary_air_moist[t].flow_mol_comp["H2O"]
            else:
                return (
                    b.primary_air[t].flow_mol_comp[j]
                    == b.primary_air_moist[t].flow_mol_comp[j]
                )

        # total combustion air mass flow rate
        @self.Expression(
            self.flowsheet().time, doc="Total combustion air mass flow rate"
        )
        def flow_mass_TCA(b, t):
            return b.primary_air[t].flow_mass + b.secondary_air[t].flow_mass

        # overall Stoichiometric ratio (SR)
        @self.Constraint(self.flowsheet().time, doc="SR equation")
        def SR_eqn(b, t):
            return (
                b.SR[t]
                * (
                    b.molflow_C_fuel[t]
                    + b.molflow_H_fuel[t] / 4
                    + b.molflow_S_fuel[t]
                    - b.molflow_O_fuel[t] / 2
                )
                == b.primary_air[t].flow_mol_comp["O2"]
                + b.secondary_air[t].flow_mol_comp["O2"]
            )

        @self.Expression(self.flowsheet().time, doc="C mole flow")
        def molflow_C_fluegas(b, t):
            return (
                b.flowrate_daf_burned[t] * b.mf_C_daf / b.atomic_mass_C
                + b.primary_air_moist[t].flow_mol_comp["CO2"]
                + b.secondary_air[t].flow_mol_comp["CO2"]
            )

        @self.Expression(self.flowsheet().time, doc="H mole flow")
        def molflow_H_fluegas(b, t):
            return (
                b.flowrate_daf_burned[t] * b.mf_H_daf / b.atomic_mass_H
                + b.primary_air_moist[t].flow_mol_comp["H2O"] * 2
                + b.secondary_air[t].flow_mol_comp["H2O"] * 2
                + b.flowrate_coal_burner[t]
                * b.mf_H2O_coal_burner[t]
                / (b.atomic_mass_H * 2 + b.atomic_mass_O)
                * 2
            )

        @self.Expression(self.flowsheet().time, doc="O mole flow")
        def molflow_O_fluegas(b, t):
            return (
                b.flowrate_daf_burned[t] * b.mf_O_daf / b.atomic_mass_O
                + b.primary_air_moist[t].flow_mol_comp["O2"] * 2
                + b.primary_air_moist[t].flow_mol_comp["CO2"] * 2
                + b.primary_air_moist[t].flow_mol_comp["H2O"]
                + b.primary_air_moist[t].flow_mol_comp["SO2"] * 2
                + b.primary_air_moist[t].flow_mol_comp["NO"]
                + b.secondary_air[t].flow_mol_comp["O2"] * 2
                + b.secondary_air[t].flow_mol_comp["CO2"] * 2
                + b.secondary_air[t].flow_mol_comp["H2O"]
                + b.secondary_air[t].flow_mol_comp["SO2"] * 2
                + b.secondary_air[t].flow_mol_comp["NO"]
                + b.flowrate_coal_burner[t]
                * b.mf_H2O_coal_burner[t]
                / (b.atomic_mass_H * 2 + b.atomic_mass_O)
            )

        @self.Expression(self.flowsheet().time, doc="N mole flow")
        def molflow_N_fluegas(b, t):
            return (
                b.flowrate_daf_burned[t] * b.mf_N_daf / b.atomic_mass_N
                + b.primary_air_moist[t].flow_mol_comp["N2"] * 2
                + b.primary_air_moist[t].flow_mol_comp["NO"]
                + b.secondary_air[t].flow_mol_comp["N2"] * 2
                + b.secondary_air[t].flow_mol_comp["NO"]
            )

        @self.Expression(self.flowsheet().time, doc="S mole flow")
        def molflow_S_fluegas(b, t):
            return (
                b.flowrate_daf_burned[t] * b.mf_S_daf / b.atomic_mass_S
                + b.primary_air_moist[t].flow_mol_comp["SO2"]
                + b.secondary_air[t].flow_mol_comp["SO2"]
            )

        # calculate flue gas flow component at flue_gas_outlet
        # NO mole fraction in flue gas is given by NOx surrogate model
        @self.Constraint(self.flowsheet().time, doc="NO at flue gas outlet")
        def NO_eqn(b, t):
            return (
                b.flue_gas[t].flow_mol_comp["NO"]
                == b.frac_mol_NOx_fluegas[t] * b.flue_gas[t].flow_mol
            )

        # N2 at outlet
        @self.Constraint(self.flowsheet().time, doc="N2 at outlet")
        def N2_eqn(b, t):
            return (
                b.flue_gas_outlet.flow_mol_comp[t, "N2"]
                == +b.molflow_N_fluegas[t] / 2
                - b.flue_gas_outlet.flow_mol_comp[t, "NO"] / 2
            )

        # SO2 at outlet
        @self.Constraint(self.flowsheet().time, doc="SO2 at outlet")
        def SO2_eqn(b, t):
            return b.flue_gas_outlet.flow_mol_comp[t, "SO2"] == b.molflow_S_fluegas[t]

        # H2O at outlet
        @self.Constraint(self.flowsheet().time, doc="H2O at outlet")
        def H2O_eqn(b, t):
            return (
                b.flue_gas_outlet.flow_mol_comp[t, "H2O"] == b.molflow_H_fluegas[t] / 2
            )

        # CO2 at outlet
        @self.Constraint(self.flowsheet().time, doc="CO2 at outlet")
        def CO2_eqn(b, t):
            return b.flue_gas_outlet.flow_mol_comp[t, "CO2"] == b.molflow_C_fluegas[t]

        # O2 at outlet
        @self.Constraint(self.flowsheet().time, doc="O2 at outlet")
        def O2_eqn(b, t):
            return (
                b.flue_gas_outlet.flow_mol_comp[t, "O2"]
                == (
                    b.molflow_O_fluegas[t]
                    - b.molflow_C_fluegas[t] * 2
                    - b.molflow_H_fluegas[t] / 2
                    - b.molflow_S_fluegas[t] * 2
                    - b.flue_gas_outlet.flow_mol_comp[t, "NO"]
                )
                / 2
            )

        # constraint for mol percent of O2 in flue gas on dry basis
        @self.Constraint(
            self.flowsheet().time, doc="Mol percent of O2 in flue gas on dry basis"
        )
        def fluegas_o2_pct_dry_eqn(b, t):
            return (
                b.fluegas_o2_pct_dry[t]
                * (
                    b.flue_gas_outlet.flow_mol_comp[t, "O2"]
                    + b.flue_gas_outlet.flow_mol_comp[t, "N2"]
                    + b.flue_gas_outlet.flow_mol_comp[t, "CO2"]
                    + b.flue_gas_outlet.flow_mol_comp[t, "SO2"]
                    + b.flue_gas_outlet.flow_mol_comp[t, "NO"]
                )
                / 100
                == b.flue_gas_outlet.flow_mol_comp[t, "O2"]
            )

        # mass flow rate of flyash containing unburned fuel
        @self.Expression(self.flowsheet().time, doc="Flyash mass flow rate  kg/s")
        def flow_mass_flyash(b, t):
            return b.flowrate_daf_flyash[t] + b.flowrate_ash[t]

    def _make_momentum_balance(self):
        # flue gas pressure is secondary air presure
        # - pressure drop through windbox
        # and burner secondary air register
        @self.Constraint(self.flowsheet().time, doc="Flue gas pressure in Pascals")
        def flue_gas_pressure_eqn(b, t):
            return (
                b.flue_gas[t].pressure * 1e-5
                == (b.secondary_air[t].pressure - b.deltaP[t]) * 1e-5
            )

        # set pressure of primary air with moisture at mill outlet
        @self.Constraint(self.flowsheet().time, doc="Mill outlet pressure")
        def primary_air_moist_pressure_eqn(b, t):
            return (
                b.primary_air[t].pressure * 1e-5
                == b.primary_air_moist[t].pressure * 1e-5
            )

    def _make_energy_balance(self):
        # temperature of primary_air_moist is equal to coal temparature
        # leaving mill, ignore the temperature of primary_air_inlet in energy
        # balance equation
        @self.Constraint(
            self.flowsheet().time, doc="Temperature of primary air leaving mill"
        )
        def primary_air_moist_temperature_eqn(b, t):
            return b.primary_air_moist[t].temperature == b.temperature_coal[t]

        # overall energy balance to calculate FEGT
        @self.Constraint(self.flowsheet().time, doc="Temperature at flue gas outlet")
        def flue_gas_temp_eqn(b, t):
            return (
                b.primary_air_moist[t].flow_mol * b.primary_air_moist[t].enth_mol
                + b.secondary_air[t].flow_mol * b.secondary_air[t].enth_mol
                + b.h_coal[t] * b.flowrate_coal_burner[t]
            ) == (
                sum(b.waterwall_heat[t, j] for j in b.zones)
                + (
                    b.platen_heat[t]
                    if self.config.has_platen_superheater is True
                    else 0
                )
                + (b.roof_heat[t] if self.config.has_roof_superheater is True else 0)
                + b.flue_gas[t].flow_mol * b.flue_gas[t].enth_mol
                # sensible heat of ash
                + b.flowrate_ash[t]
                * (
                    593 * (b.flue_gas[t].temperature - 298.15)
                    + 0.293 * (b.flue_gas[t].temperature ** 2.0 - 298.15**2.0)
                )
                # sensible heat of unburned fuel
                + b.flowrate_daf_flyash[t] * b.hs_daf_flyash[t]
            )

        # expression to calculate total heat duty of waterwall
        @self.Expression(
            self.flowsheet().time, doc="Total heat duty of all waterwall zones"
        )
        def heat_total_ww(b, t):
            return sum(b.waterwall_heat[t, j] for j in b.zones)

        # expression to calculate total heat loss through
        # waterwall, platen SH, and roof SH
        @self.Expression(self.flowsheet().time, doc="Total heat duty")
        def heat_total(b, t):
            return (
                b.heat_total_ww[t]
                + (
                    b.platen_heat[t]
                    if self.config.has_platen_superheater is True
                    else 0
                )
                + (b.roof_heat[t] if self.config.has_roof_superheater is True else 0)
            )

    def initialize_build(
        blk,
        state_args_PA=None,
        state_args_SA=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """
        Initialization routine.
        1.- initialize state blocks, using an initial guess for inlet
        primary air and secondary air.
        2.- Use PA and SA values to guess flue gas component molar flowrates,
        Temperature, and Pressure. Initialize flue gas state block.
        3.- Then, solve complete model.

        Keyword Arguments:
            state_args_PA : a dict of arguments to be passed to the property
                           package(s) for the primary air state block to
                           provide an initial state for initialization
                           (see documentation of the specific property package)
                           (default = None).
            state_args_SA : a dict of arguments to be passed to the property
                           package(s) for secondary air state block to
                           provide an initial state for initialization
                           (see documentation of the specific property package)
                           (default = None).
            outlvl : sets output level of initialisation routine
            optarg : solver options dictionary object (default=None, use
                     default solver options)
            solver : str indicating which solver to use during
                     initialization (default = None, use default solver)

        Returns:
            None
        """
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")

        # Create solver
        opt = get_solver(solver, optarg)

        # ---------------------------------------------------------------------
        # Initialize inlet property blocks
        blk.primary_air.initialize(
            outlvl=outlvl, optarg=optarg, solver=solver, state_args=state_args_PA
        )
        blk.primary_air_moist.initialize(
            outlvl=outlvl, optarg=optarg, solver=solver, state_args=state_args_PA
        )
        blk.secondary_air.initialize(
            outlvl=outlvl, optarg=optarg, solver=solver, state_args=state_args_SA
        )
        init_log.info_high("Initialization Step 1 Complete.")

        state_args = {
            "flow_mol_comp": {
                "H2O": (
                    blk.primary_air_inlet.flow_mol_comp[0, "H2O"].value
                    + blk.secondary_air_inlet.flow_mol_comp[0, "H2O"].value
                ),
                "CO2": (
                    blk.primary_air_inlet.flow_mol_comp[0, "CO2"].value
                    + blk.secondary_air_inlet.flow_mol_comp[0, "CO2"].value
                ),
                "N2": (
                    blk.primary_air_inlet.flow_mol_comp[0, "N2"].value
                    + blk.secondary_air_inlet.flow_mol_comp[0, "N2"].value
                ),
                "O2": (
                    blk.primary_air_inlet.flow_mol_comp[0, "O2"].value
                    + blk.secondary_air_inlet.flow_mol_comp[0, "O2"].value
                ),
                "SO2": (
                    blk.primary_air_inlet.flow_mol_comp[0, "SO2"].value
                    + blk.secondary_air_inlet.flow_mol_comp[0, "SO2"].value
                ),
                "NO": (
                    blk.primary_air_inlet.flow_mol_comp[0, "NO"].value
                    + blk.secondary_air_inlet.flow_mol_comp[0, "NO"].value
                ),
            },
            "temperature": 1350.00,
            "pressure": blk.primary_air_inlet.pressure[0].value,
        }
        # initialize flue gas outlet
        blk.flue_gas.initialize(state_args=state_args, outlvl=outlvl, solver=solver)
        init_log.info_high("Initialization Step 2 Complete.")

        if blk.config.calculate_PA_SA_flows is False:
            # Option 1: given PA and SA component flow rates - fixed inlets
            # fix inlet component molar flow rates
            # unfix ratio_PA2coal, SR, and fluegas_o2_pct_dry
            blk.primary_air_inlet.flow_mol_comp[...].fix()
            blk.secondary_air_inlet.flow_mol_comp[...].fix()
            blk.ratio_PA2coal.unfix()
            blk.SR.unfix()
            blk.fluegas_o2_pct_dry.unfix()
            dof = degrees_of_freedom(blk)

        else:
            # Option 2: SR, ratioPA2_coal to estimate TCA, PA, SA
            # unfix component molar flow rates, but keep T and P fixed.
            blk.primary_air_inlet.flow_mol_comp[:, :].unfix()
            blk.secondary_air_inlet.flow_mol_comp[:, :].unfix()
            dof = degrees_of_freedom(blk)

        if not dof == 0:
            raise ConfigurationError("User needs to check " "degrees of freedom")

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high("Initialization Step 3 {}.".format(idaeslog.condition(res)))
        init_log.info("Initialization Complete.")

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        # set a default waterwall zone heat scaling factor
        for v in self.waterwall_heat.values():
            if iscale.get_scaling_factor(v, warning=True) is None:
                iscale.set_scaling_factor(v, 1e-7)

        # set a default platen heat scaling factor
        if self.config.has_platen_superheater is True:
            for v in self.platen_heat.values():
                if iscale.get_scaling_factor(v, warning=True) is None:
                    iscale.set_scaling_factor(v, 1e-7)

        # set a default roof heat scaling factor
        if self.config.has_roof_superheater is True:
            for v in self.roof_heat.values():
                if iscale.get_scaling_factor(v, warning=True) is None:
                    iscale.set_scaling_factor(v, 1e-6)

        # set waterwall heat constraint scaling factor
        for t, c in self.eq_surr_waterwall_heat.items():
            sf = iscale.get_scaling_factor(
                self.waterwall_heat[t], default=1e-7, warning=True
            )
            iscale.constraint_scaling_transform(c, sf, overwrite=False)

        # set platen heat constraint scaling factor
        if self.config.has_platen_superheater is True:
            for t, c in self.eq_surr_platen_heat.items():
                sf = iscale.get_scaling_factor(
                    self.platen_heat[t], default=1e-7, warning=True
                )
                iscale.constraint_scaling_transform(c, sf, overwrite=False)

        # set roof heat constraint scaling factor
        if self.config.has_roof_superheater is True:
            for t, c in self.eq_surr_roof_heat.items():
                sf = iscale.get_scaling_factor(
                    self.roof_heat[t], default=1e-6, warning=True
                )
                iscale.constraint_scaling_transform(c, sf, overwrite=False)

        # set flue gas temperature constraint scaling factor
        for t, c in self.flue_gas_temp_eqn.items():
            sf = iscale.get_scaling_factor(
                self.platen_heat[t], default=1e-7, warning=True
            )
            iscale.constraint_scaling_transform(c, sf, overwrite=False)
