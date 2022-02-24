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
Costing package based on methods from:

    Process and Product Design Principles: Synthesis, Analysis, and
    Evaluation
    Seider, Seader, Lewin, Windagdo, 3rd Ed. John Wiley and Sons
    Chapter 22. Cost Accounting and Capital Cost Estimation
    22.2 Cost Indexes and Capital Investment

Curently, this costing package only includes methods for capital costing of
unit operations.
"""
from enum import Enum
import pyomo.environ as pyo

from idaes.generic_models.unit_models import (
    CSTR,
    Flash,
    Heater,
    HeatExchanger,
    HeatExchanger1D,
    HeatExchangerNTU,
    PFR,
    Pump,
    StoichiometricReactor)
from idaes.generic_models.unit_models.heat_exchanger \
    import HeatExchangerFlowPattern
from idaes.core.util.misc import register_units_of_measurement
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.constants import Constants
from idaes.core.util.math import smooth_max

from idaes.generic_models.costing.costing_base import CostingPackageBase

# Some more information about this module
__author__ = "Miguel Zamarripa, Andrew Lee"


class HXType(str, Enum):
    floating_head = 'floating_head'
    fixed_head = 'fixed_head'
    Utube = 'U-tube'
    kettle_vap = 'Kettle_vap'

    def __str__(self):
        return self.value


class HXMaterial(str, Enum):
    CS_CS = 'carbon steel/carbon steel'
    CS_Brass = 'carbon steel/brass'
    CS_SS = 'carbon steel/stainless steel'
    CS_Monel = 'carbon steel/monel'
    CS_Ti = 'carbon steel/titanium'
    CS_CrMo = 'carbon steel/Cr-Mo steel'
    CrMo_CrMo = 'Cr-Mo steel/Cr-Mo steel'
    SS_SS = 'stainless steel/stainless steel'
    Monel_Monel = 'monel/monel'
    Ti_Ti = 'titanium/titanium'

    def __str__(self):
        return self.value


class HXTubeLength(str, Enum):
    EightFoot = '8ft'
    TwelveFoot = '12ft'
    SixteenFoot = '16ft'
    TwentyFoot = '20ft'

    def __str__(self):
        return self.value


class VesselMaterial(str, Enum):
    CS = "carbon_steel"
    LowAlloy = "low_alloy_steel"
    SS304 = "stain_steel_304"
    SS316 = "stain_steel_316"
    Carpenter = "carpenter_20CB-3"
    Nickel200 = "nickel_200"
    Monel400 = "monel_400"
    Inconel600 = "inconel_600"
    Incoloy825 = "incoloy_825"
    Titanium = "titanium"

    def __str__(self):
        return self.value


class TrayType(str, Enum):
    Sieve = "sieve"
    Valve = "valve"
    BubbleCap = "bubble_cap"

    def __str__(self):
        return self.value


class TrayMaterial(str, Enum):
    CS = "carbon_steel"
    SS303 = "stain_steel_303"
    SS316 = "stain_steel_316"
    Carpenter = "carpenter_20CB-3"
    Monel = "monel"

    def __str__(self):
        return self.value


class HeaterMaterial(str, Enum):
    CS = "carbon_steel"
    CrMo = "Cr-Mo_alloy"
    SS = "stainless_steel"

    def __str__(self):
        return self.value


class HeaterSource(str, Enum):
    fuel = 'fuel'
    reformer = "reformer"
    pyrolysis = "pyrolysis"
    hotWater = "hot_water"
    salts = "salts"
    DowthermA = "dowtherm_a"
    steamBoiler = "steadm_boiler"

    def __str__(self):
        return self.value


class PumpMaterial(str, Enum):
    castIron = 'cast_iron'
    ductileIron = 'ductile_iron'
    castSteel = 'cast_steel'
    bronze = 'bronze'
    SS = 'stain_steel'
    HastelloyC = 'hastelloy_c'
    Monel = 'monel'
    Nickel = 'nickel'
    Titanium = 'titanium'
    NiAlBronze = 'Ni_Al_Bronze'
    CS = 'carbon_steel'

    def __str__(self):
        return self.value


class PumpType(str, Enum):
    centrifugal = 'centrifugal'
    externalGear = 'external_gear'
    reciprocating = 'reciprocating'

    def __str__(self):
        return self.value


class PumpMotorType(str, Enum):
    Open = 'open'
    Enclosed = 'enclosed'
    ExplosionProof = 'explosion_proof'

    def __str__(self):
        return self.value


# TODO : Mapping to unit models
class SSLWCosting(CostingPackageBase):

    # Register currency and conversion rates based on CE Index
    register_units_of_measurement("USD_500", "[currency]")  # base USD @ CEI 500
    register_units_of_measurement("USD2010", "500/550.8 * USD_500")
    register_units_of_measurement("USD2011", "500/585.7 * USD_500")
    register_units_of_measurement("USD2012", "500/584.6 * USD_500")
    register_units_of_measurement("USD2013", "500/567.3 * USD_500")
    register_units_of_measurement("USD2014", "500/576.1* USD_500")
    register_units_of_measurement("USD2015", "500/556.8 * USD_500")
    register_units_of_measurement("USD2016", "500/541.7 * USD_500")
    register_units_of_measurement("USD2017", "500/567.5 * USD_500")
    register_units_of_measurement("USD2018", "500/671.1 * USD_500")
    register_units_of_measurement("USD2019", "500/680.0 * USD_500")
    register_units_of_measurement("USD_394", "500/394 * USD_500")  # required for pump costing

    # Set the base year for all costs
    base_currency = pyo.units.USD2018
    # Set a base period for all operating costs
    base_period = pyo.units.year

    # TODO: Define any default flow costs of interest

    @staticmethod
    def build_global_params(blk):
        """
        This is where we can declare any global parameters we need, such as
        Lang factors, or coefficients for costing methods that should be
        shared across the process.

        You can do what you want here, so you could have e.g. sub-Blocks
        for each costing method to separate the parameters for each method.
        """
        pass

    @staticmethod
    def build_process_costs(blk):
        """
        This is where you do all your process wide costing.
        This is completely up to you, but you will have access to the
        following aggregate costs:

            1. blk.aggregate_capital_cost
            2. blk.aggregate_fixed_operating_cost
            3. blk.aggregate_variable_operating_cost
            4. blk.aggregate_flow_costs (indexed by flow type)
        """
        # TODO: Do we have any process level methods to add here?
        pass

    @staticmethod
    def initialize(blk):
        """
        Here we can add intialization steps for the things we built in
        build_process_costs.

        Note that the aggregate costs will be initialized by the framework.
        """
        # TODO: For now,  no additional process level costs to initialize
        pass

    def cost_heat_exchanger(blk,
                            hx_type=HXType.Utube,
                            material_type=HXMaterial.SS_SS,
                            tube_length=HXTubeLength.TwelveFoot,
                            integer=True):
        """
        Heat exchanger costing method.

        This method computes the purchase cost (CP) for a shell and tube heat
        exchanger (Eq. 22.43), the model computes the base cost (CB for 4 types
        of heat exchangers, such as floating head, fixed head, U-tube, and
        Kettle vaporizer), construction material factor (mat_factor),
        pressure design factor (pressure_factor), and
        tube length correction factor (length_factor),
        using Chemical Engineering base cost index of 500.

        Purchase Cost = pressure_factor *
                        material_factor *
                        length_factor *
                        base_cost_per_unit *
                        number_of_units

        Args:
            hx_type - HXType Enum indicating type of heat exchanger design,
                      default = HXType.Utube.
            material_type - HXMaterial Enum indicating material of
                            construction, default = HXMaterial.SS_SS.
            tube_length - HXTubeLength Enum indicating length of HX tubes,
                          default = HXTubeLength.TwelveFoot.
            integer - whether the number of units should be constrained to be
                      an integer or not (default = True).
        """
        # Validate arguments
        if hx_type not in HXType:
            raise ConfigurationError(
                f"{blk.unit_model.name} received invalid argument for hx_type:"
                f" {hx_type}. Argument must be a member of the HXType Enum.")
        if material_type not in HXMaterial:
            raise ConfigurationError(
                f"{blk.unit_model.name} received invalid argument for "
                f"material_type: {material_type}. Argument must be a member "
                "of the HXMaterial Enum.")
        if tube_length not in HXTubeLength:
            raise ConfigurationError(
                f"{blk.unit_model.name} received invalid argument for "
                f"tube_length: {tube_length}. Argument must be a member "
                "of the HXTubeLength Enum.")

        # Build generic costing variables
        _make_common_vars(blk, integer)

        # Length correction factor
        c_fl = {'8ft': 1.25, '12ft': 1.12, '16ft': 1.05, '20ft': 1.00}
        blk.length_factor = pyo.Param(
            mutable=True,
            initialize=c_fl[tube_length],
            doc='HX tube length correction factor')

        blk.hx_oversize = pyo.Param(mutable=True,
                                    initialize=1.1,
                                    doc='HX oversize factor (1.1 to 1.5)',
                                    units=pyo.units.ft**-2)

        # --------------------------------------------------
        # Base cost calculation - based on selected heat exchanger type:
        alpha_dict = {
            HXType.floating_head: {1: 11.9052, 2: 0.8709, 3: 0.09005},
            HXType.fixed_head: {1: 11.2927, 2: 0.8228, 3: 0.09861},
            HXType.Utube: {1: 11.3852, 2: 0.9186, 3: 0.09790},
            HXType.kettle_vap: {1: 12.2052, 2: 0.8709, 3: 0.09005}}
        alpha = alpha_dict[hx_type]

        # Convert area to square feet
        area_unit = (pyo.units.convert(blk.unit_model.area,
                                       to_units=pyo.units.ft**2) /
                     blk.number_of_units)

        def hx_cost_rule(blk):
            return (blk.base_cost_per_unit ==
                    pyo.exp(alpha[1] -
                            alpha[2]*pyo.log(area_unit*blk.hx_oversize) +
                            alpha[3]*pyo.log(area_unit*blk.hx_oversize)**2) *
                    pyo.units.USD_500)
        blk.base_cost_per_unit_eq = pyo.Constraint(rule=hx_cost_rule)

        @blk.Expression(doc="Base cost for all installed units")
        def base_cost(blk):
            return blk.base_cost_per_unit*blk.number_of_units

        # ------------------------------------------------------
        # Material of construction factor Eq. 22.44 in the reference
        blk.material_factor = pyo.Var(
            initialize=3.5,
            domain=pyo.NonNegativeReals,
            doc='Construction material correction factor')

        hx_material_factor_dict = {
            HXMaterial.CS_CS: {"A": 0.00, "B": 0.00},
            HXMaterial.CS_Brass: {"A": 1.08, "B": 0.05},
            HXMaterial.CS_SS: {"A": 1.75, "B": 0.13},
            HXMaterial.CS_Monel: {"A": 2.10, "B": 0.13},
            HXMaterial.CS_Ti: {"A": 5.20, "B": 0.16},
            HXMaterial.CS_CrMo: {"A": 1.55, "B": 0.05},
            HXMaterial.CrMo_CrMo: {"A": 1.70, "B": 0.07},
            HXMaterial.SS_SS: {"A": 2.70, "B": 0.07},
            HXMaterial.Monel_Monel: {"A": 3.30, "B": 0.08},
            HXMaterial.Ti_Ti: {"A": 9.60, "B": 0.06}}

        mf = hx_material_factor_dict[material_type]

        def hx_material_fact_rule(self):
            if material_type == HXMaterial.CS_CS:
                return blk.material_factor == 1
            else:
                return (blk.material_factor ==
                        mf["A"] + (area_unit/(100*pyo.units.ft**2))**mf["B"])
        blk.hx_material_eqn = pyo.Constraint(rule=hx_material_fact_rule)

        # ------------------------------------------------------
        # Pressure factor calculation
        t0 = blk.unit_model.flowsheet().time.first()

        # Assume higher pressure fluid is tube side
        blk.pressure_factor = pyo.Var(initialize=1,
                                      domain=pyo.NonNegativeReals,
                                      doc='Pressure design factor')

        try:
            tube_props = blk.unit_model.tube.properties_in[t0]
        except AttributeError:
            # Assume HX1D
            if (blk.unit_model.config.flow_type ==
                    HeatExchangerFlowPattern.cocurrent):
                inlet_x = 0
            else:
                inlet_x = 1
            tube_props = blk.unit_model.tube.properties[0, inlet_x]

        # Pressure units must be in psig
        pressure = (pyo.units.convert(tube_props.pressure,
                                      to_units=pyo.units.psi) -
                    pyo.units.convert(1*pyo.units.atm, to_units=pyo.units.psi))

        def hx_P_factor(blk):
            # Equation valid from 600 pisg to 3000 psig
            # return self.pressure_factor == (
            #     0.8510 + 0.1292*(pressure/600) + 0.0198*(pressure/600)**2)
            # Equation valid from 100 pisg to 2000 psig
            return blk.pressure_factor == (
                0.9803 + 0.0180 * (pressure/(100*pyo.units.psi)) +
                0.0017 * (pressure/(100*pyo.units.psi))**2)
        blk.p_factor_eq = pyo.Constraint(rule=hx_P_factor)

        # Total capital cost equation
        def hx_CP_rule(blk):
            return blk.capital_cost == (
                blk.pressure_factor*blk.material_factor *
                blk.length_factor*blk.base_cost)
        blk.capital_cost_constraint = pyo.Constraint(rule=hx_CP_rule)

    def cost_vessel(blk,
                    vertical=False,
                    material_type=VesselMaterial.CS,
                    shell_thickness=1.25*pyo.units.inch,
                    weight_limit=1,
                    aspect_ratio_range=1,
                    include_platforms_ladders=True,
                    vessel_diameter=None,
                    vessel_length=None,
                    number_of_units=1,
                    number_of_trays=None,
                    tray_material=TrayMaterial.CS,
                    tray_type=TrayType.Sieve):
        """
        Generic vessel costing method.

        Args:
            vertical - alignment of vessel; vertical if True, horizontal if
                       False (default=False).
            material_type - VesselMaterial Enum indicating material of
                            construction, default = VesselMaterial.CS.
            shell_thickness - thickness of vessel shell, including pressure
                              allowance. Default = 1.25 inches.
            weight_limit - 1: (default) 1000 to 920,000 lb, 2: 4200 to 1M lb.
                           Option 2 is only valid for vertical vessels.
            aspect_ratio_range - vertical vessels only, default = 1;
                                 1: 3 < D < 21 ft, 12 < L < 40 ft,
                                 2: 3 < D < 24 ft; 27 < L < 170 ft.
            include_platforms_ladders - whether to include platforms and
                                        ladders in costing , default = True.
            vessel_diameter - Pyomo component representing vessel diameter.
                              If not provided, assumed to be named "diameter"
            vessel_length - Pyomo component representing vessel length.
                            If not provided, assumed to be named "length".
            number_of_units - Integer or Pyomo component representing the
                               number of parallel units to be costed,
                               default = 1.
            number_of_trays - Pyomo component representing the number of
                              distillation trays in vessel (default=None)
            tray_material - Only required if number_of_trays is not None.
                            TrayMaterial Enum indicating material of
                            construction for distillation trays, default =
                            TrayMaterial.CS.
            tray_type - Only required if number_of_trays is not None.
                        TrayType Enum indicating type of distillation trays
                        to use, default = TrayMaterial.Sieve.
        """
        # Build generic costing variables
        blk.base_cost_per_unit = pyo.Var(initialize=1e5,
                                         domain=pyo.NonNegativeReals,
                                         units=pyo.units.USD_500,
                                         doc='Base cost per unit')

        blk.capital_cost = pyo.Var(initialize=1e4,
                                   domain=pyo.NonNegativeReals,
                                   bounds=(0, None),
                                   units=pyo.units.USD_500,
                                   doc='Capital cost of all units')

        # Check arguments
        if material_type not in VesselMaterial:
            raise ConfigurationError(
                f"{blk.unit_model.name} received invalid argument for "
                f"material_type: {material_type}. Argument must be a member "
                "of the VesselMaterial Enum.")

        if not vertical and number_of_trays is not None:
            raise ConfigurationError(
                f"{blk.unit_model.name} distillation trays are only supported "
                "for vertical vessels.")

        if weight_limit == 2 and not vertical:
            raise ConfigurationError(
                f"{blk.unit_model.name} weight_limit option 2 is only valid "
                "for vertical vessels.")
        elif weight_limit not in [1, 2]:
            raise ConfigurationError(
                f"{blk.unit_model.name} weight_limit argument must be 1 or 2;"
                f" recieved {weight_limit}.")

        # Check references to diameter and length
        if vessel_diameter is None:
            try:
                vessel_diameter = blk.unit_model.diameter
            except AttributeError:
                raise ConfigurationError(
                    f"{blk.unit_model.name} does not have a component named "
                    "diameter. Please provide a reference to the vessel "
                    "diameter as an argument to cost_unit.")

        if vessel_length is None:
            try:
                vessel_length = blk.unit_model.length
            except AttributeError:
                raise ConfigurationError(
                    f"{blk.unit_model.name} does not have a component named "
                    "length. Please provide a reference to the vessel "
                    "length as an argument to cost_unit.")

        D_in = pyo.units.convert(vessel_diameter, to_units=pyo.units.inch)
        L_in = pyo.units.convert(vessel_length, to_units=pyo.units.inch)

        # Material densities in lb/cubic inch
        material_factor_dict = {
            VesselMaterial.CS: {"factor": 1.0, "density": 0.284},
            VesselMaterial.LowAlloy: {"factor": 1.2, "density": 0.271},
            VesselMaterial.SS304: {"factor": 1.7, "density": 0.270},
            VesselMaterial.SS316: {"factor": 2.1, "density": 0.276},
            VesselMaterial.Carpenter: {"factor": 3.2, "density": 0.29},
            VesselMaterial.Nickel200: {"factor": 5.4, "density": 0.3216},
            VesselMaterial.Monel400: {"factor": 3.6, "density": 0.319},
            VesselMaterial.Inconel600: {"factor": 3.9, "density": 0.3071},
            VesselMaterial.Incoloy825: {"factor": 3.7, "density": 0.2903},
            VesselMaterial.Titanium: {"factor": 7.7, "density": 0.1628}}
        material_factors = material_factor_dict[material_type]

        # Users should calculate the pressure design based shell thickness
        # Pressure factor assumed to be included in thickness
        blk.shell_thickness = pyo.Param(mutable=True,
                                        doc='Shell thickness',
                                        units=pyo.units.inch)
        blk.shell_thickness.set_value(shell_thickness)
        blk.material_factor = pyo.Param(
            initialize=material_factors["factor"],
            mutable=True,
            doc='Construction material correction factor')
        blk.material_density = pyo.Param(
            initialize=material_factors["density"],
            mutable=True,
            doc='Density of the metal',
            units=pyo.units.pound/pyo.units.inch**3)

        # Calculate weight of vessel
        blk.weight = pyo.Var(initialize=1000,
                             domain=pyo.NonNegativeReals,
                             doc='Weight of vessel in lb',
                             units=pyo.units.pound)

        def weight_rule(blk):
            return blk.weight == (
                Constants.pi*(D_in + blk.shell_thickness) *
                (L_in + 0.8*D_in)*blk.shell_thickness*blk.material_density)
        blk.weight_eq = pyo.Constraint(rule=weight_rule)

        # Base Vessel cost
        # Alpha factors for correlation
        # 2nd key is weight_limit option
        alpha_dict = {
            "H": {1: {1: 8.9552, 2: -0.2330, 3: 0.04333}, 2: None},
            "V": {1: {1: 7.0132, 2: 0.18255, 3: 0.02297},
                  2: {1: 7.2756, 2: 0.18255, 3: 0.02297}}}

        if vertical:
            alpha = alpha_dict["V"][weight_limit]
        else:
            alpha = alpha_dict["H"][weight_limit]

        def base_cost_rule(blk):
            return blk.base_cost_per_unit == (
                pyo.exp(alpha[1] +
                        alpha[2]*(pyo.log(blk.weight/pyo.units.pound)) +
                        alpha[3]*(pyo.log(blk.weight/pyo.units.pound)**2)) *
                pyo.units.USD_500)
        blk.base_cost_constraint = pyo.Constraint(rule=base_cost_rule)

        # Add platform and ladder costs if required
        if include_platforms_ladders:
            SSLWCosting._cost_platforms_ladders(
                blk,
                vertical=vertical,
                aspect_ratio_range=aspect_ratio_range,
                vessel_diameter=vessel_diameter,
                vessel_length=vessel_length)

        # Add distillation trays costs if required
        if number_of_trays is not None:
            SSLWCosting._cost_distillation_trays(
                blk,
                tray_material=tray_material,
                tray_type=tray_type,
                vessel_diameter=vessel_diameter,
                number_of_trays=number_of_trays)

        # Total capital cost of vessel and ancilliary equipment
        def capital_cost_rule(blk):
            cost_expr = blk.material_factor*blk.base_cost_per_unit
            if include_platforms_ladders:
                cost_expr += blk.base_cost_platforms_ladders
            if number_of_trays is not None:
                cost_expr += blk.base_cost_trays

            return blk.capital_cost == cost_expr*number_of_units
        blk.capital_cost_constraint = pyo.Constraint(rule=capital_cost_rule)

    def _cost_platforms_ladders(blk,
                                vertical,
                                aspect_ratio_range,
                                vessel_diameter,
                                vessel_length):
        """
        Method for calculating costs of platforms and ladders
        """
        blk.base_cost_platforms_ladders = pyo.Var(
            initialize=1000,
            units=pyo.units.USD_500,
            domain=pyo.NonNegativeReals,
            doc='Base cost of platforms and ladders')

        D = pyo.units.convert(vessel_diameter, to_units=pyo.units.foot)
        L = pyo.units.convert(vessel_length, to_units=pyo.units.foot)

        if not vertical:
            def CPL_rule(blk):
                return blk.base_cost_platforms_ladders == (
                    2005*(D/pyo.units.foot)**0.20294*pyo.units.USD_500)
        else:
            def CPL_rule(blk):
                if aspect_ratio_range == 1:
                    return blk.base_cost_platforms_ladders == (
                        361.8*(D/pyo.units.foot)**0.73960 *
                        (L/pyo.units.foot)**0.70684 *
                        pyo.units.USD_500)
                elif aspect_ratio_range == 2:
                    return blk.base_cost_platforms_ladders == (
                        309.9*(D/pyo.units.foot)**0.63316 *
                        (L/pyo.units.foot)**0.80161 *
                        pyo.units.USD_500)
                else:
                    raise ConfigurationError(
                        f"{blk.unit_model.name} recieved invalid value for "
                        f"aspect_ratio_range argument: {aspect_ratio_range}. "
                        "Value must be 1 or 2.")
        blk.cost_platforms_ladders_eq = pyo.Constraint(rule=CPL_rule)

    def _cost_distillation_trays(blk,
                                 tray_material,
                                 tray_type,
                                 vessel_diameter,
                                 number_of_trays):

        # Check arguments
        if tray_material not in TrayMaterial:
            raise ConfigurationError(
                f"{blk.unit_model.name} received invalid argument for "
                f"tray_material: {tray_material}. Argument must be a member "
                "of the TrayMaterial Enum.")
        if tray_type not in TrayType:
            raise ConfigurationError(
                f"{blk.unit_model.name} received invalid argument for "
                f"tray_type: {tray_type}. Argument must be a member "
                "of the TrayType Enum.")

        D = pyo.units.convert(vessel_diameter, to_units=pyo.units.foot)

        blk.base_cost_trays = pyo.Var(initialize=1e6,
                                      units=pyo.units.USD_500,
                                      domain=pyo.NonNegativeReals,
                                      doc='Purchase cost of trays')

        tray_type_dict = {TrayType.Sieve: 1,
                          TrayType.Valve: 1.18,
                          TrayType.BubbleCap: 1.87}
        blk.tray_type_factor = pyo.Param(initialize=tray_type_dict[tray_type],
                                         mutable=True,
                                         doc='Tray type factor')

        blk.tray_material_factor = pyo.Var(
            initialize=1.0,
            doc='FTM material of construction factor for trays')
        # Alpha parameters for material factor
        alpha = {TrayMaterial.CS: {1: 1, 2: 0},
                 TrayMaterial.SS303: {1: 1.189, 2: 0.0577},
                 TrayMaterial.SS316: {1: 1.401, 2: 0.0724},
                 TrayMaterial.Carpenter: {1: 1.525, 2: 0.0788},
                 TrayMaterial.Monel: {1: 2.306, 2: 0.1120}}

        # Calculating tray factor value
        # Column diameter in ft, eqn. valid for 2 to 16 ft
        def mt_factor_rule(blk):
            return blk.tray_material_factor == (
                alpha[tray_material][1] +
                alpha[tray_material][2]*D/pyo.units.foot)
        blk.tray_material_factor_eq = pyo.Constraint(rule=mt_factor_rule)

        # Calcluate cost factor for number of trays
        blk.number_trays_factor = pyo.Var(
            initialize=1,
            units=pyo.units.dimensionless,
            doc='Cost factor for number of trays')

        def num_trays_rule(blk):
            return blk.number_trays_factor == smooth_max(
                1, 2.25/(1.0414**number_of_trays))
        blk.num_tray_factor_constraint = pyo.Constraint(rule=num_trays_rule)

        # Calcualte base cost of a single tray
        blk.base_cost_per_tray = pyo.Var(initialize=1e4,
                                         domain=pyo.NonNegativeReals,
                                         units=pyo.units.USD_500,
                                         doc='Base cost of a single tray')

        def single_tray_rule(blk):
            return blk.base_cost_per_tray == (
                468.00*pyo.exp(0.1739*D/pyo.units.foot)*pyo.units.USD_500)
        blk.single_tray_cost_constraint = pyo.Constraint(rule=single_tray_rule)

        # Capital cost of trays
        def cost_trays_rule(blk):
            return blk.base_cost_trays == (
                number_of_trays*blk.number_trays_factor*blk.tray_type_factor *
                blk.tray_material_factor*blk.base_cost_per_tray)
        blk.tray_costing_constraint = pyo.Constraint(rule=cost_trays_rule)

    def cost_vertical_vessel(blk,
                             material_type=VesselMaterial.CS,
                             shell_thickness=1.25*pyo.units.inch,
                             weight_limit=1,
                             aspect_ratio_range=1,
                             include_platforms_ladders=True,
                             vessel_diameter=None,
                             vessel_length=None,
                             number_of_units=1,
                             number_of_trays=None,
                             tray_material=TrayMaterial.CS,
                             tray_type=TrayType.Sieve):
        """
        Specific case of vessel costing method for vertical vessels.

        Args:
            material_type - VesselMaterial Enum indicating material of
                            construction, default = VesselMaterial.CS.
            shell_thickness - thickness of vessel shell, including pressure
                              allowance. Default = 1.25 inches.
            weight_limit - 1: (default) 1000 to 920,000 lb, 2: 4200 to 1M lb.
                           Option 2 is only valid for vertical vessels.
            aspect_ratio_range - vertical vessels only, default = 1;
                                 1: 3 < D < 21 ft, 12 < L < 40 ft,
                                 2: 3 < D < 24 ft; 27 < L < 170 ft.
            include_platforms_ladders - whether to include platforms and
                                        ladders in costing , default = True.
            vessel_diameter - Pyomo component representing vessel diameter.
                              If not provided, assumed to be named "diameter"
            vessel_length - Pyomo component representing vessel length.
                            If not provided, assumed to be named "length".
            number_of_units - Integer or Pyomo component representing the
                               number of parallel units to be costed,
                               default = 1.
            number_of_trays - Pyomo component representing the number of
                              distillation trays in vessel (default=None)
            tray_material - Only required if number_of_trays is not None.
                            TrayMaterial Enum indicating material of
                            construction for distillation trays, default =
                            TrayMaterial.CS.
            tray_type - Only required if number_of_trays is not None.
                        TrayType Enum indicating type of distillation trays
                        to use, default = TrayMaterial.Sieve.
        """
        SSLWCosting.cost_vessel(blk,
                                vertical=True,
                                material_type=VesselMaterial.CS,
                                shell_thickness=1.25*pyo.units.inch,
                                weight_limit=1,
                                aspect_ratio_range=1,
                                include_platforms_ladders=True,
                                vessel_diameter=None,
                                vessel_length=None,
                                number_of_units=1,
                                number_of_trays=None,
                                tray_material=TrayMaterial.CS,
                                tray_type=TrayType.Sieve)

    def cost_horizontal_vessel(blk,
                               material_type=VesselMaterial.CS,
                               shell_thickness=1.25*pyo.units.inch,
                               include_platforms_ladders=True,
                               vessel_diameter=None,
                               vessel_length=None,
                               number_of_units=1):
        """
        Specific case of vessel costing method for horizontal vessels.

        Arguments which do not apply ot horizontal vessels are excluded.

        Args:
            material_type - VesselMaterial Enum indicating material of
                            construction, default = VesselMaterial.CS.
            shell_thickness - thickness of vessel shell, including pressure
                              allowance. Default = 1.25 inches.
            include_platforms_ladders - whether to include platforms and
                                        ladders in costing , default = True.
            vessel_diameter - Pyomo component representing vessel diameter.
                              If not provided, assumed to be named "diameter"
            vessel_length - Pyomo component representing vessel length.
                            If not provided, assumed to be named "length".
            number_of_units - Integer or Pyomo component representing the
                               number of parallel units to be costed,
                               default = 1.
        """
        SSLWCosting.cost_vessel(blk,
                                vertical=False,
                                material_type=VesselMaterial.CS,
                                shell_thickness=1.25*pyo.units.inch,
                                weight_limit=1,
                                include_platforms_ladders=True,
                                vessel_diameter=None,
                                vessel_length=None,
                                number_of_units=1)

    def cost_fired_heater(blk,
                          heat_source=HeaterSource.fuel,
                          material_type=HeaterMaterial.CS,
                          integer=True):
        """
        Generic costing method for fired heaters.

        Args:
            heat_source - HeaterSource Enum indicating type of source of heat,
                          default = HeaterSource.fuel.
            material_type - HeaterMaterial Enum indicating material of
                            construction, default = HeaterMaterial.CS.
            integer - whether the number of units should be constrained to be
                      an integer or not (default = True).
        """
        # Validate arguments
        if heat_source not in HeaterSource:
            raise ConfigurationError(
                f"{blk.unit_model.name} received invalid argument for "
                f"heat_source: {heat_source}. Argument must be a member of "
                "the HeaterSource Enum.")
        if material_type not in HeaterMaterial:
            raise ConfigurationError(
                f"{blk.unit_model.name} received invalid argument for "
                f"material_type: {material_type}. Argument must be a member "
                "of the HeaterMaterial Enum.")

        # Build generic costing variables
        _make_common_vars(blk, integer)

        # Convert pressure to psi,g
        t0 = blk.unit_model.flowsheet().time.first()

        P = (pyo.units.convert(
                blk.unit_model.control_volume.properties_in[t0].pressure,
                to_units=pyo.units.psi) -
             pyo.units.convert(
                 1*pyo.units.atm,
                 to_units=pyo.units.psi))

        # Convert heat duty to BTU/hr
        Q = (pyo.units.convert(blk.unit_model.heat_duty[t0],
                               to_units=pyo.units.BTU/pyo.units.hr) /
             blk.number_of_units)

        # Material factor
        material_factor_dict = {HeaterMaterial.CS: 1.0,
                                HeaterMaterial.CrMo: 1.4,
                                HeaterMaterial.SS: 1.7}
        blk.material_factor = pyo.Param(
            initialize=material_factor_dict[material_type],
            domain=pyo.NonNegativeReals,
            doc='Construction material correction factor')

        # Pressure deisgn factor calculation
        blk.pressure_factor = pyo.Var(initialize=1.1,
                                      domain=pyo.NonNegativeReals,
                                      doc='Pressure design factor')

        def p_factor_rule(blk):
            return blk.pressure_factor == (
                0.986 - 0.0035*(P/(500.00*pyo.units.psi)) +
                0.0175*(P/(500.00*pyo.units.psi))**2)
        blk.pressure_factor_eq = pyo.Constraint(rule=p_factor_rule)

        def CB_rule(blk):
            if heat_source == HeaterSource.fuel:
                bc_expr = pyo.exp(
                    0.32325 + 0.766*pyo.log(Q/pyo.units.BTU*pyo.units.hr))
            elif heat_source == HeaterSource.reformer:
                bc_expr = (
                    0.859*(Q/pyo.units.BTU*pyo.units.hr)**0.81)
            elif heat_source == HeaterSource.pyrolysis:
                bc_expr = (
                    0.650*(Q/pyo.units.BTU*pyo.units.hr)**0.81)
            elif heat_source == HeaterSource.hotWater:
                bc_expr = pyo.exp(
                    9.593 -
                    0.3769*pyo.log((Q/pyo.units.BTU*pyo.units.hr)) +
                    0.03434*pyo.log((Q/pyo.units.BTU*pyo.units.hr))**2)
            elif heat_source == HeaterSource.salts:
                bc_expr = (
                    12.32*(Q/pyo.units.BTU*pyo.units.hr)**0.64)
            elif heat_source == HeaterSource.DowthermA:
                bc_expr = (
                    12.74*(Q/pyo.units.BTU*pyo.units.hr)**0.65)
            elif heat_source == HeaterSource.steamBoiler:
                bc_expr = (
                    0.367*(Q/pyo.units.BTU*pyo.units.hr)**0.77)

            return blk.base_cost_per_unit == bc_expr*pyo.units.USD_500

        blk.base_cost_per_unit_eq = pyo.Constraint(rule=CB_rule)

        @blk.Expression(doc="Base cost for all units installed")
        def base_cost(blk):
            return blk.base_cost_per_unit * blk.number_of_units

        # Total capital cost of heater(s)
        def CP_rule(blk):
            return blk.capital_cost == (
                blk.material_factor*blk.pressure_factor*blk.base_cost)
        blk.capital_cost_constraint = pyo.Constraint(rule=CP_rule)

    def cost_turbine(blk,
                     integer=True):
        """
        Generic costing method for turbines.

        Args:
            integer - whether the number of units should be constrained to be
                      an integer or not (default = True).
        """
        # Confirm that unit is a turbine
        if blk.unit_model.config.compressor:
            raise TypeError("cost_turbine method is only appropriate for "
                            "pressure changers with the compressor argument "
                            "equal to False.")

        # Build costing variables
        blk.capital_cost = pyo.Var(initialize=1e4,
                                   domain=pyo.NonNegativeReals,
                                   bounds=(0, None),
                                   units=pyo.units.USD_500,
                                   doc='Capital cost of all units')

        if integer is True:
            domain = pyo.Integers
        else:
            domain = pyo.NonNegativeReals
        blk.number_of_units = pyo.Var(
            initialize=1,
            domain=domain,
            bounds=(1, 100),
            doc="Number of units to install.")
        blk.number_of_units.fix(1)

        work = (-blk.unit_model.work_mechanical[
            blk.unit_model.flowsheet().time.first()])

        work_hp = pyo.units.convert(work/blk.number_of_units,
                                    to_units=pyo.units.hp)

        def CP_rule(blk):
            return blk.capital_cost == (
                blk.number_of_units*530*(work_hp/pyo.units.hp)**0.81 *
                pyo.units.USD_500)
        blk.capital_cost_constraint = pyo.Constraint(rule=CP_rule)

    def cost_pump(blk,
                  pump_type=PumpType.centrifugal,
                  material_type=PumpMaterial.SS,
                  pump_type_factor=1.4,
                  motor_type=PumpMotorType.Open,
                  integer=True):
        """
        Generic costing method for pumps.

        Args:
            pump_type - PumpType Enum indicating type of type of equipment,
                      default = PumpType.centrifugal.
            material_type - PumpMaterial Enum indicating material of
                            construction, default = PumpMaterial.SS.
                            Material type is tied to PumpType.
            pump_type_factor - empirical factor for centrigual pumps based on
                               table in source. Valid values are [1.1, 1.2, "
                               1.3, 1.4 (default), 2.1, 2.2].
            motor_type - PumpMotorType Enum indicating type of type of motor
                         to be used, default = PumpMotorType.Open.
            integer - whether the number of units should be constrained to be
                      an integer or not (default = True).
        """
        # Confirm that unit is a pump/compressor
        if not blk.unit_model.config.compressor:
            raise TypeError("cost_pump method is only appropriate for "
                            "pressure changers with the compressor argument "
                            "equal to True.")

        # Validate argument combinations
        if pump_type == PumpType.reciprocating:
            if material_type not in [PumpMaterial.ductileIron,
                                     PumpMaterial.NiAlBronze,
                                     PumpMaterial.CS,
                                     PumpMaterial.SS]:
                raise ConfigurationError(
                    f"{blk.name} invalid combination of arguments. If "
                    "pump_type == PumpType.reciprocating then material_type "
                    "must be one of PumpMaterial.ductileIron, "
                    "PumpMaterial.NiAlBronze, PumpMaterial.CS, or "
                    "PumpMaterial.SS.")
        else:
            if material_type in [PumpMaterial.NiAlBronze, PumpMaterial.CS]:
                raise ConfigurationError(
                    f"{blk.name} invalid combination of arguments. "
                    "PumpMaterial.NiAlBronze and PumpMaterial.CS material_type"
                    " are only valid for pump_type == PumpType.reciprocating.")

        # Add common variables
        blk.capital_cost = pyo.Var(initialize=1e4,
                                   domain=pyo.NonNegativeReals,
                                   bounds=(0, None),
                                   units=pyo.units.USD_500,
                                   doc='Capital cost of all units')

        if integer is True:
            domain = pyo.Integers
        else:
            domain = pyo.NonNegativeReals
        blk.number_of_units = pyo.Var(
            initialize=1,
            domain=domain,
            bounds=(1, 100),
            doc="Number of units to install.")
        blk.number_of_units.fix(1)

        t0 = blk.unit_model.flowsheet().time.first()

        work = pyo.units.convert(blk.unit_model.work_mechanical[t0] /
                                 blk.number_of_units,
                                 to_units=pyo.units.hp)

        # Pressure units required in lbf/ft^2
        deltaP = pyo.units.convert(
            blk.unit_model.deltaP[t0],
            to_units=pyo.units.psi*pyo.units.inch**2/pyo.units.foot**2)

        # Mass density units required in lb/ft^3
        dens_mass = pyo.units.convert(
            blk.unit_model.control_volume.properties_in[t0].dens_mass,
            to_units=pyo.units.pound/pyo.units.foot**3)

        # Volumetric flow units required in gpm
        Q = pyo.units.convert(
            blk.unit_model.control_volume.properties_in[t0].flow_vol /
            blk.number_of_units,
            to_units=pyo.units.gallon/pyo.units.minute)

        # Calculate pump head
        blk.pump_head = pyo.Var(
            initialize=10,
            domain=pyo.NonNegativeReals,
            doc='Pump Head in feet of fluid flowing (Pressure rise/density)',
            units=pyo.units.pound_force*pyo.units.foot/pyo.units.pound)

        def p_head_rule(blk):
            return blk.pump_head == deltaP/dens_mass
        blk.pump_head_eq = pyo.Constraint(rule=p_head_rule)

        # Pump size factor: S = Q(H)**0.5
        # Q = is the flow rate through the pump in gallons per minute
        # H = pump head in feet of flowing (pressure rise/liquid density)
        blk.size_factor = pyo.Var(
            initialize=10000,
            domain=pyo.NonNegativeReals,
            doc='Pump size factor, f(Q,pump_head)')

        def p_s_factor_rule(blk):
            return blk.size_factor == (
                Q*pyo.units.minute/pyo.units.gallon *
                (blk.pump_head *
                 pyo.units.pound/pyo.units.pound_force/pyo.units.foot)**0.5)
        blk.size_factor_eq = pyo.Constraint(rule=p_s_factor_rule)

        # Material factor
        if pump_type != PumpType.reciprocating:
            material_factor_dict = {PumpMaterial.castIron: 1.00,
                                    PumpMaterial.ductileIron: 1.15,
                                    PumpMaterial.castSteel: 1.35,
                                    PumpMaterial.bronze: 1.90,
                                    PumpMaterial.SS: 2.00,
                                    PumpMaterial.HastelloyC: 2.95,
                                    PumpMaterial.Monel: 3.30,
                                    PumpMaterial.Nickel: 3.50,
                                    PumpMaterial.Titanium: 9.70}
        else:
            material_factor_dict = {PumpMaterial.ductileIron: 1.00,
                                    PumpMaterial.SS: 2.20,
                                    PumpMaterial.NiAlBronze: 1.15,
                                    PumpMaterial.CS: 1.50}
        blk.material_factor = pyo.Param(
            initialize=material_factor_dict[material_type],
            mutable=True,
            doc='Construction material correction factor')

        # Pump type factor
        if pump_type == PumpType.centrifugal:
            pump_type_factor_dict = {1.1: 1.00,
                                     1.2: 1.50,
                                     1.3: 1.70,
                                     1.4: 2.00,
                                     2.1: 2.70,
                                     2.2: 8.90}
            try:
                ptf = pump_type_factor_dict[pump_type_factor]
            except KeyError:
                raise ConfigurationError(
                    f"{blk.name} invalid value for pump_type_factor argument: "
                    f"{pump_type_factor}. Value must be one of [1.1, 1.2, 1.3,"
                    " 1.4, 2.1, 2.2]")
        else:
            ptf = 1

        blk.FT = pyo.Param(mutable=True,
                           initialize=ptf,
                           doc='Pump type factor')

        # Base pump cost per unit
        blk.base_pump_cost_per_unit = pyo.Var(
            initialize=1e5,
            domain=pyo.NonNegativeReals,
            units=pyo.units.USD_394,
            doc='Base cost of pump (less motor) per unit')

        def base_pump_cost_rule(blk):
            if pump_type == PumpType.centrifugal:
                bpc = pyo.exp(9.7171 -
                              0.6019*pyo.log(blk.size_factor) +
                              0.0519*pyo.log(blk.size_factor)**2)
            elif pump_type == PumpType.externalGear:
                bpc = pyo.exp(
                    7.6964 +
                    0.1986*pyo.log(Q/(pyo.units.gallon/pyo.units.minute)) +
                    0.0291*pyo.log(Q/(pyo.units.gallon/pyo.units.minute))**2)
            elif pump_type == PumpType.reciprocating:
                bpc = pyo.exp(7.8103 +
                              0.26986*pyo.log(work/pyo.units.hp) +
                              0.06718*pyo.log(work/pyo.units.hp)**2)
            return blk.base_pump_cost_per_unit == bpc*pyo.units.USD_394
        blk.base_pump_cost_per_unit_eq = pyo.Constraint(
            rule=base_pump_cost_rule)

        @blk.Expression(doc="Base cost for all pumps (less motors)")
        def base_pump_cost(blk):
            return blk.base_pump_cost_per_unit * blk.number_of_units

        blk.pump_capital_cost = pyo.Var(
            initialize=100000,
            domain=pyo.NonNegativeReals,
            units=pyo.units.USD_500,
            doc='Capital cost of pumps (less motors)')

        def CP_pump_rule(blk):
            return blk.pump_capital_cost == pyo.units.convert(
                blk.FT*blk.material_factor*blk.base_pump_cost,
                to_units=pyo.units.USD_500)
        blk.pump_capital_cost_eq = pyo.Constraint(rule=CP_pump_rule)

        # Motor Costs
        pump_motor_type_dict = {PumpMotorType.Open: 1,
                                PumpMotorType.Enclosed: 1.4,
                                PumpMotorType.ExplosionProof: 1.8}
        blk.motor_FT = pyo.Param(
            mutable=True,
            initialize=pump_motor_type_dict[motor_type],
            doc='Motor type factor')

        # Efficiency of the electric motor
        eta_m = (0.80 + 0.0319*pyo.log(work/pyo.units.hp) -
                 0.00182*pyo.log(work/pyo.units.hp)**2)
        work_motor = work/eta_m

        blk.base_motor_cost_per_unit = pyo.Var(
            initialize=10000,
            domain=pyo.NonNegativeReals,
            units=pyo.units.USD_394,
            doc='Motor base purchase cost per unit')

        def base_motor_cost_rule(blk):
            return blk.base_motor_cost_per_unit == (pyo.exp(
                5.8259 + 0.13141*pyo.log(work_motor/pyo.units.hp) +
                0.053255*pyo.log(work_motor/pyo.units.hp)**2 +
                0.028628*pyo.log(work_motor/pyo.units.hp)**3 -
                0.0035549*pyo.log(work_motor/pyo.units.hp)**4) *
                pyo.units.USD_394)
        blk.base_motor_cost_eq = pyo.Constraint(rule=base_motor_cost_rule)

        @blk.Expression(doc="Base cost for all motors")
        def motor_base_cost(blk):
            return blk.base_motor_cost_per_unit * blk.number_of_units

        blk.motor_capital_cost = pyo.Var(
            initialize=100000,
            domain=pyo.NonNegativeReals,
            units=pyo.units.USD_500,
            doc='Capital cost of all motors')

        def CP_motor_rule(blk):
            return blk.motor_capital_cost == pyo.units.convert(
                blk.motor_FT*blk.motor_base_cost,
                to_units=pyo.units.USD_500)
        blk.motor_capital_cost_eq = pyo.Constraint(rule=CP_motor_rule)

        # Total capital cost (pump + electrical motor)
        def cp_cost_rule(blk):
            return blk.capital_cost == (blk.pump_capital_cost +
                                        blk.motor_capital_cost)
        blk.capital_cost_constraint = pyo.Constraint(rule=cp_cost_rule)

    # -------------------------------------------------------------------------
    # Map costing methods to unit model classes
    # Here we can provide a dict mapping unit model classes to costing methods
    # Even better, this is inheritance aware so e.g. Pump will first look for a
    # method assigned to Pump, and then fall back to PressureChanger
    unit_mapping = {CSTR: cost_vertical_vessel,
                    Flash: cost_vertical_vessel,
                    Heater: cost_fired_heater,
                    HeatExchanger: cost_heat_exchanger,
                    HeatExchanger1D: cost_heat_exchanger,
                    HeatExchangerNTU: cost_heat_exchanger,
                    PFR: cost_horizontal_vessel,
                    Pump: cost_pump,
                    StoichiometricReactor: cost_horizontal_vessel}


# -----------------------------------------------------------------------------
def _make_common_vars(blk, integer=True):
    # Build generic costing variables (most costing models need these vars)
    blk.base_cost_per_unit = pyo.Var(initialize=1e5,
                                     domain=pyo.NonNegativeReals,
                                     units=pyo.units.USD_500,
                                     doc='Base cost per unit')

    blk.capital_cost = pyo.Var(initialize=1e4,
                               domain=pyo.NonNegativeReals,
                               bounds=(0, None),
                               units=pyo.units.USD_500,
                               doc='Capital cost of all units')

    if integer is True:
        domain = pyo.Integers
    else:
        domain = pyo.NonNegativeReals
    blk.number_of_units = pyo.Var(
        initialize=1,
        domain=domain,
        bounds=(1, 100),
        doc="Number of units to install.")
    blk.number_of_units.fix(1)
