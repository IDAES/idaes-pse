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
    HeatExchanger,
    HeatExchanger1D,
    HeatExchangerNTU)
from idaes.generic_models.unit_models.heat_exchanger \
    import HeatExchangerFlowPattern
from idaes.core.util.misc import register_units_of_measurement

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


class SSLWCosting(CostingPackageBase):

    # DEfine currency and conversion rates based on CE Index
    currency_units = [
        "USD_500 = [currency]",  # base USD at CE Index of 500
        "USD2010 = 500/550.8 * USD_500",
        "USD2011 = 500/585.7 * USD_500",
        "USD2012 = 500/584.6 * USD_500",
        "USD2013 = 500/567.3 * USD_500",
        "USD2014 = 500/576.1 * USD_500",
        "USD2015 = 500/556.8 * USD_500",
        "USD2016 = 500/541.7 * USD_500",
        "USD2017 = 500/567.5 * USD_500",
        "USD2018 = 500/671.1 * USD_500",
        "USD2019 = 500/680.0 * USD_500"]
    # Need to register these units before we continue
    register_units_of_measurement(currency_units)

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
        for ach costing method to separate the parameters for each method.
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

    def hx_costing(blk,
                   hx_type=HXType.Utube,
                   material_factor=HXMaterial.SS_SS,
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
                        base_Cost

        Args:
            hx_type - HXType Enum indicating type of heatexchanger design,
                      default = HXType.Utube.
            material factor - HXMaterial Enum indicating material of
                              construction, default = HXMaterial.SS_SS.
            tube length - HXTubeLength Enum indicating length of HX tubes,
                          default = HXTubeLength.TwelveFoot.
        """
        # ------------------------ Heat Exchanger cost ------------------------
        # Heat exchanger cost
        # Build generic costing variables
        _make_common_vars(blk, integer)

        # Length correction factor
        c_fl = {'8ft': 1.25, '12ft': 1.12, '16ft': 1.05, '20ft': 1.00}
        blk.length_factor = pyo.Param(
            mutable=True,
            initialize=c_fl[tube_length],
            doc='HX tube length correction factor, FL')

        # TODO : Check that this is set correctly later
        blk.pressure_factor = pyo.Var(initialize=1,
                                      domain=pyo.NonNegativeReals,
                                      doc='Pressure design factor - FP')

        blk.material_factor = pyo.Var(
            initialize=3.5,
            domain=pyo.NonNegativeReals,
            doc='Construction material correction factor - mat_factor')

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

        blk.base_cost_per_unit

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

        mf = hx_material_factor_dict[material_factor]

        def hx_material_fact_rule(self):
            if material_factor == HXMaterial.CS_CS:
                return blk.material_factor == 1
            else:
                return (blk.material_factor ==
                        mf["A"] + (area_unit/(100*pyo.units.ft**2))**mf["B"])
        blk.hx_material_eqn = pyo.Constraint(rule=hx_material_fact_rule)

        # ------------------------------------------------------
        # Pressure factor calculation
        # assume higher pressure fluid is tube side
        try:
            tube_props = blk.unit_model.tube.properties_in[0]
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

        # ---------------------------------------------------------
        # purchase cost equation
        def hx_CP_rule(blk):
            return blk.capital_cost == (
                blk.pressure_factor*blk.material_factor *
                blk.length_factor*blk.base_cost)
        blk.capital_cost_constraint = pyo.Constraint(rule=hx_CP_rule)

    # Map costing methods to unit model classes
    # Here we can provide a dict mapping unit model classes to costing methods
    # Even better, this is inheritance aware so e.g. Pump will first look for a
    # method assigned to Pump, and then fall back to PressureChanger
    unit_mapping = {HeatExchanger: hx_costing,
                    HeatExchanger1D: hx_costing,
                    HeatExchangerNTU: hx_costing}


def _make_common_vars(blk, integer=True):
    # Build generic costing variables (all costing models need these vars)
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
