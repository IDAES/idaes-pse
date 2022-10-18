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
# =============================================================================

from pyomo.environ import (
    Integers,
    Constraint,
    Var,
    Param,
    exp,
    log,
    NonNegativeReals,
    units as pyunits,
)
from idaes.core.util.constants import Constants as const
from idaes.core.util.exceptions import ConfigurationError
from pyomo.core.base.units_container import _PyomoUnit
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from idaes.core.util import scaling as iscale

# Some more information about this module
__author__ = "Miguel Zamarripa"


def global_costing_parameters(self, year=None, integer_n_units=False):
    if year is None:
        year = "2018"
    # Cost index $/year (method argument or 2018 default)
    ce_index_dic = {
        "2019": 680,
        "2018": 671.1,
        "2017": 567.5,
        "2016": 541.7,
        "2015": 556.8,
        "2014": 576.1,
        "2013": 567.3,
        "2012": 584.6,
        "2011": 585.7,
        "2010": 550.8,
    }

    self.CE_index = Param(
        mutable=True,
        initialize=ce_index_dic[year],
        doc="Chemical Engineering Plant Cost Index $ year",
    )

    if integer_n_units not in [True, False]:
        raise ValueError(
            "{} - integer_n_units must be" " True or False".format(self.name)
        )
    self.integer = integer_n_units


def _make_vars(self):
    # build generic costing variables (all costing models need these vars)
    self.base_cost_per_unit = Var(
        initialize=1e5, domain=NonNegativeReals, doc="Unit Base Cost cost in $"
    )
    self.purchase_cost = Var(
        initialize=1e4, domain=NonNegativeReals, doc="Unit Purchase Cost in $"
    )
    # check global costing argument to build an integer
    # or a continuous variable
    try:  # first check flowsheet(), then parent_blocks()
        integer = self.parent_block().flowsheet().costing.integer
    except AttributeError:
        integer = self.parent_block().parent_block().costing.integer

    if integer is True:
        domain = Integers
    else:
        domain = NonNegativeReals
    self.number_of_units = Var(
        initialize=1,
        domain=domain,
        bounds=(1, 100),
        doc="Number of units to install - " "economics of scale",
    )
    # fix number of units by default
    self.number_of_units.fix(1)


def hx_costing(
    self,
    hx_type="U-tube",
    Mat_factor="stainless steel/stainless steel",
    length_factor="12ft",
):
    """
    Heat exchanger costing method.

    Source:
        Process and Product Design Principles: Synthesis, Analysis, and
        Evaluation
        Seider, Seader, Lewin, Windagdo, 3rd Ed. John Wiley and Sons
        Chapter 22. Cost Accounting and Capital Cost Estimation
        22.2 Cost Indexes and Capital Investment

    This method computes the purchase cost (CP) for a shell and tube heat
    exchanger (Eq. 22.43), the model computes the base cost (CB for 4 types
    of heat exchangers, such as floating head, fixed head, U-tube, and
    Kettle vaporizer), construction material factor (Mat_factor),
    pressure design factor (pressure_factor), and
    tube length correction factor (length_factor),
    using Chemical Engineering base cost index of 500.

    Purchase Cost = pressure_factor * Mat_factor * length_factor * Base Cost

    Args:
        hx_type : One of: ``'floating_head'``, ``'fixed_head'``, ``'U-tube'``
        or ``'Kettle_vap'``.
        **Default** - ``'U-tube'``

        material factor : One of: ``'stain_steel'`` or ``'carb_steel'``.
            **Default** - ``'stain_steel'``

        tube length : One of: ``'8ft'``, ``'12ft'``, ``'16ft'`` or
        ``'20ft'``.
            **Default** - ``'12ft'``
    """
    # ------------------------ Heat Exchanger cost ------------------------
    # heat exchanger cost
    # build generic costing variables
    # (base cost, purchase cost, material factor)
    _make_vars(self)

    self.L_factor = Param(
        mutable=True, initialize=1.12, doc="HX tube length correction factor, FL"
    )

    self.pressure_factor = Var(
        initialize=1, domain=NonNegativeReals, doc="Pressure design factor - FP"
    )

    self.material_factor = Var(
        initialize=3.5,
        domain=NonNegativeReals,
        doc="Construction material correction" "factor - Mat_factor",
    )

    self.hx_os = Param(
        mutable=True,
        initialize=1.1,
        doc="HX oversize factor 1.1 to 1.5",
        units=pyunits.ft**-2,
    )

    # select length correction factor
    c_fl = {"8ft": 1.25, "12ft": 1.12, "16ft": 1.05, "20ft": 1.00}
    self.L_factor = c_fl[length_factor]

    # --------------------------------------------------
    # base cost calculation
    #       # select heat exchanger type:
    alf1 = {
        "floating_head": 11.9052,
        "fixed_head": 11.2927,
        "U-tube": 11.3852,
        "Kettle_vap": 12.2052,
    }
    alf2 = {
        "floating_head": 0.8709,
        "fixed_head": 0.8228,
        "U-tube": 0.9186,
        "Kettle_vap": 0.8709,
    }
    alf3 = {
        "floating_head": 0.09005,
        "fixed_head": 0.09861,
        "U-tube": 0.09790,
        "Kettle_vap": 0.09005,
    }

    # checking units of self.parent_block().area
    area = (
        pyunits.convert(self.parent_block().area, to_units=pyunits.ft**2)
        / self.number_of_units
    )

    def hx_cost_rule(self):
        return self.base_cost_per_unit == exp(
            alf1[hx_type]
            - alf2[hx_type] * log(area * self.hx_os)
            + alf3[hx_type] * log(area * self.hx_os) ** 2
        )

    self.base_cost_per_unit_eq = Constraint(rule=hx_cost_rule)

    @self.Expression(doc="Base cost for all units installed")
    def base_cost(self):
        return self.base_cost_per_unit * self.number_of_units

    # ------------------------------------------------------
    # Material of construction factor Eq. 22.44 in the reference
    hx_material_factor_a_dic = {
        "carbon steel/carbon steel": 0.00,
        "carbon steel/brass": 1.08,
        "carbon steel/stainless steel": 1.75,
        "carbon steel/monel": 2.10,
        "carbon steel/titanium": 5.20,
        "carbon steel/Cr-Mo steel": 1.55,
        "Cr-Mo steel/Cr-Mo steel": 1.70,
        "stainless steel/stainless steel": 2.70,
        "monel/monel": 3.30,
        "titanium/titanium": 9.60,
    }

    hx_material_factor_b_dic = {
        "carbon steel/carbon steel": 0.00,
        "carbon steel/brass": 0.05,
        "carbon steel/stainless steel": 0.13,
        "carbon steel/monel": 0.13,
        "carbon steel/titanium": 0.16,
        "carbon steel/Cr-Mo steel": 0.05,
        "Cr-Mo steel/Cr-Mo steel": 0.07,
        "stainless steel/stainless steel": 0.07,
        "monel/monel": 0.08,
        "titanium/titanium": 0.06,
    }

    a = hx_material_factor_a_dic[Mat_factor]
    b = hx_material_factor_b_dic[Mat_factor]

    def hx_material_fact_rule(self):
        if Mat_factor == "carbon steel/carbon steel":
            return self.material_factor == 1
        else:
            return self.material_factor == a + (area / 100 / pyunits.ft**2) ** b

    self.hx_material_eqn = Constraint(rule=hx_material_fact_rule)

    # ------------------------------------------------------
    # Pressure factor calculation
    # doublecheck units (higher pressure fluid should be tube side)

    pressure = pyunits.convert(
        self.parent_block().tube.properties_in[0].pressure, to_units=pyunits.psi
    )

    # units must be in psig
    def hx_P_factor(self):
        # eq valid from 600 pisg to 3000 psig
        #    return self.pressure_factor == 0.8510 + 0.1292*(pressure/600)
        #                                + 0.0198*(pressure/600)**2
        # eq valid from 100 pisg to 2000 psig
        return self.pressure_factor == (
            0.9803
            + 0.0180 * (pressure / 100 / pyunits.psi)
            + 0.0017 * (pressure / 100 / pyunits.psi) ** 2
        )

    self.p_factor_eq = Constraint(rule=hx_P_factor)

    # ---------------------------------------------------------
    # purchase cost equation
    def hx_CP_rule(self):
        return self.purchase_cost == (
            self.pressure_factor
            * self.material_factor
            * self.L_factor
            * (self.parent_block().flowsheet().costing.CE_index / 500)
            * self.base_cost
        )

    self.cp_cost_eq = Constraint(rule=hx_CP_rule)


def pressure_changer_costing(
    self,
    Mat_factor="stain_steel",
    # applies for all (pump, compressor, fan, blower)
    mover_type="compressor",
    # fan, blower, compressor
    compressor_type="centrifugal",
    # only for compressor
    driver_mover_type="electrical_motor",
    # only for compressors
    pump_type="centrifugal",
    # centrifugal, external_gear, reciprocating
    pump_type_factor="1.4",
    # 1.1 to 1.4, 2.1 and 2.2
    # (needs to be wise-selected by user see table)
    pump_motor_type_factor="open",
    # centrifugal_backward, centrifugal_straight
    # vane_axial, tube_axial
    fan_type="centrifugal_backward",
    # select from table depends on fan's head
    fan_head_factor=1.45,
    # centrifugal and rotary
    blower_type="centrifugal",
):
    """
    Pressure changer costing method
    Source: Process and Product Design Principles: Synthesis, Analysis,
    and Evaluation Seider, Seader, Lewin, Windagdo, 3rd Ed. John Wiley and Sons
    Chapter 22. Cost Accounting and Capital Cost Estimation
    22.2 Cost Indexes and Capital Investment

    This method computes the purchase cost (CP) for a pressure changer unit,
    and can be used for costing of pumps, fan, blower, compressor, or Turbine.

    (created 2/21/2020)

    Arguments:
        mover_type : Only if config.compressor is True. Valid values 'fan',
                        'blower', 'compressor' (default).
        compressor_type : Only if mover_type='compressor'. Valid values
                        'centrifugal' (default), 'reciprocating', 'screw'
        driver_mover_type : Only if mover_type='compressor'. Valid values
                        'electric_motor' (default), 'steam_turbine',
                        'gas_turbine'
        pump_type : Only if config.compressor is True and
                        compressor_type='pump'. Valid values 'centrifugal',
                        'external_gear', 'reciprocating'
        pump_type_factor : Only if config.compressor is True and
                        compressor_type='pump'. Valid values 1.1 to 1.4,
                        2.1 and 2.2 (default = 1.4)
        pump_motor_type_factor : Only if config.compressor is True and
                        compressor_type='pump'. Valid values 'open' (default),
                        'enclosed', 'explosion_proof''
        Mat_factor : construction material; Valid values 'stain_steel'
                        (default), 'nickel_alloy' (compressor only)

    Returns:
        None
    """
    # build generic costing variables
    # (base cost, purchase cost, material factor)
    _make_vars(self)

    self.material_factor = Param(
        mutable=True, initialize=3, doc="Construction material correction factor"
    )

    # checking units
    if self.parent_block().config.thermodynamic_assumption.name == "pump":
        w = self.parent_block().work_fluid[self.parent_block().flowsheet().time.first()]
    else:
        w = self.parent_block().work_mechanical[
            self.parent_block().flowsheet().time.first()
        ]

    work_hp = pyunits.convert(w, to_units=pyunits.hp) / self.number_of_units

    # if compressor is == False, that means pressure changer is a Turbine
    if self.parent_block().config.compressor is False:
        #                           -1*work_hp because work is negative
        def CP_rule(self):
            return (
                self.purchase_cost
                == self.number_of_units * 530 * (-1 / pyunits.hp * work_hp) ** 0.81
            )

        self.cp_cost_eq = Constraint(rule=CP_rule)

    # if compressor is True
    # Then the pressure changer could be a Pump, Fan, Blower, or Compressor
    else:
        # if ThermodynamicAssumption == Pump (the unit is a Pump)
        if self.parent_block().config.thermodynamic_assumption.name == "pump":
            # assign work fluid  = to w
            w = self.parent_block().work_fluid[
                self.parent_block().flowsheet().time.first()
            ]

            # new variables only used by pump costing
            self.pump_head = Var(
                initialize=10,
                domain=NonNegativeReals,
                doc="Pump Head in feet of fluid flowing (Pressure rise/dens)",
                units=pyunits.pound_force * pyunits.foot / pyunits.pound,
            )

            self.size_factor = Var(
                initialize=10000,
                domain=NonNegativeReals,
                doc="Pump size factor, f(Q,pump_head)",
            )

            self.motor_base_cost_per_unit = Var(
                initialize=10000,
                domain=NonNegativeReals,
                doc="Motor base purchase" " cost in $ per unit",
            )

            self.pump_purchase_cost = Var(
                initialize=100000,
                domain=NonNegativeReals,
                doc="Pump purchase cost in $",
            )

            self.motor_purchase_cost = Var(
                initialize=100000,
                domain=NonNegativeReals,
                doc="Motor purchase cost in $",
            )

            self.FT = Param(mutable=True, initialize=1, doc="Pump-type factor")

            # Pressure units required in lbf/ft^2
            deltaP_lb_ft2 = pyunits.convert(
                self.parent_block().deltaP[0],
                to_units=pyunits.psi * pyunits.inch**2 / pyunits.foot**2,
            )

            # dens mass units required in lb/ft^3
            dens_mass_lb_ft3 = pyunits.convert(
                self.parent_block().control_volume.properties_in[0].dens_mass,
                to_units=pyunits.pound / pyunits.foot**3,
            )

            # volumetric flow units required in gpm
            Q_gpm = (
                pyunits.convert(
                    self.parent_block().control_volume.properties_in[0].flow_vol,
                    to_units=pyunits.gallon / pyunits.minute,
                )
                / self.number_of_units
            )

            # build pump_head equation
            def p_head_rule(self):
                return self.pump_head == deltaP_lb_ft2 / dens_mass_lb_ft3

            self.p_head_eq = Constraint(rule=p_head_rule)

            # S  = Q(H)**0.5 (S = Size factor for pump)
            # Q = is the flow rate through the pump in gallons per minute
            # H = pump head in feet of flowing (pressure rise/liquid density)
            # build Size Factor equation
            def p_s_factor_rule(self):
                return self.size_factor == (
                    Q_gpm
                    * pyunits.minute
                    / pyunits.gallon
                    * (
                        self.pump_head
                        * pyunits.pound
                        / pyunits.pound_force
                        / pyunits.foot
                    )
                    ** 0.5
                )

            self.s_factor_eq = Constraint(rule=p_s_factor_rule)

            # Base cost and Purchase cost for centrifugal pump
            # material factor dictionary for
            #                       centrifugal pumps and external gear pumps
            material_factor_dic = {
                "cast_iron": 1.00,
                "ductile_iron": 1.15,
                "cast_steel": 1.35,
                "bronze": 1.90,
                "stain_steel": 2.00,
                "hastelloy_c": 2.95,
                "monel": 3.30,
                "nickel": 3.50,
                "titanium": 9.70,
            }
            # material factor dictionary for Reciprocating Plunger pumps
            material_factor_reciprocal_pumps = {
                "ductile_iron": 1.00,
                "Ni_Al_Bronze": 1.15,
                "carbon_steel": 1.50,
                "stain_steel": 2.20,
            }
            # pump type factor dictionary (only used by centrifugal pumps)
            pump_type_factor_dic = {
                "1.1": 1.00,
                "1.2": 1.50,
                "1.3": 1.70,
                "1.4": 2.00,
                "2.1": 2.70,
                "2.2": 8.90,
            }

            if pump_type == "centrifugal":
                self.material_factor = material_factor_dic[Mat_factor]
                self.FT = pump_type_factor_dic[pump_type_factor]

            elif pump_type == "external_gear":
                self.material_factor = material_factor_dic[Mat_factor]
                self.FT = 1

            elif pump_type == "reciprocating":
                self.material_factor = material_factor_reciprocal_pumps[Mat_factor]
                self.FT = 1
            else:
                # PYLINT-TODO-FIX fix exception message with correct number of arguments
                raise ValueError(
                    "{} - pump type not supported. "  # pylint: disable=E1305
                    "Please see documentation for "
                    "supported options.".format(self.name, pump_type)
                )

            # pump cost correlations ---------------------------------
            def base_pump_rule(self):
                if pump_type == "centrifugal":
                    return self.base_cost_per_unit == exp(
                        9.7171
                        - 0.6019 * log(self.size_factor)
                        + 0.0519 * log(self.size_factor) ** 2
                    )
                elif pump_type == "external_gear":
                    return self.base_cost_per_unit == exp(
                        7.6964 + 0.1986 * log(Q_gpm) + 0.0291 * log(Q_gpm) ** 2
                    )

                elif pump_type == "reciprocating":
                    # brake horsepower with efficiency np typically = 90%
                    PB = (Q_gpm * self.pump_head * dens_mass_lb_ft3 / 7.48052) / (
                        33000 * 0.90
                    )
                    return self.base_cost_per_unit == exp(
                        7.8103 + 0.26986 * log(PB) + 0.06718 * log(PB) ** 2
                    )
                else:
                    # PYLINT-TODO-FIX fix exception message with correct number of arguments
                    raise ValueError(
                        "{} - pump type not supported. "  # pylint: disable=E1305
                        "Please see documentation for "
                        "supported options.".format(self.name, pump_type)
                    )

            self.base_pump_cost_per_unit_eq = Constraint(rule=base_pump_rule)

            @self.Expression(doc="Base cost for all units installed")
            def base_cost(self):
                return self.base_cost_per_unit * self.number_of_units

            def CP_pump_rule(self):
                return self.pump_purchase_cost == (
                    self.FT
                    * self.material_factor
                    * (self.parent_block().flowsheet().costing.CE_index / 394)
                    * self.base_cost
                )

            self.cp_pump_cost_eq = Constraint(rule=CP_pump_rule)

            # electric motor cost correlations ------------------------------
            pump_motor_type_dic = {"open": 1, "enclosed": 1.4, "explosion_proof": 1.8}
            self.motor_FT = Param(
                mutable=True,
                initialize=pump_motor_type_dic[pump_motor_type_factor],
                doc="Motor type factor",
            )

            # pump fractional efficiency
            np = (
                -0.316
                + 0.24015 * log(Q_gpm * pyunits.minute / pyunits.gallon)
                - 0.01199 * log(Q_gpm * pyunits.minute / pyunits.gallon) ** 2
            )
            # fractional efficiency of the electric motor
            nm = (
                0.80
                + 0.0319 * log(work_hp / pyunits.hp)
                - 0.00182 * log(work_hp / pyunits.hp) ** 2
            )

            # power consumption in horsepower
            @self.Expression()
            def power_consumption_hp(self):
                return (Q_gpm * self.pump_head * dens_mass_lb_ft3 / 7.48052) / (
                    33000 * np * nm
                )

            def base_motor_cost_rule(self):
                pc_hp = self.power_consumption_hp / pyunits.get_units(
                    self.power_consumption_hp
                )
                return self.motor_base_cost_per_unit == exp(
                    5.8259
                    + 0.13141 * log(pc_hp)
                    + 0.053255 * log(pc_hp) ** 2
                    + 0.028628 * log(pc_hp) ** 3
                    - 0.0035549 * log(pc_hp) ** 4
                )

            self.base_motor_cost_eq = Constraint(rule=base_motor_cost_rule)

            @self.Expression(doc="Base cost for all units installed")
            def motor_base_cost(self):
                return self.motor_base_cost_per_unit * self.number_of_units

            def CP_motor_rule(self):
                return self.motor_purchase_cost == (
                    self.motor_FT
                    * self.parent_block().flowsheet().costing.CE_index
                    / 394
                    * self.motor_base_cost
                )

            self.cp_motor_cost_eq = Constraint(rule=CP_motor_rule)

            # Total pump cost (pump + electrical motor)
            def cp_cost_rule(self):
                return (
                    self.purchase_cost
                    == self.motor_purchase_cost + self.pump_purchase_cost
                )

            self.cp_cost_eq = Constraint(rule=cp_cost_rule)
        # ends pump costing code

        # compressor = True, and using isothermal assumption
        # (costing not needed)
        elif (self.parent_block().config.thermodynamic_assumption.name) == "isothermal":
            # PYLINT-TODO-FIX fix exception message with correct number of arguments
            raise ValueError(
                "{} - pressure changers with isothermal "  # pylint: disable=E1305
                "assumption are too simple to be costed. ".format(self.name, mover_type)
            )
        # if config.compressor is = True
        # if thermodynamic_assumption is not = Pump
        # (pressure changer could be a Fan, Blower, or Compressor)
        else:
            # The user has to select mover_type [compressor or Fan or Blower]
            if mover_type == "compressor":
                # Compressor Purchase Cost Correlation
                FD_param = {
                    "electrical_motor": 1,
                    "steam_turbine": 1.15,
                    "gas_turbine": 1.25,
                }
                self.FD = Param(
                    mutable=True,
                    initialize=FD_param[driver_mover_type],
                    doc="Mover drive factor",
                )
                # compressor material factor dictionary
                material_factor_dic = {
                    "carbon_steel": 1,
                    "stain_steel": 2.5,
                    "nickel_alloy": 5.0,
                }

                self.material_factor = material_factor_dic[Mat_factor]

                c_alf1 = {"centrifugal": 7.58, "reciprocating": 7.9661, "screw": 8.1238}

                c_alf2 = {"centrifugal": 0.8, "reciprocating": 0.8, "screw": 0.7243}

                # Purchase cost rule
                def CB_rule(self):
                    return self.base_cost_per_unit == exp(
                        c_alf1[compressor_type]
                        + c_alf2[compressor_type] * log(work_hp / pyunits.hp)
                    )

                self.base_cost_per_unit_eq = Constraint(rule=CB_rule)

                @self.Expression(doc="Base cost for all units installed")
                def base_cost(self):
                    return self.base_cost_per_unit * self.number_of_units

                def CP_rule(self):
                    return self.purchase_cost == (
                        self.FD
                        * self.material_factor
                        * (self.parent_block().flowsheet().costing.CE_index / 500)
                        * self.base_cost
                    )

                self.cp_cost_eq = Constraint(rule=CP_rule)

            # FAN Costing Model --------------------------------------------
            elif mover_type == "fan":

                # volumetric flow units required in cfm
                if (
                    self.parent_block()
                    .config.property_package.get_metadata()
                    .properties["flow_vol"]["units"]
                ) == "m^3/s":
                    # 1 m3/s = 2118.88 cfm
                    Q_cfm = (
                        self.parent_block().control_volume.properties_in[0].flow_vol
                        * 2118.88
                    )
                elif (
                    self.parent_block()
                    .config.property_package.get_metadata()
                    .properties["flow_vol"]["units"]
                ) == "m^3/hr":
                    Q_cfm = (
                        self.parent_block().control_volume.properties_in[0].flow_vol
                        * 0.589
                    )  # 1 m3/hr = 0.589cfm
                elif (
                    self.parent_block()
                    .config.property_package.get_metadata()
                    .properties["flow_vol"]["units"]
                ) == "ft^3/s":
                    Q_cfm = (
                        self.parent_block().control_volume.properties_in[0].flow_vol
                        * 60
                    )  # 1ft3/s*60s/min = ft3/m
                # end volumetric flow units
                else:
                    # PYLINT-TODO-FIX fix exception message with correct number of arguments
                    raise ValueError(
                        "{} - volumetric flowrate units "  # pylint: disable=E1305
                        "not supported. "
                        "Please see documentation for list of "
                        "supported units.".format(self.name, mover_type)
                    )

                # fan cost correlation
                c_alf1 = {
                    "centrifugal_backward": 11.0757,
                    "centrifugal_straight": 12.1678,
                    "vane_axial": 9.5229,
                    "tube_axial": 6.12905,
                }

                c_alf2 = {
                    "centrifugal_backward": 1.12906,
                    "centrifugal_straight": 1.31363,
                    "vane_axial": 0.97566,
                    "tube_axial": 0.40254,
                }

                c_alf3 = {
                    "centrifugal_backward": 0.08860,
                    "centrifugal_straight": 0.09974,
                    "vane_axial": 0.08532,
                    "tube_axial": 0.05787,
                }

                self.head_factor = fan_head_factor

                material_factor_dic = {
                    "carbon_steel": 1.0,
                    "fiberglass": 1.8,
                    "stain_steel": 2.5,
                    "nickel_alloy": 5.0,
                }

                self.material_factor = material_factor_dic[Mat_factor]

                # Base cost
                def CB_rule(self):
                    return self.base_cost_per_unit == exp(
                        c_alf1[fan_type]
                        - c_alf2[fan_type] * (log(Q_cfm))
                        + c_alf3[fan_type] * (log(Q_cfm) ** 2)
                    )

                self.base_cost_per_unit_eq = Constraint(rule=CB_rule)

                @self.Expression(doc="Base cost for all units installed")
                def base_cost(self):
                    return self.base_cost_per_unit * self.number_of_units

                def CP_rule(self):
                    return self.purchase_cost == (
                        self.material_factor
                        * self.head_factor
                        * (self.parent_block().flowsheet().costing.CE_index / 500)
                        * self.base_cost
                    )

                self.cp_cost_eq = Constraint(rule=CP_rule)

            # Blower Costing -------------------------------------------------
            elif mover_type == "blower":
                # Blower Cost Correlation
                c_alf1 = {"centrifugal": 6.8929, "rotary": 7.59176}

                c_alf2 = {"centrifugal": 0.7900, "rotary": 0.7932}

                c_alf3 = {"centrifugal": 0.0, "rotary": 0.012900}

                material_factor_dic = {
                    "carbon_steel": 1.0,
                    "aluminum": 0.60,
                    "fiberglass": 1.8,
                    "stain_steel": 2.5,
                    "nickel_alloy": 5.0,
                }
                self.material_factor = material_factor_dic[Mat_factor]

                # Base cost
                def CB_rule(self):
                    return self.base_cost_per_unit == exp(
                        c_alf1[blower_type]
                        + c_alf2[blower_type] * (log(work_hp))
                        - c_alf3[blower_type] * (log(work_hp) ** 2)
                    )

                self.base_cost_per_unit_eq = Constraint(rule=CB_rule)

                @self.Expression(doc="Base cost for all units installed")
                def base_cost(self):
                    return self.base_cost_per_unit * self.number_of_units

                def CP_rule(self):
                    return self.purchase_cost == (
                        self.material_factor
                        * (self.parent_block().flowsheet().costing.CE_index / 500)
                        * self.base_cost
                    )

                self.cp_cost_eq = Constraint(rule=CP_rule)
            else:
                # PYLINT-TODO-FIX fix exception message with correct number of arguments
                raise ValueError(
                    "{} - mover type not supported. "  # pylint: disable=E1305
                    "Please see documentation for list of "
                    "supported mover types.".format(self.name, mover_type)
                )


def initialize(self):
    """
    Section to initialize costing blocks for all unit models
    """
    try:
        # only heat exchanger and Fired heater calculate pressure factor
        calculate_variable_from_constraint(self.pressure_factor, self.p_factor_eq)
    except AttributeError:
        pass
    if hasattr(self, "hx_material_eqn"):
        # section to initialize heat exchanger costing block
        calculate_variable_from_constraint(self.material_factor, self.hx_material_eqn)

    elif hasattr(self, "pump_head"):
        # section to initialize pressure changer costing block
        # only for pump and electric motor costing initialization
        # all the other pressure changers only include purchase cost and/or
        # base_cost_per_unit (compressors, turbines, fans, and blowers)
        calculate_variable_from_constraint(self.pump_head, self.p_head_eq)

        calculate_variable_from_constraint(self.size_factor, self.s_factor_eq)

        calculate_variable_from_constraint(
            self.base_cost_per_unit, self.base_pump_cost_per_unit_eq
        )

        calculate_variable_from_constraint(
            self.pump_purchase_cost, self.cp_pump_cost_eq
        )

        calculate_variable_from_constraint(
            self.motor_base_cost_per_unit, self.base_motor_cost_eq
        )
        calculate_variable_from_constraint(
            self.motor_purchase_cost, self.cp_motor_cost_eq
        )

        calculate_variable_from_constraint(
            self.motor_purchase_cost, self.cp_motor_cost_eq
        )

        calculate_variable_from_constraint(self.purchase_cost, self.cp_cost_eq)

    elif hasattr(self, "vessel_purchase_cost"):
        # section to initialize vessel costing block
        calculate_variable_from_constraint(self.weight, self.weight_eq)

        calculate_variable_from_constraint(self.base_cost_per_unit, self.cv_cost_eq)

        if hasattr(self, "base_cost_platf_ladders "):
            calculate_variable_from_constraint(
                self.base_cost_platf_ladders, self.CPL_eq
            )

        if hasattr(self, "base_cost_trays"):
            calculate_variable_from_constraint(
                self.tray_material_factor, self.mt_factor_eq
            )
            calculate_variable_from_constraint(self.base_cost_trays, self.bc_tray_eq)
            calculate_variable_from_constraint(
                self.purchase_cost_trays, self.cp_tray_eq
            )
        # calculate CPL cost before vessel purchase cost eq
        calculate_variable_from_constraint(self.vessel_purchase_cost, self.cp_vessel_eq)

    try:
        # base_cost_units variable doesn't exist in turbine costing
        calculate_variable_from_constraint(
            self.base_cost_per_unit, self.base_cost_per_unit_eq
        )
    except AttributeError:
        pass

    # all costing methods use purchase cost
    calculate_variable_from_constraint(self.purchase_cost, self.cp_cost_eq)


def calculate_scaling_factors(self):
    try:
        # only turbine model doesn't include base cost
        # scaling base_cost_per_unit
        sf_a = iscale.get_scaling_factor(
            self.base_cost_per_unit, default=1e-4, warning=True
        )
        iscale.constraint_scaling_transform(
            self.base_cost_per_unit_eq, sf_a, overwrite=False
        )
    except AttributeError:
        pass
    # scaling purchase cost (both var/eqn exist in all costing methods)
    s_purchase_cost = iscale.get_scaling_factor(
        self.purchase_cost, default=1e-4, warning=True
    )

    iscale.constraint_scaling_transform(
        self.cp_cost_eq, s_purchase_cost, overwrite=False
    )
