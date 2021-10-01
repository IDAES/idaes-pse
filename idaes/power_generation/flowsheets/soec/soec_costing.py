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

__author__ = "Alex Noring"

import pyomo.environ as pyo
from pyomo.util.calc_var_value import calculate_variable_from_constraint

from idaes.power_generation.costing.power_plant_costing import (
    get_PP_costing,
    get_total_TPC,
    costing_initialization,
    get_fixed_OM_costs,
    get_variable_OM_costs,
    initialize_fixed_OM_costs,
    initialize_variable_OM_costs,
)

from idaes.core.util.unit_costing import initialize as cost_init
from idaes.core.util.unit_costing import fired_heater_costing


def add_total_plant_cost(b, project_contingency, process_contingency):

    # b.costing.total_plant_cost = pyo.Var(bounds=(0, 1e4))
    try:
        b.parent_block().generic_costing_units.append(b)
    except AttributeError:
        b.parent_block().generic_costing_units = []
        b.parent_block().generic_costing_units.append(b)

    @b.costing.Expression()
    def total_plant_cost(c):
        return c.purchase_cost * project_contingency * process_contingency / 1e6


def get_soec_capital_costing(m):

    m.fs.get_costing(year="2018")

    # soec stack cost - fixed price per kW
    m.fs.soec.costing = pyo.Block()
    m.fs.soec.costing.total_plant_cost = pyo.Var(
        initialize=130,
        # bounds=(0, 1e4),
        doc="total plant cost in $MM",
    )

    @m.fs.soec.costing.Constraint()
    def total_plant_cost_eq(c):
        soec_cost = 421  # $/kW
        return c.total_plant_cost * 1e6 == soec_cost * m.fs.soec.total_power[0] / 1000

    # air preheater, ng preheaters, and bhx1 - U-tube HXs
    # costed with IDAES generic heat exchanger correlation
    m.fs.air_preheater.get_costing(hx_type="U-tube")
    add_total_plant_cost(m.fs.air_preheater, 1.15, 1.15)

    m.fs.ng_preheater.get_costing(hx_type="U-tube")
    add_total_plant_cost(m.fs.ng_preheater, 1.15, 1.15)

    m.fs.bhx1.get_costing(hx_type="U-tube")
    add_total_plant_cost(m.fs.bhx1, 1.15, 1.15)

    # bhx2
    # costed with IDAES generic fired heater correlation
    m.fs.bhx2.costing = pyo.Block()
    fired_heater_costing(
        m.fs.bhx2.costing,
        fired_type="steam_boiler",
        ref_parameter_pressure=m.fs.bhx2.tube_inlet.pressure[0],
    )
    add_total_plant_cost(m.fs.bhx2, 1.15, 1.15)

    # H2 compresso
    # costed with IDAES generic compressor correlations
    # this is a placeholder until H2 compressor scaling is released
    for unit in [m.fs.hcmp01, m.fs.hcmp02, m.fs.hcmp03, m.fs.hcmp04]:
        unit.get_costing()
        add_total_plant_cost(unit, 1.15, 1.15)

    # all the following equipment is scaled based on the bit baseline report

    # hxo2 and hxh2 - HRSGs
    HRSG_accounts = ["7.1", "7.2"]

    # get process parameters in correct units
    MMBtu_per_hr = pyo.units.MBtu / pyo.units.hr
    hxo2_duty = pyo.units.convert(m.fs.hxo2.heat_duty[0], MMBtu_per_hr)
    hxh2_duty = pyo.units.convert(m.fs.hxh2.heat_duty[0], MMBtu_per_hr)

    get_PP_costing(m.fs.hxo2, HRSG_accounts, hxo2_duty, "MMBtu/hr", 6)
    get_PP_costing(m.fs.hxh2, HRSG_accounts, hxh2_duty, "MMBtu/hr", 6)

    # carbon capture system

    # accounts that scale on CO2 flowrate in lb/hr
    # 5.1.a is the portion of Cansolv CO2 removal system scaled on CO2 flow
    # 5.12 is gas cleanup foundations
    CO2_removal_accounts_A = ["5.1.a", "5.12"]
    m.fs.CO2_removal_A = pyo.Block()

    lb_per_hr = pyo.units.lb / pyo.units.hr
    CO2_molar_mass = 44.01 * pyo.units.g / pyo.units.mol
    CO2_flow = pyo.units.convert(
        (
            0.9
            * CO2_molar_mass
            * (  # assume a 90% capture rate
                m.fs.ng_preheater.shell_outlet.flow_mol[0]
                * m.fs.ng_preheater.shell_outlet.mole_frac_comp[0, "CO2"]
                + m.fs.air_preheater.shell_outlet.flow_mol[0]
                * m.fs.air_preheater.shell_outlet.mole_frac_comp[0, "CO2"]
            )
        ),
        lb_per_hr,
    )

    get_PP_costing(
        m.fs.CO2_removal_A, CO2_removal_accounts_A, CO2_flow, "lb/hr", 6, ccs="B"
    )

    # accounts that scale on CO2 flowrate in lb/hr
    # 5.1.b is the portion of Cansolv CO2 removal system scaled on acfm
    CO2_removal_accounts_B = ["5.1.b"]
    m.fs.CO2_removal_B = pyo.Block()

    acfm = pyo.units.ft ** 3 / pyo.units.min
    absorber_inlet_flow = pyo.units.convert(
        (
            m.fs.ng_preheater.shell.properties_out[0].flow_vol
            + m.fs.ng_preheater.tube.properties_out[0].flow_vol
        ),
        acfm,
    )

    get_PP_costing(
        m.fs.CO2_removal_B,
        CO2_removal_accounts_B,
        absorber_inlet_flow,
        "acfm",
        6,
        ccs="B",
    )

    # water treatment
    WT_accounts = ["3.2", "3.4"]
    m.fs.WT_block = pyo.Block()

    # scaled on water flowrate to soec
    gal_per_min = pyo.units.gallon / pyo.units.min
    raw_water_withdrawal = pyo.units.convert(
        m.fs.aux_boiler_feed_pump.control_volume.properties_in[0].flow_vol, gal_per_min
    )

    get_PP_costing(m.fs.WT_block, WT_accounts, raw_water_withdrawal, "gpm", 6)

    # accounts scaled based on aux load - accessory electrical equipment &
    # instrumentation and control
    auxilliary_load_accounts = [
        "11.2",
        "11.3",
        "11.4",
        "11.5",
        "11.6",
        "12.4",
        "12.5",
        "12.6",
        "12.7",
        "12.8",
        "12.9",
    ]
    m.fs.aux_load_block = pyo.Block()

    aux_load = (
        m.fs.soec.total_power[0] / 1000
        + m.fs.aux_boiler_feed_pump.work_mechanical[0] / 1000
    )

    get_PP_costing(m.fs.aux_load_block, auxilliary_load_accounts, aux_load, "kW", 6)

    # build constraint summing total plant costs
    get_total_TPC(m)

    # costing initialization
    calculate_variable_from_constraint(
        m.fs.soec.costing.total_plant_cost, m.fs.soec.costing.total_plant_cost_eq
    )

    """
    generic_costing_units = [
        m.fs.air_preheater,
        m.fs.ng_preheater,
        m.fs.bhx1,
        m.fs.bhx2,
        m.fs.hcmp01,
        m.fs.hcmp02,
        m.fs.hcmp03,
        m.fs.hcmp04
        ]

    for u in generic_costing_units:
        cost_init(u.costing)

        calculate_variable_from_constraint(
            u.costing.total_plant_cost,
            u.costing.total_plant_cost_eq)
    """
    costing_initialization(m.fs)
    calculate_variable_from_constraint(
        m.fs.costing.total_TPC, m.fs.costing.total_TPC_eq
    )


def lock_capital_cost(m):
    for b in m.fs.generic_costing_units:
        for v in b.costing.total_plant_cost.values():
            print(b)
            v.display()
            v.set_value(pyo.value(v))

        b.costing.deactivate()
    m.fs.costing.deactivate()
    m.fs.costing.total_TPC.fix()


def get_soec_OM_costing(m, design_h2_production=2.5 * pyo.units.kg / pyo.units.s):
    # fixed O&M costs
    get_fixed_OM_costs(m, design_h2_production, tech=6)

    @m.fs.Constraint()
    def stack_replacement_cost(fs):
        # stack replacement cost
        stack_replacement_cost = 24.29  # $/yr/kW
        return fs.costing.other_fixed_costs * 1e6 == (
            stack_replacement_cost * m.fs.soec.total_power[0] / 1000
        )

    m.fs.costing.other_fixed_costs.unfix()

    # initialize fixed costs
    # for t in m.fs.time:
    #    calculate_variable_from_constraint(m.fs.H2_product[t],
    #                                       m.fs.H2_product_rule[t])
    calculate_variable_from_constraint(
        m.fs.costing.other_fixed_costs, m.fs.stack_replacement_cost
    )
    initialize_fixed_OM_costs(m)

    # variable O&M costs

    # electricity
    @m.fs.Expression(m.fs.time)
    def plant_load(fs, t):
        return (
            m.fs.soec.total_power[t]
            + m.fs.h2_compressor_power[t]
            + m.fs.aux_boiler_feed_pump.work_mechanical[t]
        )

    # natural gas
    @m.fs.Expression(m.fs.time)
    def NG_energy_rate(fs, t):
        NG_HHV = 908839.23 * pyo.units.J / pyo.units.mol
        return m.fs.ng_preheater.tube_inlet.flow_mol[t] * NG_HHV

    # water
    @m.fs.Expression(m.fs.time)
    def water_withdrawal(fs, t):  # cm^3/s
        molar_mass = 18.015 * pyo.units.g / pyo.units.mol
        density = 0.997 * pyo.units.g / pyo.units.cm ** 3
        return (
            molar_mass
            / density
            * (
                m.fs.aux_boiler_feed_pump.inlet.flow_mol[t]
                - m.fs.hxo2.shell_outlet.flow_mol[t]
                * m.fs.hxo2.shell_outlet.mole_frac_comp[t, "H2O"]
                - m.fs.hxh2.shell_outlet.flow_mol[t]
                * m.fs.hxh2.shell_outlet.mole_frac_comp[t, "H2O"]
            )
        )

    # water treatment
    @m.fs.Expression(m.fs.time)
    def water_treatment_use(fs, t):
        use_rate = 0.00297886 * pyo.units.lb / pyo.units.gal
        return m.fs.water_withdrawal[t] * use_rate

    resources = ["electricity", "natural gas", "water", "water treatment chemicals"]
    rates = [
        m.fs.plant_load,
        m.fs.NG_energy_rate,
        m.fs.water_withdrawal,
        m.fs.water_treatment_use,
    ]
    prices = {"electricity": 30 * pyo.units.USD / pyo.units.MWh}

    print(pyo.units.convert(m.fs.h2_product_rate_mass[0], pyo.units.g / pyo.units.s))
    m.fs.h2_product_rate_mass.display()
    get_variable_OM_costs(m, m.fs.h2_product_rate_mass, resources, rates, prices=prices)

    # initialize variable costs
    initialize_variable_OM_costs(m)


def display_soec_costing(m):
    print("Capital cost: ${:.0f}M".format(pyo.value(m.fs.costing.total_TPC)))
    print(
        "Fixed O&M cost: ${:.1f}M/yr".format(
            pyo.value(m.fs.costing.total_fixed_OM_cost)
        )
    )
    print(
        "Electricity cost: ${:.2f}/kg H2".format(
            pyo.value(m.fs.H2_costing.variable_operating_costs[0, "electricity"])
        )
    )
    print(
        "Fuel cost: ${:.2f}/kg H2".format(
            pyo.value(m.fs.H2_costing.variable_operating_costs[0, "natural gas"])
        )
    )
    print(
        "Total variable O&M cost: ${:.2f}/kg H2".format(
            pyo.value(m.fs.H2_costing.total_variable_OM_cost[0])
        )
    )
