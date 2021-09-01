##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
"""
Cost evaluation using IDAES costing framework.

Author: A. Noring
"""

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


def get_soec_costing(m, evaluate_cost=False):

    m.fs.get_costing(year="2018")

    # soec stack cost - fixed price per kW
    m.fs.soec.costing = pyo.Block()
    m.fs.soec.costing.total_plant_cost = pyo.Var(initialize=30,
                                                 bounds=(0, 1e4),
                                                 doc='total plant cost in $MM')

    @m.fs.soec.costing.Constraint()
    def total_plant_cost_eq(c):
        soec_cost = 421  # $/kW
        return c.total_plant_cost*1e6 == soec_cost*m.fs.soec.power[0]/1000

    # hxa1 and hxf1 - costed with IDAES generic heat exchanger correlation
    m.fs.hxa1.get_costing()
    m.fs.hxf1.get_costing()

    # total plant cost is not included in generic costed, it needs to be built
    # so it can be added to total TPC
    m.fs.hxa1.costing.total_plant_cost = pyo.Var(bounds=(0, 1e4))
    m.fs.hxf1.costing.total_plant_cost = pyo.Var(bounds=(0, 1e4))

    @m.fs.hxa1.costing.Constraint()
    def total_plant_cost_eq(c):
        project_contingency = 1.15
        process_contingency = 1.15
        return c.total_plant_cost*1e6 == (c.purchase_cost *
                                          project_contingency *
                                          process_contingency)

    @m.fs.hxf1.costing.Constraint()
    def total_plant_cost_eq(c):
        project_contingency = 1.15
        process_contingency = 1.15
        return c.total_plant_cost*1e6 == (c.purchase_cost *
                                          project_contingency *
                                          process_contingency)

    # all other accounts are scaled based on the bit baseline report
    # hxa2 and hxf2 - incorporated into HRSG system
    HRSG_duty_accounts = ["7.1", "7.2"]
    m.fs.HRSG_block = pyo.Block()

    # scaled on HRSG duty
    MMbtu_per_hr = pyo.units.Mbtu/pyo.units.hr
    HRSG_duty = pyo.units.convert(
        (m.fs.hxa2.heat_duty[0] + m.fs.hxf2.heat_duty[0]),
        MMbtu_per_hr)

    if evaluate_cost:
        m.fs.HRSG_block.HRSG_duty = pyo.Var(initialize=pyo.value(HRSG_duty))
        m.fs.HRSG_block.HRSG_duty.fix()
        get_PP_costing(m.fs.HRSG_block, HRSG_duty_accounts,
                       m.fs.HRSG_block.HRSG_duty, "MMBtu/hr", 6)
    else:
        get_PP_costing(m.fs.HRSG_block, HRSG_duty_accounts,
                       HRSG_duty, "MMBtu/hr", 6)

    # water treatment
    WT_accounts = ["3.2"]
    m.fs.WT_block = pyo.Block()

    # scaled on water flowrate to soec
    gal_per_min = pyo.units.gallon/pyo.units.min
    raw_water_withdrawal = pyo.units.convert(
        (m.fs.hxa1.tube.properties_in[0].flow_vol +
         m.fs.hxf1.tube.properties_in[0].flow_vol),
        gal_per_min)

    if evaluate_cost:
        m.fs.WT_block.raw_water_withdrawal = pyo.Var(initialize=pyo.value(raw_water_withdrawal))
        m.fs.WT_block.raw_water_withdrawal.fix()
        get_PP_costing(m.fs.WT_block, WT_accounts,
                       m.fs.WT_block.raw_water_withdrawal, "gpm", 6)
    else:
        get_PP_costing(m.fs.WT_block, WT_accounts,
                       raw_water_withdrawal, "gpm", 6)

    # accounts scaled based on aux load - accessory electrical equipment &
    # instrumentation and control
    auxilliary_load_accounts = [
        "11.2",
        "11.3",
        "11.4",
        "11.5",
        "11.6",
        "12.1",
        "12.2",
        "12.3",
        "12.4",
        "12.5",
        "12.6",
        "12.7",
        "12.8",
        "12.9",
    ]
    m.fs.aux_load_block = pyo.Block()

    # in this version aux load is only the soec power
    aux_load = m.fs.soec.power[0]/1000

    if evaluate_cost:
        m.fs.aux_load_block.aux_load = pyo.Var(initialize=pyo.value(aux_load))
        m.fs.aux_load_block.aux_load.fix()
        get_PP_costing(m.fs.aux_load_block, auxilliary_load_accounts,
                       m.fs.aux_load_block.aux_load, "kW", 6)
    else:
        get_PP_costing(m.fs.aux_load_block, auxilliary_load_accounts,
                       aux_load, "kW", 6)

    # build constraint summing total plant costs
    get_total_TPC(m)

    # costing initialization
    calculate_variable_from_constraint(
        m.fs.soec.costing.total_plant_cost,
        m.fs.soec.costing.total_plant_cost_eq)
    m.fs.soec.costing.total_plant_cost.fix()

    cost_init(m.fs.hxa1.costing)
    cost_init(m.fs.hxf1.costing)

    calculate_variable_from_constraint(
        m.fs.hxa1.costing.total_plant_cost,
        m.fs.hxa1.costing.total_plant_cost_eq)

    calculate_variable_from_constraint(
        m.fs.hxf1.costing.total_plant_cost,
        m.fs.hxf1.costing.total_plant_cost_eq)

    costing_initialization(m.fs)

    calculate_variable_from_constraint(m.fs.costing.total_TPC,
                                       m.fs.costing.total_TPC_eq)


