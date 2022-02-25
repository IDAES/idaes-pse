###############################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University,
# West Virginia University Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
###############################################################################

"""
Costing methods for the soec and sofc operating modes of the reversible sofc

"""
__author__ = "Chinedu Okoli", "Alex Noring"

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

    b.costing.total_plant_cost = pyo.Var(initialize=1e-6,
                                         # bounds=(0, None)
                                         doc='Unit total_TPC in $MM'
                                         )
    try:
        b.parent_block().generic_costing_units.append(b)
    except AttributeError:
        b.parent_block().generic_costing_units = []
        b.parent_block().generic_costing_units.append(b)

    @b.costing.Constraint()
    def total_plant_cost_eq(c):
        return c.total_plant_cost == (c.purchase_cost * project_contingency *
                                      process_contingency / 1e6)


def get_rsofc_sofc_capital_cost(fs):

    fs.get_costing(year="2018")

    # TPC of 650 MW NG based sofc plant
    # SOEC-only cap costs for water-side equipment and H2 compr. will be added
    # NGFC_TPC = 693.019  # MM$

    # fs.costing = pyo.Block()  # This is created in the costing module
    get_total_TPC(fs)
    fs.costing.total_plant_cost = pyo.Var(
        initialize=1e-6,
        # bounds=(0, 1e4),
        doc="total plant cost in $MM",
    )

    @fs.costing.Constraint()
    def total_plant_cost_eq(c):
        NGFC_TPC = 693.019  # MM$
        return c.total_plant_cost == NGFC_TPC

    calculate_variable_from_constraint(
        fs.costing.total_plant_cost,
        fs.costing.total_plant_cost_eq)

    # TODO - Not sure why, but total_plant_cost is not assigned to total_TPC
    # TODO - ie total_TPC returns as zero
    get_total_TPC(fs)

    costing_initialization(fs)

    calculate_variable_from_constraint(
        fs.costing.total_TPC,
        fs.costing.total_TPC_eq
    )


def get_rsofc_soec_capital_cost(fs):

    fs.get_costing(year="2018")

    # TPC of 650 MW NG based sofc plant
    # SOEC-only cap costs for water-side equipment and H2 compr. will be added
    # For rsofc_soec mode - O2 and ng preheater costs are captured in NGFC_TPC
    # Water treatment and CO2 processing capcosts included in NGFC_TPC

    # bhx1 - U-tube HXs
    # costed with IDAES generic heat exchanger correlation

    fs.bhx1.get_costing(hx_type="U-tube")
    add_total_plant_cost(fs.bhx1, 1.15, 1.15)

    # bhx2
    # TODO - is this the right costing for oxycombustor with in-bed HEX tubes?
    # TODO - try costing bhx2 as hrsg and compare the numbers to the steam_boiler cost
    # costed with IDAES generic fired heater correlation
    fs.bhx2.costing = pyo.Block()
    fired_heater_costing(
        fs.bhx2.costing,
        fired_type="steam_boiler",
        ref_parameter_pressure=fs.bhx2.inlet.pressure[0],
    )
    add_total_plant_cost(fs.bhx2, 1.15, 1.15)

    # H2 compressor
    # costed with IDAES generic compressor correlations
    # TODO - this is a placeholder until H2 compressor scaling is released
    for unit in [fs.hcmp01,
                 fs.hcmp02,
                 fs.hcmp03,
                 fs.hcmp04]:
        unit.get_costing()
        add_total_plant_cost(unit, 1.15, 1.15)

    get_total_TPC(fs)

    fs.generic_costing_units = [
        fs.bhx1,
        fs.bhx2,
        fs.hcmp01,
        fs.hcmp02,
        fs.hcmp03,
        fs.hcmp04
        ]

    for u in fs.generic_costing_units:
        cost_init(u.costing)
        calculate_variable_from_constraint(
            u.costing.total_plant_cost,
            u.costing.total_plant_cost_eq)

    costing_initialization(fs)

    calculate_variable_from_constraint(
        fs.costing.total_TPC,
        fs.costing.total_TPC_eq
    )


def lock_rsofc_soec_capital_cost(fs):
    for b in fs.generic_costing_units:
        for v in b.costing.total_plant_cost.values():
            print("TPC of generic units")
            print(b)
            v.display()
            v.set_value(pyo.value(v))
        b.costing.deactivate()
    fs.costing.deactivate()
    fs.costing.total_TPC.fix()


def lock_rsofc_sofc_capital_cost(fs):
    fs.costing.deactivate()
    fs.costing.total_TPC.fix()


def get_rsofc_sofc_fixed_OM_costing(fs,
                                    design_rsofc_netpower=650 * pyo.units.MW,
                                    fixed_TPC=None):

    # TODO - deleting this Var/Constraint pair so non-zero TPC value is used
    # del fs.costing.total_TPC
    # del fs.costing.total_TPC_eq
    fs.costing.del_component(fs.costing.total_TPC)
    fs.costing.del_component(fs.costing.total_TPC_eq)

    # fixed O&M costs
    get_fixed_OM_costs(fs, design_rsofc_netpower, tech=6, fixed_TPC=fixed_TPC)
    # labor_rate=38.50
    # get_fixed_OM_costs(m, design_rsofc_netpower, labor_rate=38.50,
    #                    labor_burden=30, operators_per_shift=6, tech=6)

    @fs.costing.Constraint()
    def stack_replacement_cost(costing):
        # stack replacement cost
        stack_replacement_cost = 13.26  # $MM/yr
        return costing.other_fixed_costs == stack_replacement_cost

    fs.costing.other_fixed_costs.unfix()

    # initialize fixed costs
    calculate_variable_from_constraint(
        fs.costing.other_fixed_costs,
        fs.costing.stack_replacement_cost
    )
    initialize_fixed_OM_costs(fs)


def get_rsofc_soec_fixed_OM_costing(fs,
                                    design_rsofc_netpower=650 * pyo.units.MW,
                                    fixed_TPC=None):
    # fixed O&M costs
    get_fixed_OM_costs(fs, design_rsofc_netpower, tech=6)
    # labor_rate=38.50
    # get_fixed_OM_costs(m, design_rsofc_netpower, labor_rate=38.50,
    #                    labor_burden=30, operators_per_shift=6, tech=6)
    initialize_fixed_OM_costs(fs)


def get_rsofc_soec_variable_OM_costing(fs):
    # variable O&M costs
    # electricity
    @fs.Expression(fs.time)
    def plant_load(fs, t):
        return (
            fs.net_power[t]  # net power
        )

    # natural gas
    @fs.Expression(fs.time)
    def NG_energy_rate(fs, t):
        NG_HHV = 908839.23 * pyo.units.J / pyo.units.mol
        return fs.ng_preheater.tube_inlet.flow_mol[t] * NG_HHV

    # water
    @fs.Expression(fs.time)
    def water_withdrawal(fs, t):  # cm^3/s
        molar_mass = 18.015 * pyo.units.g / pyo.units.mol
        density = 0.997 * pyo.units.g / pyo.units.cm ** 3
        return (
            molar_mass
            / density
            * (
                fs.aux_boiler_feed_pump.inlet.flow_mol[t]
                - fs.bhx1.shell_outlet.flow_mol[t]
                * fs.bhx1.shell_outlet.mole_frac_comp[t, "H2O"]
            )
        )

    # water treatment
    @fs.Expression(fs.time)
    def water_treatment_use(fs, t):
        use_rate = 0.00297886 * pyo.units.lb / pyo.units.gal
        return fs.water_withdrawal[t] * use_rate

    resources = ["electricity", "natural gas",
                 "water", "water treatment chemicals"]
    rates = [
        fs.plant_load,
        fs.NG_energy_rate,
        fs.water_withdrawal,
        fs.water_treatment_use,
    ]
    prices = {"electricity": 30 * pyo.units.USD / pyo.units.MWh}

    print(pyo.units.convert(fs.h2_product_rate_mass[0],
                            pyo.units.g / pyo.units.s))
    fs.h2_product_rate_mass.display()
    get_variable_OM_costs(fs, fs.h2_product_rate_mass,
                          resources, rates, prices=prices)

    # initialize variable costs
    initialize_variable_OM_costs(fs)


def get_rsofc_sofc_variable_OM_costing(fs):
    # variable O&M costs
    # Natural Gas Fuel
    NG_HHV = 908839.23*pyo.units.J/pyo.units.mol

    @fs.Expression(fs.time)
    def NG_energy_rate(fs, t):
        return fs.anode_mix.feed_inlet.flow_mol[t]*NG_HHV

    # Water
    fs.condensate_flowrate = pyo.Var(fs.time, initialize=1e-3,
                                     units=pyo.units.lb/pyo.units.s)

    # scaled based on BB Case 31A steam turbine power
    @fs.Constraint(fs.time)
    def condensate_flowrate_constraint(fs, t):
        ref_flowrate = 388.586*pyo.units.lb/pyo.units.s
        ref_power = 262.8*pyo.units.MW
        return (fs.condensate_flowrate[t] == ref_flowrate *
                fs.steam_cycle_power[t]/ref_power)

    fs.BFW_circulation_rate = pyo.Var(fs.time, initialize=1e-3,
                                      units=pyo.units.lb/pyo.units.s)

    @fs.Constraint(fs.time)
    def BFW_circulation_rate_constraint(fs, t):
        return fs.BFW_circulation_rate[t] == fs.condensate_flowrate[t]

    fs.BFW_makeup_flow = pyo.Var(fs.time, initialize=1e-3,
                                 units=pyo.units.lb/pyo.units.s)

    @fs.Constraint(fs.time)
    def BFW_makeup_flow_constraint(fs, t):
        return fs.BFW_makeup_flow[t] == 0.01*fs.BFW_circulation_rate[t]

    fs.evaporation_losses = pyo.Var(fs.time, initialize=1e-3,
                                    units=pyo.units.lb/pyo.units.s)

    @fs.Constraint(fs.time)
    def evaporation_losses_constraint(fs, t):
        return fs.evaporation_losses[t] == 0.016*fs.cooling_water_flowrate[t]

    fs.blowdown_losses = pyo.Var(fs.time, initialize=1e-3,
                                 units=pyo.units.lb/pyo.units.s)

    @fs.Constraint(fs.time)
    def blowdown_losses_constraint(fs, t):
        cycles_of_concentration = 4
        return (fs.blowdown_losses[t] ==
                fs.evaporation_losses[t]/(cycles_of_concentration - 1))

    fs.total_wet_losses = pyo.Var(fs.time, initialize=1e-3,
                                  units=pyo.units.lb/pyo.units.s)

    @fs.Constraint(fs.time)
    def wet_losses_constraint(fs, t):
        return fs.total_wet_losses[t] == (fs.evaporation_losses[t] +
                                          fs.blowdown_losses[t])

    fs.water_demand = pyo.Var(fs.time, initialize=1e-3,
                              units=pyo.units.lb/pyo.units.s)

    @fs.Constraint(fs.time)
    def water_demand_constraint(fs, t):
        return fs.water_demand[t] == (fs.BFW_makeup_flow[t] +
                                      fs.total_wet_losses[t])

    fs.water_recycle = pyo.Var(fs.time, initialize=1e-3,
                               units=pyo.units.lb/pyo.units.s)

    @fs.Constraint(fs.time)
    def water_recycle_constraint(fs, t):
        molar_mass = 18*pyo.units.g/pyo.units.mol
        CPU_water = pyo.units.convert(
            fs.CPU.water.flow_mol[t]*molar_mass,
            pyo.units.lb/pyo.units.s)
        flash_water = pyo.units.convert(
            fs.flash.liq_outlet.flow_mol[t]*molar_mass,
            pyo.units.lb/pyo.units.s)
        return fs.water_recycle[t] == CPU_water + flash_water

    fs.raw_water_withdrawal = pyo.Var(fs.time, initialize=1e-3,
                                      units=pyo.units.lb/pyo.units.s)

    @fs.Constraint(fs.time)
    def raw_water_withdrawal_constraint(fs, t):
        return fs.raw_water_withdrawal[t] == (fs.water_demand[t] -
                                              fs.water_recycle[t])

    fs.process_water_discharge = pyo.Var(fs.time, initialize=9,
                                         units=pyo.units.lb/pyo.units.s)

    @fs.Constraint(fs.time)
    def water_discharge_constraint(fs, t):
        return fs.process_water_discharge[t] == 0.9*fs.blowdown_losses[t]

    fs.raw_water_consumption = pyo.Var(fs.time, initialize=1e-3,
                                       units=pyo.units.lb/pyo.units.s)

    @fs.Constraint(fs.time)
    def raw_water_consumption_constraint(fs, t):
        return fs.raw_water_consumption[t] == (fs.raw_water_withdrawal[t] -
                                               fs.process_water_discharge[t])

    # MU & WT Chem
    @fs.Expression(fs.time)
    def waste_treatment_use(fs, t):
        use_rate = 0.00297886*pyo.units.lb/pyo.units.gal
        density = 8.34*pyo.units.lb/pyo.units.gal
        return fs.water_demand[t]/density*use_rate

    # NG desulfur TDA Adsorbent
    @fs.Expression(fs.time)
    def desulfur_adsorbent_use(fs, t):
        NG_sulfur_ppm = 5/1e6
        sulfur_MW = 32*pyo.units.g/pyo.units.mol
        sulfur_capacity = 0.03
        return (fs.anode_mix.feed_inlet.flow_mol[t] *
                NG_sulfur_ppm * sulfur_MW / sulfur_capacity)

    # Methanation Catalyst
    @fs.Expression(fs.time)
    def methanation_catalyst_use(fs, t):
        syngas_flow = pyo.units.convert(
            fs.anode_mix.feed_inlet_state[0].flow_mass,
            pyo.units.lb/pyo.units.hr)
        return (3.2*syngas_flow/611167 *
                pyo.units.m**3/pyo.units.day*pyo.units.hr/pyo.units.lb)

    resources = ["natural gas",
                 "water treatment chemicals",
                 "desulfur adsorbent",
                 "methanation catalyst"]

    rates = [fs.NG_energy_rate,
             fs.waste_treatment_use,
             fs.desulfur_adsorbent_use,
             fs.methanation_catalyst_use]

    prices = {"desulfur adsorbent": 6.0297*pyo.units.USD/pyo.units.lb,
              "methanation catalyst": 601.765*pyo.units.USD/pyo.units.m**3}

    get_variable_OM_costs(fs, fs.net_power, resources, rates, prices)

    # initialize variable costs
    initialize_variable_OM_costs(fs)


def display_rsofc_costing(m):
    print("Capital cost: ${:.0f}M".format(pyo.value(m.fs.costing.total_TPC)))
    print(
        "Fixed O&M cost: ${:.1f}M/yr".format(
            pyo.value(m.fs.costing.total_fixed_OM_cost)
        )
    )
    print(
        "Electricity cost: ${:.2f}/kg H2".format(
            pyo.value(
                m.soec_fs.H2_costing.variable_operating_costs[0,
                                                              "electricity"])
        )
    )
    print(
        "Fuel cost: ${:.2f}/kg H2".format(
            pyo.value(
                m.soec_fs.H2_costing.variable_operating_costs[0,
                                                              "natural gas"])
        )
    )
    print(
        "Total variable O&M cost: ${:.2f}/kg H2".format(
            pyo.value(m.soec_fs.H2_costing.total_variable_OM_cost[0])
        )
    )
