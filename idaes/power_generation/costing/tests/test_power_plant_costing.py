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

__author__ = "Costing Team (A. Noring and M. Zamarripa)"
__version__ = "1.0.0"

import pytest

from idaes.power_generation.costing.power_plant_costing import \
    (get_sCO2_unit_cost,
     get_PP_costing,
     get_ASU_cost,
     build_flowsheet_cost_constraint,
     costing_initialization)
from idaes.core.util.model_statistics import (degrees_of_freedom)
import pyomo.environ as pyo
from idaes.generic_models.properties import iapws95
from idaes.generic_models.properties import swco2
from idaes.core import FlowsheetBlock
from idaes.core.util import get_solver


@pytest.mark.component
def test_PP_costing():
    # Create a Concrete Model as the top level object
    m = pyo.ConcreteModel()

    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.get_costing(year='2018')

    # check that the model solved properly and has 0 degrees of freedom
    assert(degrees_of_freedom(m) == 0)

    ###########################################################################
    #  Create costing constraints                                             #
    ###########################################################################

    # coal flow rate
    # accounts 1.x and 2.x are coal handling, preparation and feed
    # accounts 4.x are for boiler BOP and foundations
    coal_accounts = ['1.1', '1.2', '1.3', '1.4',
                     '2.1', '2.2', '4.11', '4.15', '4.16']
    m.fs.boiler = pyo.Block()
    m.fs.boiler.coal_mass_flow = pyo.Var(initialize=7238.95)  # tpd
    m.fs.boiler.coal_mass_flow.fix()
    get_PP_costing(m.fs.boiler, coal_accounts, m.fs.boiler.coal_mass_flow,
                   'tpd', 2)

    # total fuel feed
    # accounts 3.x are for start up systems and miscellaneous plant equipment
    # accounts 7.x are for ductwork and stack foundations
    fuel_accounts = ['3.6', '3.9', '7.3', '7.5']
    m.fs.fuel_feed = pyo.Block()
    m.fs.fuel_feed.total_fuel_feed = pyo.Var(initialize=603246)  # lb/hr
    m.fs.fuel_feed.total_fuel_feed.fix()
    get_PP_costing(m.fs.fuel_feed, fuel_accounts,
                   m.fs.fuel_feed.total_fuel_feed, 'lb/hr', 2)

    # HP BFW flow rate
    # accounts 3.x are for feedwater systems
    # account 4.9 is for the boiler
    # account 8.4 is steam piping
    BFW_accounts = ['3.1', '3.3', '3.5', '4.9', '8.4']
    m.fs.bfp = pyo.Block()
    m.fs.bfp.BFW_mass_flow = pyo.Var(initialize=5316158)  # lb/hr
    m.fs.bfp.BFW_mass_flow.fix()
    get_PP_costing(m.fs.bfp, BFW_accounts, m.fs.bfp.BFW_mass_flow, 'lb/hr', 2)

    # Steam turbine power
    # accounts 8.x are for the steam turbine and its foundations
    power_accounts = ['8.1']
    m.fs.turb = pyo.Block()
    m.fs.turb.power = pyo.Var(initialize=769600)  # kW
    m.fs.turb.power.fix()
    get_PP_costing(m.fs.turb, power_accounts, m.fs.turb.power, 'kW', 2)

    # Condernser duty
    cond_accounts = ['8.3']
    m.fs.condenser = pyo.Block()
    m.fs.condenser.duty_MMBtu = pyo.Var(initialize=2016)  # MMBtu/hr
    m.fs.condenser.duty_MMBtu.fix()
    get_PP_costing(m.fs.condenser, cond_accounts, m.fs.condenser.duty_MMBtu,
                   "MMBtu/hr", 2)

    # Circulating water flow rate
    # accounts 9.x are for circulating water systems
    # account 14.5 is for the pumphouse
    circ_accounts = ['9.2', '9.3', '9.4', '9.6', '9.7', '14.5']
    m.fs.circulating_water = pyo.Block()
    m.fs.circulating_water.vol_flow = pyo.Var(initialize=463371)  # gpm
    m.fs.circulating_water.vol_flow.fix()
    get_PP_costing(m.fs.circulating_water,
                   circ_accounts, m.fs.circulating_water.vol_flow, 'gpm', 2)

    # Ash flow rate
    # accounts are for ash storage and handling
    ash_accounts = ['10.6', '10.7', '10.9']
    m.fs.ash_handling = pyo.Block()
    m.fs.ash_handling.ash_mass_flow = pyo.Var(initialize=66903)  # lb/hr
    m.fs.ash_handling.ash_mass_flow.fix()
    get_PP_costing(m.fs.ash_handling, ash_accounts,
                   m.fs.ash_handling.ash_mass_flow, 'lb/hr', 2)

    # add total cost
    build_flowsheet_cost_constraint(m)

    # add initialize
    costing_initialization(m.fs)

    # try solving
    solver = get_solver()
    results = solver.solve(m, tee=True)

    assert results.solver.termination_condition == \
        pyo.TerminationCondition.optimal

    #  all numbers come from the NETL excel file:
    # "201.001.001_BBR4 COE Spreadsheet_Rev0U_20190919_njk.xlsm"
    assert pytest.approx(pyo.value(
        m.fs.boiler.costing.total_plant_cost['1.1']), abs=1e-1) == 2306/1e3
    assert pytest.approx(pyo.value(
        m.fs.boiler.costing.total_plant_cost['1.2']), abs=1e-1) == 6385/1e3
    assert pytest.approx(pyo.value(
        m.fs.boiler.costing.total_plant_cost['1.3']), abs=1e-1) == 59527/1e3
    assert pytest.approx(pyo.value(
        m.fs.boiler.costing.total_plant_cost['1.4']), abs=1e-1) == 8086/1e3
    assert pytest.approx(pyo.value(
        m.fs.boiler.costing.total_plant_cost['2.1']), abs=1e-1) == 4073/1e3
    assert pytest.approx(pyo.value(
        m.fs.boiler.costing.total_plant_cost['2.2']), abs=1e-1) == 13976/1e3
    assert pytest.approx(pyo.value(
        m.fs.boiler.costing.total_plant_cost['4.11']), abs=1e-1) == 3751/1e3
    assert pytest.approx(pyo.value(
        m.fs.boiler.costing.total_plant_cost['4.15']), abs=1e-1) == 197/1e3
    assert pytest.approx(pyo.value(
        m.fs.boiler.costing.total_plant_cost['4.16']), abs=1e-1) == 1014/1e3

    assert pytest.approx(pyo.value(
        m.fs.fuel_feed.costing.total_plant_cost['3.6']), abs=1e-1) == 4864/1e3
    assert pytest.approx(pyo.value(
        m.fs.fuel_feed.costing.total_plant_cost['3.9']), abs=1e-1) == 522/1e3
    assert pytest.approx(pyo.value(
        m.fs.fuel_feed.costing.total_plant_cost['7.3']), abs=1e-1) == 1710/1e3
    assert pytest.approx(pyo.value(
        m.fs.fuel_feed.costing.total_plant_cost['7.5']), abs=1e-1) == 647/1e3

    assert pytest.approx(pyo.value(
        m.fs.bfp.costing.total_plant_cost['3.1']), abs=1e-1) == 19233/1e3
    assert pytest.approx(pyo.value(
        m.fs.bfp.costing.total_plant_cost['3.3']), abs=1e-1) == 6897/1e3
    assert pytest.approx(pyo.value(
        m.fs.bfp.costing.total_plant_cost['3.5']), abs=1e-1) == 2366/1e3
    assert pytest.approx(pyo.value(
        m.fs.bfp.costing.total_plant_cost['4.9']), abs=1e-1) == 570418/1e3
    assert pytest.approx(pyo.value(
        m.fs.bfp.costing.total_plant_cost['8.4']), abs=1e-1) == 81916/1e3

    assert pytest.approx(pyo.value(
        m.fs.turb.costing.total_plant_cost['8.1']), abs=1e-1) == 110166/1e3

    assert pytest.approx(pyo.value(
        m.fs.condenser.costing.total_plant_cost['8.3']), abs=1e-1) == 20447/1e3

    assert pytest.approx(pyo.value(
        m.fs.circulating_water.costing.total_plant_cost['9.2']),
        abs=1e-1) == 4133/1e3
    assert pytest.approx(pyo.value(
        m.fs.circulating_water.costing.total_plant_cost['9.3']),
        abs=1e-1) == 25518/1e3
    assert pytest.approx(pyo.value(
        m.fs.circulating_water.costing.total_plant_cost['9.4']),
        abs=1e-1) == 19859/1e3
    assert pytest.approx(pyo.value(
        m.fs.circulating_water.costing.total_plant_cost['9.6']),
        abs=1e-1) == 2870/1e3
    assert pytest.approx(pyo.value(
        m.fs.circulating_water.costing.total_plant_cost['9.7']),
        abs=1e-1) == 2690/1e3
    assert pytest.approx(pyo.value(
        m.fs.circulating_water.costing.total_plant_cost['14.5']),
        abs=1e-1) == 464/1e3

    assert pytest.approx(pyo.value(
        m.fs.ash_handling.costing.total_plant_cost['10.6']),
        abs=1e-1) == 6429/1e3
    assert pytest.approx(pyo.value(
        m.fs.ash_handling.costing.total_plant_cost['10.7']),
        abs=1e-1) == 10725/1e3
    assert pytest.approx(pyo.value(
        m.fs.ash_handling.costing.total_plant_cost['10.9']),
        abs=1e-1) == 2564/1e3

    assert pytest.approx(pyo.value(
        m.fs.flowsheet_cost), abs=1e-1) == 993753/1e3

    return m


@pytest.mark.component
def test_power_plant_costing():
    # Create a Concrete Model as the top level object
    m = pyo.ConcreteModel()

    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.get_costing(year='2018')

    ###########################################################################
    #  Create costing constraints                                             #
    ###########################################################################

    # subcritical PC
    coal_accounts = ['1.1', '1.2', '1.3']
    m.fs.subcritical_PC = pyo.Block()
    m.fs.subcritical_PC.coal_feed_rate = pyo.Var(initialize=7613.37)  # tpd
    m.fs.subcritical_PC.coal_feed_rate.fix()
    get_PP_costing(m.fs.subcritical_PC, coal_accounts,
                   m.fs.subcritical_PC.coal_feed_rate,
                   'tpd', 1)

    # two-stage, slurry-feed IGCC
    feedwater_accounts = ['3.1', '3.3', '3.5']
    m.fs.IGCC_1 = pyo.Block()
    m.fs.IGCC_1.feedwater_flow_rate = pyo.Var(initialize=1576062.15)  # lb/hr
    m.fs.IGCC_1.feedwater_flow_rate.fix()
    get_PP_costing(m.fs.IGCC_1, feedwater_accounts,
                   m.fs.IGCC_1.feedwater_flow_rate,
                   'lb/hr', 3)

    # single-stage, slurry-feed, IGCC
    syngas_accounts = ['6.1', '6.2', '6.3']
    m.fs.IGCC_2 = pyo.Block()
    m.fs.IGCC_2.syngas_flow_rate = pyo.Var(initialize=182335.921)  # lb/hr
    m.fs.IGCC_2.syngas_flow_rate.fix()
    get_PP_costing(m.fs.IGCC_2, syngas_accounts,
                   m.fs.IGCC_2.syngas_flow_rate,
                   'lb/hr', 4)

    # single-stage, dry-feed, IGCC
    HRSG_accounts = ['7.1', '7.2']
    m.fs.IGCC_3 = pyo.Block()
    m.fs.IGCC_3.HRSG_duty = pyo.Var(initialize=1777.86)  # MMBtu/hr
    m.fs.IGCC_3.HRSG_duty.fix()
    get_PP_costing(m.fs.IGCC_3, HRSG_accounts,
                   m.fs.IGCC_3.HRSG_duty,
                   'MMBtu/hr', 5)

    # NGCC
    steam_turbine_accounts = ['8.1', '8.2', '8.5']
    m.fs.NGCC = pyo.Block()
    m.fs.NGCC.turbine_power = pyo.Var(initialize=212500)  # kW
    m.fs.NGCC.turbine_power.fix()
    get_PP_costing(m.fs.NGCC, steam_turbine_accounts,
                   m.fs.NGCC.turbine_power,
                   'kW', 6)

    # AUSC PC
    AUSC_accounts = ['4.9', '8.4']
    m.fs.AUSC = pyo.Block()
    m.fs.AUSC.feedwater_flow = pyo.Var(initialize=3298815.58)  # lb/hr
    m.fs.AUSC.feedwater_flow.fix()
    get_PP_costing(m.fs.AUSC, AUSC_accounts,
                   m.fs.AUSC.feedwater_flow,
                   'lb/hr', 7)

    # add total cost
    build_flowsheet_cost_constraint(m)

    # add initialize
    costing_initialization(m.fs)

    # try solving
    solver = get_solver()
    results = solver.solve(m, tee=True)

    assert results.solver.termination_condition == \
        pyo.TerminationCondition.optimal
    #  all numbers come from the NETL excel file
    # "201.001.001_BBR4 COE Spreadsheet_Rev0U_20190919_njk.xlsm"
    assert pytest.approx(pyo.value(
        m.fs.subcritical_PC.costing.total_plant_cost['1.1']),
        abs=1e-1) == 2379/1e3
    assert pytest.approx(pyo.value(
        m.fs.subcritical_PC.costing.total_plant_cost['1.2']),
        abs=1e-1) == 6588/1e3
    assert pytest.approx(pyo.value(
        m.fs.subcritical_PC.costing.total_plant_cost['1.3']),
        abs=1e-1) == 61409/1e3

    assert pytest.approx(pyo.value(
        m.fs.IGCC_1.costing.total_plant_cost['3.1']), abs=1e-1) == 10807/1e3
    assert pytest.approx(pyo.value(
        m.fs.IGCC_1.costing.total_plant_cost['3.3']), abs=1e-1) == 2564/1e3
    assert pytest.approx(pyo.value(
        m.fs.IGCC_1.costing.total_plant_cost['3.5']), abs=1e-1) == 923/1e3

    assert pytest.approx(pyo.value(
        m.fs.IGCC_2.costing.total_plant_cost['6.1']), abs=1e-1) == 117850/1e3
    assert pytest.approx(pyo.value(
        m.fs.IGCC_2.costing.total_plant_cost['6.2']), abs=1e-1) == 3207/1e3
    assert pytest.approx(pyo.value(
        m.fs.IGCC_2.costing.total_plant_cost['6.3']), abs=1e-1) == 3770/1e3

    assert pytest.approx(pyo.value(
        m.fs.IGCC_3.costing.total_plant_cost['7.1']), abs=1e-1) == 53530/1e3
    assert pytest.approx(pyo.value(
        m.fs.IGCC_3.costing.total_plant_cost['7.2']), abs=1e-1) == 19113/1e3

    assert pytest.approx(pyo.value(
        m.fs.NGCC.costing.total_plant_cost['8.1']), abs=1e-1) == 49468/1e3
    assert pytest.approx(pyo.value(
        m.fs.NGCC.costing.total_plant_cost['8.2']), abs=1e-1) == 565/1e3
    assert pytest.approx(pyo.value(
        m.fs.NGCC.costing.total_plant_cost['8.5']), abs=1e-1) == 4094/1e3

    assert pytest.approx(pyo.value(
        m.fs.AUSC.costing.bare_erected_cost['4.9']), abs=1e-1) == 295509/1e3
    assert pytest.approx(pyo.value(
        m.fs.AUSC.costing.bare_erected_cost['8.4']), abs=1e-1) == 57265/1e3

    return m


@pytest.mark.component
def test_sCO2_costing():
    # Create a Concrete Model as the top level object
    m = pyo.ConcreteModel()

    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.get_costing(year='2017')

    # ######################################################
    # Primary Heater
    m.fs.boiler = pyo.Block()
    m.fs.boiler.heat_duty = pyo.Var(initialize=1461.5e6)
    m.fs.boiler.heat_duty.fix()
    m.fs.boiler.temp = pyo.Var(initialize=620)  # C
    m.fs.boiler.temp.fix()
    get_sCO2_unit_cost(m.fs.boiler, 'Coal-fired heater',
                       m.fs.boiler.heat_duty*(1e-6),
                       temp_C=m.fs.boiler.temp)

    # ######################################################
    # CO2 Turbine
    m.fs.turbine = pyo.Block()
    m.fs.turbine.work_isentropic = pyo.Var(initialize=1006.2e6)
    m.fs.turbine.work_isentropic.fix()
    m.fs.turbine.temp = pyo.Var(initialize=620)
    m.fs.turbine.temp.fix()
    get_sCO2_unit_cost(m.fs.turbine, 'Axial turbine',
                       m.fs.turbine.work_isentropic*(1e-6),
                       temp_C=m.fs.turbine.temp,
                       n_equip=1)

    # ######################################################
    # Generator
    m.fs.generator = pyo.Block()
    m.fs.generator.work_isentropic = pyo.Var(initialize=1006.2e6)
    get_sCO2_unit_cost(m.fs.generator, 'Generator',
                       m.fs.turbine.work_isentropic*(1e-6),
                       n_equip=1)

    # ######################################################
    # High Temperature Recuperator
    m.fs.HTR = pyo.Block()
    m.fs.HTR.heat_duty = pyo.Var(initialize=1461e6)  # W
    m.fs.HTR.heat_duty.fix()
    m.fs.HTR.LMTD = pyo.Var(initialize=21.45)
    m.fs.HTR.LMTD.fix()
    m.fs.HTR.temp = pyo.Var(initialize=453)
    m.fs.HTR.temp.fix()

    m.fs.HTR.UA = pyo.Var(initialize=1e8)

    # gives units of W/K
    @m.fs.Constraint()
    def HTR_UA_rule(b):
        return (b.HTR.UA*b.HTR.LMTD == b.HTR.heat_duty)

    get_sCO2_unit_cost(m.fs.HTR, 'Recuperator', m.fs.HTR.UA,
                       temp_C=m.fs.HTR.temp)

    # ######################################################
    # Low Temperature Recuperator
    m.fs.LTR = pyo.Block()
    m.fs.LTR.heat_duty = pyo.Var(initialize=911.7e6)  # W
    m.fs.LTR.heat_duty.fix()
    m.fs.LTR.LMTD = pyo.Var(initialize=5.21)
    m.fs.LTR.LMTD.fix()
    m.fs.LTR.temp = pyo.Var(initialize=216)
    m.fs.LTR.temp.fix()
    m.fs.LTR.UA = pyo.Var(initialize=1e8)

    @m.fs.Constraint()
    def LTR_UA_rule(b):
        return (b.LTR.UA*b.LTR.LMTD == b.LTR.heat_duty)

    get_sCO2_unit_cost(m.fs.LTR, 'Recuperator', m.fs.LTR.UA,
                       temp_C=m.fs.LTR.temp)

    # ######################################################
    # CO2 Cooler, costed using the recouperator not dry cooler
    m.fs.co2_cooler = pyo.Block()
    m.fs.co2_cooler.heat_duty = pyo.Var(initialize=739421217)
    m.fs.co2_cooler.heat_duty.fix()
    m.fs.co2_cooler.temp = pyo.Var(initialize=81)
    m.fs.co2_cooler.temp.fix()

    # Estimating LMTD
    # Cost from report: $27,780 thousand
    # Back-calculated UA: 41819213 W/K
    # Heat duty from report: 2523 MMBTu/hr --> 739421217 W
    # Estimated LMTD: 17.68 K
    m.fs.co2_cooler.LMTD = pyo.Var(initialize=5)
    m.fs.co2_cooler.UA = pyo.Var(initialize=1e5)
    m.fs.co2_cooler.LMTD.fix(17.68)

    @m.fs.Constraint()
    def co2_cooler_UA_rule(b):
        return (b.co2_cooler.UA * b.co2_cooler.LMTD ==
                b.co2_cooler.heat_duty)

    get_sCO2_unit_cost(m.fs.co2_cooler, 'Recuperator',
                       m.fs.co2_cooler.UA,
                       temp_C=m.fs.co2_cooler.temp)

    # ######################################################
    # Main Compressor - 5.99 m^3/s in Baseline620
    m.fs.main_compressor = pyo.Block()
    m.fs.main_compressor.flow_vol = pyo.Var(initialize=5.99)
    m.fs.main_compressor.flow_vol.fix()
    get_sCO2_unit_cost(m.fs.main_compressor, 'Barrel type compressor',
                       m.fs.main_compressor.flow_vol,
                       n_equip=5.0)

    # ######################################################
    # Main Compressor Motor
    m.fs.main_compressor_motor = pyo.Block()
    m.fs.main_compressor_motor.work_isentropic = pyo.Var(initialize=159.7e6)
    m.fs.main_compressor_motor.work_isentropic.fix()
    get_sCO2_unit_cost(m.fs.main_compressor_motor, 'Open drip-proof motor',
                       m.fs.main_compressor_motor.work_isentropic*(1e-6),
                       n_equip=5.0)

    # ######################################################
    # Recompressor - 6.89 m^3/s in Baseline620
    m.fs.bypass_compressor = pyo.Block()
    m.fs.bypass_compressor.flow_vol = pyo.Var(initialize=6.89)
    m.fs.bypass_compressor.flow_vol.fix()
    get_sCO2_unit_cost(m.fs.bypass_compressor, 'Barrel type compressor',
                       m.fs.bypass_compressor.flow_vol,
                       n_equip=4.0)

    # ######################################################
    # Recompressor Motor
    m.fs.bypass_compressor_motor = pyo.Block()
    m.fs.bypass_compressor_motor.work_isentropic = pyo.Var(initialize=124.3e6)
    m.fs.bypass_compressor_motor.work_isentropic.fix()

    get_sCO2_unit_cost(m.fs.bypass_compressor_motor, 'Open drip-proof motor',
                       m.fs.bypass_compressor_motor.work_isentropic*(1e-6),
                       n_equip=4.0)

    # add total cost
    build_flowsheet_cost_constraint(m)

    # add initialize
    costing_initialization(m.fs)

    # try solving
    solver = get_solver()
    results = solver.solve(m, tee=True)

    assert results.solver.termination_condition == \
        pyo.TerminationCondition.optimal

    assert pytest.approx(pyo.value(
        m.fs.boiler.costing.equipment_cost), abs=1e-1) == 216300/1e3
    assert pytest.approx(pyo.value(
        m.fs.turbine.costing.equipment_cost), abs=1e-1) == 13160/1e3
    assert pytest.approx(pyo.value(
        m.fs.generator.costing.equipment_cost), abs=1e-1) == 4756/1e3
    assert pytest.approx(pyo.value(
        m.fs.HTR.costing.equipment_cost), abs=1e-1) == 40150/1e3
    assert pytest.approx(pyo.value(
        m.fs.LTR.costing.equipment_cost), abs=1e-1) == 81860/1e3
    assert pytest.approx(pyo.value(
        m.fs.co2_cooler.costing.equipment_cost), abs=1e-1) == 27780/1e3
    assert pytest.approx(pyo.value(
        m.fs.main_compressor.costing.equipment_cost), abs=1e-1) == 31640/1e3
    assert pytest.approx(pyo.value(
        m.fs.bypass_compressor.costing.equipment_cost), abs=1e-1) == 26360/1e3
    assert pytest.approx(
        pyo.value(m.fs.main_compressor_motor.costing.equipment_cost) +
        pyo.value(m.fs.bypass_compressor_motor.costing.equipment_cost),
        abs=1e-1) == 29130/1e3
    assert pytest.approx(pyo.value(
        m.fs.bypass_compressor.costing.equipment_cost), abs=1e-1) == 26360/1e3

    return m


@pytest.mark.component
def test_ASU_costing():
    # Create a Concrete Model as the top level object
    m = pyo.ConcreteModel()

    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.get_costing(year='2017')

    m.fs.ASU = pyo.Block()
    m.fs.ASU.O2_flow = pyo.Var()
    m.fs.ASU.O2_flow.fix(13078)  # TPD

    get_ASU_cost(m.fs.ASU, m.fs.ASU.O2_flow)

    # try solving
    solver = get_solver()
    results = solver.solve(m, tee=True)

    assert results.solver.termination_condition == \
        pyo.TerminationCondition.optimal

    m.fs.ASU.costing.bare_erected_cost.display()

    assert pytest.approx(pyo.value(
        m.fs.ASU.costing.bare_erected_cost), abs=1) == 3.2675e6/1e3

    return m
