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
'''
Test code: power plant costing of subaccounts 6.1. 
Subaccounts 6.1-6.8 consist of costing correlations for MEA solvent-based
carbon capture system, Reference: 
Kangkang Li, Ashleigh Cousins, Hai You, Paul Feron,
Weilang Luo, Jian Chen (2016). Systematic study of aqueous monoethanolamine
(MEA)-based CO2 capture process: Techno-economic assessment of the MEA
process and its improvements. Applied Energy, 165, 648-659.
'''

import pytest
from idaes.models_extra.power_generation.costing.power_plant_costing import (
    get_PP_costing,
    get_total_TPC,
    costing_initialization,
)
from idaes.core.util.model_statistics import degrees_of_freedom
import pyomo.environ as pyo
from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver

# Get default solver for testing
solver = get_solver()


@pytest.fixture(scope="module")
def build_costing():
    # Create a Concrete Model as the top level object
    m = pyo.ConcreteModel()

    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.get_costing(year="2018")

    return m


@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_ccs_units_costing():
    m = pyo.ConcreteModel()

    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.get_costing(year="2018")

    # Accounts with carbon capture system units
    # Absorber
    absorber_accounts = ["6.1.ccs"]
    m.fs.b17 = pyo.Block()
    # Obtain absorber volume in m3
    absorber_volume = 4999  # m3

    m.fs.b17.absorber_volume = pyo.Var(initialize=absorber_volume)
    m.fs.b17.absorber_volume.fix()
    get_PP_costing(m.fs.b17, 
        absorber_accounts, 
        m.fs.b17.absorber_volume, 
        "m**3", 
        6, 
        CE_index_base=567.3
    )

    # Absorber packing
    absorber_packing_accounts = ["6.2.ccs"]
    m.fs.b18 = pyo.Block()
    # Obtain absorber packing volume in m3
    absorber_packing_volume = 3773  # m3

    m.fs.b18.absorber_packing_volume = pyo.Var(initialize=absorber_packing_volume)
    m.fs.b18.absorber_packing_volume.fix()
    get_PP_costing(
        m.fs.b18, 
        absorber_packing_accounts, 
        m.fs.b18.absorber_packing_volume, 
        "m**3", 
        6, 
        CE_index_base=567.3
    )

    # Stripper
    stripper_accounts = ["6.3.ccs"]
    m.fs.b19 = pyo.Block()
    # Obtain stripper volume in m3
    stripper_volume = 386  # m3

    m.fs.b19.stripper_volume = pyo.Var(initialize=stripper_volume)
    m.fs.b19.stripper_volume.fix()
    get_PP_costing(m.fs.b19, 
        stripper_accounts, 
        m.fs.b19.stripper_volume, 
        "m**3", 
        6, 
        CE_index_base=567.3
    )

    # Stripper packing
    stripper_packing_accounts = ["6.4.ccs"]
    m.fs.b20 = pyo.Block()
    # Obtain stripper packing volume in m3
    stripper_packing_volume = 291  # m3

    m.fs.b20.stripper_packing_volume = pyo.Var(initialize=stripper_packing_volume)
    m.fs.b20.stripper_packing_volume.fix()
    get_PP_costing(
        m.fs.b20, 
        stripper_packing_accounts, 
        m.fs.b20.stripper_packing_volume, 
        "m**3", 
        6, 
        CE_index_base=567.3
    )

    # Stripper condenser
    stripper_condenser_accounts = ["6.5.ccs"]
    m.fs.b21 = pyo.Block()
    # Obtain stripper condenser area in m2
    stripper_condenser_area = 1500  # m2

    m.fs.b21.stripper_condenser_area = pyo.Var(initialize=stripper_condenser_area)
    m.fs.b21.stripper_condenser_area.fix()
    get_PP_costing(
        m.fs.b21,
        stripper_condenser_accounts,
        m.fs.b21.stripper_condenser_area,
        "m**2",
        6,
        CE_index_base=567.3
    )

    # Stripper reboiler
    stripper_reboiler_accounts = ["6.6.ccs"]
    m.fs.b22 = pyo.Block()
    # Obtain stripper reboiler area in m2
    stripper_reboiler_area = 1500  # m2

    m.fs.b22.stripper_reboiler_area = pyo.Var(initialize=stripper_reboiler_area)
    m.fs.b22.stripper_reboiler_area.fix()
    get_PP_costing(
        m.fs.b22, 
        stripper_reboiler_accounts, 
        m.fs.b22.stripper_reboiler_area, 
        "m**2", 
        6, 
        CE_index_base=567.3
    )

    # Lean rich heat exchanger
    lean_rich_hex_accounts = ["6.7.ccs"]
    m.fs.b23 = pyo.Block()
    # Obtain lean rich heat exchanger area in m2
    lean_rich_hex_area = 40700  # m2

    m.fs.b23.lean_rich_hex_area = pyo.Var(initialize=lean_rich_hex_area)
    m.fs.b23.lean_rich_hex_area.fix()
    get_PP_costing(
        m.fs.b23, 
        lean_rich_hex_accounts, 
        m.fs.b23.lean_rich_hex_area, 
        "m**2", 
        6, 
        CE_index_base=567.3
    )

    # Lean solvent cooler
    lean_solvent_cooler_accounts = ["6.8.ccs"]
    m.fs.b24 = pyo.Block()
    # Obtain lean solvent cooler area in m2
    lean_solvent_cooler_area = 600  # m2

    m.fs.b24.lean_solvent_cooler_area = pyo.Var(initialize=lean_solvent_cooler_area)
    m.fs.b24.lean_solvent_cooler_area.fix()
    get_PP_costing(
        m.fs.b24,
        lean_solvent_cooler_accounts,
        m.fs.b24.lean_solvent_cooler_area,
        "m**2",
        6,
        CE_index_base=567.3
    )
    
    # Flue gas blower
    flue_gas_blower_accounts = ["6.9.ccs"]
    m.fs.b25 = pyo.Block()
    # Obtain the amount of CO2 captured in lb/hr
    amt_co2_capture = 246552  # lb/hr

    m.fs.b25.amt_co2_capture = pyo.Var(initialize=amt_co2_capture)
    m.fs.b25.amt_co2_capture.fix()
    get_PP_costing(
        m.fs.b25,
        flue_gas_blower_accounts,
        m.fs.b25.amt_co2_capture,
        "lb/hr",
        6,
        CE_index_base=567.3
    )
    
    # Flue gas direct contact cooler
    flue_gas_dcc_accounts = ["6.10.ccs"]
    m.fs.b26 = pyo.Block()

    m.fs.b26.amt_co2_capture = pyo.Var(initialize=amt_co2_capture)
    m.fs.b26.amt_co2_capture.fix()
    get_PP_costing(
        m.fs.b26,
        flue_gas_dcc_accounts,
        m.fs.b26.amt_co2_capture,
        "lb/hr",
        6,
        CE_index_base=567.3
    )
    
    # Flue gas direct contact cooler packing
    flue_gas_dcc_packing_accounts = ["6.11.ccs"]
    m.fs.b27 = pyo.Block()

    m.fs.b27.amt_co2_capture = pyo.Var(initialize=amt_co2_capture)
    m.fs.b27.amt_co2_capture.fix()
    get_PP_costing(
        m.fs.b27,
        flue_gas_dcc_packing_accounts,
        m.fs.b27.amt_co2_capture,
        "lb/hr",
        6,
        CE_index_base=567.3
    )
    
    # Pretreatment pump
    pretreatment_pump_accounts = ["6.12.ccs"]
    m.fs.b28 = pyo.Block()

    m.fs.b28.amt_co2_capture = pyo.Var(initialize=amt_co2_capture)
    m.fs.b28.amt_co2_capture.fix()
    get_PP_costing(
        m.fs.b28,
        pretreatment_pump_accounts,
        m.fs.b28.amt_co2_capture,
        "lb/hr",
        6,
        CE_index_base=567.3
    )
    
    # Pretreatment cooler
    pretreatment_cooler_accounts = ["6.13.ccs"]
    m.fs.b29 = pyo.Block()

    m.fs.b29.amt_co2_capture = pyo.Var(initialize=amt_co2_capture)
    m.fs.b29.amt_co2_capture.fix()
    get_PP_costing(
        m.fs.b29,
        pretreatment_cooler_accounts,
        m.fs.b29.amt_co2_capture,
        "lb/hr",
        6,
        CE_index_base=567.3
    )
    
    # Pretreatment tank
    pretreatment_tank_accounts = ["6.14.ccs"]
    m.fs.b30 = pyo.Block()

    m.fs.b30.amt_co2_capture = pyo.Var(initialize=amt_co2_capture)
    m.fs.b30.amt_co2_capture.fix()
    get_PP_costing(
        m.fs.b30,
        pretreatment_tank_accounts,
        m.fs.b30.amt_co2_capture,
        "lb/hr",
        6,
        CE_index_base=567.3
    )
    
    # Washing column
    washing_column_accounts = ["6.15.ccs"]
    m.fs.b31 = pyo.Block()

    m.fs.b31.amt_co2_capture = pyo.Var(initialize=amt_co2_capture)
    m.fs.b31.amt_co2_capture.fix()
    get_PP_costing(
        m.fs.b31,
        washing_column_accounts,
        m.fs.b31.amt_co2_capture,
        "lb/hr",
        6,
        CE_index_base=567.3
    )
    
    # Washing column packing
    washing_column_packing_accounts = ["6.16.ccs"]
    m.fs.b32 = pyo.Block()

    m.fs.b32.amt_co2_capture = pyo.Var(initialize=amt_co2_capture)
    m.fs.b32.amt_co2_capture.fix()
    get_PP_costing(
        m.fs.b32,
        washing_column_packing_accounts,
        m.fs.b32.amt_co2_capture,
        "lb/hr",
        6,
        CE_index_base=567.3
    )
    
    # Washing solvent cooler
    washing_solvent_cooler_accounts = ["6.17.ccs"]
    m.fs.b33 = pyo.Block()

    m.fs.b33.amt_co2_capture = pyo.Var(initialize=amt_co2_capture)
    m.fs.b33.amt_co2_capture.fix()
    get_PP_costing(
        m.fs.b33,
        washing_solvent_cooler_accounts,
        m.fs.b33.amt_co2_capture,
        "lb/hr",
        6,
        CE_index_base=567.3
    )
    
    # Washing solvent pump
    washing_solvent_pump_accounts = ["6.18.ccs"]
    m.fs.b34 = pyo.Block()

    m.fs.b34.amt_co2_capture = pyo.Var(initialize=amt_co2_capture)
    m.fs.b34.amt_co2_capture.fix()
    get_PP_costing(
        m.fs.b34,
        washing_solvent_pump_accounts,
        m.fs.b34.amt_co2_capture,
        "lb/hr",
        6,
        CE_index_base=567.3
    )
    
    # Condenser pump
    condenser_pump_accounts = ["6.19.ccs"]
    m.fs.b35 = pyo.Block()

    m.fs.b35.amt_co2_capture = pyo.Var(initialize=amt_co2_capture)
    m.fs.b35.amt_co2_capture.fix()
    get_PP_costing(
        m.fs.b35,
        condenser_pump_accounts,
        m.fs.b35.amt_co2_capture,
        "lb/hr",
        6,
        CE_index_base=567.3
    )
    
    # Stripper reflux drum
    stripper_reflux_drum_accounts = ["6.20.ccs"]
    m.fs.b36 = pyo.Block()

    m.fs.b36.amt_co2_capture = pyo.Var(initialize=amt_co2_capture)
    m.fs.b36.amt_co2_capture.fix()
    get_PP_costing(
        m.fs.b36,
        stripper_reflux_drum_accounts,
        m.fs.b36.amt_co2_capture,
        "lb/hr",
        6,
        CE_index_base=567.3
    )
    
    # Lean solvent pump
    lean_solvent_pump_accounts = ["6.21.ccs"]
    m.fs.b37 = pyo.Block()

    m.fs.b37.amt_co2_capture = pyo.Var(initialize=amt_co2_capture)
    m.fs.b37.amt_co2_capture.fix()
    get_PP_costing(
        m.fs.b37,
        lean_solvent_pump_accounts,
        m.fs.b37.amt_co2_capture,
        "lb/hr",
        6,
        CE_index_base=567.3
    )
    
    # Solvent storage tank
    solvent_storage_tank_accounts = ["6.22.ccs"]
    m.fs.b38 = pyo.Block()

    m.fs.b38.amt_co2_capture = pyo.Var(initialize=amt_co2_capture)
    m.fs.b38.amt_co2_capture.fix()
    get_PP_costing(
        m.fs.b38,
        solvent_storage_tank_accounts,
        m.fs.b38.amt_co2_capture,
        "lb/hr",
        6,
        CE_index_base=567.3
    )
    
    # Washing solvent tank
    washing_solvent_tank_accounts = ["6.23.ccs"]
    m.fs.b39 = pyo.Block()

    m.fs.b39.amt_co2_capture = pyo.Var(initialize=amt_co2_capture)
    m.fs.b39.amt_co2_capture.fix()
    get_PP_costing(
        m.fs.b39,
        washing_solvent_tank_accounts,
        m.fs.b39.amt_co2_capture,
        "lb/hr",
        6,
        CE_index_base=567.3
    )
    
    # Solvent stripper reclaimer
    solvent_stripper_reclaimer_accounts = ["6.24.ccs"]
    m.fs.b40 = pyo.Block()

    m.fs.b40.amt_co2_capture = pyo.Var(initialize=amt_co2_capture)
    m.fs.b40.amt_co2_capture.fix()
    get_PP_costing(
        m.fs.b40,
        solvent_stripper_reclaimer_accounts,
        m.fs.b40.amt_co2_capture,
        "lb/hr",
        6,
        CE_index_base=567.3
    )
    
    # Solvent reclaimer cooler
    solvent_reclaimer_cooler_accounts = ["6.25.ccs"]
    m.fs.b41 = pyo.Block()

    m.fs.b41.amt_co2_capture = pyo.Var(initialize=amt_co2_capture)
    m.fs.b41.amt_co2_capture.fix()
    get_PP_costing(
        m.fs.b41,
        solvent_reclaimer_cooler_accounts,
        m.fs.b41.amt_co2_capture,
        "lb/hr",
        6,
        CE_index_base=567.3
    )
    
    # Solvent filtration
    solvent_filtration_accounts = ["6.26.ccs"]
    m.fs.b42 = pyo.Block()

    m.fs.b42.amt_co2_capture = pyo.Var(initialize=amt_co2_capture)
    m.fs.b42.amt_co2_capture.fix()
    get_PP_costing(
        m.fs.b42,
        solvent_filtration_accounts,
        m.fs.b42.amt_co2_capture,
        "lb/hr",
        6,
        CE_index_base=567.3
    )

    # Initialize costing
    costing_initialization(m.fs)
    assert degrees_of_freedom(m) == 0

    # Solve the model
    results = solver.solve(m, tee=True)
    assert results.solver.termination_condition == pyo.TerminationCondition.optimal
    assert results.solver.status == pyo.SolverStatus.ok

    # Absorber TPC
    assert pytest.approx(28.885, abs=0.5) == sum(
        pyo.value(m.fs.b17.costing.total_plant_cost[ac]) for ac in absorber_accounts
    )

    # Absorber packing TPC
    assert pytest.approx(28.8834, abs=0.5) == sum(
        pyo.value(m.fs.b18.costing.total_plant_cost[ac])
        for ac in absorber_packing_accounts
    )

    # Stripper TPC
    assert pytest.approx(3.539, abs=0.5) == sum(
        pyo.value(m.fs.b19.costing.total_plant_cost[ac]) for ac in stripper_accounts
    )

    # Stripper packing TPC
    assert pytest.approx(2.598, abs=0.5) == sum(
        pyo.value(m.fs.b20.costing.total_plant_cost[ac])
        for ac in stripper_packing_accounts
    )

    # Stripper condenser TPC
    assert pytest.approx(0.7104, abs=0.5) == sum(
        pyo.value(m.fs.b21.costing.total_plant_cost[ac])
        for ac in stripper_condenser_accounts
    )

    # Stripper reboiler TPC
    assert pytest.approx(3.10464, abs=0.5) == sum(
        pyo.value(m.fs.b22.costing.total_plant_cost[ac])
        for ac in stripper_reboiler_accounts
    )

    # Lean rich heat exchanger TPC
    assert pytest.approx(5.2983, abs=0.5) == sum(
        pyo.value(m.fs.b23.costing.total_plant_cost[ac])
        for ac in lean_rich_hex_accounts
    )

    # Lean solvent cooler TPC
    assert pytest.approx(0.5748, abs=0.5) == sum(
        pyo.value(m.fs.b24.costing.total_plant_cost[ac])
        for ac in lean_solvent_cooler_accounts
    )
    return m
