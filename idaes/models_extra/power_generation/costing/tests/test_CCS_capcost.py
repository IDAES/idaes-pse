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
Test code: power plant costing of subaccounts 6.1. 
Subaccounts 6.1-6.8 consist of costing correlations for MEA solvent-based
carbon capture system, Reference: 
Kangkang Li, Ashleigh Cousins, Hai You, Paul Feron,
Weilang Luo, Jian Chen (2016). Systematic study of aqueous monoethanolamine
(MEA)-based CO2 capture process: Techno-economic assessment of the MEA
process and its improvements. Applied Energy, 165, 648-659.
"""

import pytest
from idaes.models_extra.power_generation.costing.power_plant_capcost import (
    QGESSCostingData,
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
    QGESSCostingData.get_PP_costing(
        m.fs.b17,
        absorber_accounts,
        m.fs.b17.absorber_volume,
        "m**3",
        6,
        CE_index_base=567.3,
    )

    # Absorber packing
    absorber_packing_accounts = ["6.2.ccs"]
    m.fs.b18 = pyo.Block()
    # Obtain absorber packing volume in m3
    absorber_packing_volume = 3773  # m3

    m.fs.b18.absorber_packing_volume = pyo.Var(initialize=absorber_packing_volume)
    m.fs.b18.absorber_packing_volume.fix()
    QGESSCostingData.get_PP_costing(
        m.fs.b18,
        absorber_packing_accounts,
        m.fs.b18.absorber_packing_volume,
        "m**3",
        6,
        CE_index_base=567.3,
    )

    # Stripper
    stripper_accounts = ["6.3.ccs"]
    m.fs.b19 = pyo.Block()
    # Obtain stripper volume in m3
    stripper_volume = 386  # m3

    m.fs.b19.stripper_volume = pyo.Var(initialize=stripper_volume)
    m.fs.b19.stripper_volume.fix()
    QGESSCostingData.get_PP_costing(
        m.fs.b19,
        stripper_accounts,
        m.fs.b19.stripper_volume,
        "m**3",
        6,
        CE_index_base=567.3,
    )

    # Stripper packing
    stripper_packing_accounts = ["6.4.ccs"]
    m.fs.b20 = pyo.Block()
    # Obtain stripper packing volume in m3
    stripper_packing_volume = 291  # m3

    m.fs.b20.stripper_packing_volume = pyo.Var(initialize=stripper_packing_volume)
    m.fs.b20.stripper_packing_volume.fix()
    QGESSCostingData.get_PP_costing(
        m.fs.b20,
        stripper_packing_accounts,
        m.fs.b20.stripper_packing_volume,
        "m**3",
        6,
        CE_index_base=567.3,
    )

    # Stripper condenser
    stripper_condenser_accounts = ["6.5.ccs"]
    m.fs.b21 = pyo.Block()
    # Obtain stripper condenser area in m2
    stripper_condenser_area = 1500  # m2

    m.fs.b21.stripper_condenser_area = pyo.Var(initialize=stripper_condenser_area)
    m.fs.b21.stripper_condenser_area.fix()
    QGESSCostingData.get_PP_costing(
        m.fs.b21,
        stripper_condenser_accounts,
        m.fs.b21.stripper_condenser_area,
        "m**2",
        6,
        CE_index_base=567.3,
    )

    # Stripper reboiler
    stripper_reboiler_accounts = ["6.6.ccs"]
    m.fs.b22 = pyo.Block()
    # Obtain stripper reboiler area in m2
    stripper_reboiler_area = 1500  # m2

    m.fs.b22.stripper_reboiler_area = pyo.Var(initialize=stripper_reboiler_area)
    m.fs.b22.stripper_reboiler_area.fix()
    QGESSCostingData.get_PP_costing(
        m.fs.b22,
        stripper_reboiler_accounts,
        m.fs.b22.stripper_reboiler_area,
        "m**2",
        6,
        CE_index_base=567.3,
    )

    # Lean rich heat exchanger
    lean_rich_hex_accounts = ["6.7.ccs"]
    m.fs.b23 = pyo.Block()
    # Obtain lean rich heat exchanger area in m2
    lean_rich_hex_area = 40700  # m2

    m.fs.b23.lean_rich_hex_area = pyo.Var(initialize=lean_rich_hex_area)
    m.fs.b23.lean_rich_hex_area.fix()
    QGESSCostingData.get_PP_costing(
        m.fs.b23,
        lean_rich_hex_accounts,
        m.fs.b23.lean_rich_hex_area,
        "m**2",
        6,
        CE_index_base=567.3,
    )

    # Lean solvent cooler
    lean_solvent_cooler_accounts = ["6.8.ccs"]
    m.fs.b24 = pyo.Block()
    # Obtain lean solvent cooler area in m2
    lean_solvent_cooler_area = 600  # m2

    m.fs.b24.lean_solvent_cooler_area = pyo.Var(initialize=lean_solvent_cooler_area)
    m.fs.b24.lean_solvent_cooler_area.fix()
    QGESSCostingData.get_PP_costing(
        m.fs.b24,
        lean_solvent_cooler_accounts,
        m.fs.b24.lean_solvent_cooler_area,
        "m**2",
        6,
        CE_index_base=567.3,
    )

    # Flue gas blower
    flue_gas_blower_account_1 = ["6.9.1.ccs"]
    flue_gas_blower_account_2 = ["6.9.2.ccs"]

    m.fs.b25 = pyo.Block()
    m.fs.b26 = pyo.Block()
    # Obtain the CO2 product flow in lb/hr
    co2_product_flow = 245000  # lb/hr
    fg_flow = 2e6  # m3/hr

    m.fs.b25.co2_product_flow = pyo.Var(initialize=co2_product_flow)
    m.fs.b25.co2_product_flow.fix()
    m.fs.b26.fg_flow = pyo.Var(initialize=fg_flow)
    m.fs.b26.fg_flow.fix()
    QGESSCostingData.get_PP_costing(
        m.fs.b25,
        flue_gas_blower_account_1,
        m.fs.b25.co2_product_flow,
        "lb/hr",
        6,
        CE_index_base=567.3,
    )
    QGESSCostingData.get_PP_costing(
        m.fs.b26,
        flue_gas_blower_account_2,
        m.fs.b26.fg_flow,
        "m3/hr",
        6,
        CE_index_base=567.3,
    )

    # Flue gas direct contact cooler
    flue_gas_dcc_account_1 = ["6.10.1.ccs"]
    flue_gas_dcc_account_2 = ["6.10.2.ccs"]
    m.fs.b27 = pyo.Block()
    m.fs.b28 = pyo.Block()

    m.fs.b27.co2_product_flow = pyo.Var(initialize=co2_product_flow)
    m.fs.b27.co2_product_flow.fix()
    m.fs.b28.fg_flow = pyo.Var(initialize=fg_flow)
    m.fs.b28.fg_flow.fix()
    QGESSCostingData.get_PP_costing(
        m.fs.b27,
        flue_gas_dcc_account_1,
        m.fs.b27.co2_product_flow,
        "lb/hr",
        6,
        CE_index_base=567.3,
    )
    QGESSCostingData.get_PP_costing(
        m.fs.b28,
        flue_gas_dcc_account_2,
        m.fs.b28.fg_flow,
        "m3/hr",
        6,
        CE_index_base=567.3,
    )

    # Flue gas direct contact cooler packing
    flue_gas_dcc_packing_account_1 = ["6.11.1.ccs"]
    flue_gas_dcc_packing_account_2 = ["6.11.2.ccs"]
    m.fs.b29 = pyo.Block()
    m.fs.b30 = pyo.Block()

    m.fs.b29.co2_product_flow = pyo.Var(initialize=co2_product_flow)
    m.fs.b29.co2_product_flow.fix()
    m.fs.b30.fg_flow = pyo.Var(initialize=fg_flow)
    m.fs.b30.fg_flow.fix()
    QGESSCostingData.get_PP_costing(
        m.fs.b29,
        flue_gas_dcc_packing_account_1,
        m.fs.b29.co2_product_flow,
        "lb/hr",
        6,
        CE_index_base=567.3,
    )
    QGESSCostingData.get_PP_costing(
        m.fs.b30,
        flue_gas_dcc_packing_account_2,
        m.fs.b30.fg_flow,
        "m3/hr",
        6,
        CE_index_base=567.3,
    )

    # Pretreatment pump
    pretreatment_pump_account_1 = ["6.12.1.ccs"]
    pretreatment_pump_account_2 = ["6.12.2.ccs"]
    m.fs.b31 = pyo.Block()
    m.fs.b32 = pyo.Block()

    m.fs.b31.co2_product_flow = pyo.Var(initialize=co2_product_flow)
    m.fs.b31.co2_product_flow.fix()
    m.fs.b32.fg_flow = pyo.Var(initialize=fg_flow)
    m.fs.b32.fg_flow.fix()
    QGESSCostingData.get_PP_costing(
        m.fs.b31,
        pretreatment_pump_account_1,
        m.fs.b31.co2_product_flow,
        "lb/hr",
        6,
        CE_index_base=567.3,
    )
    QGESSCostingData.get_PP_costing(
        m.fs.b32,
        pretreatment_pump_account_2,
        m.fs.b32.fg_flow,
        "m3/hr",
        6,
        CE_index_base=567.3,
    )

    # Pretreatment cooler
    pretreatment_cooler_account_1 = ["6.13.1.ccs"]
    pretreatment_cooler_account_2 = ["6.13.2.ccs"]
    m.fs.b33 = pyo.Block()
    m.fs.b34 = pyo.Block()

    m.fs.b33.co2_product_flow = pyo.Var(initialize=co2_product_flow)
    m.fs.b33.co2_product_flow.fix()
    m.fs.b34.fg_flow = pyo.Var(initialize=fg_flow)
    m.fs.b34.fg_flow.fix()
    QGESSCostingData.get_PP_costing(
        m.fs.b33,
        pretreatment_cooler_account_1,
        m.fs.b33.co2_product_flow,
        "lb/hr",
        6,
        CE_index_base=567.3,
    )
    QGESSCostingData.get_PP_costing(
        m.fs.b34,
        pretreatment_cooler_account_2,
        m.fs.b34.fg_flow,
        "m3/hr",
        6,
        CE_index_base=567.3,
    )

    # Pretreatment tank
    pretreatment_tank_account_1 = ["6.14.1.ccs"]
    pretreatment_tank_account_2 = ["6.14.2.ccs"]
    m.fs.b35 = pyo.Block()
    m.fs.b36 = pyo.Block()

    m.fs.b35.co2_product_flow = pyo.Var(initialize=co2_product_flow)
    m.fs.b35.co2_product_flow.fix()
    m.fs.b36.fg_flow = pyo.Var(initialize=fg_flow)
    m.fs.b36.fg_flow.fix()
    QGESSCostingData.get_PP_costing(
        m.fs.b35,
        pretreatment_tank_account_1,
        m.fs.b35.co2_product_flow,
        "lb/hr",
        6,
        CE_index_base=567.3,
    )
    QGESSCostingData.get_PP_costing(
        m.fs.b36,
        pretreatment_tank_account_2,
        m.fs.b36.fg_flow,
        "m3/hr",
        6,
        CE_index_base=567.3,
    )

    # Washing column
    washing_column_account_1 = ["6.15.1.ccs"]
    washing_column_account_2 = ["6.15.2.ccs"]
    m.fs.b37 = pyo.Block()
    m.fs.b38 = pyo.Block()

    m.fs.b37.co2_product_flow = pyo.Var(initialize=co2_product_flow)
    m.fs.b37.co2_product_flow.fix()
    m.fs.b38.fg_flow = pyo.Var(initialize=fg_flow)
    m.fs.b38.fg_flow.fix()
    QGESSCostingData.get_PP_costing(
        m.fs.b37,
        washing_column_account_1,
        m.fs.b37.co2_product_flow,
        "lb/hr",
        6,
        CE_index_base=567.3,
    )
    QGESSCostingData.get_PP_costing(
        m.fs.b38,
        washing_column_account_2,
        m.fs.b38.fg_flow,
        "m3/hr",
        6,
        CE_index_base=567.3,
    )

    # Washing column packing
    washing_column_packing_account_1 = ["6.16.1.ccs"]
    washing_column_packing_account_2 = ["6.16.2.ccs"]
    m.fs.b39 = pyo.Block()
    m.fs.b40 = pyo.Block()

    m.fs.b39.co2_product_flow = pyo.Var(initialize=co2_product_flow)
    m.fs.b39.co2_product_flow.fix()
    m.fs.b40.fg_flow = pyo.Var(initialize=fg_flow)
    m.fs.b40.fg_flow.fix()
    QGESSCostingData.get_PP_costing(
        m.fs.b39,
        washing_column_packing_account_1,
        m.fs.b39.co2_product_flow,
        "lb/hr",
        6,
        CE_index_base=567.3,
    )
    QGESSCostingData.get_PP_costing(
        m.fs.b40,
        washing_column_packing_account_2,
        m.fs.b40.fg_flow,
        "m3/hr",
        6,
        CE_index_base=567.3,
    )

    # Washing solvent cooler
    washing_solvent_cooler_account_1 = ["6.17.1.ccs"]
    washing_solvent_cooler_account_2 = ["6.17.2.ccs"]
    m.fs.b41 = pyo.Block()
    m.fs.b42 = pyo.Block()

    m.fs.b41.co2_product_flow = pyo.Var(initialize=co2_product_flow)
    m.fs.b41.co2_product_flow.fix()
    m.fs.b42.fg_flow = pyo.Var(initialize=fg_flow)
    m.fs.b42.fg_flow.fix()
    QGESSCostingData.get_PP_costing(
        m.fs.b41,
        washing_solvent_cooler_account_1,
        m.fs.b41.co2_product_flow,
        "lb/hr",
        6,
        CE_index_base=567.3,
    )
    QGESSCostingData.get_PP_costing(
        m.fs.b42,
        washing_solvent_cooler_account_2,
        m.fs.b42.fg_flow,
        "m3/hr",
        6,
        CE_index_base=567.3,
    )

    # Washing solvent pump
    washing_solvent_pump_account_1 = ["6.18.1.ccs"]
    washing_solvent_pump_account_2 = ["6.18.2.ccs"]
    m.fs.b43 = pyo.Block()
    m.fs.b44 = pyo.Block()

    m.fs.b43.co2_product_flow = pyo.Var(initialize=co2_product_flow)
    m.fs.b43.co2_product_flow.fix()
    m.fs.b44.fg_flow = pyo.Var(initialize=fg_flow)
    m.fs.b44.fg_flow.fix()
    QGESSCostingData.get_PP_costing(
        m.fs.b43,
        washing_solvent_pump_account_1,
        m.fs.b43.co2_product_flow,
        "lb/hr",
        6,
        CE_index_base=567.3,
    )
    QGESSCostingData.get_PP_costing(
        m.fs.b44,
        washing_solvent_pump_account_2,
        m.fs.b44.fg_flow,
        "m3/hr",
        6,
        CE_index_base=567.3,
    )

    # Condenser pump
    condenser_pump_account_1 = ["6.19.1.ccs"]
    condenser_pump_account_2 = ["6.19.2.ccs"]
    m.fs.b45 = pyo.Block()
    m.fs.b46 = pyo.Block()

    m.fs.b45.co2_product_flow = pyo.Var(initialize=co2_product_flow)
    m.fs.b45.co2_product_flow.fix()
    m.fs.b46.fg_flow = pyo.Var(initialize=fg_flow)
    m.fs.b46.fg_flow.fix()
    QGESSCostingData.get_PP_costing(
        m.fs.b45,
        condenser_pump_account_1,
        m.fs.b45.co2_product_flow,
        "lb/hr",
        6,
        CE_index_base=567.3,
    )
    QGESSCostingData.get_PP_costing(
        m.fs.b46,
        condenser_pump_account_2,
        m.fs.b46.fg_flow,
        "m3/hr",
        6,
        CE_index_base=567.3,
    )

    # Stripper reflux drum
    stripper_reflux_drum_account_1 = ["6.20.1.ccs"]
    stripper_reflux_drum_account_2 = ["6.20.2.ccs"]
    m.fs.b47 = pyo.Block()
    m.fs.b48 = pyo.Block()

    m.fs.b47.co2_product_flow = pyo.Var(initialize=co2_product_flow)
    m.fs.b47.co2_product_flow.fix()
    m.fs.b48.fg_flow = pyo.Var(initialize=fg_flow)
    m.fs.b48.fg_flow.fix()
    QGESSCostingData.get_PP_costing(
        m.fs.b47,
        stripper_reflux_drum_account_1,
        m.fs.b47.co2_product_flow,
        "lb/hr",
        6,
        CE_index_base=567.3,
    )
    QGESSCostingData.get_PP_costing(
        m.fs.b48,
        stripper_reflux_drum_account_2,
        m.fs.b48.fg_flow,
        "m3/hr",
        6,
        CE_index_base=567.3,
    )

    # Lean solvent pump
    lean_solvent_pump_account_1 = ["6.21.1.ccs"]
    lean_solvent_pump_account_2 = ["6.21.2.ccs"]
    m.fs.b49 = pyo.Block()
    m.fs.b50 = pyo.Block()

    m.fs.b49.co2_product_flow = pyo.Var(initialize=co2_product_flow)
    m.fs.b49.co2_product_flow.fix()
    m.fs.b50.fg_flow = pyo.Var(initialize=fg_flow)
    m.fs.b50.fg_flow.fix()
    QGESSCostingData.get_PP_costing(
        m.fs.b49,
        lean_solvent_pump_account_1,
        m.fs.b49.co2_product_flow,
        "lb/hr",
        6,
        CE_index_base=567.3,
    )
    QGESSCostingData.get_PP_costing(
        m.fs.b50,
        lean_solvent_pump_account_2,
        m.fs.b50.fg_flow,
        "m3/hr",
        6,
        CE_index_base=567.3,
    )

    # Solvent storage tank
    solvent_storage_tank_account_1 = ["6.22.1.ccs"]
    solvent_storage_tank_account_2 = ["6.22.2.ccs"]
    m.fs.b51 = pyo.Block()
    m.fs.b52 = pyo.Block()

    m.fs.b51.co2_product_flow = pyo.Var(initialize=co2_product_flow)
    m.fs.b51.co2_product_flow.fix()
    m.fs.b52.fg_flow = pyo.Var(initialize=fg_flow)
    m.fs.b52.fg_flow.fix()
    QGESSCostingData.get_PP_costing(
        m.fs.b51,
        solvent_storage_tank_account_1,
        m.fs.b51.co2_product_flow,
        "lb/hr",
        6,
        CE_index_base=567.3,
    )
    QGESSCostingData.get_PP_costing(
        m.fs.b52,
        solvent_storage_tank_account_2,
        m.fs.b52.fg_flow,
        "m3/hr",
        6,
        CE_index_base=567.3,
    )

    # Washing solvent tank
    washing_solvent_tank_account_1 = ["6.23.1.ccs"]
    washing_solvent_tank_account_2 = ["6.23.2.ccs"]
    m.fs.b53 = pyo.Block()
    m.fs.b54 = pyo.Block()

    m.fs.b53.co2_product_flow = pyo.Var(initialize=co2_product_flow)
    m.fs.b53.co2_product_flow.fix()
    m.fs.b54.fg_flow = pyo.Var(initialize=fg_flow)
    m.fs.b54.fg_flow.fix()
    QGESSCostingData.get_PP_costing(
        m.fs.b53,
        washing_solvent_tank_account_1,
        m.fs.b53.co2_product_flow,
        "lb/hr",
        6,
        CE_index_base=567.3,
    )
    QGESSCostingData.get_PP_costing(
        m.fs.b54,
        washing_solvent_tank_account_2,
        m.fs.b54.fg_flow,
        "m3/hr",
        6,
        CE_index_base=567.3,
    )

    # Solvent stripper reclaimer
    solvent_stripper_reclaimer_account_1 = ["6.24.1.ccs"]
    solvent_stripper_reclaimer_account_2 = ["6.24.2.ccs"]
    m.fs.b55 = pyo.Block()
    m.fs.b56 = pyo.Block()

    m.fs.b55.co2_product_flow = pyo.Var(initialize=co2_product_flow)
    m.fs.b55.co2_product_flow.fix()
    m.fs.b56.fg_flow = pyo.Var(initialize=fg_flow)
    m.fs.b56.fg_flow.fix()
    QGESSCostingData.get_PP_costing(
        m.fs.b55,
        solvent_stripper_reclaimer_account_1,
        m.fs.b55.co2_product_flow,
        "lb/hr",
        6,
        CE_index_base=567.3,
    )
    QGESSCostingData.get_PP_costing(
        m.fs.b56,
        solvent_stripper_reclaimer_account_2,
        m.fs.b56.fg_flow,
        "m3/hr",
        6,
        CE_index_base=567.3,
    )

    # Solvent reclaimer cooler
    solvent_reclaimer_cooler_account_1 = ["6.25.1.ccs"]
    solvent_reclaimer_cooler_account_2 = ["6.25.2.ccs"]
    m.fs.b57 = pyo.Block()
    m.fs.b58 = pyo.Block()

    m.fs.b57.co2_product_flow = pyo.Var(initialize=co2_product_flow)
    m.fs.b57.co2_product_flow.fix()
    m.fs.b58.fg_flow = pyo.Var(initialize=fg_flow)
    m.fs.b58.fg_flow.fix()
    QGESSCostingData.get_PP_costing(
        m.fs.b57,
        solvent_reclaimer_cooler_account_1,
        m.fs.b57.co2_product_flow,
        "lb/hr",
        6,
        CE_index_base=567.3,
    )
    QGESSCostingData.get_PP_costing(
        m.fs.b58,
        solvent_reclaimer_cooler_account_2,
        m.fs.b58.fg_flow,
        "m3/hr",
        6,
        CE_index_base=567.3,
    )

    # Solvent filtration
    solvent_filtration_account_1 = ["6.26.1.ccs"]
    solvent_filtration_account_2 = ["6.26.2.ccs"]
    m.fs.b59 = pyo.Block()
    m.fs.b60 = pyo.Block()

    m.fs.b59.co2_product_flow = pyo.Var(initialize=co2_product_flow)
    m.fs.b59.co2_product_flow.fix()
    m.fs.b60.fg_flow = pyo.Var(initialize=fg_flow)
    m.fs.b60.fg_flow.fix()
    QGESSCostingData.get_PP_costing(
        m.fs.b59,
        solvent_filtration_account_1,
        m.fs.b59.co2_product_flow,
        "lb/hr",
        6,
        CE_index_base=567.3,
    )
    QGESSCostingData.get_PP_costing(
        m.fs.b60,
        solvent_filtration_account_2,
        m.fs.b60.fg_flow,
        "m3/hr",
        6,
        CE_index_base=567.3,
    )

    QGESSCostingData.get_total_TPC(m.fs)

    # Initialize costing
    QGESSCostingData.costing_initialization(m.fs)
    assert degrees_of_freedom(m) == 0

    # Solve the model
    results = solver.solve(m, tee=True)
    assert results.solver.termination_condition == pyo.TerminationCondition.optimal
    assert results.solver.status == pyo.SolverStatus.ok

    # Absorber TPC
    assert pytest.approx(24.417, rel=1e-3) == sum(
        pyo.value(m.fs.b17.total_plant_cost[ac]) for ac in absorber_accounts
    )

    # Absorber packing TPC
    assert pytest.approx(24.4160, rel=1e-3) == sum(
        pyo.value(m.fs.b18.total_plant_cost[ac]) for ac in absorber_packing_accounts
    )

    # Stripper TPC
    assert pytest.approx(2.992, rel=1e-3) == sum(
        pyo.value(m.fs.b19.total_plant_cost[ac]) for ac in stripper_accounts
    )

    # Stripper packing TPC
    assert pytest.approx(2.196, rel=1e-3) == sum(
        pyo.value(m.fs.b20.total_plant_cost[ac]) for ac in stripper_packing_accounts
    )

    # Stripper condenser TPC
    assert pytest.approx(0.601, rel=1e-3) == sum(
        pyo.value(m.fs.b21.total_plant_cost[ac]) for ac in stripper_condenser_accounts
    )

    # Stripper reboiler TPC
    assert pytest.approx(2.624, rel=1e-3) == sum(
        pyo.value(m.fs.b22.total_plant_cost[ac]) for ac in stripper_reboiler_accounts
    )

    # Lean rich heat exchanger TPC
    assert pytest.approx(4.479, rel=1e-3) == sum(
        pyo.value(m.fs.b23.total_plant_cost[ac]) for ac in lean_rich_hex_accounts
    )

    # Lean solvent cooler TPC
    assert pytest.approx(0.486, rel=1e-3) == sum(
        pyo.value(m.fs.b24.total_plant_cost[ac]) for ac in lean_solvent_cooler_accounts
    )

    # Flue gas blower TPC
    assert pytest.approx(0.6668, rel=1e-3) == sum(
        pyo.value(m.fs.b25.total_plant_cost[ac]) for ac in flue_gas_blower_account_1
    )
    assert pytest.approx(0.8907, rel=1e-3) == sum(
        pyo.value(m.fs.b26.total_plant_cost[ac]) for ac in flue_gas_blower_account_2
    )

    # Flue gas direct contact cooler TPC
    assert pytest.approx(1.8991, rel=1e-3) == sum(
        pyo.value(m.fs.b27.total_plant_cost[ac]) for ac in flue_gas_dcc_account_1
    )
    assert pytest.approx(2.5367, rel=1e-3) == sum(
        pyo.value(m.fs.b28.total_plant_cost[ac]) for ac in flue_gas_dcc_account_2
    )

    # Flue gas direct contact cooler packing TPC
    assert pytest.approx(1.6921, rel=1e-3) == sum(
        pyo.value(m.fs.b29.total_plant_cost[ac]) for ac in flue_gas_dcc_packing_account_1
    )
    assert pytest.approx(2.2603, rel=1e-3) == sum(
        pyo.value(m.fs.b30.total_plant_cost[ac]) for ac in flue_gas_dcc_packing_account_2
    )

    # Pretreatment pump TPC
    assert pytest.approx(0.08118, rel=1e-3) == sum(
        pyo.value(m.fs.b31.total_plant_cost[ac]) for ac in pretreatment_pump_account_1
    )
    assert pytest.approx(0.1084, rel=1e-3) == sum(
        pyo.value(m.fs.b32.total_plant_cost[ac]) for ac in pretreatment_pump_account_2
    )

    # Pretreatment cooler TPC
    assert pytest.approx(0.1505, rel=1e-3) == sum(
        pyo.value(m.fs.b33.total_plant_cost[ac]) for ac in pretreatment_cooler_account_1
    )
    assert pytest.approx(0.2010, rel=1e-3) == sum(
        pyo.value(m.fs.b34.total_plant_cost[ac]) for ac in pretreatment_cooler_account_2
    )

    # Pretreatment tank TPC
    assert pytest.approx(0.0675, rel=1e-3) == sum(
        pyo.value(m.fs.b35.total_plant_cost[ac]) for ac in pretreatment_tank_account_1
    )
    assert pytest.approx(0.0902, rel=1e-3) == sum(
        pyo.value(m.fs.b36.total_plant_cost[ac]) for ac in pretreatment_tank_account_2
    )

    # Washing column TPC
    assert pytest.approx(1.8171, rel=1e-3) == sum(
        pyo.value(m.fs.b37.total_plant_cost[ac]) for ac in washing_column_account_1
    )
    assert pytest.approx(2.4272, rel=1e-3) == sum(
        pyo.value(m.fs.b38.total_plant_cost[ac]) for ac in washing_column_account_2
    )

    # Washing column packing TPC
    assert pytest.approx(1.8572, rel=1e-3) == sum(
        pyo.value(m.fs.b39.total_plant_cost[ac])
        for ac in washing_column_packing_account_1
    )
    assert pytest.approx(2.4808, rel=1e-3) == sum(
        pyo.value(m.fs.b40.total_plant_cost[ac])
        for ac in washing_column_packing_account_2
    )

    # Washing solvent cooler TPC
    assert pytest.approx(0.04196, rel=1e-3) == sum(
        pyo.value(m.fs.b41.total_plant_cost[ac])
        for ac in washing_solvent_cooler_account_1
    )
    assert pytest.approx(0.05605, rel=1e-3) == sum(
        pyo.value(m.fs.b42.total_plant_cost[ac])
        for ac in washing_solvent_cooler_account_2
    )

    # Washing solvent pump TPC
    assert pytest.approx(0.008210, rel=1e-3) == sum(
        pyo.value(m.fs.b43.total_plant_cost[ac]) for ac in washing_solvent_pump_account_1
    )
    assert pytest.approx(0.01097, rel=1e-3) == sum(
        pyo.value(m.fs.b44.total_plant_cost[ac]) for ac in washing_solvent_pump_account_2
    )

    # Condenser pump TPC
    assert pytest.approx(0.02280, rel=1e-3) == sum(
        pyo.value(m.fs.b45.total_plant_cost[ac]) for ac in condenser_pump_account_1
    )
    assert pytest.approx(0.03046, rel=1e-3) == sum(
        pyo.value(m.fs.b46.total_plant_cost[ac]) for ac in condenser_pump_account_2
    )

    # Stripper reflux drum TPC
    assert pytest.approx(0.0310, rel=1e-3) == sum(
        pyo.value(m.fs.b47.total_plant_cost[ac]) for ac in stripper_reflux_drum_account_1
    )
    assert pytest.approx(0.0414, rel=1e-3) == sum(
        pyo.value(m.fs.b48.total_plant_cost[ac]) for ac in stripper_reflux_drum_account_2
    )

    # Lean solvent pump TPC
    assert pytest.approx(0.2372, rel=1e-3) == sum(
        pyo.value(m.fs.b49.total_plant_cost[ac]) for ac in lean_solvent_pump_account_1
    )
    assert pytest.approx(0.3168, rel=1e-3) == sum(
        pyo.value(m.fs.b50.total_plant_cost[ac]) for ac in lean_solvent_pump_account_2
    )

    # Solvent storage tank TPC
    assert pytest.approx(0.2700, rel=1e-3) == sum(
        pyo.value(m.fs.b51.total_plant_cost[ac]) for ac in solvent_storage_tank_account_1
    )
    assert pytest.approx(0.3607, rel=1e-3) == sum(
        pyo.value(m.fs.b52.total_plant_cost[ac]) for ac in solvent_storage_tank_account_2
    )

    # Washing solvent tank TPC
    assert pytest.approx(0.03101, rel=1e-3) == sum(
        pyo.value(m.fs.b53.total_plant_cost[ac]) for ac in washing_solvent_tank_account_1
    )
    assert pytest.approx(0.04143, rel=1e-3) == sum(
        pyo.value(m.fs.b54.total_plant_cost[ac]) for ac in washing_solvent_tank_account_2
    )

    # Solvent stripper reclaimer TPC
    assert pytest.approx(0.1314, rel=1e-3) == sum(
        pyo.value(m.fs.b55.total_plant_cost[ac])
        for ac in solvent_stripper_reclaimer_account_1
    )
    assert pytest.approx(0.1755, rel=1e-3) == sum(
        pyo.value(m.fs.b56.total_plant_cost[ac])
        for ac in solvent_stripper_reclaimer_account_2
    )

    # Solvent reclaimer cooler TPC
    assert pytest.approx(0.1231, rel=1e-3) == sum(
        pyo.value(m.fs.b57.total_plant_cost[ac])
        for ac in solvent_reclaimer_cooler_account_1
    )
    assert pytest.approx(0.1645, rel=1e-3) == sum(
        pyo.value(m.fs.b58.total_plant_cost[ac])
        for ac in solvent_reclaimer_cooler_account_2
    )

    # Solvent filtration TPC
    assert pytest.approx(0.7215, rel=1e-3) == sum(
        pyo.value(m.fs.b59.total_plant_cost[ac]) for ac in solvent_filtration_account_1
    )
    assert pytest.approx(0.9638, rel=1e-3) == sum(
        pyo.value(m.fs.b60.total_plant_cost[ac]) for ac in solvent_filtration_account_2
    )

    # Total TPC
    TPC = 0
    for idx in range(17, 61):
        block = getattr(m.fs, 'b'+str(idx))
        for ac in block.total_plant_cost.keys():
            TPC += pyo.value(block.total_plant_cost[ac])


    assert pytest.approx(85.218, rel=1e-3) == TPC

    return m
