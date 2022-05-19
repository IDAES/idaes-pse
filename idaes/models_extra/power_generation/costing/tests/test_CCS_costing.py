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
