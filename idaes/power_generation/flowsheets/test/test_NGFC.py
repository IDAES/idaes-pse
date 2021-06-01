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
Tests to make sure the NGFC example builds and solves correctly.
"""

__author__ = "Alex Noring"

import os

import pytest
import pyomo.environ as pyo
from pyomo.common.fileutils import this_file_dir

from idaes.power_generation.flowsheets.NGFC.NGFC_flowsheet import (
    build_power_island,
    build_reformer,
    set_power_island_inputs,
    set_reformer_inputs,
    scale_flowsheet,
    initialize_power_island,
    initialize_reformer,
    connect_reformer_to_power_island,
    SOFC_ROM_setup,
    add_SOFC_energy_balance,
    add_result_constraints)

from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util import get_solver
from idaes.core.util import model_serializer as ms

solver_available = pyo.SolverFactory('ipopt').available()
solver = get_solver()


@pytest.fixture(scope="module")
def m():
    m = pyo.ConcreteModel(name='NGFC without carbon capture')
    m.fs = FlowsheetBlock(default={"dynamic": False})

    build_power_island(m)
    build_reformer(m)
    set_power_island_inputs(m)
    set_reformer_inputs(m)
    connect_reformer_to_power_island(m)

    return m


@pytest.mark.component
def test_build(m):
    """Build NGFC and check for unit ops"""

    assert degrees_of_freedom(m) == 0

    assert hasattr(m.fs, 'anode')
    assert hasattr(m.fs.anode, 'heat_duty')
    assert hasattr(m.fs, 'cathode_heat')
    assert hasattr(m.fs.cathode_heat, 'heat_duty')
    assert hasattr(m.fs, 'reformer')
    assert hasattr(m.fs.reformer, 'heat_duty')
    assert hasattr(m.fs.reformer, 'deltaP')

    assert not m.fs.anode_mix.feed.flow_mol[0].fixed


@pytest.mark.integration
def test_initialize(m):
    scale_flowsheet(m)
    initialize_power_island(m)
    initialize_reformer(m)

    assert (pyo.value(m.fs.reformer.lagrange_mult[(0, 'H')]) ==
            pytest.approx(81235, 1e-3))
    assert (pyo.value(m.fs.prereformer.lagrange_mult[(0, 'H')]) ==
            pytest.approx(62743, 1e-3))
    assert (pyo.value(m.fs.anode.lagrange_mult[(0, 'H')]) ==
            pytest.approx(78300, 1e-3))

    assert (pyo.value(m.fs.bypass_rejoin.outlet.temperature[0]) ==
            pytest.approx(621.6, 1e-3))
    assert (pyo.value(m.fs.anode_hx.shell_outlet.temperature[0]) ==
            pytest.approx(826.5, 1e-3))
    assert (pyo.value(m.fs.cathode_hx.shell_outlet.temperature[0]) ==
            pytest.approx(475.8, 1e-3))


@pytest.mark.integration
def test_ROM(m):
    SOFC_ROM_setup(m)
    add_SOFC_energy_balance(m)
    add_result_constraints(m)

    assert degrees_of_freedom(m) == 0

    assert hasattr(m.fs, 'ROM_fuel_inlet_temperature')
    assert hasattr(m.fs, 'ROM_internal_reformation')
    assert m.fs.SOFC.current_density.fixed
    assert not m.fs.air_blower.inlet.flow_mol[0].fixed
    assert m.fs.SOFC.deltaT_cell.fixed
    assert hasattr(m.fs, 'stack_power')
    assert hasattr(m.fs, 'SOFC_energy_balance')
    assert not m.fs.anode.outlet.temperature[0].fixed
    assert hasattr(m.fs, 'CO2_emissions')


@pytest.mark.integration
def test_json_load(m):
    fname = os.path.join(os.path.join(os.path.dirname(this_file_dir()),
                                      "NGFC"), "NGFC_flowsheet_init.json.gz")

    ms.from_json(m, fname=fname)

    assert (pyo.value(m.fs.cathode.ion_outlet.flow_mol[0]) ==
            pytest.approx(1670.093, 1e-5))
    assert (pyo.value(m.fs.reformer_recuperator.area) ==
            pytest.approx(4512.56, 1e-3))
    assert (pyo.value(m.fs.anode.heat_duty[0]) ==
            pytest.approx(-672918626, 1e-5))
    assert (pyo.value(m.fs.CO2_emissions) ==
            pytest.approx(291.169, 1e-5))
    assert (pyo.value(m.fs.net_power) ==
            pytest.approx(659.879, 1e-5))


