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
Tests for demonstration flowsheet.

"""

import pytest

from idaes.models.flowsheets.demo_flowsheet import (
    build_flowsheet,
    set_dof,
    initialize_flowsheet,
    solve_flowsheet,
)

from pyomo.environ import value
from pyomo.network import Arc
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.models.properties.activity_coeff_models.BTX_activity_coeff_VLE import (
    BTXParameterBlock,
)
from idaes.models.unit_models import Mixer, Heater, Flash
from idaes.core.util.model_statistics import degrees_of_freedom


@pytest.fixture(scope="module")
def model():
    m = build_flowsheet()

    return m


@pytest.mark.unit
def test_build_flowsheet(model):
    assert isinstance(model.fs, FlowsheetBlock)

    assert isinstance(model.fs.BT_props, BTXParameterBlock)

    assert isinstance(model.fs.M01, Mixer)
    assert isinstance(model.fs.H02, Heater)
    assert isinstance(model.fs.F03, Flash)

    assert isinstance(model.fs.s01, Arc)
    assert isinstance(model.fs.s02, Arc)

    assert degrees_of_freedom(model) == 13


@pytest.mark.unit
def test_set_dof(model):
    set_dof(model)

    assert degrees_of_freedom(model) == 0


@pytest.mark.unit
def test_initialize_flowsheet(model):
    initialize_flowsheet(model)

    assert degrees_of_freedom(model) == 0

    assert model.fs.M01.outlet.flow_mol[0].value == pytest.approx(2.0, 1e-3)
    assert model.fs.H02.outlet.flow_mol[0].value == pytest.approx(2.0, 1e-3)
    assert model.fs.F03.vap_outlet.flow_mol[0].expr.value == pytest.approx(1.367, 1e-3)
    assert model.fs.F03.liq_outlet.flow_mol[0].expr.value == pytest.approx(0.633, 1e-3)


@pytest.mark.integration
def test_unit_consistency(model):
    assert_units_consistent(model)


@pytest.mark.unit
def test_solve_flowsheet(model):
    solve_flowsheet(model)

    assert model.fs.M01.outlet.flow_mol[0].value == pytest.approx(2.0, 1e-4)
    assert model.fs.M01.outlet.mole_frac_comp[0, "benzene"].value == pytest.approx(
        0.5, 1e-4
    )
    assert model.fs.M01.outlet.mole_frac_comp[0, "toluene"].value == pytest.approx(
        0.5, 1e-4
    )
    assert model.fs.M01.outlet.pressure[0].value == pytest.approx(101325, 1e-4)
    assert model.fs.M01.outlet.temperature[0].value == pytest.approx(368.2, 1e-4)

    assert model.fs.H02.outlet.flow_mol[0].value == pytest.approx(2.0, 1e-4)
    assert model.fs.H02.outlet.mole_frac_comp[0, "benzene"].value == pytest.approx(
        0.5, 1e-4
    )
    assert model.fs.H02.outlet.mole_frac_comp[0, "toluene"].value == pytest.approx(
        0.5, 1e-4
    )
    assert model.fs.H02.outlet.pressure[0].value == pytest.approx(101325, 1e-4)
    assert model.fs.H02.outlet.temperature[0].value == pytest.approx(370, 1e-4)

    assert value(model.fs.F03.vap_outlet.flow_mol[0]) == pytest.approx(1.3673, 1e-4)
    assert value(model.fs.F03.vap_outlet.mole_frac_comp[0, "benzene"]) == pytest.approx(
        0.5694, 1e-4
    )
    assert value(model.fs.F03.vap_outlet.mole_frac_comp[0, "toluene"]) == pytest.approx(
        0.4306, 1e-4
    )
    assert value(model.fs.F03.vap_outlet.pressure[0]) == pytest.approx(101325, 1e-4)
    assert value(model.fs.F03.vap_outlet.temperature[0]) == pytest.approx(370, 1e-4)

    assert value(model.fs.F03.liq_outlet.flow_mol[0]) == pytest.approx(0.6327, 1e-4)
    assert value(model.fs.F03.liq_outlet.mole_frac_comp[0, "benzene"]) == pytest.approx(
        0.3501, 1e-4
    )
    assert value(model.fs.F03.liq_outlet.mole_frac_comp[0, "toluene"]) == pytest.approx(
        0.6499, 1e-4
    )
    assert value(model.fs.F03.liq_outlet.pressure[0]) == pytest.approx(101325, 1e-4)
    assert value(model.fs.F03.liq_outlet.temperature[0]) == pytest.approx(370, 1e-4)
