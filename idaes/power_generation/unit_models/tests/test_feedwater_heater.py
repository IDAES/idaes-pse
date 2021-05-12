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

import pytest
import pyomo.environ as pyo
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.generic_models.properties import iapws95
from idaes.power_generation.unit_models import FWH0D
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util import get_solver

# Set up solver
solver = get_solver()


@pytest.mark.unit
def test_fwh_model():
    model = pyo.ConcreteModel()
    model.fs = FlowsheetBlock(default={
        "dynamic": False,
        "default_property_package": iapws95.Iapws95ParameterBlock()})
    model.fs.properties = model.fs.config.default_property_package
    model.fs.fwh = FWH0D(default={
        "has_desuperheat": True,
        "has_drain_cooling": True,
        "has_drain_mixer": True,
        "property_package": model.fs.properties})

    model.fs.fwh.desuperheat.inlet_1.flow_mol.fix(100)
    model.fs.fwh.desuperheat.inlet_1.flow_mol.unfix()
    model.fs.fwh.desuperheat.inlet_1.pressure.fix(201325)
    model.fs.fwh.desuperheat.inlet_1.enth_mol.fix(60000)
    model.fs.fwh.drain_mix.drain.flow_mol.fix(1)
    model.fs.fwh.drain_mix.drain.pressure.fix(201325)
    model.fs.fwh.drain_mix.drain.enth_mol.fix(20000)
    model.fs.fwh.cooling.inlet_2.flow_mol.fix(400)
    model.fs.fwh.cooling.inlet_2.pressure.fix(101325)
    model.fs.fwh.cooling.inlet_2.enth_mol.fix(3000)
    model.fs.fwh.condense.area.fix(1000)
    model.fs.fwh.condense.overall_heat_transfer_coefficient.fix(100)
    model.fs.fwh.desuperheat.area.fix(1000)
    model.fs.fwh.desuperheat.overall_heat_transfer_coefficient.fix(10)
    model.fs.fwh.cooling.area.fix(1000)
    model.fs.fwh.cooling.overall_heat_transfer_coefficient.fix(10)
    model.fs.fwh.initialize(optarg={"max_iter": 50})

    assert(degrees_of_freedom(model) == 0)
    assert(abs(pyo.value(model.fs.fwh.desuperheat.inlet_1.flow_mol[0]) -
               98.335) < 0.01)


@pytest.mark.integration
def test_fwh_units():
    model = pyo.ConcreteModel()
    model.fs = FlowsheetBlock(default={
        "dynamic": False,
        "default_property_package": iapws95.Iapws95ParameterBlock()})
    model.fs.properties = model.fs.config.default_property_package
    model.fs.fwh = FWH0D(default={
        "has_desuperheat": True,
        "has_drain_cooling": True,
        "has_drain_mixer": True,
        "property_package": model.fs.properties})

    assert_units_consistent(model)
