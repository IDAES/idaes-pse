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

from idaes.core import FlowsheetBlock
from idaes.generic_models.properties import iapws95
from idaes.power_generation.unit_models import FWH0DDynamic
from idaes.core.util.model_statistics import degrees_of_freedom

prop_available = iapws95.iapws95_available()

# See if ipopt is available and set up solver
if pyo.SolverFactory('ipopt').available():
    solver = pyo.SolverFactory('ipopt')
else:
    solver = None


@pytest.mark.integration
@pytest.mark.skipif(not prop_available, reason="IAPWS not available")
@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_fwh_model():
    model = pyo.ConcreteModel()
    model.fs = FlowsheetBlock(default={
        "dynamic": True,
        "time_set": [0, 60],
        "time_units": pyo.units.s,
        "default_property_package": iapws95.Iapws95ParameterBlock()})
    model.fs.properties = model.fs.config.default_property_package
    model.fs.fwh = FWH0DDynamic(default={
        "has_desuperheat": True,
        "has_drain_cooling": True,
        "has_drain_mixer": True,
        "condense": {"tube": {"has_pressure_change": True},
                     "shell": {"has_pressure_change": True},
                     "has_holdup": True},
        "desuperheat": {"dynamic": False},
        "cooling": {"dynamic": False, "has_holdup": False},
        "property_package": model.fs.properties})
    model.discretizer = pyo.TransformationFactory('dae.finite_difference')
    model.discretizer.apply_to(model,
                               nfe=2,
                               wrt=model.fs.time,
                               scheme="BACKWARD")
    model.fs.fwh.set_initial_condition()
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
    model.fs.fwh.condense.area.fix(600)
    model.fs.fwh.condense.overall_heat_transfer_coefficient.fix(1010)
    model.fs.fwh.desuperheat.area.fix(85)
    model.fs.fwh.desuperheat.overall_heat_transfer_coefficient.fix(145)
    model.fs.fwh.cooling.area.fix(100)
    model.fs.fwh.cooling.overall_heat_transfer_coefficient.fix(675)
    model.fs.fwh.condense.tube.deltaP[:].fix(0)
    model.fs.fwh.condense.level.fix(0.275)
    model.fs.fwh.condense.heater_diameter.fix(1.4)
    model.fs.fwh.condense.vol_frac_shell.fix(0.675)
    model.fs.fwh.condense.cond_sect_length.fix(6.4)
    model.fs.fwh.condense.tube.volume.fix(1.3)
    model.fs.fwh.initialize(optarg={"max_iter": 50})
    solver = pyo.SolverFactory("ipopt")
    solver.options = {
            "tol": 1e-7,
            "linear_solver": "ma27",
            "max_iter": 50,
    }
    # initial drain flow rate
    drain_flow0 = model.fs.fwh.cooling.outlet_1.flow_mol[0].value
    assert(abs(drain_flow0 - 67.2061) < 0.0001)
    # change drain flow rate to increase water level
    for t in model.fs.time:
        if t >= 30:
            model.fs.fwh.cooling.outlet_1.flow_mol[t].fix(drain_flow0*0.95)
    model.fs.fwh.condense.level.unfix()
    model.fs.fwh.condense.level[0].fix()
    solver.solve(model, tee=True)
    assert(degrees_of_freedom(model) == 0)
    assert(abs(model.fs.fwh.condense.level[60].value - 0.2759) < 0.0001)
