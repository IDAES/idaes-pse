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
Tests for turbine multistage model.

Author: John Eslick
"""
import pytest

from pyomo.environ import (ConcreteModel, SolverFactory, TransformationFactory,
                           Constraint, value, units as pyunits)
from pyomo.network import Arc
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.generic_models.unit_models import Heater
from idaes.power_generation.unit_models import TurbineMultistage
from idaes.generic_models.properties import iapws95
from idaes.core.util.model_statistics import (
        degrees_of_freedom,
        activated_equalities_generator)

prop_available = iapws95.iapws95_available()

# See if ipopt is available and set up solver
if SolverFactory('ipopt').available():
    solver = SolverFactory('ipopt')
    solver.options = {'tol': 1e-6}
else:
    solver = None


@pytest.fixture(scope="module")
def model():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = iapws95.Iapws95ParameterBlock()
    # roughly based on NETL baseline studies
    m.fs.turb = TurbineMultistage(default={
        "property_package": m.fs.properties,
        "num_hp": 7,
        "num_ip": 14,
        "num_lp": 11,
        "hp_split_locations": [4, 7],
        "ip_split_locations": [5, 14],
        "lp_split_locations": [4, 7, 9, 11],
        "hp_disconnect": [7],
        "ip_split_num_outlets": {14: 3}})

    # Add reheater
    m.fs.reheat = Heater(default={"property_package": m.fs.properties})
    m.fs.hp_to_reheat = Arc(source=m.fs.turb.hp_split[7].outlet_1,
                            destination=m.fs.reheat.inlet)
    m.fs.reheat_to_ip = Arc(source=m.fs.reheat.outlet,
                            destination=m.fs.turb.ip_stages[1].inlet)

    return m


@pytest.mark.integration
def test_unit_consistency(model):
    assert_units_consistent(model)


@pytest.mark.component
@pytest.mark.skipif(not prop_available, reason="IAPWS not available")
@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_initialize(model):
    """Make a turbine model and make sure it doesn't throw exception"""
    turb = model.fs.turb

    # Set the inlet of the turbine
    p = 2.4233e7
    hin = value(iapws95.htpx(T=880*pyunits.K, P=p*pyunits.Pa))
    model.fs.turb.inlet_split.inlet.enth_mol[0].fix(hin)
    model.fs.turb.inlet_split.inlet.flow_mol[0].fix(26000)
    model.fs.turb.inlet_split.inlet.pressure[0].fix(p)

    # Set the inlet of the ip section, which is disconnected
    # here to insert reheater
    p = 7.802e+06
    hin = value(iapws95.htpx(T=880*pyunits.K, P=p*pyunits.Pa))
    model.fs.turb.ip_stages[1].inlet.enth_mol[0].value = hin
    model.fs.turb.ip_stages[1].inlet.flow_mol[0].value = 25220.0
    model.fs.turb.ip_stages[1].inlet.pressure[0].value = p

    for i, s in turb.hp_stages.items():
        s.ratioP[:] = 0.88
        s.efficiency_isentropic[:] = 0.9
    for i, s in turb.ip_stages.items():
        s.ratioP[:] = 0.85
        s.efficiency_isentropic[:] = 0.9
    for i, s in turb.lp_stages.items():
        s.ratioP[:] = 0.82
        s.efficiency_isentropic[:] = 0.9

    turb.hp_split[4].split_fraction[0, "outlet_2"].fix(0.03)
    turb.hp_split[7].split_fraction[0, "outlet_2"].fix(0.03)
    turb.ip_split[5].split_fraction[0, "outlet_2"].fix(0.04)
    turb.ip_split[14].split_fraction[0, "outlet_2"].fix(0.04)
    turb.ip_split[14].split_fraction[0, "outlet_3"].fix(0.15)
    turb.lp_split[4].split_fraction[0, "outlet_2"].fix(0.04)
    turb.lp_split[7].split_fraction[0, "outlet_2"].fix(0.04)
    turb.lp_split[9].split_fraction[0, "outlet_2"].fix(0.04)
    turb.lp_split[11].split_fraction[0, "outlet_2"].fix(0.04)

    # Congiure with reheater for a full test
    turb.ip_stages[1].inlet.fix()
    turb.inlet_split.inlet.flow_mol.unfix()
    turb.inlet_mix.use_equal_pressure_constraint()
    turb.initialize(outlvl=1)
    turb.ip_stages[1].inlet.unfix()

    for t in model.fs.time:
        model.fs.reheat.inlet.flow_mol[t].value = \
            value(turb.hp_split[7].outlet_1_state[t].flow_mol)
        model.fs.reheat.inlet.enth_mol[t].value = \
            value(turb.hp_split[7].outlet_1_state[t].enth_mol)
        model.fs.reheat.inlet.pressure[t].value = \
            value(turb.hp_split[7].outlet_1_state[t].pressure)
    model.fs.reheat.initialize(outlvl=4)

    def reheat_T_rule(b, t):
        return model.fs.reheat.control_volume.properties_out[t].temperature == 880
    model.fs.reheat.temperature_out_equation = Constraint(
            model.fs.reheat.flowsheet().config.time,
            rule=reheat_T_rule)

    TransformationFactory("network.expand_arcs").apply_to(model)
    model.fs.turb.outlet_stage.control_volume.properties_out[0].pressure.fix()

    assert degrees_of_freedom(model) == 0
    solver.solve(model, tee=True)

    eq_cons = activated_equalities_generator(model)
    for c in eq_cons:
        assert(abs(c.body() - c.lower) < 1e-4)
