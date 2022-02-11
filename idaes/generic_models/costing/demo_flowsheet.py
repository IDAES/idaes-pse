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
Demonstration flowsheet for testing costing proposal
"""

import pyomo.environ as pyo
from pyomo.network import Arc
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.util import get_solver
from idaes.core.util.initialization import propagate_state
from idaes.generic_models.properties.activity_coeff_models.BTX_activity_coeff_VLE \
    import BTXParameterBlock
from idaes.generic_models.unit_models import Mixer, Heater, Flash

import idaes.logger as idaeslog

from idaes.costing.costing_base import \
    FlowsheetCostingBlock, register_unit

from demo_costing_package import MyCosting


def build_flowsheet():
    m = pyo.ConcreteModel()

    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.BT_props = BTXParameterBlock()

    m.fs.M01 = Mixer(default={"property_package": m.fs.BT_props})

    m.fs.H02 = Heater(default={"property_package": m.fs.BT_props})

    m.fs.F03 = Flash(default={"property_package": m.fs.BT_props})

    m.fs.s01 = Arc(source=m.fs.M01.outlet, destination=m.fs.H02.inlet)
    m.fs.s02 = Arc(source=m.fs.H02.outlet, destination=m.fs.F03.inlet)

    pyo.TransformationFactory("network.expand_arcs").apply_to(m.fs)

    m.fs.M01.inlet_1.flow_mol.fix(1.0)
    m.fs.M01.inlet_1.mole_frac_comp[:, "benzene"].fix(1.0)
    m.fs.M01.inlet_1.mole_frac_comp[:, "toluene"].fix(1e-5)
    m.fs.M01.inlet_1.pressure.fix(101325)
    m.fs.M01.inlet_1.temperature.fix(370)

    m.fs.M01.inlet_2.flow_mol.fix(1.0)
    m.fs.M01.inlet_2.mole_frac_comp[:, "benzene"].fix(1e-5)
    m.fs.M01.inlet_2.mole_frac_comp[:, "toluene"].fix(1.0)
    m.fs.M01.inlet_2.pressure.fix(1.3e5)
    m.fs.M01.inlet_2.temperature.fix(380)

    m.fs.H02.outlet.temperature[0].fix(370)

    m.fs.F03.heat_duty.fix(0)
    m.fs.F03.deltaP.fix(0)

    m.fs.M01.initialize(outlvl=idaeslog.WARNING)

    propagate_state(m.fs.s01)
    m.fs.H02.initialize(outlvl=idaeslog.WARNING)

    propagate_state(m.fs.s02)
    m.fs.F03.initialize(outlvl=idaeslog.WARNING)

    solver = get_solver()
    solver.solve(m, tee=False)

    m.fs.F03.vap_outlet.display()
    m.fs.F03.liq_outlet.display()

    # Add costing here
    # Add cost for heater - lets say cost = k1*heat_duty
    # Add cost for benzene feed
    # These should exercise basic code at least

    m.fs.costing = FlowsheetCostingBlock(
        default={"costing_package": MyCosting})

    assert m.fs.costing.k1.value == 42

    # TODO : For now, go with explicitly defining costing method
    m.fs.costing.cost_unit(m.fs.H02)

    # # Costing of material flows
    # First, need to put a lower bound on heat duty
    m.fs.H02.heat_duty[0].setlb(0)
    m.fs.costing.cost_flow(m.fs.H02.heat_duty[0], "steam")

    # We can also register new units and flows
    register_unit(["USD2015 = 0.75 * USD2010"])
    m.fs.costing.register_flow_type(
        "benzene_feed", 5*pyo.units.USD2015/pyo.units.mol)

    # Cost full process (incl. aggregating unit costs)
    m.fs.costing.cost_process()

    # Initialize costing
    m.fs.costing.initialize()

    # Check unit consistency
    assert_units_consistent(m)

    # Check some outputs
    m.fs.costing.aggregate_capital_cost.display()
    m.fs.costing.total_capital_cost.display()

    m.fs.costing.aggregate_fixed_operating_cost.display()
    m.fs.costing.aggregate_variable_operating_cost.display()

    m.fs.costing.aggregate_flow_costs.display()

    # We can also convert the units of the results as needed
    print(pyo.value(pyo.units.convert(m.fs.costing.total_capital_cost,
                                      to_units=pyo.units.USD2015)))


if __name__ == "__main__":
    m = build_flowsheet()
