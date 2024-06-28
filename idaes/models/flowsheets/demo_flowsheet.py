#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Demonstration flowsheet for testing purposes.

Constructs a basic flowsheet for a benzene-toluene system with a mixer, heater
and flash unit.
"""
from pyomo.environ import ConcreteModel, TransformationFactory
from pyomo.network import Arc

from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
import idaes.core.util.scaling as iscale

from idaes.core.util.initialization import propagate_state
from idaes.models.properties.activity_coeff_models.BTX_activity_coeff_VLE import (
    BTXParameterBlock,
)
from idaes.models.unit_models import Mixer, Heater, Flash

import idaes.logger as idaeslog


def build_flowsheet():
    """Build demo flowsheet"""
    m = ConcreteModel()

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.BT_props = BTXParameterBlock()

    m.fs.M01 = Mixer(property_package=m.fs.BT_props)

    m.fs.H02 = Heater(property_package=m.fs.BT_props)

    m.fs.F03 = Flash(property_package=m.fs.BT_props)

    m.fs.s01 = Arc(source=m.fs.M01.outlet, destination=m.fs.H02.inlet)
    m.fs.s02 = Arc(source=m.fs.H02.outlet, destination=m.fs.F03.inlet)

    TransformationFactory("network.expand_arcs").apply_to(m.fs)

    return m


def set_scaling(m):
    """Set scaling for demo flowsheet"""
    params = m.fs.BT_props
    params.set_default_scaling("flow_mol", 1)
    params.set_default_scaling("flow_mol_phase", 1)
    params.set_default_scaling("flow_mol_phase_comp", 1)
    params.set_default_scaling("mole_frac_comp", 1)
    params.set_default_scaling("mole_frac_phase_comp", 1)
    params.set_default_scaling("temperature", 1e-1)
    params.set_default_scaling("pressure", 1e-5)
    params.set_default_scaling("pressure_sat_comp", 1e-5, index="benzene")
    params.set_default_scaling("pressure_sat_comp", 1e-4, index="toluene")
    params.set_default_scaling("enth_mol_phase_comp", 1e-4)

    # Mixer M01
    iscale.set_scaling_factor(m.fs.M01.inlet_1_state[0].flow_mol_phase["Liq"], 1e6)
    iscale.set_scaling_factor(m.fs.M01.inlet_1_state[0].mole_frac_comp["toluene"], 1e5)
    iscale.set_scaling_factor(
        m.fs.M01.inlet_1_state[0].mole_frac_phase_comp["Liq", "toluene"], 1e5
    )
    iscale.set_scaling_factor(
        m.fs.M01.inlet_1_state[0].mole_frac_phase_comp["Vap", "toluene"], 1e6
    )
    iscale.set_scaling_factor(m.fs.M01.inlet_2_state[0].mole_frac_comp["benzene"], 1e5)
    iscale.set_scaling_factor(
        m.fs.M01.inlet_2_state[0].mole_frac_phase_comp["Liq", "benzene"], 1e5
    )
    iscale.set_scaling_factor(
        m.fs.M01.inlet_2_state[0].mole_frac_phase_comp["Vap", "benzene"], 1e5
    )

    # Heater H02
    iscale.set_scaling_factor(m.fs.H02.control_volume.heat[0], 1e-5)

    # F03
    iscale.set_scaling_factor(m.fs.F03.control_volume.heat[0], 1)

    iscale.calculate_scaling_factors(m)


def set_dof(m):
    """Set degrees of freedom for demo flowsheet"""
    eps = 1e-5
    m.fs.M01.inlet_1.flow_mol.fix(1.0)
    m.fs.M01.inlet_1.mole_frac_comp[:, "benzene"].fix(1 - eps)
    m.fs.M01.inlet_1.mole_frac_comp[:, "toluene"].fix(eps)
    m.fs.M01.inlet_1.pressure.fix(101325)
    m.fs.M01.inlet_1.temperature.fix(370)

    m.fs.M01.inlet_2.flow_mol.fix(1.0)
    m.fs.M01.inlet_2.mole_frac_comp[:, "benzene"].fix(eps)
    m.fs.M01.inlet_2.mole_frac_comp[:, "toluene"].fix(1 - eps)
    m.fs.M01.inlet_2.pressure.fix(1.3e5)
    m.fs.M01.inlet_2.temperature.fix(380)

    m.fs.H02.outlet.temperature[0].fix(370)

    m.fs.F03.heat_duty.fix(1e-6)
    m.fs.F03.deltaP.fix(1e-6)


def initialize_flowsheet(m):
    """Initialize demo flowsheet"""
    m.fs.M01.initialize(outlvl=idaeslog.WARNING)

    propagate_state(m.fs.s01)
    m.fs.H02.initialize(outlvl=idaeslog.WARNING)
    propagate_state(m.fs.s02)
    m.fs.F03.initialize(outlvl=idaeslog.WARNING)


def solve_flowsheet(m, stee=False):
    """Solve demo flowsheet"""
    solver = get_solver("ipopt_v2")
    solver.solve(m, tee=stee)


def display_results(m):
    """Display results for flowsheet"""
    m.fs.M01.outlet.display()
    m.fs.H02.outlet.display()
    m.fs.F03.vap_outlet.display()
    m.fs.F03.liq_outlet.display()


if __name__ == "__main__":
    m = build_flowsheet()

    set_scaling(m)

    set_dof(m)

    initialize_flowsheet(m)

    solve_flowsheet(m, stee=True)

    display_results(m)
