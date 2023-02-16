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
    m.fs.BT_props.set_default_scaling(
        "flow_mol_phase_comp", 1e2, index=("Liq", "benzene")
    )
    m.fs.BT_props.set_default_scaling(
        "flow_mol_phase_comp", 1e2, index=("Liq", "toluene")
    )

    iscale.set_scaling_factor(m.fs.M01.inlet_1_state[0].mole_frac_comp["benzene"], 1)
    iscale.set_scaling_factor(m.fs.M01.inlet_1_state[0].mole_frac_comp["toluene"], 1e4)
    iscale.set_scaling_factor(
        m.fs.M01.inlet_1_state[0].enth_mol_phase_comp["Liq", "benzene"], 1e-4
    )
    iscale.set_scaling_factor(
        m.fs.M01.inlet_1_state[0].enth_mol_phase_comp["Vap", "benzene"], 1e-4
    )
    iscale.set_scaling_factor(
        m.fs.M01.inlet_1_state[0].enth_mol_phase_comp["Liq", "toluene"], 1e-4
    )
    iscale.set_scaling_factor(
        m.fs.M01.inlet_1_state[0].enth_mol_phase_comp["Vap", "toluene"], 1e-4
    )
    iscale.set_scaling_factor(m.fs.M01.inlet_1_state[0].eq_total, 1)
    iscale.set_scaling_factor(m.fs.M01.inlet_1_state[0].eq_sum_mol_frac, 1)
    iscale.set_scaling_factor(m.fs.M01.inlet_1_state[0].eq_comp["benzene"], 1)
    iscale.set_scaling_factor(m.fs.M01.inlet_1_state[0].eq_comp["toluene"], 1)
    iscale.set_scaling_factor(
        m.fs.M01.inlet_1_state[0].eq_phase_equilibrium["benzene"], 1
    )
    iscale.set_scaling_factor(
        m.fs.M01.inlet_1_state[0].eq_phase_equilibrium["toluene"], 1
    )
    iscale.set_scaling_factor(m.fs.M01.inlet_1_state[0].eq_P_vap["benzene"], 1)
    iscale.set_scaling_factor(m.fs.M01.inlet_1_state[0].eq_P_vap["toluene"], 1)

    iscale.set_scaling_factor(m.fs.M01.inlet_2_state[0].mole_frac_comp["benzene"], 1e4)
    iscale.set_scaling_factor(m.fs.M01.inlet_2_state[0].mole_frac_comp["toluene"], 1)
    iscale.set_scaling_factor(
        m.fs.M01.inlet_2_state[0].enth_mol_phase_comp["Liq", "benzene"], 1e-4
    )
    iscale.set_scaling_factor(
        m.fs.M01.inlet_2_state[0].enth_mol_phase_comp["Vap", "benzene"], 1e-4
    )
    iscale.set_scaling_factor(
        m.fs.M01.inlet_2_state[0].enth_mol_phase_comp["Liq", "toluene"], 1e-4
    )
    iscale.set_scaling_factor(
        m.fs.M01.inlet_2_state[0].enth_mol_phase_comp["Vap", "toluene"], 1e-4
    )
    iscale.set_scaling_factor(m.fs.M01.inlet_2_state[0].eq_total, 1)
    iscale.set_scaling_factor(m.fs.M01.inlet_2_state[0].eq_sum_mol_frac, 1)
    iscale.set_scaling_factor(m.fs.M01.inlet_2_state[0].eq_comp["benzene"], 1)
    iscale.set_scaling_factor(m.fs.M01.inlet_2_state[0].eq_comp["toluene"], 1)
    iscale.set_scaling_factor(
        m.fs.M01.inlet_2_state[0].eq_phase_equilibrium["benzene"], 1
    )
    iscale.set_scaling_factor(
        m.fs.M01.inlet_2_state[0].eq_phase_equilibrium["toluene"], 1
    )
    iscale.set_scaling_factor(m.fs.M01.inlet_2_state[0].eq_P_vap["benzene"], 1)
    iscale.set_scaling_factor(m.fs.M01.inlet_2_state[0].eq_P_vap["toluene"], 1)

    iscale.set_scaling_factor(
        m.fs.M01.mixed_state[0].enth_mol_phase_comp["Liq", "benzene"], 1e-4
    )
    iscale.set_scaling_factor(
        m.fs.M01.mixed_state[0].enth_mol_phase_comp["Vap", "benzene"], 1e-4
    )
    iscale.set_scaling_factor(
        m.fs.M01.mixed_state[0].enth_mol_phase_comp["Liq", "toluene"], 1e-4
    )
    iscale.set_scaling_factor(
        m.fs.M01.mixed_state[0].enth_mol_phase_comp["Vap", "toluene"], 1e-4
    )
    iscale.set_scaling_factor(m.fs.M01.mixed_state[0].enthalpy_flow_terms["Vap"], 1)
    iscale.set_scaling_factor(m.fs.M01.mixed_state[0].eq_total, 1)
    iscale.set_scaling_factor(m.fs.M01.mixed_state[0].eq_sum_mol_frac, 1)
    iscale.set_scaling_factor(m.fs.M01.mixed_state[0].eq_mol_frac_out, 1)
    iscale.set_scaling_factor(m.fs.M01.mixed_state[0].eq_comp["benzene"], 1)
    iscale.set_scaling_factor(m.fs.M01.mixed_state[0].eq_comp["toluene"], 1)
    iscale.set_scaling_factor(
        m.fs.M01.mixed_state[0].eq_phase_equilibrium["benzene"], 1
    )
    iscale.set_scaling_factor(
        m.fs.M01.mixed_state[0].eq_phase_equilibrium["toluene"], 1
    )
    iscale.set_scaling_factor(m.fs.M01.mixed_state[0].eq_P_vap["benzene"], 1)
    iscale.set_scaling_factor(m.fs.M01.mixed_state[0].eq_P_vap["toluene"], 1)

    iscale.set_scaling_factor(m.fs.M01.mixed_state[0].pressure, 1e-5)
    iscale.set_scaling_factor(m.fs.M01.mixed_state[0].temperature, 1e-2)
    iscale.set_scaling_factor(m.fs.M01.mixed_state[0].flow_mol, 1)
    iscale.set_scaling_factor(m.fs.M01.mixed_state[0].mole_frac_comp["benzene"], 1e1)
    iscale.set_scaling_factor(m.fs.M01.mixed_state[0].mole_frac_comp["toluene"], 1e1)

    iscale.set_scaling_factor(m.fs.H02.control_volume.heat, 1e-5)
    iscale.set_scaling_factor(m.fs.H02.control_volume.properties_in[0].pressure, 1e-5)
    iscale.set_scaling_factor(
        m.fs.H02.control_volume.properties_in[0].temperature, 1e-2
    )
    iscale.set_scaling_factor(m.fs.H02.control_volume.properties_in[0].flow_mol, 1)
    iscale.set_scaling_factor(
        m.fs.H02.control_volume.properties_in[0].mole_frac_comp["benzene"], 1e1
    )
    iscale.set_scaling_factor(
        m.fs.H02.control_volume.properties_in[0].mole_frac_comp["toluene"], 1e1
    )
    iscale.set_scaling_factor(m.fs.H02.control_volume.properties_out[0].flow_mol, 1)
    iscale.set_scaling_factor(
        m.fs.H02.control_volume.properties_out[0].mole_frac_comp["benzene"], 1e1
    )
    iscale.set_scaling_factor(
        m.fs.H02.control_volume.properties_out[0].mole_frac_comp["toluene"], 1e1
    )
    iscale.set_scaling_factor(
        m.fs.H02.control_volume.properties_out[0].temperature, 1e-2
    )
    iscale.set_scaling_factor(m.fs.H02.control_volume.properties_out[0].pressure, 1e-5)
    iscale.set_scaling_factor(
        m.fs.H02.control_volume.properties_in[0].enth_mol_phase_comp["Liq", "benzene"],
        1e-4,
    )
    iscale.set_scaling_factor(
        m.fs.H02.control_volume.properties_in[0].enth_mol_phase_comp["Vap", "benzene"],
        1e-4,
    )
    iscale.set_scaling_factor(
        m.fs.H02.control_volume.properties_in[0].enth_mol_phase_comp["Liq", "toluene"],
        1e-4,
    )
    iscale.set_scaling_factor(
        m.fs.H02.control_volume.properties_in[0].enth_mol_phase_comp["Vap", "toluene"],
        1e-4,
    )
    iscale.set_scaling_factor(
        m.fs.H02.control_volume.properties_in[0].eq_comp["benzene"], 1
    )
    iscale.set_scaling_factor(
        m.fs.H02.control_volume.properties_in[0].eq_comp["toluene"], 1
    )
    iscale.set_scaling_factor(
        m.fs.H02.control_volume.properties_out[0].eq_comp["benzene"], 1
    )
    iscale.set_scaling_factor(
        m.fs.H02.control_volume.properties_out[0].eq_comp["toluene"], 1
    )
    iscale.set_scaling_factor(
        m.fs.H02.control_volume.properties_in[0].eq_phase_equilibrium["benzene"], 1
    )
    iscale.set_scaling_factor(
        m.fs.H02.control_volume.properties_in[0].eq_phase_equilibrium["toluene"], 1
    )
    iscale.set_scaling_factor(
        m.fs.H02.control_volume.properties_out[0].eq_phase_equilibrium["benzene"], 1
    )
    iscale.set_scaling_factor(
        m.fs.H02.control_volume.properties_out[0].eq_phase_equilibrium["toluene"], 1
    )
    iscale.set_scaling_factor(
        m.fs.H02.control_volume.properties_in[0.0].eq_P_vap["benzene"], 1
    )
    iscale.set_scaling_factor(
        m.fs.H02.control_volume.properties_in[0.0].eq_P_vap["toluene"], 1
    )
    iscale.set_scaling_factor(
        m.fs.H02.control_volume.properties_out[0.0].eq_P_vap["benzene"], 1
    )
    iscale.set_scaling_factor(
        m.fs.H02.control_volume.properties_out[0.0].eq_P_vap["toluene"], 1
    )
    iscale.set_scaling_factor(
        m.fs.H02.control_volume.properties_out[0].enth_mol_phase_comp["Liq", "benzene"],
        1e-4,
    )
    iscale.set_scaling_factor(
        m.fs.H02.control_volume.properties_out[0].enth_mol_phase_comp["Vap", "benzene"],
        1e-4,
    )
    iscale.set_scaling_factor(
        m.fs.H02.control_volume.properties_out[0].enth_mol_phase_comp["Liq", "toluene"],
        1e-4,
    )
    iscale.set_scaling_factor(
        m.fs.H02.control_volume.properties_out[0].enth_mol_phase_comp["Vap", "toluene"],
        1e-4,
    )
    iscale.set_scaling_factor(m.fs.H02.control_volume.properties_in[0.0].eq_total, 1)
    iscale.set_scaling_factor(
        m.fs.H02.control_volume.properties_in[0.0].eq_sum_mol_frac, 1
    )
    iscale.set_scaling_factor(m.fs.H02.control_volume.properties_out[0.0].eq_total, 1)
    iscale.set_scaling_factor(
        m.fs.H02.control_volume.properties_out[0.0].eq_sum_mol_frac, 1
    )
    iscale.set_scaling_factor(
        m.fs.H02.control_volume.properties_out[0.0].eq_mol_frac_out, 1
    )

    iscale.set_scaling_factor(
        m.fs.H02.control_volume.properties_in[0.0].material_flow_terms[
            "Liq", "benzene"
        ],
        1,
    )
    iscale.set_scaling_factor(
        m.fs.H02.control_volume.properties_in[0.0].material_flow_terms[
            "Vap", "benzene"
        ],
        1,
    )
    iscale.set_scaling_factor(
        m.fs.H02.control_volume.properties_in[0.0].material_flow_terms[
            "Liq", "toluene"
        ],
        1,
    )
    iscale.set_scaling_factor(
        m.fs.H02.control_volume.properties_in[0.0].material_flow_terms[
            "Vap", "toluene"
        ],
        1,
    )
    iscale.set_scaling_factor(
        m.fs.H02.control_volume.properties_in[0.0].enthalpy_flow_terms["Liq"], 1
    )
    iscale.set_scaling_factor(
        m.fs.H02.control_volume.properties_in[0.0].enthalpy_flow_terms["Vap"], 1
    )

    iscale.set_scaling_factor(m.fs.F03.control_volume.heat, 1)
    iscale.set_scaling_factor(m.fs.F03.control_volume.properties_in[0].pressure, 1e-5)
    iscale.set_scaling_factor(
        m.fs.F03.control_volume.properties_in[0].temperature, 1e-2
    )
    iscale.set_scaling_factor(m.fs.F03.control_volume.properties_in[0].flow_mol, 1)
    iscale.set_scaling_factor(
        m.fs.F03.control_volume.properties_in[0].mole_frac_comp["benzene"], 1e1
    )
    iscale.set_scaling_factor(
        m.fs.F03.control_volume.properties_in[0].mole_frac_comp["toluene"], 1e1
    )

    iscale.set_scaling_factor(m.fs.F03.control_volume.properties_in[0.0].eq_total, 1)
    iscale.set_scaling_factor(
        m.fs.F03.control_volume.properties_in[0.0].eq_sum_mol_frac, 1
    )
    iscale.set_scaling_factor(m.fs.F03.control_volume.properties_out[0.0].eq_total, 1)
    iscale.set_scaling_factor(
        m.fs.F03.control_volume.properties_out[0.0].eq_sum_mol_frac, 1
    )
    iscale.set_scaling_factor(
        m.fs.F03.control_volume.properties_out[0.0].eq_mol_frac_out, 1
    )

    iscale.set_scaling_factor(
        m.fs.F03.control_volume.properties_in[0].enth_mol_phase_comp["Liq", "benzene"],
        1e-4,
    )
    iscale.set_scaling_factor(
        m.fs.F03.control_volume.properties_in[0].enth_mol_phase_comp["Vap", "benzene"],
        1e-4,
    )
    iscale.set_scaling_factor(
        m.fs.F03.control_volume.properties_in[0].enth_mol_phase_comp["Liq", "toluene"],
        1e-4,
    )
    iscale.set_scaling_factor(
        m.fs.F03.control_volume.properties_in[0].enth_mol_phase_comp["Vap", "toluene"],
        1e-4,
    )
    iscale.set_scaling_factor(
        m.fs.F03.control_volume.properties_out[0].enth_mol_phase_comp["Liq", "benzene"],
        1e-4,
    )
    iscale.set_scaling_factor(
        m.fs.F03.control_volume.properties_out[0].enth_mol_phase_comp["Vap", "benzene"],
        1e-4,
    )
    iscale.set_scaling_factor(
        m.fs.F03.control_volume.properties_out[0].enth_mol_phase_comp["Liq", "toluene"],
        1e-4,
    )
    iscale.set_scaling_factor(
        m.fs.F03.control_volume.properties_out[0].enth_mol_phase_comp["Vap", "toluene"],
        1e-4,
    )

    iscale.set_scaling_factor(
        m.fs.F03.control_volume.properties_in[0].eq_comp["benzene"], 1
    )
    iscale.set_scaling_factor(
        m.fs.F03.control_volume.properties_in[0].eq_comp["toluene"], 1
    )
    iscale.set_scaling_factor(
        m.fs.F03.control_volume.properties_out[0].eq_comp["benzene"], 1
    )
    iscale.set_scaling_factor(
        m.fs.F03.control_volume.properties_out[0].eq_comp["toluene"], 1
    )
    iscale.set_scaling_factor(
        m.fs.F03.control_volume.properties_in[0].eq_phase_equilibrium["benzene"], 1
    )
    iscale.set_scaling_factor(
        m.fs.F03.control_volume.properties_in[0].eq_phase_equilibrium["toluene"], 1
    )
    iscale.set_scaling_factor(
        m.fs.F03.control_volume.properties_out[0].eq_phase_equilibrium["benzene"], 1
    )
    iscale.set_scaling_factor(
        m.fs.F03.control_volume.properties_out[0].eq_phase_equilibrium["toluene"], 1
    )
    iscale.set_scaling_factor(
        m.fs.F03.control_volume.properties_in[0.0].eq_P_vap["benzene"], 1
    )
    iscale.set_scaling_factor(
        m.fs.F03.control_volume.properties_in[0.0].eq_P_vap["toluene"], 1
    )
    iscale.set_scaling_factor(
        m.fs.F03.control_volume.properties_out[0.0].eq_P_vap["benzene"], 1
    )
    iscale.set_scaling_factor(
        m.fs.F03.control_volume.properties_out[0.0].eq_P_vap["toluene"], 1
    )
    iscale.set_scaling_factor(
        m.fs.F03.control_volume.properties_in[0.0].material_flow_terms[
            "Liq", "benzene"
        ],
        1,
    )
    iscale.set_scaling_factor(
        m.fs.F03.control_volume.properties_in[0.0].material_flow_terms[
            "Vap", "benzene"
        ],
        1,
    )
    iscale.set_scaling_factor(
        m.fs.F03.control_volume.properties_in[0.0].material_flow_terms[
            "Liq", "toluene"
        ],
        1,
    )
    iscale.set_scaling_factor(
        m.fs.F03.control_volume.properties_in[0.0].material_flow_terms[
            "Vap", "toluene"
        ],
        1,
    )
    iscale.set_scaling_factor(
        m.fs.F03.control_volume.properties_in[0.0].enthalpy_flow_terms["Liq"], 1
    )
    iscale.set_scaling_factor(
        m.fs.F03.control_volume.properties_in[0.0].enthalpy_flow_terms["Vap"], 1
    )

    iscale.calculate_scaling_factors(m)


def set_dof(m):
    m.fs.M01.inlet_1.flow_mol.fix(1.0)
    m.fs.M01.inlet_1.mole_frac_comp[:, "benzene"].fix(1)
    m.fs.M01.inlet_1.mole_frac_comp[:, "toluene"].fix(1e-5)
    m.fs.M01.inlet_1.pressure.fix(101325)
    m.fs.M01.inlet_1.temperature.fix(370)

    m.fs.M01.inlet_2.flow_mol.fix(1.0)
    m.fs.M01.inlet_2.mole_frac_comp[:, "benzene"].fix(1e-5)
    m.fs.M01.inlet_2.mole_frac_comp[:, "toluene"].fix(1)
    m.fs.M01.inlet_2.pressure.fix(1.3e5)
    m.fs.M01.inlet_2.temperature.fix(380)

    m.fs.H02.outlet.temperature[0].fix(370)

    m.fs.F03.heat_duty.fix(0)
    m.fs.F03.deltaP.fix(0)


def initialize_flowsheet(m):
    m.fs.M01.initialize(outlvl=idaeslog.WARNING)

    propagate_state(m.fs.s01)
    m.fs.H02.initialize(outlvl=idaeslog.WARNING)
    propagate_state(m.fs.s02)
    m.fs.F03.initialize(outlvl=idaeslog.WARNING)


def solve_flowsheet(m, stee=False):
    solver = get_solver()
    solver.solve(m, tee=stee)


def display_results(m):
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
