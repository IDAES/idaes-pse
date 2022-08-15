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
Test for Cappresse's module for NMPC.
"""

from pyomo.environ import (
    ConcreteModel,
    SolverFactory,
    units as pyunits,
    TransformationFactory,
)
from pyomo.network import Arc

from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
)
from idaes.core.util.tests.test_initialization import (
    AqueousEnzymeParameterBlock,
    EnzymeReactionParameterBlock,
)
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.unit_models import CSTR, Mixer, MomentumMixingType
from idaes.apps.caprese import nmpc
from idaes.apps.caprese.nmpc import *

__author__ = "Robert Parker"


# See if ipopt is available and set up solver
if SolverFactory("ipopt").available():
    solver = SolverFactory("ipopt")
    solver.options = {"tol": 1e-6, "mu_init": 1e-8, "bound_push": 1e-8}
else:
    solver = None


def make_model(
    horizon=6, ntfe=60, ntcp=2, inlet_E=11.91, inlet_S=12.92, steady=False, bounds=False
):
    time_set = [0, horizon]

    m = ConcreteModel(name="CSTR model for testing")
    if steady:
        m.fs = FlowsheetBlock(dynamic=False)
    else:
        m.fs = FlowsheetBlock(
            dynamic=True,
            time_set=time_set,
            time_units=pyunits.minute,
        )

    m.fs.properties = AqueousEnzymeParameterBlock()
    m.fs.reactions = EnzymeReactionParameterBlock(property_package=m.fs.properties)
    m.fs.cstr = CSTR(
        has_holdup=True,
        property_package=m.fs.properties,
        reaction_package=m.fs.reactions,
        material_balance_type=MaterialBalanceType.componentTotal,
        energy_balance_type=EnergyBalanceType.enthalpyTotal,
        momentum_balance_type=MomentumBalanceType.none,
        has_heat_of_reaction=True,
    )

    m.fs.mixer = Mixer(
        property_package=m.fs.properties,
        material_balance_type=MaterialBalanceType.componentTotal,
        momentum_mixing_type=MomentumMixingType.none,
        num_inlets=2,
        inlet_list=["S_inlet", "E_inlet"],
    )
    # Allegedly the proper energy balance is being used...

    # Time discretization
    if not steady:
        disc = TransformationFactory("dae.collocation")
        disc.apply_to(m, wrt=m.fs.time, nfe=ntfe, ncp=ntcp, scheme="LAGRANGE-RADAU")

    # Fix geometry variables
    m.fs.cstr.volume[0].fix(1.0)

    # Fix initial conditions:
    if not steady:
        for p, j in m.fs.properties.phase_list * m.fs.properties.component_list:
            if j == "Solvent":
                continue
            m.fs.cstr.control_volume.material_holdup[0, p, j].fix(0.001)
    # Note: Model does not solve when initial conditions are empty tank

    m.fs.mixer.E_inlet.conc_mol.fix(0)
    m.fs.mixer.S_inlet.conc_mol.fix(0)

    for t, j in m.fs.time * m.fs.properties.component_list:
        if j == "E":
            m.fs.mixer.E_inlet.conc_mol[t, j].fix(inlet_E)
        elif j == "S":
            m.fs.mixer.S_inlet.conc_mol[t, j].fix(inlet_S)

    m.fs.mixer.E_inlet.flow_vol.fix(0.1)
    m.fs.mixer.S_inlet.flow_vol.fix(2.1)
    m.fs.mixer.E_inlet.conc_mol[:, "Solvent"].fix(1.0)
    m.fs.mixer.S_inlet.conc_mol[:, "Solvent"].fix(1.0)

    m.fs.mixer.E_inlet.temperature.fix(290)
    m.fs.mixer.S_inlet.temperature.fix(310)

    m.fs.inlet = Arc(source=m.fs.mixer.outlet, destination=m.fs.cstr.inlet)

    if not steady:
        # This constraint is in lieu of tracking the CSTR's level and allowing
        # the outlet flow rate to be another degree of freedom.
        # ^ Not sure how to do this in IDAES.
        @m.fs.cstr.Constraint(m.fs.time, doc="Total flow rate balance")
        def total_flow_balance(cstr, t):
            return cstr.inlet.flow_vol[t] == cstr.outlet.flow_vol[t]

    # Specify initial condition for energy
    if not steady:
        m.fs.cstr.control_volume.energy_holdup[m.fs.time.first(), "aq"].fix(300)

    TransformationFactory("network.expand_arcs").apply_to(m.fs)

    if bounds:
        m.fs.mixer.E_inlet.flow_vol.setlb(0.01)
        m.fs.mixer.E_inlet.flow_vol.setub(1.0)
        m.fs.mixer.S_inlet.flow_vol.setlb(0.5)
        m.fs.mixer.S_inlet.flow_vol.setub(5.0)

        m.fs.cstr.control_volume.material_holdup.setlb(0)
        holdup = m.fs.cstr.control_volume.material_holdup
        for t in m.fs.time:
            holdup[t, "aq", "S"].setub(20)
            holdup[t, "aq", "E"].setub(1)
            holdup[t, "aq", "P"].setub(5)
            holdup[t, "aq", "C"].setub(5)

    return m


if __name__ == "__main__":
    m_plant = make_model()
    m_controller = make_model(bounds=True)
    m_ss = make_model(steady=True)

    assert degrees_of_freedom(m_plant) == 0
    assert degrees_of_freedom(m_ss) == 0
