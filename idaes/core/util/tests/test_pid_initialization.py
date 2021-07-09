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
Test for element-by-element initialization on PID controller. This is important
as PID controllers have additional time-linking variables beyond derivative
and differential variables.
"""

import pytest
from pyomo.environ import (Block, ConcreteModel, Constraint, Expression,
                           Set, SolverFactory, Var, value, Param, Reals,
                           TransformationFactory, TerminationCondition,
                           exp)
from pyomo.network import Arc, Port
from pyomo.dae import DerivativeVar
from pyomo.common.collections import ComponentMap

from idaes.core import (FlowsheetBlock, 
                        MaterialBalanceType, 
                        EnergyBalanceType,
                        MomentumBalanceType, 
                        declare_process_block_class,
                        PhysicalParameterBlock,
                        StateBlock,
                        StateBlockData,
                        ReactionParameterBlock,
                        ReactionBlockBase,
                        ReactionBlockDataBase,
                        MaterialFlowBasis)
from idaes.core.util.testing import PhysicalParameterTestBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.generic_models.unit_models import CSTR, Mixer, MomentumMixingType
from idaes.generic_models.control import PIDBlock, PIDForm
from idaes.core.util.initialization import initialize_by_time_element
from idaes.core.util.tests.test_initialization import \
        AqueousEnzymeParameterBlock, EnzymeReactionParameterBlock
from idaes.core.util import get_solver

import idaes.logger as idaeslog

__author__ = "Robert Parker"


solver = get_solver()


def make_model(horizon=6, ntfe=60, ntcp=2, inlet_E=11.91, inlet_S=12.92):

    time_set = [0, horizon]

    m = ConcreteModel(name='CSTR with level control')
    m.fs = FlowsheetBlock(default={'dynamic': True,
                                   'time_set': time_set})

    m.fs.properties = AqueousEnzymeParameterBlock()
    m.fs.reactions = EnzymeReactionParameterBlock(
            default={'property_package': m.fs.properties})
    m.fs.cstr = CSTR(default={'has_holdup': True,
                              'property_package': m.fs.properties,
                              'reaction_package': m.fs.reactions,
                              'material_balance_type': MaterialBalanceType.componentTotal,
                              'energy_balance_type': EnergyBalanceType.enthalpyTotal,
                              'momentum_balance_type': MomentumBalanceType.none,
                              'has_heat_of_reaction': True})
    # MomentumBalanceType.none used because the property package doesn't
    # include pressure.

    m.fs.mixer = Mixer(default={
        'property_package': m.fs.properties,
        'material_balance_type': MaterialBalanceType.componentTotal,
        'momentum_mixing_type': MomentumMixingType.none,
        # MomentumMixingType.none used because the property package doesn't
        # include pressure.
        'num_inlets': 2,
       'inlet_list': ['S_inlet', 'E_inlet']})
    # Allegedly the proper energy balance is being used...

    # Time discretization
    disc = TransformationFactory('dae.collocation')
    disc.apply_to(m, wrt=m.fs.time, nfe=ntfe, ncp=ntcp, scheme='LAGRANGE-RADAU')

    m.fs.pid = PIDBlock(default={'pv': m.fs.cstr.volume,
                                 'output': m.fs.cstr.outlet.flow_vol,
                                 'upper': 5.0,
                                 'lower': 0.5,
                                 'calculate_initial_integral': True,
                                 # ^ Why would initial integral be calculated
                                 # to be nonzero?
                                 'pid_form': PIDForm.velocity})

    m.fs.pid.gain.fix(-1.0)
    m.fs.pid.time_i.fix(0.1)
    m.fs.pid.time_d.fix(0.0)
    m.fs.pid.setpoint.fix(1.0)

    # Fix initial condition for volume:
    m.fs.cstr.volume.unfix()
    m.fs.cstr.volume[m.fs.time.first()].fix(1.0)

    # Fix initial conditions for other variables:
    for p, j in m.fs.properties.phase_list*m.fs.properties.component_list:
        if j == 'Solvent':
            continue
        m.fs.cstr.control_volume.material_holdup[0, p, j].fix(0.001)
    # Note: Model does not solve when initial conditions are empty tank
    m.fs.cstr.control_volume.energy_holdup[m.fs.time.first(), 'aq'].fix(300)

    m.fs.mixer.E_inlet.conc_mol.fix(0)
    m.fs.mixer.S_inlet.conc_mol.fix(0)
    m.fs.mixer.E_inlet.conc_mol[:,'Solvent'].fix(1.)
    m.fs.mixer.S_inlet.conc_mol[:,'Solvent'].fix(1.)

    for t, j in m.fs.time*m.fs.properties.component_list:
        if j == 'E':
            m.fs.mixer.E_inlet.conc_mol[t, j].fix(inlet_E)
        elif j == 'S':
            m.fs.mixer.S_inlet.conc_mol[t, j].fix(inlet_S)

    m.fs.mixer.E_inlet.flow_vol.fix(0.1)
    m.fs.mixer.S_inlet.flow_vol.fix(2.1)

    # Specify a perturbation to substrate flow rate:
    for t in m.fs.time:
        if t < horizon/4:
            continue
        else:
            m.fs.mixer.S_inlet.flow_vol[t].fix(3.0)

    m.fs.mixer.E_inlet.temperature.fix(290)
    m.fs.mixer.S_inlet.temperature.fix(310)

    m.fs.inlet = Arc(source=m.fs.mixer.outlet, destination=m.fs.cstr.inlet)

    # Fix "initial condition" for outlet flow rate, as here it cannot be
    # specified by the PID controller
    m.fs.cstr.outlet.flow_vol[m.fs.time.first()].fix(2.2)

    TransformationFactory('network.expand_arcs').apply_to(m.fs)

    return m


@pytest.mark.component
@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_initialize():
    '''Very rough test, just to make sure degrees of freedom are not violated.
    '''
    mod = make_model(horizon=2, ntfe=20, ntcp=1, inlet_E=11.91, inlet_S=12.92)

    assert degrees_of_freedom(mod) == 0

    originally_active = ComponentMap([(comp, comp.active) 
        for comp in mod.component_data_objects((Block, Constraint))])
    originally_fixed = ComponentMap([(var, var.fixed)
        for var in mod.component_data_objects(Var)])

    initialize_by_time_element(mod.fs, mod.fs.time, solver=solver, 
            outlvl=idaeslog.DEBUG,
            fix_diff_only=False)
    assert degrees_of_freedom(mod) == 0

    for comp in mod.component_data_objects((Block, Constraint)):
        assert comp.active == originally_active[comp]
    for var in mod.component_data_objects(Var):
        assert var.fixed == originally_fixed[var]


if __name__ == '__main__':
    test_initialize()
