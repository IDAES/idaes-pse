##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
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
Test for Cappresse's module for NMPC.
"""

import pytest
from pyomo.environ import (Block, ConcreteModel,  Constraint, Expression,
                           Set, SolverFactory, Var, value, 
                           TransformationFactory, TerminationCondition)

from idaes.core import (FlowsheetBlock, MaterialBalanceType, EnergyBalanceType,
        MomentumBalanceType)
from idaes.core.util.testing import (PhysicalParameterTestBlock,
        AqueousEnzymeParameterBlock, EnzymeReactionParameterBlock,
        EnzymeReactionBlock)
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.generic_models.unit_models import CSTR, Mixer, MomentumMixingType
from idaes.core.util.exceptions import ConfigurationError
from idaes.dynamic.cappresse import nmpc
import pdb

__author__ = "Robert Parker"


# See if ipopt is available and set up solver
if SolverFactory('ipopt').available():
    solver = SolverFactory('ipopt')
    solver.options = {'tol': 1e-6,
                      'mu_init': 1e-8,
                      'bound_push': 1e-8}
else:
    solver = None


def make_model():
    horizon = 6 # minutes
    time_set = [0, horizon]
    ntfe = 60 # fe spacing: 6 seconds
    ntcp = 2
    sample_time = 0.5 # 5 fe per sample

    m = ConcreteModel(name='CSTR model for testing')
    m.fs = FlowsheetBlock(default={'dynamic': True,
                                   'time_set': time_set})

    m.fs.properties = AqueousEnzymeParameterBlock()
    m.fs.reactions = EnzymeReactionParameterBlock(
            default={'property_package': m.fs.properties})
    m.fs.cstr = CSTR(default={'property_package': m.fs.properties,
                              'reaction_package': m.fs.reactions,
                              'material_balance_type': MaterialBalanceType.componentTotal,
                              'energy_balance_type': EnergyBalanceType.enthalpyTotal,
                              'momentum_balance_type': MomentumBalanceType.none,
                              'has_heat_of_reaction': True})

    m.fs.mixer = Mixer(default={
        'property_package': m.fs.properties,
        'material_balance_type': MaterialBalanceType.componentTotal,
        'momentum_mixing_type': MomentumMixingType.none,
        'num_inlets': 2,
        'inlet_list': ['S_inlet', 'E_inlet']})
    # Allegedly the proper energy balance is being used...

    # Time discretization
    disc = TransformationFactory('dae.collocation')
    disc.apply_to(m, wrt=m.fs.time, nfe=ntfe, ncp=ntcp, scheme='LAGRANGE-RADAU')

    # Fix geometry variables
    m.fs.cstr.volume.fix(1.0)

    # Fix initial conditions:
    for p, j in m.fs.properties.phase_list*m.fs.properties.component_list:
        m.fs.cstr.control_volume.material_holdup[0, p, j].fix(0)

    m.fs.mixer.E_inlet.conc_mol.fix(0)
    m.fs.mixer.S_inlet.conc_mol.fix(0)

    for t, j in m.fs.time*m.fs.properties.component_list:
        if t <= 2:
            if j == 'E':
                m.fs.mixer.E_inlet.conc_mol[t, j].fix(11.91)
            elif j == 'S':
                m.fs.mixer.S_inlet.conc_mol[t, j].fix(12.92)
        elif t <= 4:
            if j == 'E':
                m.fs.mixer.E_inlet.conc_mol[t, j].fix(5.95)
            elif j == 'S':
                m.fs.mixer.S_inlet.conc_mol[t, j].fix(12.92)
        else:
            if j == 'E':
                m.fs.mixer.E_inlet.conc_mol[t, j].fix(8.95)
            elif j == 'S':
                m.fs.mixer.S_inlet.conc_mol[t, j].fix(16.75)

    m.fs.mixer.E_inlet.fix(0.1)
    m.fs.mixer.S_inlet.fix(2.1)

    m.fs.mixer.E_inlet.temperature.fix(290)
    m.fs.mixer.S_inlet.temperature.fix(310)

    pdb.set_trace()


@pytest.fixture
def model():
    pass


@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_something():
    pass


if __name__ == '__main__':
    make_model()
