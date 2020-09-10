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
from pyomo.network import Arc
from pyomo.common.collections import ComponentSet

from idaes.core import (FlowsheetBlock, MaterialBalanceType, EnergyBalanceType,
        MomentumBalanceType)
from idaes.core.util.model_statistics import (degrees_of_freedom, 
        activated_equalities_generator)
from idaes.core.util.initialization import initialize_by_time_element
from idaes.core.util.exceptions import ConfigurationError
from idaes.generic_models.unit_models import CSTR, Mixer, MomentumMixingType
from idaes.apps.caprese import nmpc
from idaes.apps.caprese.nmpc import *
from idaes.apps.caprese.examples.cstr_model import make_model
import idaes.logger as idaeslog

__author__ = "Robert Parker"


# See if ipopt is available and set up solver
if SolverFactory('ipopt').available():
    solver = SolverFactory('ipopt')
    solver.options = {'tol': 1e-6,
                      'bound_push': 1e-8}
else:
    solver = None


def assert_categorization(model):
    init_input_set = ComponentSet([model.mixer.S_inlet.flow_vol[0],
                                   model.mixer.E_inlet.flow_vol[0]])

    init_deriv_list = [model.cstr.control_volume.energy_accumulation[0, 'aq']]
    init_diff_list = [model.cstr.control_volume.energy_holdup[0, 'aq']]
    init_fixed_list = [#model.cstr.control_volume.volume[0],
                       model.mixer.E_inlet.temperature[0],
                       model.mixer.S_inlet.temperature[0]]

    init_ic_list = [
            model.cstr.outlet.temperature[0],
            model.cstr.control_volume.material_holdup[0, 'aq', 'S'],
            model.cstr.control_volume.material_holdup[0, 'aq', 'E'],
            model.cstr.control_volume.material_accumulation[0, 'aq', 'C'],
            model.cstr.outlet.conc_mol[0, 'P'],
            model.cstr.control_volume.volume[0],
            ]

    init_alg_list = [
        model.cstr.control_volume.volume[0],
        model.cstr.outlet.flow_vol[0],
        model.cstr.outlet.temperature[0],
        model.cstr.inlet.flow_vol[0],
        model.cstr.inlet.temperature[0],
        model.mixer.outlet.flow_vol[0],
        model.mixer.outlet.temperature[0],
        ]

    for j in model.properties.component_list:
        init_deriv_list.append(
                model.cstr.control_volume.material_accumulation[0, 'aq', j])
        init_diff_list.append(
                model.cstr.control_volume.material_holdup[0, 'aq', j])
        
        init_fixed_list.append(model.mixer.E_inlet.conc_mol[0, j])
        init_fixed_list.append(model.mixer.S_inlet.conc_mol[0, j])

        init_alg_list.extend([
            model.cstr.control_volume.properties_out[0].flow_mol_comp[j],
            model.cstr.inlet.conc_mol[0, j],
            model.cstr.control_volume.properties_in[0].flow_mol_comp[j],
            model.cstr.control_volume.rate_reaction_generation[0, 'aq', j],
            model.mixer.mixed_state[0].flow_mol_comp[j],
            model.mixer.E_inlet_state[0].flow_mol_comp[j],
            model.mixer.S_inlet_state[0].flow_mol_comp[j]
            ])
        if j != 'Solvent':
            init_alg_list.append(model.cstr.outlet.conc_mol[0, j])
            init_alg_list.append(model.mixer.outlet.conc_mol[0, j])
        else:
            init_fixed_list.append(model.cstr.outlet.conc_mol[0, j])
            init_fixed_list.append(model.mixer.outlet.conc_mol[0, j])

    for r in model.reactions.rate_reaction_idx:
        init_alg_list.extend([
            model.cstr.control_volume.reactions[0].reaction_coef[r],
            model.cstr.control_volume.reactions[0].reaction_rate[r],
            model.cstr.control_volume.rate_reaction_extent[0, r]
            ])

    init_deriv_set = ComponentSet(init_deriv_list)
    init_diff_set = ComponentSet(init_diff_list)
    init_fixed_set = ComponentSet(init_fixed_list)
    init_ic_set = ComponentSet(init_ic_list)
    init_alg_set = ComponentSet(init_alg_list)

    assert model._NMPC_NAMESPACE.input_vars.n_vars == len(init_input_set)
    for v in model._NMPC_NAMESPACE.input_vars:
        assert v[0] in init_input_set

    assert model._NMPC_NAMESPACE.deriv_vars.n_vars == len(init_deriv_set)
    for v in model._NMPC_NAMESPACE.deriv_vars:
        assert v[0] in init_deriv_set

    assert len(model._NMPC_NAMESPACE.diff_vars) == len(init_deriv_set)
    for v in model._NMPC_NAMESPACE.diff_vars:
        assert v[0] in init_diff_set

    assert len(model._NMPC_NAMESPACE.fixed_vars) == len(init_fixed_set)
    for v in model._NMPC_NAMESPACE.fixed_vars:
        assert v[0] in init_fixed_set

    assert len(model._NMPC_NAMESPACE.alg_vars) == len(init_alg_set)
    for v in model._NMPC_NAMESPACE.alg_vars:
        assert v[0] in init_alg_set

    assert len(model._NMPC_NAMESPACE.ic_vars) == len(init_ic_set)
    for v in model._NMPC_NAMESPACE.ic_vars:
        assert v[0] in init_ic_set

    assert len(model._NMPC_NAMESPACE.scalar_vars) == 0

    for var in model._NMPC_NAMESPACE.deriv_vars:
        assert len(var) == len(model._NMPC_NAMESPACE.get_time())
        assert var.index_set() is model._NMPC_NAMESPACE.get_time()
    for var in model._NMPC_NAMESPACE.alg_vars:
        assert len(var) == len(model._NMPC_NAMESPACE.get_time())
        assert var.index_set() is model._NMPC_NAMESPACE.get_time()


@pytest.mark.component
def test_constructor_2():
    m_plant = make_model(horizon=6, ntfe=60, ntcp=2)
    m_controller = make_model(horizon=3, ntfe=30, ntcp=2)
    sample_time = 0.5
    # Six samples per horizon, five elements per sample

    initial_plant_inputs = [m_plant.fs.mixer.S_inlet.flow_vol[0],
                            m_plant.fs.mixer.E_inlet.flow_vol[0]]

    # Change some initial conditions - in controller model only

    # Specify temperature instead of energy holdup
    m_controller.fs.cstr.control_volume.energy_holdup\
            [m_controller.fs.time.first(), 'aq'].unfix()
    m_controller.fs.cstr.outlet.temperature[0].fix(300)

    # Specify C_P instead of holdup
    m_controller.fs.cstr.control_volume.material_holdup\
            [m_controller.fs.time.first(), 'aq', 'P'].unfix()
    m_controller.fs.cstr.outlet.conc_mol[m_controller.fs.time.first(), 'P'].fix(0)

    # specify \dot M_C insteady of M_C 
    m_controller.fs.cstr.control_volume.material_holdup\
            [m_controller.fs.time.first(), 'aq', 'C'].unfix()
    m_controller.fs.cstr.control_volume.material_accumulation\
            [m_controller.fs.time.first(), 'aq', 'C'].fix(0)

    nmpc = NMPCSim(m_plant.fs, m_plant.fs.time, m_controller.fs, 
            m_controller.fs.time,
            inputs_at_t0=initial_plant_inputs,
            solver=solver, outlvl=idaeslog.DEBUG, 
            sample_time=sample_time)
    # IPOPT output looks a little weird solving for initial conditions here...
    # has non-zero dual infeasibility, iteration 1 has a non-zero
    # regularization coefficient. (Would love to debug this with a
    # transparent NLP solver...)

    # Check that variables have been categorized properly
    assert_categorization(m_controller.fs)


if __name__ == '__main__':
    test_constructor_2()
