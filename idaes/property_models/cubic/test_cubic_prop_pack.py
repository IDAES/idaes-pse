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

from idaes.core import FlowsheetBlock
from idaes.property_models.cubic import BT_PR_VLE
from idaes.property_models.ideal import BTX_ideal_VLE

from pyomo.environ import ConcreteModel, Objective, SolverFactory


def test_P_sweep():
    m = ConcreteModel()

    m.fs = FlowsheetBlock(default={'dynamic': False})

    m.fs.props = BT_PR_VLE.BTParameterBlock(
            default={'valid_phase': ('Vap', 'Liq')})
    m.fs.props_ideal = BTX_ideal_VLE.BTXParameterBlock(
            default={'valid_phase': ('Vap', 'Liq')})

    m.fs.state = m.fs.props.state_block_class(
            default={'parameters': m.fs.props,
                     'defined_state': True})

    m.fs.obj = Objective(expr=(m.fs.state.temperature - 510)**2)

    for logP in range(8, 13, 1):
        m.fs.obj.deactivate()

        m.fs.state.flow_mol.fix(100)
        m.fs.state.mole_frac["benzene"].fix(0.5)
        m.fs.state.mole_frac["toluene"].fix(0.5)
        m.fs.state.temperature.fix(300)
        m.fs.state.pressure.fix(10**(0.5*logP))

        m.fs.state.initialize(outlvl=0)

        m.fs.state.temperature.unfix()
        m.fs.obj.activate()

        solver = SolverFactory('ipopt')
        solver.solve(m, tee=True)

        print(10**(0.5*logP))

        assert m.fs.state.flow_mol_phase["Liq"].value <= 1e-5
