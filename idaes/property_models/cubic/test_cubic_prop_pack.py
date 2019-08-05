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

from pyomo.environ import ConcreteModel, Objective, Param, SolverFactory, value


def test_P_sweep():
    m = ConcreteModel()

    m.fs = FlowsheetBlock(default={'dynamic': False})

    m.fs.props = BT_PR_VLE.BTParameterBlock(
            default={'valid_phase': ('Vap', 'Liq')})
    m.fs.props_ideal = BTX_ideal_VLE.BTXParameterBlock(
            default={'valid_phase': ('Vap', 'Liq')})

    m.fs.state = m.fs.props.state_block_class(
            default={'parameters': m.fs.props})

    P_dict = {5e4: 300, 1e5: 314, 2e5: 340, 5e5: 380, 8e5: 406, 1e6: 418}

    m.fs.T_obj = Param(initialize=300, mutable=True)
    m.fs.obj = Objective(expr=(m.fs.state.temperature - m.fs.T_obj)**2)

    for P, T in P_dict.items():
        m.fs.obj.deactivate()

        m.fs.state.flow_mol.fix(100)
        m.fs.state.mole_frac["benzene"].fix(0.5)
        m.fs.state.temperature.fix(300)
        m.fs.state.pressure.fix(P)

        m.fs.state.mole_frac["toluene"].value = (
                1 - value(m.fs.state.mole_frac["benzene"]))

        m.fs.state.initialize(outlvl=0)

        m.fs.state.temperature.unfix()
        m.fs.obj.activate()
        m.fs.T_obj.value = T + 110

        solver = SolverFactory('ipopt')
        solver.solve(m, tee=True)

        print(P, T)

        assert m.fs.state.flow_mol_phase["Liq"].value <= 1e-5
