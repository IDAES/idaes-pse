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
Author: Andrew Lee
"""


import pytest

from idaes.core import FlowsheetBlock
from idaes.generic_models.properties.core.eos.ceos import \
    cubic_roots_available
from idaes.generic_models.properties.core.examples.BT_PR import configuration
from idaes.generic_models.properties.core.generic.generic_property import (
        GenericParameterBlock)

from pyomo.environ import (ConcreteModel,
                           Objective,
                           SolverFactory,
                           TerminationCondition,
                           value)

from idaes.generic_models.properties.tests.test_harness import PropertyTestHarness
from idaes.core.util.testing import get_default_solver

import idaes.logger as idaeslog
SOUT = idaeslog.INFO

# Set module level pyest marker
pytestmark = pytest.mark.cubic_root
prop_available = cubic_roots_available()


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_default_solver()


# -----------------------------------------------------------------------------
# # Run test harness
# class TestBasicLV(PropertyTestHarness):
#     def configure(self):
#         self.prop_pack = BT_PR.BTParameterBlock
#         self.param_args = {"valid_phase": ("Liq", "Vap")}
#         self.prop_args = {}
#         self.has_density_terms = True


# class TestBasicL(PropertyTestHarness):
#     def configure(self):
#         self.prop_pack = BT_PR.BTParameterBlock
#         self.param_args = {"valid_phase": "Liq"}
#         self.prop_args = {}
#         self.has_density_terms = True


# class TestBasicV(PropertyTestHarness):
#     def configure(self):
#         self.prop_pack = BT_PR.BTParameterBlock
#         self.param_args = {"valid_phase": "Vap"}
#         self.prop_args = {}
#         self.has_density_terms = True


# -----------------------------------------------------------------------------
# Test robustness and some outputs
@pytest.mark.solver
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.skipif(not prop_available,
                    reason="Cubic root finder not available")
class TestBTExample(object):
    # def test_T_sweep(self):
    #     m = ConcreteModel()

    #     m.fs = FlowsheetBlock(default={'dynamic': False})

    #     m.fs.props = GenericParameterBlock(default=configuration)

    #     m.fs.state = m.fs.props.build_state_block(
    #             [1],
    #             default={"defined_state": True})

    #     m.fs.obj = Objective(expr=(m.fs.state[1].temperature - 510)**2)

    #     for logP in range(8, 13, 1):
    #         m.fs.obj.deactivate()

    #         m.fs.state[1].flow_mol.fix(100)
    #         m.fs.state[1].mole_frac_comp["benzene"].fix(0.5)
    #         m.fs.state[1].mole_frac_comp["toluene"].fix(0.5)
    #         m.fs.state[1].temperature.fix(300)
    #         m.fs.state[1].pressure.fix(10**(0.5*logP))

    #         m.fs.state.initialize(outlvl=SOUT)

    #         m.fs.state[1].temperature.unfix()
    #         m.fs.obj.activate()

    #         solver = SolverFactory('ipopt')
    #         results = solver.solve(m, tee=True)

    #         assert results.solver.termination_condition == \
    #             TerminationCondition.optimal
    #         assert m.fs.state[1].flow_mol_phase["Liq"].value <= 1e-5

    # def test_P_sweep(self):
    #     m = ConcreteModel()

    #     m.fs = FlowsheetBlock(default={'dynamic': False})

    #     m.fs.props = BT_PR.BTParameterBlock(
    #             default={'valid_phase': ('Vap', 'Liq')})

    #     m.fs.state = m.fs.props.build_state_block(
    #             default={'defined_state': True})

    #     for T in range(370, 500, 25):
    #         m.fs.state[1].flow_mol.fix(100)
    #         m.fs.state[1].mole_frac_comp["benzene"].fix(0.5)
    #         m.fs.state[1].mole_frac_comp["toluene"].fix(0.5)
    #         m.fs.state[1].temperature.fix(T)
    #         m.fs.state[1].pressure.fix(1e5)

    #         m.fs.state[1].initialize(outlvl=SOUT)

    #         solver = SolverFactory('ipopt')
    #         results = solver.solve(m)

    #         assert results.solver.termination_condition == \
    #             TerminationCondition.optimal

    #         while m.fs.state[1].pressure.value <= 1e6:
    #             m.fs.state[1].pressure.value = m.fs.state[1].pressure.value + 1e5
    #             solver = SolverFactory('ipopt')
    #             results = solver.solve(m)
    #             assert results.solver.termination_condition == \
    #                 TerminationCondition.optimal
    #             print(T, m.fs.state[1].pressure.value)

    # def test_T350_P1_x5(self):
    #     m = ConcreteModel()

    #     m.fs = FlowsheetBlock(default={'dynamic': False})

    #     m.fs.props = GenericParameterBlock(default=configuration)

    #     m.fs.state = m.fs.props.build_state_block(
    #             [1],
    #             default={"defined_state": True})

    #     m.fs.state[1].flow_mol.fix(100)
    #     m.fs.state[1].mole_frac_comp["benzene"].fix(0.5)
    #     m.fs.state[1].mole_frac_comp["toluene"].fix(0.5)
    #     m.fs.state[1].temperature.fix(350)
    #     m.fs.state[1].pressure.fix(1e5)

    #     # Trigger build of enthalpy and entropy
    #     m.fs.state[1].enth_mol_phase
    #     m.fs.state[1].entr_mol_phase

    #     m.fs.state.initialize(outlvl=SOUT)

    #     solver = SolverFactory('ipopt')
    #     solver.solve(m)

    #     assert pytest.approx(value(m.fs.state[1]._teq), abs=1e-1) == 365
    #     assert 0.0035346 == pytest.approx(
    #             value(m.fs.state[1].compress_fact_phase["Liq"]), 1e-5)
    #     assert 0.966749 == pytest.approx(
    #             value(m.fs.state[1].compress_fact_phase["Vap"]), 1e-5)
    #     assert pytest.approx(
    #             value(m.fs.state[1].fug_coeff_phase_comp["Liq", "benzene"]),
    #             1e-5) == 0.894676
    #     assert pytest.approx(
    #             value(m.fs.state[1].fug_coeff_phase_comp["Liq", "toluene"]),
    #             1e-5) == 0.347566
    #     assert pytest.approx(
    #             value(m.fs.state[1].fug_coeff_phase_comp["Vap", "benzene"]),
    #             1e-5) == 0.971072
    #     assert pytest.approx(
    #             value(m.fs.state[1].fug_coeff_phase_comp["Vap", "toluene"]),
    #             1e-5) == 0.959791

    #     assert pytest.approx(
    #             value(m.fs.state[1].mole_frac_phase_comp["Liq", "benzene"]),
    #             1e-5) == 0.5
    #     assert pytest.approx(
    #             value(m.fs.state[1].mole_frac_phase_comp["Liq", "toluene"]),
    #             1e-5) == 0.5
    #     assert pytest.approx(
    #             value(m.fs.state[1].mole_frac_phase_comp["Vap", "benzene"]),
    #             1e-5) == 0.70584
    #     assert pytest.approx(
    #             value(m.fs.state[1].mole_frac_phase_comp["Vap", "toluene"]),
    #             1e-5) == 0.29416

    #     assert pytest.approx(
    #             value(m.fs.state[1].enth_mol_phase["Liq"]), 1e-5) == 38942.8
    #     assert pytest.approx(
    #             value(m.fs.state[1].enth_mol_phase["Vap"]), 1e-5) == 78048.7
    #     assert pytest.approx(
    #             value(m.fs.state[1].entr_mol_phase["Liq"]), 1e-5) == -367.558
    #     assert pytest.approx(
    #             value(m.fs.state[1].entr_mol_phase["Vap"]), 1e-5) == -269.0553

    # def test_T350_P5_x5(self):
    #     m = ConcreteModel()

    #     m.fs = FlowsheetBlock(default={'dynamic': False})

    #     m.fs.props = BT_PR.BTParameterBlock(
    #             default={'valid_phase': ('Vap', 'Liq')})

    #     m.fs.state = m.fs.props.build_state_block(
    #             default={'defined_state': True})

    #     m.fs.state[1].flow_mol.fix(100)
    #     m.fs.state[1].mole_frac_comp["benzene"].fix(0.5)
    #     m.fs.state[1].mole_frac_comp["toluene"].fix(0.5)
    #     m.fs.state[1].temperature.fix(350)
    #     m.fs.state[1].pressure.fix(5e5)

    #     # Trigger build of enthalpy and entropy
    #     m.fs.state[1].enth_mol_phase
    #     m.fs.state[1].entr_mol_phase

    #     m.fs.state.initialize(outlvl=SOUT)

    #     solver = SolverFactory('ipopt')
    #     solver.solve(m)

    #     assert pytest.approx(value(m.fs.state[1]._teq), 1e-5) == 431.47
    #     assert pytest.approx(
    #             value(m.fs.state[1].compress_fact_phase["Liq"]), 1e-5) == 0.01766
    #     assert pytest.approx(
    #             value(m.fs.state[1].compress_fact_phase["Vap"]), 1e-5) == 0.80245
    #     assert pytest.approx(
    #             value(m.fs.state[1].fug_coeff_phase_comp["Liq", "benzene"]),
    #             1e-5) == 0.181229
    #     assert pytest.approx(
    #             value(m.fs.state[1].fug_coeff_phase_comp["Liq", "toluene"]),
    #             1e-5) == 0.070601
    #     assert pytest.approx(
    #             value(m.fs.state[1].fug_coeff_phase_comp["Vap", "benzene"]),
    #             1e-5) == 0.856523
    #     assert pytest.approx(
    #             value(m.fs.state[1].fug_coeff_phase_comp["Vap", "toluene"]),
    #             1e-5) == 0.799237

    #     assert pytest.approx(
    #             value(m.fs.state[1].mole_frac_phase_comp["Liq", "benzene"]),
    #             1e-5) == 0.5
    #     assert pytest.approx(
    #             value(m.fs.state[1].mole_frac_phase_comp["Liq", "toluene"]),
    #             1e-5) == 0.5
    #     assert pytest.approx(
    #             value(m.fs.state[1].mole_frac_phase_comp["Vap", "benzene"]),
    #             1e-5) == 0.65415
    #     assert pytest.approx(
    #             value(m.fs.state[1].mole_frac_phase_comp["Vap", "toluene"]),
    #             1e-5) == 0.34585

    #     assert pytest.approx(
    #             value(m.fs.state[1].enth_mol_phase["Liq"]), 1e-5) == 38966.9
    #     assert pytest.approx(
    #             value(m.fs.state[1].enth_mol_phase["Vap"]), 1e-5) == 75150.7
    #     assert pytest.approx(
    #             value(m.fs.state[1].entr_mol_phase["Liq"]), 1e-5) == -367.6064
    #     assert pytest.approx(
    #             value(m.fs.state[1].entr_mol_phase["Vap"]), 1e-5) == -287.3318

    # def test_T450_P1_x5(self):
    #     m = ConcreteModel()

    #     m.fs = FlowsheetBlock(default={'dynamic': False})

    #     m.fs.props = BT_PR.BTParameterBlock(
    #             default={'valid_phase': ('Vap', 'Liq')})

    #     m.fs.state = m.fs.props.build_state_block(
    #             default={'defined_state': True})

    #     m.fs.state[1].flow_mol.fix(100)
    #     m.fs.state[1].mole_frac_comp["benzene"].fix(0.5)
    #     m.fs.state[1].mole_frac_comp["toluene"].fix(0.5)
    #     m.fs.state[1].temperature.fix(450)
    #     m.fs.state[1].pressure.fix(1e5)

    #     # Trigger build of enthalpy and entropy
    #     m.fs.state[1].enth_mol_phase
    #     m.fs.state[1].entr_mol_phase

    #     m.fs.state.initialize(outlvl=SOUT)

    #     solver = SolverFactory('ipopt')
    #     solver.solve(m)

    #     assert pytest.approx(value(m.fs.state[1]._teq), 1e-5) == 371.4
    #     assert 0.0033583 == pytest.approx(
    #             value(m.fs.state[1].compress_fact_phase["Liq"]), 1e-5)
    #     assert 0.9821368 == pytest.approx(
    #             value(m.fs.state[1].compress_fact_phase["Vap"]), 1e-5)
    #     assert pytest.approx(
    #             value(m.fs.state[1].fug_coeff_phase_comp["Liq", "benzene"]),
    #             1e-5) == 8.069323
    #     assert pytest.approx(
    #             value(m.fs.state[1].fug_coeff_phase_comp["Liq", "toluene"]),
    #             1e-5) == 4.304955
    #     assert pytest.approx(
    #             value(m.fs.state[1].fug_coeff_phase_comp["Vap", "benzene"]),
    #             1e-5) == 0.985365
    #     assert pytest.approx(
    #             value(m.fs.state[1].fug_coeff_phase_comp["Vap", "toluene"]),
    #             1e-5) == 0.979457

    #     assert pytest.approx(
    #             value(m.fs.state[1].mole_frac_phase_comp["Liq", "benzene"]),
    #             1e-5) == 0.29861
    #     assert pytest.approx(
    #             value(m.fs.state[1].mole_frac_phase_comp["Liq", "toluene"]),
    #             1e-5) == 0.70139
    #     assert pytest.approx(
    #             value(m.fs.state[1].mole_frac_phase_comp["Vap", "benzene"]),
    #             1e-5) == 0.5
    #     assert pytest.approx(
    #             value(m.fs.state[1].mole_frac_phase_comp["Vap", "toluene"]),
    #             1e-5) == 0.5

    #     assert pytest.approx(
    #             value(m.fs.state[1].enth_mol_phase["Liq"]), 1e-5) == 49441.2
    #     assert pytest.approx(
    #             value(m.fs.state[1].enth_mol_phase["Vap"]), 1e-5) == 84175.1
    #     assert pytest.approx(
    #             value(m.fs.state[1].entr_mol_phase["Liq"]), 1e-5) == -333.836
    #     assert pytest.approx(
    #             value(m.fs.state[1].entr_mol_phase["Vap"]), 1e-5) == -247.385

    # def test_T450_P5_x5(self):
    #     m = ConcreteModel()

    #     m.fs = FlowsheetBlock(default={'dynamic': False})

    #     m.fs.props = BT_PR.BTParameterBlock(
    #             default={'valid_phase': ('Vap', 'Liq')})

    #     m.fs.state = m.fs.props.build_state_block(
    #             default={'defined_state': True})

    #     m.fs.state[1].flow_mol.fix(100)
    #     m.fs.state[1].mole_frac_comp["benzene"].fix(0.5)
    #     m.fs.state[1].mole_frac_comp["toluene"].fix(0.5)
    #     m.fs.state[1].temperature.fix(450)
    #     m.fs.state[1].pressure.fix(5e5)

    #     # Trigger build of enthalpy and entropy
    #     m.fs.state[1].enth_mol_phase
    #     m.fs.state[1].entr_mol_phase

    #     m.fs.state.initialize(outlvl=SOUT)

    #     solver = SolverFactory('ipopt')
    #     solver.solve(m)

    #     assert pytest.approx(value(m.fs.state[1]._teq), 1e-5) == 436.93
    #     assert 0.0166181 == pytest.approx(
    #             value(m.fs.state[1].compress_fact_phase["Liq"]), 1e-5)
    #     assert 0.9053766 == pytest.approx(
    #             value(m.fs.state[1].compress_fact_phase["Vap"]), 1e-5)
    #     assert pytest.approx(
    #             value(m.fs.state[1].fug_coeff_phase_comp["Liq", "benzene"]),
    #             1e-5) == 1.63308
    #     assert pytest.approx(
    #             value(m.fs.state[1].fug_coeff_phase_comp["Liq", "toluene"]),
    #             1e-5) == 0.873213
    #     assert pytest.approx(
    #             value(m.fs.state[1].fug_coeff_phase_comp["Vap", "benzene"]),
    #             1e-5) == 0.927534
    #     assert pytest.approx(
    #             value(m.fs.state[1].fug_coeff_phase_comp["Vap", "toluene"]),
    #             1e-5) == 0.898324

    #     assert pytest.approx(
    #             value(m.fs.state[1].mole_frac_phase_comp["Liq", "benzene"]),
    #             1e-5) == 0.3488737
    #     assert pytest.approx(
    #             value(m.fs.state[1].mole_frac_phase_comp["Liq", "toluene"]),
    #             1e-5) == 0.6511263
    #     assert pytest.approx(
    #             value(m.fs.state[1].mole_frac_phase_comp["Vap", "benzene"]),
    #             1e-5) == 0.5
    #     assert pytest.approx(
    #             value(m.fs.state[1].mole_frac_phase_comp["Vap", "toluene"]),
    #             1e-5) == 0.5

    #     assert pytest.approx(
    #             value(m.fs.state[1].enth_mol_phase["Liq"]), 1e-5) == 51095.2
    #     assert pytest.approx(
    #             value(m.fs.state[1].enth_mol_phase["Vap"]), 1e-5) == 83362.3
    #     assert pytest.approx(
    #             value(m.fs.state[1].entr_mol_phase["Liq"]), 1e-5) == -331.676
    #     assert pytest.approx(
    #             value(m.fs.state[1].entr_mol_phase["Vap"]), 1e-5) == -261.961

    def test_T368_P1_x5(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(default={'dynamic': False})

        m.fs.props = GenericParameterBlock(default=configuration)

        m.fs.state = m.fs.props.build_state_block(
                [1],
                default={"defined_state": True})

        m.fs.state[1].flow_mol.fix(100)
        m.fs.state[1].mole_frac_comp["benzene"].fix(0.5)
        m.fs.state[1].mole_frac_comp["toluene"].fix(0.5)
        m.fs.state[1].temperature.fix(368)
        m.fs.state[1].pressure.fix(1e5)

        # Trigger build of enthalpy and entropy
        m.fs.state[1].enth_mol_phase
        m.fs.state[1].entr_mol_phase

        m.fs.state.initialize(outlvl=SOUT)

        solver = SolverFactory('ipopt')
        solver.solve(m)

        assert pytest.approx(value(
            m.fs.state[1]._teq[("Vap", "Liq")]), 1e-5) == 368
        assert 0.003504 == pytest.approx(
                value(m.fs.state[1].compress_fact_phase["Liq"]), 1e-5)
        assert 0.97 == pytest.approx(
                value(m.fs.state[1].compress_fact_phase["Vap"]), 1e-5)
        assert pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Liq", "benzene"]),
                1e-5) == 1.492049
        assert pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Liq", "toluene"]),
                1e-5) == 0.621563
        assert pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Vap", "benzene"]),
                1e-5) == 0.97469
        assert pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Vap", "toluene"]),
                1e-5) == 0.964642

        assert pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Liq", "benzene"]),
                1e-5) == 0.4012128
        assert pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Liq", "toluene"]),
                1e-5) == 0.5987872
        assert pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Vap", "benzene"]),
                1e-5) == 0.6141738
        assert pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Vap", "toluene"]),
                1e-5) == 0.3858262

        m.fs.state[1].enth_mol_phase_comp.display()
        m.fs.state[1].entr_mol_phase_comp.display()

        assert pytest.approx(
                value(m.fs.state[1].enth_mol_phase["Liq"]), 1e-5) == 38235.1
        assert pytest.approx(
                value(m.fs.state[1].enth_mol_phase["Vap"]), 1e-5) == 77155.4
        assert pytest.approx(
                value(m.fs.state[1].entr_mol_phase["Liq"]), 1e-5) == -364.856
        assert pytest.approx(
                value(m.fs.state[1].entr_mol_phase["Vap"]), 1e-5) == -267.892

    # def test_T376_P1_x2(self):
    #     m = ConcreteModel()

    #     m.fs = FlowsheetBlock(default={'dynamic': False})

    #     m.fs.props = BT_PR.BTParameterBlock(
    #             default={'valid_phase': ('Vap', 'Liq')})

    #     m.fs.state = m.fs.props.build_state_block(
    #             default={'defined_state': True})

    #     m.fs.state[1].flow_mol.fix(100)
    #     m.fs.state[1].mole_frac_comp["benzene"].fix(0.2)
    #     m.fs.state[1].mole_frac_comp["toluene"].fix(0.8)
    #     m.fs.state[1].temperature.fix(376)
    #     m.fs.state[1].pressure.fix(1e5)

    #     # Trigger build of enthalpy and entropy
    #     m.fs.state[1].enth_mol_phase
    #     m.fs.state[1].entr_mol_phase

    #     m.fs.state.initialize(outlvl=SOUT)

    #     solver = SolverFactory('ipopt')
    #     solver.solve(m)

    #     assert pytest.approx(value(m.fs.state[1]._teq), 1e-5) == 376
    #     assert 0.00361333 == pytest.approx(
    #             value(m.fs.state[1].compress_fact_phase["Liq"]), 1e-5)
    #     assert 0.968749 == pytest.approx(
    #             value(m.fs.state[1].compress_fact_phase["Vap"]), 1e-5)
    #     assert pytest.approx(
    #             value(m.fs.state[1].fug_coeff_phase_comp["Liq", "benzene"]),
    #             1e-5) == 1.8394188
    #     assert pytest.approx(
    #             value(m.fs.state[1].fug_coeff_phase_comp["Liq", "toluene"]),
    #             1e-5) == 0.7871415
    #     assert pytest.approx(
    #             value(m.fs.state[1].fug_coeff_phase_comp["Vap", "benzene"]),
    #             1e-5) == 0.9763608
    #     assert pytest.approx(
    #             value(m.fs.state[1].fug_coeff_phase_comp["Vap", "toluene"]),
    #             1e-5) == 0.9663611

    #     assert pytest.approx(
    #             value(m.fs.state[1].mole_frac_phase_comp["Liq", "benzene"]),
    #             1e-5) == 0.17342
    #     assert pytest.approx(
    #             value(m.fs.state[1].mole_frac_phase_comp["Liq", "toluene"]),
    #             1e-5) == 0.82658
    #     assert pytest.approx(
    #             value(m.fs.state[1].mole_frac_phase_comp["Vap", "benzene"]),
    #             1e-5) == 0.3267155
    #     assert pytest.approx(
    #             value(m.fs.state[1].mole_frac_phase_comp["Vap", "toluene"]),
    #             1e-5) == 0.6732845

    #     assert pytest.approx(
    #             value(m.fs.state[1].enth_mol_phase["Liq"]), 1e-5) == 31535.8
    #     assert pytest.approx(
    #             value(m.fs.state[1].enth_mol_phase["Vap"]), 1e-5) == 69175.3
    #     assert pytest.approx(
    #             value(m.fs.state[1].entr_mol_phase["Liq"]), 1e-5) == -372.869
    #     assert pytest.approx(
    #             value(m.fs.state[1].entr_mol_phase["Vap"]), 1e-5) == -278.766



















# import pytest
# from pyomo.environ import (ConcreteModel,
#                            Constraint,
#                            Expression,
#                            ExternalFunction,
#                            Param,
#                            Set,
#                            SolverStatus,
#                            sqrt,
#                            TerminationCondition,
#                            value,
#                            Var)

# from idaes.core import (MaterialBalanceType,
#                         EnergyBalanceType,
#                         MaterialFlowBasis,
#                         Component)
# from idaes.core.util.model_statistics import (degrees_of_freedom,
#                                               fixed_variables_set,
#                                               activated_constraints_set)
# from idaes.core.util.testing import get_default_solver
# from idaes.core.util.constants import Constants as const

# from idaes.generic_models.properties.core.generic.generic_property import (
#         GenericParameterBlock)

# from idaes.generic_models.properties.core.state_definitions import FTPx
# from idaes.generic_models.properties.core.phase_equil import smooth_VLE

# from idaes.generic_models.properties.core.examples.BT_PR \
#     import configuration


# # -----------------------------------------------------------------------------
# # Get default solver for testing
# solver = get_default_solver()


# class TestParamBlock(object):
#     def test_build(self):
#         model = ConcreteModel()
#         model.params = GenericParameterBlock(default=configuration)

#         assert isinstance(model.params.phase_list, Set)
#         assert len(model.params.phase_list) == 2
#         for i in model.params.phase_list:
#             assert i in ["Liq", "Vap"]
#         assert model.params.Liq.is_liquid_phase()
#         assert model.params.Vap.is_vapor_phase()

#         assert isinstance(model.params.component_list, Set)
#         assert len(model.params.component_list) == 2
#         for i in model.params.component_list:
#             assert i in ['benzene',
#                          'toluene']
#             assert isinstance(model.params.get_component(i), Component)

#         assert isinstance(model.params._phase_component_set, Set)
#         assert len(model.params._phase_component_set) == 4
#         for i in model.params._phase_component_set:
#             assert i in [("Liq", "benzene"), ("Liq", "toluene"),
#                          ("Vap", "benzene"), ("Vap", "toluene")]

#         assert model.params.config.state_definition == FTPx

#         assert model.params.config.state_bounds == {
#                 "flow_mol": (0, 1000),
#                 "temperature": (273.15, 450),
#                 "pressure": (5e4, 1e6)}

#         assert model.params.config.phase_equilibrium_state == {
#             ("Vap", "Liq"): smooth_VLE}

#         assert isinstance(model.params.phase_equilibrium_idx, Set)
#         assert len(model.params.phase_equilibrium_idx) == 2
#         for i in model.params.phase_equilibrium_idx:
#             assert i in ["PE1", "PE2"]

#         assert model.params.phase_equilibrium_list == {
#             "PE1": {"benzene": ("Vap", "Liq")},
#             "PE2": {"toluene": ("Vap", "Liq")}}

#         assert model.params.pressure_ref.value == 1e5
#         assert model.params.temperature_ref.value == 300

#         assert isinstance(model.params.PR_kappa, Var)
#         for i in model.params.PR_kappa:
#             assert model.params.PR_kappa[i].fixed


# class TestStateBlock(object):
#     @pytest.fixture(scope="class")
#     def model(self):
#         model = ConcreteModel()
#         model.params = GenericParameterBlock(default=configuration)

#         model.props = model.params.build_state_block(
#                 [1],
#                 default={"defined_state": True})

#         return model

#     def test_build(self, model):
#         # Check state variable values and bounds
#         assert isinstance(model.props[1].flow_mol, Var)
#         assert value(model.props[1].flow_mol) == 500
#         assert model.props[1].flow_mol.ub == 1000
#         assert model.props[1].flow_mol.lb == 0

#         assert isinstance(model.props[1].pressure, Var)
#         assert value(model.props[1].pressure) == 5.25e5
#         assert model.props[1].pressure.ub == 1e6
#         assert model.props[1].pressure.lb == 5e4

#         assert isinstance(model.props[1].temperature, Var)
#         assert value(model.props[1].temperature) == 361.575
#         assert model.props[1].temperature.ub == 450
#         assert model.props[1].temperature.lb == 273.15

#         assert isinstance(model.props[1].mole_frac_comp, Var)
#         assert len(model.props[1].mole_frac_comp) == 2
#         for i in model.props[1].mole_frac_comp:
#             assert value(model.props[1].mole_frac_comp[i]) == 0.5

#         # Check supporting variables
#         assert isinstance(model.props[1].flow_mol_phase, Var)
#         assert len(model.props[1].flow_mol_phase) == 2

#         assert isinstance(model.props[1].mole_frac_phase_comp, Var)
#         assert len(model.props[1].mole_frac_phase_comp) == 4

#         assert isinstance(model.props[1].phase_frac, Var)
#         assert len(model.props[1].phase_frac) == 2

#         assert isinstance(model.props[1].total_flow_balance, Constraint)
#         assert len(model.props[1].total_flow_balance) == 1

#         assert isinstance(model.props[1].component_flow_balances, Constraint)
#         assert len(model.props[1].component_flow_balances) == 2

#         assert isinstance(model.props[1].sum_mole_frac, Constraint)
#         assert len(model.props[1].sum_mole_frac) == 1

#         assert not hasattr(model.props[1], "sum_mole_frac_out")

#         assert isinstance(model.props[1].phase_fraction_constraint, Constraint)
#         assert len(model.props[1].phase_fraction_constraint) == 2

#         # Test cubic components
#         assert isinstance(model.props[1].PR_fw, Expression)
#         assert len(model.props[1].PR_fw) == 2
#         for i in model.params.component_list:
#             omega = model.params.get_component(i).omega
#             assert value(model.props[1].PR_fw[i]) == value(
#                 0.37464 + 1.54226*omega - 0.26992*omega**2)

#         assert isinstance(model.props[1].PR_a, Expression)
#         assert len(model.props[1].PR_a) == 2
#         for i in model.params.component_list:
#             Tc = model.params.get_component(i).temperature_crit
#             Pc = model.params.get_component(i).pressure_crit
#             assert pytest.approx(value(model.props[1].PR_a[i]),
#                                  rel=1e-5) == value(
#                 0.45724*((const.gas_constant*Tc)**2/Pc) *
#                 (((1+model.props[1].PR_fw[i] *
#                    (1-sqrt(model.props[1].temperature/Tc)))**2)))

#         assert isinstance(model.props[1].PR_b, Expression)
#         assert len(model.props[1].PR_b) == 2
#         for i in model.params.component_list:
#             Tc = model.params.get_component(i).temperature_crit
#             Pc = model.params.get_component(i).pressure_crit
#             assert value(model.props[1].PR_b[i]) == value(
#                 0.07780*const.gas_constant*Tc/Pc)

#         assert isinstance(model.props[1].PR_am, Expression)
#         assert len(model.props[1].PR_am) == 2
#         for p in model.params.phase_list:
#             assert pytest.approx(value(model.props[1].PR_am[p]),
#                                  rel=1e-5) == value(
#                 sum(sum(model.props[1].mole_frac_phase_comp[p, i] *
#                         model.props[1].mole_frac_phase_comp[p, j] *
#                         sqrt(model.props[1].PR_a[i]*model.props[1].PR_a[j]) *
#                         (1-model.params.PR_kappa[i, j])
#                         for j in model.params.component_list)
#                     for i in model.params.component_list))

#         assert isinstance(model.props[1].PR_bm, Expression)
#         assert len(model.props[1].PR_bm) == 2
#         for p in model.params.phase_list:
#             assert pytest.approx(value(model.props[1].PR_bm[p]),
#                                  rel=1e-5) == value(
#                 sum(model.props[1].mole_frac_phase_comp[p, i] *
#                     model.props[1].PR_b[i]
#                     for i in model.props[1].params.component_list))

#         assert isinstance(model.props[1].PR_A, Expression)
#         assert len(model.props[1].PR_A) == 2
#         for p in model.params.phase_list:
#             assert pytest.approx(value(model.props[1].PR_A[p]),
#                                  rel=1e-5) == value(
#                     model.props[1].PR_am[p]*model.props[1].pressure /
#                     (const.gas_constant*model.props[1].temperature)**2)

#         assert isinstance(model.props[1].PR_B, Expression)
#         assert len(model.props[1].PR_B) == 2
#         for p in model.params.phase_list:
#             assert pytest.approx(value(model.props[1].PR_B[p]),
#                                  rel=1e-5) == value(
#                     model.props[1].PR_bm[p]*model.props[1].pressure /
#                     (const.gas_constant*model.props[1].temperature))

#         assert isinstance(model.props[1].PR_delta, Expression)
#         assert len(model.props[1].PR_delta) == 4
#         for p in model.params.phase_list:
#             for i in model.params.component_list:
#                 assert pytest.approx(value(model.props[1].PR_delta[p, i]),
#                                      rel=1e-5) == value(
#                     2*sqrt(model.props[1].PR_a[i])/model.props[1].PR_am[p] *
#                     sum(model.props[1].mole_frac_phase_comp[p, j] *
#                         sqrt(model.props[1].PR_a[j]) *
#                         (1-model.params.PR_kappa[i, j])
#                         for j in model.params.component_list))

#         assert isinstance(model.props[1].PR_dadT, Expression)
#         assert len(model.props[1].PR_dadT) == 2
#         for p in model.params.phase_list:
#             assert pytest.approx(value(model.props[1].PR_dadT[p]),
#                                  rel=1e-5) == value(
#                 -((const.gas_constant/2)*sqrt(0.45724) *
#                   sum(sum(model.props[1].mole_frac_phase_comp[p, i] *
#                           model.props[1].mole_frac_phase_comp[p, j] *
#                           (1-model.params.PR_kappa[i, j]) *
#                           (model.props[1].PR_fw[j] *
#                            sqrt(model.props[1].PR_a[i] *
#                                 model.params.get_component(j).temperature_crit /
#                                 model.params.get_component(j).pressure_crit) +
#                            model.props[1].PR_fw[i] *
#                            sqrt(model.props[1].PR_a[j] *
#                                 model.params.get_component(i).temperature_crit /
#                                 model.params.get_component(i).pressure_crit))
#                           for j in model.params.component_list)
#                       for i in model.params.component_list) /
#                   sqrt(model.props[1].temperature)))

#         # Test equilibrium state Expressions
#         assert isinstance(model.props[1]._PR_a_eq, Expression)
#         assert len(model.props[1]._PR_a_eq) == 2
#         for i in model.props[1]._PR_a_eq:
#             Tc = model.params.get_component(i[2]).temperature_crit
#             Pc = model.params.get_component(i[2]).pressure_crit
#             assert pytest.approx(value(model.props[1]._PR_a_eq[i]),
#                                  rel=1e-5) == value(
#                 0.45724*((const.gas_constant*Tc)**2/Pc) *
#                 (((1+model.props[1].PR_fw[i[2]] *
#                    (1-sqrt(model.props[1]._teq[i[0], i[1]]/Tc)))**2)))

#         assert isinstance(model.props[1]._PR_am_eq, Expression)
#         assert len(model.props[1]._PR_am_eq) == 2
#         for idx in model.props[1]._PR_am_eq:
#             assert pytest.approx(value(model.props[1]._PR_am_eq[idx]),
#                                  rel=1e-5) == value(
#                 sum(sum(model.props[1].mole_frac_phase_comp[idx[2], i] *
#                         model.props[1].mole_frac_phase_comp[idx[2], j] *
#                         sqrt(model.props[1]._PR_a_eq[idx[0], idx[1], i] *
#                              model.props[1]._PR_a_eq[idx[0], idx[1], j]) *
#                         (1-model.params.PR_kappa[i, j])
#                         for j in model.params.component_list)
#                     for i in model.params.component_list))

#         assert isinstance(model.props[1]._PR_A_eq, Expression)
#         assert len(model.props[1]._PR_A_eq) == 2
#         for i in model.props[1]._PR_A_eq:
#             assert pytest.approx(value(model.props[1]._PR_A_eq[i]),
#                                  rel=1e-5) == value(
#                         model.props[1]._PR_am_eq[i] *
#                         model.props[1].pressure /
#                         (const.gas_constant *
#                          model.props[1]._teq[i[0], i[1]])**2)

#         assert isinstance(model.props[1]._PR_B_eq, Expression)
#         assert len(model.props[1]._PR_B_eq) == 2
#         for i in model.props[1]._PR_B_eq:
#             assert pytest.approx(value(model.props[1]._PR_B_eq[i]),
#                                  rel=1e-5) == value(
#                         model.props[1].PR_bm[i[2]] *
#                         model.props[1].pressure /
#                         (const.gas_constant *
#                          model.props[1]._teq[i[0], i[1]]))

#         assert isinstance(model.props[1]._PR_delta_eq, Expression)
#         assert len(model.props[1]._PR_delta_eq) == 4
#         for idx in model.props[1]._PR_delta_eq:
#             assert pytest.approx(value(model.props[1]._PR_delta_eq[idx]),
#                                  rel=1e-5) == value(
#                 2*sqrt(model.props[1]._PR_a_eq[idx[0], idx[1], idx[3]]) /
#                 model.props[1]._PR_am_eq[idx[0], idx[1], idx[2]] *
#                 sum(model.props[1].mole_frac_phase_comp[idx[2], j] *
#                     sqrt(model.props[1]._PR_a_eq[idx[0], idx[1], j]) *
#                     (1-model.params.PR_kappa[idx[3], j])
#                     for j in model.params.component_list))

#         # Check for external function components
#         assert isinstance(model.props[1]._PR_ext_func_param, Param)
#         assert model.props[1]._PR_ext_func_param.value == 0  # 0 == PR
#         assert isinstance(model.props[1]._PR_proc_Z_liq, ExternalFunction)
#         assert isinstance(model.props[1]._PR_proc_Z_vap, ExternalFunction)

#     def test_get_material_flow_terms(self, model):
#         for p in model.params.phase_list:
#             for j in model.params.component_list:
#                 assert model.props[1].get_material_flow_terms(p, j) == (
#                     model.props[1].flow_mol_phase[p] *
#                     model.props[1].mole_frac_phase_comp[p, j])

#     def test_get_enthalpy_flow_terms(self, model):
#         for p in model.params.phase_list:
#             assert model.props[1].get_enthalpy_flow_terms(p) == (
#                 model.props[1].flow_mol_phase[p] *
#                 model.props[1].enth_mol_phase[p])

#     def test_get_material_density_terms(self, model):
#         for p in model.params.phase_list:
#             for j in model.params.component_list:
#                 assert model.props[1].get_material_density_terms(p, j) == (
#                     model.props[1].dens_mol_phase[p] *
#                     model.props[1].mole_frac_phase_comp[p, j])

#     def test_get_energy_density_terms(self, model):
#         for p in model.params.phase_list:
#             assert model.props[1].get_energy_density_terms(p) == (
#                 model.props[1].dens_mol_phase[p] *
#                 model.props[1].enth_mol_phase[p])

#     def test_default_material_balance_type(self, model):
#         assert model.props[1].default_material_balance_type() == \
#             MaterialBalanceType.componentTotal

#     def test_default_energy_balance_type(self, model):
#         assert model.props[1].default_energy_balance_type() == \
#             EnergyBalanceType.enthalpyTotal

#     def test_get_material_flow_basis(self, model):
#         assert model.props[1].get_material_flow_basis() == \
#             MaterialFlowBasis.molar

#     def test_define_state_vars(self, model):
#         sv = model.props[1].define_state_vars()

#         assert len(sv) == 4
#         for i in sv:
#             assert i in ["flow_mol",
#                          "mole_frac_comp",
#                          "temperature",
#                          "pressure"]

#     def test_define_port_members(self, model):
#         sv = model.props[1].define_state_vars()

#         assert len(sv) == 4
#         for i in sv:
#             assert i in ["flow_mol",
#                          "mole_frac_comp",
#                          "temperature",
#                          "pressure"]

#     def test_define_display_vars(self, model):
#         sv = model.props[1].define_display_vars()

#         assert len(sv) == 4
#         for i in sv:
#             assert i in ["flow_mol",
#                          "mole_frac_comp",
#                          "temperature",
#                          "pressure"]

#     def test_dof(self, model):
#         # Fix state
#         model.props[1].flow_mol.fix(1)
#         model.props[1].temperature.fix(368)
#         model.props[1].pressure.fix(101325)
#         model.props[1].mole_frac_comp["benzene"].fix(0.5)
#         model.props[1].mole_frac_comp["toluene"].fix(0.5)

#         assert degrees_of_freedom(model.props[1]) == 0

#     @pytest.mark.initialize
#     @pytest.mark.solver
#     @pytest.mark.skipif(solver is None, reason="Solver not available")
#     def test_initialize(self, model):
#         orig_fixed_vars = fixed_variables_set(model)
#         orig_act_consts = activated_constraints_set(model)

#         model.props.initialize(optarg={'tol': 1e-6})

#         assert degrees_of_freedom(model) == 0

#         fin_fixed_vars = fixed_variables_set(model)
#         fin_act_consts = activated_constraints_set(model)

#         assert len(fin_act_consts) == len(orig_act_consts)
#         assert len(fin_fixed_vars) == len(orig_fixed_vars)

#         for c in fin_act_consts:
#             assert c in orig_act_consts
#         for v in fin_fixed_vars:
#             assert v in orig_fixed_vars

#     @pytest.mark.solver
#     @pytest.mark.skipif(solver is None, reason="Solver not available")
#     def test_solve(self, model):
#         results = solver.solve(model)

#         # Check for optimal solution
#         assert results.solver.termination_condition == \
#             TerminationCondition.optimal
#         assert results.solver.status == SolverStatus.ok

#     @pytest.mark.initialize
#     @pytest.mark.solver
#     @pytest.mark.skipif(solver is None, reason="Solver not available")
#     def test_solution(self, model):
#         # Check phase equilibrium results
#         assert model.props[1].mole_frac_phase_comp["Liq", "benzene"].value == \
#             pytest.approx(0.4121, abs=1e-4)
#         assert model.props[1].mole_frac_phase_comp["Vap", "benzene"].value == \
#             pytest.approx(0.6339, abs=1e-4)
#         assert model.props[1].phase_frac["Vap"].value == \
#             pytest.approx(0.3961, abs=1e-4)

#     @pytest.mark.ui
#     def test_report(self, model):
#         model.props[1].report()
