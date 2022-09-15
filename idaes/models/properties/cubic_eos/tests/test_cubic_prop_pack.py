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
import pytest

from pyomo.environ import (
    ConcreteModel,
    Constraint,
    ExternalFunction,
    Expression,
    Param,
    Set,
    sqrt,
    Var,
)

from idaes.core import FlowsheetBlock, Component
from idaes.models.properties.cubic_eos.cubic_prop_pack import (
    cubic_roots_available,
    CubicParameterBlock,
    CubicStateBlock,
    CubicEoS,
    EoS_param,
)


# Set module level pyest marker
pytestmark = pytest.mark.cubic_root
prop_available = cubic_roots_available()


@pytest.mark.unit
def test_CubicEoS():
    assert len(CubicEoS) == 2
    assert CubicEoS.PR
    assert CubicEoS.SRK


class TestParameterBlock(object):
    @pytest.mark.unit
    def test_build_default(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.params = CubicParameterBlock()

        assert m.fs.params.state_block_class is CubicStateBlock
        assert m.fs.params.config.valid_phase == ("Vap", "Liq")

        assert isinstance(m.fs.params.phase_list, Set)
        assert len(m.fs.params.phase_list) == 2
        for p in m.fs.params.phase_list:
            assert p in ["Vap", "Liq"]

    @pytest.mark.unit
    def test_build_VL(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.params = CubicParameterBlock(valid_phase=("Vap", "Liq"))

        assert m.fs.params.state_block_class is CubicStateBlock
        assert m.fs.params.config.valid_phase == ("Vap", "Liq")

        assert isinstance(m.fs.params.phase_list, Set)
        assert len(m.fs.params.phase_list) == 2
        for p in m.fs.params.phase_list:
            assert p in ["Vap", "Liq"]

    @pytest.mark.unit
    def test_build_LV(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.params = CubicParameterBlock(valid_phase=("Liq", "Vap"))

        assert m.fs.params.state_block_class is CubicStateBlock
        assert m.fs.params.config.valid_phase == ("Liq", "Vap")

        assert isinstance(m.fs.params.phase_list, Set)
        assert len(m.fs.params.phase_list) == 2
        for p in m.fs.params.phase_list:
            assert p in ["Vap", "Liq"]

    @pytest.mark.unit
    def test_build_L(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.params = CubicParameterBlock(valid_phase="Liq")

        assert m.fs.params.state_block_class is CubicStateBlock
        assert m.fs.params.config.valid_phase == ("Liq")

        assert isinstance(m.fs.params.phase_list, Set)
        assert len(m.fs.params.phase_list) == 1
        for p in m.fs.params.phase_list:
            assert p in ["Liq"]

    @pytest.mark.unit
    def test_build_V(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.params = CubicParameterBlock(valid_phase="Vap")

        assert m.fs.params.state_block_class is CubicStateBlock
        assert m.fs.params.config.valid_phase == ("Vap")

        assert isinstance(m.fs.params.phase_list, Set)
        assert len(m.fs.params.phase_list) == 1
        for p in m.fs.params.phase_list:
            assert p in ["Vap"]


class TestStateBlock_LV_PR(object):
    @pytest.fixture()
    def model(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.params = CubicParameterBlock()

        m.fs.params.a = Component()
        m.fs.params.b = Component()

        m.fs.params.cubic_type = CubicEoS.PR

        m.fs.params.gas_const = Param(default=8.314462618)

        m.fs.params.pressure_crit = Param(
            m.fs.params.component_list, initialize={"a": 5e6, "b": 4e6}
        )
        m.fs.params.temperature_crit = Param(
            m.fs.params.component_list, initialize={"a": 500, "b": 600}
        )

        m.fs.params.omega = Param(
            m.fs.params.component_list, initialize={"a": 0.2, "b": 0.2}
        )

        m.fs.params.kappa = Param(
            m.fs.params.component_list,
            m.fs.params.component_list,
            initialize={
                ("a", "a"): 0.0,
                ("a", "b"): 0.0,
                ("b", "a"): 0.0,
                ("b", "b"): 0.0,
            },
        )

        return m

    @pytest.mark.unit
    def test_build_default(self, model):
        model.fs.props = model.fs.params.build_state_block([1])

        assert isinstance(model.fs.props[1].flow_mol, Var)
        assert len(model.fs.props[1].flow_mol) == 1

        assert isinstance(model.fs.props[1].mole_frac_comp, Var)
        assert len(model.fs.props[1].mole_frac_comp) == 2
        for j in model.fs.props[1].mole_frac_comp:
            assert j in ["a", "b"]

        assert isinstance(model.fs.props[1].pressure, Var)
        assert len(model.fs.props[1].pressure) == 1
        assert isinstance(model.fs.props[1].temperature, Var)
        assert len(model.fs.props[1].temperature) == 1

        assert isinstance(model.fs.props[1].flow_mol_phase, Var)
        assert len(model.fs.props[1].flow_mol_phase) == 2
        for j in model.fs.props[1].flow_mol_phase:
            assert j in ["Vap", "Liq"]

        assert isinstance(model.fs.props[1].mole_frac_phase_comp, Var)
        assert len(model.fs.props[1].mole_frac_phase_comp) == 4
        for j in model.fs.props[1].mole_frac_phase_comp:
            assert j in [("Liq", "a"), ("Liq", "b"), ("Vap", "a"), ("Vap", "b")]

        assert isinstance(model.fs.props[1].total_flow_balance, Constraint)
        assert len(model.fs.props[1].total_flow_balance) == 1
        assert str(model.fs.props[1].total_flow_balance.body) == str(
            model.fs.props[1].flow_mol_phase["Liq"]
            + model.fs.props[1].flow_mol_phase["Vap"]
            - model.fs.props[1].flow_mol
        )

        assert isinstance(model.fs.props[1].component_flow_balances, Constraint)
        assert len(model.fs.props[1].component_flow_balances) == 2
        for j in model.fs.props[1].component_flow_balances:
            assert j in ["a", "b"]
            assert str(model.fs.props[1].component_flow_balances[j].body) == str(
                model.fs.props[1].flow_mol * model.fs.props[1].mole_frac_comp[j]
                - (
                    model.fs.props[1].flow_mol_phase["Liq"]
                    * model.fs.props[1].mole_frac_phase_comp["Liq", j]
                    + model.fs.props[1].flow_mol_phase["Vap"]
                    * model.fs.props[1].mole_frac_phase_comp["Vap", j]
                )
            )

        assert isinstance(model.fs.props[1].sum_mole_frac, Constraint)
        assert len(model.fs.props[1].sum_mole_frac) == 1
        assert str(model.fs.props[1].sum_mole_frac.body) == str(
            sum(
                model.fs.props[1].mole_frac_phase_comp["Liq", i]
                for i in model.fs.params.component_list
            )
            - sum(
                model.fs.props[1].mole_frac_phase_comp["Vap", i]
                for i in model.fs.params.component_list
            )
        )

        assert isinstance(model.fs.props[1].sum_mole_frac_out, Constraint)
        assert len(model.fs.props[1].sum_mole_frac_out) == 1
        assert str(model.fs.props[1].sum_mole_frac_out.body) == str(
            sum(
                model.fs.props[1].mole_frac_comp[i]
                for i in model.fs.params.component_list
            )
        )

        assert isinstance(model.fs.props[1]._teq, Var)
        assert len(model.fs.props[1]._teq) == 1

        assert isinstance(model.fs.props[1]._t1, Var)
        assert len(model.fs.props[1]._t1) == 1

        assert isinstance(model.fs.props[1].eps_1, Param)
        assert model.fs.props[1].eps_1.value == 0.01
        assert isinstance(model.fs.props[1].eps_2, Param)
        assert model.fs.props[1].eps_2.value == 0.0005

        assert isinstance(model.fs.props[1]._t1_constraint, Constraint)
        assert len(model.fs.props[1]._t1_constraint) == 1
        assert str(model.fs.props[1]._t1_constraint.body) == str(
            model.fs.props[1]._t1
            - 0.5
            * (
                model.fs.props[1].temperature
                + model.fs.props[1].temperature_bubble
                + sqrt(
                    (
                        model.fs.props[1].temperature
                        - model.fs.props[1].temperature_bubble
                    )
                    ** 2
                    + model.fs.props[1].eps_1 ** 2
                )
            )
        )

        assert isinstance(model.fs.props[1]._teq_constraint, Constraint)
        assert len(model.fs.props[1]._teq_constraint) == 1
        assert str(model.fs.props[1]._teq_constraint.body) == str(
            model.fs.props[1]._teq
            - 0.5
            * (
                model.fs.props[1]._t1
                + model.fs.props[1].temperature_dew
                - sqrt(
                    (model.fs.props[1]._t1 - model.fs.props[1].temperature_dew) ** 2
                    + model.fs.props[1].eps_2 ** 2
                )
            )
        )

        assert isinstance(model.fs.props[1]._tr_eq, Expression)
        assert len(model.fs.props[1]._tr_eq) == 2
        for j in model.fs.props[1]._tr_eq:
            assert j in ["a", "b"]
        assert str(model.fs.props[1]._tr_eq[j].expr) == str(
            model.fs.props[1]._teq / model.fs.props[1].params.temperature_crit[j]
        )

        assert isinstance(model.fs.props[1].equilibrium_constraint, Constraint)
        assert len(model.fs.props[1].equilibrium_constraint) == 2
        for j in model.fs.props[1].equilibrium_constraint:
            assert j in ["a", "b"]
            assert str(model.fs.props[1].equilibrium_constraint[j].body) == str(
                model.fs.props[1]._log_equilibrium_cubic("Vap", j)
                - model.fs.props[1]._log_equilibrium_cubic("Liq", j)
            )

    @pytest.mark.unit
    def test_common_cubic(self, model):
        model.fs.props = model.fs.params.build_state_block([1])

        eos = CubicEoS.PR

        assert model.fs.props[1].omegaA == EoS_param[eos]["omegaA"]
        assert model.fs.props[1].EoS_Bc == EoS_param[eos]["coeff_b"]
        assert model.fs.props[1].EoS_u == EoS_param[eos]["u"]
        assert model.fs.props[1].EoS_w == EoS_param[eos]["w"]
        assert model.fs.props[1].EoS_p == sqrt(
            EoS_param[eos]["u"] ** 2 - 4 * EoS_param[eos]["w"]
        )

        assert isinstance(model.fs.props[1].fw, Param)
        assert isinstance(model.fs.props[1].b, Param)

        assert isinstance(model.fs.props[1].compress_fact_liq_func, ExternalFunction)
        assert isinstance(model.fs.props[1].compress_fact_vap_func, ExternalFunction)


class TestStateBlock_L_PR(object):
    @pytest.fixture()
    def model(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.params = CubicParameterBlock(valid_phase="Liq")

        m.fs.params.a = Component()
        m.fs.params.b = Component()

        m.fs.params.cubic_type = CubicEoS.PR

        m.fs.params.gas_const = Param(default=8.314462618)

        m.fs.params.pressure_crit = Param(
            m.fs.params.component_list, initialize={"a": 5e6, "b": 4e6}
        )
        m.fs.params.temperature_crit = Param(
            m.fs.params.component_list, initialize={"a": 500, "b": 600}
        )

        m.fs.params.omega = Param(
            m.fs.params.component_list, initialize={"a": 0.2, "b": 0.2}
        )

        m.fs.params.kappa = Param(
            m.fs.params.component_list,
            m.fs.params.component_list,
            initialize={
                ("a", "a"): 0.0,
                ("a", "b"): 0.0,
                ("b", "a"): 0.0,
                ("b", "b"): 0.0,
            },
        )

        return m

    @pytest.mark.unit
    def test_build_default(self, model):
        model.fs.props = model.fs.params.build_state_block([1])

        assert isinstance(model.fs.props[1].flow_mol, Var)
        assert len(model.fs.props[1].flow_mol) == 1

        assert isinstance(model.fs.props[1].mole_frac_comp, Var)
        assert len(model.fs.props[1].mole_frac_comp) == 2
        for j in model.fs.props[1].mole_frac_comp:
            assert j in ["a", "b"]

        assert isinstance(model.fs.props[1].pressure, Var)
        assert len(model.fs.props[1].pressure) == 1
        assert isinstance(model.fs.props[1].temperature, Var)
        assert len(model.fs.props[1].temperature) == 1

        assert isinstance(model.fs.props[1].flow_mol_phase, Var)
        assert len(model.fs.props[1].flow_mol_phase) == 1
        for j in model.fs.props[1].flow_mol_phase:
            assert j in ["Liq"]

        assert isinstance(model.fs.props[1].mole_frac_phase_comp, Var)
        assert len(model.fs.props[1].mole_frac_phase_comp) == 2
        for j in model.fs.props[1].mole_frac_phase_comp:
            assert j in [("Liq", "a"), ("Liq", "b")]

        assert isinstance(model.fs.props[1].total_flow_balance, Constraint)
        assert len(model.fs.props[1].total_flow_balance) == 1
        assert str(model.fs.props[1].total_flow_balance.body) == str(
            model.fs.props[1].flow_mol - model.fs.props[1].flow_mol_phase["Liq"]
        )

        assert isinstance(model.fs.props[1].component_flow_balances, Constraint)
        assert len(model.fs.props[1].component_flow_balances) == 2
        for j in model.fs.props[1].component_flow_balances:
            assert j in ["a", "b"]
            assert str(model.fs.props[1].component_flow_balances[j].body) == str(
                model.fs.props[1].mole_frac_comp[j]
                - model.fs.props[1].mole_frac_phase_comp["Liq", j]
            )

        assert isinstance(model.fs.props[1].sum_mole_frac_out, Constraint)
        assert len(model.fs.props[1].sum_mole_frac_out) == 1
        assert str(model.fs.props[1].sum_mole_frac_out.body) == str(
            sum(
                model.fs.props[1].mole_frac_comp[i]
                for i in model.fs.params.component_list
            )
        )

        assert isinstance(model.fs.props[1]._teq, Expression)
        assert len(model.fs.props[1]._teq) == 1
        assert str(model.fs.props[1]._teq.expr) == str(model.fs.props[1].temperature)

    @pytest.mark.unit
    def test_common_cubic(self, model):
        model.fs.props = model.fs.params.build_state_block([1])

        eos = CubicEoS.PR

        assert model.fs.props[1].omegaA == EoS_param[eos]["omegaA"]
        assert model.fs.props[1].EoS_Bc == EoS_param[eos]["coeff_b"]
        assert model.fs.props[1].EoS_u == EoS_param[eos]["u"]
        assert model.fs.props[1].EoS_w == EoS_param[eos]["w"]
        assert model.fs.props[1].EoS_p == sqrt(
            EoS_param[eos]["u"] ** 2 - 4 * EoS_param[eos]["w"]
        )

        assert isinstance(model.fs.props[1].fw, Param)
        assert isinstance(model.fs.props[1].b, Param)

        assert isinstance(model.fs.props[1].compress_fact_liq_func, ExternalFunction)


class TestStateBlock_V_PR(object):
    @pytest.fixture()
    def model(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.params = CubicParameterBlock(valid_phase="Vap")

        m.fs.params.a = Component()
        m.fs.params.b = Component()

        m.fs.params.cubic_type = CubicEoS.PR

        m.fs.params.gas_const = Param(default=8.314462618)

        m.fs.params.pressure_crit = Param(
            m.fs.params.component_list, initialize={"a": 5e6, "b": 4e6}
        )
        m.fs.params.temperature_crit = Param(
            m.fs.params.component_list, initialize={"a": 500, "b": 600}
        )

        m.fs.params.omega = Param(
            m.fs.params.component_list, initialize={"a": 0.2, "b": 0.2}
        )

        m.fs.params.kappa = Param(
            m.fs.params.component_list,
            m.fs.params.component_list,
            initialize={
                ("a", "a"): 0.0,
                ("a", "b"): 0.0,
                ("b", "a"): 0.0,
                ("b", "b"): 0.0,
            },
        )

        return m

    @pytest.mark.unit
    def test_build_default(self, model):
        model.fs.props = model.fs.params.build_state_block([1])

        assert isinstance(model.fs.props[1].flow_mol, Var)
        assert len(model.fs.props[1].flow_mol) == 1

        assert isinstance(model.fs.props[1].mole_frac_comp, Var)
        assert len(model.fs.props[1].mole_frac_comp) == 2
        for j in model.fs.props[1].mole_frac_comp:
            assert j in ["a", "b"]

        assert isinstance(model.fs.props[1].pressure, Var)
        assert len(model.fs.props[1].pressure) == 1
        assert isinstance(model.fs.props[1].temperature, Var)
        assert len(model.fs.props[1].temperature) == 1

        assert isinstance(model.fs.props[1].flow_mol_phase, Var)
        assert len(model.fs.props[1].flow_mol_phase) == 1
        for j in model.fs.props[1].flow_mol_phase:
            assert j in ["Vap"]

        assert isinstance(model.fs.props[1].mole_frac_phase_comp, Var)
        assert len(model.fs.props[1].mole_frac_phase_comp) == 2
        for j in model.fs.props[1].mole_frac_phase_comp:
            assert j in [("Vap", "a"), ("Vap", "b")]

        assert isinstance(model.fs.props[1].total_flow_balance, Constraint)
        assert len(model.fs.props[1].total_flow_balance) == 1
        assert str(model.fs.props[1].total_flow_balance.body) == str(
            model.fs.props[1].flow_mol - model.fs.props[1].flow_mol_phase["Vap"]
        )

        assert isinstance(model.fs.props[1].component_flow_balances, Constraint)
        assert len(model.fs.props[1].component_flow_balances) == 2
        for j in model.fs.props[1].component_flow_balances:
            assert j in ["a", "b"]
            assert str(model.fs.props[1].component_flow_balances[j].body) == str(
                model.fs.props[1].mole_frac_comp[j]
                - model.fs.props[1].mole_frac_phase_comp["Vap", j]
            )

        assert isinstance(model.fs.props[1].sum_mole_frac_out, Constraint)
        assert len(model.fs.props[1].sum_mole_frac_out) == 1
        assert str(model.fs.props[1].sum_mole_frac_out.body) == str(
            sum(
                model.fs.props[1].mole_frac_comp[i]
                for i in model.fs.params.component_list
            )
        )

        assert isinstance(model.fs.props[1]._teq, Expression)
        assert len(model.fs.props[1]._teq) == 1
        assert str(model.fs.props[1]._teq.expr) == str(model.fs.props[1].temperature)

    @pytest.mark.unit
    def test_common_cubic(self, model):
        model.fs.props = model.fs.params.build_state_block([1])

        eos = CubicEoS.PR

        assert model.fs.props[1].omegaA == EoS_param[eos]["omegaA"]
        assert model.fs.props[1].EoS_Bc == EoS_param[eos]["coeff_b"]
        assert model.fs.props[1].EoS_u == EoS_param[eos]["u"]
        assert model.fs.props[1].EoS_w == EoS_param[eos]["w"]
        assert model.fs.props[1].EoS_p == sqrt(
            EoS_param[eos]["u"] ** 2 - 4 * EoS_param[eos]["w"]
        )

        assert isinstance(model.fs.props[1].fw, Param)
        assert isinstance(model.fs.props[1].b, Param)

        assert isinstance(model.fs.props[1].compress_fact_vap_func, ExternalFunction)


class TestStateBlock_LV_SRK(object):
    @pytest.fixture()
    def model(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.params = CubicParameterBlock()

        m.fs.params.a = Component()
        m.fs.params.b = Component()

        m.fs.params.cubic_type = CubicEoS.SRK

        m.fs.params.gas_const = Param(default=8.314462618)

        m.fs.params.pressure_crit = Param(
            m.fs.params.component_list, initialize={"a": 5e6, "b": 4e6}
        )
        m.fs.params.temperature_crit = Param(
            m.fs.params.component_list, initialize={"a": 500, "b": 600}
        )

        m.fs.params.omega = Param(
            m.fs.params.component_list, initialize={"a": 0.2, "b": 0.2}
        )

        m.fs.params.kappa = Param(
            m.fs.params.component_list,
            m.fs.params.component_list,
            initialize={
                ("a", "a"): 0.0,
                ("a", "b"): 0.0,
                ("b", "a"): 0.0,
                ("b", "b"): 0.0,
            },
        )

        return m

    @pytest.mark.unit
    def test_build_default(self, model):
        model.fs.props = model.fs.params.build_state_block([1])

        assert isinstance(model.fs.props[1].flow_mol, Var)
        assert len(model.fs.props[1].flow_mol) == 1

        assert isinstance(model.fs.props[1].mole_frac_comp, Var)
        assert len(model.fs.props[1].mole_frac_comp) == 2
        for j in model.fs.props[1].mole_frac_comp:
            assert j in ["a", "b"]

        assert isinstance(model.fs.props[1].pressure, Var)
        assert len(model.fs.props[1].pressure) == 1
        assert isinstance(model.fs.props[1].temperature, Var)
        assert len(model.fs.props[1].temperature) == 1

        assert isinstance(model.fs.props[1].flow_mol_phase, Var)
        assert len(model.fs.props[1].flow_mol_phase) == 2
        for j in model.fs.props[1].flow_mol_phase:
            assert j in ["Vap", "Liq"]

        assert isinstance(model.fs.props[1].mole_frac_phase_comp, Var)
        assert len(model.fs.props[1].mole_frac_phase_comp) == 4
        for j in model.fs.props[1].mole_frac_phase_comp:
            assert j in [("Liq", "a"), ("Liq", "b"), ("Vap", "a"), ("Vap", "b")]

        assert isinstance(model.fs.props[1].total_flow_balance, Constraint)
        assert len(model.fs.props[1].total_flow_balance) == 1
        assert str(model.fs.props[1].total_flow_balance.body) == str(
            model.fs.props[1].flow_mol_phase["Liq"]
            + model.fs.props[1].flow_mol_phase["Vap"]
            - model.fs.props[1].flow_mol
        )

        assert isinstance(model.fs.props[1].component_flow_balances, Constraint)
        assert len(model.fs.props[1].component_flow_balances) == 2
        for j in model.fs.props[1].component_flow_balances:
            assert j in ["a", "b"]
            assert str(model.fs.props[1].component_flow_balances[j].body) == str(
                model.fs.props[1].flow_mol * model.fs.props[1].mole_frac_comp[j]
                - (
                    model.fs.props[1].flow_mol_phase["Liq"]
                    * model.fs.props[1].mole_frac_phase_comp["Liq", j]
                    + model.fs.props[1].flow_mol_phase["Vap"]
                    * model.fs.props[1].mole_frac_phase_comp["Vap", j]
                )
            )

        assert isinstance(model.fs.props[1].sum_mole_frac, Constraint)
        assert len(model.fs.props[1].sum_mole_frac) == 1
        assert str(model.fs.props[1].sum_mole_frac.body) == str(
            sum(
                model.fs.props[1].mole_frac_phase_comp["Liq", i]
                for i in model.fs.params.component_list
            )
            - sum(
                model.fs.props[1].mole_frac_phase_comp["Vap", i]
                for i in model.fs.params.component_list
            )
        )

        assert isinstance(model.fs.props[1].sum_mole_frac_out, Constraint)
        assert len(model.fs.props[1].sum_mole_frac_out) == 1
        assert str(model.fs.props[1].sum_mole_frac_out.body) == str(
            sum(
                model.fs.props[1].mole_frac_comp[i]
                for i in model.fs.params.component_list
            )
        )

        assert isinstance(model.fs.props[1]._teq, Var)
        assert len(model.fs.props[1]._teq) == 1

        assert isinstance(model.fs.props[1]._t1, Var)
        assert len(model.fs.props[1]._t1) == 1

        assert isinstance(model.fs.props[1].eps_1, Param)
        assert model.fs.props[1].eps_1.value == 0.01
        assert isinstance(model.fs.props[1].eps_2, Param)
        assert model.fs.props[1].eps_2.value == 0.0005

        assert isinstance(model.fs.props[1]._t1_constraint, Constraint)
        assert len(model.fs.props[1]._t1_constraint) == 1
        assert str(model.fs.props[1]._t1_constraint.body) == str(
            model.fs.props[1]._t1
            - 0.5
            * (
                model.fs.props[1].temperature
                + model.fs.props[1].temperature_bubble
                + sqrt(
                    (
                        model.fs.props[1].temperature
                        - model.fs.props[1].temperature_bubble
                    )
                    ** 2
                    + model.fs.props[1].eps_1 ** 2
                )
            )
        )

        assert isinstance(model.fs.props[1]._teq_constraint, Constraint)
        assert len(model.fs.props[1]._teq_constraint) == 1
        assert str(model.fs.props[1]._teq_constraint.body) == str(
            model.fs.props[1]._teq
            - 0.5
            * (
                model.fs.props[1]._t1
                + model.fs.props[1].temperature_dew
                - sqrt(
                    (model.fs.props[1]._t1 - model.fs.props[1].temperature_dew) ** 2
                    + model.fs.props[1].eps_2 ** 2
                )
            )
        )

        assert isinstance(model.fs.props[1]._tr_eq, Expression)
        assert len(model.fs.props[1]._tr_eq) == 2
        for j in model.fs.props[1]._tr_eq:
            assert j in ["a", "b"]
        assert str(model.fs.props[1]._tr_eq[j].expr) == str(
            model.fs.props[1]._teq / model.fs.props[1].params.temperature_crit[j]
        )

        assert isinstance(model.fs.props[1].equilibrium_constraint, Constraint)
        assert len(model.fs.props[1].equilibrium_constraint) == 2
        for j in model.fs.props[1].equilibrium_constraint:
            assert j in ["a", "b"]
            assert str(model.fs.props[1].equilibrium_constraint[j].body) == str(
                model.fs.props[1]._log_equilibrium_cubic("Vap", j)
                - model.fs.props[1]._log_equilibrium_cubic("Liq", j)
            )

    @pytest.mark.unit
    def test_common_cubic(self, model):
        model.fs.props = model.fs.params.build_state_block([1])

        eos = CubicEoS.SRK

        assert model.fs.props[1].omegaA == EoS_param[eos]["omegaA"]
        assert model.fs.props[1].EoS_Bc == EoS_param[eos]["coeff_b"]
        assert model.fs.props[1].EoS_u == EoS_param[eos]["u"]
        assert model.fs.props[1].EoS_w == EoS_param[eos]["w"]
        assert model.fs.props[1].EoS_p == sqrt(
            EoS_param[eos]["u"] ** 2 - 4 * EoS_param[eos]["w"]
        )

        assert isinstance(model.fs.props[1].fw, Param)
        assert isinstance(model.fs.props[1].b, Param)

        assert isinstance(model.fs.props[1].compress_fact_liq_func, ExternalFunction)
        assert isinstance(model.fs.props[1].compress_fact_vap_func, ExternalFunction)


class TestStateBlock_L_SRK(object):
    @pytest.fixture()
    def model(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.params = CubicParameterBlock(valid_phase="Liq")

        m.fs.params.a = Component()
        m.fs.params.b = Component()

        m.fs.params.cubic_type = CubicEoS.SRK

        m.fs.params.gas_const = Param(default=8.314462618)

        m.fs.params.pressure_crit = Param(
            m.fs.params.component_list, initialize={"a": 5e6, "b": 4e6}
        )
        m.fs.params.temperature_crit = Param(
            m.fs.params.component_list, initialize={"a": 500, "b": 600}
        )

        m.fs.params.omega = Param(
            m.fs.params.component_list, initialize={"a": 0.2, "b": 0.2}
        )

        m.fs.params.kappa = Param(
            m.fs.params.component_list,
            m.fs.params.component_list,
            initialize={
                ("a", "a"): 0.0,
                ("a", "b"): 0.0,
                ("b", "a"): 0.0,
                ("b", "b"): 0.0,
            },
        )

        return m

    @pytest.mark.unit
    def test_build_default(self, model):
        model.fs.props = model.fs.params.build_state_block([1])

        assert isinstance(model.fs.props[1].flow_mol, Var)
        assert len(model.fs.props[1].flow_mol) == 1

        assert isinstance(model.fs.props[1].mole_frac_comp, Var)
        assert len(model.fs.props[1].mole_frac_comp) == 2
        for j in model.fs.props[1].mole_frac_comp:
            assert j in ["a", "b"]

        assert isinstance(model.fs.props[1].pressure, Var)
        assert len(model.fs.props[1].pressure) == 1
        assert isinstance(model.fs.props[1].temperature, Var)
        assert len(model.fs.props[1].temperature) == 1

        assert isinstance(model.fs.props[1].flow_mol_phase, Var)
        assert len(model.fs.props[1].flow_mol_phase) == 1
        for j in model.fs.props[1].flow_mol_phase:
            assert j in ["Liq"]

        assert isinstance(model.fs.props[1].mole_frac_phase_comp, Var)
        assert len(model.fs.props[1].mole_frac_phase_comp) == 2
        for j in model.fs.props[1].mole_frac_phase_comp:
            assert j in [("Liq", "a"), ("Liq", "b")]

        assert isinstance(model.fs.props[1].total_flow_balance, Constraint)
        assert len(model.fs.props[1].total_flow_balance) == 1
        assert str(model.fs.props[1].total_flow_balance.body) == str(
            model.fs.props[1].flow_mol - model.fs.props[1].flow_mol_phase["Liq"]
        )

        assert isinstance(model.fs.props[1].component_flow_balances, Constraint)
        assert len(model.fs.props[1].component_flow_balances) == 2
        for j in model.fs.props[1].component_flow_balances:
            assert j in ["a", "b"]
            assert str(model.fs.props[1].component_flow_balances[j].body) == str(
                model.fs.props[1].mole_frac_comp[j]
                - model.fs.props[1].mole_frac_phase_comp["Liq", j]
            )

        assert isinstance(model.fs.props[1].sum_mole_frac_out, Constraint)
        assert len(model.fs.props[1].sum_mole_frac_out) == 1
        assert str(model.fs.props[1].sum_mole_frac_out.body) == str(
            sum(
                model.fs.props[1].mole_frac_comp[i]
                for i in model.fs.params.component_list
            )
        )

        assert isinstance(model.fs.props[1]._teq, Expression)
        assert len(model.fs.props[1]._teq) == 1
        assert str(model.fs.props[1]._teq.expr) == str(model.fs.props[1].temperature)

    @pytest.mark.unit
    def test_common_cubic(self, model):
        model.fs.props = model.fs.params.build_state_block([1])

        eos = CubicEoS.SRK

        assert model.fs.props[1].omegaA == EoS_param[eos]["omegaA"]
        assert model.fs.props[1].EoS_Bc == EoS_param[eos]["coeff_b"]
        assert model.fs.props[1].EoS_u == EoS_param[eos]["u"]
        assert model.fs.props[1].EoS_w == EoS_param[eos]["w"]
        assert model.fs.props[1].EoS_p == sqrt(
            EoS_param[eos]["u"] ** 2 - 4 * EoS_param[eos]["w"]
        )

        assert isinstance(model.fs.props[1].fw, Param)
        assert isinstance(model.fs.props[1].b, Param)

        assert isinstance(model.fs.props[1].compress_fact_liq_func, ExternalFunction)


class TestStateBlock_V_SRK(object):
    @pytest.fixture()
    def model(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.params = CubicParameterBlock(valid_phase="Vap")

        m.fs.params.a = Component()
        m.fs.params.b = Component()

        m.fs.params.cubic_type = CubicEoS.SRK

        m.fs.params.gas_const = Param(default=8.314462618)

        m.fs.params.pressure_crit = Param(
            m.fs.params.component_list, initialize={"a": 5e6, "b": 4e6}
        )
        m.fs.params.temperature_crit = Param(
            m.fs.params.component_list, initialize={"a": 500, "b": 600}
        )

        m.fs.params.omega = Param(
            m.fs.params.component_list, initialize={"a": 0.2, "b": 0.2}
        )

        m.fs.params.kappa = Param(
            m.fs.params.component_list,
            m.fs.params.component_list,
            initialize={
                ("a", "a"): 0.0,
                ("a", "b"): 0.0,
                ("b", "a"): 0.0,
                ("b", "b"): 0.0,
            },
        )

        return m

    @pytest.mark.unit
    def test_build_default(self, model):
        model.fs.props = model.fs.params.build_state_block([1])

        assert isinstance(model.fs.props[1].flow_mol, Var)
        assert len(model.fs.props[1].flow_mol) == 1

        assert isinstance(model.fs.props[1].mole_frac_comp, Var)
        assert len(model.fs.props[1].mole_frac_comp) == 2
        for j in model.fs.props[1].mole_frac_comp:
            assert j in ["a", "b"]

        assert isinstance(model.fs.props[1].pressure, Var)
        assert len(model.fs.props[1].pressure) == 1
        assert isinstance(model.fs.props[1].temperature, Var)
        assert len(model.fs.props[1].temperature) == 1

        assert isinstance(model.fs.props[1].flow_mol_phase, Var)
        assert len(model.fs.props[1].flow_mol_phase) == 1
        for j in model.fs.props[1].flow_mol_phase:
            assert j in ["Vap"]

        assert isinstance(model.fs.props[1].mole_frac_phase_comp, Var)
        assert len(model.fs.props[1].mole_frac_phase_comp) == 2
        for j in model.fs.props[1].mole_frac_phase_comp:
            assert j in [("Vap", "a"), ("Vap", "b")]

        assert isinstance(model.fs.props[1].total_flow_balance, Constraint)
        assert len(model.fs.props[1].total_flow_balance) == 1
        assert str(model.fs.props[1].total_flow_balance.body) == str(
            model.fs.props[1].flow_mol - model.fs.props[1].flow_mol_phase["Vap"]
        )

        assert isinstance(model.fs.props[1].component_flow_balances, Constraint)
        assert len(model.fs.props[1].component_flow_balances) == 2
        for j in model.fs.props[1].component_flow_balances:
            assert j in ["a", "b"]
            assert str(model.fs.props[1].component_flow_balances[j].body) == str(
                model.fs.props[1].mole_frac_comp[j]
                - model.fs.props[1].mole_frac_phase_comp["Vap", j]
            )

        assert isinstance(model.fs.props[1].sum_mole_frac_out, Constraint)
        assert len(model.fs.props[1].sum_mole_frac_out) == 1
        assert str(model.fs.props[1].sum_mole_frac_out.body) == str(
            sum(
                model.fs.props[1].mole_frac_comp[i]
                for i in model.fs.params.component_list
            )
        )

        assert isinstance(model.fs.props[1]._teq, Expression)
        assert len(model.fs.props[1]._teq) == 1
        assert str(model.fs.props[1]._teq.expr) == str(model.fs.props[1].temperature)

    @pytest.mark.unit
    def test_common_cubic(self, model):
        model.fs.props = model.fs.params.build_state_block([1])

        eos = CubicEoS.SRK

        assert model.fs.props[1].omegaA == EoS_param[eos]["omegaA"]
        assert model.fs.props[1].EoS_Bc == EoS_param[eos]["coeff_b"]
        assert model.fs.props[1].EoS_u == EoS_param[eos]["u"]
        assert model.fs.props[1].EoS_w == EoS_param[eos]["w"]
        assert model.fs.props[1].EoS_p == sqrt(
            EoS_param[eos]["u"] ** 2 - 4 * EoS_param[eos]["w"]
        )

        assert isinstance(model.fs.props[1].fw, Param)
        assert isinstance(model.fs.props[1].b, Param)

        assert isinstance(model.fs.props[1].compress_fact_vap_func, ExternalFunction)
