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
Tests for Cubic equation of state methods

Author: Andrew Lee
"""
import pytest
from sys import modules
import os

from pyomo.environ import (ConcreteModel,
                           Expression,
                           ExternalFunction,
                           Param,
                           sqrt,
                           value,
                           Var,
                           units as pyunits)
from pyomo.core.expr.numeric_expr import ExternalFunctionExpression

from idaes.core import (declare_process_block_class,
                        LiquidPhase, VaporPhase, SolidPhase)
from idaes.generic_models.properties.core.eos.ceos import Cubic, CubicType
from idaes.generic_models.properties.core.generic.generic_property import (
        GenericParameterData)
from idaes.core.util.exceptions import PropertyNotSupportedError
from idaes.core.util.constants import Constants as const
from idaes import bin_directory


# Dummy method for property method calls
def dummy_call(b, j, T):
    return 42


# Dummy method to avoid errors when setting metadata dict
def set_metadata(b):
    pass


# Dummy methods for dummied submodules
class dummy_pe():
    def return_expression(b, *args):
        # Return a dummy expression for the constraint
        return b.temperature == 100


def phase_equil(b, *args):
    pass


@declare_process_block_class("DummyParameterBlock")
class DummyParameterData(GenericParameterData):
    def configure(self):
        self.configured = True

    def parameters(self):
        self.parameters_set = True


def define_state(b):
    # Add common variables
    b.pressure = Var(initialize=101325)
    b.temperature = Var(initialize=300)
    b.mole_frac_phase_comp = Var(b.params.phase_list,
                                 b.params.component_list,
                                 initialize=0.5)

    # Set mole fractions for each phase so they can be distinguished
    if "Sol" not in b.params.phase_list:
        b.mole_frac_phase_comp["Vap", "a"].value = 0.1
        b.mole_frac_phase_comp["Vap", "b"].value = 0.2
        b.mole_frac_phase_comp["Vap", "c"].value = 0.7
        b.mole_frac_phase_comp["Liq", "a"].value = 0.6
        b.mole_frac_phase_comp["Liq", "b"].value = 0.3
        b.mole_frac_phase_comp["Liq", "c"].value = 0.1


@pytest.fixture()
def m():
    m = ConcreteModel()

    # Dummy params block
    m.params = DummyParameterBlock(default={
                "components": {
                    "a": {"phase_equilibrium_form": {("Vap", "Liq"): dummy_pe},
                          "parameter_data": {
                              "omega": 0.1,
                              "pressure_crit": 1e5,
                              "temperature_crit": 100}},
                    "b": {"phase_equilibrium_form": {("Vap", "Liq"): dummy_pe},
                          "parameter_data": {
                              "omega": 0.2,
                              "pressure_crit": 2e5,
                              "temperature_crit": 200}},
                    "c": {"phase_equilibrium_form": {("Vap", "Liq"): dummy_pe},
                          "parameter_data": {
                              "omega": 0.3,
                              "pressure_crit": 3e5,
                              "temperature_crit": 300}}},
                "phases": {
                    "Vap": {"type": VaporPhase,
                            "equation_of_state": Cubic,
                            "equation_of_state_options": {
                                "type": CubicType.PR}},
                    "Liq": {"type": LiquidPhase,
                            "equation_of_state": Cubic,
                            "equation_of_state_options": {
                                "type": CubicType.PR}}},
                "base_units": {"time": pyunits.s,
                               "length": pyunits.m,
                               "mass": pyunits.kg,
                               "amount": pyunits.mol,
                               "temperature": pyunits.K},
                "state_definition": modules[__name__],
                "pressure_ref": 1e5,
                "temperature_ref": 300,
                "phases_in_equilibrium": [("Vap", "Liq")],
                "phase_equilibrium_state": {("Vap", "Liq"): modules[__name__]},
                "parameter_data": {"PR_kappa": {("a", "a"): 0.000,
                                                ("a", "b"): 0.000,
                                                ("a", "c"): 0.000,
                                                ("b", "a"): 0.000,
                                                ("b", "b"): 0.000,
                                                ("b", "c"): 0.000,
                                                ("c", "a"): 0.000,
                                                ("c", "b"): 0.000,
                                                ("c", "c"): 0.000}}})

    m.props = m.params.state_block_class(
        [1], default={"defined_state": False,
                      "parameters": m.params,
                      "has_phase_equilibrium": True})

    # Set a distinct value for _teq so it can be distinguished from temperature
    m.props[1]._teq[("Vap", "Liq")].value = 100

    return m


@pytest.fixture()
def m_sol():
    m = ConcreteModel()

    # Dummy params block
    m.params = DummyParameterBlock(default={
                "components": {
                    "a": {"phase_equilibrium_form": {("Sol", "Liq"): dummy_pe},
                          "parameter_data": {
                              "omega": 0.1,
                              "pressure_crit": 1e5,
                              "temperature_crit": 100}},
                    "b": {"phase_equilibrium_form": {("Sol", "Liq"): dummy_pe},
                          "parameter_data": {
                              "omega": 0.2,
                              "pressure_crit": 2e5,
                              "temperature_crit": 200}},
                    "c": {"phase_equilibrium_form": {("Sol", "Liq"): dummy_pe},
                          "parameter_data": {
                              "omega": 0.3,
                              "pressure_crit": 3e5,
                              "temperature_crit": 300}}},
                "phases": {
                    "Sol": {"type": SolidPhase,
                            "equation_of_state": Cubic,
                            "equation_of_state_options": {
                                "type": CubicType.PR}},
                    "Liq": {"type": LiquidPhase,
                            "equation_of_state": Cubic,
                            "equation_of_state_options": {
                                "type": CubicType.PR}}},
                "base_units": {"time": pyunits.s,
                               "length": pyunits.m,
                               "mass": pyunits.kg,
                               "amount": pyunits.mol,
                               "temperature": pyunits.K},
                "state_definition": modules[__name__],
                "pressure_ref": 1e5,
                "temperature_ref": 300,
                "phases_in_equilibrium": [("Sol", "Liq")],
                "phase_equilibrium_state": {("Sol", "Liq"): modules[__name__]},
                "parameter_data": {"PR_kappa": {("a", "a"): 0.000,
                                                ("a", "b"): 0.000,
                                                ("a", "c"): 0.000,
                                                ("b", "a"): 0.000,
                                                ("b", "b"): 0.000,
                                                ("b", "c"): 0.000,
                                                ("c", "a"): 0.000,
                                                ("c", "b"): 0.000,
                                                ("c", "c"): 0.000}}})

    m.props = m.params.state_block_class(
        [1], default={"defined_state": False,
                      "parameters": m.params,
                      "has_phase_equilibrium": True})

    return m


# Expected values for A and B
Al = 0.0419062
Bl = 0.0262770
Av = 0.1180005
Bv = 0.0262769
Zl = 0.9870125
Zv = 0.9067390

Al_eq = 0.8325554
Bl_eq = 0.0788309
Av_eq = 1.9805792
Bv_eq = 0.0788309

Hl = -212.211
Hv = -779.305
Ul = Hl - 8.314*300*(Zl-1)
Uv = Hv - 8.314*300*(Zv-1)

# Set path to root finder .so file
_so = os.path.join(bin_directory, "cubic_roots.so")
f_Zl = ExternalFunction(library=_so, function="ceos_z_liq")
f_Zv = ExternalFunction(library=_so, function="ceos_z_vap")


@pytest.mark.unit
def test_common(m):
    # Test cubic components
    assert isinstance(m.props[1].PR_fw, Expression)
    assert len(m.props[1].PR_fw) == len(m.params.component_list)
    for i in m.params.component_list:
        omega = m.params.get_component(i).omega
        assert value(m.props[1].PR_fw[i]) == value(
            0.37464 + 1.54226*omega - 0.26992*omega**2)

    assert isinstance(m.props[1].PR_a, Expression)
    assert len(m.props[1].PR_a) == len(m.params.component_list)
    for i in m.params.component_list:
        Tc = m.params.get_component(i).temperature_crit
        Pc = m.params.get_component(i).pressure_crit
        assert pytest.approx(value(m.props[1].PR_a[i]),
                             rel=1e-5) == value(
            0.45724*((const.gas_constant*Tc)**2/Pc) *
            (((1+m.props[1].PR_fw[i] *
               (1-sqrt(m.props[1].temperature/Tc)))**2)))

    assert isinstance(m.props[1].PR_b, Expression)
    assert len(m.props[1].PR_b) == len(m.params.component_list)
    for i in m.params.component_list:
        Tc = m.params.get_component(i).temperature_crit
        Pc = m.params.get_component(i).pressure_crit
        assert value(m.props[1].PR_b[i]) == value(
            0.07780*const.gas_constant*Tc/Pc)

    assert isinstance(m.props[1].PR_am, Expression)
    assert len(m.props[1].PR_am) == len(m.params.phase_list)
    for p in m.params.phase_list:
        assert pytest.approx(value(m.props[1].PR_am[p]),
                             rel=1e-5) == value(
            sum(sum(m.props[1].mole_frac_phase_comp[p, i] *
                    m.props[1].mole_frac_phase_comp[p, j] *
                    sqrt(m.props[1].PR_a[i]*m.props[1].PR_a[j]) *
                    (1-m.params.PR_kappa[i, j])
                    for j in m.params.component_list)
                for i in m.params.component_list))

    assert isinstance(m.props[1].PR_bm, Expression)
    assert len(m.props[1].PR_bm) == len(m.params.phase_list)
    for p in m.params.phase_list:
        assert pytest.approx(value(m.props[1].PR_bm[p]),
                             rel=1e-5) == value(
            sum(m.props[1].mole_frac_phase_comp[p, i] *
                m.props[1].PR_b[i]
                for i in m.props[1].params.component_list))

    assert isinstance(m.props[1].PR_A, Expression)
    assert len(m.props[1].PR_A) == len(m.params.phase_list)
    for p in m.params.phase_list:
        assert pytest.approx(value(m.props[1].PR_A[p]),
                             rel=1e-5) == value(
                m.props[1].PR_am[p]*m.props[1].pressure /
                (const.gas_constant*m.props[1].temperature)**2)
    assert pytest.approx(value(m.props[1].PR_A["Liq"]), rel=1e-5) == Al
    assert pytest.approx(value(m.props[1].PR_A["Vap"]), rel=1e-5) == Av

    assert isinstance(m.props[1].PR_B, Expression)
    assert len(m.props[1].PR_B) == len(m.params.phase_list)
    for p in m.params.phase_list:
        assert pytest.approx(value(m.props[1].PR_B[p]),
                             rel=1e-5) == value(
                m.props[1].PR_bm[p]*m.props[1].pressure /
                (const.gas_constant*m.props[1].temperature))
    assert pytest.approx(value(m.props[1].PR_B["Liq"]), rel=1e-5) == Bl
    assert pytest.approx(value(m.props[1].PR_B["Vap"]), rel=1e-5) == Bv

    assert isinstance(m.props[1].PR_delta, Expression)
    assert len(m.props[1].PR_delta) == 6
    for p in m.params.phase_list:
        for i in m.params.component_list:
            assert pytest.approx(value(m.props[1].PR_delta[p, i]),
                                 rel=1e-5) == value(
                2*sqrt(m.props[1].PR_a[i])/m.props[1].PR_am[p] *
                sum(m.props[1].mole_frac_phase_comp[p, j] *
                    sqrt(m.props[1].PR_a[j]) *
                    (1-m.params.PR_kappa[i, j])
                    for j in m.params.component_list))

    assert isinstance(m.props[1].PR_dadT, Expression)
    assert len(m.props[1].PR_dadT) == len(m.params.phase_list)
    for p in m.params.phase_list:
        assert pytest.approx(value(m.props[1].PR_dadT[p]),
                             rel=1e-5) == value(
            -((const.gas_constant/2)*sqrt(0.45724) *
              sum(sum(m.props[1].mole_frac_phase_comp[p, i] *
                      m.props[1].mole_frac_phase_comp[p, j] *
                      (1-m.params.PR_kappa[i, j]) *
                      (m.props[1].PR_fw[j] *
                       sqrt(m.props[1].PR_a[i] *
                            m.params.get_component(j).temperature_crit /
                            m.params.get_component(j).pressure_crit) +
                       m.props[1].PR_fw[i] *
                       sqrt(m.props[1].PR_a[j] *
                            m.params.get_component(i).temperature_crit /
                            m.params.get_component(i).pressure_crit))
                      for j in m.params.component_list)
                  for i in m.params.component_list) /
              sqrt(m.props[1].temperature)))

    # Test equilibrium state Expressions
    assert isinstance(m.props[1]._PR_a_eq, Expression)
    assert len(m.props[1]._PR_a_eq) == len(m.params.component_list)
    for i in m.props[1]._PR_a_eq:
        Tc = m.params.get_component(i[2]).temperature_crit
        Pc = m.params.get_component(i[2]).pressure_crit
        assert pytest.approx(value(m.props[1]._PR_a_eq[i]),
                             rel=1e-5) == value(
            0.45724*((const.gas_constant*Tc)**2/Pc) *
            (((1+m.props[1].PR_fw[i[2]] *
               (1-sqrt(m.props[1]._teq[i[0], i[1]]/Tc)))**2)))

    assert isinstance(m.props[1]._PR_am_eq, Expression)
    assert len(m.props[1]._PR_am_eq) == len(m.params.phase_list)
    for idx in m.props[1]._PR_am_eq:
        assert pytest.approx(value(m.props[1]._PR_am_eq[idx]),
                             rel=1e-5) == value(
            sum(sum(m.props[1].mole_frac_phase_comp[idx[2], i] *
                    m.props[1].mole_frac_phase_comp[idx[2], j] *
                    sqrt(m.props[1]._PR_a_eq[idx[0], idx[1], i] *
                         m.props[1]._PR_a_eq[idx[0], idx[1], j]) *
                    (1-m.params.PR_kappa[i, j])
                    for j in m.params.component_list)
                for i in m.params.component_list))

    assert isinstance(m.props[1]._PR_A_eq, Expression)
    assert len(m.props[1]._PR_A_eq) == len(m.params.phase_list)
    for i in m.props[1]._PR_A_eq:
        assert pytest.approx(value(m.props[1]._PR_A_eq[i]),
                             rel=1e-5) == value(
                    m.props[1]._PR_am_eq[i] *
                    m.props[1].pressure /
                    (const.gas_constant *
                     m.props[1]._teq[i[0], i[1]])**2)
    assert pytest.approx(value(m.props[1]._PR_A_eq[("Vap", "Liq"), "Liq"]),
                         rel=1e-5) == Al_eq
    assert pytest.approx(value(m.props[1]._PR_A_eq[("Vap", "Liq"), "Vap"]),
                         rel=1e-5) == Av_eq

    assert isinstance(m.props[1]._PR_B_eq, Expression)
    assert len(m.props[1]._PR_B_eq) == len(m.params.phase_list)
    for i in m.props[1]._PR_B_eq:
        assert pytest.approx(value(m.props[1]._PR_B_eq[i]),
                             rel=1e-5) == value(
                    m.props[1].PR_bm[i[2]] *
                    m.props[1].pressure /
                    (const.gas_constant *
                     m.props[1]._teq[i[0], i[1]]))
    assert pytest.approx(value(m.props[1]._PR_B_eq[("Vap", "Liq"), "Liq"]),
                         rel=1e-5) == Bl_eq
    assert pytest.approx(value(m.props[1]._PR_B_eq[("Vap", "Liq"), "Vap"]),
                         rel=1e-5) == Bv_eq

    assert isinstance(m.props[1]._PR_delta_eq, Expression)
    assert len(m.props[1]._PR_delta_eq) == 6
    for idx in m.props[1]._PR_delta_eq:
        assert pytest.approx(value(m.props[1]._PR_delta_eq[idx]),
                             rel=1e-5) == value(
            2*sqrt(m.props[1]._PR_a_eq[idx[0], idx[1], idx[3]]) /
            m.props[1]._PR_am_eq[idx[0], idx[1], idx[2]] *
            sum(m.props[1].mole_frac_phase_comp[idx[2], j] *
                sqrt(m.props[1]._PR_a_eq[idx[0], idx[1], j]) *
                (1-m.params.PR_kappa[idx[3], j])
                for j in m.params.component_list))

    # Check for external function components
    assert isinstance(m.props[1]._PR_ext_func_param, Param)
    assert m.props[1]._PR_ext_func_param.value == 0  # 0 == PR
    assert isinstance(m.props[1]._PR_proc_Z_liq, ExternalFunction)
    assert isinstance(m.props[1]._PR_proc_Z_vap, ExternalFunction)


@pytest.mark.unit
def test_compress_fact_phase_Liq(m):
    assert isinstance(Cubic.compress_fact_phase(m.props[1], "Liq"),
                      ExternalFunctionExpression)
    assert pytest.approx(value(f_Zl(0, Al, Bl)), rel=1e-5) == value(
        Cubic.compress_fact_phase(m.props[1], "Liq"))
    assert pytest.approx(value(
        Cubic.compress_fact_phase(m.props[1], "Liq")), rel=1e-5) == Zl


@pytest.mark.unit
def test_compress_fact_phase_Vap(m):
    assert isinstance(Cubic.compress_fact_phase(m.props[1], "Vap"),
                      ExternalFunctionExpression)
    assert pytest.approx(value(f_Zv(0, Av, Bv)), rel=1e-5) == value(
        Cubic.compress_fact_phase(m.props[1], "Vap"))
    assert pytest.approx(value(
        Cubic.compress_fact_phase(m.props[1], "Vap")), rel=1e-5) == Zv


@pytest.mark.unit
def test_compress_fact_phase_invalid_phase(m_sol):
    with pytest.raises(PropertyNotSupportedError):
        Cubic.compress_fact_phase(m_sol.props[1], "Sol")


@pytest.mark.unit
def test_dens_mass_phase(m):
    m.props[1].dens_mol_phase = Var(m.params.phase_list)
    m.props[1].mw_phase = Var(m.params.phase_list)

    for p in m.params.phase_list:
        assert str(Cubic.dens_mass_phase(m.props[1], p)) == str(
                m.props[1].dens_mol_phase[p]*m.props[1].mw_phase[p])


@pytest.mark.unit
def test_dens_mol_phase(m):
    assert value(Cubic.dens_mol_phase(m.props[1], "Vap")) == pytest.approx(
            44.800, rel=1e-3)
    assert value(Cubic.dens_mol_phase(m.props[1], "Liq")) == pytest.approx(
            41.157, rel=1e-3)


@pytest.mark.unit
def test_dens_mol_phase_invalid_phase(m_sol):
    with pytest.raises(PropertyNotSupportedError):
        Cubic.dens_mol_phase(m_sol.props[1], "Sol")


@pytest.mark.unit
def test_energy_internal_mol_phase(m):
    for j in m.params.component_list:
        m.params.config.include_enthalpy_of_formation = False
        m.params.get_component(j).config.enth_mol_liq_comp = dummy_call
        m.params.get_component(j).config.enth_mol_ig_comp = dummy_call

    m.props[1].energy_internal_mol_phase_comp = Var(
        m.params.phase_list, m.params.component_list, initialize=1)

    assert pytest.approx(Uv, rel=1e-4) == value(
        Cubic.energy_internal_mol_phase(m.props[1], "Vap"))
    assert pytest.approx(Ul, rel=1e-4) == value(
        Cubic.energy_internal_mol_phase(m.props[1], "Liq"))


@pytest.mark.unit
def test_energy_internal_mol_phase_comp(m):
    for j in m.params.component_list:
        m.params.config.include_enthalpy_of_formation = False
        m.params.get_component(j).config.enth_mol_liq_comp = dummy_call
        m.params.get_component(j).config.enth_mol_ig_comp = dummy_call

        assert pytest.approx(Ul, rel=1e-4) == value(
            Cubic.energy_internal_mol_phase_comp(m.props[1], "Liq", j))
        assert pytest.approx(Uv, rel=1e-4) == value(
            Cubic.energy_internal_mol_phase_comp(m.props[1], "Vap", j))


@pytest.mark.unit
def test_energy_internal_mol_phase_invalid_phase(m_sol):
    with pytest.raises(PropertyNotSupportedError):
        Cubic.energy_internal_mol_phase_comp(m_sol.props[1], "Sol", "foo")


@pytest.mark.unit
def test_enth_mol_phase(m):
    for j in m.params.component_list:
        m.params.get_component(j).config.enth_mol_liq_comp = dummy_call
        m.params.get_component(j).config.enth_mol_ig_comp = dummy_call

    m.props[1].enth_mol_phase_comp = Var(m.params.phase_list,
                                         m.params.component_list,
                                         initialize=1)

    assert pytest.approx(value(
        Cubic.enth_mol_phase(m.props[1], "Vap")), rel=1e-5) == Hv
    assert pytest.approx(value(
        Cubic.enth_mol_phase(m.props[1], "Liq")), rel=1e-5) == Hl


@pytest.mark.unit
def test_enth_mol_phase_comp(m):
    for j in m.params.component_list:
        m.params.get_component(j).config.enth_mol_liq_comp = dummy_call
        m.params.get_component(j).config.enth_mol_ig_comp = dummy_call

        assert pytest.approx(Hl, rel=1e-5) == value(
            Cubic.enth_mol_phase_comp(m.props[1], "Liq", j))
        assert pytest.approx(Hv, rel=1e-5) == value(
            Cubic.enth_mol_phase_comp(m.props[1], "Vap", j))


@pytest.mark.unit
def test_enth_mol_phase_invalid_phase(m_sol):
    with pytest.raises(PropertyNotSupportedError):
        Cubic.enth_mol_phase_comp(m_sol.props[1], "Sol", "foo")


@pytest.mark.unit
def test_entr_mol_phase(m):
    for j in m.params.component_list:
        m.params.get_component(j).config.entr_mol_liq_comp = dummy_call
        m.params.get_component(j).config.entr_mol_ig_comp = dummy_call

    m.props[1].entr_mol_phase_comp = Var(m.params.phase_list,
                                         m.params.component_list,
                                         initialize=1)

    assert pytest.approx(value(
        Cubic.entr_mol_phase(m.props[1], "Vap")), rel=1e-5) == 39.9218
    assert pytest.approx(value(
        Cubic.entr_mol_phase(m.props[1], "Liq")), rel=1e-5) == 41.1621


@pytest.mark.unit
def test_entr_mol_phase_comp(m):
    entr = {("Liq", "a"): 45.4093,
            ("Vap", "a"): 59.0666,
            ("Liq", "b"): 51.1725,
            ("Vap", "b"): 53.3035,
            ("Liq", "c"): 60.3069,
            ("Vap", "c"): 42.8875}
    for j in m.params.component_list:
        m.params.get_component(j).config.entr_mol_liq_comp = dummy_call
        m.params.get_component(j).config.entr_mol_ig_comp = dummy_call

        assert pytest.approx(entr[("Liq", j)], rel=1e-5) == value(
            Cubic.entr_mol_phase_comp(m.props[1], "Liq", j))
        assert pytest.approx(entr[("Vap", j)], rel=1e-5) == value(
            Cubic.entr_mol_phase_comp(m.props[1], "Vap", j))


@pytest.mark.unit
def test_entr_mol_phase_invalid_phase(m_sol):
    with pytest.raises(PropertyNotSupportedError):
        Cubic.entr_mol_phase_comp(m_sol.props[1], "Sol", "foo")


@pytest.mark.component
def test_fug_phase_comp(m):
    for p in m.params.phase_list:
        for j in m.params.component_list:
            assert str(Cubic.fug_phase_comp(
                            m.props[1], p, j)) == str(
                m.props[1].mole_frac_phase_comp[p, j] *
                m.props[1].pressure *
                m.props[1].fug_coeff_phase_comp[p, j])


@pytest.mark.unit
def test_fug_phase_comp_invalid_phase(m_sol):
    with pytest.raises(PropertyNotSupportedError):
        Cubic.fug_phase_comp(m_sol.props[1], "Sol", "foo")


@pytest.mark.component
def test_fug_phase_comp_eq(m):
    for p in m.params.phase_list:
        for j in m.params.component_list:
            assert str(Cubic.fug_phase_comp_eq(
                            m.props[1], p, j, ("Vap", "Liq"))) == str(
                m.props[1].mole_frac_phase_comp[p, j] *
                m.props[1].pressure *
                Cubic.fug_coeff_phase_comp_eq(
                    m.props[1], p, j, ("Vap", "Liq")))


@pytest.mark.unit
def test_fug_phase_comp_eq_invalid_phase(m_sol):
    with pytest.raises(PropertyNotSupportedError):
        Cubic.fug_phase_comp_eq(m_sol.props[1], "Sol", "foo", ("Vap", "Liq"))


@pytest.mark.unit
def test_fug_coeff_phase_comp_Liq(m):
    assert pytest.approx(1.01213, rel=1e-5) == value(
        Cubic.fug_coeff_phase_comp(m.props[1], "Liq", "a"))
    assert pytest.approx(0.95919, rel=1e-5) == value(
        Cubic.fug_coeff_phase_comp(m.props[1], "Liq", "b"))
    assert pytest.approx(0.91356, rel=1e-5) == value(
        Cubic.fug_coeff_phase_comp(m.props[1], "Liq", "c"))


@pytest.mark.unit
def test_fug_coeff_phase_comp_Vap(m):
    assert pytest.approx(1.05952, rel=1e-5) == value(
            Cubic.fug_coeff_phase_comp(m.props[1], "Vap", "a"))
    assert pytest.approx(0.96070, rel=1e-5) == value(
        Cubic.fug_coeff_phase_comp(m.props[1], "Vap", "b"))
    assert pytest.approx(0.87903, rel=1e-5) == value(
        Cubic.fug_coeff_phase_comp(m.props[1], "Vap", "c"))


@pytest.mark.unit
def test_fug_coeff_phase_comp_invalid_phase(m_sol):
    with pytest.raises(PropertyNotSupportedError):
        Cubic.fug_coeff_phase_comp(m_sol.props[1], "Sol", "foo")


@pytest.mark.unit
def test_fug_coeff_phase_comp_eq_Liq(m):
    assert pytest.approx(1.22431, rel=1e-5) == value(
        Cubic.fug_coeff_phase_comp_eq(m.props[1], "Liq", "a", ("Vap", "Liq")))
    assert pytest.approx(0.0049604, rel=1e-5) == value(
        Cubic.fug_coeff_phase_comp_eq(m.props[1], "Liq", "b", ("Vap", "Liq")))
    assert pytest.approx(3.19128e-05, rel=1e-5) == value(
        Cubic.fug_coeff_phase_comp_eq(m.props[1], "Liq", "c", ("Vap", "Liq")))


@pytest.mark.unit
def test_fug_coeff_phase_comp_eq_Vap(m):
    assert pytest.approx(86.9140, rel=1e-5) == value(
        Cubic.fug_coeff_phase_comp_eq(m.props[1], "Vap", "a", ("Vap", "Liq")))
    assert pytest.approx(0.0049878, rel=1e-5) == value(
        Cubic.fug_coeff_phase_comp_eq(m.props[1], "Vap", "b", ("Vap", "Liq")))
    assert pytest.approx(6.49776e-07, rel=1e-5) == value(
        Cubic.fug_coeff_phase_comp_eq(m.props[1], "Vap", "c", ("Vap", "Liq")))


@pytest.mark.unit
def test_fug_coeff_phase_comp_eq_invalid_phase(m_sol):
    with pytest.raises(PropertyNotSupportedError):
        Cubic.fug_coeff_phase_comp_eq(
            m_sol.props[1], "Sol", "foo", ("Vap", "Liq"))


@pytest.mark.unit
def test_gibbs_mol_phase(m):
    m.props[1].enth_mol_phase = Var(m.params.phase_list)
    m.props[1].entr_mol_phase = Var(m.params.phase_list)

    for p in m.params.phase_list:
        assert str(Cubic.gibbs_mol_phase(m.props[1], p)) == str(
            m.props[1].enth_mol_phase[p] -
            m.props[1].entr_mol_phase[p]*m.props[1].temperature)


@pytest.mark.unit
def test_gibbs_mol_phase_comp(m):
    m.props[1].enth_mol_phase_comp = Var(m.params.phase_list,
                                         m.params.component_list)
    m.props[1].entr_mol_phase_comp = Var(m.params.phase_list,
                                         m.params.component_list)

    for p in m.params.phase_list:
        for j in m.params.component_list:
            assert str(Cubic.gibbs_mol_phase_comp(m.props[1], p, j)) == str(
                    m.props[1].enth_mol_phase_comp[p, j] -
                    m.props[1].entr_mol_phase_comp[p, j] *
                    m.props[1].temperature)


@pytest.mark.unit
def test_vol_mol_phase(m):
    assert value(Cubic.vol_mol_phase(m.props[1], "Vap")) == pytest.approx(
            1/44.800, rel=1e-3)
    assert value(Cubic.vol_mol_phase(m.props[1], "Liq")) == pytest.approx(
            1/41.157, rel=1e-3)


@pytest.mark.unit
def test_vol_mol_phase_invalid_phase(m_sol):
    with pytest.raises(PropertyNotSupportedError):
        Cubic.vol_mol_phase(m_sol.props[1], "Sol")
