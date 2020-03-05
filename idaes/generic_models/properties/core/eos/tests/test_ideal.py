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
Tests for ideal equation of state methods

Author: Andrew Lee
"""
import pytest
from sys import modules

from pyomo.environ import ConcreteModel, Var

from idaes.core import declare_process_block_class
from idaes.generic_models.properties.core.eos import ideal
from idaes.generic_models.properties.core.generic.generic_property import (
        GenericParameterData)
from idaes.generic_models.properties.core.generic.tests import dummy_eos
from idaes.core.util.constants import Constants as const
from idaes.core.util.exceptions import PropertyNotSupportedError


# Dummy method for property method calls
def dummy_call(b, j, T):
    return 42


@declare_process_block_class("DummyParameterBlock")
class DummyParameterData(GenericParameterData):
    def configure(self):
        self.configured = True

    def parameters(self):
        self.parameters_set = True


def define_state(b):
    b.state_defined = True


@pytest.fixture()
def m():
    m = ConcreteModel()

    # Dummy params block
    m.params = DummyParameterBlock(default={
                "component_list": ["a", "b", "c"],
                "phase_list": ["Vap", "Liq"],
                "state_definition": modules[__name__],
                "equation_of_state": {"Vap": dummy_eos,
                                      "Liq": dummy_eos}})

    m.props = m.params.state_block_class([1],
                                         default={"defined_state": False,
                                         "parameters": m.params})

    # Add common variables
    m.props[1].pressure = Var(initialize=101325)
    m.props[1].temperature = Var(initialize=300)
    m.props[1]._teq = Var(initialize=300)
    m.props[1].mole_frac_phase_comp = Var(m.params.phase_list,
                                          m.params.component_list,
                                          initialize=0.5)

    return m


def test_common(m):
    assert ideal.common(m.props) is None


def test_dens_mass_phase(m):
    m.props[1].dens_mol_phase = Var(m.params.phase_list)
    m.props[1].mw_phase = Var(m.params.phase_list)

    for p in m.params.phase_list:
        assert str(ideal.dens_mass_phase(m.props[1], p)) == str(
                m.props[1].dens_mol_phase[p]*m.props[1].mw_phase[p])


def test_dens_mol_phase_liq(m):
    m.params.config.dens_mol_liq_comp = dummy_call

    assert str(ideal.dens_mol_phase(m.props[1], "Liq")) == str(
        sum(m.props[1].mole_frac_phase_comp["Liq", j]*42
            for j in m.params.component_list))


def test_dens_mol_phase_vap(m):
    m.params.gas_const = Var()

    assert str(ideal.dens_mol_phase(m.props[1], "Vap")) == (
            str(m.props[1].pressure)+"/(8.314462618*J/mol/K*" +
            str(m.props[1].temperature)+")")


def test_dens_mol_phase_invalid_phase(m):
    with pytest.raises(PropertyNotSupportedError):
        ideal.dens_mol_phase(m.props[1], "foo")


def test_enth_mol_phase(m):
    m.props[1].enth_mol_phase_comp = Var(m.params.phase_list,
                                         m.params.component_list)

    for p in m.params.phase_list:
        assert str(ideal.enth_mol_phase(m.props[1], p)) == str(
            sum(m.props[1].mole_frac_phase_comp[p, j] *
                m.props[1].enth_mol_phase_comp[p, j]
                for j in m.params.component_list))


def test_enth_mol_phase_comp(m):
    m.params.config.enth_mol_liq_comp = dummy_call
    m.params.config.enth_mol_ig_comp = dummy_call

    for p in m.params.phase_list:
        assert str(ideal.enth_mol_phase_comp(m.props[1], p, "foo")) == str(42)


def test_enth_mol_phase_invalid_phase(m):
    with pytest.raises(PropertyNotSupportedError):
        ideal.enth_mol_phase_comp(m.props[1], "foo", "bar")


def test_entr_mol_phase(m):
    m.props[1].entr_mol_phase_comp = Var(m.params.phase_list,
                                         m.params.component_list)

    for p in m.params.phase_list:
        assert str(ideal.entr_mol_phase(m.props[1], p)) == str(
            sum(m.props[1].mole_frac_phase_comp[p, j] *
                m.props[1].entr_mol_phase_comp[p, j]
                for j in m.params.component_list))


def test_entr_mol_phase_comp(m):
    m.params.config.entr_mol_liq_comp = dummy_call
    m.params.config.entr_mol_ig_comp = dummy_call

    for p in m.params.phase_list:
        assert str(ideal.entr_mol_phase_comp(m.props[1], p, "foo")) == str(42)


def test_entr_mol_phase_invalid_phase(m):
    with pytest.raises(PropertyNotSupportedError):
        ideal.entr_mol_phase_comp(m.props[1], "foo", "bar")


def test_fug_phase_comp_liq(m):
    m.params.config.pressure_sat_comp = dummy_call

    for j in m.params.component_list:
        assert str(ideal.fug_phase_comp(m.props[1], "Liq", j)) == str(
            m.props[1].mole_frac_phase_comp["Liq", j]*42)


def test_fug_phase_comp_vap(m):
    for j in m.params.component_list:
        assert str(ideal.fug_phase_comp(m.props[1], "Vap", j)) == str(
            m.props[1].mole_frac_phase_comp["Vap", j]*m.props[1].pressure)


def test_fug_phase_comp_invalid_phase(m):
    with pytest.raises(PropertyNotSupportedError):
        ideal.fug_phase_comp(m.props[1], "foo", "bar")


def test_fug_coeff_phase_comp(m):
    for p in m.params.phase_list:
        for j in m.params.component_list:
            assert ideal.fug_coeff_phase_comp(m.props[1], p, j) == 1


def test_fug_coeff_phase_comp_invalid_phase(m):
    with pytest.raises(PropertyNotSupportedError):
        ideal.fug_coeff_phase_comp(m.props[1], "foo", "bar")


def test_gibbs_mol_phase(m):
    m.props[1].gibbs_mol_phase_comp = Var(m.params.phase_list,
                                          m.params.component_list)

    for p in m.params.phase_list:
        assert str(ideal.gibbs_mol_phase(m.props[1], p)) == str(
            sum(m.props[1].mole_frac_phase_comp[p, j] *
                m.props[1].gibbs_mol_phase_comp[p, j]
                for j in m.params.component_list))


def test_gibbs_mol_phase_comp(m):
    m.props[1].enth_mol_phase_comp = Var(m.params.phase_list,
                                         m.params.component_list)
    m.props[1].entr_mol_phase_comp = Var(m.params.phase_list,
                                         m.params.component_list)

    for p in m.params.phase_list:
        for j in m.params.component_list:
            assert str(ideal.gibbs_mol_phase_comp(m.props[1], p, j)) == str(
                    m.props[1].enth_mol_phase_comp[p, j] -
                    m.props[1].entr_mol_phase_comp[p, j] *
                    m.props[1].temperature)
