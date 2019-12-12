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

from pyomo.environ import Block, ConcreteModel, Set, Var
from pyomo.common.config import ConfigBlock

from idaes.property_models.core.eos import ideal

from idaes.core.util.exceptions import PropertyNotSupportedError


# Dummy method for property method calls
def dummy_call(b, j, T):
    return 42


@pytest.fixture()
def m():
    m = ConcreteModel()

    # Dummy params block
    m._params = Block()
    m._params.config = ConfigBlock(implicit=True)

    m._params.component_list = Set(initialize=["a", "b", "c"])
    m._params.phase_list = Set(initialize=["Vap", "Liq"])

    # Add common variables
    m.pressure = Var()
    m.temperature = Var()
    m.mole_frac_phase_comp = Var(m._params.phase_list,
                                 m._params.component_list)

    return m


def test_common(m):
    assert ideal.common(m) is None


def test_dens_mass_phase(m):
    m.dens_mol_phase = Var(m._params.phase_list)
    m.mw_phase = Var(m._params.phase_list)

    for p in m._params.phase_list:
        assert str(ideal.dens_mass_phase(m, p)) == str(
                m.dens_mol_phase[p]*m.mw_phase[p])


def test_dens_mol_phase_liq(m):
    m._params.config.dens_mol_liq_comp = dummy_call

    assert str(ideal.dens_mol_phase(m, "Liq")) == str(
        sum(m.mole_frac_phase_comp["Liq", j]*42
            for j in m._params.component_list))


def test_dens_mol_phase_vap(m):
    m._params.gas_const = Var()

    assert str(ideal.dens_mol_phase(m, "Vap")) == str(
        m.pressure/(m._params.gas_const*m.temperature))


def test_dens_mol_phase_invalid_phase(m):
    with pytest.raises(PropertyNotSupportedError):
        ideal.dens_mol_phase(m, "foo")


def test_enth_mol_phase(m):
    m.enth_mol_phase_comp = Var(m._params.phase_list,
                                m._params.component_list)

    for p in m._params.phase_list:
        assert str(ideal.enth_mol_phase(m, p)) == str(
            sum(m.mole_frac_phase_comp[p, j]*m.enth_mol_phase_comp[p, j]
                for j in m._params.component_list))


def test_enth_mol_phase_comp(m):
    m._params.config.enth_mol_liq_comp = dummy_call
    m._params.config.enth_mol_ig_comp = dummy_call

    for p in m._params.phase_list:
        assert str(ideal.enth_mol_phase_comp(m, p, "foo")) == str(42)


def test_enth_mol_phase_invalid_phase(m):
    with pytest.raises(PropertyNotSupportedError):
        ideal.enth_mol_phase_comp(m, "foo", "bar")


def test_entr_mol_phase(m):
    m.entr_mol_phase_comp = Var(m._params.phase_list,
                                m._params.component_list)

    for p in m._params.phase_list:
        assert str(ideal.entr_mol_phase(m, p)) == str(
            sum(m.mole_frac_phase_comp[p, j]*m.entr_mol_phase_comp[p, j]
                for j in m._params.component_list))


def test_entr_mol_phase_comp(m):
    m._params.config.entr_mol_liq_comp = dummy_call
    m._params.config.entr_mol_ig_comp = dummy_call

    for p in m._params.phase_list:
        assert str(ideal.entr_mol_phase_comp(m, p, "foo")) == str(42)


def test_entr_mol_phase_invalid_phase(m):
    with pytest.raises(PropertyNotSupportedError):
        ideal.entr_mol_phase_comp(m, "foo", "bar")


def test_fug_phase_comp_liq(m):
    m._params.config.pressure_sat_comp = dummy_call

    for j in m._params.component_list:
        assert str(ideal.fug_phase_comp(m, "Liq", j)) == str(
            m.mole_frac_phase_comp["Liq", j]*42)


def test_fug_phase_comp_vap(m):
    for j in m._params.component_list:
        assert str(ideal.fug_phase_comp(m, "Vap", j)) == str(
            m.mole_frac_phase_comp["Vap", j]*m.pressure)


def test_fug_phase_comp_invalid_phase(m):
    with pytest.raises(PropertyNotSupportedError):
        ideal.fug_phase_comp(m, "foo", "bar")


def test_fug_coeff_phase_comp(m):
    for p in m._params.phase_list:
        for j in m._params.component_list:
            assert ideal.fug_coeff_phase_comp(m, p, j) == 1


def test_fug_coeff_phase_comp_invalid_phase(m):
    with pytest.raises(PropertyNotSupportedError):
        ideal.fug_coeff_phase_comp(m, "foo", "bar")


def test_gibbs_mol_phase(m):
    m.gibbs_mol_phase_comp = Var(m._params.phase_list,
                                 m._params.component_list)

    for p in m._params.phase_list:
        assert str(ideal.gibbs_mol_phase(m, p)) == str(
            sum(m.mole_frac_phase_comp[p, j]*m.gibbs_mol_phase_comp[p, j]
                for j in m._params.component_list))


def test_gibbs_mol_phase_comp(m):
    m.enth_mol_phase_comp = Var(m._params.phase_list,
                                m._params.component_list)
    m.entr_mol_phase_comp = Var(m._params.phase_list,
                                m._params.component_list)

    for p in m._params.phase_list:
        for j in m._params.component_list:
            assert str(ideal.gibbs_mol_phase_comp(m, p, j)) == str(
                    m.enth_mol_phase_comp[p, j] -
                    m.entr_mol_phase_comp[p, j] *
                    m.temperature)
