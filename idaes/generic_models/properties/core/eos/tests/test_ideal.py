##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
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

from pyomo.environ import ConcreteModel, log, Var, units as pyunits

from idaes.core import (declare_process_block_class,
                        LiquidPhase, VaporPhase, SolidPhase)
from idaes.generic_models.properties.core.eos.ideal import Ideal
from idaes.generic_models.properties.core.generic.generic_property import (
        GenericParameterData)
from idaes.core.util.exceptions import PropertyNotSupportedError
from idaes.core.util.constants import Constants as const


# Dummy method for property method calls
def dummy_call(b, j, T):
    return 42


# Dummy method to avoid errors when setting metadata dict
def set_metadata(b):
    pass


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
                "components": {"a": {}, "b": {}, "c": {}},
                "phases": {
                    "Vap": {"type": VaporPhase,
                            "equation_of_state": Ideal},
                    "Liq": {"type": LiquidPhase,
                            "equation_of_state": Ideal}},
                "base_units": {"time": pyunits.s,
                               "length": pyunits.m,
                               "mass": pyunits.kg,
                               "amount": pyunits.mol,
                               "temperature": pyunits.K},
                "state_definition": modules[__name__],
                "pressure_ref": 1e5,
                "temperature_ref": 300})

    m.props = m.params.state_block_class([1],
                                         default={"defined_state": False,
                                         "parameters": m.params})

    # Add common variables
    m.props[1].pressure = Var(initialize=101325)
    m.props[1].temperature = Var(initialize=300)
    m.props[1]._teq = Var([("Vap", "Liq")], initialize=300)
    m.props[1].mole_frac_phase_comp = Var(m.params.phase_list,
                                          m.params.component_list,
                                          initialize=0.5)

    return m


@pytest.fixture()
def m_sol():
    m = ConcreteModel()

    # Dummy params block with a Solid Phase to check phase typing
    m.params = DummyParameterBlock(default={
                "components": {"a": {}, "b": {}, "c": {}},
                "phases": {
                    "Sol": {"type": SolidPhase,
                            "equation_of_state": Ideal},
                    "Liq": {"type": LiquidPhase,
                            "equation_of_state": Ideal}},
                "base_units": {"time": pyunits.s,
                               "length": pyunits.m,
                               "mass": pyunits.kg,
                               "amount": pyunits.mol,
                               "temperature": pyunits.K},
                "state_definition": modules[__name__],
                "pressure_ref": 1e5,
                "temperature_ref": 300})

    m.props = m.params.build_state_block([1],
                                         default={"defined_state": False})

    # Add common variables
    m.props[1].pressure = Var(initialize=101325)
    m.props[1].temperature = Var(initialize=300)
    m.props[1]._teq = Var([("Vap", "Liq")], initialize=300)
    m.props[1].mole_frac_phase_comp = Var(m.params.phase_list,
                                          m.params.component_list,
                                          initialize=0.5)

    return m


@pytest.mark.unit
def test_common(m):
    assert Ideal.common(m.props, "foo") is None


@pytest.mark.unit
def test_compress_fact_phase_Liq(m):
    assert Ideal.compress_fact_phase(m.props[1], "Liq") == 1


@pytest.mark.unit
def test_compress_fact_phase_Vap(m):
    assert Ideal.compress_fact_phase(m.props[1], "Vap") == 1


@pytest.mark.unit
def test_compress_fact_phase_invalid_phase(m_sol):
    with pytest.raises(PropertyNotSupportedError):
        Ideal.compress_fact_phase(m_sol.props[1], "Sol")


@pytest.mark.unit
def test_dens_mass_phase(m):
    m.props[1].dens_mol_phase = Var(m.params.phase_list)
    m.props[1].mw_phase = Var(m.params.phase_list)

    for p in m.params.phase_list:
        assert str(Ideal.dens_mass_phase(m.props[1], p)) == str(
                m.props[1].dens_mol_phase[p]*m.props[1].mw_phase[p])


@pytest.mark.unit
def test_dens_mol_phase_liq(m):
    for j in m.params.component_list:
        m.params.get_component(j).config.dens_mol_liq_comp = dummy_call

    assert str(Ideal.dens_mol_phase(m.props[1], "Liq")) == str(
        sum(m.props[1].mole_frac_phase_comp["Liq", j]*42
            for j in m.params.component_list))


@pytest.mark.unit
def test_dens_mol_phase_vap(m):
    assert str(Ideal.dens_mol_phase(m.props[1], "Vap")) == (
            'props[1].pressure/(kg * m ** 2 / K / mol / s ** 2/J / K / '
            'mol*(8.314462618*J/mol/K)*props[1].temperature)')


@pytest.mark.unit
def test_dens_mol_phase_invalid_phase(m_sol):
    with pytest.raises(PropertyNotSupportedError):
        Ideal.dens_mol_phase(m_sol.props[1], "Sol")


@pytest.mark.unit
def test_enth_mol_phase(m):
    m.props[1].enth_mol_phase_comp = Var(m.params.phase_list,
                                         m.params.component_list)

    for p in m.params.phase_list:
        assert str(Ideal.enth_mol_phase(m.props[1], p)) == str(
            sum(m.props[1].mole_frac_phase_comp[p, j] *
                m.props[1].enth_mol_phase_comp[p, j]
                for j in m.params.component_list))


@pytest.mark.unit
def test_enth_mol_phase_comp(m):
    for j in m.params.component_list:
        m.params.get_component(j).config.enth_mol_liq_comp = dummy_call
        m.params.get_component(j).config.enth_mol_ig_comp = dummy_call

        assert str(Ideal.enth_mol_phase_comp(m.props[1], "Liq", j)) == str(42)
        assert str(Ideal.enth_mol_phase_comp(m.props[1], "Vap", j)) == str(42)


@pytest.mark.unit
def test_enth_mol_phase_invalid_phase(m_sol):
    with pytest.raises(PropertyNotSupportedError):
        Ideal.enth_mol_phase_comp(m_sol.props[1], "Sol", "foo")


@pytest.mark.unit
def test_entr_mol_phase(m):
    m.props[1].entr_mol_phase_comp = Var(m.params.phase_list,
                                         m.params.component_list)

    for p in m.params.phase_list:
        assert str(Ideal.entr_mol_phase(m.props[1], p)) == str(
            sum(m.props[1].mole_frac_phase_comp[p, j] *
                m.props[1].entr_mol_phase_comp[p, j]
                for j in m.params.component_list))


@pytest.mark.unit
def test_entr_mol_phase_comp(m):
    for j in m.params.component_list:
        m.params.get_component(j).config.entr_mol_liq_comp = dummy_call
        m.params.get_component(j).config.entr_mol_ig_comp = dummy_call

        assert str(Ideal.entr_mol_phase_comp(m.props[1], "Liq", j)) == str(42)
        assert str(Ideal.entr_mol_phase_comp(m.props[1], "Vap", j)) == (
            '42 - kg * m ** 2 / K / mol / s ** 2/J / K / '
            'mol*(8.314462618*J/mol/K)*log(props[1].mole_frac_phase_comp'
            '[Vap,{}]*props[1].pressure/params.pressure_ref)'.format(j))


@pytest.mark.unit
def test_entr_mol_phase_invalid_phase(m_sol):
    with pytest.raises(PropertyNotSupportedError):
        Ideal.entr_mol_phase_comp(m_sol.props[1], "Sol", "foo")


@pytest.mark.unit
def test_fug_phase_comp_liq(m):
    for j in m.params.component_list:
        m.params.get_component(j).config.pressure_sat_comp = dummy_call

        assert (str(Ideal.fug_phase_comp(m.props[1], "Liq", j)) ==
                str(m.props[1].mole_frac_phase_comp["Liq", j] * 42))


@pytest.mark.unit
def test_fug_phase_comp_vap(m):
    for j in m.params.component_list:
        assert (str(Ideal.fug_phase_comp(m.props[1], "Vap", j)) ==
                str(m.props[1].mole_frac_phase_comp["Vap", j] *
                    m.props[1].pressure))


@pytest.mark.unit
def test_fug_phase_comp_invalid_phase(m_sol):
    with pytest.raises(PropertyNotSupportedError):
        Ideal.fug_phase_comp(m_sol.props[1], "Sol", "foo")


@pytest.mark.unit
def test_fug_phase_comp_liq_eq(m):
    for j in m.params.component_list:
        m.params.get_component(j).config.pressure_sat_comp = dummy_call

        assert (str(Ideal.fug_phase_comp_eq(
                        m.props[1], "Liq", j, ("Vap", "Liq"))) ==
                str(m.props[1].mole_frac_phase_comp["Liq", j] * 42))


@pytest.mark.unit
def test_fug_phase_comp_vap_eq(m):
    for j in m.params.component_list:
        assert (str(Ideal.fug_phase_comp_eq(
                        m.props[1], "Vap", j, ("Vap", "Liq"))) ==
                str(m.props[1].mole_frac_phase_comp["Vap", j] *
                    m.props[1].pressure))


@pytest.mark.unit
def test_fug_phase_comp_invalid_phase_eq(m_sol):
    with pytest.raises(PropertyNotSupportedError):
        Ideal.fug_phase_comp_eq(m_sol.props[1], "Sol", "foo", ("Vap", "Liq"))


@pytest.mark.unit
def test_fug_coeff_phase_comp(m):
    for p in m.params.phase_list:
        for j in m.params.component_list:
            assert Ideal.fug_coeff_phase_comp(m.props[1], p, j) == 1


@pytest.mark.unit
def test_fug_coeff_phase_comp_invalid_phase(m_sol):
    with pytest.raises(PropertyNotSupportedError):
        Ideal.fug_coeff_phase_comp(m_sol.props[1], "Sol", "foo")


@pytest.mark.unit
def test_gibbs_mol_phase(m):
    m.props[1].gibbs_mol_phase_comp = Var(m.params.phase_list,
                                          m.params.component_list)

    for p in m.params.phase_list:
        assert str(Ideal.gibbs_mol_phase(m.props[1], p)) == str(
            sum(m.props[1].mole_frac_phase_comp[p, j] *
                m.props[1].gibbs_mol_phase_comp[p, j]
                for j in m.params.component_list))


@pytest.mark.unit
def test_gibbs_mol_phase_comp(m):
    m.props[1].enth_mol_phase_comp = Var(m.params.phase_list,
                                         m.params.component_list)
    m.props[1].entr_mol_phase_comp = Var(m.params.phase_list,
                                         m.params.component_list)

    for p in m.params.phase_list:
        for j in m.params.component_list:
            assert str(Ideal.gibbs_mol_phase_comp(m.props[1], p, j)) == str(
                    m.props[1].enth_mol_phase_comp[p, j] -
                    m.props[1].entr_mol_phase_comp[p, j] *
                    m.props[1].temperature)
