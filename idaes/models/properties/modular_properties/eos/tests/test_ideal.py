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
Tests for ideal equation of state methods

Author: Andrew Lee
"""
import pytest
from sys import modules

from pyomo.environ import ConcreteModel, Var, units as pyunits, value
from pyomo.util.check_units import assert_units_equivalent

from idaes.core import (
    declare_process_block_class,
    LiquidPhase,
    VaporPhase,
    SolidPhase,
    Solute,
    Solvent,
    Apparent,
)
from idaes.models.properties.modular_properties.eos.ideal import Ideal
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterData,
)
from idaes.core.util.exceptions import ConfigurationError, PropertyNotSupportedError
from idaes.core.util.constants import Constants as const


# Dummy method for property method calls
def dummy_call(b, j, T):
    return 42


def dummy_call2(b, j, T):
    return 7


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
    m.params = DummyParameterBlock(
        components={"a": {}, "b": {}, "c": {}},
        phases={
            "Vap": {"type": VaporPhase, "equation_of_state": Ideal},
            "Liq": {"type": LiquidPhase, "equation_of_state": Ideal},
        },
        base_units={
            "time": pyunits.s,
            "length": pyunits.m,
            "mass": pyunits.kg,
            "amount": pyunits.mol,
            "temperature": pyunits.K,
        },
        state_definition=modules[__name__],
        pressure_ref=100000.0,
        temperature_ref=300,
    )

    m.props = m.params.state_block_class([1], defined_state=False, parameters=m.params)

    # Add common variables
    m.props[1].pressure = Var(initialize=101325)
    m.props[1].temperature = Var(initialize=300, units=pyunits.K)
    m.props[1]._teq = Var([("Vap", "Liq")], initialize=300)
    m.props[1].mole_frac_phase_comp = Var(
        m.params.phase_list, m.params.component_list, initialize=0.5
    )

    return m


@pytest.fixture()
def m_sol():
    m = ConcreteModel()

    # Dummy params block with a Solid Phase to check phase typing
    m.params = DummyParameterBlock(
        components={"a": {}, "b": {}, "c": {}},
        phases={
            "Sol": {"type": SolidPhase, "equation_of_state": Ideal},
            "Liq": {"type": LiquidPhase, "equation_of_state": Ideal},
        },
        base_units={
            "time": pyunits.s,
            "length": pyunits.m,
            "mass": pyunits.kg,
            "amount": pyunits.mol,
            "temperature": pyunits.K,
        },
        state_definition=modules[__name__],
        pressure_ref=100000.0,
        temperature_ref=300,
    )

    m.props = m.params.build_state_block([1], defined_state=False)

    # Add common variables
    m.props[1].pressure = Var(initialize=101325)
    m.props[1].temperature = Var(initialize=300)
    m.props[1]._teq = Var([("Vap", "Liq")], initialize=300)
    m.props[1].mole_frac_phase_comp = Var(
        m.params.phase_list, m.params.component_list, initialize=0.5
    )

    return m


@pytest.mark.unit
def test_common(m):
    assert Ideal.common(m.props, "foo") is None


@pytest.mark.unit
def test_compress_fact_phase_Liq(m):
    assert Ideal.compress_fact_phase(m.props[1], "Liq") == 0


@pytest.mark.unit
def test_compress_fact_phase_Vap(m):
    assert Ideal.compress_fact_phase(m.props[1], "Vap") == 1


@pytest.mark.unit
def test_compress_fact_phase_sol(m_sol):
    assert Ideal.compress_fact_phase(m_sol.props[1], "Sol") == 0


@pytest.mark.unit
def test_cp_mol_phase(m):
    m.props[1].cp_mol_phase_comp = Var(m.params.phase_list, m.params.component_list)

    for p in m.params.phase_list:
        assert str(Ideal.cp_mol_phase(m.props[1], p)) == str(
            sum(
                m.props[1].mole_frac_phase_comp[p, j]
                * m.props[1].cp_mol_phase_comp[p, j]
                for j in m.params.component_list
            )
        )


@pytest.mark.unit
def test_cp_mol_phase_comp(m):
    for j in m.params.component_list:
        m.params.get_component(j).config.cp_mol_liq_comp = dummy_call
        m.params.get_component(j).config.cp_mol_ig_comp = dummy_call

        assert str(Ideal.cp_mol_phase_comp(m.props[1], "Liq", j)) == str(42)
        assert str(Ideal.cp_mol_phase_comp(m.props[1], "Vap", j)) == str(42)


@pytest.mark.unit
def test_cv_mol_phase(m):
    m.props[1].cv_mol_phase_comp = Var(m.params.phase_list, m.params.component_list)

    for p in m.params.phase_list:
        assert str(Ideal.cv_mol_phase(m.props[1], p)) == str(
            sum(
                m.props[1].mole_frac_phase_comp[p, j]
                * m.props[1].cv_mol_phase_comp[p, j]
                for j in m.params.component_list
            )
        )


@pytest.mark.unit
def test_cv_mol_phase_comp(m):
    for j in m.params.component_list:
        m.params.get_component(j).config.cp_mol_liq_comp = dummy_call
        m.params.get_component(j).config.cp_mol_ig_comp = dummy_call

        assert str(Ideal.cv_mol_phase_comp(m.props[1], "Liq", j)) == str(42)
        assert str(Ideal.cv_mol_phase_comp(m.props[1], "Vap", j)) == str(
            42
            - pyunits.convert(
                const.gas_constant,
                to_units=pyunits.kg
                * pyunits.m**2
                / pyunits.s**2
                / pyunits.mol
                / pyunits.K,
            )
        )


@pytest.mark.unit
def test_heat_capacity_ratio_phase(m):

    m.props[1].cp_mol_phase = Var(m.params.phase_list)
    m.props[1].cv_mol_phase = Var(m.params.phase_list)

    for p in m.params.phase_list:
        assert str(Ideal.heat_capacity_ratio_phase(m.props[1], p)) == (
            str(m.props[1].cp_mol_phase[p] / m.props[1].cv_mol_phase[p])
        )


@pytest.mark.unit
def test_dens_mass_phase(m):
    m.props[1].dens_mol_phase = Var(m.params.phase_list)
    m.props[1].mw_phase = Var(m.params.phase_list)

    for p in m.params.phase_list:
        assert str(Ideal.dens_mass_phase(m.props[1], p)) == str(
            m.props[1].dens_mol_phase[p] * m.props[1].mw_phase[p]
        )


@pytest.mark.unit
def test_dens_mol_phase_liq(m):
    for j in m.params.component_list:
        m.params.get_component(j).config.dens_mol_liq_comp = dummy_call

    assert str(Ideal.dens_mol_phase(m.props[1], "Liq")) == str(
        1
        / (
            1 / 42 * m.props[1].mole_frac_phase_comp["Liq", "a"]
            + 1 / 42 * m.props[1].mole_frac_phase_comp["Liq", "b"]
            + 1 / 42 * m.props[1].mole_frac_phase_comp["Liq", "c"]
        )
    )


@pytest.mark.unit
def test_dens_mol_phase_vap(m):
    assert str(Ideal.dens_mol_phase(m.props[1], "Vap")) == (
        "props[1].pressure/("
        + str(Ideal.gas_constant(m.props[1]))
        + "*props[1].temperature)"
    )


@pytest.mark.unit
def test_dens_mol_phase_sol(m_sol):
    for j in m_sol.params.component_list:
        m_sol.params.get_component(j).config.dens_mol_sol_comp = dummy_call
        m_sol.params.get_component(j).config.dens_mol_liq_comp = dummy_call

    assert str(Ideal.dens_mol_phase(m_sol.props[1], "Sol")) == str(
        1
        / sum(
            m_sol.props[1].mole_frac_phase_comp["Sol", j] / 42
            for j in m_sol.params.component_list
        )
    )


@pytest.mark.unit
def test_energy_internal_mol_phase(m):
    m.props[1].energy_internal_mol_phase_comp = Var(
        m.params.phase_list, m.params.component_list
    )

    for p in m.params.phase_list:
        assert str(Ideal.energy_internal_mol_phase(m.props[1], p)) == str(
            sum(
                m.props[1].mole_frac_phase_comp[p, j]
                * m.props[1].energy_internal_mol_phase_comp[p, j]
                for j in m.params.component_list
            )
        )


@pytest.mark.unit
def test_energy_internal_mol_phase_comp_no_h_form(m):
    for j in m.params.component_list:
        m.params.config.include_enthalpy_of_formation = False
        m.params.get_component(j).config.enth_mol_liq_comp = dummy_call
        m.params.get_component(j).config.enth_mol_ig_comp = dummy_call

        assert str(Ideal.energy_internal_mol_phase_comp(m.props[1], "Liq", j)) == str(
            42
        )
        assert str(Ideal.energy_internal_mol_phase_comp(m.props[1], "Vap", j)) == str(
            42
            - pyunits.convert(
                const.gas_constant,
                to_units=pyunits.kg
                * pyunits.m**2
                / pyunits.s**2
                / pyunits.mol
                / pyunits.K,
            )
            * (m.props[1].temperature - m.params.temperature_ref)
        )


@pytest.mark.unit
def test_energy_internal_mol_phase_comp_with_h_form(m):
    for j in m.params.component_list:
        m.params.config.include_enthalpy_of_formation = True
        m.params.get_component(j).config.elemental_composition = {
            "C": 1,
            "He": 2,
            "O": 4,
        }
        m.params.get_component(j).config.enth_mol_liq_comp = dummy_call
        m.params.get_component(j).config.enth_mol_ig_comp = dummy_call

        # For liquid phase, delta(n) should be 4 (2*He + 2*O2)
        assert (
            str(Ideal.energy_internal_mol_phase_comp(m.props[1], "Liq", j))
            == "42 -4.0*("
            + str(Ideal.gas_constant(m.props[1]))
            + ")*params.temperature_ref"
        )
        # For vapor phase, delta(n) should be 3 (2*He + 2*O2 - 1*component)
        assert str(Ideal.energy_internal_mol_phase_comp(m.props[1], "Vap", j)) == str(
            42
            - pyunits.convert(
                const.gas_constant,
                to_units=pyunits.kg
                * pyunits.m**2
                / pyunits.s**2
                / pyunits.mol
                / pyunits.K,
            )
            * (m.props[1].temperature - m.params.temperature_ref)
            - 3.0
            * pyunits.convert(
                const.gas_constant,
                to_units=pyunits.kg
                * pyunits.m**2
                / pyunits.s**2
                / pyunits.mol
                / pyunits.K,
            )
            * m.params.temperature_ref
        )


@pytest.mark.unit
def test_enth_mol_phase(m):
    for j in m.params.component_list:
        m.params.get_component(j).config.enth_mol_liq_comp = dummy_call
        m.params.get_component(j).config.enth_mol_ig_comp = dummy_call
        m.params.get_component(j).config.dens_mol_liq_comp = dummy_call2

    assert str(Ideal.enth_mol_phase(m.props[1], "Vap")) == (
        "props[1].mole_frac_phase_comp[Vap,a]*42.0 + "
        "props[1].mole_frac_phase_comp[Vap,b]*42.0 + "
        "props[1].mole_frac_phase_comp[Vap,c]*42.0"
    )
    assert str(Ideal.enth_mol_phase(m.props[1], "Liq")) == str(
        sum(
            m.props[1].mole_frac_phase_comp["Liq", j] * 42
            for j in m.params.component_list
        )
        + (m.props[1].pressure - m.params.pressure_ref)
        / m.props[1].dens_mol_phase["Liq"]
    )


@pytest.mark.unit
def test_enth_mol_phase_comp(m):
    for j in m.params.component_list:
        m.params.get_component(j).config.enth_mol_liq_comp = dummy_call
        m.params.get_component(j).config.enth_mol_ig_comp = dummy_call
        m.params.get_component(j).config.dens_mol_liq_comp = dummy_call2

    for j in m.params.component_list:
        assert str(Ideal.enth_mol_phase_comp(m.props[1], "Liq", j)) == str(
            42
            + (m.props[1].pressure - m.params.pressure_ref)
            / m.props[1].dens_mol_phase["Liq"]
        )
        assert str(Ideal.enth_mol_phase_comp(m.props[1], "Vap", j)) == str(42)


@pytest.mark.unit
def test_enth_mol_phase_sol(m_sol):
    for j in m_sol.params.component_list:
        m_sol.params.get_component(j).config.enth_mol_sol_comp = dummy_call
        m_sol.params.get_component(j).config.dens_mol_liq_comp = dummy_call2
        m_sol.params.get_component(j).config.dens_mol_sol_comp = dummy_call2

    for j in m_sol.params.component_list:
        assert str(Ideal.enth_mol_phase_comp(m_sol.props[1], "Sol", j)) == str(
            42
            + (m_sol.props[1].pressure - m_sol.params.pressure_ref)
            / m_sol.props[1].dens_mol_phase["Sol"]
        )


@pytest.mark.unit
def test_entr_mol_phase(m):
    m.props[1].entr_mol_phase_comp = Var(m.params.phase_list, m.params.component_list)

    for p in m.params.phase_list:
        assert str(Ideal.entr_mol_phase(m.props[1], p)) == str(
            sum(
                m.props[1].mole_frac_phase_comp[p, j]
                * m.props[1].entr_mol_phase_comp[p, j]
                for j in m.params.component_list
            )
        )


@pytest.mark.unit
def test_entr_mol_phase_comp(m):
    for j in m.params.component_list:
        m.params.get_component(j).config.entr_mol_liq_comp = dummy_call
        m.params.get_component(j).config.entr_mol_ig_comp = dummy_call

        assert str(Ideal.entr_mol_phase_comp(m.props[1], "Liq", j)) == str(42)
        assert str(Ideal.entr_mol_phase_comp(m.props[1], "Vap", j)) == (
            "42 - "
            + str(Ideal.gas_constant(m.props[1]))
            + "*log(props[1].mole_frac_phase_comp"
            "[Vap,{}]*props[1].pressure/params.pressure_ref)".format(j)
        )


@pytest.mark.unit
def test_entr_mol_phase_sol(m_sol):
    for j in m_sol.params.component_list:
        m_sol.params.get_component(j).config.entr_mol_sol_comp = dummy_call

        assert str(Ideal.entr_mol_phase_comp(m_sol.props[1], "Sol", j)) == str(42)


@pytest.mark.unit
def test_fug_phase_comp_liq(m):
    for j in m.params.component_list:
        m.params.get_component(j).config.pressure_sat_comp = dummy_call

        assert str(Ideal.fug_phase_comp(m.props[1], "Liq", j)) == str(
            m.props[1].mole_frac_phase_comp["Liq", j] * 42
        )


@pytest.mark.unit
def test_fug_phase_comp_vap(m):
    for j in m.params.component_list:
        assert str(Ideal.fug_phase_comp(m.props[1], "Vap", j)) == str(
            m.props[1].mole_frac_phase_comp["Vap", j] * m.props[1].pressure
        )


@pytest.mark.unit
def test_fug_phase_comp_invalid_phase(m_sol):
    with pytest.raises(PropertyNotSupportedError):
        Ideal.fug_phase_comp(m_sol.props[1], "Sol", "a")


@pytest.mark.unit
def test_fug_phase_comp_liq_eq(m):
    for j in m.params.component_list:
        m.params.get_component(j).config.pressure_sat_comp = dummy_call

        assert str(
            Ideal.fug_phase_comp_eq(m.props[1], "Liq", j, ("Vap", "Liq"))
        ) == str(m.props[1].mole_frac_phase_comp["Liq", j] * 42)


@pytest.mark.unit
def test_fug_phase_comp_vap_eq(m):
    for j in m.params.component_list:
        assert str(
            Ideal.fug_phase_comp_eq(m.props[1], "Vap", j, ("Vap", "Liq"))
        ) == str(m.props[1].mole_frac_phase_comp["Vap", j] * m.props[1].pressure)


@pytest.mark.unit
def test_fug_phase_comp_invalid_phase_eq(m_sol):
    with pytest.raises(PropertyNotSupportedError):
        Ideal.fug_phase_comp_eq(m_sol.props[1], "Sol", "a", ("Vap", "Liq"))


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
    m.props[1].gibbs_mol_phase_comp = Var(m.params.phase_list, m.params.component_list)

    for p in m.params.phase_list:
        assert str(Ideal.gibbs_mol_phase(m.props[1], p)) == str(
            sum(
                m.props[1].mole_frac_phase_comp[p, j]
                * m.props[1].gibbs_mol_phase_comp[p, j]
                for j in m.params.component_list
            )
        )


@pytest.mark.unit
def test_gibbs_mol_phase_comp(m):
    m.props[1].enth_mol_phase_comp = Var(m.params.phase_list, m.params.component_list)
    m.props[1].entr_mol_phase_comp = Var(m.params.phase_list, m.params.component_list)

    for p in m.params.phase_list:
        for j in m.params.component_list:
            assert str(Ideal.gibbs_mol_phase_comp(m.props[1], p, j)) == str(
                m.props[1].enth_mol_phase_comp[p, j]
                - m.props[1].entr_mol_phase_comp[p, j] * m.props[1].temperature
            )


@pytest.mark.unit
def test_pressure_osm_phase_no_solvent(m):
    with pytest.raises(
        ConfigurationError,
        match="called for pressure_osm, but no solvents were "
        "defined. Osmotic pressure requires at least one component to be "
        "declared as a solvent.",
    ):
        m.props[1].pressure_osm_phase.display()


@pytest.mark.unit
def test_pressure_osm_phase(m):
    # Create a dummy solvent_set
    m.params.solvent_set = ["a"]
    m.props[1].dens_mol_phase = Var(
        m.params.phase_list, initialize=55e3, units=pyunits.mol / pyunits.m**3
    )

    assert_units_equivalent(m.props[1].pressure_osm_phase["Liq"], pyunits.Pa)
    assert pytest.approx(value(const.gas_constant * 300 * 55e3), rel=1e-6) == value(
        m.props[1].pressure_osm_phase["Liq"]
    )
    assert str(m.props[1].pressure_osm_phase["Liq"]._expr) == (
        str(Ideal.gas_constant(m.props[1]))
        + "*"
        + str(
            m.props[1].temperature
            * (
                m.props[1].conc_mol_phase_comp["Liq", "b"]
                + m.props[1].conc_mol_phase_comp["Liq", "c"]
            )
        )
    )
    assert len(m.props[1].pressure_osm_phase) == 1


@pytest.mark.unit
def test_pressure_osm_phase_w_apparent_component():
    m = ConcreteModel()

    # Dummy params block
    m.params = DummyParameterBlock(
        components={
            "a": {"type": Solvent},
            "b": {"type": Solute},
            "c": {"type": Apparent, "dissociation_species": {"foo": 2}},
        },
        phases={
            "Vap": {"type": VaporPhase, "equation_of_state": Ideal},
            "Liq": {"type": LiquidPhase, "equation_of_state": Ideal},
        },
        base_units={
            "time": pyunits.s,
            "length": pyunits.m,
            "mass": pyunits.kg,
            "amount": pyunits.mol,
            "temperature": pyunits.K,
        },
        state_definition=modules[__name__],
        pressure_ref=100000.0,
        temperature_ref=300,
    )

    m.props = m.params.state_block_class([1], defined_state=False, parameters=m.params)

    # Add common variables
    m.props[1].pressure = Var(initialize=101325)
    m.props[1].temperature = Var(initialize=300, units=pyunits.K)
    m.props[1]._teq = Var([("Vap", "Liq")], initialize=300)
    m.props[1].mole_frac_phase_comp = Var(
        m.params.phase_list, m.params.component_list, initialize=0.5
    )

    m.props[1].dens_mol_phase = Var(
        m.params.phase_list, initialize=55e3, units=pyunits.mol / pyunits.m**3
    )

    assert_units_equivalent(m.props[1].pressure_osm_phase["Liq"], pyunits.Pa)
    # Factor of 1.5 accounts for the fact the c is doubled.
    assert pytest.approx(
        value(const.gas_constant * 300 * 1.5 * 55e3), rel=1e-6
    ) == value(m.props[1].pressure_osm_phase["Liq"])
    assert str(m.props[1].pressure_osm_phase["Liq"]._expr) == (
        str(Ideal.gas_constant(m.props[1]))
        + "*"
        + str(
            m.props[1].temperature
            * (
                m.props[1].conc_mol_phase_comp["Liq", "b"]
                + 2 * m.props[1].conc_mol_phase_comp["Liq", "c"]
            )
        )
    )
    assert len(m.props[1].pressure_osm_phase) == 1


@pytest.mark.unit
def test_vol_mol_phase_no_methods(m):
    with pytest.raises(
        ConfigurationError,
        match="does not have a method defined to use when "
        "calculating molar volume and density for component a "
        "in phase liq. Each component must define a method for "
        "either vol_mol_liq_comp or dens_mol_liq_comp.",
    ):
        Ideal.vol_mol_phase(m.props[1], "Liq")


@pytest.mark.unit
def test_vol_mol_phase():
    m = ConcreteModel()

    # Dummy params block
    m.params = DummyParameterBlock(
        components={
            "a": {"dens_mol_liq_comp": dummy_call},
            "b": {"dens_mol_liq_comp": dummy_call},
            "c": {"vol_mol_liq_comp": dummy_call},
        },
        phases={
            "Vap": {"type": VaporPhase, "equation_of_state": Ideal},
            "Liq": {"type": LiquidPhase, "equation_of_state": Ideal},
        },
        base_units={
            "time": pyunits.s,
            "length": pyunits.m,
            "mass": pyunits.kg,
            "amount": pyunits.mol,
            "temperature": pyunits.K,
        },
        state_definition=modules[__name__],
        pressure_ref=100000.0,
        temperature_ref=300,
    )

    m.props = m.params.state_block_class([1], defined_state=False, parameters=m.params)

    # Add common variables
    m.props[1].pressure = Var(initialize=101325)
    m.props[1].temperature = Var(initialize=300, units=pyunits.K)
    m.props[1]._teq = Var([("Vap", "Liq")], initialize=300)
    m.props[1].mole_frac_phase_comp = Var(
        m.params.phase_list, m.params.component_list, initialize=0.5
    )

    for p in m.params.phase_list:
        if p == "Vap":
            assert str(Ideal.vol_mol_phase(m.props[1], p)) == (
                str(Ideal.gas_constant(m.props[1]))
                + "*props[1].temperature/props[1].pressure"
            )
        else:
            assert str(Ideal.vol_mol_phase(m.props[1], p)) == str(
                1 / 42 * m.props[1].mole_frac_phase_comp["Liq", "a"]
                + 1 / 42 * m.props[1].mole_frac_phase_comp["Liq", "b"]
                + 42 * m.props[1].mole_frac_phase_comp["Liq", "c"]
            )
