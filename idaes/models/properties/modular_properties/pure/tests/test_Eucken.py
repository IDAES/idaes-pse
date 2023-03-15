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
Tests for Eucken equation for thermal conductivity from
Properties of Gases and Liquids, 5th Ed., 10-3-1

Authors: Douglas Allan
"""

import pytest
import types

from pyomo.environ import ConcreteModel, Block, value, Var, units as pyunits
from pyomo.common.config import ConfigBlock
from pyomo.util.check_units import assert_units_equivalent

from idaes.models.properties.modular_properties.pure.Eucken import Eucken
from idaes.core.util.misc import add_object_reference
from idaes.core.base.property_meta import PropertyClassMetadata
from idaes.core.util.exceptions import ConfigurationError


def construct_dummy_model(param_dict):
    m = ConcreteModel()

    # Create a dummy parameter block
    m.params = Block()

    m.params.config = ConfigBlock(implicit=True)
    # Parameters for sulfur dioxide, from Properties of Gases and Liquids 5th ed. Appendix B.
    m.params.config.parameter_data = {
        "f_int_eucken": (1.61, pyunits.dimensionless),
    }

    # Also need to dummy configblock on the model for the test
    m.config = ConfigBlock(implicit=True)

    m.meta_object = PropertyClassMetadata()
    m.meta_object._default_units.set_units(
        temperature=pyunits.K,
        mass=pyunits.kg,
        length=pyunits.m,
        time=pyunits.s,
        amount=pyunits.mol,
    )

    def get_metadata(self):
        return m.meta_object

    m.get_metadata = types.MethodType(get_metadata, m)
    m.params.get_metadata = types.MethodType(get_metadata, m.params)

    # Create variables that should exist on param block
    m.params.mw = Var(initialize=param_dict["mw"], units=pyunits.kg / pyunits.mol)
    # Hacking cp on the fake component object
    m.params.cp_mol_pure = Var(initialize=1, units=pyunits.J / pyunits.mol / pyunits.K)
    m.params.config.cp_mol_ig_comp = Block()

    def return_expression(b, cobj, T):
        return cobj.cp_mol_pure

    m.params.config.cp_mol_ig_comp.return_expression = return_expression

    # Create a dummy state block
    m.props = Block([1])
    add_object_reference(m.props[1], "params", m.params)

    m.props[1].temperature = Var(initialize=298, units=pyunits.K)
    # Have to trick the block into considering "params" as the component name
    m.props[1]._visc_d_phase_comp = Var(
        ["Vap"], ["params"], initialize=1e-7, units=pyunits.Pa * pyunits.s
    )
    return m


@pytest.mark.unit
def test_therm_cond_phase_comp_acetone():
    m = construct_dummy_model({"mw": 58.080e-3})

    Eucken.therm_cond_phase_comp.build_parameters(m.params, "Vap")

    assert isinstance(m.params.f_int_eucken, Var)
    assert value(m.params.f_int_eucken) == pytest.approx(1.61, rel=1e-12)

    expr = Eucken.therm_cond_phase_comp.return_expression(
        m.props[1], m.params, "Vap", m.props[1].temperature
    )

    # Pulled from Table 10-2, Properties of Gases and Liquids 5th Ed.
    T_list = [353, 393, 457]  # Gas temperature in K
    therm_cond_list = [
        15.7e-3,
        19.4e-3,
        24.7e-3,
    ]  # Experimental thermal conductivity in W/(m K)
    visc_list = [90.0e-7, 100e-7, 114.5e-7]  # Viscosities in Pascal seconds
    Cv_list = [77.9, 84.2, 96.1]  # Constant *volume* heat capacity in J/(mol K)
    f_int_list = [1, 1.32, 1.15]
    err_list = [
        [-4.7, -8.7, -8.4],  # Eucken percent error
        [16, 12, 13],  # Modified Eucken percent error
        [5.1, 0.9, 1.7],
    ]  # Stiel and Thodos percent error
    for i in range(3):
        m.props[1].temperature.value = T_list[i]
        m.props[1]._visc_d_phase_comp["Vap", "params"].value = visc_list[i]
        m.params.cp_mol_pure.value = Cv_list[i] + 8.314
        for j in range(3):
            m.params.f_int_eucken.value = f_int_list[j]
            err = 100 * value(expr / therm_cond_list[i] - 1)
            # Obtain percent error, then compare to that tabulated
            # in Table 10-2
            assert err == pytest.approx(err_list[j][i], abs=0.3)

    assert_units_equivalent(expr, pyunits.W / pyunits.m / pyunits.K)


@pytest.mark.unit
def test_therm_cond_phase_comp_CO2():
    m = construct_dummy_model({"mw": 44.009e-3})

    Eucken.therm_cond_phase_comp.build_parameters(m.params, "Vap")

    assert isinstance(m.params.f_int_eucken, Var)
    assert value(m.params.f_int_eucken) == pytest.approx(1.61, rel=1e-12)

    expr = Eucken.therm_cond_phase_comp.return_expression(
        m.props[1], m.params, "Vap", m.props[1].temperature
    )

    # Pulled from Table 10-2, Properties of Gases and Liquids 5th Ed.
    T_list = [200, 300, 473, 598, 1273]  # Gas temperature in K
    therm_cond_list = [
        9.51e-3,
        16.7e-3,
        28.4e-3,
        37.9e-3,
        81.7e-3,
    ]  # Experimental thermal conductivity in W/(m K)
    visc_list = [
        101.5e-7,
        149.5e-7,
        225.0e-7,
        272.8e-7,
        465.1e-7,
    ]  # Viscosities in Pascal seconds
    Cv_list = [
        24.7,
        28.9,
        35.6,
        39.5,
        48.8,
    ]  # Constant *volume* heat capacity in J/(mol K)
    f_int_list = [1, 1.32, 1.15]
    err_list = [
        [5.3, -3.2, -2.2, -4.8, -13],  # Eucken percent error
        [15, 7.5, 11, 9.3, 2.4],  # Modified Eucken percent error
        [9.8, 1.9, 4.1, 1.9, -5.6],
    ]  # Stiel and Thodos percent error
    for i in range(5):
        m.props[1].temperature.value = T_list[i]
        m.props[1]._visc_d_phase_comp["Vap", "params"].value = visc_list[i]
        m.params.cp_mol_pure.value = Cv_list[i] + 8.314
        for j in range(3):
            m.params.f_int_eucken.value = f_int_list[j]
            err = 100 * value(expr / therm_cond_list[i] - 1)
            # Obtain percent error, then compare to that tabulated
            # in Table 10-2
            assert err == pytest.approx(err_list[j][i], abs=0.5)

    assert_units_equivalent(expr, pyunits.W / pyunits.m / pyunits.K)


@pytest.mark.unit
def test_therm_cond_phase_comp_ammonia():
    m = construct_dummy_model({"mw": 17.031e-3})

    Eucken.therm_cond_phase_comp.build_parameters(m.params, "Vap")

    assert isinstance(m.params.f_int_eucken, Var)
    assert value(m.params.f_int_eucken) == pytest.approx(1.61, rel=1e-12)

    expr = Eucken.therm_cond_phase_comp.return_expression(
        m.props[1], m.params, "Vap", m.props[1].temperature
    )

    # Pulled from Table 10-2, Properties of Gases and Liquids 5th Ed.
    T_list = [213, 273]  # Gas temperature in K
    therm_cond_list = [16.5e-3, 21.9e-3]  # Experimental thermal conductivity in W/(m K)
    visc_list = [73.2e-7, 90.6e-7]  # Viscosities in Pascal seconds
    Cv_list = [25.4, 26.7]  # Constant *volume* heat capacity in J/(mol K)
    f_int_list = [1, 1.32, 1.15]
    err_list = [
        [15, 10],  # Eucken percent error
        [26, 21],  # Modified Eucken percent error
        [20, 16],
    ]  # Stiel and Thodos percent error
    for i in range(2):
        m.props[1].temperature.value = T_list[i]
        m.props[1]._visc_d_phase_comp["Vap", "params"].value = visc_list[i]
        m.params.cp_mol_pure.value = Cv_list[i] + 8.314
        for j in range(3):
            m.params.f_int_eucken.value = f_int_list[j]
            err = 100 * value(expr / therm_cond_list[i] - 1)
            # Obtain percent error, then compare to that tabulated
            # in Table 10-2
            assert err == pytest.approx(err_list[j][i], abs=0.6)

    assert_units_equivalent(expr, pyunits.W / pyunits.m / pyunits.K)


@pytest.mark.unit
def test_no_cp_ig_error():
    m = ConcreteModel()

    # Create a dummy parameter block
    m.params = Block()

    m.params.config = ConfigBlock(implicit=True)
    # Parameters for sulfur dioxide, from Properties of Gases and Liquids 5th ed. Appendix B.
    m.params.config.parameter_data = {
        "f_int_eucken": (1.61, pyunits.dimensionless),
    }

    # Also need to dummy configblock on the model for the test
    m.config = ConfigBlock(implicit=True)

    m.meta_object = PropertyClassMetadata()
    m.meta_object._default_units.set_units(
        temperature=pyunits.K,
        mass=pyunits.kg,
        length=pyunits.m,
        time=pyunits.s,
        amount=pyunits.mol,
    )

    def get_metadata(self):
        return m.meta_object

    m.get_metadata = types.MethodType(get_metadata, m)
    m.params.get_metadata = types.MethodType(get_metadata, m.params)

    # Create variables that should exist on param block
    m.params.mw = Var(initialize=2.0, units=pyunits.kg / pyunits.mol)
    # Hacking cp on the fake component object

    # Create a dummy state block
    m.props = Block([1])
    add_object_reference(m.props[1], "params", m.params)

    m.props[1].temperature = Var(initialize=298, units=pyunits.K)
    # Have to trick the block into considering "params" as the component name
    m.props[1].visc_d_phase_comp = Var(
        ["Vap"], ["params"], initialize=1e-7, units=pyunits.Pa * pyunits.s
    )

    Eucken.therm_cond_phase_comp.build_parameters(m.params, "Vap")

    with pytest.raises(
        ConfigurationError,
        match="Cannot find method to calculate cp_mol_ig_comp for component params.",
    ):
        expr = Eucken.therm_cond_phase_comp.return_expression(
            m.props[1], m.params, "Vap", m.props[1].temperature
        )
