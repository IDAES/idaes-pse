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
Tests for constant pure component properties

Liquid properties are tested with with data for water
Ideal gas properties are tested with data for Air

All parameter indicies based on conventions used by the source

Authors: Andres J Calderon, Andrew Lee
"""

import pytest
import types

from pyomo.environ import ConcreteModel, Block, value, Var, units as pyunits
from pyomo.common.config import ConfigBlock
from pyomo.util.check_units import assert_units_equivalent

from idaes.models.properties.modular_properties.pure.ConstantProperties import *
from idaes.core.util.misc import add_object_reference
from idaes.core.base.property_meta import PropertyClassMetadata, UnitSet


@pytest.fixture()
def frame():
    m = ConcreteModel()

    # Create a dummy parameter block
    m.params = Block()

    m.params.config = ConfigBlock(implicit=True)
    m.params.config.parameter_data = {
        "dens_mol_liq_comp_coeff": (55.2046332734557, pyunits.kmol / pyunits.m**3),
        "cp_mol_liq_comp_coeff": (
            75704.2953333398,
            pyunits.J / pyunits.kmol / pyunits.K,
        ),
        "enth_mol_form_liq_comp_ref": (-285.83, pyunits.kJ / pyunits.mol),
        "entr_mol_form_liq_comp_ref": (69.95, pyunits.J / pyunits.K / pyunits.mol),
        "cp_mol_ig_comp_coeff": (29114.850, pyunits.J / pyunits.kmol / pyunits.K),
        "enth_mol_form_ig_comp_ref": (0.0, pyunits.kJ / pyunits.mol),
        "entr_mol_form_ig_comp_ref": (0.0, pyunits.J / pyunits.K / pyunits.mol),
        "dens_mol_sol_comp_coeff": (100, pyunits.kmol / pyunits.m**3),
        "cp_mol_sol_comp_coeff": (100000, pyunits.J / pyunits.kmol / pyunits.K),
        "enth_mol_form_sol_comp_ref": (-300, pyunits.kJ / pyunits.mol),
        "entr_mol_form_sol_comp_ref": (50, pyunits.J / pyunits.K / pyunits.mol),
    }
    m.params.config.include_enthalpy_of_formation = True

    # Also need to dummy configblock on the model for the test
    m.config = ConfigBlock(implicit=True)
    m.config.include_enthalpy_of_formation = True

    m.meta_object = PropertyClassMetadata()
    m.meta_object._default_units = UnitSet(
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

    # Add necessary parameters to parameter block
    m.params.temperature_ref = Var(initialize=300, units=pyunits.K)

    # Create a dummy state block
    m.props = Block([1])
    add_object_reference(m.props[1], "params", m.params)

    m.props[1].temperature = Var(initialize=300, units=pyunits.K)

    return m


@pytest.mark.unit
def test_cp_mol_liq_comp(frame):
    Constant.cp_mol_liq_comp.build_parameters(frame.params)

    assert isinstance(frame.params.cp_mol_liq_comp_coeff, Var)
    assert value(frame.params.cp_mol_liq_comp_coeff) == pytest.approx(
        75704.3 / 1000, rel=1e-5
    )

    expr = Constant.cp_mol_liq_comp.return_expression(
        frame.props[1], frame.params, frame.props[1].temperature
    )
    assert value(expr) == pytest.approx(75704.3 / 1000, rel=1e-5)

    frame.props[1].temperature.value = 400
    assert value(expr) == pytest.approx(75704.3 / 1000, rel=1e-5)

    assert_units_equivalent(expr, pyunits.J / pyunits.mol / pyunits.K)


@pytest.mark.unit
def test_enth_mol_liq_comp(frame):

    Constant.enth_mol_liq_comp.build_parameters(frame.params)

    assert isinstance(frame.params.enth_mol_form_liq_comp_ref, Var)
    assert value(frame.params.enth_mol_form_liq_comp_ref) == -285.83e3

    expr = Constant.enth_mol_liq_comp.return_expression(
        frame.props[1], frame.params, frame.props[1].temperature
    )
    assert value(expr) == value(frame.params.enth_mol_form_liq_comp_ref)

    frame.props[1].temperature.value = 301
    assert value(expr) == pytest.approx(-285754.29570466664, rel=1e-3)

    frame.props[1].temperature.value = 400
    assert value(expr) == pytest.approx(-278259.570466, rel=1e-3)

    assert_units_equivalent(expr, pyunits.J / pyunits.mol)


@pytest.mark.unit
def test_enth_mol_liq_comp_no_formation(frame):
    frame.config.include_enthalpy_of_formation = False
    frame.params.config.include_enthalpy_of_formation = False

    Constant.enth_mol_liq_comp.build_parameters(frame.params)

    assert not hasattr(frame.params, "enth_mol_form_liq_comp_ref")

    expr = Constant.enth_mol_liq_comp.return_expression(
        frame.props[1], frame.params, frame.props[1].temperature
    )
    assert value(expr) == 0

    frame.props[1].temperature.value = 400
    assert value(expr) == pytest.approx(7570.42, rel=1e-3)

    assert_units_equivalent(expr, pyunits.J / pyunits.mol)


@pytest.mark.unit
def test_entr_mol_liq_comp(frame):
    Constant.entr_mol_liq_comp.build_parameters(frame.params)

    assert isinstance(frame.params.entr_mol_form_liq_comp_ref, Var)
    assert value(frame.params.entr_mol_form_liq_comp_ref) == 69.95

    expr = Constant.entr_mol_liq_comp.return_expression(
        frame.props[1], frame.params, frame.props[1].temperature
    )
    assert value(expr) == value(frame.params.entr_mol_form_liq_comp_ref)

    frame.props[1].temperature.value = 400
    assert value(expr) == pytest.approx(91.72, rel=1e-3)

    assert_units_equivalent(expr, pyunits.J / pyunits.mol / pyunits.K)


@pytest.mark.unit
def test_dens_mol_liq_comp(frame):
    Constant.dens_mol_liq_comp.build_parameters(frame.params)

    assert isinstance(frame.params.dens_mol_liq_comp_coeff, Var)
    assert value(frame.params.dens_mol_liq_comp_coeff) == pytest.approx(
        55.204e3, rel=1e-4
    )

    expr = Constant.dens_mol_liq_comp.return_expression(
        frame.props[1], frame.params, frame.props[1].temperature
    )
    assert value(expr) == pytest.approx(55.204e3, rel=1e-4)

    frame.props[1].temperature.value = 400
    assert value(expr) == pytest.approx(55.204e3, rel=1e-4)

    assert_units_equivalent(expr, pyunits.mol / pyunits.m**3)


@pytest.mark.unit
def test_cp_mol_ig_comp(frame):
    Constant.cp_mol_ig_comp.build_parameters(frame.params)

    assert isinstance(frame.params.cp_mol_ig_comp_coeff, Var)
    assert value(frame.params.cp_mol_ig_comp_coeff) == pytest.approx(
        29114.850 / 1000, rel=1e-5
    )

    expr = Constant.cp_mol_ig_comp.return_expression(
        frame.props[1], frame.params, frame.props[1].temperature
    )
    assert value(expr) == pytest.approx(29114.850 / 1000, rel=1e-5)

    frame.props[1].temperature.value = 400
    assert value(expr) == pytest.approx(29114.850 / 1000, rel=1e-5)

    assert_units_equivalent(expr, pyunits.J / pyunits.mol / pyunits.K)


@pytest.mark.unit
def test_enth_mol_ig_comp(frame):

    Constant.enth_mol_ig_comp.build_parameters(frame.params)

    assert isinstance(frame.params.enth_mol_form_ig_comp_ref, Var)
    assert value(frame.params.enth_mol_form_ig_comp_ref) == -0.0

    expr = Constant.enth_mol_ig_comp.return_expression(
        frame.props[1], frame.params, frame.props[1].temperature
    )
    assert value(expr) == value(frame.params.enth_mol_form_ig_comp_ref)

    frame.props[1].temperature.value = 301
    assert value(expr) == pytest.approx(29.114850, rel=1e-3)

    frame.props[1].temperature.value = 400
    assert value(expr) == pytest.approx(2911.485, rel=1e-3)

    assert_units_equivalent(expr, pyunits.J / pyunits.mol)


@pytest.mark.unit
def test_enth_mol_ig_comp_no_formation(frame):

    frame.config.include_enthalpy_of_formation = False
    frame.params.config.include_enthalpy_of_formation = False

    Constant.enth_mol_ig_comp.build_parameters(frame.params)

    assert not hasattr(frame.params, "enth_mol_form_ig_comp_ref")

    expr = Constant.enth_mol_ig_comp.return_expression(
        frame.props[1], frame.params, frame.props[1].temperature
    )
    assert value(expr) == 0

    frame.props[1].temperature.value = 400
    assert value(expr) == pytest.approx(2911.485, rel=1e-3)

    assert_units_equivalent(expr, pyunits.J / pyunits.mol)


@pytest.mark.unit
def test_entr_mol_ig_comp(frame):
    Constant.entr_mol_ig_comp.build_parameters(frame.params)

    assert isinstance(frame.params.entr_mol_form_ig_comp_ref, Var)
    assert value(frame.params.entr_mol_form_ig_comp_ref) == 0.0

    expr = Constant.entr_mol_ig_comp.return_expression(
        frame.props[1], frame.params, frame.props[1].temperature
    )
    assert value(expr) == value(frame.params.entr_mol_form_ig_comp_ref)

    frame.props[1].temperature.value = 400
    assert value(expr) == pytest.approx(8.37582, rel=1e-3)

    assert_units_equivalent(expr, pyunits.J / pyunits.mol / pyunits.K)


@pytest.mark.unit
def test_cp_mol_sol_comp(frame):
    Constant.cp_mol_sol_comp.build_parameters(frame.params)

    assert isinstance(frame.params.cp_mol_sol_comp_coeff, Var)
    assert value(frame.params.cp_mol_sol_comp_coeff) == pytest.approx(
        100000 / 1000, rel=1e-5
    )

    expr = Constant.cp_mol_sol_comp.return_expression(
        frame.props[1], frame.params, frame.props[1].temperature
    )
    assert value(expr) == pytest.approx(100000 / 1000, rel=1e-5)

    frame.props[1].temperature.value = 400
    assert value(expr) == pytest.approx(100000 / 1000, rel=1e-5)

    assert_units_equivalent(expr, pyunits.J / pyunits.mol / pyunits.K)


@pytest.mark.unit
def test_enth_mol_sol_comp(frame):

    Constant.enth_mol_sol_comp.build_parameters(frame.params)

    assert isinstance(frame.params.enth_mol_form_sol_comp_ref, Var)
    assert value(frame.params.enth_mol_form_sol_comp_ref) == -300e3

    expr = Constant.enth_mol_sol_comp.return_expression(
        frame.props[1], frame.params, frame.props[1].temperature
    )
    assert value(expr) == value(frame.params.enth_mol_form_sol_comp_ref)

    frame.props[1].temperature.value = 301
    assert value(expr) == pytest.approx(-299900.0, rel=1e-3)

    frame.props[1].temperature.value = 400
    assert value(expr) == pytest.approx(-290000.0, rel=1e-3)

    assert_units_equivalent(expr, pyunits.J / pyunits.mol)


@pytest.mark.unit
def test_enth_mol_sol_comp_no_formation(frame):
    frame.config.include_enthalpy_of_formation = False
    frame.params.config.include_enthalpy_of_formation = False

    Constant.enth_mol_sol_comp.build_parameters(frame.params)

    assert not hasattr(frame.params, "enth_mol_form_sol_comp_ref")

    expr = Constant.enth_mol_sol_comp.return_expression(
        frame.props[1], frame.params, frame.props[1].temperature
    )
    assert value(expr) == 0

    frame.props[1].temperature.value = 400
    assert value(expr) == pytest.approx(10000, rel=1e-3)

    assert_units_equivalent(expr, pyunits.J / pyunits.mol)


@pytest.mark.unit
def test_entr_mol_sol_comp(frame):
    Constant.entr_mol_sol_comp.build_parameters(frame.params)

    assert isinstance(frame.params.entr_mol_form_sol_comp_ref, Var)
    assert value(frame.params.entr_mol_form_sol_comp_ref) == 50

    expr = Constant.entr_mol_sol_comp.return_expression(
        frame.props[1], frame.params, frame.props[1].temperature
    )
    assert value(expr) == value(frame.params.entr_mol_form_sol_comp_ref)

    frame.props[1].temperature.value = 400
    assert value(expr) == pytest.approx(78.77, rel=1e-3)

    assert_units_equivalent(expr, pyunits.J / pyunits.mol / pyunits.K)


@pytest.mark.unit
def test_dens_mol_sol_comp(frame):
    Constant.dens_mol_sol_comp.build_parameters(frame.params)

    assert isinstance(frame.params.dens_mol_sol_comp_coeff, Var)
    assert value(frame.params.dens_mol_sol_comp_coeff) == pytest.approx(100e3, rel=1e-4)

    expr = Constant.dens_mol_sol_comp.return_expression(
        frame.props[1], frame.params, frame.props[1].temperature
    )
    assert value(expr) == pytest.approx(100e3, rel=1e-4)

    frame.props[1].temperature.value = 400
    assert value(expr) == pytest.approx(100e3, rel=1e-4)

    assert_units_equivalent(expr, pyunits.mol / pyunits.m**3)
