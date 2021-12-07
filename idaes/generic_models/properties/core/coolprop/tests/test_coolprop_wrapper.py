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
Tests for methods from CoolProp wrapper

Authors: Andrew Lee
"""

import pytest

try:
    from CoolProp.CoolProp import PropsSI
    coolprop_present = True
except ModuleNotFoundError:
    coolprop_present = False

from pyomo.environ import ConcreteModel, Param, units as pyunits, value, Var
from pyomo.util.check_units import assert_units_equivalent

from idaes.core import FlowsheetBlock, Component, LiquidPhase
from idaes.generic_models.properties.core.eos.ceos import \
    cubic_roots_available
from idaes.generic_models.properties.core.generic.generic_property import (
        GenericParameterBlock)
from idaes.generic_models.properties.core.state_definitions import FTPx
from idaes.generic_models.properties.core.eos.ceos import Cubic, CubicType


from idaes.generic_models.properties.core.coolprop.coolprop_wrapper import \
    CoolPropWrapper


@pytest.mark.skipif(not coolprop_present, reason="CoolProp not installed")
class TestWrapper:
    @pytest.mark.unit
    def test_load_component(self):
        # Clear cached components to be sure
        CoolPropWrapper.cached_components = {}

        # Load parameters for oxygen
        prop_dict = CoolPropWrapper._load_component_data("Oxygen")

        assert prop_dict is not None
        assert "Oxygen" in CoolPropWrapper.cached_components
        assert CoolPropWrapper.cached_components["Oxygen"] is prop_dict
        assert isinstance(prop_dict, dict)
        for k in ["STATES", "ANCILLARIES"]:
            assert k in prop_dict

    @pytest.mark.unit
    def test_load_component_invalid(self):
        with pytest.raises(RuntimeError,
                           match="Failed to find component foo in CoolProp "
                           "JSON database."):
            CoolPropWrapper._load_component_data("foo")

    @pytest.mark.unit
    def test_get_component(self):
        prop_dict = CoolPropWrapper._get_component_data("Oxygen")

        assert CoolPropWrapper.cached_components["Oxygen"] is prop_dict

    @pytest.mark.unit
    def test_get_component_alias(self):
        # Load oxygen properties using aliases
        prop_dict = CoolPropWrapper._get_component_data("R732")

        assert CoolPropWrapper.cached_components["Oxygen"] is prop_dict
        assert "R732" not in CoolPropWrapper.cached_components

    @pytest.mark.unit
    def test_load_component_by_alias(self):
        # Load CO2 data using one of its aliases
        prop_dict = CoolPropWrapper._load_component_data("R744")

        assert CoolPropWrapper.cached_components["R744"] is prop_dict
        assert "R744" in CoolPropWrapper.cached_components

        # Retrieve CO2 data using its normal name
        prop_dict2 = CoolPropWrapper._get_component_data("CO2")

        assert prop_dict2 is prop_dict
        assert "CO2" not in CoolPropWrapper.cached_components

    @pytest.mark.unit
    def test_get_component_invalid(self):
        with pytest.raises(RuntimeError,
                           match="Failed to find component foo in CoolProp "
                           "JSON database."):
            CoolPropWrapper._get_component_data("foo")

    @pytest.mark.unit
    def test_flush_cached_components(self):
        CoolPropWrapper.flush_cached_components()

        assert CoolPropWrapper.cached_components == {}

    @pytest.mark.unit
    def test_get_critical_properties(self):
        Tc = CoolPropWrapper._get_critical_property(
            "Acetone", "T")
        assert Tc == (508.1, pyunits.K)

        Pc = CoolPropWrapper._get_critical_property(
            "Acetone", "p")
        assert Pc == (4700000.0, pyunits.Pa)

        rhoc = CoolPropWrapper._get_critical_property(
            "Acetone", "rhomolar")
        assert rhoc[0] == 4699.999999999999
        assert_units_equivalent(rhoc[1], pyunits.mol/pyunits.m**3)

        hc = CoolPropWrapper._get_critical_property(
            "Acetone", "hmolar")
        assert hc[0] == 31614.73051047263
        assert_units_equivalent(hc[1], pyunits.J/pyunits.mol)

        sc = CoolPropWrapper._get_critical_property(
            "Acetone", "smolar")
        assert sc[0] == 72.97112978635582
        assert_units_equivalent(sc[1], pyunits.J/pyunits.mol/pyunits.K)

    @pytest.mark.unit
    def test_get_eos_properties(self):
        omega = CoolPropWrapper._get_eos_property("Acetone", "acentric")
        assert omega == (0.3071, pyunits.dimensionless)

        mw = CoolPropWrapper._get_eos_property("Acetone", "molar_mass")
        assert mw[0] == 0.05807914
        assert_units_equivalent(mw[1], pyunits.kg/pyunits.mol)

    @pytest.mark.unit
    def test_get_parameter_value(self):
        Tc = CoolPropWrapper.get_parameter_value(
            "Acetone", "temperature_crit")
        assert Tc == (508.1, pyunits.K)

        Pc = CoolPropWrapper.get_parameter_value(
            "Acetone", "pressure_crit")
        assert Pc == (4700000.0, pyunits.Pa)

        rhoc = CoolPropWrapper.get_parameter_value(
            "Acetone", "dens_mol_crit")
        assert rhoc[0] == 4699.999999999999
        assert_units_equivalent(rhoc[1], pyunits.mol/pyunits.m**3)

        hc = CoolPropWrapper.get_parameter_value(
            "Acetone", "enth_mol_crit")
        assert hc[0] == 31614.73051047263
        assert_units_equivalent(hc[1], pyunits.J/pyunits.mol)

        sc = CoolPropWrapper.get_parameter_value(
            "Acetone", "entr_mol_crit")
        assert sc[0] == 72.97112978635582
        assert_units_equivalent(sc[1], pyunits.J/pyunits.mol/pyunits.K)

        omega = CoolPropWrapper.get_parameter_value("Acetone", "omega")
        assert omega == (0.3071, pyunits.dimensionless)

        mw = CoolPropWrapper.get_parameter_value("Acetone", "mw")
        assert mw[0] == 0.05807914
        assert_units_equivalent(mw[1], pyunits.kg/pyunits.mol)


@pytest.mark.skipif(not coolprop_present, reason="CoolProp not installed")
class TestCoolPropIntegration(object):
    @pytest.fixture(scope="class")
    def m(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(default={'dynamic': False})

        configuration = {
            # Specifying components
            "components": {
                'benzene': {"type": Component,
                            # "enth_mol_ig_comp": RPP4,
                            # "entr_mol_ig_comp": RPP4,
                            "pressure_sat_comp": CoolPropWrapper,
                            "parameter_data": {
                                "mw": CoolPropWrapper,
                                "pressure_crit": CoolPropWrapper,
                                "temperature_crit": CoolPropWrapper,
                                "omega": CoolPropWrapper}}},
            # Specifying phases
            "phases":  {'Liq': {"type": LiquidPhase,
                                "equation_of_state": Cubic,
                                "equation_of_state_options": {
                                    "type": CubicType.PR}}},

            # Set base units of measurement
            "base_units": {"time": pyunits.s,
                           "length": pyunits.m,
                           "mass": pyunits.kg,
                           "amount": pyunits.mol,
                           "temperature": pyunits.K},

            # Specifying state definition
            "state_definition": FTPx,
            "state_bounds": {"flow_mol": (0, 100, 1000, pyunits.mol/pyunits.s),
                             "temperature": (273.15, 300, 500, pyunits.K),
                             "pressure": (5e4, 1e5, 1e6, pyunits.Pa)},
            "pressure_ref": (101325, pyunits.Pa),
            "temperature_ref": (298.15, pyunits.K),

            "parameter_data": {"PR_kappa": {("benzene", "benzene"): 0.000}}}

        m.fs.props = GenericParameterBlock(default=configuration)

        m.fs.state = m.fs.props.build_state_block(
            [0], default={"defined_state": True})

        m.fs.state[0].flow_mol.fix(1)
        m.fs.state[0].pressure.fix(101325)
        m.fs.state[0].temperature.fix(300)
        m.fs.state[0].mole_frac_comp["benzene"].fix(1)

        return m

    @pytest.mark.unit
    def test_physical_constants(self, m):
        # Benzene parameters
        assert isinstance(m.fs.props.benzene.temperature_crit, Var)
        assert m.fs.props.benzene.temperature_crit.fixed
        assert value(m.fs.props.benzene.temperature_crit) == PropsSI(
            "TCRIT", "Benzene")

        assert isinstance(m.fs.props.benzene.pressure_crit, Var)
        assert m.fs.props.benzene.pressure_crit.fixed
        assert value(m.fs.props.benzene.pressure_crit) == PropsSI(
            "PCRIT", "Benzene")

        assert isinstance(m.fs.props.benzene.mw, Param)
        assert value(m.fs.props.benzene.mw) == PropsSI(
            "molarmass", "Benzene")

        assert isinstance(m.fs.props.benzene.omega, Var)
        assert m.fs.props.benzene.omega.fixed
        assert value(m.fs.props.benzene.omega) == PropsSI(
            "acentric", "Benzene")

    @pytest.mark.unit
    def test_psat(self, m):
        assert isinstance(m.fs.props.benzene.pressure_sat_coeff_n1, Var)
        assert isinstance(m.fs.props.benzene.pressure_sat_coeff_n2, Var)
        assert isinstance(m.fs.props.benzene.pressure_sat_coeff_n3, Var)
        assert isinstance(m.fs.props.benzene.pressure_sat_coeff_n4, Var)
        assert isinstance(m.fs.props.benzene.pressure_sat_coeff_n5, Var)
        assert isinstance(m.fs.props.benzene.pressure_sat_coeff_n6, Var)

        assert isinstance(m.fs.props.benzene.pressure_sat_coeff_t1, Var)
        assert isinstance(m.fs.props.benzene.pressure_sat_coeff_t2, Var)
        assert isinstance(m.fs.props.benzene.pressure_sat_coeff_t3, Var)
        assert isinstance(m.fs.props.benzene.pressure_sat_coeff_t4, Var)
        assert isinstance(m.fs.props.benzene.pressure_sat_coeff_t5, Var)
        assert isinstance(m.fs.props.benzene.pressure_sat_coeff_t6, Var)

        assert m.fs.props.benzene.pressure_sat_coeff_n1.value == \
            0.005561906558935796
        assert m.fs.props.benzene.pressure_sat_coeff_n1.fixed
        assert m.fs.props.benzene.pressure_sat_coeff_n2.value == \
            -0.08662136922915314
        assert m.fs.props.benzene.pressure_sat_coeff_n2.fixed
        assert m.fs.props.benzene.pressure_sat_coeff_n3.value == \
            -6.964182734154488
        assert m.fs.props.benzene.pressure_sat_coeff_n3.fixed
        assert m.fs.props.benzene.pressure_sat_coeff_n4.value == \
            1.1249288132278856
        assert m.fs.props.benzene.pressure_sat_coeff_n4.fixed
        assert m.fs.props.benzene.pressure_sat_coeff_n5.value == \
            -3.961859460597414
        assert m.fs.props.benzene.pressure_sat_coeff_n5.fixed
        assert m.fs.props.benzene.pressure_sat_coeff_n6.value == \
            -13.106880507410812
        assert m.fs.props.benzene.pressure_sat_coeff_n6.fixed

        assert m.fs.props.benzene.pressure_sat_coeff_t1.value == 0.037
        assert m.fs.props.benzene.pressure_sat_coeff_t1.fixed
        assert m.fs.props.benzene.pressure_sat_coeff_t2.value == 0.505
        assert m.fs.props.benzene.pressure_sat_coeff_t2.fixed
        assert m.fs.props.benzene.pressure_sat_coeff_t3.value == 1.014
        assert m.fs.props.benzene.pressure_sat_coeff_t3.fixed
        assert m.fs.props.benzene.pressure_sat_coeff_t4.value == 1.469
        assert m.fs.props.benzene.pressure_sat_coeff_t4.fixed
        assert m.fs.props.benzene.pressure_sat_coeff_t5.value == 3.711
        assert m.fs.props.benzene.pressure_sat_coeff_t5.fixed
        assert m.fs.props.benzene.pressure_sat_coeff_t6.value == 12.647
        assert m.fs.props.benzene.pressure_sat_coeff_t6.fixed

        for T in range(300, 401, 10):
            m.fs.state[0].temperature.fix(T)
            assert pytest.approx(
                PropsSI("P", "T", T, "Q", 0.5, "benzene"), rel=5e-4) == value(
                    m.fs.state[0].pressure_sat_comp["benzene"])
