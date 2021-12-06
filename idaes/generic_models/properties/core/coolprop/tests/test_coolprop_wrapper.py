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
    from CoolProp import CoolProp
    coolprop_present = True
except ModuleNotFoundError:
    coolprop_present = False

from pyomo.environ import units as pyunits
from pyomo.util.check_units import assert_units_equivalent

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
        prop_dict = CoolPropWrapper.get_component_data("Oxygen")

        assert CoolPropWrapper.cached_components["Oxygen"] is prop_dict

    @pytest.mark.unit
    def test_get_component_alias(self):
        # Load oxygen properties using aliases
        prop_dict = CoolPropWrapper.get_component_data("R732")

        assert CoolPropWrapper.cached_components["Oxygen"] is prop_dict
        assert "R732" not in CoolPropWrapper.cached_components

    @pytest.mark.unit
    def test_get_component_invalid(self):
        with pytest.raises(RuntimeError,
                           match="Failed to find component foo in CoolProp "
                           "JSON database."):
            CoolPropWrapper.get_component_data("foo")

    @pytest.mark.unit
    def test_flush_cached_components(self):
        CoolPropWrapper.flush_cached_components()

        assert CoolPropWrapper.cached_components == {}

    @pytest.mark.unit
    def test_get_critical_properties(self):
        Tc = CoolPropWrapper.get_critical_property(
            "Acetone", "temperature_crit")
        assert Tc == (508.1, pyunits.K)

        Pc = CoolPropWrapper.get_critical_property(
            "Acetone", "pressure_crit")
        assert Pc == (4700000.0, pyunits.Pa)

        rhoc = CoolPropWrapper.get_critical_property(
            "Acetone", "dens_mol_crit")
        assert rhoc[0] == 4699.999999999999
        assert_units_equivalent(rhoc[1], pyunits.mol/pyunits.m**3)

        hc = CoolPropWrapper.get_critical_property(
            "Acetone", "enth_mol_crit")
        assert hc[0] == 31614.73051047263
        assert_units_equivalent(hc[1], pyunits.J/pyunits.mol)

        sc = CoolPropWrapper.get_critical_property(
            "Acetone", "entr_mol_crit")
        assert sc[0] == 72.97112978635582
        assert_units_equivalent(sc[1], pyunits.J/pyunits.mol/pyunits.K)

    @pytest.mark.unit
    def test_get_eos_properties(self):
        omega = CoolPropWrapper.get_eos_property("Acetone", "omega")
        assert omega == (0.3071, pyunits.dimensionless)

        mw = CoolPropWrapper.get_eos_property("Acetone", "mw")
        assert mw[0] == 0.05807914
        assert_units_equivalent(mw[1], pyunits.kg/pyunits.mol)
