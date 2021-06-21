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
from pyomo.environ import ConcreteModel, value, Var, units as pyunits
from pyomo.common.collections import ComponentSet
from idaes.generic_models.properties import iapws95

from idaes.core import MaterialBalanceType, EnergyBalanceType
from idaes.core.util.exceptions import ConfigurationError

# Set module level pyest marker
pytestmark = pytest.mark.iapws
prop_available = iapws95.iapws95_available()


# -----------------------------------------------------------------------------
# Test Enums and common functions
@pytest.mark.unit
def test_htpx_invalid_args():
    with pytest.raises(ConfigurationError):
        iapws95.htpx(300*pyunits.K, P=101325*pyunits.Pa, x=0.5)

    with pytest.raises(ConfigurationError):
        iapws95.htpx(100*pyunits.K, x=0.5)

    with pytest.raises(ConfigurationError):
        iapws95.htpx(5e3*pyunits.K, x=0.5)

    with pytest.raises(ConfigurationError):
        iapws95.htpx(300*pyunits.K, P=1e-10*pyunits.Pa)

    with pytest.raises(ConfigurationError):
        iapws95.htpx(300*pyunits.K, P=1e10*pyunits.Pa)

    with pytest.raises(ConfigurationError):
        iapws95.htpx(300*pyunits.K, x=-1)

    with pytest.raises(ConfigurationError):
        iapws95.htpx(300*pyunits.K, x=2)


@pytest.mark.skipif(not prop_available, reason="IAPWS not available")
@pytest.mark.unit
def test_htpx():
    # Compariosn of IAPWS95 results to data from Spirax-Sarco steam tables
    # https://www.spiraxsarco.com/resources-and-design-tools/steam-tables
    # Retrieved 20 Nov 2019

    # MW for unit conversions
    mw = 0.01801528

    # There appears to be a small difference between the reference states
    # for IAPWS package and the reference source.
    offset = 9.22

    # Subcooled liquid
    assert iapws95.htpx(300*pyunits.K, P=101325*pyunits.Pa) == pytest.approx(
        112143 * mw + offset, 1e-5)

    # Saturated liquid
    assert iapws95.htpx(500*pyunits.K, x=0) == pytest.approx(
        974919 * mw + offset, 1e-5)

    # Wet steam
    assert iapws95.htpx(550*pyunits.K, x=0.5) == pytest.approx(
        2.00138e06 * mw + offset, 1e-5)

    # Saturated vapor
    assert iapws95.htpx(600*pyunits.K, x=1) == pytest.approx(
        2.67730e06 * mw + offset, 1e-5)

    # Superheated steam
    assert iapws95.htpx(400*pyunits.K, P=101325*pyunits.Pa) == pytest.approx(
        2.72979e06 * mw + offset, 1e-5)


@pytest.mark.unit
def test_PhaseType():
    assert len(iapws95.PhaseType) == 4

    # Check that expected values do not raise Exceptions
    iapws95.PhaseType.MIX
    iapws95.PhaseType.LG
    iapws95.PhaseType.L
    iapws95.PhaseType.G


@pytest.mark.unit
def test_StateVars():
    assert len(iapws95.StateVars) == 2

    # Check that expected values do not raise Exceptions
    iapws95.StateVars.PH
    iapws95.StateVars.TPX


# -----------------------------------------------------------------------------
# Test builds with different phase presentations and state vars
@pytest.mark.skipif(not prop_available, reason="IAPWS not available")
class TestMixPh(object):
    # This should be the default option, so test with no arguments
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.params = iapws95.Iapws95ParameterBlock(default={
                "temperature_bounds": (249*pyunits.K, 2501*pyunits.K),
                "pressure_bounds": (0.11*pyunits.Pa, 1.1e9*pyunits.Pa),
                "enthalpy_mol_bounds": (
                    0.1*pyunits.J/pyunits.mol, 1.1e5*pyunits.J/pyunits.mol),
            }
        )
        return model

    @pytest.mark.unit
    def test_config_enth_mass(self, model):
        model = ConcreteModel()
        l = 0.1/0.018
        u = 1.1e5/0.018
        model.params = iapws95.Iapws95ParameterBlock(default={
                "enthalpy_mass_bounds": (
                    l*pyunits.J/pyunits.kg,
                    u*pyunits.J/pyunits.kg)
            }
        )
        assert pytest.approx(0.1, rel=1e-3) == value(
            model.params.default_enthalpy_bounds[0])
        assert pytest.approx(1.1e5, rel=1e-3) == value(
            model.params.default_enthalpy_bounds[1])

    @pytest.mark.unit
    def test_config(self, model):
        assert model.params.config.phase_presentation == iapws95.PhaseType.MIX
        assert model.params.config.state_vars == iapws95.StateVars.PH

        assert len(model.params.phase_list) == 1
        for i in model.params.phase_list:
            assert i in ["Mix"]

        assert pytest.approx(249) == value(
            model.params.default_temperature_bounds[0])
        assert pytest.approx(2501) == value(
            model.params.default_temperature_bounds[1])
        assert pytest.approx(0.11) == value(
            model.params.default_pressure_bounds[0])
        assert pytest.approx(1.1e9) == value(
            model.params.default_pressure_bounds[1])
        assert pytest.approx(0.1) == value(
            model.params.default_enthalpy_bounds[0])
        assert pytest.approx(1.1e5) == value(
            model.params.default_enthalpy_bounds[1])

    @pytest.mark.unit
    def test_build(self, model):
        model.prop = iapws95.Iapws95StateBlock(
            [1], default={"parameters": model.params}
        )

        assert isinstance(model.prop[1].flow_mol, Var)
        assert isinstance(model.prop[1].pressure, Var)
        assert isinstance(model.prop[1].enth_mol, Var)

        assert isinstance(model.prop[1].extensive_set, ComponentSet)
        assert isinstance(model.prop[1].intensive_set, ComponentSet)

    @pytest.mark.unit
    def test_get_material_flow_terms(self, model):
        for p in model.params.phase_list:
            for j in model.params.component_list:
                assert value(
                    model.prop[1].get_material_flow_terms(p, j)) == value(
                    model.prop[1].flow_mol
                )

    @pytest.mark.unit
    def test_get_enthalpy_flow_terms(self, model):
        for p in model.params.phase_list:
            assert (
                value(model.prop[1].get_enthalpy_flow_terms(p))
                == value(model.prop[1].enth_mol * model.prop[1].flow_mol)
            )

    @pytest.mark.unit
    def test_get_material_density_terms(self, model):
        for p in model.params.phase_list:
            for j in model.params.component_list:
                assert (
                    model.prop[1].get_material_density_terms(p, j)
                    is model.prop[1].dens_mol
                )

    @pytest.mark.unit
    def test_get_energy_density_terms(self, model):
        for p in model.params.phase_list:
            assert value(
                model.prop[1].get_energy_density_terms(p)) == value(
                model.prop[1].dens_mol * model.prop[1].energy_internal_mol
            )

    @pytest.mark.unit
    def test_default_material_balance_type(self, model):
        assert (
            model.prop[1].default_material_balance_type()
            is MaterialBalanceType.componentTotal
        )

    @pytest.mark.unit
    def test_default_energy_balance_type(self, model):
        assert (
            model.prop[1].default_energy_balance_type()
            is EnergyBalanceType.enthalpyTotal
        )

    @pytest.mark.unit
    def test_define_state_vars(self, model):
        sv = model.prop[1].define_state_vars()
        assert len(sv) == 3
        for i in sv:
            assert i in ["flow_mol", "enth_mol", "pressure"]

    @pytest.mark.unit
    def test_define_port_members(self, model):
        sv = model.prop[1].define_port_members()
        assert len(sv) == 3
        for i in sv:
            assert i in ["flow_mol", "enth_mol", "pressure"]

    @pytest.mark.unit
    def test_define_display_vars(self, model):
        dv = model.prop[1].define_display_vars()
        assert len(dv) == 6
        for i in dv:
            assert i in [
                "Molar Flow (mol/s)",
                "Mass Flow (kg/s)",
                "T (K)",
                "P (Pa)",
                "Vapor Fraction",
                "Molar Enthalpy (J/mol)",
            ]

    @pytest.mark.unit
    def test_extensive_state_vars(self, model):
        for i in model.prop[1].extensive_set:
            assert i is model.prop[1].flow_mol

    @pytest.mark.unit
    def test_intensive_state_vars(self, model):
        for i in model.prop[1].intensive_set:
            assert i.name in ["prop[1].enth_mol", "prop[1].pressure"]

    @pytest.mark.unit
    def test_model_check(self, model):
        assert model.prop[1].model_check() is None

    @pytest.mark.unit
    def test_data_initialize(self, model):
        assert model.prop[1].initialize() is None

    @pytest.mark.unit
    def test_initialize(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

        flags = model.prop.initialize(hold_state=True)

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].enth_mol.fixed
        assert model.prop[1].pressure.fixed

        model.prop.release_state(flags)

        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

    @pytest.mark.unit
    def test_initialize_w_state(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

        flags = model.prop.initialize(
            hold_state=True,
            state_args={"flow_mol": 20, "enth_mol": 200, "pressure": 2000},
        )

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].enth_mol.fixed
        assert model.prop[1].pressure.fixed

        assert model.prop[1].flow_mol.value == 20
        assert model.prop[1].enth_mol.value == 200
        assert model.prop[1].pressure.value == 2000

        model.prop.release_state(flags)

        assert model.prop[1].flow_mol.value == 20
        assert model.prop[1].enth_mol.value == 200
        assert model.prop[1].pressure.value == 2000

        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

    @pytest.mark.unit
    def test_initialize_w_partial_state_1(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

        flags = model.prop.initialize(hold_state=True,
                                      state_args={"flow_mol": 30})

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].enth_mol.fixed
        assert model.prop[1].pressure.fixed

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].enth_mol.value == 200
        assert model.prop[1].pressure.value == 2000

        model.prop.release_state(flags)

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].enth_mol.value == 200
        assert model.prop[1].pressure.value == 2000

        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

    @pytest.mark.unit
    def test_initialize_w_partial_state_2(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

        flags = model.prop.initialize(hold_state=True,
                                      state_args={"enth_mol": 300})

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].enth_mol.fixed
        assert model.prop[1].pressure.fixed

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].enth_mol.value == 300
        assert model.prop[1].pressure.value == 2000

        model.prop.release_state(flags)

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].enth_mol.value == 300
        assert model.prop[1].pressure.value == 2000

        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

    @pytest.mark.unit
    def test_initialize_w_partial_state_3(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

        flags = model.prop.initialize(hold_state=True,
                                      state_args={"pressure": 3000})

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].enth_mol.fixed
        assert model.prop[1].pressure.fixed

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].enth_mol.value == 300
        assert model.prop[1].pressure.value == 3000

        model.prop.release_state(flags)

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].enth_mol.value == 300
        assert model.prop[1].pressure.value == 3000

        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

    @pytest.mark.unit
    def test_initialize_w_state_ignore_fixed(self, model):
        model.prop[1].flow_mol.fixed = True
        model.prop[1].enth_mol.fixed = True
        model.prop[1].pressure.fixed = True

        flags = model.prop.initialize(
            hold_state=True,
            state_args={"flow_mol": 20, "enth_mol": 200, "pressure": 2000},
        )

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].enth_mol.fixed
        assert model.prop[1].pressure.fixed

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].enth_mol.value == 300
        assert model.prop[1].pressure.value == 3000

        model.prop.release_state(flags)

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].enth_mol.value == 300
        assert model.prop[1].pressure.value == 3000

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].enth_mol.fixed
        assert model.prop[1].pressure.fixed


@pytest.mark.skipif(not prop_available, reason="IAPWS not available")
class TestLGPh(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.params = iapws95.Iapws95ParameterBlock(
            default={"phase_presentation": iapws95.PhaseType.LG}
        )

        return model

    @pytest.mark.unit
    def test_config(self, model):
        assert model.params.config.phase_presentation == iapws95.PhaseType.LG
        assert model.params.config.state_vars == iapws95.StateVars.PH

        assert len(model.params.phase_list) == 2
        for i in model.params.phase_list:
            assert i in ["Liq", "Vap"]

    @pytest.mark.unit
    def test_build(self, model):
        model.prop = iapws95.Iapws95StateBlock(
            [1], default={"parameters": model.params}
        )

        assert isinstance(model.prop[1].flow_mol, Var)
        assert isinstance(model.prop[1].pressure, Var)
        assert isinstance(model.prop[1].enth_mol, Var)

        assert isinstance(model.prop[1].extensive_set, ComponentSet)
        assert isinstance(model.prop[1].intensive_set, ComponentSet)

    @pytest.mark.unit
    def test_get_material_flow_terms(self, model):
        for p in model.params.phase_list:
            for j in model.params.component_list:
                assert value(
                    model.prop[1].get_material_flow_terms(p, j)) == value(
                    model.prop[1].flow_mol * model.prop[1].phase_frac[p]
                )

    @pytest.mark.unit
    def test_get_enthalpy_flow_terms(self, model):
        for p in model.params.phase_list:
            assert value(model.prop[1].get_enthalpy_flow_terms(p)) == value(
                model.prop[1].enth_mol_phase[p]
                * model.prop[1].phase_frac[p]
                * model.prop[1].flow_mol
            )

    @pytest.mark.unit
    def test_get_material_density_terms(self, model):
        for p in model.params.phase_list:
            for j in model.params.component_list:
                assert (
                    model.prop[1].get_material_density_terms(p, j)
                    is model.prop[1].dens_mol_phase[p]
                )

    @pytest.mark.unit
    def test_get_energy_density_terms(self, model):
        for p in model.params.phase_list:
            assert value(model.prop[1].get_energy_density_terms(p)) == value(
                model.prop[1].dens_mol_phase[p]
                * model.prop[1].energy_internal_mol_phase[p]
            )

    @pytest.mark.unit
    def test_default_material_balance_type(self, model):
        assert (
            model.prop[1].default_material_balance_type()
            is MaterialBalanceType.componentTotal
        )

    @pytest.mark.unit
    def test_default_energy_balance_type(self, model):
        assert (
            model.prop[1].default_energy_balance_type()
            is EnergyBalanceType.enthalpyTotal
        )

    @pytest.mark.unit
    def test_define_state_vars(self, model):
        sv = model.prop[1].define_state_vars()
        assert len(sv) == 3
        for i in sv:
            assert i in ["flow_mol", "enth_mol", "pressure"]

    @pytest.mark.unit
    def test_define_port_members(self, model):
        sv = model.prop[1].define_port_members()
        assert len(sv) == 3
        for i in sv:
            assert i in ["flow_mol", "enth_mol", "pressure"]

    @pytest.mark.unit
    def test_define_display_vars(self, model):
        dv = model.prop[1].define_display_vars()
        assert len(dv) == 6
        for i in dv:
            assert i in [
                "Molar Flow (mol/s)",
                "Mass Flow (kg/s)",
                "T (K)",
                "P (Pa)",
                "Vapor Fraction",
                "Molar Enthalpy (J/mol)",
            ]

    @pytest.mark.unit
    def test_extensive_state_vars(self, model):
        for i in model.prop[1].extensive_set:
            assert i is model.prop[1].flow_mol

    @pytest.mark.unit
    def test_intensive_state_vars(self, model):
        for i in model.prop[1].intensive_set:
            assert i.name in ["prop[1].enth_mol", "prop[1].pressure"]

    @pytest.mark.unit
    def test_model_check(self, model):
        assert model.prop[1].model_check() is None

    @pytest.mark.unit
    def test_data_initialize(self, model):
        assert model.prop[1].initialize() is None

    @pytest.mark.unit
    def test_initialize(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

        flags = model.prop.initialize(hold_state=True)

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].enth_mol.fixed
        assert model.prop[1].pressure.fixed

        model.prop.release_state(flags)

        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

    @pytest.mark.unit
    def test_initialize_w_state(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

        flags = model.prop.initialize(
            hold_state=True,
            state_args={"flow_mol": 20, "enth_mol": 200, "pressure": 2000},
        )

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].enth_mol.fixed
        assert model.prop[1].pressure.fixed

        assert model.prop[1].flow_mol.value == 20
        assert model.prop[1].enth_mol.value == 200
        assert model.prop[1].pressure.value == 2000

        model.prop.release_state(flags)

        assert model.prop[1].flow_mol.value == 20
        assert model.prop[1].enth_mol.value == 200
        assert model.prop[1].pressure.value == 2000

        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

    @pytest.mark.unit
    def test_initialize_w_partial_state_1(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

        flags = model.prop.initialize(hold_state=True,
                                      state_args={"flow_mol": 30})

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].enth_mol.fixed
        assert model.prop[1].pressure.fixed

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].enth_mol.value == 200
        assert model.prop[1].pressure.value == 2000

        model.prop.release_state(flags)

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].enth_mol.value == 200
        assert model.prop[1].pressure.value == 2000

        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

    @pytest.mark.unit
    def test_initialize_w_partial_state_2(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

        flags = model.prop.initialize(hold_state=True,
                                      state_args={"enth_mol": 300})

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].enth_mol.fixed
        assert model.prop[1].pressure.fixed

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].enth_mol.value == 300
        assert model.prop[1].pressure.value == 2000

        model.prop.release_state(flags)

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].enth_mol.value == 300
        assert model.prop[1].pressure.value == 2000

        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

    @pytest.mark.unit
    def test_initialize_w_partial_state_3(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

        flags = model.prop.initialize(hold_state=True,
                                      state_args={"pressure": 3000})

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].enth_mol.fixed
        assert model.prop[1].pressure.fixed

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].enth_mol.value == 300
        assert model.prop[1].pressure.value == 3000

        model.prop.release_state(flags)

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].enth_mol.value == 300
        assert model.prop[1].pressure.value == 3000

        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

    @pytest.mark.unit
    def test_initialize_w_state_ignore_fixed(self, model):
        model.prop[1].flow_mol.fixed = True
        model.prop[1].enth_mol.fixed = True
        model.prop[1].pressure.fixed = True

        flags = model.prop.initialize(
            hold_state=True,
            state_args={"flow_mol": 20, "enth_mol": 200, "pressure": 2000},
        )

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].enth_mol.fixed
        assert model.prop[1].pressure.fixed

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].enth_mol.value == 300
        assert model.prop[1].pressure.value == 3000

        model.prop.release_state(flags)

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].enth_mol.value == 300
        assert model.prop[1].pressure.value == 3000

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].enth_mol.fixed
        assert model.prop[1].pressure.fixed


@pytest.mark.skipif(not prop_available, reason="IAPWS not available")
class TestLPh(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.params = iapws95.Iapws95ParameterBlock(
            default={"phase_presentation": iapws95.PhaseType.L}
        )

        return model

    @pytest.mark.unit
    def test_config(self, model):
        assert model.params.config.phase_presentation == iapws95.PhaseType.L
        assert model.params.config.state_vars == iapws95.StateVars.PH

        assert len(model.params.phase_list) == 1
        for i in model.params.phase_list:
            assert i in ["Liq"]

    @pytest.mark.unit
    def test_build(self, model):
        model.prop = iapws95.Iapws95StateBlock(
            [1], default={"parameters": model.params}
        )

        assert isinstance(model.prop[1].flow_mol, Var)
        assert isinstance(model.prop[1].pressure, Var)
        assert isinstance(model.prop[1].enth_mol, Var)

        assert isinstance(model.prop[1].extensive_set, ComponentSet)
        assert isinstance(model.prop[1].intensive_set, ComponentSet)

    @pytest.mark.unit
    def test_get_material_flow_terms(self, model):
        for p in model.params.phase_list:
            for j in model.params.component_list:
                assert value(
                    model.prop[1].get_material_flow_terms(p, j)) == value(
                    model.prop[1].flow_mol * model.prop[1].phase_frac[p]
                )

    @pytest.mark.unit
    def test_get_enthalpy_flow_terms(self, model):
        for p in model.params.phase_list:
            assert value(model.prop[1].get_enthalpy_flow_terms(p)) == value(
                model.prop[1].enth_mol_phase[p]
                * model.prop[1].phase_frac[p]
                * model.prop[1].flow_mol
            )

    @pytest.mark.unit
    def test_get_material_density_terms(self, model):
        for p in model.params.phase_list:
            for j in model.params.component_list:
                assert (
                    model.prop[1].get_material_density_terms(p, j)
                    is model.prop[1].dens_mol_phase[p]
                )

    @pytest.mark.unit
    def test_get_energy_density_terms(self, model):
        for p in model.params.phase_list:
            assert value(model.prop[1].get_energy_density_terms(p)) == value(
                model.prop[1].dens_mol_phase[p]
                * model.prop[1].energy_internal_mol_phase[p]
            )

    @pytest.mark.unit
    def test_default_material_balance_type(self, model):
        assert (
            model.prop[1].default_material_balance_type()
            is MaterialBalanceType.componentTotal
        )

    @pytest.mark.unit
    def test_default_energy_balance_type(self, model):
        assert (
            model.prop[1].default_energy_balance_type()
            is EnergyBalanceType.enthalpyTotal
        )

    @pytest.mark.unit
    def test_define_state_vars(self, model):
        sv = model.prop[1].define_state_vars()
        assert len(sv) == 3
        for i in sv:
            assert i in ["flow_mol", "enth_mol", "pressure"]

    @pytest.mark.unit
    def test_define_port_members(self, model):
        sv = model.prop[1].define_port_members()
        assert len(sv) == 3
        for i in sv:
            assert i in ["flow_mol", "enth_mol", "pressure"]

    @pytest.mark.unit
    def test_define_display_vars(self, model):
        dv = model.prop[1].define_display_vars()
        assert len(dv) == 6
        for i in dv:
            assert i in [
                "Molar Flow (mol/s)",
                "Mass Flow (kg/s)",
                "T (K)",
                "P (Pa)",
                "Vapor Fraction",
                "Molar Enthalpy (J/mol)",
            ]

    @pytest.mark.unit
    def test_extensive_state_vars(self, model):
        for i in model.prop[1].extensive_set:
            assert i is model.prop[1].flow_mol

    @pytest.mark.unit
    def test_intensive_state_vars(self, model):
        for i in model.prop[1].intensive_set:
            assert i.name in ["prop[1].enth_mol", "prop[1].pressure"]

    @pytest.mark.unit
    def test_model_check(self, model):
        assert model.prop[1].model_check() is None

    @pytest.mark.unit
    def test_data_initialize(self, model):
        assert model.prop[1].initialize() is None

    @pytest.mark.unit
    def test_initialize(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

        flags = model.prop.initialize(hold_state=True)

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].enth_mol.fixed
        assert model.prop[1].pressure.fixed

        model.prop.release_state(flags)

        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

    @pytest.mark.unit
    def test_initialize_w_state(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

        flags = model.prop.initialize(
            hold_state=True,
            state_args={"flow_mol": 20, "enth_mol": 200, "pressure": 2000},
        )

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].enth_mol.fixed
        assert model.prop[1].pressure.fixed

        assert model.prop[1].flow_mol.value == 20
        assert model.prop[1].enth_mol.value == 200
        assert model.prop[1].pressure.value == 2000

        model.prop.release_state(flags)

        assert model.prop[1].flow_mol.value == 20
        assert model.prop[1].enth_mol.value == 200
        assert model.prop[1].pressure.value == 2000

        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

    @pytest.mark.unit
    def test_initialize_w_partial_state_1(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

        flags = model.prop.initialize(hold_state=True,
                                      state_args={"flow_mol": 30})

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].enth_mol.fixed
        assert model.prop[1].pressure.fixed

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].enth_mol.value == 200
        assert model.prop[1].pressure.value == 2000

        model.prop.release_state(flags)

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].enth_mol.value == 200
        assert model.prop[1].pressure.value == 2000

        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

    @pytest.mark.unit
    def test_initialize_w_partial_state_2(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

        flags = model.prop.initialize(hold_state=True,
                                      state_args={"enth_mol": 300})

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].enth_mol.fixed
        assert model.prop[1].pressure.fixed

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].enth_mol.value == 300
        assert model.prop[1].pressure.value == 2000

        model.prop.release_state(flags)

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].enth_mol.value == 300
        assert model.prop[1].pressure.value == 2000

        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

    @pytest.mark.unit
    def test_initialize_w_partial_state_3(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

        flags = model.prop.initialize(hold_state=True,
                                      state_args={"pressure": 3000})

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].enth_mol.fixed
        assert model.prop[1].pressure.fixed

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].enth_mol.value == 300
        assert model.prop[1].pressure.value == 3000

        model.prop.release_state(flags)

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].enth_mol.value == 300
        assert model.prop[1].pressure.value == 3000

        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

    @pytest.mark.unit
    def test_initialize_w_state_ignore_fixed(self, model):
        model.prop[1].flow_mol.fixed = True
        model.prop[1].enth_mol.fixed = True
        model.prop[1].pressure.fixed = True

        flags = model.prop.initialize(
            hold_state=True,
            state_args={"flow_mol": 20, "enth_mol": 200, "pressure": 2000},
        )

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].enth_mol.fixed
        assert model.prop[1].pressure.fixed

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].enth_mol.value == 300
        assert model.prop[1].pressure.value == 3000

        model.prop.release_state(flags)

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].enth_mol.value == 300
        assert model.prop[1].pressure.value == 3000

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].enth_mol.fixed
        assert model.prop[1].pressure.fixed


@pytest.mark.skipif(not prop_available, reason="IAPWS not available")
class TestGPh(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.params = iapws95.Iapws95ParameterBlock(
            default={"phase_presentation": iapws95.PhaseType.G}
        )

        return model

    @pytest.mark.unit
    def test_config(self, model):
        assert model.params.config.phase_presentation == iapws95.PhaseType.G
        assert model.params.config.state_vars == iapws95.StateVars.PH

        assert len(model.params.phase_list) == 1
        for i in model.params.phase_list:
            assert i in ["Vap"]

    @pytest.mark.unit
    def test_build(self, model):
        model.prop = iapws95.Iapws95StateBlock(
            [1], default={"parameters": model.params}
        )

        assert isinstance(model.prop[1].flow_mol, Var)
        assert isinstance(model.prop[1].pressure, Var)
        assert isinstance(model.prop[1].enth_mol, Var)

        assert isinstance(model.prop[1].extensive_set, ComponentSet)
        assert isinstance(model.prop[1].intensive_set, ComponentSet)

    @pytest.mark.unit
    def test_get_material_flow_terms(self, model):
        for p in model.params.phase_list:
            for j in model.params.component_list:
                assert value(
                    model.prop[1].get_material_flow_terms(p, j)) == value(
                    model.prop[1].flow_mol * model.prop[1].phase_frac[p]
                )

    @pytest.mark.unit
    def test_get_enthalpy_flow_terms(self, model):
        for p in model.params.phase_list:
            assert value(model.prop[1].get_enthalpy_flow_terms(p)) == value(
                model.prop[1].enth_mol_phase[p]
                * model.prop[1].phase_frac[p]
                * model.prop[1].flow_mol
            )

    @pytest.mark.unit
    def test_get_material_density_terms(self, model):
        for p in model.params.phase_list:
            for j in model.params.component_list:
                assert (
                    model.prop[1].get_material_density_terms(p, j)
                    is model.prop[1].dens_mol_phase[p]
                )

    @pytest.mark.unit
    def test_get_energy_density_terms(self, model):
        for p in model.params.phase_list:
            assert value(model.prop[1].get_energy_density_terms(p)) == value(
                model.prop[1].dens_mol_phase[p]
                * model.prop[1].energy_internal_mol_phase[p]
            )

    @pytest.mark.unit
    def test_default_material_balance_type(self, model):
        assert (
            model.prop[1].default_material_balance_type()
            is MaterialBalanceType.componentTotal
        )

    @pytest.mark.unit
    def test_default_energy_balance_type(self, model):
        assert (
            model.prop[1].default_energy_balance_type()
            is EnergyBalanceType.enthalpyTotal
        )

    @pytest.mark.unit
    def test_define_state_vars(self, model):
        sv = model.prop[1].define_state_vars()
        assert len(sv) == 3
        for i in sv:
            assert i in ["flow_mol", "enth_mol", "pressure"]

    @pytest.mark.unit
    def test_define_port_members(self, model):
        sv = model.prop[1].define_port_members()
        assert len(sv) == 3
        for i in sv:
            assert i in ["flow_mol", "enth_mol", "pressure"]

    @pytest.mark.unit
    def test_define_display_vars(self, model):
        dv = model.prop[1].define_display_vars()
        assert len(dv) == 6
        for i in dv:
            assert i in [
                "Molar Flow (mol/s)",
                "Mass Flow (kg/s)",
                "T (K)",
                "P (Pa)",
                "Vapor Fraction",
                "Molar Enthalpy (J/mol)",
            ]

    @pytest.mark.unit
    def test_extensive_state_vars(self, model):
        for i in model.prop[1].extensive_set:
            assert i is model.prop[1].flow_mol

    @pytest.mark.unit
    def test_intensive_state_vars(self, model):
        for i in model.prop[1].intensive_set:
            assert i.name in ["prop[1].enth_mol", "prop[1].pressure"]

    @pytest.mark.unit
    def test_model_check(self, model):
        assert model.prop[1].model_check() is None

    @pytest.mark.unit
    def test_data_initialize(self, model):
        assert model.prop[1].initialize() is None

    @pytest.mark.unit
    def test_initialize(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

        flags = model.prop.initialize(hold_state=True)

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].enth_mol.fixed
        assert model.prop[1].pressure.fixed

        model.prop.release_state(flags)

        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

    @pytest.mark.unit
    def test_initialize_w_state(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

        flags = model.prop.initialize(
            hold_state=True,
            state_args={"flow_mol": 20, "enth_mol": 200, "pressure": 2000},
        )

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].enth_mol.fixed
        assert model.prop[1].pressure.fixed

        assert model.prop[1].flow_mol.value == 20
        assert model.prop[1].enth_mol.value == 200
        assert model.prop[1].pressure.value == 2000

        model.prop.release_state(flags)

        assert model.prop[1].flow_mol.value == 20
        assert model.prop[1].enth_mol.value == 200
        assert model.prop[1].pressure.value == 2000

        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

    @pytest.mark.unit
    def test_initialize_w_partial_state_1(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

        flags = model.prop.initialize(hold_state=True,
                                      state_args={"flow_mol": 30})

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].enth_mol.fixed
        assert model.prop[1].pressure.fixed

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].enth_mol.value == 200
        assert model.prop[1].pressure.value == 2000

        model.prop.release_state(flags)

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].enth_mol.value == 200
        assert model.prop[1].pressure.value == 2000

        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

    @pytest.mark.unit
    def test_initialize_w_partial_state_2(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

        flags = model.prop.initialize(hold_state=True,
                                      state_args={"enth_mol": 300})

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].enth_mol.fixed
        assert model.prop[1].pressure.fixed

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].enth_mol.value == 300
        assert model.prop[1].pressure.value == 2000

        model.prop.release_state(flags)

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].enth_mol.value == 300
        assert model.prop[1].pressure.value == 2000

        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

    @pytest.mark.unit
    def test_initialize_w_partial_state_3(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

        flags = model.prop.initialize(hold_state=True,
                                      state_args={"pressure": 3000})

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].enth_mol.fixed
        assert model.prop[1].pressure.fixed

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].enth_mol.value == 300
        assert model.prop[1].pressure.value == 3000

        model.prop.release_state(flags)

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].enth_mol.value == 300
        assert model.prop[1].pressure.value == 3000

        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

    @pytest.mark.unit
    def test_initialize_w_state_ignore_fixed(self, model):
        model.prop[1].flow_mol.fixed = True
        model.prop[1].enth_mol.fixed = True
        model.prop[1].pressure.fixed = True

        flags = model.prop.initialize(
            hold_state=True,
            state_args={"flow_mol": 20, "enth_mol": 200, "pressure": 2000},
        )

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].enth_mol.fixed
        assert model.prop[1].pressure.fixed

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].enth_mol.value == 300
        assert model.prop[1].pressure.value == 3000

        model.prop.release_state(flags)

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].enth_mol.value == 300
        assert model.prop[1].pressure.value == 3000

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].enth_mol.fixed
        assert model.prop[1].pressure.fixed


@pytest.mark.skipif(not prop_available, reason="IAPWS not available")
class TestMixTpx(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.params = iapws95.Iapws95ParameterBlock(
            default={"state_vars": iapws95.StateVars.TPX}
        )

        return model

    @pytest.mark.unit
    def test_config(self, model):
        assert model.params.config.phase_presentation == iapws95.PhaseType.MIX
        assert model.params.config.state_vars == iapws95.StateVars.TPX

        assert len(model.params.phase_list) == 1
        for i in model.params.phase_list:
            assert i in ["Mix"]

    @pytest.mark.unit
    def test_build(self, model):
        model.prop = iapws95.Iapws95StateBlock(
            [1], default={"parameters": model.params}
        )

        assert isinstance(model.prop[1].flow_mol, Var)
        assert isinstance(model.prop[1].pressure, Var)
        assert isinstance(model.prop[1].temperature, Var)
        assert isinstance(model.prop[1].vapor_frac, Var)

        assert isinstance(model.prop[1].extensive_set, ComponentSet)
        assert isinstance(model.prop[1].intensive_set, ComponentSet)

    @pytest.mark.unit
    def test_get_material_flow_terms(self, model):
        for p in model.params.phase_list:
            for j in model.params.component_list:
                assert value(
                    model.prop[1].get_material_flow_terms(p, j)) == value(
                    model.prop[1].flow_mol
                )

    @pytest.mark.unit
    def test_get_enthalpy_flow_terms(self, model):
        for p in model.params.phase_list:
            assert value(
                model.prop[1].get_enthalpy_flow_terms(p)) == value(
                model.prop[1].enth_mol * model.prop[1].flow_mol
            )

    @pytest.mark.unit
    def test_get_material_density_terms(self, model):
        for p in model.params.phase_list:
            for j in model.params.component_list:
                assert (
                    model.prop[1].get_material_density_terms(p, j)
                    is model.prop[1].dens_mol
                )

    @pytest.mark.unit
    def test_get_energy_density_terms(self, model):
        for p in model.params.phase_list:
            assert value(
                model.prop[1].get_energy_density_terms(p)) == value(
                model.prop[1].dens_mol * model.prop[1].energy_internal_mol
            )

    @pytest.mark.unit
    def test_default_material_balance_type(self, model):
        assert (
            model.prop[1].default_material_balance_type()
            is MaterialBalanceType.componentTotal
        )

    @pytest.mark.unit
    def test_default_energy_balance_type(self, model):
        assert (
            model.prop[1].default_energy_balance_type()
            is EnergyBalanceType.enthalpyTotal
        )

    @pytest.mark.unit
    def test_define_state_vars(self, model):
        sv = model.prop[1].define_state_vars()
        assert len(sv) == 4
        for i in sv:
            assert i in ["flow_mol", "temperature", "pressure", "vapor_frac"]

    @pytest.mark.unit
    def test_define_port_members(self, model):
        sv = model.prop[1].define_port_members()
        assert len(sv) == 4
        for i in sv:
            assert i in ["flow_mol", "temperature", "pressure", "vapor_frac"]

    @pytest.mark.unit
    def test_define_display_vars(self, model):
        dv = model.prop[1].define_display_vars()
        assert len(dv) == 6
        for i in dv:
            assert i in [
                "Molar Flow (mol/s)",
                "Mass Flow (kg/s)",
                "T (K)",
                "P (Pa)",
                "Vapor Fraction",
                "Molar Enthalpy (J/mol)",
            ]

    @pytest.mark.unit
    def test_extensive_state_vars(self, model):
        for i in model.prop[1].extensive_set:
            assert i is model.prop[1].flow_mol

    @pytest.mark.unit
    def test_intensive_state_vars(self, model):
        for i in model.prop[1].intensive_set:
            assert i.name in [
                "prop[1].temperature",
                "prop[1].pressure",
                "prop[1].vapor_frac",
            ]

    @pytest.mark.unit
    def test_model_check(self, model):
        assert model.prop[1].model_check() is None

    @pytest.mark.unit
    def test_data_initialize(self, model):
        assert model.prop[1].initialize() is None

    @pytest.mark.unit
    def test_initialize(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

        flags = model.prop.initialize(hold_state=True)

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].temperature.fixed
        assert model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

        model.prop.release_state(flags)

        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

    @pytest.mark.unit
    def test_initialize_w_state(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

        flags = model.prop.initialize(
            hold_state=True,
            state_args={
                "flow_mol": 20,
                "temperature": 200,
                "pressure": 2000,
                "vapor_frac": 0.2,
            },
        )

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].temperature.fixed
        assert model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

        assert model.prop[1].flow_mol.value == 20
        assert model.prop[1].temperature.value == 200
        assert model.prop[1].pressure.value == 2000
        assert model.prop[1].vapor_frac.value == 0.0

        model.prop.release_state(flags)

        assert model.prop[1].flow_mol.value == 20
        assert model.prop[1].temperature.value == 200
        assert model.prop[1].pressure.value == 2000
        assert model.prop[1].vapor_frac.value == 0.0

        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

    @pytest.mark.unit
    def test_initialize_w_partial_state_1(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

        flags = model.prop.initialize(hold_state=True,
                                      state_args={"flow_mol": 30})

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].temperature.fixed
        assert model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].temperature.value == 200
        assert model.prop[1].pressure.value == 2000
        assert model.prop[1].vapor_frac.value == 0.0

        model.prop.release_state(flags)

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].temperature.value == 200
        assert model.prop[1].pressure.value == 2000
        assert model.prop[1].vapor_frac.value == 0.0

        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

    @pytest.mark.unit
    def test_initialize_w_partial_state_2(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

        flags = model.prop.initialize(hold_state=True,
                                      state_args={"temperature": 300})

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].temperature.fixed
        assert model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].temperature.value == 300
        assert model.prop[1].pressure.value == 2000
        assert model.prop[1].vapor_frac.value == 0.0

        model.prop.release_state(flags)

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].temperature.value == 300
        assert model.prop[1].pressure.value == 2000
        assert model.prop[1].vapor_frac.value == 0.0

        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

    @pytest.mark.unit
    def test_initialize_w_partial_state_3(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

        flags = model.prop.initialize(hold_state=True,
                                      state_args={"pressure": 3000})

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].temperature.fixed
        assert model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].temperature.value == 300
        assert model.prop[1].pressure.value == 3000
        assert model.prop[1].vapor_frac.value == 0.0

        model.prop.release_state(flags)

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].temperature.value == 300
        assert model.prop[1].pressure.value == 3000
        assert model.prop[1].vapor_frac.value == 0.0

        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

    @pytest.mark.unit
    def test_initialize_w_state_ignore_fixed(self, model):
        model.prop[1].flow_mol.fixed = True
        model.prop[1].temperature.fixed = True
        model.prop[1].pressure.fixed = True

        flags = model.prop.initialize(
            hold_state=True,
            state_args={"flow_mol": 20, "temperature": 200, "pressure": 2000},
        )

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].temperature.fixed
        assert model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].temperature.value == 300
        assert model.prop[1].pressure.value == 3000
        assert model.prop[1].vapor_frac.value == 0.0

        model.prop.release_state(flags)

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].temperature.value == 300
        assert model.prop[1].pressure.value == 3000
        assert model.prop[1].vapor_frac.value == 0.0

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].temperature.fixed
        assert model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed


@pytest.mark.skipif(not prop_available, reason="IAPWS not available")
class TestLgTpx(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.params = iapws95.Iapws95ParameterBlock(
            default={
                "phase_presentation": iapws95.PhaseType.LG,
                "state_vars": iapws95.StateVars.TPX,
            }
        )

        return model

    @pytest.mark.unit
    def test_config(self, model):
        assert model.params.config.phase_presentation == iapws95.PhaseType.LG
        assert model.params.config.state_vars == iapws95.StateVars.TPX

        assert len(model.params.phase_list) == 2
        for i in model.params.phase_list:
            assert i in ["Liq", "Vap"]

    @pytest.mark.unit
    def test_build(self, model):
        model.prop = iapws95.Iapws95StateBlock(
            [1], default={"parameters": model.params}
        )

        assert isinstance(model.prop[1].flow_mol, Var)
        assert isinstance(model.prop[1].pressure, Var)
        assert isinstance(model.prop[1].temperature, Var)
        assert isinstance(model.prop[1].vapor_frac, Var)

        assert isinstance(model.prop[1].extensive_set, ComponentSet)
        assert isinstance(model.prop[1].intensive_set, ComponentSet)

    @pytest.mark.unit
    def test_get_material_flow_terms(self, model):
        for p in model.params.phase_list:
            for j in model.params.component_list:
                assert value(
                    model.prop[1].get_material_flow_terms(p, j)) == value(
                    model.prop[1].flow_mol * model.prop[1].phase_frac[p]
                )

    @pytest.mark.unit
    def test_get_enthalpy_flow_terms(self, model):
        for p in model.params.phase_list:
            assert value(model.prop[1].get_enthalpy_flow_terms(p)) == value(
                model.prop[1].enth_mol_phase[p]
                * model.prop[1].phase_frac[p]
                * model.prop[1].flow_mol
            )

    @pytest.mark.unit
    def test_get_material_density_terms(self, model):
        for p in model.params.phase_list:
            for j in model.params.component_list:
                assert (
                    model.prop[1].get_material_density_terms(p, j)
                    is model.prop[1].dens_mol_phase[p]
                )

    @pytest.mark.unit
    def test_get_energy_density_terms(self, model):
        for p in model.params.phase_list:
            assert value(model.prop[1].get_energy_density_terms(p)) == value(
                model.prop[1].dens_mol_phase[p]
                * model.prop[1].energy_internal_mol_phase[p]
            )

    @pytest.mark.unit
    def test_default_material_balance_type(self, model):
        assert (
            model.prop[1].default_material_balance_type()
            is MaterialBalanceType.componentTotal
        )

    @pytest.mark.unit
    def test_default_energy_balance_type(self, model):
        assert (
            model.prop[1].default_energy_balance_type()
            is EnergyBalanceType.enthalpyTotal
        )

    @pytest.mark.unit
    def test_define_state_vars(self, model):
        sv = model.prop[1].define_state_vars()
        assert len(sv) == 4
        for i in sv:
            assert i in ["flow_mol", "temperature", "pressure", "vapor_frac"]

    @pytest.mark.unit
    def test_define_port_members(self, model):
        sv = model.prop[1].define_port_members()
        assert len(sv) == 4
        for i in sv:
            assert i in ["flow_mol", "temperature", "pressure", "vapor_frac"]

    @pytest.mark.unit
    def test_define_display_vars(self, model):
        dv = model.prop[1].define_display_vars()
        assert len(dv) == 6
        for i in dv:
            assert i in [
                "Molar Flow (mol/s)",
                "Mass Flow (kg/s)",
                "T (K)",
                "P (Pa)",
                "Vapor Fraction",
                "Molar Enthalpy (J/mol)",
            ]

    @pytest.mark.unit
    def test_extensive_state_vars(self, model):
        for i in model.prop[1].extensive_set:
            assert i is model.prop[1].flow_mol

    @pytest.mark.unit
    def test_intensive_state_vars(self, model):
        for i in model.prop[1].intensive_set:
            assert i.name in [
                "prop[1].temperature",
                "prop[1].pressure",
                "prop[1].vapor_frac",
            ]

    @pytest.mark.unit
    def test_model_check(self, model):
        assert model.prop[1].model_check() is None

    @pytest.mark.unit
    def test_data_initialize(self, model):
        assert model.prop[1].initialize() is None

    @pytest.mark.unit
    def test_initialize(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert not model.prop[1].vapor_frac.fixed

        flags = model.prop.initialize(hold_state=True)

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].temperature.fixed
        assert model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

        model.prop.release_state(flags)

        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert not model.prop[1].vapor_frac.fixed

    @pytest.mark.unit
    def test_initialize_w_state(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert not model.prop[1].vapor_frac.fixed

        flags = model.prop.initialize(
            hold_state=True,
            state_args={"flow_mol": 20, "temperature": 200, "pressure": 2000},
        )

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].temperature.fixed
        assert model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

        assert model.prop[1].flow_mol.value == 20
        assert model.prop[1].temperature.value == 200
        assert model.prop[1].pressure.value == 2000

        model.prop.release_state(flags)

        assert model.prop[1].flow_mol.value == 20
        assert model.prop[1].temperature.value == 200
        assert model.prop[1].pressure.value == 2000

        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert not model.prop[1].vapor_frac.fixed

    @pytest.mark.unit
    def test_initialize_w_partial_state_1(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert not model.prop[1].vapor_frac.fixed

        flags = model.prop.initialize(hold_state=True,
                                      state_args={"flow_mol": 30})

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].temperature.fixed
        assert model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].temperature.value == 200
        assert model.prop[1].pressure.value == 2000

        model.prop.release_state(flags)

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].temperature.value == 200
        assert model.prop[1].pressure.value == 2000

        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert not model.prop[1].vapor_frac.fixed

    @pytest.mark.unit
    def test_initialize_w_partial_state_2(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert not model.prop[1].vapor_frac.fixed

        flags = model.prop.initialize(hold_state=True,
                                      state_args={"temperature": 300})

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].temperature.fixed
        assert model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].temperature.value == 300
        assert model.prop[1].pressure.value == 2000

        model.prop.release_state(flags)

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].temperature.value == 300
        assert model.prop[1].pressure.value == 2000

        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert not model.prop[1].vapor_frac.fixed

    @pytest.mark.unit
    def test_initialize_w_partial_state_3(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert not model.prop[1].vapor_frac.fixed

        flags = model.prop.initialize(hold_state=True,
                                      state_args={"pressure": 3000})

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].temperature.fixed
        assert model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].temperature.value == 300
        assert model.prop[1].pressure.value == 3000

        model.prop.release_state(flags)

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].temperature.value == 300
        assert model.prop[1].pressure.value == 3000

        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert not model.prop[1].vapor_frac.fixed

    @pytest.mark.unit
    def test_initialize_w_state_ignore_fixed(self, model):
        model.prop[1].flow_mol.fixed = True
        model.prop[1].temperature.fixed = True
        model.prop[1].pressure.fixed = True

        flags = model.prop.initialize(
            hold_state=True,
            state_args={"flow_mol": 20, "temperature": 200, "pressure": 2000},
        )

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].temperature.fixed
        assert model.prop[1].pressure.fixed

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].temperature.value == 300
        assert model.prop[1].pressure.value == 3000

        model.prop.release_state(flags)

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].temperature.value == 300
        assert model.prop[1].pressure.value == 3000

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].temperature.fixed
        assert model.prop[1].pressure.fixed


@pytest.mark.skipif(not prop_available, reason="IAPWS not available")
class TestLTpx(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.params = iapws95.Iapws95ParameterBlock(
            default={
                "phase_presentation": iapws95.PhaseType.L,
                "state_vars": iapws95.StateVars.TPX,
            }
        )

        return model

    @pytest.mark.unit
    def test_config(self, model):
        assert model.params.config.phase_presentation == iapws95.PhaseType.L
        assert model.params.config.state_vars == iapws95.StateVars.TPX

        assert len(model.params.phase_list) == 1
        for i in model.params.phase_list:
            assert i in ["Liq"]

    @pytest.mark.unit
    def test_build(self, model):
        model.prop = iapws95.Iapws95StateBlock(
            [1], default={"parameters": model.params}
        )

        assert isinstance(model.prop[1].flow_mol, Var)
        assert isinstance(model.prop[1].pressure, Var)
        assert isinstance(model.prop[1].temperature, Var)
        assert isinstance(model.prop[1].vapor_frac, Var)

        assert isinstance(model.prop[1].extensive_set, ComponentSet)
        assert isinstance(model.prop[1].intensive_set, ComponentSet)

    @pytest.mark.unit
    def test_get_material_flow_terms(self, model):
        for p in model.params.phase_list:
            for j in model.params.component_list:
                assert value(
                    model.prop[1].get_material_flow_terms(p, j)) == value(
                    model.prop[1].flow_mol * model.prop[1].phase_frac[p]
                )

    @pytest.mark.unit
    def test_get_enthalpy_flow_terms(self, model):
        for p in model.params.phase_list:
            assert value(model.prop[1].get_enthalpy_flow_terms(p)) == value(
                model.prop[1].enth_mol_phase[p]
                * model.prop[1].phase_frac[p]
                * model.prop[1].flow_mol
            )

    @pytest.mark.unit
    def test_get_material_density_terms(self, model):
        for p in model.params.phase_list:
            for j in model.params.component_list:
                assert (
                    model.prop[1].get_material_density_terms(p, j)
                    is model.prop[1].dens_mol_phase[p]
                )

    @pytest.mark.unit
    def test_get_energy_density_terms(self, model):
        for p in model.params.phase_list:
            assert value(model.prop[1].get_energy_density_terms(p)) == value(
                model.prop[1].dens_mol_phase[p]
                * model.prop[1].energy_internal_mol_phase[p]
            )

    @pytest.mark.unit
    def test_default_material_balance_type(self, model):
        assert (
            model.prop[1].default_material_balance_type()
            is MaterialBalanceType.componentTotal
        )

    @pytest.mark.unit
    def test_default_energy_balance_type(self, model):
        assert (
            model.prop[1].default_energy_balance_type()
            is EnergyBalanceType.enthalpyTotal
        )

    @pytest.mark.unit
    def test_define_state_vars(self, model):
        sv = model.prop[1].define_state_vars()
        assert len(sv) == 3
        for i in sv:
            assert i in ["flow_mol", "temperature", "pressure"]

    @pytest.mark.unit
    def test_define_port_members(self, model):
        sv = model.prop[1].define_port_members()
        assert len(sv) == 3
        for i in sv:
            assert i in ["flow_mol", "temperature", "pressure"]

    @pytest.mark.unit
    def test_define_display_vars(self, model):
        dv = model.prop[1].define_display_vars()
        assert len(dv) == 6
        for i in dv:
            assert i in [
                "Molar Flow (mol/s)",
                "Mass Flow (kg/s)",
                "T (K)",
                "P (Pa)",
                "Vapor Fraction",
                "Molar Enthalpy (J/mol)",
            ]

    @pytest.mark.unit
    def test_extensive_state_vars(self, model):
        for i in model.prop[1].extensive_set:
            assert i is model.prop[1].flow_mol

    @pytest.mark.unit
    def test_intensive_state_vars(self, model):
        for i in model.prop[1].intensive_set:
            assert i.name in [
                "prop[1].temperature",
                "prop[1].pressure",
                "prop[1].vapor_frac",
            ]

    @pytest.mark.unit
    def test_model_check(self, model):
        assert model.prop[1].model_check() is None

    @pytest.mark.unit
    def test_data_initialize(self, model):
        assert model.prop[1].initialize() is None

    @pytest.mark.unit
    def test_initialize(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

        flags = model.prop.initialize(hold_state=True)

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].temperature.fixed
        assert model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

        model.prop.release_state(flags)

        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

    @pytest.mark.unit
    def test_initialize_w_state(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

        flags = model.prop.initialize(
            hold_state=True,
            state_args={"flow_mol": 20, "temperature": 200, "pressure": 2000},
        )

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].temperature.fixed
        assert model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

        assert model.prop[1].flow_mol.value == 20
        assert model.prop[1].temperature.value == 200
        assert model.prop[1].pressure.value == 2000

        model.prop.release_state(flags)

        assert model.prop[1].flow_mol.value == 20
        assert model.prop[1].temperature.value == 200
        assert model.prop[1].pressure.value == 2000

        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

    @pytest.mark.unit
    def test_initialize_w_partial_state_1(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

        flags = model.prop.initialize(hold_state=True,
                                      state_args={"flow_mol": 30})

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].temperature.fixed
        assert model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].temperature.value == 200
        assert model.prop[1].pressure.value == 2000

        model.prop.release_state(flags)

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].temperature.value == 200
        assert model.prop[1].pressure.value == 2000

        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

    @pytest.mark.unit
    def test_initialize_w_partial_state_2(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

        flags = model.prop.initialize(hold_state=True,
                                      state_args={"temperature": 300})

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].temperature.fixed
        assert model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].temperature.value == 300
        assert model.prop[1].pressure.value == 2000

        model.prop.release_state(flags)

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].temperature.value == 300
        assert model.prop[1].pressure.value == 2000

        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

    @pytest.mark.unit
    def test_initialize_w_partial_state_3(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

        flags = model.prop.initialize(hold_state=True,
                                      state_args={"pressure": 3000})

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].temperature.fixed
        assert model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].temperature.value == 300
        assert model.prop[1].pressure.value == 3000

        model.prop.release_state(flags)

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].temperature.value == 300
        assert model.prop[1].pressure.value == 3000

        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

    @pytest.mark.unit
    def test_initialize_w_state_ignore_fixed(self, model):
        model.prop[1].flow_mol.fixed = True
        model.prop[1].temperature.fixed = True
        model.prop[1].pressure.fixed = True

        flags = model.prop.initialize(
            hold_state=True,
            state_args={"flow_mol": 20, "temperature": 200, "pressure": 2000},
        )

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].temperature.fixed
        assert model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].temperature.value == 300
        assert model.prop[1].pressure.value == 3000

        model.prop.release_state(flags)

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].temperature.value == 300
        assert model.prop[1].pressure.value == 3000

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].temperature.fixed
        assert model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed


@pytest.mark.skipif(not prop_available, reason="IAPWS not available")
class TestGTpx(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.params = iapws95.Iapws95ParameterBlock(
            default={
                "phase_presentation": iapws95.PhaseType.G,
                "state_vars": iapws95.StateVars.TPX,
            }
        )
        return model

    @pytest.mark.unit
    def test_config(self, model):
        assert model.params.config.phase_presentation == iapws95.PhaseType.G
        assert model.params.config.state_vars == iapws95.StateVars.TPX

        assert len(model.params.phase_list) == 1
        for i in model.params.phase_list:
            assert i in ["Vap"]

    @pytest.mark.unit
    def test_build(self, model):
        model.prop = iapws95.Iapws95StateBlock(
            [1], default={"parameters": model.params}
        )

        assert isinstance(model.prop[1].flow_mol, Var)
        assert isinstance(model.prop[1].pressure, Var)
        assert isinstance(model.prop[1].temperature, Var)
        assert isinstance(model.prop[1].vapor_frac, Var)

        assert isinstance(model.prop[1].extensive_set, ComponentSet)
        assert isinstance(model.prop[1].intensive_set, ComponentSet)

    @pytest.mark.unit
    def test_get_material_flow_terms(self, model):
        for p in model.params.phase_list:
            for j in model.params.component_list:
                assert value(
                    model.prop[1].get_material_flow_terms(p, j)) == value(
                    model.prop[1].flow_mol * model.prop[1].phase_frac[p]
                )

    @pytest.mark.unit
    def test_get_enthalpy_flow_terms(self, model):
        for p in model.params.phase_list:
            assert value(model.prop[1].get_enthalpy_flow_terms(p)) == value(
                model.prop[1].enth_mol_phase[p]
                * model.prop[1].phase_frac[p]
                * model.prop[1].flow_mol
            )

    @pytest.mark.unit
    def test_get_material_density_terms(self, model):
        for p in model.params.phase_list:
            for j in model.params.component_list:
                assert (
                    model.prop[1].get_material_density_terms(p, j)
                    is model.prop[1].dens_mol_phase[p]
                )

    @pytest.mark.unit
    def test_get_energy_density_terms(self, model):
        for p in model.params.phase_list:
            assert value(model.prop[1].get_energy_density_terms(p)) == value(
                model.prop[1].dens_mol_phase[p]
                * model.prop[1].energy_internal_mol_phase[p]
            )

    @pytest.mark.unit
    def test_default_material_balance_type(self, model):
        assert (
            model.prop[1].default_material_balance_type()
            is MaterialBalanceType.componentTotal
        )

    @pytest.mark.unit
    def test_default_energy_balance_type(self, model):
        assert (
            model.prop[1].default_energy_balance_type()
            is EnergyBalanceType.enthalpyTotal
        )

    @pytest.mark.unit
    def test_define_state_vars(self, model):
        sv = model.prop[1].define_state_vars()
        assert len(sv) == 3
        for i in sv:
            assert i in ["flow_mol", "temperature", "pressure"]

    @pytest.mark.unit
    def test_define_port_members(self, model):
        sv = model.prop[1].define_port_members()
        assert len(sv) == 3
        for i in sv:
            assert i in ["flow_mol", "temperature", "pressure"]

    @pytest.mark.unit
    def test_define_display_vars(self, model):
        dv = model.prop[1].define_display_vars()
        assert len(dv) == 6
        for i in dv:
            assert i in [
                "Molar Flow (mol/s)",
                "Mass Flow (kg/s)",
                "T (K)",
                "P (Pa)",
                "Vapor Fraction",
                "Molar Enthalpy (J/mol)",
            ]

    @pytest.mark.unit
    def test_extensive_state_vars(self, model):
        for i in model.prop[1].extensive_set:
            assert i is model.prop[1].flow_mol

    @pytest.mark.unit
    def test_intensive_state_vars(self, model):
        for i in model.prop[1].intensive_set:
            assert i.name in [
                "prop[1].temperature",
                "prop[1].pressure",
                "prop[1].vapor_frac",
            ]

    @pytest.mark.unit
    def test_model_check(self, model):
        assert model.prop[1].model_check() is None

    @pytest.mark.unit
    def test_data_initialize(self, model):
        assert model.prop[1].initialize() is None

    @pytest.mark.unit
    def test_initialize(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

        flags = model.prop.initialize(hold_state=True)

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].temperature.fixed
        assert model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

        model.prop.release_state(flags)

        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

    @pytest.mark.unit
    def test_initialize_w_state(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

        flags = model.prop.initialize(
            hold_state=True,
            state_args={"flow_mol": 20, "temperature": 200, "pressure": 2000},
        )

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].temperature.fixed
        assert model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

        assert model.prop[1].flow_mol.value == 20
        assert model.prop[1].temperature.value == 200
        assert model.prop[1].pressure.value == 2000

        model.prop.release_state(flags)

        assert model.prop[1].flow_mol.value == 20
        assert model.prop[1].temperature.value == 200
        assert model.prop[1].pressure.value == 2000

        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

    @pytest.mark.unit
    def test_initialize_w_partial_state_1(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

        flags = model.prop.initialize(hold_state=True,
                                      state_args={"flow_mol": 30})

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].temperature.fixed
        assert model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].temperature.value == 200
        assert model.prop[1].pressure.value == 2000

        model.prop.release_state(flags)

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].temperature.value == 200
        assert model.prop[1].pressure.value == 2000

        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

    @pytest.mark.unit
    def test_initialize_w_partial_state_2(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

        flags = model.prop.initialize(hold_state=True,
                                      state_args={"temperature": 300})

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].temperature.fixed
        assert model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].temperature.value == 300
        assert model.prop[1].pressure.value == 2000

        model.prop.release_state(flags)

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].temperature.value == 300
        assert model.prop[1].pressure.value == 2000

        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

    @pytest.mark.unit
    def test_initialize_w_partial_state_3(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

        flags = model.prop.initialize(hold_state=True,
                                      state_args={"pressure": 3000})

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].temperature.fixed
        assert model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].temperature.value == 300
        assert model.prop[1].pressure.value == 3000

        model.prop.release_state(flags)

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].temperature.value == 300
        assert model.prop[1].pressure.value == 3000

        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

    @pytest.mark.unit
    def test_initialize_w_state_ignore_fixed(self, model):
        model.prop[1].flow_mol.fixed = True
        model.prop[1].temperature.fixed = True
        model.prop[1].pressure.fixed = True

        flags = model.prop.initialize(
            hold_state=True,
            state_args={"flow_mol": 20, "temperature": 200, "pressure": 2000},
        )

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].temperature.fixed
        assert model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].temperature.value == 300
        assert model.prop[1].pressure.value == 3000

        model.prop.release_state(flags)

        assert model.prop[1].flow_mol.value == 30
        assert model.prop[1].temperature.value == 300
        assert model.prop[1].pressure.value == 3000

        assert model.prop[1].flow_mol.fixed
        assert model.prop[1].temperature.fixed
        assert model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed
