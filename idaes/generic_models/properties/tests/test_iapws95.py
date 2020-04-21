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

import pytest
from pyomo.environ import ConcreteModel, value, Var
from pyomo.core.kernel.component_set import ComponentSet
from pyomo.common.fileutils import this_file_dir
from idaes.generic_models.properties import iapws95
import csv
import os

from idaes.core import MaterialBalanceType, EnergyBalanceType
from idaes.core.util.exceptions import ConfigurationError

# Set module level pyest marker
pytestmark = pytest.mark.iapws
prop_available = iapws95.iapws95_available()


# -----------------------------------------------------------------------------
# Test Enums and common functions
def test_htpx_invalid_args():
    with pytest.raises(ConfigurationError):
        iapws95.htpx(300, P=101325, x=0.5)

    with pytest.raises(ConfigurationError):
        iapws95.htpx(100, x=0.5)

    with pytest.raises(ConfigurationError):
        iapws95.htpx(5e3, x=0.5)

    with pytest.raises(ConfigurationError):
        iapws95.htpx(300, P=0)

    with pytest.raises(ConfigurationError):
        iapws95.htpx(300, P=1e10)

    with pytest.raises(ConfigurationError):
        iapws95.htpx(300, x=-1)

    with pytest.raises(ConfigurationError):
        iapws95.htpx(300, x=2)


@pytest.mark.skipif(not prop_available, reason="IAPWS not available")
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
    assert iapws95.htpx(300, P=101325) == pytest.approx(112143 * mw + offset, 1e-5)

    # Saturated liquid
    assert iapws95.htpx(500, x=0) == pytest.approx(974919 * mw + offset, 1e-5)

    # Wet steam
    assert iapws95.htpx(550, x=0.5) == pytest.approx(2.00138e06 * mw + offset, 1e-5)

    # Saturated vapor
    assert iapws95.htpx(600, x=1) == pytest.approx(2.67730e06 * mw + offset, 1e-5)

    # Superheated steam
    assert iapws95.htpx(400, P=101325) == pytest.approx(2.72979e06 * mw + offset, 1e-5)


def test_PhaseType():
    assert len(iapws95.PhaseType) == 4

    # Check that expected values do not raise Exceptions
    iapws95.PhaseType.MIX
    iapws95.PhaseType.LG
    iapws95.PhaseType.L
    iapws95.PhaseType.G


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
        model.params = iapws95.Iapws95ParameterBlock()

        return model

    def test_config(self, model):
        assert model.params.config.phase_presentation == iapws95.PhaseType.MIX
        assert model.params.config.state_vars == iapws95.StateVars.PH

        assert len(model.params.phase_list) == 1
        for i in model.params.phase_list:
            assert i in ["Mix"]

    def test_build(self, model):
        model.prop = iapws95.Iapws95StateBlock(
            [1], default={"parameters": model.params}
        )

        assert isinstance(model.prop[1].flow_mol, Var)
        assert isinstance(model.prop[1].pressure, Var)
        assert isinstance(model.prop[1].enth_mol, Var)

        assert isinstance(model.prop[1].extensive_set, ComponentSet)
        assert isinstance(model.prop[1].intensive_set, ComponentSet)

    def test_get_material_flow_terms(self, model):
        for p in model.params.phase_list:
            for j in model.params.component_list:
                assert value(model.prop[1].get_material_flow_terms(p, j)) == value(
                    model.prop[1].flow_mol
                )

    def test_get_enthalpy_flow_terms(self, model):
        for p in model.params.phase_list:
            assert (
                model.prop[1].get_enthalpy_flow_terms(p)
                == model.prop[1].enth_mol * model.prop[1].flow_mol
            )

    def test_get_material_density_terms(self, model):
        for p in model.params.phase_list:
            for j in model.params.component_list:
                assert (
                    model.prop[1].get_material_density_terms(p, j)
                    is model.prop[1].dens_mol
                )

    def test_get_energy_density_terms(self, model):
        for p in model.params.phase_list:
            assert (
                model.prop[1].get_energy_density_terms(p)
                == model.prop[1].dens_mol * model.prop[1].energy_internal_mol
            )

    def test_default_material_balance_type(self, model):
        assert (
            model.prop[1].default_material_balance_type()
            is MaterialBalanceType.componentTotal
        )

    def test_default_energy_balance_type(self, model):
        assert (
            model.prop[1].default_energy_balance_type()
            is EnergyBalanceType.enthalpyTotal
        )

    def test_define_state_vars(self, model):
        sv = model.prop[1].define_state_vars()
        assert len(sv) == 3
        for i in sv:
            assert i in ["flow_mol", "enth_mol", "pressure"]

    def test_define_port_members(self, model):
        sv = model.prop[1].define_port_members()
        assert len(sv) == 3
        for i in sv:
            assert i in ["flow_mol", "enth_mol", "pressure"]

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

    def test_extensive_state_vars(self, model):
        for i in model.prop[1].extensive_set:
            assert i is model.prop[1].flow_mol

    def test_intensive_state_vars(self, model):
        for i in model.prop[1].intensive_set:
            assert i in [model.prop[1].enth_mol, model.prop[1].pressure]

    def test_model_check(self, model):
        assert model.prop[1].model_check() is None

    def test_data_initialize(self, model):
        assert model.prop[1].initialize() is None

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

    def test_initialize_w_partial_state_1(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

        flags = model.prop.initialize(hold_state=True, state_args={"flow_mol": 30})

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

    def test_initialize_w_partial_state_2(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

        flags = model.prop.initialize(hold_state=True, state_args={"enth_mol": 300})

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

    def test_initialize_w_partial_state_3(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

        flags = model.prop.initialize(hold_state=True, state_args={"pressure": 3000})

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

    def test_config(self, model):
        assert model.params.config.phase_presentation == iapws95.PhaseType.LG
        assert model.params.config.state_vars == iapws95.StateVars.PH

        assert len(model.params.phase_list) == 2
        for i in model.params.phase_list:
            assert i in ["Liq", "Vap"]

    def test_build(self, model):
        model.prop = iapws95.Iapws95StateBlock(
            [1], default={"parameters": model.params}
        )

        assert isinstance(model.prop[1].flow_mol, Var)
        assert isinstance(model.prop[1].pressure, Var)
        assert isinstance(model.prop[1].enth_mol, Var)

        assert isinstance(model.prop[1].extensive_set, ComponentSet)
        assert isinstance(model.prop[1].intensive_set, ComponentSet)

    def test_get_material_flow_terms(self, model):
        for p in model.params.phase_list:
            for j in model.params.component_list:
                assert (
                    model.prop[1].get_material_flow_terms(p, j)
                    == model.prop[1].flow_mol * model.prop[1].phase_frac[p]
                )

    def test_get_enthalpy_flow_terms(self, model):
        for p in model.params.phase_list:
            assert model.prop[1].get_enthalpy_flow_terms(p) == (
                model.prop[1].enth_mol_phase[p]
                * model.prop[1].phase_frac[p]
                * model.prop[1].flow_mol
            )

    def test_get_material_density_terms(self, model):
        for p in model.params.phase_list:
            for j in model.params.component_list:
                assert (
                    model.prop[1].get_material_density_terms(p, j)
                    is model.prop[1].dens_mol_phase[p]
                )

    def test_get_energy_density_terms(self, model):
        for p in model.params.phase_list:
            assert model.prop[1].get_energy_density_terms(p) == (
                model.prop[1].dens_mol_phase[p]
                * model.prop[1].energy_internal_mol_phase[p]
            )

    def test_default_material_balance_type(self, model):
        assert (
            model.prop[1].default_material_balance_type()
            is MaterialBalanceType.componentTotal
        )

    def test_default_energy_balance_type(self, model):
        assert (
            model.prop[1].default_energy_balance_type()
            is EnergyBalanceType.enthalpyTotal
        )

    def test_define_state_vars(self, model):
        sv = model.prop[1].define_state_vars()
        assert len(sv) == 3
        for i in sv:
            assert i in ["flow_mol", "enth_mol", "pressure"]

    def test_define_port_members(self, model):
        sv = model.prop[1].define_port_members()
        assert len(sv) == 3
        for i in sv:
            assert i in ["flow_mol", "enth_mol", "pressure"]

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

    def test_extensive_state_vars(self, model):
        for i in model.prop[1].extensive_set:
            assert i is model.prop[1].flow_mol

    def test_intensive_state_vars(self, model):
        for i in model.prop[1].intensive_set:
            assert i in [model.prop[1].enth_mol, model.prop[1].pressure]

    def test_model_check(self, model):
        assert model.prop[1].model_check() is None

    def test_data_initialize(self, model):
        assert model.prop[1].initialize() is None

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

    def test_initialize_w_partial_state_1(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

        flags = model.prop.initialize(hold_state=True, state_args={"flow_mol": 30})

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

    def test_initialize_w_partial_state_2(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

        flags = model.prop.initialize(hold_state=True, state_args={"enth_mol": 300})

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

    def test_initialize_w_partial_state_3(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

        flags = model.prop.initialize(hold_state=True, state_args={"pressure": 3000})

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

    def test_config(self, model):
        assert model.params.config.phase_presentation == iapws95.PhaseType.L
        assert model.params.config.state_vars == iapws95.StateVars.PH

        assert len(model.params.phase_list) == 1
        for i in model.params.phase_list:
            assert i in ["Liq"]

    def test_build(self, model):
        model.prop = iapws95.Iapws95StateBlock(
            [1], default={"parameters": model.params}
        )

        assert isinstance(model.prop[1].flow_mol, Var)
        assert isinstance(model.prop[1].pressure, Var)
        assert isinstance(model.prop[1].enth_mol, Var)

        assert isinstance(model.prop[1].extensive_set, ComponentSet)
        assert isinstance(model.prop[1].intensive_set, ComponentSet)

    def test_get_material_flow_terms(self, model):
        for p in model.params.phase_list:
            for j in model.params.component_list:
                assert (
                    model.prop[1].get_material_flow_terms(p, j)
                    == model.prop[1].flow_mol * model.prop[1].phase_frac[p]
                )

    def test_get_enthalpy_flow_terms(self, model):
        for p in model.params.phase_list:
            assert model.prop[1].get_enthalpy_flow_terms(p) == (
                model.prop[1].enth_mol_phase[p]
                * model.prop[1].phase_frac[p]
                * model.prop[1].flow_mol
            )

    def test_get_material_density_terms(self, model):
        for p in model.params.phase_list:
            for j in model.params.component_list:
                assert (
                    model.prop[1].get_material_density_terms(p, j)
                    is model.prop[1].dens_mol_phase[p]
                )

    def test_get_energy_density_terms(self, model):
        for p in model.params.phase_list:
            assert model.prop[1].get_energy_density_terms(p) == (
                model.prop[1].dens_mol_phase[p]
                * model.prop[1].energy_internal_mol_phase[p]
            )

    def test_default_material_balance_type(self, model):
        assert (
            model.prop[1].default_material_balance_type()
            is MaterialBalanceType.componentTotal
        )

    def test_default_energy_balance_type(self, model):
        assert (
            model.prop[1].default_energy_balance_type()
            is EnergyBalanceType.enthalpyTotal
        )

    def test_define_state_vars(self, model):
        sv = model.prop[1].define_state_vars()
        assert len(sv) == 3
        for i in sv:
            assert i in ["flow_mol", "enth_mol", "pressure"]

    def test_define_port_members(self, model):
        sv = model.prop[1].define_port_members()
        assert len(sv) == 3
        for i in sv:
            assert i in ["flow_mol", "enth_mol", "pressure"]

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

    def test_extensive_state_vars(self, model):
        for i in model.prop[1].extensive_set:
            assert i is model.prop[1].flow_mol

    def test_intensive_state_vars(self, model):
        for i in model.prop[1].intensive_set:
            assert i in [model.prop[1].enth_mol, model.prop[1].pressure]

    def test_model_check(self, model):
        assert model.prop[1].model_check() is None

    def test_data_initialize(self, model):
        assert model.prop[1].initialize() is None

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

    def test_initialize_w_partial_state_1(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

        flags = model.prop.initialize(hold_state=True, state_args={"flow_mol": 30})

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

    def test_initialize_w_partial_state_2(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

        flags = model.prop.initialize(hold_state=True, state_args={"enth_mol": 300})

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

    def test_initialize_w_partial_state_3(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

        flags = model.prop.initialize(hold_state=True, state_args={"pressure": 3000})

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

    def test_config(self, model):
        assert model.params.config.phase_presentation == iapws95.PhaseType.G
        assert model.params.config.state_vars == iapws95.StateVars.PH

        assert len(model.params.phase_list) == 1
        for i in model.params.phase_list:
            assert i in ["Vap"]

    def test_build(self, model):
        model.prop = iapws95.Iapws95StateBlock(
            [1], default={"parameters": model.params}
        )

        assert isinstance(model.prop[1].flow_mol, Var)
        assert isinstance(model.prop[1].pressure, Var)
        assert isinstance(model.prop[1].enth_mol, Var)

        assert isinstance(model.prop[1].extensive_set, ComponentSet)
        assert isinstance(model.prop[1].intensive_set, ComponentSet)

    def test_get_material_flow_terms(self, model):
        for p in model.params.phase_list:
            for j in model.params.component_list:
                assert (
                    model.prop[1].get_material_flow_terms(p, j)
                    == model.prop[1].flow_mol * model.prop[1].phase_frac[p]
                )

    def test_get_enthalpy_flow_terms(self, model):
        for p in model.params.phase_list:
            assert model.prop[1].get_enthalpy_flow_terms(p) == (
                model.prop[1].enth_mol_phase[p]
                * model.prop[1].phase_frac[p]
                * model.prop[1].flow_mol
            )

    def test_get_material_density_terms(self, model):
        for p in model.params.phase_list:
            for j in model.params.component_list:
                assert (
                    model.prop[1].get_material_density_terms(p, j)
                    is model.prop[1].dens_mol_phase[p]
                )

    def test_get_energy_density_terms(self, model):
        for p in model.params.phase_list:
            assert model.prop[1].get_energy_density_terms(p) == (
                model.prop[1].dens_mol_phase[p]
                * model.prop[1].energy_internal_mol_phase[p]
            )

    def test_default_material_balance_type(self, model):
        assert (
            model.prop[1].default_material_balance_type()
            is MaterialBalanceType.componentTotal
        )

    def test_default_energy_balance_type(self, model):
        assert (
            model.prop[1].default_energy_balance_type()
            is EnergyBalanceType.enthalpyTotal
        )

    def test_define_state_vars(self, model):
        sv = model.prop[1].define_state_vars()
        assert len(sv) == 3
        for i in sv:
            assert i in ["flow_mol", "enth_mol", "pressure"]

    def test_define_port_members(self, model):
        sv = model.prop[1].define_port_members()
        assert len(sv) == 3
        for i in sv:
            assert i in ["flow_mol", "enth_mol", "pressure"]

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

    def test_extensive_state_vars(self, model):
        for i in model.prop[1].extensive_set:
            assert i is model.prop[1].flow_mol

    def test_intensive_state_vars(self, model):
        for i in model.prop[1].intensive_set:
            assert i in [model.prop[1].enth_mol, model.prop[1].pressure]

    def test_model_check(self, model):
        assert model.prop[1].model_check() is None

    def test_data_initialize(self, model):
        assert model.prop[1].initialize() is None

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

    def test_initialize_w_partial_state_1(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

        flags = model.prop.initialize(hold_state=True, state_args={"flow_mol": 30})

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

    def test_initialize_w_partial_state_2(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

        flags = model.prop.initialize(hold_state=True, state_args={"enth_mol": 300})

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

    def test_initialize_w_partial_state_3(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].enth_mol.fixed
        assert not model.prop[1].pressure.fixed

        flags = model.prop.initialize(hold_state=True, state_args={"pressure": 3000})

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

    def test_config(self, model):
        assert model.params.config.phase_presentation == iapws95.PhaseType.MIX
        assert model.params.config.state_vars == iapws95.StateVars.TPX

        assert len(model.params.phase_list) == 1
        for i in model.params.phase_list:
            assert i in ["Mix"]

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

    def test_get_material_flow_terms(self, model):
        for p in model.params.phase_list:
            for j in model.params.component_list:
                assert value(model.prop[1].get_material_flow_terms(p, j)) == value(
                    model.prop[1].flow_mol
                )

    def test_get_enthalpy_flow_terms(self, model):
        for p in model.params.phase_list:
            assert (
                model.prop[1].get_enthalpy_flow_terms(p)
                == model.prop[1].enth_mol * model.prop[1].flow_mol
            )

    def test_get_material_density_terms(self, model):
        for p in model.params.phase_list:
            for j in model.params.component_list:
                assert (
                    model.prop[1].get_material_density_terms(p, j)
                    is model.prop[1].dens_mol
                )

    def test_get_energy_density_terms(self, model):
        for p in model.params.phase_list:
            assert (
                model.prop[1].get_energy_density_terms(p)
                == model.prop[1].dens_mol * model.prop[1].energy_internal_mol
            )

    def test_default_material_balance_type(self, model):
        assert (
            model.prop[1].default_material_balance_type()
            is MaterialBalanceType.componentTotal
        )

    def test_default_energy_balance_type(self, model):
        assert (
            model.prop[1].default_energy_balance_type()
            is EnergyBalanceType.enthalpyTotal
        )

    def test_define_state_vars(self, model):
        sv = model.prop[1].define_state_vars()
        assert len(sv) == 4
        for i in sv:
            assert i in ["flow_mol", "temperature", "pressure", "vapor_frac"]

    def test_define_port_members(self, model):
        sv = model.prop[1].define_port_members()
        assert len(sv) == 4
        for i in sv:
            assert i in ["flow_mol", "temperature", "pressure", "vapor_frac"]

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

    def test_extensive_state_vars(self, model):
        for i in model.prop[1].extensive_set:
            assert i is model.prop[1].flow_mol

    def test_intensive_state_vars(self, model):
        for i in model.prop[1].intensive_set:
            assert i in [
                model.prop[1].temperature,
                model.prop[1].pressure,
                model.prop[1].vapor_frac,
            ]

    def test_model_check(self, model):
        assert model.prop[1].model_check() is None

    def test_data_initialize(self, model):
        assert model.prop[1].initialize() is None

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

    def test_initialize_w_partial_state_1(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

        flags = model.prop.initialize(hold_state=True, state_args={"flow_mol": 30})

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

    def test_initialize_w_partial_state_2(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

        flags = model.prop.initialize(hold_state=True, state_args={"temperature": 300})

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

    def test_initialize_w_partial_state_3(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

        flags = model.prop.initialize(hold_state=True, state_args={"pressure": 3000})

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

    def test_config(self, model):
        assert model.params.config.phase_presentation == iapws95.PhaseType.LG
        assert model.params.config.state_vars == iapws95.StateVars.TPX

        assert len(model.params.phase_list) == 2
        for i in model.params.phase_list:
            assert i in ["Liq", "Vap"]

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

    def test_get_material_flow_terms(self, model):
        for p in model.params.phase_list:
            for j in model.params.component_list:
                assert (
                    model.prop[1].get_material_flow_terms(p, j)
                    == model.prop[1].flow_mol * model.prop[1].phase_frac[p]
                )

    def test_get_enthalpy_flow_terms(self, model):
        for p in model.params.phase_list:
            assert model.prop[1].get_enthalpy_flow_terms(p) == (
                model.prop[1].enth_mol_phase[p]
                * model.prop[1].phase_frac[p]
                * model.prop[1].flow_mol
            )

    def test_get_material_density_terms(self, model):
        for p in model.params.phase_list:
            for j in model.params.component_list:
                assert (
                    model.prop[1].get_material_density_terms(p, j)
                    is model.prop[1].dens_mol_phase[p]
                )

    def test_get_energy_density_terms(self, model):
        for p in model.params.phase_list:
            assert model.prop[1].get_energy_density_terms(p) == (
                model.prop[1].dens_mol_phase[p]
                * model.prop[1].energy_internal_mol_phase[p]
            )

    def test_default_material_balance_type(self, model):
        assert (
            model.prop[1].default_material_balance_type()
            is MaterialBalanceType.componentTotal
        )

    def test_default_energy_balance_type(self, model):
        assert (
            model.prop[1].default_energy_balance_type()
            is EnergyBalanceType.enthalpyTotal
        )

    def test_define_state_vars(self, model):
        sv = model.prop[1].define_state_vars()
        assert len(sv) == 4
        for i in sv:
            assert i in ["flow_mol", "temperature", "pressure", "vapor_frac"]

    def test_define_port_members(self, model):
        sv = model.prop[1].define_port_members()
        assert len(sv) == 4
        for i in sv:
            assert i in ["flow_mol", "temperature", "pressure", "vapor_frac"]

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

    def test_extensive_state_vars(self, model):
        for i in model.prop[1].extensive_set:
            assert i is model.prop[1].flow_mol

    def test_intensive_state_vars(self, model):
        for i in model.prop[1].intensive_set:
            assert i in [
                model.prop[1].temperature,
                model.prop[1].pressure,
                model.prop[1].vapor_frac,
            ]

    def test_model_check(self, model):
        assert model.prop[1].model_check() is None

    def test_data_initialize(self, model):
        assert model.prop[1].initialize() is None

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

    def test_initialize_w_partial_state_1(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert not model.prop[1].vapor_frac.fixed

        flags = model.prop.initialize(hold_state=True, state_args={"flow_mol": 30})

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

    def test_initialize_w_partial_state_2(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert not model.prop[1].vapor_frac.fixed

        flags = model.prop.initialize(hold_state=True, state_args={"temperature": 300})

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

    def test_initialize_w_partial_state_3(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert not model.prop[1].vapor_frac.fixed

        flags = model.prop.initialize(hold_state=True, state_args={"pressure": 3000})

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

    def test_config(self, model):
        assert model.params.config.phase_presentation == iapws95.PhaseType.L
        assert model.params.config.state_vars == iapws95.StateVars.TPX

        assert len(model.params.phase_list) == 1
        for i in model.params.phase_list:
            assert i in ["Liq"]

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

    def test_get_material_flow_terms(self, model):
        for p in model.params.phase_list:
            for j in model.params.component_list:
                assert (
                    model.prop[1].get_material_flow_terms(p, j)
                    == model.prop[1].flow_mol * model.prop[1].phase_frac[p]
                )

    def test_get_enthalpy_flow_terms(self, model):
        for p in model.params.phase_list:
            assert model.prop[1].get_enthalpy_flow_terms(p) == (
                model.prop[1].enth_mol_phase[p]
                * model.prop[1].phase_frac[p]
                * model.prop[1].flow_mol
            )

    def test_get_material_density_terms(self, model):
        for p in model.params.phase_list:
            for j in model.params.component_list:
                assert (
                    model.prop[1].get_material_density_terms(p, j)
                    is model.prop[1].dens_mol_phase[p]
                )

    def test_get_energy_density_terms(self, model):
        for p in model.params.phase_list:
            assert model.prop[1].get_energy_density_terms(p) == (
                model.prop[1].dens_mol_phase[p]
                * model.prop[1].energy_internal_mol_phase[p]
            )

    def test_default_material_balance_type(self, model):
        assert (
            model.prop[1].default_material_balance_type()
            is MaterialBalanceType.componentTotal
        )

    def test_default_energy_balance_type(self, model):
        assert (
            model.prop[1].default_energy_balance_type()
            is EnergyBalanceType.enthalpyTotal
        )

    def test_define_state_vars(self, model):
        sv = model.prop[1].define_state_vars()
        assert len(sv) == 3
        for i in sv:
            assert i in ["flow_mol", "temperature", "pressure"]

    def test_define_port_members(self, model):
        sv = model.prop[1].define_port_members()
        assert len(sv) == 3
        for i in sv:
            assert i in ["flow_mol", "temperature", "pressure"]

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

    def test_extensive_state_vars(self, model):
        for i in model.prop[1].extensive_set:
            assert i is model.prop[1].flow_mol

    def test_intensive_state_vars(self, model):
        for i in model.prop[1].intensive_set:
            assert i in [
                model.prop[1].temperature,
                model.prop[1].pressure,
                model.prop[1].vapor_frac,
            ]

    def test_model_check(self, model):
        assert model.prop[1].model_check() is None

    def test_data_initialize(self, model):
        assert model.prop[1].initialize() is None

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

    def test_initialize_w_partial_state_1(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

        flags = model.prop.initialize(hold_state=True, state_args={"flow_mol": 30})

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

    def test_initialize_w_partial_state_2(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

        flags = model.prop.initialize(hold_state=True, state_args={"temperature": 300})

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

    def test_initialize_w_partial_state_3(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

        flags = model.prop.initialize(hold_state=True, state_args={"pressure": 3000})

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

    def test_config(self, model):
        assert model.params.config.phase_presentation == iapws95.PhaseType.G
        assert model.params.config.state_vars == iapws95.StateVars.TPX

        assert len(model.params.phase_list) == 1
        for i in model.params.phase_list:
            assert i in ["Vap"]

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

    def test_get_material_flow_terms(self, model):
        for p in model.params.phase_list:
            for j in model.params.component_list:
                assert (
                    model.prop[1].get_material_flow_terms(p, j)
                    == model.prop[1].flow_mol * model.prop[1].phase_frac[p]
                )

    def test_get_enthalpy_flow_terms(self, model):
        for p in model.params.phase_list:
            assert model.prop[1].get_enthalpy_flow_terms(p) == (
                model.prop[1].enth_mol_phase[p]
                * model.prop[1].phase_frac[p]
                * model.prop[1].flow_mol
            )

    def test_get_material_density_terms(self, model):
        for p in model.params.phase_list:
            for j in model.params.component_list:
                assert (
                    model.prop[1].get_material_density_terms(p, j)
                    is model.prop[1].dens_mol_phase[p]
                )

    def test_get_energy_density_terms(self, model):
        for p in model.params.phase_list:
            assert model.prop[1].get_energy_density_terms(p) == (
                model.prop[1].dens_mol_phase[p]
                * model.prop[1].energy_internal_mol_phase[p]
            )

    def test_default_material_balance_type(self, model):
        assert (
            model.prop[1].default_material_balance_type()
            is MaterialBalanceType.componentTotal
        )

    def test_default_energy_balance_type(self, model):
        assert (
            model.prop[1].default_energy_balance_type()
            is EnergyBalanceType.enthalpyTotal
        )

    def test_define_state_vars(self, model):
        sv = model.prop[1].define_state_vars()
        assert len(sv) == 3
        for i in sv:
            assert i in ["flow_mol", "temperature", "pressure"]

    def test_define_port_members(self, model):
        sv = model.prop[1].define_port_members()
        assert len(sv) == 3
        for i in sv:
            assert i in ["flow_mol", "temperature", "pressure"]

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

    def test_extensive_state_vars(self, model):
        for i in model.prop[1].extensive_set:
            assert i is model.prop[1].flow_mol

    def test_intensive_state_vars(self, model):
        for i in model.prop[1].intensive_set:
            assert i in [
                model.prop[1].temperature,
                model.prop[1].pressure,
                model.prop[1].vapor_frac,
            ]

    def test_model_check(self, model):
        assert model.prop[1].model_check() is None

    def test_data_initialize(self, model):
        assert model.prop[1].initialize() is None

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

    def test_initialize_w_partial_state_1(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

        flags = model.prop.initialize(hold_state=True, state_args={"flow_mol": 30})

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

    def test_initialize_w_partial_state_2(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

        flags = model.prop.initialize(hold_state=True, state_args={"temperature": 300})

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

    def test_initialize_w_partial_state_3(self, model):
        assert not model.prop[1].flow_mol.fixed
        assert not model.prop[1].temperature.fixed
        assert not model.prop[1].pressure.fixed
        assert model.prop[1].vapor_frac.fixed

        flags = model.prop.initialize(hold_state=True, state_args={"pressure": 3000})

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


# -----------------------------------------------------------------------------
# Test values against data tables
def read_data(fname, col):
    dfile = os.path.join(this_file_dir(), fname)
    cond = []  # Tuple (T [K],P [Pa], data) pressure in file is MPa
    with open(dfile, "r") as csvfile:
        dat = csv.reader(csvfile, delimiter="\t", quotechar='"')
        for i in range(7):
            next(dat) # skip header
        for row in dat:
            try:
                x = float(row[col])
            except:
                x = row[col]
            cond.append((float(row[0]), float(row[1]) * 1e6, x))
    return cond


@pytest.mark.slow
@pytest.mark.skipif(not prop_available, reason="IAPWS not available")
@pytest.mark.nocircleci()
class TestPHvalues(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.prop_param = iapws95.Iapws95ParameterBlock()
        model.prop_in = iapws95.Iapws95StateBlock(
            default={"parameters": model.prop_param}
        )

        return model
    
    def test_tau_sat(self, model):
        cond = read_data("sat_prop_iapws95_nist_webbook.txt", col=2)
        for c in cond:
            tau = value(model.prop_in.func_tau_sat(c[1] / 1000.0))
            T = 647.096 / tau
            print("{}, {}, {}".format(c[1], c[0], T))
            assert abs(T - c[0]) < 0.1

    def test_liquid_density_sat(self, model):
        cond = read_data("sat_prop_iapws95_nist_webbook.txt", col=2)
        for c in cond:
            if c[0] > 645:  # getting very close to critical point
                tol = 0.05
            else:
                tol = 0.001
            model.prop_in.temperature.set_value(c[0])
            model.prop_in.pressure = c[1]
            rho = value(model.prop_in.dens_mass_phase["Liq"])
            assert rho == pytest.approx(c[2], rel=tol)

    def test_vapor_density_sat(self, model):
        cond = read_data("sat_prop_iapws95_nist_webbook.txt", col=14)
        for c in cond:
            if c[0] > 645:  # getting very close to critical point
                tol = 0.01
            else:
                tol = 0.001
            model.prop_in.temperature.set_value(c[0])
            model.prop_in.pressure = c[1]
            rho = value(model.prop_in.dens_mass_phase["Vap"])
            assert rho == pytest.approx(c[2], rel=tol)

    def test_liquid_enthalpy_sat(self, model):
        cond = read_data("sat_prop_iapws95_nist_webbook.txt", col=5)
        for c in cond:
            if c[0] > 645:  # getting very close to critical point
                tol = 0.01
            else:
                tol = 0.001
            model.prop_in.pressure = c[1]
            enth = value(
                model.prop_in.enth_mol_sat_phase["Liq"] / model.prop_in.mw / 1000.0
            )
            assert enth == pytest.approx(c[2], rel=tol)

    def test_vapor_enthalpy_sat(self, model):
        cond = read_data("sat_prop_iapws95_nist_webbook.txt", col=17)
        for c in cond:
            if c[0] > 645:  # getting very close to critical point
                tol = 0.01
            else:
                tol = 0.001
            model.prop_in.pressure = c[1]
            enth = value(
                model.prop_in.enth_mol_sat_phase["Vap"] / model.prop_in.mw / 1000.0
            )
            assert enth == pytest.approx(c[2], rel=tol)

    def test_enthalpy_of_vaporization(self, model):
        cond_liq = read_data("sat_prop_iapws95_nist_webbook.txt", col=5)
        cond_vap = read_data("sat_prop_iapws95_nist_webbook.txt", col=17)
        for i, c in enumerate(cond_liq):
            if c[0] > 645:  # getting very close to critical point
                tol = 0.05
            else:
                tol = 0.001
            model.prop_in.pressure.value = c[1]
            enth = value(model.prop_in.dh_vap_mol / model.prop_in.mw / 1000.0)
            enth_dat = cond_vap[i][2] - c[2]
            if enth < 1e-1:
                assert enth == pytest.approx(enth_dat, abs=1e-1)
            else:
                assert enth == pytest.approx(enth_dat, rel=tol)
        # Over Critical Pressure
        model.prop_in.pressure = model.prop_in.config.parameters.pressure_crit * 1.1
        enth = value(model.prop_in.dh_vap_mol / model.prop_in.mw / 1000.0)
        assert abs(enth) < 0.001

    def test_density(self, model):
        cond = read_data("prop_iapws95_nist_webbook.txt", col=2)
        phase = read_data("prop_iapws95_nist_webbook.txt", col=13)
        for i, c in enumerate(cond):
            if phase[i][2] in ["liquid", "supercritical"]:
                p = "Liq"
            else:
                p = "Vap"
            model.prop_in.temperature.set_value(c[0])
            model.prop_in.pressure = c[1]
            rho = value(model.prop_in.dens_mass_phase[p])
            if rho > 250 and rho < 420 and c[0] < 700 and c[0] > 645:
                tol = 0.03  # steep part in sc region
            elif c[1] < 20:
                tol = 0.01  # very low pressure < 20 Pa
            else:
                tol = 0.001
            print(c[0], c[1], p)
            assert rho == pytest.approx(c[2], rel=tol)

    def test_enthalpy(self, model):
        cond = read_data("prop_iapws95_nist_webbook.txt", col=5)
        phase = read_data("prop_iapws95_nist_webbook.txt", col=13)
        for i, c in enumerate(cond):
            if phase[i][2] in ["liquid", "supercritical"]:
                p = "Liq"
            else:
                p = "Vap"
            model.prop_in.temperature.set_value(c[0])
            model.prop_in.pressure = c[1]
            h = value(model.prop_in.enth_mol_phase[p] / model.prop_in.mw / 1000)
            rho = value(model.prop_in.dens_mass_phase[p])
            if rho > 250 and rho < 420 and c[0] < 700 and c[0] > 640:
                tol = 0.03  # steep part in sc region
            elif c[1] < 20:
                tol = 0.005  # very low pressure < 20 Pa
            else:
                tol = 0.0015
            assert h == pytest.approx(c[2], rel=tol)

    def test_enthalpy_vapor_as_function_of_p_and_tau(self, model):
        cond = read_data("prop_iapws95_nist_webbook.txt", col=5)
        phase = read_data("prop_iapws95_nist_webbook.txt", col=13)
        for i, c in enumerate(cond):
            if phase[i][2] in ["liquid", "supercritical"]:
                continue
            model.prop_in.temperature.set_value(c[0])
            h = value(model.prop_in.func_hvpt(c[1] / 1000, 647.096 / c[0]))
            rho = value(model.prop_in.dens_mass_phase["Vap"])
            if rho > 250 and rho < 420 and c[0] < 700 and c[0] > 640:
                tol = 0.03  # steep part in sc region
            else:
                tol = 0.003
            assert abs(h - c[2]) / c[2] < tol

    def test_enthalpy_liquid_as_function_of_p_and_tau(self, model):
        cond = read_data("prop_iapws95_nist_webbook.txt", col=5)
        phase = read_data("prop_iapws95_nist_webbook.txt", col=13)
        for i, c in enumerate(cond):
            p = phase[i][2]
            if p in ["vapor"]:
                continue
            model.prop_in.temperature.set_value(c[0])
            h = value(model.prop_in.func_hlpt(c[1] / 1000, 647.096 / c[0]))
            rho = value(model.prop_in.dens_mass_phase["Liq"])
            if rho > 250 and rho < 420 and c[0] < 700 and c[0] > 640:
                tol = 0.03  # steep part in sc region
            else:
                tol = 0.006
            assert abs(h - c[2]) / c[2] < tol

    def test_entropy_vapor_as_function_of_p_and_tau(self, model):
        cond = read_data("prop_iapws95_nist_webbook.txt", col=6)
        phase = read_data("prop_iapws95_nist_webbook.txt", col=13)
        for i, c in enumerate(cond):
            if phase[i][2] in ["liquid", "supercritical"]:
                continue
            model.prop_in.temperature.set_value(c[0])
            s = value(model.prop_in.func_svpt(c[1] / 1000, 647.096 / c[0]))
            rho = value(model.prop_in.dens_mass_phase["Vap"])
            if rho > 250 and rho < 420 and c[0] < 700 and c[0] > 640:
                tol = 0.03  # steep part in sc region
            else:
                tol = 0.003
            assert abs(s - c[2]) / c[2] < tol

    def test_entropy_liquid_as_function_of_p_and_tau(self, model):
        cond = read_data("prop_iapws95_nist_webbook.txt", col=6)
        phase = read_data("prop_iapws95_nist_webbook.txt", col=13)
        for i, c in enumerate(cond):
            p = phase[i][2]
            if p in ["vapor"]:
                continue
            model.prop_in.temperature.set_value(c[0])
            s = value(model.prop_in.func_slpt(c[1] / 1000, 647.096 / c[0]))
            rho = value(model.prop_in.dens_mass_phase["Liq"])
            if rho > 250 and rho < 420 and c[0] < 700 and c[0] > 640:
                tol = 0.03  # steep part in sc region
            else:
                tol = 0.006
            assert abs(s - c[2]) / c[2] < tol

    def test_entropy(self, model):
        cond = read_data("prop_iapws95_nist_webbook.txt", col=6)
        phase = read_data("prop_iapws95_nist_webbook.txt", col=13)
        for i, c in enumerate(cond):
            if phase[i][2] in ["liquid", "supercritical"]:
                p = "Liq"
            else:
                p = "Vap"
            model.prop_in.temperature.set_value(c[0])
            model.prop_in.pressure = c[1]
            s = value(model.prop_in.entr_mol_phase[p] / model.prop_in.mw / 1000)
            rho = value(model.prop_in.dens_mass_phase[p])
            if rho > 250 and rho < 420 and c[0] < 700 and c[0] > 640:
                tol = 0.03  # steep part in sc region
            else:
                tol = 0.003
            assert abs(s - c[2]) / c[2] < tol

    def test_internal_energy(self, model):
        cond = read_data("prop_iapws95_nist_webbook.txt", col=4)
        phase = read_data("prop_iapws95_nist_webbook.txt", col=13)
        for i, c in enumerate(cond):
            if phase[i][2] in ["liquid", "supercritical"]:
                p = "Liq"
            else:
                p = "Vap"
            model.prop_in.temperature.set_value(c[0])
            model.prop_in.pressure = c[1]
            u = value(
                model.prop_in.energy_internal_mol_phase[p] / model.prop_in.mw / 1000
            )
            rho = value(model.prop_in.dens_mass_phase[p])
            if rho > 250 and rho < 420 and c[0] < 700 and c[0] > 640:
                tol = 0.02  # steep part in sc region
            else:
                tol = 0.002
            assert abs(u - c[2]) / c[2] < tol

    def test_speed_of_sound(self, model):
        cond = read_data("prop_iapws95_nist_webbook.txt", col=9)
        phase = read_data("prop_iapws95_nist_webbook.txt", col=13)
        for i, c in enumerate(cond):
            if c[2] == "undefined":
                continue
            if (c[0] > 640 and c[0] < 659) and (c[1] > 2.1e7 and c[1] < 2.5e7):
                # near critical and non-analytic terms were omitted
                continue
            if phase[i][2] in ["liquid", "supercritical"]:
                p = "Liq"
            else:
                p = "Vap"
            model.prop_in.temperature.set_value(c[0])
            model.prop_in.pressure = c[1]
            w = value(model.prop_in.speed_sound_phase[p])
            rho = value(model.prop_in.dens_mass_phase[p])
            if rho > 250 and rho < 420 and c[0] < 700 and c[0] > 640:
                tol = 0.03  # steep part in sc region
            else:
                tol = 0.005
            assert abs(w - c[2]) / c[2] < tol

    def test_cp(self, model):
        cond = read_data("prop_iapws95_nist_webbook.txt", col=8)
        phase = read_data("prop_iapws95_nist_webbook.txt", col=13)
        for i, c in enumerate(cond):
            if c[2] == "undefined":
                continue
            if (c[0] > 640 and c[0] < 680) and (c[1] > 2.1e7 and c[1] < 3.1e7):
                # near critical and non-analytic terms were omitted
                continue
            if phase[i][2] in ["liquid", "supercritical"]:
                p = "Liq"
            else:
                p = "Vap"
            model.prop_in.temperature.set_value(c[0])
            model.prop_in.pressure = c[1]
            cp = value(model.prop_in.cp_mol_phase[p] / model.prop_in.mw / 1000)
            rho = value(model.prop_in.dens_mass_phase[p])
            if rho > 250 and rho < 420 and c[0] < 700 and c[0] > 640:
                tol = 0.03  # steep part in sc region
            else:
                tol = 0.005
            assert abs(cp - c[2]) / c[2] < tol

    def test_cv(self, model):
        cond = read_data("prop_iapws95_nist_webbook.txt", col=7)
        phase = read_data("prop_iapws95_nist_webbook.txt", col=13)
        for i, c in enumerate(cond):
            if c[2] == "undefined":
                continue
            if (c[0] > 640 and c[0] < 680) and (c[1] > 2.1e7 and c[1] < 3.1e7):
                # near critical and non-analytic terms were omitted
                continue
            if phase[i][2] in ["liquid", "supercritical"]:
                p = "Liq"
            else:
                p = "Vap"
            model.prop_in.temperature.set_value(c[0])
            model.prop_in.pressure = c[1]
            cv = value(model.prop_in.cv_mol_phase[p] / model.prop_in.mw / 1000)
            rho = value(model.prop_in.dens_mass_phase[p])
            if rho > 250 and rho < 420 and c[0] < 700 and c[0] > 640:
                tol = 0.03  # steep part in sc region
            else:
                tol = 0.003
            assert abs(cv - c[2]) / c[2] < tol
