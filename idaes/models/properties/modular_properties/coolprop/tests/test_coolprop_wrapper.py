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
from numpy import arange

from pyomo.environ import (
    Block,
    ConcreteModel,
    Param,
    SolverStatus,
    TerminationCondition,
    units as pyunits,
    value,
    Var,
)
from pyomo.util.check_units import assert_units_equivalent

from idaes.core import FlowsheetBlock, Component, LiquidPhase, VaporPhase
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.eos.ceos import Cubic, CubicType
from idaes.core.solvers import get_solver
from idaes.models.properties.modular_properties.pure.ConstantProperties import Constant


CoolProp = pytest.importorskip("CoolProp.CoolProp", reason="CoolProp not installed")

from idaes.models.properties.modular_properties.coolprop.coolprop_wrapper import (
    CoolPropWrapper,
    CoolPropExpressionError,
    CoolPropPropertyError,
)

solver = get_solver()


class TestWrapper:
    @pytest.mark.unit
    def test_load_component(self):
        # Clear cached components to be sure
        CoolPropWrapper._cached_components = {}

        # Load parameters for oxygen
        prop_dict = CoolPropWrapper._load_component_data("Oxygen")

        assert prop_dict is not None
        assert "Oxygen" in CoolPropWrapper._cached_components
        assert CoolPropWrapper._cached_components["Oxygen"] is prop_dict
        assert isinstance(prop_dict, dict)
        for k in ["STATES", "ANCILLARIES"]:
            assert k in prop_dict

    @pytest.mark.unit
    def test_load_component_invalid(self):
        with pytest.raises(
            RuntimeError,
            match="Failed to find component foo in CoolProp " "JSON database.",
        ):
            CoolPropWrapper._load_component_data("foo")

    @pytest.mark.unit
    def test_get_component(self):
        prop_dict = CoolPropWrapper._get_component_data("Oxygen")

        assert CoolPropWrapper._cached_components["Oxygen"] is prop_dict

    @pytest.mark.unit
    def test_get_component_alias(self):
        # Load oxygen properties using aliases
        prop_dict = CoolPropWrapper._get_component_data("R732")

        assert CoolPropWrapper._cached_components["Oxygen"] is prop_dict
        assert "R732" in CoolPropWrapper._cached_components
        assert (
            CoolPropWrapper._cached_components["R732"]
            is CoolPropWrapper._cached_components["Oxygen"]
        )

    @pytest.mark.unit
    def test_load_component_by_alias(self):
        # Load CO2 data using one of its aliases
        prop_dict = CoolPropWrapper._load_component_data("R744")

        assert CoolPropWrapper._cached_components["R744"] is prop_dict
        assert "R744" in CoolPropWrapper._cached_components

        # Retrieve CO2 data using its normal name
        prop_dict2 = CoolPropWrapper._get_component_data("CO2")

        assert prop_dict2 is prop_dict
        assert "CO2" in CoolPropWrapper._cached_components
        assert (
            CoolPropWrapper._cached_components["CO2"]
            is CoolPropWrapper._cached_components["R744"]
        )

    @pytest.mark.unit
    def test_get_component_invalid(self):
        with pytest.raises(
            RuntimeError,
            match="Failed to find component foo in CoolProp " "JSON database.",
        ):
            CoolPropWrapper._get_component_data("foo")

    @pytest.mark.unit
    def test_flush_cached_components(self):
        CoolPropWrapper.flush_cached_components()

        assert CoolPropWrapper._cached_components == {}

    @pytest.mark.unit
    def test_get_critical_properties(self):
        Tc = CoolPropWrapper._get_critical_property("Acetone", "T")
        assert Tc == (508.1, pyunits.K)

        Pc = CoolPropWrapper._get_critical_property("Acetone", "p")
        assert Pc == (4700000.0, pyunits.Pa)

        rhoc = CoolPropWrapper._get_critical_property("Acetone", "rhomolar")
        assert rhoc[0] == 4699.999999999999
        assert_units_equivalent(rhoc[1], pyunits.mol / pyunits.m**3)

        hc = CoolPropWrapper._get_critical_property("Acetone", "hmolar")
        assert hc[0] == 31614.73051047263
        assert_units_equivalent(hc[1], pyunits.J / pyunits.mol)

        sc = CoolPropWrapper._get_critical_property("Acetone", "smolar")
        assert sc[0] == 72.97112978635582
        assert_units_equivalent(sc[1], pyunits.J / pyunits.mol / pyunits.K)

    @pytest.mark.unit
    def test_get_eos_properties(self):
        omega = CoolPropWrapper._get_eos_property("Acetone", "acentric")
        assert omega == (0.3071, pyunits.dimensionless)

        mw = CoolPropWrapper._get_eos_property("Acetone", "molar_mass")
        assert mw[0] == 0.05807914
        assert_units_equivalent(mw[1], pyunits.kg / pyunits.mol)

    @pytest.mark.unit
    def test_get_parameter_value(self):
        Tc = CoolPropWrapper.get_parameter_value("Acetone", "temperature_crit")
        assert Tc == (508.1, pyunits.K)

        Pc = CoolPropWrapper.get_parameter_value("Acetone", "pressure_crit")
        assert Pc == (4700000.0, pyunits.Pa)

        rhoc = CoolPropWrapper.get_parameter_value("Acetone", "dens_mol_crit")
        assert rhoc[0] == 4699.999999999999
        assert_units_equivalent(rhoc[1], pyunits.mol / pyunits.m**3)

        hc = CoolPropWrapper.get_parameter_value("Acetone", "enth_mol_crit")
        assert hc[0] == 31614.73051047263
        assert_units_equivalent(hc[1], pyunits.J / pyunits.mol)

        sc = CoolPropWrapper.get_parameter_value("Acetone", "entr_mol_crit")
        assert sc[0] == 72.97112978635582
        assert_units_equivalent(sc[1], pyunits.J / pyunits.mol / pyunits.K)

        omega = CoolPropWrapper.get_parameter_value("Acetone", "omega")
        assert omega == (0.3071, pyunits.dimensionless)

        mw = CoolPropWrapper.get_parameter_value("Acetone", "mw")
        assert mw[0] == 0.05807914
        assert_units_equivalent(mw[1], pyunits.kg / pyunits.mol)

    @pytest.mark.unit
    def test_pressure_sat_no_params(self):
        m = ConcreteModel()
        # Dummy a Component object with the name of a property which does not
        # have data for saturation pressure
        m.Air = Block()

        with pytest.raises(
            CoolPropPropertyError,
            match="Could not retrieve parameters for "
            "pressure_sat of component Air from CoolProp. This "
            "likely indicates that CoolProp does not have "
            "values for the necessary parameters.",
        ):
            CoolPropWrapper.pressure_sat_comp.build_parameters(m.Air)

    @pytest.mark.unit
    def test_pressure_sat_unrecognised_form(self):
        m = ConcreteModel()
        # First, add a dummy component to the cached components which uses
        # unsupoorted form
        CoolPropWrapper._cached_components["TestComp"] = {
            "ANCILLARIES": {"pS": {"type": "foo", "using_tau_r": True}},
            "INFO": {"ALIASES": ["testcomp"], "NAME": "TestComp"},
        }
        m.TestComp = Block()

        with pytest.raises(
            CoolPropExpressionError,
            match="Found unsupported expression form for "
            "pressure_sat of component TestComp. This likely "
            "occured due to changes in CoolProp and the "
            "interface should be updated.",
        ):
            CoolPropWrapper.pressure_sat_comp.build_parameters(m.TestComp)

        # Pressure_sat uses has two parts to form. Set type to supported form
        # and using_tau_r to False (unsupported)
        CoolPropWrapper._cached_components["TestComp"]["ANCILLARIES"]["pS"][
            "type"
        ] = "pL"
        CoolPropWrapper._cached_components["TestComp"]["ANCILLARIES"]["pS"][
            "using_tau_r"
        ] = False
        with pytest.raises(
            CoolPropExpressionError,
            match="Found unsupported expression form for "
            "pressure_sat of component TestComp. This likely "
            "occured due to changes in CoolProp and the "
            "interface should be updated.",
        ):
            CoolPropWrapper.pressure_sat_comp.build_parameters(m.TestComp)


class TestCoolPropProperties(object):
    @pytest.fixture(scope="class")
    def m(self):
        # Clear cached components to ensure clean slate
        CoolPropWrapper.flush_cached_components()

        m = ConcreteModel()

        m.fs = FlowsheetBlock(dynamic=False)

        configuration = {
            # Specifying components
            "components": {
                "benzene": {
                    "type": Component,
                    "dens_mol_liq_comp": CoolPropWrapper,
                    "enth_mol_liq_comp": CoolPropWrapper,
                    "enth_mol_ig_comp": CoolPropWrapper,
                    "entr_mol_liq_comp": CoolPropWrapper,
                    "entr_mol_ig_comp": CoolPropWrapper,
                    "pressure_sat_comp": CoolPropWrapper,
                    "parameter_data": {
                        "mw": CoolPropWrapper,
                        "dens_mol_crit": CoolPropWrapper,
                        "pressure_crit": CoolPropWrapper,
                        "temperature_crit": CoolPropWrapper,
                        "omega": CoolPropWrapper,
                    },
                }
            },
            # Specifying phases
            "phases": {
                "Liq": {
                    "type": LiquidPhase,
                    "equation_of_state": Cubic,
                    "equation_of_state_options": {"type": CubicType.PR},
                }
            },
            # Set base units of measurement
            "base_units": {
                "time": pyunits.s,
                "length": pyunits.m,
                "mass": pyunits.kg,
                "amount": pyunits.mol,
                "temperature": pyunits.K,
            },
            # Specifying state definition
            "state_definition": FTPx,
            "state_bounds": {
                "flow_mol": (0, 100, 1000, pyunits.mol / pyunits.s),
                "temperature": (273.15, 300, 500, pyunits.K),
                "pressure": (5e4, 1e5, 1e6, pyunits.Pa),
            },
            "pressure_ref": (101325, pyunits.Pa),
            "temperature_ref": (298.15, pyunits.K),
            "parameter_data": {"PR_kappa": {("benzene", "benzene"): 0.000}},
        }

        m.fs.props = GenericParameterBlock(**configuration)

        m.fs.state = m.fs.props.build_state_block([0], defined_state=True)

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
        assert value(m.fs.props.benzene.temperature_crit) == CoolProp.PropsSI(
            "TCRIT", "Benzene"
        )

        assert isinstance(m.fs.props.benzene.pressure_crit, Var)
        assert m.fs.props.benzene.pressure_crit.fixed
        assert value(m.fs.props.benzene.pressure_crit) == CoolProp.PropsSI(
            "PCRIT", "Benzene"
        )

        assert isinstance(m.fs.props.benzene.dens_mol_crit, Var)
        assert m.fs.props.benzene.dens_mol_crit.fixed
        assert value(m.fs.props.benzene.dens_mol_crit) == CoolProp.PropsSI(
            "RHOMOLAR_CRITICAL", "Benzene"
        )

        assert isinstance(m.fs.props.benzene.mw, Param)
        assert value(m.fs.props.benzene.mw) == CoolProp.PropsSI("molarmass", "Benzene")

        assert isinstance(m.fs.props.benzene.omega, Var)
        assert m.fs.props.benzene.omega.fixed
        assert value(m.fs.props.benzene.omega) == CoolProp.PropsSI(
            "acentric", "Benzene"
        )

    @pytest.mark.unit
    def test_dens_mol_liq_comp(self, m):
        assert isinstance(m.fs.props.benzene.dens_mol_liq_comp_coeff_n1, Var)
        assert isinstance(m.fs.props.benzene.dens_mol_liq_comp_coeff_n2, Var)
        assert isinstance(m.fs.props.benzene.dens_mol_liq_comp_coeff_n3, Var)
        assert isinstance(m.fs.props.benzene.dens_mol_liq_comp_coeff_n4, Var)
        assert isinstance(m.fs.props.benzene.dens_mol_liq_comp_coeff_n5, Var)
        assert isinstance(m.fs.props.benzene.dens_mol_liq_comp_coeff_n6, Var)

        assert isinstance(m.fs.props.benzene.dens_mol_liq_comp_coeff_t1, Var)
        assert isinstance(m.fs.props.benzene.dens_mol_liq_comp_coeff_t2, Var)
        assert isinstance(m.fs.props.benzene.dens_mol_liq_comp_coeff_t3, Var)
        assert isinstance(m.fs.props.benzene.dens_mol_liq_comp_coeff_t4, Var)
        assert isinstance(m.fs.props.benzene.dens_mol_liq_comp_coeff_t5, Var)
        assert isinstance(m.fs.props.benzene.dens_mol_liq_comp_coeff_t6, Var)

        assert m.fs.props.benzene.dens_mol_liq_comp_coeff_n1.value == pytest.approx(
            2.852587673922022, rel=1e-10
        )
        assert m.fs.props.benzene.dens_mol_liq_comp_coeff_n1.fixed
        assert m.fs.props.benzene.dens_mol_liq_comp_coeff_n2.value == pytest.approx(
            -0.5596547795188646, rel=1e-10
        )
        assert m.fs.props.benzene.dens_mol_liq_comp_coeff_n2.fixed
        assert m.fs.props.benzene.dens_mol_liq_comp_coeff_n3.value == pytest.approx(
            14.872052666571532, rel=1e-10
        )
        assert m.fs.props.benzene.dens_mol_liq_comp_coeff_n3.fixed
        assert m.fs.props.benzene.dens_mol_liq_comp_coeff_n4.value == pytest.approx(
            -66.42959110979461, rel=1e-10
        )
        assert m.fs.props.benzene.dens_mol_liq_comp_coeff_n4.fixed
        assert m.fs.props.benzene.dens_mol_liq_comp_coeff_n5.value == pytest.approx(
            1158.1329856052375, rel=1e-10
        )
        assert m.fs.props.benzene.dens_mol_liq_comp_coeff_n5.fixed
        assert m.fs.props.benzene.dens_mol_liq_comp_coeff_n6.value == pytest.approx(
            -3128.774352224071, rel=1e-10
        )
        assert m.fs.props.benzene.dens_mol_liq_comp_coeff_n6.fixed

        assert m.fs.props.benzene.dens_mol_liq_comp_coeff_t1.value == pytest.approx(
            0.407, rel=1e-10
        )
        assert m.fs.props.benzene.dens_mol_liq_comp_coeff_t1.fixed
        assert m.fs.props.benzene.dens_mol_liq_comp_coeff_t2.value == pytest.approx(
            0.565, rel=1e-10
        )
        assert m.fs.props.benzene.dens_mol_liq_comp_coeff_t2.fixed
        assert m.fs.props.benzene.dens_mol_liq_comp_coeff_t3.value == pytest.approx(
            4.029, rel=1e-10
        )
        assert m.fs.props.benzene.dens_mol_liq_comp_coeff_t3.fixed
        assert m.fs.props.benzene.dens_mol_liq_comp_coeff_t4.value == pytest.approx(
            5.699, rel=1e-10
        )
        assert m.fs.props.benzene.dens_mol_liq_comp_coeff_t4.fixed
        assert m.fs.props.benzene.dens_mol_liq_comp_coeff_t5.value == pytest.approx(
            9.989, rel=1e-10
        )
        assert m.fs.props.benzene.dens_mol_liq_comp_coeff_t5.fixed
        assert m.fs.props.benzene.dens_mol_liq_comp_coeff_t6.value == pytest.approx(
            12.299, rel=1e-10
        )
        assert m.fs.props.benzene.dens_mol_liq_comp_coeff_t6.fixed

        # CoolProp results are non-ideal, which results in deviations
        for T in range(300, 401, 10):
            m.fs.state[0].temperature.fix(T)
            assert pytest.approx(
                CoolProp.PropsSI("DMOLAR", "T", T, "Q", 0, "benzene"), rel=5e-3
            ) == value(
                CoolPropWrapper.dens_mol_liq_comp.return_expression(
                    m.fs.state[0], m.fs.props.benzene, T * pyunits.K
                )
            )

    @pytest.mark.unit
    def test_enth_mol_liq(self, m):
        assert isinstance(m.fs.props.benzene.enth_mol_liq_comp_coeff_A0, Var)
        assert isinstance(m.fs.props.benzene.enth_mol_liq_comp_coeff_A1, Var)
        assert isinstance(m.fs.props.benzene.enth_mol_liq_comp_coeff_A2, Var)
        assert isinstance(m.fs.props.benzene.enth_mol_liq_comp_coeff_A3, Var)
        assert isinstance(m.fs.props.benzene.enth_mol_liq_comp_coeff_A4, Var)
        assert isinstance(m.fs.props.benzene.enth_mol_liq_comp_coeff_A5, Var)
        assert isinstance(m.fs.props.benzene.enth_mol_liq_comp_coeff_A6, Var)
        assert isinstance(m.fs.props.benzene.enth_mol_liq_comp_coeff_A7, Var)

        assert isinstance(m.fs.props.benzene.enth_mol_liq_comp_coeff_B0, Var)
        assert isinstance(m.fs.props.benzene.enth_mol_liq_comp_coeff_B1, Var)

        assert m.fs.props.benzene.enth_mol_liq_comp_coeff_A0.value == pytest.approx(
            -59920.08975429371, rel=1e-10
        )
        assert m.fs.props.benzene.enth_mol_liq_comp_coeff_A0.fixed
        assert m.fs.props.benzene.enth_mol_liq_comp_coeff_A1.value == pytest.approx(
            -435.3004255527287, rel=1e-10
        )
        assert m.fs.props.benzene.enth_mol_liq_comp_coeff_A1.fixed
        assert m.fs.props.benzene.enth_mol_liq_comp_coeff_A2.value == pytest.approx(
            5.686269274827653, rel=1e-10
        )
        assert m.fs.props.benzene.enth_mol_liq_comp_coeff_A2.fixed
        assert m.fs.props.benzene.enth_mol_liq_comp_coeff_A3.value == pytest.approx(
            -0.026551688405615677, rel=1e-10
        )
        assert m.fs.props.benzene.enth_mol_liq_comp_coeff_A3.fixed
        assert m.fs.props.benzene.enth_mol_liq_comp_coeff_A4.value == pytest.approx(
            7.158706209734695e-05, rel=1e-10
        )
        assert m.fs.props.benzene.enth_mol_liq_comp_coeff_A4.fixed
        assert m.fs.props.benzene.enth_mol_liq_comp_coeff_A5.value == pytest.approx(
            -1.1450245890590754e-07, rel=1e-10
        )
        assert m.fs.props.benzene.enth_mol_liq_comp_coeff_A5.fixed
        assert m.fs.props.benzene.enth_mol_liq_comp_coeff_A6.value == pytest.approx(
            1.0020646255928023e-10, rel=1e-10
        )
        assert m.fs.props.benzene.enth_mol_liq_comp_coeff_A6.fixed
        assert m.fs.props.benzene.enth_mol_liq_comp_coeff_A7.value == pytest.approx(
            -3.714155743547148e-14, rel=1e-10
        )
        assert m.fs.props.benzene.enth_mol_liq_comp_coeff_A7.fixed

        assert m.fs.props.benzene.enth_mol_liq_comp_coeff_B0.value == pytest.approx(
            1, rel=1e-10
        )
        assert m.fs.props.benzene.enth_mol_liq_comp_coeff_B0.fixed
        assert m.fs.props.benzene.enth_mol_liq_comp_coeff_B1.value == pytest.approx(
            -0.0017636924804910123, rel=1e-10
        )
        assert m.fs.props.benzene.enth_mol_liq_comp_coeff_B1.fixed

        assert m.fs.props.benzene.enth_mol_liq_comp_anchor.value == pytest.approx(
            54228.40275249632, rel=1e-10
        )
        assert m.fs.props.benzene.enth_mol_liq_comp_anchor.fixed

        for T in range(300, 401, 10):
            m.fs.state[0].temperature.fix(T)
            assert pytest.approx(
                CoolProp.PropsSI("Hmolar", "T", T, "Q", 0, "benzene"), rel=5e-4
            ) == value(
                CoolPropWrapper.enth_mol_liq_comp.return_expression(
                    m.fs.state[0], m.fs.props.benzene, T * pyunits.K
                )
            )

    @pytest.mark.unit
    def test_enth_mol_ig(self, m):
        assert isinstance(m.fs.props.benzene.enth_mol_ig_comp_coeff_A0, Var)
        assert isinstance(m.fs.props.benzene.enth_mol_ig_comp_coeff_A1, Var)
        assert isinstance(m.fs.props.benzene.enth_mol_ig_comp_coeff_A2, Var)
        assert isinstance(m.fs.props.benzene.enth_mol_ig_comp_coeff_A3, Var)
        assert isinstance(m.fs.props.benzene.enth_mol_ig_comp_coeff_A4, Var)
        assert isinstance(m.fs.props.benzene.enth_mol_ig_comp_coeff_A5, Var)
        assert isinstance(m.fs.props.benzene.enth_mol_ig_comp_coeff_A6, Var)
        assert isinstance(m.fs.props.benzene.enth_mol_ig_comp_coeff_A7, Var)

        assert isinstance(m.fs.props.benzene.enth_mol_ig_comp_coeff_B0, Var)
        assert isinstance(m.fs.props.benzene.enth_mol_ig_comp_coeff_B1, Var)

        assert m.fs.props.benzene.enth_mol_ig_comp_coeff_A0.value == pytest.approx(
            -72325.6503937073, rel=1e-10
        )
        assert m.fs.props.benzene.enth_mol_ig_comp_coeff_A0.fixed
        assert m.fs.props.benzene.enth_mol_ig_comp_coeff_A1.value == pytest.approx(
            2230.356360119971, rel=1e-10
        )
        assert m.fs.props.benzene.enth_mol_ig_comp_coeff_A1.fixed
        assert m.fs.props.benzene.enth_mol_ig_comp_coeff_A2.value == pytest.approx(
            -19.220224456239638, rel=1e-10
        )
        assert m.fs.props.benzene.enth_mol_ig_comp_coeff_A2.fixed
        assert m.fs.props.benzene.enth_mol_ig_comp_coeff_A3.value == pytest.approx(
            0.08579743322690633, rel=1e-10
        )
        assert m.fs.props.benzene.enth_mol_ig_comp_coeff_A3.fixed
        assert m.fs.props.benzene.enth_mol_ig_comp_coeff_A4.value == pytest.approx(
            -0.00022509135257932873, rel=1e-10
        )
        assert m.fs.props.benzene.enth_mol_ig_comp_coeff_A4.fixed
        assert m.fs.props.benzene.enth_mol_ig_comp_coeff_A5.value == pytest.approx(
            3.4964051978046773e-07, rel=1e-10
        )
        assert m.fs.props.benzene.enth_mol_ig_comp_coeff_A5.fixed
        assert m.fs.props.benzene.enth_mol_ig_comp_coeff_A6.value == pytest.approx(
            -2.9861697218437715e-10, rel=1e-10
        )
        assert m.fs.props.benzene.enth_mol_ig_comp_coeff_A6.fixed
        assert m.fs.props.benzene.enth_mol_ig_comp_coeff_A7.value == pytest.approx(
            1.0849858354638855e-13, rel=1e-10
        )
        assert m.fs.props.benzene.enth_mol_ig_comp_coeff_A7.fixed

        assert m.fs.props.benzene.enth_mol_ig_comp_coeff_B0.value == pytest.approx(
            1, rel=1e-10
        )
        assert m.fs.props.benzene.enth_mol_ig_comp_coeff_B0.fixed
        assert m.fs.props.benzene.enth_mol_ig_comp_coeff_B1.value == pytest.approx(
            -0.0017620752364379978, rel=1e-10
        )
        assert m.fs.props.benzene.enth_mol_ig_comp_coeff_B1.fixed

        for T in range(300, 401, 10):
            m.fs.state[0].temperature.fix(T)
            assert pytest.approx(
                CoolProp.PropsSI("Hmolar", "T", T, "Q", 1, "benzene"), rel=5e-4
            ) == value(
                CoolPropWrapper.enth_mol_ig_comp.return_expression(
                    m.fs.state[0], m.fs.props.benzene, T * pyunits.K
                )
            )

    @pytest.mark.unit
    def test_entr_mol_liq(self, m):
        assert isinstance(m.fs.props.benzene.entr_mol_liq_comp_coeff_A0, Var)
        assert isinstance(m.fs.props.benzene.entr_mol_liq_comp_coeff_A1, Var)
        assert isinstance(m.fs.props.benzene.entr_mol_liq_comp_coeff_A2, Var)
        assert isinstance(m.fs.props.benzene.entr_mol_liq_comp_coeff_A3, Var)
        assert isinstance(m.fs.props.benzene.entr_mol_liq_comp_coeff_A4, Var)
        assert isinstance(m.fs.props.benzene.entr_mol_liq_comp_coeff_A5, Var)
        assert isinstance(m.fs.props.benzene.entr_mol_liq_comp_coeff_A6, Var)
        assert isinstance(m.fs.props.benzene.entr_mol_liq_comp_coeff_A7, Var)

        assert isinstance(m.fs.props.benzene.entr_mol_liq_comp_coeff_B0, Var)
        assert isinstance(m.fs.props.benzene.entr_mol_liq_comp_coeff_B1, Var)

        assert m.fs.props.benzene.entr_mol_liq_comp_coeff_A0.value == pytest.approx(
            -367.28736141903704, rel=1e-10
        )
        assert m.fs.props.benzene.entr_mol_liq_comp_coeff_A0.fixed
        assert m.fs.props.benzene.entr_mol_liq_comp_coeff_A1.value == pytest.approx(
            1.9742974534821924, rel=1e-10
        )
        assert m.fs.props.benzene.entr_mol_liq_comp_coeff_A1.fixed
        assert m.fs.props.benzene.entr_mol_liq_comp_coeff_A2.value == pytest.approx(
            -0.0044315658996985936, rel=1e-10
        )
        assert m.fs.props.benzene.entr_mol_liq_comp_coeff_A2.fixed
        assert m.fs.props.benzene.entr_mol_liq_comp_coeff_A3.value == pytest.approx(
            1.3891733323299664e-06, rel=1e-10
        )
        assert m.fs.props.benzene.entr_mol_liq_comp_coeff_A3.fixed
        assert m.fs.props.benzene.entr_mol_liq_comp_coeff_A4.value == pytest.approx(
            2.2861650201840276e-08, rel=1e-10
        )
        assert m.fs.props.benzene.entr_mol_liq_comp_coeff_A4.fixed
        assert m.fs.props.benzene.entr_mol_liq_comp_coeff_A5.value == pytest.approx(
            -6.533821280556607e-11, rel=1e-10
        )
        assert m.fs.props.benzene.entr_mol_liq_comp_coeff_A5.fixed
        assert m.fs.props.benzene.entr_mol_liq_comp_coeff_A6.value == pytest.approx(
            7.560967575590022e-14, rel=1e-10
        )
        assert m.fs.props.benzene.entr_mol_liq_comp_coeff_A6.fixed
        assert m.fs.props.benzene.entr_mol_liq_comp_coeff_A7.value == pytest.approx(
            -3.327152569338899e-17, rel=1e-10
        )
        assert m.fs.props.benzene.entr_mol_liq_comp_coeff_A7.fixed

        assert m.fs.props.benzene.entr_mol_liq_comp_coeff_B0.value == pytest.approx(
            1, rel=1e-10
        )
        assert m.fs.props.benzene.entr_mol_liq_comp_coeff_B0.fixed
        assert m.fs.props.benzene.entr_mol_liq_comp_coeff_B1.value == pytest.approx(
            -0.0017637184112409214, rel=1e-10
        )
        assert m.fs.props.benzene.entr_mol_liq_comp_coeff_B1.fixed

        assert m.fs.props.benzene.entr_mol_liq_comp_anchor.value == pytest.approx(
            108.66827371730366, rel=1e-10
        )
        assert m.fs.props.benzene.entr_mol_liq_comp_anchor.fixed

        for T in range(300, 401, 10):
            m.fs.state[0].temperature.fix(T)
            assert pytest.approx(
                CoolProp.PropsSI("Smolar", "T", T, "Q", 0, "benzene"), rel=5e-4
            ) == value(
                CoolPropWrapper.entr_mol_liq_comp.return_expression(
                    m.fs.state[0], m.fs.props.benzene, T * pyunits.K
                )
            )

    @pytest.mark.unit
    def test_entr_mol_ig(self, m):
        assert isinstance(m.fs.props.benzene.entr_mol_ig_comp_coeff_A0, Var)
        assert isinstance(m.fs.props.benzene.entr_mol_ig_comp_coeff_A1, Var)
        assert isinstance(m.fs.props.benzene.entr_mol_ig_comp_coeff_A2, Var)
        assert isinstance(m.fs.props.benzene.entr_mol_ig_comp_coeff_A3, Var)
        assert isinstance(m.fs.props.benzene.entr_mol_ig_comp_coeff_A4, Var)
        assert isinstance(m.fs.props.benzene.entr_mol_ig_comp_coeff_A5, Var)
        assert isinstance(m.fs.props.benzene.entr_mol_ig_comp_coeff_A6, Var)
        assert isinstance(m.fs.props.benzene.entr_mol_ig_comp_coeff_A7, Var)

        assert isinstance(m.fs.props.benzene.entr_mol_ig_comp_coeff_B0, Var)
        assert isinstance(m.fs.props.benzene.entr_mol_ig_comp_coeff_B1, Var)

        assert m.fs.props.benzene.entr_mol_ig_comp_coeff_A0.value == pytest.approx(
            805.7903359684354, rel=1e-10
        )
        assert m.fs.props.benzene.entr_mol_ig_comp_coeff_A0.fixed
        assert m.fs.props.benzene.entr_mol_ig_comp_coeff_A1.value == pytest.approx(
            -7.054023595110686, rel=1e-10
        )
        assert m.fs.props.benzene.entr_mol_ig_comp_coeff_A1.fixed
        assert m.fs.props.benzene.entr_mol_ig_comp_coeff_A2.value == pytest.approx(
            0.026948845504297693, rel=1e-10
        )
        assert m.fs.props.benzene.entr_mol_ig_comp_coeff_A2.fixed
        assert m.fs.props.benzene.entr_mol_ig_comp_coeff_A3.value == pytest.approx(
            -4.872969272569453e-05, rel=1e-10
        )
        assert m.fs.props.benzene.entr_mol_ig_comp_coeff_A3.fixed
        assert m.fs.props.benzene.entr_mol_ig_comp_coeff_A4.value == pytest.approx(
            1.4263950361676005e-08, rel=1e-10
        )
        assert m.fs.props.benzene.entr_mol_ig_comp_coeff_A4.fixed
        assert m.fs.props.benzene.entr_mol_ig_comp_coeff_A5.value == pytest.approx(
            9.666423804354728e-11, rel=1e-10
        )
        assert m.fs.props.benzene.entr_mol_ig_comp_coeff_A5.fixed
        assert m.fs.props.benzene.entr_mol_ig_comp_coeff_A6.value == pytest.approx(
            -1.5555968001310246e-13, rel=1e-10
        )
        assert m.fs.props.benzene.entr_mol_ig_comp_coeff_A6.fixed
        assert m.fs.props.benzene.entr_mol_ig_comp_coeff_A7.value == pytest.approx(
            7.656482336490157e-17, rel=1e-10
        )
        assert m.fs.props.benzene.entr_mol_ig_comp_coeff_A7.fixed

        assert m.fs.props.benzene.entr_mol_ig_comp_coeff_B0.value == pytest.approx(
            1, rel=1e-10
        )
        assert m.fs.props.benzene.entr_mol_ig_comp_coeff_B0.fixed
        assert m.fs.props.benzene.entr_mol_ig_comp_coeff_B1.value == pytest.approx(
            -0.0017620273519612887, rel=1e-10
        )
        assert m.fs.props.benzene.entr_mol_ig_comp_coeff_B1.fixed

        for T in range(300, 401, 10):
            m.fs.state[0].temperature.fix(T)
            assert pytest.approx(
                CoolProp.PropsSI("Smolar", "T", T, "Q", 1, "benzene"), rel=5e-4
            ) == value(
                CoolPropWrapper.entr_mol_ig_comp.return_expression(
                    m.fs.state[0], m.fs.props.benzene, T * pyunits.K
                )
            )

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

        assert m.fs.props.benzene.pressure_sat_coeff_n1.value == pytest.approx(
            0.005561906558935796, rel=1e-10
        )
        assert m.fs.props.benzene.pressure_sat_coeff_n1.fixed
        assert m.fs.props.benzene.pressure_sat_coeff_n2.value == pytest.approx(
            -0.08662136922915314, rel=1e-10
        )
        assert m.fs.props.benzene.pressure_sat_coeff_n2.fixed
        assert m.fs.props.benzene.pressure_sat_coeff_n3.value == pytest.approx(
            -6.964182734154488, rel=1e-10
        )
        assert m.fs.props.benzene.pressure_sat_coeff_n3.fixed
        assert m.fs.props.benzene.pressure_sat_coeff_n4.value == pytest.approx(
            1.1249288132278856, rel=1e-10
        )
        assert m.fs.props.benzene.pressure_sat_coeff_n4.fixed
        assert m.fs.props.benzene.pressure_sat_coeff_n5.value == pytest.approx(
            -3.961859460597414, rel=1e-10
        )
        assert m.fs.props.benzene.pressure_sat_coeff_n5.fixed
        assert m.fs.props.benzene.pressure_sat_coeff_n6.value == pytest.approx(
            -13.106880507410812, rel=1e-10
        )
        assert m.fs.props.benzene.pressure_sat_coeff_n6.fixed

        assert m.fs.props.benzene.pressure_sat_coeff_t1.value == pytest.approx(
            0.037, rel=1e-10
        )
        assert m.fs.props.benzene.pressure_sat_coeff_t1.fixed
        assert m.fs.props.benzene.pressure_sat_coeff_t2.value == pytest.approx(
            0.505, rel=1e-10
        )
        assert m.fs.props.benzene.pressure_sat_coeff_t2.fixed
        assert m.fs.props.benzene.pressure_sat_coeff_t3.value == pytest.approx(
            1.014, rel=1e-10
        )
        assert m.fs.props.benzene.pressure_sat_coeff_t3.fixed
        assert m.fs.props.benzene.pressure_sat_coeff_t4.value == pytest.approx(
            1.469, rel=1e-10
        )
        assert m.fs.props.benzene.pressure_sat_coeff_t4.fixed
        assert m.fs.props.benzene.pressure_sat_coeff_t5.value == pytest.approx(
            3.711, rel=1e-10
        )
        assert m.fs.props.benzene.pressure_sat_coeff_t5.fixed
        assert m.fs.props.benzene.pressure_sat_coeff_t6.value == pytest.approx(
            12.647, rel=1e-10
        )
        assert m.fs.props.benzene.pressure_sat_coeff_t6.fixed

        for T in range(300, 401, 10):
            m.fs.state[0].temperature.fix(T)
            assert pytest.approx(
                CoolProp.PropsSI("P", "T", T, "Q", 0.5, "benzene"), rel=5e-4
            ) == value(m.fs.state[0].pressure_sat_comp["benzene"])


class TestVerifyExcessLiq(object):
    @pytest.fixture(scope="class")
    def m(self):
        # Clear cached components to ensure clean slate
        CoolPropWrapper.flush_cached_components()

        m = ConcreteModel()

        m.fs = FlowsheetBlock(dynamic=False)

        configuration = {
            # Specifying components
            "components": {
                "benzene": {
                    "type": Component,
                    "elemental_composition": {"H": 6, "C": 6},
                    "dens_mol_liq_comp": CoolPropWrapper,
                    "enth_mol_ig_comp": Constant,
                    "entr_mol_ig_comp": Constant,
                    "parameter_data": {
                        "mw": CoolPropWrapper,
                        "dens_mol_crit": CoolPropWrapper,
                        "pressure_crit": CoolPropWrapper,
                        "temperature_crit": CoolPropWrapper,
                        "omega": CoolPropWrapper,
                        "cp_mol_ig_comp_coeff": 0,
                        "enth_mol_form_ig_comp_ref": 0,
                        "entr_mol_form_ig_comp_ref": 0,
                    },
                }
            },
            # Specifying phases
            "phases": {
                "Liq": {
                    "type": LiquidPhase,
                    "equation_of_state": Cubic,
                    "equation_of_state_options": {"type": CubicType.PR},
                }
            },
            # Set base units of measurement
            "base_units": {
                "time": pyunits.s,
                "length": pyunits.m,
                "mass": pyunits.kg,
                "amount": pyunits.mol,
                "temperature": pyunits.K,
            },
            # Specifying state definition
            "state_definition": FTPx,
            "state_bounds": {
                "flow_mol": (0, 100, 1000, pyunits.mol / pyunits.s),
                "temperature": (273.15, 300, 500, pyunits.K),
                "pressure": (5e4, 1e5, 1e6, pyunits.Pa),
            },
            "pressure_ref": (101325, pyunits.Pa),
            "temperature_ref": (298.15, pyunits.K),
            "parameter_data": {"PR_kappa": {("benzene", "benzene"): 0.000}},
        }

        m.fs.props = GenericParameterBlock(**configuration)

        m.fs.state = m.fs.props.build_state_block([0], defined_state=True)

        m.fs.state[0].flow_mol.fix(1)
        m.fs.state[0].mole_frac_comp["benzene"].fix(1)

        return m

    @pytest.mark.integration
    def test_cubic_liquid(self, m):
        for P in range(1, 11):
            Td = CoolProp.PropsSI("T", "P", P * 1e5, "Q", 0.5, "PR::benzene")
            Tmin = CoolProp.PropsSI("TMIN", "benzene")

            m.fs.state[0].pressure.fix(P * 1e5)
            m.fs.state[0].temperature.fix(Tmin)

            m.fs.state.initialize()

            for T in arange(Tmin, Td, 10):
                print(P, T)
                m.fs.state[0].temperature.fix(T)

                results = solver.solve(m.fs)

                assert (
                    results.solver.termination_condition == TerminationCondition.optimal
                )
                assert results.solver.status == SolverStatus.ok

                # Check results
                assert pytest.approx(
                    CoolProp.PropsSI("Z", "T", T, "P", P * 1e5, "PR::benzene"), rel=1e-8
                ) == value(m.fs.state[0].compress_fact_phase["Liq"])

                assert pytest.approx(
                    CoolProp.PropsSI("DMOLAR", "T", T, "P", P * 1e5, "PR::benzene"),
                    rel=1e-6,
                ) == value(m.fs.state[0].dens_mol_phase["Liq"])

                assert pytest.approx(
                    CoolProp.PropsSI(
                        "HMOLAR_RESIDUAL", "T", T, "P", P * 1e5, "PR::benzene"
                    ),
                    rel=1e-6,
                ) == value(m.fs.state[0].enth_mol_phase["Liq"])

    @pytest.mark.integration
    def test_cubic_liquid_entr(self, m):
        Td = CoolProp.PropsSI("T", "P", 101325, "Q", 0.5, "PR::benzene")
        Tmin = CoolProp.PropsSI("TMIN", "benzene")

        for T in arange(Tmin, Td, 10):
            m.fs.state[0].pressure.fix(101325)
            m.fs.state[0].temperature.fix(T)

            m.fs.state.initialize()

            S0_CP = CoolProp.PropsSI("SMOLAR", "T", T, "P", 101325, "PR::benzene")
            S0_I = value(m.fs.state[0].entr_mol_phase["Liq"])

            for P in range(1, 11):
                print(T, P)
                m.fs.state[0].pressure.fix(P * 1e5)

                results = solver.solve(m.fs)

                assert (
                    results.solver.termination_condition == TerminationCondition.optimal
                )
                assert results.solver.status == SolverStatus.ok

                assert pytest.approx(
                    CoolProp.PropsSI("SMOLAR", "T", T, "P", P * 1e5, "PR::benzene")
                    - S0_CP,
                    rel=1e-4,
                ) == value(m.fs.state[0].entr_mol_phase["Liq"] - S0_I)


class TestVerifyExcessVap(object):
    @pytest.fixture(scope="class")
    def m(self):
        # Clear cached components to ensure clean slate
        CoolPropWrapper.flush_cached_components()

        m = ConcreteModel()

        m.fs = FlowsheetBlock(dynamic=False)

        # NBP reference state: h=0, s=0 for saturated liquid at 1 atmosphere
        CoolProp.set_reference_state("benzene", "NBP")

        # Get Tsat for benezene at 1 atm
        Tref = CoolProp.PropsSI("T", "P", 101325, "Q", 0.5, "PR::benzene")

        configuration = {
            # Specifying components
            "components": {
                "benzene": {
                    "type": Component,
                    "elemental_composition": {"H": 6, "C": 6},
                    "dens_mol_liq_comp": CoolPropWrapper,
                    "enth_mol_ig_comp": Constant,
                    "entr_mol_ig_comp": Constant,
                    "parameter_data": {
                        "mw": CoolPropWrapper,
                        "dens_mol_crit": CoolPropWrapper,
                        "pressure_crit": CoolPropWrapper,
                        "temperature_crit": CoolPropWrapper,
                        "omega": CoolPropWrapper,
                        "cp_mol_ig_comp_coeff": 0,
                        "enth_mol_form_ig_comp_ref": 0,
                        "entr_mol_form_ig_comp_ref": 0,
                    },
                }
            },
            # Specifying phases
            "phases": {
                "Vap": {
                    "type": VaporPhase,
                    "equation_of_state": Cubic,
                    "equation_of_state_options": {"type": CubicType.PR},
                }
            },
            # Set base units of measurement
            "base_units": {
                "time": pyunits.s,
                "length": pyunits.m,
                "mass": pyunits.kg,
                "amount": pyunits.mol,
                "temperature": pyunits.K,
            },
            # Specifying state definition
            "state_definition": FTPx,
            "state_bounds": {
                "flow_mol": (0, 100, 1000, pyunits.mol / pyunits.s),
                "temperature": (273.15, 300, 730, pyunits.K),
                "pressure": (5e4, 1e5, 1e6, pyunits.Pa),
            },
            "pressure_ref": (101325, pyunits.Pa),
            "temperature_ref": (Tref, pyunits.K),
            "parameter_data": {"PR_kappa": {("benzene", "benzene"): 0.000}},
        }

        m.fs.props = GenericParameterBlock(**configuration)

        m.fs.state = m.fs.props.build_state_block([0], defined_state=True)

        m.fs.state[0].flow_mol.fix(1)
        m.fs.state[0].mole_frac_comp["benzene"].fix(1)

        return m

    @pytest.mark.integration
    def test_cubic_vapor(self, m):
        for P in range(1, 11):
            Tb = CoolProp.PropsSI("T", "P", P * 1e5, "Q", 0.5, "PR::benzene")
            Tmax = CoolProp.PropsSI("TMAX", "benzene")

            m.fs.state[0].pressure.fix(P * 1e5)
            m.fs.state[0].temperature.fix(Tb + 1)

            m.fs.state.initialize()

            for T in arange(Tb + 1, Tmax, 10):
                m.fs.state[0].temperature.fix(T)

                results = solver.solve(m.fs)

                assert (
                    results.solver.termination_condition == TerminationCondition.optimal
                )
                assert results.solver.status == SolverStatus.ok

                # Check results
                assert pytest.approx(
                    CoolProp.PropsSI("Z", "T", T, "P", P * 1e5, "PR::benzene"), rel=1e-8
                ) == value(m.fs.state[0].compress_fact_phase["Vap"])

                assert pytest.approx(
                    CoolProp.PropsSI("DMOLAR", "T", T, "P", P * 1e5, "PR::benzene"),
                    rel=1e-6,
                ) == value(m.fs.state[0].dens_mol_phase["Vap"])

                assert pytest.approx(
                    CoolProp.PropsSI(
                        "HMOLAR_RESIDUAL", "T", T, "P", P * 1e5, "PR::benzene"
                    ),
                    rel=1e-6,
                ) == value(m.fs.state[0].enth_mol_phase["Vap"])

    @pytest.mark.integration
    def test_cubic_vapor_entr(self, m):
        Tb = CoolProp.PropsSI("T", "P", 1e6, "Q", 0.5, "PR::benzene")
        Tmax = CoolProp.PropsSI("TMAX", "benzene")

        for T in arange(Tb + 1, Tmax, 10):
            m.fs.state[0].pressure.fix(101325)
            m.fs.state[0].temperature.fix(T)

            m.fs.state.initialize()

            S0_CP = CoolProp.PropsSI("SMOLAR", "T", T, "P", 101325, "PR::benzene")
            S0_I = value(m.fs.state[0].entr_mol_phase["Vap"])

            for P in range(1, 11):
                print(T, P)
                m.fs.state[0].pressure.fix(P * 1e5)

                results = solver.solve(m.fs)

                assert (
                    results.solver.termination_condition == TerminationCondition.optimal
                )
                assert results.solver.status == SolverStatus.ok

                assert pytest.approx(
                    CoolProp.PropsSI("SMOLAR", "T", T, "P", P * 1e5, "PR::benzene")
                    - S0_CP,
                    rel=1e-4,
                ) == value(m.fs.state[0].entr_mol_phase["Vap"] - S0_I)
