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
Authors: Anuja Deshpande, Andrew Lee
"""
import pytest
import numpy as np

from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    Set,
    value,
    Var,
    units as pyunits,
)
from pyomo.util.check_units import assert_units_consistent
from pyomo.common.unittest import assertStructuredAlmostEqual

from idaes.core import Component

from idaes.core import FlowsheetBlock

from idaes.models.unit_models import Flash

from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    fixed_variables_set,
    activated_constraints_set,
)
from idaes.core.solvers import get_solver

from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)

from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.phase_equil import SmoothVLE

from idaes.models.properties.modular_properties.examples.CO2_H2O_Ideal_VLE import (
    configuration,
)


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


def _as_quantity(x):
    unit = pyunits.get_units(x)
    if unit is None:
        unit = pyunits.dimensionless
    return value(x) * unit._get_pint_unit()


class TestParamBlock(object):
    @pytest.mark.unit
    def test_build(self):
        model = ConcreteModel()
        model.params = GenericParameterBlock(**configuration)

        assert isinstance(model.params.phase_list, Set)
        assert len(model.params.phase_list) == 2
        for i in model.params.phase_list:
            assert i in ["Liq", "Vap"]
        assert model.params.Liq.is_liquid_phase()
        assert model.params.Vap.is_vapor_phase()

        assert isinstance(model.params.component_list, Set)
        assert len(model.params.component_list) == 2
        for i in model.params.component_list:
            assert i in ["H2O", "CO2"]
            assert isinstance(model.params.get_component(i), Component)

        assert isinstance(model.params._phase_component_set, Set)
        assert len(model.params._phase_component_set) == 3
        for i in model.params._phase_component_set:
            assert i in [("Liq", "H2O"), ("Vap", "H2O"), ("Vap", "CO2")]

        assert model.params.config.state_definition == FTPx

        assertStructuredAlmostEqual(
            model.params.config.state_bounds,
            {
                "flow_mol": (0, 10, 20, pyunits.mol / pyunits.s),
                "temperature": (273.15, 323.15, 1000, pyunits.K),
                "pressure": (5e4, 108900, 1e7, pyunits.Pa),
            },
            item_callback=_as_quantity,
        )

        assert model.params.config.phase_equilibrium_state == {
            ("Vap", "Liq"): SmoothVLE
        }

        assert isinstance(model.params.phase_equilibrium_idx, Set)
        assert len(model.params.phase_equilibrium_idx) == 1
        for i in model.params.phase_equilibrium_idx:
            assert i in ["PE1"]

        assert model.params.phase_equilibrium_list == {"PE1": {"H2O": ("Vap", "Liq")}}

        assert model.params.pressure_ref.value == 101325
        assert model.params.temperature_ref.value == 298.15

        assert model.params.H2O.mw.value == 18.0153e-3
        assert model.params.H2O.pressure_crit.value == 220.64e5
        assert model.params.H2O.temperature_crit.value == 647

        assert model.params.CO2.mw.value == 44.0095e-3
        assert model.params.CO2.pressure_crit.value == 73.825e5
        assert model.params.CO2.temperature_crit.value == 304.23

        assert_units_consistent(model)


class TestStateBlock(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.params = GenericParameterBlock(**configuration)

        model.props = model.params.build_state_block([1], defined_state=True)

        model.props[1].flow_mol.fix(10)
        model.props[1].temperature.fix(323.15)
        model.props[1].pressure.fix(108900)
        model.props[1].mole_frac_comp["H2O"].fix(0.5)
        model.props[1].mole_frac_comp["CO2"].fix(0.5)

        return model

    @pytest.mark.unit
    def test_build(self, model):
        # Check state variable values and bounds
        assert isinstance(model.props[1].flow_mol, Var)
        assert value(model.props[1].flow_mol) == 10
        assert model.props[1].flow_mol.ub == 20
        assert model.props[1].flow_mol.lb == 0

        assert isinstance(model.props[1].pressure, Var)
        assert value(model.props[1].pressure) == 108900
        assert model.props[1].pressure.ub == 1e7
        assert model.props[1].pressure.lb == 5e4

        assert isinstance(model.props[1].temperature, Var)
        assert value(model.props[1].temperature) == 323.15
        assert model.props[1].temperature.ub == 1000
        assert model.props[1].temperature.lb == 273.15

        assert isinstance(model.props[1].mole_frac_comp, Var)
        assert len(model.props[1].mole_frac_comp) == 2
        assert value(model.props[1].mole_frac_comp["H2O"]) == 0.5
        assert value(model.props[1].mole_frac_comp["CO2"]) == 0.5

        assert_units_consistent(model)

    @pytest.mark.unit
    def test_define_state_vars(self, model):
        sv = model.props[1].define_state_vars()

        assert len(sv) == 4
        for i in sv:
            assert i in ["flow_mol", "mole_frac_comp", "temperature", "pressure"]

    @pytest.mark.unit
    def test_define_port_members(self, model):
        sv = model.props[1].define_state_vars()

        assert len(sv) == 4
        for i in sv:
            assert i in ["flow_mol", "mole_frac_comp", "temperature", "pressure"]

    @pytest.mark.unit
    def test_define_display_vars(self, model):
        sv = model.props[1].define_display_vars()

        assert len(sv) == 4
        for i in sv:
            assert i in [
                "Total Molar Flowrate",
                "Total Mole Fraction",
                "Temperature",
                "Pressure",
            ]

    @pytest.mark.unit
    def test_dof(self, model):
        assert degrees_of_freedom(model.props[1]) == 0

    @pytest.mark.component
    def test_initialize(self, model):
        orig_fixed_vars = fixed_variables_set(model)
        orig_act_consts = activated_constraints_set(model)

        model.props.initialize(optarg={"tol": 1e-6})

        assert degrees_of_freedom(model) == 0

        fin_fixed_vars = fixed_variables_set(model)
        fin_act_consts = activated_constraints_set(model)

        assert len(fin_act_consts) == len(orig_act_consts)
        assert len(fin_fixed_vars) == len(orig_fixed_vars)

        for c in fin_act_consts:
            assert c in orig_act_consts
        for v in fin_fixed_vars:
            assert v in orig_fixed_vars

    @pytest.mark.component
    def test_solve(self, model):
        results = solver.solve(model)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, model):
        # Check phase equilibrium results
        assert model.props[1].mole_frac_phase_comp["Liq", "H2O"].value == pytest.approx(
            1, abs=1e-4
        )
        assert model.props[1].mole_frac_phase_comp["Vap", "H2O"].value == pytest.approx(
            0.11501, abs=1e-4
        )
        assert model.props[1].mole_frac_phase_comp["Vap", "CO2"].value == pytest.approx(
            0.88499, abs=1e-4
        )
        assert model.props[1].phase_frac["Vap"].value == pytest.approx(
            0.56498, abs=1e-4
        )

    @pytest.mark.integration
    def test_temp_swing(self):
        # Create a flash model with the CO2-H2O property package
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = GenericParameterBlock(**configuration)
        m.fs.flash = Flash(property_package=m.fs.properties)

        # Fix inlet stream state variables
        m.fs.flash.inlet.flow_mol.fix(9.89433124673833)  # mol/s
        m.fs.flash.inlet.mole_frac_comp[0, "CO2"].fix(0.13805801934749645)
        m.fs.flash.inlet.mole_frac_comp[0, "H2O"].fix(0.8619419806525035)
        m.fs.flash.inlet.pressure.fix(183430)  # Pa
        m.fs.flash.inlet.temperature.fix(396.79057912844183)  # K

        # Fix flash and its outlet conditions
        m.fs.flash.deltaP.fix(0)
        m.fs.flash.vap_outlet.temperature.fix(313.15)

        # Initialize the flash model
        m.fs.flash.initialize()

        # Create a dictionary of expected solution for flash outlet temperature
        # sweep
        temp_range = list(np.linspace(313, 396))

        expected_vapor_frac = [
            0.14388,
            0.14445,
            0.14508,
            0.14576,
            0.1465,
            0.14731,
            0.14818,
            0.14913,
            0.15017,
            0.15129,
            0.15251,
            0.15384,
            0.15528,
            0.15685,
            0.15856,
            0.16042,
            0.16245,
            0.16467,
            0.16709,
            0.16974,
            0.17265,
            0.17584,
            0.17935,
            0.18323,
            0.18751,
            0.19226,
            0.19755,
            0.20346,
            0.21008,
            0.21755,
            0.22601,
            0.23565,
            0.24673,
            0.25956,
            0.27456,
            0.29229,
            0.31354,
            0.33942,
            0.37157,
            0.4125,
            0.46628,
            0.53993,
            0.64678,
            0.81547,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
        ]

        expected_heat_duty = [
            -396276.55508,
            -394859.42896,
            -393421.28226,
            -391960.40602,
            -390474.94512,
            -388962.88119,
            -387422.01297,
            -385849.93371,
            -384244.00483,
            -382601.32549,
            -380918.69691,
            -379192.58068,
            -377419.04976,
            -375593.73054,
            -373711.73419,
            -371767.57480,
            -369755.07135,
            -367667.22968,
            -365496.09940,
            -363232.59958,
            -360866.30485,
            -358385.18097,
            -355775.25563,
            -353020.20490,
            -350100.82935,
            -346994.38367,
            -343673.70980,
            -340106.10300,
            -336251.80960,
            -332062.00898,
            -327476.06061,
            -322417.68382,
            -316789.55435,
            -310465.49473,
            -303278.90949,
            -295005.17406,
            -285333.93764,
            -273823.89005,
            -259825.51107,
            -242341.82570,
            -219760.16114,
            -189290.02362,
            -145647.55666,
            -77469.59283,
            -3219.88910,
            -2631.66067,
            -2043.09220,
            -1454.17760,
            -864.91585,
            -275.30623,
        ]

        outvals = zip(expected_vapor_frac, expected_heat_duty)
        expected_sol = dict(zip(temp_range, outvals))

        # Solve the model for a range of flash outlet temperatures
        # Perform flash outlet temperature sweep and test the solution
        for t in temp_range:
            m.fs.flash.vap_outlet.temperature.fix(t)
            res = solver.solve(m)
            assert res.solver.termination_condition == "optimal"
            frac = value(m.fs.flash.vap_outlet.flow_mol[0]) / value(
                m.fs.flash.inlet.flow_mol[0]
            )
            assert frac == pytest.approx(expected_sol[t][0], abs=1e-4)
            hduty = value(m.fs.flash.heat_duty[0])
            assert hduty == pytest.approx(expected_sol[t][1], rel=1e-4)

    @pytest.mark.unit
    def test_report(self, model):
        model.props[1].report()
