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
Author: Andrew Lee, Alejandro Garciadiego
"""
import pytest
from pyomo.environ import (ConcreteModel,
                           Set,
                           SolverStatus,
                           TerminationCondition,
                           value,
                           Var,
                           units as pyunits)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import Component, FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester
from idaes.core.util import get_default_solver


from idaes.generic_models.properties.core.generic.generic_property import (
        GenericParameterBlock)

from idaes.generic_models.unit_models import Flash
from idaes.generic_models.properties.core.state_definitions import FTPx
from idaes.generic_models.properties.core.phase_equil import smooth_VLE

from idaes.generic_models.properties.core.examples.CO2_bmimPF6_PR \
    import configuration


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_default_solver()


# Test for configuration dictionaries with parameters from Properties of Gases
# and liquids 4th edition
class TestParamBlock(object):
    @pytest.mark.unit
    def test_build(self):
        model = ConcreteModel()

        model.param = GenericParameterBlock(default=configuration)

        assert isinstance(model.param.phase_list, Set)
        assert len(model.param.phase_list) == 2
        for i in model.param.phase_list:
            assert i in ["Liq", "Vap"]
        assert model.param.Liq.is_liquid_phase()
        assert model.param.Vap.is_vapor_phase()

        assert isinstance(model.param.component_list, Set)
        assert len(model.param.component_list) == 2
        for i in model.param.component_list:
            assert i in ['bmimPF6',
                         'carbon_dioxide']
            assert isinstance(model.param.get_component(i), Component)

        assert isinstance(model.param._phase_component_set, Set)
        assert len(model.param._phase_component_set) == 3
        for i in model.param._phase_component_set:
            assert i in [("Liq", "bmimPF6"), ("Liq", "carbon_dioxide"),
                         ("Vap", "carbon_dioxide")]

        assert model.param.config.state_definition == FTPx

        assert model.param.config.state_bounds == {
            "flow_mol": (0, 100, 1000, pyunits.mol/pyunits.s),
            "temperature": (10, 300, 500, pyunits.K),
            "pressure": (5e-4, 1e5, 1e10, pyunits.Pa)}

        assert model.param.config.phase_equilibrium_state == {
            ("Vap", "Liq"): smooth_VLE}

        assert isinstance(model.param.phase_equilibrium_idx, Set)
        assert len(model.param.phase_equilibrium_idx) == 1
        for i in model.param.phase_equilibrium_idx:
            assert i in ["PE1"]

        assert model.param.phase_equilibrium_list == {
            "PE1": {"carbon_dioxide": ("Vap", "Liq")}}

        assert model.param.pressure_ref.value == 101325
        assert model.param.temperature_ref.value == 298.15

        assert model.param.bmimPF6.mw.value == 284.18E-3
        assert model.param.bmimPF6.pressure_crit.value == 24e5
        assert model.param.bmimPF6.temperature_crit.value == 860

        assert model.param.carbon_dioxide.mw.value == 44.010E-3
        assert model.param.carbon_dioxide.pressure_crit.value == 71.8e5
        assert model.param.carbon_dioxide.temperature_crit.value == 304.1

        assert_units_consistent(model)


class TestStateBlock(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()

        model.fs = FlowsheetBlock(default={"dynamic": False})

        model.fs.param = GenericParameterBlock(default=configuration)

        model.fs.props = model.fs.param.build_state_block(
                [1],
                default={"defined_state": True})

        # Fix state
        model.fs.props[1].flow_mol.fix(1)
        model.fs.props[1].temperature.fix(298.15)
        model.fs.props[1].pressure.fix(1214713.75)
        model.fs.props[1].mole_frac_comp["carbon_dioxide"].fix(0.2)
        model.fs.props[1].mole_frac_comp["bmimPF6"].fix(0.8)

        assert degrees_of_freedom(model.fs.props[1]) == 0

        return model

    @pytest.mark.unit
    def test_build(self, model):
        # Check state variable values and bounds
        assert isinstance(model.fs.props[1].flow_mol, Var)
        assert value(model.fs.props[1].flow_mol) == 1
        assert model.fs.props[1].flow_mol.ub == 1000
        assert model.fs.props[1].flow_mol.lb == 0

        assert isinstance(model.fs.props[1].pressure, Var)
        assert value(model.fs.props[1].pressure) == 1214713.75
        assert model.fs.props[1].pressure.ub == 1e10
        assert model.fs.props[1].pressure.lb == 5e-4

        assert isinstance(model.fs.props[1].temperature, Var)
        assert value(model.fs.props[1].temperature) == 298.15
        assert model.fs.props[1].temperature.ub == 500
        assert model.fs.props[1].temperature.lb == 10

        assert isinstance(model.fs.props[1].mole_frac_comp, Var)
        assert len(model.fs.props[1].mole_frac_comp) == 2
        assert value(model.fs.props[1].mole_frac_comp["carbon_dioxide"]) == 0.2
        assert value(model.fs.props[1].mole_frac_comp["bmimPF6"]) == 0.8

        assert_units_consistent(model)

    @pytest.mark.unit
    def test_define_state_vars(self, model):
        sv = model.fs.props[1].define_state_vars()

        assert len(sv) == 4
        for i in sv:
            assert i in ["flow_mol",
                         "mole_frac_comp",
                         "temperature",
                         "pressure"]

    @pytest.mark.unit
    def test_define_port_members(self, model):
        sv = model.fs.props[1].define_state_vars()

        assert len(sv) == 4
        for i in sv:
            assert i in ["flow_mol",
                         "mole_frac_comp",
                         "temperature",
                         "pressure"]

    @pytest.mark.unit
    def test_define_display_vars(self, model):
        sv = model.fs.props[1].define_display_vars()

        assert len(sv) == 4
        for i in sv:
            assert i in ["Total Molar Flowrate",
                         "Total Mole Fraction",
                         "Temperature",
                         "Pressure"]

    @pytest.mark.unit
    def test_unit_dof(self, model):
        model.fs.unit = Flash(default={"property_package": model.fs.param,
                                       "has_heat_transfer": False,
                                       "has_pressure_change": False})
        # Fix state
        model.fs.unit.inlet.flow_mol.fix(1)
        model.fs.unit.inlet.temperature.fix(200.00)
        model.fs.unit.inlet.pressure.fix(101325)
        model.fs.unit.inlet.mole_frac_comp[0, "carbon_dioxide"].fix(1/2)
        model.fs.unit.inlet.mole_frac_comp[0, "bmimPF6"].fix(1/2)

        assert degrees_of_freedom(model.fs.unit) == 0

    @pytest.mark.component
    def test_initialize(self, model):

        initialization_tester(model, dof=0)

    @pytest.mark.component
    def test_solve(self, model):
        results = solver.solve(model)

        # Check for optimal solution
        assert results.solver.termination_condition == \
            TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.component
    def test_solution(self, model):
        # Check phase equilibrium results
        assert model.fs.unit.liq_outlet.mole_frac_comp[
            0, "carbon_dioxide"].value == \
            pytest.approx(0.3119, abs=1e-4)
        assert model.fs.unit.vap_outlet.mole_frac_comp[
            0, "carbon_dioxide"].value == \
            pytest.approx(1.0000, abs=1e-4)
        assert (model.fs.unit.vap_outlet.flow_mol[0].value /
                model.fs.unit.liq_outlet.flow_mol[0].value) == \
            pytest.approx(0.37619, abs=1e-4)

    @pytest.mark.ui
    @pytest.mark.unit
    def test_report(self, model):
        model.fs.props[1].report()
