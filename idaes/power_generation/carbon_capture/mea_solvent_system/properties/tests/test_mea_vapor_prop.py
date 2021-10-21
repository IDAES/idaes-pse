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
from pyomo.environ import (ConcreteModel,
                           Set,
                           SolverStatus,
                           TerminationCondition,
                           value,
                           Var,
                           units as pyunits)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import Component
from idaes.core.util import get_solver
from idaes.generic_models.properties.core.generic.generic_property import (
    GenericParameterBlock)
from idaes.generic_models.properties.core.state_definitions import FTPx
from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              fixed_variables_set,
                                              activated_constraints_set)

from idaes.power_generation.carbon_capture.mea_solvent_system.properties.mea_vapor_prop \
 import configuration_vapor_abs, configuration_vapor_reg

from idaes.core.components import Component
from pyomo.common.config import ConfigBlock,ConfigValue
from idaes.core import (declare_process_block_class, ProcessBlockData)


# Set up logger
#_log = idaeslog.getLogger(__name__)

@declare_process_block_class("Component", block_class=Component)
class ComponentData(ProcessBlockData):
    """
    Component type for additional properties.
    """
    CONFIG = ConfigBlock()

    CONFIG.declare("visc_d_comp", ConfigValue(
        description="Method to use to calculate vapor phase viscosity",
        doc="Method to use to vapor viscosity. Users "))
    def build(self):
        super().build()

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


# Test for configuration dictionaries
class TestParamBlock(object):
    @pytest.mark.unit
    def test_abs_build(self):
        model = ConcreteModel()
        model.params = GenericParameterBlock(default=configuration_vapor_abs)

        assert isinstance(model.params.phase_list, Set)
        assert len(model.params.phase_list) == 1
        for i in model.params.phase_list:
            assert i in ["Vap"]
        assert model.params.Vap.is_vapor_phase()

        assert isinstance(model.params.component_list, Set)
        assert len(model.params.component_list) == 4
        for i in model.params.component_list:
            assert i in ['N2',
                         'O2',
                         'H2O',
                         'CO2']
        assert isinstance(model.params._phase_component_set, Set)
        assert len(model.params._phase_component_set) == 4
        for i in model.params._phase_component_set:
            assert i in [
                ("Vap", "N2"),
                ("Vap", "O2"), ("Vap", "H2O"), ("Vap", "CO2")]

        assert model.params.config.state_definition == FTPx

        assert model.params.pressure_ref.value == 101325
        assert model.params.temperature_ref.value == 298.15

        assert model.params.N2.mw.value == 28.014E-3
        assert model.params.N2.pressure_crit.value == 33.958E5
        assert model.params.N2.temperature_crit.value == 126.192
        assert model.params.N2.diffus_volume.value == 18.5

        assert model.params.O2.mw.value == 31.999E-3
        assert model.params.O2.pressure_crit.value == 50.464e5
        assert model.params.O2.temperature_crit.value == 154.599
        assert model.params.O2.diffus_volume.value == 16.3

        assert model.params.H2O.mw.value == 18.015E-3
        assert model.params.H2O.pressure_crit.value == 220.64e5
        assert model.params.H2O.temperature_crit.value == 647.096

        assert model.params.CO2.mw.value == 44.009E-3
        assert model.params.CO2.pressure_crit.value == 73.773e5
        assert model.params.CO2.temperature_crit.value == 304.128

        assert_units_consistent(model)

    @pytest.mark.unit
    def test_reg_build(self):
        model = ConcreteModel()
        model.params = GenericParameterBlock(default=configuration_vapor_reg)

        assert isinstance(model.params.phase_list, Set)
        assert len(model.params.phase_list) == 1
        for i in model.params.phase_list:
            assert i in ["Vap"]
        assert model.params.Vap.is_vapor_phase()

        assert isinstance(model.params.component_list, Set)
        assert len(model.params.component_list) == 2
        assert isinstance(model.params._phase_component_set, Set)
        assert len(model.params._phase_component_set) == 2
        for i in model.params._phase_component_set:
            assert i in [("Vap", "H2O"), ("Vap", "CO2")]

        assert model.params.config.state_definition == FTPx

        assert model.params.pressure_ref.value == 101325
        assert model.params.temperature_ref.value == 298.15

        assert model.params.H2O.mw.value == 18.015E-3
        assert model.params.H2O.pressure_crit.value == 220.64e5
        assert model.params.H2O.temperature_crit.value == 647.096
        assert model.params.H2O.diffus_volume.value == 13.1

        assert model.params.CO2.mw.value == 44.009E-3
        assert model.params.CO2.pressure_crit.value == 73.773e5
        assert model.params.CO2.temperature_crit.value == 304.128
        assert model.params.CO2.diffus_volume.value == 26.7

        assert_units_consistent(model)

class TestStateBlock(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.params = GenericParameterBlock(default=configuration_vapor_abs)

        model.props = model.params.build_state_block(
                [1],
                default={"defined_state": True})
        model.props[1].flow_mol.fix(22)
        model.props[1].temperature.fix(313)
        model.props[1].pressure.fix(101325)
        model.props[1].mole_frac_comp["CO2"].fix(0.09183)
        model.props[1].mole_frac_comp["H2O"].fix(0.08495)
        model.props[1].mole_frac_comp["N2"].fix(0.73900)
        model.props[1].mole_frac_comp["O2"].fix(0.08422)

        return model

    @pytest.mark.unit
    def test_build(self, model):
        # Check state variable values and bounds
        assert isinstance(model.props[1].flow_mol, Var)
        assert value(model.props[1].flow_mol) == 22
        assert model.props[1].flow_mol.ub == 1000
        assert model.props[1].flow_mol.lb == 0

        assert isinstance(model.props[1].pressure, Var)
        assert value(model.props[1].pressure) == 101325
        assert model.props[1].pressure.ub == 1e6
        assert model.props[1].pressure.lb == 5e4

        assert isinstance(model.props[1].temperature, Var)
        assert value(model.props[1].temperature) == 313
        assert model.props[1].temperature.ub == 450
        assert model.props[1].temperature.lb == 273.15

        assert isinstance(model.props[1].mole_frac_comp, Var)
        assert len(model.props[1].mole_frac_comp) == 4

        assert_units_consistent(model)

    @pytest.mark.component
    def test_initialize(self):
        model = ConcreteModel()
        model.params = GenericParameterBlock(default=configuration_vapor_abs)

        model.props = model.params.build_state_block(
                [1],
                default={"defined_state": True})
        model.props[1].flow_mol.fix(22)
        model.props[1].temperature.fix(313)
        model.props[1].pressure.fix(101325)
        model.props[1].mole_frac_comp["CO2"].fix(0.09183)
        model.props[1].mole_frac_comp["H2O"].fix(0.08495)
        model.props[1].mole_frac_comp["N2"].fix(0.73900)
        model.props[1].mole_frac_comp["O2"].fix(0.08422)

        assert degrees_of_freedom(model.props[1]) == 0

        orig_fixed_vars = fixed_variables_set(model)
        orig_act_consts = activated_constraints_set(model)

        model.props.initialize(optarg={'tol': 1e-6})

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
        assert results.solver.termination_condition == \
            TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.component
    def test_solution(self, model):
        # Check Transport transporties
        model = ConcreteModel()
        model.params = GenericParameterBlock(default=configuration_vapor_abs)

        model.props = model.params.build_state_block(
                [1],
                default={"defined_state": True})
        model.props[1].flow_mol.fix(22)
        model.props[1].temperature.fix(313)
        model.props[1].pressure.fix(101325)
        model.props[1].mole_frac_comp["CO2"].fix(0.09183)
        model.props[1].mole_frac_comp["H2O"].fix(0.08495)
        model.props[1].mole_frac_comp["N2"].fix(0.73900)
        model.props[1].mole_frac_comp["O2"].fix(0.08422)

        results = solver.solve(model)

        # viscosity
        assert value(model.props[1].visc_d_comp["Vap", "N2"]) == \
            pytest.approx(1.8372e-05, abs=1e-4)
        assert value(model.props[1].visc_d_comp["Vap", "O2"]) == \
            pytest.approx(2.1312e-05, abs=1e-4)
        assert value(model.props[1].visc_d_comp["Vap", "H2O"]) == \
            pytest.approx(1.03378e-05, abs=1e-4)
        assert value(model.props[1].visc_d_comp["Vap", "CO2"]) == \
            pytest.approx(1.5675e-05, abs=1e-4)
        assert value(model.props[1].visc_d_phase["Vap"]) == \
            pytest.approx(1.761358e-05, abs=1e-4)

        # Thermal conductivity
        assert value(model.props[1].therm_cond_comp["Vap", "N2"]) == \
            pytest.approx(0.0265, abs=1e-4)
        assert value(model.props[1].therm_cond_comp["Vap", "O2"]) == \
            pytest.approx(0.02764, abs=1e-4)
        assert value(model.props[1].therm_cond_comp["Vap", "H2O"]) == \
            pytest.approx(0.01904, abs=1e-4)
        assert value(model.props[1].therm_cond_comp["Vap", "CO2"]) == \
            pytest.approx(0.01763, abs=1e-4)
        assert value(model.props[1].therm_cond_phase["Vap"]) == \
            pytest.approx(0.0249, abs=1e-4)

        # Diffusivity
        assert value(model.props[1].diffus_phase_comp["Vap", "N2"]) == \
            pytest.approx(2.18375e-05, abs=1e-4)
        assert value(model.props[1].diffus_phase_comp["Vap", "O2"]) == \
            pytest.approx(2.23035e-05, abs=1e-4)
        assert value(model.props[1].diffus_phase_comp["Vap", "H2O"]) == \
            pytest.approx(2.75097e-05, abs=1e-4)
        assert value(model.props[1].diffus_phase_comp["Vap", "CO2"]) == \
            pytest.approx(1.81086e-05, abs=1e-4)

    @pytest.mark.unit
    def test_report(self, model):
        model.props[1].report()

