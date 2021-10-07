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
from idaes.power_generation.carbon_capture.mea_solvent_system.properties.MEA_PCC \
 import configuration_vapor_abs, configuration_vapor_reg


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
            assert isinstance(model.params.get_component(i), Component)

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

        assert model.params.O2.mw.value == 31.999E-3
        assert model.params.O2.pressure_crit.value == 50.464e5
        assert model.params.O2.temperature_crit.value == 154.599

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
        for i in model.params.component_list:
            assert i in ['H2O','CO2']
            assert isinstance(model.params.get_component(i), Component)

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

        assert model.params.CO2.mw.value == 44.009E-3
        assert model.params.CO2.pressure_crit.value == 73.773e5
        assert model.params.CO2.temperature_crit.value == 304.128

        assert_units_consistent(model)

