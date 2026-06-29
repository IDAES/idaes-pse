#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Tests for solvent reboiler unit model.
Authors: Andrew Lee
"""

import pytest
from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    Constraint,
    Param,
    units,
    value,
)
from pyomo.util.check_units import assert_units_consistent, assert_units_equivalent

from idaes.core import FlowsheetBlock
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.testing import initialization_tester
from idaes.core.util import scaling as iscale
from idaes.core.solvers import get_solver

from idaes.models_extra.column_models.solvent_reboiler import SolventReboiler
from idaes.models_extra.column_models.properties.MEA_solvent import (
    configuration as aqueous_mea,
)
from idaes.models_extra.column_models.properties.MEA_vapor import flue_gas
from idaes.core.util.testing import assert_solution_equivalent

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


def call_initialization_old_API(model):
    initialization_tester(
        model,
        liquid_state_args={
            "pressure": 183700,
            "temperature": 393.8,
            "flow_mol": 74.33,
            "mole_frac_comp": {"CO2": 0.0285, "H2O": 0.8491, "MEA": 0.1224},
        },
    )


def call_initialization_new_API(model):
    prefix = "liquid_phase.properties_out[0]."
    liquid_guess = {
        prefix + "pressure": 183700,
        prefix + "temperature": 393.8,
        prefix + "flow_mol": 74.33,
    }
    mole_frac_guesses = {"CO2": 0.0285, "H2O": 0.8491, "MEA": 0.1224}
    for j, val in mole_frac_guesses.items():
        liquid_guess[prefix + f"mole_frac_comp[{j}]"] = val

    initializer_object = model.fs.unit.default_initializer()
    initializer_object.initialize(model.fs.unit, initial_guesses=liquid_guess)


def vapor_flow_specify_model(model):
    model.fs.unit.vapor_reboil.flow_mol[0].fix(9.56)


def heat_duty_specify_model(model):
    model.fs.unit.heat_duty.fix(430.61e3)


_rel = 1e-5
_abs = None

vapor_flow_expected_results = {
    "bottoms.flow_mol": {
        0: (74.33, _rel, _abs),
    },
    "bottoms.mole_frac_comp": {
        (0, "CO2"): (0.0285221, _rel, _abs),
        (0, "MEA"): (0.122455, _rel, _abs),
        (0, "H2O"): (0.849023, _rel, _abs),
    },
    "bottoms.pressure": {0: (183700, _rel, _abs)},
    "bottoms.temperature": {0: (393.773, _rel, _abs)},
    "vapor_reboil.flow_mol": {
        0: (9.56, _rel, _abs),
    },
    "vapor_reboil.mole_frac_comp": {
        (0, "CO2"): (0.0643063, _rel, _abs),
        (0, "H2O"): (0.935693, _rel, _abs),
        (0, "N2"): (0, None, 1e-8),
        (0, "O2"): (0, None, 1e-8),
    },
    "vapor_reboil.pressure": {
        0: (183700, _rel, _abs),
    },
    "vapor_reboil.temperature": {
        0: (393.773, _rel, _abs),
    },
    "heat_duty": {0: (420983, _rel, _abs)},
}

heat_duty_expected_results = {
    "bottoms.flow_mol": {
        0: (74.1048, _rel, _abs),
    },
    "bottoms.mole_frac_comp": {
        (0, "CO2"): (0.0285059, _rel, _abs),
        (0, "MEA"): (0.122827, _rel, _abs),
        (0, "H2O"): (0.848667, _rel, _abs),
    },
    "bottoms.pressure": {
        0: (183700, _rel, _abs),
    },
    "bottoms.temperature": {
        0: (393.810, _rel, _abs),
    },
    "vapor_reboil.flow_mol": {
        0: (9.7852, _rel, _abs),
    },
    "vapor_reboil.mole_frac_comp": {
        (0, "CO2"): (0.063605, _rel, _abs),
        (0, "H2O"): (0.936395, _rel, _abs),
        (0, "N2"): (0, None, 1e-8),
        (0, "O2"): (0, None, 1e-8),
    },
    "vapor_reboil.pressure": {
        0: (183700, _rel, _abs),
    },
    "vapor_reboil.temperature": {
        0: (393.810, _rel, _abs),
    },
}


vapor_flow_specification = {
    "specify_model": vapor_flow_specify_model,
    "expected_results": vapor_flow_expected_results,
}
heat_duty_specification = {
    "specify_model": heat_duty_specify_model,
    "expected_results": heat_duty_expected_results,
}


# -----------------------------------------------------------------------------
@pytest.mark.parametrize(
    "initialization_method",
    [call_initialization_old_API, call_initialization_new_API],
    scope="class",
)
@pytest.mark.parametrize(
    "specification_type",
    [vapor_flow_specification, heat_duty_specification],
    scope="class",
)
class TestAbsorber(object):
    @pytest.fixture(scope="class")
    def model(self, initialization_method, specification_type):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.liquid_properties = GenericParameterBlock(**aqueous_mea)
        m.fs.vapor_properties = GenericParameterBlock(**flue_gas)

        m.fs.unit = SolventReboiler(
            liquid_property_package=m.fs.liquid_properties,
            vapor_property_package=m.fs.vapor_properties,
        )

        m.fs.unit.inlet.flow_mol[0].fix(83.89)
        m.fs.unit.inlet.temperature[0].fix(392.5)
        m.fs.unit.inlet.pressure[0].fix(183700)
        m.fs.unit.inlet.mole_frac_comp[0, "CO2"].fix(0.0326)
        m.fs.unit.inlet.mole_frac_comp[0, "H2O"].fix(0.8589)
        m.fs.unit.inlet.mole_frac_comp[0, "MEA"].fix(0.1085)

        specification_type["specify_model"](m)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, initialization_method, specification_type, model):

        assert hasattr(model.fs.unit, "inlet")
        assert len(model.fs.unit.inlet.vars) == 4
        assert hasattr(model.fs.unit.inlet, "flow_mol")
        assert hasattr(model.fs.unit.inlet, "mole_frac_comp")
        assert hasattr(model.fs.unit.inlet, "temperature")
        assert hasattr(model.fs.unit.inlet, "pressure")

        assert hasattr(model.fs.unit, "bottoms")
        assert len(model.fs.unit.bottoms.vars) == 4
        assert hasattr(model.fs.unit.bottoms, "flow_mol")
        assert hasattr(model.fs.unit.bottoms, "mole_frac_comp")
        assert hasattr(model.fs.unit.bottoms, "temperature")
        assert hasattr(model.fs.unit.bottoms, "pressure")

        assert hasattr(model.fs.unit, "vapor_reboil")
        assert len(model.fs.unit.vapor_reboil.vars) == 4
        assert hasattr(model.fs.unit.vapor_reboil, "flow_mol")
        assert hasattr(model.fs.unit.vapor_reboil, "mole_frac_comp")
        assert hasattr(model.fs.unit.vapor_reboil, "temperature")
        assert hasattr(model.fs.unit.vapor_reboil, "pressure")

        assert isinstance(model.fs.unit.unit_material_balance, Constraint)
        assert isinstance(model.fs.unit.unit_enthalpy_balance, Constraint)
        assert isinstance(model.fs.unit.unit_temperature_equality, Constraint)
        assert isinstance(model.fs.unit.unit_pressure_balance, Constraint)
        assert isinstance(model.fs.unit.zero_flow_param, Param)

        assert number_variables(model.fs.unit) == 84
        assert number_total_constraints(model.fs.unit) == 77
        assert number_unused_variables(model.fs.unit) == 0

    @pytest.mark.component
    def test_units(self, initialization_method, specification_type, model):
        assert_units_consistent(model)
        assert_units_equivalent(model.fs.unit.heat_duty[0], units.W)

    @pytest.mark.unit
    def test_dof(self, initialization_method, specification_type, model):
        assert degrees_of_freedom(model) == 0

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, initialization_method, specification_type, model):
        initialization_method(model)

    # @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, initialization_method, specification_type, model):
        results = solver.solve(model)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, initialization_method, specification_type, model):
        assert_solution_equivalent(
            model.fs.unit, specification_type["expected_results"]
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, initialization_method, specification_type, model):
        assert (
            abs(
                value(
                    model.fs.unit.inlet.flow_mol[0]
                    - model.fs.unit.bottoms.flow_mol[0]
                    - model.fs.unit.vapor_reboil.flow_mol[0]
                )
            )
            <= 1e-6
        )

        assert (
            abs(
                value(
                    model.fs.unit.inlet.flow_mol[0]
                    * model.fs.unit.inlet.mole_frac_comp[0, "CO2"]
                    - model.fs.unit.bottoms.flow_mol[0]
                    * model.fs.unit.bottoms.mole_frac_comp[0, "CO2"]
                    - model.fs.unit.vapor_reboil.flow_mol[0]
                    * model.fs.unit.vapor_reboil.mole_frac_comp[0, "CO2"]
                )
            )
            <= 1e-6
        )
        assert (
            abs(
                value(
                    model.fs.unit.inlet.flow_mol[0]
                    * model.fs.unit.inlet.mole_frac_comp[0, "H2O"]
                    - model.fs.unit.bottoms.flow_mol[0]
                    * model.fs.unit.bottoms.mole_frac_comp[0, "H2O"]
                    - model.fs.unit.vapor_reboil.flow_mol[0]
                    * model.fs.unit.vapor_reboil.mole_frac_comp[0, "H2O"]
                )
            )
            <= 1e-6
        )
        assert (
            abs(
                value(
                    model.fs.unit.inlet.flow_mol[0]
                    * model.fs.unit.inlet.mole_frac_comp[0, "MEA"]
                    - model.fs.unit.bottoms.flow_mol[0]
                    * model.fs.unit.bottoms.mole_frac_comp[0, "MEA"]
                )
            )
            <= 1e-6
        )

        assert (
            abs(
                value(
                    model.fs.unit.liquid_phase.properties_in[0]._enthalpy_flow_term[
                        "Liq"
                    ]
                    - model.fs.unit.liquid_phase.properties_out[0]._enthalpy_flow_term[
                        "Liq"
                    ]
                    - model.fs.unit.vapor_phase[0]._enthalpy_flow_term["Vap"]
                    + model.fs.unit.heat_duty[0]
                )
            )
            <= 1e-6
        )
