#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Integration tests for constraint scaling by expected magnitude in Custom Scaling.

Tests will use the Gibbs Reactor example as the test case as it is very
poorly scaled in native form.

Author: Andrew Lee
"""

import os
import pytest

from pyomo.environ import (
    assert_optimal_termination,
    ConcreteModel,
    Constraint,
    Suffix,
    TransformationFactory,
    value,
    Var,
    units,
)

from idaes.core import FlowsheetBlock
from idaes.models.unit_models.gibbs_reactor import GibbsReactor
from idaes.models.properties.activity_coeff_models.methane_combustion_ideal import (
    MethaneParameterBlock as MethaneCombustionParameterBlock,
)
from idaes.core.util.testing import PhysicalParameterTestBlock, initialization_tester
from idaes.core.solvers import get_solver
from idaes.core.util import to_json, from_json, StoreSpec
from idaes.core.util.scaling import jacobian_cond
from idaes.core.scaling import AutoScaler, CustomScalerBase, set_scaling_factor


FILENAME = "gibbs_solution.json"
local_path = os.path.dirname(os.path.realpath(__file__))
fname = os.path.join(local_path, FILENAME)


def build_model():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = MethaneCombustionParameterBlock()

    m.fs.unit = GibbsReactor(
        property_package=m.fs.properties,
        has_heat_transfer=True,
        has_pressure_change=True,
    )

    m.fs.unit.inlet.flow_mol[0].fix(230.0)
    m.fs.unit.inlet.mole_frac_comp[0, "H2"].fix(0.0435)
    m.fs.unit.inlet.mole_frac_comp[0, "N2"].fix(0.6522)
    m.fs.unit.inlet.mole_frac_comp[0, "O2"].fix(0.1739)
    m.fs.unit.inlet.mole_frac_comp[0, "CO2"].fix(1e-5)
    m.fs.unit.inlet.mole_frac_comp[0, "CH4"].fix(0.1304)
    m.fs.unit.inlet.mole_frac_comp[0, "CO"].fix(1e-5)
    m.fs.unit.inlet.mole_frac_comp[0, "H2O"].fix(1e-5)
    m.fs.unit.inlet.mole_frac_comp[0, "NH3"].fix(1e-5)
    m.fs.unit.inlet.temperature[0].fix(1500.0)
    m.fs.unit.inlet.pressure[0].fix(101325.0)

    m.fs.unit.outlet.temperature[0].fix(2844.38)
    m.fs.unit.deltaP.fix(0)

    return m


@pytest.fixture
def gibbs():
    m = build_model()

    # Load solution
    from_json(m, fname=fname, wts=StoreSpec().value())
    # Make sure we have no suffixes loaded
    assert not hasattr(m.fs.unit, "scaling_factor")

    # Autoscale variables by magnitude
    scaler = AutoScaler()
    scaler.scale_variables_by_magnitude(m)
    assert hasattr(m.fs.unit, "scaling_factor")

    return m


@pytest.mark.integration
def test_verify_model_load(gibbs):
    for v in gibbs.component_data_objects(ctype=Var, descend_into=True):
        assert v.value is not None


@pytest.mark.integration
def test_autoscale_L2_norm(gibbs):
    scaler = AutoScaler()
    scaler.scale_constraints_by_jacobian_norm(gibbs, norm=2)

    scaled = jacobian_cond(gibbs, scaled=True)

    assert scaled == pytest.approx(2510.945, rel=1e-5)


@pytest.mark.integration
def test_autoscale_L1_norm(gibbs):
    scaler = AutoScaler()
    scaler.scale_constraints_by_jacobian_norm(gibbs, norm=1)

    scaled = jacobian_cond(gibbs, scaled=True)

    assert scaled == pytest.approx(2986.994, rel=1e-5)


@pytest.mark.integration
def test_nominal_magnitude_harmonic(gibbs):
    scaler = CustomScalerBase()

    for c in gibbs.component_data_objects(ctype=Constraint, descend_into=True):
        scaler.scale_constraint_by_nominal_value(c, scheme="harmonic_mean")

    assert jacobian_cond(gibbs, scaled=True) == pytest.approx(2.83944e12, rel=1e-5)


@pytest.mark.integration
def test_nominal_magnitude_inv_max(gibbs):
    scaler = CustomScalerBase()

    for c in gibbs.component_data_objects(ctype=Constraint, descend_into=True):
        scaler.scale_constraint_by_nominal_value(c, scheme="inverse_maximum")

    assert jacobian_cond(gibbs, scaled=True) == pytest.approx(784576, rel=1e-5)


@pytest.mark.integration
def test_nominal_magnitude_inv_min(gibbs):
    scaler = CustomScalerBase()

    for c in gibbs.component_data_objects(ctype=Constraint, descend_into=True):
        scaler.scale_constraint_by_nominal_value(c, scheme="inverse_minimum")

    assert jacobian_cond(gibbs, scaled=True) == pytest.approx(5.601e12, rel=1e-5)


@pytest.mark.integration
def test_nominal_magnitude_inv_sum(gibbs):
    scaler = CustomScalerBase()

    for c in gibbs.component_data_objects(ctype=Constraint, descend_into=True):
        scaler.scale_constraint_by_nominal_value(c, scheme="inverse_sum")

    assert jacobian_cond(gibbs, scaled=True) == pytest.approx(1501632, rel=1e-5)


@pytest.mark.integration
def test_nominal_magnitude_inv_rss(gibbs):
    scaler = CustomScalerBase()

    for c in gibbs.component_data_objects(ctype=Constraint, descend_into=True):
        scaler.scale_constraint_by_nominal_value(c, scheme="inverse_root_sum_squared")

    assert jacobian_cond(gibbs, scaled=True) == pytest.approx(959994, rel=1e-5)


@pytest.mark.integration
def test_scale_constraint_by_nominal_derivative_2norm_perfect_information(gibbs):
    scaler = CustomScalerBase()

    for c in gibbs.component_data_objects(ctype=Constraint, descend_into=True):
        scaler.scale_constraint_by_nominal_derivative_norm(c)

    scaled = jacobian_cond(gibbs, scaled=True)
    assert scaled == pytest.approx(3.07419e06, rel=1e-5)


@pytest.mark.integration
def test_scale_constraint_by_nominal_derivative_2norm_imperfect_information():
    # Build a fresh model with no scaling factors
    model = build_model()

    # Set imperfect scaling factors for all variables, representing an initial "best-guess"
    # Feed states are known exactly - set scaling based on these
    set_scaling_factor(
        model.fs.unit.control_volume.properties_in[0.0].flow_mol, 1 / 230
    )
    set_scaling_factor(
        model.fs.unit.control_volume.properties_in[0.0].flow_mol_phase, 1 / 230
    )  # Only 1 phase, so we "know" this
    set_scaling_factor(
        model.fs.unit.control_volume.properties_in[0.0].mole_frac_comp["H2"], 1 / 0.0435
    )
    set_scaling_factor(
        model.fs.unit.control_volume.properties_in[0.0].mole_frac_comp["N2"], 1 / 0.6522
    )
    set_scaling_factor(
        model.fs.unit.control_volume.properties_in[0.0].mole_frac_comp["O2"], 1 / 0.1739
    )
    set_scaling_factor(
        model.fs.unit.control_volume.properties_in[0.0].mole_frac_comp["CO2"], 1e5
    )
    set_scaling_factor(
        model.fs.unit.control_volume.properties_in[0.0].mole_frac_comp["CH4"],
        1 / 0.1304,
    )
    set_scaling_factor(
        model.fs.unit.control_volume.properties_in[0.0].mole_frac_comp["CO"], 1e5
    )
    set_scaling_factor(
        model.fs.unit.control_volume.properties_in[0.0].mole_frac_comp["H2O"], 1e5
    )
    set_scaling_factor(
        model.fs.unit.control_volume.properties_in[0.0].mole_frac_comp["NH3"], 1e5
    )
    set_scaling_factor(
        model.fs.unit.control_volume.properties_in[0.0].temperature, 1 / 1500
    )
    set_scaling_factor(model.fs.unit.control_volume.properties_in[0.0].pressure, 1e-5)
    # Assume user does not know anything about enthalpy

    # Best guesses for unit model and outlet state conditions
    set_scaling_factor(model.fs.unit.control_volume.heat[0.0], 1e-6)

    set_scaling_factor(model.fs.unit.control_volume.properties_out[0.0].flow_mol, 1e-2)
    set_scaling_factor(
        model.fs.unit.control_volume.properties_out[0.0].flow_mol_phase, 1e-2
    )  # Only 1 phase, so we "know" this
    # N2 is inert, so will be order 0.1, assume CH4 and H2 are near-totally consumed, assume most O2 consumed
    # Assume moderate amounts of CO2 and H2O, small amounts of CO, trace NH3 NH3
    set_scaling_factor(
        model.fs.unit.control_volume.properties_out[0.0].mole_frac_comp["H2"], 1e4
    )
    set_scaling_factor(
        model.fs.unit.control_volume.properties_out[0.0].mole_frac_comp["N2"], 10
    )
    set_scaling_factor(
        model.fs.unit.control_volume.properties_out[0.0].mole_frac_comp["O2"], 1e2
    )
    set_scaling_factor(
        model.fs.unit.control_volume.properties_out[0.0].mole_frac_comp["CO2"], 10
    )
    set_scaling_factor(
        model.fs.unit.control_volume.properties_out[0.0].mole_frac_comp["CH4"], 1e4
    )
    set_scaling_factor(
        model.fs.unit.control_volume.properties_out[0.0].mole_frac_comp["CO"], 1e3
    )
    set_scaling_factor(
        model.fs.unit.control_volume.properties_out[0.0].mole_frac_comp["H2O"], 10
    )
    set_scaling_factor(
        model.fs.unit.control_volume.properties_out[0.0].mole_frac_comp["NH3"], 1e4
    )
    set_scaling_factor(
        model.fs.unit.control_volume.properties_out[0.0].temperature, 1e-3
    )
    set_scaling_factor(model.fs.unit.control_volume.properties_out[0.0].pressure, 1e-5)

    scaler = CustomScalerBase()
    for c in model.component_data_objects(ctype=Constraint, descend_into=True):
        scaler.scale_constraint_by_nominal_derivative_norm(c)

    scaled = jacobian_cond(model, scaled=True)
    assert scaled == pytest.approx(1.06128e11, rel=1e-5)


@pytest.mark.integration
def test_scale_constraint_by_nominal_derivative_1norm_perfect_information(gibbs):
    scaler = CustomScalerBase()

    for c in gibbs.component_data_objects(ctype=Constraint, descend_into=True):
        scaler.scale_constraint_by_nominal_derivative_norm(c, norm=1)

    scaled = jacobian_cond(gibbs, scaled=True)
    assert scaled == pytest.approx(2.060153e06, rel=1e-5)


@pytest.mark.integration
def test_scale_constraint_by_nominal_derivative_clean_up(gibbs):
    # Confirm that the scaler did not change any values in the model
    # Build a fresh model with no scaling factors or initial values
    model = build_model()
    # Clone the model to use as a reference case
    refmodel = model.clone()

    # Set imperfect scaling factors for all variables, representing an initial "best-guess"
    # Feed states are known exactly - set scaling based on these
    set_scaling_factor(
        model.fs.unit.control_volume.properties_in[0.0].flow_mol, 1 / 230
    )
    set_scaling_factor(
        model.fs.unit.control_volume.properties_in[0.0].flow_mol_phase, 1 / 230
    )  # Only 1 phase, so we "know" this
    set_scaling_factor(
        model.fs.unit.control_volume.properties_in[0.0].mole_frac_comp["H2"], 1 / 0.0435
    )
    set_scaling_factor(
        model.fs.unit.control_volume.properties_in[0.0].mole_frac_comp["N2"], 1 / 0.6522
    )
    set_scaling_factor(
        model.fs.unit.control_volume.properties_in[0.0].mole_frac_comp["O2"], 1 / 0.1739
    )
    set_scaling_factor(
        model.fs.unit.control_volume.properties_in[0.0].mole_frac_comp["CO2"], 1e5
    )
    set_scaling_factor(
        model.fs.unit.control_volume.properties_in[0.0].mole_frac_comp["CH4"],
        1 / 0.1304,
    )
    set_scaling_factor(
        model.fs.unit.control_volume.properties_in[0.0].mole_frac_comp["CO"], 1e5
    )
    set_scaling_factor(
        model.fs.unit.control_volume.properties_in[0.0].mole_frac_comp["H2O"], 1e5
    )
    set_scaling_factor(
        model.fs.unit.control_volume.properties_in[0.0].mole_frac_comp["NH3"], 1e5
    )
    set_scaling_factor(
        model.fs.unit.control_volume.properties_in[0.0].temperature, 1 / 1500
    )
    set_scaling_factor(model.fs.unit.control_volume.properties_in[0.0].pressure, 1e-5)
    # Assume user does not know anything about enthalpy

    # Best guesses for unit model and outlet state conditions
    set_scaling_factor(model.fs.unit.control_volume.heat[0.0], 1e-6)

    set_scaling_factor(model.fs.unit.control_volume.properties_out[0.0].flow_mol, 1e-2)
    set_scaling_factor(
        model.fs.unit.control_volume.properties_out[0.0].flow_mol_phase, 1e-2
    )  # Only 1 phase, so we "know" this
    # N2 is inert, so will be order 0.1, assume CH4 and H2 are near-totally consumed, assume most O2 consumed
    # Assume moderate amounts of CO2 and H2O, small amounts of CO, trace NH3
    set_scaling_factor(
        model.fs.unit.control_volume.properties_out[0.0].mole_frac_comp["H2"], 1e4
    )
    set_scaling_factor(
        model.fs.unit.control_volume.properties_out[0.0].mole_frac_comp["N2"], 10
    )
    set_scaling_factor(
        model.fs.unit.control_volume.properties_out[0.0].mole_frac_comp["O2"], 1e2
    )
    set_scaling_factor(
        model.fs.unit.control_volume.properties_out[0.0].mole_frac_comp["CO2"], 10
    )
    set_scaling_factor(
        model.fs.unit.control_volume.properties_out[0.0].mole_frac_comp["CH4"], 1e4
    )
    set_scaling_factor(
        model.fs.unit.control_volume.properties_out[0.0].mole_frac_comp["CO"], 1e3
    )
    set_scaling_factor(
        model.fs.unit.control_volume.properties_out[0.0].mole_frac_comp["H2O"], 10
    )
    set_scaling_factor(
        model.fs.unit.control_volume.properties_out[0.0].mole_frac_comp["NH3"], 1e4
    )
    set_scaling_factor(
        model.fs.unit.control_volume.properties_out[0.0].temperature, 1e-3
    )
    set_scaling_factor(model.fs.unit.control_volume.properties_out[0.0].pressure, 1e-5)

    scaler = CustomScalerBase()
    for c in model.component_data_objects(ctype=Constraint, descend_into=True):
        scaler.scale_constraint_by_nominal_derivative_norm(c)

    # Iterate overall all variables in model and compare value to reference model
    for v in model.component_data_objects(ctype=Var, descend_into=True):
        refv = refmodel.find_component(v.name)

        assert v.value == refv.value


if __name__ == "__main__":
    # Run file to regenerate solution json
    m = build_model()

    initialization_tester(
        m,
        optarg={"tol": 1e-6},
        state_args={
            "temperature": 2844.38,
            "pressure": 101325.0,
            "flow_mol": 251.05,
            "mole_frac_comp": {
                "CH4": 1e-5,
                "CO": 0.0916,
                "CO2": 0.0281,
                "H2": 0.1155,
                "H2O": 0.1633,
                "N2": 0.5975,
                "NH3": 1e-5,
                "O2": 0.0067,
            },
        },
    )

    scaler = AutoScaler()
    scaler.scale_variables_by_magnitude(m)
    scaler.scale_constraints_by_jacobian_norm(m)

    solver = get_solver("ipopt_v2", writer_config={"scale_model": True})
    results = solver.solve(m, tee=True)

    # Check for optimal solution
    assert_optimal_termination(results)

    to_json(m, fname=fname, human_read=True, wts=StoreSpec().value())
