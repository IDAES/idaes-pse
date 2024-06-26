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
from idaes.core.util import to_json, from_json
from idaes.core.util.scaling import jacobian_cond
from idaes.core.scaling import AutoScaler
from idaes.core.scaling import CustomScalerBase


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
    from_json(m, fname=fname)
    # Make sure we have no suffixes loaded
    assert not hasattr(m.fs.unit, "scaling_factor")

    # Autoscale variables by magnitude
    scaler = AutoScaler()
    scaler.variables_by_magnitude(m)
    assert hasattr(m.fs.unit, "scaling_factor")

    return m


@pytest.mark.integration
def test_verify_model_load(gibbs):
    for v in gibbs.component_data_objects(ctype=Var, descend_into=True):
        assert v.value is not None


@pytest.mark.integration
def test_autoscale_L2_norm(gibbs):
    scaler = AutoScaler()
    scaler.constraints_by_jacobian_norm(gibbs, norm=2)

    scaled = jacobian_cond(gibbs, scaled=True)

    assert scaled == pytest.approx(2510.945, rel=1e-5)


@pytest.mark.integration
def test_autoscale_L1_norm(gibbs):
    scaler = AutoScaler()
    scaler.constraints_by_jacobian_norm(gibbs, norm=1)

    scaled = jacobian_cond(gibbs, scaled=True)

    assert scaled == pytest.approx(2986.994, rel=1e-5)


@pytest.mark.integration
def test_nominal_magnitude_harmonic(gibbs):
    scaler = CustomScalerBase()

    for c in gibbs.component_data_objects(ctype=Constraint, descend_into=True):
        scaler.scale_constraint_by_nominal_value(c, scheme="harmonic_mean")

    scaled = jacobian_cond(gibbs, scaled=True)

    assert scaled == pytest.approx(2.83719e12, rel=1e-5)


@pytest.mark.integration
def test_nominal_magnitude_max(gibbs):
    scaler = CustomScalerBase()

    for c in gibbs.component_data_objects(ctype=Constraint, descend_into=True):
        scaler.scale_constraint_by_nominal_value(c, scheme="inverse_maximum")

    scaled = jacobian_cond(gibbs, scaled=True)

    assert scaled == pytest.approx(1.77413e11, rel=1e-5)


@pytest.mark.integration
def test_nominal_magnitude_min(gibbs):
    scaler = CustomScalerBase()

    for c in gibbs.component_data_objects(ctype=Constraint, descend_into=True):
        scaler.scale_constraint_by_nominal_value(c, scheme="inverse_minimum")

    scaled = jacobian_cond(gibbs, scaled=True)

    assert scaled == pytest.approx(5.59871e12, rel=1e-5)


@pytest.mark.integration
def test_nominal_magnitude_inv_mean(gibbs):
    scaler = CustomScalerBase()

    for c in gibbs.component_data_objects(ctype=Constraint, descend_into=True):
        scaler.scale_constraint_by_nominal_value(c, scheme="inverse_sum")

    scaled = jacobian_cond(gibbs, scaled=True)

    assert scaled == pytest.approx(8.76894e10, rel=1e-5)


@pytest.mark.integration
def test_nominal_magnitude_inv_rss(gibbs):
    scaler = CustomScalerBase()

    for c in gibbs.component_data_objects(ctype=Constraint, descend_into=True):
        scaler.scale_constraint_by_nominal_value(c, scheme="inverse_root_sum_squared")

    scaled = jacobian_cond(gibbs, scaled=True)

    assert scaled == pytest.approx(1.30134e11, rel=1e-5)


@pytest.mark.integration
def test_scale_constraint_by_nominal_jacobian_norm(gibbs):
    scaler = CustomScalerBase()

    for c in gibbs.component_data_objects(ctype=Constraint, descend_into=True):
        scaler.scale_constraint_by_nominal_jacobian_norm(c)

    scaled = jacobian_cond(gibbs, scaled=True)
    print(scaled)
    assert scaled == pytest.approx(1.73937e16, rel=1e-5)
    assert False


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
    scaler.variables_by_magnitude(m)
    scaler.constraints_by_jacobian_norm(m)

    solver = get_solver("ipopt_v2", writer_config={"scale_model": True})
    results = solver.solve(m, tee=True)

    # Check for optimal solution
    assert_optimal_termination(results)

    to_json(m, fname=fname, human_read=True)
