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
Tests for Scaling Profiler

Author: Andrew Lee
"""

from io import StringIO
import os
import pytest

from pyomo.environ import ConcreteModel, Constraint, value, Var

from idaes.core import FlowsheetBlock
from idaes.models.unit_models.gibbs_reactor import GibbsReactor
from idaes.models.properties.activity_coeff_models.methane_combustion_ideal import (
    MethaneParameterBlock as MethaneCombustionParameterBlock,
)
from idaes.core.util import from_json, StoreSpec
from idaes.core.scaling import set_scaling_factor
from idaes.core.scaling.scaler_profiling import ScalingProfiler


def demo_model():
    m = ConcreteModel()

    m.v1 = Var(initialize=2)
    m.v1.fix()
    m.v2 = Var(initialize=4)
    m.v3 = Var(initialize=-6)

    m.c1 = Constraint(expr=m.v2 == m.v1**2)
    m.c2 = Constraint(expr=0 == m.v1 + m.v2 + m.v3)

    return m


def demo_scaling(model):
    set_scaling_factor(model.v2, 0.5)


def demo_pertubration(model):
    model.v1.fix(3)


class TestScalingProfiler:
    @pytest.mark.unit
    def test_init_class(self):
        sp = ScalingProfiler(
            build_model=demo_model,
            user_scaling=demo_scaling,
            perturb_state=demo_pertubration,
        )

        assert sp._build_model is demo_model
        assert sp._user_scaling is demo_scaling
        assert sp._perturb_state is demo_pertubration
        assert sp._solver is not None
        assert sp._scaling_methods is not None

    @pytest.mark.unit
    def test_init_class_w_args(self):
        sp = ScalingProfiler(
            build_model=demo_model,
            user_scaling=demo_scaling,
            perturb_state=demo_pertubration,
            solver="foo",
            scaling_methods="bar",
        )

        assert sp._build_model is demo_model
        assert sp._user_scaling is demo_scaling
        assert sp._perturb_state is demo_pertubration
        assert sp._solver == "foo"
        assert sp._scaling_methods == "bar"

    @pytest.mark.unit
    def test_scale_vars_user(self):
        sp = ScalingProfiler(
            build_model=demo_model,
            user_scaling=demo_scaling,
            perturb_state=demo_pertubration,
        )

        model = sp._build_model()
        sp._scale_vars(model)

        # Should apply user scaling which will only scale v2
        assert len(model.scaling_factor) == 1
        assert model.scaling_factor[model.v2] == 0.5

    @pytest.mark.unit
    def test_scale_vars_auto(self):
        sp = ScalingProfiler(
            build_model=demo_model,
            user_scaling=demo_scaling,
            perturb_state=demo_pertubration,
        )

        model = sp._build_model()
        sp._scale_vars(model, perfect=True)

        # Should apply auto scaling which will scale all vars
        assert len(model.scaling_factor) == 3
        assert model.scaling_factor[model.v1] == 0.5
        assert model.scaling_factor[model.v2] == 0.25
        assert model.scaling_factor[model.v3] == pytest.approx(1 / 6, rel=1e-8)

    @pytest.mark.unit
    def test_apply_scaling_no_method(self):
        sp = ScalingProfiler(
            build_model=demo_model,
            user_scaling=demo_scaling,
        )

        model = sp._build_model()
        sp._apply_scaling(model, None, block_based="foo")

        # Should not apply any constraint scaling - no scaling suffix should be present
        assert not hasattr(model, "scaling_factor")

    @pytest.mark.unit
    def test_apply_scaling_constraint_based_method(self):
        # Dummy method to apply scaling by constraint
        def dummy_scaler(constraint):
            set_scaling_factor(constraint, 2)

        sp = ScalingProfiler(
            build_model=demo_model,
            user_scaling=demo_scaling,
        )

        # Constraint based method, so use block_based=False
        model = sp._build_model()
        sp._apply_scaling(model, dummy_scaler, block_based=False)

        # Should have scaling factors of 2 for all constraints
        assert len(model.scaling_factor) == 2
        assert model.scaling_factor[model.c1] == 2
        assert model.scaling_factor[model.c2] == 2

    @pytest.mark.unit
    def test_apply_scaling_block_based_method(self):
        # Dummy method to apply scaling by block
        def dummy_scaler(block):
            set_scaling_factor(block.c1, 2)
            set_scaling_factor(block.c2, 4)

        sp = ScalingProfiler(
            build_model=demo_model,
            user_scaling=demo_scaling,
        )

        # Block based method, so use block_based=True
        model = sp._build_model()
        sp._apply_scaling(model, dummy_scaler, block_based=True)

        # Should have scaling factors of for all constraints
        assert len(model.scaling_factor) == 2
        assert model.scaling_factor[model.c1] == 2
        assert model.scaling_factor[model.c2] == 4

    @pytest.mark.unit
    def test_solve_perturbed_state_no_callback(self):
        sp = ScalingProfiler(
            build_model=demo_model,
            user_scaling=demo_scaling,
        )

        model = sp._build_model()
        res = sp._solved_perturbed_state(model)

        assert res == {}

    @pytest.mark.unit
    def test_solve_perturbed_state(self):
        sp = ScalingProfiler(
            build_model=demo_model,
            user_scaling=demo_scaling,
            perturb_state=demo_pertubration,
        )

        model = sp._build_model()
        res = sp._solved_perturbed_state(model)

        assert value(model.v3) == pytest.approx(-12, rel=1e-6)

        # Model is trivially presolvable, so expect convergence in 0 iterations
        assert res["solved"]
        assert (
            res["termination_message"]
            == "TerminationCondition.convergenceCriteriaSatisfied"
        )
        assert res["iterations"] == 0


# Case study using Gibbs reactor model
# Get solution json from scaling tests
FILENAME = "gibbs_solution.json"
local_path = os.path.dirname(os.path.realpath(__file__))
fname = os.path.join(local_path, FILENAME)


def build_model():
    model = ConcreteModel()
    model.fs = FlowsheetBlock(dynamic=False)

    model.fs.properties = MethaneCombustionParameterBlock()

    model.fs.unit = GibbsReactor(
        property_package=model.fs.properties,
        has_heat_transfer=True,
        has_pressure_change=True,
    )

    model.fs.unit.inlet.flow_mol[0].fix(230.0)
    model.fs.unit.inlet.mole_frac_comp[0, "H2"].fix(0.0435)
    model.fs.unit.inlet.mole_frac_comp[0, "N2"].fix(0.6522)
    model.fs.unit.inlet.mole_frac_comp[0, "O2"].fix(0.1739)
    model.fs.unit.inlet.mole_frac_comp[0, "CO2"].fix(1e-5)
    model.fs.unit.inlet.mole_frac_comp[0, "CH4"].fix(0.1304)
    model.fs.unit.inlet.mole_frac_comp[0, "CO"].fix(1e-5)
    model.fs.unit.inlet.mole_frac_comp[0, "H2O"].fix(1e-5)
    model.fs.unit.inlet.mole_frac_comp[0, "NH3"].fix(1e-5)
    model.fs.unit.inlet.temperature[0].fix(1500.0)
    model.fs.unit.inlet.pressure[0].fix(101325.0)

    model.fs.unit.outlet.temperature[0].fix(2844.38)
    model.fs.unit.deltaP.fix(0)

    from_json(model, fname=fname, wts=StoreSpec.value())

    return model


def scale_vars(model):
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

    return model


def perturb_solution(model):
    # Decrease O2 in feed
    model.fs.unit.inlet.mole_frac_comp[0, "N2"].fix(0.7222)
    model.fs.unit.inlet.mole_frac_comp[0, "O2"].fix(0.1039)


expected_profile = {
    "Unscaled": {
        "Manual": {
            "condition_number": 5.70342e17,
            "solved": False,
            "termination_message": "TerminationCondition.locallyInfeasible",
            "iterations": 57,
            "iters_in_restoration": 27,
            "iters_w_regularization": 21,
        }
    },
    "Vars Only": {
        "Manual": {
            "condition_number": 9.24503e16,
            "solved": False,
            "termination_message": "TerminationCondition.locallyInfeasible",
            "iterations": 82,
            "iters_in_restoration": 82,
            "iters_w_regularization": 39,
        },
        "Auto": {
            "condition_number": 6.57667e14,
            "solved": True,
            "termination_message": "TerminationCondition.convergenceCriteriaSatisfied",
            "iterations": 9,
            "iters_in_restoration": 0,
            "iters_w_regularization": 0,
        },
    },
    "Harmonic": {
        "Manual": {
            "condition_number": 3.73643e18,
            "solved": False,
            "termination_message": "TerminationCondition.locallyInfeasible",
            "iterations": 39,
            "iters_in_restoration": 39,
            "iters_w_regularization": 0,
        },
        "Auto": {
            "condition_number": 2.83944e12,
            "solved": True,
            "termination_message": "TerminationCondition.convergenceCriteriaSatisfied",
            "iterations": 17,
            "iters_in_restoration": 0,
            "iters_w_regularization": 0,
        },
    },
    "Inverse Sum": {
        "Manual": {
            "condition_number": 9.31670e15,
            "solved": False,
            "termination_message": "TerminationCondition.iterationLimit",
            "iterations": 200,
            "iters_in_restoration": 200,
            "iters_w_regularization": 75,
        },
        "Auto": {
            "condition_number": 1.50163e6,
            "solved": True,
            "termination_message": "TerminationCondition.convergenceCriteriaSatisfied",
            "iterations": 6,
            "iters_in_restoration": 0,
            "iters_w_regularization": 0,
        },
    },
    "Inverse Root Sum Squares": {
        "Manual": {
            "condition_number": 1.15511e16,
            "solved": False,
            "termination_message": "TerminationCondition.iterationLimit",
            "iterations": 200,
            "iters_in_restoration": 201,
            "iters_w_regularization": 107,
        },
        "Auto": {
            "condition_number": 9.59994e5,
            "solved": True,
            "termination_message": "TerminationCondition.convergenceCriteriaSatisfied",
            "iterations": 6,
            "iters_in_restoration": 0,
            "iters_w_regularization": 0,
        },
    },
    "Inverse Maximum": {
        "Manual": {
            "condition_number": 1.094304e16,
            "solved": False,
            "termination_message": "TerminationCondition.iterationLimit",
            "iterations": 200,
            "iters_in_restoration": 197,
            "iters_w_regularization": 75,
        },
        "Auto": {
            "condition_number": 7.84576e5,
            "solved": True,
            "termination_message": "TerminationCondition.convergenceCriteriaSatisfied",
            "iterations": 6,
            "iters_in_restoration": 0,
            "iters_w_regularization": 0,
        },
    },
    "Inverse Minimum": {
        "Manual": {
            "condition_number": 7.34636e18,
            "solved": False,
            "termination_message": "TerminationCondition.locallyInfeasible",
            "iterations": 49,
            "iters_in_restoration": 49,
            "iters_w_regularization": 1,
        },
        "Auto": {
            "condition_number": 5.600998e12,
            "solved": True,
            "termination_message": "TerminationCondition.convergenceCriteriaSatisfied",
            "iterations": 16,
            "iters_in_restoration": 0,
            "iters_w_regularization": 0,
        },
    },
    "Nominal L1 Norm": {
        "Manual": {
            "condition_number": 1.18925e16,
            "solved": False,
            "termination_message": "TerminationCondition.locallyInfeasible",
            "iterations": 61,
            "iters_in_restoration": 60,
            "iters_w_regularization": 15,
        },
        "Auto": {
            "condition_number": 2.06015e6,
            "solved": True,
            "termination_message": "TerminationCondition.convergenceCriteriaSatisfied",
            "iterations": 4,
            "iters_in_restoration": 0,
            "iters_w_regularization": 0,
        },
    },
    "Nominal L2 Norm": {
        "Manual": {
            "condition_number": 1.18824e16,
            "solved": False,
            "termination_message": "TerminationCondition.locallyInfeasible",
            "iterations": 53,
            "iters_in_restoration": 50,
            "iters_w_regularization": 7,
        },
        "Auto": {
            "condition_number": 3.07419e6,
            "solved": True,
            "termination_message": "TerminationCondition.convergenceCriteriaSatisfied",
            "iterations": 4,
            "iters_in_restoration": 0,
            "iters_w_regularization": 0,
        },
    },
    "Actual L1 Norm": {
        "Manual": {
            "condition_number": 1.46059e9,
            "solved": False,
            "termination_message": "TerminationCondition.locallyInfeasible",
            "iterations": 29,
            "iters_in_restoration": 29,
            "iters_w_regularization": 0,
        },
        "Auto": {
            "condition_number": 2986.99,
            "solved": True,
            "termination_message": "TerminationCondition.convergenceCriteriaSatisfied",
            "iterations": 6,
            "iters_in_restoration": 0,
            "iters_w_regularization": 0,
        },
    },
    "Actual L2 Norm": {
        "Manual": {
            "condition_number": 6.61297e8,
            "solved": False,
            "termination_message": "TerminationCondition.locallyInfeasible",
            "iterations": 29,
            "iters_in_restoration": 29,
            "iters_w_regularization": 0,
        },
        "Auto": {
            "condition_number": 2510.95,
            "solved": True,
            "termination_message": "TerminationCondition.convergenceCriteriaSatisfied",
            "iterations": 6,
            "iters_in_restoration": 0,
            "iters_w_regularization": 0,
        },
    },
}


@pytest.mark.unit
def test_write_profile_report():
    sp = ScalingProfiler(
        build_model=build_model,
        user_scaling=scale_vars,
        perturb_state=perturb_solution,
    )

    stream = StringIO()

    sp.write_profile_report(results=expected_profile, stream=stream)

    expected = """
============================================================================
Scaling Profile Report
----------------------------------------------------------------------------
Scaling Method           || User Scaling           || Perfect Scaling
Unscaled                 || 5.703E+17 | Failed 57  ||
Vars Only                || 9.245E+16 | Failed 82  || 6.577E+14 | Solved 9  
Harmonic                 || 3.736E+18 | Failed 39  || 2.839E+12 | Solved 17 
Inverse Sum              || 9.317E+15 | Failed 200 || 1.502E+06 | Solved 6  
Inverse Root Sum Squares || 1.155E+16 | Failed 200 || 9.600E+05 | Solved 6  
Inverse Maximum          || 1.094E+16 | Failed 200 || 7.846E+05 | Solved 6  
Inverse Minimum          || 7.346E+18 | Failed 49  || 5.601E+12 | Solved 16 
Nominal L1 Norm          || 1.189E+16 | Failed 61  || 2.060E+06 | Solved 4  
Nominal L2 Norm          || 1.188E+16 | Failed 53  || 3.074E+06 | Solved 4  
Actual L1 Norm           || 1.461E+09 | Failed 29  || 2.987E+03 | Solved 6  
Actual L2 Norm           || 6.613E+08 | Failed 29  || 2.511E+03 | Solved 6  
============================================================================
"""

    assert stream.getvalue() == expected


@pytest.mark.integration
def test_case_study_profiling():
    sp = ScalingProfiler(
        build_model=build_model,
        user_scaling=scale_vars,
        perturb_state=perturb_solution,
    )

    results = sp.profile_scaling_methods()

    for cmeth, stats in results.items():
        for vmeth in ["Manual", "Auto"]:
            if cmeth == "Unscaled" and vmeth == "Auto":
                # Unscaled does not have data for auto
                continue
            rstats = stats[vmeth]
            xstats = expected_profile[cmeth][vmeth]
            assert rstats["condition_number"] == pytest.approx(
                xstats["condition_number"], rel=1e-5
            )
            assert rstats["solved"] == xstats["solved"]
            assert rstats["termination_message"] == xstats["termination_message"]
            for iters in [
                "iterations",
                "iters_in_restoration",
                "iters_w_regularization",
            ]:
                # We will allow a variance of 2 iteration in this test to avoid being overly fragile
                assert rstats[iters] == pytest.approx(xstats[iters], abs=2)


@pytest.mark.integration
def test_report_scaling_profiles():
    sp = ScalingProfiler(
        build_model=build_model,
        user_scaling=scale_vars,
        perturb_state=perturb_solution,
    )

    stream = StringIO()

    sp.report_scaling_profiles(stream=stream)

    expected = """
============================================================================
Scaling Profile Report
----------------------------------------------------------------------------
Scaling Method           || User Scaling           || Perfect Scaling
Unscaled                 || 5.703E+17 | Failed 57  ||
Vars Only                || 9.245E+16 | Failed 82  || 6.577E+14 | Solved 9  
Harmonic                 || 3.736E+18 | Failed 39  || 2.839E+12 | Solved 17 
Inverse Sum              || 9.317E+15 | Failed 200 || 1.502E+06 | Solved 6  
Inverse Root Sum Squares || 1.155E+16 | Failed 200 || 9.600E+05 | Solved 6  
Inverse Maximum          || 1.094E+16 | Failed 200 || 7.846E+05 | Solved 6  
Inverse Minimum          || 7.346E+18 | Failed 49  || 5.601E+12 | Solved 16 
Nominal L1 Norm          || 1.189E+16 | Failed 61  || 2.060E+06 | Solved 4  
Nominal L2 Norm          || 1.188E+16 | Failed 53  || 3.074E+06 | Solved 4  
Actual L1 Norm           || 1.461E+09 | Failed 29  || 2.987E+03 | Solved 6  
Actual L2 Norm           || 6.613E+08 | Failed 29  || 2.511E+03 | Solved 6  
============================================================================
"""

    assert stream.getvalue() == expected
