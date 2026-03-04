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
This module contains tests for the convergence analysis tools.
"""

from copy import deepcopy
from io import StringIO
import os
from copy import deepcopy
from pyomo.contrib.pynumero.interfaces.external_grey_box import (
    ExternalGreyBoxBlock,
    ExternalGreyBoxModel,
)

from pandas import DataFrame
import pytest

from pyomo.environ import (
    ConcreteModel,
    Constraint,
    SolverFactory,
    value,
    Var,
    assert_optimal_termination,
)
from pyomo.common.fileutils import this_file_dir
from pyomo.common.tempfiles import TempfileManager

from idaes.core.util.diagnostics_tools.convergence_analysis import (
    IpoptConvergenceAnalysis,
)
from idaes.core.util.parameter_sweep import (
    SequentialSweepRunner,
    ParameterSweepSpecification,
)
from idaes.core.surrogate.pysmo.sampling import (
    UniformSampling,
)

__author__ = "Alex Dowling, Douglas Allan, Andrew Lee"


currdir = this_file_dir()
PSWEEP_JSON_PATH = os.path.join(currdir, "..", "..", "tests", "load_psweep.json")


ca_dict = {
    "specification": {
        "inputs": {
            "v2": {
                "pyomo_path": "v2",
                "lower": 2,
                "upper": 6,
            },
        },
        "sampling_method": "UniformSampling",
        "sample_size": [2],
        "samples": {
            "index": [0, 1],
            "columns": ["v2"],
            "data": [[2.0], [6.0]],
            "index_names": [None],
            "column_names": [None],
        },
    },
    "results": {
        0: {
            "success": True,
            "results": 2,
        },
        1: {
            "success": True,
            "results": 6,
        },
    },
}

ca_res = {
    "specification": {
        "inputs": {"v2": {"pyomo_path": "v2", "lower": 2, "upper": 6}},
        "sampling_method": "UniformSampling",
        "sample_size": [2],
        "samples": {
            "index": [0, 1],
            "columns": ["v2"],
            "data": [[2.0], [6.0]],
            "index_names": [None],
            "column_names": [None],
        },
    },
    "results": {
        0: {
            "success": False,
            "results": {
                "iters": 7,
                "iters_in_restoration": 4,
                "iters_w_regularization": 0,
                "time": 0.0,
                "numerical_issues": True,
            },
        },
        1: {
            "success": False,
            "results": {
                "iters": 7,
                "iters_in_restoration": 4,
                "iters_w_regularization": 0,
                "time": 0.0,
                "numerical_issues": True,
            },
        },
    },
}


class TestIpoptConvergenceAnalysis:
    @pytest.fixture
    def model(self):
        m = ConcreteModel()

        m.v1 = Var(bounds=(None, 1))
        m.v2 = Var()
        m.c = Constraint(expr=m.v1 == m.v2)

        m.v2.fix(0)

        return m

    @pytest.mark.component
    def test_with_grey_box(self):

        class BasicGrayBox(ExternalGreyBoxModel):
            def input_names(self):
                return ["a1", "a2", "a3"]

            def output_names(self):
                return ["o1", "o2"]

            def equality_constraint_names(self):
                return ["a_sum"]

            def evaluate_equality_constraints(self):
                a1 = self._input_values[0]
                a2 = self._input_values[1]
                return [a1 * 0.5 + a2]

        m = ConcreteModel()

        m.gb = ExternalGreyBoxBlock(external_model=BasicGrayBox())
        with pytest.raises(NotImplementedError):
            IpoptConvergenceAnalysis(model=m)

    @pytest.mark.unit
    def test_init(self, model):
        ca = IpoptConvergenceAnalysis(model)

        assert ca._model is model
        assert isinstance(ca._psweep, SequentialSweepRunner)
        assert isinstance(ca.results, dict)
        assert ca.config.input_specification is None
        assert ca.config.solver_options is None

    @pytest.mark.unit
    def test_init_solver_options(self, model):
        solver_obj = SolverFactory("ipopt")
        options = solver_obj.options["max_iter"] = 100
        ca = IpoptConvergenceAnalysis(
            model, solver_obj=solver_obj, solver_options=options
        )

        assert ca._model is model
        assert isinstance(ca._psweep, SequentialSweepRunner)
        assert isinstance(ca.results, dict)
        assert ca.config.input_specification is None
        assert ca.config.solver_options is not None

    @pytest.mark.unit
    def test_build_model(self, model):
        ca = IpoptConvergenceAnalysis(model)

        clone = ca._build_model()

        assert clone is not model
        assert isinstance(clone.v1, Var)
        assert isinstance(clone.v2, Var)
        assert isinstance(clone.c, Constraint)

    @pytest.mark.unit
    def test_parse_ipopt_output(self, model):
        ca = IpoptConvergenceAnalysis(model)

        fname = os.path.join(currdir, "ipopt_output.txt")
        iters, restoration, regularization, time = ca._parse_ipopt_output(fname)

        assert iters == 43
        assert restoration == 39
        assert regularization == 4
        assert time == 0.016 + 0.035

    @pytest.mark.component
    @pytest.mark.solver
    def test_run_ipopt_with_stats(self):
        m = ConcreteModel()
        m.v1 = Var(initialize=1)
        m.c1 = Constraint(expr=m.v1 == 4)

        ca = IpoptConvergenceAnalysis(m)
        solver = SolverFactory("ipopt")

        (
            status,
            iters,
            iters_in_restoration,
            iters_w_regularization,
            time,
        ) = ca._run_ipopt_with_stats(m, solver)

        assert_optimal_termination(status)
        assert iters == 1
        assert iters_in_restoration == 0
        assert iters_w_regularization == 0
        assert isinstance(time, float)

    @pytest.mark.component
    @pytest.mark.solver
    def test_run_model(self, model):
        ca = IpoptConvergenceAnalysis(model)

        model.v2.fix(0.5)

        solver = SolverFactory("ipopt")

        success, run_stats = ca._run_model(model, solver)

        assert success
        assert value(model.v1) == pytest.approx(0.5, rel=1e-8)

        assert len(run_stats) == 4
        assert run_stats[0] == 1
        assert run_stats[1] == 0
        assert run_stats[2] == 0

    @pytest.mark.unit
    def test_build_outputs(self, model):
        ca = IpoptConvergenceAnalysis(model)

        model.v1.set_value(0.5)
        model.v2.fix(0.5)

        results = ca._build_outputs(model, (1, 2, 3, 4))

        assert results == {
            "iters": 1,
            "iters_in_restoration": 2,
            "iters_w_regularization": 3,
            "time": 4,
            "numerical_issues": False,
        }

    @pytest.mark.unit
    def test_build_outputs_with_warnings(self, model):
        ca = IpoptConvergenceAnalysis(model)

        model.v1.set_value(4)
        model.v2.fix(0.5)

        results = ca._build_outputs(model, (1, 2, 3, 4))

        assert results == {
            "iters": 1,
            "iters_in_restoration": 2,
            "iters_w_regularization": 3,
            "time": 4,
            "numerical_issues": True,
        }

    @pytest.mark.unit
    def test_recourse(self, model):
        ca = IpoptConvergenceAnalysis(model)

        assert ca._recourse(model) == {
            "iters": -1,
            "iters_in_restoration": -1,
            "iters_w_regularization": -1,
            "time": -1,
            "numerical_issues": -1,
        }

    @pytest.mark.integration
    @pytest.mark.solver
    def test_run_convergence_analysis(self, model):
        spec = ParameterSweepSpecification()
        spec.add_sampled_input("v2", lower=0, upper=3)
        spec.set_sampling_method(UniformSampling)
        spec.set_sample_size([4])

        ca = IpoptConvergenceAnalysis(model, input_specification=spec)

        ca.run_convergence_analysis()

        assert isinstance(ca.results, dict)
        assert len(ca.results) == 4

        # Ignore time, as it is too noisy to test
        # Sample 0 should solve cleanly
        assert ca.results[0]["success"]
        assert ca.results[0]["results"]["iters"] == 0
        assert ca.results[0]["results"]["iters_in_restoration"] == 0
        assert ca.results[0]["results"]["iters_w_regularization"] == 0
        assert not ca.results[0]["results"]["numerical_issues"]

        # Sample 1 should solve, but have issues due to bound on v1
        assert ca.results[1]["success"]
        assert ca.results[1]["results"]["iters"] == pytest.approx(3, abs=1)
        assert ca.results[1]["results"]["iters_in_restoration"] == 0
        assert ca.results[1]["results"]["iters_w_regularization"] == 0
        assert ca.results[1]["results"]["numerical_issues"]

        # Other iterations should fail due to bound
        assert not ca.results[2]["success"]
        assert ca.results[2]["results"]["iters"] == pytest.approx(7, abs=1)
        assert ca.results[2]["results"]["iters_in_restoration"] == pytest.approx(
            4, abs=1
        )
        assert ca.results[2]["results"]["iters_w_regularization"] == 0
        assert ca.results[2]["results"]["numerical_issues"]

        assert not ca.results[3]["success"]
        assert ca.results[3]["results"]["iters"] == pytest.approx(8, abs=1)
        assert ca.results[3]["results"]["iters_in_restoration"] == pytest.approx(
            5, abs=1
        )
        assert ca.results[3]["results"]["iters_w_regularization"] == 0
        assert ca.results[3]["results"]["numerical_issues"]

    @pytest.fixture(scope="class")
    def ca_with_results(self):
        spec = ParameterSweepSpecification()
        spec.set_sampling_method(UniformSampling)
        spec.add_sampled_input("v2", 2, 6)
        spec.set_sample_size([2])
        spec.generate_samples()

        ca = IpoptConvergenceAnalysis(
            model=ConcreteModel(),
            input_specification=spec,
        )

        ca._psweep._results = {
            0: {"success": True, "results": 2},
            1: {"success": True, "results": 6},
        }

        return ca

    @pytest.mark.unit
    def test_report_convergence_summary(self):
        stream = StringIO()

        ca = IpoptConvergenceAnalysis(
            model=ConcreteModel(),
        )

        ca._psweep._results = {
            0: {
                "success": True,
                "results": {
                    "iters_in_restoration": 1,
                    "iters_w_regularization": 0,
                    "numerical_issues": 10,
                },
            },
            1: {
                "success": True,
                "results": {
                    "iters_in_restoration": 0,
                    "iters_w_regularization": 5,
                    "numerical_issues": 5,
                },
            },
            2: {
                "success": False,
                "results": {
                    "iters_in_restoration": 0,
                    "iters_w_regularization": 0,
                    "numerical_issues": 0,
                },
            },
        }

        ca.report_convergence_summary(stream)

        expected = """Successes: 2, Failures 1 (66.66666666666667%)
Runs with Restoration: 1
Runs with Regularization: 1
Runs with Numerical Issues: 2
"""

        assert stream.getvalue() == expected

    @pytest.mark.component
    def test_to_dict(self, ca_with_results):
        outdict = ca_with_results.to_dict()
        assert outdict == ca_dict

    @pytest.mark.unit
    def test_from_dict(self):
        ca = IpoptConvergenceAnalysis(
            model=ConcreteModel(),
        )

        ca.from_dict(ca_dict)

        input_spec = ca._psweep.get_input_specification()

        assert isinstance(input_spec, ParameterSweepSpecification)
        assert len(input_spec.inputs) == 1

        assert input_spec.sampling_method is UniformSampling
        assert isinstance(input_spec.samples, DataFrame)
        assert input_spec.sample_size == [2]

        assert isinstance(ca.results, dict)
        assert len(ca.results) == 2

        for i in [0, 1]:
            assert ca.results[i]["success"]
            assert ca.results[i]["results"] == 2 + i * 4

    @pytest.mark.component
    def test_to_json_file(self, ca_with_results):
        temp_context = TempfileManager.new_context()
        tmpfile = temp_context.create_tempfile(suffix=".json")

        ca_with_results.to_json_file(tmpfile)

        with open(tmpfile, "r") as f:
            lines = f.read()
        f.close()

        expected = """{
   "specification": {
      "inputs": {
         "v2": {
            "pyomo_path": "v2",
            "lower": 2,
            "upper": 6
         }
      },
      "sampling_method": "UniformSampling",
      "sample_size": [
         2
      ],
      "samples": {
         "index": [
            0,
            1
         ],
         "columns": [
            "v2"
         ],
         "data": [
            [
               2.0
            ],
            [
               6.0
            ]
         ],
         "index_names": [
            null
         ],
         "column_names": [
            null
         ]
      }
   },
   "results": {
      "0": {
         "success": true,
         "results": 2
      },
      "1": {
         "success": true,
         "results": 6
      }
   }
}"""

        assert lines == expected

        # Check for clean up
        temp_context.release(remove=True)
        assert not os.path.exists(tmpfile)

    @pytest.mark.unit
    def test_load_from_json_file(self):
        fname = PSWEEP_JSON_PATH

        ca = IpoptConvergenceAnalysis(
            model=ConcreteModel(),
        )
        ca.from_json_file(fname)

        input_spec = ca._psweep.get_input_specification()

        assert isinstance(input_spec, ParameterSweepSpecification)
        assert len(input_spec.inputs) == 1

        assert input_spec.sampling_method is UniformSampling
        assert isinstance(input_spec.samples, DataFrame)
        assert input_spec.sample_size == [2]

        assert isinstance(ca.results, dict)
        assert len(ca.results) == 2

        for i in [0, 1]:
            assert ca.results[i]["success"]
            assert ca.results[i]["results"] == 2 + i * 4

    @pytest.mark.integration
    @pytest.mark.solver
    def test_run_convergence_analysis_from_dict(self, model):
        ca = IpoptConvergenceAnalysis(
            model=model,
        )
        ca.run_convergence_analysis_from_dict(ca_dict)

        input_spec = ca._psweep.get_input_specification()

        assert isinstance(input_spec, ParameterSweepSpecification)
        assert len(input_spec.inputs) == 1

        assert input_spec.sampling_method is UniformSampling
        assert isinstance(input_spec.samples, DataFrame)
        assert input_spec.sample_size == [2]

        assert isinstance(ca.results, dict)
        assert len(ca.results) == 2

        assert not ca.results[0]["success"]
        assert ca.results[0]["results"]["iters"] == pytest.approx(7, abs=1)
        assert ca.results[0]["results"]["iters_in_restoration"] == pytest.approx(
            4, abs=1
        )
        assert ca.results[0]["results"]["iters_w_regularization"] == 0
        assert ca.results[0]["results"]["numerical_issues"]

        assert not ca.results[1]["success"]
        assert ca.results[1]["results"]["iters"] == pytest.approx(7, abs=1)
        assert ca.results[1]["results"]["iters_in_restoration"] == pytest.approx(
            4, abs=1
        )
        assert ca.results[1]["results"]["iters_w_regularization"] == 0
        assert ca.results[1]["results"]["numerical_issues"]

    @pytest.mark.integration
    @pytest.mark.solver
    def test_run_convergence_analysis_from_file(self, model):
        fname = PSWEEP_JSON_PATH

        ca = IpoptConvergenceAnalysis(
            model=model,
        )
        ca.run_convergence_analysis_from_file(fname)

        input_spec = ca._psweep.get_input_specification()

        assert isinstance(input_spec, ParameterSweepSpecification)
        assert len(input_spec.inputs) == 1

        assert input_spec.sampling_method is UniformSampling
        assert isinstance(input_spec.samples, DataFrame)
        assert input_spec.sample_size == [2]

        assert isinstance(ca.results, dict)
        assert len(ca.results) == 2

        assert not ca.results[0]["success"]
        assert ca.results[0]["results"]["iters"] == pytest.approx(7, abs=1)
        assert ca.results[0]["results"]["iters_in_restoration"] == pytest.approx(
            4, abs=1
        )
        assert ca.results[0]["results"]["iters_w_regularization"] == 0
        assert ca.results[0]["results"]["numerical_issues"]

        assert not ca.results[1]["success"]
        assert ca.results[1]["results"]["iters"] == pytest.approx(7, abs=1)
        assert ca.results[1]["results"]["iters_in_restoration"] == pytest.approx(
            4, abs=1
        )
        assert ca.results[1]["results"]["iters_w_regularization"] == 0
        assert ca.results[1]["results"]["numerical_issues"]

    @pytest.fixture(scope="class")
    def conv_anal(self):
        ca = IpoptConvergenceAnalysis(
            model=ConcreteModel(),
        )
        ca.from_dict(ca_res)

        return ca

    @pytest.mark.unit
    def test_compare_results_to_dict_ok(self, conv_anal):
        diffs = conv_anal._compare_results_to_dict(ca_res)

        assert diffs["success"] == []
        assert diffs["iters"] == []
        assert diffs["iters_in_restoration"] == []
        assert diffs["iters_w_regularization"] == []
        assert diffs["numerical_issues"] == []

    @pytest.mark.unit
    def test_compare_results_to_dict_success(self, conv_anal):
        ca_copy = deepcopy(ca_res)
        ca_copy["results"][0]["success"] = True

        diffs = conv_anal._compare_results_to_dict(ca_copy)

        assert diffs["success"] == [0]
        assert diffs["iters"] == []
        assert diffs["iters_in_restoration"] == []
        assert diffs["iters_w_regularization"] == []
        assert diffs["numerical_issues"] == []

    @pytest.mark.unit
    def test_compare_results_to_dict_iters(self, conv_anal):
        ca_copy = deepcopy(ca_res)
        ca_copy["results"][0]["results"]["iters"] = 8
        ca_copy["results"][1]["results"]["iters"] = 9

        diffs = conv_anal._compare_results_to_dict(ca_copy)

        assert diffs["success"] == []
        assert diffs["iters"] == [1]
        assert diffs["iters_in_restoration"] == []
        assert diffs["iters_w_regularization"] == []
        assert diffs["numerical_issues"] == []

        diffs = conv_anal._compare_results_to_dict(ca_copy, abs_tol=0, rel_tol=0)

        assert diffs["success"] == []
        assert diffs["iters"] == [0, 1]
        assert diffs["iters_in_restoration"] == []
        assert diffs["iters_w_regularization"] == []
        assert diffs["numerical_issues"] == []

    @pytest.mark.unit
    def test_compare_results_to_dict_restoration(self, conv_anal):
        ca_copy = deepcopy(ca_res)
        ca_copy["results"][0]["results"]["iters_in_restoration"] = 5
        ca_copy["results"][1]["results"]["iters_in_restoration"] = 6

        diffs = conv_anal._compare_results_to_dict(ca_copy)

        assert diffs["success"] == []
        assert diffs["iters"] == []
        assert diffs["iters_in_restoration"] == [1]
        assert diffs["iters_w_regularization"] == []
        assert diffs["numerical_issues"] == []

        diffs = conv_anal._compare_results_to_dict(ca_copy, abs_tol=0, rel_tol=0)

        assert diffs["success"] == []
        assert diffs["iters"] == []
        assert diffs["iters_in_restoration"] == [0, 1]
        assert diffs["iters_w_regularization"] == []
        assert diffs["numerical_issues"] == []

    @pytest.mark.unit
    def test_compare_results_to_dict_regularization(self, conv_anal):
        ca_copy = deepcopy(ca_res)
        ca_copy["results"][0]["results"]["iters_w_regularization"] = 1
        ca_copy["results"][1]["results"]["iters_w_regularization"] = 2

        diffs = conv_anal._compare_results_to_dict(ca_copy)

        assert diffs["success"] == []
        assert diffs["iters"] == []
        assert diffs["iters_in_restoration"] == []
        assert diffs["iters_w_regularization"] == [1]
        assert diffs["numerical_issues"] == []

        diffs = conv_anal._compare_results_to_dict(ca_copy, abs_tol=0, rel_tol=0)

        assert diffs["success"] == []
        assert diffs["iters"] == []
        assert diffs["iters_in_restoration"] == []
        assert diffs["iters_w_regularization"] == [0, 1]
        assert diffs["numerical_issues"] == []

    @pytest.mark.unit
    def test_compare_results_to_dict_numerical_issues(self, conv_anal):
        ca_copy = deepcopy(ca_res)
        ca_copy["results"][1]["results"]["numerical_issues"] = False

        diffs = conv_anal._compare_results_to_dict(ca_copy)

        assert diffs["success"] == []
        assert diffs["iters"] == []
        assert diffs["iters_in_restoration"] == []
        assert diffs["iters_w_regularization"] == []
        assert diffs["numerical_issues"] == [1]

    @pytest.mark.integration
    @pytest.mark.solver
    def test_compare_convergence_to_baseline(self, model):
        fname = os.path.join(currdir, "convergence_baseline.json")

        ca = IpoptConvergenceAnalysis(
            model=model,
        )

        diffs = ca.compare_convergence_to_baseline(fname)

        # Baseline has incorrect values
        assert diffs == {
            "success": [0],
            "iters": [],
            "iters_in_restoration": [],
            "iters_w_regularization": [1],
            "numerical_issues": [],
        }

    @pytest.mark.integration
    @pytest.mark.solver
    def test_assert_baseline_comparison(self, model):
        fname = os.path.join(currdir, "convergence_baseline.json")

        ca = IpoptConvergenceAnalysis(
            model=model,
        )

        # Baseline has incorrect values
        with pytest.raises(
            AssertionError,
            match="Convergence analysis does not match baseline",
        ):
            ca.assert_baseline_comparison(fname)
