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


import pytest
import os
import os.path
from collections import OrderedDict

from pandas import DataFrame
from pandas.testing import assert_frame_equal

from pyomo.environ import ConcreteModel, Constraint, Var, Param, value, SolverFactory
from pyomo.common.fileutils import this_file_dir
from pyomo.solvers.plugins.solvers.IPOPT import IPOPT
from pyomo.common.tempfiles import TempfileManager

from idaes.core.util.parameter_sweep import ParameterSweepSpecification, ParameterSweep
from idaes.core.surrogate.pysmo.sampling import (
    LatinHypercubeSampling,
    UniformSampling,
    HaltonSampling,
    HammersleySampling,
    CVTSampling,
)
from idaes.core.util.exceptions import ConfigurationError

currdir = this_file_dir()


expected_todict = {
    "inputs": OrderedDict(
        [
            (
                "foo",
                OrderedDict(
                    [
                        ("pyomo_path", "model.foo"),
                        ("lower", 0),
                        ("upper", 10),
                    ]
                ),
            ),
            (
                "bar",
                OrderedDict(
                    [
                        ("pyomo_path", "model.bar"),
                        ("lower", 20),
                        ("upper", 40),
                    ]
                ),
            ),
        ]
    ),
    "sampling_method": "UniformSampling",
    "sample_size": [4, 3],
    "samples": {
        "index": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
        "columns": ["foo", "bar"],
        "data": [
            [0.0, 20.0],
            [0.0, 30.0],
            [0.0, 40.0],
            [3.333333333333333, 20.0],
            [3.333333333333333, 30.0],
            [3.333333333333333, 40.0],
            [6.666666666666666, 20.0],
            [6.666666666666666, 30.0],
            [6.666666666666666, 40.0],
            [10.0, 20.0],
            [10.0, 30.0],
            [10.0, 40.0],
        ],
        "index_names": [None],
        "column_names": [None],
    },
}


class TestParameterSweepSpecification:
    @pytest.mark.unit
    def test_init(self):
        spec = ParameterSweepSpecification()

        assert isinstance(spec.inputs, OrderedDict)
        assert len(spec.inputs) == 0

        assert spec.sampling_method is None
        assert spec.samples is None
        assert spec.sample_size is None

    @pytest.mark.unit
    def test_add_sampled_input(self):
        spec = ParameterSweepSpecification()

        spec.add_sampled_input("foo", "model.foo", 0, 10)
        spec.add_sampled_input("bar", "model.bar", 20, 40)

        assert len(spec.inputs) == 2

        assert isinstance(spec.inputs["foo"], OrderedDict)
        assert spec.inputs["foo"]["pyomo_path"] == "model.foo"
        assert spec.inputs["foo"]["lower"] == 0
        assert spec.inputs["foo"]["upper"] == 10

        assert isinstance(spec.inputs["bar"], OrderedDict)
        assert spec.inputs["bar"]["pyomo_path"] == "model.bar"
        assert spec.inputs["bar"]["lower"] == 20
        assert spec.inputs["bar"]["upper"] == 40

    @pytest.mark.unit
    def test_generate_pysmo_data_input(self):
        spec = ParameterSweepSpecification()

        spec.add_sampled_input("foo", "model.foo", 0, 10)
        spec.add_sampled_input("bar", "model.bar", 20, 40)

        data_spec = spec._generate_pysmo_data_input()

        assert isinstance(data_spec, list)
        assert len(data_spec) == 2

        assert data_spec[0] == [0, 20]
        assert data_spec[1] == [10, 40]

    @pytest.mark.unit
    def test_set_sampling_method(self):
        spec = ParameterSweepSpecification()

        spec.set_sampling_method(LatinHypercubeSampling)
        assert spec.sampling_method is LatinHypercubeSampling

        spec.set_sampling_method(UniformSampling)
        assert spec.sampling_method is UniformSampling

        spec.set_sampling_method(HaltonSampling)
        assert spec.sampling_method is HaltonSampling

        spec.set_sampling_method(HammersleySampling)
        assert spec.sampling_method is HammersleySampling

        spec.set_sampling_method(CVTSampling)
        assert spec.sampling_method is CVTSampling

    @pytest.mark.unit
    def test_set_sample_size(self):
        spec = ParameterSweepSpecification()

        spec.set_sample_size(10)

        assert spec.sample_size == 10

    @pytest.mark.unit
    def test_set_sampling_method_invalid(self):
        spec = ParameterSweepSpecification()

        with pytest.raises(
            TypeError,
            match="Sampling method must be an instance of a Pysmo SamplingMethod "
            "\(received foo\)",
        ):
            spec.set_sampling_method("foo")

    @pytest.mark.unit
    def test_generate_samples_no_inputs(self):
        spec = ParameterSweepSpecification()
        spec.set_sampling_method(LatinHypercubeSampling)

        with pytest.raises(
            ValueError, match="Please identify at least on input variable to sample."
        ):
            spec.generate_samples()

    @pytest.mark.unit
    def test_generate_samples_none(self):
        spec = ParameterSweepSpecification()
        spec.set_sampling_method(LatinHypercubeSampling)

        spec.add_sampled_input("foo", "model.foo", 0, 10)
        spec.add_sampled_input("bar", "model.bar", 20, 40)

        with pytest.raises(ValueError, match="Please set a sample size."):
            spec.generate_samples()

    @pytest.mark.unit
    def test_generate_samples_zero(self):
        spec = ParameterSweepSpecification()
        spec.set_sampling_method(LatinHypercubeSampling)

        spec.add_sampled_input("foo", "model.foo", 0, 10)
        spec.add_sampled_input("bar", "model.bar", 20, 40)

        spec.set_sample_size(0)

        with pytest.raises(
            ValueError, match="sample_size must be an integer greater than 1."
        ):
            spec.generate_samples()

    @pytest.mark.unit
    def test_generate_samples_not_int(self):
        spec = ParameterSweepSpecification()
        spec.set_sampling_method(LatinHypercubeSampling)

        spec.add_sampled_input("foo", "model.foo", 0, 10)
        spec.add_sampled_input("bar", "model.bar", 20, 40)

        spec.set_sample_size(0.5)

        with pytest.raises(ValueError, match="sample_size must be an integer."):
            spec.generate_samples()

    @pytest.mark.unit
    def test_generate_samples_uniform_not_list(self):
        spec = ParameterSweepSpecification()
        spec.set_sampling_method(UniformSampling)

        spec.add_sampled_input("foo", "model.foo", 0, 10)
        spec.add_sampled_input("bar", "model.bar", 20, 40)

        spec.set_sample_size(0.5)

        with pytest.raises(
            TypeError,
            match="For UniformSampling, sample_size must be list of integers.",
        ):
            spec.generate_samples()

    @pytest.mark.component
    def test_generate_samples_FF(self):
        spec = ParameterSweepSpecification()
        spec.set_sampling_method(UniformSampling)

        spec.add_sampled_input("foo", "model.foo", 0, 10)
        spec.add_sampled_input("bar", "model.bar", 20, 40)

        spec.set_sample_size([4, 3])

        samples = spec.generate_samples()
        assert isinstance(spec.samples, DataFrame)
        assert spec._samples.shape == (12, 2)

        expected = DataFrame(
            data=[
                [0.0, 20.0],
                [0.0, 30.0],
                [0.0, 40.0],
                [10 / 3, 20.0],
                [10 / 3, 30.0],
                [10 / 3, 40.0],
                [20 / 3, 20.0],
                [20 / 3, 30.0],
                [20 / 3, 40.0],
                [10.0, 20.0],
                [10.0, 30.0],
                [10.0, 40.0],
            ],
            columns=["foo", "bar"],
        )

        assert_frame_equal(samples, expected)

    @pytest.mark.integration
    def test_generate_samples_LHC(self):
        spec = ParameterSweepSpecification()
        spec.set_sampling_method(LatinHypercubeSampling)

        spec.add_sampled_input("foo", "model.foo", 0, 10)
        spec.add_sampled_input("bar", "model.bar", 20, 40)

        spec.set_sample_size(10)

        spec.generate_samples()
        assert spec._samples.shape == (10, 2)

    @pytest.mark.integration
    def test_generate_samples_CVT(self):
        spec = ParameterSweepSpecification()
        spec.set_sampling_method(CVTSampling)

        spec.add_sampled_input("foo", "model.foo", 0, 10)
        spec.add_sampled_input("bar", "model.bar", 20, 40)

        spec.set_sample_size(10)

        spec.generate_samples()
        assert spec._samples.shape == (10, 2)

    @pytest.mark.integration
    def test_generate_samples_Halton(self):
        spec = ParameterSweepSpecification()
        spec.set_sampling_method(HaltonSampling)

        spec.add_sampled_input("foo", "model.foo", 0, 10)
        spec.add_sampled_input("bar", "model.bar", 20, 40)

        spec.set_sample_size(10)

        spec.generate_samples()
        assert spec._samples.shape == (10, 2)

    @pytest.mark.integration
    def test_generate_samples_Hammersley(self):
        spec = ParameterSweepSpecification()
        spec.set_sampling_method(HammersleySampling)

        spec.add_sampled_input("foo", "model.foo", 0, 10)
        spec.add_sampled_input("bar", "model.bar", 20, 40)

        spec.set_sample_size(10)

        spec.generate_samples()
        assert spec._samples.shape == (10, 2)

    @pytest.mark.unit
    def test_to_dict(self):
        spec = ParameterSweepSpecification()
        spec.set_sampling_method(UniformSampling)
        spec.add_sampled_input("foo", "model.foo", 0, 10)
        spec.add_sampled_input("bar", "model.bar", 20, 40)
        spec.set_sample_size([4, 3])
        spec.generate_samples()

        outdict = spec.to_dict()

        assert outdict == expected_todict

    @pytest.mark.unit
    def test_from_dict(self):
        spec = ParameterSweepSpecification()

        spec.from_dict(expected_todict)

        assert isinstance(spec.inputs, OrderedDict)
        assert len(spec.inputs) == 2

        assert isinstance(spec.inputs["foo"], OrderedDict)
        assert spec.inputs["foo"]["pyomo_path"] == "model.foo"
        assert spec.inputs["foo"]["lower"] == 0
        assert spec.inputs["foo"]["upper"] == 10

        assert isinstance(spec.inputs["bar"], OrderedDict)
        assert spec.inputs["bar"]["pyomo_path"] == "model.bar"
        assert spec.inputs["bar"]["lower"] == 20
        assert spec.inputs["bar"]["upper"] == 40

        assert spec.sampling_method is UniformSampling
        assert spec.sample_size == [4, 3]

        assert isinstance(spec.samples, DataFrame)
        assert spec.samples.shape == (12, 2)

        expected = DataFrame(
            data=[
                [0.0, 20.0],
                [0.0, 30.0],
                [0.0, 40.0],
                [10 / 3, 20.0],
                [10 / 3, 30.0],
                [10 / 3, 40.0],
                [20 / 3, 20.0],
                [20 / 3, 30.0],
                [20 / 3, 40.0],
                [10.0, 20.0],
                [10.0, 30.0],
                [10.0, 40.0],
            ],
            columns=["foo", "bar"],
        )

        assert_frame_equal(spec.samples, expected)

    @pytest.mark.component
    def test_to_json_file(self):
        temp_context = TempfileManager.new_context()
        tmpfile = temp_context.create_tempfile(suffix=".json")

        spec = ParameterSweepSpecification()
        spec.set_sampling_method(UniformSampling)
        spec.add_sampled_input("foo", "model.foo", 0, 10)
        spec.add_sampled_input("bar", "model.bar", 20, 40)
        spec.set_sample_size([4, 3])
        spec.generate_samples()

        spec.to_json_file(tmpfile)

        with open(tmpfile, "r") as f:
            lines = f.read()
        f.close()

        expected = """{
   "inputs": {
      "foo": {
         "pyomo_path": "model.foo",
         "lower": 0,
         "upper": 10
      },
      "bar": {
         "pyomo_path": "model.bar",
         "lower": 20,
         "upper": 40
      }
   },
   "sampling_method": "UniformSampling",
   "sample_size": [
      4,
      3
   ],
   "samples": {
      "index": [
         0,
         1,
         2,
         3,
         4,
         5,
         6,
         7,
         8,
         9,
         10,
         11
      ],
      "columns": [
         "foo",
         "bar"
      ],
      "data": [
         [
            0.0,
            20.0
         ],
         [
            0.0,
            30.0
         ],
         [
            0.0,
            40.0
         ],
         [
            3.333333333333333,
            20.0
         ],
         [
            3.333333333333333,
            30.0
         ],
         [
            3.333333333333333,
            40.0
         ],
         [
            6.666666666666666,
            20.0
         ],
         [
            6.666666666666666,
            30.0
         ],
         [
            6.666666666666666,
            40.0
         ],
         [
            10.0,
            20.0
         ],
         [
            10.0,
            30.0
         ],
         [
            10.0,
            40.0
         ]
      ],
      "index_names": [
         null
      ],
      "column_names": [
         null
      ]
   }
}"""

        assert lines == expected

        # Check for clean up
        temp_context.release(remove=True)
        assert not os.path.exists(tmpfile)

    @pytest.mark.unit
    def test_from_json_file(self):
        fname = os.path.join(currdir, "load_spec.json")

        spec = ParameterSweepSpecification()

        spec.from_json_file(fname)

        assert isinstance(spec.inputs, OrderedDict)
        assert len(spec.inputs) == 2

        assert isinstance(spec.inputs["foo"], OrderedDict)
        assert spec.inputs["foo"]["pyomo_path"] == "model.foo"
        assert spec.inputs["foo"]["lower"] == 0
        assert spec.inputs["foo"]["upper"] == 10

        assert isinstance(spec.inputs["bar"], OrderedDict)
        assert spec.inputs["bar"]["pyomo_path"] == "model.bar"
        assert spec.inputs["bar"]["lower"] == 20
        assert spec.inputs["bar"]["upper"] == 40

        assert spec.sampling_method is UniformSampling
        assert spec.sample_size == [4, 3]

        assert isinstance(spec.samples, DataFrame)
        assert spec.samples.shape == (12, 2)

        expected = DataFrame(
            data=[
                [0.0, 20.0],
                [0.0, 30.0],
                [0.0, 40.0],
                [10 / 3, 20.0],
                [10 / 3, 30.0],
                [10 / 3, 40.0],
                [20 / 3, 20.0],
                [20 / 3, 30.0],
                [20 / 3, 40.0],
                [10.0, 20.0],
                [10.0, 30.0],
                [10.0, 40.0],
            ],
            columns=["foo", "bar"],
        )

        assert_frame_equal(spec.samples, expected)


# spec = ParameterSweepSpecification()
# spec.set_sampling_method(LatinHypercubeSampling)
# spec.add_sampled_input("v1", "v1", 0, 10)
# spec.add_sampled_input("v2", "v2", 20, 40)
# spec.add_sampled_input("p2", "p2", -10, 10)
# spec.set_sample_size(5)
# spec.generate_samples()
#
#
# ceval_dict = {
#     "specification": {
#         "inputs": OrderedDict(
#             [
#                 (
#                     "v2",
#                     OrderedDict(
#                         [
#                             ("pyomo_path", "v2"),
#                             ("lower", 2),
#                             ("upper", 6),
#                             ("mean", None),
#                             ("std", None),
#                             ("distribution", "normal"),
#                         ]
#                     ),
#                 )
#             ]
#         ),
#         "sampling_method": "UniformSampling",
#         "sample_size": [2],
#         "samples": {
#             "index": [0, 1],
#             "columns": ["v2"],
#             "data": [[2.0], [6.0]],
#             "index_names": [None],
#             "column_names": [None],
#         },
#     },
#     "results": OrderedDict(
#         [
#             (
#                 0,
#                 {
#                     "solved": True,
#                     "iters": 1,
#                     "iters_in_restoration": 0,
#                     "iters_w_regularization": 0,
#                     "time": 0.0,
#                     "numerical_issues": 0,
#                 },
#             ),
#             (
#                 1,
#                 {
#                     "solved": True,
#                     "iters": 1,
#                     "iters_in_restoration": 0,
#                     "iters_w_regularization": 0,
#                     "time": 0.0,
#                     "numerical_issues": 0,
#                 },
#             ),
#         ]
#     ),
# }
#
#
# class TestConvergenceEvaluation:
#     def build_model(self):
#         m = ConcreteModel()
#         m.v1 = Var(initialize=1)
#         m.v2 = Var(initialize=1)
#         m.v3 = Var(initialize=1)
#         m.p1 = Param(initialize=1)
#         m.p2 = Param(initialize=1, mutable=True)
#         m.p3 = Param(initialize=1, mutable=True)
#
#         m.v1.fix(1)
#         m.v2.fix(1)
#
#         return m
#
#     @pytest.mark.unit
#     def test_init(self):
#         ceval = cb.ConvergenceEvaluation()
#
#         assert ceval.config.build_model is None
#         assert ceval.config.input_specification is None
#         assert ceval.config.ipopt_options == {}
#         assert ceval._input_spec is None
#         assert ceval._results is None
#
#     @pytest.mark.unit
#     def test_get_initialized_model_none(self):
#         ceval = cb.ConvergenceEvaluation()
#
#         with pytest.raises(
#             ConfigurationError,
#             match="Please specify a method to construct the model of interest.",
#         ):
#             ceval.get_initialized_model()
#
#     @pytest.mark.unit
#     def test_get_initialized_model(self):
#         ceval = cb.ConvergenceEvaluation(
#             build_model=self.build_model,
#         )
#
#         m2 = ceval.get_initialized_model()
#
#         assert isinstance(m2, ConcreteModel)
#         assert isinstance(m2.v1, Var)
#         assert isinstance(m2.v2, Var)
#         assert isinstance(m2.p1, Param)
#         assert isinstance(m2.p2, Param)
#         assert isinstance(m2.p3, Param)
#
#     @pytest.mark.unit
#     def test_get_specification(self):
#         ceval = cb.ConvergenceEvaluation(
#             build_model=self.build_model,
#             input_specification=spec,
#         )
#
#         s2 = ceval.get_specification()
#
#         assert s2 is spec
#         assert ceval._input_spec is spec
#
#     @pytest.mark.unit
#     def test_get_specification_none(self):
#         ceval = cb.ConvergenceEvaluation(
#             build_model=self.build_model,
#         )
#
#         with pytest.raises(
#             ConfigurationError,
#             match="Please specify an input specification to use for sampling.",
#         ):
#             ceval.get_specification()
#
#     @pytest.mark.unit
#     def test_get_input_samples(self):
#         ceval = cb.ConvergenceEvaluation(
#             build_model=self.build_model,
#             input_specification=spec,
#         )
#
#         s = ceval.get_input_samples()
#
#         assert s is spec.samples
#
#     @pytest.mark.unit
#     def test_get_input_samples_none(self):
#         spec2 = cb.ConvergenceEvaluationSpecification()
#         spec2.set_sampling_method(LatinHypercubeSampling)
#         spec2.add_sampled_input("v1", "v1", 0, 10)
#         spec2.add_sampled_input("v2", "v2", 20, 40)
#         spec2.add_sampled_input("p2", "p2", -10, 10)
#
#         ceval = cb.ConvergenceEvaluation(
#             build_model=self.build_model,
#             input_specification=spec2,
#         )
#
#         with pytest.raises(
#             ValueError,
#             match="Please set a sample size.",
#         ):
#             ceval.get_input_samples()
#
#     @pytest.mark.unit
#     def test_set_input_values(self):
#         ceval = cb.ConvergenceEvaluation(
#             build_model=self.build_model,
#             input_specification=spec,
#         )
#
#         model = ceval.get_initialized_model()
#
#         for i in range(5):
#             ceval._set_input_values(model, i)
#
#             assert value(model.v1) == pytest.approx(spec.samples["v1"][i], abs=1e-10)
#             assert value(model.v2) == pytest.approx(spec.samples["v2"][i], abs=1e-10)
#             assert value(model.p2) == pytest.approx(spec.samples["p2"][i], abs=1e-10)
#
#             assert value(model.v3) == pytest.approx(1, abs=1e-10)
#             assert value(model.p1) == pytest.approx(1, abs=1e-10)
#             assert value(model.p3) == pytest.approx(1, abs=1e-10)
#
#     @pytest.mark.unit
#     def test_set_input_immutable_param(self):
#         def bm():
#             m = ConcreteModel()
#             m.p1 = Param()
#             return m
#
#         spec2 = cb.ConvergenceEvaluationSpecification()
#         spec2.set_sampling_method(LatinHypercubeSampling)
#         spec2.add_sampled_input("p1", "p1", 0, 10)
#         spec2.set_sample_size(2)
#
#         ceval = cb.ConvergenceEvaluation(
#             build_model=bm,
#             input_specification=spec2,
#         )
#
#         model = ceval.get_initialized_model()
#
#         with pytest.raises(
#             ValueError,
#             match="Convergence testing found an input of type Param that "
#             "was not mutable \(p1\). Please make sure all "
#             "sampled inputs are either mutable params or fixed vars.",
#         ):
#             ceval._set_input_values(model, 0)
#
#     @pytest.mark.unit
#     def test_set_input_unfixed_var(self):
#         def bm():
#             m = ConcreteModel()
#             m.v1 = Var()
#             return m
#
#         spec2 = cb.ConvergenceEvaluationSpecification()
#         spec2.set_sampling_method(LatinHypercubeSampling)
#         spec2.add_sampled_input("v1", "v1", 0, 10)
#         spec2.set_sample_size(2)
#
#         ceval = cb.ConvergenceEvaluation(
#             build_model=bm,
#             input_specification=spec2,
#         )
#
#         model = ceval.get_initialized_model()
#
#         with pytest.raises(
#             ValueError,
#             match="Convergence testing found an input of type Var that "
#             "was not fixed \(v1\). Please make sure all "
#             "sampled inputs are either mutable params or fixed vars.",
#         ):
#             ceval._set_input_values(model, 0)
#
#     @pytest.mark.unit
#     def test_set_input_invalid_ctype(self):
#         def bm():
#             m = ConcreteModel()
#             m.v1 = Var()
#             m.c1 = Constraint(expr=m.v1 == 2)
#             return m
#
#         spec2 = cb.ConvergenceEvaluationSpecification()
#         spec2.set_sampling_method(LatinHypercubeSampling)
#         spec2.add_sampled_input("c1", "c1", 0, 10)
#         spec2.set_sample_size(2)
#
#         ceval = cb.ConvergenceEvaluation(
#             build_model=bm,
#             input_specification=spec2,
#         )
#
#         model = ceval.get_initialized_model()
#
#         with pytest.raises(
#             ValueError,
#             match="Failed to find a valid input component \(must be "
#             "a fixed Var or a mutable Param\). Instead, "
#             "pyomo_path: c1 returned: c1.",
#         ):
#             ceval._set_input_values(model, 0)
#
#     @pytest.mark.unit
#     def test_get_solver(self):
#         ceval = cb.ConvergenceEvaluation(
#             build_model=self.build_model,
#             input_specification=spec,
#         )
#
#         solver = ceval._get_solver()
#
#         assert isinstance(solver, IPOPT)
#         assert solver.options == {}
#
#     @pytest.mark.unit
#     def test_get_solver_w_options(self):
#         ceval = cb.ConvergenceEvaluation(
#             build_model=self.build_model,
#             input_specification=spec,
#             ipopt_options={"bound_push": 1e-5},
#         )
#
#         solver = ceval._get_solver()
#
#         assert isinstance(solver, IPOPT)
#         assert solver.options == {"bound_push": 1e-5}
#
#     @pytest.mark.component
#     @pytest.mark.solver
#     def test_execute_single_sample(self):
#         def build_model():
#             m = ConcreteModel()
#             m.v1 = Var(initialize=1)
#             m.v2 = Var(initialize=4)
#             m.c1 = Constraint(expr=m.v1 == m.v2)
#
#             m.v2.fix()
#
#             return m
#
#         spec2 = cb.ConvergenceEvaluationSpecification()
#         spec2.set_sampling_method(LatinHypercubeSampling)
#         spec2.add_sampled_input("v2", "v2", 2, 6)
#         spec2.set_sample_size(2)
#
#         ceval = cb.ConvergenceEvaluation(
#             build_model=build_model,
#             input_specification=spec2,
#         )
#
#         results_dict = ceval.execute_single_sample(1)
#
#         assert results_dict["solved"]
#         assert results_dict["iters"] == 1
#         assert results_dict["iters_in_restoration"] == 0
#         assert results_dict["iters_w_regularization"] == 0
#         assert results_dict["time"] < 0.01
#         assert results_dict["numerical_issues"] == 0
#
#     @pytest.mark.component
#     @pytest.mark.solver
#     def test_run_convergence_evaluation(self):
#         def build_model():
#             m = ConcreteModel()
#             m.v1 = Var(initialize=1)
#             m.v2 = Var(initialize=4)
#             m.c1 = Constraint(expr=m.v1 == m.v2)
#
#             m.v2.fix()
#
#             return m
#
#         spec2 = cb.ConvergenceEvaluationSpecification()
#         spec2.set_sampling_method(UniformSampling)
#         spec2.add_sampled_input("v2", "v2", 2, 6)
#         spec2.set_sample_size([2])
#
#         ceval = cb.ConvergenceEvaluation(
#             build_model=build_model,
#             input_specification=spec2,
#         )
#
#         results = ceval.run_convergence_evaluation()
#
#         assert isinstance(results, dict)
#         assert len(results) == 2
#
#         for i in [0, 1]:
#             assert results[i]["solved"]
#             assert results[i]["iters"] == 1
#             assert results[i]["iters_in_restoration"] == 0
#             assert results[i]["iters_w_regularization"] == 0
#             assert results[i]["time"] < 0.01
#             assert results[i]["numerical_issues"] == 0
#
#         assert results is ceval.results
#
#     @pytest.fixture(scope="class")
#     def ceval_executed(self):
#         def build_model():
#             m = ConcreteModel()
#             m.v1 = Var(initialize=1)
#             m.v2 = Var(initialize=4)
#             m.c1 = Constraint(expr=m.v1 == m.v2)
#
#             m.v2.fix()
#
#             return m
#
#         spec2 = cb.ConvergenceEvaluationSpecification()
#         spec2.set_sampling_method(UniformSampling)
#         spec2.add_sampled_input("v2", "v2", 2, 6)
#         spec2.set_sample_size([2])
#
#         ceval = cb.ConvergenceEvaluation(
#             build_model=build_model,
#             input_specification=spec2,
#         )
#
#         ceval.run_convergence_evaluation()
#
#         return ceval
#
#     @pytest.mark.component
#     @pytest.mark.solver
#     def test_to_dict(self, ceval_executed):
#         outdict = ceval_executed.to_dict()
#         assert outdict == ceval_dict
#
#     @pytest.mark.unit
#     def test_from_dict(self):
#         def build_model():
#             m = ConcreteModel()
#             m.v1 = Var(initialize=1)
#             m.v2 = Var(initialize=4)
#             m.c1 = Constraint(expr=m.v1 == m.v2)
#
#             m.v2.fix()
#
#             return m
#
#         ceval = cb.ConvergenceEvaluation(
#             build_model=build_model,
#         )
#
#         ceval.from_dict(ceval_dict)
#
#         assert isinstance(ceval._input_spec, cb.ConvergenceEvaluationSpecification)
#         assert len(ceval._input_spec.inputs) == 1
#
#         assert ceval._input_spec.sampling_method is UniformSampling
#         assert isinstance(ceval._input_spec.samples, DataFrame)
#         assert ceval._input_spec.sample_size == [2]
#
#         assert isinstance(ceval._results, dict)
#         assert len(ceval._results) == 2
#
#         for i in [0, 1]:
#             assert ceval._results[i]["solved"]
#             assert ceval._results[i]["iters"] == 1
#             assert ceval._results[i]["iters_in_restoration"] == 0
#             assert ceval._results[i]["iters_w_regularization"] == 0
#             assert ceval._results[i]["time"] < 0.01
#             assert ceval._results[i]["numerical_issues"] == 0
#
#     @pytest.mark.component
#     def test_write_to_json_file(self, ceval_executed):
#         temp_context = TempfileManager.new_context()
#         tmpfile = temp_context.create_tempfile(suffix=".json")
#
#         ceval_executed.write_to_json_file(tmpfile)
#
#         with open(tmpfile, "r") as f:
#             lines = f.read()
#         f.close()
#
#         expected = """{
#    "specification": {
#       "inputs": {
#          "v2": {
#             "pyomo_path": "v2",
#             "lower": 2,
#             "upper": 6,
#             "mean": null,
#             "std": null,
#             "distribution": "normal"
#          }
#       },
#       "sampling_method": "UniformSampling",
#       "sample_size": [
#          2
#       ],
#       "samples": {
#          "index": [
#             0,
#             1
#          ],
#          "columns": [
#             "v2"
#          ],
#          "data": [
#             [
#                2.0
#             ],
#             [
#                6.0
#             ]
#          ],
#          "index_names": [
#             null
#          ],
#          "column_names": [
#             null
#          ]
#       }
#    },
#    "results": {
#       "0": {
#          "solved": true,
#          "iters": 1,
#          "iters_in_restoration": 0,
#          "iters_w_regularization": 0,
#          "time": 0.0,
#          "numerical_issues": 0
#       },
#       "1": {
#          "solved": true,
#          "iters": 1,
#          "iters_in_restoration": 0,
#          "iters_w_regularization": 0,
#          "time": 0.0,
#          "numerical_issues": 0
#       }
#    }
# }"""
#
#         assert lines == expected
#
#         # Check for clean up
#         temp_context.release(remove=True)
#         assert not os.path.exists(tmpfile)
#
#     @pytest.mark.unit
#     def test_load_from_json_file(self):
#         fname = os.path.join(currdir, "load_ceval.json")
#
#         ceval = cb.ConvergenceEvaluation()
#         ceval.load_from_json_file(fname)
#
#         assert ceval.config.build_model is None
#
#         assert isinstance(ceval._input_spec, cb.ConvergenceEvaluationSpecification)
#         assert len(ceval._input_spec.inputs) == 1
#
#         assert ceval._input_spec.sampling_method is UniformSampling
#         assert isinstance(ceval._input_spec.samples, DataFrame)
#         assert ceval._input_spec.sample_size == [2]
#
#         assert isinstance(ceval._results, dict)
#         assert len(ceval._results) == 2
#
#         for i in [0, 1]:
#             assert ceval._results[i]["solved"]
#             assert ceval._results[i]["iters"] == 1
#             assert ceval._results[i]["iters_in_restoration"] == 0
#             assert ceval._results[i]["iters_w_regularization"] == 0
#             assert ceval._results[i]["time"] < 0.01
#             assert ceval._results[i]["numerical_issues"] == 0
#
#     @pytest.mark.component
#     @pytest.mark.solver
#     def test_verify_samples_from_dict(self, ceval_executed):
#         ceval_executed._verify_samples_from_dict(ceval_dict)
#
#         with pytest.raises(
#             ValueError,
#             match="Samples in comparison evaluation do not match current evaluation",
#         ):
#             fail_dict = {
#                 "specification": {
#                     "samples": {
#                         "index": [0, 1],
#                         "columns": ["v2"],
#                         "data": [[1.0], [6.0]],  # Changed sample value from 2 to 1
#                         "index_names": [None],
#                         "column_names": [None],
#                     },
#                 },
#             }
#             ceval_executed._verify_samples_from_dict(fail_dict)
#
#     @pytest.mark.component
#     @pytest.mark.solver
#     def test_compare_results_to_dict_pass(self, ceval_executed):
#         (
#             solved_diff,
#             iters_diff,
#             restore_diff,
#             reg_diff,
#             num_iss_diff,
#         ) = ceval_executed.compare_results_to_dict(ceval_dict)
#
#         assert solved_diff == []
#         assert iters_diff == []
#         assert restore_diff == []
#         assert reg_diff == []
#         assert num_iss_diff == []
#
#     @pytest.mark.component
#     @pytest.mark.solver
#     def test_compare_results_to_dict_solved(self, ceval_executed):
#         fail_dict = {
#             "results": OrderedDict(
#                 [
#                     (
#                         0,
#                         {
#                             "solved": False,
#                             "iters": 1,
#                             "iters_in_restoration": 0,
#                             "iters_w_regularization": 0,
#                             "time": 0.0,
#                             "numerical_issues": 0,
#                         },
#                     ),
#                     (
#                         1,
#                         {
#                             "solved": True,
#                             "iters": 1,
#                             "iters_in_restoration": 0,
#                             "iters_w_regularization": 0,
#                             "time": 0.0,
#                             "numerical_issues": 0,
#                         },
#                     ),
#                 ]
#             ),
#         }
#
#         (
#             solved_diff,
#             iters_diff,
#             restore_diff,
#             reg_diff,
#             num_iss_diff,
#         ) = ceval_executed.compare_results_to_dict(fail_dict, verify_samples=False)
#
#         assert solved_diff == [0]
#         assert iters_diff == []
#         assert restore_diff == []
#         assert reg_diff == []
#         assert num_iss_diff == []
#
#     @pytest.mark.component
#     @pytest.mark.solver
#     def test_compare_results_to_dict_iters(self, ceval_executed):
#         fail_dict = {
#             "results": OrderedDict(
#                 [
#                     (
#                         0,
#                         {
#                             "solved": True,
#                             "iters": 2,
#                             "iters_in_restoration": 0,
#                             "iters_w_regularization": 0,
#                             "time": 0.0,
#                             "numerical_issues": 0,
#                         },
#                     ),
#                     (
#                         1,
#                         {
#                             "solved": True,
#                             "iters": 3,
#                             "iters_in_restoration": 0,
#                             "iters_w_regularization": 0,
#                             "time": 0.0,
#                             "numerical_issues": 0,
#                         },
#                     ),
#                 ]
#             ),
#         }
#
#         (
#             solved_diff,
#             iters_diff,
#             restore_diff,
#             reg_diff,
#             num_iss_diff,
#         ) = ceval_executed.compare_results_to_dict(fail_dict, verify_samples=False)
#
#         assert solved_diff == []
#         assert iters_diff == [1]
#         assert restore_diff == []
#         assert reg_diff == []
#         assert num_iss_diff == []
#
#     @pytest.mark.component
#     @pytest.mark.solver
#     def test_compare_results_to_dict_restoration(self, ceval_executed):
#         fail_dict = {
#             "results": OrderedDict(
#                 [
#                     (
#                         0,
#                         {
#                             "solved": True,
#                             "iters": 1,
#                             "iters_in_restoration": 1,
#                             "iters_w_regularization": 0,
#                             "time": 0.0,
#                             "numerical_issues": 0,
#                         },
#                     ),
#                     (
#                         1,
#                         {
#                             "solved": True,
#                             "iters": 1,
#                             "iters_in_restoration": 2,
#                             "iters_w_regularization": 0,
#                             "time": 0.0,
#                             "numerical_issues": 0,
#                         },
#                     ),
#                 ]
#             ),
#         }
#
#         (
#             solved_diff,
#             iters_diff,
#             restore_diff,
#             reg_diff,
#             num_iss_diff,
#         ) = ceval_executed.compare_results_to_dict(fail_dict, verify_samples=False)
#
#         assert solved_diff == []
#         assert iters_diff == []
#         assert restore_diff == [1]
#         assert reg_diff == []
#         assert num_iss_diff == []
#
#     @pytest.mark.component
#     @pytest.mark.solver
#     def test_compare_results_to_dict_regularization(self, ceval_executed):
#         fail_dict = {
#             "results": OrderedDict(
#                 [
#                     (
#                         0,
#                         {
#                             "solved": True,
#                             "iters": 1,
#                             "iters_in_restoration": 0,
#                             "iters_w_regularization": 2,
#                             "time": 0.0,
#                             "numerical_issues": 0,
#                         },
#                     ),
#                     (
#                         1,
#                         {
#                             "solved": True,
#                             "iters": 1,
#                             "iters_in_restoration": 0,
#                             "iters_w_regularization": 1,
#                             "time": 0.0,
#                             "numerical_issues": 0,
#                         },
#                     ),
#                 ]
#             ),
#         }
#
#         (
#             solved_diff,
#             iters_diff,
#             restore_diff,
#             reg_diff,
#             num_iss_diff,
#         ) = ceval_executed.compare_results_to_dict(fail_dict, verify_samples=False)
#
#         assert solved_diff == []
#         assert iters_diff == []
#         assert restore_diff == []
#         assert reg_diff == [0]
#         assert num_iss_diff == []
#
#     @pytest.mark.component
#     @pytest.mark.solver
#     def test_compare_results_to_dict_numerical_issues(self, ceval_executed):
#         fail_dict = {
#             "results": OrderedDict(
#                 [
#                     (
#                         0,
#                         {
#                             "solved": True,
#                             "iters": 1,
#                             "iters_in_restoration": 0,
#                             "iters_w_regularization": 0,
#                             "time": 0.0,
#                             "numerical_issues": 1,
#                         },
#                     ),
#                     (
#                         1,
#                         {
#                             "solved": True,
#                             "iters": 1,
#                             "iters_in_restoration": 0,
#                             "iters_w_regularization": 0,
#                             "time": 0.0,
#                             "numerical_issues": 0,
#                         },
#                     ),
#                 ]
#             ),
#         }
#
#         (
#             solved_diff,
#             iters_diff,
#             restore_diff,
#             reg_diff,
#             num_iss_diff,
#         ) = ceval_executed.compare_results_to_dict(fail_dict, verify_samples=False)
#
#         assert solved_diff == []
#         assert iters_diff == []
#         assert restore_diff == []
#         assert reg_diff == []
#         assert num_iss_diff == [0]
#
#     @pytest.mark.component
#     @pytest.mark.solver
#     def test_compare_results_to_json_file(self, ceval_executed):
#         fname = os.path.join(currdir, "load_ceval.json")
#
#         (
#             solved_diff,
#             iters_diff,
#             restore_diff,
#             reg_diff,
#             num_iss_diff,
#         ) = ceval_executed.compare_results_to_json_file(fname)
#
#         assert solved_diff == []
#         assert iters_diff == []
#         assert restore_diff == []
#         assert reg_diff == []
#         assert num_iss_diff == []
#
#     @pytest.mark.component
#     @pytest.mark.solver
#     def test_run_baseline_comparison(self):
#         fname = os.path.join(currdir, "load_ceval.json")
#
#         def build_model():
#             m = ConcreteModel()
#             m.v1 = Var(initialize=1)
#             m.v2 = Var(initialize=4)
#             m.c1 = Constraint(expr=m.v1 == m.v2)
#
#             m.v2.fix()
#
#             return m
#
#         ceval = cb.ConvergenceEvaluation(
#             build_model=build_model,
#         )
#         (
#             solved_diff,
#             iters_diff,
#             restore_diff,
#             reg_diff,
#             num_iss_diff,
#         ) = ceval.run_baseline_comparison(fname)
#
#         assert solved_diff == []
#         assert iters_diff == []
#         assert restore_diff == []
#         assert reg_diff == []
#         assert num_iss_diff == []
