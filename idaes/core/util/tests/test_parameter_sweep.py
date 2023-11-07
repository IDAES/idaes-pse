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

from pyomo.environ import (
    ConcreteModel,
    Constraint,
    Var,
    Param,
    value,
    SolverFactory,
    assert_optimal_termination,
)
from pyomo.common.fileutils import this_file_dir
from pyomo.solvers.plugins.solvers.IPOPT import IPOPT
from pyomo.common.tempfiles import TempfileManager

from idaes.core.util.parameter_sweep import (
    ParameterSweepSpecification,
    ParameterSweepBase,
    SequentialSweepRunner,
)
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


spec = ParameterSweepSpecification()
spec.set_sampling_method(LatinHypercubeSampling)
spec.add_sampled_input("v1", "v1", 0, 10)
spec.add_sampled_input("v2", "v2", 20, 40)
spec.add_sampled_input("p2", "p2", -10, 10)
spec.set_sample_size(5)
spec.generate_samples()


psweep_dict = {
    "specification": {
        "inputs": OrderedDict(
            [
                (
                    "v2",
                    OrderedDict(
                        [
                            ("pyomo_path", "v2"),
                            ("lower", 2),
                            ("upper", 6),
                        ]
                    ),
                )
            ]
        ),
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
    "results": OrderedDict(
        [
            (
                0,
                {
                    "solved": True,
                    "results": 2,
                },
            ),
            (
                1,
                {
                    "solved": True,
                    "results": 6,
                },
            ),
        ]
    ),
}


class TestParameterSweepBase:
    def build_model(self):
        m = ConcreteModel()
        m.v1 = Var(initialize=1)
        m.v2 = Var(initialize=1)
        m.v3 = Var(initialize=1)
        m.p1 = Param(initialize=1)
        m.p2 = Param(initialize=1, mutable=True)
        m.p3 = Param(initialize=1, mutable=True)

        m.v1.fix(1)
        m.v2.fix(1)

        return m

    @pytest.mark.unit
    def test_init(self):
        psweep = ParameterSweepBase()

        assert psweep.config.build_model is None
        assert psweep.config.run_model is None
        assert psweep.config.collect_results is None
        assert psweep.config.failure_recourse is None
        assert psweep.config.input_specification is None
        assert psweep.config.solver is None
        assert psweep.config.solver_options is None
        assert isinstance(psweep.results, OrderedDict)
        assert len(psweep.results) == 0

    @pytest.mark.unit
    def test_get_initialized_model_none(self):
        psweep = ParameterSweepBase()

        with pytest.raises(
            ConfigurationError,
            match="Please specify a method to construct the model of interest.",
        ):
            psweep.get_initialized_model()

    @pytest.mark.unit
    def test_get_initialized_model(self):
        psweep = ParameterSweepBase(
            build_model=self.build_model,
        )

        m2 = psweep.get_initialized_model()

        assert isinstance(m2, ConcreteModel)
        assert isinstance(m2.v1, Var)
        assert isinstance(m2.v2, Var)
        assert isinstance(m2.p1, Param)
        assert isinstance(m2.p2, Param)
        assert isinstance(m2.p3, Param)

    @pytest.mark.unit
    def test_get_specification(self):
        psweep = ParameterSweepBase(
            build_model=self.build_model,
            input_specification=spec,
        )

        s2 = psweep.get_input_specification()

        assert s2 is spec

    @pytest.mark.unit
    def test_get_specification_none(self):
        psweep = ParameterSweepBase(
            build_model=self.build_model,
        )

        with pytest.raises(
            ConfigurationError,
            match="Please specify an input specification to use for sampling.",
        ):
            psweep.get_input_specification()

    @pytest.mark.unit
    def test_get_input_samples(self):
        psweep = ParameterSweepBase(
            build_model=self.build_model,
            input_specification=spec,
        )

        s = psweep.get_input_samples()

        assert s is spec.samples

    @pytest.mark.unit
    def test_get_input_samples_none(self):
        spec2 = ParameterSweepSpecification()
        spec2.set_sampling_method(LatinHypercubeSampling)
        spec2.add_sampled_input("v1", "v1", 0, 10)
        spec2.add_sampled_input("v2", "v2", 20, 40)
        spec2.add_sampled_input("p2", "p2", -10, 10)

        psweep = ParameterSweepBase(
            build_model=self.build_model,
            input_specification=spec2,
        )

        with pytest.raises(
            ValueError,
            match="Please set a sample size.",
        ):
            psweep.get_input_samples()

    @pytest.mark.unit
    def test_set_input_values(self):
        psweep = ParameterSweepBase(
            build_model=self.build_model,
            input_specification=spec,
        )

        model = psweep.get_initialized_model()

        for i in range(5):
            psweep.set_input_values(model, i)

            assert value(model.v1) == pytest.approx(spec.samples["v1"][i], abs=1e-10)
            assert value(model.v2) == pytest.approx(spec.samples["v2"][i], abs=1e-10)
            assert value(model.p2) == pytest.approx(spec.samples["p2"][i], abs=1e-10)

            assert value(model.v3) == pytest.approx(1, abs=1e-10)
            assert value(model.p1) == pytest.approx(1, abs=1e-10)
            assert value(model.p3) == pytest.approx(1, abs=1e-10)

    @pytest.mark.unit
    def test_set_input_immutable_param(self):
        def bm():
            m = ConcreteModel()
            m.p1 = Param()
            return m

        spec2 = ParameterSweepSpecification()
        spec2.set_sampling_method(LatinHypercubeSampling)
        spec2.add_sampled_input("p1", "p1", 0, 10)
        spec2.set_sample_size(2)

        psweep = ParameterSweepBase(
            build_model=bm,
            input_specification=spec2,
        )

        model = psweep.get_initialized_model()

        with pytest.raises(
            ValueError,
            match="Convergence testing found an input of type Param that "
            "was not mutable \(p1\). Please make sure all "
            "sampled inputs are either mutable params or fixed vars.",
        ):
            psweep.set_input_values(model, 0)

    @pytest.mark.unit
    def test_set_input_unfixed_var(self):
        def bm():
            m = ConcreteModel()
            m.v1 = Var()
            return m

        spec2 = ParameterSweepSpecification()
        spec2.set_sampling_method(LatinHypercubeSampling)
        spec2.add_sampled_input("v1", "v1", 0, 10)
        spec2.set_sample_size(2)

        psweep = ParameterSweepBase(
            build_model=bm,
            input_specification=spec2,
        )

        model = psweep.get_initialized_model()

        with pytest.raises(
            ValueError,
            match="Convergence testing found an input of type Var that "
            "was not fixed \(v1\). Please make sure all "
            "sampled inputs are either mutable params or fixed vars.",
        ):
            psweep.set_input_values(model, 0)

    @pytest.mark.unit
    def test_set_input_invalid_ctype(self):
        def bm():
            m = ConcreteModel()
            m.v1 = Var()
            m.c1 = Constraint(expr=m.v1 == 2)
            return m

        spec2 = ParameterSweepSpecification()
        spec2.set_sampling_method(LatinHypercubeSampling)
        spec2.add_sampled_input("c1", "c1", 0, 10)
        spec2.set_sample_size(2)

        psweep = ParameterSweepBase(
            build_model=bm,
            input_specification=spec2,
        )

        model = psweep.get_initialized_model()

        with pytest.raises(
            ValueError,
            match="Failed to find a valid input component \(must be "
            "a fixed Var or a mutable Param\). Instead, "
            "pyomo_path: c1 returned: c1.",
        ):
            psweep.set_input_values(model, 0)

    @pytest.mark.unit
    def test_get_solver(self):
        psweep = ParameterSweepBase(
            build_model=self.build_model,
            input_specification=spec,
        )

        with pytest.raises(
            ConfigurationError,
            match="Please specify a solver to use.",
        ):
            psweep._get_solver()

    @pytest.mark.unit
    def test_get_solver_w_options(self):
        psweep = ParameterSweepBase(
            build_model=self.build_model,
            input_specification=spec,
            solver="ipopt",
            solver_options={"bound_push": 1e-5},
        )

        solver = psweep._get_solver()

        assert isinstance(solver, IPOPT)
        assert solver.options == {"bound_push": 1e-5}

    @pytest.mark.component
    @pytest.mark.solver
    def test_run_model_none(self):
        psweep = ParameterSweepBase()

        model = ConcreteModel()
        model.v = Var()
        model.c = Constraint(expr=model.v == 4)

        solver = SolverFactory("ipopt")

        status = psweep.run_model(model, solver)

        assert_optimal_termination(status)
        assert value(model.v) == pytest.approx(4, rel=1e-8)

    @pytest.mark.unit
    def test_run_model_dummy(self):
        def dummy_run(model, solver):
            return "foo"

        psweep = ParameterSweepBase(run_model=dummy_run)

        model = ConcreteModel()
        status = psweep.run_model(model, "bar")

        assert status == "foo"

    @pytest.mark.unit
    def test_collect_results_none(self):
        psweep = ParameterSweepBase()

        with pytest.raises(
            ConfigurationError,
            match="Please provide a method to collect results from sample run.",
        ):
            psweep.collect_results("foo", "bar")

    @pytest.mark.unit
    def test_collect_results(self):
        def dummy_collect(model, status):
            return value(model.v)

        model = ConcreteModel()
        model.v = Var(initialize=1)
        model.c = Constraint(expr=model.v == 4)

        psweep = ParameterSweepBase(
            collect_results=dummy_collect,
        )

        results = psweep.collect_results(model, "foo")

        assert results == 1

    @pytest.mark.component
    @pytest.mark.solver
    def test_execute_single_sample(self):
        def build_model():
            m = ConcreteModel()
            m.v1 = Var(initialize=1)
            m.v2 = Var(initialize=4)
            m.c1 = Constraint(expr=m.v1 == m.v2)

            m.v2.fix()

            return m

        def collect_results(model, status):
            return value(model.v1)

        spec2 = ParameterSweepSpecification()
        spec2.set_sampling_method(UniformSampling)
        spec2.add_sampled_input("v2", "v2", 2, 6)
        spec2.set_sample_size([2])

        psweep = ParameterSweepBase(
            build_model=build_model,
            input_specification=spec2,
            collect_results=collect_results,
            solver="ipopt",
        )

        results, solved = psweep.execute_single_sample(1)

        assert results == pytest.approx(6, rel=1e-8)
        assert solved

    @pytest.mark.component
    def test_execute_single_sample_recourse_none(self):
        def build_model():
            m = ConcreteModel()
            m.v1 = Var(initialize=1)
            m.v2 = Var(initialize=4)
            m.c1 = Constraint(expr=m.v1 == m.v2)

            m.v2.fix()

            return m

        def dummy_solver(model):
            raise Exception

        def dummy_get_solver():
            return dummy_solver

        spec2 = ParameterSweepSpecification()
        spec2.set_sampling_method(UniformSampling)
        spec2.add_sampled_input("v2", "v2", 2, 6)
        spec2.set_sample_size([2])

        psweep = ParameterSweepBase(
            build_model=build_model,
            input_specification=spec2,
        )

        psweep._get_solver = dummy_get_solver

        results, solved = psweep.execute_single_sample(1)

        assert results is None
        assert not solved

    @pytest.mark.component
    def test_execute_single_sample_recourse_custom(self):
        def build_model():
            m = ConcreteModel()
            m.v1 = Var(initialize=1)
            m.v2 = Var(initialize=4)
            m.c1 = Constraint(expr=m.v1 == m.v2)

            m.v2.fix()

            return m

        def dummy_solver(model):
            raise Exception

        def dummy_get_solver():
            return dummy_solver

        def recourse(model):
            return "foo", "bar"

        spec2 = ParameterSweepSpecification()
        spec2.set_sampling_method(UniformSampling)
        spec2.add_sampled_input("v2", "v2", 2, 6)
        spec2.set_sample_size([2])

        psweep = ParameterSweepBase(
            build_model=build_model,
            input_specification=spec2,
            failure_recourse=recourse,
        )

        psweep._get_solver = dummy_get_solver

        results, solved = psweep.execute_single_sample(1)

        assert results == "foo"
        assert solved == "bar"

    @pytest.fixture(scope="class")
    def psweep_with_results(self):
        spec2 = ParameterSweepSpecification()
        spec2.set_sampling_method(UniformSampling)
        spec2.add_sampled_input("v2", "v2", 2, 6)
        spec2.set_sample_size([2])
        spec2.generate_samples()

        psweep = ParameterSweepBase(
            input_specification=spec2,
        )

        psweep._results = OrderedDict(
            {
                0: {"solved": True, "results": 2},
                1: {"solved": True, "results": 6},
            }
        )

        return psweep

    @pytest.mark.component
    @pytest.mark.solver
    def test_to_dict(self, psweep_with_results):
        outdict = psweep_with_results.to_dict()
        assert outdict == psweep_dict

    @pytest.mark.unit
    def test_from_dict(self):
        def build_model():
            m = ConcreteModel()
            m.v1 = Var(initialize=1)
            m.v2 = Var(initialize=4)
            m.c1 = Constraint(expr=m.v1 == m.v2)

            m.v2.fix()

            return m

        psweep = ParameterSweepBase(
            build_model=build_model,
        )

        psweep.from_dict(psweep_dict)

        assert isinstance(psweep._input_spec, ParameterSweepSpecification)
        assert len(psweep._input_spec.inputs) == 1

        assert psweep._input_spec.sampling_method is UniformSampling
        assert isinstance(psweep._input_spec.samples, DataFrame)
        assert psweep._input_spec.sample_size == [2]

        assert isinstance(psweep._results, dict)
        assert len(psweep._results) == 2

        for i in [0, 1]:
            assert psweep._results[i]["solved"]
            assert psweep._results[i]["results"] == 2 + i * 4

    @pytest.mark.component
    def test_to_json_file(self, psweep_with_results):
        temp_context = TempfileManager.new_context()
        tmpfile = temp_context.create_tempfile(suffix=".json")

        psweep_with_results.to_json_file(tmpfile)

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
         "solved": true,
         "results": 2
      },
      "1": {
         "solved": true,
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
        fname = os.path.join(currdir, "load_psweep.json")

        psweep = ParameterSweepBase()
        psweep.from_json_file(fname)

        assert psweep.config.build_model is None

        assert isinstance(psweep._input_spec, ParameterSweepSpecification)
        assert len(psweep._input_spec.inputs) == 1

        assert psweep._input_spec.sampling_method is UniformSampling
        assert isinstance(psweep._input_spec.samples, DataFrame)
        assert psweep._input_spec.sample_size == [2]

        assert isinstance(psweep._results, dict)
        assert len(psweep._results) == 2

        for i in [0, 1]:
            assert psweep._results[i]["solved"]
            assert psweep._results[i]["results"] == 2 + i * 4

    @pytest.mark.unit
    def test_execute_parameter_sweep(self):
        psweep = ParameterSweepBase()

        with pytest.raises(NotImplementedError):
            psweep.execute_parameter_sweep()


class TestSequentialSweepRunner:
    @pytest.mark.component
    @pytest.mark.solver
    def test_sequential_runner(self):
        def build_model():
            m = ConcreteModel()
            m.v1 = Var(initialize=1)
            m.v2 = Var(initialize=4)
            m.c1 = Constraint(expr=m.v1 == m.v2 - 1e-3)

            m.v2.fix()

            return m

        def collect_results(model, status):
            return value(model.v1)

        spec2 = ParameterSweepSpecification()
        spec2.set_sampling_method(UniformSampling)
        spec2.add_sampled_input("v2", "v2", 2, 6)
        spec2.set_sample_size([2])
        spec2.generate_samples()

        psweep = SequentialSweepRunner(
            build_model=build_model,
            input_specification=spec2,
            solver="ipopt",
            collect_results=collect_results,
        )

        psweep.execute_parameter_sweep()

        assert psweep.results[0]["solved"]
        assert psweep.results[0]["results"] == pytest.approx(2 - 1e-3, rel=1e-8)

        assert psweep.results[1]["solved"]
        assert psweep.results[1]["results"] == pytest.approx(6 - 1e-3, rel=1e-8)
