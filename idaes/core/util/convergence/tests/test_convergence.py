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
Tests for the convergence testing module

Author: Carl Laird
"""
import io
import json
import pytest
import os
import os.path
from collections import OrderedDict

import pyomo.environ as pe
from pyomo.common.fileutils import this_file_dir
from pyomo.common.unittest import assertStructuredAlmostEqual
import idaes.core.util.convergence.convergence_base as cb
import idaes

# See if ipopt is available and set up solver
ipopt_available = pe.SolverFactory("ipopt").available()
ceval_fixedvar_mutableparam_str = (
    "idaes.core.util.convergence.tests."
    "conv_eval_classes.ConvEvalFixedVarMutableParam"
)
ceval_fixedvar_immutableparam_str = (
    "idaes.core.util.convergence.tests."
    "conv_eval_classes.ConvEvalFixedVarImmutableParam"
)
ceval_unfixedvar_mutableparam_str = (
    "idaes.core.util.convergence.tests."
    "conv_eval_classes.ConvEvalUnfixedVarMutableParam"
)

currdir = this_file_dir()
wrtdir = idaes.testing_directory


@pytest.mark.unit
def test_convergence_evaluation_specification_file_fixedvar_mutableparam():
    ceval_class = cb._class_import(ceval_fixedvar_mutableparam_str)
    ceval = ceval_class()
    import idaes.core.util.convergence.tests.conv_eval_classes as cev

    assert ceval.__class__ == cev.ConvEvalFixedVarMutableParam().__class__

    spec = ceval.get_specification()
    fname = os.path.join(wrtdir, "ceval_fixedvar_mutableparam.3.42.json")
    cb.write_sample_file(
        spec, fname, ceval_fixedvar_mutableparam_str, n_points=3, seed=42
    )

    baseline_fname = os.path.join(
        currdir, "ceval_fixedvar_mutableparam.3.42.baseline.json"
    )

    with open(fname) as FILE:
        test = json.load(FILE)
    with open(baseline_fname) as FILE:
        baseline = json.load(FILE)
    assertStructuredAlmostEqual(baseline, test, abstol=1e-8)

    if os.path.exists(fname):
        os.remove(fname)


@pytest.mark.unit
def test_convergence_evaluation_specification_file_unfixedvar_mutableparam():
    ceval_class = cb._class_import(ceval_unfixedvar_mutableparam_str)
    ceval = ceval_class()
    import idaes.core.util.convergence.tests.conv_eval_classes as cev

    assert ceval.__class__ == cev.ConvEvalUnfixedVarMutableParam().__class__

    spec = ceval.get_specification()
    fname = os.path.join(wrtdir, "ceval_unfixedvar_mutableparam.3.42.json")
    cb.write_sample_file(
        spec, fname, ceval_unfixedvar_mutableparam_str, n_points=3, seed=42
    )

    baseline_fname = os.path.join(
        currdir, "ceval_unfixedvar_mutableparam.3.42.baseline.json"
    )

    with open(fname) as FILE:
        test = json.load(FILE)
    with open(baseline_fname) as FILE:
        baseline = json.load(FILE)
    assertStructuredAlmostEqual(baseline, test, abstol=1e-8)

    # expect an exception because var is not fixed
    with pytest.raises(ValueError):
        (
            inputs,
            samples,
            global_results,
        ) = cb.run_convergence_evaluation_from_sample_file(fname)

    if os.path.exists(fname):
        os.remove(fname)


@pytest.mark.unit
def test_convergence_evaluation_stats_from_to():
    ceval_class = cb._class_import(ceval_fixedvar_mutableparam_str)
    ceval = ceval_class()
    import idaes.core.util.convergence.tests.conv_eval_classes as cev

    assert ceval.__class__ == cev.ConvEvalFixedVarMutableParam().__class__

    spec = ceval.get_specification()
    fname = os.path.join(wrtdir, "ceval_fixedvar_mutableparam.3.43.json")
    cb.write_sample_file(
        spec, fname, ceval_fixedvar_mutableparam_str, n_points=3, seed=43
    )

    inputs, samples, results = cb.run_convergence_evaluation_from_sample_file(fname)

    s = cb.save_convergence_statistics(inputs, results)
    d = s.to_dict()
    s2 = cb.Stats(from_dict=d)
    fname1 = os.path.join(wrtdir, "stats1.json")

    with open(fname1, "w") as f:
        s.to_json(f)
    with open(fname1) as FILE:
        baseline = json.load(FILE)

    buf = io.StringIO()
    s2.to_json(buf)
    test = json.loads(buf.getvalue())
    assertStructuredAlmostEqual(baseline, test, abstol=1e-8)

    s3 = cb.Stats(from_json=fname1)
    buf = io.StringIO()
    s3.to_json(buf)
    test = json.loads(buf.getvalue())
    assertStructuredAlmostEqual(baseline, test, abstol=1e-8)

    if os.path.exists(fname):
        os.remove(fname)
    if os.path.exists(fname1):
        os.remove(fname1)


@pytest.mark.skipif(not ipopt_available, reason="Ipopt solver not available")
@pytest.mark.unit
def test_convergence_evaluation_single():
    ceval_class = cb._class_import(ceval_fixedvar_mutableparam_str)
    ceval = ceval_class()
    import idaes.core.util.convergence.tests.conv_eval_classes as cev

    assert ceval.__class__ == cev.ConvEvalFixedVarMutableParam().__class__

    spec = ceval.get_specification()
    fname = os.path.join(wrtdir, "ceval_fixedvar_mutableparam.3.43.json")
    cb.write_sample_file(
        spec, fname, ceval_fixedvar_mutableparam_str, n_points=3, seed=43
    )
    s, c, i, rest, reg, i = cb.run_single_sample_from_sample_file(
        fname, name="Sample-1"
    )
    assert c == True


@pytest.mark.unit
def test_convergence_evaluation_specification_file_fixedvar_immutableparam():
    ceval_class = cb._class_import(ceval_fixedvar_immutableparam_str)
    ceval = ceval_class()
    import idaes.core.util.convergence.tests.conv_eval_classes as cev

    assert ceval.__class__ == cev.ConvEvalFixedVarImmutableParam().__class__

    spec = ceval.get_specification()
    fname = os.path.join(wrtdir, "ceval_fixedvar_immutableparam.3.42.json")
    cb.write_sample_file(
        spec, fname, ceval_fixedvar_immutableparam_str, n_points=3, seed=42
    )

    baseline_fname = os.path.join(
        currdir, "ceval_fixedvar_immutableparam.3.42.baseline.json"
    )

    with open(fname) as FILE:
        test = json.load(FILE)
    with open(baseline_fname) as FILE:
        baseline = json.load(FILE)
    assertStructuredAlmostEqual(baseline, test, abstol=1e-8)

    # expect an exception because param is not mutable
    with pytest.raises(ValueError):
        (
            inputs,
            samples,
            global_results,
        ) = cb.run_convergence_evaluation_from_sample_file(fname)

    if os.path.exists(fname):
        os.remove(fname)


@pytest.mark.skipif(not ipopt_available, reason="Ipopt solver not available")
@pytest.mark.unit
def test_convergence_evaluation_fixedvar_mutableparam():
    ceval_class = cb._class_import(ceval_fixedvar_mutableparam_str)
    ceval = ceval_class()
    import idaes.core.util.convergence.tests.conv_eval_classes as cev

    assert ceval.__class__ == cev.ConvEvalFixedVarMutableParam().__class__

    spec = ceval.get_specification()
    fname = os.path.join(currdir, "ceval_fixedvar_mutableparam.3.43.json")
    cb.write_sample_file(
        spec, fname, ceval_fixedvar_mutableparam_str, n_points=3, seed=43
    )

    inputs, samples, global_results = cb.run_convergence_evaluation_from_sample_file(
        fname
    )

    # put the results into a json file for comparison
    # jsondict = dict(inputs=inputs, samples=samples,
    #                 global_results=global_results)
    # results_fname = os.path.join(
    #        currdir, 'ceval_fixedvar_mutableparam.3.43.results.json')
    # with open(results_fname, 'w') as fd:
    #     json.dump(jsondict, fd, indent=3)

    # compare results
    assert global_results[0]["name"] == "Sample-1"
    assert global_results[0]["solved"]
    # This should take 14 iterations to converge with the IDAES solver
    # distribution, but due to various solver factors the number of iterations
    # could vary.
    assert global_results[0]["iters"] == pytest.approx(14, abs=2)

    assert global_results[1]["name"] == "Sample-2"
    assert global_results[1]["solved"]
    # This should take 15 iterations to converge with the IDAES solver
    # distribution, but due to various solver factors the number of iterations
    # could vary.
    assert global_results[1]["iters"] == pytest.approx(15, abs=2)

    assert global_results[2]["name"] == "Sample-3"
    assert global_results[2]["solved"]
    # This should take 12 iterations to converge with the IDAES solver
    # distribution, but due to various solver factors the number of iterations
    # could vary.
    assert global_results[2]["iters"] == pytest.approx(12, abs=2)

    if os.path.exists(fname):
        os.remove(fname)
    # if os.path.exists(results_fname):
    #     os.remove(results_fname)


@pytest.mark.unit
def test_parse_ipopt_output():
    fname = os.path.join(currdir, "ipopt_output.txt")
    iters, restoration, regularization, time = cb._parse_ipopt_output(fname)

    assert iters == 43
    assert restoration == 39
    assert regularization == 4
    assert time == 0.016 + 0.035


@pytest.mark.unit
def test_compare_convergence_runs_all_same():
    run1 = [
        OrderedDict(
            [
                ("name", "Sample_1"),
                ("solved", True),
                ("iters", 11),
                ("iters_in_restoration", 3),
                ("iters_w_regularization", 7),
                ("time", 21),
            ]
        ),
        OrderedDict(
            [
                ("name", "Sample_2"),
                ("solved", False),
                ("iters", 12),
                ("iters_in_restoration", 4),
                ("iters_w_regularization", 8),
                ("time", 22),
            ]
        ),
        OrderedDict(
            [
                ("name", "Sample_3"),
                ("solved", True),
                ("iters", 13),
                ("iters_in_restoration", 5),
                ("iters_w_regularization", 9),
                ("time", 23),
            ]
        ),
    ]

    baseline = OrderedDict(
        [
            (
                "Sample_1",
                OrderedDict(
                    [
                        ("solved", True),
                        ("iters", 11),
                        ("iters_in_restoration", 3),
                        ("iters_w_regularization", 7),
                        ("time", 21),
                    ]
                ),
            ),
            (
                "Sample_2",
                OrderedDict(
                    [
                        ("solved", False),
                        ("iters", 12),
                        ("iters_in_restoration", 4),
                        ("iters_w_regularization", 8),
                        ("time", 22),
                    ]
                ),
            ),
            (
                "Sample_3",
                OrderedDict(
                    [
                        ("solved", True),
                        ("iters", 13),
                        ("iters_in_restoration", 5),
                        ("iters_w_regularization", 9),
                        ("time", 23),
                    ]
                ),
            ),
        ]
    )

    diff_solves, diff_iters, diff_rest, diff_reg = cb._compare_run_to_baseline(
        run1, baseline
    )

    assert diff_solves == []
    assert diff_iters == []
    assert diff_rest == []
    assert diff_reg == []


@pytest.mark.unit
def test_compare_convergence_runs_invalid_sample_name():
    run1 = [
        OrderedDict(
            [
                ("name", "Sample_X"),
                ("solved", True),
                ("iters", 11),
                ("iters_in_restoration", 3),
                ("iters_w_regularization", 7),
                ("time", 21),
            ]
        ),
        OrderedDict(
            [
                ("name", "Sample_2"),
                ("solved", False),
                ("iters", 12),
                ("iters_in_restoration", 4),
                ("iters_w_regularization", 8),
                ("time", 22),
            ]
        ),
        OrderedDict(
            [
                ("name", "Sample_3"),
                ("solved", True),
                ("iters", 13),
                ("iters_in_restoration", 5),
                ("iters_w_regularization", 9),
                ("time", 23),
            ]
        ),
    ]

    baseline = OrderedDict(
        [
            (
                "Sample_1",
                OrderedDict(
                    [
                        ("solved", True),
                        ("iters", 11),
                        ("iters_in_restoration", 3),
                        ("iters_w_regularization", 7),
                        ("time", 21),
                    ]
                ),
            ),
            (
                "Sample_2",
                OrderedDict(
                    [
                        ("solved", False),
                        ("iters", 12),
                        ("iters_in_restoration", 4),
                        ("iters_w_regularization", 8),
                        ("time", 22),
                    ]
                ),
            ),
            (
                "Sample_3",
                OrderedDict(
                    [
                        ("solved", True),
                        ("iters", 13),
                        ("iters_in_restoration", 5),
                        ("iters_w_regularization", 9),
                        ("time", 23),
                    ]
                ),
            ),
        ]
    )

    with pytest.raises(
        KeyError,
        match="baseline does not contain a sample with the name Sample_X. Please check that "
        "both convergence evaluation runs used the same set of samples.",
    ):
        diff_solves, diff_iters, diff_rest, diff_reg = cb._compare_run_to_baseline(
            run1, baseline
        )


@pytest.mark.unit
def test_compare_convergence_runs_differences_abs():
    run1 = [
        OrderedDict(
            [
                ("name", "Sample_1"),
                ("solved", True),  # diff
                ("iters", 11),  # same
                ("iters_in_restoration", 5),  # diff
                ("iters_w_regularization", 8),  # in tol
                ("time", 21),
            ]
        ),
        OrderedDict(
            [
                ("name", "Sample_2"),
                ("solved", False),  # diff
                ("iters", 13),  # in tol
                ("iters_in_restoration", 3),  # in tol
                ("iters_w_regularization", 6),  # diff
                ("time", 22),
            ]
        ),
        OrderedDict(
            [
                ("name", "Sample_3"),
                ("solved", True),  # same
                ("iters", 15),  # diff
                ("iters_in_restoration", 5),  # same
                ("iters_w_regularization", 9),  # same
                ("time", 23),
            ]
        ),
    ]

    baseline = OrderedDict(
        [
            (
                "Sample_1",
                OrderedDict(
                    [
                        ("solved", False),
                        ("iters", 11),
                        ("iters_in_restoration", 3),
                        ("iters_w_regularization", 7),
                        ("time", 21),
                    ]
                ),
            ),
            (
                "Sample_2",
                OrderedDict(
                    [
                        ("solved", True),
                        ("iters", 12),
                        ("iters_in_restoration", 4),
                        ("iters_w_regularization", 8),
                        ("time", 22),
                    ]
                ),
            ),
            (
                "Sample_3",
                OrderedDict(
                    [
                        ("solved", True),
                        ("iters", 13),
                        ("iters_in_restoration", 5),
                        ("iters_w_regularization", 9),
                        ("time", 23),
                    ]
                ),
            ),
        ]
    )

    diff_solves, diff_iters, diff_rest, diff_reg = cb._compare_run_to_baseline(
        run1, baseline, rel_tol=1e-5, abs_tol=1
    )

    assert diff_solves == ["Sample_1", "Sample_2"]
    assert diff_iters == ["Sample_3"]
    assert diff_rest == ["Sample_1"]
    assert diff_reg == ["Sample_2"]


@pytest.mark.unit
def test_compare_convergence_runs_differences_rel():
    run1 = [
        OrderedDict(
            [
                ("name", "Sample_1"),
                ("solved", True),  # diff
                ("iters", 11),  # same
                ("iters_in_restoration", 5),  # diff
                ("iters_w_regularization", 8),  # in tol
                ("time", 21),
            ]
        ),
        OrderedDict(
            [
                ("name", "Sample_2"),
                ("solved", True),  # diff
                ("iters", 13),  # in tol
                ("iters_in_restoration", 4),  # same
                ("iters_w_regularization", 4),  # diff
                ("time", 22),
            ]
        ),
        OrderedDict(
            [
                ("name", "Sample_3"),
                ("solved", False),  # diff
                ("iters", 20),  # diff
                ("iters_in_restoration", 4),  # in tol
                ("iters_w_regularization", 9),  # same
                ("time", 23),
            ]
        ),
    ]

    baseline = OrderedDict(
        [
            (
                "Sample_1",
                OrderedDict(
                    [
                        ("solved", False),
                        ("iters", 11),
                        ("iters_in_restoration", 3),
                        ("iters_w_regularization", 7),
                        ("time", 21),
                    ]
                ),
            ),
            (
                "Sample_2",
                OrderedDict(
                    [
                        ("solved", True),
                        ("iters", 12),
                        ("iters_in_restoration", 4),
                        ("iters_w_regularization", 8),
                        ("time", 22),
                    ]
                ),
            ),
            (
                "Sample_3",
                OrderedDict(
                    [
                        ("solved", True),
                        ("iters", 13),
                        ("iters_in_restoration", 5),
                        ("iters_w_regularization", 9),
                        ("time", 23),
                    ]
                ),
            ),
        ]
    )

    diff_solves, diff_iters, diff_rest, diff_reg = cb._compare_run_to_baseline(
        run1, baseline, rel_tol=0.25, abs_tol=0
    )

    assert diff_solves == ["Sample_1", "Sample_3"]
    assert diff_iters == ["Sample_3"]
    assert diff_rest == ["Sample_1"]
    assert diff_reg == ["Sample_2"]
