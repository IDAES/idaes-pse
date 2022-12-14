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
"""Commandline interface for convergence testing tools"""

__author__ = "John Eslick"

import importlib
import re
import click
import json
import logging
from pyomo.common.dependencies import attempt_import
from idaes.commands import cb

cnv = attempt_import("idaes.core.util.convergence.convergence_base")[0]
dmf = attempt_import("idaes.core.dmf")[0]
dmf_error = attempt_import("idaes.core.dmf.errors")[0]

_log = logging.getLogger("idaes.commands.convergence")


@cb.command(name="convergence-sample", help="Create a convergence sample file.")
@click.option(
    "-e",
    "--evaluation-class",
    default=None,
    type=str,
    required=True,
    help="Convergence evaluation class",
)
@click.option(
    "-s",
    "--sample-file",
    default=None,
    type=str,
    required=True,
    help="Output sample file",
)
@click.option(
    "-N",
    "--number-samples",
    default=None,
    type=int,
    required=True,
    help="Number of samples",
)
@click.option("--seed", default=None, type=int, help="Random number generator seed")
@click.option(
    "-m",
    "--convergence_module",
    default=None,
    type=str,
    required=False,
    help="Addtional module that registers ConvergenceEvaluation classes",
)
def convergence_sample(
    evaluation_class, sample_file, number_samples, seed, convergence_module
):
    import idaes.models.convergence
    import idaes.models_extra.convergence

    if convergence_module is not None:
        mod = importlib.import_module(convergence_module)
    if evaluation_class in cnv.convergence_classes:
        evaluation_class = cnv.convergence_classes[evaluation_class]
    try:
        conv_eval_class = cnv._class_import(evaluation_class)
        conv_eval = conv_eval_class()
    except Exception as e:
        click.echo(
            "Failed to find the specified convergence_evaluation_class "
            "with error: {}".format(str(e))
        )
        raise ValueError("Invalid convergence_evaluation_class specified (-e).")

    spec = conv_eval.get_specification()
    cnv.write_sample_file(
        eval_spec=spec,
        filename=sample_file,
        convergence_evaluation_class_str=evaluation_class,
        n_points=number_samples,
        seed=seed,
    )


@cb.command(name="convergence-eval", help="Run convergence sample evaluation.")
@click.option(
    "-s",
    "--sample-file",
    default=None,
    type=str,
    required=True,
    help="Path of sample file to run",
)
@click.option(
    "-D",
    "--dmf",
    default=None,
    type=str,
    required=False,
    help="Path to Data Managment Framwork (DMF) workspace",
)
@click.option(
    "-r",
    "--report-file",
    default=None,
    type=str,
    required=False,
    help="Optional text report file path",
)
@click.option(
    "-j",
    "--json-file",
    default=None,
    type=str,
    required=False,
    help="Optional json file to save results and stats",
)
@click.option(
    "-m",
    "--convergence_module",
    default=None,
    type=str,
    required=False,
    help="Addtional module that registers ConvergenceEvaluation classes",
)
@click.option(
    "--single-sample",
    default=None,
    type=str,
    help="Run only a single sample with given name",
)
def convergence_eval(
    sample_file, dmf, report_file, json_file, convergence_module, single_sample
):
    import idaes.models.convergence
    import idaes.models_extra.convergence

    if convergence_module is not None:
        mod = importlib.import_module(convergence_module)
    if dmf is not None:
        try:
            dmf = dmf.DMF(dmf)
        except dmf_error.DMFError as err:
            _log.error("Unable to init DMF: {}".format(err))
            return -1
    if single_sample is None:
        (inputs, samples, results) = cnv.run_convergence_evaluation_from_sample_file(
            sample_file=sample_file
        )
        if results is not None:
            cnv.save_convergence_statistics(
                inputs, results, dmf=dmf, report_path=report_file, json_path=json_file
            )
    else:
        results = cnv.run_single_sample_from_sample_file(
            sample_file=sample_file, name=single_sample
        )
        click.echo(
            json.dumps(
                {"solved": results[1], "iters": results[2], "time": results[3]},
                indent=4,
            )
        )


@cb.command(name="convergence-search", help="Search for convergence test classes.")
@click.option(
    "-r",
    "--regex",
    default=None,
    type=str,
    help="Optional regular expression to filter registered convergence classes",
)
@click.option(
    "-m",
    "--convergence_module",
    default=None,
    type=str,
    required=False,
    help="Optional additional module that registers convergence classes",
)
def convergence_search(regex, convergence_module):
    import idaes.models.convergence
    import idaes.models_extra.convergence

    if convergence_module is not None:
        mod = importlib.import_module(convergence_module)
    if regex is not None:
        pat = re.compile(regex)
    else:
        pat = None
    l = []
    for k, v in cnv.convergence_classes.items():
        if pat is None or pat.match(k) or pat.match(v):
            l.append(k)
    for k in sorted(l):
        click.echo(f"{k}:\n   {cnv.convergence_classes[k]}")
