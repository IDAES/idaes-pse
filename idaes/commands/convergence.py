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
"""Commandline interface for convergence testing tools"""
# This code is deprecated
# pylint: disable=missing-function-docstring

__author__ = "John Eslick"

import click
from pyomo.common.dependencies import attempt_import
from idaes.commands import cb

cnv = attempt_import("idaes.core.util.convergence.convergence_base")[0]


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
    help="Additional module that registers ConvergenceEvaluation classes",
)
def convergence_sample(
    evaluation_class, sample_file, number_samples, seed, convergence_module
):
    click.echo(
        "The command line interface for convergence testing has been deprecated and no longer works. "
        "A new API has been added to run convergence tests directly from the ConvergenceEvaluation object."
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
    help="Path to Data Management Framework (DMF) workspace",
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
    help="Additional module that registers ConvergenceEvaluation classes",
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
    click.echo(
        "The command line interface for convergence testing has been deprecated and no longer works. "
        "A new API has been added to run convergence tests directly from the ConvergenceEvaluation object."
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
    click.echo(
        "The command line interface for convergence testing has been deprecated and no longer works. "
        "A new API has been added to run convergence tests directly from the ConvergenceEvaluation object."
    )
