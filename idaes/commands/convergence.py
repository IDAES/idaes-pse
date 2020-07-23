##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
"""Commandline interface for convergence testing tools"""

__author__ = "John Eslick"

import inspect
import re
import os
import click
import logging
from pyomo.common.dependencies import attempt_import
from idaes.commands import cb
import idaes

cnv = attempt_import('idaes.core.util.convergence.convergence_base')[0]
dmf = attempt_import('idaes.dmf')[0]
dmf_error = attempt_import('idaes.dmf.errors')[0]
pkgutil = attempt_import('pkgutil')[0]


_log = logging.getLogger("idaes.commands.convergence")

@cb.command(name="convergence-sample", help="Create a convergence sample file.")
@click.option('-e', '--evaluation-class', default=None, type=str, required=True,
    help="Convergence evaluation class")
@click.option('-s', '--sample-file', default=None, type=str, required=True,
    help="Output sample file")
@click.option('-N', '--number-samples', default=None, type=int, required=True,
    help="Number of samples")
@click.option('--seed', default=None, type=int)
def convergence_sample(
    evaluation_class, sample_file, number_samples, seed):
    try:
        conv_eval_class = cnv._class_import(evaluation_class)
        conv_eval = conv_eval_class()
    except Exception as e:
        click.echo('Failed to find the specified convergence_evaluation_class '
              'with error: {}'.format(str(e)))
        raise ValueError('Invalid convergence_evaluation_class specified (-e).')

    spec = conv_eval.get_specification()
    cnv.write_sample_file(
        eval_spec=spec,
        filename=sample_file,
        convergence_evaluation_class_str=evaluation_class,
        n_points=number_samples,
        seed=seed,
    )

@cb.command(name="convergence-eval", help="Run convergence sample evaluation.")
@click.option('-s', '--sample-file', default=None, type=str, required=True)
@click.option('-D', '--dmf', default=None, type=str)
def convergence_eval(sample_file, dmf):
    if dmf is not None:
        try:
            dmf = dmf.DMF(dmf)
        except dmf_error.DMFError as err:
            _log.error('Unable to init DMF: {}'.format(err))
            return -1
    (inputs, samples, results) = cnv.run_convergence_evaluation_from_sample_file(
        sample_file=sample_file
    )
    if results is not None:
        cnv.save_convergence_statistics(inputs, results, dmf=dmf)


@cb.command(name="convergence-search", help="Search for convergence test classes.")
@click.option('-r', '--regex', default=".+ConvergenceEvaluation$", type=str)
def convergence_search(regex):
    pat = re.compile(regex)
    for loader, module_name, is_pkg in pkgutil.walk_packages(idaes.__path__):
        try:
            m = loader.find_module(module_name).load_module(module_name)
            c = inspect.getmembers(m, inspect.isclass)
        except:
            continue
        if c:
            for i in c:
                cname = ".".join(["idaes", m.__name__, i[0]])
                if pat.match(i[0]):
                    click.echo(cname)
