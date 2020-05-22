##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
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
import pkgutil
import re
import os
import click
import logging
import idaes.solvers
from idaes.commands import cb

import idaes.core.util.convergence.convergence_base as cnv
from idaes.dmf import DMF
from idaes.dmf.errors import DMFError

_log = logging.getLogger("idaes.commands.convergence")

@cb.command(name="convergence-sample", help="Create a convergence sample file.")
@click.option('-e', '--evaluation-class', default=None, type=str, required=True)
@click.option('-s', '--sample-file', default=None, type=str, required=True)
@click.option('-N', '--number-samples', default=None, type=int, required=True)
@click.option('--seed', default=None, type=int)
def convergence_sample(evaluation_class, sample_file, number_samples, seed):
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
            dmf = DMF(dmf)
        except DMFError as err:
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
