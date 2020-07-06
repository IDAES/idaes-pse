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
"""
This module is a command-line script for executing convergence evaluation
testing on IDAES models.

Convergence evaluation testing is used to verify reliable convergence of a
model over a range of conditions for inputs and parameters. The developer of
the test must create a ConvergenceEvaluation class prior to executing any
convergence testing (see convergence_base.py for documentation).

Convergence evaluation testing is a two step process. In the first step, a json
file is created that contains a set of points sampled from the provided inputs.
This step only needs to be done once - up front. The second step, which should
be executed any time there is a major code change that could impact the model,
takes that set of sampled points and solves the model at each of the points,
collecting convergence statistics (success/failure, iterations, and solution
time).

To find help on convergence.py::

   $ python convergence.py --help

You will see that there are some subcommands. To find help on a particular
subcommand::

   $ python convergence.py  <subcommand> --help

To create a sample file, you can use a command-line like the following (this
should be done once by the  model developer for a few different sample sizes)::

   $ python ../../../core/util/convergence/convergence.py create-sample-file
         -s PressureChanger-10.json
         -N 10 --seed=42
         -e idaes.models.convergence.pressure_changer.
             pressure_changer_conv_eval.PressureChangerConvergenceEvaluation

More commonly, to run the convergence evaluation::

   $ python ../../../core/util/convergence/convergence.py run-eval
         -s PressureChanger-10.json

Note that the convergence evaluation can also be run in parallel if you have
installed MPI and mpi4py using a command line like the following::

   $ mpirun -np 4 python ../../../core/util/convergence/convergence.py run-eval
         -s PressureChanger-10.json

"""
import argparse
import logging
import os
import sys
#
import idaes.core.util.convergence.convergence_base as cb
from idaes.dmf import DMF
from idaes.dmf.errors import DMFError

_log = logging.getLogger(__name__)
_f = logging.Formatter(fmt='%(asctime)s %(name)s [%(levelname)s] %(message)s')
_h = logging.StreamHandler()
_h.setFormatter(_f)
_log.addHandler(_h)


def _parse_arguments():
    # register the command line arguments
    parser = argparse.ArgumentParser(prog='convergence')
    subparsers = parser.add_subparsers()
    subparsers.required = True
    subparsers.help = "Set the subcommand:"
    subparsers.metavar = 'cmd'
    subparsers.dest = 'command'
    run_report_subparser = \
        subparsers.add_parser('run-eval',
                              help="runs the convergence evaluation"
                                   " defined in the sample file given by (-s)")

    create_sample_file_subparser = \
        subparsers.add_parser('create-sample-file',
                              help="create a json file with a specified number"
                                   " of sample points (-N) based on the input"
                                   " parameter distributions defined in the"
                                   " convergence evaluation class (-e)")

    # create validation function for file existence
    def file_exists_validate(parser_obj, arg):
        if not os.path.exists(arg) or not os.path.isfile(arg):
            parser_obj.error('File: {} does not exist.'.format(arg))
        else:
            return arg

    run_report_subparser.add_argument('-s', '--sample-file',
                                      # nargs=1,
                                      default=None,
                                      type=lambda arg: file_exists_validate(
                                              parser, arg),
                                      # choices=
                                      required=True,
                                      help="Path to the json file that "
                                           "contains the points that should be"
                                           " executed in the convergence"
                                           " evaluation.",
                                      metavar='filepath',
                                      dest='sample_file')
    run_report_subparser.add_argument('-D', '--dmf', dest='dmfcfg',
                                      metavar='DIR', default=None,
                                      help='Use DMF configuration at DIR '
                                           '(default=do not use DMF)')
    run_report_subparser.add_argument('-v', '--verbose', dest='vb',
                                      action='count', default=0,
                                      help='Increase output verbosity')

    create_sample_file_subparser.add_argument(
                '-e',
                '--convergence-evaluation-class',
                # nargs=1,
                default=None,
                type=str,
                # choices=
                required=True,
                help="Specify the module.class that identifies the convergence"
                     " evaluation object that should be used to identify the"
                     " input parameters for sampling.",
                metavar='class_str',
                dest='convergence_evaluation_class_str')

    create_sample_file_subparser.add_argument(
                '-s', '--sample-file',
                # nargs=1,
                default=None,
                type=str,
                # choices=
                required=True,
                help="Path to the points file that should be created.",
                metavar='filepath',
                dest='sample_file')

    create_sample_file_subparser.add_argument(
                '-N', '--number-samples',
                # nargs=1,
                default=None,
                type=int,
                # choices=
                required=False,
                help="Used with subcommands: 'create-sample-file'. The number"
                     " of samples that should be generated when creating"
                     " the sample points file.",
                metavar='N',
                dest='number_samples')

    create_sample_file_subparser.add_argument(
                '--seed',
                # nargs=1,
                default=None,
                type=int,
                # choices=
                required=False,
                help="Used with subcommands: 'create-sample-file'. The seed"
                     " value that should be used when generating the random"
                     " samples. If None, a default seeding will be used.",
                metavar='seed',
                dest='seed')

    return parser.parse_args()


def main():
    args = _parse_arguments()

    if args.vb > 1:
        _log.setLevel(logging.DEBUG)
    elif args.vb > 0:
        _log.setLevel(logging.INFO)
    else:
        _log.setLevel(logging.WARN)

    if args.command == 'create-sample-file':
        assert args.number_samples is not None
        assert args.sample_file is not None
        assert args.convergence_evaluation_class_str is not None

        try:
            conv_eval_class = cb._class_import(
                                    args.convergence_evaluation_class_str)
            conv_eval = conv_eval_class()
        except Exception as e:
            print('Failed to find the specified convergence_evaluation_class '
                  'with error: {}'.format(str(e)))
            raise ValueError(
                    'Invalid convergence_evaluation_class specified (-e).')

        spec = conv_eval.get_specification()
        cb.write_sample_file(eval_spec=spec,
                             filename=args.sample_file,
                             convergence_evaluation_class_str=
                                 args.convergence_evaluation_class_str,
                             n_points=args.number_samples,
                             seed=args.seed)
    else:
        if args.dmfcfg is None:
            dmf = None
        else:
            try:
                dmf = DMF(args.dmfcfg)
            except DMFError as err:
                _log.error('Unable to init DMF: {}'.format(err))
                return -1
        (inputs, samples, results) = \
            cb.run_convergence_evaluation_from_sample_file(
                    sample_file=args.sample_file)
        if results is not None:
            cb.save_convergence_statistics(inputs, results, dmf=dmf)
    return 0


if __name__ == '__main__':
    sys.exit(main())
