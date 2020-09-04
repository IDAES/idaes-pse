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
This module provides the base classes and methods for running convergence
evaluations on IDAES models. The convergence evaluation runs a given model over
a set of sample points to ensure reliable convergence over the parameter space.

The module requires the user to provide:
  - a set of inputs along with their lower bound, upper bound, mean,
and standard deviation.
  - an initialized Pyomo model
  - a Pyomo solver with appropriate options

The module executes convergence evaluation in two steps. In the first step, a
json file is created that containsa set of points sampled from the provided
inputs. This step only needs to be done once - up front. The second step, which
should be executed any time there is a major code change that could impact the
model, takes that set of sampled points and solves the model at each of the
points, collecting convergence statistics (success/failure, iterations, and
solution time).

This can be used as a tool to evaluate model convergence reliability over the
defined input space, or to verify that convergence performance is not
decreasing with framework and/or model changes.

In order to write a convergence evaluation for your model, you must inherit a
class from ConvergenceEvaluation, and implement three methods:

- get_specification: This method should create and return a
    ConvergenceEvaluationSpecification object. There are methods on
    ConvergenceEvaluationSpecification to add inputs. These inputs contain a
    string that identifies a Pyomo Param or Var object, the lower and upper
    bounds, and the mean and standard deviation to be used for sampling. When
    samples are generated, they are drawn from a normal distribution, and then
    truncated by the lower or upper bounds.
- get_initialized_model: This method should create and return a Pyomo model
    object that is already initialized and ready to be solved. This model will
    be modified according to the sampled inputs, and then it will be solved.
- get_solver: This method should return an instance of the Pyomo solver that
    will be used for the analysis.

There are methods to create the sample points file (on
ConvergenceEvaluationSpecification), to run a convergence evaluation
(run_convergence_evaluation), and print the results in table form
(print_convergence_statistics).

However, this package can also be executed using the command-line interface.
See the documentation in convergence.py for more information.
"""
# stdlib
from collections import OrderedDict
import getpass
import importlib as il
import json
import logging
import numpy as np
import sys
from io import StringIO
# pyomo
import pyutilib.services
from pyutilib.misc import capture_output
from pyomo.core import Param, Var
from pyomo.opt import TerminationCondition
from pyomo.common.log import LoggingIntercept
# idaes
import idaes.core.util.convergence.mpi_utils as mpiu
from idaes.dmf import resource


class ConvergenceEvaluationSpecification(object):
    def __init__(self):
        self._inputs = OrderedDict()

    def add_sampled_input(
        self,
        name,
        pyomo_path,
        lower,
        upper,
        mean=None,
        std=None,
        distribution="normal"
    ):
        """
        Add an input that should be sampled when forming the set of
        specific points that need to be run for the convergence evaluation

        The input will be sampled assuming a normal distribution
        (with given mean and standard devation)
        truncated to the values given by lower and upper bounds

        Parameters
        ----------
        name : str
           The name of the input.
        pyomo_path : str
           A string representation of the path to the variable or parameter to
           be sampled. This string will be executed to retrieve the Pyomo
           component.
        lower : float
           A lower bound on the input variable or parameter.
        upper : float
           An upper bound on the input variable or parameter
        mean : float
           The mean value to use when generating normal distribution samples
        std : float
           The standard deviation to use when generating normal distribution samples
        distribution : str
            The Distribution type {"normal", "uniform"}

        Returns
        -------
           N/A
        """
        # ToDo: put some error checking here ... Maybe we should have the model
        # ToDo: already? Can use to check if the pyomo_path is valid? check if
        # ToDo: bounds will be violated?
        spec = OrderedDict()
        spec['pyomo_path'] = pyomo_path
        spec['lower'] = lower
        spec['upper'] = upper
        spec['mean'] = mean
        spec['std'] = std
        spec['distribution'] = distribution
        self._inputs[name] = spec

    @property
    def inputs(self):
        return self._inputs


class ConvergenceEvaluation(object):
    def __init__(self):
        self._sampling_specifications = list()

    def get_specification(self):
        """
        User should override this method to return an instance of the
        ConvergenceEvaluationSpecification for this particular model and test
        set.

        The basic flow for this method is:
           - Create a ConvergenceEvaluationSpecification
           - Call add_sampled_input for every input that should be varied.
           - return the ConvergenceEvaluationSpecification

        Returns
        -------
           ConvergenceEvaluationSpecification
        """
        raise NotImplementedError('Not implemented in the base class. This'
                                  ' should be overridden in the derived class')

    def get_initialized_model(self):
        """
        User should override this method to return an initialized model that is
        ready to solve. The convergence evaluation methods will change the
        values of parameters or variables according to the sampling
        specifications.

        Returns
        -------
           Pyomo model : return a Pyomo model object that is initialized and
                        ready to solve. This is the model object that will be
                        used in the evaluation.
        """
        raise NotImplementedError('Not implemented in the base class. This'
                                  ' should be overridden in the derived class')

    def get_solver(self):
        """
        User should create and return the solver that will be used for the
        convergence evaluation (including any necessary options)

        Returns
        -------
           Pyomo solver

        """
        # ToDo: We may want this to be standardized across all the models
        # within IDAES (should not need different solver options for different
        # unit models)
        raise NotImplementedError('Not implemented in the base class. This'
                                  ' should be overridden in the derived class')


def _class_import(class_path):
    # Note, this method assumes that everything in front
    # of the last dot is a module, followed by one class
    # I don't think this will work for nested classes
    tokens = class_path.split('.')
    modpath = '.'.join(tokens[0:-1])
    mod = il.import_module(modpath)
    ret_class = getattr(mod, tokens[-1])
    return ret_class


def _run_ipopt_with_stats(model, solver, max_iter=500, max_cpu_time=120):
    """
    Run the solver (must be ipopt) and return the convergence statistics

    Parameters
    ----------
    model : Pyomo model
       The pyomo model to be solved

    solver : Pyomo solver
       The pyomo solver to use - it must be ipopt, but with whichever options
       are preferred

    max_iter : int
       The maximum number of iterations to allow for ipopt

    max_cpu_time : int
       The maximum cpu time to allow for ipopt (in seconds)

    Returns
    -------
       Returns a tuple with (solve status object, bool (solve successful or
       not), number of iters, solve time)
    """
    # ToDo: Check that the "solver" is, in fact, IPOPT

    pyutilib.services.TempfileManager.push()
    tempfile = pyutilib.services.TempfileManager.create_tempfile(
                        suffix='ipopt_out',
                        text=True)
    opts = {'output_file': tempfile,
            'max_iter': max_iter,
            'max_cpu_time': max_cpu_time}

    status_obj = solver.solve(model, options=opts, tee=True)
    solved = True
    if status_obj.solver.termination_condition != TerminationCondition.optimal:
        solved = False

    iters = 0
    time = 0
    # parse the output file to get the iteration count, solver times, etc.
    with open(tempfile, 'r') as f:
        for line in f:
            if line.startswith('Number of Iterations....:'):
                tokens = line.split()
                iters = int(tokens[3])
            elif line.startswith(
                    'Total CPU secs in IPOPT (w/o function evaluations)   ='):
                tokens = line.split()
                time += float(tokens[9])
            elif line.startswith(
                    'Total CPU secs in NLP function evaluations           ='):
                tokens = line.split()
                time += float(tokens[8])

    pyutilib.services.TempfileManager.pop(remove=True)
    return status_obj, solved, iters, time


def _progress_bar(fraction, msg, length=20):
    length = length - 2
    n_complete = int(length*fraction)
    n_remaining = length - n_complete
    characters = ['*']*n_complete         # ['*' for i in range(n_complete)]
    characters.extend(['-']*n_remaining)  # ['-' for i in range(n_remaining)])
    sys.stdout.write('%5.1f%s [%s] %s\n' % (
                            fraction*100.0, '%', ''.join(characters), msg))
    sys.stdout.flush()


def _set_model_parameters_from_sample(model, inputs, sample_point):
    """
    This method takes the parameter values in sample_points and sets them on
    the model to prepare it for the new solve

    Parameters
    ----------
    model : Pyomo model
       The pyomo model that will be modified with values from the new
       sample_points

    inputs : dict
       The inputs dictionary from the sample-file

    sample_point : dict
       The dictionary of var and/or parameter values for the sample point

    Returns
    -------
       N/A
    """
    for k, v in sample_point.items():
        if k == '_name':
            # because parallel task manager does not handle dictionaries,
            # we add an _name entry to the sample point dictionary
            # ToDo: remove this when parallel task manager is fixed
            continue
        # k stores the "name" of the input
        # need to get the pyomo path of the input from the inputs structure
        pyomo_path = inputs[k]['pyomo_path']

        comp = model.find_component(pyomo_path)
        try:
            ctype = comp.ctype
        except AttributeError:
            ctype = None

        if ctype is Param:
            if comp.is_constant():
                raise ValueError(
                        'Convergence testing found an input of type Param that'
                        'was not mutable. Please make sure all sampled inputs'
                        ' are either mutable params or fixed vars.')
            comp.set_value(v)
            # print('Just set', comp,'to', float(v))

        elif ctype is Var:
            if not comp.is_fixed():
                raise ValueError(
                        'Convergence testing found an input of type Var that'
                        'was not fixed. Please make sure all sampled inputs'
                        ' are either mutable params or fixed vars.')
            comp.set_value(float(v))
            # print('Just set', comp,'to', float(v))
        else:
            raise ValueError('Failed to find a valid input component (must be'
                             ' a fixed Var or a mutable Param). Instead,'
                             ' pyomo_path: {} returned: {}'
                             .format(pyomo_path, comp))


def write_sample_file(eval_spec,
                      filename,
                      convergence_evaluation_class_str,
                      n_points,
                      seed=None):
    """
    Samples the space of the inputs defined in the eval_spec, and creates a
    json file with all the points to be used in executing a convergence
    evaluation

    Parameters
    ----------
    filename : str
       The filename for the json file that will be created containing all the
       points to be run
    eval_spec : ConvergenceEvaluationSpecification
       The convergence evaluation specification object that we would like to
       sample
    convergence_evaluation_class_str : str
       Python string that identifies the convergence evaluation class for this
       specific evaluation. This is usually in the form of module.class_name.
    n_points : int
       The total number of points that should be created
    seed : int or None
       The seed to be used when generating samples. If set to None, then the
       seed is not set
    Returns
    -------
       N/A
    """
    if seed is not None:
        np.random.seed(seed)

    # build the samples
    samples = OrderedDict()
    for i in range(n_points):
        sample = samples['Sample-{}'.format(i+1)] = OrderedDict()
        for k, v in eval_spec.inputs.items():
            if v["distribution"] == "normal":
                s = np.random.normal(loc=v['mean'], scale=v['std'])
                s = v["lower"] if s < v["lower"] else s
                s = v["upper"] if s > v["upper"] else s
            elif v["distribution"] == "uniform":
                s = np.random.uniform(low=v['lower'], high=v['upper'])
            sample[k] = s

    # create the dictionary storing all the necessary information
    jsondict = OrderedDict()
    jsondict['inputs'] = OrderedDict(eval_spec.inputs)
    jsondict['n_points'] = len(samples)
    jsondict['seed'] = seed
    jsondict['samples'] = samples
    jsondict['convergence_evaluation_class_str'] = \
        convergence_evaluation_class_str

    with open(filename, 'w') as fd:
        json.dump(jsondict, fd, indent=3)


def run_convergence_evaluation_from_sample_file(sample_file):
    # load the sample file
    try:
        with open(sample_file, 'r') as fd:
            jsondict = json.load(fd, object_pairs_hook=OrderedDict)
    except Exception as e:
        print('Error reading json file: {}'.format(str(e)))
        raise ValueError('Problem reading json file: {}'.format(sample_file))

    # create the convergence evaluation object
    convergence_evaluation_class_str = jsondict[
            'convergence_evaluation_class_str']
    try:
        convergence_evaluation_class = _class_import(
                                            convergence_evaluation_class_str)
        conv_eval = convergence_evaluation_class()
    except Exception as e:
        print('Error creating the class: {}. Error returned: {}'.format(
                    convergence_evaluation_class_str,
                    str(e)))
        raise ValueError(
                'Invalid value specified for convergence_evaluation_class_str:'
                '{} in sample file: {}'.format(
                        convergence_evaluation_class_str, sample_file))

    return run_convergence_evaluation(jsondict, conv_eval)


def run_convergence_evaluation(sample_file_dict, conv_eval):
    """
    Run convergence evaluation and generate the statistics based on information
    in the sample_file.

    Parameters
    ----------
    sample_file_dict : dict
        Dictionary created by ConvergenceEvaluationSpecification that contains
        the input and sample point information

    conv_eval : ConvergenceEvaluation
        The ConvergenceEvaluation object that should be used

    Returns
    -------
       N/A
    """
    inputs = sample_file_dict['inputs']
    samples = sample_file_dict['samples']

    # current parallel task manager code does not work with dictionaries, so
    # convert samples to a list
    # ToDo: fix and test parallel task manager with dictionaries and change
    # this
    samples_list = list()
    for k, v in samples.items():
        v['_name'] = k
        samples_list.append(v)
    n_samples = len(samples_list)

    task_mgr = mpiu.ParallelTaskManager(n_samples)
    local_samples_list = task_mgr.global_to_local_data(samples_list)

    results = list()
    for (si, ss) in enumerate(local_samples_list):
        sample_name = ss['_name']
        # print progress on the rank-0 process
        if task_mgr.is_root():
            _progress_bar(float(si) / float(len(local_samples_list)),
                          'Root Process: {}'.format(sample_name))

        # capture the output
        # ToDo: make this an option and turn off for single sample execution
        output_buffer = StringIO()
        with LoggingIntercept(output_buffer, 'idaes', logging.ERROR):
            with capture_output():  # as str_out:
                model = conv_eval.get_initialized_model()
                _set_model_parameters_from_sample(model, inputs, ss)
                solver = conv_eval.get_solver()
                (status_obj, solved, iters, time) = \
                    _run_ipopt_with_stats(model, solver)

        # run without output capture
        # model = conv_eval.get_initialized_model()
        # _set_model_parameters_from_sample(model, inputs, ss)
        # solver = conv_eval.get_solver()
        # (status_obj, solved, iters, time) = \
        #        _run_ipopt_with_stats(model, solver)

        if not solved:
            print('Sample: {} failed to converge.'.format(sample_name))

        results_dict = OrderedDict()
        results_dict['name'] = sample_name
        results_dict['sample_point'] = ss
        results_dict['solved'] = solved
        results_dict['iters'] = iters
        results_dict['time'] = time
        results.append(results_dict)

    global_results = task_mgr.gather_global_data(results)
    return inputs, samples, global_results


def save_convergence_statistics(inputs, results, dmf=None):
    s = Stats(results)
    if dmf is None:
        print_convergence_statistics(inputs, results, s)
    else:
        save_results_to_dmf(dmf, inputs, results, s)


class Stats(object):
    def __init__(self, results):
        self.notable_cases, self.failed_cases = [], []
        self.iters_successful = list()
        self.time_successful = list()

        # loop through and gather some data
        for r in results:
            if r['solved'] is True:
                self.iters_successful.append(r['iters'])
                self.time_successful.append(r['time'])
            else:
                self.failed_cases.append(r)
        # data for summary table
        self.iters_min = 0
        self.iters_mean = 0
        self.iters_std = 0
        self.iters_max = 0
        if len(self.iters_successful) > 0:
            self.iters_min = int(np.min(self.iters_successful))
            self.iters_mean = int(np.mean(self.iters_successful))
            self.iters_std = int(np.std(self.iters_successful))
            self.iters_max = int(np.max(self.iters_successful))
        self.time_min = 0
        self.time_mean = 0
        self.time_std = 0
        self.time_max = 0
        if len(self.time_successful) > 0:
            self.time_min = float(np.min(self.time_successful))
            self.time_mean = float(np.mean(self.time_successful))
            self.time_std = float(np.std(self.time_successful))
            self.time_max = float(np.max(self.time_successful))


def print_convergence_statistics(inputs, results, s):
    """
    Print the statistics returned from run_convergence_evaluation in a set of
    tables

    Parameters
    ----------
    inputs : dict
       The inputs dictionary returned by run_convergence_evaluation
    results : dict
       The results dictionary returned by run_convergence_evaluation

    Returns
    -------
       N/A
    """

    print()
    print('==== Scenario Statistics ====')
    print('%20s %10s %10s %10s %10s' % ('Parameter', 'Min', 'Mean',
                                        'StdDev', 'Max'))
    print('-'*64)
    for k, v in inputs.items():
        values = [results[i]['sample_point'][k] for i in range(len(results))]
        print('%20s %10g %10g %10g %10g' % (k, np.min(values),
                                            float(np.mean(values)),
                                            float(np.std(values)),
                                            float(np.max(values))))
    print('-'*64)

    print()
    print('==== Summary ====')
    print('Number of Successful Cases (solved=True): %.0f/%.0f' % (
            len(results) - len(s.failed_cases), len(results)))
    print('... Iterations      (min, -1std, mean, +1std, max):'
          '%5.0f, %5f, %5f, %5f, %5.0f'
          % (s.iters_min, s.iters_mean - s.iters_std,
             s.iters_mean, s.iters_mean + s.iters_std,
             s.iters_max))

    print('... Solver Time (s) (min, -1std, mean, +1std, max):'
          '%5f, %5f, %5f, %5f, %5f'
          % (s.time_min, s.time_mean - s.time_std,
             s.time_mean, s.time_mean + s.time_std,
             s.time_max))

    # print the detailed table
    print()
    print('==== Table of Results ====')
    print()
    print('%4s %20s %10s %10s %10s' % ('Flag', 'Name', 'Solved',
                                       'Iters', 'Time'))
    print('-' * 58)
    for r in results:
        flag = ''
        if r['solved'] is not True:
            flag = 'F'
        else:
            if r['time'] > s.time_mean + 2.0 * s.time_std and \
                    r['time'] > s.time_mean + 5.0:
                # add a more absolute check when s.time_std is small
                flag += 'T'
            if r['iters'] > s.iters_mean + 2.0 * s.iters_std and \
                    r['iters'] > s.iters_mean + 5:
                # add a more absolute check when s.iters_std is small
                flag += 'I'

        if flag != '':
            r['flag'] = flag
            s.notable_cases.append(r)

        print('%4s %20s %10s %10.0f %10.2f' % (
                flag, r['name'], r['solved'], r['iters'], r['time']))
        if r['solved'] is not True:
            s.failed_cases.append(r)
        else:
            s.iters_successful.append(r['iters'])
            s.time_successful.append(r['time'])

    print('-' * 58)

    print()
    print('==== Notable Cases ====')
    if len(s.notable_cases) == 0:
        print('... None')
    for c in s.notable_cases:
        msg = ''
        if c['flag'] == 'F':
            msg = 'failed solve'
        elif c['flag'] == 'T':
            msg = 'long runtime'
        elif c['flag'] == 'I':
            msg = 'high iteration count'
        elif c['flag'] == 'TI':
            msg = 'high iteration count; long runtime'

        print(c['name'], ':', msg)
        for k, v in c['sample_point'].items():
            if k == '_name':
                # we need this because the parallel task manager does not
                # handle dicts  and we add _name to identify the sample name
                # ToDo: fix this when the parallel task manager is extended to
                # handle dicts
                continue
            if inputs[k]['distribution'] == "normal":
                mean = inputs[k]['mean']
                std = inputs[k]['std']
                alpha = (v - mean) / std
                print('       %20s: %10g (%5.2f standard deviations'
                      ' above/below the mean)' % (k, v, alpha))


def save_results_to_dmf(dmf, inputs, results, stats):
    """Save results of run, along with stats, to DMF.

    Args:
        dmf (DMF): Data management framework object
        inputs (dict): Run inputs
        results (dict): Run results
        stats (Stats): Calculated result statistics

    Returns:
        None
    """
    d = {}  # data
    x = {}
    for k, v in inputs.items():
        values = [results[i]['sample_point'][k] for i in range(len(results))]
        x[k] = {
            'min': float(np.min(values)),
            'mean': float(np.mean(values)),
            'stdev': float(np.std(values)),
            'max': float(np.max(values))
        }
    d['scenarios'] = x
    d['summary'] = {
        'cases_success': len(results) - len(stats.failed_cases),
        'cases_total': len(results),
        'iters': {
            'min': stats.iters_min,
            '-1std': stats.iters_mean - stats.iters_std,
            'mean': stats.iters_mean,
            '+1std': stats.iters_mean + stats.iters_std,
            'max': stats.iters_max
        },
        'time': {
            'min': stats.time_min,
            '-1std': stats.time_mean - stats.time_std,
            'mean': stats.time_mean,
            '+1std': stats.time_mean + stats.time_std,
            'max': stats.time_max
        }
    }
    tbl = []
    for r in results:
        notable = False
        solved = r['solved']
        item = {
            'solved': solved,
            'name': r['name'],
            'iters': r['iters'],
            'time': r['time']
        }
        for vtype in 'time', 'iters':
            mean, std = [getattr(stats, a) for a in ('{}_mean'.format(vtype),
                                                     '{}_std'.format(vtype))]
            if r[vtype] > mean + 2. * std and r[vtype] > mean + 5:
                item['outlier_{}'.format(vtype)] = True
                notable = True
        if not solved:
            stats.failed_cases.append(r)
            notable = True
        else:
            stats.iters_successful.append(r['iters'])
            stats.time_successful.append(r['time'])
        tbl.append(item)
        if notable:
            stats.notable_cases.append(r)
    d['results_table'] = tbl
    # notable cases
    cases = []
    for r in stats.notable_cases:
        msg, is_outlier = [], False
        if not r['solved']:
            msg.append('failed solve')
        else:
            is_outlier = True
            if r['outlier_time']:
                msg.append('long runtime')
            if r['outlier_iters']:
                msg.append('high iteration count')
        case = {'messages': msg}
        if is_outlier:
            outliers = {}
            for k, v in r['sample_point'].items():
                if k == '_name':
                    # we need this because the parallel task manager does not
                    # handle dicts and we add _name to identify the sample name
                    # TODO: Fix this when the parallel task manager is extended
                    # TODO: to handle dicts
                    continue
                mean = inputs[k]['mean']
                std = inputs[k]['std']
                alpha = (v - mean) / std
                outliers[k] = (v, alpha)
            case['outliers'] = outliers
        else:
            case['outliers'] = []
        cases.append(case)
    d['notable_cases'] = cases
    # Build and save resource object
    rsrc = resource.Resource(value={
        'name': 'convergence_results',
        'desc': 'statistics returned from run_convergence_evaluation',
        'creator': {'name': getpass.getuser()},
        'data': d}, type_=resource.ResourceTypes.data)
    dmf.add(rsrc)
