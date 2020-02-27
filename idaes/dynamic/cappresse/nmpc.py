# -*- coding: UTF-8 -*-
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
"""
Class for performing NMPC simulations of IDAES flowsheets
"""

from pyomo.environ import (Block, Constraint, Var, TerminationCondition,
        SolverFactory)
from pyomo.kernel import ComponentSet
from pyomo.dae import ContinuousSet, DerivativeVar
from pyomo.dae.flatten import flatten_dae_variables

from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.dyn_utils import get_activity_dict, deactivate_model_at
import idaes.logger as idaeslog
import pdb

__author__ = "Robert Parker and David Thierry"


class NMPCSim(object):
    """
    Main class for NMPC simulations of IDAES flowsheets. 
    """

    def __init__(self, plant_model, controller_model, steady_model,
            initial_inputs, **kwargs):
        # Maybe include a kwarg for require_steady - if False, set-point is not
        # forced to be a steady state

        solver = kwargs.pop('solver', SolverFactory('ipopt'))

        # Logger properties:
        outlvl = kwargs.pop('outlvl', idaeslog.NOTSET)

        # Set up attributes
        self.p_mod = plant_model
        self.c_mod = controller_model

        self.p_mod._originally_active = get_activity_dict(self.p_mod)
        self.solve_initial_conditions(self.p_mod, initial_inputs,
                activity_dict=self.p_mod._originally_active,
                solver=solver)

        self.categorize_variables(self.p_mod, initial_inputs)

        pdb.set_trace()

        pass

    def validate_models(self, plant_model, controller_model, steady_model):
        if not (isinstance(plant_model, FlowsheetBlock) and 
                isinstance(controller_model, FlowsheetBlock) and
                isinstance(steady_model, FlowsheetBlock)):
            raise ValueError(
                    'Provided models must be FlowsheetBlocks')
        if (plant_model.model() is controller_model.model() or
            plant_model.model() is steady_model.model() or
            steady_model.model() is controller_model.model()):
            raise ValueError(
                    'Provided models must not live in the same top-level'
                    'ConcreteModel')

    def solve_initial_conditions(self, model, initial_inputs, **kwargs):
        # Record which Constraints/Variables are initially inactive
        # deactivate model except at t=0
        # fix initial inputs - raise error if no value
        # raise error if not square
        # re-activate model
        #
        # Later include option to skip solve for consistent initial conditions

        was_originally_active = kwargs.pop('activity_dict', None)
        solver = kwargs.pop('solver', SolverFactory('ipopt'))
        outlvl = kwargs.pop('outlvl', idaeslog.NOTSET)
        solver_log = idaeslog.getSolveLogger(__name__, level=outlvl)
        init_log = idaeslog.getInitLogger(__name__, level=outlvl)

        toplevel = model.model()

        non_initial_time = [t for t in model.time]
        non_initial_time.remove(model.time.first())
        deactivated = deactivate_model_at(model, model.time, non_initial_time, 
                outlvl=idaeslog.ERROR)

        for vardata in initial_inputs:
            if vardata.model() is not toplevel:
                raise ValueError(
                        f"Trying to fix an input that does not belong to model"
                        " {toplevel.name}. Are 'initial_inputs' arguments"
                        " contained in the plant model?")
            if vardata.value == None:
                raise ValueError(
                        "Trying to solve for consistent initial conditions with "
                        "an input of value None. Inputs must have valid values "
                        "to solve for consistent initial conditions")
            else:
                vardata.fix()

        if degrees_of_freedom(model) != 0:
            raise ValueError(
                    f"Non-zero degrees of freedom in model {toplevel.name}"
                    " when solving for consistent initial conditions. Have the "
                    "right number of initial conditions in the plant model been"
                    " fixed?")

        with idaeslog.solver_log(solver_log, level=idaeslog.DEBUG) as slc:
            results = solver.solve(model, tee=slc.tee)

        if results.solver.termination_condition == TerminationCondition.optimal:
            init_log.info(
                    'Successfully solved for consistent initial conditions')
        else:
            init_log.error('Failed to solve for consistent initial conditions')
            raise ValueError(
            f'Falied to solve {toplevel.name} for consistent initial conditions')

        if was_originally_active is not None:
            for t in non_initial_time:
                for comp in deactivated[t]:
                    if was_originally_active[id(comp)]:
                        comp.activate()

    def categorize_variables(self, model, initial_inputs):
        # This is the caveman implementation
        # Would be much faster to do this during flattening
        # Seriously. This is sooo slow
        time = model.time
        t0 = time.first()

        # Shouldn't search for algebraic variables at time.first()
        # because some could be fixed as initial conditions instead
        # of differential variables. 
        # At any other time, however, anything that is not an input/dv/deriv
        # and that is not fixed is algebraic. This assumes that the user
        # has not fixed any algebraic variables (at non-initial time points),
        # which is fine as doing so will give a degree of freedom error in
        # the simulation/initialization.
        t1 = time.get_finite_elements()[1]

        deriv_vars = []
        diff_vars = []
        input_vars = []
        alg_vars = []
        not_vars = []

        scalar_vars, dae_vars = flatten_dae_variables(model, time)

        # Copy before enumerating so I can safely pop without skipping elements
        # ^No. Still changes indices. 
        for i, var in reversed(list(enumerate(dae_vars.copy()))):
            for ini in initial_inputs:
                if var[t0] is ini:
                    input_vars.append(dae_vars.pop(i))
                    break

        #other_copy = dae_vars.copy()
        for i, var in reversed(list(enumerate(dae_vars.copy()))):
            parent = var[t0].parent_component()
            index = var[t0].index()
            if not isinstance(parent, DerivativeVar):
                continue
            if time not in ComponentSet(parent.get_continuousset_list()):
                continue
            deriv_vars.append(dae_vars.pop(i))

            # now need to locate the state variable slice in dae_vars
            # state var:
            #     var[t0].parent_component().get_state_var()[var[t0].index()]

            found = False
            # Don't need to copy here as only expect to pop one entry in this
            # loop
            for j, diffvar in enumerate(dae_vars):
                if not diffvar[t0] is parent.get_state_var()[index]:
                    continue
                found = True
                diff_vars.append(dae_vars.pop(j))
                break

            if not found:
                #raise ValueError(
                print(f'Could not find state for derivative variable {var[t0].name}')

        for i, var in reversed(list(enumerate(dae_vars.copy()))):
            # If the variable is still in the list of time-indexed vars,
            # it must either be fixed (not a var) or be an algebraic var
            if var[t1].fixed:
                not_vars.append(dae_vars.pop(i))
            else:
                alg_vars.append(dae_vars.pop(i))

        model.deriv_vars = deriv_vars
        model.diff_vars = diff_vars
        model.input_vars = input_vars
        model.alg_vars = alg_vars
        model.not_vars = not_vars

