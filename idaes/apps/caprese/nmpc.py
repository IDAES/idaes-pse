# -*- coding: utf-8 -*-
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

from pyomo.environ import (
        Block, 
        Constraint, 
        Var, 
        TerminationCondition,
        SolverFactory, 
        Objective, 
        NonNegativeReals, 
        Reals, 
        TransformationFactory, 
        Reference, 
        value,
        )
from pyomo.core.base.range import remainder
from pyomo.common.collections import ComponentMap
from pyomo.dae.initialization import (
        solve_consistent_initial_conditions,
        get_inconsistent_initial_conditions,
        )
from pyutilib.misc.config import ConfigDict, ConfigValue

from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import (
        degrees_of_freedom, 
        activated_equalities_generator,
        )
from idaes.core.util.dyn_utils import (
        deactivate_model_at,
        path_from_block, 
        find_comp_in_block, 
        find_comp_in_block_at_time,
        )
from idaes.apps.caprese.common.config import (
        ControlInitOption,
        ElementInitializationInputOption,
        TimeResolutionOption,
        ControlPenaltyType,
        VariableCategory)
from idaes.apps.caprese.model import (
        DynamicModelHelper,
        ControllerHelper,
        )
#from idaes.apps.caprese.util import (
#        initialize_by_element_in_range,
#        find_slices_in_model, 
#        NMPCVarGroup, 
#        NMPCVarLocator, 
#        copy_values_at_time, 
#        validate_list_of_vardata, 
#        validate_list_of_vardata_value_tuples, 
#        validate_solver,
#        find_point_in_continuousset,
#        get_violated_bounds_at_time)
#from idaes.apps.caprese.base_class import DynamicBase
import idaes.logger as idaeslog

__author__ = "Robert Parker and David Thierry"


class NMPCSim(object):
    """
    Main class for NMPC simulations of Pyomo models.
    """
    # pyomo.common.config.add_docstring_list

    def __init__(self, plant_model=None, plant_time_set=None, 
        controller_model=None, controller_time_set=None, inputs_at_t0=None,
        sample_time=None, **kwargs):
        """Constructor method. Accepts plant and controller models needed for 
        NMPC simulation, as well as time sets (Pyomo Sets) in each model
        Inputs at the first time point in the plant model are also required.
        Models provided are added to the NMPCSim instance as attributes.
        This constructor solves for consistent initial conditions 
        in the plant and controller and performs categorization into lists of
        differential, derivative, algebraic, input, fixed, and scalar variables,
        which are added as attributes to a _NMPC_NAMESPACE Block on each model.

        Args:
            plant_model : Plant Pyomo model, NMPC of which will be 
                          simulated. Currently this must contain the entire 
                          timespan it is desired to simulate.
            plant_time_set : Set to treat as time in the plant model
            controller_model : Model to be used to calculate control inputs
                               for the plant. Control inputs in controller
                               must exist in the plant, and initial condition
                               variables in the plant must exist in the 
                               controller.
            controller_time_set : Set to treat as time in the controller model
            inputs_at_t0 : List of VarData objects containing the variables
                             to be treated as control inputs, at time.first().
            solver : Solver to be used for verification of consistent initial 
                     conditions, will also be used as the default solver if
                     another is not provided for initializing or solving the 
                     optimal control problem.
            outlvl : IDAES logger output level. Default is idaes.logger.INFO.
                     To see solver output, use idaes.logger.DEBUG.
            sample_time : Length of time each control input will be held for.
                          This must be an integer multiple of the (finite
                          element) discretization spacing in both the plant
                          and controller models. Default is to use the 
                          controller model's discretization spacing.

        """
        self.plant = DynamicModelHelper(
                plant_model,
                plant_time_set,
                inputs_at_t0,
                )

        init_controller_inputs = [
                find_comp_in_block(controller_model, plant_model, comp)
                for comp in inputs_at_t0
                ]
        self.controller = ControllerHelper(
                controller_model,
                controller_time_set,
                init_controller_inputs,
                )

        controller_measurements = self.controller.measured_vars
        plant_measurements = self.plant.find_components(
                self.controller.model,
                controller_measurements,
                self.controller.time,
                )

        plant_inputs = self.plant.category_dict[VariableCategory.INPUT]
        controller_inputs = self.controller.find_components(
                self.plant.model,
                plant_inputs,
                self.plant.time,
                )

        self.controller_measurement_map = ComponentMap(
                zip(plant_measurements, controller_measurements))

        self.plant_input_map = ComponentMap(
                zip(controller_inputs, plant_inputs))

        self.controller.set_sample_time(sample_time)
        self.plant.set_sample_time(sample_time)
        self.sample_time = sample_time

        self.current_plant_time = 0

    # TODO: put next two methods on a model wrapper.
    def has_consistent_initial_conditions(self, model, **kwargs):
        """
        Finds constraints at time.first() that are violated by more than
        tolerance. Returns True if any are found.
        """
        # This will raise an error if any constraints at t0 cannot be
        # evaluated, i.e. contain a variable of value None.
        namespace = getattr(model, self.get_namespace_name())
        time = namespace.get_time()
        config = self.config(kwargs)
        tolerance = config.tolerance
        inconsistent = get_inconsistent_initial_conditions(
                model, 
                time, 
                tol=tolerance,
                suppress_warnings=True)
        return not inconsistent

    def solve_consistent_initial_conditions(self, model, **kwargs):
        """
        Uses pyomo.dae.initialization solve_consistent_initial_conditions
        function to solve for consistent initial conditions. Inputs are
        fixed at time.first() in attempt to eliminate degrees of freedom.
        """
        namespace = getattr(model, self.get_namespace_name())
        time = namespace.get_time()
        strip_bounds = kwargs.pop('strip_bounds', True)
        config = self.config(kwargs)
        outlvl = config.outlvl
        solver = config.solver
        solver_log = idaeslog.getSolveLogger('nmpc', level=outlvl)
        t0 = time.first()

        previously_fixed = ComponentMap()
        for var in namespace.input_vars:
            var0 = var[t0]
            previously_fixed[var0] = var0.fixed
            var0.fix()

        if strip_bounds:
            strip_var_bounds = TransformationFactory(
                                           'contrib.strip_var_bounds')
            strip_var_bounds.apply_to(model, reversible=True)

        with idaeslog.solver_log(solver_log, level=idaeslog.DEBUG) as slc:
            result = solve_consistent_initial_conditions(model, time, solver,
                    tee=slc.tee)

        if strip_bounds:
            strip_var_bounds.revert(model)

        for var, was_fixed in previously_fixed.items():
            if not was_fixed:
                var.unfix()

        return result

# TODO: This functionality is useful, but I need to rethink how
# I want to do it.
#    def calculate_error_between_states(self, mod1, mod2, t1, t2, 
#            Q_matrix=[],
#            categories=[VariableCategory.DIFFERENTIAL],
#            **kwargs):
#        pass
