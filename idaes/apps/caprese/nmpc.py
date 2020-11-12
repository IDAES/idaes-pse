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
        DynamicBlock,
        )
from idaes.apps.caprese.controller import (
        ControllerBlock,
        )
import idaes.logger as idaeslog

__author__ = "Robert Parker and David Thierry"


class NMPCSim(object):
    """
    Main class for NMPC simulations of Pyomo models.
    """
    # pyomo.common.config.add_docstring_list

    def __init__(self, plant_model=None, plant_time_set=None, 
        controller_model=None, controller_time_set=None, inputs_at_t0=None,
        measurements=None, sample_time=None, **kwargs):
        """
        """
        init_plant_measurements = [
                find_comp_in_block(plant_model, controller_model, comp)
                for comp in measurements
                ]
        self.plant = DynamicBlock(
                model=plant_model,
                time=plant_time_set,
                inputs=inputs_at_t0,
                measurements=init_plant_measurements,
                )
        self.plant.construct()

        init_controller_inputs = [
                find_comp_in_block(controller_model, plant_model, comp)
                for comp in inputs_at_t0
                ]
        self.controller = ControllerBlock(
                model=controller_model,
                time=controller_time_set,
                inputs=init_controller_inputs,
                measurements=measurements,
                )
        self.controller.construct()

        controller_measurements = self.controller.measurement_vars
        plant_measurements = self.plant.find_components(
                self.controller.mod,
                controller_measurements,
                self.controller.time,
                )

        plant_inputs = self.plant.category_dict[VariableCategory.INPUT]
        controller_inputs = self.controller.find_components(
                self.plant.mod,
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
    def has_consistent_initial_conditions(self, block, **kwargs):
        """
        Finds constraints at time.first() that are violated by more than
        tolerance. Returns True if any are found.
        """
        # This will raise an error if any constraints at t0 cannot be
        # evaluated, i.e. contain a variable of value None.
        time = block.time
        inconsistent = get_inconsistent_initial_conditions(
                block, 
                time, 
                tol=1e-8,
                suppress_warnings=True)
        return not inconsistent

    def solve_consistent_initial_conditions(self, block, **kwargs):
        """
        Uses pyomo.dae.initialization solve_consistent_initial_conditions
        function to solve for consistent initial conditions. Inputs are
        fixed at time.first() in attempt to eliminate degrees of freedom.
        """
        time = block.time
        strip_bounds = kwargs.pop('strip_bounds', True)
        solver = SolverFactory('ipopt')
        solver.set_options({
            'linear_solver': 'ma57',
            })
        solver_log = idaeslog.getSolveLogger('nmpc', level=idaeslog.DEBUG)
        t0 = time.first()

        previously_fixed = ComponentMap()
        for var in block.vectors.input[:,t0]:
            previously_fixed[var] = var.fixed
            var.fix()

        if strip_bounds:
            strip_var_bounds = TransformationFactory(
                                           'contrib.strip_var_bounds')
            strip_var_bounds.apply_to(block, reversible=True)

        with idaeslog.solver_log(solver_log, level=idaeslog.DEBUG) as slc:
            result = solve_consistent_initial_conditions(block, time, solver,
                    tee=slc.tee)

        if strip_bounds:
            strip_var_bounds.revert(block)

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
