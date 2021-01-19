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

import random
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
        ComponentUID,
        )
from pyomo.core.base.range import remainder
from pyomo.common.collections import ComponentMap
from pyomo.dae.initialization import (
        solve_consistent_initial_conditions,
        get_inconsistent_initial_conditions,
        )
from pyomo.util.slices import slice_component_along_sets
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
    # TODO: pyomo.common.config.add_docstring_list

    def __init__(self, plant_model=None, plant_time_set=None, 
        controller_model=None, controller_time_set=None, inputs_at_t0=None,
        measurements=None, sample_time=None, **kwargs):
        """
        Measurements must be defined in the controller model.
        Inputs must be defined in the plant model.
        """
        # To find components in a model given a name,
        # modulo the index of some set:
        # i.   slice the component along the set
        # ii.  create a cuid from that slice
        # iii. get a (any) component in the new model from the cuid
        # iv.  slice the new component along the corresponding set
        # v.   create a reference to that slice
        # vi.  access the reference at the index you want (optional)
        self.measurement_cuids = [
                ComponentUID(
                slice_component_along_sets(comp, (controller_time_set,)))
                for comp in measurements
                ]
        self.input_cuids = [
                ComponentUID(
                slice_component_along_sets(comp, (plant_time_set,)))
                for comp in inputs_at_t0
                ]

        init_plant_measurements = []
        for cuid in self.measurement_cuids:
            # Here I perform steps iii. and iv. of the above.
            # With the CUID rewrite, I should be able to combine these steps
            # and get the slice directly from the CUID.
            for comp in cuid.list_components(plant_model):
                break
            _slice = slice_component_along_sets(comp, (plant_time_set,))
            ref = Reference(_slice)
            t0 = plant_time_set.first()
            init_plant_measurements.append(ref[t0])

        self.plant = DynamicBlock(
                model=plant_model,
                time=plant_time_set,
                inputs=inputs_at_t0,
                measurements=init_plant_measurements,
                )
        self.plant.construct()

        # Here we repeat essentially the same "find component"
        # procedure as above.
        init_controller_inputs = []
        for cuid in self.input_cuids:
            for comp in cuid.list_components(controller_model):
                break
            _slice = slice_component_along_sets(comp, (controller_time_set,))
            ref = Reference(_slice)
            t0 = controller_time_set.first()
            init_controller_inputs.append(ref[t0])
        self.controller = ControllerBlock(
                model=controller_model,
                time=controller_time_set,
                inputs=init_controller_inputs,
                measurements=measurements,
                )
        self.controller.construct()

        if sample_time is not None:
            self.controller.set_sample_time(sample_time)
            self.plant.set_sample_time(sample_time)
            self.sample_time = sample_time

        t0 = controller_time_set.first()
        self.noise_bounds = [(0.0, var[t0].ub) for var in 
                self.controller.measurement_vars]
        self.noise_function = random.gauss
