# -*- coding: utf-8 -*-
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
"""
Class for performing NMPC simulations of IDAES flowsheets
"""

import random
from pyomo.environ import (
    ComponentUID,
)
from pyomo.util.slices import slice_component_along_sets
from pyomo.common.config import ConfigDict, ConfigValue

from idaes.apps.caprese.dynamic_block import (
    DynamicBlock,
)
from idaes.apps.caprese.controller import (
    ControllerBlock,
)

__author__ = "Robert Parker and David Thierry"


class NMPCSim(object):
    """
    This is a user-facing class to perform NMPC simulations with Pyomo
    models for both plant and controller. The user must provide the
    models to use for each, along with sets to treat as "time,"
    inputs in the plant model, and measurements in the controller
    model. Its functionality is primarily to ensure that these components
    (as defined by the names relative to the corresponding provided models)
    exist on both models.
    """

    # TODO: pyomo.common.config.add_docstring_list

    def __init__(
        self,
        plant_model=None,
        plant_time_set=None,
        controller_model=None,
        controller_time_set=None,
        inputs_at_t0=None,
        measurements=None,
        sample_time=None,
        **kwargs
    ):
        """
        Measurements must be defined in the controller model.
        Inputs must be defined in the plant model.
        """
        # To find components in a model given a name,
        # modulo the index of some set:
        # i.   slice the component along the set
        # ii.  create a cuid from that slice
        # iii. get a reference to the slice from the cuid on the new model
        # iv.  access the reference at the index you want (optional)
        self.measurement_cuids = [
            ComponentUID(slice_component_along_sets(comp, (controller_time_set,)))
            for comp in measurements
        ]
        self.input_cuids = [
            ComponentUID(slice_component_along_sets(comp, (plant_time_set,)))
            for comp in inputs_at_t0
        ]

        p_t0 = plant_time_set.first()
        init_plant_measurements = [
            cuid.find_component_on(plant_model)[p_t0] for cuid in self.measurement_cuids
        ]

        self.plant = DynamicBlock(
            model=plant_model,
            time=plant_time_set,
            inputs=inputs_at_t0,
            measurements=init_plant_measurements,
        )
        self.plant.construct()

        # Here we repeat essentially the same "find component"
        # procedure as above.
        c_t0 = controller_time_set.first()
        init_controller_inputs = [
            cuid.find_component_on(controller_model)[c_t0] for cuid in self.input_cuids
        ]
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
