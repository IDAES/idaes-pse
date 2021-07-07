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
Class for performing MHE simulations of IDAES flowsheets
This is just a temporary class. I expect this will merge with nmpc.py,
and then have something like "DynamicSim", which can provide either
(nmpc + plant), or (mhe + plant), or (mhe + nmpc + plant). 
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
from idaes.apps.caprese.estimator import (
        EstimatorBlock,
        )

__author__ = "Kuan-Han Lin"


class MHESim(object):
    """
    This is a user-facing class to perform MHE simulations with Pyomo
    models for both plant and estimator. The user must provide the
    models to use for each, along with sets to treat as "time,"
    inputs in the plant model, and measurements in the estimator
    model. Its functionality is primarily to ensure that these components
    (as defined by the names relative to the corresponding provided models)
    exist on both models.
    """
    # TODO: pyomo.common.config.add_docstring_list

    def __init__(self, plant_model=None, plant_time_set=None, 
        estimator_model=None, estimator_time_set=None, inputs_at_t0=None,
        measurements=None, sample_time=None, **kwargs):
        """
        Measurements must be defined in the estimator model.
        Inputs must be defined in the plant model.
        """
        # To find components in a model given a name,
        # modulo the index of some set:
        # i.   slice the component along the set
        # ii.  create a cuid from that slice
        # iii. get a reference to the slice from the cuid on the new model
        # iv.  access the reference at the index you want (optional)
        self.measurement_cuids = [
                ComponentUID(
                slice_component_along_sets(comp, (estimator_time_set,)))
                for comp in measurements
                ]
        self.input_cuids = [
                ComponentUID(
                slice_component_along_sets(comp, (plant_time_set,)))
                for comp in inputs_at_t0
                ]

        p_t0 = plant_time_set.first()
        init_plant_measurements = [
                cuid.find_component_on(plant_model)[p_t0]
                for cuid in self.measurement_cuids
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
        c_t0 = estimator_time_set.first()
        init_estimator_inputs = [
                cuid.find_component_on(estimator_model)[c_t0]
                for cuid in self.input_cuids
                ]
        self.estimator = EstimatorBlock(
                model=estimator_model,
                time=estimator_time_set,
                inputs=init_estimator_inputs,
                measurements=measurements,
                )
        self.estimator.construct()

        if sample_time is not None:
            self.estimator.set_sample_time(sample_time)
            self.plant.set_sample_time(sample_time)
            self.sample_time = sample_time
            
        #The following should be organized in a function
        self.estimator._add_sample_point_set()
        self.estimator._add_actual_measurement_param()
        self.estimator._add_measurement_error()
        self.estimator._add_measurement_constraint()
        self.estimator._add_model_disturbance()
        self.estimator._add_disturbance_to_differential_cons()
        self.estimator._add_MHE_variables_into_vectors()

