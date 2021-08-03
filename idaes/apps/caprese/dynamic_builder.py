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
Class for performing plant, MHE, NMPC simulations of IDAES flowsheets.
"""

import random
from pyomo.environ import (
        ComponentUID,
        )
from pyomo.util.slices import slice_component_along_sets
from pyomo.common.config import ConfigDict, ConfigValue

from idaes.apps.caprese.dynamic_block import DynamicBlock
from idaes.apps.caprese.controller import ControllerBlock
from idaes.apps.caprese.estimator import EstimatorBlock
from idaes.apps.caprese.common.config import VariableCategory as VC

__author__ = "Robert Parker, David Thierry, and Kuan-Han Lin"


def use_CUID_to_capture_vars_t0_in_given_model(cuid_list, 
                                               target_t0,
                                               target_mod):
    varlist_in_target_model = [cuid.find_component_on(target_mod)[target_t0]
                               for cuid in cuid_list]
            
    return varlist_in_target_model
            
class DynamicSim(object):
    def __init__(self,
                 plant_model = None,
                 plant_time_set = None,
                 estimator_model = None,
                 estimator_time_set = None,
                 controller_model = None, 
                 controller_time_set = None,
                 inputs_at_t0 = None,
                 measurements_at_t0 = None,
                 sample_time = None,):
        
        self.input_cuids = [
                            ComponentUID(
                            slice_component_along_sets(comp, (plant_time_set,)))
                            for comp in inputs_at_t0
                            ]
        self.measurement_cuids = [
                                    ComponentUID(
                                    slice_component_along_sets(comp, (plant_time_set,)))
                                    for comp in measurements_at_t0
                                    ]
        
        self.plant_is_existing = False
        self.estimator_is_existing = False
        self.controller_is_existing = False
        
        #----------------------------------------------------------------------                
        if plant_model:
            self.plant_is_existing = True
            
            # Capture vars in plant model
            plant_inputs_t0 = inputs_at_t0
            plant_measurements_t0 = measurements_at_t0
            
            self.plant = DynamicBlock(
                    model=plant_model,
                    time=plant_time_set,
                    inputs=plant_inputs_t0,
                    measurements=plant_measurements_t0,
                    )
            self.plant.construct()
            
            if sample_time is not None:
                self.plant.set_sample_time(sample_time)
                
                
                
        #----------------------------------------------------------------------                
        if estimator_model:
            self.estimator_is_existing = True
            
            # Capture vars in estimator model
            estimator_inputs_t0 = use_CUID_to_capture_vars_t0_in_given_model(self.input_cuids,
                                                                              estimator_time_set.first(),
                                                                              estimator_model,)
            estimator_measurements_t0 = use_CUID_to_capture_vars_t0_in_given_model(self.measurement_cuids,
                                                                                    estimator_time_set.first(),
                                                                                    estimator_model,)
            
            self.estimator = EstimatorBlock(
                    model=estimator_model,
                    time=estimator_time_set,
                    inputs=estimator_inputs_t0,
                    measurements=estimator_measurements_t0,
                    sample_time = sample_time, #MHE requires the sample time.
                    )
            self.estimator.construct()
    
            if sample_time is not None:
                self.estimator.set_sample_time(sample_time)
        
        
        #----------------------------------------------------------------------
        if controller_model:
            self.controller_is_existing = True
            
            # Capture vars in controller model
            controller_inputs_t0 = use_CUID_to_capture_vars_t0_in_given_model(self.input_cuids,
                                                                              controller_time_set.first(),
                                                                              controller_model,)
            controller_measurements_t0 = use_CUID_to_capture_vars_t0_in_given_model(self.measurement_cuids,
                                                                                    controller_time_set.first(),
                                                                                    controller_model,)
            
            self.controller = ControllerBlock(
                    model=controller_model,
                    time=controller_time_set,
                    inputs=controller_inputs_t0,
                    measurements=controller_measurements_t0,
                    )
            self.controller.construct()        
            
            if sample_time is not None:
                self.controller.set_sample_time(sample_time)
             
            # Controller should know whether the estimator exists or not    
            if self.estimator_is_existing:
                self.controller.estimator_is_existing = True
            else:
                self.controller.estimator_is_existing = False
    
    
        if sample_time is not None:
            self.sample_time = sample_time
            

        #check number of measurements and differential vars if there is no mhe
        if self.controller_is_existing and not self.estimator_is_existing:
            if len(measurements_at_t0) != \
                len(self.controller.differential_vars):
            
                raise RuntimeError(
                    "The controller is declared but there is no estimator. \n"
                    "Therefore, the number of declared measurements should equal to "
                    "the number of differential variables.")
                
        if self.controller_is_existing and self.estimator_is_existing:
            # In this case, controller should take the estimation results from MHE, 
            # instead of the measurements directly from the plant for initial conditions.       
            
            # Initial conditions are defined by the differential variables, 
            # not by measurements.
            controller_t0 = self.controller.time.first()
            self.controller.vectors.measurement[:, controller_t0].unfix()
            self.controller.vectors.differential[:, controller_t0].fix()
