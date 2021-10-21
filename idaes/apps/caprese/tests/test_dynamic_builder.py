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
Test for dynamic builder.
"""


import pytest
import pyomo.environ as pyo
from idaes.apps.caprese.dynamic_builder import (
            capture_vars_t0_w_CUID,
            DynamicSim,)
from idaes.apps.caprese.tests.test_simple_model import (
            make_model,)

__author__ = "Kuan-Han Lin"


solver_available = pyo.SolverFactory('ipopt').available()
if solver_available:
    solver = pyo.SolverFactory('ipopt')
else:
    solver = None


class TestDynamicBuilder(object):
    
    @pytest.mark.unit
    def test_capture_vars_t0_w_CUID(self):
        mod1 = make_model(horizon = 0.5, nfe = 2)
        mod2 = make_model(horizon = 1., nfe = 2)
        
        mod1t0 = mod1.time.first()
        mod1_varlist = [mod1.conc[:, "A"], mod1.conc[:, "B"]]
        mod1_cuids = [pyo.ComponentUID(var) for var in mod1_varlist]
        
        mod2t0 = mod2.time.first()
        mod2_varlist = capture_vars_t0_w_CUID(mod1_cuids, mod2t0, mod2)

        # Make sure return variables are in mod2
        assert all(var.model() is mod2 for var in mod2_varlist)
        pred_str_cuids = ["conc[0,A]", "conc[0,B]"]
        assert all(str(pyo.ComponentUID(var)) == pred_str for var, pred_str in zip(mod2_varlist, pred_str_cuids))
        
        
    @pytest.mark.unit
    def test_DynamicSim_case1(self):
        '''
        Case 1: plant only.

        '''
        sample_time = 0.5
        p_mod = make_model(horizon = sample_time, nfe = 2)
        p_time = p_mod.time
        t0 = p_time.first()
        inputs = [p_mod.flow_in[t0]]
        measurements = [p_mod.conc[t0,'A']]                

        dyna = DynamicSim(plant_model = p_mod,
                          plant_time_set = p_time,
                          inputs_at_t0 = inputs,
                          measurements_at_t0 = measurements,
                          sample_time = sample_time)
        
        assert type(dyna) is DynamicSim
        assert dyna.has_plant
        assert not dyna.has_controller
        assert not dyna.has_estimator
        assert hasattr(dyna, "plant")
        assert hasattr(dyna.plant, "sample_points")
        assert dyna.plant.sample_points == [0, 0.5]
        
        
    @pytest.mark.unit
    def test_DynamicSim_case2(self):
        '''
        Case 2: plant and controller exist.
        
        '''
        try:
            sample_time = 0.5
            p_mod = make_model(horizon = sample_time, nfe = 2)
            p_time = p_mod.time
            t0 = p_time.first()
            inputs = [p_mod.flow_in[t0]]
            measurements = [p_mod.conc[t0,'A']]                
            
            c_mod = make_model(horizon = 2, nfe = 8)
            c_time = c_mod.time
        
            dyna_fail = DynamicSim(plant_model = p_mod,
                                   plant_time_set = p_time,
                                   controller_model = c_mod,
                                   controller_time_set = c_time,
                                   inputs_at_t0 = inputs,
                                   measurements_at_t0 = measurements,
                                   sample_time = sample_time)
        except RuntimeError:
            sample_time = 0.5
            p_mod = make_model(horizon = sample_time, nfe = 2)
            p_time = p_mod.time
            t0 = p_time.first()
            inputs = [p_mod.flow_in[t0]]
            measurements = [p_mod.conc[t0,'A'], p_mod.conc[t0, "B"],]
            
            c_mod = make_model(horizon = 2, nfe = 8)
            c_time = c_mod.time
            
            dyna = DynamicSim(plant_model = p_mod,
                              plant_time_set = p_time,
                              controller_model = c_mod,
                              controller_time_set = c_time,
                              inputs_at_t0 = inputs,
                              measurements_at_t0 = measurements,
                              sample_time = sample_time)
            
            
        assert type(dyna) is DynamicSim
        assert dyna.has_plant
        assert dyna.has_controller
        assert not dyna.has_estimator
        assert not dyna.controller.has_estimator
        assert hasattr(dyna, "plant")
        assert hasattr(dyna, "controller")
        assert hasattr(dyna.plant, "sample_points")
        assert dyna.plant.sample_points == [0, 0.5]
        assert hasattr(dyna.controller, "sample_points")
        assert dyna.controller.sample_points == [0.5*i for i in range(0, 5)]
        

    @pytest.mark.unit
    def test_DynamicSim_case3(self):
        '''
        Case 3: plant and estimator exist.
        
        '''
        sample_time = 0.5
        p_mod = make_model(horizon = sample_time, nfe = 2)
        p_time = p_mod.time
        t0 = p_time.first()
        inputs = [p_mod.flow_in[t0]]
        measurements = [p_mod.conc[t0,'A']]                
        
        e_mod = make_model(horizon = 2.5, nfe = 5)
        e_time = e_mod.time
    
        dyna = DynamicSim(plant_model = p_mod,
                          plant_time_set = p_time,
                          estimator_model = e_mod,
                          estimator_time_set = e_time,
                          inputs_at_t0 = inputs,
                          measurements_at_t0 = measurements,
                          sample_time = sample_time)
    
        assert type(dyna) is DynamicSim
        assert dyna.has_plant
        assert not dyna.has_controller
        assert dyna.has_estimator
        assert hasattr(dyna, "plant")
        assert hasattr(dyna, "estimator")
        assert hasattr(dyna.plant, "sample_points")
        assert dyna.plant.sample_points == [0, 0.5]
        assert hasattr(dyna.estimator, "sample_points")
        assert dyna.estimator.sample_points == [0.5*i for i in range(0, 6)]

                
    @pytest.mark.unit
    def test_DynamicSim_case4(self):
        '''
        Case 4: plant, controller, and estimator exist.
        
        '''
        sample_time = 0.5
        p_mod = make_model(horizon = sample_time, nfe = 2)
        p_time = p_mod.time
        t0 = p_time.first()
        inputs = [p_mod.flow_in[t0]]
        measurements = [p_mod.flow_out[t0]]                

        c_mod = make_model(horizon = 2, nfe = 8)
        c_time = c_mod.time    
        
        e_mod = make_model(horizon = 2.5, nfe = 5)
        e_time = e_mod.time
    
        dyna = DynamicSim(plant_model = p_mod,
                          plant_time_set = p_time,
                          controller_model = c_mod,
                          controller_time_set = c_time,
                          estimator_model = e_mod,
                          estimator_time_set = e_time,
                          inputs_at_t0 = inputs,
                          measurements_at_t0 = measurements,
                          sample_time = sample_time)
    
        assert type(dyna) is DynamicSim
        assert dyna.has_plant
        assert dyna.has_controller
        assert dyna.controller.has_estimator
        assert dyna.has_estimator
        assert hasattr(dyna, "plant")
        assert hasattr(dyna, "controller")
        assert hasattr(dyna, "estimator")
        assert hasattr(dyna.plant, "sample_points")
        assert dyna.plant.sample_points == [0, 0.5]
        assert hasattr(dyna.estimator, "sample_points")
        assert dyna.estimator.sample_points == [0.5*i for i in range(0, 6)]
        assert hasattr(dyna.controller, "sample_points")
        assert dyna.controller.sample_points == [0.5*i for i in range(0, 5)]
        
        # check fixed initial conditions in the controller are differential 
        # varisables instead of measurements.
        c_t0 = c_time.first()
        assert all(not vart0.fixed 
                   for vart0 in dyna.controller.vectors.measurement[:, c_t0])
        assert all(dyna.controller.vectors.differential[:, c_t0].fixed)
