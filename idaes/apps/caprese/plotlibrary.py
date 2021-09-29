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

import pandas as pd
import matplotlib.pyplot as plt
from pyomo.environ import (
    ComponentUID,
    Reference,
    )

__author__ = "Kuan-Han Lin"

        
def plot_plant_state_evolution(varslist,
                               plant_dataframe):
    
    for ind, var in enumerate(varslist):
        str_cuid = str(ComponentUID(var.referent))
        plant_xaxis = plant_dataframe.index
        if str_cuid not in plant_dataframe.columns:
            print("Given state ", var.name, " is not saved in the dataframe. ")
            continue
        plant_state = plant_dataframe[str_cuid]
        
        title = "Simulation result: "+ str_cuid
        plt.figure(title)
        plt.plot(plant_xaxis, plant_state)#, label = "Real state")
        # plt.legend()
        plt.xlabel("time")
        plt.title(title)
    plt.show()
                
def plot_control_input(inputs, 
                       plant_dataframe):
    
    for ind, var in enumerate(inputs):
        str_cuid = str(ComponentUID(var.referent))
        plant_xaxis = plant_dataframe.index
        if str_cuid not in plant_dataframe.columns:
            print("Given input ", var.name, " is not saved in the dataframe. ")
            continue
        input_data = plant_dataframe[str_cuid]
        
        title = "Control input: " + str_cuid
        plt.figure(title)
        plt.plot(plant_xaxis, input_data, "r")
        plt.xlabel("time")
        plt.title(title)
    plt.show()
        
def plot_setpoint_tracking_results(varslist, 
                                   plant_dataframe, 
                                   setpoint_dataframe):        
    
    for ind, var in enumerate(varslist):
        str_cuid = str(ComponentUID(var.referent))
        plant_xaxis = plant_dataframe.index
        if str_cuid not in plant_dataframe.columns:
            print("Given state ", var.name, " is not saved in the dataframe. ")
            continue
        real_state = plant_dataframe[str_cuid]
        setpoint_df_xaxis = setpoint_dataframe.index
        state_setpoint = setpoint_dataframe[str_cuid]
        
        
        title = "Setpoint tracking result: "+ str_cuid
        plt.figure(title)
        plt.plot(plant_xaxis, real_state, label = "Real state")
        plt.plot(setpoint_df_xaxis, state_setpoint, label = "Setpoint")
        plt.legend()
        plt.xlabel("time")
        plt.title(title)
    plt.show()
    
def plot_estimation_results(varslist, 
                            plant_dataframe,
                            estimator_dataframe):
        
    for ind, var in enumerate(varslist):
        str_cuid = str(ComponentUID(var.referent))
        plant_xaxis = plant_dataframe.index
        real_state = plant_dataframe[str_cuid]
        estimator_xaxis = estimator_dataframe.index
        estimate = estimator_dataframe[str_cuid]
        
        title = "Estimation result: "+ str_cuid
        plt.figure(title)
        plt.plot(plant_xaxis, real_state, label = "Real state")
        plt.plot(estimator_xaxis, estimate, "o",label = "estimate")
        plt.legend()
        plt.xlabel("time")
        plt.title(title)
    plt.show()