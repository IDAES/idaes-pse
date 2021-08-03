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


class PLANT_PlotLibrary(object):
    def get_plant_dataframe(self):
        raise NotImplementedError("Must implement a method to get plant data")
        
    def plot_plant_state_evolution(self,
                                   varslist = None,
                                   plant_dataframe = None,):
        
        if plant_dataframe is None:
            plant_dataframe = self.get_plant_dataframe()
            
        if varslist is None:
            raise RuntimeError("No specific state is declared. "
                               "Please declare some states for plotting.")

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
                    
    def plot_control_input(self, 
                           inputs = None, 
                           dataframe=None):
        if dataframe is None:
            dataframe = self.get_plant_dataframe()
            
        if inputs is None:
            raise RuntimeError("No specific input is declared."
                               "Please declare some control inputs for plotting.")
        
        for ind, var in enumerate(inputs):
            str_cuid = str(ComponentUID(var.referent))
            input_data = dataframe[str_cuid]
            
            title = "Control input: " + str_cuid
            plt.figure(title)
            plt.plot(dataframe.index, input_data, "r")
            plt.xlabel("time")
            plt.title(title)
        plt.show()


class NMPC_PlotLibrary(PLANT_PlotLibrary):    
    def get_controller_dataframe(self):
        raise NotImplementedError("Must implement a method to get controller data")
        
    def plot_setpoint_tracking_results(self, 
                                       varslist = None, 
                                       plant_dataframe = None, 
                                       controller_dataframe = None):        
        if plant_dataframe is None:
            plant_dataframe = self.get_plant_dataframe()
        
        if controller_dataframe is None:
            controller_dataframe = self.get_controller_dataframe()
        
        if varslist is None:
            raise RuntimeError("No specific state is declared. "
                               "Please declare some states for plotting.")
        
        for ind, var in enumerate(varslist):
            str_cuid = str(ComponentUID(var.referent))
            plant_xaxis = plant_dataframe.index
            if str_cuid not in plant_dataframe.columns:
                print("Given state ", var.name, " is not saved in the dataframe. ")
                continue
            real_state = plant_dataframe[str_cuid]
            controller_xaxis = controller_dataframe.index
            state_setpoint = controller_dataframe[str_cuid + "_setpoint"]
            
            
            title = "Setpoint tracking result: "+ str_cuid
            plt.figure(title)
            plt.plot(plant_xaxis, real_state, label = "Real state")
            plt.plot(controller_xaxis, state_setpoint, label = "Setpoint")
            plt.legend()
            plt.xlabel("time")
            plt.title(title)
        plt.show()
        

class MHE_PlotLibrary(PLANT_PlotLibrary):
    def get_estimator_dataframe(self):
        raise NotImplementedError("Must implement a method to get estimator data")
    
    def plot_estimation_results(self, 
                                varslist = None, 
                                plant_dataframe = None,
                                estimator_dataframe = None,):
        
        if plant_dataframe is None:
            plant_dataframe = self.get_plant_dataframe()
            
        if estimator_dataframe is None:
            estimator_dataframe = self.get_estimator_dataframe()
            
        if varslist is None:
            raise RuntimeError("No specific state is declared."
                               "Please declare some states for plotting.")
            
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