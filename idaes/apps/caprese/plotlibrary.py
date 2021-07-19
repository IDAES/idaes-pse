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
from pyomo.environ import (ComponentUID,)
from pyomo.util.slices import slice_component_along_sets

__author__ = "Kuan-Han Lin"

def convert_input_in_piecewise_form(xaxis, input_data):
    delta = xaxis[1] - xaxis[0]
    piecewise_ind = []
    for ind in xaxis:
        #drop the first index
        if ind != xaxis[0]:
            piecewise_ind += [ind-delta*0.999, ind]
            
    piecewise_data = []
    for inx, item in enumerate(input_data):
        if inx != 0:
            piecewise_data += [item, item]
            
    return piecewise_ind, piecewise_data

class MHE_PlotLibrary(object):
    def load_data(self):
        self.plant_data = self.plant.dataframe
        self.estimator_data = self.estimator.dataframe
    
    def plot_estimation_result(self, vars_soi = None, xlabel = "time"):
        self.load_data()
        
        if xlabel == "time":
            xaxis = self.plant_data["time"]
        elif xlabel == "iterations":
            xaxis = self.plant_data.index
        else:
            raise RuntimeError("Cannot plot figures, "
                               "xlabel should be either 'time' or 'iterations'")
            
        n_figures = len(vars_soi)
        CUID_list = self.plant.diffvar_cuid
        for figind in range(n_figures):
            curr_comp = vars_soi[figind]
            curr_cuid = ComponentUID(curr_comp)
            if curr_cuid not in CUID_list:
                print(str(curr_cuid), " is not correct.")
                print("Please give the differential variable as 'IndexedComponent_slice', "
                        "e.g. Tall[:, Tj]")
                continue
            
            real_state = self.plant_data[str(curr_cuid)]
            estimated_state = self.estimator_data[str(curr_cuid)]
            plt.figure(figind)
            plt.plot(xaxis, real_state, label = "real state")
            plt.plot(xaxis, estimated_state, label = "estimated state")
            plt.legend()
            plt.xlabel(xlabel)
            plt.title("Estimation result: "+ str(curr_cuid))
        plt.show()
        
class NMPC_PlotLibrary(object):
    def load_data(self):
        self.plant_data = self.plant.dataframe
        self.controller_data = self.controller.dataframe
        
    def plot_control_input(self, vars_soi = None, xlabel = "time"):
        self.load_data()
        
        if xlabel == "time":
            xaxis = self.plant_data["time"]
        elif xlabel == "iterations":
            xaxis = self.plant_data.index
        else:
            raise RuntimeError("Cannot plot figures, "
                               "xlabel should be either 'time' or 'iterations'")
        
        if vars_soi is None:
            print("No specific control input is declared.")
            print("Plot all control inputs.")
            soi_cuids = self.plant.inputvar_cuid
        else:
            #check whether given inputs are inputs and in the correct form
            soi_cuids = []
            for item in vars_soi:
                item_cuid = ComponentUID(item)
                if item_cuid not in self.plant.inputvar_cuid:
                    print(str(item_cuid), " is not correct.")
                    print("Please give the input variable as 'IndexedComponent_slice', "
                            "e.g. Tjinb[:]")
                    continue
                else:
                    soi_cuids.append(item_cuid)
        
        n_figures = len(soi_cuids)
        for figind in range(n_figures):
            curr_cuid = soi_cuids[figind]
            input_data = self.plant_data[str(curr_cuid)]
            piecewise_ind, piecewise_data = convert_input_in_piecewise_form(xaxis, 
                                                                            input_data)
            plt.figure(figind+100)
            plt.plot(piecewise_ind, piecewise_data, "r")
            # plt.ylim([min(piecewise_data)*0.99, max(piecewise_data)]*1.01)
            plt.xlabel(xlabel)
            plt.title("Control input: " + str(curr_cuid))
        plt.show()
        
    def plot_setpoint_tracking_result(self, vars_soi = None, xlabel = "time"):
        self.load_data()
        
        if xlabel == "time":
            xaxis = self.plant_data["time"]
        elif xlabel == "iterations":
            xaxis = self.plant_data.index
        else:
            raise RuntimeError("Cannot plot figures, "
                               "xlabel should be either 'time' or 'iterations'")
            
        n_figures = len(vars_soi)
        CUID_list = self.plant.diffvar_cuid
        for figind in range(n_figures):
            curr_comp = vars_soi[figind]
            curr_cuid = ComponentUID(curr_comp)
            if curr_cuid not in CUID_list:
                print(str(curr_cuid), " is not correct.")
                print("Please give the differential variable as 'IndexedComponent_slice', "
                        "e.g. Tall[:, Tj]")
                continue
            
            real_state = self.plant_data[str(curr_cuid)]
            state_setpoint = self.controller_data[str(curr_cuid) + "_setpoint"]
            if str(state_setpoint[0]) == "nan":
                state_setpoint[0] = state_setpoint[1]
            plt.figure(figind)
            plt.plot(xaxis, real_state, label = "real state")
            plt.plot(xaxis, state_setpoint, label = "setpoint")
            plt.legend()
            plt.xlabel(xlabel)
            plt.title("Setpoint tracking result: "+ str(curr_cuid))
        plt.show()
            
        
        
        
