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
import numpy as np
from pyomo.environ import (ComponentUID,)
from idaes.apps.caprese.common.config import (
        VariableCategory as VC,
        ConstraintCategory as CC,
        )
from pyomo.util.slices import slice_component_along_sets

__author__ = "Kuan-Han Lin"

def get_cuid_from_nmpcvar_list(var_list, cs):
    cuid_list = []
    str_list = []
    for var in var_list:
        t0 = cs.first()
        sliced_set = slice_component_along_sets(var[cs.first()], (cs,))
        cuid = ComponentUID(sliced_set)
        cuid_list.append(cuid)
        str_list.append(str(cuid))
    return cuid_list, str_list


class DynamicDataManger(object):
    def __init__(self):
        diffvar_cuid_list, diffvar_str_list = get_cuid_from_nmpcvar_list(self.category_dict[VC.DIFFERENTIAL],
                                                                         self.time)
        inputvar_cuid_list, inputvar_str_list = get_cuid_from_nmpcvar_list(self.category_dict[VC.INPUT],
                                                                           self.time)
        meavar_cuid_list, meavar_str_list = get_cuid_from_nmpcvar_list(self.measurement_vars,
                                                                       self.time)
        self.diffvar_cuid = diffvar_cuid_list
        self.diffvar_str = diffvar_str_list
        self.inputvar_cuid = inputvar_cuid_list
        self.inputvar_str = inputvar_str_list
        self.meavar_cuid = meavar_cuid_list
        self.meavar_str = meavar_str_list
        
    def save_data_in_dataframe(self, dynamic_data):
        rowind = len(self.dataframe.index)
        curr_time = rowind * self.sample_time
        curr_data = [curr_time] + dynamic_data
        self.dataframe.loc[rowind] = curr_data
        
        
class PlantDataManger(DynamicDataManger):
    def __init__(self):
        print("Save plant states and control inputs in plant's attribute 'dataframe'")
        
    def construct_plant_dataframe(self):
        super(PlantDataManger, self).__init__()

        columns = ["time"] + self.diffvar_str + self.inputvar_str
        self.dataframe = pd.DataFrame(columns = columns)
        
    def initialize_plant_dataframe(self, skip_t0 = False):
        self.construct_plant_dataframe()
        
        if not skip_t0:
            t0 = self.time.first()
            diffvar_list = [self.DIFFERENTIAL_BLOCK[bind].var[t0].value 
                            for bind in self.DIFFERENTIAL_SET]
            inputvar_list = [np.nan for var in range(len(self.INPUT_SET))]
            row0 = [t0] + diffvar_list + inputvar_list
            
            rowind = len(self.dataframe.index)
            self.dataframe.loc[rowind] = row0
        
    def record_plant_data(self):
        t_sp = self.sample_points[1]
        inputvar_list = [self.INPUT_BLOCK[bind].var[t_sp].value
                         for bind in self.INPUT_SET]
        diffvar_list = [self.DIFFERENTIAL_BLOCK[bind].var[t_sp].value 
                        for bind in self.DIFFERENTIAL_SET]
        plant_data = diffvar_list + inputvar_list
        self.save_data_in_dataframe(plant_data)
    
    
class ControllerDataManger(DynamicDataManger):
    def __init__(self):
        print("Save initial states, setpoints, and control inputs in controller's attribute 'dataframe'")
        
    def construct_controller_dataframe(self):
        super(ControllerDataManger, self).__init__()
        
        diffvar_ic_columns = [item+"_ic" for item in self.diffvar_str]
        diffvar_setpoint_columns = [item+"_setpoint" for item in self.diffvar_str]
        
        columns = ["time"] + diffvar_ic_columns +  diffvar_setpoint_columns + self.inputvar_str
        self.dataframe = pd.DataFrame(columns = columns)
        
    def initialize_controller_dataframe(self, skip_t0 = False):
        self.construct_controller_dataframe()
        
        if not skip_t0:
            row0 = [np.nan for ind in range(len(self.dataframe.columns))]
            rowind = len(self.dataframe.index)
            self.dataframe.loc[rowind] = row0
        
    def record_controller_data(self):
        t0 = self.time.first()
        t_sp = self.sample_points[1]
        
        diffvar_ic_list = [self.DIFFERENTIAL_BLOCK[bind].var[t0].value
                           for bind in self.DIFFERENTIAL_SET]
        diffvar_sp_list = [self.DIFFERENTIAL_BLOCK[bind].var.setpoint
                           for bind in self.DIFFERENTIAL_SET]
        inputvar_list = [self.INPUT_BLOCK[bind].var[t_sp].value
                         for bind in self.INPUT_SET]
        nmpc_data = diffvar_ic_list + diffvar_sp_list + inputvar_list
        self.save_data_in_dataframe(nmpc_data)
        
    
class EstimatorDataManger(DynamicDataManger):
    def __init__(self):
        print("Save estimated states and actual measurements in estimator's attribute 'dataframe'")
        
    def construct_estimator_dataframe(self):
        super(EstimatorDataManger, self).__init__()

        actualmea_columns = [item+"_actualmeas" for item in self.meavar_str]
        
        columns = ["time"] + self.diffvar_str + actualmea_columns
        self.dataframe = pd.DataFrame(columns = columns)
        
    def initialize_estimator_dataframe(self, skip_t0 = False):
        self.construct_estimator_dataframe()
        
        if not skip_t0:
            row0 = [np.nan for ind in range(len(self.dataframe.columns))]
            rowind = len(self.dataframe.index)
            self.dataframe.loc[rowind] = row0
        
    def record_estimator_data(self):
        tlast = self.time.last()
        
        diffvar_tlast_list = [self.DIFFERENTIAL_BLOCK[bind].var[tlast].value
                             for bind in self.DIFFERENTIAL_SET]
        actualmeavar_tlast_list = [self.ACTUALMEASUREMENT_BLOCK[bind].var[tlast].value
                                  for bind in self.ACTUALMEASUREMENT_SET]
        mhe_data = diffvar_tlast_list + actualmeavar_tlast_list
        self.save_data_in_dataframe(mhe_data)
        
        