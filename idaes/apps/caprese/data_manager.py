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

# import pandas as pd
# import numpy as np
# from pyomo.environ import (ComponentUID,)
# from idaes.apps.caprese.common.config import (
#         VariableCategory as VC,
#         ConstraintCategory as CC,
#         )
# from pyomo.util.slices import slice_component_along_sets

import pyomo.environ as pyo

from pyomo.core.base.componentuid import ComponentUID
from pyomo.common.collections import ComponentMap
import pandas as pd
from collections import OrderedDict
from idaes.apps.caprese.plotlibrary import (
                        NMPC_PlotLibrary,
                        MHE_PlotLibrary,)

__author__ = "Robert Parker and Kuan-Han Lin"


def empty_dataframe_from_variables(variables, rename_map=None):
    if rename_map is None:
        # Sometimes we construct a variable with a programatically
        # generated name that is not that meaningful to th user, e.g.
        # ACTUALMEASUREMENTBLOCK[1].var. This argument lets us rename
        # to something more meaningful, such as temperature_actual.
        rename_map = ComponentMap()
    var_dict = OrderedDict()
    hash_set = set()
    var_dict["iteration"] = []
    for var in variables:
        # TODO: In this case, do we need to make sure we have the
        # same indexing sets?
        #
        #idx_set = var.index_set()
        #set_hash = hash(tuple(idx_set))
        #if len(hash_set) == 0:
        #    hash_set.add(set_hash)
        #elif set_hash not in hash_set:
        #    raise ValueError()
        if var in rename_map:
            key = rename_map[var]
        else:
            if var.is_reference():
                # NOTE: This will fail if var is a reference-to-list
                # or reference-to-dict.
                key = str(ComponentUID(var.referent))
            else:
                key = str(ComponentUID(var))
        var_dict[key] = []
        dataframe = pd.DataFrame(var_dict)
        # Change the type of iteration column to integer
        dataframe["iteration"] = dataframe["iteration"].astype(int)
    return dataframe


def add_variable_values_to_dataframe(
        dataframe,
        variables,
        iteration,
        time_subset=None,
        rename_map=None,
        initial_time=None,
        time_map=None,
        ):
    if rename_map is None:
        rename_map = ComponentMap()
    
    if initial_time is None:
        if len(dataframe.index) == 0:
            initial_time = 0.0
        else:
            initial_time = dataframe.index[-1]

    df_map = OrderedDict()
    hash_set = set()
    for var in variables:
        idx_set = var.index_set()
        set_hash = hash(tuple(idx_set))
        if len(hash_set) == 0:
            hash_set.add(set_hash)

            if time_subset is None:
                time_subset = idx_set
                
            df_map["iteration"] = len(time_subset)*[iteration]

        elif set_hash not in hash_set:
            raise ValueError()

        if var in rename_map:
            key = rename_map[var]
        else:
            if var.is_reference():
                # NOTE: This will fail if var is a reference-to-list
                # or reference-to-dict.
                key = str(ComponentUID(var.referent))                
            else:
                key = str(ComponentUID(var))

        # append only the desired values
        df_map[key] = [var[t].value for t in time_subset] # Values to append for each variable

    if time_map:
        time_points = time_map.values()
    else:
        time_points = [initial_time + t for t in time_subset]
    df = pd.DataFrame(df_map, index=time_points)
    
    return dataframe.append(df)

def add_setpoint_column_into_dataframe(dataframe, variables):
    for var in variables:
        column_name = str(ComponentUID(var.referent)) + "_setpoint"
        dataframe[column_name] = []
    return dataframe

def add_variable_setpoints_to_dataframe(dataframe, variables, time_subset, map_for_user_given_vars = None):
    for var in variables:
        column_name = str(ComponentUID(var.referent)) + "_setpoint"
        column_ind = dataframe.columns.get_loc(column_name)
        start_row_ind = len(dataframe.index) - len(time_subset)
        for pt in range(len(time_subset)):
            row_ind = start_row_ind + pt
            if var in map_for_user_given_vars:
                dataframe.iat[row_ind, column_ind] = map_for_user_given_vars[var].setpoint
            else:
                dataframe.iat[row_ind, column_ind] = var.setpoint
    return dataframe
        

class PlantDataManager(object):
    def __init__(self, 
                 plantblock,
                 states_of_interest = [],
                 inputs_of_interest = []):
                
        # vardatamap = self.plantblock.vardata_map
        # t0 = self.plantblock.time.first()
        if states_of_interest != []:
            cuid_states = [ComponentUID(var.referent) for var in states_of_interest]
            self.plant_states_of_interest = [cuid.find_component_on(plantblock)
                                             for cuid in cuid_states]
        else:
            self.plant_states_of_interest = []
            
        if inputs_of_interest != []:
            cuid_inputs = [ComponentUID(var.referent) for var in inputs_of_interest]
            self.plant_inputs_of_interest = [cuid.find_component_on(plantblock)
                                             for cuid in cuid_inputs]
        else:
            self.plant_inputs_of_interest = []
        
        self.plantblock = plantblock
        self.plant_vars_of_interest = self.plant_states_of_interest + \
                                        self.plant_inputs_of_interest + \
                                            self.plantblock.differential_vars + \
                                                self.plantblock.input_vars

        self.plant_df = empty_dataframe_from_variables(self.plant_vars_of_interest)
        
    def get_plant_dataframe(self):
        return self.plant_df
        
    def save_initial_plant_data(self):
        plant_states_to_save = self.plant_states_of_interest + self.plantblock.differential_vars
        self.plant_df = add_variable_values_to_dataframe(self.plant_df,
                                                         plant_states_to_save, #no inputs
                                                         iteration = 0,
                                                         time_subset = [self.plantblock.time.first()])
        

    def save_plant_data(self, iteration):
        #skip time.first()
        time_subset = self.plantblock.time.ordered_data()[1:]
        self.plant_df = add_variable_values_to_dataframe(self.plant_df,
                                                         self.plant_vars_of_interest,
                                                         iteration,
                                                         time_subset = time_subset,)
        
class ControllerDataManager(PlantDataManager, NMPC_PlotLibrary):
    def __init__(self, 
                 plantblock, 
                 controllerblock, 
                 states_of_interest = [],
                 inputs_of_interest = []):
        super(ControllerDataManager, self).__init__(plantblock, 
                                                   states_of_interest, 
                                                   inputs_of_interest,)
        
        self.controllerblock = controllerblock
        # Convert vars in plant to vars in controller
        if states_of_interest != []:
            cuid_states = [ComponentUID(var.referent) for var in states_of_interest]
            controller_states_of_interest = [cuid.find_component_on(controllerblock)
                                             for cuid in cuid_states]
        else:
            controller_states_of_interest = []
            
        if inputs_of_interest != []:
            cuid_inputs = [ComponentUID(var.referent) for var in inputs_of_interest]
            controller_inputs_of_interest = [cuid.find_component_on(controllerblock)
                                             for cuid in cuid_inputs]
        else:
            controller_inputs_of_interest = []
        
        
        self.controller_vars_of_interest = controller_inputs_of_interest + controllerblock.input_vars   
        self.controller_df = empty_dataframe_from_variables(self.controller_vars_of_interest)
        
        self.states_need_setpoints = controller_states_of_interest + controllerblock.differential_vars
        self.controller_df = add_setpoint_column_into_dataframe(self.controller_df, 
                                                                self.states_need_setpoints)
        vardata_map = self.controllerblock.vardata_map
        t0 = self.controllerblock.time.first()
        self.user_given_vars_map_nmpcvar = ComponentMap((var, vardata_map[var[t0]]) for var in controller_states_of_interest)

  
    def get_controller_dataframe(self):
        return self.controller_df
    
    def save_controller_data(self, iteration):
        #skip time.first()
        time = self.controllerblock.time
        #Save values in the first sample time
        time_subset = [t for t in time if t <= self.controllerblock.sample_points[1] 
                                           and t != time.first()]
        self.controller_df = add_variable_values_to_dataframe(self.controller_df,
                                                              self.controller_vars_of_interest,
                                                              iteration,
                                                              time_subset = time_subset,)

        self.controller_df = add_variable_setpoints_to_dataframe(self.controller_df,
                                                                 self.states_need_setpoints,
                                                                 time_subset,
                                                                 self.user_given_vars_map_nmpcvar,)
        
class EstimatorDataManager(PlantDataManager, MHE_PlotLibrary):
    def __init__(self, 
                 plantblock, 
                 estimatorblock, 
                 states_of_interest = [],):
        super(EstimatorDataManager, self).__init__(plantblock, 
                                                   states_of_interest,)
        
        self.estimatorblock = estimatorblock
        # Convert vars in plant to vars in controller
        if states_of_interest != []:
            cuid_states = [ComponentUID(var.referent) for var in states_of_interest]
            estimator_states_of_interest = [cuid.find_component_on(estimatorblock)
                                            for cuid in cuid_states]
            
        self.estimator_vars_of_interest = estimator_states_of_interest + estimatorblock.differential_vars
        self.estimator_df = empty_dataframe_from_variables(self.estimator_vars_of_interest)
        
    def get_estimator_dataframe(self):
        return self.estimator_df
    
    def save_estimator_data(self, iteration):
        time = self.estimatorblock.time
        t_last = time.last()
        time_map = OrderedDict([(t_last, (iteration+1)*self.estimatorblock.sample_time),])
        self.estimator_df = add_variable_values_to_dataframe(self.estimator_df,
                                                             self.estimator_vars_of_interest,
                                                             iteration,
                                                             time_subset = [t_last],
                                                             time_map = time_map,)
        
        
#-----------------------------------------------------------------------------------------------------
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
        
        
# class PlantDataManger(DynamicDataManger):
#     def __init__(self):
#         print("Save plant states and control inputs in plant's attribute 'dataframe'")
        
#     def construct_plant_dataframe(self):
#         super(PlantDataManger, self).__init__()

#         columns = ["time"] + self.diffvar_str + self.inputvar_str
#         self.dataframe = pd.DataFrame(columns = columns)
        
#     def initialize_plant_dataframe(self, skip_t0 = False):
#         self.construct_plant_dataframe()
        
#         if not skip_t0:
#             t0 = self.time.first()
#             diffvar_list = [self.DIFFERENTIAL_BLOCK[bind].var[t0].value 
#                             for bind in self.DIFFERENTIAL_SET]
#             inputvar_list = [np.nan for var in range(len(self.INPUT_SET))]
#             row0 = [t0] + diffvar_list + inputvar_list
            
#             rowind = len(self.dataframe.index)
#             self.dataframe.loc[rowind] = row0
        
#     def record_plant_data(self):
#         t_sp = self.sample_points[1]
#         inputvar_list = [self.INPUT_BLOCK[bind].var[t_sp].value
#                          for bind in self.INPUT_SET]
#         diffvar_list = [self.DIFFERENTIAL_BLOCK[bind].var[t_sp].value 
#                         for bind in self.DIFFERENTIAL_SET]
#         plant_data = diffvar_list + inputvar_list
#         self.save_data_in_dataframe(plant_data)
    
    
# class ControllerDataManger(DynamicDataManger):
#     def __init__(self):
#         print("Save initial states, setpoints, and control inputs in controller's attribute 'dataframe'")
        
#     def construct_controller_dataframe(self):
#         super(ControllerDataManger, self).__init__()
        
#         diffvar_ic_columns = [item+"_ic" for item in self.diffvar_str]
#         diffvar_setpoint_columns = [item+"_setpoint" for item in self.diffvar_str]
        
#         columns = ["time"] + diffvar_ic_columns +  diffvar_setpoint_columns + self.inputvar_str
#         self.dataframe = pd.DataFrame(columns = columns)
        
#     def initialize_controller_dataframe(self, skip_t0 = False):
#         self.construct_controller_dataframe()
        
#         if not skip_t0:
#             row0 = [np.nan for ind in range(len(self.dataframe.columns))]
#             rowind = len(self.dataframe.index)
#             self.dataframe.loc[rowind] = row0
        
#     def record_controller_data(self):
#         t0 = self.time.first()
#         t_sp = self.sample_points[1]
        
#         diffvar_ic_list = [self.DIFFERENTIAL_BLOCK[bind].var[t0].value
#                            for bind in self.DIFFERENTIAL_SET]
#         diffvar_sp_list = [self.DIFFERENTIAL_BLOCK[bind].var.setpoint
#                            for bind in self.DIFFERENTIAL_SET]
#         inputvar_list = [self.INPUT_BLOCK[bind].var[t_sp].value
#                          for bind in self.INPUT_SET]
#         nmpc_data = diffvar_ic_list + diffvar_sp_list + inputvar_list
#         self.save_data_in_dataframe(nmpc_data)
        
    
# class EstimatorDataManger(DynamicDataManger):
#     def __init__(self):
#         print("Save estimated states and actual measurements in estimator's attribute 'dataframe'")
        
#     def construct_estimator_dataframe(self):
#         super(EstimatorDataManger, self).__init__()

#         actualmeas_columns = [item+"_actualmeas" for item in self.meavar_str]
        
#         columns = ["time"] + self.diffvar_str + actualmeas_columns
#         self.dataframe = pd.DataFrame(columns = columns)
        
#     def initialize_estimator_dataframe(self, skip_t0 = False):
#         self.construct_estimator_dataframe()
        
#         if not skip_t0:
#             row0 = [np.nan for ind in range(len(self.dataframe.columns))]
#             rowind = len(self.dataframe.index)
#             self.dataframe.loc[rowind] = row0
        
#     def record_estimator_data(self):
#         tlast = self.time.last()
        
#         diffvar_tlast_list = [self.DIFFERENTIAL_BLOCK[bind].var[tlast].value
#                              for bind in self.DIFFERENTIAL_SET]
#         actualmeavar_tlast_list = [self.ACTUALMEASUREMENT_BLOCK[bind].var[tlast].value
#                                   for bind in self.ACTUALMEASUREMENT_SET]
#         mhe_data = diffvar_tlast_list + actualmeavar_tlast_list
#         self.save_data_in_dataframe(mhe_data)
        
        