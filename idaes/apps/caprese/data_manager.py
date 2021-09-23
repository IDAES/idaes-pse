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

import pyomo.environ as pyo

from pyomo.core.base.componentuid import ComponentUID
from pyomo.common.collections import ComponentMap
import pandas as pd
from collections import OrderedDict
from idaes.apps.caprese.plotlibrary import (
                        PLANT_PlotLibrary,
                        NMPC_PlotLibrary,
                        MHE_PlotLibrary,)

__author__ = "Robert Parker and Kuan-Han Lin"


def empty_dataframe_from_variables(variables, rename_map=None):
    '''
    This function creates a pandas dataframe with variables in columns.

    Parameters
    ----------
    variables : list of variables of interest
    rename_map : dictionary or componentmap that maps the variable to a 
                    specific name string.
        
    Returns
    -------
    dataframe : pandas dataframe

    '''
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
        # KHL:I think we can leave the check alone here and do it when adding values. 
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
        # Change the type of iteration column from float to integer
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
    '''
    This function appends data to the dataframe. 

    Parameters
    ----------
    dataframe : pandas dataframe that we want to append the data to.
    variables : list of variables of interest.
    iteration : current iteration.
    time_subset : time indices of interset in the variables.
    rename_map : dictionary or componentmap that maps the variable to a 
                    specific name string.
    initial_time : initial time for this data saving process.
    time_map : map the time in time_subset to real time points
    '''
    
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
        time_points = [time_map[t] for t in time_subset]
    else:
        time_points = [initial_time + t for t in time_subset]
        
    df = pd.DataFrame(df_map, index=time_points)
    
    return dataframe.append(df)


def add_setpoint_column_into_dataframe(dataframe, variables):
    '''
    This function adds additional columns to save the setpoint of states for NMPC 
    right after the dataframe is constructed.
    
    parameters
    -------------
    dataframe: constructed dataframe
    variables: list of variables that their setpoints will be saved in the dataframe
    '''
    
    for var in variables:
        column_name = str(ComponentUID(var.referent)) + "_setpoint"
        dataframe[column_name] = []
    return dataframe

def add_variable_setpoints_to_dataframe(dataframe, 
                                        variables, 
                                        time_subset, 
                                        map_for_user_given_vars = None):
    '''
    Save the setpoints for states of interest in nmpc's dataframe.
    '''
    
    for var in variables:
        column_name = str(ComponentUID(var.referent)) + "_setpoint"
        # Get column index for "iat", which is used later
        column_ind = dataframe.columns.get_loc(column_name)
        start_row_ind = len(dataframe.index) - len(time_subset)
        for pt in range(len(time_subset)):
            row_ind = start_row_ind + pt
            if var in map_for_user_given_vars:
                # User given variables are not nmpc_var, so they don't have setpoint attribute.
                # User a componentmap to get corresponding nmpc_var.
                dataframe.iat[row_ind, column_ind] = map_for_user_given_vars[var].setpoint
            else:
                dataframe.iat[row_ind, column_ind] = var.setpoint
    return dataframe
        

class PlantDataManager(PLANT_PlotLibrary):
    def __init__(self, 
                 plantblock,
                 user_interested_states = None,
                 user_interested_inputs = None):
                
        # Make sure given variables are all in plantblock
        if user_interested_states is not None:
            cuid_states = [ComponentUID(var.referent) for var in user_interested_states]
            self.plant_user_interested_states = [cuid.find_component_on(plantblock)
                                                 for cuid in cuid_states]
        else:
            self.plant_user_interested_states = []
            
        if user_interested_inputs is not None:
            cuid_inputs = [ComponentUID(var.referent) for var in user_interested_inputs]
            self.plant_user_interested_inputs = [cuid.find_component_on(plantblock)
                                                 for cuid in cuid_inputs]
        else:
            self.plant_user_interested_inputs = []
        
        self.plantblock = plantblock
        self.plant_vars_of_interest = self.plant_user_interested_states + \
                                        self.plant_user_interested_inputs + \
                                            self.plantblock.differential_vars + \
                                                self.plantblock.input_vars

        self.plant_df = empty_dataframe_from_variables(self.plant_vars_of_interest)
        
    def get_plant_dataframe(self):
        return self.plant_df
        
    def save_initial_plant_data(self):
        plant_states_to_save = self.plant_user_interested_states + self.plantblock.differential_vars
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
                 user_interested_states = None,
                 user_interested_inputs = None):
        
        if not hasattr(self, "plant_df"):
            super(ControllerDataManager, self).__init__(plantblock, 
                                                        user_interested_states, 
                                                        user_interested_inputs,)
        
        self.controllerblock = controllerblock
        # Convert vars in plant to vars in controller
        if user_interested_states is not None:
            cuid_states = [ComponentUID(var.referent) for var in user_interested_states]
            self.controller_user_interested_states = [cuid.find_component_on(controllerblock)
                                                      for cuid in cuid_states]
        else:
            self.controller_user_interested_states = []
            
        if user_interested_inputs is not None:
            cuid_inputs = [ComponentUID(var.referent) for var in user_interested_inputs]
            self.controller_user_interested_inputs = [cuid.find_component_on(controllerblock)
                                                      for cuid in cuid_inputs]
        else:
            self.controller_user_interested_inputs = []
        
        
        self.controller_vars_of_interest = self.controller_user_interested_inputs + \
                                                controllerblock.input_vars   
        self.controller_df = empty_dataframe_from_variables(self.controller_vars_of_interest)
        
        self.states_need_setpoints = self.controller_user_interested_states + \
                                        controllerblock.differential_vars
        self.controller_df = add_setpoint_column_into_dataframe(self.controller_df, 
                                                                self.states_need_setpoints)
        
        # Important!!!
        # User given variables are not nmpc_var, so they don't have setpoint attribute.
        # Create this componentmap to map the given vars to corresponding nmpc_vars.
        vardata_map = self.controllerblock.vardata_map
        t0 = self.controllerblock.time.first()
        self.user_given_vars_map_nmpcvar = ComponentMap((var, vardata_map[var[t0]]) 
                                                        for var in self.controller_user_interested_states)

  
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
                 user_interested_states = None,):
        
        if not hasattr(self, "plant_df"):
            super(EstimatorDataManager, self).__init__(plantblock, 
                                                       user_interested_states,)
        
        self.estimatorblock = estimatorblock
        # Convert vars in plant to vars in estimator
        if user_interested_states is not None:
            cuid_states = [ComponentUID(var.referent) for var in user_interested_states]
            self.estimator_user_interested_states = [cuid.find_component_on(estimatorblock)
                                                     for cuid in cuid_states]
        else:
            self.estimator_user_interested_states = []
            
        self.estimator_vars_of_interest = self.estimator_user_interested_states + \
                                                estimatorblock.differential_vars
        self.estimator_df = empty_dataframe_from_variables(self.estimator_vars_of_interest)
        
    def get_estimator_dataframe(self):
        return self.estimator_df
    
    def save_estimator_data(self, iteration):
        time = self.estimatorblock.time
        t_last = time.last()
        # Estimation starts from the very first sample time.
        time_map = OrderedDict([(t_last, (iteration+1)*self.estimatorblock.sample_time),])
        self.estimator_df = add_variable_values_to_dataframe(self.estimator_df,
                                                             self.estimator_vars_of_interest,
                                                             iteration,
                                                             time_subset = [t_last],
                                                             time_map = time_map,)
        
class DynamicDataManager(ControllerDataManager, EstimatorDataManager):
    def __init__(self,
                 plantblock,
                 controllerblock,
                 estimatorblock,
                 user_interested_states = None,
                 user_interested_inputs = None,):
        
        # Create plant dataframe
        super(EstimatorDataManager, self).__init__(plantblock,
                                                   user_interested_states,
                                                   user_interested_inputs,)
        # Create estimator dataframe
        super(ControllerDataManager, self).__init__(plantblock = None, # plant dataframe is already there
                                                    estimatorblock = estimatorblock,
                                                    user_interested_states = user_interested_states,)
        # Create controller dataframe
        super(DynamicDataManager, self).__init__(plantblock = None, # plant dataframe is already there
                                                 controllerblock = controllerblock,
                                                 user_interested_states = user_interested_states,
                                                 user_interested_inputs = user_interested_inputs,)
