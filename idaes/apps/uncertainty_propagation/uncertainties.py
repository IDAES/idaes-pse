##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################

import idaes
import os
import pandas as pd
import numpy as np
from scipy import sparse
import pyomo.contrib.parmest.parmest as parmest
from pyomo.environ import *
from pyomo.opt import SolverFactory
import shutil
import logging
from collections import namedtuple
from idaes.apps.uncertainty_propagation.sens import SensitivityInterface, get_dsdp, get_dfds_dcds # will replace with pyomo (# Pyomo PR 1613: https://github.com/pyomo/pyomo/pull/1613/) 

logger = logging.getLogger('idaes.apps.uncertainty_propagation')

def quantify_propagate_uncertainty(model_function, model_uncertain,  data, theta_names, obj_function=None, tee=False, diagnostic_mode=False, solver_options=None):
    """
    This function calculates error propagation of the objective function and constraints. The parmest uses 'model_function' to estimate uncertain parameters. The uncertain parameters in 'model_uncertain' are fixed with the estimated values. The function 'quantify_propagate_uncertainty' calculates error propagation of objective function and constraints in the 'model_uncertain'.
    
    Args:
        model_function(function)                                                  : A python Function that generates an instance of the Pyomo model using 'data' as the input argument
        model_uncertain(function or Pyomo ConcreteModel)                          : Function is a a python Function that generates an instance of the Pyomo model
        data(pandas DataFrame, list of dictionaries, or list of json file names)  : Data that is used to build an instance of the Pyomo model and build the objective function
        theta_names(list of strings)                                              : List of Var names to estimate
        obj_function(function, optional)                                          : Function used to formulate parameter estimation objective, generally sum of squared error between measurements and model variables.    
        tee(bool, optional)                                                       : Indicates that ef solver output should be teed
        diagnostic_mode(bool, optional)                                           : If True, print diagnostics from the solver
        solver_options(dict, optional)                                            : Provides options to the solver (also the name of an attribute)
    
    Returns:
        tuple   : results object containing the all information including:
            - results.obj(float)                    : Objective function value for the given obj_function 
            - results.theta(dict)                   : Estimated parameters
            - results.theta_names(dict)             : Names of parameters
            - results.cov(numpy.ndarray)            : Covariance of theta
            - results.gradient_f(numpy.ndarray)     : Gradient vector of the objective function with respect to the (decision variables, parameters) at the optimal solution
            - results.gradient_c(scipy.sparse.coo.coo_matrix)     : Gradient vector of the constraints with respect to the (decision variables, parameters) at the optimal solution. Each row contains [column number, row number, and value], colum order follows variable order in col. If no constraint exists, return []. Note: Python sparse matrix index starts from 0.  To match with k_aug sparse form, need to add 1 for each row and col index.
            - results.dsdp(numpy.ndarray)           : Gradient vector of the (decision variables, parameters) with respect to paramerters (=theta_name). number of rows = len(theta_name), number of columns= len(col)
            - results.propagation_c(numpy.ndarray)  : Error propagation in the constraints (only constraints have theta_name)  e.g) [[constraint iumber1, error1],[constraint number2, error2]]             
            - results.propagation_f(numpy.float64)  : Error propagation in the objective function, df/dp*cov_p*df/dp + (df/dx*dx/dp)*cov_p*(df/dx*dx/dp)
            - results.col(list)                     : List of variable names  
    
    Raises:
        TypeError:  When tee entry is not Boolean
        TypeError:  When diagnostic_mode entry is not Boolean
        TypeError:  When solver_options entry is not None and a Dictionary
        Warnings:   When an element of theta_names includes a space
    
    """

    if not isinstance(tee, bool):
        raise TypeError('tee  must be boolean.')
    if not isinstance(diagnostic_mode, bool):
        raise TypeError('diagnostic_mode  must be boolean.')    
    if not solver_options==None:
        if not isinstance(solver_options, dict):
            raise TypeError('solver_options must be dictionary.')
    # Remove all "'" and " " in theta_names
    theta_names, var_dic,variable_clean = clean_variable_name(theta_names)
    parmest_class = parmest.Estimator(model_function, data, theta_names, obj_function,
                 tee, diagnostic_mode, solver_options)
    obj, theta,cov = parmest_class.theta_est(calc_cov=True)
    # Convert theta keys to the original name 
    # Revert theta_names to be original
    if variable_clean:
        theta_out = {}
        for i in var_dic.keys():
            theta_out[var_dic[i]] = theta[i]
        theta_names = [var_dic[v] for v in theta_names] 
    else:
        theta_out = theta
    propagate_results  =  propagate_uncertainty(model_uncertain, theta, cov, theta_names, tee)
   
    if len(propagate_results.gradient_c)>0:
        I,J,V =  propagate_results.gradient_c[:,1].flatten().astype(int)-1, propagate_results.gradient_c[:,0].flatten().astype(int)-1, propagate_results.gradient_c[:,2].flatten()
        coo_ = sparse.coo_matrix((V,(I,J)))
    else:
        coo_ = propagate_results.gradient_c

    Output = namedtuple('Output',['obj', 'theta', 'theta_names', 'cov','gradient_f', 'gradient_c', 'dsdp', 'propagation_c', 'propagation_f','col'])
    results= Output(obj, theta_out, theta_names, cov, 
                     propagate_results.gradient_f,
                     propagate_results.gradient_c, 
                     propagate_results.dsdp,
                     coo_,
                     propagate_results.propagation_f,
                     propagate_results.col)
    return results

def propagate_uncertainty(model_uncertain, theta, cov, theta_names, tee=False, solver_options=None):
    """
    This function calculates gradient vector, expectation, and variance of the objective function and constraints  of the model for given estimated optimal parameters and covariance matrix of parameters. It calculates error propagation of the objective function and constraints by using gradient vector and covariance matrix.
    
    Args:
        model_uncertain(function or Pyomo ConcreteModel) : Function is a a python Function that generates an instance of the Pyomo model
        theta(dict)                    : Estimated parameters 
        cov(numpy.ndarray)             : Covariance matrix of parameters 
        theta_names(list of strings)   : List of estimated Var names
        var_dic(dictionary)            : If any original variable contains apostrophes, need an auxiliary dictionary with keys theta_namess without apostrophes, values with apostrophes. e.g) var_dic: {'fs.properties.tau[benzene,toluene]': "fs.properties.tau['benzene','toluene']", 'fs.properties.tau[toluene,benzene]': "fs.properties.tau['toluene','benzene']"} 
        tee(bool, optional)            : Indicates that ef solver output should be teed
        solver_options(dict, optional) : Provides options to the solver (also the name of an attribute)

     Returns:
         tuple   : results object containing the all information including
            - results.gradient_f(numpy.ndarray)    : Gradient vector of the objective function with respect to the (decision variables, parameters) at the optimal solution
            - results.gradient_c(numpy.ndarray)    : Gradient vector of the constraints with respect to the (decision variables, parameters) at the optimal solution. Each row contains [column number, row number, and value], colum order follows variable order in col. If no constraint exists, return []
            - results.dsdp(numpy.ndarray)          : Gradient vector of the (decision variables, parameters) with respect to paramerters (=theta_name). number of rows = len(theta_name), number of columns= len(col)
            - results.propagation_c(numpy.ndarray) : Error propagation in the constraints (only constraints have theta_name)  e.g) [[constraint iumber1, error1],[constraint number2, error2]]             
            - results.propagation_f(numpy.float64) : Error propagation in the objective function, df/dp*cov_p*df/dp + (df/dx*dx/dp)*cov_p*(df/dx*dx/dp)
            - results.col(list)                    : List of variable names  
    Raises:
        Exception:  When model_uncertain is neither 'ConcreteModel' nor 'function'.
    """

    # define Pyomo model
    if type(model_uncertain).__name__ == 'ConcreteModel':
        model = model_uncertain
    elif type(model_uncertain).__name__ == 'function':
        model = model_uncertain()
    else:
        raise TypeError('model_uncertain must be either python function or Pyomo ConcreteModel.')
    # Remove all "'" in theta_names     
    theta_names, var_dic,variable_clean = clean_variable_name(theta_names)
    for v in theta_names:
        eval('model.'+var_dic[v]).setlb(theta[v])
        eval('model.'+var_dic[v]).setub(theta[v])
    # get gradient of the objective function, constraints, and the column number of each theta
    gradient_f, gradient_c, col, line_dic = get_dfds_dcds(model, theta_names, tee)
    dsdp, col  = get_dsdp(model, theta_names, theta, var_dic,tee)      
    
    gradient_f_dic = {}
    for i in range(len(col)):
        gradient_f_dic["d(f)/d("+col[i]+")"] = gradient_f[i]
    num_constraints = len(list(model.component_data_objects(Constraint,
                                                            active=True,
                                                            descend_into=True)))
    gradient_c_dic = {}
    if num_constraints > 0 :
        gradient_c = np.array([i for i in gradient_c if not np.isclose(i[2],0)])
        row_number, col_number = np.shape(gradient_c)
        gradient_c_dic = {}
        for i in range(row_number):
            gradient_c_dic["d(c"+ str(int(gradient_c[i,1]))+")/d("+col[int(gradient_c[i,0]-1)]+")"] = gradient_c[i,2]
    dsdp_dic = {}
    for i in range(len(theta_names)):
        for j in range(len(col)):
            if SensitivityInterface.get_default_block_name() not in col[j]:
                dsdp_dic["d("+col[j] +")/d("+theta_names[i]+")"] =  -dsdp[i, j]
    x_list = [x for x in col if x not in var_dic.keys()] #x_list includes only decision variables (exclude theta_names) 
    p_list = theta_names #This makes cov and fxxp have the same order 
    fxxp = []
    for i in p_list:
        fxxp_tmp = 0
        for j in x_list:
            fxxp_tmp = fxxp_tmp + gradient_f_dic['d(f)/d('+j+')']*dsdp_dic['d('+j+')/d('+i+')']
        fxxp.append(fxxp_tmp)
    fxxp = np.array(fxxp)
    # calculate error propagation df/dp*cov_p*df/dp + (df/dx*dx/dp)*cov_p*(df/dx*dx/dp)
    # objective function: 
    for cc in list(model.component_data_objects(Objective,
                                            active=True,
                                            descend_into=True)):
        # save gradient of the objective with respect to only theta_names
        if any(ele in str(cc.expr) for ele in theta_names):
            gradient_f_theta = gradient_f[[int(line_dic[v]) -1 for v in theta_names]]
    propagation_f = {}
    # calculate error propagation of the objective fuction
    if 'gradient_f_theta' in locals():
        propagation_f = np.dot(gradient_f_theta,np.dot(cov,np.transpose(gradient_f_theta))) + np.dot(fxxp,np.dot(cov,np.transpose(fxxp)))
    else:
        propagation_f = np.dot(fxxp,np.dot(cov,np.transpose(fxxp)))
    # constraints: 
    # calculate error propagation dc/dp*cov_p*dc/dp + (dc/dx*dx/dp)*cov_p*(dc/dx*dx/dp)
    num_constraints = len(list(model.component_data_objects(Constraint,
                                                            active=True,
                                                            descend_into=True))) 
    if num_constraints > 0:
        # save constraints into gradient_cc only if constraints includes theta 
        # gradient_c[i] is a sparse matrix with column number, row number, and value
        gradient_cc = []
        for i in range(gradient_c.shape[0]):
            if gradient_c[i,[0]] in line_dic.values():
                gradient_cc.append(gradient_c[i])
        gradient_cc = np.array(gradient_cc)
        # save unique constraints numbers that contain theta, index starts from 1. 
        constriant_number = list(set(gradient_cc[:,[1]].flatten().astype(int)))       
        # convert sparse matrix to dense matrix with (value, (row number, column number))
        I,J,V =  gradient_cc[:,[1]].flatten().astype(int)-1, gradient_cc[:,[0]].flatten().astype(int)-1,  gradient_cc[:,[2]].flatten()
        gradient_cc = sparse.coo_matrix((V,(I,J))).todense()
        # propagation_c includes `dc/ds*cov_p*dc/ds`
        propagation_c = []
        # calculate error propagation of constraints
        for r in range(1,num_constraints+1):
            if r in constriant_number:
                propagation_c.append([int(r), float(np.dot(gradient_cc[r-1],np.dot(cov,np.transpose(gradient_cc[r-1]))))])
    else:
        propagation_c = []
    propagation_c = np.array(propagation_c)
    Output = namedtuple('Output',['gradient_f', 'gradient_c', 'dsdp', 'propagation_c','propagation_f','col'])
    results= Output(gradient_f, gradient_c, dsdp, propagation_c, propagation_f,col)
    return results


def clean_variable_name(theta_names):
    """
    This function removes all ' and spaces in theta_names. Note: The  current theta_est(calc_cov=True) of parmest in Pyomo doesn't allow ' and spaces in the variable names.
       
    Args:
        theta_names(list of strings) : List of Var names
    
    Returns:
            : It returns the following variables:
            
            - theta_names_out(list of strings) : List of Var names after removing  all ' and spaces
            - var_dic(dict)                    : Dictionary with keys converted theta_names and values origianl theta_names 
    """

    # Variable names cannot have "'" for parmest_class.theta_est(calc_cov=True)
    # Save original variables name in to var_dic
    # Remove all "'" and " " in theta_names
    var_dic = {}
    theta_names_out = []
    clean = False
    for i in range(len(theta_names)):
        theta_tmp = theta_names[i].replace("'", '')
        theta_tmp = theta_tmp.replace(" ", '')
        theta_names_out.append(theta_tmp)
        var_dic[theta_tmp] = theta_names[i]
        if "'" in theta_names[i] or " " in theta_names[i] :
            logger.warning(theta_names[i]+ " includes ' or space.")
            logger.warning("The cleaned name: " + theta_tmp)
            clean = True
    if clean:
       logger.warning("All ' and spaces in theta_names are removed.")
    return theta_names_out, var_dic, clean
