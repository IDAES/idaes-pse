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
import pandas as pd
import numpy as np
from scipy import sparse
import pyomo.contrib.parmest.parmest as parmest
from pyomo.environ import *
from pyomo.opt import SolverFactory
import shutil
import logging
from collections import namedtuple

logger = logging.getLogger('idaes.apps.uncertainty_propagation')

def quantify_propagate_uncertainty(model_function, model_uncertain,  data, theta_names, obj_function=None, 
                 tee=False, diagnostic_mode=False, solver_options=None):
    """This function calculates error propagation of the objective function and constraints. 
    The parmest uses 'model_function' to estimate uncertain parameters. The uncertain parameters in 
    'model_uncertain' are fixed with the estimated values. The function 'quantify_propagate_uncertainty' 
    calculates error propagation of objective function and constraints in the 'model_uncertain'.    
    
    Parameters
    ----------
    model_function: function
        A python Function that generates an instance of the Pyomo model using
        'data' as the input argument
    model_uncertain: function or Pyomo ConcreteModel
        function is a a python Function that generates an instance of the Pyomo model
    data: pandas DataFrame, list of dictionaries, or list of json file names
        Data that is used to build an instance of the Pyomo model and 
        build the objective function
    theta_names: list of strings
        List of Var names to estimate
    obj_function: function, optional
        Function used to formulate parameter estimation objective, 
        generally sum of squared error between measurements and model variables.    
    tee: bool, optional
        Indicates that ef solver output should be teed
    diagnostic_mode: bool, optional
        If True, print diagnostics from the solver
    solver_options: dict, optional
        Provides options to the solver (also the name of an attribute)
    
    Returns
    -------
    results : namedtuple
        results.obj: float
            objective function value for the given obj_function 
        results.theta_out: dictionary
            Estimated parameters
        results.cov: numpy.ndarray
            Covariance of theta_out
        results.propagation_f: dict
            error propagation in the objective function with the dictionary key, 'objective' 
            if the objective function doesn't include any uncertain parameters, return an empty dictionary        
        results.propagation_c: dict
            error propagation in the constraints with the dictionary keys, 'constraints l'
            where l is the line number. 
            if no constraint includes uncertain parameters, return an empty dictionary      
        
    Raises
    ------
    ValueError
        When tee entry is not Boolean
    ValueError
        When diagnostic_mode entry is not Boolean
    ValueError
        When solver_options entry is not None and a Dictionary
    Warnings
        When an element of theta_names includes a space
    """
    if not isinstance(tee, bool):
        raise ValueError('tee  must be boolean.')
    if not isinstance(diagnostic_mode, bool):
        raise ValueError('diagnostic_mode  must be boolean.')    
    if not solver_options==None:
        if not isinstance(solver_options, dict):
            raise ValueError('solver_options must be dictionary.')

    
    # Remove all "'" and " " in theta_names
    variable_clean = False
    for _ in theta_names:
        if " " in _ or "'" in _:
            variable_clean = True    
    if variable_clean:
        theta_names, var_dic = clean_variable_name(theta_names)
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
    propagation_f, propagation_c =  propagate_uncertainty(model_uncertain, theta, cov, theta_names, tee)
    
    Output = namedtuple('Output',['obj', 'theta_out', 'cov', 'propagation_f', 'propagation_c'])
    results= Output(obj, theta_out, cov, propagation_f, propagation_c)
    
    return results

def propagate_uncertainty(model_uncertain, theta, cov, theta_names, tee=False, solver_options=None):
    """This function calculates gradient vector, expectation, and variance of the objective function and constraints
    of the model for given estimated optimal parameters and covariance matrix of parameters. It calculates 
    error propagation of the objective function and constraints by using gradient vector and covariance matrix. 
    
    Parameters
    ----------
    model_uncertain: function or Pyomo ConcreteModel
        function is a a python Function that generates an instance of the Pyomo model
    theta: dict
        Estimated parameters 
    cov: numpy.ndarray
        Covariance matrix of parameters 
    theta_names: list of strings
        List of estimated Var names
    var_dic: dictionary
        If any original variable contains "'", need an auxiliary dictionary  with keys theta_namess without "'", values with "'".
        e.g) var_dic: {'fs.properties.tau[benzene,toluene]': "fs.properties.tau['benzene','toluene']", 'fs.properties.tau[toluene,benzene]': "fs.properties.tau['toluene','benzene']"} 
    tee: bool, optional
        Indicates that ef solver output should be teed
    solver_options: dict, optional
        Provides options to the solver (also the name of an attribute)
    
    Returns
    -------
    propagation_f: dictionary
        error propagation in the objective function with the dictionary key, 'objective' 
        if the objective function doesn't include any uncertain parameters, return an empty dictionary        
    propagation_c: dictionary
        error propagation in the constraints with the dictionary keys, 'constraints r'
        where r is the line number. 
        if no constraint includes uncertain parameters, return an empty dictionary      
 
    Raises
    ------
    Exception
        When model_uncertain is neither 'ConcreteModel' nor 'function'.
    """
    # define Pyomo model
    if type(model_uncertain).__name__ == 'ConcreteModel':
        model = model_uncertain
    elif type(model_uncertain).__name__ == 'function':
        model = model_uncertain()
    else:
        raise Exception('model_uncertain must be either python function or Pyomo ConcreteModel.')

    # Remove all "'" in theta_names     
    variable_clean = False
    for _ in theta_names:
        if " " in _ or "'" in _:
            variable_clean = True
    if variable_clean:
        theta_names, var_dic = clean_variable_name(theta_names)
        # Fix the estiamted paramerters 
        for v in theta_names:
            eval('model.'+var_dic[v]).setlb(theta[v])
            eval('model.'+var_dic[v]).setub(theta[v])
    else:
        for v in theta_names:
            eval('model.'+v).setlb(theta[v])
            eval('model.'+v).setub(theta[v])
    # get gradient of the objective function, constraints, and the column number of each theta
    gradient_f, gradient_c,line_dic = get_sensitivity(model, theta_names, tee)
    
    # calculate error propagation g*cov*g 
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
        propagation_f['objective'] = np.dot(gradient_f_theta,np.dot(cov,np.transpose(gradient_f_theta)))
   

    # constraints: 
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
        gradient_dic = {}
        propagation_c = {}
        
        # calculate error propagation of constraints
        for r in constriant_number:
            gradient_dic['constraints '+str(r)] = gradient_cc[r-1]
            propagation_c['constraints '+str(r)] = float(np.dot(gradient_cc[r-1],np.dot(cov,np.transpose(gradient_cc[r-1]))))
    else:
        propagation_c = {}

    return propagation_f, propagation_c

def get_sensitivity(model, theta_names, tee=False, solver_options=None):
    """This function calculates gradient vector of the objective function and constraints with respect to the variables in theta_names.
    
    Parameters
    ----------
    model: Pyomo ConcreteModel
        model should includes an objective function 
    theta_names: list of strings
        List of Var names
    tee: bool, optional
        Indicates that ef solver output should be teed
    solver_options: dict, optional
        Provides options to the solver (also the name of an attribute)
    
    Returns
    -------
    gradient_f: numpy.ndarray
        gradient vector of the objective function with respect to all decision variables at the optimal solution
    gradient_c: numpy.ndarray
        gradient vector of the constraints with respect to all decision variables at the optimal solution
        Each row contains column number, row number, and value
        If no constraint exists, return []
    line_dic: dict
        column numbers of the theta_names in the model. Index starts from 1
    Raises
    ------
    RuntimeError
        When ipopt or kaug or dotsens is not available
    Exception
        When ipopt fails 
    """
    #Create the solver plugin using the ASL interface
    ipopt = SolverFactory('ipopt',solver_io='nl')
    if solver_options is not None:
        ipopt.options = solver_options
    kaug = SolverFactory('k_aug',solver_io='nl')
    dotsens = SolverFactory('dot_sens',solver_io='nl')
    if not ipopt.available(False):
        raise RuntimeError('ipopt is not available')
    if not kaug.available(False):
        raise RuntimeError('k_aug is not available')
    if not dotsens.available(False):
        raise RuntimeError('dotsens is not available')
    
    # Declare Suffixes
    model.dual = Suffix(direction = Suffix.IMPORT)
    model.ipopt_zL_out = Suffix(direction=Suffix.IMPORT)
    model.ipopt_zU_out = Suffix(direction=Suffix.IMPORT)
    model.ipopt_zL_in = Suffix(direction=Suffix.EXPORT)
    model.ipopt_zU_in = Suffix(direction=Suffix.EXPORT)

    # K_AUG SUFFIXES
    model.dof_v = Suffix(direction=Suffix.EXPORT)  #: SUFFIX FOR K_AUG
    model.rh_name = Suffix(direction=Suffix.IMPORT)  #: SUFFIX FOR K_AUG AS WELL
    kaug.options["print_kkt"] = ""
    results = ipopt.solve(model,tee=tee)

    # Rasie Exception if ipopt fails 
    if (results.solver.status == pyomo.opt.SolverStatus.warning):
        raise Exception(results.solver.Message)
    
    for o in model.component_objects(Objective, active=True):
        f_mean = value(o)
    model.ipopt_zL_in.update(model.ipopt_zL_out)
    model.ipopt_zU_in.update(model.ipopt_zU_out)
    #: run k_aug
    kaug.solve(model, tee=tee)  #: always call k_aug AFTER ipopt.
    model.write('col_row.nl', format='nl', io_options={'symbolic_solver_labels':True})
    # get the column numbers of theta
    line_dic = {}
    for v in theta_names:
        line_dic[v] = line_num('col_row.col', v)
    # load gradient of the objective function
    gradient_f = np.loadtxt("./GJH/gradient_f_print.txt")
    # load gradient of all constraints (sparse)
    # If no constraint exists, return []
    num_constraints = len(list(model.component_data_objects(Constraint,
                                                            active=True,
                                                            descend_into=True)))
    if num_constraints > 0 :
        gradient_c = np.loadtxt("./GJH/A_print.txt")
    else:
        gradient_c = np.array([])

    # remove all generated files
    shutil.move("col_row.nl", "./GJH/")
    shutil.move("col_row.col", "./GJH/")
    shutil.move("col_row.row", "./GJH/")
    shutil.rmtree('GJH', ignore_errors=True)
    return gradient_f, gradient_c, line_dic

def line_num(file_name, target):
    """This function returns the line number contains 'target' in the file_name
    This function identities constraints that have variables in theta_names.     

    Parameters
    ----------
    file_name: string
        file name includes information of variabe order (col_row.col)
    target: string   
        variable name to check  
    Returns
    -------
    count: int
        line number of target in the file_name
        
    Raises
    ------
    Exception
        When col_row.col doesnot include target
    """
    with open(file_name) as f:
        count = int(1)
        for line in f:
            if line.strip() == target:
                return int(count)
            count += 1
    raise Exception("col_row.col should includes target")

def clean_variable_name(theta_names):
    """This eunction removes all ' and spaces in theta_names.
       Note: The  current theta_est(calc_cov=True) of parmest in Pyomo
       doesn't allow ' and spaces in the variable names    
       
    Parameters
    ----------
    theta_names: list of strings
        List of Var names
    
    Returns
    -------
    theta_names_out: list of strings
        List of Var names after removing  all ' and spaces
    var_dic: dict
       dictionary with keys converted theta_names and values origianl theta_names 
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
    return theta_names_out, var_dic

