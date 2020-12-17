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


def quantify_propagate_unucertainty(model_function, model_uncertain,  data, theta_names, obj_function=None, 
                 tee=False, diagnostic_mode=False, solver_options=None):
    """This function calculates error propagation of the objective function and constraints.  
    
    Parameters
    ----------
    model_function: function
        Function that generates an instance of the Pyomo model using
        'data' as the input argument
    model_uncertain: function
        Function that generates an instance of the Pyomo model using
        'theta' and 'theta_names'  as the input arguments
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
    obj: float
        objective function value for the given obj_function 
    theta: dictionary
        Estimated parameters
    cov: numpy.ndarray
        Covariance of theta
    propagation_f: dict
        error propagation in the objective function with the dictionary key, 'objective' 
        if the objective function doesn't include any uncertain parameters, return an empty dictionary        
    propagation_c: dict
        error propagation in the constraints with the dictionary keys, 'constraints l'
        where l is the line number. 
        if no constraint includes uncertain parameters, return an empty dictionary      
        
    Raises
    ------
    Exception
        When tee entry is not Boolean
    Exception
        When diagnostic_mode entry is not Boolean
    Exception
        When solver_options entry is not None and a Dictionary
    """
    if not isinstance(tee, bool):
        raise Exception('tee  must be boolean.')
    if not isinstance(diagnostic_mode, bool):
        raise Exception('diagnostic_mode  must be boolean.')    
    if not solver_options==None:
        if not isinstance(solver_options, dict):
            raise Exception('solver_options must be dictionary.')

    parmest_class = parmest.Estimator(model_function, data, theta_names, obj_function,
                 tee, diagnostic_mode, solver_options)
    obj, theta, cov = parmest_class.theta_est(calc_cov=True)

    propagation_f, propagation_c =  propagate_uncertainty(model_uncertain, theta, cov, theta_names)
    return obj, theta, cov, propagation_f, propagation_c

def propagate_uncertainty(model_uncertain, theta, cov, theta_names, tee=False, solver_options=None):
    """This function calculates gradient vector, expectation, and variance of the objective function 
    of the model for given estimated optimal parameters and covariance matrix of parameters. It calculates 
    error propagation of the objective function and constraints by using gradient vector and covariance matrix. 
    
    Parameters
    ----------
    model_uncertain: function
        Function that generates an instance of the Pyomo model using
        'theta' and 'theta_names;  as the input argument
    theta: dict
        Estimated parameters 
    cov: numpy.ndarray
        Covariance matrix of parameters 
    theta_names: list of strings
        List of estimated Var names
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
    RuntimeError
        When ipopt or kaug or dotsens is not available
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

    model = model_uncertain(theta,theta_names)
    
    model.dual = Suffix(direction=Suffix.IMPORT_EXPORT)
    model.ipopt_zL_out = Suffix(direction=Suffix.IMPORT)
    model.ipopt_zU_out = Suffix(direction=Suffix.IMPORT)
    model.ipopt_zL_in = Suffix(direction=Suffix.EXPORT)
    model.ipopt_zU_in = Suffix(direction=Suffix.EXPORT)

    #: K_AUG SUFFIXES
    model.dof_v = Suffix(direction=Suffix.EXPORT)  #: SUFFIX FOR K_AUG
    model.rh_name = Suffix(direction=Suffix.IMPORT)  #: SUFFIX FOR K_AUG AS WELL
    kaug.options["print_kkt"] = ""
    ipopt.solve(model,tee=tee)
    for o in model.component_objects(Objective, active=True):
        f_mean = value(o)
    model.ipopt_zL_in.update(model.ipopt_zL_out)
    model.ipopt_zU_in.update(model.ipopt_zU_out)
    #: k_aug
    kaug.solve(model, tee=tee)  #: always call k_aug AFTER ipopt.
    model.write('col_row.nl', format='nl', io_options={'symbolic_solver_labels':True})
    # get the column numbers of theta
    line_dic = {} 
    for v in theta_names:
        line_dic[v] = line_num('col_row.col', v)
    # load gradient of the objective function
    gradient_f = np.loadtxt("./GJH/gradient_f_print.txt")
    # load gradient of all constraints (sparse)
    gradient_c = np.loadtxt("./GJH/A_print.txt")
    # remove all generated files
    shutil.move("col_row.nl", "./GJH/")
    shutil.move("col_row.col", "./GJH/")
    shutil.move("col_row.row", "./GJH/")
    shutil.rmtree('GJH', ignore_errors=True)
    
    # objective function: get gradients of theta and calculate error propagation g*cov*g 
    for cc in list(model.component_data_objects(Objective,
                                            active=True,
                                            descend_into=True)):
        if any(ele in str(cc.expr) for ele in theta_names):
            gradient_f_theta = gradient_f[[int(line_dic[v]) -1 for v in theta_names]]
    propagation_f = {}
    if 'gradient_f_theta' in locals():
        propagation_f['objective'] = np.dot(gradient_f_theta,np.dot(cov,np.transpose(gradient_f_theta)))
   

    # constraintsn: get gradients of theta and calculate error propagation g*cov*g
    num_constraints = len(list(model.component_data_objects(Constraint,
                                                            active=True,
                                                            descend_into=True))) 
    if num_constraints > 0:
        gradient_cc = []
        for i in range(gradient_c.shape[0]):
            if gradient_c[i,[0]] in line_dic.values():
                gradient_cc.append(gradient_c[i])
        gradient_cc = np.array(gradient_cc)
        I,J,V =  gradient_cc[:,[1]].flatten().astype(int)-1, gradient_cc[:,[0]].flatten().astype(int)-1,  gradient_cc[:,[2]].flatten()
        constriant_number = list(set(gradient_cc[:,[1]].flatten().astype(int)))

        gradient_cc = sparse.coo_matrix((V,(I,J))).todense()
        gradient_dic = {}
        propagation_c = {}
        for r in constriant_number:
            gradient_dic['constraints '+str(r)] = gradient_cc[r-1]
            propagation_c['constraints '+str(r)] = float(np.dot(gradient_cc[r-1],np.dot(cov,np.transpose(gradient_cc[r-1]))))
    else:
        propagation_c = {}

    return propagation_f, propagation_c


def line_num(file_name, target):
    """This function returns the line number contains 'target' in the file_name
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


if __name__ == '__main__':
    # min (y-f(x,theta))^2
    from rooney_biegler import rooney_biegler_model,rooney_biegler_model_opt
    variable_name = ['asymptote', 'rate_constant']
    data = pd.DataFrame(data=[[1,8.3],[2,10.3],[3,19.0],
                              [4,16.0],[5,15.6],[7,19.8]],
                        columns=['hour', 'y'])

    def SSE(model, data):
        expr = sum((data.y[i] - model.response_function[data.hour[i]])**2 for i in data.index)
        return expr

    obj, theta, cov, propagation_f, propagation_c =  quantify_propagate_unucertainty(rooney_biegler_model,rooney_biegler_model_opt, data, variable_name, SSE)
