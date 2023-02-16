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
from idaes.apps.uncertainty_propagation.sens import (
    SensitivityInterface,
    get_dsdp,
    get_dfds_dcds,
)

# will replace with pyomo
# (Pyomo PR 1613: https://github.com/pyomo/pyomo/pull/1613/)

logger = logging.getLogger("idaes.apps.uncertainty_propagation")


def quantify_propagate_uncertainty(
    model_function,
    model_uncertain,
    data,
    theta_names,
    obj_function=None,
    tee=False,
    diagnostic_mode=False,
    solver_options=None,
    covariance_n=None,
):
    """This function calculates error propagation of the objective function and
    constraints. The parmest uses 'model_function' to estimate uncertain
    parameters. The uncertain parameters in 'model_uncertain' are fixed with
    the estimated values. The function 'quantify_propagate_uncertainty'
    calculates error propagation of objective function and constraints
    in the 'model_uncertain'.

    The following terms are used to define the output dimensions:
    Ncon   = number of constraints
    Nvar   = number of variables (Nx + Ntheta)
    Nx     = the number of decision (primal) variables
    Ntheta = number of uncertain parameters.

    Parameters
    ----------
    model_function : function
        A python Function that generates an instance of the Pyomo model using
        'data' as the input argument
    model_uncertain : function or Pyomo ConcreteModel
        Function is a python/ Function that generates an instance of the
        Pyomo model
    data : pandas DataFrame, list of dictionary, or list of json file names
        Data that is used to build an instance of the Pyomo model and build the
        objective function
    theta_names : list of strings
        List of Var names to estimate
    obj_function : function, optional
        Function used to formulate parameter estimation objective, generally
        sum of squared error between measurements and model variables,
        by default None
    tee : bool, optional
        Indicates that ef solver output should be teed, by default False
    diagnostic_mode : bool, optional
        If True, print diagnostics from the solver, by default False
    solver_options : dict, optional
        Provides options to the solver (also the name of an attribute),
        by default None
    covariance_n : int, optional
        Number of datapoints to use in the objective function to
        calculate the covariance matrix.  If omitted, defaults to
        len(data)

    Returns
    -------
    tuple
        results object containing the all information including

        - results.obj: float
            Real number. Objective function value for the given obj_function
        - results.theta: dict
            Size Ntheta python dictionary. Estimated parameters
        - results.theta_names: list
            Size Ntheta list. Names of parameters
        - results.cov: numpy.ndarray
            Ntheta by Ntheta matrix. Covariance of theta
        - results.gradient_f: numpy.ndarray
            Length Nvar array. Gradient vector of the objective function with
            respect to the (decision variables, parameters) at the optimal
            solution
        - results.gradient_c: scipy.sparse.csr.csr_matrix
            Ncon by Nvar size sparse matrix. Gradient vector of the constraints
            with respect to the (decision variables, parameters) at the optimal
            solution.
        - results.dsdp: scipy.sparse.csr.csr_matrix
            Ntheta by Nvar size sparse matrix. Gradient vector of the
            (decision variables, parameters) with respect to paramerters
            (=theta_name). number of rows = len(theta_name),
            number of columns= len(col)
        - results.propagation_c: numpy.ndarray
            Length Ncon array. Error propagation in the constraints,
            dc/dp*cov_p*dc/dp + (dc/dx*dx/dp)*cov_p*(dc/dx*dx/dp)
        - results.propagation_f: numpy.float64
            Real number. Error propagation in the objective function,
            df/dp*cov_p*df/dp + (df/dx*dx/dp)*cov_p*(df/dx*dx/dp)
        - results.col: list
            Size Nvar. List of variable names. Note that variables names
            includes both decision variable and uncertain parameter names.
            The order can be mixed.
        - results.row: list
            Size Ncon+1. List of constraints and objective function names

    Raises
    ------
    TypeError
        When tee entry is not Boolean
    TypeError
        When diagnostic_mode entry is not Boolean
    TypeError
        When solver_options entry is not None and a Dictionary
    Warnings
        When an element of theta_names includes a space

    """

    if not isinstance(tee, bool):
        raise TypeError("tee  must be boolean.")
    if not isinstance(diagnostic_mode, bool):
        raise TypeError("diagnostic_mode  must be boolean.")
    if not solver_options == None:
        if not isinstance(solver_options, dict):
            raise TypeError("solver_options must be dictionary.")
    if covariance_n is None:
        if isinstance(data, pd.DataFrame):
            covariance_n = len(data.index)
        else:
            covariance_n = len(data)
        logger.info(
            "covariance_n omitted from quantify_propagate_uncertainty().  "
            f"Assuming {covariance_n}"
        )
    # Remove all "'" and " " in theta_names
    theta_names, var_dic, variable_clean = clean_variable_name(theta_names)
    parmest_class = parmest.Estimator(
        model_function,
        data,
        theta_names,
        obj_function,
        tee,
        diagnostic_mode,
        solver_options,
    )
    obj, theta, cov = parmest_class.theta_est(calc_cov=True, cov_n=covariance_n)
    # Convert theta keys to the original name
    # Revert theta_names to be original
    if variable_clean:
        theta_out = {}
        for i in var_dic.keys():
            theta_out[var_dic[i]] = theta[i]
        theta_names = [var_dic[v] for v in theta_names]
    else:
        theta_out = theta
    propagate_results = propagate_uncertainty(
        model_uncertain, theta, cov, theta_names, tee
    )

    Output = namedtuple(
        "Output",
        [
            "obj",
            "theta",
            "theta_names",
            "cov",
            "gradient_f",
            "gradient_c",
            "dsdp",
            "propagation_c",
            "propagation_f",
            "col",
            "row",
        ],
    )
    results = Output(
        obj,
        theta_out,
        theta_names,
        cov,
        propagate_results.gradient_f,
        propagate_results.gradient_c,
        propagate_results.dsdp,
        propagate_results.propagation_c,
        propagate_results.propagation_f,
        propagate_results.col,
        propagate_results.row,
    )
    return results


def propagate_uncertainty(
    model_uncertain, theta, cov, theta_names, tee=False, solver_options=None
):

    """This function calculates gradient vector, expectation, and variance of
    the objective function and constraints  of the model for given estimated
    optimal parameters and covariance matrix of parameters.
    It calculates error propagation of the objective function and constraints
    by using gradient vector and covariance matrix.

    The following terms are used to define the output dimensions:
    Ncon   = number of constraints
    Nvar   = number of variables (Nx + Ntheta)
    Nx     = the number of decision (primal) variables
    Np = number of uncertain parameters.

    Parameters
    ----------
    model_uncertain : function or Pyomo ConcreteModel
        Function is a python/ Function that generates an instance of the
        Pyomo model
    theta : dict
        Size Ntheta python dictionary. Estimated parameters
    cov : numpy.ndarray
        Ntheta by Ntheta matrix. Covariance matrix of parameters
    theta_names : list of strings
        Size Ntheta. List of estimated l theta names
    tee : bool, optional
        Indicates that ef solver output should be teed, by default False
    solver_options : dict, optional
        Provides options to the solver (also the name of an attribute),
        by default None

    Returns
    -------
    tuple
        results object containing the all information including

        - results.gradient_f: numpy.ndarray
            Length Nvar array. Gradient vector of the objective function with
            respect to the (decision variables, parameters) at the optimal
            solution
        - results.gradient_c: scipy.sparse.csr.csr_matrix
            Ncon by Nvar size sparse matrix. Gradient vector of the constraints
            with respect to the (decision variables, parameters) at the optimal
            solution.
        - results.dsdp: scipy.sparse.csr.csr_matrix
            Ntheta by Nvar size sparse matrix. Gradient vector of the
            (decision variables, parameters) with respect to paramerters
            (=theta_name). number of rows = len(theta_name),
            number of columns= len(col)
        - results.propagation_c: numpy.ndarray
            Length Ncon array. Error propagation in the constraints,
            dc/dp*cov_p*dc/dp + (dc/dx*dx/dp)*cov_p*(dc/dx*dx/dp)
        - results.propagation_f: numpy.float64
            Real number. Error propagation in the objective function,
            df/dp*cov_p*df/dp + (df/dx*dx/dp)*cov_p*(df/dx*dx/dp)
        - results.col: list
            Size Nvar. List of variable names. Note that variables names
            includes both decision variable and uncertain parameter names.
            The order can be mixed.
        - results.row: list
            Size Ncon+1. List of constraints and objective function names

    Raises
    ------
    Exception
        if model_uncertain is neither 'ConcreteModel' nor 'function'.
    """
    # define Pyomo model
    try:
        if isinstance(model_uncertain, Block):
            model = model_uncertain
        else:
            model = model_uncertain()
    except TypeError as e:
        raise """
        model_uncertain must be either python function or Pyomo ConcreteModel"""

    # check and convert (optional) covariance matrix
    if isinstance(cov, np.ndarray):
        cov_ = cov
    else:
        cov_ = cov.to_numpy()

    if len(cov_.shape) != 2:
        raise ValueError("cov must be a 2-dimensional matrix or dataframe")

    if cov_.shape[0] != cov_.shape[1]:
        raise ValueError("cov must be square")

    if cov_.shape[0] != len(theta_names):
        raise ValueError(
            """cov must be a n x n matrix or dataframe where 
                         n = len(theta_names)"""
        )

    if len(theta_names) != len(theta):
        raise ValueError(
            """theta_names and theta must have the same number 
                          of elements"""
        )

    # Remove all "'" in theta_names
    theta_names, var_dic, variable_clean = clean_variable_name(theta_names)
    for v in theta_names:
        model.find_component(var_dic[v]).setlb(theta[v])
        model.find_component(var_dic[v]).setub(theta[v])
    # get gradient of the objective function, constraints,
    # and the column,row names
    dsdp, col = get_dsdp(model, theta_names, theta, var_dic, tee)
    dsdp = dsdp.toarray().T  # change shape, Nvar by Ntheta
    gradient_f, gradient_c, col, row, line_dic = get_dfds_dcds(model, theta_names, tee)
    num_constraints = len(
        list(model.component_data_objects(Constraint, active=True, descend_into=True))
    )

    # calculate error propagation of the objective fuction
    # = df/dp*cov_p*df/dp + (df/dx*dx/dp)*cov_p*(df/dx*dx/dp)
    # = (df/ds*ds/dp)*cov*(df/ds*ds/dp)
    # step 1. df/ds*ds/dp

    # [1 x (Nx + Np)] matrix
    df_ds = np.reshape(gradient_f, (1, len(gradient_f)))

    # [(Nx + Np) x (Np) ] matrix
    ds_dp = np.reshape(dsdp, (len(dsdp), len(theta_names)))

    # [1 x Np ] matrix
    fssp = np.matmul(df_ds, ds_dp)

    # step 2. (df/ds*ds/dp)*cov*(df/ds*ds/dp).transpose()
    # [1 x Np] x [Np x Np ] x [Np x 1] = [1 x 1]
    propagation_f = fssp @ cov_ @ fssp.transpose()
    #
    assert propagation_f.shape == (
        1,
        1,
    ), """propagation_f should be 
    a 1 x 1 matrix. Something is wrong if you are seeing this..."""

    # convert to scalar
    propagation_f = propagation_f[0, 0]

    # calculate error propagation of constraints
    # = dc/dp*cov_p*dc/dp + (dc/dx*dx/dp)*cov_p*(dc/dx*dx/dp)
    # = (dc/ds*ds/dp)*cov*(dc/ds*ds/dp)
    num_constraints = len(
        list(model.component_data_objects(Constraint, active=True, descend_into=True))
    )
    if num_constraints > 0:
        # gradient_c rearrange.
        # k_aug sparse form is [col_idx, row_idx, val] with index starts from 1
        # python sparse form is [row_idx, col_idx, val] with index starts from 0
        # note: vairable 'row' from get_dfds_dcds includes objecive function
        # name. i.e, Ncon = len(row)-1
        """
        row_idx = gradient_c[:,1]-1
        col_idx = gradient_c[:,0]-1
        data = gradient_c[:,2]
        gradient_c
        = sparse.csr_matrix((data, (row_idx, col_idx)),
                             shape=(len(row)-1, len(col)))
        """
        # step 1. dc/ds*ds/dp
        # [Ncon x (Nx + Np)] x [(Nx + Np] x Np] = [Ncon x Np] matrix
        cssp = np.matmul(gradient_c.toarray(), dsdp)

        # step 2. (dc/ds*ds/dp)*cov*(dc/ds*ds/dp).transpose()
        # [Ncon x Np] x [Np x Np] x [Np x Ncon]

        # Not sure if this is correct.
        # propagation_c = np.sum(np.matmul(cssp,cov)*cssp,axis=1)

        # Updated to match propgation_f
        propagation_c = cssp @ cov_ @ cssp.transpose()
    else:
        propagation_c = np.array([])

    Output = namedtuple(
        "Output",
        [
            "gradient_f",
            "gradient_c",
            "dsdp",
            "propagation_c",
            "propagation_f",
            "col",
            "row",
        ],
    )
    results = Output(
        gradient_f,
        gradient_c,
        sparse.csr_matrix(dsdp),
        propagation_c,
        propagation_f,
        col,
        row,
    )
    return results


# TODO: Improve the robustness of Parmest then remove this function.
def clean_variable_name(theta_names):
    """This function removes all ' and spaces in theta_names. Note that
    the current theta_est(calc_cov=True) of parmest in Pyomo doesn't allow ' and
    spaces in the variable names. Once a future version of Parmest fixes this issue,
    this function can be depreciated.

    Parameters
    ----------
    theta_names : list of strings
        List of Var names

    Returns
    -------
        It returns the following variables

        - theta_names_out: list of strings
            List of Var names after removing  all ' and spaces
        - var_dic: dict
            Dictionary with keys converted theta_names and values origianl
            theta_names
    """
    # Variable names cannot have "'" for parmest_class.theta_est(calc_cov=True)
    # Save original variables name in to var_dic
    # Remove all "'" and " " in theta_names
    var_dic = {}
    theta_names_out = []
    clean = False
    for i in range(len(theta_names)):
        theta_tmp = theta_names[i].replace("'", "")
        theta_tmp = theta_tmp.replace(" ", "")
        theta_names_out.append(theta_tmp)
        var_dic[theta_tmp] = theta_names[i]
        if "'" in theta_names[i] or " " in theta_names[i]:
            logger.warning(theta_names[i] + " includes ' or space.")
            logger.warning("The cleaned name: " + theta_tmp)
            clean = True
    if clean:
        logger.warning("All ' and spaces in theta_names are removed.")
    return theta_names_out, var_dic, clean
