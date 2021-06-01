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
"""
This module contains model diagnostic utility functions for use in IDAES (Pyomo) models.
"""

import pytest

# Need to update
import pyomo.environ as pyo
                           
# TODO: Add pyomo.dae test case
'''
from pyomo.environ import TransformationFactory
from pyomo.dae import ContinuousSet, DerivativeVar
'''

# Need to update
from idaes.core.util.model_diagnostics import *

# Author: Alex Dowling

# This was from 
# @pytest.fixture()
def problem1():
    m = pyo.ConcreteModel()

    m.I = pyo.Set(initialize=[i for i in range(5)])

    m.x = pyo.Var(m.I,bounds=(-10,10),initialize=1.0)

    m.con1 = pyo.Constraint(expr=m.x[0] + m.x[1] - m.x[3] >= 10)
    m.con2 = pyo.Constraint(expr=m.x[0]*m.x[3] + m.x[1] >= 0)
    m.con3 = pyo.Constraint(expr=m.x[4]*m.x[3] + m.x[0]*m.x[3] - m.x[4] == 0)

    m.obj = pyo.Objective(expr=sum(m.x[i]**2 for i in m.I))

    return m

def example2(with_degenerate_constraint=True):
    ''' Create the Pyomo model for Example 2
    
    Arguments:
        with_degenerate_constraint: Boolean, if True, include the redundant linear constraint
    
    Returns:
        m2: Pyomo model
    '''
    
    m2 = pyo.ConcreteModel()

    m2.I = pyo.Set(initialize=[i for i in range(1,4)])

    m2.x = pyo.Var(m2.I,bounds=(0,5),initialize=1.0)

    m2.con1 = pyo.Constraint(expr=m2.x[1] + m2.x[2] >= 1)
    m2.con2 = pyo.Constraint(expr=m2.x[1] + m2.x[2] + m2.x[3] == 1)
    m2.con3 = pyo.Constraint(expr=m2.x[2] - 2*m2.x[3] <= 1)
    m2.con4 = pyo.Constraint(expr=m2.x[1] + m2.x[3] >= 1)
    
    if with_degenerate_constraint:
        m2.con5 = pyo.Constraint(expr=m2.x[1] + m2.x[2] + m2.x[3] == 1)

    m2.obj = pyo.Objective(expr=sum(m2.x[i] for i in m2.I))
    
    return m2

def extract_constraint_names(cs):
    ''' Get constraint names from ComponentSet
    
    Arguments:
        cs: ComponentSet object
        
    Return:
        constraint_names: list of constraint names (strings)
            
    '''
    
    constraint_names = []
    for i in cs:
        constraint_names.append(i.name)
    return constraint_names

# Problem 1
@pytest.mark.skipif(not pyo.SolverFactory('ipopt').available(False), reason="no Ipopt")
@pytest.mark.unit
def test_problem1():
    # Create test problem
    m = problem1()
    
    # Specify Ipopt as the solver
    opt = pyo.SolverFactory('ipopt')

    # Specifying an iteration limit of 0 allows us to inspect the initial point
    opt.options['max_iter'] = 0

    # "Solving" the model with an iteration limit of 0 load the initial point and applies
    # any preprocessors (e.g., enforces bounds)
    opt.solve(m, tee=True)

    # Create Degeneracy Hunter object
    dh = DegeneracyHunter(m) 
    
    # Find constraints with residuals > 0.1 
    initial_point_constraints = dh.check_residuals(tol=0.1)
    
    # Check there are 2 constraints with large residuals
    assert len(initial_point_constraints) == 2
    
    initial_point_constraint_names = extract_constraint_names(initial_point_constraints)
    
    # Check first constraint
    assert initial_point_constraint_names[0] == 'con1'
    
    # Check second constraint
    assert initial_point_constraint_names[1] == 'con3'
    
    opt.options['max_iter'] = 50
    
    # Solve
    opt.solve(m, tee=True)
    
    # Find constraints with residuals > 0.1 
    solution_constraints = dh.check_residuals(tol=1E-6)
    
    # Check at the solution no constraints are violated
    assert len(solution_constraints) == 0
    
    # Check no constraints are near their bounds
    solution_bounds = dh.check_variable_bounds(tol=0.1)
    
    # Check at the solution no constraints are violated
    assert len(solution_bounds) == 0
    

# Problem 2 without degenerate constraint
@pytest.mark.skipif(not pyo.SolverFactory('ipopt').available(False), reason="no Ipopt")
@pytest.mark.unit
def test_problem2_without_degenerate_constraint():

    # Create test problem instance
    m2 = example2(with_degenerate_constraint=False)
    
    # Specify Ipopt as the solver
    opt = pyo.SolverFactory('ipopt')
    
    # Specifying an iteration limit of 0 allows us to inspect the initial point
    opt.options['max_iter'] = 0

    # "Solving" the model with an iteration limit of 0 load the initial point and applies
    # any preprocessors (e.g., enforces bounds)
    opt.solve(m2, tee=True)
    
    # Create Degeneracy Hunter object
    dh2 = DegeneracyHunter(m2)
    
    # Check for violated constraints at the initial point
    initial_point_constraints = dh2.check_residuals(tol=0.1)
    
    # Check there are 1 constraints with large residuals
    assert len(initial_point_constraints) == 1
    
    initial_point_constraint_names = extract_constraint_names(initial_point_constraints)
    
    # Check first constraint
    assert initial_point_constraint_names[0] == 'con2'
    
    # Resolve
    opt.options['max_iter'] = 500
    opt.solve(m2, tee=True)
    
    # Check solution
    x_sln = []
    
    for i in m2.I:
        x_sln.append(m2.x[i]())
    
    assert pytest.approx(x_sln[0], abs=1E-6) == 1.0
    assert pytest.approx(x_sln[1], abs=1E-6) == 0.0
    assert pytest.approx(x_sln[2], abs=1E-6) == 0.0
    



# Problem 2 with degenerate constraint
@pytest.mark.skipif(not pyo.SolverFactory('ipopt').available(False), reason="no Ipopt")
@pytest.mark.unit
def test_problem2_with_degenerate_constraint():

    # Create test problem instance
    m2 = example2(with_degenerate_constraint=True)
    
    # Specify Ipopt as the solver
    opt = pyo.SolverFactory('ipopt')

    # Specifying an iteration limit of 0 allows us to inspect the initial point
    opt.options['max_iter'] = 0

    # "Solving" the model with an iteration limit of 0 load the initial point and applies
    # any preprocessors (e.g., enforces bounds)
    opt.solve(m2, tee=True)
    
    # Create Degeneracy Hunter object
    dh2 = DegeneracyHunter(m2)
    
    # Check for violated constraints at the initial point
    initial_point_constraints = dh2.check_residuals(tol=0.1)
    
    # Check there are 2 constraints with large residuals
    assert len(initial_point_constraints) == 2
    
    initial_point_constraint_names = extract_constraint_names(initial_point_constraints)
    
    # Check first constraint
    assert initial_point_constraint_names[0] == 'con2'
    
    # Check first constraint
    assert initial_point_constraint_names[1] == 'con5'
    
    # Resolve
    opt.options['max_iter'] = 500
    opt.solve(m2, tee=True)
    
    # Check solution
    x_sln = []
    
    for i in m2.I:
        x_sln.append(m2.x[i]())
    
    assert pytest.approx(x_sln[0], abs=1E-6) == 1.0
    assert pytest.approx(x_sln[1], abs=1E-6) == 0.0
    assert pytest.approx(x_sln[2], abs=1E-6) == 0.0
    
    # Check the rank
    n_rank_deficient = dh2.check_rank_equality_constraints()
    
    assert n_rank_deficient == 1
    
    # TODO: Add MILP solver to idaes get-extensions and add more tests