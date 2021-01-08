''' Degeneracy Hunter is a collection of utility functions to assist in mathematical
modeling in Pyomo.

Alex Dowling, University of Notre Dame

'''

from pyomo.environ import *
from pyomo.core.kernel.component_set import ComponentSet
from pyomo.core.expr.visitor import identify_variables
from pyomo.contrib.pynumero.interfaces.pyomo_nlp import PyomoNLP
import numpy as np
from scipy.sparse.linalg import svds
from scipy.sparse import issparse, find

from pyomo.opt import SolverStatus, TerminationCondition

import matplotlib.pyplot as plt


# This library is already used in Pyomo
import operator

'''
This code is from Robby Parker.
'''

def get_con_eq_idx_map(interface, constraint_names):
    # This is super janky. TODO: Add functionality
    # to PyNumero.
    # constraint_names is an argument because generating
    # constraint_names list from PyNumero takes forever.
    # TODO: Find way to make AslNLP::constraint_names()
    # not ridiculously slow.
    mask = interface._con_full_eq_mask
    temp = np.array(constraint_names)
    eq_names = list(temp[mask])
    return {name: i for i, name in enumerate(eq_names)}

def get_eq_con_list(interface, constraint_names):
    mask = interface._con_full_eq_mask
    temp = np.array(constraint_names)
    eq_names = list(temp[mask])
    return eq_names

class DegeneracyHunter():

    def __init__(self,block_or_jac,solver=None):
        ''' Initialize Degeneracy Hunter Object
    
        Arguments:
            block_or_jac: Pyomo model or Jacobian
        
        '''
        
        if type(block_or_jac) is ConcreteModel:
        
            # Add Pyomo model to the object
            self.block = block_or_jac
        
            # setup pynumero interface
            self.nlp = PyomoNLP(self.block)

            # calculate Jacobian of equality constraints in COO sparse matrix format
            jac_eq = self.nlp.evaluate_jacobian_eq()

            # save the Jacobian
            self.jac_eq = jac_eq
        
            # Create a list of equality constraint names
            pyomo_constraints = self.nlp.get_pyomo_constraints()
            constraint_names = [c.name for c in pyomo_constraints]
            self.name2eq_idx = get_con_eq_idx_map(self.nlp, constraint_names)
            self.eq_con_list = get_eq_con_list(self.nlp, constraint_names)
        
        else:
            self.jac_eq = block_or_jac
            
            self.eq_con_list = None
        
        # number of equality constraints, variables
        self.n_eq = self.jac_eq.shape[0]
        self.n_var = self.jac_eq.shape[1]
        
        
        # Define default candidate equations (enumerate)
        candidate_eqns = range(self.n_eq)
                
        # Initialize solver
        if solver is None:
            self.solver = SolverFactory('gurobi')
            self.solver.options = {'NumericFocus':3}
            
        else:
            self.solver = solver
            
        # Create spot to store singular values
        self.s = None
        
        # Set constants for MILPs
        self.max_nu = 1E5
        self.min_nonzero_nu = 1E-5
        
        
    # TODO: migrate these improves to model_statistics.large_residuals_set
    def check_residuals(self,tol=1e-5, print_level=2, sort=True):
        """
        Method to return a ComponentSet of all Constraint components with a
        residual greater than a given threshold which appear in a model.
        Args:
            block : model to be studied
            tol : residual threshold for inclusion in ComponentSet
            print_level: controls to extend of output to the screen
                0: nothing printed
                1: only name of constraint printed
                2: each constraint is pretty printed
                3: pretty print each constraint, then print value for included variable
            sort: sort residuals in descending order for printing
        Returns:
            A ComponentSet including all Constraint components with a residual
            greater than tol which appear in block
        """
        large_residuals = dict()
        large_residuals_set = ComponentSet()
        for c in self.block.component_data_objects(
                ctype=Constraint, active=True, descend_into=True):
            r = 0.0 # residual
        
            # check the lower bound
            # skip if inequality constraint
            if c.lower is None:
                r_temp = 0
            else:
                r_temp = value(c.lower - c.body())
            # update the residual
            if c.active and r_temp > r:
                r = r_temp
        
            # check the upper bound
            # skip if inequality constraint
            if c.upper is None:
                r_temp = 0
            else:
                r_temp = value(c.body() - c.upper)

            # update the residual
            if c.active and r_temp > r:
                r = r_temp
            
            # save residual if it is above threshold
            if r > tol:
                large_residuals_set.add(c)
            
                # add to dictionary for sorting
                # key: constraint, value: residual
                if print_level > 0:
                    large_residuals[c] = r
    
        if print_level > 0:
            print(" ")
            if len(large_residuals) > 0:
                print("All constraints with residuals larger than",tol,":")
                if print_level == 1:
                    print("Count\tName\t|residual|")
                
                if sort:
                    large_residuals = dict(sorted(large_residuals.items(), key=operator.itemgetter(1),reverse=True))
        
                for i, (c,r) in enumerate(large_residuals.items()):
                    if print_level == 0:
                        # Basic print statement. count, constraint, residual
                        print(i,"\t",c,"\t",r)
                    else:
                        # Pretty print constraint
                        print("\ncount =",i,"\t|residual| =",r)
                        c.pprint()
                
                    if print_level == 2:
                        # print values and bounds for each variable in the constraint
                        print("variable\tlower\tvalue\tupper")
                        for v in identify_variables(c.body):
                            self.print_variable_bounds(v)
            else:
                print("No constraints with residuals larger than",tol,"!")
                
        return large_residuals_set

    # Migrate these improvements to model_statistics.variables_near_bounds_generator
    def check_variable_bounds(self,tol=1e-5, skip_lb = False, skip_ub = False, LOUD=True):
        """
        Method to return a ComponentSet of all variables within a tolerance
        of their bounds.
        Args:
            block : model to be studied
            tol : residual threshold for inclusion in ComponentSet
        Returns:
            A ComponentSet including all Constraint components with a residual
            greater than tol which appear in block
        """
        variables_near_bounds_set = ComponentSet()
        for c in self.block.component_data_objects(
                ctype=Var):
        
            added = False
        
            if not c.lb is None and not skip_lb:
                if c.value - c.lb < tol:
                    variables_near_bounds_set.add(c)
                    added = True
                
            if not c.ub is None and not skip_ub and not added:
                if c.ub - c.value < tol:
                    variables_near_bounds_set.add(c)
    
        if LOUD:
            print(" ")
            if len(variables_near_bounds_set) > 0:
                print("Variables within",tol,"of their bounds:")
                print("variable\tlower\tvalue\tupper")
                for v in variables_near_bounds_set:
                    self.print_variable_bounds(v)
            else:
                print("No variables within",tol,"of their bounds.")
            
        return variables_near_bounds_set
    
    def check_rank_equality_constraints(self):
        """
        Method to check the rank of the Jacobian of the equality constraints
        """
        
        print("\nChecking rank of Jacobian of equality constraints...")
        
        if self.s is None:
            self.svd_analysis()
        
        n = len(self.s)
        
        print("Smallest singular value(s):")
        for i in range(n):
            print("%.3E" % self.s[i])

    @staticmethod
    def _prepare_ids_milp(jac_eq, M=1E5):
        '''
        Prepare MILP to compute the irreducible degenerate set
        
        Argument:
            jac_eq Jacobian of equality constraints [matrix]
            M: largest value for nu
            
        Returns:
            m_dh: Pyomo model to calculate irreducible degenerate sets
        '''
    
        n_eq = jac_eq.shape[0]
        n_var = jac_eq.shape[1]
    
        # Create Pyomo model for irreducible degenerate set
        m_dh = ConcreteModel()

        # Create index for constraints
        m_dh.C = Set(initialize=[i for i in range(n_eq)])

        m_dh.V = Set(initialize=[i for i in range(n_var)])

        # Add variables
        m_dh.nu = Var(m_dh.C, bounds=(-M,M),initialize=1.0)
        m_dh.y = Var(m_dh.C, domain=Binary)

        # Constraint to enforce set is degenerate
        if issparse(jac_eq):
            m_dh.J = jac_eq.tocsc()
            
            def eq_degenerate(m_dh, v):
            
                # Find the columns with non-zero entries
                C_ = find(m_dh.J[:,v])[0]
                return sum(m_dh.J[c,v] * m_dh.nu[c] for c in C_) == 0
            
        else:
            m_dh.J = jac_eq
        
            def eq_degenerate(m_dh, v):
                return sum(m_dh.J[c,v] * m_dh.nu[c] for c in m_dh.C) == 0
            
        m_dh.degenerate = Constraint(m_dh.V, rule=eq_degenerate)

        def eq_lower(m_dh, c):
            return -M * m_dh.y[c] <= m_dh.nu[c]
        m_dh.lower = Constraint(m_dh.C, rule=eq_lower)

        def eq_upper(m_dh, c):
            return  m_dh.nu[c] <= M * m_dh.y[c]
        m_dh.upper= Constraint(m_dh.C, rule=eq_upper)

        m_dh.obj = Objective(expr=sum(m_dh.y[c] for c in m_dh.C))
        
        return m_dh
    
    @staticmethod
    def _prepare_find_candidates_milp(jac_eq, M = 1E5, m_small = 1E-5):
        '''
        Prepare MILP to find candidate equations for consider for IDS
        
        Argument:
            jac_eq Jacobian of equality constraints [matrix]
            M: maximum value for nu
            m_small: smallest value for nu to be considered non-zero
            
        Returns:
            m_fc: Pyomo model to find candidates
        '''
    
        n_eq = jac_eq.shape[0]
        n_var = jac_eq.shape[1]
    
        # Create Pyomo model for irreducible degenerate set
        m_dh = ConcreteModel()

        # Create index for constraints
        m_dh.C = Set(initialize=[i for i in range(n_eq)])

        m_dh.V = Set(initialize=[i for i in range(n_var)])

        # Specify minimum size for nu to be considered non-zero
        m_dh.m_small = m_small

        # Add variables
        m_dh.nu = Var(m_dh.C, bounds=(-M-m_small,M+m_small),initialize=1.0)
        m_dh.y_pos = Var(m_dh.C, domain=Binary)
        m_dh.y_neg = Var(m_dh.C, domain=Binary)
        m_dh.abs_nu = Var(m_dh.C, bounds=(0, M + m_small))
        
        # Positive exclusive or negative
        def eq_pos_xor_negative(m, c):
            return m.y_pos[c] + m.y_neg[c] <= 1
        
        m_dh.pos_xor_neg = Constraint(m_dh.C)

        # Constraint to enforce set is degenerate
        if issparse(jac_eq):
            m_dh.J = jac_eq.tocsc()
            
            def eq_degenerate(m_dh, v):
            
                # Find the columns with non-zero entries
                C_ = find(m_dh.J[:,v])[0]
                return sum(m_dh.J[c,v] * m_dh.nu[c] for c in C_) == 0
            
        else:
            m_dh.J = jac_eq
        
            def eq_degenerate(m_dh, v):
                return sum(m_dh.J[c,v] * m_dh.nu[c] for c in m_dh.C) == 0
            
        m_dh.degenerate = Constraint(m_dh.V, rule=eq_degenerate)

        # When y_pos = 1, nu >= m_small
        # When y_pos = 0, nu >= - m_small
        def eq_pos_lower(m_dh, c):
            return m_dh.nu[c] >= -m_small + 2*m_small * m_dh.y_pos[c]
        m_dh.pos_lower = Constraint(m_dh.C, rule=eq_pos_lower)

        # When y_pos = 1, nu <= M + m_small
        # When y_pos = 0, nu <= m_small
        def eq_pos_upper(m_dh, c):
            return m_dh.nu[c] <= M * m_dh.y_pos[c] + m_small
        m_dh.pos_upper = Constraint(m_dh.C, rule=eq_pos_upper)
        
        # When y_neg = 1, nu <= -m_small
        # When y_neg = 0, nu <= m_small
        def eq_neg_upper(m_dh, c):
            return m_dh.nu[c] <= m_small - 2*m_small * m_dh.y_neg[c]
        m_dh.neg_upper = Constraint(m_dh.C, rule=eq_neg_upper)

        # When y_neg = 1, nu >= -M - m_small
        # When y_neg = 0, nu >= - m_small
        def eq_neg_lower(m_dh, c):
            return m_dh.nu[c] >= -M * m_dh.y_neg[c] - m_small
        m_dh.neg_lower = Constraint(m_dh.C, rule=eq_neg_lower)
        
        # Absolute value
        def eq_abs_lower(m_dh,c):
            return -m_dh.abs_nu[c] <= m_dh.nu[c]
        m_dh.abs_lower = Constraint(m_dh.C, rule=eq_abs_lower)
        
        def eq_abs_upper(m_dh,c):
            return m_dh.nu[c] <= m_dh.abs_nu[c]
        m_dh.abs_upper = Constraint(m_dh.C, rule=eq_abs_upper)
        
        # At least one constraint must be in the degenerate set
        m_dh.degenerate_set_nonempty = Constraint(expr=sum(m_dh.y_pos[c] + m_dh.y_neg[c] for c in m_dh.C) >= 1)

        # Minimize the L1-norm of nu
        m_dh.obj = Objective(expr=sum(m_dh.abs_nu[c] for c in m_dh.C))
        
        return m_dh

    
    @staticmethod
    def _check_candidate_ids(ids_milp,solver,c,eq_con_list=None,tee=False):
        ''' Solve MILP to check if equation 'c' is a main component in an irreducible
            degenerate set
            
            Arguments:
                ids_milp: Pyomo model to calculate IDS
                solver: Pyomo solver (must support MILP)
                c: index for the constraint to consider [integer]

            Returns:
                ids: either None or dictionary containing the IDS
        '''
        
        # Fix weight on candidate equation
        ids_milp.nu[c].fix(1.0)
        
        # Solve MILP
        results = solver.solve(ids_milp, tee=tee)

        ids_milp.nu[c].unfix()

        if results.solver.Status == SolverStatus.ok and results.solver.termination_condition == TerminationCondition.optimal:
            # We found an irreducible degenerate set
            
            # Create empty dictionary
            ids_ = {}
            
            for i in ids_milp.C:
                # Check if constraint is included
                if ids_milp.y[i]() > 0.9:
                    # If it is, save the value of nu
                    if eq_con_list is None:
                        name = i
                    else:
                        name = eq_con_list[i]
                    ids_[name] = ids_milp.nu[i]()
            return ids_
        else:
            return None
    
    @staticmethod
    def _find_candidate_eqs(candidates_milp, solver, eq_con_list=None,tee=False):
        ''' Solve MILP to check if equation 'c' is a main component in an irreducible
            degenerate set
            
            Arguments:
                candidates_milp: Pyomo model to calculate IDS
                solver: Pyomo solver (must support MILP)
                c: index for the constraint to consider [integer]

            Returns:
                candidate_eqns: either None or list of indicies
                degenerate_set: either None or dictionary containing the degenerate_set
        '''
        
        results = solver.solve(candidates_milp, tee=tee)
        
        if results.solver.Status == SolverStatus.ok and results.solver.termination_condition == TerminationCondition.optimal:
            # We found a degenerate set
            
            # Create empty dictionary
            ds_ = {}
            
            # Create empty list
            candidate_eqns = []
            
            for i in candidates_milp.C:
                # Check if constraint is included
                if candidates_milp.abs_nu[i]() > candidates_milp.m_small * 0.99:
                    # If it is, save the value of nu
                    if eq_con_list is None:
                        name = i
                    else:
                        name = eq_con_list[i]
                    ds_[name] = candidates_milp.nu[i]()
                    candidate_eqns.append(i)
                    
            return candidate_eqns, ds_
        else:
            return None, None
        
    
    def svd_analysis(self,n_smallest_sv = 10):
        '''
        Perform SVD analysis of the constraint Jacobian
        '''
        
        # Determine the number of singular values to compute
        n_sv = min(n_smallest_sv, min(self.n_eq, self.n_var) - 1)
        
        # Perform SVD
        # Recall J is a n_eq x n_var matrix
        # Thus U is a n_eq x n_eq matrix
        # And V is a n_var x n_var
        # (U or V may be smaller in economy mode)
        # Thus we really only care about U
        u, s, v = svds(self.jac_eq, k = n_sv, which='SM')
        
        # Save results
        self.u = u
        self.s = s
        self.v = v
    
    def find_candidate_equations(self,verbose=True,tee=False):
        '''
        Solve MILP to find a degenerate set and candidate equations
        '''
        
        if verbose:
            print("*** Searching for a Single Degenerate Set ***")
            print("Building MILP model...")
        self.candidates_milp = self._prepare_find_candidates_milp(self.jac_eq,
                                                                    self.max_nu,
                                                                    self.min_nonzero_nu)
        
        if verbose:
            print("Solving MILP model...")
        ce, ds = self._find_candidate_eqs(self.candidates_milp,
                                            self.solver,
                                            self.eq_con_list,
                                            tee)
        
        if ce is not None:
            self.candidate_eqns = ce
        
        return ds
        
    
    def find_irreducible_degenerate_sets(self,verbose=True,tee=False):
        """
        Compute irreducible degenerate sets
        """
        
        irreducible_degenerate_sets = []
        
        if verbose:
            print("*** Searching for Irreducible Degenerate Sets ***")
            print("Building MILP model...")
        self.dh_milp = self._prepare_ids_milp(self.jac_eq, self.max_nu)
                
        # Loop over candidate equations
        for c in self.candidate_eqns:
        
            if verbose:
                print("Solving MILP",c+1,"of",len(self.candidate_eqns),"...")
        
            # Check if equation 'c' is a major element of an IDS
            ids_ = self._check_candidate_ids(self.dh_milp,
                                                self.solver,
                                                c,
                                                self.eq_con_list,
                                                tee)
            
            if ids_ is not None:
                irreducible_degenerate_sets.append(ids_)

        if verbose:
            for i,s in enumerate(irreducible_degenerate_sets):
                print("\nIrreducible Degenerate Set",i)
                print("nu\tConstraint Name")
                for k,v in s.items():
                    print(v,"\t",k)
        
        return irreducible_degenerate_sets

    ### Helper Functions
    
    @staticmethod
    def print_variable_bounds(v):
        ''' Print variable, bounds, and value
        Argument:
            v: variable
        Return:
            nothing
        '''
        print(v,"\t\t",v.lb,"\t",v.value,"\t",v.ub)