##Fit the constant values of each tree using regression

import numpy as np
import random
import re
import os
from IPython.display import FileLink
from IPython.display import display as disp

import math

from sklearn.metrics import mean_squared_error, r2_score

from pyomo.environ import *
from pyomo.core import value
from pyomo.common.errors import ApplicationError
from pyomo.opt import TerminationCondition

np.random.seed(42)

def constant_fit(es, cst_count, X, y, opt, tree):

    c_lo = min(-100, np.min(X))
    c_up = max(100, np.max(X))

    if cst_count > 0:
        
        cm = ConcreteModel()
        
        for num in range(cst_count):
            es = re.sub("cst", f"cm.c[{num}]", es, count=1)
        cm.c = Var([i for i in range(cst_count)], domain=Reals)
        # cm.c = Var([i for i in range(cst_count)], domain=Reals, bounds=(c_lo, c_up))
        for i in range(cst_count):
            cm.c[i] = 1 + random.random()
        n_data_set = [i for i in range(len(y))]
        cm.y_pred = Var(n_data_set, domain=Reals)
        # cm.y_pred = Var(n_data_set, domain=Reals, initialize=np.mean(y))

        pyomo_es = es.replace(":","i")
        pyomo_es_eval = [0]*len(y)

        # import pickle
        # with open("pyomo_es_eval.pkl", "wb") as f:
        #     pickle.dump(pyomo_es_eval, f)

        for i in range(len(y)):
            try:
                pyomo_es_eval[i] = eval(pyomo_es)
            except (ZeroDivisionError, ValueError, TypeError):
                print('Tree is not valid')
                return None, 999999, 0
            
            if pyomo_es_eval is None or 'inf' in str(pyomo_es_eval[i]) or 'nan' in str(pyomo_es_eval[i]):
                print('Tree is not valid')
                return None, 999999, 0
            if '/(0)' in str(pyomo_es_eval[i]) or '/(0*' in str(pyomo_es_eval[i]) or '/(0/' in str(pyomo_es_eval[i]):
                print('Tree is not valid')
                return None, 999999, 0
            # if re.search(r'/(0)', str(pyomo_es_eval[i])) or re.search(r'/\(0\*', str(pyomo_es_eval[i])) or re.search(r'/\(0/', str(pyomo_es_eval[i])):
            #     return None, 999999, 0
            # if re.search(r'/(0)', str(pyomo_es_eval[i])):
            #     return None, 999999, 0
            # if re.search(r'/\(0/', str(pyomo_es_eval[i])):
            #     return None, 999999, 0
            # if re.search(r'/\(0\*', str(pyomo_es_eval[i])):
            #     return None, 999999, 0

        cm.y_pred_constr = ConstraintList()
        sumup = 0
        for i in range(len(y)):
            cm.y_pred_constr.add(cm.y_pred[i] == pyomo_es_eval[i])
            sumup += (y[i] - cm.y_pred[i]) ** 2

        cm.obj = Objective(expr=(sumup/len(y))**0.5)
        try:
            results = opt.solve(cm, tee=True, keepfiles=True, symbolic_solver_labels=True)
            print('Running BARON to optimize the tree expression')
            
            # The solver output file path
            solver_output_file = f'summary_tree_{tree}'
            disp(FileLink(solver_output_file))

            if results.solver.status == SolverStatus.ok:
                print('Solver Status:', results.solver.status)   
            else:
                print('Solver Status:', results.solver.status)

                # Open the file
                with open(solver_output_file, 'r') as file:
                    content = file.read()
                # Print the contents on the screen
                print(content)
            print('Termination Condition:', results.solver.termination_condition)

            if results.solver.termination_condition != TerminationCondition.optimal:
                print('Rerunning BARON with IPOPT allowable')
                opt.options['MaxIter'] = 10  
                opt.options['AllowIPOPT'] = 1
                results = opt.solve(cm, tee=True, keepfiles=True, symbolic_solver_labels=True)
                opt.options['AllowIPOPT'] = 0

                if results.solver.status == SolverStatus.ok:
                    print('Solver Status:', results.solver.status)   
                else:
                    print('Solver Status:', results.solver.status)

                    # Open the file
                    with open(solver_output_file, 'r') as file:
                        content = file.read()
                    # Print the contents on the screen
                    print(content)

                    # Open the file
                    #os.startfile(solver_output_file) 
                print('Termination Condition:', results.solver.termination_condition)

        except ApplicationError:
            print('Tree is not valid')
            return None, 999999, 0        
        # try:
        #     results = opt.solve(cm, tee=True, keepfiles=True, symbolic_solver_labels=True, timeout=600)
        #     if results.solver.termination_condition == TerminationCondition.maxTimeLimit:
        #         return None, 999999, 0
        # except ApplicationError:
        #     return None, 999999, 0
        

        try:
            es2 = es[:]            
            for i in range(cst_count):
                es2 = es2.replace("cm.c[{}]".format(i), str(value(cm.c[i])))            
            obj = value(cm.obj)
            for i in range(len(y)):
                y_pred = eval(es2)
            r_squared = r2_score(y, y_pred)
            print('Tree is valid')
        except: # If the expression is not valid, just
            es2 = es[:]
            obj = 999999 
            r_squared = 0
            print('Tree is not valid')

    else:
        es2 = es[:]
        try:
            obj = EvaluateExpressionString(es2, X, y)[1]
            r_squared = EvaluateExpressionString(es2, X, y)[2]
            print('Tree does not contain a constant')
        except: # If the expression is not valid, just
            obj = 999999
            r_squared = 0
            print('Tree is not valid')
            pass

    #To regularize RMSE
    max_min = np.max(y) - np.min(y)
    obj = obj / max_min

    return es2, obj, r_squared


def EvaluateExpressionString(expression_string, Xt, yt):
    try:
        X = Xt
        y_pred = eval(expression_string)
    except:
        y_pred = []
        for i in range(len(yt)):
            X = np.reshape(Xt[i, :], (1, Xt.shape[1]))
            y_pred.append(eval(expression_string))
        y_pred = np.array(y_pred)

    RMSE = mean_squared_error(yt, y_pred) ** 0.5
    r_squared = r2_score(yt, y_pred)

    return y_pred, RMSE, r_squared