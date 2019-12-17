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
"""
Script to test functionality of key components for use in workshops.

Created on Mon Apr  1 13:32:37 2019

@author: alee
"""
check_count=0
# Check imports
try:
    from pyomo.environ import *
    from pyomo.opt import SolverStatus, TerminationCondition
    from pyomo.network import Arc, SequentialDecomposition
    print("Pyomo Import Checks:        Passed")
    check_count += 1
except:
    print("Pyomo Import Checks:        FAILED")

try:
    from idaes.core import *
    from idaes.unit_models import (PressureChanger,
                               CSTR,
                               Flash,
                               Heater,
                               Mixer,
                               Separator)
    from idaes.unit_models.pressure_changer import ThermodynamicAssumption

    from idaes.core.util.model_statistics import degrees_of_freedom

    print("IDAES Import Checks:        Passed")
    check_count += 1
except:
    print("IDAES Import Checks:        FAILED")

# -----------------------------------------------------------------------------
# Test available solvers
if SolverFactory('ipopt').available():
    print("Solver Availability Check:  Passed")
    check_count += 1
else:
    print("Solver Availability Check:  FAILED")

# -----------------------------------------------------------------------------
# Check model construction and solving
m = ConcreteModel()

m.fs = FlowsheetBlock(default={"dynamic": False})

m.fs.v = Var(m.fs.time)

def cons_rule(b, t):
    return b.v[t] == 1

m.fs.c = Constraint(m.fs.time, rule=cons_rule)

# Create a solver
solver = SolverFactory('ipopt')
results = solver.solve(m.fs)

if (results.solver.termination_condition == TerminationCondition.optimal
        and results.solver.status == SolverStatus.ok):
    print("Simple Model Check:         Passed")
    check_count += 1
else:
    print("Simple Model Check:         FAILED")

# -----------------------------------------------------------------------------
# Summary
print()
if check_count == 4:
    print("All Good!")
else:
    print("Something is not right. Please contact someone for assistance.")
