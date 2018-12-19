from pyomo.environ import (ConcreteModel, SolverFactory, TerminationCondition,
                           SolverStatus)

from idaes.core import FlowsheetBlock
from idaes.unit_models.flash import Flash as FL
from idaes.property_models import BFW_properties as properties
from idaes.ui.report import degrees_of_freedom


# -----------------------------------------------------------------------------
# Create a flowsheet for test
m = ConcreteModel()

m.param_block = properties.PhysicalParameterBlock()
m.prop_block = properties.StateBlock(parameters=m.param_block)

m.fs.display()

solver = SolverFactory('ipopt')
solver.options = {'tol': 1e-6,
                  'mu_init': 1e-8,
                  'bound_push': 1e-8}
