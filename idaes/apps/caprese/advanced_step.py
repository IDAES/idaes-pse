from pyomo.core.base.constraint import Var, Constraint, Param, value

class AdvancedStepManager(object):
    def __init__(self, block, dof_vars, index):
        self.block = block
        self.dof_vars = dof_vars
        self.n_dof_vars = len(dof_vars)
        self.index = index

        # Remove this line if I find some reason to not construct
        # Constraints right away. If this happens, I will need 
        # to have a "constructed" flag that __enter__ can check.
        self.make_dof_constraints()

    def make_dof_constraints(self):
        block = self.block
        dof_vars = self.dof_vars
        n_dof_vars = self.n_dof_vars
        index = self.index

        # TODO: There could be multiple "dof" components depending
        # on perspective. User should provide a name, and we should
        # use unique_component_name.
        block.dof_set = Set(range(n_dof_vars))

        dof_param_dict = {
                i: value(dof_vars[i][index])
                for i in range(n_dof_vars)
                }
        block.dof_param = Param(dof_set, initialize=dof_param_dict)

        def dof_constraint_rule(b, i):
            return dof_vars[i][index] == b.dof_param[i]
        block.dof_constraint = Constraint(
                self.dof_set, 
                rule=dof_constraint_rule,
                )
        block.dof_constraint.deactivate()

    def __enter__(self):
        index = self.index
        dof_set = self.block.dof_set
        self.block.dof_constraints.activate()
        for i in dof_set:
            # TODO: Should I check that dof_vars were previously fixed?
            dof_vars[i][index].unfix()
        return self

    def __exit__(self, e_type, e_val, e_tb):
        index = self.index
        dof_set = sself.block.dof_set
        self.block.dof_constraints.deactivate()
        for i in dof_set:
            dof_vars[i][index].fix()
