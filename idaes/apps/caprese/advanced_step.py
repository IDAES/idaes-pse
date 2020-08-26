from pyomo.environ import Var, Constraint, Param, value, Suffix

K_AUG_SUFFIXES = [
        ("ipopt_zL_out", Suffix.IMPORT),
        ("ipopt_zU_out", Suffix.IMPORT),
        ("ipopt_zL_in", Suffix.EXPORT),
        ("ipopt_zU_in", Suffix.EXPORT),
        ("dual", Suffix.IMPORT_EXPORT),
        ("npdp", Suffix.EXPORT),
        ("dof_v", Suffix.EXPORT),
        ]

class AdvancedStepManager(object):
    def __init__(self, block, wrt_vars, index):
        self.block = block
        self.wrt_vars = wrt_vars
        self.n_wrt_vars = len(wrt_vars)
        self.index = index

        # Remove this line if I find some reason to not construct
        # Constraints right away. If this happens, I will need 
        # to have a "constructed" flag that __enter__ can check.
        self.make_wrt_constraints()

    def make_param_constraints(self):
        block = self.block
        param_vars = self.wrt_vars
        n_wrt_vars = self.n_wrt_vars
        index = self.index

        # TODO: There could be multiple "wrt" components depending
        # on perspective. User should provide a name, and we should
        # use unique_component_name.
        block.wrt_set = Set(range(n_wrt_vars))

        wrt_param_dict = {
                i: value(wrt_vars[i][index])
                for i in range(n_wrt_vars)
                }
        block.wrt_param = Param(
                wrt_set,
                initialize=wrt_param_dict,
                mutable=True,
                )

        def wrt_constraint_rule(b, i):
            return wrt_vars[i][index] == b.wrt_param[i]
        block.wrt_constraint = Constraint(
                self.wrt_set, 
                rule=wrt_constraint_rule,
                )
        block.wrt_constraint.deactivate()

    def __enter__(self):
        index = self.index
        wrt_set = self.block.wrt_set
        self.block.wrt_constraints.activate()
        for i in wrt_set:
            # TODO: Should I check that wrt_vars were previously fixed?
            wrt_vars[i][index].unfix()
        return self

    def __exit__(self, e_type, e_val, e_tb):
        index = self.index
        wrt_set = sself.block.wrt_set
        self.block.wrt_constraints.deactivate()
        for i in wrt_set:
            wrt_vars[i][index].fix()
