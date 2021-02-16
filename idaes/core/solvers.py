from pyomo.environ import SolverFactory
import idaes

class SolverWrapper(object):
    def __init__(self, name):
        self.name = name
        if name == 'default':
            self.solver = None
            doc = "IDAES Configured Default Solver"
        else:
            self.solver = SolverFactory.get_class(name)
            doc = SolverFactory.doc(name)
        SolverFactory.unregister(name)
        # We should really make a public method to register a class without using a decorator
        SolverFactory.register(name, doc)(self)

    def __call__(self, *args, **kwargs):
        if self.name == "default":
            name = idaes.cfg.default_solver
            solver = SolverFactory.get_class(name)
        else:
            name = self.name
            solver = self.solver
        if name in idaes.cfg:
            for k, v in idaes.cfg[name].items():
                if k not in kwargs:
                    kwargs[k] = v
                elif k == "options":
                    # options is in ConfigBlock and in kwargs, treat "options"
                    # special so individual options can have defaults not jut
                    # the whole options block
                    for opk, opv in v.items():
                        if opk not in kwargs["options"]:
                            kwargs["options"][opk] = opv
        return solver(*args, **kwargs)

for c in list(SolverFactory):
    SolverWrapper(c)
SolverWrapper('default')
