from pyomo.environ import SolverFactory as _SolverFactory
from pyomo.opt.base.solvers import SolverFactoryClass as _SolverFactoryClass
import idaes

class SolverFactoryClass(_SolverFactoryClass):
    def __call__(self, _name, **kwargs):
        s = super().__call__(_name, **kwargs)
        if _name in idaes.cfg and "options" in idaes.cfg[_name]:
            for k, v in idaes.cfg[_name]["options"].items():
                if k not in s.options:
                    s.options[k] = v
        return s

SolverFactory = SolverFactoryClass("idaes solver factory config wrapper")
# since this is a new object, link to Pyomo SolverFactory registered solvers
SolverFactory._cls = _SolverFactory._cls
