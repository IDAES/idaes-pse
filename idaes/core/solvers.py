from pyomo.environ import SolverFactory as _SolverFactory
from pyomo.opt.base.solvers import SolverFactoryClass as _SolverFactoryClass
import idaes

class SolverFactoryClass(_SolverFactoryClass):
    def __call__(self, _name=None, **kwargs):
        if _name is None:
            _name = idaes.cfg["default_solver"]
        s = super().__call__(_name, **kwargs)
        if _name in idaes.cfg and "options" in idaes.cfg[_name]:
            if hasattr(idaes.cfg[_name]["options"], "value"): #ConfigBlock
                s.options.update(idaes.cfg[_name]["options"].value())
            else: # not ConfigBlock
                s.options.update(idaes.cfg[_name][options])
        return s

SolverFactory = SolverFactoryClass("idaes solver factory config wrapper")
# since this is a new object, link to Pyomo SolverFactory registered solvers
SolverFactory._cls = _SolverFactory._cls
