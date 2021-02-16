import pyomo.environ as pyo
from pyomo.opt.base.solvers import SolverFactoryClass as _SolverFactoryClass
import idaes

class SolverFactoryClass(_SolverFactoryClass):
    def __call__(self, _name=None, **kwargs):
        if _name is None:
            _name = idaes.cfg["default_solver"]
        if _name in idaes.cfg:
            for k, v in idaes.cfg[_name].items():
                if k not in kwargs:
                    kwargs[k] = v
                elif k == "options":
                    # options is in ConfigBlock and in kwargs, treat "options"
                    # special so individual options can have defaults not jut
                    # the whole options block
                    for opk, opv in v.items():
                        if opk not in kwargs["options"]:
                            kwargs["options"][opk] = opv
        return super().__call__(_name, **kwargs)

SolverFactory = SolverFactoryClass("idaes solver factory config wrapper")
pyo.PyomoSolverFactory = pyo.SolverFactory # Keep two available
pyo.IDAESSolverFactory = SolverFactory

def _update_solver_factory_registry():
    """Make sure the SolverFactory objects are in syc with the solver registry"""
    pyo.PyomoSolverFactory._cls = pyo.SolverFactory._cls
    pyo.IDAESSolverFactory._cls = pyo.SolverFactory._cls

def _use_pyomo_solver_factory():
    """Use the standard Pyomo SolverFactory object"""
    pyo.SolverFactory = pyo.PyomoSolverFactory

def _use_idaes_solver_factory():
    """Use the IDAES solver factory object that checks the IDAES config for
    default args"""
    pyo.SolverFactory = pyo.IDAESSolverFactory

# If this module is imported assume should use the IDAES version of SolverFactory
_update_solver_factory_registry()
_use_idaes_solver_factory()
