from pyomo.solvers.plugins.solvers.IPOPT import IPOPT
from pyomo.environ import SolverFactory

import idaes

@SolverFactory.register('ipopt-idaes', doc='The Ipopt NLP solver, with IDAES options')
class IpoptIdaes(IPOPT):
    def __init__(self, **kwds):
        super().__init__(**kwds)
        for k, v in idaes.cfg["ipopt-idaes"]["options"].items():
            self.options[k] = v
