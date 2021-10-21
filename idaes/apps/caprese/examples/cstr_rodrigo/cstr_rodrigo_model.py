#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
"""
Caprese example for NMPC with Rodrigo's CSTR model
"""

from pyomo.environ import (Block, ConcreteModel,  Constraint, Expression,
                            Set, SolverFactory, Var, value, 
                            TransformationFactory, TerminationCondition, Param,
                            exp)
from pyomo.dae.contset import ContinuousSet
from pyomo.dae.diffvar import DerivativeVar

from idaes.core.util.model_statistics import degrees_of_freedom


__author__ = "Kuan-Han Lin"


# See if ipopt is available and set up solver
if SolverFactory('ipopt').available():
    solver = SolverFactory('ipopt')
    solver.options = {'tol': 1e-6,
                      'mu_init': 1e-8,
                      'bound_push': 1e-8}
else:
    solver = None
    

def make_model(horizon=10, ntfe=5, ntcp=3, steady=False, bounds=False):
    #sampling time = 2
    
    m = ConcreteModel(name = "CSTR_rodrigo")
    
    if steady:
        m.t = Set(initialize = [1.0])
    else:
        m.t = ContinuousSet(bounds=(0, horizon))
        
    m.T_ind = Set(initialize = ["T", "Tj"])
    
    #Parameters
    m.Cainb = Param(default=1.0)
    m.Tinb = Param(default=275.0)
        
    m.V = Param(initialize=100)
    m.UA = Param(initialize=20000 * 60)
    m.rho = Param(initialize=1000)
    m.Cp = Param(initialize=4.2)
    m.Vw = Param(initialize=10)
    m.rhow = Param(initialize=1000)
    m.Cpw = Param(initialize=4.2)
    m.k0 = Param(initialize=4.11e13)
    m.E = Param(initialize=76534.704)
    m.R = Param(initialize=8.314472)
    m.Er = Param(initialize=lambda m: (value(m.E) / value(m.R)))
    m.dH = Param(initialize=596619.)

    m.F = Param(m.t, mutable=True, default=1.20E+02)
    m.Fw = Param(m.t, mutable=True, default=3.00E+01)
    
    #States
    m.Ca = Var(m.t,  initialize=1.61E-02)
    m.Tall = Var(m.t, m.T_ind)
    
    #Algebraic var
    m.k = Var(m.t,  initialize=4.71E+02)
    
    #Controls
    m.Tjinb = Var(m.t, initialize=250)
    
    #Derivative variables
    if steady:
        m.Cadot = Var(m.t, initialize = 0.0)
        m.Talldot = Var(m.t, m.T_ind, initialize = 0.0)
        m.Cadot.fix()
        m.Talldot.fix()
    else:
        m.Cadot = DerivativeVar(m.Ca, initialize=-3.59E+01)
        m.Talldot = DerivativeVar(m.Tall, wrt = m.t)
        
    #Constraint rules
    def _rule_k(m, i):
        return m.k[i] == m.k0 * exp(-m.Er / m.Tall[i, "T"])

    def _rule_ca(m, i):
        rule = (m.Cadot[i] == (m.F[i] / m.V) * 
                              (m.Cainb - m.Ca[i]) - 2 * 
                              m.k[i] * 
                              m.Ca[i] ** 2)
        return rule

    def _rule_t(m, i):
        return m.Talldot[i, "T"] == ((m.F[i] / m.V) * (m.Tinb - m.Tall[i, "T"]) + 
                                     2.0 * m.dH / (m.rho * m.Cp) * m.k[i] * m.Ca[i] ** 2 -
                                     m.UA / (m.V * m.rho * m.Cp) * (m.Tall[i, "T"] - m.Tall[i, "Tj"]))

    def _rule_tj(m, i):
        return m.Talldot[i, "Tj"] == ((m.Fw[i] / m.Vw) * (m.Tjinb[i] - m.Tall[i, "Tj"]) + 
                                      m.UA/(m.Vw * m.rhow * m.Cpw) * (m.Tall[i, "T"] - m.Tall[i, "Tj"]))
    
    #Constraints
    m.kdef = Constraint(m.t, rule = _rule_k)
    m.de_ca = Constraint(m.t, rule = _rule_ca)
    m.de_T = Constraint(m.t, rule = _rule_t)
    m.de_Tj = Constraint(m.t, rule = _rule_tj)
        
    # Time discretization
    if not steady:
        disc = TransformationFactory('dae.collocation')
        disc.apply_to(m, wrt=m.t, nfe=ntfe, ncp=ntcp, scheme='LAGRANGE-RADAU')
        
    # Set initialize values for IndexedVar with string index
    m.Tall[:, "T"] = 3.92E+02
    m.Tall[:, "Tj"] = 3.78E+02
    m.Talldot[:, "T"] = 5.19E+03
    m.Talldot[:, "Tj"] = -9.7E+02
    
    #Set initial conditions for states
    if not steady:
        m.Ca[0].fix(1.92E-02)
        m.Tall[0, "T"].fix(3.84E+02)
        m.Tall[0, "Tj"].fix(3.71E+02)
    
    #Fix input here to pass dof check 
    #(not necessary for running nmpc; it will be fixed if it's not fix here)
    m.Tjinb.fix(250.0)
    
        
    #bounds
    if bounds:
        m.Ca.setlb(0.0)
        m.Ca.setub(0.04)
        m.Tall[:,"T"].setlb(3.5E+02)
        m.Tall[:,"T"].setub(4.5E+02)
        m.Tall[:,"Tj"].setlb(3.5E+02)
        m.Tall[:,"Tj"].setub(4.5E+02)
        m.Tjinb.setlb(200.)
        m.Tjinb.setub(300.)
        
    return m

# states = ["Ca", "T", "Tj"]
# controls = ["Tjinb"]

#first ss
# Ca: 0.019193793974995824
# T: 384.0072426119906
# Tj: 371.2735227257833
# Tjinb(control input): 250

#second ss
# Ca: 0.017999999838289406
# T: 386.0954447399273
# Tj: 374.1675457210488
# Tjinb(control input): 260.5685074460155

if __name__ == '__main__':
    m_plant = make_model()
    m_ss = make_model(steady = True)
    
    assert degrees_of_freedom(m_plant) == 0
    assert degrees_of_freedom(m_ss) == 0
    
