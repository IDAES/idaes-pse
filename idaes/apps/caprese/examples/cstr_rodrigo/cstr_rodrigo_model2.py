# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 19:45:54 2020

@author: greg6
"""

# from pyomo.environ import (Block, ConcreteModel,  Constraint, Expression,
#                            Set, SolverFactory, Var, value, 
#                            TransformationFactory, TerminationCondition)
from pyomo.environ import*
from pyomo.dae import*

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

    m.F = Param(m.t, mutable=True, default=1.2000000000000000E+02)
    m.Fw = Param(m.t, mutable=True, default=3.0000000000000000E+01)
    
    #States
    m.Ca = Var(m.t,  initialize=1.60659680385930765667001907104350E-02)
    # m.T = Var(m.t,  initialize=3.92336059452774350120307644829154E+02)
    # m.Tj = Var(m.t,  initialize=3.77995395658401662331016268581152E+02)
    m.Tall = Var(m.t, m.T_ind)
    
    #Algebraic var
    m.k = Var(m.t,  initialize=4.70706140E+02)
    
    #Controls
    m.Tjinb = Var(m.t, initialize=250)
    
    #: These guys have to be zero at the steady-state (steady).
    zero0 = dict.fromkeys(m.t)
    for key in zero0.keys():
        zero0[key] = 0.0
    if steady:
        m.Cadot = Var(m.t, initialize = 0.0)
        m.Tdot = Var(m.t, initialize = 0.0)
        m.Tjdot = Var(m.t, initialize = 0.0)
        m.Cadot.fix()
        m.Tdot.fix()
        m.Tjdot.fix()
    else:
        m.Cadot = DerivativeVar(m.Ca, initialize=-3.58709135E+01)
        # m.Tdot = DerivativeVar(m.T, initialize=5.19191848E+03)
        # m.Tjdot = DerivativeVar(m.Tj, initialize=-9.70467399E+02)
        m.Talldot = DerivativeVar(m.Tall, wrt = m.t)
        
    #Constraint rules

    # Using Constraint.Skip on equations at t0 is not supported right now.
    # This is only because I will attempt to solve for consistent initial
    # conditions, but at the initial conditions, this model will not exist.
    # TODO: Add as option ability to solve/not solve for consistent ICs
    def _rule_k(m, i):
        #if i == 0:
        #    return Constraint.Skip
        #else:
        return m.k[i] == m.k0 * exp(-m.Er / m.Tall[i, "T"])

    def _rule_ca(m, i):
        #if i == 0:
        #    return Constraint.Skip
        #else:
        rule = (m.Cadot[i] == (m.F[i] / m.V) * 
                              (m.Cainb - m.Ca[i]) - 2 * 
                              m.k[i] * 
                              m.Ca[i] ** 2)
        return rule

    def _rule_t(m, i):
        #if i == 0:
        #    return Constraint.Skip
        #else:
        return m.Talldot[i, "T"] == (m.F[i] / m.V) * (m.Tinb - m.Tall[i, "T"]) + \
        2.0 * m.dH / (m.rho * m.Cp) * m.k[i] * m.Ca[i] ** 2 -\
        m.UA / (m.V * m.rho * m.Cp) * (m.Tall[i, "T"] - m.Tall[i, "Tj"])

    def _rule_tj(m, i):
        #if i == 0:
        #    return Constraint.Skip
        #else:
        return (m.Talldot[i, "Tj"] ==
               (m.Fw[i] / m.Vw) * 
               (m.Tjinb[i] - m.Tall[i, "Tj"]) + 
               m.UA / 
               (m.Vw * m.rhow * m.Cpw) * 
               (m.Tall[i, "T"] - m.Tall[i, "Tj"]))
    
    #Constraints
    m.kdef = Constraint(m.t, rule = _rule_k)
    m.de_ca = Constraint(m.t, rule = _rule_ca)
    m.de_T = Constraint(m.t, rule = _rule_t)
    m.de_Tj = Constraint(m.t, rule = _rule_tj)
    
    # #ics rule
    # m.Ca_ic = Param( default=1.9193793974995963E-02, mutable=True)
    # m.T_ic = Param( default=3.8400724261199036E+02, mutable=True)
    # m.Tj_ic = Param( default=3.7127352272578315E+02, mutable=True)
    
    # # # #
    # # The following are not used:
    # def _rule_ca0(m):
    #     return m.Ca[0] == m.Ca_ic

    # def _rule_t0(m):
    #     return m.T[0] == m.T_ic

    # def _rule_tj0(m):
    #     return m.Tj[0] == m.Tj_ic
    # # # #
        
    # Time discretization
    if not steady:
        disc = TransformationFactory('dae.collocation')
        disc.apply_to(m, wrt=m.t, nfe=ntfe, ncp=ntcp, scheme='LAGRANGE-RADAU')
    m.Tall[:, "T"] = 3.92336059452774350120307644829154E+02
    m.Tall[:, "Tj"] = 3.77995395658401662331016268581152E+02
    m.Talldot[:, "T"] = 5.19191848E+03
    m.Talldot[:, "Tj"] = -9.70467399E+02
    
    if not steady:
        m.Ca[0].fix(1.9193793974995963E-02)
        m.Tall[0, "T"].fix(3.8400724261199036E+02)
        m.Tall[0, "Tj"].fix(3.7127352272578315E+02)
        
    # Fix control inputs (why????
    m.Tjinb.fix(250.0)
    
    #m.not_use = Var(m.t)
    #m.not_use.fix(1.0)
        
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
