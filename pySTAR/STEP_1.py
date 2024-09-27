#Build NLP to be solved by BARON
#NLP is formed by relaxing MINLP constraints of symbolic regression tree

from pyomo.environ import *
import numpy as np

#Symbolic regression tree (m.y relaxed)
def SRT(
X, #input data (n_data = X.shape[0], p = X.shape[1])
y, #target values
):
    #Hyperparameters
    depth = 3 ###
    eps = 1e-6
    c_lo = min(-100, np.min(X))
    c_up = max(100, np.max(X))
    v_lo = -100
    v_up = 100

    #Parameters given from dataset
    n_data = X.shape[0]
    n_data_set = [i+1 for i in range(n_data)] #to index of instances start with 1 and end with n_data
    p = X.shape[1]
    Xposi, Xnega, Xzero = XsetGenerator(X) #set of positive, negative, and zero values of X
    
    #Operators
    VARS = ["x_" + str(i+1) for i in range(p)]
    L = ["cst"] + VARS
    B = ["+", "-", "*", "/"]
    U = ["**0.5","log","exp"]
    Opair = []
    Opair.append(["log","exp"])
    O = B + U + L

    #Nodes and their allowed operators
    N, T = node_definer(depth) #Index of the nodes
    NnotT = np.setdiff1d(N,T)
    Nperfect = NnotT
    Y = Y_define(O, NnotT, L, T) #set of allowed pair of (operator, node index)


    ###### Model Set Up ############################################################################
    m = ConcreteModel()
    m.y = Var(Y, domain=UnitInterval)
    m.c = Var(N, domain=NonNegativeReals, bounds=(c_lo, c_up))
    m.v = Var(n_data_set, N, domain=Reals, bounds=(v_lo, v_up))
    #m.eps = Var(n_data_set, N, domain=NonNegativeReals, bounds=(eps,0.5), initialize=eps)

    #Objective
    def sse(m):
        #return sum((y[i-1] - m.v[i, 1]) ** 2 for i in range(1, n_data+1)) ##Scale Down
        return ((sum((y[i-1] - m.v[i, 1]) ** 2 for i in range(1, n_data+1)))/n_data)**0.5
    m.obj = Objective(rule=sse, sense=minimize)

    #Tree-defining constraints
    def tdc23a_rule(m,n): #if n is operator, 2*n+1 must exist
        return sum(m.y[(n,o)] for o in B + U) == sum(m.y[(2*n+1,o)] for o in O if (2*n+1,o) in Y)
    m.tdc23a = Constraint(NnotT, rule=tdc23a_rule)
    def tdc23b_rule(m,n): #If n is a binary operator, 2*n must exist
        return sum(m.y[(n,o)] for o in B) == sum(m.y[(2*n,o)] for o in O if (2*n,o) in Y)
    m.tdc23b = Constraint(NnotT, rule=tdc23b_rule)
    def tdcA1a_rule(m,n): #Only one operator or operand can be assigned or not assigned at all for each node
        return sum(m.y[(n,o)] for o in O if (n,o) in Y) <= 1
    m.tdcA1a = Constraint(N, rule=tdcA1a_rule)
    def tdcA1b_rule(m,n): #At least one variable must appear in the expression
        return sum(sum(m.y[(n,o)] for o in VARS if (n,o) in Y) for n in N) >= 1    
    m.tdcA1b = Constraint(rule=tdcA1b_rule)

    ### Value-defining constraints ################################################################
    def vdc24a_rule(m,i,n): #if n is variable, return the value of the variable. If it is not a variable let it to range between v_lo and v_up
        return m.v[i,n] <= sum(X[i-1, j-1] * m.y[(n, "x_" + str(j))] for j in range(1, p+1)) + v_up*sum(m.y[(n, o)] for o in B+U+["cst"] if (n,o) in Y)
    m.vdc24a = Constraint(n_data_set, N, rule=vdc24a_rule)
    def vdc24b_rule(m,i,n):
        return m.v[i,n] >= sum(X[i-1, j-1] * m.y[(n, "x_" + str(j))] for j in range(1, p+1)) + v_lo*sum(m.y[(n, o)] for o in B+U+["cst"] if (n,o) in Y)
    m.vdc24b = Constraint(n_data_set, N, rule=vdc24b_rule)

    def consta_rule(m,i,n): #If n is constant, return the value of the constant. If it is not a constant let it to range between gap between v bound and c bound
        return m.v[i,n] - m.c[n] <= (v_up - c_lo)*(1-m.y[(n,"cst")])
    m.consta = Constraint(n_data_set, N, rule=consta_rule)
    def constb_rule(m,i,n):
        return m.v[i,n] - m.c[n] >= (v_lo - c_up)*(1-m.y[(n,"cst")])
    m.constb = Constraint(n_data_set, N, rule=constb_rule)

    def adda_rule(m,i,n): #If n is addition, return the sum of the two children. If it is not an addition let it to range between gap between v bound
        return m.v[i,n] - (m.v[i, 2*n] + m.v[i,2*n+1]) <= (v_up - 2*v_lo)*(1-m.y[(n,"+")])
    m.adda = Constraint(n_data_set, NnotT, rule=adda_rule)
    def addb_rule(m,i,n):
        return m.v[i,n] - (m.v[i, 2*n] + m.v[i,2*n+1]) >= (v_lo - 2*v_up)*(1-m.y[(n,"+")])
    m.addb = Constraint(n_data_set, NnotT, rule=addb_rule)

    def suba_rule(m,i,n): #If n is subtraction, return the difference of the two children. If it is not a subtraction let it to range between gap between v bound
        return m.v[i,n] - (m.v[i, 2*n] - m.v[i,2*n+1]) <= (2*v_up - v_lo)*(1-m.y[(n,"-")])
    m.suba = Constraint(n_data_set, NnotT, rule=suba_rule)
    def subb_rule(m,i,n):
        return m.v[i,n] - (m.v[i, 2*n] - m.v[i,2*n+1]) >= (2*v_lo - v_up)*(1-m.y[(n,"-")])
    m.subb = Constraint(n_data_set, NnotT, rule=subb_rule)

    def multa_rule(m,i,n): #If n is multiplication, return the product of the two children. If it is not a multiplication let it to range between gap between v bound
        return m.v[i,n] - (m.v[i, 2*n] * m.v[i,2*n+1]) <= v_up - min(v_lo**2, v_lo*v_up, v_up**2)*(1-m.y[(n,"*")])
    m.multa = Constraint(n_data_set, NnotT, rule=multa_rule)
    def multb_rule(m,i,n):
        return m.v[i,n] - (m.v[i, 2*n] * m.v[i,2*n+1]) >= v_lo - max(v_lo**2, v_up**2)*(1-m.y[(n,"*")])
    m.multb = Constraint(n_data_set, NnotT, rule=multb_rule)

    def diva_rule(m, i, n): #If n is division, return the division of the two children. Do not let it get divided by zero
        return m.v[i, n] * m.v[i, 2 * n + 1] - (m.v[i, 2 * n]) <= (
            max(v_lo ** 2, v_up ** 2) - v_lo
        ) * (1 - m.y[(n, "/")])
    m.diva = Constraint(n_data_set, NnotT, rule=diva_rule)
    def divb_rule(m, i, n):
        return m.v[i, n] * m.v[i, 2 * n + 1] - (m.v[i, 2 * n]) >= (
            min(v_lo ** 2, v_lo * v_up, v_up ** 2) - v_up
        ) * (1 - m.y[(n, "/")])
    m.divb = Constraint(n_data_set, NnotT, rule=divb_rule)
    def divc_rule(m, i, n):
        #return eps * m.y[(n, "/")] <= m.v[i, 2 * n] ** 2
        return ((m.v[i, 2*n+1])**2)*m.y[(n,"/")] >= m.y[(n,"/")]-1 + eps
    m.divc = Constraint(n_data_set, NnotT, rule=divc_rule)
    # def divd_rule(m, i, n):
    #     return eps * m.y[(n, "/")] <= m.v[i, 2 * n + 1] ** 2
    # m.divd = Constraint(n_data_set, NnotT, rule=divd_rule)

    def sqra_rule(m, i, n): #If n is square root, return the square of the child. Do not let it get square rooted by negative values
        return m.v[i, n] ** 2 - (m.v[i, 2 * n + 1]) <= (
            max(v_lo ** 2, v_up ** 2) - v_lo
        ) * (1 - m.y[(n, "**0.5")])
    m.sqra = Constraint(n_data_set, NnotT, rule=sqra_rule)
    def sqrb_rule(m, i, n):
        return m.v[i, n] ** 2 - (m.v[i, 2 * n + 1]) >= (-v_up) * (1 - m.y[(n, "**0.5")])
    m.sqrb = Constraint(n_data_set, NnotT, rule=sqrb_rule)
    def sqrc_rule(m, i, n):
        # return eps - m.v[i, 2 * n + 1] <= (eps - v_lo) * (1 - m.y[(n, "**0.5")])
        return m.v[i, 2*n+1] * m.y[n, "**0.5"] >= 0
    m.sqrc = Constraint(n_data_set, NnotT, rule=sqrc_rule)
    
    def expa_rule(m, i, n): #If n is exponential, return the exponential of the child. 
        return m.v[i, n] - exp(m.v[i, 2 * n + 1]) <= v_up * (1 - m.y[(n, "exp")])
    m.expa = Constraint(n_data_set, NnotT, rule=expa_rule)
    def expb_rule(m, i, n):
        if v_up >= 10:
            e_up = 1e5
        else:
            e_up = exp(v_up)
        return m.v[i, n] - exp(m.v[i, 2 * n + 1]) >= (v_lo - e_up) * (
            1 - m.y[(n, "exp")]
        )
    m.expb = Constraint(n_data_set, NnotT, rule=expb_rule)

    def loga_rule(m, i, n): #If n is logarithm, return the logarithm of the child and do not let it get logarithmed by negative values
        if v_up - v_lo >= 10:
            e_up_lo = 1e5
        elif exp(v_up - v_lo) < -1e-5:
            e_up_lo = -1e5
        else:
            e_up_lo = exp(v_up - v_lo)
        return exp(m.v[i, n]) - m.v[i, 2 * n + 1] <= e_up_lo * (1 - m.y[(n, "log")])
    m.loga = Constraint(n_data_set, NnotT, rule=loga_rule)
    def logb_rule(m, i, n):
        return exp(m.v[i, n]) - m.v[i, 2 * n + 1] >= -v_up * (1 - m.y[(n, "log")])
    m.logb = Constraint(n_data_set, NnotT, rule=logb_rule)
    def logc_rule(m, i, n):
        #return eps - m.v[i, 2 * n + 1] <= (eps - v_lo) * (1 - m.y[(n, "log")])
        return m.v[i, 2*n+1] * m.y[n, "log"] >= m.y[n, "log"]-1 + eps
    m.logc = Constraint(n_data_set, NnotT, rule=logc_rule)

    # Redundancy-eliminating constraints
    def rec28a_rule(m, n):
            if (2 * n + 1, "-") in Y and (n, "+") in Y:
                return m.y[(n, "+")] + m.y[(2 * n + 1, "-")] <= 1
            else:
                return Constraint.Skip
    m.rec28a = Constraint(np.setdiff1d(N, Nperfect), rule=rec28a_rule)
    def rec28b_rule(m, n):
        if (2 * n + 1, "/") in Y and (n, "*") in Y:
            return m.y[(n, "*")] + m.y[(2 * n + 1, "/")] <= 1
        else:
            return Constraint.Skip
    m.rec28b = Constraint(np.setdiff1d(N, Nperfect), rule=rec28b_rule)
    def rec28c_rule(m, n):
        return m.y[(2 * n + 1, "cst")] <= m.y[(n, "+")] + m.y[(n, "*")]
    m.rec28c = Constraint(NnotT, rule=rec28c_rule)
    def recA12d_rule(m, n):
        return m.y[(2 * n, "cst")] + m.y[(2 * n + 1, "cst")] <= 1
    m.recA12d = Constraint(NnotT, rule=recA12d_rule)
    def recA12e_rule(m, n, o1, o2):
        if (2 * n + 1, o2) in Y:
            return m.y[(n, o1)] + m.y[(2 * n + 1, o2)] <= 1
        else:
            return Constraint.Skip
    m.recA12e = Constraint(NnotT, Opair, rule=recA12e_rule)
    def recA12f_rule(m, n, o1, o2):
        if (2 * n + 1, o1) in Y:
            return m.y[(n, o2)] + m.y[(2 * n + 1, o1)] <= 1
        else:
            return Constraint.Skip
    m.recA12f = Constraint(NnotT, Opair, rule=recA12f_rule)

    #Implication cuts
    def ic12a_rule(m, j, n):
        return m.y[(n, "/")] + m.y[(2 * n + 1, "x_" + str(j))] <= 1
    m.ic12a = Constraint(Xzero, NnotT, rule=ic12a_rule)
    def ic12b_rule(m, j, n):
        return m.y[(n, "**0.5")] + m.y[(2 * n + 1, "x_" + str(j))] <= 1
    m.ic12b = Constraint(Xnega, NnotT, rule=ic12b_rule)
    def ic12c_rule(m, j, n):
        return m.y[(n, "log")] + m.y[(2 * n + 1, "x_" + str(j))] <= 1
    m.ic12c = Constraint(Xnega+Xzero, NnotT, rule=ic12c_rule)

    #symmetry breaking constraints
    def sbc14_rule(m, n):
        return m.v[1, 2 * n] - m.v[1, 2 * n + 1] >= (v_lo - v_up) * (
            1 - m.y[(n, "+")] - m.y[(n, "*")]
        )
    m.sbc14 = Constraint(Nperfect, rule=sbc14_rule)

    return m, Y, B, U, L, NnotT, T, N, c_lo, c_up



def XsetGenerator(X): #list of positive, negative, and zero values of X
    Xposi = []
    Xnega = []
    Xzero = []
    for var, x in enumerate(X.T):
        if np.any(x > 0):
            Xposi.append(var+1)
        if np.any(x < 0):
            Xnega.append(var+1)
        if np.any(x == 0):
            Xzero.append(var+1)
    return Xposi, Xnega, Xzero

def node_definer(depth): #Index of the nodes
    Nodes = [1]
    TerminalNodes = []
    Old_Level_Nodes = [1]
    for level in range(1, depth + 1):
        New_Level_Nodes = []
        for node in Old_Level_Nodes:
            New_Level_Nodes.append(2*node)  
            New_Level_Nodes.append(2*node + 1)
        Nodes += New_Level_Nodes  
        if level != depth:
            Old_Level_Nodes = New_Level_Nodes
        else:
            TerminalNodes += New_Level_Nodes  
    return Nodes, TerminalNodes

def Y_define(O, NnotT, L, T): #set of allowed pair of (operator, node index)
    Y = []
    for o in O:
        for n in NnotT:
            Y.append((n, o))
    for l in L:
        for t in T:
            Y.append((t, l))
    return Y

