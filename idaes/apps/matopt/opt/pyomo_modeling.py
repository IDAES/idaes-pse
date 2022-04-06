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
from idaes.logger import getModelLogger

logging = getModelLogger("MatOptModel")

from pyomo.environ import *
from pyomo.core.base.var import _GeneralVarData
from pyomo.core.expr.numeric_expr import (
    MonomialTermExpression,
    SumExpression,
    NegationExpression,
)
from ..util.util import isZero, areEqual

DBL_TOL = 1e-5
DEFAULT_BIG_M = 9999
DEFAULT_EPS = DBL_TOL


# ================================================
# ==========      UTILITY FUNCTIONS     ==========
# ================================================
def getLB(e):
    """Calculates the lower bound (if possible) of a linear expression.

    Args:
        e: A Pyomo expression, Pyomo variable, float, or int.

    Returns:
        (float) lower bound from interval arithmetic. None if negative infinity.

    """
    # Future work: use Pyomo function to achieve this functionality
    # return compute_bounds_on_expr(e)[0]
    if isinstance(e, _GeneralVarData):
        if e.is_fixed():
            return value(e)
        else:
            return e.lb
    elif isinstance(e, MonomialTermExpression):
        assert isinstance(e.args[1], pyomo.core.base.var._GeneralVarData), (
            "This code relies on the assumption that the only variable "
            "in a monomial expression is the second argument"
        )
        if e.args[0] > 0:
            sub_expr_result = getLB(e.args[1])
            return e.args[0] * sub_expr_result if sub_expr_result is not None else None
        else:
            sub_expr_result = getUB(e.args[1])
            return e.args[0] * sub_expr_result if sub_expr_result is not None else None
    elif isinstance(e, SumExpression):
        result = 0.0
        for elem in e.args:
            elem_lb = getLB(elem)
            if elem_lb is not None:
                result += elem_lb
            else:
                return None
        return result
    elif isinstance(e, NegationExpression):
        elem_ub = getUB(e.args[0])
        return -elem_ub if elem_ub is not None else None
    elif isinstance(e, (float, int)):
        return e
    else:
        print(type(e))
        raise NotImplementedError(
            "Unsupported expression type for lower bound calculation. Please check getLB() function in "
            "pyomo_modeling.py for supported types."
        )


def getUB(e):
    """Calculates the upper bound (if possible) of a linear expression.

    Args:
        e: A Pyomo expression, Pyomo variable, float, or int.

    Returns:
        (float) upper bound from interval arithmetic. None if infinity.

    """
    # Future work: use Pyomo function to achieve this functionality
    # return compute_bounds_on_expr(e)[1]
    if isinstance(e, _GeneralVarData):
        if e.is_fixed():
            return value(e)
        else:
            return e.ub
    elif isinstance(e, MonomialTermExpression):
        assert isinstance(e.args[1], pyomo.core.base.var._GeneralVarData), (
            "This code relies on the assumption that the only variable "
            "in a monomial expression is the second argument"
        )
        if e.args[0] > 0:
            sub_expr_result = getUB(e.args[1])
            return e.args[0] * sub_expr_result if sub_expr_result is not None else None
        else:
            sub_expr_result = getLB(e.args[1])
            return e.args[0] * sub_expr_result if sub_expr_result is not None else None
    elif isinstance(e, SumExpression):
        result = 0.0
        for elem in e.args:
            elem_ub = getUB(elem)
            if elem_ub is not None:
                result += elem_ub
            else:
                return None
        return result
    elif isinstance(e, NegationExpression):
        elem_lb = getLB(e.args[0])
        return -elem_lb if elem_lb is not None else None
    elif isinstance(e, (float, int)):
        return e
    else:
        print(type(e))
        raise NotImplementedError(
            "Unsupported expression type for upper bound calculation. Please check getUB() function in "
            "pyomo_modeling.py for supported types."
        )


# ================================================
# ==========   GENERIC MODEL ELEMENTS   ==========
# ================================================
def makeMyPyomoBaseModel(C, Atoms=None, Confs=None):
    """
    Make the Pyomo model for a basic materials design problem.

    Creates the basic sets and variables that make up the problem.
    All variables are created sparsely, so they are initialized
    only after they are actually referenced in a constraint. In
    this way, we smartly eliminate unnecessary variables and
    constraints.

    Basic Variables:
        Yi: Presence of building block at site i
        Xij: Presence of building blocks at both sites i and j
        Ci: Count of building block bonds to neighbors to site i

    Common Descriptors:
        Zi: Presence of target site at i

    Atom-Specific Variables:
        Yik: Presence of building block of type k at site i
        Xijkl: Presence of type k and l at sites i and j, respectively
        Cikl: Count of neighbors of type l next to site i with type k

    Conformation-Specific Descriptors:
        Zic: Presence of conformation c at site i

    The basic variables and atom-specific variables above are
    automatically encoded by calling addConsForGeneralVars.
    The descriptor variables Zi and Zic must be explicitly constrained.
    Some standard approaches are formalized in:

    * addConsBoundDescriptorsWithImpl
    * addConsIndFromDescriptors
    * addConsZicFromYi
    * addConsZicFromYiLifted
    * addConsZicFromYik
    * addConsZicFromYikLifted

    Additional variables and constraints  can be created and managed on
    a model-specific level.

    Args:
        C (Canvas): The design space over which to model the problem.
        Atoms (list<Atom/Any>): Optional, the set of building blocks.
            If present, generates a meaningful set for K building blocks.
            Else, only the general variables for presence/absence are
            meaningful. (Default value = None)
        Confs (list<Design>): Optional, the set of conformations to
            potentialy use for indicator variables.
            (Default value = None)

    Returns:
        (ConcreteModel): Pyomo model with basic sets and variables initialized.

    """
    m = ConcreteModel()
    # Adding Model formulation information
    m.Canvas = C
    m.Ni = C.NeighborhoodIndexes
    m.Atoms = Atoms
    m.Confs = Confs
    # Adding Basic Sets
    m.nI = len(C)
    m.I = Set(initialize=range(m.nI), ordered=True)
    # Adding Basic Variables
    m.Yi = Var(m.I, domain=Binary, dense=False)
    m.Xij = Var(m.I, m.I, domain=Binary, dense=False)

    def _ruleCiBounds(m, i):
        return 0, len(m.Ni[i])

    m.Ci = Var(m.I, domain=NonNegativeIntegers, bounds=_ruleCiBounds, dense=False)
    # Adding Common Descriptors
    m.Zi = Var(m.I, domain=Binary, dense=False)
    # Adding Atoms Set
    m.nK = len(Atoms) if Atoms is not None else 0
    m.K = Set(initialize=Atoms)
    m.Yik = Var(m.I, m.K, domain=Binary, dense=False)
    m.Xijkl = Var(m.I, m.I, m.K, m.K, domain=Binary, dense=False)

    def _ruleCiklBounds(m, i, k, l):
        return 0, len(m.Ni[i])

    m.Cikl = Var(
        m.I, m.K, m.K, domain=NonNegativeIntegers, bounds=_ruleCiklBounds, dense=False
    )
    # Adding Confs Set
    m.nC = len(Confs) if Confs is not None else 0
    m.C = Set(initialize=range(m.nC))
    m.Zic = Var(m.I, m.C, domain=Binary, dense=False)
    return m


def _addConsCiFromCikl(m):
    m.AssignCiFromCikl = Constraint(m.I)
    for i in m.Ci.keys():
        m.AssignCiFromCikl.add(
            index=i,
            expr=(
                m.Ci[i]
                == sum(
                    m.Cikl[i, k, l]
                    for k in m.K
                    for l in m.K
                    if k is not None and l is not None
                )
            ),
        )


def _addConsCiFromXij(m):
    m.AssignCiFromXij = Constraint(m.I)
    for i in m.Ci.keys():
        m.AssignCiFromXij.add(
            index=i,
            expr=(m.Ci[i] == sum(m.Xij[i, j] for j in m.Ni[i] if j is not None)),
        )


def _addConsCiklFromXijkl(m):
    m.AssignCiklFromXijkl = Constraint(m.I, m.K, m.K)
    for i, k, l in m.Cikl.keys():
        m.AssignCiklFromXijkl.add(
            index=(i, k, l),
            expr=(
                m.Cikl[i, k, l]
                == sum(m.Xijkl[i, j, k, l] for j in m.Ni[i] if j is not None)
            ),
        )


def _addConsXijFromXijkl(m):
    m.AssignXijFromXijkl = Constraint(m.I, m.I)
    for i, j in m.Xij.keys():
        m.AssignXijFromXijkl.add(
            index=(i, j),
            expr=(
                m.Xij[i, j]
                == sum(
                    m.Xijkl[i, j, k, l]
                    for k in m.K
                    for l in m.K
                    if k is not None and l is not None
                )
            ),
        )


def _addConsXijFromYi(m):
    m.AssignXijFromYi1 = Constraint(m.I, m.I)
    m.AssignXijFromYi2 = Constraint(m.I, m.I)
    m.AssignXijFromYi3 = Constraint(m.I, m.I)
    for i, j in m.Xij.keys():
        m.AssignXijFromYi1.add(index=(i, j), expr=(m.Xij[i, j] <= m.Yi[i]))
        m.AssignXijFromYi2.add(index=(i, j), expr=(m.Xij[i, j] <= m.Yi[j]))
        m.AssignXijFromYi3.add(
            index=(i, j), expr=(m.Xij[i, j] >= m.Yi[i] + m.Yi[j] - 1)
        )


def _addConsXijklFromYik(m):
    m.AssignXijklFromYik1 = Constraint(m.I, m.I, m.K, m.K)
    m.AssignXijklFromYik2 = Constraint(m.I, m.I, m.K, m.K)
    m.AssignXijklFromYik3 = Constraint(m.I, m.I, m.K, m.K)
    for i, j, k, l in m.Xijkl.keys():
        m.AssignXijklFromYik1.add(
            index=(i, j, k, l), expr=(m.Xijkl[i, j, k, l] <= m.Yik[i, k])
        )
        m.AssignXijklFromYik2.add(
            index=(i, j, k, l), expr=(m.Xijkl[i, j, k, l] <= m.Yik[j, l])
        )
        m.AssignXijklFromYik3.add(
            index=(i, j, k, l),
            expr=(m.Xijkl[i, j, k, l] >= m.Yik[i, k] + m.Yik[j, l] - 1),
        )


def _addConsYiFromYik(m):
    m.AssignYiFromYik = Constraint(m.I)
    for i in m.Yi.keys():
        m.AssignYiFromYik.add(
            index=i, expr=(m.Yi[i] == sum(m.Yik[i, k] for k in m.K if k is not None))
        )


def _addConsYikSOS1(m):
    m.AssignYikSOS1 = Constraint(m.I)
    for i in m.Yik.index_set().set_tuple[0]:
        m.AssignYikSOS1.add(index=i, expr=(sum(m.Yik[i, k] for k in m.K) <= 1))


def addConsForGeneralVars(m):
    """Scan over the model and encode constraints for all basic variables present.

    Encodes variables in a chain from derived to most basc:
        Ci   <- Xij   <- Yi
        Cikl <- Xijkl <- Yik

    Also can bridge from type-specific to type-agnostic variables:
        Ci  <- Cikl
        Xij <- Xijkl
        Yi  <- Yik

    Args:
        m (ConcreteModel): Model to encode basic constraints.

    Returns:
        None.

    Raises:
        NotImplementedErorr: Several competing encodings exist and make
            the set of constraints not clearly defined.

    """
    # Defining Ci
    if len(m.Ci) > 0:
        if len(m.Xij) > 0 and len(m.Cikl) > 0:
            _addConsCiFromCikl(m)
            # NOTE: Chose to add only one formulation because the second is redundant
        elif len(m.Cikl) > 0 or len(m.Xijkl) > 0 or len(m.Yik) > 0:
            _addConsCiFromCikl(m)
        elif len(m.Xij) > 0 or len(m.Yi) > 0 or len(m.Yik) == 0:
            _addConsCiFromXij(m)
        else:
            raise NotImplementedError(
                "Ci lacks a proper way to be defined. User should define MaterialDescriptor Ci explicitly using valid "
                "DescriptorRule"
            )
            # Defining Cikl
    if len(m.Cikl) > 0:
        _addConsCiklFromXijkl(m)
    # Defining Xij
    if len(m.Xij) > 0:
        if len(m.Xijkl) > 0 and len(m.Yi) > 0:
            # NOTE: Chose to add both formulation because we believe neither is dominating
            _addConsXijFromXijkl(m)
            _addConsXijFromYi(m)
        elif len(m.Xijkl) > 0 or len(m.Yik) > 0:
            _addConsXijFromXijkl(m)
        elif len(m.Yi) > 0 or len(m.Yik) == 0:
            _addConsXijFromYi(m)
        else:
            raise NotImplementedError(
                "Xij lacks a propper way to be defined. User should define MaterialDescriptor Xij explicitly using "
                "valid DescriptorRule"
            )
            # Defining Xijkl
    if len(m.Xijkl) > 0:
        _addConsXijklFromYik(m)
    # Defining Zic
    if len(m.Zic) > 0:
        if m.nK > 1:
            addConsZicFromYikLifted(m)
        else:
            addConsZicFromYiLifted(m)
    # Defining Yi
    if len(m.Yi) > 0:
        if len(m.Yik) > 0:
            _addConsYiFromYik(m)
        else:
            # Yi can serve as a suitable basic variable if no 'k' index is in the model
            pass
    # Define Yik
    # Still need to restrict Yik to be SOS1
    if len(m.K) > 1:
        _addConsYikSOS1(m)


def fixYik(m, i, k, val):
    if isZero(val, DBL_TOL):
        fixYikDown(m, i, k)
    else:
        assert areEqual(val, 1.0, DBL_TOL)
        fixYikUp(m, i, k)


def fixYi(m, i, val):
    if isZero(val, DBL_TOL):
        fixYiDown(m, i)
    else:
        assert areEqual(val, 1.0, DBL_TOL)
        fixYiUp(m, i)


def fixXijkl(m, i, j, k, l, val):
    if isZero(val, DBL_TOL):
        fixXijklDown(m, i, j, k, l)
    else:
        assert areEqual(val, 1.0, DBL_TOL)
        fixXijklUp(m, i, j, k, l)


def fixXij(m, i, j, val):
    if isZero(val, DBL_TOL):
        fixXijDown(m, i, j)
    else:
        assert areEqual(val, 1.0, DBL_TOL)
        fixXijUp(m, i, j)


def fixCikl(m, i, k, l, val):
    m.Cikl[i, k, l].fix(val)
    checkCiklBounds(m, i, k, l)


def fixCi(m, i, val):
    m.Ci[i].fix(val)
    checkCiBounds(m, i)


def fixZic(m, i, c, val):
    if isZero(val, DBL_TOL):
        fixZicDown(m, i, c)
    else:
        assert areEqual(val, 1.0, DBL_TOL)
        fixZicUp(m, i, c)


def fixYikUp(m, i, k):
    if m.Yik[i, k].is_fixed():
        assert value(m.Yik[i, k]) > 0.5
        return  # unwind the stack if this has already been fixed
    for kp in m.K:
        if kp != k and m.Yik[i, kp].is_fixed() and value(m.Yik[i, kp]) > 0.5:
            raise ValueError("Tried to set multiple (k) in Yik to 1")
    m.Yik[i, k].fix(1)
    fixYiUp(m, i)
    for j in m.Ni[i]:
        if j is not None:
            for l in m.K:
                if getLB(m.Yik[j, l]) > 0.5:
                    fixXijklUp(m, i, j, k, l)
    for ip in m.I:
        if i in m.Ni[ip]:
            for l in m.K:
                if getLB(m.Yik[ip, l]) > 0.5:
                    fixXijklUp(m, ip, i, l, k)
    for kp in m.K:
        if kp != k:
            fixYikDown(m, i, kp)


def fixYiUp(m, i):
    if m.Yi[i].is_fixed():
        assert value(m.Yi[i]) > 0.5
        return
    m.Yi[i].fix(1)
    for j in m.Ni[i]:
        if j is not None:
            if getLB(m.Yi[j]) > 0.5:
                fixXijUp(m, i, j)
    for ip in m.I:
        if i in m.Ni[ip] and getLB(m.Yi[ip]) > 0.5:
            fixXijUp(m, ip, i)
    checkYikSum(m, i)


def checkYikSum(m, i):
    if len(m.Yik) == 0:
        return  # Don't introduce type-dependent variables
    sumYikLB = sum(getLB(m.Yik[i, k]) for k in m.K)
    sumYikUB = sum(getUB(m.Yik[i, k]) for k in m.K)
    assert sumYikLB <= sumYikUB
    if areEqual(sumYikUB, 1.0, 0.1) and getLB(m.Yi[i]) > 0.5:
        assert sumYikUB > 0.5
        # Can fix the one remaining Yik to one
        for k in m.K:
            if not m.Yik[i, k].is_fixed():
                fixYikUp(m, i, k)
                break
    elif sumYikUB < 0.5:
        fixYiDown(m, i)


def fixXijklUp(m, i, j, k, l, checkCikl=True):
    if m.Xijkl[i, j, k, l].is_fixed():
        assert value(m.Xijkl[i, j, k, l]) > 0.5
        return
    for kp in m.K:
        for lp in m.K:
            if (
                kp != k
                and lp != l
                and m.Xijkl[i, j, kp, lp].is_fixed()
                and value(m.Xijkl[i, j, kp, lp]) > 0.5
            ):
                raise ValueError("Tried to set multiple (k,l) in Xijkl to 1")
    m.Xijkl[i, j, k, l].fix(1)
    fixYikUp(m, i, k)
    fixYikUp(m, j, l)
    fixXijUp(m, i, j)
    for kp in m.K:
        for lp in m.K:
            if kp != k and lp != l:
                fixXijklDown(m, i, j, kp, lp)
    if checkCikl:
        checkCiklBounds(m, i, k, l)


def fixXijUp(m, i, j, checkCi=True):
    if m.Xij[i, j].is_fixed():
        assert value(m.Xij[i, j]) > 0.5
        return
    m.Xij[i, j].fix(1)
    fixYiUp(m, i)
    fixYiUp(m, j)
    if checkCi:
        checkCiBounds(m, i)


def checkCiklBounds(m, i, k, l):
    sumXijklLB = sum(
        1
        for j in m.Ni[i]
        if (
            j is not None
            and m.Xijkl[i, j, k, l].is_fixed()
            and value(m.Xijkl[i, j, k, l]) > 0.5
        )
    )
    sumXijklUB = sum(1 for j in m.Ni[i] if j is not None) - sum(
        1
        for j in m.Ni[i]
        if (
            j is not None
            and m.Xijkl[i, j, k, l].is_fixed()
            and value(m.Xijkl[i, j, k, l]) < 0.5
        )
    )
    assert sumXijklLB <= sumXijklUB
    if sumXijklUB < getLB(m.Cikl[i, k, l]):
        raise ValueError("Fixing basic vars implied an infeasible problem")
    if sumXijklLB > getUB(m.Cikl[i, k, l]):
        raise ValueError("Fixing basic vars implied an infeasible problem")
    if getUB(m.Cikl[i, k, l]) > sumXijklUB:
        m.Cikl[i, k, l].setub(sumXijklUB)
    if getLB(m.Cikl[i, k, l]) < sumXijklLB:  # Tighten Cikl bounds
        m.Cikl[i, k, l].setlb(sumXijklLB)  # Tighten Cikl bounds
    if areEqual(m.Cikl[i, k, l].lb, m.Cikl[i, k, l].ub, DBL_TOL):
        m.Cikl[i, k, l].fix(m.Cikl[i, k, l].lb)  # Fix if bounds are closed
    if areEqual(getLB(m.Cikl[i, k, l]), sumXijklUB, 0.1):
        # fix all unfixed Xijkl to one
        for j in m.Ni[i]:
            if j is not None and not m.Xijkl[i, j, k, l].is_fixed():
                fixXijklUp(m, i, j, k, l, checkCikl=False)
    if areEqual(getUB(m.Cikl[i, k, l]), sumXijklLB, 0.1):
        # fix all unfixed Xijkl to zero
        for j in m.Ni[i]:
            if j is not None and not m.Xijkl[i, j, k, l].is_fixed():
                fixXijklDown(m, i, j, k, l, checkCikl=False)


def checkCiBounds(m, i):
    sumXijLB = sum(
        1
        for j in m.Ni[i]
        if (j is not None and m.Xij[i, j].is_fixed() and value(m.Xij[i, j]) > 0.5)
    )
    sumXijUB = sum(1 for j in m.Ni[i] if j is not None) - sum(
        1
        for j in m.Ni[i]
        if (j is not None and m.Xij[i, j].is_fixed() and value(m.Xij[i, j]) < 0.5)
    )
    assert sumXijLB <= sumXijUB
    if sumXijUB < getLB(m.Ci[i]):
        raise ValueError("Fixing basic vars implied an infeasible problem")
    if sumXijLB > getUB(m.Ci[i]):
        raise ValueError("Fixing basic vars implied an infeasible problem")
    if getUB(m.Ci[i]) > sumXijUB:
        m.Ci[i].setub(sumXijUB)
    if getLB(m.Ci[i]) < sumXijLB:  # Tighten Ci bounds
        m.Ci[i].setlb(sumXijLB)  # Tighten Ci bounds
    if areEqual(m.Ci[i].lb, m.Ci[i].ub, DBL_TOL):
        m.Ci[i].fix(m.Ci[i].lb)  # Fix if bounds are closed
    if areEqual(getLB(m.Ci[i]), sumXijUB, 0.1):
        # fix all unfixed Xij to one
        for j in m.Ni[i]:
            if j is not None and not m.Xij[i, j].is_fixed():
                fixXijUp(m, i, j, checkCi=False)
    if areEqual(getUB(m.Ci[i]), sumXijLB, 0.1):
        # fix all unfixed Xij to zero
        for j in m.Ni[i]:
            if j is not None and not m.Xij[i, j].is_fixed():
                fixXijDown(m, i, j, checkCi=False)


def fixYikDown(m, i, k):
    if m.Yik[i, k].is_fixed():
        assert value(m.Yik[i, k]) < 0.5
        return
    m.Yik[i, k].fix(0)
    for j in m.Ni[i]:
        if j is not None:
            for l in m.K:
                fixXijklDown(m, i, j, k, l)
    checkYikSum(m, i)
    for l in m.K:
        fixCiklDown(m, i, k, l)


def fixYiDown(m, i):
    if m.Yi[i].is_fixed():
        assert value(m.Yi[i]) < 0.5
        return
    m.Yi[i].fix(0)
    for k in m.K:
        fixYikDown(m, i, k)
    for j in m.Ni[i]:
        if j is not None:
            fixXijDown(m, i, j)
    fixCiDown(m, i)


def fixXijklDown(m, i, j, k, l, checkCikl=True):
    if m.Xijkl[i, j, k, l].is_fixed():
        assert value(m.Xijkl[i, j, k, l]) < 0.5
        return
    m.Xijkl[i, j, k, l].fix(0)
    if checkCikl:
        checkCiklBounds(m, i, k, l)


def fixXijDown(m, i, j, checkCi=True):
    if m.Xij[i, j].is_fixed():
        assert value(m.Xij[i, j]) < 0.5
        return
    m.Xij[i, j].fix(0)
    if checkCi:
        checkCiBounds(m, i)


def fixCiklDown(m, i, k, l):
    if m.Cikl[i, k, l].is_fixed():
        assert value(m.Cikl[i, k, l]) < 0.5
        return
    m.Cikl[i, k, l].fix(0)
    checkCiklBounds(m, i, k, l)


def fixCiDown(m, i):
    if m.Ci[i].is_fixed():
        assert value(m.Ci[i]) < 0.5
        return
    m.Ci[i].fix(0)
    for k in m.K:
        for l in m.L:
            fixCiklDown(m, i, k, l)
    checkCiBounds(m, i)


def fixZicUp(m, i, c):
    if m.Zic[i, c].is_fixed():
        assert value(m.Zic[i, c]) > 0.5
        return
    m.Zic[i, c].fix(1)
    for cp in m.C:
        if cp != c:
            fixZicDown(m, i, cp)
    for l, j in enumerate(m.Ni[i]):
        if j is not None:
            for k in m.K:
                if m.Confs[c][l] == k:
                    fixYikUp(m, j, k)
                else:
                    fixYikDown(m, j, k)


def fixZicDown(m, i, c):
    if m.Zic[i, c].is_fixed():
        assert value(m.Zic[i, c]) < 0.5
        return
    m.Zic[i, c].fix(0)


def setDesignFromYik(D, m, blnSetNoneOtherwise=True):
    """Set the elements of a Design according to variables.

    Uses the Pyomo Set m.Atoms to determine which building blocks
    correspond to variable values.

    Args:
        D (Design): Design that matches the Canvas that the model
            was constructed with.
        m (ConcreteModel): Pyomo model with a solution avialable and an
            Atoms set defined.
        blnSetNoneOtherwise (bool): Flag to control whether the contents
            should be erased before assignment. (Default value = True)

    Returns:
        None.

    """
    if blnSetNoneOtherwise:
        D.setContents(None)
    for i, k in m.Yik.keys():
        if value(m.Yik[i, k]) > 0.5:
            D.setContent(i, k)


def setDesignFromYi(D, m, blnSetNoneOtherwise=True):
    """Set the elements of a Design according to variables.

    Uses the Pyomo Set m.Atoms to determine which building block
    corresponds to variable values. Assumes that m.Atoms[0]
    corresponds to 1 and None corresponds to 0.

    Args:
        D (Design): Design that matches the Canvas that the model
            was constructed with.
        m (ConcreteModel): Pyomo model with a solution avialable and an
            Atoms set defined.
        blnSetNoneOtherwise (bool): Flag to control whether the contents
            should be erased before assignment. (Default value = True)

    Returns:
        None.

    """
    if blnSetNoneOtherwise:
        D.setContents(None)
    for i in m.Yi.keys():
        if value(m.Yi[i]) > 0.5:
            D.setContent(i, m.Atoms[0])
        else:
            D.setContent(i, None)


def setDesignFromModel(D, m, blnSetNoneOtherwise=True):
    """Set the elements of a Design according to variables.

    Chooses between setDesignFromYik and setDesignFromYi.
    See those functions for more information.

    Args:
        D (Design): Design that matches the Canvas that the model
            was constructed with.
        m (ConcreteModel): Pyomo model with a solution avialable and an
            Atoms set defined.
        blnSetNoneOtherwise (bool): Flag to control whether the contents
            should be erased before assignment. (Default value = True)

    Returns:
        None.

    Raises:
        NotImplementedError: If it is not clear how to decide between
            competing choices.

    """
    if len(m.Yik) > 0:
        setDesignFromYik(D, m, blnSetNoneOtherwise=blnSetNoneOtherwise)
    elif len(m.Yi) > 0:
        setDesignFromYi(D, m, blnSetNoneOtherwise=blnSetNoneOtherwise)
    else:
        raise NotImplementedError(
            "User should define MaterialDescriptor Yi or Yik explicitly using valid DescriptorRule"
        )


def _validYiGivenD(m, D):
    for i in m.Yi.keys():
        if value(m.Yi[i]) > 0.5:
            if D.Contents[i] is None:
                # CASE: Atom should be present in Design but was missing
                logging.debug(
                    "_validYiGivenD value(m.Yi[{}])={} D.Contents[{}]={}".format(
                        i, value(m.Yi[i]), i, D.Contents[i]
                    )
                )
                return False
        else:
            if D.Contents[i] is not None:
                # CASE: Atom should not present in Design but was found
                logging.debug(
                    "_validYiGivenD value(m.Yi[{}])={} D.Contents[{}]={}".format(
                        i, value(m.Yi[i]), i, D.Contents[i]
                    )
                )
                return False
    return True


def _validYikGivenD(m, D):
    for i, k in m.Yik.keys():
        if value(m.Yik[i, k]) > 0.5:
            if D.Contents[i] != k:
                # CASE: Atom of type k should be present, but was not found in Design
                logging.debug(
                    "_validYikGivenD value(m.Yik[{}])={} Contents[i]={}".format(
                        (i, k), value(m.Yik[i, k]), i, D.Contents[i]
                    )
                )
                return False
        else:
            if D.Contents[i] == k:
                # CASE: Atom of type k should not be present, but was found in Design
                logging.debug(
                    "_validYikGivenD value(m.Yik[{}])={} Contents[i]={}".format(
                        (i, k), value(m.Yik[i, k]), i, D.Contents[i]
                    )
                )
                return False
    return True


def _validXijGivenD(m, D):
    for i, j in m.Xij.keys():
        if value(m.Xij[i, j]) > 0.5:
            if D.Contents[i] is None or D.Contents[j] is None:
                # CASE: Bond should be present but missing in Design
                logging.debug(
                    "_validXijGivenD value(m.Xij[{}])={} D.Contents[{}]={} D.Contents[{}]={}".format(
                        (i, j), value(m.Xij[i, j]), i, D.Contents[i], j, D.Contents[j]
                    )
                )
                return False
        else:
            if D.Contents[i] is not None and D.Contents[j] is not None:
                # CASE: No bond should be present, but was found in Design
                logging.debug(
                    "_validXijGivenD value(m.Xij[{}])={} D.Contents[{}]={} D.Contents[{}]={}".format(
                        (i, j), value(m.Xij[i, j]), i, D.Contents[i], j, D.Contents[j]
                    )
                )
                return False
    return True


def _validXijklGivenD(m, D):
    for i, j, k, l in m.Xijkl.keys():
        if value(m.Xijkl[i, j, k, l]) > 0.5:
            if D.Contents[i] != k or D.Contents[j] != l:
                # CASE: Bond of type k,l should be present, but missing in Design
                logging.debug(
                    "_validXijklGivenD value(m.Xijkl[{}])={} D.Contents[{}]={} D.Contents[{}]={}".format(
                        (i, j, k, l),
                        value(m.Xijkl[i, j, k, l]),
                        i,
                        D.Contents[{}],
                        j,
                        D.Contents[{}],
                    )
                )
                return False
        else:
            if D.Contents[i] == k and D.Contents[j] == l:
                # CASE: Bond of type k,l should not be prese,t but was found in Design
                logging.debug(
                    "_validXijklGivenD value(m.Xijkl[{}])={} D.Contents[{}]={} D.Contents[{}]={}".format(
                        (i, j, k, l),
                        value(m.Xijkl[i, j, k, l]),
                        i,
                        D.Contents[{}],
                        j,
                        D.Contents[{}],
                    )
                )
                return False
    return True


def _validCiGivenD(m, D, DBL_TOL=1e-5):
    for i in m.Ci.keys():
        Count = 0
        if D.Contents[i] is not None:
            for j in m.Ni[i]:
                if j is not None and D.Contents[j] is not None:
                    Count += 1
        if not areEqual(Count, value(m.Ci[i]), DBL_TOL):
            # CASE: Count of non-void neighbors did not match
            logging.debug("_validCiGivenD value(m.Ci[{}])={}".format(i, value(m.Ci[i])))
            return False
    return True


def _validCiklGivenD(m, D, DBL_TOL=1e-5):
    for i, k, l in m.Cikl.keys():
        Count = 0
        if D.Contents[i] == k:
            for j in m.Ni[i]:
                if j is not None and D.Contents[j] == l:
                    Count += 1
        if not areEqual(Count, value(m.Cikl[i, k, l]), DBL_TOL):
            # CASE: Count of bonds of type k,l did not mach
            logging.debug(
                "_validCiklGivenD value(m.Cikl[{}])={}".format(
                    (i, k, l), value(m.Cikl[i, k, l])
                )
            )
            return False
    return True


def _validVarsGivenD(m, D):
    if len(m.Yi) > 0 and not _validYiGivenD(m, D):
        return False
    if len(m.Yik) > 0 and not _validYikGivenD(m, D):
        return False
    if len(m.Xij) > 0 and not _validXijGivenD(m, D):
        return False
    if len(m.Xijkl) > 0 and not _validXijklGivenD(m, D):
        return False
    if len(m.Ci) > 0 and not _validCiGivenD(m, D):
        return False
    if len(m.Cikl) > 0 and not _validCiklGivenD(m, D):
        return False
    return True


def validModelSoln(m, D):
    """Assess if a Design matches variable values.

    This function checks that a building blocks in a Design
    are consistent with information from a model. Importantly,
    it is one-directional in that we do *not* check if the
    model has all of the corresponding variables present. It
    is not generally clear how to validate formulated models, but it
    is clear how to validate Designs given a model.

    Args:
        m (ConcreteModel): Pyomo model with a solution avialable.
        D (Design): Design with Contents set.

    Returns:
        (bool) True if Design matches the information in the model
            variables.

    """
    if not _validVarsGivenD(m, D):
        return False
    return True


# ================================================
# ==========   COMMON MODEL ELEMENTS    ==========
# ================================================
# --- COMMON CONSTRAINTS
def addConsBoundDescriptorsWithImpl(
    m,
    ConName,
    Descs,
    Impls,
    LBs=None,
    UBs=None,
    DefaultBigM=DEFAULT_BIG_M,
    blnWarnDefaultBigM=True,
):
    """Add implication constraints to a Pyomo model.

    The form of the constraint is:
    LB - M*(1-Impl) <= Desc <= UB + M*(1-Impl)

    The descriptors, implication variables, and bounds are assumed to
    be indexed over the Canvas. Big-M values are identified by
    interval arithmetic. In case this is not possible, default big-M
    values can be provided.

    Args:
        m (ConcreteModel): Pyomo model to attach constraint to.
        ConName (str): Name of new constraint to attach to the model.
        Descs (list<Expr>): Pyomo expressions to bound with implication.
        Impls (list<Var>/IndexedVar): Pyomo implication variables.
        LBs (list<float/Param>): Optional, lower bounds to enforce if
            implication is active. (Default value = None)
        UBs (list<float/Param>): Optional, upper bounds to enforce if
            implication is active. (Default value = None)
        DefaultBigM (float): Optional, big-M value to use if one can
            not be found from interval arithmetic.
            (Default value = DEFAULT_BIG_M)
        blnWarnDefaultBigM (bool): Optional, flag to enable warnings
            if the default big-M value is used. (Default value = True)

    Returns:
        None.

    """
    if LBs is not None:
        ConNameLB = ConName + "_ImplLB"
        setattr(m, ConNameLB, Constraint(m.I))
        ConLB = getattr(m, ConNameLB)
        for i in m.I:
            ExprLB = getLB(Descs[i])
            if ExprLB is None:
                if blnWarnDefaultBigM:
                    logging.warning(
                        "addConsBoundDescriptorsWithImpl used a DefaultBigM"
                    )
                MLB = -DefaultBigM
            else:
                MLB = ExprLB - LBs[i]
            ConLB.add(index=i, expr=(MLB * (1 - Impls[i]) <= Descs[i] - LBs[i]))

    if UBs is not None:
        ConNameUB = ConName + "_ImplUB"
        setattr(m, ConNameUB, Constraint(m.I))
        ConUB = getattr(m, ConNameUB)
        for i in m.I:
            ExprUB = getUB(Descs[i])
            if ExprUB is None:
                if blnWarnDefaultBigM:
                    logging.warning(
                        "addConsBoundDescriptorsWithImpl used a DefaultBigM"
                    )
                MUB = DefaultBigM
            else:
                MUB = ExprUB - UBs[i]
            ConUB.add(index=i, expr=(MUB * (1 - Impls[i]) >= Descs[i] - UBs[i]))


def addConsIndFromDescriptors(
    m,
    ConName,
    Descs,
    Inds,
    LBs=None,
    UBs=None,
    Eps=DEFAULT_EPS,
    DefaultBigM=DEFAULT_BIG_M,
    blnWarnDefaultBigM=True,
):
    """Add indicator constraints to a Pyomo model.

    The form of the constraint are:
        Desc >= UB + Eps - M*(UBInd)
        Desc <= LB - Eps + M*(LBInd)
        Ind >= LBInd + UBInd - 1

    The descriptors, indication variables, and bounds are assumed to
    be indexed over the Canvas. Big-M values are identified by
    interval arithmetic. In case this is not possible, default big-M
    values can be provided.

    Args:
        m (ConcreteModel): Pyomo model to attach constraint to.
        ConName (str): Name of new constraint to attach to the model.
        Descs (list<Expr>): Pyomo expressions to bound with implication.
        Inds (list<Var>/IndexedVar): Pyomo indicator variables.
        LBs (list<float/Param>): Optional, lower bounds to enforce if
            implication is active. (Default value = None)
        UBs (list<float/Param>): Optional, upper bounds to enforce if
            implication is active. (Default value = None)
        Eps (float): Optional, tolerance to use for modeling
            constraint equality.
        DefaultBigM (float): Optional, big-M value to use if one can
            not be found from interval arithmetic.
            (Default value = DEFAULT_BIG_M)
        blnWarnDefaultBigM (bool): Optional, flag to enable warnings
            if the default big-M value is used. (Default value = True)

    Returns:
        None.

    """
    if LBs is not None and UBs is not None:
        ConNameLink = ConName + "_Link"
        VarNameLBInd = ConName + "_LBIndi"
        VarNameUBInd = ConName + "_UBIndi"
        setattr(m, ConNameLink, Constraint(m.I))
        setattr(m, VarNameLBInd, Var(m.I, bounds=Binary, dense=False))
        setattr(m, VarNameUBInd, Var(m.I, bounds=Binary, dense=False))
        ConLink = getattr(m, ConNameLink)
        IndsLB = getattr(m, VarNameLBInd)
        IndsUB = getattr(m, VarNameUBInd)
        for i in m.I:
            ConLink.add(index=i, expr=(Inds[i] >= IndsLB[i] + IndsUB[i] - 1))
    elif LBs is not None:
        IndsLB = Inds
    elif UBs is not None:
        IndsUB = Inds
    else:
        return

    if LBs is not None:
        ConNameLB = ConName + "_LB"
        setattr(m, ConNameLB, Constraint(m.I))
        ConLB = getattr(m, ConNameLB)
        for i in m.I:
            ExprUB = getUB(Descs[i])
            if ExprUB is None:
                if blnWarnDefaultBigM:
                    logging.warning("addConsIndFromDescriptors used a DefaultBigM")
                MLB = DefaultBigM
            else:
                MLB = -LBs[i] + Eps + ExprUB
            ConLB.add(index=i, expr=(Descs[i] <= LBs[i] - Eps + MLB * (IndsLB[i])))
    if UBs is not None:
        ConNameUB = ConName + "_UB"
        setattr(m, ConNameUB, Constraint(m.I))
        ConUB = getattr(m, ConNameUB)
        for i in m.I:
            ExprLB = getLB(Descs[i])
            if ExprLB is None:
                if blnWarnDefaultBigM:
                    logging.warning("addConsIndFromDescriptors used a DefaultBigM")
                MUB = -DefaultBigM
            else:
                MUB = -UBs[i] - Eps + ExprLB
            ConUB.add(index=i, expr=(Descs[i] >= UBs[i] + Eps + MUB * (IndsUB[i])))


def addConsLocalBudgets(m, ConName, Vs, Subsets=None, BudgetLBs=None, BudgetUBs=None):
    """Add budget constraints to a Pyomo model.

    The form of the constraints are:
        BudgetLB <= sum(Vs[j] for j in Subsets[i]) <= BudgetUB

    The variables, subsets, and budget bounds are assumed to be
    indexed over Canvas locations. The budget is written over
    local neighborhoods defined in Subsets. If the Subsets keyword
    is not provided, then the Canvas neighborhoods are used.

    Args:
        m (ConcreteModel): Pyomo model to attach constraint to.
        ConName (str): Name of new constraint to attach to the model.
        Vs (list<Var/Expr>): Pyomo expressions to add up in a local budget.
        Subsets: Optional, Neighborhoods to sum over. If not provided,
            the Canvas Neighborhoods are used. (Default value = None)
        BudgetLBs: Optional, lower bounds on the local summations to
            enforce. (Default value = None)
        BudgetUBs: Optional, upper bounds on the local summations to
            enforce. (Default value = None)

    Returns:
        None.

    """
    if BudgetLBs is None and BudgetUBs is None:
        return
    if BudgetLBs is None:
        BudgetLBs = [None] * len(m.I)
    elif BudgetUBs is None:
        BudgetUBs = [None] * len(m.I)
    if Subsets is None:
        Subsets = m.Ni
    else:
        assert len(Subsets) == len(m.I)
    setattr(m, ConName, Constraint(m.I))
    Con = getattr(m, ConName)
    for i in m.I:
        if Subsets[i] is not None:
            Expr = 0
            for j in Subsets[i]:
                if j is not None:
                    Expr += Vs[j]
            Con.add(index=i, expr=(BudgetLBs[i], Expr, BudgetUBs[i]))


def addConGlobalBudget(m, ConName, Vs, Subset=None, BudgetLB=None, BudgetUB=None):
    """Add budget constraints to a Pyomo model.

    The form of the constraint is:
        BudgetLB <= sum(Vs[i] for i in Subset) <= BudgetUB

    The variables are assumed to be indexed over Canvas locations.
    In contrast to addConLocalBudgets, this function applies a single
    budget constraint written over the entire Canvas. If the Subsets
    keyword is provided, then a subset of Canvas locations are added.

    Args:
        m (ConcreteModel): Pyomo model to attach constraint to.
        ConName (str): Name of new constraint to attach to the model.
        Vs (list<Var/Expr>): Pyomo expressions to add up in a budget.
        Subset: Optional, indices to sum over. If not provided, the
            Canvas indices are used. (Default value = None)
        BudgetLB: Optional, lower bound on the summation to enforce.
            (Default value = None)
        BudgetUB: Optional, upper bound on the summation to enforce.
            (Default value = None)

    Returns:
        None.

    """
    if BudgetLB is None and BudgetUB is None:
        return
    if Subset is None:
        Subset = m.I
    Expr = sum(Vs[i] for i in Subset)
    setattr(m, ConName, Constraint(expr=(BudgetLB, Expr, BudgetUB)))


def addConsLocalBudgetsByTypes(
    m, ConName, Vs, Subsets=None, Types=None, BudgetLBs=None, BudgetUBs=None
):
    """Add budget constraints to a Pyomo model.

    The form of the constraints are:
        BudgetLBs[k] <= sum(Vs[j,k] for j in Subsets[i]) <= BudgetUBs[k]

    The variables, subsets, and budget bounds are assumed to be
    indexed over Canvas locations. The budget is written over
    local neighborhoods defined in Subsets. If the Subsets keyword
    is not provided, then the Canvas neighborhoods are used.

    Args:
        m (ConcreteModel): Pyomo model to attach constraint to.
        ConName (str): Name of new constraint to attach to the model.
        Vs (dict<Var/Expr>): Pyomo expressions to add up in a local budget.
        Subsets: Optional, Neighborhoods to sum over. If not provided,
            the Canvas Neighborhoods are used. (Default value = None)
        BudgetLBs: Optional, lower bounds on the local summations to
            enforce. (Default value = None)
        BudgetUBs: Optional, upper bounds on the local summations to
            enforce. (Default value = None)

    Returns:
        None.

    """
    if BudgetLBs is None and BudgetUBs is None:
        return
    if Subsets is None:
        Subsets = m.Ni
    else:
        assert len(Subsets) == len(m.I)
    if Types is None:
        Types = m.K
    if BudgetLBs is None:
        BudgetLBs = {(i, k): None for i in m.I for k in Types}
    elif BudgetUBs is None:
        BudgetUBs = {(i, k): None for i in m.I for k in Types}
    setattr(m, ConName, Constraint(m.I, Types))
    Con = getattr(m, ConName)
    for k in Types:
        for i in m.I:
            if Subsets[i] is not None:
                Expr = 0
                for j in Subsets[i]:
                    if j is not None:
                        Expr += Vs[j, k]
                Con.add(index=(i, k), expr=(BudgetLBs[i, k], Expr, BudgetUBs[i, k]))


def addConsGlobalBudgetsByTypes(
    m, ConName, Vs, Subset=None, Types=None, BudgetLBs=None, BudgetUBs=None
):
    """Add budget constraints to a Pyomo model.

    The form of the constraint is:
        BudgetLBs[k] <= sum(Vs[i,k] for i in Subset) <= BudgetUBs[k]

    The variables are assumed to be indexed over Canvas locations
    and types.
    In contrast to addConLocalBudgetsByTypes, this function applies
    a single budget constraint written over the entire Canvas. If the
    Subsets keyword is provided, then a subset of Canvas locations
    are added.

    Args:
        m (ConcreteModel): Pyomo model to attach constraint to.
        ConName (str): Name of new constraint to attach to the model.
        Vs (dict<k:Var/Expr>): Pyomo expressions to add up in a budget.
        Subset: Optional, indices to sum over. If not provided, the
            Canvas indices are used. (Default value = None)
        BudgetLBs: Optional, lower bound on the summation to enforce.
            (Default value = None)
        BudgetUBs: Optional, upper bound on the summation to enforce.
            (Default value = None)

    Returns:
        None.

    """
    if BudgetLBs is None and BudgetUBs is None:
        return
    if Subset is None:
        Subset = m.I
    if Types is None:
        Types = m.K
    if BudgetLBs is None:
        BudgetLBs = {k: None for k in Types}
    elif BudgetUBs is None:
        BudgetUBs = {k: None for k in Types}
    setattr(m, ConName, Constraint(Types))
    Con = getattr(m, ConName)
    for k in Types:
        Expr = sum(Vs[i, k] for i in Subset)
        Con.add(index=k, expr=(BudgetLBs[k], Expr, BudgetUBs[k]))


def addConsZicMutExc(m):
    """Add constraints for mutually exclusive conformations.

    Args:
        m (ConcreteModel): Pyomo model to attach constraints to.

    Returns:
        None.

    """

    def _ruleZicMutExc(m, i):
        return sum(m.Zic[i, c] for c in m.C) <= 1

    m.AssignZicMutExc = Constraint(m.I, rule=_ruleZicMutExc)


def addConsZicColExh(m):
    """
    Add constraints for collectively exhaustive conformations.

    Args:
        m (ConcreteModel): Pyomo model to attach constraints to.

    Returns:
        None.

    """

    def _ruleZicColExh(m, i):
        return sum(m.Zic[i, c] for c in m.C) >= 1

    m.AssignZicColExh = Constraint(m.I, rule=_ruleZicColExh)


def addConsZicMutExcColExh(m):
    """
    Add constraints for mutually exclusive, collectively
    exhaustive of conformations.

    Args:
        m (ConcreteModel): Pyomo model to attach constraints to.

    Returns:
        None.

    """

    def _ruleZicMutExcColExh(m, i):
        return sum(m.Zic[i, c] for c in m.C) == 1

    m.AssignZicMutExcColExh = Constraint(m.I, rule=_ruleZicMutExcColExh)


def addConsZicFromYi(m, blnConfsAreMutExc=True, blnConfsAreColExh=False):
    """
    Add constraints for indicating conformations.

    This function produces indicator constraints for the presence of
    conformation 'c' at location 'i' if building blocks in neighboring
    locations match the pattern for conformation 'c'. Neighboring
    locations are defined by the set m.Ni. In this function, any
    types of building block are considered. The variable Yi should be
    equal to 1 when any building block is present.
    This function only writes assignment constraints for m.Zic variables
    previously referenced in the model.

    Args:
        m (ConcreteModel): Pyomo model to attach constraints to.
        blnConfsAreMutExc (bool): Flag to indicate whether to add
            constraints enforcing mutually exclusive conformations.
            (Default value = True)
        blnConfsAreColExh (bool): Flag to indicate whether to add
            constraints enforcing collectively exhaustive conformations.
            NOTE: Since this function only loops over Zic variables
            previously referenced in the model, it is important
            to be carefule when using this flag. (Default value = False)

    Returns:
        None.

    """
    m.AssignZicFromYi1 = Constraint(m.I, m.C, m.I)
    m.AssignZicFromYi2 = Constraint(m.I, m.C)
    for i, c in m.Zic.keys():
        assert len(m.Ni[i]) == len(m.Confs[c])
        RelaxExpr = 0
        for l, j in enumerate(m.Ni[i]):
            if j is not None:
                if m.Confs[c][l] is not None:
                    RelaxExpr += m.Yi[j] - 1
                    m.AssignZicFromYi1.add(
                        index=(i, c, j), expr=(m.Zic[i, c] <= m.Yi[j])
                    )
                else:
                    assert m.Confs[c][l] is None
                    RelaxExpr += -m.Yi[j]
                    m.AssignZicFromYi1.add(
                        index=(i, c, j), expr=(m.Zic[i, c] <= 1 - m.Yi[j])
                    )
        m.AssignZicFromYi2.add(index=(i, c), expr=(m.Zic[i, c] >= 1 + RelaxExpr))
    if blnConfsAreMutExc and blnConfsAreColExh:
        addConsZicMutExcColExh(m)
    elif blnConfsAreMutExc:
        addConsZicMutExc(m)
    elif blnConfsAreColExh:
        addConsZicColExh(m)


def addConsZicFromYiLifted(m, blnConfsAreMutExc=True, blnConfsAreColExh=False):
    """Add constraints for indicating conformations.

    See documentation for addConsZicFromYi for more information.
    This version uses a lifted set of constraints.

    Args:
        m (ConcreteModel): Pyomo model to attach constraints to.
        blnConfsAreMutExc (bool): Flag to indicate whether to add
            constraints enforcing mutually exclusive conformations.
            (Default value = True)
        blnConfsAreColExh (bool): Flag to indicate whether to add
            constraints enforcing collectively exhaustive conformations.
            NOTE: Since this function only loops over Zic variables
            previously referenced in the model, it is important
            to be carefule when using this flag. (Default value = False)

    Returns:
        None.

    """
    m.AssignZicFromYiLifted1 = Constraint(m.I, m.I)
    m.AssignZicFromYiLifted2 = Constraint(m.I, m.I)
    m.AssignZicFromYiLifted3 = Constraint(m.I, m.C)
    iAlreadyProcessed = set()
    for i, _ in m.Zic.keys():
        if i not in iAlreadyProcessed:
            for l, j in enumerate(m.Ni[i]):
                if j is not None:
                    PosZic = sum(
                        m.Zic[i, c]
                        for c in m.Zic[i, :].wildcard_keys()
                        if m.Confs[c][l] is not None
                    )
                    NegZic = sum(
                        m.Zic[i, c]
                        for c in m.Zic[i, :].wildcard_keys()
                        if m.Confs[c][l] is None
                    )
                    if PosZic is not 0:
                        m.AssignZicFromYiLifted1.add(
                            index=(i, j), expr=(PosZic <= m.Yi[j])
                        )
                    if NegZic is not 0:
                        m.AssignZicFromYiLifted2.add(
                            index=(i, j), expr=(NegZic <= 1 - m.Yi[j])
                        )
            iAlreadyProcessed.add(i)
    for i, c in m.Zic.keys():
        RelaxExpr = 0
        for l, j in enumerate(m.Ni[i]):
            if j is not None:
                if m.Confs[c][l] is not None:
                    RelaxExpr += m.Yi[j] - 1
                else:
                    RelaxExpr += -m.Yi[j]
        m.AssignZicFromYiLifted3.add(index=(i, c), expr=(m.Zic[i, c] >= 1 + RelaxExpr))
    if blnConfsAreMutExc and blnConfsAreColExh:
        addConsZicMutExcColExh(m)
    elif blnConfsAreMutExc:
        addConsZicMutExc(m)
    elif blnConfsAreColExh:
        addConsZicColExh(m)


def addConsZicFromYik(m, blnConfsAreMutExc=True, blnConfsAreColExh=False):
    """Add constraints for indicating conformations.

    This function produces indicator constraints for the presence of
    conformation 'c' at location 'i' if building blocks of type 'k'
    in neighboring locations match the pattern for conformation 'c'.
    Neighboring locations are defined by the set m.Ni. Building block
    types are defined in m.Atoms and correspond to the 'k' index.
    The variable Yik should be equal to 1 when a building block of type
    k is present in the conformation.
    This function only writes assignment constraints for m.Zic variables
    previously referenced in the model.

    Args:
        m (ConcreteModel): Pyomo model to attach constraints to.
        blnConfsAreMutExc (bool): Flag to indicate whether to add
            constraints enforcing mutually exclusive conformations.
            (Default value = True)
        blnConfsAreColExh (bool): Flag to indicate whether to add
            constraints enforcing collectively exhaustive conformations.
            NOTE: Since this function only loops over Zic variables
            previously referenced in the model, it is important
            to be carefule when using this flag. (Default value = False)

    Returns:
        None.

    """
    m.AssignZicFromYik1 = Constraint(m.I, m.C, m.I, m.K)
    m.AssignZicFromYik2 = Constraint(m.I, m.C)
    for i, c in m.Zic.keys():
        assert len(m.Ni[i]) == len(m.Confs[c])
        RelaxExpr = 0
        for l, j in enumerate(m.Ni[i]):
            if j is not None:
                if m.Confs[c][l] is None:
                    RelaxExpr += m.Yi[j]
                for k in m.K:
                    if m.Confs[c][l] == k:
                        RelaxExpr += 1 - m.Yik[j, k]
                        m.AssignZicFromYik1.add(
                            index=(i, c, j, k), expr=(m.Zic[i, c] <= m.Yik[j, k])
                        )
                    else:
                        m.AssignZicFromYik1.add(
                            index=(i, c, j, k), expr=(m.Zic[i, c] <= 1 - m.Yik[j, k])
                        )
        m.AssignZicFromYik2.add(index=(i, c), expr=(m.Zic[i, c] >= 1 - RelaxExpr))
    if blnConfsAreMutExc and blnConfsAreColExh:
        addConsZicMutExcColExh(m)
    elif blnConfsAreMutExc:
        addConsZicMutExc(m)
    elif blnConfsAreColExh:
        addConsZicColExh(m)


def addConsZicFromYikLifted(m, blnConfsAreMutExc=True, blnConfsAreColExh=False):
    """Add constraints for indicating conformations.

    See documentation for addConsZicFromYik for more information.
    This version uses a lifted set of constraints.

    Args:
        m (ConcreteModel): Pyomo model to attach constraints to.
        blnConfsAreMutExc (bool): Flag to indicate whether to add
            constraints enforcing mutually exclusive conformations.
            (Default value = True)
        blnConfsAreColExh (bool): Flag to indicate whether to add
            constraints enforcing collectively exhaustive conformations.
            NOTE: Since this function only loops over Zic variables
            previously referenced in the model, it is important
            to be carefule when using this flag. (Default value = False)

    Returns:
        None.

    """
    m.AssignZicFromYikLifted1 = Constraint(m.I, m.I, m.K)
    m.AssignZicFromYikLifted2 = Constraint(m.I, m.I, m.K)
    m.AssignZicFromYikLifted3 = Constraint(m.I, m.C)
    iAlreadyProcessed = set()
    for i, _ in m.Zic.keys():
        if i not in iAlreadyProcessed:
            for l, j in enumerate(m.Ni[i]):
                if j is not None:
                    for k in m.K:
                        PosZic = sum(
                            m.Zic[i, c]
                            for c in m.Zic[i, :].wildcard_keys()
                            if m.Confs[c][l] == k
                        )
                        NegZic = sum(
                            m.Zic[i, c]
                            for c in m.Zic[i, :].wildcard_keys()
                            if m.Confs[c][l] != k
                        )
                        if PosZic is not 0:
                            m.AssignZicFromYikLifted1.add(
                                index=(i, j, k), expr=(PosZic <= m.Yik[j, k])
                            )
                        if NegZic is not 0:
                            m.AssignZicFromYikLifted2.add(
                                index=(i, j, k), expr=(NegZic <= 1 - m.Yik[j, k])
                            )
            iAlreadyProcessed.add(i)
    for i, c in m.Zic.keys():
        RelaxExpr = 0
        for l, j in enumerate(m.Ni[i]):
            if j is not None:
                if m.Confs[c][l] is None:
                    RelaxExpr += m.Yi[j]
                for k in m.K:
                    if m.Confs[c][l] == k:
                        RelaxExpr += 1 - m.Yik[j, k]
        m.AssignZicFromYikLifted3.add(index=(i, c), expr=(m.Zic[i, c] >= 1 - RelaxExpr))
    if blnConfsAreMutExc and blnConfsAreColExh:
        addConsZicMutExcColExh(m)
    elif blnConfsAreMutExc:
        addConsZicMutExc(m)
    elif blnConfsAreColExh:
        addConsZicColExh(m)


# --- COMMON OBJECTIVE FUNCTIONS
def addObjMaxSumZi(m):
    """Add objective for maximizing target sites.

    Args:
        m (ConcreteModel): Pyomo model to add objective to.

    Returns:
        None.

    """
    m.Obj = Objective(expr=sum(m.Zi[i] for i in m.I), sense=maximize)


def addObjMinSumYi(m):
    """Add objective for minimizing building blocks used.

    This objective is mainly useful for debugging models.
    Fix some (supposed) feasible solution, use this objective,
    and see which building blocks were realized.

    Args:
        m (ConcreteModel): Pyomo model to add objective to.

    Returns:
        None.

    """
    m.Obj = Objective(expr=sum(m.Yi[i] for i in m.I), sense=minimize)


def addObjMaxSumZic(m, Alphas=None, ISubset=None, ConfSubset=None):
    """Add objective for maximizing target conformations.

    Args:
        m (ConcreteModel): Pyomo model to add objective to.
        Alphas (list<float/Param>): Optional, weighting to give each
            conformation. If not present, assumes a weight of 1.0.
            (Default value = None)
        ISubset (list<int>/Set): Optional, subset of locations to
            consider. If not present, assumes the entire Canvas.
            (Default value = None)
        ConfSubset: Optional, subset of conformations to consider.
            If not present, assumes all conformations.
            (Default value = None)

    Returns:
        None.

    """
    if ISubset is None:
        ISubset = m.I
    if ConfSubset is None:
        ConfSubset = m.C
    if Alphas is None:
        Alphas = [1] * len(ISubset)
    m.Obj = Objective(
        expr=sum(Alphas[c] * m.Zic[i, c] for i in ISubset for c in ConfSubset),
        sense=maximize,
    )
