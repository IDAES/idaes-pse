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

import re

from pyomo.core.expr.sympy_tools import (
    _functionMap,
    _pyomo_operator_map,
    _configure_sympy,
)

try:
    import sympy

    _configure_sympy(sympy, True)
except:
    pass

from pyomo.environ import ExternalFunction, Var, Expression, value, units as pu
from pyomo.core.base.constraint import _ConstraintData, Constraint
from pyomo.core.base.expression import _ExpressionData
from pyomo.core.base.block import _BlockData
from pyomo.core.expr.visitor import StreamBasedExpressionVisitor
from pyomo.core.expr.numeric_expr import ExternalFunctionExpression
from pyomo.core.expr import current as EXPR, native_types
from pyomo.common.collections import ComponentMap

# TODO<jce> Look into things like sum operator and template expressions


def _add_latex_subscripts(x, s):
    if not "_" in x:
        return "{}_{{ {} }}".format(x, s)
    else:
        return re.sub(r"^(.+_)({.+}|.)", "\\1{{ \\2 ,{} }}".format(s), x)


def deduplicate_symbol(x, v, used):
    """
    Check if x is a duplicated LaTeX symbol if so add incrementing Di subscript

    Args:
        x: symbol string
        v: pyomo object
        used: dictionary of pyomo objects with symbols as keys

    Returns:
        Returns a unique symbol.  If x was not in used keys, returns x,
        otherwise adds exponents to make it unique.
    """
    y = x
    k = 1
    while x in used and id(used[x]) != id(v):
        x = _add_latex_subscripts(y, "D{}".format(k))
        k += 1
        if k > 1000:
            # Either the symbol string is not updating or there are lots of dupes
            break
    used[x] = v
    return x


class PyomoSympyBimap(object):
    """
    This is based on the class of the same name in pyomo.core.base.symbolic, but
    it adds mapping latex symbols to the sympy symbols. This will get you pretty
    equations when using sympy's LaTeX writer.
    """

    def __init__(self):
        self.pyomo2sympy = ComponentMap()
        self.parent_symbol = ComponentMap()
        self.sympy2pyomo = {}
        self.sympy2latex = {}
        self.used = {}
        self.i_var = 0
        self.i_expr = 0
        self.i_func = 0
        self.i = 0

    def getPyomoSymbol(self, sympy_object, default=None):
        return self.sympy2pyomo.get(sympy_object, default)

    def getSympySymbol(self, pyomo_object):
        if pyomo_object in self.pyomo2sympy:
            return self.pyomo2sympy[pyomo_object]
        else:
            return self._add_sympy(pyomo_object)

    def sympyVars(self):
        return self.sympy2pyomo.keys()

    def _add_sympy(self, pyomo_object):
        parent_object = pyomo_object.parent_component()
        if isinstance(parent_object, Var):
            sympy_class = sympy.Symbol
            base_name = "x"
            i = self.i_var
            self.i_var += 1
        elif isinstance(parent_object, Expression):
            sympy_class = sympy.Symbol
            base_name = "u"
            i = self.i_expr
            self.i_expr += 1
        elif isinstance(parent_object, ExternalFunction):
            sympy_class = sympy.Function
            base_name = "func"
            i = self.i_func
            self.i_func += 1
        else:
            raise Exception("Should be Var, Exression, or ExternalFunction")

        if parent_object.is_indexed() and parent_object in self.parent_symbol:
            x = self.parent_symbol[parent_object][0]
            latex_symbol = self.parent_symbol[parent_object][1]
        else:
            x = "{}_{}".format(base_name, i)
            base_latex = getattr(parent_object, "latex_symbol", None)
            if base_latex is not None:
                latex_symbol = deduplicate_symbol(base_latex, parent_object, self.used)
            else:
                latex_symbol = x
        if parent_object.is_indexed():
            if parent_object not in self.parent_symbol:
                self.parent_symbol[parent_object] = (x, latex_symbol)
            idx = pyomo_object.index()
            idxl = idx
            if isinstance(idx, tuple):
                idxl = ",".join(idx)
                idx = "|".join(idx)
            x = "{}[{}]".format(x, idx)
            latex_symbol = _add_latex_subscripts(latex_symbol, idxl)
        sympy_obj = sympy_class(x)
        self.pyomo2sympy[pyomo_object] = sympy_obj
        self.sympy2pyomo[sympy_obj] = pyomo_object
        self.sympy2latex[sympy_obj] = latex_symbol
        return sympy_obj


class Pyomo2SympyVisitor(StreamBasedExpressionVisitor):
    """
    This is based on the class of the same name in pyomo.core.base.symbolic, but
    it catches ExternalFunctions and does not decend into named expressions.
    """

    def __init__(self, object_map):
        super(Pyomo2SympyVisitor, self).__init__()
        self.object_map = object_map

    def exitNode(self, node, values):
        if isinstance(node, ExternalFunctionExpression):
            # catch ExternalFunction
            _op = self.object_map.getSympySymbol(node._fcn)
        else:
            if node.__class__ is EXPR.UnaryFunctionExpression:
                return _functionMap[node._name](values[0])
            _op = _pyomo_operator_map.get(node.__class__, None)
        if _op is None:
            return node._apply_operation(values)
        else:
            return _op(*tuple(values))

    def beforeChild(self, node, child, child_idx):
        # Don't replace native or sympy types
        if type(child) in native_types:
            return False, child
        # We will not descend into named expressions...
        if child.is_expression_type():
            if child.is_named_expression_type():
                # To keep expressions from becoming too crazy complicated
                # treat names expressions like variables.
                return False, self.object_map.getSympySymbol(child)
            else:
                return True, None
        # Replace pyomo variables with sympy variables
        if child.is_potentially_variable():
            return False, self.object_map.getSympySymbol(child)
        # Everything else is a constant...
        return False, value(child)


def sympify_expression(expr):
    """
    Converts Pyomo expressions to sympy expressions.
    This is based on the function of the same name in pyomo.core.base.symbolic.
    The difference between this and the Pymomo is that this one checks if the
    expr argument is a named expression and expands it anyway. This allows named
    expressions to only be expanded if they are the top level object.
    """
    #
    # Create the visitor and call it.
    #
    object_map = PyomoSympyBimap()
    visitor = Pyomo2SympyVisitor(object_map)
    is_expr, ans = visitor.beforeChild(None, expr, None)
    try:  # If I explicitly ask for a named expression then descend into it.
        if expr.is_named_expression_type():
            is_expr = True
    except:
        pass
    if not is_expr:  # and not expr.is_named_expression_type():
        return object_map, ans

    return object_map, visitor.walk_expression(expr)


def _add_docs(object_map, docs, typ, head):
    """
    Adds documentation for a set of pyomo components to a markdown table

    Args:
        object_map (PyomoSympyBimap): Pyomo, sympy, LaTeX mapping
        docs: string containing a markdown table
        typ: the class of objects to document (Var, Expression, ExternalFunction)
        head: a string to used in the sybol table heading for this class of objects
    Returns:
        A new string markdown table with added doc rows.
    """
    docked = set()  # components already documented, mainly for indexed compoents
    whead = True  # write heading before adding first item

    if not isinstance(object_map, (list, tuple)):
        object_map = [object_map]

    for om in object_map:
        for i, sc in enumerate(om.sympyVars()):
            c = om.getPyomoSymbol(sc)
            cdat = c
            c = c.parent_component()  # Document the parent for indexed comps
            if not isinstance(c, typ):
                continue
            if whead:  # add heading if is first entry in this section
                docs += f"**{head}** | **Doc** | **Path** | **UoM**\n"
                docs += ":--- | :--- | :--- | :---\n"
                whead = False
            if id(c) not in docked:
                docked.add(id(c))  # make sure we don't get a line for every index
                try:  # just document the parent of indexed vars
                    s = om.parent_symbol[c][1]
                except KeyError:  # non-indexed vars
                    s = om.sympy2latex[sc]
                docs += "${}$|{}|{}|{}\n".format(s, c.doc, c, pu.get_units(cdat))
    return docs


def to_latex(expr):
    """Return a sympy expression for the given Pyomo expression

    Args:
        expr (Expression): Pyomo expression
    Returns:
        (dict): keys: sympy_expr, a sympy expression; where, markdown string
            with documentation table; latex_expr, a LaTeX string representation
            of the expression.
    """
    object_map, sympy_expr = sympify_expression(expr)

    # This next bit documents the expression, could use a lot of work, but
    # for now it generates markdown tables that are resonably readable in a
    # jupyter notebook.
    docs = "\nSymbol | Doc | Path | UoM\n ---: | --- | --- | ---\n"
    docs = _add_docs(object_map, docs, Var, "Variable")
    docs = _add_docs(object_map, docs, Expression, "Expression")
    docs = _add_docs(object_map, docs, ExternalFunction, "Function")

    # probably should break this up, but this will do for now.
    return {
        "sympy_expr": sympy_expr,
        "where": docs,
        "latex_expr": sympy.latex(sympy_expr, symbol_names=object_map.sympy2latex),
        "object_map": object_map,
    }


def document_constraints(
    comp, doc=True, descend_into=True, fixed_vars=False, to_doc=None
):
    """
    Provides nicely formatted constraint documetntation in markdown format,
    assuming the $$latex math$$ and $latex math$ syntax is supported.

    Args:
        comp: A Pyomo component to document in {_ConstraintData, _ExpressionData,
                _BlockData}.
        doc: True adds a documentation table for each constraint or expression.
                Due to the way symbols are semi-automatiaclly generated, the
                exact symbol definitions may be unique to each constraint or
                expression, if unique LaTeX symbols were not provided
                everywhere in a block.
        descend_into: If True, look in subblocks for constraints.

    Returns:
        A string in markdown format with equations in LaTeX form.
    """
    if to_doc is None:
        to_doc = []
    s = None
    if isinstance(comp, _ExpressionData):
        d = to_latex(comp)
        to_doc.append(d["object_map"])
        try:
            sy = comp.latex_symbol
        except:
            sy = None
        try:
            d["latex_expr"] = comp.latex_nice_expr
        except:
            pass
        if sy is not None:
            if doc:
                s = "$${} = {}$$\n{}".format(sy, d["latex_expr"], d["where"])
            else:
                s = "$${} = {}$$".format(sy, d["latex_expr"])
        else:
            if doc:
                s = "$${}$$\n{}".format(d["latex_expr"], d["where"])
            else:
                s = "$${}$$".format(d["latex_expr"])
    elif isinstance(comp, _ConstraintData):
        d = to_latex(comp.body)
        to_doc.append(d["object_map"])
        try:
            return f"$${comp.latex_nice_expr}$$"
        except:
            pass
        if comp.upper != comp.lower:
            if doc:
                s = "$${} \le {}\le {}$$\n{}".format(
                    comp.lower, d["latex_expr"], comp.upper, d["where"]
                )
            else:
                s = "$${} \le {}\le {}$$".format(
                    comp.lower, d["latex_expr"], comp.upper
                )
        else:
            if doc:
                s = "$${} = {}$$\n{}".format(comp.lower, d["latex_expr"], d["where"])
            else:
                s = "$${} = {}$$".format(comp.lower, d["latex_expr"])
    elif isinstance(comp, _BlockData):
        cs = []
        for c in comp.component_data_objects(Constraint, descend_into=descend_into):
            if not c.active:
                continue
            cs.append("\n**Constraint:** {}\n".format(c))
            cs.append(document_constraints(c, doc, to_doc=to_doc))
        for c in comp.component_data_objects(Expression, descend_into=descend_into):
            if not c.active:
                continue
            cs.append("\n**Expression:** {}\n".format(c))
            cs.append(document_constraints(c, doc, to_doc=to_doc))
        if fixed_vars:
            for c in comp.component_data_objects(Var, descend_into=descend_into):
                if not c.fixed:
                    continue
                cs.append(f"**Fixed Var:** {c} = {value(c)}")
                try:
                    sy = c.latex_symbol
                except:
                    sy = None
                if sy is not None:
                    cs.append(
                        "$${} = {} \\text{{ {} }}$$".format(
                            sy, value(c), pu.get_units(c)
                        )
                    )
        s = "\n".join(cs)
        docs = "\n"
        docs_var = _add_docs(to_doc, docs, Var, "Variable")
        docs_expr = _add_docs(to_doc, docs, Expression, "Expression")
        s = "\n".join([s, docs_var, docs_expr])
    return s
