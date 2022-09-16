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

__author__ = "John Eslick"

"""Transformation to replace variables with other variables."""
import pyomo.environ as pyo
from pyomo.core.base.transformation import TransformationFactory
from pyomo.core.plugins.transform.hierarchy import NonIsomorphicTransformation
from pyomo.core.expr import current as EXPR
from pyomo.contrib.fbbt.fbbt import compute_bounds_on_expr
from pyomo.repn import generate_standard_repn
import idaes.logger as idaeslog


_log = idaeslog.getLogger(__name__)


@TransformationFactory.register(
    "simple_equality_eliminator",
    doc="Eliminate simple equalities in the form a*x + b*y = c or a*x = b",
)
class SimpleEqualityEliminator(NonIsomorphicTransformation):
    def _get_subs(self, instance):
        subs = {}  # Substitute one var for another from a * x + b * y + c = 0
        subs_map = {}  # id -> var
        fixes = []  # fix a variable from a * x + c = 0
        cnstr = set()  # constraints to deactivate
        rset = set()  # variables used in a sub or fixed
        for c in instance.component_data_objects(pyo.Constraint, active=True):
            if (
                pyo.value(c.lower) is not None
                and pyo.value(c.lower) == pyo.value(c.upper)
                and c.body.polynomial_degree() == 1
            ):
                repn = generate_standard_repn(c.body)
                if len(repn.nonlinear_vars) != 0 or len(repn.quadratic_vars) != 0:
                    continue
                if len(repn.linear_vars) > 2:
                    continue
                elif len(repn.linear_vars) < 1:
                    _log.warning("Constraint with no vars {}: {}".format(c, c.expr))
                    continue
                b0 = repn.constant - pyo.value(c.upper)
                v0 = repn.linear_vars[0]
                a0 = repn.linear_coefs[0]
                if id(v0) in rset:
                    continue
                elif len(repn.linear_vars) == 1:
                    fixes.append((v0, -b0 / a0))
                    rset.add(id(v0))
                    cnstr.add(c)
                    continue

                v1 = repn.linear_vars[1]
                a1 = repn.linear_coefs[1]
                if id(v1) in rset:
                    continue

                subs[id(v0)] = -(b0 + a1 * v1) / a0
                cnstr.add(c)
                rset.add(id(v0))
                rset.add(id(v1))

                if self.reversible:
                    subs_map[id(v0)] = v0

                _log.debug("Sub: {} = {}".format(v0, subs[id(v0)]))

                # Use the tightest set of bounds from v0 and v1
                lb, ub = compute_bounds_on_expr(-(b0 + a0 * v0) / a1)

                if lb is not None:
                    if v1.lb is None or lb > v1.lb:
                        v1.setlb(lb)
                if ub is not None:
                    if v1.ub is None or ub < v1.ub:
                        v1.setub(ub)

        return subs, cnstr, fixes, subs_map

    def _apply_to(self, instance, max_iter=5, reversible=True):
        """
        Apply the transformation.  This is called by ``apply_to`` in the
        superclass, and should not be called directly.  ``apply_to`` takes the
        same arguments.

        Args:
            instance: A block or model to apply the transformation to

        Returns:
            None
        """
        self.reversible = reversible
        if reversible:
            self._instance = instance
            self._all_subs = []
            self._subs_map = {}
            self._expr_map = {}
            self._all_fixes = []
            self._all_deactivate = []
            self._original = {}

            # The named expressions could be changed as a side effect of the
            # constraint expression replacements, so for maximum saftey, just
            # store all the expressions for Expressions
            for c in instance.component_data_objects(
                pyo.Expression,
                descend_into=True,
            ):
                self._original[id(c)] = c.expr
                self._expr_map[id(c)] = c

        nr_tot = 0
        # repeat elimination until no more can be eliminated or hit max_iter
        for i in range(max_iter):
            subs, cnstr, fixes, subs_map = self._get_subs(instance)

            if reversible:
                self._all_fixes += fixes
                self._all_deactivate += cnstr
                self._all_subs.append(subs)
                self._subs_map.update(subs_map)

            nr = len(cnstr)
            if nr == 0:
                break
            nr_tot += nr

            for c in cnstr:  # deactivate constraints that aren't needed
                c.deactivate()
            for v in fixes:  # fix variables that can be fixed
                v[0].fix(v[1])

            # Do replacements in Expressions, Constraints, and Objectives
            # where one var is replaced with a linear expression containing
            # another
            vis = EXPR.ExpressionReplacementVisitor(
                substitute=subs,
                descend_into_named_expressions=True,
                remove_named_expressions=False,
            )
            for c in instance.component_data_objects(
                (pyo.Constraint, pyo.Objective), descend_into=True, active=True
            ):
                if id(c) not in self._original and reversible:
                    self._original[id(c)] = c.expr
                    self._expr_map[id(c)] = c
                c.set_value(expr=vis.walk_expression(c.expr))

        _log.info("Eliminated {} variables and constraints".format(nr_tot))

    @staticmethod
    def get_logger(level=None):
        """Return the logger for this, so its a little easier to get the debugging
        output, by changing the level.
        """
        if level is not None:
            _log.setLevel(level)
        return _log

    def revert(self):
        """Revert model to pretransformation state, using substitutions to
        calcualte values of varaibles that were removed from the problem. This
        applies to the last reversible transformation performed with this object.
        """
        try:
            instance = self._instance
        except AttributeError:
            _log.warning("Nothing to revert.")

        for c in self._all_deactivate:
            c.activate()
        for c in self._all_fixes:
            c[0].unfix()
        for cid in self._original:
            c = self._expr_map[cid]
            c.set_value(expr=self._original[cid])

        # The problem should be back, now fill in values for the variables that
        # were removed
        for subs in reversed(self._all_subs):
            for sid in subs:
                self._subs_map[sid].value = pyo.value(subs[sid])

        del self._instance
        del self._all_subs
        del self._all_fixes
        del self._all_deactivate
        del self._original
