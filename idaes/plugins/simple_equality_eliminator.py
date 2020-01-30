##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################

__author__ = "John Eslick"

"""Transformation to replace variables with other variables."""
import pyomo.environ as pyo
from pyomo.core.base.plugin import TransformationFactory
from pyomo.core.plugins.transform.hierarchy import NonIsomorphicTransformation
from pyomo.core.expr import current as EXPR
from pyomo.repn import generate_standard_repn
import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)


@TransformationFactory.register(
    'simple_equality_eliminator',
    doc="Eliminate simple equalities in the form a*x + b*y = c or a*x = b")
class SimpleEqualityElimninator(NonIsomorphicTransformation):

    def _get_subs(self, instance):
        subs = {} # Substitute on var for another from a * x + b * y + c = 0
        fixes = [] # fix a variable from a * x + c = 0
        cnstr = set() # constraints ro deactivate
        rset = set() # varaibles used in a sub or fixed
        for c in instance.component_data_objects(pyo.Constraint, active=True):
            if (
                pyo.value(c.lower) is not None and
                pyo.value(c.lower) == pyo.value(c.upper) and
                c.body.polynomial_degree() == 1
            ):
                repn = generate_standard_repn(c.body)
                assert len(repn.nonlinear_vars) == 0
                assert len(repn.quadratic_vars) == 0
                if len(repn.linear_vars) > 2:
                    continue
                elif len(repn.linear_vars) < 1:
                    _log.warning("Constraint with no vars {}: {}".format(c, c.expr))
                    continue
                b0 = repn.constant - pyo.value(c.upper)
                v0 = repn.linear_vars[0]
                a0 = repn.linear_coefs[0]
                if len(repn.linear_vars) == 1:
                    if id(v0) in rset:
                        continue
                    fixes.append((v0, -b0/a0))
                    rset.add(id(v0))
                    cnstr.add(c)
                    continue
                v1 = repn.linear_vars[1]
                a1 = repn.linear_coefs[1]

                if id(v0) not in rset and id(v1) not in rset:
                    k = v0
                    e = -(b0 + a1 * v1) / a0
                    v = v1
                else:
                    continue

                subs[id(k)] = e
                cnstr.add(c)
                rset.add(id(v0))
                rset.add(id(v1))

                _log.debug("Sub: {} = {}".format(k, e))
        return subs, cnstr, fixes


    def _apply_to(self, instance):
        """
        Apply the transformation.  This is called by ``apply_to`` in the
        superclass, and should not be called directly.  ``apply_to`` takes the
        same arguments.

        Args:
            instance: A block or model to apply the transformation to

        Returns:
            None
        """
        nr_tot = 0
        for i in range(100): # repeat elimination until no more can be eliminated
            subs, cnstr, fixes = self._get_subs(instance)
            nr = len(cnstr)
            if nr == 0:
                break
            nr_tot += nr

            vis = EXPR.ExpressionReplacementVisitor(
                substitute=subs,
                descend_into_named_expressions=True,
                remove_named_expressions=False,
            )

            for c in cnstr: # deactivate constraints that aren't needed
                c.deactivate()
            for v in fixes: # fix variables that can be fixed
                v[0].fix(v[1])

            # Do replacements in Expressions, Constraints, and Objectives
            # where one var is replaced with a linear expression containitng
            # another
            for c in instance.component_data_objects(
                (pyo.Constraint, pyo.Objective),
                descend_into=True,
                active=True
            ):
                c.set_value(expr=vis.dfs_postorder_stack(c.expr))
        _log.info("Eliminated {} variables and constraints".format(nr_tot))

    @staticmethod
    def get_logger(level=None):
        """Return the logger for this, so its a little easier to see the
        get the debugging output, by changing the level.
        """
        if level is not None:
            _log.setLevel(level)
        return _log
