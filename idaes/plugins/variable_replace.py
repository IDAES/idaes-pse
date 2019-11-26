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
from pyomo.common.config import ConfigBlock, ConfigValue, add_docstring_list
from pyomo.core.beta.dict_objects import VarDict


@TransformationFactory.register(
    'replace_variables',
    doc="Replace variables with other variables.")
class ReplaceVariables(NonIsomorphicTransformation):
    """Replace variables with other variables in a Block.

    Keyword arguments below are specified for the ``apply_to`` and
    ``create_using`` functions.

    """
    CONFIG = ConfigBlock()
    CONFIG.declare("substitute", ConfigValue(
        default=[],
        description="List-like of tuples where the first item in a tuple is a "
                    "Pyomo variable to be replaced and the second item in the "
                    "tuple is a Pyomo variable to replace it with. This "
                    "transformation is not reversible."
    ))

    __doc__ = add_docstring_list(__doc__, CONFIG)

    @staticmethod
    def replace(instance, substitute):
        d = {}
        for r in substitute:
            d[id(r[0])] = r[1]
        vis = EXPR.ExpressionReplacementVisitor(substitute=d)
        for c in instance.component_data_objects(pyo.Constraint, descend_into=True):
            c.set_value(expr=vis.dfs_postorder_stack(c.expr))
        for e in instance.component_data_objects(pyo.Expression, descend_into=True):
            vis.dfs_postorder_stack(e)

    def _apply_to(self, instance, **kwds):
        config = self.CONFIG(kwds)
        self.replace(instance, config.substitute)
