# -*- coding: UTF-8 -*-
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
"""
This module contains utility functions for reporting structural statistics of
IDAES models.
"""

__author__ = "Andrew Lee"

import sys

from pyomo.environ import Block, Constraint, Expression, Objective, Var, value
from pyomo.dae import DerivativeVar
from pyomo.core.expr.current import identify_variables
from pyomo.core.kernel.component_set import ComponentSet

from idaes.core.util.exceptions import ConfigurationError


class ModelStatistics():
    """
    Creates a model statistics object for quantifying and reporting model
    statistics.
    """
    def __init__(self, block, always_recalculate=False):
        self.model_object = block
        self.always_recalculate = always_recalculate

    def _activated_block_component_generator(self, ctype):
        for b in self.model_object.component_data_objects(
                    ctype=Block, active=True, descend_into=True):
            for c in b.component_data_objects(ctype=ctype,
                                              active=None,
                                              descend_into=False):
                yield c

    # -------------------------------------------------------------------------
    # Block methods
    def total_block_set(self, recalculate=None):
        if recalculate is None:
            recalculate = self.always_recalculate

        if hasattr(self._total_block_set):
            if not recalculate:
                return self._total_block_set
            else:
                self.del_component(self._total_block_set)

        self._total_block_set = ComponentSet(
                self.model_object.component_data_objects(
                        ctype=Block, active=None, descend_into=True))
        self._total_block_set.add(self)

        return self._total_block_set

    def number_total_blocks(self, recalculate=None):
        if recalculate is None:
            recalculate = self.always_recalculate

        if hasattr(self, "_total_block_set"):
            return len(self.total_block_set(recalculate=recalculate))
        else:
            b = 1  # Start at 1 to include self
            for o in self.model_object.component_data_objects(
                        ctype=Block, active=None, descend_into=True):
                b += 1
            return b

    def deactivated_block_set(self, recalculate=None):
        if recalculate is None:
            recalculate = self.always_recalculate

        if hasattr(self._deactivated_block_set):
            if not recalculate:
                return self._deactivated_block_set
            else:
                self.del_component(self._deactivated_block_set)

        self._deactivated_block_set = ComponentSet(
                self.model_object.component_data_objects(
                        ctype=Block, active=False, descend_into=True))

        return self._deactivated_block_set

    def number_deactivated_blocks(self, recalculate=None):
        if recalculate is None:
            recalculate = self.always_recalculate

        if hasattr(self, "_deactivated_block_set"):
            return len(self.deactivated_block_set(recalculate=recalculate))
        else:
            b = 0
            for o in self.model_object.component_data_objects(
                        ctype=Block, active=False, descend_into=True):
                b += 1
            return b

    # -------------------------------------------------------------------------
    # Basic Constraint methods
    def total_constraint_set(self, recalculate=None):
        if recalculate is None:
            recalculate = self.always_recalculate

        if hasattr(self, "_total_constraint_set"):
            if not recalculate:
                return self._total_constraint_set
            else:
                self.del_component(self._total_constraint_set)

        self._total_constraint_set = ComponentSet()

        self._total_constraint_set.update(
                self._activated_block_component_generator(ctype=Constraint))

        return self._total_constraint_set

    def number_total_constraints(self, recalculate=None):
        if recalculate is None:
            recalculate = self.always_recalculate

        if hasattr(self, "_total_constraint_set"):
            return len(self.total_constraint_set(recalculate=recalculate))
        else:
            tc = 0
            for c in self._activated_block_component_generator(
                    ctype=Constraint):
                tc += 1
            return tc

    # -------------------------------------------------------------------------
    # Equality Constraints
    def total_equality_generator(self):
        for c in self._activated_block_component_generator(ctype=Constraint):
            if (c.upper is not None and
                    c.lower is not None and
                    c.upper == c.lower):
                yield c

    def total_equality_set(self, recalculate=None):
        if recalculate is None:
            recalculate = self.always_recalculate

        if hasattr(self, "_total_equality_set"):
            if not recalculate:
                return self._total_equality_set
            else:
                self.del_component(self._total_equality_set)

        self._total_equality_set = ComponentSet()
        self._total_equality_set.update(self.total_equality_generator())

        return self._total_equality_set

    def number_total_equalities(self, recalculate=None):
        if recalculate is None:
            recalculate = self.always_recalculate

        if hasattr(self, "_total_equality_set"):
            return len(self.total_equality_set(recalculate=recalculate))
        else:
            tc = 0
            for c in self.total_equality_generator():
                tc += 1
            return tc

    def activated_equalities_generator(self):
        for c in self.model_object.component_data_objects(
                    Constraint, active=True, descend_into=True):
            if (c.upper is not None and c.lower is not None and
                    c.upper == c.lower):
                yield c

    def activated_equalities_set(self, recalculate=None):
        if recalculate is None:
            recalculate = self.always_recalculate

        if hasattr(self._activated_equalities_set):
            if not recalculate:
                return self._activated_equalities_set
            else:
                self.del_component(self._activated_equalities_set)

        self._activated_equalities_set = ComponentSet()

        if hasattr(self, "_equalities_set"):
            # Make set of activated equalities from this
            self._activated_equalities_set.update(
                    ec for ec in self._equalities_set if ec.active)
        elif hasattr(self, "_total_constraint_set"):
            # Make set of activated equalities from this
            for c in self._total_constraint_set:
                if (c.upper is not None and c.lower is not None and
                        c.upper == c.lower and c.active):
                    self._activated_equalities_set.add(c)
        else:
            # Build set from generator
            for c in self.activated_equalities_generator():
                self._activated_equalities_set.add(c)

        return self._activated_equalities_set

    def number_activated_equalities(self, recalculate=None):
        if recalculate is None:
            recalculate = self.always_recalculate

        if hasattr(self, "_activated_equalities_set"):
            return len(self.activated_equalities_set(recalculate=recalculate))
        else:
            c = 0
            for o in self.activated_equalities_generator():
                c += 1
            return c

    def deactivated_equality_generator(self):
        for c in self.total_equality_generator():
            if not c.active:
                yield c

    def deactivated_equality_set(self, recalculate=None):
        if recalculate is None:
            recalculate = self.always_recalculate

        if hasattr(self, "_deactivated_equality_set"):
            if not recalculate:
                return self._deactivated_equality_set
            else:
                self.del_component(self._deactivated_equality_set)

        self._deactivated_equality_set = ComponentSet()
        self._deactivated_equality_set.update(
                self.deactivated_equality_generator())

        return self._deactivated_equality_set

    def number_deactivated_equalities(self, recalculate=None):
        if recalculate is None:
            recalculate = self.always_recalculate

        if hasattr(self, "_deactivated_equality_set"):
            return len(self.deactivated_equality_set(recalculate=recalculate))
        else:
            tc = 0
            for c in self.deactivated_equality_generator():
                tc += 1
            return tc

    # -------------------------------------------------------------------------
    # Inequality Constraints
    def total_inequality_generator(self):
        for c in self._activated_block_component_generator(ctype=Constraint):
            if c.upper is None or c.lower is None:
                yield c

    def total_inequality_set(self, recalculate=None):
        if recalculate is None:
            recalculate = self.always_recalculate

        if hasattr(self, "_total_inequality_set"):
            if not recalculate:
                return self._total_inequality_set
            else:
                self.del_component(self._total_inequality_set)

        self._total_inequality_set = ComponentSet()
        self._total_inequality_set.update(self.total_inequality_generator())

        return self._total_inequality_set

    def number_total_inequalities(self, recalculate=None):
        if recalculate is None:
            recalculate = self.always_recalculate

        if hasattr(self, "_total_inequality_set"):
            return len(self.total_inequality_set(recalculate=recalculate))
        else:
            tc = 0
            for c in self.total_inequality_generator():
                tc += 1
            return tc

    def activated_inequalities_generator(self):
        for c in self.model_object.component_data_objects(
                    Constraint, active=True, descend_into=True):
            if c.upper is None or c.lower is None:
                yield c

    def activated_inequalities_set(self, recalculate=None):
        if recalculate is None:
            recalculate = self.always_recalculate

        if hasattr(self._activated_inequalities_set):
            if not recalculate:
                return self._activated_inequalities_set
            else:
                self.del_component(self._activated_inequalities_set)

        self._activated_inequalities_set = ComponentSet()

        if hasattr(self, "_inequalities_set"):
            # Make set of activated inequalities from this
            self._activated_inequalities_set.update(
                    ec for ec in self._inequalities_set if ec.active)
        elif hasattr(self, "_total_constraint_set"):
            # Make set of activated inequalities from this
            for c in self._constraint_set:
                if (c.upper is None or c.lower) is None and c.active:
                    self._activated_inequalities_set.add(c)
        else:
            # Build set from generator
            for c in self.activated_inequalities_generator():
                self._activated_inequalities_set.add(c)

        return self._activated_inequalities_set

    def number_activated_inequalities(self, recalculate=None):
        if recalculate is None:
            recalculate = self.always_recalculate

        if hasattr(self, "_activated_inequality_set"):
            return len(self.activated_inequalities_set(
                    recalculate=recalculate))
        else:
            c = 0
            for o in self.activated_inequalities_generator():
                c += 1
            return c

    def deactivated_inequality_generator(self):
        for c in self.total_inequality_generator():
            if not c.active:
                yield c

    def deactivated_inequality_set(self, recalculate=None):
        if recalculate is None:
            recalculate = self.always_recalculate

        if hasattr(self, "_deactivated_inequality_set"):
            if not recalculate:
                return self._deactivated_inequality_set
            else:
                self.del_component(self._deactivated_inequality_set)

        self._deactivated_inequality_set = ComponentSet()
        self._deactivated_inequality_set.update(
                self.deactivated_inequality_generator())

        return self._deactivated_inequality_set

    def number_deactivated_inequalities(self, recalculate=None):
        if recalculate is None:
            recalculate = self.always_recalculate

        if hasattr(self, "_deactivated_inequality_set"):
            return len(self.deactivated_inequality_set(
                    recalculate=recalculate))
        else:
            tc = 0
            for c in self.deactivated_inequality_generator():
                tc += 1
            return tc

    # -------------------------------------------------------------------------
    # Basic Variable Methods
    # Always use ComponentSets for Vars to avoid duplication of References
    # i.e. number methods should alwys use the ComponentSet, not a generator
    def variables_set(self, recalculate=None):
        if recalculate is None:
            recalculate = self.always_recalculate

        if hasattr(self, "_variables_set"):
            if not recalculate:
                return self._variables_set
            else:
                self.del_component(self._variables_set)

        self._variables_set = ComponentSet()

        self._variables_set.update(
                self.model_object.component_data_objects(
                        ctype=Var,
                        active=True,
                        descend_into=True))

        return self._variables_set

    def number_variables(self, recalculate=None):
        if recalculate is None:
            recalculate = self.always_recalculate

        return len(self.variables_set(recalculate=recalculate))

    def fixed_variables_set(self, recalculate=None):
        if recalculate is None:
            recalculate = self.always_recalculate

        if hasattr(self, "_fixed_variables_set"):
            if not recalculate:
                return self._fixed_variables_set
            else:
                self.del_component(self._fixed_variables_set)

        self._fixed_variables_set = ComponentSet()

        if hasattr(self, "_variables_set"):
            for v in self.variables_set(recalculate=recalculate):
                if v.fixed:
                    self._fixed_variables_set.add(v)
        else:
            for v in self.model_object.component_data_objects(
                        ctype=Var,
                        active=True,
                        descend_into=True):
                if v.fixed:
                    self._fixed_variables_set.add(v)

        return self._fixed_variables_set

    def number_fixed_variables(self, recalculate=None):
        if recalculate is None:
            recalculate = self.always_recalculate

        return len(self.fixed_variables_set(recalculate=recalculate))

    # -------------------------------------------------------------------------
    # Variables in Constraints
    def variables_in_activated_constraints_set(self, recalculate=None):
        if recalculate is None:
            recalculate = self.always_recalculate

        if hasattr(self, "_variables_in_activated_constraints_set"):
            if not recalculate:
                return self._variables_in_activated_constraints_set
            else:
                self.del_component(
                        self._variables_in_activated_constraints_set)

        self._variables_in_activated_constraints_set = ComponentSet()

        if hasattr(self, "_total_constraint_set"):
            cs = self.total_constraint_set(recalculate=recalculate)
        else:
            cs = self.model_object.component_data_objects(ctype=Constraint,
                                                          active=True,
                                                          descend_into=True)

        for c in cs:
            for v in identify_variables(c.body):
                self._variables_in_activated_constraints_set.add(v)

        return self._variables_in_activated_constraints_set

    def number_variables_in_activated_constraints(self, recalculate=None):
        if recalculate is None:
            recalculate = self.always_recalculate

        return len(self.variables_in_activated_constraints_set(
                recalculate=recalculate))

    def variables_in_activated_equalities_set(self, recalculate=None):
        if recalculate is None:
            recalculate = self.always_recalculate

        if hasattr(self, "_variables_in_activated_equalities_set"):
            if not recalculate:
                return self._variables_in_activated_equalities_set
            else:
                self.del_component(self._variables_in_activated_equalities_set)

        self._variables_in_activated_equalities_set = ComponentSet()

        if hasattr(self, "_activated_equalities_set"):
            ec = self.activated_equalities_set(recalculate=recalculate)
        else:
            ec = self.activated_equalities_generator()

        for c in ec:
            for v in identify_variables(c.body):
                self._variables_in_activated_equalities_set.add(v)

        return self._variables_in_activated_equalities_set

    def number_variables_in_activated_equalities(self, recalculate=None):
        if recalculate is None:
            recalculate = self.always_recalculate

        return len(self.variables_in_activated_equalities_set(
                        recalculate=recalculate))

    def variables_in_activated_inequalities_set(self, recalculate=None):
        if recalculate is None:
            recalculate = self.always_recalculate

        if hasattr(self, "_variables_in_activated_inequalities_set"):
            if not recalculate:
                return self._variables_in_activated_inequalities_set
            else:
                self.del_component(
                        self._variables_in_activated_inequalities_set)

        self._variables_in_activated_inequalities_set = ComponentSet()

        if hasattr(self, "_activated_inequalities_set"):
            ec = self.activated_inequalities_set(recalculate=recalculate)
        else:
            ec = self.activated_inequalities_generator()

        for c in ec:
            for v in identify_variables(c.body):
                self._variables_in_activated_inequalities_set.add(v)

        return self._variables_in_activated_inequalities_set

    def number_variables_in_activated_inequalities(self, recalculate=None):
        if recalculate is None:
            recalculate = self.always_recalculate

        return len(self.variables_in_activated_inequalities_set(
                        recalculate=recalculate))

    def variables_only_in_inequalities(self, recalculate=None):
        if recalculate is None:
            recalculate = self.always_recalculate

        if hasattr(self, "_variables_only_in_inequalities_set"):
            if not recalculate:
                return self._variables_only_in_inequalities_set
            else:
                self.del_component(self._variables_only_in_inequalities_set)

        self._variables_only_in_inequalities_set = (
                self.variables_in_activated_inequalities_set(
                        recalculate=recalculate) -
                self.variables_in_activated_equalities_set(
                        recalculate=recalculate))

        return self._variables_only_in_inequalities_set

    def number_variables_only_in_inequalities(self, recalculate=None):
        if recalculate is None:
            recalculate = self.always_recalculate

        return len(self.variables_only_in_inequalities(
                    recalculate=recalculate))

    def fixed_variables_only_in_inequalities(self, recalculate=None):
        if recalculate is None:
            recalculate = self.always_recalculate

        if hasattr(self, "_fixed_variables_only_in_inequalities_set"):
            if not recalculate:
                return self._fixed_variables_only_in_inequalities_set
            else:
                self.del_component(
                        self._fixed_variables_only_in_inequalities_set)

        self._fixed_variables_only_in_inequalities_set = ComponentSet()

        for v in self.variables_only_in_inequalities(recalculate=recalculate):
            if v.fixed:
                self._fixed_variables_only_in_inequalities_set.add(v)

        return self._fixed_variables_only_in_inequalities_set

    def number_fixed_variables_only_in_inequalities(self, recalculate=None):
        if recalculate is None:
            recalculate = self.always_recalculate

        return len(self.fixed_variables_only_in_inequalities(
                    recalculate=recalculate))

    # -------------------------------------------------------------------------
    # Fixed Variables in Constraints
    def fixed_variables_in_activated_equalities_set(self, recalculate=None):
        if recalculate is None:
            recalculate = self.always_recalculate

        if hasattr(self, "_fixed_variables_in_activated_equalities_set"):
            if not recalculate:
                return self._fixed_variables_in_activated_equalities_set
            else:
                self.del_component(
                        self._fixed_variables_in_activated_equalities_set)

        self._fixed_variables_in_activated_equalities_set = ComponentSet()

        for v in self.variables_in_activated_equalities_set(
                recalculate=recalculate):
            if v.fixed:
                self._fixed_variables_in_activated_equalities_set.add(v)

        return self._fixed_variables_in_activated_equalities_set

    def number_fixed_variables_in_activated_equalities(self, recalculate=None):
        if recalculate is None:
            recalculate = self.always_recalculate

        return len(self.fixed_variables_in_activated_equalities_set(
                recalculate=recalculate))

    # -------------------------------------------------------------------------
    # Unused and un-Transformed Variables
    def unused_variables_set(self, recalculate=None):
        if recalculate is None:
            recalculate = self.always_recalculate

        if hasattr(self, "_unused_variables_set"):
            if not recalculate:
                return self._unused_variables_set
            else:
                self.del_component(self._unused_variables_set)

        self._unused_variables_set = (
            self.variables_set(recalculate=recalculate) -
            self.variables_in_activated_constraints_set(
                    recalculate=recalculate))

        return self._unused_variables_set

    def number_unused_variables(self, recalculate=None):
        if recalculate is None:
            recalculate = self.always_recalculate

        return len(self.unused_variables_set(recalculate=recalculate))

    def fixed_unused_variables_set(self, recalculate=None):
        if recalculate is None:
            recalculate = self.always_recalculate

        if hasattr(self, "_fixed_unused_variables_set"):
            if not recalculate:
                return self._fixed_unused_variables_set
            else:
                self.del_component(self._fixed_unused_variables_set)

        self._fixed_unused_variables_set = ComponentSet()

        for v in self.unused_variables_set(recalculate=recalculate):
            if v.fixed:
                self._fixed_unused_variables_set.add(v)

        return self._fixed_unused_variables_set

    def number_fixed_unused_variables(self, recalculate=None):
        if recalculate is None:
            recalculate = self.always_recalculate

        return len(self.fixed_unused_variables_set(recalculate=recalculate))

    def derivative_variables_set(self, recalculate=None):
        if recalculate is None:
            recalculate = self.always_recalculate

        if hasattr(self, "_derivative_variables_set"):
            if not recalculate:
                return self._derivative_variables_set
            else:
                self.del_component(self._derivative_variables_set)

        self._derivative_variables_set = ComponentSet()

        self._derivative_variables_set.update(
                self.model_object.component_data_objects(
                        ctype=Var,
                        active=True,
                        descend_into=True))

        return self._derivative_variables_set

    def number_derivative_variables(self, recalculate=None):
        if recalculate is None:
            recalculate = self.always_recalculate

        return len(self.derivative_variables_set(recalculate=recalculate))

    # -------------------------------------------------------------------------
    # Objective methods
    def total_objective_generator(self):
        for o in self._activated_block_component_generator(ctype=Objective):
            yield o

    def total_objective_set(self, recalculate=None):
        if recalculate is None:
            recalculate = self.always_recalculate

        if hasattr(self, "_total_objective_set"):
            if not recalculate:
                return self._total_objective_set
            else:
                self.del_component(self._total_objective_set)

        self._total_objective_set = ComponentSet()
        self._total_objective_set.update(
                self.total_objective_generator())

        return self._total_objective_set

    def number_total_objectives(self, recalculate=None):
        if recalculate is None:
            recalculate = self.always_recalculate

        if hasattr(self, "_total_objective_set"):
            return len(self.total_objective_set(
                    recalculate=recalculate))
        else:
            to = 0
            for o in self.total_objective_generator():
                to += 1
            return to

    def deactivated_objective_generator(self):
        for o in self._activated_block_component_generator(ctype=Objective):
            if not o.active:
                yield o

    def deactivated_objective_set(self, recalculate=None):
        if recalculate is None:
            recalculate = self.always_recalculate

        if hasattr(self, "_deactivated_objective_set"):
            if not recalculate:
                return self._deactivated_objective_set
            else:
                self.del_component(self._deactivated_objective_set)

        self._deactivated_objective_set = ComponentSet()
        self._deactivated_objective_set.update(
                self.deactivated_objective_generator())

        return self._deactivated_objective_set

    def number_deactivated_objectives(self, recalculate=None):
        if recalculate is None:
            recalculate = self.always_recalculate

        if hasattr(self, "_deactivated_objective_set"):
            return len(self.deactivated_objective_set(
                    recalculate=recalculate))
        else:
            to = 0
            for o in self.deactivated_objective_generator():
                to += 1
            return to

    # -------------------------------------------------------------------------
    # Expression methods
    # Always use ComponentsSets here to avoid duplication of References
    def expressions_set(self, recalculate=None):
        if recalculate is None:
            recalculate = self.always_recalculate

        if hasattr(self, "_expressions_set"):
            if not recalculate:
                return self._expressions_set
            else:
                self.del_component(self._expressions_set)

        self._expressions_set = ComponentSet(
                self.model_object.component_data_objects(
                        ctype=Expression,
                        active=True,
                        descend_into=True))

        return self._expressions_set

    def number_expressions(self, recalculate=None):
        if recalculate is None:
            recalculate = self.always_recalculate

        return len(self.expressions_set(recalculate=recalculate))

    # -------------------------------------------------------------------------
    # Other model statistics
    def degrees_of_freedom(self, recalculate=None):
        if recalculate is None:
            recalculate = self.always_recalculate

        return (self.number_variables_in_activated_equalities(
                        recalculate=recalculate) -
                self.number_fixed_variables_in_activated_equalities(
                        recalculate=recalculate) -
                self.number_activated_equalities(
                        recalculate=recalculate))

    def large_residual_set(self, tol=1e-5, recalculate=None):
        if recalculate is None:
            recalculate = self.always_recalculate

        if hasattr(self, "_large_residual_set"):
            if not recalculate:
                return self._large_residual_set
            else:
                self.del_component(self._large_residual_set)

        self._large_residual_set = ComponentSet()

        for c in self.model_object.component_data_objects(ctype=Constraint,
                                                          active=True,
                                                          descend_into=True):
            if c.active and value(c.lower - c.body()) > tol:
                self._large_residual_set.add(c)
            elif c.active and value(c.body() - c.upper) > tol:
                self._large_residual_set.add(c)

        return self._large_residual_set

    def number_large_residuals(self, tol=1e-5, recalculate=None):
        if recalculate is None:
            recalculate = self.always_recalculate

        if hasattr(self, "_large_residual_set"):
            return len(self.large_residual_set(
                    recalculate=recalculate))
        else:
            lr = 0
            for c in self.model_object.component_data_objects(
                    ctype=Constraint, active=True, descend_into=True):
                if c.active and value(c.lower - c.body()) > tol:
                    lr += 1
                elif c.active and value(c.body() - c.upper) > tol:
                    lr += 1
            return lr

    # -------------------------------------------------------------------------
    # Reporting methods
    def report(self,
               ostream=None,
               recalculate=None):
        """
        Method to print a report of the model statistics for a Pyomo Block

        Args:
            block - the Block object to report statistics from
            recalculate - whether model statistics should be recalculated or
                    not. Options are True, False and None (default). Setting
                    this to False will save time, but will not capture any
                    changes made to the model state since the statistics were
                    generated. If None, will use the global
                    'always_recalcualte` setting from the ModelStatistics
                    object.

        Returns:
            Printed output of the model statistics
        """
        if ostream is None:
            ostream = sys.stdout

        if recalculate is None:
            recalculate = self.always_recalculate

        tab = " "*4
        header = '='*72

        if self.model_object.name == "unknown":
            name_str = ""
        else:
            name_str = f"-  {self.model_object.name}"

        ostream.write("\n")
        ostream.write(header+"\n")
        ostream.write(f"Model Statistics  {name_str} \n")
        ostream.write("\n")
        ostream.write(f"Degrees of Freedom: "
                      f"{self.degrees_of_freedom(recalculate=recalculate)} \n")
        ostream.write("\n")
        ostream.write(f"Total No. Variables: "
                      f"{self.number_variables(recalculate=recalculate)} \n")
        ostream.write(f"{tab}No. Fixed Variables: "
                      f"{self.number_fixed_variables(recalculate=recalculate)}"
                      f"\n")
        ostream.write(
            f"{tab}No. Unused Variables: "
            f"{self.number_unused_variables(recalculate=recalculate)} (Fixed):"
            f"{self.number_fixed_unused_variables(recalculate=recalculate)})"
            f"\n")
        nv_alias = self.number_variables_only_in_inequalities
        nfv_alias = self.number_fixed_variables_only_in_inequalities
        ostream.write(
            f"{tab}No. Variables only in Inequalities:"
            f" {nv_alias(recalculate=recalculate)}"
            f" (Fixed: {nfv_alias(recalculate=recalculate)}) \n")
        ostream.write("\n")
        ostream.write(
                f"Total No. Constraints: "
                f"{self.number_total_constraints(recalculate=recalculate)} \n")
        ostream.write(
            f"{tab}No. Equality Constraints: "
            f"{self.number_total_equalities(recalculate=recalculate)}"
            f" (Deactivated: "
            f"{self.number_deactivated_equalities(recalculate=recalculate)})"
            f"\n")
        ostream.write(
            f"{tab}No. Inequality Constraints: "
            f"{self.number_total_inequalities(recalculate=recalculate)}"
            f" (Deactivated: "
            f"{self.number_deactivated_inequalities(recalculate=recalculate)})"
            f"\n")
        ostream.write("\n")
        ostream.write(
            f"No. Objectives: "
            f"{self.number_total_objectives(recalculate=recalculate)}"
            f" (Deactivated: "
            f"{self.number_deactivated_objectives(recalculate=recalculate)})"
            f"\n")
        ostream.write("\n")
        ostream.write(
            f"No. Blocks: {self.number_total_blocks(recalculate=recalculate)}"
            f" (Deactivated: "
            f"{self.number_deactivated_blocks(recalculate=recalculate)}) \n")
        ostream.write(f"No. Expressions: "
                      f"{self.number_expressions(recalculate=recalculate)} \n")
        ostream.write(header+"\n")
        ostream.write("\n")
