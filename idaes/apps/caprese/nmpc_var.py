# -*- coding: utf-8 -*-
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
""" Variable classes for time-indexed variables that can be treated specially
by NMPC.
"""

from pyomo.core.base.var import IndexedVar
from pyomo.core.base.indexed_component_slice import IndexedComponent_slice


class NmpcVar(IndexedVar):
    """
    This class is meant to be used as a ctype for references to time-only
    slices. It provides attributes that we would like to store for each
    "variable," but not for different time-indices of the same "variable."
    Consistent with the NMPC/DAE literature, a "variable" here is
    something indexed only by time.

    Users may want to access these "time-only variables," and may do so
    by looping over `component_objects`. For example:

    >>> for var in model.component_objects(NmpcVar):
    >>>     var[:].set_value(var[0])

    """

    def __init__(self, *args, **kwargs):
        if not args:
            raise NotImplementedError(
                "%s component must be indexed by at least one set." % self.__class__
            )
        self.setpoint = kwargs.pop("setpoint", None)
        self.weight = kwargs.pop("weight", None)
        self.variance = kwargs.pop("variance", None)
        self.nominal = kwargs.pop("nominal", None)
        self.noise_bounds = kwargs.pop("noise_bounds", None)
        kwargs.setdefault("ctype", type(self))
        super(NmpcVar, self).__init__(*args, **kwargs)


"""
The following classes serve as custom ctypes for references
to time-only slices and calls to `component_objects` when
the user or developer wants to specify a particular subset
of the "NMPC variables."
"""


class DiffVar(NmpcVar):
    _attr = "differential"


class DerivVar(NmpcVar):
    _attr = "derivative"


class AlgVar(NmpcVar):
    _attr = "algebraic"


class InputVar(NmpcVar):
    _attr = "input"


class FixedVar(NmpcVar):
    _attr = "fixed"


class MeasuredVar(NmpcVar):
    _attr = "measurement"


class _NmpcVector(IndexedVar):
    """
    This is a class used as a ctype for references to `NmpcVar` objects.
    It is used to construct a "vector" of these variables with a
    "Pyomothonic" API. These "vectors" are indexed by an integer describing
    the "coordinate" of the variable within the vector and by time.
    In addition to accessing the underlying data objects with a 2-dimensional
    index, this class adds methods to iterate over the underlying `NmpcVar`
    objects. This is necessary to access the special attributes on these
    variables that would not make sense to attach to a data object.

    For example:

    >>> input_vars = Reference(INPUT_BLOCK[:].var[:], ctype=_NmpcVector)
    >>> input_vars[:, 0].fix()
    >>> input_vars.set_setpoint(input_setpoint_list)

    """

    def _generate_referenced_vars(self):
        _slice = self.referent.duplicate()
        # This class should only be instantiated as a reference...
        assert type(_slice) is IndexedComponent_slice
        # _slice looks something like BLOCK[:].var[:]
        assert _slice._call_stack[-1][0] == IndexedComponent_slice.get_item
        assert _slice._call_stack[-1][1] == slice(None)
        assert _slice._call_stack[-2][0] == IndexedComponent_slice.get_attribute
        assert _slice._call_stack[-3][0] == IndexedComponent_slice.slice_info
        _slice._call_stack.pop()
        _slice._len -= 1
        for var in _slice:
            assert isinstance(var, NmpcVar)
            yield var

    def set_setpoint(self, setpoint):
        referent_gen = self._generate_referenced_vars()
        try:
            for var, sp in zip(referent_gen, setpoint):
                var.setpoint = sp
        except TypeError:
            for var in self._generate_referenced_vars():
                var.setpoint = setpoint

    def get_setpoint(self):
        for var in self._generate_referenced_vars():
            yield var.setpoint

    @property
    def values(self):
        referent_gen = self._generate_referenced_vars()
        return list(list(var[t].value for t in var) for var in referent_gen)

    @values.setter
    def values(self, vals):
        referent_gen = self._generate_referenced_vars()
        try:
            for var, val in zip(referent_gen, vals):
                # var is time-indexed
                var[:].set_value(val)
        except TypeError:
            # A scalar value was provided. Set for all i, t
            self[...].set_value(vals)
