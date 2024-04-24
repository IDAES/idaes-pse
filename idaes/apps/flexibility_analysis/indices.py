#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
This module contains some utility classes for creating meaningful indices
when formulating the flexibility test problem.
"""
class _ConIndex(object):
    def __init__(self, con, bound):
        self._con = con
        self._bound = bound

    @property
    def con(self):
        # pylint: disable=missing-function-docstring
        return self._con

    @property
    def bound(self):
        # pylint: disable=missing-function-docstring
        return self._bound

    def __repr__(self):
        if self.bound is None:
            return str(self.con)
        else:
            return str((str(self.con), str(self.bound)))

    def __str__(self):
        return repr(self)

    def __eq__(self, other):
        if isinstance(other, _ConIndex):
            return self.con is other.con and self.bound is other.bound
        return False

    def __hash__(self):
        return hash((self.con, self.bound))


class _VarIndex(object):
    def __init__(self, var, bound):
        self._var = var
        self._bound = bound

    @property
    def var(self):
        # pylint: disable=missing-function-docstring
        return self._var

    @property
    def bound(self):
        # pylint: disable=missing-function-docstring
        return self._bound

    def __repr__(self):
        if self.bound is None:
            return str(self.var)
        else:
            return str((str(self.var), str(self.bound)))

    def __str__(self):
        return repr(self)

    def __eq__(self, other):
        if isinstance(other, _VarIndex):
            return self.var is other.var and self.bound is other.bound
        return False

    def __hash__(self):
        return hash((id(self.var), self.bound))
