#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
from pyomo.core.base.indexed_component import UnindexedComponent_set
from pyomo.core.base.global_set import UnindexedComponent_index
from pyomo.network.port import Port, PortData

__author__ = "Douglas Allan"
# Classes for inlet and outlet ports to assist other functions
# in classifying port and flowsheet topology

# This diamond class structure was copied over from Pyomo
# ____________________________________________________________________________________
#
# Pyomo: Python Optimization Modeling Objects
# Copyright (c) 2008-2026 National Technology and Engineering Solutions of Sandia, LLC
# Under the terms of Contract DE-NA0003525 with National Technology and Engineering
# Solutions of Sandia, LLC, the U.S. Government retains certain rights in this
# software.  This software is distributed under the 3-clause BSD License.
# ____________________________________________________________________________________


class InletPortData(PortData):
    """
    Element of an InletPort
    """

    pass


class InletPort(Port):
    """
    Class to identify whether a port on a model is an inlet port.
    """

    def __new__(cls, *args, **kwds):
        if cls != InletPort:
            return super(Port, cls).__new__(cls)
        if not args or (args[0] is UnindexedComponent_set and len(args) == 1):
            return ScalarInletPort.__new__(ScalarInletPort)
        else:
            return IndexedInletPort.__new__(IndexedInletPort)


class ScalarInletPort(InletPort, InletPortData):
    """
    Unindexed InletPort
    """

    def __init__(self, *args, **kwd):
        InletPortData.__init__(self, component=self)
        InletPort.__init__(self, *args, **kwd)
        self._index = UnindexedComponent_index


class IndexedInletPort(InletPort):
    """
    Indexed InletPort
    """

    pass


class OutletPortData(PortData):
    """
    Element of an OutletPort
    """

    pass


class OutletPort(Port):
    """
    Class to identify whether a port on a model is an outlet port.
    """

    def __new__(cls, *args, **kwds):
        if cls != OutletPort:
            return super(Port, cls).__new__(cls)
        if not args or (args[0] is UnindexedComponent_set and len(args) == 1):
            return ScalarOutletPort.__new__(ScalarOutletPort)
        else:
            return IndexedOutletPort.__new__(IndexedOutletPort)


class ScalarOutletPort(OutletPort, OutletPortData):
    """
    Unindexed OutletPort
    """

    def __init__(self, *args, **kwd):
        OutletPortData.__init__(self, component=self)
        OutletPort.__init__(self, *args, **kwd)
        self._index = UnindexedComponent_index


class IndexedOutletPort(OutletPort):
    """
    Indexed OutletPort
    """

    pass
