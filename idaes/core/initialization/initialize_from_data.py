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
Initializer class for implementing initialization from a data source
"""
from idaes.core.initialization.initializer_base import InitializerBase

__author__ = "Andrew Lee"


class FromDataInitializer(InitializerBase):
    """
    This is a general purpose Initializer object which attempts to initialize a
    model from user provided data.

    Data can be provided in either json format or as a dict-like structure. The loaded
    solution is then checked to ensure that it satisfies all constraints in the model.

    """

    def initialization_routine(self, model):
        # No action required as data has been loaded in the load_initial_guesses method
        pass
