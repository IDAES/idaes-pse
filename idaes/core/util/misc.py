
##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes".
##############################################################################
"""
This module contains miscalaneous utility functions for use in IDAES models.
"""

__author__ = "Andrew Lee"


def add_object_reference(self, local_name, remote_object):
    """
    Method to create a reference in the local model to a remote Pyomo object.
    This method should only be used where Pyomo Reference objects are not
    suitable (such as for referencing scalar Pyomo objects where the None
    index is undesirable).

    Args:
        local_name : name to use for local reference (str)
        remote_object : object to make a reference to

    Returns:
        None
    """
    try:
        object.__setattr__(self, local_name, remote_object)
    except AttributeError:
        raise AttributeError("{} failed to construct reference to {} - remote "
                             "object does not exist.".format(self.name,
                                                             remote_object))


__author__ = "Jaffer Ghouse"


def extract_data(data_dict):
    """
    General method that returns a rule to extract data from a python
    dictionary. This method allows the param block to have a database for
    a parameter but extract a subset of this data to initialize a Pyomo
    param object.
    """
    def _rule_initialize(m, *args):
        if len(args) > 1:
            return data_dict[args]
        else:
            return data_dict[args[0]]
    return _rule_initialize
