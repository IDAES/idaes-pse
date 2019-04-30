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

import pandas as pd
import pyomo.environ as pyo
from idaes.core.property_base import StateBlock, StateBlockData

__author__ = "John Eslick"

def stream_table(streams, attributes, heading=None):
    """
    Create a Pandas DataFrame that shows the material state in streams.

    Args:
        streams (dict): A dictionary with stream name keys and StateBlockData
            objects for values.  The stream names do not need to correspond to
            Arcs in the flowhseet. Any name can be associated with a state
            block. Use an OrderedDict to show the streams in a specific order,
            otherwise the dataframe can be sorted later.
        attributes (list or tupple of strings): Attributes to report from a
            StateBlock, can be a Var, Param, or Expression. If an attribute
            doesn't exist or doesn't have a valid value, it will be treated as
            missing data.
        heading (list or tupple of srings): A list of strings that will be used
            as column headings. If None the attribute names will be used.
    Returns:
        (DataFrame): A Pandas dataframe containting a stream table
    """
    if heading is None: heading = attributes
    st = pd.DataFrame(columns=heading)
    row = [None]*len(attributes) # not a big deal but save time on realloc
    for key, s in streams.items():
        for i, a in enumerate(attributes):
            try:
                v = getattr(s, a, None)
                v = pyo.value(v, exception=False)
            except ZeroDivisionError:
                v = None
            row[i] = v
        st.loc[key] = row
    return st

def state_table(m, attributes, heading=None):
    """
    Create a Pandas dataframe that shows the material state in every state block.

    Args:
        m (Block): Pyomo model or block from which to create a state block table
        attributes (list or tupple of strings): Attributes to report from a
            StateBlock, can be a Var, Param, or Expression. If an attribute
            doesn't exist or doesn't have a valid value, it will be treated as
            missing data.
        heading (list or tupple of srings): A list of strings that will be used
            as column headings. If None the attribute names will be used.
    Returns:
        (DataFrame): A Pandas DataFrame with a StateBlock table
    """
    streams = {} #make a dict for a stream table containing all the state blocks
    for c in m.component_objects():
        if isinstance(c, StateBlock):
            for i in c:
                streams[c[i].name] = c[i]
    return stream_table(streams, attributes=attributes, heading=heading)
