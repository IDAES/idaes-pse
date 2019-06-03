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

from pandas import DataFrame
from pyomo.environ import value
from pyomo.network import Arc, Port

from idaes.core.util.exceptions import ConfigurationError

__author__ = "John Eslick, Andrew Lee"


def stream_table_dataframe_to_string(stream_table, **kwargs):
    """
    Method to print a stream table from a dataframe. Method takes any argument
    understood by DataFrame.to_string
    """
    # Set some default values for keyword arguments
    na_rep = kwargs.pop("na_rep", "-")
    justify = kwargs.pop("justify", "center")
    float_format = kwargs.pop(
            "float_format",
            lambda x: "{:#.2f}".format(x) if x >= 1 else "{:#.2g}".format(x))

    # Print stream table
    return stream_table.to_string(na_rep=na_rep,
                                  justify=justify,
                                  float_format=float_format,
                                  **kwargs)


def create_stream_table_dataframe(streams,
                                  true_state=False,
                                  time_point=0,
                                  orient='columns'):

    stream_attributes = {}

    for n in streams.keys():
        try:
            if isinstance(streams[n], Arc):
                # Use destination of Arc, as inlets are more likely (?) to be
                # fully-defined StateBlocks
                sb = _get_state_from_port(streams[n].destination)
            elif isinstance(streams[n], Port):
                sb = _get_state_from_port(streams[n])
            else:
                sb = streams[n]

            if true_state:
                disp_dict = sb[time_point].define_state_vars()
            else:
                disp_dict = sb[time_point].define_display_vars()

            stream_attributes[n] = {}

            for k in disp_dict:
                for i in disp_dict[k]:
                    if i is None:
                        stream_attributes[n][k] = value(disp_dict[k][i])
                    else:
                        stream_attributes[n][k+" "+i] = value(disp_dict[k][i])
        except AttributeError:
            raise TypeError(
                    f"Unrecognised component provided in stream argument "
                    f"{streams[n]}. get_stream_table_attributes only "
                    f"supports Arcs, Ports or StateBlocks.")

    stream_table = DataFrame.from_dict(stream_attributes, orient=orient)

    return stream_table


def _get_state_from_port(port):
    # Check port for _state_block attribute
    try:
        return port._state_block
    except AttributeError:
        # Port was not created by IDAES add_port methods. Return exception for
        # the user to fix.
        raise ConfigurationError(
                f"Port {port.name} does not have a _state_block attribute, "
                f"thus cannot determine StateBlock to use for collecting data."
                f" Please provide the associated StateBLock instead, or use "
                f"the IDAES add_port methods to create the Port.")


def stream_table(streams, attributes, heading=None):
    """
    Create a Pandas DataFrame that shows the material state in streams.

    Args:
        streams (dict): A dictionary with stream name keys and StateBlockData
            objects for values.  The stream names do not need to correspond to
            Arcs in the flowhseet. Any name can be associated with a state
            block. Use an OrderedDict to show the streams in a specific order,
            otherwise the dataframe can be sorted later.
        attributes (list or tuple of strings): Attributes to report from a
            StateBlock, can be a Var, Param, or Expression. If an attribute
            doesn't exist or doesn't have a valid value, it will be treated as
            missing data.
        heading (list or tuple of srings): A list of strings that will be used
            as column headings. If None the attribute names will be used.
    Returns:
        (DataFrame): A Pandas dataframe containing a stream table
    """
    if heading is None:
        heading = attributes
    st = DataFrame(columns=heading)
    row = [None]*len(attributes)  # not a big deal but save time on realloc
    for key, s in streams.items():
        for i, a in enumerate(attributes):
            try:
                v = getattr(s, a, None)
                v = value(v, exception=False)
            except ZeroDivisionError:
                v = None
            row[i] = v
        st.loc[key] = row
    return st


def state_table(m, attributes, heading=None):
    """
    Create a Pandas dataframe that shows the material state in every state
    block.

    Args:
        m (Block): Pyomo model or block from which to create a state block
            table
        attributes (list or tuple of strings): Attributes to report from a
            StateBlock, can be a Var, Param, or Expression. If an attribute
            doesn't exist or doesn't have a valid value, it will be treated as
            missing data.
        heading (list or tuple of srings): A list of strings that will be used
            as column headings. If None the attribute names will be used.
    Returns:
        (DataFrame): A Pandas DataFrame with a StateBlock table
    """
    streams = {}  # make a dict for a stream table containing all state blocks
    for c in m.component_objects():
        if isinstance(c, StateBlock):
            for i in c:
                streams[c[i].name] = c[i]
    return stream_table(streams, attributes=attributes, heading=heading)
