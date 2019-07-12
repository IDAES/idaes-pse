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
from collections import OrderedDict
from pyomo.environ import value
from pyomo.network import Arc, Port

from idaes.core.util.exceptions import ConfigurationError

__author__ = "John Eslick, Andrew Lee"


def create_stream_table_dataframe(streams,
                                  true_state=False,
                                  time_point=0,
                                  orient='columns'):
    """
    Method to create a stream table in the form of a pandas dataframe. Method
    takes a dict with name keys and stream values. Use an OrderedDict to list
    the streams in a specific order, otherwise the dataframe can be sorted
    later.

    Args:
        streams : dict with name keys and stream values. Names will be used as
            display names for stream table, and streams may be Arcs, Ports or
            StateBlocks.
        true_state : indicated whether the stream table should contain the
            display variables define in the StateBlock (False, default) or the
            state variables (True).
        time_point : point in the time domain at which to generate stream table
            (default = 0)
        orient : orientation of stream table. Accepted values are 'columns'
            (default) where streams are displayed as columns, or 'index' where
            stream are displayed as rows.

    Returns:
        A pandas DataFrame containing the stream table data.
    """
    stream_attributes = OrderedDict()

    for n in streams.keys():
        try:
            if isinstance(streams[n], Arc):
                # Use destination of Arc, as inlets are more likely (?) to be
                # fully-defined StateBlocks
                sb = _get_state_from_port(streams[n].destination, time_point)
            elif isinstance(streams[n], Port):
                sb = _get_state_from_port(streams[n], time_point)
            else:
                sb = streams[n][time_point]

            if true_state:
                disp_dict = sb.define_state_vars()
            else:
                disp_dict = sb.define_display_vars()

            stream_attributes[n] = {}

            for k in disp_dict:
                for i in disp_dict[k]:
                    if i is None:
                        stream_attributes[n][k] = value(disp_dict[k][i])
                    else:
                        stream_attributes[n][k+" "+str(i)] = \
                            value(disp_dict[k][i])
        except (AttributeError, KeyError):
            raise TypeError(
                    f"Unrecognised component provided in stream argument "
                    f"{streams[n]}. get_stream_table_attributes only "
                    f"supports Arcs, Ports or StateBlocks.")

    return DataFrame.from_dict(stream_attributes, orient=orient)


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
            lambda x: "{:#.5g}".format(x))

    # Print stream table
    return stream_table.to_string(na_rep=na_rep,
                                  justify=justify,
                                  float_format=float_format,
                                  **kwargs)


def _get_state_from_port(port, time_point):
    # Check port for _state_block attribute
    try:
        if len(port._state_block) == 1:
            return port._state_block[0][time_point]
        else:
            return port._state_block[0][time_point, port._state_block[1]]
    except AttributeError:
        # Port was not created by IDAES add_port methods. Return exception for
        # the user to fix.
        raise ConfigurationError(
                f"Port {port.name} does not have a _state_block attribute, "
                f"thus cannot determine StateBlock to use for collecting data."
                f" Please provide the associated StateBlock instead, or use "
                f"the IDAES add_port methods to create the Port.")


def generate_table(blocks, attributes, heading=None):
    """
    Create a Pandas DataFrame that contains a list of user-defined attributes
    from a set of Blocks.

    Args:
        blocks (dict): A dictionary with name keys and BlockData objects for
            values. Any name can be associated with a block. Use an OrderedDict
            to show the blocks in a specific order, otherwise the dataframe can
            be sorted later.
        attributes (list or tuple of strings): Attributes to report from a
            Block, can be a Var, Param, or Expression. If an attribute doesn't
            exist or doesn't have a valid value, it will be treated as missing
            data.
        heading (list or tuple of srings): A list of strings that will be used
            as column headings. If None the attribute names will be used.
    Returns:
        (DataFrame): A Pandas dataframe containing a data table
    """
    if heading is None:
        heading = attributes
    st = DataFrame(columns=heading)
    row = [None]*len(attributes)  # not a big deal but save time on realloc
    for key, s in blocks.items():
        for i, a in enumerate(attributes):
            try:
                v = getattr(s, a, None)
                v = value(v, exception=False)
            except ZeroDivisionError:
                v = None
            row[i] = v
        st.loc[key] = row
    return st
