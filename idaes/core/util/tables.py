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

from pandas import DataFrame
from collections import OrderedDict
from pyomo.environ import value
from pyomo.network import Arc, Port
from pyomo.core.base.units_container import units

import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)

__author__ = "John Eslick, Andrew Lee"


def arcs_to_stream_dict(
    blk, additional=None, descend_into=True, sort=False, prepend=None, s=None
):
    """
    Creates a stream dictionary from the Arcs in a model, using the Arc names as
    keys. This can be used to automate the creation of the streams dictionary
    needed for the ``create_stream_table_dataframe()`` and  ``stream_states_dict()``
    functions.

    Args:
        blk (pyomo.environ._BlockData): Pyomo model to search for Arcs
        additional (dict): Additional states to add to the stream dictionary,
            which aren't represented by arcs in blk, for example feed or
            product streams without Arcs attached or states internal to a unit
            model.
        descend_into (bool): If True, search subblocks for Arcs as well. The
            default is True.
        sort (bool): If True sort keys and return an OrderedDict
        prepend (str): Prepend a string to the arc name joined with a '.'.
            This can be useful to prevent conflicting names when sub blocks
            contain Arcs that have the same names when used in combination
            with descend_into=False.
        s (dict): Add streams to an existing stream dict.

    Returns:
        Dictionary with Arc names as keys and the Arcs as values.

    """
    if s is None:
        s = {}
    for c in blk.component_objects(Arc, descend_into=descend_into):
        key = c.getname()
        if prepend is not None:
            key = ".".join([prepend, key])
        s[key] = c
    if additional is not None:
        s.update(additional)
    if sort:
        s = OrderedDict(sorted(s.items()))
    return s


def stream_states_dict(streams, time_point=0):
    """
    Method to create a dictionary of state block representing stream states.
    This takes a dict with stream name keys and stream values.

    Args:
        streams : dict with name keys and stream values. Names will be used as
            display names for stream table, and streams may be Arcs, Ports or
            StateBlocks.
        time_point : point in the time domain at which to generate stream table
            (default = 0)

    Returns:
        A pandas DataFrame containing the stream table data.
    """
    stream_dict = OrderedDict()

    def _stream_dict_add(sb, n, i=None):
        """add a line to the stream table"""
        if i is None:
            key = n
        else:
            key = "{}[{}]".format(n, i)
        stream_dict[key] = sb

    for n in streams.keys():
        if isinstance(streams[n], Arc):
            for i, a in streams[n].items():
                try:
                    # if getting the StateBlock from the destination port
                    # fails for any reason try the source port. This could
                    # happen if a port does not have an associated
                    # StateBlock. For example a surrogate model may not
                    # use state blocks, unit models may handle physical
                    # properties without state blocks, or the port could
                    # be used to serve the purpose of a translator block.
                    sb = _get_state_from_port(a.ports[1], time_point)
                except:
                    sb = _get_state_from_port(a.ports[0], time_point)
                _stream_dict_add(sb, n, i)
        elif isinstance(streams[n], Port):
            sb = _get_state_from_port(streams[n], time_point)
            _stream_dict_add(sb, n)
        else:
            # _IndexedStateBlock is a private class, so cannot directly test
            # whether  streams[n] is one or not.
            try:
                sb = streams[n][time_point]
            except KeyError as err:
                raise TypeError(
                    f"Either component type of stream argument {streams[n]} "
                    f"is unindexed or {time_point} is not a member of its "
                    f"indexing set."
                ) from err
            _stream_dict_add(sb, n)
    return stream_dict


def tag_state_quantities(blocks, attributes, labels, exception=False):
    """Take a stream states dictionary, and return a tag dictionary for stream
    quantities.  This takes a dictionary (blk) that has state block labels as
    keys and state blocks as values.  The attributes are a list of attributes to
    tag.  If an element of the attribute list is list-like, the fist element is
    the attribute and the remaining elements are indexes.  Lables provides a list
    of attribute lables to be used to create the tag.  Tags are blk_key + label
    for the attribute.

    Args:
        blocks (dict): Dictionary of state blocks.  The key is the block label to
            be used in the tag, and the value is a state block.
        attributes (list-like): A list of attriutes to tag.  It is okay if a
            particular attribute does not exist in a state bock.  This allows
            you to mix state blocks with differnt sets of attributes. If an
            attribute is indexed, the attribute can be specified as a list or
            tuple where the first element is the attribute and the remaining
            elements are indexes.
        labels (list-like): These are attribute lables.  The order corresponds to the
            attribute list.  They are used to create the tags.  Tags are in the
            form blk.key + label.
        exception (bool): If True, raise exceptions releated to invalid or
            missing indexes. If false missing or bad indexes are ignored and
            None is used for the table value.  Setting this to False allows
            tables where some state blocks have the same attributes with differnt
            indexing. (default is True)

    Return:
        (dict): Dictionary where the keys are tags and the values are model
            attributes, usually Pyomo component data objects.
    """

    tags = {}
    if labels is None:
        lables = attributes
        for a in attributes:
            if isinstance(a, (tuple, list)):
                if len(a) == 2:
                    # in case there are multiple indexes and user gives tuple
                    label = f"{a[0]}[{a[1]}]"
                if len(a) > 2:
                    label = f"{a[0]}[{a[1:]}]"
                else:
                    label = a[0]

    for key, s in blocks.items():
        for i, a in enumerate(attributes):
            j = None
            if isinstance(a, (list, tuple)):
                # if a is list or tuple, the first element should be the
                # attribute and the remaining elements should be indexes.
                if len(a) == 2:
                    j = a[1]  # catch user supplying list-like of indexes
                if len(a) > 2:
                    j = a[1:]
                # if len(a) == 1, we'll say that's fine here.  Don't know why you
                # would put the attribute in a list-like if not indexed, but I'll
                # allow it.
                a = a[0]
            v = getattr(s, a, None)
            if j is not None and v is not None:
                try:
                    v = v[j]
                except KeyError:
                    if not exception:
                        v = None
                    else:
                        _log.error(f"{j} is not a valid index of {a}")
                        raise KeyError(f"{j} is not a valid index of {a}")
            try:
                value(v, exception=False)
            except TypeError:
                if not exception:
                    v = None
                else:
                    _log.error(f"Cannot calculate value of {a} (may be subscriptable)")
                    raise TypeError(
                        f"Cannot calculate value of {a} (may be subscriptable)"
                    )
            except ZeroDivisionError:
                pass  # this one is okay
            if v is not None:
                tags[f"{key}{labels[i]}"] = v
    return tags


def create_stream_table_dataframe(
    streams, true_state=False, time_point=0, orient="columns", add_units=False
):
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
        add_units : Add a Units column to the dataframe representing the units
            of the stream values.

    Returns:
        A pandas DataFrame containing the stream table data.
    """
    stream_attributes = OrderedDict()
    stream_states = stream_states_dict(streams=streams, time_point=time_point)
    full_keys = []  # List of all rows in dataframe to fill in missing data

    if add_units and stream_states:
        stream_attributes["Units"] = {}

    for key, sb in stream_states.items():
        stream_attributes[key] = {}
        if true_state:
            disp_dict = sb.define_state_vars()
        else:
            disp_dict = sb.define_display_vars()
        for k in disp_dict:
            for i in disp_dict[k]:
                stream_key = k if i is None else f"{k} {i}"
                stream_attributes[key][stream_key] = value(disp_dict[k][i])
                if add_units:
                    pyomo_unit = units.get_units(disp_dict[k][i])
                    if pyomo_unit is not None:
                        pint_unit = pyomo_unit._get_pint_unit()
                        stream_attributes["Units"][stream_key] = {
                            "raw": str(pyomo_unit),
                            "html": "{:~H}".format(pint_unit),
                            "latex": "{:~L}".format(pint_unit),
                        }
                    else:
                        stream_attributes["Units"][stream_key] = None
                if stream_key not in full_keys:
                    full_keys.append(stream_key)

    # Check for missing rows in any stream, and fill with "-" if needed
    for k, v in stream_attributes.items():
        for r in full_keys:
            if r not in v.keys():
                # Missing row, fill with placeholder
                v[r] = "-"

    return DataFrame.from_dict(stream_attributes, orient=orient)


def stream_table_dataframe_to_string(stream_table, **kwargs):
    """
    Method to print a stream table from a dataframe. Method takes any argument
    understood by DataFrame.to_string
    """
    # Set some default values for keyword arguments
    na_rep = kwargs.pop("na_rep", "-")
    justify = kwargs.pop("justify", "center")
    float_format = kwargs.pop("float_format", lambda x: "{:#.5g}".format(x))

    # Print stream table
    return stream_table.to_string(
        na_rep=na_rep, justify=justify, float_format=float_format, **kwargs
    )


def _get_state_from_port(port, time_point):
    """
    Attempt to find a StateBlock-like object connected to a Port. If the
    object is indexed both in space and time, assume that the time index
    comes first.  If no components are assigned to the Port, raise a
    ValueError. If the first component's parent block has no index, raise an
    AttributeError. If different variables on the port appear to be connected
    to different state blocks, raise a RuntimeError.

    Args:
        port (pyomo.network.Port): a port with variables derived from some
            single StateBlock
        time_point : point in the time domain at which to index StateBlock
            (default = 0)

    Returns:
        (StateBlock-like) : an object containing all the components contained
            in the port.
    """
    vlist = list(port.iter_vars())
    states = [v.parent_block().parent_component() for v in vlist]

    if len(vlist) == 0:
        raise ValueError(
            f"No block could be retrieved from Port {port.name} "
            f"because it contains no components."
        )
    # Check the number of indices of the parent property block. If its indexed
    # both in space and time, keep the second, spatial index and throw out the
    # first, temporal index. If that ordering is changed, this method will
    # need to be changed as well.
    try:
        idx = vlist[0].parent_block().index()
    except AttributeError as err:
        raise AttributeError(
            f"No block could be retrieved from Port {port.name} "
            f"because block {vlist[0].parent_block().name} has no index."
        ) from err
    # Assuming the time index is always first and the spatial indices are all
    # the same
    if isinstance(idx, tuple):
        idx = (time_point, vlist[0].parent_block().index()[1:])

    else:
        idx = (time_point,)
    # This method also assumes that ports with different spatial indices won't
    # end up at the same port. Otherwise this check is insufficient.
    if all(states[0] is s for s in states):
        return states[0][idx]
    raise RuntimeError(
        f"No block could be retrieved from Port {port.name} "
        f"because components are derived from multiple blocks."
    )


def generate_table(blocks, attributes, heading=None, exception=True):
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
        exception (bool): If True, raise exceptions releated to invalid or
            missing indexes. If false missing or bad indexes are ignored and
            None is used for the table value.  Setting this to False allows
            tables where some state blocks have the same attributes with differnt
            indexing. (default is True)
    Returns:
        (DataFrame): A Pandas dataframe containing a data table
    """
    if heading is None:
        heading = attributes
    st = DataFrame(columns=heading)
    row = [None] * len(attributes)  # not a big deal but save time on realloc
    for key, s in blocks.items():
        for i, a in enumerate(attributes):
            j = None
            if isinstance(a, (list, tuple)):
                # if a is list or tuple, assume index supplied
                try:
                    assert len(a) > 1
                except AssertionError:
                    _log.error(f"An index must be supplided for attribute {a[0]}")
                    raise AssertionError(
                        f"An index must be supplided for attribute {a[0]}"
                    )
                j = a[1:]
                a = a[0]
            v = getattr(s, a, None)
            if j is not None and v is not None:
                try:
                    v = v[j]
                except KeyError:
                    if not exception:
                        v = None
                    else:
                        _log.error(f"{j} is not a valid index of {a}")
                        raise KeyError(f"{j} is not a valid index of {a}")
            try:
                v = value(v, exception=False)
            except TypeError:
                if not exception:
                    v = None
                else:
                    _log.error(f"Cannot calculate value of {a} (may be subscriptable)")
                    raise TypeError(
                        f"Cannot calculate value of {a} (may be subscriptable)"
                    )
            except ZeroDivisionError:
                v = None
            row[i] = v
        st.loc[key] = row
    return st
