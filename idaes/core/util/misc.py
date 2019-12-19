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
This module contains miscellaneous utility functions for use in IDAES models.
"""
import xml.dom.minidom
import pyomo.environ as pyo

# Author: Andrew Lee
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
        raise AttributeError(
            "{} failed to construct reference to {} - remote "
            "object does not exist.".format(self.name, remote_object)
        )


# Author: Jaffer Ghouse
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


# Author: John Eslick
def TagReference(s, description=""):
    """
    Create a Pyomo reference with an added description string attribute to
    describe the reference. The intended use for these references is to create a
    time-indexed reference to variables in a model corresponding to plant
    measurment tags.

    Args:
        s: Pyomo time slice of a variable or expression
        description (str): A description the measurment

    Returns:
        A Pyomo Reference object with an added doc attribute
    """
    r = pyo.Reference(s)
    r.description = description
    return r


# Author John Eslick
def svg_tag(
    tags,
    svg,
    outfile=None,
    idx=None,
    tag_map=None,
    show_tags=False,
    byte_encoding="utf-8",
):
    """
    Replace text in a SVG with tag values for the model. This works by looking
    for text elements in the SVG with IDs that match the tags or are in tag_map.

    Args:
        tags: A dictionary where the key is the tag and the value is a Pyomo
            Refernce.  The refernce could be indexed. In yypical IDAES
            applications the references would be indexed by time.
        svg: a file pointer or a string continaing svg contents
        outfile: a file name to save the results, if None don't save
        idx: if None not indexed, otherwise an index in the indexing set of the
            reference
        tag_map: dictionary with svg id keys and tag values, to map svg ids to
            tags
        show_tags: Put tag labels of the diagram instead of numbers
        byte_encoding: If svg is given as a byte-array, use this encoding to
            convert it to a string.

    Returns:
        String for SVG
    """
    if isinstance(svg, str):  # assume this is svg content string
        pass
    elif isinstance(svg, bytes):
        svg = svg.decode(byte_encoding)
    elif hasattr(svg, "read"):  # file-like object to svg
        svg = svg.read()
    else:
        raise TypeError("SVG must either be a string or a file-like object")
    # Make tag map here because the tags may not make valid XML IDs if no
    # tag_map provided we'll go ahead and handle XML @ (maybe more in future)
    if tag_map is None:
        tag_map = dict()
        for tag in tags:
            new_tag = tag.replace("@", "_")
            tag_map[new_tag] = tag
    # Search for text in the svg that has an id in tags
    doc = xml.dom.minidom.parseString(svg)
    texts = doc.getElementsByTagName("text")
    for t in texts:
        id = t.attributes["id"].value
        if id in tag_map:
            # if it's multiline change last line
            tspan = t.getElementsByTagName("tspan")[-1].childNodes[0]
            try:
                if show_tags:
                    val = tag_map[id]
                elif idx is None:
                    val = pyo.value(tags[tag_map[id]], exception=False)
                else:
                    val = pyo.value(tags[tag_map[id]][idx], exception=False)
            except ZeroDivisionError:
                val = "Divide_by_0"
            try:
                tspan.nodeValue = "{:.4e}".format(val)
            except ValueError:  # whatever it is can't be scientific notation
                tspan.nodeValue = val

    new_svg = doc.toxml()
    # If outfile is provided save to a file
    if outfile is not None:
        with open(outfile, "w") as f:
            f.write(new_svg)
    return new_svg


# Author: John Eslick
def copy_port_values(destination, source):
    """
    Copy the variable values in the source port to the destination port. The
    ports must containt the same variables.

    Args:
        (pyomo.Port): Copy values from this port
        (pyomo.Port): Copy values to this port

    Returns:
        None
    """
    for k, v in destination.vars.items():
        if isinstance(v, pyo.Var):
            for i in v:
                v[i].value = pyo.value(source.vars[k][i])
