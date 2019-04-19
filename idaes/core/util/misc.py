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
        raise AttributeError("{} failed to construct reference to {} - remote "
                             "object does not exist.".format(self.name,
                                                             remote_object))

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
def svg_tag(tags, svg, outfile=None, time=None):
    """
    Replace text in a svg with tag values for the model.

    Args:
        tags: A dictionary where the key is the tag and the value is a Pyomo
            Refernce indexed by time.
        svg: a file pointer or a string continaing svg contents
        outfile: a file name to save the results, if None return svg string
        time: if None use the first time, otherwise a time in the time domain to
            report results at.

    Returns:
        String for svg
    """

    if isinstance(svg, str):
        svg_str = svg
    elif hasattr(svg, "read"):
        svg_str = svg.read()
    else:
        raise TypeError("SVG must either be a string or a file-like object")

    doc = xml.dom.minidom.parseString(svg_str)
    name = doc.getElementsByTagName('text')
    for t in name:
        if(t.attributes['id'].value in tags):
            name = doc.getElementsByTagName('text')

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
