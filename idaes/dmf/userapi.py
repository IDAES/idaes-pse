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
"""
Data Management Framework high-level functions.
"""
# stdlib
import logging
import sys
# package
from idaes.dmf import DMF
from idaes.dmf import errors
from idaes.dmf import propindex

__author__ = 'Dan Gunter <dkgunter@lbl.gov>'

_log = logging.getLogger(__name__)


def get_workspace(path='', name=None, desc=None, create=False, errs=None,
                  **kwargs):
    """Create or load a DMF workspace.

    If the :class:`DMF` constructor, throws an exception, this catches
    it and prints the error to the provided stream (or stdout).

    See :class:`DMF` for details on arguments.

    Args:
        path (str): Path to workspace.
        name (str): Name to be used for workspace.
        desc (str): Longer description of workspace.
        create (bool): If the path to the workspace does not exist,
                       this controls whether to create it.
        errs (object): Stream for errors, stdout is used if None

    Returns:
        DMF: New instance, or None if it failed.
    """
    dmf = None
    try:
        dmf = DMF(path=path, name=name, desc=desc, create=create, **kwargs)
    except errors.DMFError as err:
        if errs is None:
            errs = sys.stdout
        msg = 'Error creating DMF workspace\n'
        if isinstance(err, errors.DMFError) and not create:
            msg += 'Directory not found, and "create" flag is False\n'
            msg += 'If you want to create the workspace, try again with ' \
                   'create=True\n'
        else:
            msg += '{}\n'.format(err)
        msg += '\npath: {}\nname: {}\ndesc: {}\n'.format(path, name, desc)
        errs.write(msg)
    return dmf


def find_property_packages(dmf, properties=None):
    """Find all property packages matching provided criteria.

    Return the matching packages as a generator.

    Args:
        dmf (DMF): Data Management Framework instance
        properties (List[str]): Names of properties that must be present in the
              returned packages.
    Returns:
        Generator[idaes.dmf.resource.Resource]: Each object has the property
            data (properties and default units) in its `.data` attribute.
    """
    query = {'tags': [propindex.DMFVisitor.INDEXED_PROPERTY_TAG]}
    for prop in properties:
        field = 'data.properties.' + prop
        query[field] = True
    return dmf.find(query)


def index_property_packages(dmf):
    propindex.index_property_metadata(dmf)
