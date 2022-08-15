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
Exception classes.
"""
import logging

__author__ = "Dan Gunter <dkgunter@lbl.gov>"

_log = logging.getLogger(__name__)


class DMFError(Exception):
    def __init__(self, detailed_error="No details"):
        msg = "{}".format(detailed_error)
        super(DMFError, self).__init__(msg)


class ParseError(Exception):
    pass


class CommandError(Exception):
    def __init__(self, command, operation, details):
        msg = 'Operation "{op}" in command "{c}" failed: {d}'.format(
            op=operation, c=command, d=details
        )
        super(CommandError, self).__init__(msg)


class WorkspaceError(DMFError):
    pass


class WorkspaceNotFoundError(WorkspaceError):
    def __init__(self, from_dir):
        msg = 'Workspace not found for path "{}" '.format(from_dir)
        super(WorkspaceNotFoundError, self).__init__(msg)


class WorkspaceCannotCreateError(WorkspaceError):
    def __init__(self, path):
        msg = 'Unable to create new workspace at "{}"'.format(path)
        super(WorkspaceCannotCreateError, self).__init__(msg)


class WorkspaceConfNotFoundError(WorkspaceError):
    def __init__(self, path):
        msg = 'Workspace config not found at path "{}" '.format(path)
        super(WorkspaceConfNotFoundError, self).__init__(msg)


class WorkspaceConfMissingField(WorkspaceError):
    def __init__(self, path, name, desc):
        msg = 'Workspace config at path "{}" missing {} field "{}"'.format(
            path, desc, name
        )
        super(WorkspaceConfMissingField, self).__init__(msg)


class FileError(Exception):
    pass


class ResourceError(Exception):
    pass


class NoSuchResourceError(ResourceError):
    def __init__(self, name=None, id_=None):
        if name is not None and id_ is None:
            msg = 'No resource of type "{}" found'.format(name)
        elif id_ is not None and name is None:
            msg = 'No resource "{}" found'.format(id_)
        elif name is None and id_ is None:
            msg = "Resource not found"
        else:
            msg = 'No resource "{}" of type "{}" found'.format(id_, name)
        super(NoSuchResourceError, self).__init__(msg)


class DuplicateResourceError(ResourceError):
    def __init__(self, op, id_):
        msg = 'While executing "{}": Duplicate resource "{}"'.format(op, id_)
        super(DuplicateResourceError, self).__init__(msg)


class BadResourceError(ResourceError):
    pass


class SearchError(Exception):
    def __init__(self, spec, problem):
        msg = 'Search "{}" failed: {}'.format(spec, problem)
        super(SearchError, self).__init__(msg)


class ModuleFormatError(Exception):
    def __init__(self, module_name, type_, what):
        msg = (
            'Python module "{}" does not conform to conventions for a '
            "{} module: {}".format(module_name, type_, what)
        )
        super(ModuleFormatError, self).__init__(msg)


class DmfError(Exception):
    pass


class InvalidRelationError(DmfError):
    def __init__(self, subj, pred, obj):
        msg = "Invalid relation: {} --({})--> {}".format(subj, pred, obj)
        super(InvalidRelationError, self).__init__(msg)


class DataFormatError(DmfError):
    def __init__(self, dtype, err):
        msg = 'Bad data format for type "{}":\n{}'.format(dtype, err)
        super(DataFormatError, self).__init__(msg)


# Alamo


class AlamoError(DmfError):
    def __init__(self, msg):
        super(AlamoError, self).__init__("ALAMO Error: {}".format(msg))


class AlamoDisabledError(AlamoError):
    def __init__(self):
        super(AlamoDisabledError, self).__init__("ALAMO is disabled")
