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
Resource representaitons.
"""
# stdlib
from collections import namedtuple
from datetime import datetime
import getpass
import logging
import os
import pprint
import re
import uuid
# third-party
import jsonschema
import pendulum
import six
# local
from .util import datetime_timestamp

__author__ = 'Dan Gunter <dkgunter@lbl.gov>'

_log = logging.getLogger(__name__)

#: Constants for relation predicates
PR_DERIVED = 'derived'  # derivedFrom
PR_CONTAINS = 'contains'
PR_USES = 'uses'
PR_VERSION = 'version'
RELATION_PREDICATES = {PR_DERIVED, PR_CONTAINS, PR_USES, PR_VERSION}

#: Constants for resource 'types'
TY_EXPERIMENT = 'experiment'
TY_TABULAR = 'tabular_data'
TY_PROPERTY = 'propertydb'
TY_FLOWSHEET = 'flowsheet'
TY_NOTEBOOK = 'notebook'
TY_CODE = 'code'
TY_SURRMOD = 'surrogate_model'
TY_DATA = 'data'
TY_OTHER = 'other'
RESOURCE_TYPES = {TY_EXPERIMENT, TY_TABULAR, TY_PROPERTY, TY_FLOWSHEET,
                  TY_NOTEBOOK, TY_CODE, TY_SURRMOD, TY_DATA, TY_OTHER}

# Constants for fields in stored relations
RR_PRED = 'predicate'
RR_SUBJ = 'subject'
RR_OBJ = 'object'
RR_ID = 'identifier'
RR_ROLE = 'role'

RESOURCE_SCHEMA = {
    "$schema": "http://json-schema.org/draft-04/schema#",
    "id": "http://idaes.org",
    "definitions": {
        "SemanticVersion": {
            "title": "Version",
            "description": "Resource version using semantic versioning",
            "type": "array",
            "items": [
                {"type": "integer"},
                {"type": "integer"},
                {"type": "integer"},
                {"type": "string"}
            ],
            "minItems": 4
        }
    },
    "type": "object",
    "properties": {
        "aliases": {
            "type": "array",
            "items": {
                "type": "string"
            }
        },
        "codes": {
            "type": "array",
            "items": {
                "type": "object",
                "properties": {
                    "type": {
                        "type": "string",
                        "enum": ["method", "function", "module", "class",
                                 "file",
                                 "package", "repository", "notebook"]
                    },
                    "desc": {"type": "string"},
                    "name": {"type": "string"},
                    "language": {"type": "string"},
                    "idhash": {"type": "string"},
                    "location": {"type": "string"},
                    "version": {"$ref": "#/definitions/SemanticVersion"}
                },
                "required": ["name"]
            }
        },
        "collaborators": {
            "type": "array",
            "items": {
                "type": "object",
                "properties": {
                    "email": {
                        "type": "string",
                        "format": "email"
                    },
                    "name": {"type": "string"}
                },
                "required": ["name"]
            }
        },
        "created": {
            "type": "number"
        },
        "creator": {
            "type": "object",
            "properties": {
                "email": {
                    "type": "string",
                    "format": "email"
                },
                "name": {"type": "string"}
            },
            "required": ["name"]
        },
        "data": {
            "type": "object"
        },
        "datafiles": {
            "type": "array",
            "items": {
                "type": "object",
                "properties": {
                    "desc": {"type": "string"},
                    "metadata": {
                        "type": "object"
                    },
                    "mimetype": {"type": "string"},
                    "path": {"type": "string"},
                    "is_copy": {"type": "boolean"}
                },
                "required": ["path"]
            }
        },
        "datafiles_dir": {
            "type": "string"
        },
        "desc": {
            "type": "string"
        },
        "id_": {
            "type": "string"
        },
        "modified": {
            "type": "number"
        },
        "relations": {
            "type": "array",
            "items": {
                "type": "object",
                "properties": {
                    RR_PRED: {
                        "type": "string",
                        "enum": list(RELATION_PREDICATES)
                    },
                    RR_ID: {"type": "string"},
                    RR_ROLE: {
                        "type": "string",
                        "enum": [RR_SUBJ, RR_OBJ]
                    }
                },
                "required": [RR_PRED, RR_ID, RR_ROLE]
            }
        },
        "sources": {
            "type": "array",
            "items": {
                "type": "object",
                "properties": {
                    "date": {"type": "number"},
                    "doi": {"type": "string"},
                    "isbn": {"type": "string"},
                    "language": {"type": "string"},
                    "source": {"type": "string"}
                }
            }
        },
        "tags": {
            "type": "array",
            "items": {
                "type": "string",
            }
        },
        "type": {
            "type": "string",
            "enum": list(RESOURCE_TYPES)
        },
        "version_info": {
            "type": "object",
            "properties": {
                "created": {"type": "number"},
                "name": {"type": "string"},
                "version": {"$ref": "#/definitions/SemanticVersion"}
            }
        }
    },
    "required": ["id_"],
    "additionalProperties": False
}


class Dict(dict):
    """Subclass of dict that has a 'dirty' bit.
    """
    def __init__(self, *args, **kwargs):
        super(Dict, self).__init__(*args, **kwargs)
        self._dirty = True

    def __setitem__(self, key, value):
        self._dirty = True
        super(Dict, self).__setitem__(key, value)

    def set_clean(self):
        self._dirty = False

    def is_dirty(self):
        return self._dirty


class Resource(object):
    """Core object for the Data Management Framework.
    """
    ID_FIELD = 'id_'     #: Identifier field name constant
    TYPE_FIELD = 'type'  #: Resource type field name constant

    def __init__(self, value=None, type_=None):
        self._set_defaults()
        if value:
            self.v.update(value)
        if type_ is not None:
            self.v[self.TYPE_FIELD] = type_
        self._validator = jsonschema.Draft4Validator(RESOURCE_SCHEMA)
        self._validations = 0  # count validations; mostly for testing
        self.do_copy = self.is_tmp = False  # flags for copying datafiles

    def _set_defaults(self):
        now = date_float(pendulum.utcnow())
        self.v = Dict({
            self.ID_FIELD: identifier_str(),
            self.TYPE_FIELD: TY_OTHER,
            'aliases': [],
            'codes': [],
            'collaborators': [],
            'created': now,
            'modified': now,
            'creator': {'name': getpass.getuser()},
            'data': {},
            'datafiles': [],
            'datafiles_dir': '',
            'desc': '',
            'relations': [],
            'sources': [],
            'tags': [],
            'version_info': {
                'created': now,
                'version': (0, 0, 0),
                'name': ''
            }
        })

    def _massage_values(self):
        try:
            # convert dates
            for item in self.v['sources']:
                if not isinstance(item['date'], float):
                    item['date'] = date_float(item['date'])
            if not isinstance(self.v['created'], float):
                self.v['created'] = date_float(self.v['created'])
            if not isinstance(self.v['modified'], float):
                self.v['modified'] = date_float(self.v['modified'])
            if not isinstance(self.v['version_info']['created'], float):
                self.v['version_info']['created'] = date_float(
                    self.v['version_info']['created'])
            # convert versions
            if not isinstance(self.v['version_info']['version'], list):
                self.v['version_info']['version'] = version_list(
                    self.v['version_info']['version'])
            for i, code in enumerate(self.v['codes']):
                if not isinstance(code['version'], list):
                    code['version'] = version_list(code['version'])
                    self.v['codes'][i] = code
        except (TypeError, ValueError, KeyError) as err:
            raise ValueError('While converting resource values: {}'
                             .format(err))
        self.v.set_clean()

    def validate(self):
        if self.v.is_dirty():
            self._massage_values()
            self._validator.validate(self.v)
            self._validations += 1

    @property
    def id(self):
        """Get resource identifier.
        """
        return self.v[self.ID_FIELD]

    def set_id(self, value=None):
        self.v[self.ID_FIELD] = identifier_str(value)

    @property
    def type(self):
        """Get resource type.
        """
        return self.v[self.TYPE_FIELD]

    @property
    def data(self):
        """Get JSON data for this resource.
        """
        return self.v['data']

    @data.setter
    def data(self, value):
        """Set JSON data for this resource.
        """
        self.v['data'] = value

    def get_datafiles(self, mode='r'):
        """Generate readable file objects for 'datafiles' in resource.

        Args:
            mode (str): Mode for ``open()``
        Returns:
            generator: Generates ``file`` objects.
        """
        dfdir = self.v['datafiles_dir']
        for datafile in self.v['datafiles']:
            if not dfdir:
                path = datafile['path']
            else:
                path = os.path.join(dfdir, datafile['path'])
            fp = open(path, mode=mode)
            yield fp

    def _repr_text_(self):
        return pprint.pformat(self.v, indent=2)

#
# Function(s) to help creating [two-way] relations
# between resources
#


#: Provide attribute access to an RDF subject, predicate, object triple
Triple = namedtuple('Triple', 'subject predicate object')


def create_relation(rel):
    """Create a relationship between two Resource instances.

    Relations are stored in both the `subject` and `object` resources, in
    the following way::

        If R = (subject)S, (predicate)P, and (object)O
        then store the following:
          In S.relations: {predicate: P, identifier:O.id, role:subject}
          In O.relations: {predicate: P, identifier:S.id, role:object}

    Args:
        rel (Triple): Relation triple. The 'subject' and 'object' parts
                    should be :class:`Resource`, and the 'predicate' should
                    be a simple string.
    Returns:
        None
    Raises:
        ValueError: if this relation already exists in the subject or
                    object resource, or the predicate is not in the list
                    of valid ones in RELATION_PREDICATES
    """
    if rel.predicate not in RELATION_PREDICATES:
        raise ValueError('Bad predicate: "{}" not in: {}'.format(
            rel.predicate, ', '.join(list(RELATION_PREDICATES))))
    rel_d = {RR_PRED: rel.predicate,
             RR_ID: rel.object.v[Resource.ID_FIELD],
             RR_ROLE: RR_SUBJ}
    if rel_d in rel.subject.v['relations']:
        raise ValueError('Duplicate relation for subject: {}'.format(rel))
    rel.subject.v['relations'].append(rel_d)
    rel_d = {RR_PRED: rel.predicate,
             RR_ID: rel.subject.v[Resource.ID_FIELD],
             RR_ROLE: RR_OBJ}
    # note: hard for this to happen unless the relation was added manually
    if rel_d in rel.object.v['relations']:
        raise ValueError('Duplicate relation for object: {}'.format(rel))
    rel.object.v['relations'].append(rel_d)


def create_relation_args(*args):
    """Syntactic sugar to take 3 args instead of a Triple.
    """
    return create_relation(Triple(*args))


def triple_from_resource_relations(id_, rrel):
    """Create a Triple from one entry in resource['relations'].

    Args:
        id_ (str): Identifier of the containing resource.
        rrel (dict): Stored relation with three keys, see `create_relation()`.
    Return:
        Triple: A triple
    """
    if rrel[RR_ROLE] == RR_SUBJ:
        rel = Triple(id_, rrel[RR_PRED], rrel[RR_ID])
    else:
        rel = Triple(rrel[RR_ID], rrel[RR_PRED], id_)
    return rel


#
# Some handy-dandy conversion functions.
#


def date_float(value):
    def bad_date(e):
        raise ValueError('Cannot convert date "{}" to float: {}'
                         .format(value, e))
    dt, usec = None, 0
    if isinstance(value, pendulum.Pendulum):
        return value.timestamp()
    elif isinstance(value, datetime):
        dt = value
    elif isinstance(value, tuple):
        try:
            dt = datetime(*value)
        except TypeError as err:
            bad_date(err)
    elif isinstance(value, six.string_types):
        try:
            dt = pendulum.parse(value)
        except pendulum.exceptions.ParserError as err:
            bad_date(err)
    elif isinstance(value, float) or isinstance(value, int):
        try:
            dt = datetime.fromtimestamp(value)
        except ValueError as err:
            bad_date(err)
    if dt is None:
        raise ValueError('Cannot convert date, value is "None"')
    return datetime_timestamp(dt) + usec  # just a float


# def isoformat(ts):
#     return datetime.fromtimestamp(ts).isoformat()


def version_list(value):
    """Semantic version.

    Three numeric identifiers, separated by a dot.
    Trailing non-numeric characters allowed.

    Inputs, string or tuple, may have less than three numeric
    identifiers, but internally the value will be padded with
    zeros to always be of length four.

    A leading dash or underscore in the trailing non-numeric characters
    is removed.
    
    Some examples of valid inputs and how they translate to 4-part versions:

    .. testcode::

        version_list('1')
        version_list('1.1')
        version_list('1a')
        version_list('1.12.1')
        version_list('1.12.13-1')

    .. testoutput::

        [1, 0, 0, '']
        [1, 1, 0, '']
        [1, 0, 0, 'a']
        [1, 12, 1, '']
        [1, 12, 13, '1']

    Some example inputs that will fail:

    .. testcode::

        for bad_input in ('rc3',      # too short
                          '1.a.1.',   # non-number in middle
                          '1.12.13.x' # too long
            ):
            try:
                version_list(bad_input)
            except ValueError:
                print(f"failed: {bad_input}")

    .. testoutput::

        failed: rc3
        failed: 1.a.1.
        failed: 1.12.13.x

    Returns:
        list: [major:int, minor:int, debug:int, release-type:str]
    """
    def bad_version(v):
        raise ValueError("Bad version: {}".format(v))
    ver = ()
    if isinstance(value, list) or isinstance(value, tuple):
        ver = value
    elif isinstance(value, str):
        ver = value.split('.', 2)
    elif isinstance(value, int):
        ver = (value, 0, 0)
    else:
        bad_version(value)
    if len(ver) < 1:
        bad_version(value)
    verlist = []
    # leading version numbers
    for i in range(len(ver) - 1):
        try:
            verlist.append(int(ver[i]))
        except ValueError:
            bad_version(value)
    # last version number
    s = ver[-1]
    extra = ''
    if isinstance(s, int):
        verlist.append(s if len(verlist) < 3 else str(s))
    elif isinstance(s, six.string_types):
        if s:
            m = re.match('([0-9]+)?(.*)', s)
            if m.group(1) is not None:
                verlist.append(int(m.group(1)))
            extra = m.group(2)
    else:  # last version must be int or str
        bad_version(value)
    # must have at least one numbered version
    if len(verlist) == 0:
        bad_version(value)
    # pad with zeros, and add non-numeric ID
    while len(verlist) < 3:
        verlist.append(0)
    if extra and extra[0] == '.':
        # cannot start extra version with '.'
        bad_version(value)
    if extra and extra[0] in ('-', '_'):
        extra = extra[1:]
    verlist.append(extra)
    return verlist


def format_version(values):
    s = '{}.{}.{}'.format(*values[:3])
    if values[3]:
        s += '-{}'.format(values[3])
    return s


def identifier_str(value=None):
    """Unique identifier.

    Args:
        value (str): If given, validate that it is a 32-byte str
                     If not given or None, set new random value.
    """
    # regular expression for identifier: hex string len=32
    id_expr = '[0-9a-f]{32}'
    if value is None:
        value = uuid.uuid4().hex
    elif not re.match(id_expr, value):
        raise ValueError('Bad format for identifier "{}", must match '
                         'regular exprresion "{}"'.format(value, id_expr))
    return value
