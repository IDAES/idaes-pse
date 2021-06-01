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
Resource representaitons.
"""
# stdlib
import abc
import argparse
from collections import namedtuple
from datetime import datetime
import getpass
import hashlib
import json
from json import JSONDecodeError
import logging
import os
from pathlib import Path
import pprint
import re
import sys
from typing import List, Union, Iterator, Optional, IO
import uuid

# third-party
import jsonschema
import pandas
import yaml

# local
from .util import datetime_timestamp, parse_datetime

__author__ = "Dan Gunter"

_log = logging.getLogger(__name__)


class ProgLangExt:
    """Helper class to map from file extensions to
    names of the programming language.
    """

    _extmap = {
        "py": "Python",
        "pyc": "Python/compiled",
        "c": "C",
        "cpp": "C++",
        "cxx": "C++",
        "f": "FORTRAN",
        "f77": "FORTRAN",
        "f90": "FORTRAN",
        "js": "JavaScript",
        "jl": "Julia",
    }

    @classmethod
    def guess(cls, ext, default=None):
        ext = ext.lower()
        return cls._extmap.get(ext, default)


class Predicates:
    """Constants for relation predicates (types of relations).

    This can be extended at runtime by simply adding non-underscore attributes
    to this class, e.g.:

        Predicates.destroys = "destroys"
    """
    derived = "derived"  #: object is derived from the subject
    contains = "contains"  #: object is contained by the subject
    uses = "uses"  #: object is used by the subject
    version = "version"  #: object is a (new) version of the subject

    @classmethod
    def valid(cls, value) -> bool:
        """Whether the 'value' is a valid predicate.
        """
        return value in cls.all()

    @classmethod
    def all(cls) -> List[str]:
        """Return all predicates, as a list of strings.
        """
        return [k for k in cls.__dict__ if k and (k != "all") and (k[0] != "_")]


class ResourceTypes:
    experiment = "experiment"               #: Experiment(s)
    tabular = "tabular_data"                #: Tabular data
    publication = "publication"             #: Published work(s)
    property = "propertydb"                 #: Property data
    flowsheet = "flowsheet"                 #: Process flowsheet
    notebook = "notebook"                   #: Jupyter Notebook
    code = "code"                           #: Source code(s)
    surrogate_model = "surrogate_model"     #: Surrogate model
    data = "data"                           #: Generic data
    other = "json"                          #: JSON data
    json = "other"                          #: User-defined type of resource
    resource_json = "resource_json"         #: JSON serialized resource

    @classmethod
    def all(cls) -> List[str]:
        """Return all resource type names, as a list of strings.
        """
        return [k for k in cls.__dict__ if k and (k != "all") and (k[0] != "_")]


# Constants for fields in stored relations
RR_PRED = "predicate"
RR_SUBJ = "subject"
RR_OBJ = "object"
RR_ID = "identifier"
RR_ROLE = "role"

#: Schema for a resource
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
                {"type": "string"},
            ],
            "minItems": 4,
        }
    },
    "type": "object",
    "properties": {
        "aliases": {"type": "array", "items": {"type": "string"}},
        "collaborators": {
            "type": "array",
            "items": {
                "type": "object",
                "properties": {
                    "email": {"type": "string", "format": "email"},
                    "name": {"type": "string"},
                },
                "required": ["name"],
            },
        },
        "created": {"type": "number"},
        "creator": {
            "type": "object",
            "properties": {
                "email": {"type": "string", "format": "email"},
                "name": {"type": "string"},
            },
            "required": ["name"],
        },
        "data": {"type": "object"},
        "datafiles": {
            "type": "array",
            "items": {
                "type": "object",
                "properties": {
                    "desc": {"type": "string"},
                    "metadata": {"type": "object"},
                    "mimetype": {"type": "string"},
                    "path": {"type": "string"},
                    "sha1": {"type": "string"},
                    "is_copy": {"type": "boolean"},
                },
                "required": ["path"],
            },
        },
        "datafiles_dir": {"type": "string"},
        "desc": {"type": "string"},
        "id_": {"type": "string"},
        "modified": {"type": "number"},
        "relations": {
            "type": "array",
            "items": {
                "type": "object",
                "properties": {
                    RR_PRED: {"type": "string", "enum": Predicates.all()},
                    RR_ID: {"type": "string"},
                    RR_ROLE: {"type": "string", "enum": [RR_SUBJ, RR_OBJ]},
                },
                "required": [RR_PRED, RR_ID, RR_ROLE],
            },
        },
        "tags": {"type": "array", "items": {"type": "string"}},
        "type": {"type": "string", "enum": ResourceTypes.all()},
        "version_info": {
            "type": "object",
            "properties": {
                "created": {"type": "number"},
                "name": {"type": "string"},
                "version": {"$ref": "#/definitions/SemanticVersion"},
            },
        },
        # These are only required for certain values of `type`
        # For type = code
        "codes": {
            "type": "array",
            "items": {
                "type": "object",
                "properties": {
                    "type": {
                        "type": "string",
                        "enum": [
                            "method",
                            "function",
                            "module",
                            "class",
                            "file",
                            "package",
                            "repository",
                            "notebook",
                            "block"
                        ],
                    },
                    "desc": {"type": "string"},
                    "name": {"type": "string"},
                    "language": {"type": "string"},
                    "idhash": {"type": "string"},
                    "location": {"type": "string"},
                    "inline": {"type": "string"},
                    "version": {"$ref": "#/definitions/SemanticVersion"},
                },
                "required": ["name"],
            },
        },
        # For type = publication
        "sources": {
            "type": "array",
            "items": {
                "type": "object",
                "properties": {
                    "date": {"type": "number"},
                    "doi": {"type": "string"},
                    "isbn": {"type": "string"},
                    "language": {"type": "string"},
                    "source": {"type": "string"},
                }
            }
        }
    },
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

    ID_FIELD = "id_"  #: Identifier field name constant
    ID_LENGTH = 32  #: Full-length of identifier
    TYPE_FIELD = "type"  #: Resource type field name constant

    def __init__(self, value: dict = None, type_: str = None):
        if type_ is None:
            type_ = ResourceTypes.other
        self._set_defaults(type_)
        if value:
            self.set_values(value)
        self._validator = jsonschema.Draft4Validator(RESOURCE_SCHEMA)
        self._validations = 0  # count validations; mostly for testing
        self.do_copy = self.is_tmp = False  # flags for copying datafiles

    def _set_defaults(self, t):
        now = datetime.now().timestamp()
        self.v = Dict(
            {
                self.ID_FIELD: identifier_str(),
                self.TYPE_FIELD: t,
                "aliases": [],
                "collaborators": [],
                "created": now,
                "modified": now,
                "creator": {"name": getpass.getuser()},
                "data": {},
                "datafiles": [],
                "datafiles_dir": "",
                "desc": "",
                "relations": [],
                "tags": [],
                "version_info": {"created": now, "version": (0, 0, 0), "name": ""},
            }
        )
        if t == ResourceTypes.code:
            self.v["codes"] = []
        elif t == ResourceTypes.publication:
            self.v["sources"] = []

    @property
    def sources(self):
        return self.v.get("sources", [])

    @property
    def codes(self):
        return self.v.get("codes", [])

    def _massage_values(self):
        try:
            # convert dates
            for item in self.sources:
                if not isinstance(item["date"], float):
                    item["date"] = date_float(item["date"])
            if not isinstance(self.v["created"], float):
                self.v["created"] = date_float(self.v["created"])
            if not isinstance(self.v["modified"], float):
                self.v["modified"] = date_float(self.v["modified"])
            if not isinstance(self.v["version_info"]["created"], float):
                self.v["version_info"]["created"] = date_float(
                    self.v["version_info"]["created"]
                )
            # convert versions
            if not isinstance(self.v["version_info"]["version"], list):
                self.v["version_info"]["version"] = version_list(
                    self.v["version_info"]["version"]
                )
            for i, code in enumerate(self.codes):
                if "version" in code:
                    if not isinstance(code["version"], list):
                        code["version"] = version_list(code["version"])
                        self.v["codes"][i] = code
        except (TypeError, ValueError, KeyError) as err:
            raise ValueError("While converting resource values: {}".format(err))
        self.v.set_clean()

    def validate(self):
        if self.v.is_dirty():
            self._massage_values()
            self._validator.validate(self.v)
            self._validations += 1

    # Some exceptions for communicating problems on import
    class InferResourceTypeError(Exception):
        pass

    class LoadResourceError(Exception):
        def __init__(self, inferred_type, msg):
            super().__init__(f"resource type '{inferred_type}': {msg}")

    @classmethod
    def from_file(
        cls, path: str, as_type: str = None, strict: bool = True, do_copy: bool = True
    ) -> "Resource":
        """Import resource from a file.

        Args:
            path: File path
            as_type: Resource type. If None/empty, then inferred from path.
            strict: If True, fail when file extension and contents don't match.
                    If False, always fall through to generic resource.
            do_copy: If True (the default), copy the files; else do not

        Raises:
            InferResourceTypeError: if resource type does not match inferred/specified
            LoadResourceError: if resource import failed
        """
        path = Path(path)
        if as_type:
            if as_type == ResourceTypes.json:  # make sure resources validate
                try:
                    parsed = json.load(path.open())
                    jsonschema.Draft4Validator(RESOURCE_SCHEMA).validate(parsed)
                except (UnicodeDecodeError, JSONDecodeError):
                    raise ValueError("Resource is not well-formed JSON")
                except jsonschema.ValidationError as err:
                    raise ValueError(f"Resource does not match schema: {err}")
            else:
                parsed = None
        else:
            as_type, parsed = cls._infer_resource_type(path, strict)
        importer = cls._get_resource_importer(
            as_type, path, parsed=parsed, do_copy=do_copy
        )
        return importer.create()

    @classmethod
    def _infer_resource_type(cls, path: Path, strict: bool):
        default_type = ResourceTypes.other
        try:
            if path.suffix == ".ipynb":
                return ResourceTypes.notebook, None
            if path.suffix == ".py":
                return ResourceTypes.code, None
            if path.suffix == ".json":
                max_bytes = 1e6  # arbitrary limit
                # over max_bytes? generic
                file_size = path.stat().st_size
                if file_size > max_bytes:
                    _log.warning(
                        f"Not attempting to parse JSON, file size "
                        f"{file_size} > {max_bytes}"
                    )
                    return default_type, None
                # see if it's a Resource
                try:
                    parsed = json.load(path.open())
                except (UnicodeDecodeError, JSONDecodeError):
                    raise ValueError("File ending in '.json' is not valid JSON")
                try:
                    jsonschema.Draft4Validator(RESOURCE_SCHEMA).validate(parsed)
                except jsonschema.ValidationError:
                    return ResourceTypes.json, parsed  # generic JSON
                return ResourceTypes.resource_json, parsed
        except ValueError as err:
            if strict:
                raise cls.InferResourceTypeError(str(err))
            _log.warning(f"{err}: treating as generic file")
        return default_type, None

    @classmethod
    def _get_resource_importer(
        cls, type_: str, path: Path, parsed=None, **kwargs
    ) -> "ResourceImporter":
        E = cls.LoadResourceError  # alias for exception raised
        if type_ == ResourceTypes.notebook:
            try:
                nb = json.load(open(str(path)))
            except (UnicodeDecodeError, JSONDecodeError):
                raise E(ResourceTypes.notebook, "not valid JSON")
            for key in "cells", "metadata", "nbformat":
                if key not in nb:
                    raise E(ResourceTypes.notebook, f"missing key: {key}")
            return JupyterNotebookImporter(path, **kwargs)
        if type_ == ResourceTypes.code:
            language = ProgLangExt.guess(path.suffix, default="unknown")
            return CodeImporter(path, language, **kwargs)
        if type_ == ResourceTypes.json:
            return JsonFileImporter(path, **kwargs)
        if type_ == ResourceTypes.resource_json:
            return SerializedResourceImporter(path, parsed, **kwargs)
        return FileImporter(path, **kwargs)

    @property
    def id(self):
        """Get resource identifier.
        """
        return self.v[self.ID_FIELD]

    def set_id(self, value=None):
        self.v[self.ID_FIELD] = identifier_str(value)

    def set_field(self, key: str, value):
        """Set field/value pairs with some special shortcuts.

        Args:
            key: Field name. Special values:
                * "name": Set as first in the list of aliases
                * "author": Set as first in the list of collaborators
            value: Field value.

        Returns:
            None
        """
        if key == "name":
            self.v["aliases"].insert(0, value)
        elif key == "author":
            self.v["collaborators"].insert(0, value)
        else:
            self.v[key] = value

    def set_values(self, values):
        for k, v in values.items():
            self.set_field(k, v)

    @property
    def name(self):
        """Get resource name (first alias).
        """
        try:
            nm = self.v["aliases"][0]
        except IndexError:
            nm = ""
        return nm

    @property
    def type(self):
        """Get resource type.
        """
        return self.v[self.TYPE_FIELD]

    @property
    def data(self):
        """Get JSON data for this resource.
        """
        return self.v["data"]

    @data.setter
    def data(self, value):
        """Set JSON data for this resource.
        """
        self.v["data"] = value

    def get_datafiles(
        self, mode: Optional[str] = None, ignore_errors: bool = False
    ) -> Iterator[Union[Path, IO]]:
        """Generate paths or file objects for 'datafiles' in resource.

        Args:
            mode: If given, call `open(mode)` on the path and return a file object. Otherwise, return a Path object.
            ignore_errors: If true, ignore failures on `open(mode)` of a path. This will only have an effect if
                  mode is not None (otherwise, no attempt is made to open the paths).
        Returns:
            generator: Generates Path or file objects, depending on mode
        """
        df_dir = self.v["datafiles_dir"]
        for datafile in self.v["datafiles"]:
            if not df_dir:
                path = Path(datafile["path"])
            else:
                df = Path(datafile["path"])
                if df.is_absolute():
                    path = df
                else:
                    path = Path(df_dir) / df
            if mode is None:
                yield path
            else:
                try:
                    fp = path.open(mode=mode)
                except FileNotFoundError:
                    if ignore_errors:
                        _log.warning("Failed to open path '{path}' in mode '{mode}'")
                    else:
                        raise
                yield fp

    def _repr_text_(self):
        return pprint.pformat(self.v, indent=2)

    def formatted_source(self) -> str:
        result = []
        for src in self.sources:
            s = f"{src['source']}"
            if src["isbn"]:
                s += f" ISBN: {src['isbn']}"
            if src["date"]:
                s += f" Date: {src['date']}"
            result.append(s)
        return "\n".join(result)


#
# Function(s) to help creating [two-way] relations
# between resources
#


#: Provide attribute access to an RDF subject, predicate, object triple
Triple = namedtuple("Triple", "subject predicate object")


def create_relation(*args):
    """Create a relationship between two Resource instances.

    Relations are stored in both the `subject` and `object` resources, in
    the following way::

        If R = (subject)S, (predicate)P, and (object)O
        then store the following:
          In S.relations: {predicate: P, identifier:O.id, role:subject}
          In O.relations: {predicate: P, identifier:S.id, role:object}

    Args:
        args: Either a `Triple`, or three arguments that can be made into one.
              The 'subject' and 'object' parts should be :class:`Resource`, and
              the 'predicate' should be a simple string.
    Returns:
        None
    Raises:
        ValueError: if this relation already exists in the subject or
                    object resource, or the predicate is not in the list
                    of valid ones in Predicates
    """
    if len(args) == 1:
        triple = args[0]
    elif len(args) == 3:
        triple = Triple(*args)
    else:
        raise ValueError(
            "Wrong number of arguments. Must be one argument (Triple) or "
            "three arguments (subject:Resource, predicate, object:resource)"
        )
    return _create_relation(triple)


def _create_relation(rel):
    if not Predicates.valid(rel.predicate):
        raise ValueError(
            'Bad predicate: "{}" not in: {}'.format(
                rel.predicate, ", ".join(Predicates.all())
            )
        )
    try:
        _ = rel.subject.v[Resource.ID_FIELD], rel.object.v[Resource.ID_FIELD]
    except AttributeError:
        raise TypeError("Relation subject and object must be Resources")

    rel_d = {
        RR_PRED: rel.predicate,
        RR_ID: rel.object.v[Resource.ID_FIELD],
        RR_ROLE: RR_SUBJ,
    }
    if rel_d in rel.subject.v["relations"]:
        raise ValueError("Duplicate relation for subject: {}".format(rel))
    rel.subject.v["relations"].append(rel_d)
    rel_d = {
        RR_PRED: rel.predicate,
        RR_ID: rel.subject.v[Resource.ID_FIELD],
        RR_ROLE: RR_OBJ,
    }
    # note: hard for this to happen unless the relation was added manually
    if rel_d in rel.object.v["relations"]:
        raise ValueError("Duplicate relation for object: {}".format(rel))
    rel.object.v["relations"].append(rel_d)


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
    """Convert a date to a floating point seconds since the UNIX epoch.
    """

    def bad_date(e):
        raise ValueError('Cannot convert date "{}" to float: {}'.format(value, e))

    dt, usec = None, 0
    if isinstance(value, datetime):
        dt = value
    elif isinstance(value, tuple):
        try:
            dt = datetime(*value)
        except TypeError as err:
            bad_date(err)
    elif isinstance(value, str):
        try:
            dt = parse_datetime(value)
        except ValueError as err:
            bad_date(err)
    elif isinstance(value, float) or isinstance(value, int):
        try:
            dt = datetime.fromtimestamp(value)
        except ValueError as err:
            bad_date(err)
    if dt is None:
        raise ValueError('Cannot convert date, value is "None"')
    return datetime_timestamp(dt) + usec  # just a float


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

    .. testsetup:: version_list

        from idaes.dmf.resource import version_list

    .. doctest:: version_list

        >>> version_list('1')
        [1, 0, 0, '']
        >>> version_list('1.1')
        [1, 1, 0, '']
        >>> version_list('1a')
        [1, 0, 0, 'a']
        >>> version_list('1.12.1')
        [1, 12, 1, '']
        >>> version_list('1.12.13-1')
        [1, 12, 13, '1']

    Some examples of invalid inputs:

    .. doctest:: version_list

        >>> for bad_input in ('rc3',      # too short
        ...                   '1.a.1.',   # non-number in middle
        ...                   '1.12.13.x' # too long
        ...     ):
        ...     try:
        ...         version_list(bad_input)
        ...     except ValueError:
        ...         print(f"failed: {bad_input}")
        ...
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
        ver = value.split(".", 2)
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
    extra = ""
    if isinstance(s, int):
        verlist.append(s if len(verlist) < 3 else str(s))
    elif isinstance(s, str):
        if s:
            m = re.match("([0-9]+)?(.*)", s)
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
    if extra and extra[0] == ".":
        # cannot start extra version with '.'
        bad_version(value)
    if extra and extra[0] in ("-", "_"):
        extra = extra[1:]
    verlist.append(extra)
    return verlist


def format_version(values):
    s = "{}.{}.{}".format(*values[:3])
    if len(values) > 3 and values[3]:
        s += "-{}".format(values[3])
    return s


def identifier_str(value=None, allow_prefix=False):
    """Generate or validate a unique identifier.

    If generating, you will get a UUID in hex format

    .. testsetup:: idstr

        from idaes.dmf.resource import identifier_str

    .. doctest:: idstr

        >>> identifier_str()  #doctest: +ELLIPSIS
        '...'

    If validating, anything that is not 32 lowercase letters
    or digits will fail.

    .. doctest:: idstr

        >>> identifier_str('A' * 32)   #doctest: +NORMALIZE_WHITESPACE
        Traceback (most recent call last):
        ValueError: Bad format for identifier "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA":
        must match regular expression "[0-9a-f]{32}"

    Args:
        value (str): If given, validate that it is a 32-byte str
                     If not given or None, set new random value.
    Raises:
        ValuError, if a value is given, and it is invalid.
    """
    # regular expression for identifier: hex string len=32
    if allow_prefix:
        id_expr = "[0-9a-f]{1,32}"
    else:
        id_expr = "[0-9a-f]{32}"
    if value is None:
        value = uuid.uuid4().hex
    elif not re.match(id_expr, value):
        raise ValueError(
            'bad format for identifier "{}": must match '
            'regular expression "{}"'.format(value, id_expr)
        )
    return value


class TidyUnitData:
    """Handle "tidy data" with per-column units.

    This can be used to convert from a simple dictionary/json
    representation like this::

            {
              "variables": ["compound", "pressure"],
              "units": [null|None, "Pa"],
              "observations": [
                ["benzene", 4890000.0],
                ...etc..
              ]
            }

    into a pandas DataFrame. A convenience method is provided for returning
    the data in a format easily dealt with when creating unit block parameters.
    Note that the keys in the preceding dictionary match the names of the
    parameters in the constructor (so you can pass this directly in as `**arg`).

    Attributes:
        units (list): Units for each column, None where no units are defined
        table (pandas.DataFrame): The observation data
    """

    def __init__(
        self,
        data: dict = None,
        variables: List = None,
        units: List = None,
        observations: List = None,
    ):
        """Constructor.

        Args:
            data: Optional, data dict (overrides variables/units/observations)
            variables: List of variables (table header)
            units: List of units, or None, same length as `variables`
            observations: Rows of the body of the table
        Raises:
            ValueError: For bad `data` (missing keys, not dict, etc.), or mismatches
                        in lengths of various pieces.
        """
        if data:
            try:
                variables, units, observations = (
                    data["variables"],
                    data["units"],
                    data["observations"],
                )
            except KeyError as err:
                raise ValueError(f"Missing expected key in `data` param: {err}")
            except TypeError as err:
                raise ValueError(f"Bad value for `data` param: {err}")
        n = len(variables)
        if n == 0:
            self.df, self.units = pandas.DataFrame(), ()
            return
        if len(units) != n:
            raise ValueError(
                f"Length of units {len(units)} " f"must match length of header ({n})"
            )
        self.units = units
        self.table = pandas.DataFrame(data=observations, columns=variables)

    @property
    def param_data(self) -> dict:
        """Data in a form easily consumed by unit block params.

        The dictionary returned is like ``{ (key1, key2, ..): value }``,
        where the keys are values from all columns except the last,
        and the value is the last column.
        """
        d = {}
        for row in self.table.itertuples():
            key = tuple(row[1:-1])
            if len(key) == 1:
                key = key[0]
            value = row[-1]
            d[key] = value
        return d


#
# Import Resource of varying types from file
#


class ResourceImporter(abc.ABC):
    """Base class for Resource importers.
    """

    def __init__(self, path: Path, do_copy: bool = None):
        self._path = path
        self._do_copy = do_copy

    def create(self) -> Resource:
        """Factory method.
        """
        r = self._create()
        r.validate()
        return r

    @abc.abstractmethod
    def _create(self) -> Resource:
        pass

    def _add_datafiles(self, r):
        abspath = str(self._path.absolute())
        file_hash = self._hash_file(abspath)
        r.v["datafiles"].append(
            {
                "desc": self._path.name,
                "path": abspath,
                "do_copy": self._do_copy,
                "sha1": file_hash,
            }
        )

    def _hash_file(self, path):
        blksz, h = 1 << 16, hashlib.sha1()
        with open(path, "rb") as f:
            blk = f.read(blksz)
            while blk:
                h.update(blk)
                blk = f.read(blksz)
        return h.hexdigest()


class JupyterNotebookImporter(ResourceImporter):
    def _create(self) -> Resource:
        r = Resource(type_=ResourceTypes.notebook)
        # XXX: add notebook 'metadata' as FilePath metadata attr
        self._add_datafiles(r)
        r.v["desc"] = self._path.name
        return r


class CodeImporter(ResourceImporter):
    def __init__(self, path, language, **kwargs):
        super().__init__(path, **kwargs)
        self.language = language

    def _create(self) -> Resource:
        r = Resource(type_=ResourceTypes.code)
        r.codes.append(
            {"name": self._path.name, "language": self.language, "type": "module"}
        )
        self._add_datafiles(r)
        r.v["desc"] = self._path.name
        return r


class FileImporter(ResourceImporter):
    def _create(self) -> Resource:
        r = Resource(type_=ResourceTypes.data)
        self._add_datafiles(r)
        r.v["desc"] = str(self._path)
        return r


class JsonFileImporter(ResourceImporter):
    def _create(self) -> Resource:
        r = Resource(type_=ResourceTypes.json)
        self._add_datafiles(r)
        r.v["desc"] = str(self._path)
        return r


class SerializedResourceImporter(ResourceImporter):
    def __init__(self, path, parsed, **kwargs):
        super().__init__(path, **kwargs)
        self.parsed = parsed

    def _create(self) -> Resource:
        r = Resource(value=self.parsed)
        return r


#
# Fill any JSON-schema-constrained instance with
# dummy values
#


def add_dummy_values(validator_class):
    validate_properties = validator_class.VALIDATORS["properties"]

    dummy_for_type = {
        "object": {},
        "boolean": False,
        "null": None,
        "number": 0,
        "integer": 0,
        "string": "",
        "any": "",
    }

    def set_defaults(validator, properties, instance, schema):
        for prop, subschema in properties.items():
            # print(f"@@ at prop={prop}")
            if "$ref" in subschema:
                refkey = subschema["$ref"].split("/")[-1]
                subschema = RESOURCE_SCHEMA["definitions"][refkey]
                # print(f"@@ ref {refkey} resolved, schema={subschema}")
            if "enum" in subschema:
                value = subschema["enum"][0]
            elif subschema["type"] == "array":
                # print("@@ array subschema")
                if "minItems" in subschema:
                    value = [
                        dummy_for_type[item["type"]] for item in subschema["items"]
                    ]
                elif subschema["items"]["type"] == "object":
                    value = [{}]
                else:
                    value = []
                # print(f"@@ array value = {value}")
            else:
                value = dummy_for_type[subschema["type"]]
            instance.setdefault(prop, value)
        for error in validate_properties(validator, properties, instance, schema):
            yield error

    return jsonschema.validators.extend(validator_class, {"properties": set_defaults})


# Create dummy resource

DummyValueValidator = add_dummy_values(jsonschema.Draft4Validator)
DUMMY_RESOURCE = {}
DummyValueValidator(RESOURCE_SCHEMA).validate(DUMMY_RESOURCE)
# print("@@ dummy resource:\n{}".format(DUMMY_RESOURCE))


def schema_as_yaml():
    """Export resource schema as YAML suitable for embedding into, e.g.,
       an OpenAPI spec.
    """
    return yaml.dump(RESOURCE_SCHEMA)


#
# Things to do if run as a script
#


if __name__ == "__main__":
    # parse command line
    ap = argparse.ArgumentParser()
    actions = {
        "json_schema": "print resource schema as JSON",
        "yaml_schema": "print resource schema as YAML",
    }
    ap.add_argument(
        "action", help="Action when run as a script", choices=list(actions.keys())
    )
    args = ap.parse_args()
    # take appropriate action
    if args.action == "json_schema":
        json.dump(RESOURCE_SCHEMA, sys.stdout, indent=2)
    elif args.action == "yaml_schema":
        print(schema_as_yaml())
    else:
        print("nothing to do")
    # exit
    sys.exit(0)
