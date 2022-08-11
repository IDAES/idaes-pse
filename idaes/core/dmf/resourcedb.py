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
Resource database.
"""
# system
from datetime import datetime
import logging
import re

# third party
from tinydb import TinyDB, Query

# local
from . import errors
from .resource import Resource
from .resource import Triple, triple_from_resource_relations

__author__ = "Dan Gunter <dkgunter@lbl.gov>"

_log = logging.getLogger(__name__)


class ResourceDB(object):
    """A database interface to all the resources within a given DMF workspace."""

    def __init__(self, dbfile=None, connection=None):
        """Initialize from DMF and given configuration field.

        Args:
            dbfile (str): DB location
            connection: If non-empty, this is an
                existing connection that should be re-used, instead of
                trying to connect to the location in `dbfile`.

        Raises:
            ValueError, if dbfile
             and connection are both None
        """
        self._db = None
        self._gr = None

        if connection is not None:
            self._db = connection
        elif dbfile is not None:
            try:
                db = TinyDB(dbfile)
            except IOError:
                raise errors.FileError('Cannot open resource DB "{}"'.format(dbfile))
            # turn off caching, otherwise update() does not work properly
            self._db = db.table("resources", cache_size=0)

    def __len__(self):
        return len(self._db)

    def find(self, filter_dict, id_only=False, flags=0):
        """Find and return records based on the provided filter.

        Args:
            filter_dict (dict): Search filter. For syntax, see docs in
                                :meth:`.dmf.DMF.find`.
            id_only (bool): If true, return only the identifier of each
                resource; otherwise a Resource object is returned.
            flags (int): Flag values for, e.g., regex searches

        Returns:
            generator of int|Resource, depending on the value of `id_only`
        """

        def as_resource(_r):
            _log.debug(f"as_resource: id_={_r['id_']}")
            rsrc = Resource(value=_r)
            rsrc.v["doc_id"] = _r.doc_id
            return rsrc

        # with no filter, do a find-all
        if not filter_dict:
            for r in self._db.all():
                if id_only:
                    yield r.doc_id
                else:
                    yield as_resource(r)
            return
        filter_expr = self._create_filter_expr(filter_dict, flags)

        # return results for query
        _log.debug("Find resources matching: {}".format(filter_expr))
        results = self._db.search(filter_expr)
        for r in results:
            if id_only:
                yield r.doc_id
            else:
                _log.debug(f"got resource: {r}")
                yield as_resource(r)

    @classmethod
    def _create_filter_expr(cls, filter_dict, flags=0):
        # convert filter to expression for TinyDB
        q = Query()
        filter_expr = None

        def fadd(c):
            """Update filter, or set to 1st condition, & return new value."""
            return c if filter_expr is None else filter_expr & c

        for k, v in filter_dict.items():
            if not k:
                continue  # XXX: Issue a warning?
            # strip off list-query operator
            qry_all = False
            if isinstance(v, list) and k.endswith("!"):
                k, qry_all = k[:-1], True
            # Get the query[field1][field2][..][fieldN] object
            qry = q
            for field in k.split("."):
                qry = qry[field]
            # A list value means find these value(s) in the
            # list value of the record. By default, any matching
            # value is a match. If the key has a '!' appended, then
            # require all values for a match.
            if isinstance(v, list):
                if len(v) > 0:
                    # if values are dicts, nest a query
                    # XXX: this only works one level deep
                    # XXX: only one nested query is used/allowed
                    if isinstance(v[0], dict):
                        list_expr = None
                        for list_k, list_v in v[0].items():
                            list_qry = Query()
                            for field in list_k.split("."):
                                list_qry = list_qry[field]
                            list_cond = cls._expr_to_query(
                                list_qry, list_v, flags=flags
                            )
                            if list_expr is None:
                                list_expr = list_cond
                            else:
                                list_expr = list_expr & list_cond
                    # otherwise, simply put the values in there
                    else:
                        list_expr = tuple(v)
                    # add nested "query" (or value tuple)
                    if qry_all:
                        cond = qry.all(list_expr)
                    else:
                        cond = qry.any(list_expr)
            else:
                cond = cls._expr_to_query(qry, v, flags=flags)
            filter_expr = fadd(cond)
        return filter_expr

    @classmethod
    def _expr_to_query(cls, qry, v, flags=None):
        """Get a query from a filter expr.

        There are two types of non-list values:
         1. equality is just {key: value}
         2. inequalities are {key: {op: value, ...}}
           operators are "$<shell test op>", such as "$lt",
           just like in MongoDB; see _op_cond() method for details.
        """
        result = None
        if isinstance(v, dict):
            for op_key, op_value in v.items():
                tv = cls._value_transform(op_value)
                cond = cls._op_cond(qry, op_key, tv)
                result = cond if result is None else result & cond
        else:
            tv = cls._value_transform(v)
            if v is True:
                result = qry.exists()
            elif tv is False:
                result = ~qry.exists()
            elif hasattr(tv, "match"):  # regex
                result = qry.matches(tv.pattern, flags=flags)
            else:
                result = qry == tv
        return result

    @staticmethod
    def _value_transform(v):
        # transform dates into timestamps
        if isinstance(v, datetime):
            return v.timestamp()
        # support for special string '@' values
        elif isinstance(v, str) and len(v) > 0 and v[0] == "@":
            if v == "@true":
                return True
            if v == "@false":
                return False
            return v
        # for a regex, return a Python regex obj
        elif isinstance(v, str) and len(v) > 0 and v[0] == "~":
            return re.compile(v[1:])
        # default is no transformation
        else:
            return v

    @staticmethod
    def _op_cond(query, op, value):
        # sanity check that operator is truthy
        if not op:
            raise ValueError(f"empty operator for value `{value}`")
        # just a clumsy switch statement..
        if op == "$gt":
            cond = query > value
        elif op == "$ge":
            cond = query >= value
        elif op == "$lt":
            cond = query < value
        elif op == "$le":
            cond = query <= value
        elif op == "$ne":
            cond = query != value
        else:
            raise ValueError("Unexpected operator: {}".format(op))
        return cond

    def find_one(self, *args, **kwargs):
        """Same as `find()`, but returning only first value or None."""
        result = None
        for value in self.find(*args, **kwargs):
            result = value
            break
        return result

    def find_related(self, id_, filter_dict=None, outgoing=True, maxdepth=0, meta=None):
        """Find all resources connected to the identified one.

        Args:
            id_ (str): Unique ID of target resource.
            filter_dict (dict): Filter to these resources
            outgoing:
            maxdepth:
            meta (List[str]): Metadata fields to extract
        Returns:
            Generator of (depth, relation, metadata)
        Raises:
            KeyError if the resource is not found.
        """
        if maxdepth <= 0:
            maxdepth = 9223372036854775807
        # Get an iterator over all the resources, optionally
        # filtered by an expression, as for find().
        if filter_dict:
            filter_expr = self._create_filter_expr(filter_dict)
            resource_list = self._db.search(filter_expr)
        else:
            resource_list = self._db.all()
        # build adjacency list representing connections between resources
        relation_map = {}
        for rsrc in resource_list:
            # print(f"@@ process resource: {rsrc[Resource.ID_FIELD]}")
            # print(f"@@ resource relations: {rsrc['relations']}")
            for rrel in rsrc["relations"]:
                uuid = rsrc[Resource.ID_FIELD]
                rel = triple_from_resource_relations(uuid, rrel)
                # skip start of edge, get metadata, etc. at end of edge
                if (outgoing and rel.subject == uuid) or (
                    not outgoing and rel.object == uuid
                ):
                    # print(f"@@ do nothing for: {uuid}")
                    continue
                # add entry in map
                key = rel.subject if outgoing else rel.object
                # print(f"@@ add to relation map, key={key}")
                meta_info = {k: rsrc[k] for k in meta}
                value = (rel.subject, rel.predicate, rel.object, meta_info)
                value_list = relation_map.get(key, None)
                if value_list is None:
                    relation_map[key] = [value]
                else:
                    value_list.append(value)
        _log.debug(f"built relation map: {relation_map}")
        # stop if there are no connections
        if id_ not in relation_map:
            return
        # Do a depth-first search through the edges, yield-ing
        # the relations as we go
        q, depth, visited = [], 0, set()
        q.extend(relation_map[id_])
        visited.add(id_)
        while len(q) > 0 and depth < maxdepth:
            depth += 1
            # visit all the nodes in the queue
            n = len(q)
            for i in range(n):
                relation = Triple(*q[i][:3])
                yield (depth, relation, q[i][3])
                if depth < maxdepth:
                    # Follow relations from subject or object, depending on
                    # the "direction" that we are searching.
                    next_id = relation.object if outgoing else relation.subject
                    # If there are relations, and we haven't already been to
                    # this node, add them at the end of the queue; we will
                    # visit them at the next depth increment.
                    if next_id in relation_map and next_id not in visited:
                        q.extend(relation_map[next_id])
                        visited.add(next_id)
            q = q[n:]  # pop off all the nodes we just visited

    def get(self, identifier):
        """Get a resource by identifier.

        Args:
          identifier: Internal identifier

        Returns:
            (Resource) A resource or None
        """

        def as_resource(_r):
            rsrc = Resource(value=_r)
            rsrc.v["doc_id"] = _r.doc_id
            return rsrc

        item = self._db.get(doc_id=identifier)
        if item is None:
            return None
        return as_resource(item)

    def put(self, resource):
        """Put this resource into the database.

        Args:
            resource (Resource): The resource to add

        Returns:
            None

        Raises:
            errors.DuplicateResourceError: If there is already a resource
                in the database with the same "id".
        """
        _log.debug(f"put resource id={resource.id}")
        # check for same id
        qry = Query()
        if self._db.contains(qry.id_ == resource.id):
            raise errors.DuplicateResourceError("put", resource.id)
        # add resource
        self._db.insert(resource.v)

    def delete(self, id_=None, idlist=None, filter_dict=None, internal_ids=False):
        """Delete one or more resources with given identifiers.

        Args:
            id_ (Union[str,int]): If given, delete this id.
            idlist (list): If given, delete ids in this list
            filter_dict (dict): If given, perform a search and
                           delete ids it finds.
            internal_ids (bool): If True, treat identifiers as numeric
                (internal) identifiers. Otherwise treat them as
                resource (string) indentifiers.
        Returns:
            (list[str]) Identifiers
        """
        if internal_ids:
            doc_ids = idlist if idlist else [id_]
            self._db.remove(doc_ids=doc_ids)
        else:
            ID = Resource.ID_FIELD
            if filter_dict:
                cond = self._create_filter_expr(filter_dict)
            elif id_:
                cond = self._create_filter_expr({ID: id_})
            elif idlist:
                cond = self._create_filter_expr({ID: [idlist]})
            else:
                return
            self._db.remove(cond=cond)

    def update(self, id_, new_dict):
        """Update the identified resource with new values.

        Args:
            id_ (int): Identifier of resource to update
            new_dict (dict): New dictionary of resource values
        Returns:
            None
        Raises:
            ValueError: If new resource is of wrong type
            KeyError: If old resource is not found
        """
        _log.debug("update.start")
        id_cond = {Resource.ID_FIELD: id_}
        old = self.find_one(id_cond)
        if old is None:
            raise errors.NoSuchResourceError(id_=id_)
        T = Resource.TYPE_FIELD
        if old.v[T] != new_dict[T]:
            raise ValueError(
                'New resource type="{}" does not '
                'match current resource type "{}"'.format(new_dict[T], old.v[T])
            )
        _log.debug(f"update old resource: {old.v}")
        changed = {}
        for k, v in new_dict.items():
            if k not in old.v:
                changed[k] = v
            elif old.v[k] != v:
                changed[k] = v
        _log.debug(f"update resource {id_} with new values: {changed}")
        self._db.update(changed, self._create_filter_expr(id_cond))
