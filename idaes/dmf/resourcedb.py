##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes".
##############################################################################
"""
Resource database.
"""
# system
from datetime import datetime
import logging
# third party
import pendulum
import six
from tinydb import TinyDB, Query
# local
from . import errors
from .resource import Resource
from .resource import Triple, triple_from_resource_relations

__author__ = 'Dan Gunter <dkgunter@lbl.gov>'

_log = logging.getLogger(__name__)


class ResourceDB(object):
    """A database interface to all the resources within a given DMF workspace.
    """
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
                raise errors.FileError('Cannot open resource DB "{}"'
                                       .format(dbfile))
            self._db = db

    def __len__(self):
        return len(self._db)

    def find(self, filter_dict, id_only=False):
        """Find and return records based on the provided filter.

        Args:
            filter_dict (dict): Search filter. For syntax, see docs in
                                :meth:`.dmf.DMF.find`.
            id_only (bool): If true, return only the identifier of each
                resource; otherwise a Resource object is returned.

        Returns:
            (list of int|Resource) Depending on the value of `id_only`
        """
        def as_resource(_r):
            rsrc = Resource(value=_r)
            rsrc.v['doc_id'] = _r.doc_id
            return rsrc
        # with no filter, do a find-all
        if not filter_dict:
            for r in self._db.all():
                if id_only:
                    yield r.eid
                else:
                    yield as_resource(r)
            return
        filter_expr = self._create_filter_expr(filter_dict)

        # return results for query
        _log.debug('Find resources matching: {}'.format(filter_expr))
        results = self._db.search(filter_expr)
        for r in results:
            if id_only:
                yield r.eid
            else:
                yield as_resource(r)

    @classmethod
    def _create_filter_expr(cls, filter_dict):
        # convert filter to expression for TinyDB
        q = Query()
        filter_expr = None

        def fadd(c):
            """Update filter, or set to 1st condition, & return new value."""
            return c if filter_expr is None else filter_expr & c

        for k, v in six.iteritems(filter_dict):
            if not k:
                continue  # XXX: Issue a warning?
            # strip off list-query operator
            qry_all = False
            if isinstance(v, list) and k.endswith('!'):
                k, qry_all = k[:-1], True
            # Get the query[field1][field2][..][fieldN] object
            qry = q
            for field in k.split('.'):
                qry = qry[field]
            # A list value means find these value(s) in the
            # list value of the record. By default, any matching
            # value is a match. If the key has a '!' appended, then
            # require all values for a match.
            if isinstance(v, list):
                if qry_all:
                    cond = qry.all(tuple(v))
                else:
                    cond = qry.any(tuple(v))
                filter_expr = fadd(cond)
            else:
                # There are two types of non-list values:
                # - equality is just {key: value}
                # - inequalities are {key: {op: value, ...}}
                #   operators are "$<shell test op>", such as "$lt",
                #   just like in MongoDB; see _op_cond() method for details.
                if isinstance(v, dict):
                    for op_key, op_value in six.iteritems(v):
                        tv = cls._value_transform(op_value)
                        cond = cls._op_cond(qry, op_key, tv)
                        filter_expr = fadd(cond)
                elif v is True:
                    cond = qry.exists()
                    filter_expr = fadd(cond)
                elif v is False:
                    cond = ~ qry.exists()
                    filter_expr = fadd(cond)
                else:
                    tv = cls._value_transform(v)
                    cond = qry == tv
                    filter_expr = fadd(cond)
        return filter_expr

    @staticmethod
    def _value_transform(v):
        if isinstance(v, datetime) or isinstance(v, pendulum.Pendulum):
            if isinstance(v, datetime):
                pv = pendulum.create(v.year, v.month, v.day, v.hour, v.minute,
                                     v.second, v.microsecond, v.tzname())
            else:
                pv = v
            return pv.timestamp()
        else:
            return v

    @staticmethod
    def _op_cond(query, op, value):
        # just a clumsy switch statement..
        if op == '$gt':
            cond = query > value
        elif op == '$ge':
            cond = query >= value
        elif op == '$lt':
            cond = query < value
        elif op == '$le':
            cond = query <= value
        elif op == '$ne':
            cond = query != value
        else:
            raise ValueError('Unexpected operator: {}'.format(op))
        return cond

    def find_one(self, *args, **kwargs):
        """Same as `find()`, but returning only first value or None.
        """
        result = None
        for value in self.find(*args, **kwargs):
            result = value
            break
        return result

    def find_related(self, id_, filter_dict=None, outgoing=True,
                     maxdepth=0, meta=None):
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
            for rrel in rsrc['relations']:
                uuid = rsrc[Resource.ID_FIELD]
                rel = triple_from_resource_relations(uuid, rrel)
                # skip start of edge, get metadata, etc. at end of edge
                if (outgoing and rel.subject == uuid) or \
                        (not outgoing and rel.object == uuid):
                    continue
                # add entry in map
                key = rel.subject if outgoing else rel.object
                meta_info = {k: rsrc[k] for k in meta}
                value = (rel.subject, rel.predicate, rel.object, meta_info)
                value_list = relation_map.get(key, None)
                if value_list is None:
                    relation_map[key] = [value]
                else:
                    value_list.append(value)
        # print('@@ built relation map: {}'.format(relation_map))
        # make sure root node exists
        if id_ not in relation_map:
            raise KeyError('Resource not found: {}'.format(id_))
        # Do a depth-first search through the edges, yield-ing
        # the relations as we go
        q, depth, visited = [], 0, set()
        q.extend(relation_map[id_])
        visited.add(id_)
        while len(q) > 0 and depth < maxdepth:
            depth += 1
            # print('@@ depth={:d}. Queue: {}'.format(depth, q))
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
            rsrc.v['doc_id'] = _r.doc_id
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
        # check for same id
        qry = Query()
        if self._db.contains(qry.id_ == resource.id):
            raise errors.DuplicateResourceError("put", resource.id)
        # add resource
        self._db.insert(resource.v)

    def delete(self, id_=None, idlist=None, filter_dict=None,
               internal_ids=False):
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
        id_cond = {Resource.ID_FIELD: id_}
        old = self.find_one(id_cond)
        if old is None:
            raise KeyError('Cannot find resource id={}'.format(id_))
        T = Resource.TYPE_FIELD
        if old.v[T] != new_dict[T]:
            raise ValueError('New resource type="{}" does not '
                             'match current resource type "{}"'
                             .format(new_dict[T], old.v[T]))

        changed = {k: v for k, v in six.iteritems(new_dict)
                   if old.v.get(k, None) != v}
        # note: above does not check for missing/new keys
        self._db.update(changed, self._create_filter_expr(id_cond))
