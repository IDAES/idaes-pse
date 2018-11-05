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
Data Management Framework
"""
# stdlib
import logging
import os
import shutil
import sys
import uuid
# third-party
import pendulum
import six
from traitlets import HasTraits, default, observe
from traitlets import Unicode
import yaml
# local
from . import errors, workspace, resourcedb, propdata
from . import resource
from .util import get_logger, mkdir_p

_log = get_logger('dmf')


class DMFConfig(object):
    """Global DMF configuration.

    Every time you create an instance of the :class:`DMF <idaes.dmf.dmf.DMF>`,
    or run a ``dmf`` command on the command-line, the library opens the
    global DMF configuration file to figure out wthe default workspace
    (and, eventually, other values).

    The default location for this configuration file is "~/.dmf", i.e.
    the file named ".dmf" in the user's home directory. This can be
    modified programmatically by changing the "filename" attribute of
    this class.

    The contents of .dmf are formatted as `YAML`_, with the following
    keys defined:

        workspace
            Path to the default workspace directory.

    An example file is shown below:

    .. code-block:: YAML

        {workspace: /tmp/newdir}

    .. _YAML: http://www.yaml.org/
    """
    # configuration file location
    filename = os.path.expanduser('~/.dmf')
    # known keys
    WORKSPACE = 'workspace'
    _keys = [WORKSPACE]
    # defaults
    DEFAULTS = {WORKSPACE: os.path.abspath('.')}

    def __init__(self, defaults=None):
        """Create new configuration.

        Args:
            defaults (dict): Default values to use if nothing is found.
                            If the value is 'True', use self.DEFAULTS
        Raises:
            IOError: If self.filename is not openable or parse-able
        """
        self.c, fp = {}, None
        try:
            fp = open(self.filename, 'rb')
        except IOError:
            if defaults is None:
                raise
            if defaults is True:
                self.c.update(self.DEFAULTS)
            else:
                self.c.update(defaults)
        if fp:
            try:
                self._parse(fp)
            except ValueError as err:
                raise IOError('Parsing configuration at "{}": {}'
                              .format(self.filename, err))

    def save(self):
        try:
            fp = open(self.filename, 'w')
        except IOError:
            raise
        try:
            self._write(fp)
        except ValueError as err:
            raise IOError('Failed to write in YAML format to <{}>: {}'
                          .format(self.filename, err))

    @property
    def workspace(self):
        return self.c[self.WORKSPACE]

    def _parse(self, fp):
        try:
            y = yaml.load(fp)
        except Exception as err:
            raise ValueError(err)
        if y:
            for k, v in six.iteritems(y):
                if k in self._keys:
                    self.c[k] = v

    def _write(self, fp):
        yaml.dump(self.c, fp)


class DMF(workspace.Workspace, HasTraits):
    """Data Management Framework (DMF).

    Expected usage is to instantiate this class, once, and then use it for
    storing, searching, and retrieve :term:`resource` \s that
    are required for the given analysis.

    For details on the configuration files used by the DMF, see
    documentation for :class:`DMFConfig` (global configuration) and
    :class:`idaes.dmf.workspace.Workspace`.
    """
    db_file = Unicode(help='Database file name')
    datafile_dir = Unicode(help='Data file directory, '
                                'relative to DMF root')

    CONF_DB_FILE = 'db_file'
    CONF_DATA_DIR = 'datafile_dir'
    CONF_HELP_PATH = workspace.Fields.DOC_HTML_PATH

    # logging should really provide this
    _levelnames = {'fatal': logging.FATAL, 'error': logging.ERROR,
                   'warn': logging.WARN, 'warning': logging.WARN,
                   'info': logging.INFO, 'debug': logging.DEBUG}

    def __init__(self, path='', name=None, desc=None, create=False,
                 **ws_kwargs):
        """Create or load DMF workspace.

        Args:
            path (str): Path to workspace. If given, will override any
                  global configuration. If not given (empty or None), then
                  global configuration will be checked first. If empty, and
                  the global configuration does not parse or exist, then
                  "." will be used.
            name (str): Name to be used for workspace.
            desc (str): Longer description of workspace.
            create (bool): If the path to the workspace does not exist,
                           this controls whether to create it or raise
                           an error.
            **ws_kwargs: Keyword arguments for :meth:`workspace.Workspace()`
                         constructor.
        Raises:
            errors.DMFWorkspaceNotFoundError: If workspace is not found, and
                 the `create` option was not given in the keyword args.
            errors.DMFBadWorkspaceError: If workspace is found, but there
                 are errors with its configuration.
        """
        # get global configuration
        conf = DMFConfig(defaults={})
        # get path, if not specified, from configuration
        if not path:
            path = conf.c.get(DMFConfig.WORKSPACE, '.')
        # set up workspace
        ws_kwargs['create'] = create
        try:
            workspace.Workspace.__init__(self, path, **ws_kwargs)
        except OSError as err:
            raise errors.DMFBadWorkspaceError(path,
                                              'Cannot create: {}'.format(err))
        except errors.WorkspaceNotFoundError:
            raise errors.DMFWorkspaceNotFoundError(path)
        except errors.WorkspaceConfNotFoundError:
            raise errors.DMFBadWorkspaceError(path, 'Configuration not found')
        except errors.WorkspaceConfMissingField as err:
            raise errors.DMFBadWorkspaceError(path, str(err))
        try:
            meta_dict = self.meta  # call associated code
        except (errors.ParseError, ValueError) as err:
            msg = 'Configuration parse error: {}'.format(err)
            raise errors.DMFBadWorkspaceError(path, msg)
        try:
            self._validate_conf(meta_dict)
        except ValueError as err:
            msg = 'Configuration validation error: {}'.format(err)
            raise errors.DMFBadWorkspaceError(path, msg)
        # set up logging
        if workspace.Fields.LOG_CONF in meta_dict:
            try:
                self._configure_logging(meta_dict[workspace.Fields.LOG_CONF])
            except ValueError as err:
                msg = 'Configuration, logging section, error: {}'.format(err)
                raise errors.DMFBadWorkspaceError(path, msg)
        # set up rest of DMF
        path = os.path.join(self.root, self.db_file)
        self._db = resourcedb.ResourceDB(path)
        self._datafile_path = os.path.join(self.root, self.datafile_dir)
        if not os.path.exists(self._datafile_path):
            os.mkdir(self._datafile_path, 0o750)
        # add create/modified date, and optional name/description
        _w = workspace.Workspace
        right_now = pendulum.now().to_datetime_string()
        meta = {_w.CONF_CREATED: right_now, _w.CONF_MODIFIED: right_now}
        if name:
            meta[_w.CONF_NAME] = name
        if desc:
            meta[_w.CONF_DESC] = desc
        self.set_meta(meta)

    def _configure_logging(self, conf):
        """Configure logging for DMF.

        Expected schema::

            <logger-name>:
                level: <levelname>              default=NOTSET
                output: file|_stdout_|_stderr_  default=_stderr_

        Args:
            conf (dict): Configuration dict
        Raises:
            ValueError: for bad configuration values
        """
        for lognm in conf.keys():
            name = lognm.lower()
            if name == 'root':
                log = get_logger()
            else:
                log = get_logger(lognm)
            subconf = conf[lognm]
            if 'output' in subconf:
                dest = subconf['output']
                if dest == '_stdout_':
                    h = logging.StreamHandler(stream=sys.stdout)
                elif dest == '_stderr_':
                    h = logging.StreamHandler(stream=sys.stderr)
                else:
                    try:
                        h = logging.FileHandler(dest)
                    except IOError:
                        raise ValueError('Cannot open output file "{}" for '
                                         'logger "{}"'.format(dest, lognm))
                log.addHandler(h)
            log.setLevel(logging.NOTSET)
            if 'level' in subconf:
                levelnm = subconf['level'].lower()
                level = self._levelnames.get(levelnm, None)
                if level is None:
                    opt = ', '.join(self._levelnames.keys())
                    raise ValueError('Bad level "{}" for logger "{}". Must be '
                                     'one of: {}'.format(levelnm, lognm, opt))
                log.setLevel(level)

    def _validate_conf(self, c):
        if self.CONF_HELP_PATH not in c:
            _log.warn('Path to built HTML documentation is not set. '
                      'The DMF "help" command will not work. To set '
                      'this path, set "{}" in the DMF configuration file.'
                      .format(self.CONF_HELP_PATH))

    @default(CONF_DB_FILE)
    def _default_db_file(self):
        return self.meta.get(self.CONF_DB_FILE, 'resourcedb.json')

    @default(CONF_DATA_DIR)
    def _default_res_dir(self):
        return self.meta.get(self.CONF_DATA_DIR, 'files')

    @observe(CONF_DB_FILE, CONF_DATA_DIR, CONF_HELP_PATH)
    def _observe_setting(self, change):
        if change['type'] != 'change':
            return
        values = {change['name']: change['new']}
        self.set_meta(values)

    def add(self, rsrc):
        """Add a resource and associated files.

        If the resource has 'datafiles', there are some special
        values that cause those files to be copied and possibly the
        original removed at this point. There are attributes `do_copy`
        and `is_tmp` on the resource, and also potentially keys of the
        same name in the datafiles themselves. If present, the datafile
        key/value pairs will override the attributes in the resource.
        For `do_copy`, the original file will be copied into the
        DMF workspace. If `do_copy` is True, then if `is_tmp` is also
        True the original file will be removed (after the copy is made,
        of course).

        Args:
            rsrc (resource.Resource): The resource
        Returns:
            (str) Resource ID
        Raises:
            DMFError, DuplicateResourceError
        """
        # Copy files as necessary
        # Note: this updates paths in the Resource, so should come first
        if 'datafiles' in rsrc.v:
            self._copy_files(rsrc)
        # Add resource
        try:
            self._db.put(rsrc)
        except errors.DuplicateResourceError as err:
            _log.error('Cannot add resource: {}'.format(err))
            raise
        return rsrc.id

    def _copy_files(self, rsrc):
        if 'datafiles_dir' in rsrc.v:
            # If there is a datafiles_dir, use it
            ddir = rsrc.v['datafiles_dir']
        else:
            # If no datafiles_dir, create a random subdir of the DMF
            # configured `_datafile_path`. The subdir prevents name
            # collisions across resources.
            random_subdir = uuid.uuid4().hex
            ddir = os.path.join(self._datafile_path, random_subdir)
        try:
            mkdir_p(ddir)
        except os.error as err:
            raise errors.DMFError('Cannot make dir "{}": {}'
                                  .format(ddir, err))
        for datafile in rsrc.v['datafiles']:
            if 'do_copy' in datafile:
                do_copy = datafile['do_copy']
            else:
                do_copy = rsrc.do_copy
            if do_copy:
                # The `do_copy` flag says do a copy of this datafile from its
                # current path, say /a/path/to/file, into the resource's
                # datafile-dir, say /a/dir/for/resources/, resulting in
                # e.g. /a/dir/for/resources/file.
                filepath = datafile['path']
                filedir, filename = os.path.split(filepath)
                copydir = os.path.join(ddir, filename)
                shutil.copy2(filepath, copydir)
                # The `is_tmp` flag means to remove the original resource file
                # after the copy is done.
                if 'is_tmp' in datafile:
                    is_tmp = datafile['is_tmp']
                else:
                    is_tmp = rsrc.is_tmp
                if is_tmp:
                    try:
                        os.unlink(filepath)
                    except OSError as err:
                        _log.error('Removing temporary datafile "{}": {}'
                                   .format(filepath, err))
                    if 'is_tmp' in datafile:  # remove this directive
                        del datafile['is_tmp']
                datafile['path'] = filename
                datafile['is_copy'] = True
                if 'do_copy' in datafile:  # remove this directive
                    del datafile['do_copy']
            else:
                datafile['is_copy'] = False
            # For idempotence, turn off these flags post-copy
            rsrc.do_copy = rsrc.is_tmp = False
        rsrc.v['datafiles_dir'] = ddir

    def count(self):
        return len(self._db)

    def fetch_one(self, rid):
        """Fetch one resource, from its identifier.

        Args:
            rid (str): Resource identifier
        Returns:
            (resource.Resource) The found resource, or None if no match
        """
        item = self._db.find_one({resource.Resource.ID_FIELD: rid})
        return self._postproc_resource(item)

    def fetch_many(self, rid_list):
        """Fetch multiple resources, by their identifiers.

        Args:
            rid_list (list): List of integer resource identifers

        Returns:
            (list of resource.Resource) List of found resources (may be empty)
        """
        for rid in rid_list:
            yield self.fetch_one(rid)

    def find(self, filter_dict=None, id_only=False):
        """Find and return resources matching the filter.

        The filter syntax is a subset of the MongoDB filter syntax.
        This means that it is represented as a dictionary, where
        each key is an attribute or nested attribute name, and each
        value is the value against which to match. There are five
        possible types of values:

        1. scalar string or number (int, float): Match resources that
           have this exact value for the given attribute.
        2. date, as datetime.datetime or pendulum.Pendulum instance: Match
           resources that have this exact date for the given attribute.
        3. list: Match resources that have a list value for this attribute,
           and for which any of the values in the provided list are in the
           resource's corresponding value. If a '!' is appended to the key
           name, then this will be interpreted as a directive to only match
           resources for which *all* values in the provided list are present.
        4. dict: This is an inequality, with one or more key/value pairs.
           The key is the type of inequality and the value is the numeric
           value for that range. All keys begin with '$'. The possible
           inequalities are:

                - "$lt": Less than (<)
                - "$le": Less than or equal (<=)
                - "$gt": Greater than (>)
                - "$ge": Greater than or equal (>=)
                - "$ne": Not equal to (!=)

        5. Boolean True means does the field exist, and False means
           does it *not* exist.

        Args:
            filter_dict (dict): Search filter.
            id_only (bool): If true, return only the identifier of each
                resource; otherwise a Resource object is returned.

        Returns:
            (list of int|Resource) Depending on the value of `id_only`.
        """
        return (self._postproc_resource(r) for r in
                self._db.find(filter_dict=filter_dict, id_only=id_only))

    def find_related(self, rsrc, filter_dict=None, maxdepth=0, meta=None,
                     outgoing=True):
        """Find related resources.

        Args:
            rsrc (resource.Resource): Resource starting point
            filter_dict (dict): See parameter of same name in :meth:`find`.
            maxdepth (int): Maximum depth of search (starts at 1)
            meta (List[str]): Metadata fields to extract for meta part
            outgoing (bool): If True, look at outgoing relations. Otherwise
                             look at incoming relations. e.g. if A 'uses' B
                             and if True, would find B starting from A.
                             If False, would find A starting from B.
        Returns:
            Generates  triples (depth, Triple, meta), where
            the depth is an integer (starting at 1), the Triple is a
            simple namedtuple wrapping (subject, object, predicate), and
            `meta` is a dict of metadata for the endpoint of the relation
            (the object if outgoing=True, the subject if outgoing=False)
            for the fields provided in the `meta` parameter.
        Raises:
            NoSuchResourceError: if the starting resource is not found
        """
        try:
            return self._db.find_related(rsrc.id, outgoing=outgoing,
                                         maxdepth=maxdepth, meta=meta)
        except KeyError:
            raise errors.NoSuchResourceError(id_=rsrc.id)

    def remove(self, identifier=None, filter_dict=None, update_relations=True):
        """Remove one or more resources, from its identifier or a filter.
        Unless told otherwise, this method will scan the DB and remove
        all relations that involve this resource.

        Args:
            identifier (str): Identifier for a resource.
            filter_dict (dict): Filter to use instead of identifier
            update_relations (bool): If True (the default), scan the DB and
                remove all relations that involve this identifier.
        """
        if not any((identifier, filter_dict)):
            return None
        if identifier:
            id_list = [self.fetch_one(identifier).id]
        else:
            id_list = self.find(filter_dict=filter_dict, id_only=True)
        if not id_list:
            return
        self._db.delete(idlist=id_list)
        # If requested, remove this resource from all the relations where it
        # was a subject or object
        if update_relations:
            for rsrc in self.find():
                keep = [rel for rel in rsrc.v['relations']
                        if rel['identifier'] != identifier]
                # if anything was removed, update the resource
                if len(keep) < len(rsrc.v['relations']):
                    rsrc.v['relations'] = keep
                    # save back to DMF
                    self.update(rsrc)

        if identifier is not None:
            rsrc = self.fetch_one(identifier)
            if rsrc is None:
                _log.error('Cannot find resource id={} to remove'.format(
                    identifier))
                return None
            self._db.delete(identifier)
        else:
            idlist = self.find(filter_dict=filter_dict, id_only=True)
            for eid in idlist:
                self.remove(identifier=eid, update_relations=update_relations)

    def update(self, rsrc, sync_relations=False, upsert=False):
        """Update/insert stored resource.

        Args:
            rsrc (resource.Resource): Resource instance
            sync_relations (bool): If True, and if resource
                exists in the DB, then the "relations" attribute of
                the provided resource will be changed to the stored value.
            upsert (bool): If true, and the resource is not in the DMF,
                           then insert it. If false, and the resource is
                           not in the DMF, then do nothing.
        Returns:
            bool: True if the resource was updated or added, False if nothing
                  was done.
        Raises:
            errors.DMFError: If the input resource was invalid.
        """
        did_update = False
        # sanity-check input
        if not isinstance(rsrc, resource.Resource):
            raise TypeError('Resource type expected, got: {}'.format(
                type(rsrc)))
        # synchronize relations
        if sync_relations:
            db_rsrc = self.fetch_one(rsrc.id)
            if db_rsrc is not None:
                # print('@@ updating relations ({}) to ({})'
                #       .format(rsrc.v['relations'], db_rsrc.v['relations']))
                rsrc.v['relations'] = db_rsrc.v['relations']
        # update or insert new values
        try:
            self._db.update(rsrc.id, rsrc.v)
        except KeyError:
            if upsert:
                self._db.put(rsrc)
                did_update = True
            else:
                raise
        except ValueError as err:
            raise errors.DMFError('Bad value for new resource: {}'
                                  .format(err))
        return did_update

    def _postproc_resource(self, r):
        """Perform any additional changes to resources retrieved
        before passing them up to the application.
        """
        return r

    def __str__(self):
        return 'DMF config="{}"'.format(self._conf)


def get_propertydb_table(rsrc):
    return propdata.PropertyTable.load(rsrc.datafiles[0].fullpath)

