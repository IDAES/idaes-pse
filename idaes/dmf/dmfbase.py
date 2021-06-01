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
Data Management Framework
"""
# stdlib
from datetime import datetime
from io import IOBase
import logging
import os
import pathlib
import re
import shutil
import sys
import uuid
from typing import Generator, Union

# third-party
from traitlets import HasTraits, default, observe
from traitlets import Unicode
import yaml

# local
from . import errors
from .resource import Resource
from . import resourcedb
from . import workspace
from .util import mkdir_p, yaml_load


__author__ = "Dan Gunter"

_log = logging.getLogger(__name__)


class DMFConfig(object):
    """Global DMF configuration.

    Every time you create an instance of the :class:`DMF`
    or run a ``dmf`` command on the command-line, the library opens the
    global DMF configuration file to figure out the default workspace
    (and, eventually, other values).

    The default location for this configuration file is "~/.dmf", i.e.
    the file named ".dmf" in the user's home directory. This can be
    modified programmatically by changing the "filename" attribute of
    this class.

    The contents of the configuration are formatted as `YAML`_ with
    the following keys defined:

        workspace
            Path to the default workspace directory.

    .. _YAML: http://www.yaml.org/
    """

    # configuration file location
    _filename = os.path.expanduser("~/.dmf")
    # known keys
    WORKSPACE = "workspace"
    _keys = [WORKSPACE]
    # defaults
    DEFAULTS = {WORKSPACE: os.path.abspath(".")}

    def __init__(self, defaults=None):
        """Create new configuration.

        Args:
            defaults (dict): Default values to use if nothing is found.
                            If not provided, use self.DEFAULTS
        Raises:
            IOError: If self.filename is not openable or readable
            ValueError: If the file can't be parsed as YAML
        """
        # init to defaults
        self.c = {}
        if defaults:
            self.c.update(defaults)
        else:
            self.c.update(self.DEFAULTS)
        # try to open the config file
        fp = None
        try:
            if not os.path.exists(self._filename):
                raise IOError("File not found: {}".format(self._filename))
            fp = open(self._filename, "rb")
        except IOError as err:
            _log.warning(
                "Unable to open global DMF configuration file "
                "for reading: {}. Using default configuration values.".format(err)
            )
        # if we got a config file, parse it
        if fp:
            try:
                self._parse(fp)
            except ValueError as err:
                raise ValueError(
                    'Parsing configuration at "{}": {}'.format(self._filename, err)
                )

    @classmethod
    def configuration_exists(cls):
        return pathlib.Path(cls._filename).exists()

    @classmethod
    def configuration_path(cls) -> pathlib.Path:
        return pathlib.Path(cls._filename)

    def save(self):
        _log.info(f"Saving configuration location to: {self._filename}")
        try:
            fp = open(self._filename, "w")
        except IOError as err:
            _log.info(f"Saving configuration to '{self._filename}' failed: {err}")
            raise
        _log.debug(f"Writing new global config to: {self._filename}")
        yaml.dump(self.c, fp)

    @property
    def workspace(self):
        return self.c[self.WORKSPACE]

    def _parse(self, fp):
        try:
            y = yaml_load(fp)
        except Exception as err:
            raise ValueError(err)
        if y:
            for k, v in y.items():
                if k in self._keys:
                    self.c[k] = v


class DMF(workspace.Workspace, HasTraits):
    """Data Management Framework (DMF).

    Expected usage is to instantiate this class, once, and then use it for
    storing, searching, and retrieving resources that
    are required for the given analysis.

    For details on the configuration files used by the DMF, see
    documentation for :class:`DMFConfig` (global configuration) and
    :class:`idaes.dmf.workspace.Workspace`.
    """

    db_file = Unicode(help="Database file name")
    datafile_dir = Unicode(help="Data file directory, " "relative to DMF root")

    CONF_DB_FILE = "db_file"
    CONF_DATA_DIR = "datafile_dir"
    CONF_HELP_PATH = workspace.Fields.DOC_HTML_PATH

    # logging should really provide this
    _levelnames = {
        "fatal": logging.FATAL,
        "error": logging.ERROR,
        "warn": logging.WARN,
        "warning": logging.WARN,
        "info": logging.INFO,
        "debug": logging.DEBUG,
    }

    def __init__(
        self, path: Union[pathlib.Path, str] = "", name=None, desc=None, create=False,
            existing_ok=True, save_path=False, **ws_kwargs
    ):
        """Create or load DMF workspace.

        Args:
            path (str|Path): Path to workspace. If given, will override any
                  global configuration. If not given (empty or None), then
                  global configuration will be checked first. If empty, and
                  the global configuration does not parse or exist, then
                  "." will be used.
            name (str): Name to be used for workspace.
            desc (str): Longer description of workspace.
            create (bool): If the path to the workspace does not exist,
                           this controls whether to create it or raise
                           an error.
            existing_ok (bool): See :class:`workspace.Workspace` constructor
            save_path: If True, save provided path globally
            **ws_kwargs: Keyword arguments for :meth:`workspace.Workspace()`
                         constructor.
        Raises:
            WorkspaceError: Any error that could be raised by the Workspace
                            superclass.
        """
        # get global configuration
        conf = DMFConfig(defaults={})
        # normalize path to a string
        path = str(path)
        # get path, if not specified, from configuration
        if not path:
            path = conf.c.get(DMFConfig.WORKSPACE, ".")
        elif save_path:
            conf.c[DMFConfig.WORKSPACE] = os.path.abspath(path)
        # set up workspace
        ws_kwargs["create"] = create
        ws_kwargs["existing_ok"] = existing_ok
        try:
            super(DMF, self).__init__(path, **ws_kwargs)
        except (errors.ParseError, ValueError) as err:
            msg = 'Configuration "{}", parse error: {}'.format(path, err)
            raise errors.WorkspaceError(msg)
        # check local config
        self._validate_conf(self.meta)
        # merge selected global config values into local
        # XXX: not done
        # set up logging
        if workspace.Fields.LOG_CONF in self.meta:
            try:
                self._configure_logging()
            except ValueError as err:
                msg = 'Configuration "{}", logging section error: {}'.format(path, err)
                raise errors.WorkspaceError(msg)
        # set up rest of DMF
        path = os.path.join(self.root, self.db_file)
        self._db = resourcedb.ResourceDB(path)
        self._datafile_path = os.path.join(self.root, self.datafile_dir)
        if not os.path.exists(self._datafile_path):
            os.mkdir(self._datafile_path, 0o750)
        # add create/modified date, and optional name/description
        _w = workspace.Workspace
        right_now = datetime.isoformat(datetime.now())
        meta = {_w.CONF_CREATED: right_now, _w.CONF_MODIFIED: right_now}
        if name:
            meta[_w.CONF_NAME] = name
        if desc:
            meta[_w.CONF_DESC] = desc
        self.set_meta(meta)
        # save global config mods
        conf.save()
        # initialize resources for this instance
        self._resources = {}

    def _configure_logging(self):
        """Configure logging for DMF.

        Expected schema::

            <logger-name>:
                level: <levelname>              default=NOTSET
                output: file|_stdout_|_stderr_  default=_stderr_

        Raises:
            ValueError: for bad configuration values
        """
        conf = self.meta[workspace.Fields.LOG_CONF]
        for lognm in conf.keys():
            name = lognm.lower()
            if name == "root":
                log = self._get_logger()
            else:
                log = self._get_logger(lognm)
            subconf = conf[lognm]
            if "output" in subconf:
                dest = subconf["output"]
                if dest == "_stdout_":
                    h = logging.StreamHandler(stream=sys.stdout)
                elif dest == "_stderr_":
                    h = logging.StreamHandler(stream=sys.stderr)
                else:
                    try:
                        h = logging.FileHandler(dest)
                    except IOError:
                        raise ValueError(
                            'Cannot open output file "{}" for '
                            'logger "{}"'.format(dest, lognm)
                        )
                log.addHandler(h)
            log.setLevel(logging.NOTSET)
            if "level" in subconf:
                levelnm = subconf["level"].lower()
                level = self._levelnames.get(levelnm, None)
                if level is None:
                    opt = ", ".join(self._levelnames.keys())
                    raise ValueError(
                        'Bad level "{}" for logger "{}". Must be '
                        "one of: {}".format(levelnm, lognm, opt)
                    )
                log.setLevel(level)
            if "format" in subconf:
                fmt = logging.Formatter(subconf["format"])
                if log.hasHandlers():
                    for h in log.handlers:
                        h.setFormatter(fmt)

    @staticmethod
    def _get_logger(name=None):
        """Get a logger by absolute or short-hand name.
        """
        if name is None:
            fullname = "idaes"
        # idaes.<whatever> just use name as-is
        elif re.match(r"idaes\.[a-zA-Z_.]+", name):
            fullname = name
        # .<whatever>, root at idaes
        elif re.match(r"\.[a-zA-Z_.]+", name):
            fullname = "idaes" + name
        # some special names
        elif name in ("dmf", "core", "models", "ui"):
            fullname = "idaes." + name
        # when in doubt just take value as-is
        else:
            fullname = name
        return logging.getLogger(fullname)

    def _validate_conf(self, c):
        if self.CONF_HELP_PATH not in c:
            # _log.info(
            #     "Path to built HTML documentation is not set. "
            #     'The DMF "help" command will not work. To set '
            #     'this path, set "{}" in the DMF configuration file.'.format(
            #         self.CONF_HELP_PATH
            #     )
            #)
            pass

    @default(CONF_DB_FILE)
    def _default_db_file(self):
        return self.meta.get(self.CONF_DB_FILE, "resourcedb.json")

    @default(CONF_DATA_DIR)
    def _default_res_dir(self):
        return self.meta.get(self.CONF_DATA_DIR, "files")

    @observe(CONF_DB_FILE, CONF_DATA_DIR, CONF_HELP_PATH)
    def _observe_setting(self, change):
        if change["type"] != "change":
            return
        values = {change["name"]: change["new"]}
        self.set_meta(values)

    @property
    def datafiles_path(self):
        return self._datafile_path

    @property
    def workspace_path(self) -> pathlib.Path:
        """Path to workspace directory.
        """
        return pathlib.Path(self.root)

    def new(
        self,
        file: Union[str, pathlib.Path] = None,
        file_kw: dict = None,
        resource=None,
        type_: str = None,
        values: dict = None,
        **kwargs,
    ) -> Resource:
        """Creates a resource and adds it to the DMF. There are a number of ways to specify how to
        create this resource, detailed below.

        Args:
            file (optional): File from which to construct the resource
            file_kw (optional): Additional keywords for `Resource.from_file`
            resource (optional): Existing resource to copy/add into this DMF instance
            type_ (optional): When constructing from kwargs, or value, this specifies the type of
                resource to create, just like the same-named argument to the Resource constructor.
            values (optional): An alternative to kwargs, useful for keys that are not valid
                Python identifiers.
            kwargs: Construct resource from key/value pairs, or add these pairs to a resource
                 constructed with other arguments.

        Returns:
             The resource object

        Raises:
            DMFError, DuplicateResourceError: Trouble building resource
            ValueError: invalid argument combination
        """
        constructed_from_kwargs = False
        if file is not None:
            if file_kw is None:
                file_kw = {}
            rsrc = Resource.from_file(str(file), **file_kw)
        elif resource is not None:
            rsrc = resource
        elif values is not None:
            rsrc = Resource(value=values, type_=type_)
        elif kwargs is not None:
            rsrc = Resource(value=kwargs, type_=type_)
            constructed_from_kwargs = True
        else:
            raise ValueError("No arguments given to constructor")
        # add kwargs, optionally
        if (not constructed_from_kwargs) and kwargs:
            rsrc.set_values(kwargs)
        # add to DMF
        rsrc_id = self.add(rsrc)
        return rsrc

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

        Resources added during the lifetime of this DMF instance are remembered,
        so that `update()` with no arguments applies to all of them.

        Args:
            rsrc (Resource): The resource
        Returns:
            (str) Resource ID
        Raises:
            DMFError, DuplicateResourceError
        """
        # Copy files as necessary
        # Note: this updates paths in the Resource, so should come first
        if "datafiles" in rsrc.v:
            self._copy_files(rsrc)
        # Add resource
        try:
            self._db.put(rsrc)
        except errors.DuplicateResourceError as err:
            _log.error("Cannot add resource: {}".format(err))
            raise
        # if that worked, remember in session store
        self._resources[rsrc.id] = rsrc
        return rsrc.id

    def attach(self, resource: Resource) -> Resource:
        """Associate this resource with this instance of the DMF, so that e.g. `update()`
        with no arguments will save updated values.

        Args:
            resource: The resource to attach

        Returns:
            input value (for chaining)
        """
        self._resources[resource.id] = resource
        return resource

    def _copy_files(self, rsrc):
        if rsrc.v.get("datafiles_dir", None):
            # If there is a datafiles_dir, use it
            ddir = rsrc.v["datafiles_dir"]
            _log.debug(f"_copy_files: use existing datafiles dir '{ddir}'")
        else:
            # If no datafiles_dir, create a random subdir of the DMF
            # configured `_datafile_path`. The subdir prevents name
            # collisions across resources.
            random_subdir = uuid.uuid4().hex
            ddir = os.path.join(self._datafile_path, random_subdir)
            _log.debug(f"_copy_files: create new datafiles dir '{ddir}'")
        try:
            mkdir_p(ddir)
        except os.error as err:
            raise errors.DMFError('Cannot make dir "{}": {}'.format(ddir, err))
        for datafile in rsrc.v["datafiles"]:
            if "do_copy" in datafile:
                do_copy = datafile["do_copy"]
            else:
                do_copy = rsrc.do_copy
            if do_copy:
                # The `do_copy` flag says do a copy of this datafile from its
                # current path, say /a/path/to/file, into the resource's
                # datafile-dir, say /a/dir/for/resources/, resulting in
                # e.g. /a/dir/for/resources/file.
                filepath = datafile["path"]
                _, filename = os.path.split(filepath)
                copydir = os.path.join(ddir, filename)
                _log.debug(
                    'Copying datafile "{}" to directory "{}"'.format(filepath, copydir)
                )
                try:
                    shutil.copy2(filepath, copydir)
                except (IOError, OSError) as err:
                    msg = (
                        'Cannot copy datafile from "{}" to DMF '
                        'directory "{}": {}'.format(filepath, copydir, err)
                    )
                    _log.error(msg)
                    raise errors.DMFError(msg)
                # The `is_tmp` flag means to remove the original resource file
                # after the copy is done.
                if "is_tmp" in datafile:
                    is_tmp = datafile["is_tmp"]
                else:
                    is_tmp = rsrc.is_tmp
                if is_tmp:
                    _log.debug(
                        "Temporary datafile flag is on, removing "
                        'original datafile "{}"'.format(filepath)
                    )
                    try:
                        os.unlink(filepath)
                    except OSError as err:
                        _log.error(
                            'Removing temporary datafile "{}": {}'.format(filepath, err)
                        )
                    if "is_tmp" in datafile:  # remove this directive
                        del datafile["is_tmp"]
                datafile["path"] = filename
                datafile["is_copy"] = True
                if "do_copy" in datafile:  # remove this directive
                    del datafile["do_copy"]
            else:
                datafile["is_copy"] = False
        # For idempotence, turn off these flags post-copy
        rsrc.do_copy = rsrc.is_tmp = False
        # Make sure datafiles dir is in sync
        rsrc.v["datafiles_dir"] = ddir

    def count(self):
        return len(self._db)

    def fetch_one(self, rid, id_only=False):
        """Fetch one resource, from its identifier.

        Args:
            rid (str): Resource identifier
            id_only (bool): If true, return only the identifier of each
                resource; otherwise a Resource object is returned.
        Returns:
            (Resource) The found resource, or None if no match
        """
        item = self._db.find_one({Resource.ID_FIELD: rid}, id_only=id_only)
        if item is None:
            return None
        elif id_only:
            return item
        else:
            return self._postproc_resource(item)

    def find(self, filter_dict=None, name=None, id_only=False, re_flags=0):
        """Find and return resources matching the filter.

        The filter syntax is a subset of the MongoDB filter syntax.
        This means that it is represented as a dictionary, where
        each key is an attribute or nested attribute name, and each
        value is the value against which to match. There are six
        possible types of values:

        1. scalar string or number (int, float): Match resources that
           have this exact value for the given attribute.
        2. special scalars "@<value>":

                - "@true"/"@false": boolean (bare True/False will test existence)

        3. date, as datetime.datetime  instance: Match
           resources that have this exact date for the given attribute.
        4. list: Match resources that have a list value for this attribute,
           and for which any of the values in the provided list are in the
           resource's corresponding value. If a '!' is appended to the key
           name, then this will be interpreted as a directive to only match
           resources for which *all* values in the provided list are present.
        5. dict: This is an inequality, with one or more key/value pairs.
           The key is the type of inequality and the value is the numeric
           value for that range. All keys begin with '$'. The possible
           inequalities are:

                - "$lt": Less than (<)
                - "$le": Less than or equal (<=)
                - "$gt": Greater than (>)
                - "$ge": Greater than or equal (>=)
                - "$ne": Not equal to (!=)

        6. Boolean True means does the field exist, and False means
           does it *not* exist.
        7. Regular expression, string "~<expr>" and `re_flags`
           for flags (understood: re.IGNORECASE)

        Args:
            filter_dict (dict): Search filter.
            name (str): If present, add {'aliases': [<name>]} to filter_dict. This
                        is syntactic sugar for a common case.
            id_only (bool): If true, return only the identifier of each
                resource; otherwise a Resource object is returned.
            re_flags (int): Flags for regex filters

        Returns:
            (list of int|Resource) Depending on the value of `id_only`.
        """
        if name is not None:
            if filter_dict is None:
                filter_dict = {"aliases": [name]}
            elif "aliases" in filter_dict:
                filter_dict["aliases"].append(name)
            else:
                filter_dict["aliases"] = [name]
        elif filter_dict is None:
            filter_dict = {}
        elif not hasattr(filter_dict, "items"):
            raise TypeError(
                "Parameter 'filter_dict' must be a dictionary, not a "
                f"{type(filter_dict)}"
            )
        return (
            self._postproc_resource(r)
            for r in self._db.find(
                filter_dict=filter_dict, id_only=id_only, flags=re_flags
            )
        )

    def find_one(self, *args, **kwargs) -> Resource:
        """Find and return one matching resource.

        This is a wrapper around `find()` that does two things:

            1. Only returns the first result (if any)

            2. By default, calls `attach()` on the result to add it to the DMF.
               You can change this by passing `attach=False` as a keyword.

        Returns:
            Resource, or None or an empty list.
        """
        attach_to_dmf = True
        if 'attach' in kwargs:
            attach_to_dmf = kwargs['attach']
            del kwargs['attach']
        results = self.find(*args, **kwargs)
        if results is None:
            return None
        result_list = list(results)
        if len(result_list) == 0:
            return []
        result = result_list[0]
        if attach_to_dmf:
            self.attach(result)
        return result

    def find_by_id(self, identifier: str, id_only=False) -> Generator:
        """Find resources by their identifier or identifier prefix.
        """
        if len(identifier) == Resource.ID_LENGTH:
            for rsrc in self._db.find(
                {Resource.ID_FIELD: identifier}, id_only=id_only
            ):
                yield rsrc
        else:
            regex, flags = f"{identifier}[a-z]*", re.IGNORECASE
            for rsrc in self._db.find(
                {Resource.ID_FIELD: "~" + regex}, flags=flags, id_only=id_only
            ):
                yield rsrc

    def find_related(
        self, rsrc, filter_dict=None, maxdepth=0, meta=None, outgoing=True
    ):
        """Find related resources.

        Args:
            rsrc (Resource): Resource starting point
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
        if meta is None:
            meta = [
                Resource.ID_FIELD,
                Resource.TYPE_FIELD,
                "desc",
                "version_info",
            ]
        else:
            if Resource.ID_FIELD not in meta:
                meta.insert(0, Resource.ID_FIELD)
        try:
            return self._db.find_related(
                rsrc.id,
                outgoing=outgoing,
                maxdepth=maxdepth,
                meta=meta,
                filter_dict=filter_dict,
            )
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
        id_list, rid_list = None, None
        if identifier:
            # sanity check identifier type
            if not hasattr(identifier, "lower"):
                raise TypeError(
                    f"identifier argument is not a string. type={type(identifier)}"
                )
            identifier = str(identifier)
            id_one = self.fetch_one(identifier, id_only=True)
            id_list = None if id_one is None else [id_one]
            rid_list = [identifier]
        else:
            id_list = list(self.find(filter_dict=filter_dict, id_only=True))
            if id_list:
                rid_list = []
                for i in id_list:
                    rsrc = self._db.get(i)
                    rid_list.append(rsrc.id)
        if not id_list:
            _log.info(
                "Cannot remove resource-id={} filter={}: Not found".format(
                    identifier, filter_dict
                )
            )
            return
        # TODO: Delete all associated data files!!
        self._db.delete(idlist=id_list, internal_ids=True)
        # delete any added during this session
        for rsrc_id in id_list:
            if rsrc_id in self._resources:
                del self._resources[rsrc_id]
        # If requested, remove deleted resources from all the relations
        # where it was a subject or object
        if update_relations:
            # look at all resources in DB
            for rsrc in self.find():
                # for each one figure out which relations to keep
                keep = []
                for rel in rsrc.v["relations"]:
                    # if none of the removed resource ids are present, keep it
                    if rel["identifier"] not in rid_list:
                        keep.append(rel)
                # if we didn't keep all the relations, update the resource
                if len(keep) < len(rsrc.v["relations"]):
                    rsrc.v["relations"] = keep
                    # save back to DMF
                    self.update(rsrc)

    def update(
        self,
        rsrc: Resource = None,
        sync_relations: bool = False,
        upsert: bool = False,
    ):
        """Update/insert stored resource.

        Args:
            rsrc (optional): Resource instance
            sync_relations (optional): If True, and if resource exists in the DB, then the "relations" attribute of
                the provided resource will be changed to the stored value.
            upsert (optional): If true, and the resource is not in the DMF, then insert it. If false, and the
                resource is not in the DMF, then do nothing.
        Returns:
            bool, int: For a single resource: True if the resource was updated or added, False if nothing was done.
                If no specific resource was given, returns the number of resources updated. You can compare this
                with `resource_count()`.
        Raises:
            errors.DMFError: If the input resource was invalid.
        """
        did_update = False
        # with no specific resource, update all resources added to this instance
        if rsrc is None:
            return self._update_resources(sync_relations)
        # sanity-check input
        if not isinstance(rsrc, Resource):
            raise TypeError("Resource type expected, got: {}".format(type(rsrc)))
        # synchronize relations
        if sync_relations:
            _log.debug("synchronize relations")
            db_rsrc = self.fetch_one(rsrc.id)
            if db_rsrc is not None:
                rsrc.v["relations"] = db_rsrc.v["relations"]
        # update or insert new values
        try:
            self._db.update(rsrc.id, rsrc.v)
            did_update = True
        except errors.NoSuchResourceError:
            if upsert:
                self._db.put(rsrc)
                did_update = True
            else:
                raise
        except ValueError as err:
            raise errors.DMFError("Bad value for new resource: {}".format(err))
        return did_update

    def _update_resources(self, sync_relations) -> int:
        valid_resources = {}
        for key, rsrc in self._resources.items():
            try:
                self.update(rsrc=rsrc, sync_relations=sync_relations, upsert=False)
            except errors.NoSuchResourceError as err:
                _log.warning(f"During update, resource not found: {err}")
            except Exception as err:
                _log.warning(f"During update, unknown error: {err}")
            else:
                # if it was OK, add to valid resources dict
                valid_resources[key] = rsrc
        # keep only valid resources
        self._resources = valid_resources
        return len(valid_resources)

    def _postproc_resource(self, r):
        """Perform any additional changes to resources retrieved
        before passing them up to the application.
        """
        return r

    def __str__(self):
        return 'DMF config="{}"'.format(self._conf)

    @property
    def resource_count(self) -> int:
        """How many resources have been added to this instance of the DMF.
        """
        return len(self._resources)


def get_propertydb_table(rsrc):
    from idaes.dmf import propdata

    return propdata.PropertyTable.load(rsrc.datafiles[0].fullpath)


# 1. scalar string or number (int, float): Match resources that
#            have this exact value for the given attribute.
#         2. special scalars "@<value>":
#
#                 - "@true"/"@false": boolean (bare True/False will test existence)
#
#         3. date, as datetime.datetime instance: Match
#            resources that have this exact date for the given attribute.
#         4. list: Match resources that have a list value for this attribute,
#            and for which any of the values in the provided list are in the
#            resource's corresponding value. If a '!' is appended to the key
#            name, then this will be interpreted as a directive to only match
#            resources for which *all* values in the provided list are present.
#         5. dict: This is an inequality, with one or more key/value pairs.
#            The key is the type of inequality and the value is the numeric
#            value for that range. All keys begin with '$'. The possible
#            inequalities are:
#
#                 - "$lt": Less than (<)
#                 - "$le": Less than or equal (<=)
#                 - "$gt": Greater than (>)
#                 - "$ge": Greater than or equal (>=)
#                 - "$ne": Not equal to (!=)
#
#         6. Boolean True means does the field exist, and False means
#            does it *not* exist.
#         7. Regular expression, string "~<expr>" and `re_flags`
#            for flags (understood: re.IGNORECASE)
