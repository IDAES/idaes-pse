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
DMF support for standard IDAES datasets.

See :mod:`idaes.core.datasets` for user-facing API.
"""
# stdlib
from collections import namedtuple
import json
import logging
from pathlib import Path
from pkg_resources import get_distribution
from typing import Tuple, Dict, Union, List

# package
from idaes.core.dmf import DMF, resource
from idaes.core.dmf.resource import Predicates
from idaes.core.dmf.tables import Table

__authors__ = ["Dan Gunter (LBNL)"]
__author__ = __authors__[0]

# Name of IDAES distribution, same as NAME in setup.py
NAME = "idaes-pse"

# Keep this the same as constant DMF_DATA_ROOT in setup.py
DMF_DATA_ROOT = "data"

_log = logging.getLogger(__name__)


class FileMissingError(Exception):
    pass


class ConfigurationError(Exception):
    pass


def get_dataset_workspace() -> Path:
    dist = get_distribution(NAME)  # XXX check for errors
    return Path(dist.location) / DMF_DATA_ROOT


class Dataset:
    """Base class for datasets."""

    CONF_NAME = "dataset.json"

    def __init__(self, workspace=None):
        if workspace is None:
            workspace = get_dataset_workspace()
        self._dmf = DMF(workspace, create=False)

    def _load_conf(self, directory):
        self._dir = Path(directory)
        if not self._dir.exists():
            raise FileMissingError(f"Dataset directory '{directory}' not found")

        path = self._dir / self.CONF_NAME
        if not path.exists():
            raise FileMissingError(f"Required configuration file '{path}' not found")
        with open(path, "r", encoding="utf-8") as f:
            try:
                conf = json.load(f)
            except json.JSONDecodeError as err:
                raise ConfigurationError(
                    f"Cannot parse configuration file '{path}':" f" {err}"
                )
        self._conf, self._conf_path = conf, path


class PublicationDataset(Dataset):
    """Load a publication and associated data from a directory that includes
    a configuration file (in JSON format, see :class:`Dataset` base class).
    """

    PUBLICATION_KEY = "text"
    NAME_KEY = "name"
    #  Allow picking a tab? SPREADSHEET_TAB_KEY = "tab"
    TABLES_KEY = "tables"

    copy_flag = True

    def __init__(self, workspace=None):
        super(PublicationDataset, self).__init__(workspace=workspace)
        self._conf, self._conf_path = None, None

    def retrieve(self, publication_name) -> Tuple[resource.Resource, Dict[str, Table]]:
        pub = self._dmf.find_one(name=publication_name)
        if not pub:
            raise KeyError(f"No publication found for name '{publication_name}'")
        tables = {}
        for r in self._dmf.find_related_resources(pub, predicate=Predicates.derived):
            tables[r.name] = r.table
        if not tables:
            _log.warning(f"No tables found for publication '{publication_name}'")
        return pub, tables

    def load(self, directory):
        """Load publication from the provided directory."""
        directory = Path(directory)
        self._load_conf(directory)

        # Check for required keys
        for key in self.PUBLICATION_KEY, self.NAME_KEY:
            if key not in self._conf:
                raise ConfigurationError(
                    f"Required key '{key}' missing "
                    f"from configuration file '{self._conf_path}'"
                )

        # Create publication resource
        _log.debug("publication_resource.create.begin")
        name = self._conf[self.NAME_KEY]
        pub_r = resource.Resource(type_=resource.ResourceTypes.publication, name=name)
        # Add publication text and associated metadata
        pub = self._conf[self.PUBLICATION_KEY]
        try:
            file_ = pub["file"]
        except KeyError:
            raise ConfigurationError(
                f"Required key 'file' in section '{self.PUBLICATION_KEY}' "
                f"missing from configuration file '{self._conf_path}'"
            )
        meta = {
            "date": pub.get("date", ""),
            "doi": pub.get("doi", ""),
            "isbn": pub.get("isbn", ""),
            "language": pub.get("language", "english"),
        }
        if "isbn" in pub:
            meta["source"] = f"{pub.get('authors', 'Anon.')}, "
            f"\"{pub.get('title', 'no title')}\". "
            f"{pub.get('publisher', '')} ({pub.get('date', '')})"
        else:
            meta["source"] = (
                f"{pub.get('authors', 'Anon.')}, "
                f"\"{pub.get('title', 'no title')}\". "
                f"{pub.get('venue', '')} ({pub.get('date', '')})"
            )
        pub_r.add_data_file(directory / file_, do_copy=self.copy_flag)
        pub_r.sources.append(meta)
        pub_r.desc = pub.get("title", name)
        self._dmf.add(pub_r)
        _log.debug(f"publication_resource.create.end id={pub_r.id}")

        # Create and add a resource for each table
        tables = self._conf.get(self.TABLES_KEY, None)
        if tables:
            num = 0
            for table in tables:
                _log.debug("publication_table_resource.create.begin")
                tbl_r = self._create_table_resource(table, directory, name, pub_r)
                num += 1
                table_id = str(tbl_r.id) if tbl_r else "<no table>"
                _log.debug(
                    f"publication_table_resource.create.end " f"id={table_id} num={num}"
                )

        # Update all the derivation relationships
        _log.debug("dmf_update.begin")
        self._dmf.update()
        _log.debug("dmf_update.end")

    _default_table_name = "unknown"  # Default name for a table

    def _create_table_resource(self, table, directory, pub_name, pub_resource):
        # get table name
        table_name = table.get("name", self._default_table_name)
        # set up datafile
        if "datafile" not in table:
            _log.warning(f"Skipping table '{table_name}': no datafile")
            return None
        table_datafile = table["datafile"]
        table_path = directory / table_datafile
        if not table_path.exists:
            raise ConfigurationError(
                f"Cannot find data file for table '{table_name}': {table_path}"
            )
        # create new resource
        tbl_r = resource.Resource(type_=resource.ResourceTypes.tabular)
        # populate resource
        table_desc = table.get("description", None)
        if table_desc is None:
            _log.warning(f"No description given for table data " f"{table_datafile}")
        _log.debug(f"Adding table: path={table_path} desc={table_desc}")
        tbl_r.add_table(table_path, desc=table_desc, do_copy=self.copy_flag)
        tbl_r.desc = table_desc
        tbl_r.name = table_name  # set name for this table
        tbl_r.add_tag(pub_name)  # add tag shared with publication
        self._dmf.add(tbl_r)
        resource.create_relation(pub_resource, resource.Predicates.derived, tbl_r)
        return tbl_r


AvailableResult = namedtuple("AvailableResult", "Class description")


class Publication:
    """Abstract superclass for the public interface to a publication-derived dataset.

    Do not instantiate directly. Instead, subclass and pass the appropriate
    name of the dataset (i.e., its name in the DMF) to this constructor.
    See :class:`idaes.core.datasets.Pitzer` for an example.
    """

    def __init__(self, name, workspace=None):
        self._ds = PublicationDataset(workspace=workspace)
        self._name = name
        self._pub, self._tables = self._ds.retrieve(name)

    def list_tables(self) -> List[str]:
        return list(self._tables.keys())

    def get_table(self, name) -> Union[Table, None]:
        return self._tables.get(name, None)
