##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
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
Storage of models for the IDAES Flowsheet Visualizer.

Currently implemented methods are trivial storage in memory and storage to a file.
"""
# stdlib
from abc import ABC, abstractmethod
import json
from pathlib import Path
from typing import Dict, Union
# package
from idaes import logger

_log = logger.getLogger(__name__)


class DataStore(ABC):

    @abstractmethod
    def save(self, data: Union[Dict, str]):
        pass

    @abstractmethod
    def load(self) -> Union[Dict, str]:
        pass

    @classmethod
    def create(cls, dest=None) -> 'DataStore':
        """Factory method to create and return the appropriate DataStore subclass
        given a destination.

        Args:
            dest: If a string or Path, return a FileDataStore

        Raises:
            ValueError if `dest` can't be matched to a DataStore subclass.
        """
        if isinstance(dest, str):
            return FileDataStore(Path(dest))
        elif isinstance(dest, Path):
            return FileDataStore(dest)
        elif dest is None:
            return MemoryDataStore()
        else:
            raise ValueError(f"Unknown destination '{dest}' for type '{type(dest)}'")

    def __str__(self):
        return "generic storage; should not be used directly"


class FileDataStore(DataStore):
    def __init__(self, path):
        self._p = path
        self._is_json = False

    def save(self, data):
        with self._p.open("w") as fp:
            if isinstance(data, dict):
                json.dump(data, fp)
                self._is_json = True
            else:
                fp.write(data)
                self._is_json = False

    def load(self):
        with self._p.open("r") as fp:
            if self._is_json:
                data = json.load(fp)
            else:
                data = fp.read()
        return data

    def __str__(self):
        return f"file storage at '{self._p}'"


class MemoryDataStore(DataStore):
    def __init__(self):
        self._data = None

    def save(self, data):
        self._data = data

    def load(self):
        return self._data

    def __str__(self):
        return "memory storage"


class DataStoreManager:
    """Manage operations on multiple id/data-store pairs.
    """
    def __init__(self):
        self._id_store = {}
        self._id_path = {}  # maps identifiers to Path objects

    def add(self, id_: str, store: DataStore):
        """Add an identifier and associated storage location.

        Args:
            id_: Identifier
            store: Where to store it
        """
        self._id_store[id_] = store

    def _check_id(self, id_: str):
        if id_ not in self._id_store:
            raise KeyError(f"No storage associated with identifier '{id_}'")

    def save(self, id_: str, data: Union[Dict, str]):
        self._check_id(id_)
        self._id_store[id_].save(data)

    def load(self, id_: str) -> Union[Dict, str]:
        self._check_id(id_)
        return self._id_store[id_].load()

