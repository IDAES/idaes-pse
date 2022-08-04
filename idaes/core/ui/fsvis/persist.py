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
from . import errors

_log = logger.getLogger(__name__)


class DataStore(ABC):
    @abstractmethod
    def save(self, data: Union[Dict, str]):
        """Save data.

        Args:
            data: Data to save.
        """
        pass

    @abstractmethod
    def load(self) -> Dict:
        """Load data.

        Returns:
            data, as a dict (even if string of JSON was passed as input to `save()`)
        """
        pass

    @classmethod
    def create(cls, dest=None) -> "DataStore":
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
            print("create memory store")
            return MemoryDataStore()
        else:
            raise ValueError(f"Unknown destination '{dest}' for type '{type(dest)}'")

    def __str__(self):
        """Return string version of self.

        Used for equality tests, so make sure this is unique as needed for equality tests
        between stores to behave correctly.
        """
        return "generic storage; should not be used directly"

    def __eq__(self, other):
        return str(other) == str(self)

    @property
    def filename(self):
        return ""


class FileDataStore(DataStore):
    def __init__(self, path: Path):
        self._p = path

    @property
    def path(self) -> Path:
        return Path(str(self._p))  # return a copy

    def save(self, data):
        """Save data to a file.

        Args:
            data: If a dict, must be serializable by `json.dump()`, otherwise
                  treated as a string of JSON (e.g. '{"name": "value"}')

        Returns:
            None

        Raises:
            DataStoreError, if serialization or I/O fails
        """
        _log.debug(f"Save to file: {self._p}")
        try:
            with self._p.open("w") as fp:
                if isinstance(data, dict):
                    try:
                        json.dump(data, fp)
                    except TypeError as err:
                        raise errors.DatastoreSerializeError(data, err, stream=fp)
                else:
                    _parse_json(data)  # validation
                    fp.write(str(data))
        except ValueError as err:
            raise errors.DatastoreError(str(err))
        except IOError as err:
            raise errors.DatastoreSaveError(f"IO error with datastore: {err}")

    def load(self):
        _log.debug(f"Load from file: {self._p}")
        try:
            with self._p.open("r") as fp:
                try:
                    data = json.load(fp)
                except json.JSONDecodeError as err:
                    raise ValueError(f"Reading JSON failed: {err}")
        except FileNotFoundError:
            # normalize errors finding stored object to ValueError
            raise ValueError(f"File '{self._p}' not found")
        return data

    def __str__(self):
        return f"file '{self._p}'"

    @property
    def filename(self):
        return str(self._p)


class MemoryDataStore(DataStore):
    def __init__(self):
        self._data = None

    def save(self, data: Union[str, Dict]):
        """Store data in memory.

        Args:
            data: If a dict, must be serializable by `json.dump()`, otherwise
                  treated as a string of JSON (e.g. '{"name": "value"}').
                  This is enforced even though, obviously, the data can be "saved"
                  in memory even if it is *not* JSON-serializable

        Returns:
            None

        Raises:
            DataStoreError, if serialization fails
        """
        if isinstance(data, dict):
            try:
                json.dumps(data)
            except TypeError as err:
                raise errors.DatastoreSerializeError(data, err)
            self._data = data
        else:
            self._data = _parse_json(data)

    def load(self) -> Dict:
        if self._data is None:
            raise ValueError("Data is empty")
        return self._data

    def __str__(self):
        return "__MEMORY__"


def _parse_json(data) -> Dict:
    """Parse string of the data to JSON, with desired exceptions raised.

    Raises:
        errors.DatastoreSerializeError: If the parse fails
    """
    s = str(data)
    try:
        data = json.loads(s)
    except json.decoder.JSONDecodeError as err:
        err2 = f"could not decode string of JSON: {err}"
        raise errors.DatastoreSerializeError(s[:256], err2)
    return data


class DataStoreManager:
    """Manage operations on multiple id/data-store pairs."""

    def __init__(self):
        self._id_store = {}
        self._id_path = {}  # maps identifiers to Path objects

    def add(self, id_: str, store: DataStore) -> bool:
        """Add an identifier and associated storage location.

        If the same datastore is already there, do nothing.

        Args:
            id_: Identifier
            store: Where to store it
        Returns:
            True if we added a new store, False if we did nothing
        """
        if id_ in self._id_store and self._id_store[id_] == store:
            _log.debug(f"Use existing store, {self._id_store[id_]}, for '{id_}'")
            added = False
        else:
            self._id_store[id_] = store
            added = True
        return added

    def save(self, id_: str, data: Union[Dict, str]):
        """Save flowsheet with given identifier.

        Args:
            id_: Flowsheet identifier
            data: Data for the flowsheet, either as serialized JSON or a Python dict

        Returns:
            None

        Raises:
            KeyError if the flowsheet is not found
            DatastoreSerializeError on JSON errors
        """
        self._find(id_).save(data)
        _log.debug(f"Flowsheet '{id_}' saved")

    def load(self, id_: str) -> Dict:
        """Load a flowhseet with a given identifier.

        Args:
            id_: Flowsheet identifier

        Returns:
            Flowsheet (always as a dict, no matter how it was saved)

        Raises:
            KeyError if the flowsheet is not found
            ValueError on JSON errors
        """
        value = self._find(id_).load()
        _log.debug(f"Flowsheet '{id_}' loaded")
        return value

    def _find(self, id_: str):
        try:
            store = self._id_store[id_]
        except KeyError:
            raise KeyError(f"Unknown flowsheet '{id_}'")
        return store
