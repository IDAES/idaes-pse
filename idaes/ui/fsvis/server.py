"""
Backend logic for flowsheet visualization server
"""
from typing import Dict


class DataStorage:
    """Trivial data storage.
    """
    def __init__(self):
        self._data = {}

    def save(self, id_: str, data: Dict):
        self._data[id_] = data

    def fetch(self, id_: str):
        return self._data.get(id_, None)


