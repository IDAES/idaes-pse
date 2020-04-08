"""
Backend logic for flowsheet visualization server
"""
import json

from typing import Dict

from idaes.ui.flowsheet_comparer import compare_models, model_jointjs_conversion


class DataStorage:
    """Trivial data storage.
    """
    def __init__(self):
        self._data = {}

    def save(self, id_: str, data: Dict):
        self._data[id_] = data

    def fetch(self, id_: str):
        data = self._data.get(id_, None)

        # If there is no data then return an empty model
        if not data:
            data = {"cells": [], "model": {"unit_models": {}, "arcs": {}, "id": id_}}
        return data

    def update(self, id_: str, new_flowsheet: str):
        old_json = self.fetch(id_)

        if not old_json:
            self.save(id_, new_flowsheet)
            return

        diff_model, model_json = compare_models(old_json, new_flowsheet)
        new_json = model_jointjs_conversion(diff_model, model_json)

        self.save(id_, new_json)
