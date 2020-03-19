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
        return self._data.get(id_, None)

    def compare(self, id_: str, new_flowsheet: str):
        old_json = self.fetch(id_)

        if not old_json:
            self.save(id_, new_flowsheet)
            return

        diff_model, model_json = compare_models(old_json, new_flowsheet)
        new_json = model_jointjs_conversion(diff_model, model_json)

        self.save(id_, new_json)
