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
