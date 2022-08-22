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
#
# Author: Abdelrahman Elbashandy
#################################################################################



from itertools import product
from idaes.core.ui.flowsheet import FlowsheetSerializer


class UnitModelsPositioning:
    """Represents icon positioning information for a given unit model."""

    def __init__(self, adj_list: dict, unit_models: dict):
        """Construct the layout of a given directed graphical flowsheet.

        Args:
            adj_list: Adjacency list, a directed graph representation for the given flowsheet
            unit_models: FlowsheetSerializer's self.unit_models to get each unit_model's info

        """
        self._adj_list = adj_list
        self._unit_models = unit_models

        self._feeds = set()
        self._products = set()

        self._identify_feeds_products()

        # Total number of units
        self._N = len(unit_models)
        # number of non feed nor product unit models
        self._n = self._N - len(feeds) - len(products)

    def _identify_feeds_products(self):
        """Identify feeds and products. This will be useful for multiple
           reasons:
            - feeds will be our entry point to the graph
            - calculate the number of non feed nor product unit models to know 
              the grid size
        """
        for _, unit_info in self._unit_models.items():
            unit_name, unit_type = unit_info['name'], unit_info['type']
            if unit_type == FlowsheetSerializer.FEED:
                self._feeds.add(unit_name)
            elif unit_type == FlowsheetSerializer.PRODUCT:
                self._products.add(unit_name)

    def _locate_positions(self):
        """Locating the position of each unit model"""
        

