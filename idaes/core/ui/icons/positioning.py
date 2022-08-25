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
# Author: Abdelrahman Elbashandy - aaelbashandy@lbl.gov
#################################################################################


from collections import deque
from idaes.core.ui.flowsheet import FlowsheetSerializer

class Node:
    """A node represents a unit model or an element in JointJs terms"""
    def __init__(self, id: str, level: int = 0, rank: int = 0) -> None:
        self.id = id
        self.level = level # equivalent to row number in a matrix
        self.rank = rank   # equivalent to col number in a matrix
        self.parent_level = -1
        self.parent_rank = -1
        self.direct_parents = set()
        self.direct_children = set()

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
        self._nodes = {}

        self._identify_feeds_products()

        # Total number of units
        self._N = len(unit_models)

        # metadata for the positioning algorithm
        self._allocated_positions = {} # Save used (x,y) positions
        self._max_depth = 0
        self._levels = {}

        # These are how much we move each unit away from each other in each direction
        self._X = 70
        self._Y = 50

        # Computing
        self._assign_default_positions()


    def set_X(self, x: int):
        """sets X to x. X is how much we move in a horizontal direction"""
        self._X = x


    def set_Y(self, y: int):
        """sets Y to y. Y is how much we move in a vertical direction"""
        self._Y = y


    def _assign_default_positions(self):
        """Old implementation of unit model positioning. It positions unit models
           on a diagonal.
        """
        x_pos = 10
        y_pos = 10
        y_starting_pos = 10

        for _, unit_info in self._unit_models.items():
            unit_name, unit_type = unit_info['name'], unit_info['type']
            # If x_pos it greater than 700 then start another diagonal line
            if x_pos >= 700:
                x_pos = 100
                y_pos = y_starting_pos
                y_starting_pos += 100
            else:
                x_pos += 100
                y_pos += 100

            self._allocated_positions[unit_name] = (x_pos, y_pos)


    def _identify_feeds_products(self):
        """Identify feeds and products. This will be useful for multiple
           reasons:
            - feeds will be our entry point to the graph
            - calculate the number of non feed nor product unit models to know 
              the grid size
            - be able to control where feeds and products should be located
        """
        for _, unit_info in self._unit_models.items():
            unit_name, unit_type = unit_info['name'], unit_info['type']
            # Create a node for this unit model and add it to our nodes map
            self._nodes[unit_name] = Node(id = unit_name)
            if unit_type == FlowsheetSerializer.FEED:
                self._feeds.add(unit_name)
            elif unit_type == FlowsheetSerializer.PRODUCT:
                self._products.add(unit_name)


    def _publish_levels_and_ranks(self):
        """This sets each node or unit_model to a level and a rank.
           
           A level is like a trunk of a tree everything will branch off. Each
           level shares the same name of a feed. The feed that shares the name
           of the level doesn't have an effect on the level itself. That feed
           just got lucky to have that level get its name as an identifier.

           A rank is the depth of that node. It's equivalent to a matrix's
           column level.
        
           We also find the max depth of the maximally stretched directed
           cyclic graph.
           This is useful for to know the maximum the number of columns we can
           stretch the graph upon.

           Implementation: BFS - Breadth First Search to know the oldest
           parent(s) each node has because that oldest parent(s) will be the
           main factor of where the node is going to be located.
        """
        queue = deque()

        visited_nodes = set()

        # Add the feed unit models to our queue first as they are our entry point
        for feed_name in self._feeds:
            feed_node = self._nodes[feed_name]
            feed_node.level = feed_name
            self._levels[feed_name] = 0
            queue.append(feed_node)
            visited_nodes.add(feed_name) # add to visited nodes
            self._max_depth = 1 # We know there are feeds

        while queue:
            node = queue.popleft()
            for child in self._adj_list[node.id]:
                assert child in self._nodes, f"child name '{child}' should have a node"

                child_node = self._nodes[child]
                if child not in visited_nodes:
                    # Prepare node and its child relationships
                    child_node.parent_level = node.level
                    child_node.parent_rank = node.rank
                    child_node.level = node.level
                    child_node.rank = node.rank + 1
                    child_node.direct_parents.add(node.id)
                    node.direct_children.add(child_node.id)

                    queue.append(child_node)
                    visited_nodes.add(child) # add to visited nodes
                    self._max_depth = max(self._max_depth, child_node.rank)

                    # This is for just feeds
                    if node.level in self._levels:
                        self._levels[node.level] += 1
        # After finishing this previous loop, all the nodes should be visited
        assert len(visited_nodes) == self._N, "Some nodes aren't accessible through"

        # Merge unused levels
        # .copy() because we are modifying the original self._levels
        for level, count in self._levels.copy().items():
            if count == 0:
                self._levels.pop(level)
                assert len(self._adj_list[level]) > 0, f"Feed '${level}' isn't connected to any unit model"
                # Merge this feed into any of its children's level.
                # self._adj_list[level] is a 'set' remember? Not indexable.
                any_child_node = self._nodes[next(iter(self._adj_list[level]))]
                self._nodes[level].level = any_child_node.level


    def _assign_positions(self):
        """Assign an (x,y) position for each unit_model"""
        
