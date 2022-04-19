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
Tests for Keras Pruning Utils
"""

import pytest
import numpy as np
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Input
from idaes.core.surrogate.keras.prune import count_nodes, prune_sequential

@pytest.mark.unit
def test_node_count_helper():
    model = Sequential()
    model.add(Input(1))
    model.add(Dense(50))
    model.add(Dense(25, activation='relu'))

    # Input layer is not counted
    assert(count_nodes(model) == 75)

@pytest.mark.unit
def test_pruning():
    model = Sequential()
    model.add(Input(1))
    model.add(Dense(5))
    model.add(Dense(5, activation='relu'))
    model.add(Dense(1))

    """Construct weight matrix where the following nodes are pruned:
        Layer 1 Node 1: All zeros
        Layer 1 Node 2: Forward layer all zeros
        Layer 1 Node 3: Forward layer all zeros once one pruning step occurs
        Layer 2 Node 4: Forward layer all zeros
    """
    l1_w = np.array([[0, 1, 2, 3, 4]])

    b = np.array([1, 2, 3, 4, 5])

    l2_w = np.array([
            [1, 0, 0, 6, 7],
            [4, 0, 0, 9, 10],
            [1, 0, 0, 7, 9],
            [2, 0, 3, 5, 9],
            [7, 0, 0, 8, 10],
        ]).transpose()

    l3_w = np.array([[1, 1, 1, 0, 1]]).transpose()
    l3_b = np.array([1])

    w = [l1_w, b, l2_w, b, l3_w, l3_b]

    model.set_weights(w)

    # Input layer is not counted
    assert(count_nodes(model) == 11)

    pruned_model = prune_sequential(model)
    w = pruned_model.get_weights()
    cfg = pruned_model.get_config()

    # Check correct nodes were removed
    assert(count_nodes(pruned_model) == 7)
    assert(w[0].shape == (1, 2))
    assert(w[2].shape == (2, 4))
    assert(w[1][0] == 4)

    # Check Model parameters
    assert(cfg['layers'][0]['class_name'] == 'InputLayer')
    assert(cfg['layers'][1]['config']['activation'] == 'linear')
    assert(cfg['layers'][2]['config']['activation'] == 'relu')
    assert (cfg['layers'][3]['config']['activation'] == 'linear')


