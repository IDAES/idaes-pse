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
Tests for Keras Sparsification Utils
"""

import pytest
from tensorflow.keras.layers import Dense, Input
from tensorflow.keras.models import Sequential
from idaes.core.surrogate.keras.sparsify import *

@pytest.mark.unit
def test_sparsification():
    model = Sequential()
    model.add(Input(1))
    model.add(Dense(5))
    model.add(Dense(5, activation='relu'))
    model.add(Dense(1))

    l1_w = np.array([[0, 1, 2, 3, 4]])

    b = np.array([1, 2, 3, 4, 5])

    l2_w = np.array([
        [-1, 5, 10, 6, 7],
        [4, -1, 15, 9, 10],
        [1, -2, 0, 7, 9],
        [2, 10, 3, 5, 9],
        [7, 11, 16, -8, 10]
    ]).transpose()

    l3_w = np.array([[1, -5, 6, -8, 3]]).transpose()
    l3_b = np.array([1])

    w = [l1_w, b, l2_w, b, l3_w, l3_b]

    model.set_weights(w)

    # Check helper function to count zero weights
    assert(count_N_zero_weights(model.get_weights()) == 2)

    # 3 weights from l1, 15 from l2, 3 from l3
    s = 3.0/5.0
    sparse = sparsify_sequential(model, s)

    l1_sparse = np.array([[0, 0, 0, 3, 4]])

    l2_sparse = np.array([
        [0, 0, 10, 0, 0],
        [0, 0, 15, 9, 10],
        [0, 0, 0, 0, 9],
        [0, 10, 0, 0, 9],
        [0, 11, 16, 0, 10]
    ]).transpose()

    l3_sparse = np.array([[0, 0, 6, -8, 0]]).transpose()
    w_sparse = model.get_weights()

    # Check weights
    np.testing.assert_allclose(w_sparse[0], l1_sparse)
    np.testing.assert_allclose(w_sparse[2], l2_sparse)
    np.testing.assert_allclose(w_sparse[4], l3_sparse)

    # Check biases are unchanged
    np.testing.assert_allclose(w_sparse[1], b)
    np.testing.assert_allclose(w_sparse[3], b)
    np.testing.assert_allclose(w_sparse[5], l3_b)