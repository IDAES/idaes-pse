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
Utility function for the creation of custom sparsification loops for keras models
"""
import numpy as np
from pyomo.common.dependencies import attempt_import

keras, keras_available = attempt_import("tensorflow.keras")

# Takes an N dimensional array and finds the kth lowest magnitude elements
def get_indices_smallest_weights(w, k):

    w = abs(w)
    idx = np.argpartition(w.ravel(), k)
    return np.array(np.unravel_index(idx, w.shape))[:, range(k)].transpose().tolist()

# Gets the total number of zero weights in a model
def count_N_zero_weights(w):

    N = 0

    # Remove biases from the weights array
    weights = [w[2*i] for i in range(len(w)//2)]

    for layer in weights:
        unique, counts = np.unique(layer, return_counts=True)
        mapping = dict(zip(unique, counts))
        if 0 in mapping:
            N += mapping[0]

    return N

# Runs one sparsification step using the desired final sparsity
def sparsify_sequential(model, sparsity, inplace=True):

    # For each layer in the NN sparsify st percent of the weights
    w = model.get_weights()

    # Loop over the weights matrices (i += 2 because w is weights + biases, don't sparsify biases)
    for i in range(0, len(w), 2):
        # Get layer weight matrix
        l_w = w[i]

        # Get the number of nodes which should be zero
        N_sparsify = int(sparsity * len(np.matrix.flatten(l_w)))

        # Get the indexes for the weights with the lowest magnitude and set them to zero
        idxs = get_indices_smallest_weights(l_w, N_sparsify)
        for idx in idxs:
            l_w[tuple(idx)] = 0

        # Update the weight matrix with the new layer weights
        w[i] = l_w

    # print(count_N_zero_weights(model.get_weights()), count_N_zero_weights(w))

    # Update the model weights with the new weight array - if inplace modify directly, else clone (requires recompiling)
    if inplace:
        model.set_weights(w)
        return model
    else:
        new_model = keras.clone_model(model)
        new_model.set_weights(w)
        return new_model
