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
Utility functions to prune nodes from keras models
"""

from tensorflow.keras.models import Sequential
import numpy as np

"""
Finds all nodes which have constant output in a layer.
These nodes may affect NN output (i.e. non-zero bias, activation(0) != 0, activation_forward(0) != 0, etc) but do not
change with the inputs of the neural network.
"""
def get_constant_nodes(l_w, l_w_forward):
    # Reshape the layer so it is in the form [node #, node weights]
    l_w = np.transpose(l_w)

    # If there is a forward layer (IE not last layer) reshape it to (forward node, prev node weight)
    if l_w_forward is not None:
        l_w_forward = np.transpose(l_w_forward)

    # Inactive nodes will be labelled by their index
    inactive_nodes = []

    # Get node number indexes of all inactive nodes
    for node_idx, node_w in enumerate(l_w):
        node_w = node_w.flatten()

        # If the weight array is filled with 0s it has a constant output and should be pruned
        if np.count_nonzero(node_w) == 0:
            inactive_nodes.append(node_idx)

        # Else check the upstream nodes. If all upstream nodes have w = 0 for the downstream node, it is inactive
        elif l_w_forward is not None:

            remove_node = True

            # Check the node_idx column in each forward node for a weight that is not 0
            for forward_node in l_w_forward:

                # If there is a node in the forward layer that doesnt have w = 0 for the current node, don't prune
                if forward_node[node_idx] != 0:
                    remove_node = False
                    break

            # All nodes in the forward layer zero out this node's output, prune
            if remove_node:
                inactive_nodes.append(node_idx)

    return inactive_nodes


# Construct weight matrix for new layer
def rebuild_weights(w, remove_nodes):
    for l, removal_indexes in enumerate(remove_nodes.values()):
        # Remove biases
        w[2 * l + 1] = np.delete(w[2 * l + 1], removal_indexes, axis=0)

        # Remove weight columns (current layer)
        w[2 * l] = np.delete(w[2 * l], removal_indexes, axis=1)

        # Remove weights of nodes from previous layer that were removed
        if l != 0:
            w[2 * l] = np.delete(w[2 * l], last_index, axis=0)

        last_index = removal_indexes

    return w


# Constructs a new compiled neural net with fewer nodes using an updated weight matrix
def build_new_NN(model_old, w_new):
    config = model_old.get_config()
    for i, layer_config in enumerate(config['layers']):

        # Check if the layer has a units key
        if 'config' in layer_config and 'units' in layer_config['config']:
            # Config has extra input layer (not in w matrix) so adjust the index to match previous functions
            l_idx = i - 1

            # Get the new number of neurons by looking at the weights matrix
            neurons = w_new[l_idx * 2].shape[-1]

            # Change the layer config
            layer_config['config']['units'] = neurons
            config['layers'][i] = layer_config

    model = Sequential.from_config(config)

    # Set the weights of the model
    model.set_weights(w_new)
    return model


# Helper function to count the number of nodes in a neural network
def count_nodes(model):
    cfg = model.get_config()
    N_nodes = 0
    for i, layer_config in enumerate(cfg['layers']):
        if 'config' in layer_config and 'units' in layer_config['config']:
            N_nodes += layer_config['config']['units']
    return N_nodes


# Prune the nodes that have zero weights or are zeroed out by future nodes
def pruning_step(model):
    w = model.get_weights()

    # Container to keep track of all nodes within a layer that must be removed, indexed by [layer][idx node to remove]
    nodes_to_remove = {}

    # Conduct one pruning step on the model
    for l in range(0, len(w), 2):
        l_w = w[l]

        # If there is a layer ahead of this layer, get the forward layer weights
        l_w_forward = w[l + 2] if l < len(w) - 2 else None

        # Get an array of indexes for the nodes that must be removed and append to the main list
        nodes_to_remove[l] = get_constant_nodes(l_w, l_w_forward)

        # Check for removal of input layer nodes (l=0, non-empty nodes_to_remove)
        """if l == 0 and nodes_to_remove:
            print(f"Warning: input nodes {nodes_to_remove} have constant output and are being removed."
                  f"The returned NN will have a different input shape.")

        # Check for the removal of output layer nodes (l=len(w)-2, non-empty nodes_to_remove)
        if l == len(w) - 2 and nodes_to_remove:
            # Get bias value for removed nodes as this is the constant value of the output
            biases_removed_nodes = [l_b[idx] for idx in nodes_to_remove]

            print(
                f"Warning: output nodes {nodes_to_remove} have constant output and are being removed from the NN."
                f"This indicates the output is constant with the given input variables and is equal to the bias.\n"
                f"The constant value for these outputs was: {biases_removed_nodes}")"""

    w_new = rebuild_weights(w, nodes_to_remove)

    return build_new_NN(model, w_new)

# Helper function for the user to compare nodes in the pruned vs unpruned models
def count_nodes(model):
    cfg = model.get_config()
    N_nodes = 0
    for i, layer_config in enumerate(cfg['layers']):
        if 'config' in layer_config and 'units' in layer_config['config']:
            N_nodes += layer_config['config']['units']
    return N_nodes


"""
Prunes nodes from a keras sequential model by removing nodes with constant outputs.
While inactive nodes have an affect on the output, they are pruned because their contribution can be shifted to active
nodes. The model should be retrained after pruning to recover accuracy after removal of nodes.
"""
def prune_sequential(model, max_steps=100, verbose=0):

    # Do while loop to continue pruning until model nodes are constant or max pruning steps occur to prevent inf loop
    for step in range(max_steps):

        initial_nodes = count_nodes(model)
        model = pruning_step(model)
        final_nodes = count_nodes(model)

        removed_nodes = initial_nodes - final_nodes

        if verbose == 1:
            print(f"Pruning step {step+1} removed {removed_nodes} nodes")

        if removed_nodes == 0:
            break

    return model
