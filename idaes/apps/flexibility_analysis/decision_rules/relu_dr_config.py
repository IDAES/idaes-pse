#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
from .dr_config import DRConfig
from pyomo.common.config import ConfigValue


class ReluDRConfig(DRConfig):
    r"""
    A class for specifying options for constructing neural network based decision rules
    for use in the flexibility test problem.

    Attributes
    ----------
    n_layers: int
        The number of layers in the neural network (default=: 4)
    n_nodes_per_layer: int
        The number of nodes in each layer of the neural network (default: 4)
    tensorflow_seed: int
        The seed to pass to tensorflow during training
    scale_inputs: bool
        If False, the inputs to the neural network (uncertain parameter values)
        will not be scaled for training (default: True)
    scale_outputs: bool
        If False, the outputs to the neural network (controls)
        will not be scaled for training (default: True)
    epochs: int
        The number of epochs to use in training the neural network (default: 2000)
    batch_size: int
        The batch size to use in training the neural network (default: 20)
    learning_rate: float
        The learning rate for training the neural network (default: None)
    plot_history: bool
        If True, the training history will be plotted (default: False)
    """

    def __init__(
        self,
        description=None,
        doc=None,
        implicit=False,
        implicit_domain=None,
        visibility=0,
    ):
        super().__init__(
            description=description,
            doc=doc,
            implicit=implicit,
            implicit_domain=implicit_domain,
            visibility=visibility,
        )
        self.n_layers: int = self.declare(
            "n_layers", ConfigValue(domain=int, default=4)
        )
        self.n_nodes_per_layer: int = self.declare(
            "n_nodes_per_layer", ConfigValue(domain=int, default=4)
        )
        self.tensorflow_seed: int = self.declare(
            "tensorflow_seed", ConfigValue(domain=int, default=0)
        )
        self.scale_inputs: bool = self.declare(
            "scale_inputs", ConfigValue(domain=bool, default=True)
        )
        self.scale_outputs: bool = self.declare(
            "scale_outputs", ConfigValue(domain=bool, default=True)
        )
        self.epochs: int = self.declare("epochs", ConfigValue(domain=int, default=2000))
        self.batch_size: int = self.declare(
            "batch_size", ConfigValue(domain=int, default=20)
        )
        self.learning_rate = self.declare("learning_rate", ConfigValue(default=None))
        self.plot_history: bool = self.declare(
            "plot_history", ConfigValue(domain=bool, default=False)
        )
