#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Interface for importing Linear-tree models into IDAES
"""
# TODO: Missing docstrings
# pylint: disable=missing-class-docstring
# pylint: disable=missing-function-docstring

import json
import os.path
import pickle

import pandas as pd

from pyomo.common.dependencies import attempt_import

from idaes.core.surrogate.sampling.scaling import OffsetScaler

from idaes.core.surrogate.omlt_base_surrogate_class import OMLTSurrogate

lt, lt_available = attempt_import("lineartree")
omlt, omlt_available = attempt_import("omlt")

if omlt_available:
    from omlt.linear_tree import (
        LinearTreeDefinition,
        LinearTreeGDPFormulation,
        LinearTreeHybridBigMFormulation,
    )


class LinearTreeSurrogate(OMLTSurrogate):
    def __init__(
        self,
        lt_model,
        input_labels,
        output_labels,
        input_bounds,
        input_scaler=None,
        output_scaler=None,
    ):
        """
        Standard SurrogateObject for surrogates based on Linear Tree models.
        Utilizes the OMLT framework for importing Linear Tree models to IDAES.

        Contains methods to both populate a Pyomo Block with constraints
        representing the surrogate and to evaluate the surrogate a set of user
        provided points.

        This constructor should only be used when first creating the surrogate within IDAES.
        Once created, this object can be stored to disk using save_to_folder and loaded
        with load_from_folder

        Args:
           lt_model: Linear-tree model
              This is the Linear-tree model that will be loaded.
           input_labels: list of str
              The ordered list of labels corresponding to the inputs in the linear-tree model
           output_labels: list of str
              The ordered list of labels corresponding to the outputs in the linear-tree model
           input_bounds: None of dict of tuples
              Keys correspond to each of the input labels and values are the tuples of
              bounds (lb, ub)
           input_scaler: None or OffsetScaler
              The scaler to be used for the inputs. If None, then no scaler is used
           output_scaler: None of OffsetScaler
              The scaler to be used for the outputs. If None, then no scaler is used
        """
        super().__init__(
            input_labels=input_labels,
            output_labels=output_labels,
            input_bounds=input_bounds,
            input_scaler=input_scaler,
            output_scaler=output_scaler,
        )
        self._lt_model = lt_model

    def populate_block(self, block, additional_options=None):
        """
        Method to populate a Pyomo Block with the linear-tree model constraints.

        Args:
           block: Pyomo Block component
              The block to be populated with variables and/or constraints.
           additional_options: dict or None
              If not None, then should be a dict with the following keys;
              'formulation': LinearTreeSurrogate.Formulation
              The formulation to use with OMLT. Possible values are LINEAR_TREE_GDP_BIGM,
              LINEAR_TREE_GDP_HULL, LINEAR_TREE_GDP_MBIGM, or LINEAR_TREE_HYBRID_BIGM
              (default is LINEAR_TREE_GDP_BIGM)
        """
        formulation = additional_options.pop(
            "formulation", LinearTreeSurrogate.Formulation.LINEAR_TREE_GDP_BIGM
        )
        omlt_scaling, scaled_input_bounds = self.generate_omlt_scaling_objecets()
        scaled_keys = list(scaled_input_bounds.keys())
        unscaled_keys = list(self.input_bounds().keys())
        unscaled_input_bounds = {
            scaled_keys[idx]: self.input_bounds()[unscaled_keys[idx]]
            for idx, _ in enumerate(scaled_keys)
        }
        lt = LinearTreeDefinition(
            self._lt_model,
            scaling_object=omlt_scaling,
            scaled_input_bounds=scaled_input_bounds,
            unscaled_input_bounds=unscaled_input_bounds,
        )

        if formulation == LinearTreeSurrogate.Formulation.LINEAR_TREE_GDP_BIGM:
            formulation_object = LinearTreeGDPFormulation(lt, transformation="bigm")
        elif formulation == LinearTreeSurrogate.Formulation.LINEAR_TREE_GDP_HULL:
            formulation_object = LinearTreeGDPFormulation(lt, transformation="hull")
        elif formulation == LinearTreeSurrogate.Formulation.LINEAR_TREE_HYBRID_BIGM:
            formulation_object = LinearTreeHybridBigMFormulation(lt)
        else:
            raise ValueError(
                'An unrecognized formulation "{}" was passed to '
                "LinearTreeSurrogate.populate_block. Please pass a valid "
                "formulation.".format(formulation)
            )
        self.populate_block_with_model(block, formulation_object)

    def evaluate_surrogate(self, inputs):
        """
        Method to evaluate linear-tree model at a set of input values.

        Args:
            inputs: numpy array of input values. First dimension of array
                must match the number of input variables.

        Returns:
            outputs: numpy array of values for all outputs evaluated at input
                points.
        """
        x = inputs
        if self._input_scaler is not None:
            x = self._input_scaler.scale(x)
        y = self._lt_model.predict(x.to_numpy())

        # y is a numpy array, make it a dataframe
        y = pd.DataFrame(
            data=y, columns=self.output_labels(), index=inputs.index, dtype="float64"
        )
        if self._output_scaler is not None:
            y = self._output_scaler.unscale(y)
        return y

    def save_to_folder(self, lt_folder_name, lt_model_name="idaes_linear_tree_model"):
        """
        Save the surrogate object to disk by providing the name of the
        folder to contain the linear-tree model and additional IDAES metadata

        Args:
           folder_name: str
              The name of the folder to contain the linear-tree model and additional
              IDAES metadata
        """
        with open(os.path.join(lt_folder_name, lt_model_name + ".pkl"), "wb") as FILE:
            pickle.dump(self._lt_model, FILE)
        info = dict()
        info["input_scaler"] = None
        if self._input_scaler is not None:
            info["input_scaler"] = self._input_scaler.to_dict()
        info["output_scaler"] = None
        if self._output_scaler is not None:
            info["output_scaler"] = self._output_scaler.to_dict()

        # serialize information from the base class
        info["input_labels"] = self.input_labels()
        info["output_labels"] = self.output_labels()
        info["input_bounds"] = self.input_bounds()

        with open(os.path.join(lt_folder_name, "idaes_info.json"), "w") as fd:
            json.dump(info, fd)

    @classmethod
    def load_from_folder(cls, lt_folder_name, lt_model_name="idaes_linear_tree_model"):
        """
        Load the surrogate object from disk by providing the name of the
        folder holding the linear-tree model

        Args:
           folder_name: str
              The name of the folder containing the Linear-tree model and additional
              IDAES metadata

        Returns: an instance of LinearTreeSurrogate
        """

        with open(os.path.join(lt_folder_name, lt_model_name + ".pkl"), "rb") as FILE:
            lt_model = pickle.load(FILE)

        with open(os.path.join(lt_folder_name, "idaes_info.json")) as fd:
            info = json.load(fd)

        input_scaler = None
        if info["input_scaler"] is not None:
            input_scaler = OffsetScaler.from_dict(info["input_scaler"])

        output_scaler = None
        if info["output_scaler"] is not None:
            output_scaler = OffsetScaler.from_dict(info["output_scaler"])

        return LinearTreeSurrogate(
            lt_model=lt_model,
            input_labels=info["input_labels"],
            output_labels=info["output_labels"],
            input_bounds=info["input_bounds"],
            input_scaler=input_scaler,
            output_scaler=output_scaler,
        )


def save_linear_tree_pickle(lt, path, name):
    with open(os.path.join(path, "{}.pickle".format(name)), "wb") as file:
        pickle.dump(lt, file)


def load_linear_tree_pickle(path, name):
    with open(os.path.join(path, "{}.pickle".format(name)), "rb") as file:
        lt = pickle.load(file)
    return lt
