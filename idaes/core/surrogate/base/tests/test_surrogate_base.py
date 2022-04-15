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
Tests for SurrogateBase classes
"""
import pytest

import pandas as pd

from idaes.core.surrogate.base.surrogate_base import SurrogateTrainer, SurrogateBase

tdata = {"x1": [1, 2, 3, 4], "x2": [10, 20, 30, 40], "z1": [11, 22, 33, 44]}
training_data = pd.DataFrame.from_dict(tdata)

vdata = {"x1": [5], "x2": [50], "z1": [55]}
validation_data = pd.DataFrame.from_dict(vdata)


class TestSurrogateTrainer:
    @pytest.fixture(scope="class")
    def trainer(self):
        return SurrogateTrainer(
            input_labels=["x1", "x2"],
            output_labels=["z1"],
            training_dataframe=training_data,
            validation_dataframe=validation_data,
        )

    @pytest.mark.unit
    def test_surrogate_trainer_init(self, trainer):
        assert trainer._input_labels == ["x1", "x2"]
        assert trainer._output_labels == ["z1"]
        assert trainer._training_dataframe is training_data
        assert trainer._validation_dataframe is validation_data
        assert trainer._input_bounds == {"x1": (1, 4), "x2": (10, 40)}

    @pytest.mark.unit
    def test_surrogate_trainer_no_validation(self):
        trainer = SurrogateTrainer(
            input_labels=["x1", "x2"],
            output_labels=["z1"],
            training_dataframe=training_data,
        )

        assert trainer._input_labels == ["x1", "x2"]
        assert trainer._output_labels == ["z1"]
        assert trainer._training_dataframe is training_data
        assert trainer._validation_dataframe is None
        assert trainer._input_bounds == {"x1": (1, 4), "x2": (10, 40)}

    @pytest.mark.unit
    def test_surrogate_trainer_no_input_labels(self):
        with pytest.raises(
            ValueError,
            match="SurrogateTrainer requires a list of "
            "input_labels and a list of output_labels which "
            "must both have a length of at least one",
        ):
            SurrogateTrainer(
                input_labels=None,
                output_labels=["z1"],
                training_dataframe=training_data,
            )

    @pytest.mark.unit
    def test_surrogate_trainer_input_labels_len_0(self):
        with pytest.raises(
            ValueError,
            match="SurrogateTrainer requires a list of "
            "input_labels and a list of output_labels which "
            "must both have a length of at least one",
        ):
            SurrogateTrainer(
                input_labels=[], output_labels=["z1"], training_dataframe=training_data
            )

    @pytest.mark.unit
    def test_surrogate_trainer_no_output_labels(self):
        with pytest.raises(
            ValueError,
            match="SurrogateTrainer requires a list of "
            "input_labels and a list of output_labels which "
            "must both have a length of at least one",
        ):
            SurrogateTrainer(
                input_labels=["x1", "x2"],
                output_labels=None,
                training_dataframe=training_data,
            )

    @pytest.mark.unit
    def test_surrogate_trainer_output_labels_len_0(self):
        with pytest.raises(
            ValueError,
            match="SurrogateTrainer requires a list of "
            "input_labels and a list of output_labels which "
            "must both have a length of at least one",
        ):
            SurrogateTrainer(
                input_labels=["x1", "x2"],
                output_labels=[],
                training_dataframe=training_data,
            )

    @pytest.mark.unit
    def test_surrogate_trainer_labels_overlap(self):
        with pytest.raises(
            ValueError,
            match="Duplicate label found in input_labels " "and/or output_labels.",
        ):
            SurrogateTrainer(
                input_labels=["x1", "x2"],
                output_labels=["x1"],
                training_dataframe=training_data,
            )

    @pytest.mark.unit
    def test_surrogate_trainer_input_label_mismatch(self):
        with pytest.raises(
            ValueError,
            match="The following input labels were not found in "
            "the training data columns: {'foo'}",
        ):
            SurrogateTrainer(
                input_labels=["x1", "foo"],
                output_labels=["z1"],
                training_dataframe=training_data,
            )

    @pytest.mark.unit
    def test_surrogate_trainer_input_label_mismatch_verification(self):
        vdata_x = {"x1": [5], "foo": [50], "z1": [55]}
        validation_data_x = pd.DataFrame.from_dict(vdata_x)

        with pytest.raises(
            ValueError,
            match="The following input labels were not found in "
            "the validation data columns: {'x2'}",
        ):
            SurrogateTrainer(
                input_labels=["x1", "x2"],
                output_labels=["z1"],
                training_dataframe=training_data,
                validation_dataframe=validation_data_x,
            )

    @pytest.mark.unit
    def test_surrogate_trainer_output_label_mismatch(self):
        with pytest.raises(
            ValueError,
            match="The following output labels were not found in "
            "the training data columns: {'foo'}.",
        ):
            SurrogateTrainer(
                input_labels=["x1", "x2"],
                output_labels=["foo"],
                training_dataframe=training_data,
            )

    @pytest.mark.unit
    def test_surrogate_trainer_output_label_mismatch_verification(self):
        vdata_x = {"x1": [5], "x2": [50], "foo": [55]}
        validation_data_x = pd.DataFrame.from_dict(vdata_x)

        with pytest.raises(
            ValueError,
            match="The following output labels were not found in "
            "the validation data columns: {'z1'}.",
        ):
            SurrogateTrainer(
                input_labels=["x1", "x2"],
                output_labels=["z1"],
                training_dataframe=training_data,
                validation_dataframe=validation_data_x,
            )

    @pytest.mark.unit
    def test_surrogate_trainer_w_bounds(self):
        trainer = SurrogateTrainer(
            input_labels=["x1", "x2"],
            output_labels=["z1"],
            input_bounds={"x1": (0, 5), "x2": (0, 50)},
            training_dataframe=training_data,
            validation_dataframe=validation_data,
        )

        assert trainer._input_bounds == {"x1": (0, 5), "x2": (0, 50)}

    @pytest.mark.unit
    def test_n_inputs(self, trainer):
        assert trainer.n_inputs() == 2

    @pytest.mark.unit
    def test_n_outputs(self, trainer):
        assert trainer.n_outputs() == 1

    @pytest.mark.unit
    def test_input_labels(self, trainer):
        assert trainer.input_labels() == trainer._input_labels

    @pytest.mark.unit
    def test_output_labels(self, trainer):
        assert trainer.output_labels() == trainer._output_labels

    @pytest.mark.unit
    def test_input_bounds(self, trainer):
        assert trainer.input_bounds() == trainer._input_bounds

    @pytest.mark.unit
    def test_train_surrogate(self, trainer):
        with pytest.raises(
            NotImplementedError,
            match="train_surrogate called, but not implemented " "on the derived class",
        ):
            trainer.train_surrogate()


class TestSurrogateBase:
    @pytest.fixture(scope="class")
    def surrogate(self):
        return SurrogateBase(
            input_labels=["x1", "x2"],
            output_labels=["z1"],
            input_bounds={"x1": (0, 5), "x2": (0, 50)},
        )

    @pytest.mark.unit
    def test_init(self, surrogate):
        assert surrogate._input_labels == ["x1", "x2"]
        assert surrogate._output_labels == ["z1"]
        assert surrogate._input_bounds == {"x1": (0, 5), "x2": (0, 50)}

    @pytest.mark.unit
    def test_n_inputs(self, surrogate):
        assert surrogate.n_inputs() == 2

    @pytest.mark.unit
    def test_n_outputs(self, surrogate):
        assert surrogate.n_outputs() == 1

    @pytest.mark.unit
    def test_input_labels(self, surrogate):
        assert surrogate.input_labels() == surrogate._input_labels

    @pytest.mark.unit
    def test_output_labels(self, surrogate):
        assert surrogate.output_labels() == surrogate._output_labels

    @pytest.mark.unit
    def test_input_bounds(self, surrogate):
        assert surrogate.input_bounds() == surrogate._input_bounds

    @pytest.mark.unit
    def test_populate_block(self, surrogate):
        with pytest.raises(
            NotImplementedError,
            match="SurrogateModel class has not implemented " "populate_block method.",
        ):
            surrogate.populate_block("foo")

    @pytest.mark.unit
    def test_evaluate_surrogate(self, surrogate):
        with pytest.raises(
            NotImplementedError,
            match="SurrogateModel class has not implemented an "
            "evaluate_surrogate method.",
        ):
            surrogate.evaluate_surrogate(dataframe=training_data)

    @pytest.mark.unit
    def test_save(self, surrogate):
        with pytest.raises(
            NotImplementedError,
            match='"save" should be implemented in the'
            " class derived from SurrogateBase",
        ):
            surrogate.save(None)

    @pytest.mark.unit
    def test_load(self):
        with pytest.raises(
            NotImplementedError,
            match='"load" should be implemented in the'
            " class derived from SurrogateBase",
        ):
            SurrogateBase.load(None)
