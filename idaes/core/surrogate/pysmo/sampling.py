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

from __future__ import division, print_function
from six import string_types

# from builtins import int, str
import numpy as np
import pandas as pd
import warnings
import itertools

__author__ = "Oluwamayowa Amusat"


class FeatureScaling:
    """

    A class for scaling and unscaling input and output data. The class contains three main functions
    """

    def __init__(self):
        pass

    @staticmethod
    def data_scaling_minmax(data):
        """

        This function performs column-wise minimax scaling on the input dataset.

            Args:
                data (NumPy Array or Pandas Dataframe): The input data set to be scaled. Must be a numpy array or dataframe.

            Returns:
                scaled_data(NumPy Array): A 2-D numpy array containing the scaled data. All array values will be between [0, 1].
                data_minimum(NumPy Array): A 2-D row vector containing the column-wise minimums of the input data
                data_maximum(NumPy Array): A 2-D row vector containing the column-wise maximums of the input data

            Raises:
                TypeError: Raised when the input data is not a numpy array or dataframe
        """
        # Confirm that data type is an array or DataFrame
        if isinstance(data, np.ndarray):
            input_data = data
        elif isinstance(data, pd.DataFrame):
            input_data = data.values
        else:
            raise TypeError(
                "original_data_input: Pandas dataframe or numpy array required."
            )

        if input_data.ndim == 1:
            input_data = input_data.reshape(len(input_data), 1)
        data_minimum = np.min(input_data, axis=0)
        data_maximum = np.max(input_data, axis=0)
        scale = data_maximum - data_minimum
        scale[scale == 0.0] = 1.0
        scaled_data = (input_data - data_minimum) / scale
        # scaled_data = (input_data - data_minimum) / (data_maximum - data_minimum)
        data_minimum = data_minimum.reshape(1, data_minimum.shape[0])
        data_maximum = data_maximum.reshape(1, data_maximum.shape[0])
        return scaled_data, data_minimum, data_maximum

    @staticmethod
    def data_unscaling_minmax(x_scaled, x_min, x_max):
        """

        This function performs column-wise un-scaling on the a minmax-scaled input dataset.

            Args:
                x_scaled(NumPy Array): The input data set to be un-scaled. Data values should be between 0 and 1.
                x_min(NumPy Array): 1-D or 2-D (n-by-1) vector containing the actual minimum value for each column. Must contain same number of elements as the number of columns in x_scaled.
                x_max(NumPy Array): 1-D or 2-D (n-by-1) vector containing the actual maximum value for each column. Must contain same number of elements as the number of columns in x_scaled.

            Returns:
                unscaled_data(NumPy Array): A 2-D numpy array containing the scaled data, unscaled_data = x_min + x_scaled * (x_max - x_min)

            Raises:
                IndexError: Function raises index error when the dimensions of the arrays are inconsistent.
        """
        # Check if it can be evaluated. Will return index error if dimensions are wrong
        if x_scaled.ndim == 1:  # Check if 1D, and convert to 2D if required.
            x_scaled = x_scaled.reshape(len(x_scaled), 1)
        if (x_scaled.shape[1] != x_min.size) or (x_scaled.shape[1] != x_max.size):
            raise IndexError("Dimensionality problems with data for un-scaling.")
        unscaled_data = x_min + x_scaled * (x_max - x_min)
        return unscaled_data

    # @staticmethod
    # def data_scaling_standardization(data):
    #     # Confirm that data type is an array or DataFrame
    #     if isinstance(data, np.ndarray):
    #         input_data = data
    #     elif isinstance(data, pd.DataFrame):
    #         input_data = data.values
    #     else:
    #         raise TypeError('original_data_input: Pandas dataframe or numpy array required.')
    #
    #     if input_data.ndim == 1:
    #         input_data = input_data.reshape(len(input_data), 1)
    #
    #     data_mean = np.mean(input_data, axis=0)
    #     data_stdev = np.std(input_data, axis=0)
    #     scaled_data = (input_data - data_mean) / data_stdev
    #     data_mean = data_mean.reshape(1, data_mean.shape[0])
    #     data_stdev = data_stdev.reshape(1, data_stdev.shape[0])
    #     return scaled_data, data_mean, data_stdev


class SamplingMethods:
    def nearest_neighbour(self, full_data, a):
        """
        Function determines the closest point to a in data_input (user provided data).
        This is done by determining the input data with the smallest L2 distance from a.

        The function:
        1. Calculates the L2 distance between all the input data points and a,
        2. Sorts the input data based on the calculated L2-distances, and
        3. Selects the sample point in the first row (after sorting) as the closest sample point.

        Args:
            self: contains, among other things, the input data.
            full_data: refers to the input dataset supplied by the user.
            a: a single row vector containing the sample point we want to find the closest sample to.

        Returns:
            closest_point: a row vector containing the closest point to a in self.x_data
        """
        no_y_vars = self.x_data.shape[1] - full_data.shape[1]
        dist = full_data[:, :no_y_vars] - a
        l2_norm = np.sqrt(np.sum((dist**2), axis=1))
        l2_norm = l2_norm.reshape(l2_norm.shape[0], 1)
        distances = np.append(full_data, l2_norm, 1)
        sorted_distances = distances[distances[:, -1].argsort()]
        closest_point = sorted_distances[0, :-1]
        return closest_point

    def points_selection(self, full_data, generated_sample_points):
        """
        Uses L2-distance evaluation (implemented in nearest_neighbour) to find closest available points in original data to those generated by the sampling technique.
        Calls the nearest_neighbour function for each row in the input data.

        Args:
            full_data: refers to the input dataset supplied by the user.
            generated_sample_points(NumPy Array): The vector of points (number_of_sample rows) for which the closest points in the original data are to be found. Each row represents a sample point.

        Returns:
            equivalent_points: Array containing the points (in rows) most similar to those in generated_sample_points
        """

        equivalent_points = np.zeros(
            (generated_sample_points.shape[0], len(self.data_headers))
        )
        for i in range(0, generated_sample_points.shape[0]):
            closest_point = self.nearest_neighbour(
                full_data, generated_sample_points[i, :]
            )
            equivalent_points[i, :] = closest_point
        return equivalent_points

    def sample_point_selection(self, full_data, sample_points, sampling_type):
        if sampling_type == "selection":
            sd = FeatureScaling()
            scaled_data, data_min, data_max = sd.data_scaling_minmax(full_data)
            points_closest_scaled = self.points_selection(scaled_data, sample_points)
            points_closest_unscaled = sd.data_unscaling_minmax(
                points_closest_scaled, data_min, data_max
            )

            unique_sample_points = np.unique(points_closest_unscaled, axis=0)
            if unique_sample_points.shape[0] < points_closest_unscaled.shape[0]:
                warnings.warn(
                    "The returned number of samples is less than the requested number due to repetitions during nearest neighbour selection."
                )
            print(
                "\nNumber of unique samples returned by sampling algorithm:",
                unique_sample_points.shape[0],
            )

        elif sampling_type == "creation":
            sd = FeatureScaling()
            unique_sample_points = sd.data_unscaling_minmax(
                sample_points, full_data[0, :], full_data[1, :]
            )

        return unique_sample_points

    def prime_number_generator(self, n):
        """
        Function generates a list of the first n prime numbers

            Args:
                n(int): Number of prime numbers required

            Returns:
                prime_list(list): A list of the first n prime numbers

        Example: Generate first three prime numbers
            >>  prime_number_generator(3)
            >> [2, 3, 5]

        """
        # Alternative way of generating primes using list generators
        # prime_list = []
        # current_no = 2
        # while len(prime_list) < n:
        #     matching_objs = next((o for o in range(2, current_no) if current_no % o == 0), 0)
        #     if matching_objs==0:
        #         prime_list.append(current_no)
        #     current_no += 1

        prime_list = []
        current_no = 2
        while len(prime_list) < n:
            for i in range(2, current_no):
                if (current_no % i) == 0:
                    break
            else:
                prime_list.append(current_no)
            current_no += 1
        return prime_list

    def base_conversion(self, a, b):
        """
        Function converts integer a from base 10 to base b

            Args:
                a(int): Number to be converted, base 10
                b(int): Base required

            Returns:
                string_representation(list): List containing strings of individual digits of "a" in the new base "b"

        Examples: Convert (i) 5 to base 2 and (ii) 57 to base 47
            >>  base_conversion(5, 2)
            >> ['1', '0', '1']

            >>  base_conversion(57, 47)
            >> ['1', '10']

        """

        string_representation = []
        if a < b:
            string_representation.append(str(a))
        else:
            while a > 0:
                a, c = (a // b, a % b)
                string_representation.append(str(c))
            string_representation = string_representation[::-1]
        return string_representation

    def prime_base_to_decimal(self, num, base):
        """
        ===============================================================================================================
        Function converts a fractional number "num" in base "base" to base 10. Reverses the process in base_conversion
        Note: The first string element is ignored, since this would be zero for a fractional number.

            Args:
                num(list): Number in base b to be converted. The number must be represented as a list containing individual digits of the base, with the first entry as zero.
                b(int): Original base

            Returns:
                decimal_equivalent(float): Fractional number in base 10

        Examples:
        Convert 0.01 (base 2) to base 10
            >>  prime_base_to_decimal(['0', '0', '1'], 2)  # Represents 0.01 in base 2
            >> 0.25

        Convert 0.01 (base 20) to base 10
            >>  prime_base_to_decimal(['0', '0', '1'], 20)  # Represents 0.01 in base 20
            >> 0.0025
        ================================================================================================================

        """
        binary = num
        decimal_equivalent = 0
        # Convert fractional part decimal equivalent
        for i in range(1, len(binary)):
            decimal_equivalent += int(binary[i]) / (base**i)
        return decimal_equivalent

    def data_sequencing(self, no_samples, prime_base):
        """
        ===============================================================================================================
        Function which generates the first no_samples elements of the Halton or Hammersley sequence based on the prime number prime_base
        The steps for generating the first no_samples of the sequence are as follows:
        1. Create a list of numbers between 0 and no_samples --- nums = [0, 1, 2, ..., no_samples]
        2. Convert each element in nums into its base form based on the prime number prime_base, reverse the base digits of each number in num
        3. Add a decimal point in front of the reversed number
        4. Convert the reversed numbers back to base 10
            Args:
                no_samples(int): Number of Halton/Hammersley sequence elements required
                prime_base(int): Current prime number to be used as base

            Returns:
                sequence_decimal(NumPy Array): 1-D array containing the first no_samples elements of the sequence based on prime_base

        Examples:
        First three elements of the Halton sequence based on base 2
            >>  data_sequencing(self, 3, 2)
            >> [0, 0.5, 0.75]
        ================================================================================================================

        """
        pure_numbers = np.arange(0, no_samples)
        bitwise_rep = []
        reversed_bitwise_rep = []
        sequence_bitwise = []
        sequence_decimal = np.zeros((no_samples, 1))
        for i in range(0, no_samples):
            base_rep = self.base_conversion(pure_numbers[i], prime_base)
            bitwise_rep.append(base_rep)
            reversed_bitwise_rep.append(base_rep[::-1])
            sequence_bitwise.append(["0."] + reversed_bitwise_rep[i])
            sequence_decimal[i, 0] = self.prime_base_to_decimal(
                sequence_bitwise[i], prime_base
            )
        sequence_decimal = sequence_decimal.reshape(
            sequence_decimal.shape[0],
        )
        return sequence_decimal

    def selection_columns_preprocessing(self, data_input, xlabels, ylabels):
        """
        Pre-processing data for multiple output selection case.

        Args:
            data_input:     data supplied by user (dataframe or numpy array)
            xlabels:        list of input variables
            ylabels:        list of output variables
        """
        self.df_flag = True

        if isinstance(data_input, pd.DataFrame):
            if xlabels is None:
                xlabels = []
            if ylabels is None:
                ylabels = []
            set_of_labels = xlabels + ylabels

            if not set_of_labels:
                # Empty list assumes one output in final column - all other variables treated as inputs
                print(
                    "No column information provided. All except last column will be considered as x variables."
                )
                self.data = data_input.values
                self.data_headers = data_input.columns.values.tolist()
                self.x_data = self.data[:, :-1]
                self.data_headers_xvars = self.data_headers[:-1]  # Last column is y
            elif not set(set_of_labels).issubset(data_input.columns):
                raise IndexError(
                    'Invalid column names in "xlabels" or "ylabels": ',
                    [k for k in set_of_labels if k not in data_input.columns],
                )
            elif not ylabels:
                # When only xlabels are supplied, assume all non-inputs are outputs
                self.x_data = data_input.filter(xlabels).values
                self.data_headers_xvars = xlabels
                non_x_vars = [k for k in data_input.columns if k not in xlabels]
                y_data = data_input.filter(non_x_vars).values
                self.data_headers = self.data_headers_xvars + non_x_vars
                self.data = np.concatenate((self.x_data, y_data), axis=1)
            elif not xlabels:
                # When only ylabels are supplied, assume all non-outputs are inputs
                xvars = [k for k in data_input.columns if k not in ylabels]
                self.x_data = data_input.filter(xvars).values
                self.data_headers_xvars = xvars
                y_data = data_input.filter(ylabels).values
                self.data_headers = self.data_headers_xvars + ylabels
                self.data = np.concatenate((self.x_data, y_data), axis=1)
            else:
                # Only provided columns are considered - any column not in xlabels and ylabels is dropped completely.
                dropped_cols = [k for k in data_input.columns if k not in set_of_labels]
                if len(dropped_cols) > 0:
                    warn_str = "The following columns were dropped: " + str(
                        dropped_cols
                    )
                    warnings.warn(warn_str)
                self.x_data = data_input.filter(xlabels).values
                self.data_headers = set_of_labels
                self.data_headers_xvars = xlabels
                self.data = data_input.filter(set_of_labels).values

        elif isinstance(data_input, np.ndarray):
            self.df_flag = False

            if xlabels is None:
                xlabels = []
            if ylabels is None:
                ylabels = []
            set_of_labels = xlabels + ylabels

            if not set_of_labels:
                # Empty list assumes one output in final column
                print(
                    "No column information provided. All except last column will be considered as x variables."
                )
                self.data = data_input
                self.data_headers = [i for i in range(0, data_input.shape[1])]
                self.x_data = self.data[:, :-1]
                self.data_headers_xvars = self.data_headers[:-1]  # Last column is y
            elif not set(set_of_labels).issubset(range(0, data_input.shape[1])):
                raise IndexError(
                    'Invalid column names in "xlabels" or "ylabels": ',
                    [
                        k
                        for k in set_of_labels
                        if k not in range(0, data_input.shape[1])
                    ],
                )
            elif not ylabels:
                # Assumes all non-inputs are outputs
                self.x_data = data_input[:, xlabels]
                self.data_headers_xvars = xlabels
                non_x_vars = [
                    k for k in range(0, data_input.shape[1]) if k not in xlabels
                ]
                y_data = data_input[:, non_x_vars]
                self.data_headers = self.data_headers_xvars + non_x_vars
                self.data = np.concatenate((self.x_data, y_data), axis=1)
            elif not xlabels:
                # Assumes all non-outputs are inputs
                xvars = [k for k in range(0, data_input.shape[1]) if k not in ylabels]
                self.x_data = data_input[:, xvars]
                self.data_headers_xvars = xvars
                y_data = data_input[:, ylabels]
                self.data_headers = self.data_headers_xvars + ylabels
                self.data = np.concatenate((self.x_data, y_data), axis=1)
            else:
                # Only provided columns are considered
                dropped_cols = [
                    k for k in range(0, data_input.shape[1]) if k not in set_of_labels
                ]
                if len(dropped_cols) > 0:
                    warn_str = "The following columns were dropped: " + str(
                        dropped_cols
                    )
                    warnings.warn(warn_str)
                self.x_data = data_input[:, xlabels]
                self.data_headers = set_of_labels
                self.data_headers_xvars = xlabels
                self.data = data_input[:, set_of_labels]


class LatinHypercubeSampling(SamplingMethods):
    """
    A class that performs Latin Hypercube Sampling. The function returns LHS samples which have been selected randomly after sample space stratification.

    It should be noted that no minimax criterion has been used in this implementation, so the LHS samples selected will not have space-filling properties.

    To use: call class with inputs, and then run ``sample_points`` method.

    **Example:**

    .. code-block:: python

        # To select 10 LHS samples from "data"
        >>> b = rbf.LatinHypercubeSampling(data, 10, sampling_type="selection")
        >>> samples = b.sample_points()

    """

    def __init__(
        self,
        data_input,
        number_of_samples=None,
        sampling_type=None,
        xlabels=None,
        ylabels=None,
    ):
        """
        Initialization of **LatinHypercubeSampling** class. Two inputs are required.

        Args:
            data_input (NumPy Array, Pandas Dataframe or list) :  The input data set or range to be sampled.

                - When the aim is to select a set of samples from an existing dataset, the dataset must be a NumPy Array or a Pandas Dataframe and **sampling_type** option must be set to "selection". A single output variable (y) is assumed to be supplied in the last column if **xlabels** and **ylabels** are not supplied.
                - When the aim is to generate a set of samples from a data range, the dataset must be a list containing two lists of equal lengths which contain the variable bounds and **sampling_type** option must be set to "creation". It is assumed that the range contains no output variable information in this case.

            number_of_samples (int): The number of samples to be generated. Should be a positive integer less than or equal to the number of entries (rows) in **data_input**.
            sampling_type (str) : Option which determines whether the algorithm selects samples from an existing dataset ("selection") or attempts to generate sample from a supplied range ("creation"). Default is "creation".

        Keyword Args:
            xlabels (list): List of column names (if **data_input** is a dataframe) or column numbers (if **data_input** is an array) for the independent/input  variables.  Only used in "selection" mode. Default is None.
            ylabels (list): List of column names (if **data_input** is a dataframe) or column numbers (if **data_input** is an array) for the dependent/output variables. Only used in "selection" mode. Default is None.

        Returns:
            **self** function containing the input information

        Raises:
            ValueError: The input data (**data_input**) is the wrong type.

            IndexError: When invalid column names are supplied in **xlabels** or **ylabels**

            Exception: When **number_of_samples** is invalid (not an integer, too large, zero, or negative)

        """

        if sampling_type is None:
            sampling_type = "creation"
            self.sampling_type = sampling_type
            print("Creation-type sampling will be used.")
        elif not isinstance(sampling_type, string_types):
            raise Exception("Invalid sampling type entry. Must be of type <str>.")
        elif (sampling_type.lower() == "creation") or (
            sampling_type.lower() == "selection"
        ):
            sampling_type = sampling_type.lower()
            self.sampling_type = sampling_type
        else:
            raise Exception(
                'Invalid sampling type requirement entered. Enter "creation" for sampling from a range or "selection" for selecting samples from a dataset.'
            )
        print("Sampling type: ", self.sampling_type, "\n")

        if self.sampling_type == "selection":
            if isinstance(data_input, (pd.DataFrame, np.ndarray)):
                self.selection_columns_preprocessing(data_input, xlabels, ylabels)
            else:
                raise ValueError(
                    'Pandas dataframe or numpy array required for sampling_type "selection."'
                )

            # Catch potential errors in number_of_samples
            if number_of_samples is None:
                print(
                    "\nNo entry for number of samples to be generated. The default value of 5 will be used."
                )
                number_of_samples = 5
            elif number_of_samples > self.data.shape[0]:
                raise Exception(
                    "LHS sample size cannot be greater than number of samples in the input data set"
                )
            elif not isinstance(number_of_samples, int):
                raise Exception("number_of_samples must be an integer.")
            elif number_of_samples <= 0:
                raise Exception("number_of_samples must a positive, non-zero integer.")
            self.number_of_samples = number_of_samples

        elif self.sampling_type == "creation":
            if not isinstance(data_input, list):
                raise ValueError(
                    'List entry of two elements expected for sampling_type "creation."'
                )
            elif len(data_input) != 2:
                raise Exception("data_input must contain two lists of equal lengths.")
            elif not isinstance(data_input[0], list) or not isinstance(
                data_input[1], list
            ):
                raise Exception("data_input must contain two lists of equal lengths.")
            elif len(data_input[0]) != len(data_input[1]):
                raise Exception("data_input must contain two lists of equal lengths.")
            elif data_input[0] == data_input[1]:
                raise Exception("Invalid entry: both lists are equal.")
            else:
                bounds_array = np.zeros(
                    (
                        2,
                        len(data_input[0]),
                    )
                )
                bounds_array[0, :] = np.array(data_input[0])
                bounds_array[1, :] = np.array(data_input[1])
                data_headers = []
            self.data = bounds_array
            self.data_headers = data_headers

            # Catch potential errors in number_of_samples
            if number_of_samples is None:
                print(
                    "\nNo entry for number of samples to be generated. The default value of 5 will be used."
                )
                number_of_samples = 5
            elif not isinstance(number_of_samples, int):
                raise Exception("number_of_samples must be an integer.")
            elif number_of_samples <= 0:
                raise Exception("number_of_samples must a positive, non-zero integer.")
            self.number_of_samples = number_of_samples
            self.x_data = bounds_array  # Only x data will be present in this case

    def variable_sample_creation(self, variable_min, variable_max):
        """

        Function that generates the required number of sample points for a given variable within a specified range using stratification.
        The function divides the variable sample space into self.number_of_samples equal strata and generates a single random sample from each strata based on its lower and upper bound.

        Args:
            self
            variable_min(float): The lower bound of the sample space region. Should be a single number.
            variable_max(float): The upper bound of the sample space region. Should be a single number.

        Returns:
            var_samples(NumPy Array): A numpy array of size (self.number_of_samples x 1) containing the randomly generated points from each strata
        """

        strata_size = 1 / self.number_of_samples
        var_samples = np.zeros((self.number_of_samples, 1))
        for i in range(self.number_of_samples):
            strata_lb = i * strata_size
            sample_point = strata_lb + (np.random.rand() * strata_size)
            var_samples[i, 0] = (
                sample_point * (variable_max - variable_min)
            ) + variable_min
        return var_samples

    def lhs_points_generation(self):
        """
        Generate points within each strata for each variable based on stratification. When invoked, it:
        1. Determines the mimumum and maximum value for each feature (column),
        2. Calls the variable_sample_creation function on each feature, passing in its mimmum and maximum
        3. Returns an array containing the points selected in each strata of each column

        Returns:
            sample_points_vector(NumPy Array): Array containing the columns of the random samples generated in each strata.
        """

        ns, nf = np.shape(self.x_data)
        sample_points_vector = np.zeros(
            (self.number_of_samples, nf)
        )  # Array containing points in each interval for each variable
        for i in range(nf):
            variable_min = 0  # np.min(self.x_data[:, i])
            variable_max = 1  # np.max(self.x_data[:, i])
            var_samples = self.variable_sample_creation(
                variable_min, variable_max
            )  # Data generation step
            sample_points_vector[:, i] = var_samples[:, 0]
        return sample_points_vector

    @staticmethod
    def random_shuffling(vector_of_points):
        """
        This function carries out random shuffling of column data to generate samples.
        Data in each of the columns  in the input array is shuffled separately, meaning that the rows of the resultant array will contain random samples from the sample space.

        Args:
            vector_of_points(NumPy Array): Array containing ordered points generated from stratification. Should usually be the output of the lhs_points_generation function. Each column self.number_of_samples elements.

        Returns:
            vector_of_points(NumPy Array): 2-D array containing the shuffled data. Should contain number_of_sample rows, with each row representing a potential random sample from within the sample space.

        """

        _, nf = np.shape(vector_of_points)
        for i in range(0, nf):
            z_col = vector_of_points[:, i]
            np.random.shuffle(z_col)
            vector_of_points[:, i] = z_col
        return vector_of_points

    def sample_points(self):
        """
        ``sample_points`` generates or selects Latin Hypercube samples from an input dataset or data range. When called, it:

            1. generates samples points from stratified regions by calling the ``lhs_points_generation``,
            2. generates potential sample points by random shuffling, and
            3. when a dataset is provided, selects the closest available samples to the theoretical sample points from within the input data.

        Returns:
            NumPy Array or Pandas Dataframe:     A numpy array or Pandas dataframe containing **number_of_samples** points selected or generated by LHS.

        """

        vector_of_points = (
            self.lhs_points_generation()
        )  # Assumes [X, Y] data is supplied.
        generated_sample_points = self.random_shuffling(vector_of_points)
        unique_sample_points = self.sample_point_selection(
            self.data, generated_sample_points, self.sampling_type
        )

        if len(self.data_headers) > 0 and self.df_flag:
            unique_sample_points = pd.DataFrame(
                unique_sample_points, columns=self.data_headers
            )
        return unique_sample_points


class UniformSampling(SamplingMethods):
    """
    A class that performs Uniform Sampling. Depending on the settings, the algorithm either returns samples from an input dataset  which have been selected using Euclidean distance minimization after the uniform samples have been generated,
    or returns samples from a supplied data range.

    Full-factorial samples are based on dividing the space of each variable randomly and then generating all possible variable combinations.

    - The number of points to be sampled per variable needs to be specified in a list.

    To use: call class with inputs, and then ``sample_points`` function

    **Example:**

    .. code-block:: python

        # To select 50 samples on a (10 x 5) grid in a 2D space:
        >>> b = rbf.UniformSampling(data, [10, 5], sampling_type="selection")
        >>> samples = b.sample_points()

    """

    def __init__(
        self,
        data_input,
        list_of_samples_per_variable,
        sampling_type=None,
        xlabels=None,
        ylabels=None,
        edges=None,
    ):
        """
        Initialization of UniformSampling class. Three inputs are required.

        Args:
            data_input (NumPy Array, Pandas Dataframe or list) :  The input data set or range to be sampled.

                - When the aim is to select a set of samples from an existing dataset, the dataset must be a NumPy Array or a Pandas Dataframe and **sampling_type** option must be set to "selection". A single output variable (y) is assumed to be supplied in the last column if **xlabels** and **ylabels** are not supplied.
                - When the aim is to generate a set of samples from a data range, the dataset must be a list containing two lists of equal lengths which contain the variable bounds and **sampling_type** option must be set to "creation". It is assumed that the range contains no output variable information in this case.

            list_of_samples_per_variable (list): The list containing the number of subdivisions for each variable. Each dimension (variable) must be represented by a positive integer variable greater than 1.
            sampling_type (str) : Option which determines whether the algorithm selects samples from an existing dataset ("selection") or attempts to generate sample from a supplied range ("creation"). Default is "creation".

        Keyword Args:
            xlabels (list): List of column names (if **data_input** is a dataframe) or column numbers (if **data_input** is an array) for the independent/input  variables.  Only used in "selection" mode. Default is None.
            ylabels (list): List of column names (if **data_input** is a dataframe) or column numbers (if **data_input** is an array) for the dependent/output variables. Only used in "selection" mode. Default is None.
            edges (bool): Boolean variable representing how the points should be selected. A value of True (default) indicates the points should be equally spaced edge to edge, otherwise they will be in the centres of the bins filling the unit cube

        Returns:
            **self** function containing the input information

        Raises:
            ValueError: The **data_input** is the wrong type

            ValueError: When **list_of_samples_per_variable** is of the wrong length, is not a list or contains elements other than integers

            IndexError: When invalid column names are supplied in **xlabels** or **ylabels**

            Exception: When **edges** entry is not Boolean

        """
        if sampling_type is None:
            sampling_type = "creation"
            self.sampling_type = sampling_type
            print("Creation-type sampling will be used.")
        elif not isinstance(sampling_type, string_types):
            raise Exception("Invalid sampling type entry. Must be of type <str>.")
        elif (sampling_type.lower() == "creation") or (
            sampling_type.lower() == "selection"
        ):
            sampling_type = sampling_type.lower()
            self.sampling_type = sampling_type
        else:
            raise Exception(
                'Invalid sampling type requirement entered. Enter "creation" for sampling from a range or "selection" for selecting samples from a dataset.'
            )
        print("Sampling type: ", self.sampling_type, "\n")

        if self.sampling_type == "selection":
            if isinstance(data_input, (pd.DataFrame, np.ndarray)):
                self.selection_columns_preprocessing(data_input, xlabels, ylabels)
            else:
                raise ValueError(
                    'Pandas dataframe or numpy array required for sampling_type "selection."'
                )

        elif self.sampling_type == "creation":
            if not isinstance(data_input, list):
                raise ValueError(
                    'List entry of two elements expected for sampling_type "creation."'
                )
            elif len(data_input) != 2:
                raise Exception("data_input must contain two lists of equal lengths.")
            elif not isinstance(data_input[0], list) or not isinstance(
                data_input[1], list
            ):
                raise Exception("data_input must contain two lists of equal lengths.")
            elif len(data_input[0]) != len(data_input[1]):
                raise Exception("data_input must contain two lists of equal lengths.")
            elif data_input[0] == data_input[1]:
                raise Exception("Invalid entry: both lists are equal.")
            else:
                bounds_array = np.zeros(
                    (
                        2,
                        len(data_input[0]),
                    )
                )
                bounds_array[0, :] = np.array(data_input[0])
                bounds_array[1, :] = np.array(data_input[1])
                data_headers = []
            self.data = bounds_array
            self.x_data = bounds_array
            self.data_headers = data_headers

        if edges is None:
            edges = True
            self.edge = edges
        elif not isinstance(edges, bool):
            raise Exception('Invalid "edges" entry. Must be boolean')
        elif (edges is True) or (edges is False):
            self.edge = edges

        # Check that list_of_samples_per_variable is a list, list length is correct, all dimensions greater than 1 and all list values are integers
        if not isinstance(list_of_samples_per_variable, list):
            raise TypeError("list_of_samples_per_variable: list required.")
        if len(list_of_samples_per_variable) != self.x_data.shape[1]:
            raise ValueError(
                "Length of list_of_samples_per_variable must equal the number of variables."
            )
        if min(list_of_samples_per_variable) < 2:
            raise ValueError(
                "All variables must have at least two points per dimension"
            )
        if all(isinstance(q, int) for q in list_of_samples_per_variable) is False:
            raise TypeError("All values in list must be integers")

        self.dim_vector = list_of_samples_per_variable
        self.number_of_samples = int(np.prod(self.dim_vector))

        if (
            self.sampling_type == "selection"
            and self.number_of_samples > self.data.shape[0]
        ):
            raise Exception(
                "Sample size cannot be greater than number of samples in the input data set"
            )

    def sample_points(self):
        """
        ``sample_points`` generates or selects full-factorial designs from an input dataset or data range.

        Returns:
            NumPy Array or Pandas Dataframe:     A numpy array or Pandas dataframe containing the sample points generated or selected by full-factorial sampling.

        """

        points_spread = []
        if self.edge is True:
            for i in self.dim_vector:
                variable_spread = np.arange(i) / (i - 1)
                points_spread.append(variable_spread)
        elif self.edge is False:
            for i in self.dim_vector:
                variable_spread = np.arange(i + 1) / i
                shifted_points = [
                    (variable_spread[i] + variable_spread[i - 1]) / 2
                    for i in range(1, len(variable_spread))
                ]
                points_spread.append(shifted_points)
        samples_list = list(itertools.product(*points_spread))
        samples_array = np.asarray(samples_list)
        unique_sample_points = self.sample_point_selection(
            self.data, samples_array, self.sampling_type
        )
        if len(self.data_headers) > 0 and self.df_flag:
            unique_sample_points = pd.DataFrame(
                unique_sample_points, columns=self.data_headers
            )
        return unique_sample_points


class HaltonSampling(SamplingMethods):
    """
    A class that performs Halton Sampling.

    Halton samples are based on the reversing/flipping the base conversion of numbers using primes.

    To generate n samples in a :math:`p`-dimensional space, the first :math:`p` prime numbers are used to generate the samples.

    Note:
        Use of this method is limited to use in low-dimensionality problems (less than 10 variables). At higher dimensions, the performance of the sampling method has been shown to degrade.

    To use: call class with inputs, and then ``sample_points`` function.

    **Example:**

    .. code-block:: python

        # For the first 10 Halton samples in a 2-D space:
        >>> b = rbf.HaltonSampling(data, 10, sampling_type="selection")
        >>> samples = b.sample_points()

    """

    def __init__(
        self,
        data_input,
        number_of_samples=None,
        sampling_type=None,
        xlabels=None,
        ylabels=None,
    ):
        """

        Initialization of **HaltonSampling** class. Two inputs are required.

        Args:
            data_input (NumPy Array, Pandas Dataframe or list) : The input data set or range to be sampled.

                - When the aim is to select a set of samples from an existing dataset, the dataset must be a NumPy Array or a Pandas Dataframe and **sampling_type** option must be set to "selection". A single output variable (y) is assumed to be supplied in the last column if **xlabels** and **ylabels** are not supplied.
                - When the aim is to generate a set of samples from a data range, the dataset must be a list containing two lists of equal lengths which contain the variable bounds and **sampling_type** option must be set to "creation". It is assumed that the range contains no output variable information  in this case.

            number_of_samples(int): The number of samples to be generated. Should be a positive integer less than or equal to the number of entries (rows) in **data_input**.
            sampling_type(str) : Option which determines whether the algorithm selects samples from an existing dataset ("selection") or attempts to generate sample from a supplied range ("creation"). Default is "creation".

        Keyword Args:
            xlabels (list): List of column names (if **data_input** is a dataframe) or column numbers (if **data_input** is an array) for the independent/input  variables.  Only used in "selection" mode. Default is None.
            ylabels (list): List of column names (if **data_input** is a dataframe) or column numbers (if **data_input** is an array) for the dependent/output variables. Only used in "selection" mode. Default is None.

        Returns:
            **self** function containing the input information.

        Raises:
            ValueError: The **data_input** is the wrong type.

            IndexError: When invalid column names are supplied in **xlabels** or **ylabels**

            Exception: When the **number_of_samples** is invalid (not an integer, too large, zero or negative.)

        """
        if sampling_type is None:
            sampling_type = "creation"
            self.sampling_type = sampling_type
            print("Creation-type sampling will be used.")
        elif not isinstance(sampling_type, string_types):
            raise Exception("Invalid sampling type entry. Must be of type <str>.")
        elif (sampling_type.lower() == "creation") or (
            sampling_type.lower() == "selection"
        ):
            sampling_type = sampling_type.lower()
            self.sampling_type = sampling_type
        else:
            raise Exception(
                'Invalid sampling type requirement entered. Enter "creation" for sampling from a range or "selection" for selecting samples from a dataset.'
            )
        print("Sampling type: ", self.sampling_type, "\n")

        if self.sampling_type == "selection":
            if isinstance(data_input, (pd.DataFrame, np.ndarray)):
                self.selection_columns_preprocessing(data_input, xlabels, ylabels)
            else:
                raise ValueError(
                    'Pandas dataframe or numpy array required for sampling_type "selection."'
                )

            # Catch potential errors in number_of_samples
            if number_of_samples is None:
                print(
                    "\nNo entry for number of samples to be generated. The default value of 5 will be used."
                )
                number_of_samples = 5
            elif number_of_samples > self.data.shape[0]:
                raise Exception(
                    "Sample size cannot be greater than number of samples in the input data set"
                )
            elif not isinstance(number_of_samples, int):
                raise Exception("number_of_samples must be an integer.")
            elif number_of_samples <= 0:
                raise Exception("number_of_samples must a positive, non-zero integer.")
            self.number_of_samples = number_of_samples

        elif self.sampling_type == "creation":
            if not isinstance(data_input, list):
                raise ValueError(
                    'List entry of two elements expected for sampling_type "creation."'
                )
            elif len(data_input) != 2:
                raise Exception("data_input must contain two lists of equal lengths.")
            elif not isinstance(data_input[0], list) or not isinstance(
                data_input[1], list
            ):
                raise Exception("data_input must contain two lists of equal lengths.")
            elif len(data_input[0]) != len(data_input[1]):
                raise Exception("data_input must contain two lists of equal lengths.")
            elif data_input[0] == data_input[1]:
                raise Exception("Invalid entry: both lists are equal.")
            else:
                bounds_array = np.zeros(
                    (
                        2,
                        len(data_input[0]),
                    )
                )
                bounds_array[0, :] = np.array(data_input[0])
                bounds_array[1, :] = np.array(data_input[1])
                data_headers = []
            self.data = bounds_array
            self.data_headers = data_headers

            # Catch potential errors in number_of_samples
            if number_of_samples is None:
                print(
                    "\nNo entry for number of samples to be generated. The default value of 5 will be used."
                )
                number_of_samples = 5
            elif not isinstance(number_of_samples, int):
                raise Exception("number_of_samples must be an integer.")
            elif number_of_samples <= 0:
                raise Exception("number_of_samples must a positive, non-zero integer.")
            self.number_of_samples = number_of_samples
            self.x_data = bounds_array  # Only x data will be present in this case

        if self.x_data.shape[1] > 10:
            raise Exception(
                "Dimensionality problem: This method is not available for problems with dimensionality > 10: the performance of the method degrades substantially at higher dimensions"
            )

    def sample_points(self):
        """
        The ``sample_points`` method generates the Halton samples. The steps followed here are:

            1. Determine the number of features in the input data.
            2. Generate the list of primes to be considered by calling ``prime_number_generator`` from the sampling superclass.
            3. Create the first **number_of_samples** elements of the Halton sequence for each prime.
            4. Create the Halton samples by combining the corresponding elements of the Halton sequences for each prime.
            5. When in "selection" mode, determine the closest corresponding point in the input dataset using Euclidean distance minimization. This is done by calling the ``nearest_neighbours`` method in the sampling superclass.

        Returns:
            NumPy Array or Pandas Dataframe:     A numpy array or Pandas dataframe containing **number_of_samples** Halton sample points.

        """
        no_features = self.x_data.shape[1]
        # Generate list of no_features prime numbers
        prime_list = self.prime_number_generator(no_features)
        sample_points = np.zeros((self.number_of_samples, no_features))
        for i in range(0, no_features):
            sample_points[:, i] = self.data_sequencing(
                self.number_of_samples, prime_list[i]
            )
        # Scale input data, then find data points closest in sample space. Unscale before returning points
        unique_sample_points = self.sample_point_selection(
            self.data, sample_points, self.sampling_type
        )
        if len(self.data_headers) > 0 and self.df_flag:
            unique_sample_points = pd.DataFrame(
                unique_sample_points, columns=self.data_headers
            )
        return unique_sample_points


class HammersleySampling(SamplingMethods):
    """
    A class that performs Hammersley Sampling.

    Hammersley samples are generated in a similar way to Halton samples - based on the reversing/flipping the base conversion of numbers using primes.

    To generate :math:`n` samples in a :math:`p`-dimensional space, the first :math:`\\left(p-1\\right)` prime numbers are used to generate the samples. The first dimension is obtained by uniformly dividing the region into **no_samples points**.

    Note:
        Use of this method is limited to use in low-dimensionality problems (less than 10 variables). At higher dimensionalities, the performance of the sampling method has been shown to degrade.

    To use: call class with inputs, and then ``sample_points`` function.

    **Example:**

    .. code-block:: python

        # For the first 10 Hammersley samples in a 2-D space:
        >>> b = rbf.HammersleySampling(data, 10, sampling_type="selection")
        >>> samples = b.sample_points()

    """

    def __init__(
        self,
        data_input,
        number_of_samples=None,
        sampling_type=None,
        xlabels=None,
        ylabels=None,
    ):
        """
        Initialization of **HammersleySampling** class. Two inputs are required.

        Args:
            data_input (NumPy Array, Pandas Dataframe or list): The input data set or range to be sampled.

               - When the aim is to select a set of samples from an existing dataset, the dataset must be a NumPy Array or a Pandas Dataframe and **sampling_type** option must be set to "selection". A single output variable (y) is assumed to be supplied in the last column if **xlabels** and **ylabels** are not supplied.
               - When the aim is to generate a set of samples from a data range, the dataset must be a list containing two lists of equal lengths which contain the variable bounds and **sampling_type** option must be set to "creation". It is assumed that the range contains no output variable information  in this case.

            number_of_samples(int): The number of samples to be generated. Should be a positive integer less than or equal to the number of entries (rows) in **data_input**.
            sampling_type(str) : Option which determines whether the algorithm selects samples from an existing dataset ("selection") or attempts to generate sample from a supplied range ("creation"). Default is "creation".

        Keyword Args:
            xlabels (list): List of column names (if **data_input** is a dataframe) or column numbers (if **data_input** is an array) for the independent/input  variables.  Only used in "selection" mode. Default is None.
            ylabels (list): List of column names (if **data_input** is a dataframe) or column numbers (if **data_input** is an array) for the dependent/output variables. Only used in "selection" mode. Default is None.

            Returns:
                **self** function containing the input information.

            Raises:
                ValueError: When **data_input** is the wrong type.

                IndexError: When invalid column names are supplied in **xlabels** or **ylabels**

                Exception: When the **number_of_samples** is invalid (not an integer, too large, zero, negative)

        """
        if sampling_type is None:
            sampling_type = "creation"
            self.sampling_type = sampling_type
            print("Creation-type sampling will be used.")
        elif not isinstance(sampling_type, string_types):
            raise Exception("Invalid sampling type entry. Must be of type <str>.")
        elif (sampling_type.lower() == "creation") or (
            sampling_type.lower() == "selection"
        ):
            sampling_type = sampling_type.lower()
            self.sampling_type = sampling_type
        else:
            raise Exception(
                'Invalid sampling type requirement entered. Enter "creation" for sampling from a range or "selection" for selecting samples from a dataset.'
            )
        print("Sampling type: ", self.sampling_type, "\n")

        if self.sampling_type == "selection":

            if isinstance(data_input, (pd.DataFrame, np.ndarray)):
                self.selection_columns_preprocessing(data_input, xlabels, ylabels)
            else:
                raise ValueError(
                    'Pandas dataframe or numpy array required for sampling_type "selection."'
                )

            # Catch potential errors in number_of_samples
            if number_of_samples is None:
                print(
                    "\nNo entry for number of samples to be generated. The default value of 5 will be used."
                )
                number_of_samples = 5
            elif number_of_samples > self.data.shape[0]:
                raise Exception(
                    "Sample size cannot be greater than number of samples in the input data set"
                )
            elif not isinstance(number_of_samples, int):
                raise Exception("number_of_samples must be an integer.")
            elif number_of_samples <= 0:
                raise Exception("number_of_samples must a positive, non-zero integer.")
            self.number_of_samples = number_of_samples

        elif self.sampling_type == "creation":
            if not isinstance(data_input, list):
                raise ValueError(
                    'List entry of two elements expected for sampling_type "creation."'
                )
            elif len(data_input) != 2:
                raise Exception("data_input must contain two lists of equal lengths.")
            elif not isinstance(data_input[0], list) or not isinstance(
                data_input[1], list
            ):
                raise Exception("data_input must contain two lists of equal lengths.")
            elif len(data_input[0]) != len(data_input[1]):
                raise Exception("data_input must contain two lists of equal lengths.")
            elif data_input[0] == data_input[1]:
                raise Exception("Invalid entry: both lists are equal.")
            else:
                bounds_array = np.zeros(
                    (
                        2,
                        len(data_input[0]),
                    )
                )
                bounds_array[0, :] = np.array(data_input[0])
                bounds_array[1, :] = np.array(data_input[1])
                data_headers = []
            self.data = bounds_array
            self.data_headers = data_headers

            # Catch potential errors in number_of_samples
            if number_of_samples is None:
                print(
                    "\nNo entry for number of samples to be generated. The default value of 5 will be used."
                )
                number_of_samples = 5
            elif not isinstance(number_of_samples, int):
                raise Exception("number_of_samples must be an integer.")
            elif number_of_samples <= 0:
                raise Exception("number_of_samples must a positive, non-zero integer.")
            self.number_of_samples = number_of_samples
            self.x_data = bounds_array  # Only x data will be present in this case

        if self.x_data.shape[1] > 10:
            raise Exception(
                "Dimensionality problem: This method is not available for problems with dimensionality > 10: the performance of the method degrades substantially at higher dimensions"
            )

    def sample_points(self):
        """
        The **sampling_type** method generates the Hammersley sample points. The steps followed here are:

            1. Determine the number of features :math:`n_{f}` in the input data.
            2. Generate the list of :math:`\\left(n_{f}-1\\right)` primes to be considered by calling prime_number_generator.
            3. Divide the space [0,**number_of_samples**-1] into **number_of_samples** places to obtain the first dimension for the Hammersley sequence.
            4. For the other :math:`\\left(n_{f}-1\\right)` dimensions, create first **number_of_samples** elements of the Hammersley sequence for each of the :math:`\\left(n_{f}-1\\right)` primes.
            5. Create the Hammersley samples by combining the corresponding elements of the Hammersley sequences created in steps 3 and 4
            6. When in "selection" mode, determine the closest corresponding point in the input dataset using Euclidean distance minimization. This is done by calling the ``nearest_neighbours`` method in the sampling superclass.

        Returns:
            NumPy Array or Pandas Dataframe:     A numpy array or Pandas dataframe containing **number_of_samples** Hammersley sample points.

        """
        no_features = self.x_data.shape[1]
        if no_features == 1:
            prime_list = []
        else:
            prime_list = self.prime_number_generator(no_features - 1)
        sample_points = np.zeros((self.number_of_samples, no_features))
        sample_points[:, 0] = (
            np.arange(0, self.number_of_samples)
        ) / self.number_of_samples
        for i in range(0, len(prime_list)):
            sample_points[:, i + 1] = self.data_sequencing(
                self.number_of_samples, prime_list[i]
            )

        unique_sample_points = self.sample_point_selection(
            self.data, sample_points, self.sampling_type
        )
        if len(self.data_headers) > 0 and self.df_flag:
            unique_sample_points = pd.DataFrame(
                unique_sample_points, columns=self.data_headers
            )
        return unique_sample_points


class CVTSampling(SamplingMethods):
    """
    A class that constructs Centroidal Voronoi Tessellation (CVT) samples.

    CVT sampling is based on the generation of samples in which the generators of the Voronoi tessellations and the mass centroids coincide.

    To use: call class with inputs, and then ``sample_points`` function.

    **Example:**

    .. code-block:: python

        # For the first 10 CVT samples in a 2-D space:
        >>> b = rbf.CVTSampling(data_bounds, 10, tolerance = 1e-5, sampling_type="creation")
        >>> samples = b.sample_points()

    """

    def __init__(
        self,
        data_input,
        number_of_samples=None,
        tolerance=None,
        sampling_type=None,
        xlabels=None,
        ylabels=None,
    ):
        """
        Initialization of CVTSampling class. Two inputs are required, while an optional option to control the solution accuracy may be specified.

        Args:
            data_input (NumPy Array, Pandas Dataframe or list): The input data set or range to be sampled.
                - When the aim is to select a set of samples from an existing dataset, the dataset must be a NumPy Array or a Pandas Dataframe and **sampling_type** option must be set to "selection". A single output variable (y) is assumed to be supplied in the last column if **xlabels** and **ylabels** are not supplied.
                - When the aim is to generate a set of samples from a data range, the dataset must be a list containing two lists of equal lengths which contain the variable bounds and **sampling_type** option must be set to "creation". It is assumed that the range contains no output variable information  in this case.


            number_of_samples(int): The number of samples to be generated. Should be a positive integer less than or equal to the number of entries (rows) in **data_input**.
            sampling_type(str) : Option which determines whether the algorithm selects samples from an existing dataset ("selection") or attempts to generate sample from a supplied range ("creation"). Default is "creation".

        Keyword Args:
            xlabels (list): List of column names (if **data_input** is a dataframe) or column numbers (if **data_input** is an array) for the independent/input  variables.  Only used in "selection" mode. Default is None.
            ylabels (list): List of column names (if **data_input** is a dataframe) or column numbers (if **data_input** is an array) for the dependent/output variables. Only used in "selection" mode. Default is None.
            tolerance(float): Maximum allowable Euclidean distance between centres from consectutive iterations of the algorithm. Termination condition for algorithm.

                - The smaller the value of tolerance, the better the solution but the longer the algorithm requires to converge. Default value is :math:`10^{-7}`.

        Returns:
                **self** function containing the input information.


        Raises:
                ValueError: When **data_input** is the wrong type.

                IndexError: When invalid column names are supplied in **xlabels** or **ylabels**

                Exception: When the **number_of_samples** is invalid (not an integer, too large, zero, negative)

                Exception: When the tolerance specified is too loose (tolerance > 0.1) or invalid

                warnings.warn: when the tolerance specified by the user is too tight (tolerance < :math:`10^{-9}`)

        """
        if sampling_type is None:
            sampling_type = "creation"
            self.sampling_type = sampling_type
            print("Creation-type sampling will be used.")
        elif not isinstance(sampling_type, string_types):
            raise Exception("Invalid sampling type entry. Must be of type <str>.")
        elif (sampling_type.lower() == "creation") or (
            sampling_type.lower() == "selection"
        ):
            sampling_type = sampling_type.lower()
            self.sampling_type = sampling_type
        else:
            raise Exception(
                'Invalid sampling type requirement entered. Enter "creation" for sampling from a range or "selection" for selecting samples from a dataset.'
            )
        print("Sampling type: ", self.sampling_type, "\n")

        if self.sampling_type == "selection":
            if isinstance(data_input, (pd.DataFrame, np.ndarray)):
                self.selection_columns_preprocessing(data_input, xlabels, ylabels)
            else:
                raise ValueError(
                    'Pandas dataframe or numpy array required for sampling_type "selection."'
                )

            # Catch potential errors in number_of_samples
            if number_of_samples is None:
                print(
                    "\nNo entry for number of samples to be generated. The default value of 5 will be used."
                )
                number_of_samples = 5
            elif number_of_samples > self.data.shape[0]:
                raise Exception(
                    "CVT sample size cannot be greater than number of samples in the input data set"
                )
            elif not isinstance(number_of_samples, int):
                raise Exception("number_of_samples must be an integer.")
            elif number_of_samples <= 0:
                raise Exception("number_of_samples must a positive, non-zero integer.")
            self.number_of_centres = number_of_samples

        elif self.sampling_type == "creation":
            if not isinstance(data_input, list):
                raise ValueError(
                    'List entry of two elements expected for sampling_type "creation."'
                )
            elif len(data_input) != 2:
                raise Exception("data_input must contain two lists of equal lengths.")
            elif not isinstance(data_input[0], list) or not isinstance(
                data_input[1], list
            ):
                raise Exception("data_input must contain two lists of equal lengths.")
            elif len(data_input[0]) != len(data_input[1]):
                raise Exception("data_input must contain two lists of equal lengths.")
            elif data_input[0] == data_input[1]:
                raise Exception("Invalid entry: both lists are equal.")
            else:
                bounds_array = np.zeros(
                    (
                        2,
                        len(data_input[0]),
                    )
                )
                bounds_array[0, :] = np.array(data_input[0])
                bounds_array[1, :] = np.array(data_input[1])
                data_headers = []
            self.data = bounds_array
            self.data_headers = data_headers

            # Catch potential errors in number_of_samples
            if number_of_samples is None:
                print(
                    "\nNo entry for number of samples to be generated. The default value of 5 will be used."
                )
                number_of_samples = 5
            elif not isinstance(number_of_samples, int):
                raise Exception("number_of_samples must be an integer.")
            elif number_of_samples <= 0:
                raise Exception("number_of_samples must a positive, non-zero integer.")
            self.number_of_centres = number_of_samples

            x_data = bounds_array  # Only x data will be present in this case
            if x_data.ndim == 1:
                x_data = x_data.reshape(len(x_data), 1)
            self.x_data = x_data
            # self.y_data = []

        if tolerance is None:
            tolerance = 1e-7
        elif tolerance > 0.1:
            raise Exception("Tolerance must be less than 0.1 to achieve good results")
        elif tolerance < 1e-9:
            warnings.warn(
                "Tolerance too tight. CVT algorithm may take long time to converge."
            )
        elif (tolerance < 0.1) and (tolerance > 1e-9):
            tolerance = tolerance
        else:
            raise Exception("Invalid tolerance input")
        self.eps = tolerance

    @staticmethod
    def random_sample_selection(no_samples, no_features):
        """
        Function generates a the required number of samples (no_samples) within an no_features-dimensional space.
        This is achieved by generating an m x n 2-D array using numpy's random.rand function, where

            - m = number of training samples to be generated, and'
            - n = number of design features/variables (dimensionality of the problem).

        Args:
            no_samples(int): The number of samples to be generated.
            no_features(int): Number of design features/variables in the input data.

        Returns:
            random_points(NumPy Array): 2-D array of size no_samples x no_features generated from a uniform distribution.

        Example: Generate three samples for a two-dimensional problem
            >>  rbf.CVTSampling.random_sample_selection(3, 2)
            >> array([[0.03149075, 0.70566624], [0.48319597, 0.03810093], [0.19962214, 0.57641408]])

        """
        random_points = np.random.rand(no_samples, no_features)
        return random_points

    @staticmethod
    def eucl_distance(u, v):
        """
        The function eucl_distance(u,v) calculates Euclidean distance between two points or arrays u and v.

        Args:
            u, v (NumPy Array): Two points or arrays with the same number of features (same second dimension)

        Returns:
            euc_d(NumPy Array): Array of size (u.shape[0] x 1) containing Euclidean distances.

        """
        d = u - v
        d_sq = d**2
        euc_d = np.sqrt(np.sum(d_sq, axis=1))
        return euc_d

    @staticmethod
    def create_centres(
        initial_centres, current_random_points, current_centres, counter
    ):
        """
        The function create_centres generates new mass centroids for the design space based on McQueen's method.
        The mass centroids are created based on the previous mass centroids and the mean of random data sampling the design space.

            Args:
                initial_centres(NumPy Array): A 2-D array containing the current mass centroids, size no_samples x no_features.
                current_random_points(NumPy Array): A 2-D array containing several points generated randomly from within the design space.
                current_centres(NumPy Array): Array containing the index number of the closest mass centroid of each point in current_random_points, representing its class.
                counter(int): current iteration number

            Returns:
                centres(NumPy Array): A 2-D array containing the new mass centroids, size no_samples x no_features.

        The steps carried out in the function at each iteration are:
        (1) Classify the current random points in current_random_points based on their centres
        (2) Evaluate the mean of the random points in each class
        (3) Create the new centres as the weighted average of the current centres (initial_centres) and the mean data calculated in the second step. The weighting is done based on the number of iterations (counter).

        """
        centres = np.zeros((initial_centres.shape[0], initial_centres.shape[1]))
        current_centres = current_centres.reshape(current_centres.shape[0], 1)
        for i in range(0, initial_centres.shape[0]):
            data_matrix = current_random_points[current_centres[:, 0] == i]
            m_prime, n_prime = data_matrix.shape
            if m_prime == 0:
                centres[i, :] = np.mean(initial_centres, axis=0)
            else:
                centres[i, :] = np.mean(data_matrix, axis=0)

        # Weighted average based on previous number of iterations
        centres = ((counter * initial_centres) + centres) / (counter + 1)
        return centres

    def sample_points(self):
        """
        The ``sample_points`` method determines the best/optimal centre points (centroids) for a data set based on the minimization of the total distance between points and centres.

        Procedure based on McQueen's algorithm: iteratively minimize distance, and re-position centroids.
        Centre re-calculation done as the mean of each data cluster around each centre.

        Returns:
            NumPy Array or Pandas Dataframe:     A numpy array or Pandas dataframe containing the final **number_of_samples** centroids obtained by the CVT algorithm.

        """
        _, n = self.x_data.shape
        size_multiple = 1000
        initial_centres = self.random_sample_selection(self.number_of_centres, n)
        # Iterative optimization process
        cost_old = 0
        cost_new = 0
        cost_change = float("Inf")
        counter = 1
        while (cost_change > self.eps) and (counter <= 1000):
            cost_old = cost_new
            current_random_points = self.random_sample_selection(
                self.number_of_centres * size_multiple, n
            )
            distance_matrix = np.zeros(
                (current_random_points.shape[0], initial_centres.shape[0])
            )  # Vector to store distances from centroids
            current_centres = np.zeros(
                (current_random_points.shape[0], 1)
            )  # Vector containing the centroid each point belongs to

            # Calculate distance between random points and centres, sort and estimate new centres
            for i in range(0, self.number_of_centres):
                distance_matrix[:, i] = self.eucl_distance(
                    current_random_points, initial_centres[i, :]
                )
            current_centres = np.argmin(distance_matrix, axis=1)
            new_centres = self.create_centres(
                initial_centres, current_random_points, current_centres, counter
            )

            # Estimate distance between new and old centres
            distance_btw_centres = self.eucl_distance(new_centres, initial_centres)
            cost_new = np.sqrt(np.sum(distance_btw_centres**2))
            cost_change = np.abs(cost_old - cost_new)
            counter += 1
            # print(counter, cost_change)
            if cost_change >= self.eps:
                initial_centres = new_centres

        sample_points = new_centres

        unique_sample_points = self.sample_point_selection(
            self.data, sample_points, self.sampling_type
        )
        if len(self.data_headers) > 0 and self.df_flag:
            unique_sample_points = pd.DataFrame(
                unique_sample_points, columns=self.data_headers
            )
        return unique_sample_points
