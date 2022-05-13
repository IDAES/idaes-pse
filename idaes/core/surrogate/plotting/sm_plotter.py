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
Task: Artificial Intelligence/Machine Learning
Subtask: General Unified Surrogate Object - Plotting Methods
Author: B. Paul
"""

# Import statements
import numpy as np
import matplotlib.pyplot as plt

# TODO: Consider adding seaborn in future
from matplotlib.backends.backend_pdf import PdfPages

from itertools import combinations  # used to pick xvar pairs in scatter3D


"""
Standard plotting methods (scatter, parity and residual plots).
"""


def surrogate_scatter2D(surrogate, dataframe, filename=None, show=True):
    input_data = dataframe[surrogate.input_labels()]
    output_data = dataframe[surrogate.output_labels()]
    output_surrogate = surrogate.evaluate_surrogate(input_data)
    _scatter2D(
        xdata=input_data.values,
        zdata=output_data.values,
        zfit=output_surrogate.values,
        xlabels=surrogate.input_labels(),
        zlabels=surrogate.output_labels(),
        show=show,
        filename=filename,
    )


def _scatter2D(
    xdata, zdata, zfit, xlabels=None, zlabels=None, show=True, filename=None
):
    """
    Plots training data as 2D scatter, with overlayed model output data.

    Required: input variable data, output variable data, model fit data

    Optional: input and output labels (will be generated as needed), filename
    to save plots to PDF

    Generally, x (input) variables are indexed by 'i', and z (output) variables
    are indexed by 'j'.
    """

    numins = np.shape(xdata)[1]  # number of input variables
    numouts = np.shape(zdata)[1]  # number of output variables

    # check and add labels if missing
    if xlabels is None:
        xlabels = []
        for i in range(numins):
            xlabels.append("x" + str(i + 1))
    if zlabels is None:
        zlabels = []
        for j in range(numouts):
            zlabels.append("z" + str(j + 1))

    if filename is not None:
        if filename[-4:] != ".pdf":  # checking if extension is present
            raise Exception(
                "Filename does not end in .pdf, please amend "
                "filename with the proper extension."
            )
        pdf = PdfPages(filename)

    for j in range(numouts):  # loop over all outputs, zj
        for i in range(numins):  # plot every possible zj = f(xi)
            fig = plt.figure()
            ax = fig.add_subplot()
            ax.scatter(xdata[:, i], zdata[:, j], c="grey", marker="s", label="Data")
            ax.scatter(xdata[:, i], zfit[:, j], c="b", marker=".", label="Model")

            ax.set_xlabel(xlabels[i])
            ax.set_ylabel(zlabels[j])
            ax.set_title("2D Scatter Plot")
            ax.legend()

            if show is True:
                plt.show()
            if filename is not None:
                pdf.savefig(fig)
    if filename is not None:  # place outside loop to avoid closing/reopening
        pdf.close()


def surrogate_scatter3D(surrogate, dataframe, filename=None, show=True):
    input_data = dataframe[surrogate.input_labels()]
    output_data = dataframe[surrogate.output_labels()]
    output_surrogate = surrogate.evaluate_surrogate(input_data)
    _scatter3D(
        xdata=input_data.values,
        zdata=output_data.values,
        zfit=output_surrogate.values,
        xlabels=surrogate.input_labels(),
        zlabels=surrogate.output_labels(),
        show=show,
        filename=filename,
    )


def _scatter3D(
    xdata, zdata, zfit, xlabels=None, zlabels=None, show=True, filename=None
):
    """
    Plots training data as 3D scatter, with overlayed model output data.

    Required: input variable data, output variable data, model fit data

    Optional: input and output labels (will be generated as needed), filename
    to save plots to PDF

    Generally, x (input) variables are indexed by 'i', and z (output) variables
    are indexed by 'j'.
    """

    numins = np.shape(xdata)[1]  # number of input variables
    numouts = np.shape(zdata)[1]  # number of output variables

    # check and add labels if missing
    if xlabels is None:
        xlabels = []
        for i in range(numins):
            xlabels.append("x" + str(i + 1))
    if zlabels is None:
        zlabels = []
        for j in range(numouts):
            zlabels.append("z" + str(j + 1))

    if filename is not None:
        if filename[-4:] != ".pdf":  # checking if extension is present
            raise Exception(
                "Filename does not end in .pdf, please amend "
                "filename with the proper extension."
            )
        pdf = PdfPages(filename)

    for j in range(numouts):  # loop over all outputs, zj
        for pair in list(combinations(range(numins), 2)):  # pick two x vars
            a, b = pair[0], pair[1]  # indices for the x variables picked
            fig = plt.figure()
            ax = fig.add_subplot(projection="3d")

            ax.scatter(
                xdata[:, a],
                xdata[:, b],
                zdata[:, j],
                c="grey",
                marker="s",
                label="Data",
            )
            ax.scatter(
                xdata[:, a], xdata[:, b], zfit[:, j], c="b", marker=".", label="Model"
            )
            ax.set_xlabel(xlabels[a])
            ax.set_ylabel(xlabels[b])
            ax.set_zlabel(zlabels[j])
            ax.set_title("3D Scatter Plot")
            ax.legend()

            if show is True:
                plt.show()
            if filename is not None:
                pdf.savefig(fig)
    if filename is not None:  # place outside loop to avoid closing/reopening
        pdf.close()


def surrogate_parity(surrogate, dataframe, filename=None, show=True):
    input_data = dataframe[surrogate.input_labels()]
    output_data = dataframe[surrogate.output_labels()]
    output_surrogate = surrogate.evaluate_surrogate(input_data)
    _parity(
        zdata=output_data.values,
        zfit=output_surrogate.values,
        zlabels=surrogate.output_labels(),
        show=show,
        filename=filename,
    )


def _parity(zdata, zfit, zlabels=None, show=True, filename=None):
    """
    Plots model output against data output.

    Required: output variable data, model fit data

    Optional: output labels (will be generated as needed), filename to save
    plots to PDF

    Generally, x (input) variables are indexed by 'i', and z (output) variables
    are indexed by 'j'.
    """
    numouts = np.shape(zdata)[1]  # number of output variables

    # check and add labels if missing

    if zlabels is None:
        zlabels = []
        for j in range(numouts):
            zlabels.append("z" + str(j + 1))

    if filename is not None:
        if filename[-4:] != ".pdf":  # checking if extension is present
            raise Exception(
                "Filename does not end in .pdf, please amend "
                "filename with the proper extension."
            )
        pdf = PdfPages(filename)

    for j in range(numouts):  # loop over all outputs, zj
        fig = plt.figure()
        ax = fig.add_subplot()

        ax.plot(zdata[:, j], zdata[:, j], c="grey", label="Data")
        ax.scatter(zdata[:, j], zfit[:, j], s=3, c="b", marker=".", label="Predictions")

        ax.set_xlabel(zlabels[j])
        ax.set_ylabel("Model Output")
        ax.set_title("Parity Plot")
        ax.legend()

        if show is True:
            plt.show()
        if filename is not None:
            pdf.savefig(fig)
    if filename is not None:  # place outside loop to avoid closing/reopening
        pdf.close()


def surrogate_residual(
    surrogate, dataframe, filename=None, relative_error=False, show=True
):
    input_data = dataframe[surrogate.input_labels()]
    output_data = dataframe[surrogate.output_labels()]
    output_surrogate = surrogate.evaluate_surrogate(input_data)
    residual = np.abs(output_data - output_surrogate)
    if relative_error is True:
        residual = np.divide(residual, np.maximum(output_data, 1.0))
    _residual(
        xdata=input_data.values,
        residual=residual.values,
        xlabels=surrogate.input_labels(),
        elabels=surrogate.output_labels(),
        show=show,
        filename=filename,
    )


def _residual(xdata, residual, xlabels=None, elabels=None, show=True, filename=None):
    """
    Plots model residual against data output.

    Required: input variable data, model error data

    Optional: input and error labels (will be generated as needed), filename
    to save plots to PDF

    Generally, x (input) variables are indexed by 'i', and z (output) variables
    are indexed by 'j'.
    """
    numins = np.shape(xdata)[1]  # number of input variables
    numouts = np.shape(residual)[1]  # number of output variables

    # check and add labels if missing

    if xlabels is None:
        xlabels = []
        for j in range(numins):
            xlabels.append("x" + str(j + 1))
    if elabels is None:
        for i in range(numouts):
            elabels.append("z" + str(j + 1) + " Error")

    if filename is not None:
        if filename[-4:] != ".pdf":  # checking if extension is present
            raise Exception(
                "Filename does not end in .pdf, please amend "
                "filename with the proper extension."
            )
        pdf = PdfPages(filename)

    for i in range(numins):
        for j in range(numouts):  # loop over all outputs, zj
            fig = plt.figure()
            ax = fig.add_subplot()
            ax.scatter(xdata[:, i], residual[:, j], marker=".", label=elabels[j])
            ax.set_xlabel(xlabels[i])
            ax.set_ylabel(elabels[j])
            ax.set_title("Residual Plot")
            ax.legend()

            if show is True:
                plt.show()
            if filename is not None:
                pdf.savefig(fig)
    if filename is not None:  # place outside loop to avoid closing/reopening
        pdf.close()
