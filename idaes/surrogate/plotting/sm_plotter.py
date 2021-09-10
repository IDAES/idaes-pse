##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
"""
Task: Artificial Intelligence/Machine Learning
Subtask: General Unified Surrogate Object - Plotting Methods
Author: B. Paul
"""

# Import statements
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from itertools import combinations  # used to pick xvar pairs in scatter3D


"""
Standard plotting methods (scatter, parity and residual plots).
"""


def scatter2D(xdata, zdata, xtest, zfit, xlabels=None, zlabels=None, dz=None,
              dzfit=None, clo=None, chi=None, clabel=None, show=True,
              PDF=False, filename='results_scatter2D.pdf'):
    """
    Plots training data as 2D scatter, with overlayed model output data. If
    available, can overlay confidence intervals on scatter plot and can
    generate additional scatter plots for first derivative data.

    Required: input variable data, output variable data, model fit data

    Optional: input and output labels (will be generated as needed), first
    derivative data (if given, derivative fit data required as well),
    confidence interval data.

    Generally, x (input) variables are indexed by 'i', z (output) variables
    are indexded by 'j', and dummy variables are indexed by 'w' when needed.
    """

    numins = np.shape(xdata)[1]  # number of input variables
    numouts = np.shape(zdata)[1]  # number of output variables

    # check and add labels if missing
    if xlabels is None:
        xlabels = []
        for i in range(numins):
            xlabels.append('x' + str(i+1))
    if zlabels is None:
        zlabels = []
        for j in range(numouts):
            zlabels.append('z' + str(j+1))
    if clabel is None:
        clabel = 'Confidence Intervals'

    pointidx = -1  # dummy index for each plot produced
    fig, ax = [], []
    for j in range(numouts):  # loop over all outputs, zj
        for i in range(numins):  # plot every possible zj = f(xi)
            pointidx += 1
            fig.append(plt.figure())
            ax.append(fig[pointidx].add_subplot())
            ax[pointidx].scatter(xdata[:, i], zdata[:, j], c='b', marker='s',
                                 label='Data')
            ax[pointidx].scatter(xtest[:, i], zfit[:, j], c='r',
                                 marker='.', label='Model')

            # plot confidence intervals if given data
            if clo is not None and chi is not None:
                ax[pointidx].scatter(xdata[:, i], clo[:, j], c='g',
                                     marker='.', label=clabel)
                ax[pointidx].scatter(xdata[:, i], chi[:, j], c='g',
                                     marker='.', label=clabel)

            ax[pointidx].set_xlabel(xlabels[i])
            ax[pointidx].set_ylabel(zlabels[j])
            ax[pointidx].set_title('2D Scatter Plot')
            ax[pointidx].legend()

            if dz is not None and dzfit is not None:
                pointidx += 1
                fig.append(plt.figure())
                ax.append(fig[pointidx].add_subplot())
                ax[pointidx].scatter(xdata[:, i], dz[:, j], c='b', marker='s',
                                     label='Data')
                ax[pointidx].scatter(xdata[:, i], dzfit[:, j], c='r',
                                     marker='.', label='Model')
                ax[pointidx].set_xlabel(xlabels[i])
                ax[pointidx].set_ylabel('d[' + str(zlabels[j]) + ']/d[' +
                                        xlabels[i] + ']')
                ax[pointidx].set_title('2D Derivative Plot')
                ax[pointidx].legend()

            if show is True:
                plt.show()
    if PDF is True:
        pdfPrint(fig, filename)


def scatter3D(xdata, zdata, xtest, zfit, xlabels=None, zlabels=None, show=True,
              PDF=False, filename='results_scatter3D.pdf'):
    """
    Plots training data as 3D scatter, with overlayed model output data.

    Required: input variable data, output variable data, model fit data

    Optional: input and output labels (will be generated as needed)

    Generally, x (input) variables are indexed by 'i', z (output) variables
    are indexded by 'j', and dummy variables are indexed by 'w' when needed.
    """

    numins = np.shape(xdata)[1]  # number of input variables
    numouts = np.shape(zdata)[1]  # number of output variables

    # check and add labels if missing
    if xlabels is None:
        xlabels = []
        for i in range(numins):
            xlabels.append('x' + str(i+1))
    if zlabels is None:
        zlabels = []
        for j in range(numouts):
            zlabels.append('z' + str(j+1))

    pairidx = -1  # dummy index for each plot produced
    fig, ax = [], []
    for j in range(numouts):  # loop over all outputs, zj
        for pair in list(combinations(range(numins), 2)):  # pick two x vars
            pairidx += 1  # index for new plot
            a, b = pair[0], pair[1]  # indices for the x variables picked
            fig.append(plt.figure())
            ax.append(fig[pairidx].add_subplot(projection='3d'))

            ax[pairidx].scatter(xdata[:, a], xdata[:, b], zdata[:, j], c='b',
                                marker='s', label='Data')
            ax[pairidx].scatter(xtest[:, a], xtest[:, b], zfit[:, j], c='r',
                                marker='.', label='Model')
            ax[pairidx].set_xlabel(xlabels[a])
            ax[pairidx].set_ylabel(xlabels[b])
            ax[pairidx].set_zlabel(zlabels[j])
            ax[pairidx].set_title('3D Scatter Plot')
            ax[pairidx].legend()

            if show is True:
                plt.show()
    if PDF is True:
        pdfPrint(fig, filename)


def parity(zdata, zfit, zlabels=None, clo=None, chi=None,
           clabel=None, show=True, PDF=False, filename='results_parity.pdf'):
    """
    Plots model output against data output, with confidence interval data
    if provided to the method.
    """
    numouts = np.shape(zdata)[1]  # number of output variables

    # check and add labels if missing

    if zlabels is None:
        zlabels = []
        for j in range(numouts):
            zlabels.append('z' + str(j+1))

    fig, ax = [], []
    for j in range(numouts):  # loop over all outputs, zj
        fig.append(plt.figure())
        ax.append(fig[j].add_subplot())

        ax[j].plot(zdata[:, j], zdata[:, j], c='b', label='Data')
        ax[j].scatter(zdata[:, j], zfit[:, j], c='r', marker='s',
                      label='Predictions')

        # plot confidence intervals if given data
        if clo is not None and chi is not None:
            ax[j].scatter(zdata[:, j], clo[:, j], c='g',
                          marker='.', label=clabel)
            ax[j].scatter(zdata[:, j], chi[:, j], c='g',
                          marker='.', label=clabel)

        ax[j].set_xlabel(zlabels[j])
        ax[j].set_ylabel('Model Output')
        ax[j].set_title('Parity Plot')
        ax[j].legend()

        if show is True:
            plt.show()
    if PDF is True:
        pdfPrint(fig, filename)


def residual(zdata, e, zlabels=None, elabel=None, show=True, PDF=False,
             filename='results_residual.pdf'):
    """
    Plots model error against data output.
    """

    numouts = np.shape(zdata)[1]  # number of output variables

    # check and add labels if missing

    if zlabels is None:
        zlabels = []
        for j in range(numouts):
            zlabels.append('z' + str(j+1))
    if elabel is None:
        elabel = 'Model Error'

    fig, ax = [], []
    for j in range(numouts):  # loop over all outputs, zj

        fig.append(plt.figure())
        ax.append(fig[j].add_subplot())
        ax[j].scatter(zdata[:, j], e[:, j], c='b', marker='s',
                      label=elabel)
        ax[j].set_xlabel(zlabels[j])
        ax[j].set_ylabel(elabel)
        ax[j].set_title('Residual Plot')
        ax[j].legend()

        if show is True:
            plt.show()
    if PDF is True:
        pdfPrint(fig, filename)


def pdfPrint(fig, filename):
    """
    Print input figure list to a single PDF file, with a specified name.
    """
    pp = PdfPages(filename)
    for plot in fig:
        pp.savefig(plot)
    pp.close()


def extractData(data):  # might not need if we ensure data is all 2D a priori
    """
    Ensures any input vectors (numrows,0) or lists (numrows,) are converted to
    arrays (numrows,1) to prevent indexing errors while plotting.

    Note - the plotting methods don't depend on this, and eventually this will
    be taken care by the SurrogateObject output.
    """
    if len(np.shape(data)) == 2 or data is None:  # check formatting
        array = np.array(data)
    else:  # reformat values as an n x 1 array
        temp = np.empty((np.shape(data)[0], 1))
        count = 0
        for val in data:
            temp[count] = data[count]
            count += 1
        array = np.array(temp)

    return array
