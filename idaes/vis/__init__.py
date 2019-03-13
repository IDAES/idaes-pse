##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
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
The `idaes.vis` subpackage contains the framework and implementation
of plots that are expected to be of general utility within the IDAES
framework.

For users, an entry point is provided for IDAES classes to produce
plots with the :class:`idaes.vis.plotbase.PlotRegistry` singleton.

Plots will inherit from the interface in :class:`idaes.vis.plotbase.PlotBase`,
which provides some basic methods.

The current implementations all use the Python "bokeh" package, and can
be found in :mod:`idaes.vis.bokeh_plots`.

For more details, please refer to the visualization section of the
main IDAES documentation.
"""
__all__ = ["bokeh_plots", "plotbase", "plotutils"]
