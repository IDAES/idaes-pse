
Visualization
=============

Contents
--------

.. toctree::
    :maxdepth: 1

    modelvis
    vis_hx
    vis_prof

.. warning::
    The visualization library is still in active development and we
    hope to improve on it in future releases. Please use its
    functionality at your own discretion.

Overview
--------

The JupyterLab Model Visualization Tool is composed of a JupyterLab extension 
that displays a complete flowsheet. The tool is used in conjunction with the 
`idaes.dmf.ui` subpackage, which serializes the model as a `.idaes.vis` file 
which is read by the JupyterLab extension.

The `idaes.vis` subpackage contains the framework and implementation
of plots that are expected to be of general utility within the IDAES
framework.

For users, an entry point is provided for IDAES classes to produce
plots with the :class:`idaes.vis.plotbase.PlotRegistry` singleton.

Plots will inherit from the interface in :class:`idaes.vis.plotbase.PlotBase`,
which provides some basic methods.

The current implementations all use the Python "bokeh" package, and can
be found in :mod:`idaes.vis.bokeh_plots`.



