
Visualization
=============

Contents
--------

.. toctree::
    :maxdepth: 1

    vis_hx
    vis_prof

.. warning::
    The visualization library is still in active development and we
    hope to improve on it in future releases. Please use its
    functionality at your own discretion.

Overview
--------

The `idaes.vis` subpackage contains the framework and implementation
of plots that are expected to be of general utility within the IDAES
framework.

For users, an entry point is provided for IDAES classes to produce
plots with the :class:`idaes.vis.plotbase.PlotRegistry` singleton.

Plots will inherit from the interface in :class:`idaes.vis.plotbase.PlotBase`,
which provides some basic methods.

The current implementations all use the Python "bokeh" package, and can
be found in :mod:`idaes.vis.bokeh_plots`.



