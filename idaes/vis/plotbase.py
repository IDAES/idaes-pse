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
Base classes for visualization and plotting in IDAES.

Create new plots by inheriting from :class:`PlotBase`. See the
:mod:`idaes.vis.bokeh_plots` module for examples.
"""
# stdlib
from abc import ABC, abstractmethod
from typing import Type, List
# third party
import pandas as pd


class PlotBase(ABC):
    """Abstract base class for a plot.
    """

    def __init__(self, current_plot):
        """Constructor.

        Args:
            current_plot: Set the current plot.
        """
        self.current_plot = current_plot

    @abstractmethod
    def annotate(self, x, y, label: str):
        """Annotate a plot with a given point and a label.

        Args:
            x: Value of independent variable.
            y: Value of dependent variable.
            label: Text label.
        """
        pass

    @abstractmethod
    def resize(self, height: int = -1, width: int = -1):
        """Resize a plot's height and width.

        Args:
            height: Height in screen units.
            width: Width in screen units.

        Returns:
            None

        Raises:
            None
        """
        pass

    @abstractmethod
    def save(self, destination: str):
        """Save the current plot object to HTML in filepath provided by destination.

        Args:
            destination: Valid file path to save HTML to.

        Returns:
            filename where HTML is saved.

        Raises:
            None
        """
        pass

    @abstractmethod
    def show(self, in_notebook=True):
        """Display plot in a Jupyter notebook.

        Args:
            in_notebook: Display in Jupyter notebook or generate HTML file.

        Returns:
            None

        Raises:
            None
        """
        pass

    @classmethod
    def validate(cls, data_frame: pd.DataFrame, x: str, y: List, legend=None):
        """Validate that the plot parameters are valid.

        Args:
            data_frame: a pandas data frame of any type.
            x: Key in data-frame to use as x-axis.
            y: Keys in data-frame to use as y-axis.
            legend: List of labels to use as legend for a plot.
        Returns:
            True, '' on valid data frames (if x and y are in the data frame keys)
            False, "message" on invalid data
        """
        if data_frame is None or data_frame.empty:
            return False, "Empty data frame"
        if not legend and legend is not None:
            return False, "Bad legend labels"
        data_keys = set(data_frame.keys())
        if x not in data_keys:
            return False, f"The x-axis key {x} is not a column"
        if not set(y).issubset(data_keys):
            diff = set(y) - data_keys
            if len(diff) > 1:
                return False, f"The y-axis keys {diff} are not columns"
            else:
                return False, f"y-axis key {diff} is not a column"
        return True, ''


class PlotRegistry(object):
    """Set of associations between objects/classes + a plot name, and
    the parameters and values needed to perform the plot.

    The basic idea is to create a set of named plots associated with
    a given IDAES class, and then allow the user or other APIs to invoke that
    plot once the data is populated in an instance of the class. This keeps
    the details of how to create plots of a given type in the classes that
    will create them.

    For example::

        class MyIdaesClass(ProcessBase):
          # .. code for the class
          def plot_setup(self, plot_class):
              # .. details of creating plot_instance from object contents ..
              return plot_instance
        PlotRegistry().register(MyIdaesClass, 'basic', MyIdaesClass.plot_setup)

        # .. and, later ..
        obj = MyIdaesClass(...)
        # .. do things that fill "obj" with data ..
        # now create the plot
        plot = PlotRegistry().get(obj, 'basic')
        plot.show()

    XXX: This class is not actually used (yet) by any of the IDAES models.
    """

    __objs = {}  # shared state for all instances

    def __init__(self):
        self._objs = self.__objs

    def register(
        self, obj, name: str, plot: Type, setup_fn=None, overwrite: bool = False
    ):
        """Register an object/plot combination.

        Args:
            obj: Class or instance
            name: Name for this association
            plot: Plot class
            setup_fn: Optional setup function to call. Function should take
                      two arguments: `plot` class instance, `obj` assoc. with plot.
            overwrite: If true, allow overwrite of existing entry in the registry
        """
        if obj in self._objs:
            entry = self._objs[obj]
        else:
            entry = {}
            self._objs[obj] = entry
        if name in entry and not overwrite:
            raise KeyError(
                f"Cannot register {name}: would overwrite existing "
                f"entry for object, and 'overwrite' parameter is not set"
            )
        entry[name] = (plot, setup_fn)

    def get(self, obj, name: str):
        """Get a plot object for the given object + name.

        Args:
            obj: Object for which to get the plot
            name: Registered name of plot to get

        Returns:
            Return value of setup function given to :meth:`register()`,
            or, if that is empty, the retrieved plot object.
        """
        if obj not in self._objs:
            raise KeyError(f"No plots registered for {obj}")
        entry = self._objs[obj]
        if name not in entry:
            raise KeyError(f"No plot named '{name}' found in entry for object {obj}")
        plot, setup_fn = entry[name]
        if setup_fn is not None:
            result = setup_fn(plot, obj)
        else:
            result = plot
        return result

    def remove_all(self):
        """Remove all entries from the registry.

        Since the registry is a singleton, this removes all entries from ALL
        instances. Use with care.
        """
        for k in list(self._objs.keys()):
            del self._objs[k]
