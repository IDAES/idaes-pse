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
This module contains utility functions to generate phase equilibrium data and
plots.
"""

__author__ = "Alejandro Garciadiego"

# Import objects from pyomo package
from pyomo.environ import (
    check_optimal_termination,
    value,
    units as pyunits,
)
import idaes.logger as idaeslog
from idaes.core.solvers import get_solver

import idaes.logger as idaeslog

# Import plotting functions
import matplotlib.pyplot as plt
import numpy as np


def Txy_diagram(
    model,
    component_1,
    component_2,
    pressure,
    num_points=20,
    temperature=298.15,
    figure_name=None,
    print_legend=True,
    include_pressure=False,
    print_level=idaeslog.NOTSET,
    solver=None,
    solver_op=None,
):

    """
    This method generates T-x-y plots. Given the components, pressure and property dictionary
    this function calls Txy_data() to generate T-x-y data and once the data has
    been generated calls build_txy_diagrams() to create a plot.

    Args:
        component_1: Component which composition will be plotted in x axis
        component_2: Component which composition will decrease in x axis
        pressure: Pressure at which the bubble and drew temperatures will be calculated
        temperature: Temperature at which to initialize state block
        num_points: Number of data point to be calculated
        properties: property package which contains parameters to calculate bubble
        and dew temperatures for the mixture of the compnents specified.
        figure_name: if a figure name is included the plot will save with the name
        figure_name.png
        print_legend (bool): = If True, include legend to distinguish between
            Bubble and dew temperature. The default is True.
        include_pressure (bool) = If True, print pressure at which the plot is
        calculated in legends. The default is False.
        print_level: printing level from initialization
        solver: solver to use (default=None, use IDAES default solver)
        solver_op: solver options

    Returns:
        Plot
    """
    # Run txy_ data funtion to obtain bubble and dew twmperatures
    Txy_data_to_plot = Txy_data(
        model,
        component_1,
        component_2,
        pressure,
        num_points,
        temperature,
        print_level,
        solver,
        solver_op,
    )

    # Run diagrams function to convert t-x-y data into a plot
    build_txy_diagrams(Txy_data_to_plot, figure_name, print_legend, include_pressure)


def Txy_data(
    model,
    component_1,
    component_2,
    pressure,
    num_points=20,
    temperature=298.15,
    print_level=idaeslog.NOTSET,
    solver=None,
    solver_op=None,
):
    """
    Function to generate T-x-y data. The function builds a state block and extracts
    bubble and dew temperatures at P pressure for N number of compositions.
    As N is increased increase the time of the calculation will increase and
    create a smoother looking plot.

    Args:
        component_1: Component 1
        component_2: Component 2
        pressure: Pressure at which the bubble and drew temperatures will be calculates
        temperature: Temperature at which to initialize state block
        num_points: Number of data point to be calculated
        model: Model wit intialized Property package which contains data to calculate
        bubble and dew temperatures for  component 1 and component 2
        print_level: printing level from initialization
        solver: solver to use (default=None, use IDAES default solver)
        solver_op: solver options

    Returns:
        (Class): A class containing the T-x-y data

    """

    components = list(model.params.component_list)
    components_used = [component_1, component_2]
    components_not_used = list(set(components) - set(components_used))

    # Add properties parameter blocks to the flowsheet with specifications

    model.props = model.params.build_state_block([1], defined_state=True)

    # Set intial concentration of component 1 close to 1
    x = 0.99

    # Set conditions for flash unit model
    model.props[1].mole_frac_comp[component_1].fix(x)

    for i in components_not_used:
        model.props[1].mole_frac_comp[i].fix(1e-5)

    xs = sum(value(model.props[1].mole_frac_comp[i]) for i in components_not_used)

    model.props[1].mole_frac_comp[component_2].fix(1 - x - xs)

    model.props[1].flow_mol.fix(1)
    model.props[1].temperature.fix(temperature)
    model.props[1].pressure.fix(pressure)

    # Initialize flash unit model
    model.props[1].calculate_scaling_factors()
    model.props.initialize(solver=solver, optarg=solver_op, outlvl=print_level)

    solver = get_solver(solver, solver_op)

    # Create an array of compositions with N number of points
    x_d = np.linspace(x, 1 - x - xs, num_points)

    # Create emprty arrays for concentration, bubble temperature and dew temperature
    X = []
    Tbubb = []
    Tdew = []

    # Obtain pressure and temperature units from the unit model
    Punit = pyunits.get_units(model.props[1].pressure)
    Tunit = pyunits.get_units(model.props[1].temperature)

    count = 1
    # Create and run loop to calculate temperatures at every composition
    for i in range(len(x_d)):
        model.props[1].mole_frac_comp[component_1].fix(x_d[i])
        model.props[1].mole_frac_comp[component_2].fix(1 - x_d[i] - xs)
        # solve the model
        status = solver.solve(model, tee=False)
        # If solution is optimal store the concentration, and calculated temperatures in the created arrays
        if check_optimal_termination(status):

            print(
                "Case: ", count, " Optimal. ", component_1, "x = {:.2f}".format(x_d[i])
            )

            if hasattr(model.props[1], "_mole_frac_tdew") and hasattr(
                model.props[1], "_mole_frac_tbub"
            ):
                Tbubb.append(value(model.props[1].temperature_bubble["Vap", "Liq"]))
                Tdew.append(value(model.props[1].temperature_dew["Vap", "Liq"]))

            elif hasattr(model.props[1], "_mole_frac_tdew"):
                print("One of the components only exists in vapor phase.")
                Tdew.append(value(model.props[1].temperature_dew["Vap", "Liq"]))

            elif hasattr(model.props[1], "_mole_frac_tbub"):
                print("One of the components only exists in liquid phase.")
                Tbubb.append(value(model.props[1].temperature_bubble["Vap", "Liq"]))

            X.append(x_d[i])

        # If the solver did not solve to an optimal solution, do not store the data point
        else:
            print(
                "Case: ", count, " No Result", component_1, "x = {:.2f}".format(x_d[i])
            )
        count += 1

    # Call TXYData function and store the data in TD class
    TD = TXYDataClass(component_1, component_2, Punit, Tunit, pressure)
    TD.TBubb = Tbubb
    TD.TDew = Tdew
    TD.x = X

    # Return the data class with all the information of the calculations
    return TD


# Author: Alejandro Garciadiego
class TXYDataClass:
    """
    Write data needed for build_txy_diagrams() into a class. The class can be
    obtained by running Txy_data() or by assigining values to the class.
    """

    def __init__(self, component_1, component_2, Punits, Tunits, pressure):
        """
        Args:
            component_1: Component 1
            component_2: Component 2
            Punits: Initial value of heat of hot utility
            Tunits: Initial value of heat to be removed by cold utility
            pressure: Pressure at which the T-x-y data was evaluated

        Returns:
            (Class): A class containing the T-x-y data
        """

        # Build
        self.Component_1 = component_1
        self.Component_2 = component_2

        # Assign units of pressure and temperature
        self.Punits = Punits
        self.Tunits = Tunits

        # Assign pressure at which the data has been calculated
        self.P = pressure

        # Create arrays for data
        self.TBubb = []
        self.TDew = []
        self.x = []

    def Temp_Bubb(self, data_list):
        """
        Args:
            data_list: Bubble temperature array

        Returns:
            None
        """
        self.TBubb = data_list

    def Temp_Dew(self, data_list_2):
        """
        Args:
            data_list_2: Dew temperature array

        Returns:
            None
        """
        self.Tdew = data_list_2

    def composition(self, data_list_3):
        """
        Args:
            data_list_3: x data array

        Returns:
            None
        """
        self.x = data_list_3


# Author: Alejandro Garciadiego
def build_txy_diagrams(
    txy_data, figure_name=None, print_legend=True, include_pressure=False
):
    """
    Args:
        txy_data: Txy data class includes components bubble and dew
        temperatures, compositions, components, pressure, and units.
        figure_name: if a figure name is included the plot will save with the name
        figure_name.png
        print_legend (bool): = If True, include legend to distinguish between
            Bubble and dew temperature. The default is True.
        include_pressure (bool) = If True, print pressure at which the plot is
        calculated in legends. The default is False.
    Returns:
        t-x-y plot
    """

    # Declare a plot and it's size
    ig, ax = plt.subplots(figsize=(12, 8))

    if len(txy_data.TBubb) and len(txy_data.TDew) > 0:
        if include_pressure == True:
            # Plot results for bubble temperature

            ax.plot(
                txy_data.x,
                txy_data.TBubb,
                "r",
                label="Bubble Temp P = "
                + str(txy_data.Press)
                + " "
                + str(txy_data.Punits),
                linewidth=1.5,
            )

            # Plot results for dew temperature
            ax.plot(
                txy_data.x,
                txy_data.TDew,
                "b",
                label="Dew Temp P = "
                + str(txy_data.Press)
                + " "
                + str(txy_data.Punits),
                linewidth=1.5,
            )

        else:
            # Plot results for bubble temperature
            ax.plot(
                txy_data.x, txy_data.TBubb, "r", label="Bubble Temp ", linewidth=1.5
            )

            # Plot results for dew temperature
            ax.plot(txy_data.x, txy_data.TDew, "b", label="Dew Temp", linewidth=1.5)

    elif len(txy_data.TDew) == 0:

        if include_pressure == True:
            # Plot results for bubble temperature

            # Plot results for dew temperature
            ax.plot(
                txy_data.x,
                txy_data.TBubb,
                "b",
                label="Dew Temp P = "
                + str(txy_data.Press)
                + " "
                + str(txy_data.Punits),
                linewidth=1.5,
            )

        else:
            # Plot results for dew temperature
            ax.plot(txy_data.x, txy_data.TBubb, "b", label="Dew Temp", linewidth=1.5)

    elif len(txy_data.TBubb) == 0:

        if include_pressure == True:
            # Plot results for bubble temperature

            # Plot results for dew temperature
            ax.plot(
                txy_data.x,
                txy_data.TDew,
                "b",
                label="Bubble Temp P = "
                + str(txy_data.Press)
                + " "
                + str(txy_data.Punits),
                linewidth=1.5,
            )

        else:
            # Plot results for dew temperature
            ax.plot(txy_data.x, txy_data.TDew, "b", label="Dew Temp", linewidth=1.5)

    # Include grid
    ax.grid()

    # Declare labels and fontsize
    plt.xlabel(txy_data.Component_1 + " concentration (mol/mol)", fontsize=20)
    plt.ylabel("Temperature [" + str(txy_data.Tunits) + "]", fontsize=20)

    # Declare plot title
    plt.title(
        "T-x-y diagram " + txy_data.Component_1 + "-" + txy_data.Component_2,
        fontsize=24,
    )

    # Set limits of 0-1 mole fraction
    plt.xlim(0.0, 1)

    # Declare legend and fontsize
    if print_legend == True:
        plt.legend(fontsize=16)

    if figure_name:
        plt.savefig(str(figure_name) + ".png")

    # Show plot
    plt.show()
