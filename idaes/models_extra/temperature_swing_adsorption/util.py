#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Utility functions for temperature swing adsorption models
"""

__author__ = "Daison Yancy Caballero"

from sys import stdout
import textwrap
from pandas import DataFrame
import matplotlib.pyplot as plt

from pyomo.environ import value

from idaes.core.util.tables import stream_table_dataframe_to_string


def tsa_summary(tsa, stream=stdout, export=False):
    """
    Prints a summary of tsa variables.

    Args:
        tsa (FixedBedTSA0D): FixedBedTSA0D model to summarize.
        stream (iostream, optional): Output location. Defaults to sys.stdout.
        export (bool, optional): Flag to save as csv file. Defaults to False.

    Returns:
        None.

    """

    var_dict = tsa.get_var_dict()

    summary_dir = {}
    summary_dir["Value"] = {}
    summary_dir["pos"] = {}

    count = 1
    for k, v in var_dict.items():
        summary_dir["Value"][k] = value(v)
        summary_dir["pos"][k] = count
        count += 1

    df = DataFrame.from_dict(summary_dir, orient="columns")
    del df["pos"]
    if export:
        df.to_csv(f"{tsa.local_name}_summary.csv")

    stream.write("\n" + "=" * 84)
    stream.write("\n" + f"Summary - {tsa.local_name}")
    stream.write("\n" + "-" * 84)
    stream.write(textwrap.indent(stream_table_dataframe_to_string(df), " " * 4))
    stream.write("\n" + "=" * 84 + "\n")

    return df


def plot_tsa_profiles(tsa):
    """
    Plot method for common fixed bed TSA variables

    Variables plotted:
        T : Temperature [K]
        P : Pressure [bar]
        y : Mole fraction of CO2 in the gas phase

    """
    # heating profiles
    time_domain = tsa.heating.time_domain
    t_heating = list(time_domain)
    t_heating = [i * value(tsa.heating.time) / 60 for i in t_heating]
    y_heating = [value(tsa.heating.mole_frac[t, "CO2"]) for t in time_domain]
    T_heating = [value(tsa.heating.temperature[t]) for t in time_domain]
    P_heating = [value(tsa.pressure_adsorption) * 1e-5 for i in range(len(t_heating))]

    # cooling profiles
    time_domain = tsa.cooling.time_domain
    t_cooling = list(time_domain)
    t_cooling = [t_heating[-1] + i * value(tsa.cooling.time) / 60 for i in t_cooling]
    y_cooling = [value(tsa.cooling.mole_frac[t, "CO2"]) for t in time_domain]
    T_cooling = [value(tsa.cooling.temperature[t]) for t in time_domain]
    P_cooling = [value(tsa.cooling.pressure[t]) for t in time_domain]
    P_cooling = [i * 1e-5 for i in P_cooling]

    # pressurization profiles
    t_pressurization = [
        t_cooling[-1] * 60
        + i
        * ((t_cooling[-1] * 60 + value(tsa.pressurization.time)) - t_cooling[-1] * 60)
        / (1000 - 1)
        for i in range(1000)
    ]
    t_pressurization.pop(0)
    t_pressurization = [i / 60 for i in t_pressurization]
    n_press = len(t_pressurization)
    tf = tsa.cooling.time_domain.last()
    y_pressurization = [value(tsa.cooling.mole_frac[tf, "CO2"])] * n_press
    y_pressurization[-1] = value(tsa.pressurization.mole_frac["CO2"])
    T_pressurization = [value(tsa.cooling.temperature[tf])] * n_press
    P_pressurization = [value(tsa.cooling.pressure[tf]) * 1e-5] * n_press
    P_pressurization[-1] = value(tsa.pressure_adsorption) * 1e-5

    # adsorption profiles
    t_adsorption = [
        t_pressurization[-1] * 60
        + i
        * (
            (t_pressurization[-1] * 60 + value(tsa.adsorption.time))
            - t_pressurization[-1] * 60
        )
        / (1000 - 1)
        for i in range(1000)
    ]
    t_adsorption.pop(0)
    t_adsorption = [i / 60 for i in t_adsorption]
    n_ads = len(t_adsorption)
    y_adsorption = [value(tsa.pressurization.mole_frac["CO2"])] * n_ads
    y_adsorption[-1] = value(tsa.mole_frac_in["CO2"])
    T_adsorption = [value(tsa.temperature_adsorption)] * n_ads
    P_adsorption = [value(tsa.pressure_adsorption) * 1e-5] * n_ads

    # plotting profiles - whole cycle
    t_cycle = t_heating + t_cooling + t_pressurization + t_adsorption
    y_cycle = y_heating + y_cooling + y_pressurization + y_adsorption
    T_cycle = T_heating + T_cooling + T_pressurization + T_adsorption
    P_cycle = P_heating + P_cooling + P_pressurization + P_adsorption

    plots = plt.subplots(3, figsize=(5, 8))
    (ax1, ax2, ax3) = plots[1]
    ax1.plot(t_cycle, T_cycle)
    ax1.set_ylabel("Temperature [K]")
    ax1.set_xlim(0, max(t_cycle) + 1)
    ax1.set_ylim(value(tsa.temperature_cooling), value(tsa.temperature_heating))
    ax1.axvline(t_heating[-1], lw=0.7, color="r", ls="--")
    ax1.axvline(t_cooling[-1], lw=0.7, color="r", ls="--")
    ax1.axvline(t_pressurization[-1], lw=0.7, color="r", ls="--")
    ax1.axvline(t_adsorption[-1], lw=0.7, color="r", ls="--")
    ax1.grid()
    ax2.plot(t_cycle, P_cycle)
    ax2.set_ylabel("Pressure [bar]")
    ax2.set_xlim(0, max(t_cycle) + 1)
    ax2.set_ylim(0, 1.05)
    ax2.axvline(t_heating[-1], lw=0.7, color="r", ls="--")
    ax2.axvline(t_cooling[-1], lw=0.7, color="r", ls="--")
    ax2.axvline(t_pressurization[-1], lw=0.7, color="r", ls="--")
    ax2.axvline(t_adsorption[-1], lw=0.7, color="r", ls="--")
    ax2.grid()
    ax3.plot(t_cycle, y_cycle)
    ax3.set_ylabel(r"Mole Fraction of $\mathrm{CO_{2}}$ [-]")
    ax3.set_xlabel("Time [min]")
    ax3.set_xlim(0, max(t_cycle) + 1)
    ax3.set_ylim(0, 1.05)
    ax3.axvline(t_heating[-1], lw=0.7, color="r", ls="--")
    ax3.axvline(t_cooling[-1], lw=0.7, color="r", ls="--")
    ax3.axvline(t_pressurization[-1], lw=0.7, color="r", ls="--")
    ax3.axvline(t_adsorption[-1], lw=0.7, color="r", ls="--")
    ax3.grid()

    plt.show()
