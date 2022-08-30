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
This module contains functions to read and manage data for use in parameter
esitmation, data reconciliation, and validation.
"""

__author__ = "John Eslick"

import logging
import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pint
import os

import pyomo.environ as pyo
import warnings

try:
    import seaborn as sns
    from PyPDF2 import PdfFileMerger
except ImportError:
    sns = None

_log = logging.getLogger(__file__)


def _strip(tag):
    """
    Tag renaming function to remove whitespace, depending on the csv format
    column heading items can pick up some extra whitespace when read into Pandas
    """
    return tag.strip()


# Some common unit string conversions, these are ones we've come across that
# are not handled by pint. We can organize and refine known unit strings
# more in the future.
_unit_strings = {
    # Pressure
    "PSI": "psi",
    "PSIA": "psi",
    "psia": "psi",
    "PSIG": "psig",
    "INWC": "inH2O gauge",
    "IN WC": "inH2O gauge",
    "IN/WC": "inH2O gauge",
    '" H2O': "inH2O gauge",
    "INHG": "inHg",
    "IN HG": "inHg",
    "IN/HG": "inHg",
    "HGA": "inHg",
    "IN HGA": "inHg",
    # Fraction
    "PCT": "percent",
    "pct": "percent",
    "PERCT": "percent",
    "PERCT.": "percent",
    "PCNT": "percent",
    "PPM": "ppm",
    "PPB": "ppb",
    "% OPEN": "percent open",
    "% CLSD": "percent closed",
    "% CLOSED": "percent closed",
    # Length
    "IN": "in",
    "INS": "in",
    "INCHES": "in",
    "Inches": "in",
    "FT": "ft",
    "FEET": "ft",
    "FOOT": "ft",
    "Feet": "ft",
    "MILS": "minch",
    # Speed
    "MPH": "mile/hr",
    "IPS": "in/s",
    # Volume
    "KGAL": "kgal",
    # Vol Flow
    "GPM": "gal/min",
    "gpm": "gal/min",
    "CFM": "ft^3/min",
    "KCFM": "ft^3/mmin",
    "SCFM": "ft^3/min",
    "KSCFM": "ft^3/mmin",  # be careful with this one
    # I don't know how to indicate its
    # a volumetric flow at standard
    # conditions
    # Angle
    "DEG": "deg",
    # Angular Speed
    "RPM": "rpm",
    # Fequency
    "HZ": "hz",
    # Temperature
    "DEG F": "degF",
    "Deg F": "degF",
    "deg F": "degF",
    "DEG C": "degC",
    "Deg C": "degC",
    "deg C": "degC",
    "DEGF": "degF",
    "DegF": "degF",
    "DEGC": "degC",
    "DegC": "degC",
    # Temperature Difference
    "DELTA DEG F": "delta_degF",
    "DETLA Deg F": "delta_degF",
    "DETLA deg F": "delta_degF",
    "DETLA DEG C": "delta_degC",
    "DETLA Deg C": "delta_degC",
    "DELTA deg C": "delta_degC",
    "DELTA DEGF": "delta_degF",
    "DELTA DegF": "delta_degF",
    "DELTA degF": "delta_degF",
    "DELTA DEGC": "delta_degC",
    "DELTA DegC": "delta_degC",
    "DELTA degC": "delta_degC",
    "Delta DEG F": "delta_degF",
    "Delta Deg F": "delta_degF",
    "Delta deg F": "delta_degF",
    "Delta DEG C": "delta_degC",
    "Delta Deg C": "delta_degC",
    "Delta deg C": "delta_degC",
    "Delta DEGF": "delta_degF",
    "Delta DegF": "delta_degF",
    "Delta degF": "delta_degF",
    "Delta DEGC": "delta_degC",
    "Delta DegC": "delta_degC",
    "Delta degC": "delta_degC",
    "delta DEG F": "delta_degF",
    "delta Deg F": "delta_degF",
    "delta deg F": "delta_degF",
    "delta DEG C": "delta_degC",
    "delta Deg C": "delta_degC",
    "delta deg C": "delta_degC",
    "delta DEGF": "delta_degF",
    "delta DegF": "delta_degF",
    "delta degF": "delta_degF",
    "delta DEGC": "delta_degC",
    "delta DegC": "delta_degC",
    "delta degC": "delta_degC",
    # Energy
    "MBTU": "kbtu",
    # Mass
    "MLB": "klb",
    "K LB": "klb",
    "K LBS": "klb",
    "lb.": "lb",
    # Mass flow
    "TPH": "ton/hr",
    "tph": "ton/hr",
    "KLB/HR": "klb/hr",
    "KPPH": "klb/hr",
    # Current
    "AMP": "amp",
    "AMPS": "amp",
    "Amps": "amp",
    "Amp": "amp",
    "AMP AC": "amp",
    # pH
    "PH": "pH",
    # VARS (volt-amp reactive)
    "VARS": "VAR",
    "MVARS": "MVAR",
}

_gauge_pressures = {"psig": "psi", "inH2O gauge": "inH2O"}

_ignore_units = [
    "percent",
    "ppm",
    "ppb",
    "pH",
    "VAR",
    "MVAR",
    "H2O",
    "percent open",
    "percent closed",
]

_register_new_units = [
    "in_H2O = 248.84 Pa = inH2O",
]


def unit_convert(
    x,
    frm,
    to=None,
    system=None,
    unit_string_map={},
    ignore_units=[],
    gauge_pressures={},
    ambient_pressure=1.0,
    ambient_pressure_unit="atm",
):
    """Convert the quantity x to a different set of units. X can be a numpy array
    or pandas series. The from unit is translated into a string that pint
    can recognize by first looking in unit_string_map then looking in
    know aliases defined in this file. If it is neither place it will be given
    to pint as-is. This translation of the unit is done so that data can be read
    in with the original provided units.

    Args:
        x (float, numpy.array, pandas.series): quantity to convert
        frm (str): original unit string
        to (str): new unit string, or specify "system"
        system (str): unit system to covert to, or specify "to"
        unit_string_map (dict): keys are unit strings and values are
            corresponding strings that pint can recognize.  This only applies to
            the from string.
        ignore_units (list, or tuple): units to not convert
        gauge_pressures (dict): keys are units strings to be considered gauge
            pressures and the values are corresponding absolute pressure units
        ambient_pressure (float, numpy.array, pandas.series): pressure to add
            to gauge pressure to convert it to absolute pressure. The default
            is 1. The unit is atm by default, but can be changed with the
            ambient_pressure_unit argument.
        ambient_pressure_unit (str): Unit for ambient pressure, default is atm,
            and should be a unit recognized by pint
    Returns:
        (tuple): quantity and unit string
    """
    ureg = pint.UnitRegistry(system=system)
    for u in _register_new_units:
        ureg.define(u)
    if frm in unit_string_map:
        frm = unit_string_map[frm]
    elif frm in _unit_strings:
        frm = _unit_strings[frm]
    # Now check for gauge pressure
    gauge = False
    if frm in gauge_pressures:
        gauge = True
        frm = gauge_pressures[frm]
    elif frm in _gauge_pressures:
        gauge = True
        frm = _gauge_pressures[frm]
    q = ureg.Quantity
    if (frm in _ignore_units) or (frm in ignore_units):
        return (x, frm)
    else:
        try:
            ureg.parse_expression(frm)
        except pint.errors.UndefinedUnitError:
            warnings.warn(
                "In unit conversion, from unit '{}' is not defined."
                " No conversion.".format(frm),
                UserWarning,
            )
            return x, frm
    if to is None:
        y = q(np.array(x), ureg.parse_expression(frm)).to_base_units()
    else:
        y = q(np.array(x), ureg.parse_expression(frm)).to(to)
    if gauge:
        # convert gauge pressure to absolute
        y = y + ambient_pressure * ureg.parse_expression(ambient_pressure_unit)
    return (y.magnitude, str(y.units))


def update_metadata_model_references(model, metadata):
    """
    Create model references from refernce strings in the metadata. This updates
    the 'reference' field in the metadata.

    Args:
        model (pyomo.environ.Block): Pyomo model
        metadata (dict): Tag metadata dictionary

    Returns:
        None
    """
    for tag, md in metadata.items():
        if md["reference_string"]:
            try:
                md["reference"] = pyo.Reference(
                    eval(md["reference_string"], {"m": model})
                )
            except (KeyError, AttributeError, NameError):
                warnings.warn(
                    "Tag reference {} not found".format(md["reference_string"]),
                    UserWarning,
                )


# This prevents breakage, can remove after example updates.
upadate_metadata_model_references = update_metadata_model_references


def read_data(
    csv_file,
    csv_file_metadata,
    model=None,
    rename_mapper=None,
    unit_system=None,
    ambient_pressure=1.0,
    ambient_pressure_unit="atm",
):
    """
    Read CSV data into a Pandas DataFrame.

    The data should be in a form where the first row contains column headings
    where each column is labeled with a data tag, and the first column contains
    data point labels or time stamps. The metadata should be in a csv file where
    the first column is the tag name, the second column is the model reference (
    which can be empty), the third column is the tag description, and the fourth
    column is the unit of measure string. Any additional information can be
    added to columns after the fourth column and will be ignored. The units of
    measure should be something that is recognized by pint, or in the aliases
    defined in this file. Any tags not listed in the metadata will be dropped.

    The function returns two items a pandas.DataFrame containing process data,
    and a dictionary with tag metadata.  The metadata dictionary keys are tag name,
    and the values are dictionaries with the keys: "reference_string", "description",
    "units", and "reference".


    Args:
        csv_file (str): Path of file to read
        csv_file_metadata (str): Path of csv file to read column metadata from
        model (pyomo.environ.ConcreteModel): Optional model to map tags to
        rename_mapper (Callable): Optional function to rename tags
        unit_system (str): Optional system of units to atempt convert to
        ambient_pressure (float, numpy.array, pandas.series, str): Optional
            pressure to use to convert gauge pressure to absolute. If a string is
            supplied, the corresponding data tag is assumed to be ambient pressure.
        ambient_pressure_unit (str): Optional ambient pressure unit, should be a
            unit recognized by pint.

    Returns:
        (pandas.DataFrame, dict)
    """
    # read file
    df = pd.read_csv(csv_file, parse_dates=True, index_col=0)
    # Drop empty columns
    df.drop(df.columns[df.columns.str.contains("Unnamed")], axis=1, inplace=True)
    df.rename(mapper=_strip, axis="columns", inplace=True)
    if rename_mapper:
        # Change tag names in some systematic way with the function rename_mapper
        df.rename(mapper=rename_mapper, axis="columns", inplace=True)
    metadata = {}
    if csv_file_metadata:
        with open(csv_file_metadata, "r") as f:
            reader = csv.reader(f)
            for line in reader:
                tag = line[0].strip()
                if rename_mapper:
                    tag = rename_mapper(tag)
                metadata[tag] = {
                    "reference_string": line[1].strip(),
                    "reference": None,
                    "description": line[2].strip(),
                    "units": line[3].strip(),
                }
    # If a model was provided, map the tags with a reference string to the model
    if model:
        update_metadata_model_references(model, metadata)
    # Drop the columns with no metadata (assuming those are columns to ignore)
    for tag in df:
        if tag not in metadata:
            df.drop(tag, axis=1, inplace=True)

    # Check if a data tag was specified to use as ambient pressure in conversion
    # of gauge pressures.  If so, get the numbers and replace the tag string
    if isinstance(ambient_pressure, str):
        try:
            ambient_pressure = np.array(df[ambient_pressure])
        except KeyError:
            _log.exception(
                "Tag '{}' does not exist for ambient pressure".format(ambient_pressure)
            )
            raise

    # If unit_system is specified bulk convert everything to that system of units
    # also update the meta data
    if unit_system:
        for tag in df:
            df[tag], metadata[tag]["units"] = unit_convert(
                df[tag],
                metadata[tag]["units"],
                system=unit_system,
                ambient_pressure=ambient_pressure,
                ambient_pressure_unit=ambient_pressure_unit,
            )

    return df, metadata


def _bin_number(x, bin_size):
    return np.array(x / bin_size, dtype=int)


def bin_data(df, bin_by, bin_no, bin_nom, bin_size, min_value=None, max_value=None):
    """
    Sort data into bins by a column value.  If the min or max are given and
    the value in bin_by for a row is out of the range [min, max], the row is
    dropped from the data frame.

    Args:
        df (pandas.DataFrame): Data frame to add bin information to
        bin_by (str): A column for values to bin by
        bin_no (str): A new column for bin number
        bin_nom (str): A new column for the mid-point value of bin_by
        bin_size (float): size of a bin
        min_value (in {float, None}): Smallest value to keep or None for no lower
        max_value (in (float, None}): Largest value to keep or None for no upper

    Returns:
        (dict): returns the data frame, and a dictionary with the number of rows
            in each bin.
    """

    # Drop rows outside [min, max]
    df.drop(index=df.index[np.isnan(df[bin_by])], inplace=True)
    if min_value is not None:
        df.drop(index=df.index[df[bin_by] < min_value], inplace=True)
    else:
        min_value = min(df[bin_by])
    if max_value is not None:
        df.drop(index=df.index[df[bin_by] > max_value], inplace=True)

    # Want the bins to line up so 0 is between bins and want min_value in bin 0.
    bin_offset = _bin_number(min_value, bin_size)
    df[bin_no] = _bin_number(df[bin_by], bin_size) - bin_offset
    df[bin_nom] = bin_size * (df[bin_no] + bin_offset + 0.5)
    a, b = np.unique(df[bin_no], return_counts=True)
    hist = dict(zip(a, b))
    return hist


def bin_stdev(df, bin_no, min_data=4):
    """
    Calculate the standard deviation for each column in each bin.

    Args:
        df (pandas.DataFrame): pandas data frame that is a bin number column
        bin_no (str): Column to group by, usually contains bin number
        min_data (int): Minimum number of data points requitred to calculate
            standard deviation for a bin (default=4)

    Returns:
        dict: key is the bin number and the value is a pandas.Serries with column
            standard deviations
    """
    nos = np.unique(df[bin_no])
    res = {}
    for i in nos:
        idx = df.index[df[bin_no] == i]
        if len(idx) < min_data:
            continue
        df2 = df.loc[idx]
        res[i] = df2.std(axis=0)
    return res


def data_rec_plot_book(
    df_data,
    df_rec,
    bin_nom,
    file="data_rec_plot_book.pdf",
    tmp_dir="tmp_plots",
    xlabel=None,
    metadata=None,
    cols=None,
    skip_cols=[],
):
    """
    Make box and whisker plots from process data compared to data rec results
    based on bins from the bin_data() function.  The df_data and df_rec data
    frames should have the same index set and the df_data data frame contains
    the bin data.  This will plot the intersection of columns containg numerical
    data.

    Args:
        df_data: data frame with original data
        df_rec: data frame with reconciled data
        bin_nom: bin mid-point value column
        file: path for generated pdf
        tmp_dir: a directory to store temporary plots in
        xlabel: Label for x axis
        metadata: tag meta data dictionary
        cols: List of columns to plot, if None plot all
        skip_cols: List of columns not to plot, this overrides cols

    Return:
        None

    """
    if sns is None:
        _log.error(
            "Plotting data requires the 'seaborn' and 'PyPDF2' packages. "
            "Install the required packages before using the data_book() function. "
            "Plot terminated."
        )
        return

    if not os.path.isdir(tmp_dir):
        os.mkdir(tmp_dir)

    pdfs = []
    flierprops = dict(markerfacecolor="0.5", markersize=2, marker="o", linestyle="none")
    f = plt.figure(figsize=(16, 9))
    if cols is None:
        cols = df_data.columns

    cols = sorted(set(cols).intersection(df_rec.columns))

    f = plt.figure(figsize=(16, 9))
    ax = sns.countplot(x=bin_nom, data=df_data)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    fname = os.path.join(tmp_dir, "plot_hist.pdf")
    f.savefig(fname, bbox_inches="tight")
    pdfs.append(fname)
    plt.close(f)

    for i, col in enumerate(cols):
        if col in skip_cols:
            continue
        f = plt.figure(i, figsize=(16, 9))

        x = pd.concat([df_data[bin_nom], df_data[bin_nom]], ignore_index=True)
        y = pd.concat([df_data[col], df_rec[col]], ignore_index=True)
        h = ["Data"] * len(df_data.index) + ["Reconciled"] * len(df_data.index)
        ax = sns.boxplot(x=x, y=y, hue=h, flierprops=flierprops)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        if metadata is not None:
            md = metadata.get(col, {})
            yl = "{} {} [{}]".format(
                col, md.get("description", ""), md.get("units", "none")
            )
        else:
            yl = col
        fname = os.path.join(tmp_dir, f"plot_{i}.pdf")
        ax.set(xlabel=xlabel, ylabel=yl)
        f.savefig(fname, bbox_inches="tight")
        pdfs.append(fname)
        plt.close(f)

    # Combine pdfs into one multi-page document
    writer = PdfFileMerger()
    for pdf in pdfs:
        writer.append(pdf)
    writer.write(file)
    _log.info(f"Plot written to {file}.")


def data_plot_book(
    df,
    bin_nom,
    file="data_plot_book.pdf",
    tmp_dir="tmp_plots",
    xlabel=None,
    metadata=None,
    cols=None,
    skip_cols=[],
):
    """
    Make box and whisker plots from process data based on bins from the
    bin_data() function.

    Args:
        df: data frame
        bin_nom: bin mid-point value column
        file: path for generated pdf
        tmp_dir: a directory to store temporary plots in
        xlabel: Label for x axis
        metadata: tag meta data dictionary

    Return:
        None

    """
    if sns is None:
        _log.error(
            "Plotting data requires the 'seaborn' and 'PyPDF2' packages. "
            "Install the required packages before using the data_book() function. "
            "Plot terminated."
        )
        return

    if not os.path.isdir(tmp_dir):
        os.mkdir(tmp_dir)

    pdfs = []
    flierprops = dict(markerfacecolor="0.5", markersize=2, marker="o", linestyle="none")
    f = plt.figure(figsize=(16, 9))
    if cols is None:
        cols = sorted(df.columns)
    else:
        cols = sorted(cols)

    f = plt.figure(figsize=(16, 9))
    ax = sns.countplot(x=bin_nom, data=df)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    fname = os.path.join(tmp_dir, "plot_hist.pdf")
    f.savefig(fname, bbox_inches="tight")
    pdfs.append(fname)
    plt.close(f)

    for i, col in enumerate(cols):
        if col in skip_cols:
            continue
        f = plt.figure(i, figsize=(16, 9))
        ax = sns.boxplot(x=df[bin_nom], y=df[col], flierprops=flierprops)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        if metadata is not None:
            md = metadata.get(col, {})
            yl = "{} {} [{}]".format(
                col, md.get("description", ""), md.get("units", "none")
            )
        else:
            yl = col
        fname = os.path.join(tmp_dir, f"plot_{i}.pdf")
        ax.set(xlabel=xlabel, ylabel=yl)
        f.savefig(fname, bbox_inches="tight")
        pdfs.append(fname)
        plt.close(f)

    # Combine pdfs into one multi-page document
    writer = PdfFileMerger()
    for pdf in pdfs:
        writer.append(pdf)
    writer.write(file)
    _log.info(f"Plot written to {file}.")
