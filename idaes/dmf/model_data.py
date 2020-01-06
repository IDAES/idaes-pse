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
This module contains functions to read and manage data for use in parameter
esitmation, data reconciliation, and validation.
"""

__author__ = "John Eslick"

import logging
import csv
import pandas as pd
import numpy as np
import pint

import pyomo.environ as pyo
import warnings


def _strip(tag):
    """
    Tag renaming function to remove whitespace, depending on the csv format
    column heading items can pick up some extra whitespace when read into Pandas
    """
    return tag.strip()


_log = logging.getLogger(__file__)

# Some common unit string conversions, these are ones we've come across that
# are not handled by pint. We can organize and refine known unit strings
# more in the future.
_unit_strings = {
    # Pressure
    "PSI": "psi",
    "PSIA": "psi",
    "psia": "psi",
    "PSIG": "psig",
    "INWC": "in water",
    "IN WC": "in water",
    "IN/WC": "in water",
    '" H2O': "in water",
    "INHG": "in hg",
    "IN HG": "in hg",
    "IN/HG": "in hg",
    "HGA": "in hg",
    "IN HGA": "in hg",
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

_gauge_pressures = {"psig": "psi", "in water gauge": "in water", "in hg gauge": "in hg"}

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

    Args:
        csv_file (str): Path of file to read
        csv_file_metadata (str): Path of csv file to read column metadata from
        model (ConcreteModel): Optional model to map tags to
        rename_mapper (function): Optional function to rename tags
        unit_system (str): Optional system of units to atempt convert to
        ambient_pressure (float, numpy.array, pandas.series, str): Optional
            pressure to use to convert gauge pressure to absolute if a string is
            supplied the corresponding data tag is assumed to be ambient pressure
        ambient_pressure_unit (str): Optional ambient pressure unit, should be a
            unit recognized by pint.

    Returns:
        (DataFrame): A Pandas data frame with tags in columns and rows indexed
            by time.
        (dict): Column metadata, units of measure, description, and model
            mapping information.
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
        for tag, md in metadata.items():
            if md["reference_string"]:
                try:
                    md["reference"] = pyo.Reference(
                        eval(md["reference_string"], {"m": model})
                    )
                except KeyError:
                    warnings.warn(
                        "Tag reference {} not found".format(md["reference_string"]),
                        UserWarning,
                    )
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
