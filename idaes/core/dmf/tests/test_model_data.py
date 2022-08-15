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
This module contains data utility function tests.
"""
import os
import warnings
import pytest
import numpy as np
import pyomo.environ as pyo
import pandas as pd

from pyomo.common.fileutils import this_file_dir

import idaes.core.dmf.model_data as da

_data_dir = os.path.join(this_file_dir(), "data_files")


@pytest.mark.unit
def test_bin_data():
    def make_data_frame():
        return pd.DataFrame(
            data={
                "power": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, None],
                "x1": [2, 3, 4, 6, 7, 9, 3, 1, 3, 3, 2, 3, 4],
                "x2": [4, 5, 6, 3, 5, 8, 4, 1, 4, 5, 2, 6, 5],
            }
        )

    df = make_data_frame()
    hist = da.bin_data(
        df,
        bin_by="power",
        bin_no="pb no.",
        bin_nom="pb power",
        bin_size=4,
        min_value=2,
        max_value=10,
    )
    assert hist[0] == 2
    assert hist[1] == 4
    assert hist[2] == 3
    assert 3 not in hist
    assert "pb no." in df
    assert "pb power" in df

    df = make_data_frame()
    hist = da.bin_data(
        df, bin_by="power", bin_no="pb no.", bin_nom="pb power", bin_size=4
    )
    assert hist[0] == 3
    assert hist[1] == 4
    assert hist[2] == 4
    assert hist[3] == 1
    assert "pb no." in df
    assert "pb power" in df

    r = da.bin_stdev(df, bin_no="pb no.", min_data=3)

    assert 0 in r
    assert 1 in r
    assert 2 in r
    assert 3 not in r  # not enough points

    assert r[0]["x1"] == pytest.approx(1.0)


@pytest.mark.component
def test_map_data():
    data1 = os.path.join(_data_dir, "data1.csv")
    data1_meta = os.path.join(_data_dir, "data1_meta.csv")
    m = pyo.ConcreteModel("Data Test Model")
    m.time = pyo.Set(initialize=[1, 2, 3])
    m.pressure = pyo.Var(m.time, doc="pressure (Pa)", initialize=101325)
    m.temperature = pyo.Var(m.time, doc="temperature (K)", initialize=300)
    m.volume = pyo.Var(m.time, doc="volume (m^3)", initialize=10)

    def retag(tag):
        return tag.replace(".junk", "")

    df, df_meta = da.read_data(
        data1, data1_meta, model=m, rename_mapper=retag, unit_system="mks"
    )

    # Check for expected columns in data and meta data
    assert "T" in df
    assert "P" in df
    assert "V" in df
    assert "T" in df_meta
    assert "P" in df_meta
    assert "V" in df_meta

    # Check that the unit strings updated after conversion
    assert df_meta["T"]["units"] == "kelvin"
    # this next unit is Pa
    assert df_meta["P"]["units"] == "kilogram / meter / second ** 2"
    assert df_meta["V"]["units"] == "meter ** 3"

    # Check that the unit conversions are okay
    assert df["T"]["1901-3-3 12:00"] == pytest.approx(300, rel=1e-4)
    assert df["P"]["1901-3-3 12:00"] == pytest.approx(200000, rel=1e-4)
    assert df["V"]["1901-3-3 12:00"] == pytest.approx(5.187286689, rel=1e-4)

    # Check the mapping of the tags to the model (the 1 key is the time indexed
    # from the model, because the reference is for a time-indexed variable)
    assert pyo.value(df_meta["T"]["reference"][1]) == pytest.approx(300, rel=1e-4)
    assert pyo.value(df_meta["P"]["reference"][1]) == pytest.approx(101325, rel=1e-4)
    assert pyo.value(df_meta["V"]["reference"][1]) == pytest.approx(10, rel=1e-4)


@pytest.mark.component
def test_map_data_use_ambient_pressure():
    # Try out PSIG and barometric pressure in PSIA.
    data1 = os.path.join(_data_dir, "data1.csv")
    data1_meta = os.path.join(_data_dir, "data1_meta.csv")
    m = pyo.ConcreteModel("Data Test Model")
    m.time = pyo.Set(initialize=[1, 2, 3])
    m.pressure = pyo.Var(m.time, doc="pressure (Pa)", initialize=101325)
    m.temperature = pyo.Var(m.time, doc="temperature (K)", initialize=300)
    m.volume = pyo.Var(m.time, doc="volume (m^3)", initialize=10)

    def retag(tag):
        return tag.replace(".junk", "")

    df, df_meta = da.read_data(
        data1,
        data1_meta,
        model=m,
        rename_mapper=retag,
        unit_system="mks",
        ambient_pressure="Pamb",
        ambient_pressure_unit="psi",
    )

    # Check that the unit conversions are okay
    assert df["P"]["1901-3-3 12:00"] == pytest.approx(195886, rel=1e-4)


@pytest.mark.component
def test_map_data_use_ambient_pressure2():
    # Try out inches of water column and barometric pressure in inHg.
    data2 = os.path.join(_data_dir, "data2.csv")
    data2_meta = os.path.join(_data_dir, "data2_meta.csv")
    m = pyo.ConcreteModel("Data Test Model")
    m.time = pyo.Set(initialize=[1, 2, 3])
    m.pressure = pyo.Var(m.time, doc="pressure (Pa)", initialize=101325)
    m.temperature = pyo.Var(m.time, doc="temperature (K)", initialize=300)
    m.volume = pyo.Var(m.time, doc="volume (m^3)", initialize=10)

    def retag(tag):
        return tag.replace(".junk", "")

    df, df_meta = da.read_data(
        data2,
        data2_meta,
        model=m,
        rename_mapper=retag,
        unit_system="mks",
        ambient_pressure="Pamb",
        ambient_pressure_unit="inHg",
    )

    # Check that the unit conversions are okay, same pressures as data1
    # differnt units, so the data read in should be the same
    assert df["P"]["1901-3-3 10:00"] == pytest.approx(96526.6, rel=1e-2)
    assert df["P"]["1901-3-3 12:00"] == pytest.approx(195886, rel=1e-2)


@pytest.mark.component
def test_unit_coversion():
    # spot test some unit conversions and features
    # da.unit_convert(x, frm, to=None, system=None, unit_string_map={},
    #                 ignore_units=[], gauge_pressures={}, atm=1.0):

    p_atm = np.array([1, 2, 3])
    p_psi, unit = da.unit_convert(p_atm, "atm", "psi")

    assert p_psi[0] == pytest.approx(14.7, rel=1e-2)
    assert p_psi[1] == pytest.approx(14.7 * 2, rel=1e-2)
    assert p_psi[2] == pytest.approx(14.7 * 3, rel=1e-2)
    assert unit == "pound_force_per_square_inch"

    # ppb is on the list of units to ignore, and not attempt to convert
    p, unit = da.unit_convert(p_atm, "ppb", "psi")
    assert p[0] == pytest.approx(1, rel=1e-2)
    assert unit == "ppb"

    # psig is on the list of gauge pressures.
    p, unit = da.unit_convert(p_psi, "psig", "atm")

    assert p[0] == pytest.approx(2, rel=1e-1)
    assert p[1] == pytest.approx(3, rel=1e-1)
    assert p[2] == pytest.approx(4, rel=1e-1)

    # check the general system of units conversion
    p_pa, unit = da.unit_convert(p_psi, "psi", system="mks")

    assert p_pa[0] == pytest.approx(101325, rel=1e-1)
    assert unit == "kilogram / meter / second ** 2"  # AKA Pa

    # Test for unit conversion of gauge pressue with different atmosperic
    # pressure values
    p, unit = da.unit_convert(
        p_psi, "psig", "atm", ambient_pressure=np.array([1, 1.1, 1.2])
    )

    assert p[0] == pytest.approx(2, rel=1e-1)
    assert p[1] == pytest.approx(3.1, rel=1e-1)
    assert p[2] == pytest.approx(4.2, rel=1e-1)

    # Again but make sure it works with a scalar to
    p, unit = da.unit_convert(p_psi, "psig", "atm", ambient_pressure=1.2)

    assert p[0] == pytest.approx(2.2, rel=1e-1)
    assert p[1] == pytest.approx(3.2, rel=1e-1)
    assert p[2] == pytest.approx(4.2, rel=1e-1)

    # test custom unit string mapping
    p, unit = da.unit_convert(
        p_psi, "MYPRESSURE", "atm", unit_string_map={"MYPRESSURE": "psi"}
    )

    assert p[0] == pytest.approx(1, rel=1e-1)
    assert p[1] == pytest.approx(2, rel=1e-1)
    assert p[2] == pytest.approx(3, rel=1e-1)

    # Test that a unit that doesn't exist remains unchanged
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        p, unit = da.unit_convert(p_psi, "MYPRESSURE", "atm")
        assert len(w) > 0
        found = False
        for wa in w:
            if (
                issubclass(wa.category, UserWarning)
                and str(wa.message) == "In unit conversion, from unit 'MYPRESSURE'"
                " is not defined. No conversion."
            ):
                found = True
                break
        if not found:
            raise Exception("Expected warning about undefined unit not found.")

    assert p_psi[0] == pytest.approx(14.7, rel=1e-1)
    assert unit == "MYPRESSURE"
