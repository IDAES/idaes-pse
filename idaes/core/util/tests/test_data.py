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
This module contains data utility function tests.
"""
import os
import pytest
import pyomo.environ as pyo

from pyomo.common.fileutils import this_file_dir

import idaes.core.util.data as da

_data_dir = os.path.join(this_file_dir(), "data_files")

def test_map_data():
    data1 = os.path.join(_data_dir, "data1.csv")
    data1_meta = os.path.join(_data_dir, "data1_meta.csv")
    m = pyo.ConcreteModel("Data Test Model")
    m.time = pyo.Set(initialize=[1,2,3])
    m.pressure = pyo.Var(m.time, doc="pressure (Pa)", initialize=101325)
    m.temperature = pyo.Var(m.time, doc="temperature (K)", initialize=300)
    m.volume = pyo.Var(m.time, doc="volume (m^3)", initialize=10)

    def retag(tag):
        return tag.replace(".junk", "")

    df, df_meta = da.read_data(data1, data1_meta, model=m, rename_mapper=retag,
                               unit_system="mks")

    # Check for expected columns in data and meta data
    assert("T" in df)
    assert("P" in df)
    assert("V" in df)
    assert("T" in df_meta)
    assert("P" in df_meta)
    assert("V" in df_meta)

    # Check that the unit stings updated after conversion
    assert(df_meta["T"]["units"] == "kelvin")
    # this next unit is Pa
    assert(df_meta["P"]["units"] == "kilogram / meter / second ** 2")
    assert(df_meta["V"]["units"] == "meter ** 3")

    # Check that the unit conversions are okay
    assert(df["T"]["1901-3-3 12:00"] == pytest.approx(300, rel=1e-4))
    assert(df["P"]["1901-3-3 12:00"] == pytest.approx(200000, rel=1e-4))
    assert(df["V"]["1901-3-3 12:00"] == pytest.approx(5.187286689, rel=1e-4))

    # Check the mapping of the tags to the model (the 1 key is the time indexed
    # from the model, because the reference is for a time-indexed variable)
    assert(pyo.value(df_meta["T"]["reference"][1]) ==
        pytest.approx(300, rel=1e-4))
    assert(pyo.value(df_meta["P"]["reference"][1]) ==
        pytest.approx(101325, rel=1e-4))
    assert(pyo.value(df_meta["V"]["reference"][1]) ==
        pytest.approx(10, rel=1e-4))
