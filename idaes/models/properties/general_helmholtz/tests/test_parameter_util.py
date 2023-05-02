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
"""Test the Helmholtz EoS parameter writing utility.  Property calculations 
tested here do not use IDAES properties or the external functions. This directly
tests the expressions created by the WriteParameters class that will be written 
to parameter and expression files to verify they are generated correctly.

The tests her use three points from the saturation curve neat critical, near triple,
and in between.  This data comes from the original equation of state papers.  Comments
include citations.

On some test points, I had to infer extra significant figures for density.  For example,
this density of a liquid near the triple point is insensitive to pressure, so to calculate
the pressure from density of a liquid at the triple point, a lot of significant figures 
are needed. As long a I could add significant figures to the density and not change the 
reported value when rounded, I assumed it was okay. 
"""

import pytest
from idaes.models.properties.general_helmholtz.components.parameters.h2o import (
    main as h2o_main,
)
from idaes.models.properties.general_helmholtz.components.parameters.co2 import (
    main as co2_main,
)
from idaes.models.properties.general_helmholtz.components.parameters.propane import (
    main as propane_main,
)


def _common_1(sat_thermo_data, we):
    # Check the EoS expressions, the tolerance is a little loose
    # due to lack of sig. figs. in reported data
    for pnt in sat_thermo_data.values():
        assert we.calculate_pressure(rho=pnt["rhov"], T=pnt["T"]) == pytest.approx(
            pnt["p"], rel=1e-2, abs=1e-3
        )
        assert we.calculate_pressure(rho=pnt["rhol"], T=pnt["T"]) == pytest.approx(
            pnt["p"], rel=1e-2, abs=1e-3
        )
        assert we.calculate_enthalpy(rho=pnt["rhov"], T=pnt["T"]) == pytest.approx(
            pnt["hv"], rel=1e-2, abs=1e-3
        )
        assert we.calculate_enthalpy(rho=pnt["rhol"], T=pnt["T"]) == pytest.approx(
            pnt["hl"], rel=1e-2, abs=1e-3
        )
        assert we.calculate_entropy(rho=pnt["rhov"], T=pnt["T"]) == pytest.approx(
            pnt["sv"], rel=1e-2, abs=1e-3
        )
        assert we.calculate_entropy(rho=pnt["rhol"], T=pnt["T"]) == pytest.approx(
            pnt["sl"], rel=1e-2, abs=1e-3
        )

    # Check the approximate sat density curves
    rhol, rhov = we.approx_sat_curves([sat_thermo_data[2]["T"]])
    assert rhol[0] == pytest.approx(sat_thermo_data[2]["rhol"], rel=1e-1)
    assert rhov[0] == pytest.approx(sat_thermo_data[2]["rhov"], rel=1e-1)


@pytest.mark.unit
def test_h2o():
    # Some test data from:
    #
    # Wagner, W.,  A. Pruss (2002). The IAPWS Formulation 1995 for the
    #   Thermodynamic Properties of Ordinary Water Substance for General and
    #   Scientific Use. J. Phys. Chem. Ref. Data, 31, 387-535.
    sat_thermo_data = {
        1: {  # near critical
            "T": 273.16,
            "p": 0.612,
            "rhol": 999.79252,
            "hl": 0.001,
            "sl": 0.0000,
            "rhov": 0.00485,
            "hv": 2500.92,
            "sv": 9.1555,
        },
        2: {  # between critical and triple point
            "T": 400,
            "p": 245.77,
            "rhol": 937.486,
            "hl": 532.953,
            "sl": 1.6013,
            "rhov": 1.3694,
            "hv": 2715.70,
            "sv": 7.0581,
        },
        3: {  # near triple point
            "T": 646,
            "p": 21775,
            "rhol": 402.96,
            "hl": 1963.49,
            "sl": 4.2214,
            "rhov": 243.46,
            "hv": 2238.06,
            "sv": 4.6465,
        },
    }

    we = h2o_main(dry_run=True)
    _common_1(sat_thermo_data, we)


@pytest.mark.unit
def test_co2():
    # Some test data from:
    #
    # Span, R., and W. Wagner (1996). A New Equation of State for Carbon Dioxide
    #    Covering the Fluid Region from the Triple-Point Temperature to 1100 K as
    #    Pressures up to 800 MPa. Journal of Physical and Chemical Reference Data,
    #    25, 1509.
    sat_thermo_data = {
        1: {  # near critical
            "T": 216.592,
            "p": 517.96,
            "rhol": 1178.46,
            "hl": -426.74,
            "sl": -2.2177,
            "rhov": 13.761,
            "hv": -76.364,
            "sv": -0.59999,
        },
        2: {  # between critical and triple point
            "T": 250,
            "p": 1785.0,
            "rhol": 1045.97,
            "hl": -359.07,
            "sl": -1.9323,
            "rhov": 46.644,
            "hv": -69.376,
            "sv": -0.77492,
        },
        3: {  # near triple point
            "T": 303,
            "p": 7189.0,
            "rhol": 599.86,
            "hl": -203.73,
            "sl": -1.4004,
            "rhov": 339.00,
            "hv": -139.91,
            "sv": -1.1897,
        },
    }

    we = co2_main(dry_run=True)
    _common_1(sat_thermo_data, we)


@pytest.mark.unit
def test_propane():
    # Some test data from:
    #
    # Lemmon, E. , McLinden, M. and Wagner, W. (2009), Thermodynamic Properties of
    #   Propane. III. A Reference Equation of State for Temperatures from the Melting
    #   Line to 650 K and Pressures up to 1000 MPa, Journal of Chemical and Engineering
    #   Data, [online], https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=832438
    #   (Accessed April 14, 2023)
    sat_thermo_data = {
        1: {  #
            "T": -187.625 + 273.15,
            "p": 0.17203e-6,
            "rhol": 733.125204,
            "hl": -196.64,
            "sl": -1.396,
            "rhov": 0.1071e-7,
            "hv": 366.26,
            "sv": 5.186,
        },
        2: {  # between critical and triple point
            "T": 273.15,
            "p": 474.46,
            "rhol": 528.59,
            "hl": 200.00,
            "sl": 1.000,
            "rhov": 10.351,
            "hv": 574.87,
            "sv": 2.372,
        },
        3: {  # near triple point
            "T": 95 + 273.15,
            "p": 4119.5,
            "rhol": 286.51,
            "hl": 516.33,
            "sl": 1.948,
            "rhov": 156.31,
            "hv": 595.81,
            "sv": 2.164,
        },
    }

    we = propane_main(dry_run=True)
    _common_1(sat_thermo_data, we)
