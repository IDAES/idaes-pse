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

The tests here use three points from the saturation curve near critical, near triple,
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
from idaes.models.properties.general_helmholtz.components.parameters.r1234ze import (
    main as r1234ze_main,
)
from idaes.models.properties.general_helmholtz.components.parameters.r125 import (
    main as r125_main,
)
from idaes.models.properties.general_helmholtz.components.parameters.r134a import (
    main as r134a_main,
)
from idaes.models.properties.general_helmholtz.components.parameters.r227ea import (
    main as r227ea_main,
)
from idaes.models.properties.general_helmholtz.components.parameters.r32 import (
    main as r32_main,
)


def _common_sat(sat_thermo_data, we):
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


def _common_pressure(thermo_data, we):
    # Check the EoS expressions, the tolerance is a little loose
    # due to lack of sig. figs. in reported data
    for pnt in thermo_data.values():
        assert we.calculate_pressure(rho=pnt["rho"], T=pnt["T"]) == pytest.approx(
            pnt["p"], rel=1e-2, abs=1e-3
        )


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
    _common_sat(sat_thermo_data, we)


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
    _common_sat(sat_thermo_data, we)


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
        1: {  # near triple
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
        3: {  # near critical point
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
    _common_sat(sat_thermo_data, we)


@pytest.mark.unit
def test_r1234ze():
    # Some test data from:
    #
    # Monika Thol and Eric W. Lemmon. Equation of State for the Thermodynamic
    #    Properties of trans-1,3,3,3-Tetrafluoropropene [R-1234ze(E)]. Int. J.
    #    Thermophys, 37(3):1–16, 2016. doi:10.1007/s10765-016-2040-6.
    thermo_data = {
        1: {
            "T": 200,
            "p": 2161.4,
            "rho": 12.6 * 114.0416,
        },
        2: {
            "T": 200,
            "p": 0,
            "rho": 1e-12,  # Need a small number here to avoid 0/0
        },
        3: {
            "T": 350,
            "p": 95041.56,
            "rho": 11.4 * 114.0416,
        },
        4: {
            "T": 360,
            "p": 2039.103,
            "rho": 1.0 * 114.0416,
        },
        4: {
            "T": 383,
            "p": 3670.105,
            "rho": 4.29 * 114.0416,
        },
        5: {
            "T": 420,
            "p": 19954.47,
            "rho": 8.0 * 114.0416,
        },
    }

    we = r1234ze_main(dry_run=True)
    _common_pressure(thermo_data, we)


@pytest.mark.unit
def test_r125():
    # Some test data from:
    #
    # Lemmon, E.W., R.T. Jacobsen. A New Functional Form and New Fitting Techniques
    #    for Equations of State with Application to Pentafluoroethane (HFC-125),
    #    J. Phys. Chem. Ref. Data, Vol. 34, No. 1, 2005
    sat_thermo_data = {
        1: {  # near triple point
            "T": -100.63 + 273.15,
            "p": 2.91,
            "rhol": 1690.6804,
            "hl": 87.130,
            "sl": 0.49022,
            "rhov": 0.24462,
            "hv": 277.39,
            "sv": 1.5931,
        },
        2: {  # between critical and triple point
            "T": 273.15,
            "p": 670.52,
            "rhol": 1319.8,
            "hl": 200.0,
            "sl": 1.000,
            "rhov": 42.070,
            "hv": 333.16,
            "sv": 1.4875,
        },
        3: {  # near critical point
            "T": 65 + 273.15,
            "p": 3536.97,
            "rhol": 735.11,
            "hl": 304.88,
            "sl": 1.3311,
            "rhov": 416.57,
            "hv": 332.24,
            "sv": 1.412,
        },
    }

    we = r125_main(dry_run=True)
    _common_sat(sat_thermo_data, we)


@pytest.mark.unit
def test_r134a():
    # Some test data from:
    #
    # Tillner-Roth, R.; Baehr, H.D., An International Standard Formulation for the",
    #    Thermodynamic Properties of 1,1,1,2-Tetrafluoroethane (HFC-134a) for",
    #    Temperatures from 170 K to 455 K and Pressures up to 70 MPa, J. Phys. Chem.",
    #    Ref. Data, 1994, 23, 5, 657-729, https://doi.org/10.1063/1.555958"
    sat_thermo_data = {
        1: {  # near triple point
            "T": 169.85,
            "p": 0.39,
            "rhol": 1591.10745,
            "hl": 71.454,
            "sl": 0.4126,
            "rhov": 0.02817,
            "hv": 334.94,
            "sv": 1.9639,
        },
        2: {  # between critical and triple point
            "T": 310.0,
            "p": 933.40,
            "rhol": 1159.9,
            "hl": 251.73,
            "sl": 1.1756,
            "rhov": 45.785,
            "hv": 418.03,
            "sv": 1.7121,
        },
        3: {  # near critical point
            "T": 372,
            "p": 3881.1,
            "rhol": 693.10,
            "hl": 367.71,
            "sl": 1.5041,
            "rhov": 334.83,
            "hv": 412.67,
            "sv": 1.6250,
        },
    }

    we = r134a_main(dry_run=True)
    _common_sat(sat_thermo_data, we)


@pytest.mark.unit
def test_r32():
    # Some test data from:
    #
    # Tillner-Roth, R., A. Yokozeki. An International Standard Equation of State for",
    #    Difluoromethane (R-32) for Temperatures from the Triple Point at 136.34 K to",
    #    435 K and Pressures up to 70 MPa, Journal of Physical and Chemical Reference",
    #    Data 26, 1273 (1997); https://doi.org/10.1063/1.556002"
    sat_thermo_data = {
        1: {  # near triple point
            "T": -136.81 + 273.15,
            "p": 0.05446079,
            "rhol": 1429.2733,
            "hl": -19.07,
            "sl": -0.104,
            "rhov": 0.0025,
            "hv": 444.31,
            "sv": 3.2937,
        },
        2: {  # between critical and triple point
            "T": 273.15,
            "p": 813.10,
            "rhol": 1055.25,
            "hl": 200.00,
            "sl": 1.0,
            "rhov": 22.091,
            "hv": 515.30,
            "sv": 2.1543,
        },
        3: {  # near critical point
            "T": 76 + 273.15,
            "p": 5531.5,
            "rhol": 583.33,
            "hl": 378.03,
            "sl": 1.5470,
            "rhov": 275.00,
            "hv": 455.86,
            "sv": 1.7699,
        },
    }

    we = r32_main(dry_run=True)
    _common_sat(sat_thermo_data, we)
