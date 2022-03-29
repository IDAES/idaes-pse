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
Tests for the cubic root finder external functions
"""

import pytest
import os

from numpy import logspace

from pyomo.environ import (ConcreteModel,
                           ExternalFunction,
                           value)
from pyomo.core.base.external import AMPLExternalFunction

from idaes import bin_directory

__author__ = "Andrew Lee"

# Set path to root finder .so file
_so = os.path.join(bin_directory, "cubic_roots.so")

# Set module level pyest marker
pytestmark = pytest.mark.cubic_root
prop_available = os.path.isfile(_so)

# Define parameters for different supported equations of state
EoS_param = {
        0: {'u': 2, 'w': -1},
        1: {'u': 1, 'w': 0}
        }

# Set paramter for number of points to test in reduced pressure and temperature
# Use Tr and Pr as A and B are linked
SAMPLES = 250
t_set = logspace(-1, 2, num=SAMPLES, base=10, endpoint=True)
p_set = logspace(-2, 2, num=SAMPLES, base=10, endpoint=True)

# Absolute tolerance for root finder checks
TOL = 1e-5
# Parameter indicating whether to test partial derivatives
TEST_DERS = True
# Relative finite difference step for partial derivatives
DEL = 1e-4
# Relative tolerance for accepting partial derivatives
FD_TOL = 5e-2


def between(y, x1, x2):
    return 0 > (y-x1)*(y-x2)


@pytest.fixture()
def root_finder():
    m = ConcreteModel()

    # Define external function methods
    m.proc_Z_liq = ExternalFunction(library=_so,
                                    function="ceos_z_liq")
    m.proc_Z_vap = ExternalFunction(library=_so,
                                    function="ceos_z_vap")
    m.proc_Z_liq_x = ExternalFunction(library=_so,
                                      function="ceos_z_liq_extend")
    m.proc_Z_vap_x = ExternalFunction(library=_so,
                                      function="ceos_z_vap_extend")

    return m


# TODO - Need tests for how external function behaves when given an invalid
# TODO - eos type. Currently seems to return a value, which probably should not
# TODO - happen.


@pytest.mark.integration
@pytest.mark.skipif(not prop_available,
                    reason="Cubic root finder not available")
def test_roots_Z_liq(root_finder):
    for eos_type in [0, 1]:
        u = EoS_param[eos_type]["u"]
        w = EoS_param[eos_type]["w"]

        for T in t_set:
            for P in p_set:
                # Calculate A and B parameters from Tr and Pr
                A = 0.5*P/T**2
                B = 0.1*P/T

                # Get results of external function call
                f = root_finder.proc_Z_liq
                assert(isinstance(f, AMPLExternalFunction))
                Z, g, h = f.evaluate_fgh(args=(eos_type, A, B))

                # Calculate parameters of cubic
                c1 = 1
                c2 = -(1+B-u*B)
                c3 = (A-u*B-(u-w)*B**2)
                c4 = -(A*B+w*B**2+w*B**3)

                # Calculate residual and derivatives w.r.t. Z
                res = c1*Z**3 + c2*Z**2 + c3*Z + c4
                dz = 3*c1*Z**2 + 2*c2*Z + c3
                dz2 = 6*c1*Z + 2*c2

                try:
                    # Residual can be extemely sensitive to value of Z, so
                    # if residual fails tolerance, test size of Newton step to
                    # converge to root
                    assert (pytest.approx(0, abs=TOL) == res or
                            pytest.approx(0, abs=TOL) == value(res/dz))

                    # Check derivative signs to confirm correct root
                    assert dz >= 0  # Should always have non-negative slope

                    # Determine number of roots - calculate discriminant
                    dis = (18*c1*c2*c3*c4 - 4*c2**3*c4 + c2**2*c3**2 -
                           4*c1*c3**3 - 27*c1**2*c4**2)

                    if dis >= 0:
                        # Cubic has 2 or 3 real roots
                        # Second derivative should be non-positive
                        assert dz2 <= 0
                    # otherwise no need to check 2nd derivative

                    if TEST_DERS:
                        # Perform finite differences on A and B
                        ZAp, gAp, hAp = f.evaluate_fgh(
                                args=(eos_type, A*(1+DEL), B))
                        ZAm, gAm, hAm = f.evaluate_fgh(
                                args=(eos_type, A*(1-DEL), B))

                        ZBp, gBp, hBp = f.evaluate_fgh(
                                    args=(eos_type, A, B*(1+DEL)))
                        ZBm, gBm, hBm = f.evaluate_fgh(
                                    args=(eos_type, A, B*(1-DEL)))

                        # Check variance in Z values. A very large difference
                        # indicates a transition between single and multiple
                        # root regions, and hence that the partial derivatvies
                        # will be very sensitive.
                        # In these cases, skip the derivative tests.
                        if abs(ZAp - Z) > 1e-3 or abs(ZAm - Z) > 1e-3:
                            A_skip = True
                        else:
                            A_skip = False
                        if (abs(ZBp - Z) > 1e-3 or
                                abs(ZBm - Z) > 1e-3 or
                                abs(dis) < 1e-7):
                            B_skip = True
                        else:
                            B_skip = False

                        # Test gradient terms
                        # Calculate numerical first partial derivative
                        if not A_skip:
                            dZdA_p = (ZAp-Z)/(A*DEL)
                            dZdA_m = (Z-ZAm)/(A*DEL)

                        if not B_skip:
                            dZdB_p = (ZBp-Z)/(B*DEL)
                            dZdB_m = (Z-ZBm)/(B*DEL)

                        # Partial derivative w.r.t. EoS identifier
                        assert g[0] == 0

                        # Check that external function value lies within TOL of
                        # least one of the numerical values (delta+ or delta-),
                        # OR lies between the two numerical values and within
                        # 10*TOL of one of the numerical values
                        if not A_skip:
                            assert (
                                pytest.approx(dZdA_p, FD_TOL) == g[1] or
                                pytest.approx(dZdA_m, FD_TOL) == g[1] or
                                (between(g[1], dZdA_p, dZdA_m) and (
                                 pytest.approx(dZdA_p, 10*FD_TOL) == g[1] or
                                 pytest.approx(dZdA_m, 10*FD_TOL) == g[1])))

                        if not B_skip:
                            assert (
                                pytest.approx(dZdB_p, FD_TOL) == g[2] or
                                pytest.approx(dZdB_m, FD_TOL) == g[2] or
                                (between(g[2], dZdB_p, dZdB_m) and (
                                 pytest.approx(dZdB_p, 10*FD_TOL) == g[2] or
                                 pytest.approx(dZdB_m, 10*FD_TOL) == g[2])))

                        # Test hessian terms
                        # Calculate numerical second partial derivatives
                        if not A_skip:
                            d2ZdA2_p = (gAp[1]-g[1])/(A*DEL)
                            d2ZdA2_m = (g[1]-gAm[1])/(A*DEL)

                        if not B_skip:
                            d2ZdB2_p = (gBp[2]-g[2])/(B*DEL)
                            d2ZdB2_m = (g[2]-gBm[2])/(B*DEL)

                        if not A_skip and not B_skip:
                            d2ZdAB_p = (gBp[1]-g[1])/(B*DEL)
                            d2ZdAB_m = (g[1]-gBm[1])/(B*DEL)

                        # Partial derivatives w.r.t eos_type
                        assert h[0] == 0
                        assert h[1] == 0
                        assert h[3] == 0

                        if not A_skip:
                            assert (pytest.approx(d2ZdA2_p, FD_TOL) == h[2] or
                                    pytest.approx(d2ZdA2_m, FD_TOL) == h[2] or
                                    between(h[2], d2ZdA2_p, d2ZdA2_m))
                        if not A_skip and not B_skip:
                            assert (pytest.approx(d2ZdAB_p, FD_TOL) == h[4] or
                                    pytest.approx(d2ZdAB_m, FD_TOL) == h[4] or
                                    between(h[4], d2ZdAB_p, d2ZdAB_m))
                        if not B_skip:
                            # Second derivative w.r.t. B is very sensitive near
                            # point that roots disappear, and is at a
                            # maximum (or minimum) so skip tests if close to
                            # this point
                            assert (pytest.approx(d2ZdB2_p, FD_TOL) == h[5] or
                                    pytest.approx(d2ZdB2_m, FD_TOL) == h[5] or
                                    between(h[5], d2ZdB2_p, d2ZdB2_m))

                except AssertionError:
                    # Print values at failure and raise exception
                    print(eos_type, T, P, A, B, Z)
                    raise


@pytest.mark.integration
@pytest.mark.skipif(not prop_available,
                    reason="Cubic root finder not available")
def test_roots_Z_vap(root_finder):
    for eos_type in [0, 1]:
        u = EoS_param[eos_type]["u"]
        w = EoS_param[eos_type]["w"]

        for T in t_set:
            for P in p_set:
                # Calculate A and B parameters from Tr and Pr
                A = 0.5*P/T**2
                B = 0.1*P/T

                # Get results of external function call
                f = root_finder.proc_Z_vap
                assert(isinstance(f, AMPLExternalFunction))
                Z, g, h = f.evaluate_fgh(args=(eos_type, A, B))

                # Calculate parameters of cubic
                c1 = 1
                c2 = -(1+B-u*B)
                c3 = (A-u*B-(u-w)*B**2)
                c4 = -(A*B+w*B**2+w*B**3)

                # Calculate residual and derivatives w.r.t. Z
                res = c1*Z**3 + c2*Z**2 + c3*Z + c4
                dz = 3*c1*Z**2 + 2*c2*Z + c3
                dz2 = 6*c1*Z + 2*c2

                try:
                    # Residual can be extemely sensitive to value of Z, so
                    # if residual fails tolerance, test size of Newton step to
                    # converge to root
                    assert (pytest.approx(0, abs=TOL) == res or
                            pytest.approx(0, abs=TOL) == value(res/dz))

                    # Check derivative signs to confirm correct root
                    assert dz >= 0  # Should always have non-negative slope

                    # Determine number of roots - calculate discriminant
                    dis = (18*c1*c2*c3*c4 - 4*c2**3*c4 + c2**2*c3**2 -
                           4*c1*c3**3 - 27*c1**2*c4**2)

                    if dis >= 0:
                        # Cubic has 2 or 3 real roots
                        # Second derivative should be non-negative
                        assert dz2 >= 0
                    # otherwise no need to check 2nd derivative

                    if TEST_DERS:
                        # Perform finite differences on A and B
                        ZAp, gAp, hAp = f.evaluate_fgh(
                                args=(eos_type, A*(1+DEL), B))
                        ZAm, gAm, hAm = f.evaluate_fgh(
                                args=(eos_type, A*(1-DEL), B))

                        ZBp, gBp, hBp = f.evaluate_fgh(
                                    args=(eos_type, A, B*(1+DEL)))
                        ZBm, gBm, hBm = f.evaluate_fgh(
                                    args=(eos_type, A, B*(1-DEL)))

                        # Check variance in Z values. A very large difference
                        # indicates a transition between single and multiple
                        # root regions, and hence that the partial derivatvies
                        # will be very sensitive.
                        # In these cases, skip the derivative tests.
                        if abs(ZAp - Z) > 1e-3 or abs(ZAm - Z) > 1e-3:
                            A_skip = True
                        else:
                            A_skip = False
                        if (abs(ZBp - Z) > 1e-3 or
                                abs(ZBm - Z) > 1e-3 or
                                abs(dis) < 1e-7):
                            B_skip = True
                        else:
                            B_skip = False

                        # Test gradient terms
                        # Calculate numerical first partial derivative
                        if not A_skip:
                            dZdA_p = (ZAp-Z)/(A*DEL)
                            dZdA_m = (Z-ZAm)/(A*DEL)

                        if not B_skip:
                            dZdB_p = (ZBp-Z)/(B*DEL)
                            dZdB_m = (Z-ZBm)/(B*DEL)

                        # Partial derivative w.r.t. EoS identifier
                        assert g[0] == 0

                        # Check that external function value lies within TOL of
                        # least one of the numerical values (delta+ or delta-),
                        # OR lies between the two numerical values and within
                        # 10*TOL of one of the numerical values
                        if not A_skip:
                            assert (
                                pytest.approx(dZdA_p, FD_TOL) == g[1] or
                                pytest.approx(dZdA_m, FD_TOL) == g[1] or
                                (between(g[1], dZdA_p, dZdA_m) and (
                                 pytest.approx(dZdA_p, 10*FD_TOL) == g[1] or
                                 pytest.approx(dZdA_m, 10*FD_TOL) == g[1])))

                        if not B_skip:
                            assert (
                                pytest.approx(dZdB_p, FD_TOL) == g[2] or
                                pytest.approx(dZdB_m, FD_TOL) == g[2] or
                                (between(g[2], dZdB_p, dZdB_m) and (
                                 pytest.approx(dZdB_p, 10*FD_TOL) == g[2] or
                                 pytest.approx(dZdB_m, 10*FD_TOL) == g[2])))

                        # Test hessian terms
                        # Calculate numerical second partial derivatives
                        if not A_skip:
                            d2ZdA2_p = (gAp[1]-g[1])/(A*DEL)
                            d2ZdA2_m = (g[1]-gAm[1])/(A*DEL)

                        if not B_skip:
                            d2ZdB2_p = (gBp[2]-g[2])/(B*DEL)
                            d2ZdB2_m = (g[2]-gBm[2])/(B*DEL)

                        if not A_skip and not B_skip:
                            d2ZdAB_p = (gBp[1]-g[1])/(B*DEL)
                            d2ZdAB_m = (g[1]-gBm[1])/(B*DEL)

                        # Partial derivatives w.r.t eos_type
                        assert h[0] == 0
                        assert h[1] == 0
                        assert h[3] == 0

                        if not A_skip:
                            assert (pytest.approx(d2ZdA2_p, FD_TOL) == h[2] or
                                    pytest.approx(d2ZdA2_m, FD_TOL) == h[2] or
                                    between(h[2], d2ZdA2_p, d2ZdA2_m))
                        if not A_skip and not B_skip:
                            assert (pytest.approx(d2ZdAB_p, FD_TOL) == h[4] or
                                    pytest.approx(d2ZdAB_m, FD_TOL) == h[4] or
                                    between(h[4], d2ZdAB_p, d2ZdAB_m))
                        if not B_skip:
                            # Second derivative w.r.t. B is very sensitive near
                            # point that roots disappear, and is at a
                            # maximum (or minimum) so skip tests if close to
                            # this point
                            assert (pytest.approx(d2ZdB2_p, FD_TOL) == h[5] or
                                    pytest.approx(d2ZdB2_m, FD_TOL) == h[5] or
                                    between(h[5], d2ZdB2_p, d2ZdB2_m))

                except AssertionError:
                    # Print values at failure and raise exception
                    print(eos_type, T, P, A, B, Z)
                    raise


@pytest.mark.integration
@pytest.mark.skipif(not prop_available,
                    reason="Cubic root finder not available")
def test_roots_Z_liq_ext(root_finder):
    for eos_type in [0, 1]:
        u = EoS_param[eos_type]["u"]
        w = EoS_param[eos_type]["w"]

        for T in t_set:
            for P in p_set:
                # Calculate A and B parameters from Tr and Pr
                A = 0.5*P/T**2
                B = 0.1*P/T

                # Get results of external function call
                f = root_finder.proc_Z_liq_x
                assert(isinstance(f, AMPLExternalFunction))
                Z, g, h = f.evaluate_fgh(args=(eos_type, A, B))

                # Calculate parameters of cubic
                c1 = 1
                c2 = -(1+B-u*B)
                c3 = (A-u*B-(u-w)*B**2)
                c4 = -(A*B+w*B**2+w*B**3)

                det = c2**2 - 3*c3
                a = -(1.0/3.0)*(c2 + det**0.5)

                # Check to see if extension is triggered
                if det <= 0 or (a**3 + c2*a**2 + c3*a + c4) >= 0:
                    # Extension is not used
                    # Calculate residual and derivatives w.r.t. Z
                    res = c1*Z**3 + c2*Z**2 + c3*Z + c4
                    dz = 3*c1*Z**2 + 2*c2*Z + c3
                    dz2 = 6*c1*Z + 2*c2

                    try:
                        # Residual can be extemely sensitive to value of Z, so
                        # if residual fails tolerance, test size of Newton step
                        # to converge to root
                        assert (pytest.approx(0, abs=TOL) == res or
                                pytest.approx(0, abs=TOL) == value(res/dz))

                        # Check derivative signs to confirm correct root
                        assert dz >= 0  # Should always have non-negative slope

                        # Determine number of roots - calculate discriminant
                        dis = (18*c1*c2*c3*c4 - 4*c2**3*c4 + c2**2*c3**2 -
                               4*c1*c3**3 - 27*c1**2*c4**2)

                        if dis >= 0:
                            # Cubic has 2 or 3 real roots
                            # Second derivative should be non-positive
                            assert dz2 <= 0
                        # otherwise no need to check 2nd derivative

                    except AssertionError:
                        # Print values at failure and raise exception
                        print(eos_type, T, P, A, B, Z, det, a)
                        raise

                else:
                    # Extention is used, calculate extended root
                    c1x = 2*a
                    c2x = -c2 - 3.0*c1x
                    c3x = 3*c1x**2 + 2*c2*c1x + c3
                    c4x = c4 - 0.75*c1x**3 - 0.5*c2*c1x**2

                    # Calculate residual and derivatives w.r.t. Z_ext
                    res = c1*Z**3 + c2x*Z**2 + c3x*Z + c4x
                    dz = 3*c1*Z**2 + 2*c2x*Z + c3x

                    try:
                        # Residual can be extemely sensitive to value of Z, so
                        # if residual fails tolerance, test size of Newton step
                        # to converge to root
                        assert (pytest.approx(0, abs=TOL) == res or
                                pytest.approx(0, abs=TOL) == value(res/dz))

                        # Check derivative signs to confirm correct root
                        assert dz >= 0  # Should always have non-negative slope

                        # Determine number of roots - calculate discriminant
                        dis = (18*c1*c2x*c3x*c4x - 4*c2x**3*c4x +
                               c2x**2*c3x**2 - 4*c1*c3x**3 - 27*c1**2*c4x**2)

                        # Second derivative could be anything, don't check

                    except AssertionError:
                        # Print values at failure and raise exception
                        print(eos_type, T, P, A, B, Z)
                        raise

                if TEST_DERS:
                    try:
                        # Perform finite differences on A and B
                        ZAp, gAp, hAp = f.evaluate_fgh(
                                args=(eos_type, A*(1+DEL), B))
                        ZAm, gAm, hAm = f.evaluate_fgh(
                                args=(eos_type, A*(1-DEL), B))

                        ZBp, gBp, hBp = f.evaluate_fgh(
                                    args=(eos_type, A, B*(1+DEL)))
                        ZBm, gBm, hBm = f.evaluate_fgh(
                                    args=(eos_type, A, B*(1-DEL)))

                        # Check variance in Z values. A very large
                        # difference indicates a transition between
                        # single and multiple root regions, and hence that
                        # the partial derivatvies will be very sensitive.
                        # In these cases, skip the derivative tests.
                        if (abs(ZAp - Z) > 1e-3 or
                                abs(ZAm - Z) > 1e-3 or
                                abs(a) < 1e-1):
                            A_skip = True
                        else:
                            A_skip = False
                        if (abs(ZBp - Z) > 1e-3 or
                                abs(ZBm - Z) > 1e-3 or
                                abs(dis) < 1e-7 or
                                abs(a) < 1e-1):
                            B_skip = True
                        else:
                            B_skip = False

                        # Test gradient terms
                        # Calculate numerical first partial derivative
                        if not A_skip:
                            dZdA_p = (ZAp-Z)/(A*DEL)
                            dZdA_m = (Z-ZAm)/(A*DEL)

                        if not B_skip:
                            dZdB_p = (ZBp-Z)/(B*DEL)
                            dZdB_m = (Z-ZBm)/(B*DEL)

                        # Partial derivative w.r.t. EoS identifier
                        assert g[0] == 0

                        # Check that external function value lies within
                        # TOL of least one of the numerical values (delta+
                        # or delta-), OR lies between the two numerical
                        # values and within 10*TOL of one of the numerical
                        # values
                        if not A_skip:
                            assert (
                                pytest.approx(dZdA_p, FD_TOL) == g[1] or
                                pytest.approx(dZdA_m, FD_TOL) == g[1] or
                                (between(g[1], dZdA_p, dZdA_m) and (
                                 pytest.approx(dZdA_p, 10*FD_TOL) ==
                                 g[1] or
                                 pytest.approx(dZdA_m, 10*FD_TOL) ==
                                 g[1])))

                        if not B_skip:
                            assert (
                                pytest.approx(dZdB_p, FD_TOL) == g[2] or
                                pytest.approx(dZdB_m, FD_TOL) == g[2] or
                                (between(g[2], dZdB_p, dZdB_m) and (
                                 pytest.approx(dZdB_p, 10*FD_TOL) ==
                                 g[2] or
                                 pytest.approx(dZdB_m, 10*FD_TOL) ==
                                 g[2])))

                        # Test hessian terms
                        # Calculate numerical second partial derivatives
                        if not A_skip:
                            d2ZdA2_p = (gAp[1]-g[1])/(A*DEL)
                            d2ZdA2_m = (g[1]-gAm[1])/(A*DEL)

                        if not B_skip:
                            d2ZdB2_p = (gBp[2]-g[2])/(B*DEL)
                            d2ZdB2_m = (g[2]-gBm[2])/(B*DEL)

                        if not A_skip and not B_skip:
                            d2ZdAB_p = (gBp[1]-g[1])/(B*DEL)
                            d2ZdAB_m = (g[1]-gBm[1])/(B*DEL)

                        # Partial derivatives w.r.t eos_type
                        assert h[0] == 0
                        assert h[1] == 0
                        assert h[3] == 0

                        if not A_skip:
                            assert (
                                pytest.approx(d2ZdA2_p, FD_TOL) == h[2] or
                                pytest.approx(d2ZdA2_m, FD_TOL) == h[2] or
                                between(h[2], d2ZdA2_p, d2ZdA2_m))
                        if not A_skip and not B_skip:
                            assert (
                                pytest.approx(d2ZdAB_p, FD_TOL) == h[4] or
                                pytest.approx(d2ZdAB_m, FD_TOL) == h[4] or
                                between(h[4], d2ZdAB_p, d2ZdAB_m))
                        if not B_skip:
                            # Second derivative w.r.t. B is very sensitive
                            # near point that roots disappear, and is at a
                            # maximum (or minimum) so skip tests if close
                            # to this point
                            assert (
                                pytest.approx(d2ZdB2_p, FD_TOL) == h[5] or
                                pytest.approx(d2ZdB2_m, FD_TOL) == h[5] or
                                between(h[5], d2ZdB2_p, d2ZdB2_m))

                    except AssertionError:
                        # Print values at failure and raise exception
                        print(eos_type, T, P, A, B, Z)
                        raise


@pytest.mark.integration
@pytest.mark.skipif(not prop_available,
                    reason="Cubic root finder not available")
def test_roots_Z_vap_ext(root_finder):
    for eos_type in [0, 1]:
        u = EoS_param[eos_type]["u"]
        w = EoS_param[eos_type]["w"]

        for T in t_set:
            for P in p_set:
                # Calculate A and B parameters from Tr and Pr
                A = 0.5*P/T**2
                B = 0.1*P/T

                # Get results of external function call
                f = root_finder.proc_Z_vap_x
                assert(isinstance(f, AMPLExternalFunction))
                Z, g, h = f.evaluate_fgh(args=(eos_type, A, B))

                # Calculate parameters of cubic
                c1 = 1
                c2 = -(1+B-u*B)
                c3 = (A-u*B-(u-w)*B**2)
                c4 = -(A*B+w*B**2+w*B**3)

                det = c2**2 - 3*c3
                a = -(1.0/3.0)*(c2 - det**0.5)

                # Check to see if extension is triggered
                if det <= 0 or (a**3 + c2*a**2 + c3*a + c4) <= 0:
                    # Extension is not used
                    # Calculate residual and derivatives w.r.t. Z
                    res = c1*Z**3 + c2*Z**2 + c3*Z + c4
                    dz = 3*c1*Z**2 + 2*c2*Z + c3
                    dz2 = 6*c1*Z + 2*c2

                    try:
                        # Residual can be extemely sensitive to value of Z, so
                        # if residual fails tolerance, test size of Newton step
                        # to converge to root
                        assert (pytest.approx(0, abs=TOL) == res or
                                pytest.approx(0, abs=TOL) == value(res/dz))

                        # Check derivative signs to confirm correct root
                        assert dz >= 0  # Should always have non-negative slope

                        # Determine number of roots - calculate discriminant
                        dis = (18*c1*c2*c3*c4 - 4*c2**3*c4 + c2**2*c3**2 -
                               4*c1*c3**3 - 27*c1**2*c4**2)

                        if dis >= 0:
                            # Cubic has 2 or 3 real roots
                            # Second derivative should be non-negative
                            assert dz2 >= 0
                        # otherwise no need to check 2nd derivative

                    except AssertionError:
                        # Print values at failure and raise exception
                        print(eos_type, T, P, A, B, Z)
                        raise

                else:
                    # Extention is used, calculate extended root
                    c1x = 2*a
                    c2x = -c2 - 3.0*c1x
                    c3x = 3*c1x**2 + 2*c2*c1x + c3
                    c4x = c4 - 0.75*c1x**3 - 0.5*c2*c1x**2

                    # Calculate residual and derivatives w.r.t. Z_ext
                    res = c1*Z**3 + c2x*Z**2 + c3x*Z + c4x
                    dz = 3*c1*Z**2 + 2*c2x*Z + c3x

                    try:
                        # Residual can be extemely sensitive to value of Z, so
                        # if residual fails tolerance, test size of Newton step
                        # to converge to root
                        assert (pytest.approx(0, abs=TOL) == res or
                                pytest.approx(0, abs=TOL) == value(res/dz))

                        # Check derivative signs to confirm correct root
                        assert dz >= 0  # Should always have non-negative slope

                        # Determine number of roots - calculate discriminant
                        dis = (18*c1*c2x*c3x*c4x - 4*c2x**3*c4x +
                               c2x**2*c3x**2 - 4*c1*c3x**3 - 27*c1**2*c4x**2)

                        # Second derivative could be anything, don't check

                    except AssertionError:
                        # Print values at failure and raise exception
                        print(eos_type, T, P, A, B, Z)
                        raise

                if TEST_DERS:
                    try:
                        # Perform finite differences on A and B
                        ZAp, gAp, hAp = f.evaluate_fgh(
                                args=(eos_type, A*(1+DEL), B))
                        ZAm, gAm, hAm = f.evaluate_fgh(
                                args=(eos_type, A*(1-DEL), B))

                        ZBp, gBp, hBp = f.evaluate_fgh(
                                    args=(eos_type, A, B*(1+DEL)))
                        ZBm, gBm, hBm = f.evaluate_fgh(
                                    args=(eos_type, A, B*(1-DEL)))

                        # Check variance in Z values. A very large
                        # difference indicates a transition between
                        # single and multiple root regions, and hence that
                        # the partial derivatvies will be very sensitive.
                        # In these cases, skip the derivative tests.
                        if (abs(ZAp - Z) > 1e-3 or
                                abs(ZAm - Z) > 1e-3 or
                                abs(a) < 0.5):
                            A_skip = True
                        else:
                            A_skip = False
                        if (abs(ZBp - Z) > 1e-3 or
                                abs(ZBm - Z) > 1e-3 or
                                abs(dis) < 1e-7 or
                                abs(a) < 0.5):
                            B_skip = True
                        else:
                            B_skip = False

                        # Test gradient terms
                        # Calculate numerical first partial derivative
                        if not A_skip:
                            dZdA_p = (ZAp-Z)/(A*DEL)
                            dZdA_m = (Z-ZAm)/(A*DEL)

                        if not B_skip:
                            dZdB_p = (ZBp-Z)/(B*DEL)
                            dZdB_m = (Z-ZBm)/(B*DEL)

                        # Partial derivative w.r.t. EoS identifier
                        assert g[0] == 0

                        # Check that external function value lies within
                        # TOL of least one of the numerical values (delta+
                        # or delta-), OR lies between the two numerical
                        # values and within 10*TOL of one of the numerical
                        # values
                        if not A_skip:
                            assert (
                                pytest.approx(dZdA_p, FD_TOL) == g[1] or
                                pytest.approx(dZdA_m, FD_TOL) == g[1] or
                                (between(g[1], dZdA_p, dZdA_m) and (
                                 pytest.approx(dZdA_p, 10*FD_TOL) ==
                                 g[1] or
                                 pytest.approx(dZdA_m, 10*FD_TOL) ==
                                 g[1])))

                        if not B_skip:
                            assert (
                                pytest.approx(dZdB_p, FD_TOL) == g[2] or
                                pytest.approx(dZdB_m, FD_TOL) == g[2] or
                                (between(g[2], dZdB_p, dZdB_m) and (
                                 pytest.approx(dZdB_p, 10*FD_TOL) ==
                                 g[2] or
                                 pytest.approx(dZdB_m, 10*FD_TOL) ==
                                 g[2])))

                        # Test hessian terms
                        # Calculate numerical second partial derivatives
                        if not A_skip:
                            d2ZdA2_p = (gAp[1]-g[1])/(A*DEL)
                            d2ZdA2_m = (g[1]-gAm[1])/(A*DEL)

                        if not B_skip:
                            d2ZdB2_p = (gBp[2]-g[2])/(B*DEL)
                            d2ZdB2_m = (g[2]-gBm[2])/(B*DEL)

                        if not A_skip and not B_skip:
                            d2ZdAB_p = (gBp[1]-g[1])/(B*DEL)
                            d2ZdAB_m = (g[1]-gBm[1])/(B*DEL)

                        # Partial derivatives w.r.t eos_type
                        assert h[0] == 0
                        assert h[1] == 0
                        assert h[3] == 0

                        if not A_skip:
                            assert (
                                pytest.approx(d2ZdA2_p, FD_TOL) == h[2] or
                                pytest.approx(d2ZdA2_m, FD_TOL) == h[2] or
                                between(h[2], d2ZdA2_p, d2ZdA2_m))
                        if not A_skip and not B_skip:
                            assert (
                                pytest.approx(d2ZdAB_p, FD_TOL) == h[4] or
                                pytest.approx(d2ZdAB_m, FD_TOL) == h[4] or
                                between(h[4], d2ZdAB_p, d2ZdAB_m))
                        if not B_skip:
                            # Second derivative w.r.t. B is very sensitive
                            # near point that roots disappear, and is at a
                            # maximum (or minimum) so skip tests if close
                            # to this point
                            assert (
                                pytest.approx(d2ZdB2_p, FD_TOL) == h[5] or
                                pytest.approx(d2ZdB2_m, FD_TOL) == h[5] or
                                between(h[5], d2ZdB2_p, d2ZdB2_m))

                    except AssertionError:
                        # Print values at failure and raise exception
                        print(eos_type, T, P, A, B, Z)
                        raise
