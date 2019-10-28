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
Tests for the cubic root finder external functions
"""

import pytest
import os

from numpy import logspace

from pyomo.environ import (ConcreteModel,
                           ExternalFunction,
                           value)
from pyomo.core.base.external import AMPLExternalFunction

__author__ = "Andrew Lee"

# Set path to root finder .so file
_so = os.path.join(os.path.dirname(__file__), "../cubic_roots.so")

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
                        ZERO_CUT = 1e-10  # Aany derivative less than this is 0

                        # Partial derivatives w.r.t. EoS identifier
                        assert g[0] == 0
                        assert h[0] == 0

                        # Partial derivatives w.r.t. B
                        # Use small finite difference step and get new values
                        ZAp, gAp, hAp = f.evaluate_fgh(
                                args=(eos_type, A*(1+DEL), B))
                        ZAm, gAm, hAm = f.evaluate_fgh(
                                args=(eos_type, A*(1-DEL), B))

                        # Calculate numerical first partial derivative
                        dZdA_p = (ZAp-Z)/(A*DEL)
                        dZdA_m = (Z-ZAm)/(A*DEL)

                        # Calculate numerical second partial derivative
                        # For this, use finite differences on the external
                        # function gradients.
                        d2ZdA2_p = (gAp[1]-g[1])/(A*DEL)
                        d2ZdA2_m = (g[1]-gAm[1])/(A*DEL)

                        # Check that external function value agrees with one of
                        # the numerical values (delta+ or delta-), or lies
                        # between the two numerical values.
                        if abs(g[1]) > ZERO_CUT:
                            assert (dZdA_p == pytest.approx(g[1], FD_TOL) or
                                    dZdA_m == pytest.approx(g[1], FD_TOL) or
                                    between(g[1], dZdA_p, dZdA_m))

                        if abs(h[1]) > ZERO_CUT:
                            assert (d2ZdA2_p == pytest.approx(h[1], FD_TOL) or
                                    d2ZdA2_m == pytest.approx(h[1], FD_TOL) or
                                    between(h[1], d2ZdA2_p, d2ZdA2_m))

                        # Partial derivatives w.r.t. B
                        # Use small finite difference step and get new values
                        ZBp, gBp, hBp = f.evaluate_fgh(
                                args=(eos_type, A, B*(1+DEL)))
                        ZBm, gBm, hBm = f.evaluate_fgh(
                                args=(eos_type, A, B*(1-DEL)))

                        # Calculate numerical first partial derivative
                        dZdB_p = (ZBp-Z)/(B*DEL)
                        dZdB_m = (Z-ZBm)/(B*DEL)

                        # Calculate numerical second partial derivative
                        # For this, use finite differences on the external
                        # function gradients.
                        d2ZdB2_p = (gBp[2]-g[2])/(B*DEL)
                        d2ZdB2_m = (g[2]-gBm[2])/(B*DEL)

                        # Check that external function value agrees with one of
                        # the numerical values (delta+ or delta-), or lies
                        # between the two numerical values.
                        if abs(g[2]) > ZERO_CUT:
                            assert (dZdB_p == pytest.approx(g[2], FD_TOL) or
                                    dZdB_m == pytest.approx(g[2], FD_TOL) or
                                    between(g[2], dZdB_p, dZdB_m))

                        if abs(h[3]) > ZERO_CUT:
                            assert (d2ZdB2_p == pytest.approx(h[3], FD_TOL) or
                                    d2ZdB2_m == pytest.approx(h[3], FD_TOL) or
                                    between(h[3], d2ZdB2_p, d2ZdB2_m))

                except AssertionError:
                    # Print values at failure and raise exception
                    print(eos_type, T, P, A, B, Z)
                    raise


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
                        # Second derivative should be non-positive
                        assert dz2 >= 0
                    # otherwise no need to check 2nd derivative

                    if TEST_DERS:
                        ZERO_CUT = 1e-10  # Aany derivative less than this is 0

                        # Partial derivatives w.r.t. EoS identifier
                        assert g[0] == 0
                        assert h[0] == 0

                        # Partial derivatives w.r.t. B
                        # Use small finite difference step and get new values
                        ZAp, gAp, hAp = f.evaluate_fgh(
                                args=(eos_type, A*(1+DEL), B))
                        ZAm, gAm, hAm = f.evaluate_fgh(
                                args=(eos_type, A*(1-DEL), B))

                        # Calculate numerical first partial derivative
                        dZdA_p = (ZAp-Z)/(A*DEL)
                        dZdA_m = (Z-ZAm)/(A*DEL)

                        # Calculate numerical second partial derivative
                        # For this, use finite differences on the external
                        # function gradients.
                        d2ZdA2_p = (gAp[1]-g[1])/(A*DEL)
                        d2ZdA2_m = (g[1]-gAm[1])/(A*DEL)

                        # Check that external function value agrees with one of
                        # the numerical values (delta+ or delta-), or lies
                        # between the two numerical values.
                        if abs(g[1]) > ZERO_CUT:
                            assert (dZdA_p == pytest.approx(g[1], FD_TOL) or
                                    dZdA_m == pytest.approx(g[1], FD_TOL) or
                                    between(g[1], dZdA_p, dZdA_m))

                        if abs(h[1]) > ZERO_CUT:
                            assert (d2ZdA2_p == pytest.approx(h[1], FD_TOL) or
                                    d2ZdA2_m == pytest.approx(h[1], FD_TOL) or
                                    between(h[1], d2ZdA2_p, d2ZdA2_m))

                        # Partial derivatives w.r.t. B
                        # Use small finite difference step and get new values
                        ZBp, gBp, hBp = f.evaluate_fgh(
                                args=(eos_type, A, B*(1+DEL)))
                        ZBm, gBm, hBm = f.evaluate_fgh(
                                args=(eos_type, A, B*(1-DEL)))

                        # Calculate numerical first partial derivative
                        dZdB_p = (ZBp-Z)/(B*DEL)
                        dZdB_m = (Z-ZBm)/(B*DEL)

                        # Calculate numerical second partial derivative
                        # For this, use finite differences on the external
                        # function gradients.
                        d2ZdB2_p = (gBp[2]-g[2])/(B*DEL)
                        d2ZdB2_m = (g[2]-gBm[2])/(B*DEL)

                        # Check that external function value agrees with one of
                        # the numerical values (delta+ or delta-), or lies
                        # between the two numerical values.
                        if abs(g[2]) > ZERO_CUT:
                            assert (dZdB_p == pytest.approx(g[2], FD_TOL) or
                                    dZdB_m == pytest.approx(g[2], FD_TOL) or
                                    between(g[2], dZdB_p, dZdB_m))

                        if abs(h[3]) > ZERO_CUT:
                            assert (d2ZdB2_p == pytest.approx(h[3], FD_TOL) or
                                    d2ZdB2_m == pytest.approx(h[3], FD_TOL) or
                                    between(h[3], d2ZdB2_p, d2ZdB2_m))

                except AssertionError:
                    # Print values at failure and raise exception
                    print(eos_type, T, P, A, B, Z)
                    raise


# TODO : Add tests for extended functions. These are harder, and I need to look
# into what these do where a root is absent.
