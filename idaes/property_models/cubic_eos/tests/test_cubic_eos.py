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
                           Expression,
                           ExternalFunction,
                           Param,
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


# Set paramter for number of points to test in pressure and temperature space
SAMPLES = 250
t_set = logspace(-3, 3, num=SAMPLES, base=10, endpoint=True)
p_set = logspace(3, -3, num=SAMPLES, base=10, endpoint=True)

# Set tolerance for root finder checks
TOL = 1e-8

# Parameter indicating whether to test partial derivatives
TEST_DERS = True
# Relative finite difference step for partial derivatives
DEL = 1e-6
# Tolerance for accepting partial derivatives based on finite differences
FD_TOL = 1e-2


def between(y, x1, x2):
    return 0 > (y-x1)*(y-x2)


@pytest.mark.skipif(not prop_available,
                    reason="Cubic root finder not available")
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


def test_roots_Z_liq(root_finder):
    for eos_type in [0, 1]:
        u = EoS_param[eos_type]["u"]
        w = EoS_param[eos_type]["w"]

        for T in t_set:
            for P in p_set:
                A = P**2/T**2
                B = P/T

                f = root_finder.proc_Z_liq
                assert(isinstance(f, AMPLExternalFunction))
                Z, g, h = f.evaluate_fgh(args=(eos_type, A, B))

                # Calculate parameters of cubic
                c1 = 1
                c2 = (1+B-u*B)
                c3 = (A-u*B-(u-w)*B**2)
                c4 = (A*B+w*B**2+w*B**3)

                res = (Z**3 - (1+B-u*B)*Z**2 + (A-u*B-(u-w)*B**2)*Z -
                       (A*B+w*B**2+w*B**3))
                dz = (3*Z**2 - 2*(1+B-u*B)*Z + (A-u*B-(u-w)*B**2))
                dz2 = (6*Z - 2*(1+B-u*B))

                # Residual can be extemely sensitive to value of Z, so instead
                # test size of Newton step to converge to root
                assert pytest.approx(0, abs=TOL) == value(
                        res/dz)

                # Check derivative signs to confirm correct root
                assert dz >= 0  # Should always have non-negative slope

                # Determine number of roots - calculate discriminant
                dis = (18*c1*c2*c3*c4 - 4*c2**3*c4 + c2**2*c3**2 -
                       4*c1*c3**3 + 27*c1**2*c4**2)

                if dis >= 0:
                    # Cubic has 2 or 3 real roots
                    # Second derivative should be non-positive
                    assert dz2 <= 0
                # otherwise, only one root and no need to check 2nd derivative

                if TEST_DERS:
                    ZERO_CUT = 1e-10 # Consider any derivative less than this 0

                    # Partial derivatives w.r.t. EoS identifier - should be 0
                    assert g[0] == 0
                    assert h[0] == 0

                    # Partial derivatives w.r.t. B
                    ZAp, gAp, hAp = f.evaluate_fgh(
                            args=(eos_type, A*(1+DEL), B))
                    ZAm, gAm, hAm = f.evaluate_fgh(
                            args=(eos_type, A*(1-DEL), B))

                    dZdA_p = (ZAp-Z)/(A*DEL)
                    dZdA_m = (Z-ZAm)/(A*DEL)

                    d2ZdA2_p = (gAp[1]-g[1])/(A*DEL)
                    d2ZdA2_m = (g[1]-gAm[1])/(A*DEL)

                    if abs(g[1]) > ZERO_CUT:
                        assert (dZdA_p == pytest.approx(g[1], FD_TOL) or
                                dZdA_m == pytest.approx(g[1], FD_TOL) or
                                between(g[1], dZdA_p, dZdA_m))

                    if abs(h[1]) > ZERO_CUT:
                        assert (d2ZdA2_p == pytest.approx(h[1], FD_TOL) or
                                d2ZdA2_m == pytest.approx(h[1], FD_TOL) or
                                between(h[1], d2ZdA2_p, d2ZdA2_m))

                    # Partial derivatives w.r.t. B
                    ZBp, gBp, hBp = f.evaluate_fgh(
                            args=(eos_type, A, B*(1+DEL)))
                    ZBm, gBm, hBm = f.evaluate_fgh(
                            args=(eos_type, A, B*(1-DEL)))

                    dZdB_p = (ZBp-Z)/(B*DEL)
                    dZdB_m = (Z-ZBm)/(B*DEL)

                    d2ZdB2_p = (gBp[2]-g[2])/(B*DEL)
                    d2ZdB2_m = (g[2]-gBm[2])/(B*DEL)

                    if abs(g[2]) > ZERO_CUT:
                        assert (dZdB_p == pytest.approx(g[2], FD_TOL) or
                                dZdB_m == pytest.approx(g[2], FD_TOL) or
                                between(g[2], dZdB_p, dZdB_m))

                    if abs(h[2]) > ZERO_CUT:
                        assert (d2ZdB2_p == pytest.approx(h[2], FD_TOL) or
                                d2ZdB2_m == pytest.approx(h[2], FD_TOL) or
                                between(h[2], d2ZdB2_p, d2ZdB2_m))


# TODO : Add tests for extended functions. These are harder, and I need to look
# into hwat these do where a root is absent.
