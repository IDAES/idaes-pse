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
"""Test the basic 0 = z^3 + bz^2 + cz + d calculator."""

import os
import pytest
import pyomo.environ as pyo
import idaes


def cubic_function(z, b, c, d):
    return z**3 + b * z**2 + c * z + d


def derivs(z, b, c, d):
    # Use the tested derivative formulas to make sure the external function
    # keeps returning the right answers. These formulas were confirmed and
    # documented in a jupyter notebook.
    g = [None] * 3
    h = [None] * 6
    g[0] = 1.0 / (-1 + c / z**2 + 2 * d / z**3)
    g[1] = 1.0 / (-2 * z - b + d / z**2)
    g[2] = 1.0 / (-3 * z**2 - 2 * b * z - c)
    """ the hessian structure is:
            b c d
          +------
        b | 0 1 3
        c | - 2 4
        d | - - 5
    """
    h[0] = g[0] ** 3 * (2 * c / z**3 + 6 * d / z**4)
    h[1] = g[0] ** 2 * (2 * c / z**3 * g[1] + 6 * d / z**4 * g[1] - 1 / z**2)
    h[2] = g[1] ** 3 * (2 + 2 * d / z**3)
    h[3] = g[0] ** 2 * (2 * c / z**3 * g[2] + 6 * d / z**4 * g[2] - 2 / z**3)
    h[4] = g[1] ** 2 * (2 * g[2] + 2 * d / z**3 * g[2] - 1 / z**2)
    h[5] = g[2] ** 3 * (6 * z + 2 * b)
    return g, h


@pytest.mark.unit
def test_general_cubic_root_finder():
    plib = os.path.join(idaes.bin_directory, "cubic_roots.so")
    m = pyo.ConcreteModel()
    m.croot_l = pyo.ExternalFunction(library=plib, function="cubic_root_l")
    m.croot_h = pyo.ExternalFunction(library=plib, function="cubic_root_h")
    param_dict = {
        1: {"b": -3, "c": 0.5, "d": -1, "three": False},
        2: {"b": -3, "c": 0.5, "d": 1, "three": True},
        3: {"b": 1, "c": 0.5, "d": -1, "three": False},
        4: {"b": 1, "c": 0.5, "d": 2, "three": False},
    }

    for k, v in param_dict.items():
        b = v["b"]
        c = v["c"]
        d = v["d"]
        fl, gl, hl = m.croot_l.evaluate_fgh(args=(b, c, d))
        fh, gh, hh = m.croot_h.evaluate_fgh(args=(b, c, d))
        if v["three"]:
            assert fl < fh
        else:
            assert pytest.approx(fl, abs=1e-5) == fh
        assert pytest.approx(0, abs=1e-6) == cubic_function(fl, b, c, d)
        assert pytest.approx(0, abs=1e-6) == cubic_function(fh, b, c, d)
        g, h = derivs(fl, b, c, d)
        for i, ge in enumerate(g):
            assert pytest.approx(ge, abs=1e-5) == gl[i]
        for i, he in enumerate(h):
            assert pytest.approx(he, abs=1e-5) == hl[i]
        g, h = derivs(fh, b, c, d)
        for i, ge in enumerate(g):
            assert pytest.approx(ge, abs=1e-5) == gh[i]
        for i, he in enumerate(h):
            assert pytest.approx(he, abs=1e-5) == hh[i]
