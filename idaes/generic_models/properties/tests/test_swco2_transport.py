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

__author__ = "John Eslick"

import pytest
from pyomo.environ import ConcreteModel, value
import idaes.generic_models.properties.swco2 as swco2
import idaes


class TestSWCO2(object):

    @pytest.fixture(scope="class")
    def model2(self):
        model = ConcreteModel()
        model.prop_param = swco2.SWCO2ParameterBlock()
        model.prop_in = swco2.SWCO2StateBlock(
            default={"parameters": model.prop_param}
        )
        return model

    @pytest.mark.unit
    def test_transport(self, model2):
        """Transport property tests
        """

        data = (
           #(T,      P,        mu,          tc,         tol_mu,  tol_tc)
            (200,    0.1e6,    10.06e-6,    9.63e-3,    0.02,    0.02),
            (240,    0.1e6,    12.07e-6,    12.23e-3,   0.02,    0.02),
            (240,    10.0e6,   188.91e-6,   159.24e-3,  0.02,    0.02),
            (240,    25e6,     214.16e-6,   171.59e-3,  0.02,    0.02),
            (280,    0.1e6,    14.05e-6,    15.19e-3,   0.02,    0.02),
            (280,    2.5e6,    14.51e-6,    17.46e-3,   0.02,    0.04),
            (280,    30.0e6,   134.98e-6,   136.63e-3,  0.02,    0.02),
            (280,    50.0e6,   160.73e-6,   152.67e-3,  0.02,    0.02),
            (306,    7.5e6,    24.06e-6,    23.93e-3,   0.02,    0.40), # supercritcal near critical
            (400,    5.0e6,    20.37e-6,    27.46e-3,   0.02,    0.02),
            (400,    50.0e6,   65.67e-6,    86.48e-3,   0.02,    0.02),
            (400,    100e6,    101.41e-6,   121.41e-3,  0.02,    0.02),
            (600,    50e6,     42.20e-6,    64.11e-3,   0.02,    0.02),
            (1000,   5e6,      41.42e-6,    71.18e-3,   0.02,    0.02),
            (1000,   50e6,     46.16e-6,    80.33e-3,   0.02,    0.02),
            (1000,   100e6,    54.79e-6,    92.16e-3,   0.02,    0.02),
        )

        for d in data:
            tol_mu = d[4]
            tol_tc = d[5]
            P = d[1]
            T = d[0]
            mu_data = d[2]
            tc_data = d[3]
            model2.prop_in.temperature.set_value(T)
            model2.prop_in.pressure = P
            Tsat = value(model2.prop_in.temperature_sat)
            if P >= 7.377e6 and T >= 304.1282:
                # super critical, which we clasify as liquid
                ph = "Liq"
            if P <= 0.5179e6:
                # below triple point
                ph = "Vap"
            elif T > Tsat + 0.5:
                ph = "Vap"
            elif T < Tsat - 0.5:
                ph = "Liq"
            else:
                # if too close to sat, don't want to get into a situation where
                # I don't know if it's liquid or vapor, so using extreme
                # caution
                # If we're woried about it, we can add tests on the sat curve,
                # and test both liquid and vapor
                continue

            rho = value(model2.prop_in.dens_mass_phase[ph])
            mu = value(model2.prop_in.visc_d_phase[ph])
            tc = value(model2.prop_in.therm_cond_phase[ph])
            print("T = {}, P = {}, mu = {}, tc = {}, rho = {}, phase = {}".
                  format(T, P, mu, tc, rho, ph)
                  )
            print("Data: mu = {}, tc = {}".format(mu_data, tc_data))
            assert tc == pytest.approx(tc_data, rel=tol_tc)
            assert mu == pytest.approx(mu_data, rel=tol_mu)
