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
from pyomo.environ import ConcreteModel, value, units as pyunits
import idaes.generic_models.properties.swco2 as swco2
from idaes.generic_models.unit_models import Compressor
from idaes.core import FlowsheetBlock
import idaes
from idaes.core.util import get_solver

solver = get_solver()


class TestIntegration(object):
    @pytest.fixture(scope="class")
    def compressor_model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.properties = swco2.SWCO2ParameterBlock()
        m.fs.unit = Compressor(default={"property_package": m.fs.properties})
        return m

    @pytest.mark.unit
    def test_verify(self, compressor_model):
        model = compressor_model
        # Verify the turbine results against 3 known test cases

        # Case Data (90% isentropic efficency)
        cases = {
            "F": (1000, 1000),  # mol/s
            "Tin": (500, 300),  # K
            "Pin": (10, 100),  # kPa
            "W": (3414.29266, 2796.30966),  # kW
            "Tout": (574.744119, 372.6675676),  # K
            "Pout": (20, 250),  # kPa
            "xout": (1.0, 1.0),  # vapor fraction
            "Tisen": (567.418852, 365.7680891),
        }

        for i, F in enumerate(cases["F"]):
            Tin = cases["Tin"][i]
            Tout = cases["Tout"][i]
            Pin = cases["Pin"][i]*1000
            Pout = cases["Pout"][i]*1000
            hin = swco2.htpx(T=Tin*pyunits.K, P=Pin*pyunits.Pa)
            W = cases["W"][i]*1000
            Tis = cases["Tisen"][i]
            xout = cases["xout"][i]

            model.fs.unit.inlet.flow_mol[0].fix(F)
            model.fs.unit.inlet.enth_mol[0].fix(hin)
            model.fs.unit.inlet.pressure[0].fix(Pin)
            model.fs.unit.deltaP.fix(Pout - Pin)
            model.fs.unit.efficiency_isentropic.fix(0.9)
            model.fs.unit.initialize(optarg={'tol': 1e-6})
            solver.solve(model)

            Tout = pytest.approx(cases["Tout"][i], rel=1e-2)
            Pout = pytest.approx(cases["Pout"][i]*1000, rel=1e-2)
            Pout = pytest.approx(cases["Pout"][i]*1000, rel=1e-2)
            W = pytest.approx(cases["W"][i]*1000, rel=1e-2)
            xout = pytest.approx(xout, rel=1e-2)
            prop_out = model.fs.unit.control_volume.properties_out[0]
            prop_in = model.fs.unit.control_volume.properties_in[0]
            prop_is = model.fs.unit.properties_isentropic[0]

            assert value(prop_in.temperature) == pytest.approx(Tin, rel=1e-3)
            assert value(prop_is.temperature) == pytest.approx(Tis, rel=1e-3)
            assert value(model.fs.unit.control_volume.work[0]) == W
            assert value(prop_out.pressure) == Pout
            assert value(prop_out.temperature) == Tout
            assert value(prop_out.vapor_frac) == xout
