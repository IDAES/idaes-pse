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
import pyomo.common.unittest as unittest
import pytest

from pyomo.contrib.incidence_analysis import IncidenceGraphInterface

from idaes.gas_distribution.flowsheets.simple_network_model import (
    make_simple_model,
)
from idaes.gas_distribution.flowsheets.run_simple_network_optimization import (
    run_dynamic_optimization,
)

class TestSimpleNetwork(unittest.TestCase):

    @pytest.mark.component
    def test_simple_network_model(self):
        model = make_simple_model()

        # Fix supply pressure, boost pressure, and demand flow rate
        model.fs.compressor.boost_pressure[:].fix()
        model.fs.nodes[0].state[:].pressure.fix()
        model.fs.nodes[1].demands[0].flow_mol[:].fix()

        # Model has zero degrees of freedom and a perfect matching
        # between constraints and variables.
        igraph = IncidenceGraphInterface(model)
        M = len(igraph.constraints)
        N = len(igraph.variables)
        matching = igraph.maximum_matching()
        self.assertEqual(N, M)
        self.assertEqual(len(matching), N)

    @pytest.mark.component
    def test_dynamic_optimization_result(self):
        """
        Test that dynamic optimization solve gives results we expect.
        These values are taken from the solve of this dynamic optimization
        problem itself, not a non-IDAES reproduction of the same solve,
        so this test only asserts that the values produced by the flowsheet
        script don't change.

        """
        sim_data = run_dynamic_optimization()

        time_list, value_map = sim_data

        pred_values = {}
        pred_values["fs.nodes[0].state[*].flow_mol"] = [
            1.6666e4, 1.7659e4, 1.9825e4, 2.1226e4, 2.1218e4, 2.0301e4,
            1.9687e4, 1.9957e4, 2.0113e4, 1.9948e4, 1.9991e4, 2.0030e4,
            1.9977e4, 2.0006e4, 2.0001e4, 1.9998e4, 1.9999e4, 2.0000e4,
            2.0000e4, 1.9999e4, 1.9999e4,
        ]
        pred_values["fs.nodes[1].state[*].flow_mol"] = [
            1.6666e4, 1.6666e4, 1.6666e4, 1.6666e4, 1.6666e4, 1.9999e4,
            1.9999e4, 1.9999e4, 1.9999e4, 1.9999e4, 1.9999e4, 1.9999e4,
            1.9999e4, 1.9999e4, 1.9999e4, 1.9999e4, 1.9999e4, 1.9999e4,
            1.9999e4, 1.9999e4, 1.9999e4,
        ]
        pred_values["fs.compressor.boost_pressure[*]"] = [
            7.000, 7.264, 7.994, 8.795, 9.320, 9.313, 9.174, 9.193, 9.237,
            9.209, 9.212, 9.222, 9.212, 9.216, 9.216, 9.216, 9.216, 9.216,
            9.216, 9.216, 9.216,
        ]
        pred_values["fs.nodes[1].state[*].pressure"] = [
            50.806, 50.812, 50.852, 50.979, 51.232, 50.599, 50.500, 50.500,
            50.504, 50.501, 50.500, 50.500, 50.500, 50.500, 50.500, 50.500,
            50.500, 50.500, 50.500, 50.500, 50.500,
        ]

        # Only test values for the four variables above
        value_map = {name: value_map[name] for name in pred_values}

        self.assertStructuredAlmostEqual(pred_values, value_map, reltol=0.01)


if __name__ == "__main__":
    unittest.main()
