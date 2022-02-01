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

from idaes.apps.nmpc.examples.pipeline.run_pipeline_nmpc import run_nmpc

class TestPipelineNMPC(unittest.TestCase):

    def test_small_nmpc_simulation(self):
        simulation_data, applied_inputs = run_nmpc(
            simulation_horizon=4.0,
            sample_period=2.0,
        )

        # These values were not reproduced with another model, just taken
        # from this NMPC simulation. They are here to make sure the result
        # doesn't change.
        pred_inlet_flow = [
            3.000e5, 4.927e5, 4.446e5, 4.460e5, 4.557e5, 4.653e5,
            4.722e5, 4.785e5, 4.842e5,
        ]
        pred_outlet_pressure = [
            50.91, 47.82, 46.19, 45.06, 43.58, 43.08, 42.62, 42.21, 41.83,
        ]

        pred_outlet_flow = [3.00000e5, 5.53003e5, 5.35472e5]
        pred_inlet_pressure = [57.000, 58.905, 58.920]

        actual_inlet_flow = simulation_data[1][
            "fs.pipeline.control_volume.flow_mass[*,0.0]"
        ]
        actual_outlet_pressure = simulation_data[1][
            "fs.pipeline.control_volume.pressure[*,1.0]"
        ]

        actual_outlet_flow = applied_inputs[1][
            "fs.pipeline.control_volume.flow_mass[*,1.0]"
        ]
        actual_inlet_pressure = applied_inputs[1][
            "fs.pipeline.control_volume.pressure[*,0.0]"
        ]

        self.assertStructuredAlmostEqual(
            pred_inlet_flow, actual_inlet_flow, reltol=1e-3
        )
        self.assertStructuredAlmostEqual(
            pred_outlet_pressure, actual_outlet_pressure, reltol=1e-3
        )
        self.assertStructuredAlmostEqual(
            pred_outlet_flow, actual_outlet_flow, reltol=1e-3
        )
        self.assertStructuredAlmostEqual(
            pred_inlet_pressure, actual_inlet_pressure, reltol=1e-3
        )


if __name__ == "__main__":
    unittest.main()
