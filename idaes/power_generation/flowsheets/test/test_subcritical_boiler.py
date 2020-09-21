##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
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
Make sure the subcritical boiler recirculation system example solves.
"""

__author__ = "Miguel Zamarripa"

import pytest
import pyomo.environ as pyo
from idaes.power_generation.flowsheets.subcritical_power_plant.\
    subcritical_boiler import main
from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              activated_equalities_generator)
from idaes.generic_models.properties import iapws95
from idaes.core.util.testing import get_default_solver
import idaes.core.util.scaling as iscale

solver_available = pyo.SolverFactory('ipopt').available()
prop_available = iapws95.iapws95_available()
solver = get_default_solver()


@pytest.fixture(scope="module")
def build_flowsheet():
    return main()


@pytest.mark.unit
def test_basic_build(build_flowsheet):
    """Make a turbine model and make sure it doesn't throw exception"""
    m = build_flowsheet
    # Check unit config arguments
    assert m.fs.Waterwalls[1].config.has_heat_transfer
    assert m.fs.Waterwalls[1].config.has_pressure_change
    assert m.fs.Waterwalls[1].config.property_package is m.fs.prop_water

    assert m.fs.downcomer.config.has_heat_transfer
    assert isinstance(m.fs.downcomer.heat_duty, pyo.Var)
    assert isinstance(m.fs.downcomer.deltaP, pyo.Var)
    assert isinstance(m.fs.drum.drum_level, pyo.Var)


@pytest.mark.skipif(not prop_available, reason="IAPWS not available")
@pytest.mark.skipif(not solver_available, reason="Solver not available")
@pytest.mark.component
def test_init(build_flowsheet):
    m = build_flowsheet
    # check that the model solved properly and has 0 degrees of freedom
    assert(degrees_of_freedom(m) == 0)

    m.fs.drum.feedwater_inlet.flow_mol[:].fix()
    m.fs.drum.feedwater_inlet.pressure[:].unfix()
    m.fs.drum.feedwater_inlet.enth_mol[:].fix()
    optarg = {
            "tol": 1e-6,
            "max_iter": 20}
    solver.options = optarg
    # set scaling parameters
    for i in m.fs.ww_zones:
        iscale.set_scaling_factor(m.fs.Waterwalls[i].heat_flux_conv[0], 1e-5)
    iscale.calculate_scaling_factors(m)
    results = solver.solve(m, tee=True)
    assert results.solver.termination_condition == \
        pyo.TerminationCondition.optimal
    assert results.solver.status == pyo.SolverStatus.ok

    assert (pyo.value(m.fs.downcomer.deltaP[0]) > 0)

    assert (pytest.approx(0, abs=1e-3) ==
            pyo.value(m.fs.Waterwalls[10].control_volume.
                      properties_out[0].flow_mol - m.fs.Waterwalls[1].
                      control_volume.properties_in[0].flow_mol))

    assert (pyo.value(m.fs.drum.feedwater_inlet.flow_mol[0]
                      - m.fs.drum.steam_outlet.flow_mol[0]) ==
            pytest.approx(0, abs=1e-3))
