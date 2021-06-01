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
Tests for SOFC ROM.

Author: Alex Noring
"""

import pytest

from idaes.power_generation.properties.NGFC.ROM.SOFC_ROM import \
    build_SOFC_ROM, initialize_SOFC_ROM
from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              number_variables,
                                              number_total_constraints,
                                              number_unused_variables)
from idaes.core.util import get_solver

from pyomo.environ import (ConcreteModel,
                           Block,
                           SolverStatus,
                           TerminationCondition,
                           value)
from pyomo.util.check_units import assert_units_consistent


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


# -----------------------------------------------------------------------------
class TestSOFCROM(object):
    @pytest.fixture()
    def m(self):
        m = ConcreteModel()
        build_SOFC_ROM(m)
        return m

    @pytest.mark.build
    @pytest.mark.integration
    def test_build(self, m):
        assert degrees_of_freedom(m) == 9
        assert number_variables(m) == 13562
        assert number_total_constraints(m) == 13552
        assert number_unused_variables(m) == 0
        assert_units_consistent(m)

        assert hasattr(m.SOFC, 'ROM_input')
        assert hasattr(m.SOFC, 'norm_input')
        assert hasattr(m.SOFC, 'F')
        assert hasattr(m.SOFC, 'R')
        assert hasattr(m.SOFC, 'norm_output')
        assert hasattr(m.SOFC, 'ROM_output')

        assert hasattr(m.SOFC, 'current_density')
        assert hasattr(m.SOFC, 'fuel_temperature')
        assert hasattr(m.SOFC, 'internal_reforming')
        assert hasattr(m.SOFC, 'air_temperature')
        assert hasattr(m.SOFC, 'air_recirculation')
        assert hasattr(m.SOFC, 'OTC')
        assert hasattr(m.SOFC, 'fuel_util')
        assert hasattr(m.SOFC, 'air_util')
        assert hasattr(m.SOFC, 'pressure')

        assert hasattr(m.SOFC, 'anode_outlet_temperature')
        assert hasattr(m.SOFC, 'cathode_outlet_temperature')
        assert hasattr(m.SOFC, 'stack_voltage')
        assert hasattr(m.SOFC, 'max_cell_temperature')
        assert hasattr(m.SOFC, 'deltaT_cell')

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.integration
    def test_initialize(self, m):
        initialize_SOFC_ROM(m.SOFC)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.integration
    def test_solve(self, m):
        m.SOFC.current_density.fix(4000)
        m.SOFC.fuel_temperature.fix(348.3)
        m.SOFC.internal_reforming.fix(0.6)
        m.SOFC.air_temperature.fix(617.3)
        m.SOFC.air_recirculation.fix(0.5)
        m.SOFC.OTC.fix(2.1)
        m.SOFC.fuel_util.fix(0.8)
        m.SOFC.air_util.fix(0.449)
        m.SOFC.pressure.fix(1)

        assert degrees_of_freedom(m) == 0

        results = solver.solve(m)

        assert results.solver.termination_condition == \
            TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

        assert (pytest.approx(705.5, abs=1e-1) ==
                value(m.SOFC.anode_outlet_temperature))
        assert (pytest.approx(726.3, abs=1e-1) ==
                value(m.SOFC.cathode_outlet_temperature))
        assert (pytest.approx(0.867, abs=1e-3) ==
                value(m.SOFC.stack_voltage))
        assert (pytest.approx(749.1, abs=1e-1) ==
                value(m.SOFC.max_cell_temperature))
        assert (pytest.approx(99.9, abs=1e-1) ==
                value(m.SOFC.deltaT_cell))
