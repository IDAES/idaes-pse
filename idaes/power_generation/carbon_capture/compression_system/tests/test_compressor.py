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
Compressor model
The model consists equations and dimensionless numbers without
using geometry information.

"""
import pytest
# Import Pyomo libraries
import pyomo.environ as pyo
from pyomo.util.check_units import assert_units_consistent

# Import IDAES core
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom

# Import Unit Model Modules
# from idaes.generic_models.properties import iapws95
import idaes.generic_models.properties.swco2 as swco2

from idaes.power_generation.carbon_capture.compression_system.compressor import InletStage

from idaes.core.util.testing import get_default_solver, initialization_tester
# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_default_solver()

# -----------------------------------------------------------------------------


@pytest.fixture(scope="module")
def build_unit():
    # Create a Concrete Model as the top level object
    m = pyo.ConcreteModel()
    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties_co2 = swco2.SWCO2ParameterBlock() 
    # fs.properties_h2o = iapws95.Iapws95ParameterBlock()  
    m.fs.unit = InletStage(default={"property_package": m.fs.properties_co2})
        
    # Set the compressor inlet conditions and an initial flow guess
    # compressor stage 1
    p = 1.01353 * 1e5  # Pa
    t = 40.0113 + 273.15  # K
    fin = 1615.36  # mol/s
    hin_co2 = swco2.htpx(T=t*pyo.units.K, P=p*pyo.units.Pa)

    # inlet stream
    m.fs.unit.inlet.flow_mol[:].fix(fin)
    m.fs.unit.inlet.enth_mol[:].fix(hin_co2)
    m.fs.unit.inlet.pressure[:].fix(p)
    m.fs.unit.mass_flow_coeff.fix(0.09)

    # inlet specifications
    m.fs.unit.ratioP[0].value = 2.41692
    m.fs.unit.r2.value = 0.67654
    return m


@pytest.mark.unit
def test_basic_build(build_unit):
    """Make a model and make sure it doesn't throw exception"""
    m = build_unit
    assert degrees_of_freedom(m) == 0
    # Check unit config arguments
    assert len(m.fs.unit.config) == 10
    # assert m.fs.unit.config.has_heat_transfer
    assert m.fs.unit.config.thermodynamic_assumption
    assert m.fs.unit.config.property_package is m.fs.properties_co2


# @pytest.mark.integration
# def test_units(build_unit):
#     assert_units_consistent(build_unit)


# @pytest.mark.skipif(not iapws95.iapws95_available(),
#                     reason="IAPWS not available")
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_initialize_drum(build_unit):    
    initialization_tester(build_unit, dof=0)


# @pytest.mark.skipif(not iapws95.iapws95_available(),
#                     reason="IAPWS not available")
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_run_drum(build_unit):
    m = build_unit  

    optarg = {"tol": 1e-7,
              "linear_solver": "ma27",
              "max_iter": 40}
    solver.options = optarg
    # solve model
    results = solver.solve(m, tee=True)
    # Check for optimal solution
    assert results.solver.termination_condition == \
        pyo.TerminationCondition.optimal
    assert results.solver.status == pyo.SolverStatus.ok
    assert degrees_of_freedom(m) == 0
    
    # energy balance
    assert (pytest.approx(0, abs=1e-3) ==
            pyo.value(m.fs.unit.inlet.flow_mol[0]
                      * m.fs.unit.inlet.enth_mol[0]
                      - m.fs.unit.outlet.flow_mol[0]
                      * m.fs.unit.outlet.enth_mol[0]
                      + m.fs.unit.work_mechanical[0]))
    # pressure change
    assert (pytest.approx(439456013.08821166, abs=1e-3) ==
            pyo.value(m.fs.unit.deltaP[0]))
    # mass balance
    assert (pytest.approx(0, abs=1e-3) ==
            pyo.value(m.fs.unit.inlet.flow_mol[0]                      
                      - m.fs.unit.outlet.flow_mol[0]
                      ))
