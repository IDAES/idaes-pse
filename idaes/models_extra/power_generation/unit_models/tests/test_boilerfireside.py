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
"""
Unit operation model for a steam heater applicable to platen superheater
and roof superheater, model main equations:

* Heat is given by fire-side boiler model
* Calculate pressure change due to friction
* Calculate slag layer wall temperature
* Consider a layer of metal and a layer of slag

Created on Aug 27, 2020 by Boiler Team (J. Ma, M. Zamarripa)
"""
import pytest

# Import Pyomo libraries
import pyomo.environ as pyo

# Import IDAES core
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom

# Import Unit Model Modules
from idaes.models.properties import iapws95

from idaes.core.util.testing import initialization_tester
from idaes.core.solvers import get_solver
from idaes.models_extra.power_generation.properties import FlueGasParameterBlock
from idaes.models_extra.power_generation.unit_models.boiler_fireside import (
    BoilerFireside,
)
from idaes.models_extra.power_generation.unit_models.tests.datadictionary import (
    data_dic,
)

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()

# -----------------------------------------------------------------------------


@pytest.fixture(scope="module")
def build_unit():
    # Create a Concrete Model as the top level object
    m = pyo.ConcreteModel()
    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(dynamic=False)
    # Add property packages to flowsheet library
    m.fs.prop_fluegas = FlueGasParameterBlock()
    # boiler based on surrogate
    m.fs.unit = BoilerFireside(
        dynamic=False,
        property_package=m.fs.prop_fluegas,
        calculate_PA_SA_flows=False,
        number_of_zones=12,
        has_platen_superheater=True,
        has_roof_superheater=True,
        surrogate_dictionary=data_dic,
    )
    return m


@pytest.mark.unit
def test_basic_build(build_unit):
    """Make a boiler fire side model
    and make sure it doesn't throw exception"""
    m = build_unit
    assert degrees_of_freedom(m) == 29
    # Check unit config arguments
    assert len(m.fs.unit.config) == 9
    assert m.fs.unit.config.number_of_zones == 12
    assert isinstance(m.fs.unit.waterwall_heat, pyo.Var)
    assert isinstance(m.fs.unit.platen_heat, pyo.Var)
    assert isinstance(m.fs.unit.roof_heat, pyo.Var)
    assert m.fs.unit.config.property_package is m.fs.prop_fluegas


@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS not available")
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_initialize_unit(build_unit):
    m = build_unit
    # platen superheater
    # fix coal ultimate analysis on dry basis
    # using Coal Illinois No. 6 bituminous coal
    m.fs.unit.mf_C_coal_dry.fix(0.6300)
    m.fs.unit.mf_H_coal_dry.fix(0.0412)
    m.fs.unit.mf_O_coal_dry.fix(0.1300)
    m.fs.unit.mf_N_coal_dry.fix(0.0127)
    m.fs.unit.mf_S_coal_dry.fix(0.0284)
    m.fs.unit.mf_Ash_coal_dry.fix(0.1577)
    m.fs.unit.hhv_coal_dry.fix(2.581e007)
    m.fs.unit.frac_moisture_vaporized[:].fix(0.6)
    m.fs.unit.mf_H2O_coal_raw[:].fix(0.2)  # moisture content
    m.fs.unit.flowrate_coal_raw[:].fix(44.0)  # kg/s

    m.fs.unit.wall_temperature_waterwall[:, :].fix(690)
    m.fs.unit.wall_temperature_platen[:].fix(750)
    m.fs.unit.wall_temperature_roof[:].fix(650)
    m.fs.unit.fcorrection_heat_ww.fix(0.95)
    m.fs.unit.fcorrection_heat_platen.fix(0.95)
    # SCPC simulation approx flue gas = 21290.6999  # mol/s
    flow_mol_pa = 21290.6999 * 0.34  # approx. 1/3 as Primary air
    flow_mol_sa = 21290.6999 * 0.66  # approx. 2/3 as Secondary air

    state_args_PA = {
        "flow_mol_comp": {
            "H2O": 0.0078267 * flow_mol_pa,
            "CO2": 0.000337339 * flow_mol_pa,
            "N2": 0.783994 * flow_mol_pa,
            "O2": 0.20784 * flow_mol_pa,
            "SO2": 1e-5 * flow_mol_pa,
            "NO": 1e-5 * flow_mol_pa,
        },
        "temperature": 333.15,
        "pressure": 101325.00,
    }
    state_args_SA = {
        "flow_mol_comp": {
            "H2O": 0.0078267 * flow_mol_sa,
            "CO2": 0.000337339 * flow_mol_sa,
            "N2": 0.783994 * flow_mol_sa,
            "O2": 0.20784 * flow_mol_sa,
            "SO2": 1e-5 * flow_mol_sa,
            "NO": 1e-5 * flow_mol_sa,
        },
        "temperature": 650.15,
        "pressure": 101325.00,
    }

    for i in m.fs.prop_fluegas.component_list:
        m.fs.unit.primary_air_inlet.flow_mol_comp[:, i].fix(
            state_args_PA["flow_mol_comp"][i]
        )
        m.fs.unit.secondary_air_inlet.flow_mol_comp[:, i].fix(
            state_args_SA["flow_mol_comp"][i]
        )
    m.fs.unit.primary_air_inlet.pressure[:].fix(101325.00)
    m.fs.unit.secondary_air_inlet.pressure[:].fix(101325.00)
    m.fs.unit.primary_air_inlet.temperature[:].fix(333.15)
    m.fs.unit.secondary_air_inlet.temperature[:].fix(650.15)
    m.fs.unit.temperature_coal[:].fix(335.15)
    m.fs.unit.SR_lf.fix(1.0)
    m.fs.unit.deltaP.fix(1000)
    assert degrees_of_freedom(m) == 0

    initialization_tester(
        build_unit, dof=0, state_args_PA=state_args_PA, state_args_SA=state_args_SA
    )


@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS not available")
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_solve_unit(build_unit):
    m = build_unit
    assert degrees_of_freedom(m) == 0

    results = solver.solve(m, tee=True)
    # Check for optimal solution
    assert pyo.check_optimal_termination(results)

    # flue gas outlet temperature
    assert pytest.approx(1236.780228, abs=1e-3) == pyo.value(
        m.fs.unit.flue_gas_outlet.temperature[0]
    )

    # check material balance
    assert (
        pytest.approx(
            pyo.value(
                m.fs.unit.primary_air_moist[0].flow_mass
                + m.fs.unit.secondary_air[0].flow_mass
                + m.fs.unit.flowrate_coal_burner[0]
                - m.fs.unit.flue_gas[0].flow_mass
                - m.fs.unit.flow_mass_flyash[0]
            ),
            abs=1e-2,
        )
        == 0
    )

    #     # check energy balance
    assert (
        pytest.approx(
            pyo.value(
                (
                    m.fs.unit.primary_air_moist[0].flow_mol
                    * m.fs.unit.primary_air_moist[0].enth_mol
                    + m.fs.unit.secondary_air[0].flow_mol
                    * m.fs.unit.secondary_air[0].enth_mol
                    + m.fs.unit.h_coal[0] * m.fs.unit.flowrate_coal_burner[0]
                )
                - (
                    sum(m.fs.unit.waterwall_heat[0, j] for j in m.fs.unit.zones)
                    + (
                        m.fs.unit.platen_heat[0]
                        if m.fs.unit.config.has_platen_superheater is True
                        else 0
                    )
                    + (
                        m.fs.unit.roof_heat[0]
                        if m.fs.unit.config.has_roof_superheater is True
                        else 0
                    )
                    + m.fs.unit.flue_gas[0].flow_mol * m.fs.unit.flue_gas[0].enth_mol
                    + m.fs.unit.flowrate_ash[0]
                    * (
                        593 * (m.fs.unit.flue_gas[0].temperature - 298.15)
                        + 0.293
                        * (m.fs.unit.flue_gas[0].temperature ** 2.0 - 298.15**2.0)
                    )
                    + m.fs.unit.flowrate_daf_flyash[0] * m.fs.unit.hs_daf_flyash[0]
                )
            ),
            abs=1e-3,
        )
        == 0
    )


@pytest.fixture(scope="module")
def build_unit_option2():
    # Create a Concrete Model as the top level object
    m = pyo.ConcreteModel()
    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(dynamic=False)
    # Add property packages to flowsheet library
    m.fs.prop_fluegas = FlueGasParameterBlock()
    # boiler based on surrogate
    m.fs.unit = BoilerFireside(
        dynamic=False,
        property_package=m.fs.prop_fluegas,
        calculate_PA_SA_flows=True,
        number_of_zones=12,
        has_platen_superheater=True,
        has_roof_superheater=True,
        surrogate_dictionary=data_dic,
    )
    return m


@pytest.mark.unit
def test_build_option2(build_unit_option2):
    """Make a boiler fire side model with option 2
    and make sure it doesn't throw exception"""
    m = build_unit_option2
    assert degrees_of_freedom(m) == 19
    # Check unit config arguments
    assert len(m.fs.unit.config) == 9
    assert m.fs.unit.config.number_of_zones == 12
    assert isinstance(m.fs.unit.waterwall_heat, pyo.Var)
    assert isinstance(m.fs.unit.platen_heat, pyo.Var)
    assert isinstance(m.fs.unit.SR, pyo.Var)
    assert isinstance(m.fs.unit.ratio_PA2coal, pyo.Var)
    assert isinstance(m.fs.unit.roof_heat, pyo.Var)
    assert m.fs.unit.config.property_package is m.fs.prop_fluegas
