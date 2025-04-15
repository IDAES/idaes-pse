#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Tests for MEA solvent column model
Author: Anuja Deshpande, Andrew Lee
"""
# Import Python libraries
import copy
import pytest

# Import Pyomo libraries
from pyomo.environ import (
    ConcreteModel,
    value,
    assert_optimal_termination,
    units as pyunits,
)
import pyomo.common.unittest as unittest

# Import IDAES Libraries
import idaes
from idaes.core import FlowsheetBlock
from idaes.models_extra.column_models.MEAsolvent_column import MEAColumn
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models_extra.column_models.properties.MEA_vapor import (
    flue_gas as vaporconfig,
)
from idaes.models_extra.column_models.properties.MEA_solvent import (
    configuration as liquidconfig,
)

from idaes.core.util.model_statistics import degrees_of_freedom
from pyomo.util.check_units import assert_units_consistent
from idaes.core.util.testing import initialization_tester
from idaes.core.solvers import get_solver
from idaes.core.util.performance import PerformanceBaseClass
import idaes.core.util.scaling as iscale

import idaes.logger as idaeslog

solver = get_solver()


def _scale_mea_liquid_params(params, scaling_factor_flow_mol=3e-4):
    params.set_default_scaling("enth_mol_phase", 3e-4)
    params.set_default_scaling("pressure", 1e-5)
    params.set_default_scaling("temperature", 1)
    params.set_default_scaling("flow_mol", scaling_factor_flow_mol)
    params.set_default_scaling("flow_mol_phase", scaling_factor_flow_mol)

    params.set_default_scaling(
        "flow_mass_phase", scaling_factor_flow_mol / 24e-3
    )  # MW mixture ~= 24 g/Mol
    params.set_default_scaling("dens_mol_phase", 1 / 43000)
    params.set_default_scaling("visc_d_phase", 700)
    params.set_default_scaling("log_k_eq", 1)

    mole_frac_scaling_factors = {
        "H2O": 1,
        "MEA": 10,
        "CO2": 20,
    }
    mole_frac_true_scaling_factors = {
        "CO2": 1e4,  # Could go to 1e4 or 3e4
        "H2O": 1,
        "HCO3_-": 1000,
        "MEA": 30,
        "MEACOO_-": 30,
        "MEA_+": 30,
    }
    for comp, sf_x in mole_frac_scaling_factors.items():
        params.set_default_scaling("mole_frac_comp", sf_x, index=comp)
        params.set_default_scaling("mole_frac_phase_comp", sf_x, index=("Liq", comp))
        params.set_default_scaling(
            "flow_mol_phase_comp", sf_x * scaling_factor_flow_mol, index=("Liq", comp)
        )

    for comp, sf_x in mole_frac_true_scaling_factors.items():
        params.set_default_scaling(
            "mole_frac_phase_comp_true", sf_x, index=("Liq", comp)
        )
        params.set_default_scaling(
            "flow_mol_phase_comp_true",
            sf_x * scaling_factor_flow_mol,
            index=("Liq", comp),
        )

    params.set_default_scaling(
        "apparent_inherent_reaction_extent", 2e-2, index="bicarbonate"
    )
    params.set_default_scaling(
        "apparent_inherent_reaction_extent", 1e-3, index="carbamate"
    )


def _scale_mea_vapor_params(params, scaling_factor_flow_mol=3e-4):
    params.set_default_scaling("enth_mol_phase", 3e-4)
    params.set_default_scaling("pressure", 1e-5)
    params.set_default_scaling("temperature", 1)
    params.set_default_scaling("flow_mol", scaling_factor_flow_mol)
    params.set_default_scaling("flow_mol_phase", scaling_factor_flow_mol)

    params.set_default_scaling("flow_mol_phase_comp", 2 * scaling_factor_flow_mol)
    params.set_default_scaling(
        "flow_mass_phase", scaling_factor_flow_mol / 24e-3
    )  # Say MW ~=24 g/mol
    params.set_default_scaling("visc_d_phase", 6e4)


# -----------------------------------------------------------------------------
def build_model():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    # Set up molar flow bounds
    vaporconfig_copy = copy.deepcopy(vaporconfig)
    liquidconfig_copy = copy.deepcopy(liquidconfig)
    vaporconfig_copy["state_bounds"]["flow_mol"] = (
        0,
        1,
        1e6,
        pyunits.mol / pyunits.s,
    )
    liquidconfig_copy["state_bounds"]["flow_mol"] = (
        0,
        1,
        1e6,
        pyunits.mol / pyunits.s,
    )

    # Set up property package
    m.fs.vapor_properties = GenericParameterBlock(**vaporconfig_copy)
    m.fs.liquid_properties = GenericParameterBlock(**liquidconfig_copy)
    # Number of finite elements and finite element list in the spatial domain
    x_nfe = 40

    # Create an instance of the column in the flowsheet
    m.fs.unit = MEAColumn(
        finite_elements=x_nfe,
        vapor_phase={"property_package": m.fs.vapor_properties},
        liquid_phase={"property_package": m.fs.liquid_properties},
    )

    # Fix column design variables
    # Absorber diameter
    m.fs.unit.diameter_column.fix(18)

    # Absorber length
    m.fs.unit.length_column.fix(21.6)  # meter

    # Fix operating conditions

    # Flue gas inlet
    m.fs.unit.vapor_inlet.flow_mol.fix(19000)  # mol/sec
    m.fs.unit.vapor_inlet.temperature.fix(313.15)  # K
    m.fs.unit.vapor_inlet.pressure.fix(100000)  # Pa
    m.fs.unit.vapor_inlet.mole_frac_comp[0, "CO2"].fix(0.041)
    m.fs.unit.vapor_inlet.mole_frac_comp[0, "H2O"].fix(0.074)
    m.fs.unit.vapor_inlet.mole_frac_comp[0, "N2"].fix(0.763)
    m.fs.unit.vapor_inlet.mole_frac_comp[0, "O2"].fix(0.122)

    # Lean solvent inlet
    m.fs.unit.liquid_inlet.flow_mol.fix(24000)  # mol/sec
    m.fs.unit.liquid_inlet.temperature.fix(313.15)  # K
    m.fs.unit.liquid_inlet.pressure.fix(100000)  # Pa
    m.fs.unit.liquid_inlet.mole_frac_comp[0, "CO2"].fix(0.022)
    m.fs.unit.liquid_inlet.mole_frac_comp[0, "H2O"].fix(0.868)
    m.fs.unit.liquid_inlet.mole_frac_comp[0, "MEA"].fix(0.11)

    _scale_mea_liquid_params(m.fs.liquid_properties)
    _scale_mea_vapor_params(m.fs.vapor_properties)
    for t in m.fs.time:
        for x in m.fs.unit.liquid_phase.length_domain:
            iscale.set_scaling_factor(m.fs.unit.velocity_liq[t, x], 20)
            iscale.set_scaling_factor(
                m.fs.unit.interphase_mass_transfer[t, x, "CO2"], 1 / 20
            )
            iscale.set_scaling_factor(
                m.fs.unit.interphase_mass_transfer[t, x, "H2O"], 1 / 100
            )

    iscale.calculate_scaling_factors(m)

    return m


@pytest.mark.performance
class Test_MEAColumn_Performance(PerformanceBaseClass, unittest.TestCase):
    def build_model(self):
        return build_model()

    def initialize_model(self, model):
        model.fs.unit.initialize(
            optarg={
                "bound_push": 1e-6,
                "max_iter": 400,
            }
        )

    def solve_model(self, model):
        with idaes.temporary_config_ctx():
            # Get default solver for testing
            solver.options["bound_push"] = 1e-10

            res = solver.solve(model)

            assert_optimal_termination(res)


class TestAbsorber:
    @pytest.fixture(scope="class")
    def model(self):
        return build_model()

    @pytest.mark.unit
    def test_degrees_of_freedom(self, model):
        assert degrees_of_freedom(model) == 0

    @pytest.mark.component
    def test_unit_consistency(self, model):
        assert_units_consistent(model)

    @pytest.mark.component
    # @pytest.mark.xfail  # TODO: Remove once model is fixed
    def test_initialize(self, model):
        initialization_tester(
            model,
            outlvl=idaeslog.DEBUG,
            optarg={
                "bound_push": 1e-6,
                "max_iter": 400,
            },
        )

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    # @pytest.mark.xfail  # TODO: Remove once model is fixed
    def test_solve(self, model):
        with idaes.temporary_config_ctx():
            # Get default solver for testing
            solver.options["bound_push"] = 1e-10

            res = solver.solve(model)

        # Solver status/condition
        assert_optimal_termination(res)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    # @pytest.mark.xfail  # TODO: Remove once model is fixed
    def test_solution(self, model):
        assert pytest.approx(19436.66052, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.flow_mol[0]
        )
        assert pytest.approx(0.00126265, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.mole_frac_comp[0, "CO2"]
        )
        assert pytest.approx(0.13361959, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.mole_frac_comp[0, "H2O"]
        )
        assert pytest.approx(0.74585857, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.mole_frac_comp[0, "N2"]
        )
        assert pytest.approx(0.11925917, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.mole_frac_comp[0, "O2"]
        )
        assert pytest.approx(100000.0, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.pressure[0]
        )
        assert pytest.approx(327.04361, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.temperature[0]
        )

        assert pytest.approx(23563.339470, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.flow_mol[0]
        )
        assert pytest.approx(0.05442600, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.mole_frac_comp[0, "CO2"]
        )
        assert pytest.approx(0.83353555, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.mole_frac_comp[0, "H2O"]
        )
        assert pytest.approx(0.11203845, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.mole_frac_comp[0, "MEA"]
        )
        assert pytest.approx(100000.0, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.pressure[0]
        )
        assert pytest.approx(315.2712896, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.temperature[0]
        )

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, model):
        vap_in = model.fs.unit.vapor_phase.properties[0, 0]
        vap_out = model.fs.unit.vapor_phase.properties[0, 1]
        liq_in = model.fs.unit.liquid_phase.properties[0, 1]
        liq_out = model.fs.unit.liquid_phase.properties[0, 0]

        # Material conservation
        for j in ["CO2", "H2O"]:
            normalization = value(
                vap_in.get_material_flow_terms("Vap", j)
                + liq_in.get_material_flow_terms("Liq", j)
            )
            assert 1e-6 >= abs(
                value(
                    vap_in.get_material_flow_terms("Vap", j)
                    + liq_in.get_material_flow_terms("Liq", j)
                    - vap_out.get_material_flow_terms("Vap", j)
                    - liq_out.get_material_flow_terms("Liq", j)
                )
                / normalization
            )
        for j in ["N2", "O2"]:
            normalization = value(vap_in.get_material_flow_terms("Vap", j))
            assert 1e-6 >= abs(
                value(
                    vap_in.get_material_flow_terms("Vap", j)
                    - vap_out.get_material_flow_terms("Vap", j)
                )
                / normalization
            )
        for j in ["MEA"]:
            normalization = value(liq_in.get_material_flow_terms("Liq", j))
            assert 1e-6 >= abs(
                value(
                    liq_in.get_material_flow_terms("Liq", j)
                    - liq_out.get_material_flow_terms("Liq", j)
                )
                / normalization
            )

        # Energy conservation
        normalization = abs(
            value(
                vap_out.get_enthalpy_flow_terms("Vap")
                - vap_in.get_enthalpy_flow_terms("Vap")
            )
        )
        assert 1e-6 >= abs(
            value(
                vap_in.get_enthalpy_flow_terms("Vap")
                + liq_in.get_enthalpy_flow_terms("Liq")
                - vap_out.get_enthalpy_flow_terms("Vap")
                - liq_out.get_enthalpy_flow_terms("Liq")
            )
            / normalization
        )
