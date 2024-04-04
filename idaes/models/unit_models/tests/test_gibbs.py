#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Tests for Gibbs reactor.

Author: Andrew Lee
"""
import pytest

from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    Constraint,
    Suffix,
    TransformationFactory,
    value,
    units,
)

from idaes.core import FlowsheetBlock, EnergyBalanceType, MomentumBalanceType
from idaes.models.unit_models.gibbs_reactor import GibbsReactor
from idaes.models.properties.activity_coeff_models.methane_combustion_ideal import (
    MethaneParameterBlock as MethaneCombustionParameterBlock,
)
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    fixed_variables_set,
    activated_constraints_set,
    number_unused_variables,
)
from idaes.core.util.testing import PhysicalParameterTestBlock, initialization_tester
from idaes.core.solvers import get_solver
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.initialization import (
    BlockTriangularizationInitializer,
    SingleControlVolumeUnitInitializer,
    InitializationStatus,
)
from idaes.core.util import DiagnosticsToolbox


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


# -----------------------------------------------------------------------------
@pytest.mark.unit
def test_config():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = PhysicalParameterTestBlock()

    m.fs.unit = GibbsReactor(property_package=m.fs.properties)

    # Check unit config arguments
    assert len(m.fs.unit.config) == 9

    assert not m.fs.unit.config.dynamic
    assert not m.fs.unit.config.has_holdup
    assert m.fs.unit.config.energy_balance_type == EnergyBalanceType.useDefault
    assert m.fs.unit.config.momentum_balance_type == MomentumBalanceType.pressureTotal
    assert not m.fs.unit.config.has_heat_transfer
    assert not m.fs.unit.config.has_pressure_change
    assert m.fs.unit.config.property_package is m.fs.properties
    assert m.fs.unit.config.inert_species == []

    assert isinstance(m.fs.unit.gibbs_minimization, Constraint)
    assert len(m.fs.unit.gibbs_minimization) == 4

    assert not hasattr(m.fs.unit, "inert_species_balance")


@pytest.mark.unit
def test_inerts():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = PhysicalParameterTestBlock()

    m.fs.unit = GibbsReactor(property_package=m.fs.properties, inert_species=["c1"])

    assert isinstance(m.fs.unit.inert_species_balance, Constraint)
    assert len(m.fs.unit.inert_species_balance) == 2
    assert m.fs.unit.inert_species_balance[0, "p1", "c1"] != Constraint.Skip
    assert m.fs.unit.inert_species_balance[0, "p2", "c1"] != Constraint.Skip

    assert isinstance(m.fs.unit.gibbs_minimization, Constraint)
    assert len(m.fs.unit.gibbs_minimization) == 2


@pytest.mark.unit
def test_inerts_dependent_w_multi_phase():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = PhysicalParameterTestBlock()
    # Change elemental composition to introduce dependency
    m.fs.properties.element_comp = {
        "c1": {"H": 0, "He": 0, "Li": 3},
        "c2": {"H": 4, "He": 5, "Li": 0},
    }

    m.fs.unit = GibbsReactor(property_package=m.fs.properties, inert_species=["c1"])

    assert isinstance(m.fs.unit.inert_species_balance, Constraint)
    assert len(m.fs.unit.inert_species_balance) == 2
    assert m.fs.unit.inert_species_balance[0, "p1", "c1"] != Constraint.Skip
    assert m.fs.unit.inert_species_balance[0, "p2", "c1"] != Constraint.Skip

    assert isinstance(m.fs.unit.gibbs_minimization, Constraint)
    assert len(m.fs.unit.gibbs_minimization) == 2


@pytest.mark.unit
def test_inerts_dependent_w_single_phase():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = PhysicalParameterTestBlock()
    # Set phase list to only have 1 phase
    m.fs.properties.phase_list = ["p1"]
    # Change elemental composition to introduce dependency
    m.fs.properties.element_comp = {
        "c1": {"H": 0, "He": 0, "Li": 3},
        "c2": {"H": 4, "He": 5, "Li": 0},
    }

    m.fs.unit = GibbsReactor(property_package=m.fs.properties, inert_species=["c1"])

    assert isinstance(m.fs.unit.inert_species_balance, Constraint)
    assert len(m.fs.unit.inert_species_balance) == 0
    assert (0, "p1", "c1") not in m.fs.unit.inert_species_balance

    assert isinstance(m.fs.unit.gibbs_minimization, Constraint)
    assert len(m.fs.unit.gibbs_minimization) == 1


@pytest.mark.unit
def test_invalid_inert():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = PhysicalParameterTestBlock()

    with pytest.raises(
        ConfigurationError,
        match="fs.unit invalid component in inert_species "
        "argument. foo is not in the property package "
        "component list.",
    ):
        m.fs.unit = GibbsReactor(
            property_package=m.fs.properties, inert_species=["foo"]
        )


# -----------------------------------------------------------------------------
class TestMethane(object):
    @pytest.fixture(scope="class")
    def methane(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = MethaneCombustionParameterBlock()

        m.fs.unit = GibbsReactor(
            property_package=m.fs.properties,
            has_heat_transfer=True,
            has_pressure_change=True,
        )

        m.fs.unit.inlet.flow_mol[0].fix(230.0)
        m.fs.unit.inlet.mole_frac_comp[0, "H2"].fix(0.0435)
        m.fs.unit.inlet.mole_frac_comp[0, "N2"].fix(0.6522)
        m.fs.unit.inlet.mole_frac_comp[0, "O2"].fix(0.1739)
        m.fs.unit.inlet.mole_frac_comp[0, "CO2"].fix(1e-5)
        m.fs.unit.inlet.mole_frac_comp[0, "CH4"].fix(0.1304)
        m.fs.unit.inlet.mole_frac_comp[0, "CO"].fix(1e-5)
        m.fs.unit.inlet.mole_frac_comp[0, "H2O"].fix(1e-5)
        m.fs.unit.inlet.mole_frac_comp[0, "NH3"].fix(1e-5)
        m.fs.unit.inlet.temperature[0].fix(1500.0)
        m.fs.unit.inlet.pressure[0].fix(101325.0)

        m.fs.unit.outlet.temperature[0].fix(2844.38)
        m.fs.unit.deltaP.fix(0)

        # Fix some bounds to avoid potential log(0)
        # TODO: This really should be fixed in the property package, but breaks other tests
        m.fs.unit.control_volume.properties_out[0].pressure.setlb(1000)
        m.fs.unit.control_volume.properties_out[0].mole_frac_phase_comp.setlb(1e-12)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, methane):
        assert hasattr(methane.fs.unit, "inlet")
        assert len(methane.fs.unit.inlet.vars) == 4
        assert hasattr(methane.fs.unit.inlet, "flow_mol")
        assert hasattr(methane.fs.unit.inlet, "mole_frac_comp")
        assert hasattr(methane.fs.unit.inlet, "temperature")
        assert hasattr(methane.fs.unit.inlet, "pressure")

        assert hasattr(methane.fs.unit, "outlet")
        assert len(methane.fs.unit.outlet.vars) == 4
        assert hasattr(methane.fs.unit.outlet, "flow_mol")
        assert hasattr(methane.fs.unit.outlet, "mole_frac_comp")
        assert hasattr(methane.fs.unit.outlet, "temperature")
        assert hasattr(methane.fs.unit.outlet, "pressure")

        assert hasattr(methane.fs.unit, "gibbs_minimization")
        assert hasattr(methane.fs.unit, "heat_duty")
        assert hasattr(methane.fs.unit, "deltaP")

        assert number_variables(methane) == 80
        assert number_total_constraints(methane) == 67
        assert number_unused_variables(methane) == 0

    @pytest.mark.component
    def test_structural_issues(self, methane):
        dt = DiagnosticsToolbox(methane)
        dt.assert_no_structural_warnings()

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize_temperature(self, methane):
        initialization_tester(
            methane,
            optarg={"tol": 1e-6},
            state_args={
                "temperature": 2844.38,
                "pressure": 101325.0,
                "flow_mol": 251.05,
                "mole_frac_comp": {
                    "CH4": 1e-5,
                    "CO": 0.0916,
                    "CO2": 0.0281,
                    "H2": 0.1155,
                    "H2O": 0.1633,
                    "N2": 0.5975,
                    "NH3": 1e-5,
                    "O2": 0.0067,
                },
            },
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve_temperature(self, methane):
        methane.scaling_factor = Suffix(direction=Suffix.EXPORT)

        methane.scaling_factor[
            methane.fs.unit.control_volume.element_balances[0.0, "C"]
        ] = 0.0038968315684515787
        methane.scaling_factor[
            methane.fs.unit.control_volume.element_balances[0.0, "H"]
        ] = 0.0009690314543471861
        methane.scaling_factor[
            methane.fs.unit.control_volume.element_balances[0.0, "N"]
        ] = 0.0016665906198716563
        methane.scaling_factor[
            methane.fs.unit.control_volume.element_balances[0.0, "O"]
        ] = 0.0067608566657646
        methane.scaling_factor[
            methane.fs.unit.control_volume.enthalpy_balances[0.0]
        ] = 6.343688225967796e-08
        methane.scaling_factor[methane.fs.unit.control_volume.heat[0.0]] = (
            1.3415588575040103e-07
        )
        methane.scaling_factor[methane.fs.unit.control_volume.pressure_balance[0.0]] = (
            9.869232667160129e-06
        )
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].enth_mol_phase["Vap"]
        ] = 0.00010271414106049353
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].enth_mol_phase_comp[
                "Vap", "CH4"
            ]
        ] = 1.3404825737265415e-05
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].enth_mol_phase_comp[
                "Vap", "CO2"
            ]
        ] = 2.5411669038422445e-06
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].enth_mol_phase_comp[
                "Vap", "CO"
            ]
        ] = 9.047317470370035e-06
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].enth_mol_phase_comp[
                "Vap", "H2O"
            ]
        ] = 4.135136252739528e-06
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].enth_mol_phase_comp[
                "Vap", "H2"
            ]
        ] = 1
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].enth_mol_phase_comp[
                "Vap", "N2"
            ]
        ] = 1
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].enth_mol_phase_comp[
                "Vap", "NH3"
            ]
        ] = 2.1786492374727668e-05
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].enth_mol_phase_comp[
                "Vap", "O2"
            ]
        ] = 1
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].eq_comp["CH4"]
        ] = 0.03334222459322486
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].eq_comp["CO2"]
        ] = 434.782608695652
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].eq_comp["CO"]
        ] = 434.782608695652
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].eq_comp["H2O"]
        ] = 434.782608695652
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].eq_comp["H2"]
        ] = 0.09995002498750626
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].eq_comp["N2"]
        ] = 0.00666640001066624
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].eq_comp["NH3"]
        ] = 434.782608695652
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].eq_comp["O2"]
        ] = 0.025001875140635548
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].eq_enth_mol_phase["Vap"]
        ] = 5.9334197643529735e-05
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].eq_enth_mol_phase_comp[
                "Vap", "CH4"
            ]
        ] = 1.3404825737265415e-05
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].eq_enth_mol_phase_comp[
                "Vap", "CO2"
            ]
        ] = 2.5411669038422445e-06
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].eq_enth_mol_phase_comp[
                "Vap", "CO"
            ]
        ] = 9.047317470370035e-06
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].eq_enth_mol_phase_comp[
                "Vap", "H2O"
            ]
        ] = 4.135136252739528e-06
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].eq_enth_mol_phase_comp[
                "Vap", "H2"
            ]
        ] = 1.0
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].eq_enth_mol_phase_comp[
                "Vap", "N2"
            ]
        ] = 1.0
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].eq_enth_mol_phase_comp[
                "Vap", "NH3"
            ]
        ] = 2.1786492374727668e-05
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].eq_enth_mol_phase_comp[
                "Vap", "O2"
            ]
        ] = 1.0
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].eq_total
        ] = 0.004347826086956522
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].flow_mol_phase["Vap"]
        ] = 0.004347826086956522
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].mole_frac_phase_comp[
                "Vap", "CH4"
            ]
        ] = 7.668711656441719
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].mole_frac_phase_comp[
                "Vap", "CO2"
            ]
        ] = 99999.99999999997
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].mole_frac_phase_comp[
                "Vap", "CO"
            ]
        ] = 99999.99999999997
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].mole_frac_phase_comp[
                "Vap", "H2O"
            ]
        ] = 99999.99999999997
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].mole_frac_phase_comp[
                "Vap", "H2"
            ]
        ] = 22.98850574712644
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].mole_frac_phase_comp[
                "Vap", "N2"
            ]
        ] = 1.5332720024532351
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].mole_frac_phase_comp[
                "Vap", "NH3"
            ]
        ] = 99999.99999999997
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].mole_frac_phase_comp[
                "Vap", "O2"
            ]
        ] = 5.750431282346176
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].enth_mol_phase["Vap"]
        ] = 2.5797225634987406e-05
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].enth_mol_phase_comp[
                "Vap", "CH4"
            ]
        ] = 9.131356578608373e-05
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].enth_mol_phase_comp[
                "Vap", "CO2"
            ]
        ] = 3.8432408761344524e-06
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].enth_mol_phase_comp[
                "Vap", "CO"
            ]
        ] = 8.617482264231038e-06
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].enth_mol_phase_comp[
                "Vap", "H2O"
            ]
        ] = 5.614692038136825e-06
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].enth_mol_phase_comp[
                "Vap", "H2"
            ]
        ] = 1.1286872760356127e-05
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].enth_mol_phase_comp[
                "Vap", "N2"
            ]
        ] = 0.0007529656318755848
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].enth_mol_phase_comp[
                "Vap", "NH3"
            ]
        ] = 0.000147611121632822
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].enth_mol_phase_comp[
                "Vap", "O2"
            ]
        ] = 0.00011456899776931075
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].entr_mol_phase_comp[
                "Vap", "CH4"
            ]
        ] = 0.002273113949091129
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].entr_mol_phase_comp[
                "Vap", "CO2"
            ]
        ] = 0.0032746835561418483
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].entr_mol_phase_comp[
                "Vap", "CO"
            ]
        ] = 0.004593745058364088
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].entr_mol_phase_comp[
                "Vap", "H2O"
            ]
        ] = 0.004278478461249582
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].entr_mol_phase_comp[
                "Vap", "H2"
            ]
        ] = 0.005287458822573278
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].entr_mol_phase_comp[
                "Vap", "N2"
            ]
        ] = 0.005016168207707185
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].entr_mol_phase_comp[
                "Vap", "NH3"
            ]
        ] = 0.0029937456289480776
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].entr_mol_phase_comp[
                "Vap", "O2"
            ]
        ] = 0.003654125560818891
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_comp["CH4"]
        ] = 0.002827763345315764
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_comp["CO2"]
        ] = 0.08851222463945824
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_comp["CO"]
        ] = 0.020535852025143013
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_comp["H2O"]
        ] = 0.011302753862264373
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_comp["H2"]
        ] = 0.019411801141636542
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_comp["N2"]
        ] = 0.0033331760499149495
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_comp["NH3"]
        ] = 3476.4619148053275
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_comp["O2"]
        ] = 8.288029770353534
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_enth_mol_phase["Vap"]
        ] = 1.5616764299049252e-05
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_enth_mol_phase_comp[
                "Vap", "CH4"
            ]
        ] = 9.131356578608373e-05
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_enth_mol_phase_comp[
                "Vap", "CO2"
            ]
        ] = 3.8432408761344524e-06
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_enth_mol_phase_comp[
                "Vap", "CO"
            ]
        ] = 8.617482264231038e-06
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_enth_mol_phase_comp[
                "Vap", "H2O"
            ]
        ] = 5.614692038136825e-06
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_enth_mol_phase_comp[
                "Vap", "H2"
            ]
        ] = 1.1286872760356127e-05
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_enth_mol_phase_comp[
                "Vap", "N2"
            ]
        ] = 0.0007529656318755848
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_enth_mol_phase_comp[
                "Vap", "NH3"
            ]
        ] = 0.000147611121632822
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_enth_mol_phase_comp[
                "Vap", "O2"
            ]
        ] = 0.00011456899776931075
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_entr_mol_phase_comp[
                "Vap", "CH4"
            ]
        ] = 1.2281784483868407e-12
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_entr_mol_phase_comp[
                "Vap", "CO2"
            ]
        ] = 0.0032722586564032413
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_entr_mol_phase_comp[
                "Vap", "CO"
            ]
        ] = 0.004587058233184441
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_entr_mol_phase_comp[
                "Vap", "H2O"
            ]
        ] = 0.004273074490321421
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_entr_mol_phase_comp[
                "Vap", "H2"
            ]
        ] = 0.005277269330062221
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_entr_mol_phase_comp[
                "Vap", "N2"
            ]
        ] = 0.005007465522090417
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_entr_mol_phase_comp[
                "Vap", "NH3"
            ]
        ] = 0.0029918924828130663
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_entr_mol_phase_comp[
                "Vap", "O2"
            ]
        ] = 0.003650757205096906
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_gibbs_mol_phase_comp[
                "Vap", "CH4"
            ]
        ] = 5.675584743853875e-07
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_gibbs_mol_phase_comp[
                "Vap", "CO2"
            ]
        ] = 6.906678667100022e-07
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_gibbs_mol_phase_comp[
                "Vap", "CO"
            ]
        ] = 1.0328408888379269e-06
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_gibbs_mol_phase_comp[
                "Vap", "H2O"
            ]
        ] = 9.189402632308884e-07
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_gibbs_mol_phase_comp[
                "Vap", "H2"
            ]
        ] = 1.415412212742754e-06
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_gibbs_mol_phase_comp[
                "Vap", "N2"
            ]
        ] = 1.248468142622135e-06
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_gibbs_mol_phase_comp[
                "Vap", "NH3"
            ]
        ] = 7.41580819400193e-07
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_gibbs_mol_phase_comp[
                "Vap", "O2"
            ]
        ] = 9.033010042075156e-07
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_mol_frac_out
        ] = 0.8416262137210224
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_total
        ] = 0.002827763345315764
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].flow_mol
        ] = 0.003999061274127067
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].flow_mol_phase["Vap"]
        ] = 0.003999061274127067
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].gibbs_mol_phase_comp[
                "Vap", "CH4"
            ]
        ] = 8.062155808993221e-07
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].gibbs_mol_phase_comp[
                "Vap", "CO2"
            ]
        ] = 8.859012593383004e-07
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].gibbs_mol_phase_comp[
                "Vap", "CO"
            ]
        ] = 1.3601211872489997e-06
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].gibbs_mol_phase_comp[
                "Vap", "H2O"
            ]
        ] = 1.1863588686277076e-06
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].gibbs_mol_phase_comp[
                "Vap", "H2"
            ]
        ] = 2.225437251662786e-06
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].gibbs_mol_phase_comp[
                "Vap", "N2"
            ]
        ] = 1.767676703310131e-06
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].gibbs_mol_phase_comp[
                "Vap", "NH3"
            ]
        ] = 1.0450609388337302e-06
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].gibbs_mol_phase_comp[
                "Vap", "O2"
            ]
        ] = 1.2704369868483197e-06
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].mole_frac_comp["CH4"]
        ] = 1
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].mole_frac_comp["CO2"]
        ] = 44.26650084714147
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].mole_frac_comp["CO"]
        ] = 10.270336270166638
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].mole_frac_comp["H2O"]
        ] = 5.6527035158927825
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].mole_frac_comp["H2"]
        ] = 9.70817890049651
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].mole_frac_comp["N2"]
        ] = 1.6669792340916454
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].mole_frac_comp["NH3"]
        ] = 1738638.9837520984
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].mole_frac_comp["O2"]
        ] = 4144.987636961217
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[
                "Vap", "CH4"
            ]
        ] = 1
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[
                "Vap", "CO2"
            ]
        ] = 44.26650084714147
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[
                "Vap", "CO"
            ]
        ] = 10.270336270166638
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[
                "Vap", "H2O"
            ]
        ] = 5.6527035158927825
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[
                "Vap", "H2"
            ]
        ] = 9.70817890049651
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[
                "Vap", "N2"
            ]
        ] = 1.6669792340916454
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[
                "Vap", "NH3"
            ]
        ] = 1738638.9837520984
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[
                "Vap", "O2"
            ]
        ] = 4144.987636961217
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].pressure
        ] = 9.869232667160129e-06
        methane.scaling_factor[
            methane.fs.unit.gibbs_minimization[0.0, "Vap", "CH4"]
        ] = 6.37201815837067e-07
        methane.scaling_factor[
            methane.fs.unit.gibbs_minimization[0.0, "Vap", "CO2"]
        ] = 7.052606407804767e-07
        methane.scaling_factor[methane.fs.unit.gibbs_minimization[0.0, "Vap", "CO"]] = (
            1.1096130740800528e-06
        )
        methane.scaling_factor[
            methane.fs.unit.gibbs_minimization[0.0, "Vap", "H2O"]
        ] = 9.679516364316312e-07
        methane.scaling_factor[methane.fs.unit.gibbs_minimization[0.0, "Vap", "H2"]] = (
            1.5736217717559093e-06
        )
        methane.scaling_factor[methane.fs.unit.gibbs_minimization[0.0, "Vap", "N2"]] = (
            1.2499361838560746e-06
        )
        methane.scaling_factor[
            methane.fs.unit.gibbs_minimization[0.0, "Vap", "NH3"]
        ] = 8.304717457164266e-07
        methane.scaling_factor[methane.fs.unit.gibbs_minimization[0.0, "Vap", "O2"]] = (
            8.983346084706515e-07
        )
        methane.scaling_factor[methane.fs.unit.lagrange_mult[0.0, "C"]] = (
            2.926858666284934e-06
        )
        methane.scaling_factor[methane.fs.unit.lagrange_mult[0.0, "H"]] = (
            4.450874503325572e-06
        )
        methane.scaling_factor[methane.fs.unit.lagrange_mult[0.0, "N"]] = (
            3.535353406620262e-06
        )
        methane.scaling_factor[methane.fs.unit.lagrange_mult[0.0, "O"]] = (
            2.5408739736966393e-06
        )

        scaling = TransformationFactory("core.scale_model")
        sm = scaling.create_using(methane, rename=False)

        results = solver.solve(sm)

        scaling.propagate_solution(sm, methane)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution_temperature(self, methane):
        assert pytest.approx(250.06, abs=1e-2) == value(
            methane.fs.unit.outlet.flow_mol[0]
        )
        assert pytest.approx(0.0, abs=1e-4) == value(
            methane.fs.unit.outlet.mole_frac_comp[0, "CH4"]
        )
        assert pytest.approx(0.0974, abs=1e-4) == value(
            methane.fs.unit.outlet.mole_frac_comp[0, "CO"]
        )
        assert pytest.approx(0.0226, abs=1e-4) == value(
            methane.fs.unit.outlet.mole_frac_comp[0, "CO2"]
        )
        assert pytest.approx(0.1030, abs=1e-4) == value(
            methane.fs.unit.outlet.mole_frac_comp[0, "H2"]
        )
        assert pytest.approx(0.1769, abs=1e-4) == value(
            methane.fs.unit.outlet.mole_frac_comp[0, "H2O"]
        )
        assert pytest.approx(0.5999, abs=1e-4) == value(
            methane.fs.unit.outlet.mole_frac_comp[0, "N2"]
        )
        assert pytest.approx(0.0, abs=1e-5) == value(
            methane.fs.unit.outlet.mole_frac_comp[0, "NH3"]
        )
        assert pytest.approx(0.0002, abs=1e-4) == value(
            methane.fs.unit.outlet.mole_frac_comp[0, "O2"]
        )
        assert pytest.approx(-7454077, abs=1e2) == value(methane.fs.unit.heat_duty[0])
        assert pytest.approx(101325.0, abs=1e-2) == value(
            methane.fs.unit.outlet.pressure[0]
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation_temperature(self, methane):

        assert (
            abs(
                value(
                    (
                        methane.fs.unit.inlet.flow_mol[0]
                        * methane.fs.unit.control_volume.properties_in[
                            0
                        ].enth_mol_phase["Vap"]
                        - methane.fs.unit.outlet.flow_mol[0]
                        * methane.fs.unit.control_volume.properties_out[
                            0
                        ].enth_mol_phase["Vap"]
                        + methane.fs.unit.heat_duty[0]
                    )
                    / methane.fs.unit.heat_duty[0]
                )
            )
            <= 1e-6
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_numerical_issues(self, methane):
        scaling = TransformationFactory("core.scale_model")
        sm = scaling.create_using(methane, rename=False)

        dt = DiagnosticsToolbox(sm)
        dt.assert_no_numerical_warnings()

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize_duty(self, methane):
        methane.fs.unit.outlet.temperature[0].unfix()
        methane.fs.unit.heat_duty.fix(-7454077)
        assert degrees_of_freedom(methane) == 0

        orig_fixed_vars = fixed_variables_set(methane)
        orig_act_consts = activated_constraints_set(methane)

        methane.fs.unit.initialize(
            optarg={"tol": 1e-6},
            state_args={
                "temperature": 2844.38,
                "pressure": 101325.0,
                "flow_mol": 251.05,
                "mole_frac_comp": {
                    "CH4": 1e-5,
                    "CO": 0.0916,
                    "CO2": 0.0281,
                    "H2": 0.1155,
                    "H2O": 0.1633,
                    "N2": 0.5975,
                    "NH3": 1e-5,
                    "O2": 0.0067,
                },
            },
        )

        assert degrees_of_freedom(methane) == 0

        fin_fixed_vars = fixed_variables_set(methane)
        fin_act_consts = activated_constraints_set(methane)

        assert len(fin_act_consts) == len(orig_act_consts)
        assert len(fin_fixed_vars) == len(orig_fixed_vars)

        for c in fin_act_consts:
            assert c in orig_act_consts
        for v in fin_fixed_vars:
            assert v in orig_fixed_vars

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve_heat_duty(self, methane):
        # Remove previous scaling factors
        methane.del_component(methane.scaling_factor)

        methane.scaling_factor = Suffix(direction=Suffix.EXPORT)

        methane.scaling_factor[
            methane.fs.unit.control_volume.element_balances[0.0, "C"]
        ] = 0.003895771416286546
        methane.scaling_factor[
            methane.fs.unit.control_volume.element_balances[0.0, "H"]
        ] = 0.0009691213173623994
        methane.scaling_factor[
            methane.fs.unit.control_volume.element_balances[0.0, "N"]
        ] = 0.0016665908743918427
        methane.scaling_factor[
            methane.fs.unit.control_volume.element_balances[0.0, "O"]
        ] = 0.006733906435857895
        methane.scaling_factor[
            methane.fs.unit.control_volume.enthalpy_balances[0.0]
        ] = 7.199382898409605e-08
        methane.scaling_factor[methane.fs.unit.control_volume.pressure_balance[0.0]] = (
            9.869232667160129e-06
        )
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].enth_mol_phase["Vap"]
        ] = 0.00010271414106049353
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].enth_mol_phase_comp[
                "Vap", "CH4"
            ]
        ] = 1.3404825737265415e-05
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].enth_mol_phase_comp[
                "Vap", "CO2"
            ]
        ] = 2.5411669038422445e-06
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].enth_mol_phase_comp[
                "Vap", "CO"
            ]
        ] = 9.047317470370035e-06
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].enth_mol_phase_comp[
                "Vap", "H2O"
            ]
        ] = 4.135136252739528e-06
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].enth_mol_phase_comp[
                "Vap", "H2"
            ]
        ] = 1
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].enth_mol_phase_comp[
                "Vap", "N2"
            ]
        ] = 1
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].enth_mol_phase_comp[
                "Vap", "NH3"
            ]
        ] = 2.1786492374727668e-05
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].enth_mol_phase_comp[
                "Vap", "O2"
            ]
        ] = 1
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].eq_comp["CH4"]
        ] = 0.03334222459322486
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].eq_comp["CO2"]
        ] = 434.782608695652
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].eq_comp["CO"]
        ] = 434.782608695652
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].eq_comp["H2O"]
        ] = 434.782608695652
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].eq_comp["H2"]
        ] = 0.09995002498750626
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].eq_comp["N2"]
        ] = 0.00666640001066624
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].eq_comp["NH3"]
        ] = 434.782608695652
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].eq_comp["O2"]
        ] = 0.025001875140635548
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].eq_enth_mol_phase["Vap"]
        ] = 5.9334197643529735e-05
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].eq_enth_mol_phase_comp[
                "Vap", "CH4"
            ]
        ] = 1.3404825737265415e-05
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].eq_enth_mol_phase_comp[
                "Vap", "CO2"
            ]
        ] = 2.5411669038422445e-06
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].eq_enth_mol_phase_comp[
                "Vap", "CO"
            ]
        ] = 9.047317470370035e-06
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].eq_enth_mol_phase_comp[
                "Vap", "H2O"
            ]
        ] = 4.135136252739528e-06
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].eq_enth_mol_phase_comp[
                "Vap", "H2"
            ]
        ] = 1.0
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].eq_enth_mol_phase_comp[
                "Vap", "N2"
            ]
        ] = 1.0
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].eq_enth_mol_phase_comp[
                "Vap", "NH3"
            ]
        ] = 2.1786492374727668e-05
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].eq_enth_mol_phase_comp[
                "Vap", "O2"
            ]
        ] = 1.0
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].eq_total
        ] = 0.004347826086956522
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].flow_mol_phase["Vap"]
        ] = 0.004347826086956522
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].mole_frac_phase_comp[
                "Vap", "CH4"
            ]
        ] = 7.668711656441719
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].mole_frac_phase_comp[
                "Vap", "CO2"
            ]
        ] = 99999.99999999997
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].mole_frac_phase_comp[
                "Vap", "CO"
            ]
        ] = 99999.99999999997
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].mole_frac_phase_comp[
                "Vap", "H2O"
            ]
        ] = 99999.99999999997
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].mole_frac_phase_comp[
                "Vap", "H2"
            ]
        ] = 22.98850574712644
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].mole_frac_phase_comp[
                "Vap", "N2"
            ]
        ] = 1.5332720024532351
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].mole_frac_phase_comp[
                "Vap", "NH3"
            ]
        ] = 99999.99999999997
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_in[0.0].mole_frac_phase_comp[
                "Vap", "O2"
            ]
        ] = 5.750431282346176
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].enth_mol_phase["Vap"]
        ] = 2.579089235408562e-05
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].enth_mol_phase_comp[
                "Vap", "CH4"
            ]
        ] = 0.00013694447144485197
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].enth_mol_phase_comp[
                "Vap", "CO2"
            ]
        ] = 3.048997164495024e-06
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].enth_mol_phase_comp[
                "Vap", "CO"
            ]
        ] = 1.032900097747103e-05
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].enth_mol_phase_comp[
                "Vap", "H2O"
            ]
        ] = 4.977043447527836e-06
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].enth_mol_phase_comp[
                "Vap", "H2"
            ]
        ] = 2.3899280715507977e-05
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].enth_mol_phase_comp[
                "Vap", "N2"
            ]
        ] = 6.31659002014246e-05
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].enth_mol_phase_comp[
                "Vap", "NH3"
            ]
        ] = 0.0002135671616883595
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].enth_mol_phase_comp[
                "Vap", "O2"
            ]
        ] = 8.194560075137228e-05
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].entr_mol_phase_comp[
                "Vap", "CH4"
            ]
        ] = 0.002416655313418302
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].entr_mol_phase_comp[
                "Vap", "CO2"
            ]
        ] = 0.0035442516743102712
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].entr_mol_phase_comp[
                "Vap", "CO"
            ]
        ] = 0.004456473678792036
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].entr_mol_phase_comp[
                "Vap", "H2O"
            ]
        ] = 0.004453410785746496
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].entr_mol_phase_comp[
                "Vap", "H2"
            ]
        ] = 0.005822282444612275
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].entr_mol_phase_comp[
                "Vap", "N2"
            ]
        ] = 0.004883049120795377
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].entr_mol_phase_comp[
                "Vap", "NH3"
            ]
        ] = 0.0030506224994020967
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].entr_mol_phase_comp[
                "Vap", "O2"
            ]
        ] = 0.0031405744684195823
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_comp["CH4"]
        ] = 0.0028284396804244337
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_comp["CO2"]
        ] = 0.12404165785988744
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_comp["CO"]
        ] = 0.019256175462170257
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_comp["H2O"]
        ] = 0.010875650523130753
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_comp["H2"]
        ] = 0.020815988487732813
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_comp["N2"]
        ] = 0.003333178086069927
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_comp["NH3"]
        ] = 1528.6054813156613
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_comp["O2"]
        ] = 696.9862666112366
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_enth_mol_phase["Vap"]
        ] = 1.4497108873652463e-05
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_enth_mol_phase_comp[
                "Vap", "CH4"
            ]
        ] = 6.9059089891022035e-06
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_enth_mol_phase_comp[
                "Vap", "CO2"
            ]
        ] = 2.444674689765088e-06
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_enth_mol_phase_comp[
                "Vap", "CO"
            ]
        ] = 1.0075547572052885e-05
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_enth_mol_phase_comp[
                "Vap", "H2O"
            ]
        ] = 4.335306573370538e-06
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_enth_mol_phase_comp[
                "Vap", "H2"
            ]
        ] = 5.8333326983598445e-06
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_enth_mol_phase_comp[
                "Vap", "N2"
            ]
        ] = 5.6424839623995114e-05
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_enth_mol_phase_comp[
                "Vap", "NH3"
            ]
        ] = 1.7166243852082998e-05
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_enth_mol_phase_comp[
                "Vap", "O2"
            ]
        ] = 3.12933823458493e-05
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_entr_mol_phase_comp[
                "Vap", "CH4"
            ]
        ] = 1.1897376657481283e-11
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_entr_mol_phase_comp[
                "Vap", "CO2"
            ]
        ] = 0.003320984391726235
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_entr_mol_phase_comp[
                "Vap", "CO"
            ]
        ] = 0.004446525169209352
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_entr_mol_phase_comp[
                "Vap", "H2O"
            ]
        ] = 0.004347558430163033
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_entr_mol_phase_comp[
                "Vap", "H2"
            ]
        ] = 0.00536900022840844
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_entr_mol_phase_comp[
                "Vap", "N2"
            ]
        ] = 0.004874347940539925
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_entr_mol_phase_comp[
                "Vap", "NH3"
            ]
        ] = 0.003039965534465019
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_entr_mol_phase_comp[
                "Vap", "O2"
            ]
        ] = 0.0031359732136764075
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_gibbs_mol_phase_comp[
                "Vap", "CH4"
            ]
        ] = 5.950956432532686e-07
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_gibbs_mol_phase_comp[
                "Vap", "CO2"
            ]
        ] = 7.15349973612882e-07
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_gibbs_mol_phase_comp[
                "Vap", "CO"
            ]
        ] = 1.0278492215315099e-06
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_gibbs_mol_phase_comp[
                "Vap", "H2O"
            ]
        ] = 9.452962659471516e-07
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_gibbs_mol_phase_comp[
                "Vap", "H2"
            ]
        ] = 1.484262777482121e-06
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_gibbs_mol_phase_comp[
                "Vap", "N2"
            ]
        ] = 1.2185350436147026e-06
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_gibbs_mol_phase_comp[
                "Vap", "NH3"
            ]
        ] = 7.51571959805243e-07
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_gibbs_mol_phase_comp[
                "Vap", "O2"
            ]
        ] = 7.795160552447015e-07
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_mol_frac_out
        ] = 0.8409200257976712
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].eq_total
        ] = 0.0028284396804244332
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].flow_mol
        ] = 0.004000017756410457
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].flow_mol_phase["Vap"]
        ] = 0.004000017756410457
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].gibbs_mol_phase_comp[
                "Vap", "CH4"
            ]
        ] = 1.0256057645532026e-06
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].gibbs_mol_phase_comp[
                "Vap", "CO2"
            ]
        ] = 1.0123237586232911e-06
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].gibbs_mol_phase_comp[
                "Vap", "CO"
            ]
        ] = 1.6087623724494755e-06
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].gibbs_mol_phase_comp[
                "Vap", "H2O"
            ]
        ] = 1.377288144399382e-06
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].gibbs_mol_phase_comp[
                "Vap", "H2"
            ]
        ] = 2.77905882129248e-06
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].gibbs_mol_phase_comp[
                "Vap", "N2"
            ]
        ] = 2.159334308222618e-06
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].gibbs_mol_phase_comp[
                "Vap", "NH3"
            ]
        ] = 1.2965057859756173e-06
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].gibbs_mol_phase_comp[
                "Vap", "O2"
            ]
        ] = 1.3652606771737144e-06
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].mole_frac_comp["CH4"]
        ] = 1
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].mole_frac_comp["CO2"]
        ] = 62.020553614342035
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].mole_frac_comp["CO"]
        ] = 9.628044991205439
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].mole_frac_comp["H2O"]
        ] = 5.437801122608198
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].mole_frac_comp["H2"]
        ] = 10.407948041917045
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].mole_frac_comp["N2"]
        ] = 1.6665816449080268
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].mole_frac_comp["NH3"]
        ] = 764299.3478545976
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].mole_frac_comp["O2"]
        ] = 348491.5863157065
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[
                "Vap", "CH4"
            ]
        ] = 1
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[
                "Vap", "CO2"
            ]
        ] = 62.020553614342035
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[
                "Vap", "CO"
            ]
        ] = 9.628044991205439
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[
                "Vap", "H2O"
            ]
        ] = 5.437801122608198
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[
                "Vap", "H2"
            ]
        ] = 10.407948041917045
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[
                "Vap", "N2"
            ]
        ] = 1.6665816449080268
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[
                "Vap", "NH3"
            ]
        ] = 764299.3478545976
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[
                "Vap", "O2"
            ]
        ] = 348491.5863157065
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].pressure
        ] = 9.869232667160129e-06
        methane.scaling_factor[
            methane.fs.unit.control_volume.properties_out[0.0].temperature
        ] = 0.0004275929186958249
        methane.scaling_factor[
            methane.fs.unit.gibbs_minimization[0.0, "Vap", "CH4"]
        ] = 8.074436602893497e-07
        methane.scaling_factor[
            methane.fs.unit.gibbs_minimization[0.0, "Vap", "CO2"]
        ] = 7.961846197368761e-07
        methane.scaling_factor[methane.fs.unit.gibbs_minimization[0.0, "Vap", "CO"]] = (
            1.306639729794169e-06
        )
        methane.scaling_factor[
            methane.fs.unit.gibbs_minimization[0.0, "Vap", "H2O"]
        ] = 1.1245365150847992e-06
        methane.scaling_factor[methane.fs.unit.gibbs_minimization[0.0, "Vap", "H2"]] = (
            1.9650913378522065e-06
        )
        methane.scaling_factor[methane.fs.unit.gibbs_minimization[0.0, "Vap", "N2"]] = (
            1.5268799321929755e-06
        )
        methane.scaling_factor[
            methane.fs.unit.gibbs_minimization[0.0, "Vap", "NH3"]
        ] = 1.0314993300171925e-06
        methane.scaling_factor[methane.fs.unit.gibbs_minimization[0.0, "Vap", "O2"]] = (
            9.653850829168712e-07
        )
        methane.scaling_factor[methane.fs.unit.lagrange_mult[0.0, "C"]] = (
            3.915957066473961e-06
        )
        methane.scaling_factor[methane.fs.unit.lagrange_mult[0.0, "H"]] = (
            5.55811764258496e-06
        )
        methane.scaling_factor[methane.fs.unit.lagrange_mult[0.0, "N"]] = (
            4.318668616445236e-06
        )
        methane.scaling_factor[methane.fs.unit.lagrange_mult[0.0, "O"]] = (
            2.730521354347429e-06
        )

        scaling = TransformationFactory("core.scale_model")
        sm = scaling.create_using(methane, rename=False)

        results = solver.solve(sm)

        scaling.propagate_solution(sm, methane)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution_duty(self, methane):
        methane.fs.unit.outlet.display()
        assert pytest.approx(250.06, abs=1e-1) == value(
            methane.fs.unit.outlet.flow_mol[0]
        )
        assert pytest.approx(0.0, abs=1e-4) == value(
            methane.fs.unit.outlet.mole_frac_comp[0, "CH4"]
        )
        assert pytest.approx(0.097367, abs=1e-4) == value(
            methane.fs.unit.outlet.mole_frac_comp[0, "CO"]
        )
        assert pytest.approx(0.022591, abs=1e-4) == value(
            methane.fs.unit.outlet.mole_frac_comp[0, "CO2"]
        )
        assert pytest.approx(0.10301, abs=1e-4) == value(
            methane.fs.unit.outlet.mole_frac_comp[0, "H2"]
        )
        assert pytest.approx(0.17691, abs=1e-4) == value(
            methane.fs.unit.outlet.mole_frac_comp[0, "H2O"]
        )
        assert pytest.approx(0.59989, abs=1e-4) == value(
            methane.fs.unit.outlet.mole_frac_comp[0, "N2"]
        )
        assert pytest.approx(0.0, abs=1e-5) == value(
            methane.fs.unit.outlet.mole_frac_comp[0, "NH3"]
        )
        assert pytest.approx(0.00024128, abs=1e-4) == value(
            methane.fs.unit.outlet.mole_frac_comp[0, "O2"]
        )
        assert pytest.approx(-7454077, abs=1e2) == value(methane.fs.unit.heat_duty[0])
        assert pytest.approx(101325.0, abs=1e-2) == value(
            methane.fs.unit.outlet.pressure[0]
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation_duty(self, methane):
        assert (
            abs(
                value(
                    methane.fs.unit.inlet.flow_mol[0]
                    * methane.fs.unit.control_volume.properties_in[0].enth_mol_phase[
                        "Vap"
                    ]
                    - methane.fs.unit.outlet.flow_mol[0]
                    * methane.fs.unit.control_volume.properties_out[0].enth_mol_phase[
                        "Vap"
                    ]
                    + methane.fs.unit.heat_duty[0]
                )
            )
            <= 1e-4
        )

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_performance_contents(self, methane):
        perf_dict = methane.fs.unit._get_performance_contents()

        assert perf_dict == {
            "vars": {
                "Heat Duty": methane.fs.unit.heat_duty[0],
                "Pressure Change": methane.fs.unit.deltaP[0],
            }
        }


class TestInitializers:
    @pytest.fixture
    def model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = MethaneCombustionParameterBlock()

        m.fs.unit = GibbsReactor(
            property_package=m.fs.properties,
            has_heat_transfer=True,
            has_pressure_change=True,
        )

        m.fs.unit.inlet.flow_mol[0].set_value(230.0)
        m.fs.unit.inlet.mole_frac_comp[0, "H2"].set_value(0.0435)
        m.fs.unit.inlet.mole_frac_comp[0, "N2"].set_value(0.6522)
        m.fs.unit.inlet.mole_frac_comp[0, "O2"].set_value(0.1739)
        m.fs.unit.inlet.mole_frac_comp[0, "CO2"].set_value(1e-5)
        m.fs.unit.inlet.mole_frac_comp[0, "CH4"].set_value(0.1304)
        m.fs.unit.inlet.mole_frac_comp[0, "CO"].set_value(1e-5)
        m.fs.unit.inlet.mole_frac_comp[0, "H2O"].set_value(1e-5)
        m.fs.unit.inlet.mole_frac_comp[0, "NH3"].set_value(1e-5)
        m.fs.unit.inlet.temperature[0].set_value(1500.0)
        m.fs.unit.inlet.pressure[0].set_value(101325.0)

        m.fs.unit.outlet.temperature[0].fix(2844.38)
        m.fs.unit.deltaP.fix(0)

        return m

    @pytest.mark.component
    def test_general_hierarchical(self, model):
        initializer = SingleControlVolumeUnitInitializer()
        initializer.initialize(
            model.fs.unit,
            initial_guesses={
                "control_volume.properties_out[0].pressure": 101325.0,
                "control_volume.properties_out[0].flow_mol": 251.05,
                "control_volume.properties_out[0].mole_frac_comp[CH4]": 1e-5,
                "control_volume.properties_out[0].mole_frac_comp[CO]": 0.0916,
                "control_volume.properties_out[0].mole_frac_comp[CO2]": 0.0281,
                "control_volume.properties_out[0].mole_frac_comp[H2]": 0.1155,
                "control_volume.properties_out[0].mole_frac_comp[H2O]": 0.1633,
                "control_volume.properties_out[0].mole_frac_comp[N2]": 0.59478,
                "control_volume.properties_out[0].mole_frac_comp[NH3]": 1e-5,
                "control_volume.properties_out[0].mole_frac_comp[O2]": 0.0067,
            },
        )

        assert initializer.summary[model.fs.unit]["status"] == InitializationStatus.Ok

        assert pytest.approx(250.06, abs=1e-2) == value(
            model.fs.unit.outlet.flow_mol[0]
        )
        assert pytest.approx(0.0, abs=1e-4) == value(
            model.fs.unit.outlet.mole_frac_comp[0, "CH4"]
        )
        assert pytest.approx(0.0974, abs=1e-4) == value(
            model.fs.unit.outlet.mole_frac_comp[0, "CO"]
        )
        assert pytest.approx(0.0226, abs=1e-4) == value(
            model.fs.unit.outlet.mole_frac_comp[0, "CO2"]
        )
        assert pytest.approx(0.1030, abs=1e-4) == value(
            model.fs.unit.outlet.mole_frac_comp[0, "H2"]
        )
        assert pytest.approx(0.1769, abs=1e-4) == value(
            model.fs.unit.outlet.mole_frac_comp[0, "H2O"]
        )
        assert pytest.approx(0.5999, abs=1e-4) == value(
            model.fs.unit.outlet.mole_frac_comp[0, "N2"]
        )
        assert pytest.approx(0.0, abs=1e-5) == value(
            model.fs.unit.outlet.mole_frac_comp[0, "NH3"]
        )
        assert pytest.approx(0.0002, abs=1e-4) == value(
            model.fs.unit.outlet.mole_frac_comp[0, "O2"]
        )
        assert pytest.approx(-7454077, abs=1e2) == value(model.fs.unit.heat_duty[0])
        assert pytest.approx(101325.0, abs=1e-2) == value(
            model.fs.unit.outlet.pressure[0]
        )

        assert not model.fs.unit.inlet.flow_mol[0].fixed
        assert not model.fs.unit.inlet.mole_frac_comp[0, "H2"].fixed
        assert not model.fs.unit.inlet.mole_frac_comp[0, "N2"].fixed
        assert not model.fs.unit.inlet.mole_frac_comp[0, "O2"].fixed
        assert not model.fs.unit.inlet.mole_frac_comp[0, "CO2"].fixed
        assert not model.fs.unit.inlet.mole_frac_comp[0, "CH4"].fixed
        assert not model.fs.unit.inlet.mole_frac_comp[0, "CO"].fixed
        assert not model.fs.unit.inlet.mole_frac_comp[0, "H2O"].fixed
        assert not model.fs.unit.inlet.mole_frac_comp[0, "NH3"].fixed
        assert not model.fs.unit.inlet.temperature[0].fixed
        assert not model.fs.unit.inlet.pressure[0].fixed

    @pytest.mark.component
    def test_block_triangularization(self, model):
        initializer = BlockTriangularizationInitializer(constraint_tolerance=2e-5)
        initializer.initialize(
            model.fs.unit,
            initial_guesses={
                "control_volume.properties_out[0].pressure": 101325.0,
                "control_volume.properties_out[0].flow_mol": 251.05,
                "control_volume.properties_out[0].mole_frac_comp[CH4]": 1e-5,
                "control_volume.properties_out[0].mole_frac_comp[CO]": 0.0916,
                "control_volume.properties_out[0].mole_frac_comp[CO2]": 0.0281,
                "control_volume.properties_out[0].mole_frac_comp[H2]": 0.1155,
                "control_volume.properties_out[0].mole_frac_comp[H2O]": 0.1633,
                "control_volume.properties_out[0].mole_frac_comp[N2]": 0.59478,
                "control_volume.properties_out[0].mole_frac_comp[NH3]": 1e-5,
                "control_volume.properties_out[0].mole_frac_comp[O2]": 0.0067,
            },
        )

        assert initializer.summary[model.fs.unit]["status"] == InitializationStatus.Ok

        assert pytest.approx(250.06, abs=1e-2) == value(
            model.fs.unit.outlet.flow_mol[0]
        )
        assert pytest.approx(0.0, abs=1e-4) == value(
            model.fs.unit.outlet.mole_frac_comp[0, "CH4"]
        )
        assert pytest.approx(0.0974, abs=1e-4) == value(
            model.fs.unit.outlet.mole_frac_comp[0, "CO"]
        )
        assert pytest.approx(0.0226, abs=1e-4) == value(
            model.fs.unit.outlet.mole_frac_comp[0, "CO2"]
        )
        assert pytest.approx(0.1030, abs=1e-4) == value(
            model.fs.unit.outlet.mole_frac_comp[0, "H2"]
        )
        assert pytest.approx(0.1769, abs=1e-4) == value(
            model.fs.unit.outlet.mole_frac_comp[0, "H2O"]
        )
        assert pytest.approx(0.5999, abs=1e-4) == value(
            model.fs.unit.outlet.mole_frac_comp[0, "N2"]
        )
        assert pytest.approx(0.0, abs=1e-5) == value(
            model.fs.unit.outlet.mole_frac_comp[0, "NH3"]
        )
        assert pytest.approx(0.0002, abs=1e-4) == value(
            model.fs.unit.outlet.mole_frac_comp[0, "O2"]
        )
        assert pytest.approx(-7454077, abs=1e2) == value(model.fs.unit.heat_duty[0])
        assert pytest.approx(101325.0, abs=1e-2) == value(
            model.fs.unit.outlet.pressure[0]
        )

        assert not model.fs.unit.inlet.flow_mol[0].fixed
        assert not model.fs.unit.inlet.mole_frac_comp[0, "H2"].fixed
        assert not model.fs.unit.inlet.mole_frac_comp[0, "N2"].fixed
        assert not model.fs.unit.inlet.mole_frac_comp[0, "O2"].fixed
        assert not model.fs.unit.inlet.mole_frac_comp[0, "CO2"].fixed
        assert not model.fs.unit.inlet.mole_frac_comp[0, "CH4"].fixed
        assert not model.fs.unit.inlet.mole_frac_comp[0, "CO"].fixed
        assert not model.fs.unit.inlet.mole_frac_comp[0, "H2O"].fixed
        assert not model.fs.unit.inlet.mole_frac_comp[0, "NH3"].fixed
        assert not model.fs.unit.inlet.temperature[0].fixed
        assert not model.fs.unit.inlet.pressure[0].fixed
