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
Tests for Extractor unit model.
Authors: Andrew Lee
"""

import pytest

from pyomo.environ import (
    assert_optimal_termination,
    ConcreteModel,
    Constraint,
    RangeSet,
    Set,
    units,
    value,
    Var,
)
from pyomo.network import Port
from pyomo.common.config import ConfigBlock
from pyomo.util.check_units import assert_units_consistent, assert_units_equivalent

from idaes.core import (
    FlowsheetBlock,
    FlowDirection,
    declare_process_block_class,
    PhysicalParameterBlock,
    StateBlock,
    StateBlockData,
    Component,
    Phase,
    MaterialFlowBasis,
)
from idaes.models.unit_models.extractor import ExtractorCascade, ExtractorCascadeData
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.solvers import get_solver
from idaes.core.util.exceptions import ConfigurationError


# -----------------------------------------------------------------------------
# Property packages for testing
@declare_process_block_class("Parameters1")
class Parameter1Data(PhysicalParameterBlock):
    def build(self):
        super().build()

        self.phase1 = Phase()

        self.solvent1 = Component()
        self.solute1 = Component()
        self.solute2 = Component()
        self.solute3 = Component()

        self._state_block_class = StateBlock1

    @classmethod
    def define_metadata(cls, obj):
        obj.add_default_units(
            {
                "time": units.s,
                "length": units.m,
                "mass": units.kg,
                "amount": units.mol,
                "temperature": units.K,
            }
        )


class SBlock1Base(StateBlock):
    def initialize(blk, **kwargs):
        pass

    def release_state(blk, **kwargs):
        pass


@declare_process_block_class("StateBlock1", block_class=SBlock1Base)
class StateBlock1Data(StateBlockData):
    CONFIG = ConfigBlock(implicit=True)

    def build(self):
        super().build()

        self.flow_mol_phase_comp = Var(
            self.phase_component_set,
            units=units.mol / units.s,
        )
        self.enth_flow = Var(
            units=units.J / units.s,
        )
        self.pressure = Var(units=units.Pa)

    def get_material_flow_terms(self, p, j):
        return self.flow_mol_phase_comp[p, j]

    def get_enthalpy_flow_terms(self, p):
        return self.enth_flow

    def get_material_flow_basis(self):
        return MaterialFlowBasis.molar

    def define_state_vars(self):
        return {
            "flow_mol_phase_comp": self.flow_mol_phase_comp,
            "enth_flow": self.enth_flow,
            "pressure": self.pressure,
        }


@declare_process_block_class("Parameters2")
class Parameter2Data(PhysicalParameterBlock):
    def build(self):
        super().build()

        self.phase1 = Phase()

        self.solvent2 = Component()
        self.solute1 = Component()
        self.solute2 = Component()

        self._state_block_class = StateBlock2

    @classmethod
    def define_metadata(cls, obj):
        obj.add_default_units(
            {
                "time": units.s,
                "length": units.m,
                "mass": units.kg,
                "amount": units.mol,
                "temperature": units.K,
            }
        )


class SBlock2Base(StateBlock):
    def initialize(blk, **kwargs):
        pass

    def release_state(blk, **kwargs):
        pass


@declare_process_block_class("StateBlock2", block_class=SBlock1Base)
class StateBlock2Data(StateBlockData):
    CONFIG = ConfigBlock(implicit=True)

    def build(self):
        super().build()

        self.flow_mol_phase_comp = Var(
            self.phase_component_set,
            units=units.mol / units.s,
        )
        self.enth_flow = Var(
            units=units.J / units.s,
        )
        self.pressure = Var(units=units.Pa)

    def get_material_flow_terms(self, p, j):
        return self.flow_mol_phase_comp[p, j]

    def get_enthalpy_flow_terms(self, p):
        return self.enth_flow

    def get_material_flow_basis(self):
        return MaterialFlowBasis.molar

    def define_state_vars(self):
        return {
            "flow_mol_phase_comp": self.flow_mol_phase_comp,
            "enth_flow": self.enth_flow,
            "pressure": self.pressure,
        }


@declare_process_block_class("Parameters3")
class Parameter3Data(PhysicalParameterBlock):
    def build(self):
        super().build()

        self.phase1 = Phase()

        self.solvent4 = Component()

        self._state_block_class = StateBlock3

    @classmethod
    def define_metadata(cls, obj):
        obj.add_default_units(
            {
                "time": units.s,
                "length": units.m,
                "mass": units.kg,
                "amount": units.mol,
                "temperature": units.K,
            }
        )


class SBlock3Base(StateBlock):
    def initialize(blk, **kwargs):
        pass

    def release_state(blk, **kwargs):
        pass


@declare_process_block_class("StateBlock3", block_class=SBlock1Base)
class StateBlock3Data(StateBlockData):
    CONFIG = ConfigBlock(implicit=True)

    def build(self):
        super().build()

    def get_material_flow_basis(self):
        return MaterialFlowBasis.mass


# -----------------------------------------------------------------------------
# Frame class for unit testing
@declare_process_block_class("ECFrame")
class ECFrameData(ExtractorCascadeData):
    def build(self):
        super(ExtractorCascadeData, self).build()


# -----------------------------------------------------------------------------
class TestBuild:
    @pytest.fixture
    def model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties1 = Parameters1()
        m.fs.properties2 = Parameters2()

        m.fs.unit = ECFrame(
            number_of_stages=2,
            streams={
                "stream1": {"property_package": m.fs.properties1},
                "stream2": {
                    "property_package": m.fs.properties2,
                    "flow_direction": FlowDirection.backward,
                },
            },
        )

        return m

    @pytest.mark.unit
    def test_config(self, model):
        assert not model.fs.unit.config.dynamic
        assert not model.fs.unit.config.has_holdup
        assert model.fs.unit.config.number_of_stages == 2
        assert "stream1" in model.fs.unit.config.streams
        assert "stream2" in model.fs.unit.config.streams
        assert model.fs.unit.config.interacting_streams is None

    @pytest.mark.unit
    def test_verify_inputs_too_few_streams(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties1 = Parameters1()

        m.fs.unit = ECFrame(
            number_of_stages=2,
            streams={
                "stream1": {"property_package": m.fs.properties1},
            },
        )

        with pytest.raises(
            ConfigurationError,
            match="ExtractorCascade models must define at least two streams; received "
            "\['stream1'\]",
        ):
            m.fs.unit._verify_inputs()

    @pytest.mark.unit
    def test_verify_inputs_dynamic(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=True, time_units=units.s)

        m.fs.properties1 = Parameters1()

        m.fs.unit = ECFrame(
            number_of_stages=2,
            streams={
                "stream1": {"property_package": m.fs.properties1},
                "stream2": {"property_package": m.fs.properties1},
            },
        )

        with pytest.raises(
            NotImplementedError,
            match="ExtractorCascade model does not support dynamics yet.",
        ):
            m.fs.unit._verify_inputs()

    @pytest.mark.unit
    def test_verify_inputs_no_common_components(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties1 = Parameters1()
        m.fs.properties3 = Parameters3()

        m.fs.unit = ECFrame(
            number_of_stages=2,
            streams={
                "stream1": {"property_package": m.fs.properties1},
                "stream2": {"property_package": m.fs.properties3},
            },
        )

        with pytest.raises(
            ConfigurationError,
            match="No common components found in property packages. Extractor "
            "model assumes mass transfer occurs between components with the "
            "same name in different streams.",
        ):
            m.fs.unit._verify_inputs()

    @pytest.mark.unit
    def test_verify_inputs_construct_components(self, model):
        model.fs.unit._verify_inputs()

        assert isinstance(model.fs.unit.stages, RangeSet)
        assert len(model.fs.unit.stages) == 2

        assert isinstance(model.fs.unit.stream_component_interactions, Set)
        # One stream pair with two common components
        assert len(model.fs.unit.stream_component_interactions) == 2
        for k in model.fs.unit.stream_component_interactions:
            assert k in [
                ("stream1", "stream2", "solute1"),
                ("stream1", "stream2", "solute2"),
            ]

    @pytest.mark.unit
    def test_build_state_blocks(self, model):
        model.fs.unit._verify_inputs()
        flow_basis, uom = model.fs.unit._build_state_blocks()

        assert flow_basis == MaterialFlowBasis.molar
        assert uom == model.fs.properties1.get_metadata().derived_units

        assert isinstance(model.fs.unit.stream1, StateBlock1)
        assert len(model.fs.unit.stream1) == 2
        assert not model.fs.unit.stream1[0, 1].config.defined_state
        assert not model.fs.unit.stream1[0, 2].config.defined_state

        assert isinstance(model.fs.unit.stream1_inlet_state, StateBlock1)
        assert len(model.fs.unit.stream1_inlet_state) == 1
        assert model.fs.unit.stream1_inlet_state[0].config.defined_state

        assert isinstance(model.fs.unit.stream2, StateBlock2)
        assert len(model.fs.unit.stream2) == 2
        assert not model.fs.unit.stream2[0, 1].config.defined_state
        assert not model.fs.unit.stream2[0, 2].config.defined_state

        assert isinstance(model.fs.unit.stream2_inlet_state, StateBlock2)
        assert len(model.fs.unit.stream2_inlet_state) == 1
        assert model.fs.unit.stream2_inlet_state[0].config.defined_state

    @pytest.mark.unit
    def test_build_state_blocks_different_flow_basis(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties1 = Parameters1()
        m.fs.properties3 = Parameters3()

        m.fs.unit = ECFrame(
            number_of_stages=2,
            streams={
                "stream1": {"property_package": m.fs.properties1},
                # Properties3 has a different flow basis
                "stream2": {"property_package": m.fs.properties3},
            },
        )

        m.fs.unit.stages = Set(initialize=[1, 2])

        with pytest.raises(
            ConfigurationError,
            match="Property packages use different flow bases: ExtractionCascade "
            "requires all property packages to use the same basis. stream2 uses "
            "MaterialFlowBasis.mass, whilst first stream uses "
            "MaterialFlowBasis.molar.",
        ):
            m.fs.unit._build_state_blocks()

    @pytest.mark.unit
    def test_material_balances(self, model):
        model.fs.unit._verify_inputs()
        flow_basis, uom = model.fs.unit._build_state_blocks()
        model.fs.unit._build_material_balance_constraints(flow_basis, uom)

        assert isinstance(model.fs.unit.material_transfer_term, Var)
        # One stream pair with two common components over two stages and 1 time point
        assert len(model.fs.unit.material_transfer_term) == 4
        assert_units_equivalent(
            model.fs.unit.material_transfer_term._units, units.mol / units.s
        )

        assert isinstance(model.fs.unit.stream1_material_balance, Constraint)
        # 1 time point, 2 stages, 4 components
        assert len(model.fs.unit.stream1_material_balance) == 8

        for j in ["solvent1", "solute3"]:  # no mass transfer, forward flow
            assert str(model.fs.unit.stream1_material_balance[0, 1, j]._expr) == str(
                0
                == model.fs.unit.stream1_inlet_state[0].flow_mol_phase_comp["phase1", j]
                - model.fs.unit.stream1[0, 1].flow_mol_phase_comp["phase1", j]
            )
            assert str(model.fs.unit.stream1_material_balance[0, 2, j]._expr) == str(
                0
                == model.fs.unit.stream1[0, 1].flow_mol_phase_comp["phase1", j]
                - model.fs.unit.stream1[0, 2].flow_mol_phase_comp["phase1", j]
            )
        for j in ["solute1", "solute2"]:  # has +ve mass transfer, forward flow
            assert str(model.fs.unit.stream1_material_balance[0, 1, j]._expr) == str(
                0
                == model.fs.unit.stream1_inlet_state[0].flow_mol_phase_comp["phase1", j]
                - model.fs.unit.stream1[0, 1].flow_mol_phase_comp["phase1", j]
                + model.fs.unit.material_transfer_term[0, 1, "stream1", "stream2", j]
            )
            assert str(model.fs.unit.stream1_material_balance[0, 2, j]._expr) == str(
                0
                == model.fs.unit.stream1[0, 1].flow_mol_phase_comp["phase1", j]
                - model.fs.unit.stream1[0, 2].flow_mol_phase_comp["phase1", j]
                + model.fs.unit.material_transfer_term[0, 2, "stream1", "stream2", j]
            )

        assert isinstance(model.fs.unit.stream2_material_balance, Constraint)
        # 1 time point, 2 stages, 3 components
        assert len(model.fs.unit.stream2_material_balance) == 6
        for j in ["solvent2"]:  # no mass transfer, reverse flow
            assert str(model.fs.unit.stream2_material_balance[0, 2, j]._expr) == str(
                0
                == model.fs.unit.stream2_inlet_state[0].flow_mol_phase_comp["phase1", j]
                - model.fs.unit.stream2[0, 2].flow_mol_phase_comp["phase1", j]
            )
            assert str(model.fs.unit.stream2_material_balance[0, 1, j]._expr) == str(
                0
                == model.fs.unit.stream2[0, 2].flow_mol_phase_comp["phase1", j]
                - model.fs.unit.stream2[0, 1].flow_mol_phase_comp["phase1", j]
            )
        for j in ["solute1", "solute2"]:  # has -ve mass transfer, reverse flow
            assert str(model.fs.unit.stream2_material_balance[0, 2, j]._expr) == str(
                0
                == model.fs.unit.stream2_inlet_state[0].flow_mol_phase_comp["phase1", j]
                - model.fs.unit.stream2[0, 2].flow_mol_phase_comp["phase1", j]
                - model.fs.unit.material_transfer_term[0, 2, "stream1", "stream2", j]
            )
            assert str(model.fs.unit.stream2_material_balance[0, 1, j]._expr) == str(
                0
                == model.fs.unit.stream2[0, 2].flow_mol_phase_comp["phase1", j]
                - model.fs.unit.stream2[0, 1].flow_mol_phase_comp["phase1", j]
                - model.fs.unit.material_transfer_term[0, 1, "stream1", "stream2", j]
            )

    @pytest.mark.unit
    def test_energy_balances(self, model):
        model.fs.unit._verify_inputs()
        _, uom = model.fs.unit._build_state_blocks()
        model.fs.unit._build_energy_balance_constraints(uom)
        # TODO: Check transfer terms

        assert isinstance(model.fs.unit.stream1_energy_balance, Constraint)
        # 1 time point, 2 stages
        assert len(model.fs.unit.stream1_energy_balance) == 2

        assert str(model.fs.unit.stream1_energy_balance[0, 1]._expr) == str(
            0
            == units.convert(
                model.fs.unit.stream1_inlet_state[0].enth_flow
                - model.fs.unit.stream1[0, 1].enth_flow,
                units.kg * units.m**2 / units.s**3,
            )
        )
        assert str(model.fs.unit.stream1_energy_balance[0, 2]._expr) == str(
            0
            == units.convert(
                model.fs.unit.stream1[0, 1].enth_flow
                - model.fs.unit.stream1[0, 2].enth_flow,
                units.kg * units.m**2 / units.s**3,
            )
        )

        assert isinstance(model.fs.unit.stream2_energy_balance, Constraint)
        # 1 time point, 2 stages
        assert len(model.fs.unit.stream2_energy_balance) == 2

        assert str(model.fs.unit.stream2_energy_balance[0, 2]._expr) == str(
            0
            == units.convert(
                model.fs.unit.stream2_inlet_state[0].enth_flow
                - model.fs.unit.stream2[0, 2].enth_flow,
                units.kg * units.m**2 / units.s**3,
            )
        )
        assert str(model.fs.unit.stream2_energy_balance[0, 1]._expr) == str(
            0
            == units.convert(
                model.fs.unit.stream2[0, 2].enth_flow
                - model.fs.unit.stream2[0, 1].enth_flow,
                units.kg * units.m**2 / units.s**3,
            )
        )

    @pytest.mark.unit
    def test_pressure_balances(self, model):
        model.fs.unit._verify_inputs()
        _, uom = model.fs.unit._build_state_blocks()
        model.fs.unit._build_pressure_balance_constraints(uom)
        # TODO: Check deltaP terms

        assert isinstance(model.fs.unit.stream1_pressure_balance, Constraint)
        # 1 time point, 2 stages
        assert len(model.fs.unit.stream1_pressure_balance) == 2

        assert str(model.fs.unit.stream1_pressure_balance[0, 1]._expr) == str(
            0
            == units.convert(
                model.fs.unit.stream1_inlet_state[0].pressure
                - model.fs.unit.stream1[0, 1].pressure,
                units.kg / units.m / units.s**2,
            )
        )
        assert str(model.fs.unit.stream1_pressure_balance[0, 2]._expr) == str(
            0
            == units.convert(
                model.fs.unit.stream1[0, 1].pressure
                - model.fs.unit.stream1[0, 2].pressure,
                units.kg / units.m / units.s**2,
            )
        )

        assert isinstance(model.fs.unit.stream2_pressure_balance, Constraint)
        # 1 time point, 2 stages
        assert len(model.fs.unit.stream2_pressure_balance) == 2

        assert str(model.fs.unit.stream2_pressure_balance[0, 2]._expr) == str(
            0
            == units.convert(
                model.fs.unit.stream2_inlet_state[0].pressure
                - model.fs.unit.stream2[0, 2].pressure,
                units.kg / units.m / units.s**2,
            )
        )
        assert str(model.fs.unit.stream2_pressure_balance[0, 1]._expr) == str(
            0
            == units.convert(
                model.fs.unit.stream2[0, 2].pressure
                - model.fs.unit.stream2[0, 1].pressure,
                units.kg / units.m / units.s**2,
            )
        )

    @pytest.mark.unit
    def test_ports(self, model):
        model.fs.unit._verify_inputs()
        model.fs.unit._build_state_blocks()
        model.fs.unit._build_ports()

        assert isinstance(model.fs.unit.stream1_inlet, Port)
        for p, j in model.fs.unit.stream1.phase_component_set:
            assert (
                model.fs.unit.stream1_inlet.flow_mol_phase_comp[0, p, j]
                is model.fs.unit.stream1_inlet_state[0].flow_mol_phase_comp[p, j]
            )
        assert (
            model.fs.unit.stream1_inlet.enth_flow[0]
            is model.fs.unit.stream1_inlet_state[0].enth_flow
        )

        assert isinstance(model.fs.unit.stream1_outlet, Port)
        for p, j in model.fs.unit.stream1.phase_component_set:
            assert (
                model.fs.unit.stream1_outlet.flow_mol_phase_comp[0, p, j]
                is model.fs.unit.stream1[0, 2].flow_mol_phase_comp[p, j]
            )
        assert (
            model.fs.unit.stream1_outlet.enth_flow[0]
            is model.fs.unit.stream1[0, 2].enth_flow
        )

        assert isinstance(model.fs.unit.stream2_inlet, Port)
        for p, j in model.fs.unit.stream2.phase_component_set:
            assert (
                model.fs.unit.stream2_inlet.flow_mol_phase_comp[0, p, j]
                is model.fs.unit.stream2_inlet_state[0].flow_mol_phase_comp[p, j]
            )
        assert (
            model.fs.unit.stream2_inlet.enth_flow[0]
            is model.fs.unit.stream2_inlet_state[0].enth_flow
        )

        assert isinstance(model.fs.unit.stream2_outlet, Port)
        for p, j in model.fs.unit.stream2.phase_component_set:
            assert (
                model.fs.unit.stream2_outlet.flow_mol_phase_comp[0, p, j]
                is model.fs.unit.stream2[0, 1].flow_mol_phase_comp[p, j]
            )
        assert (
            model.fs.unit.stream2_outlet.enth_flow[0]
            is model.fs.unit.stream2[0, 1].enth_flow
        )


class TestToyProblem:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties1 = Parameters1()
        m.fs.properties2 = Parameters2()

        m.fs.unit = ExtractorCascade(
            number_of_stages=2,
            streams={
                "stream1": {"property_package": m.fs.properties1},
                "stream2": {
                    "property_package": m.fs.properties2,
                    "flow_direction": FlowDirection.backward,
                },
            },
        )

        return m

    @pytest.mark.unit
    def test_degrees_of_freedom(self, model):
        # Expect 15 DoF: 6 stream1 inlets, 5 stream2 inlets and 2x2 mass transfer terms
        assert (degrees_of_freedom(model)) == 15

    @pytest.mark.component
    def test_unit_consistency(self, model):
        assert_units_consistent(model)

    @pytest.mark.component
    def test_toy_problem(self, model):
        model.fs.unit.stream1_inlet.flow_mol_phase_comp[0, "phase1", "solvent1"].fix(2)
        model.fs.unit.stream1_inlet.flow_mol_phase_comp[0, "phase1", "solute1"].fix(3)
        model.fs.unit.stream1_inlet.flow_mol_phase_comp[0, "phase1", "solute2"].fix(4)
        model.fs.unit.stream1_inlet.flow_mol_phase_comp[0, "phase1", "solute3"].fix(5)
        model.fs.unit.stream1_inlet.enth_flow[0].fix(5000)
        model.fs.unit.stream1_inlet.pressure[0].fix(1.2e5)

        model.fs.unit.stream2_inlet.flow_mol_phase_comp[0, "phase1", "solvent2"].fix(11)
        model.fs.unit.stream2_inlet.flow_mol_phase_comp[0, "phase1", "solute1"].fix(12)
        model.fs.unit.stream2_inlet.flow_mol_phase_comp[0, "phase1", "solute2"].fix(13)
        model.fs.unit.stream2_inlet.enth_flow[0].fix(7000)
        model.fs.unit.stream2_inlet.pressure[0].fix(2e5)

        model.fs.unit.material_transfer_term[0, :, "stream1", "stream2", "solute1"].fix(
            0.5
        )
        model.fs.unit.material_transfer_term[0, :, "stream1", "stream2", "solute2"].fix(
            -0.5
        )

        assert (degrees_of_freedom(model)) == 0

        solver = get_solver()
        results = solver.solve(model)

        assert_optimal_termination(results)

        assert value(
            model.fs.unit.stream1_outlet.flow_mol_phase_comp[0, "phase1", "solvent1"]
        ) == pytest.approx(2, rel=1e-5)
        assert value(
            model.fs.unit.stream1_outlet.flow_mol_phase_comp[0, "phase1", "solute1"]
        ) == pytest.approx(
            4, rel=1e-5
        )  # 3 + 0.5 + 0.5
        assert value(
            model.fs.unit.stream1_outlet.flow_mol_phase_comp[0, "phase1", "solute2"]
        ) == pytest.approx(
            3, rel=1e-5
        )  # 3 - 0.5 - 0.5
        assert value(
            model.fs.unit.stream1_outlet.flow_mol_phase_comp[0, "phase1", "solute3"]
        ) == pytest.approx(5, rel=1e-5)

        assert value(
            model.fs.unit.stream2_outlet.flow_mol_phase_comp[0, "phase1", "solvent2"]
        ) == pytest.approx(11, rel=1e-5)
        assert value(
            model.fs.unit.stream2_outlet.flow_mol_phase_comp[0, "phase1", "solute1"]
        ) == pytest.approx(
            11, rel=1e-5
        )  # 12 - 0.5 - 0.5
        assert value(
            model.fs.unit.stream2_outlet.flow_mol_phase_comp[0, "phase1", "solute2"]
        ) == pytest.approx(
            14, rel=1e-5
        )  # 13 + 0.5 + 0.5

        assert value(model.fs.unit.stream1_outlet.enth_flow[0]) == pytest.approx(
            5000, rel=1e-5
        )
        assert value(model.fs.unit.stream2_outlet.enth_flow[0]) == pytest.approx(
            7000, rel=1e-5
        )

        assert value(model.fs.unit.stream1_outlet.pressure[0]) == pytest.approx(
            1.2e5, rel=1e-5
        )
        assert value(model.fs.unit.stream2_outlet.pressure[0]) == pytest.approx(
            2e5, rel=1e-5
        )
