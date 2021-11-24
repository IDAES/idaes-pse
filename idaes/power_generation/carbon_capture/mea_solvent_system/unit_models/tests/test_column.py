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
Author: Paul Akula, Anuja Deshpande
"""
# Import Python libraries
import pytest

# Import Pyomo libraries
from pyomo.environ import \
    ConcreteModel, value, SolverStatus, TerminationCondition
from pyomo.util.check_units import assert_units_consistent

# Import IDAES Libraries
from idaes.core import FlowsheetBlock
from idaes.power_generation.carbon_capture.mea_solvent_system.unit_models.column \
    import PackedColumn, ProcessType
from idaes.power_generation.carbon_capture.mea_solvent_system.properties.vapor_prop \
    import VaporParameterBlock
from idaes.power_generation.carbon_capture.mea_solvent_system.properties.liquid_prop \
    import LiquidParameterBlock

from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester
from idaes.core.util import get_solver


# -----------------------------------------------------------------------------
solver = get_solver()


class TestAbsorberSS:
    @pytest.fixture(scope="class")
    def model(self):
        # Spacial domain finite elemets and finite element list
        x_nfe = 10
        x_nfe_list = [i / x_nfe for i in range(x_nfe + 1)]

        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        # Set up property package
        m.fs.vapor_properties = VaporParameterBlock()
        m.fs.liquid_properties = LiquidParameterBlock()

        # Create instance of column on flowsheet
        m.fs.unit = PackedColumn(default={
            "process_type": ProcessType.absorber,
            "finite_elements": x_nfe,
            "length_domain_set": x_nfe_list,
            "transformation_method": "dae.finite_difference",
            "vapor_side": {
                "transformation_scheme": "BACKWARD",
                "property_package": m.fs.vapor_properties,
                "has_pressure_change": False},
            "liquid_side":
            {
                "transformation_scheme": "FORWARD",
                "property_package": m.fs.liquid_properties
            }})

        # Fix  input variables
        m.fs.unit.diameter_column.fix(0.64135)
        m.fs.unit.length_column.fix(18.15)

        # Vapor
        m.fs.unit.vapor_inlet.flow_mol.fix(21.48)
        m.fs.unit.vapor_inlet.temperature.fix(317.88)
        m.fs.unit.vapor_inlet.pressure.fix(107650)
        m.fs.unit.vapor_inlet.mole_frac_comp[0, "CO2"].fix(0.11453)
        m.fs.unit.vapor_inlet.mole_frac_comp[0, "H2O"].fix(0.08526)
        m.fs.unit.vapor_inlet.mole_frac_comp[0, "N2"].fix(0.73821)
        m.fs.unit.vapor_inlet.mole_frac_comp[0, "O2"].fix(0.06200)
        # Liquid
        m.fs.unit.liquid_inlet.flow_mol.fix(37.55)
        m.fs.unit.liquid_inlet.temperature.fix(319.87)
        m.fs.unit.liquid_inlet.mole_frac_comp[0, "CO2"].fix(0.00963)
        m.fs.unit.liquid_inlet.mole_frac_comp[0, "H2O"].fix(0.87435)
        m.fs.unit.liquid_inlet.mole_frac_comp[0, "MEA"].fix(0.11602)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.dynamic is False
        assert model.fs.unit.config.liquid_side.transformation_scheme == \
            'FORWARD'
        assert model.fs.unit.config.vapor_side.transformation_scheme == \
            'BACKWARD'

        assert hasattr(model.fs.unit, "vapor_inlet")
        assert hasattr(model.fs.unit, "vapor_outlet")
        assert hasattr(model.fs.unit, "liquid_inlet")
        assert hasattr(model.fs.unit, "liquid_outlet")

    @pytest.mark.unit
    def test_degrees_of_freedom(self, model):
        assert degrees_of_freedom(model) == 0

    @pytest.mark.component
    def test_units(self, model):
        assert_units_consistent(model)

    @pytest.mark.component
    def test_initialization(self, model):
        initialization_tester(model)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, model):
        res = solver.solve(model)

        # Solver status/condition
        assert res.solver.termination_condition == TerminationCondition.optimal
        assert res.solver.status == SolverStatus.ok

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, model):
        # model.fs.unit.velocity_vap.display()
        # model.fs.unit.velocity_liq.display()
        # model.fs.unit.area_interfacial.display()
        # model.fs.unit.holdup_liq.display()
        # model.fs.unit.holdup_vap.display()
        model.fs.unit.h_v.display()
        # model.fs.unit.k_l_CO2.display()
        # model.fs.unit.phi.display()
        # model.fs.unit.pressure_equil.display()
        assert False

        assert pytest.approx(24.5379, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.flow_mol[0])
        assert pytest.approx(0.0275158, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.mole_frac_comp[0, "CO2"])
        assert pytest.approx(0.271997, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.mole_frac_comp[0, "H2O"])
        assert pytest.approx(0.646214, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.mole_frac_comp[0, "N2"])
        assert pytest.approx(0.0542735, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.mole_frac_comp[0, "O2"])
        assert pytest.approx(107650, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.pressure[0])
        assert pytest.approx(346.177, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.temperature[0])

        assert pytest.approx(34.4921, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.flow_mol[0])
        assert pytest.approx(0.0622325, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.mole_frac_comp[0, "CO2"])
        assert pytest.approx(0.811462, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.mole_frac_comp[0, "H2O"])
        assert pytest.approx(0.126306, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.mole_frac_comp[0, "MEA"])
        assert pytest.approx(107650, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.pressure[0])
        assert pytest.approx(322.953, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.temperature[0])

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, model):
        vap_in = model.fs.unit.vapor_phase.properties[0, 0]
        vap_out = model.fs.unit.vapor_phase.properties[0, 1]
        liq_in = model.fs.unit.liquid_phase.properties[0, 1]
        liq_out = model.fs.unit.liquid_phase.properties[0, 0]

        # Material conservation
        for j in ["CO2", "H2O"]:
            assert 1e-6 >= abs(value(
                vap_in.get_material_flow_terms("Vap", j) +
                liq_in.get_material_flow_terms("Liq", j) -
                vap_out.get_material_flow_terms("Vap", j) -
                liq_out.get_material_flow_terms("Liq", j)))
        for j in ["N2", "O2"]:
            assert 1e-6 >= abs(value(
                vap_in.get_material_flow_terms("Vap", j) -
                vap_out.get_material_flow_terms("Vap", j)))
        for j in ["MEA"]:
            assert 1e-6 >= abs(value(
                liq_in.get_material_flow_terms("Liq", j) -
                liq_out.get_material_flow_terms("Liq", j)))

        # Energy conservation
        # assert 1e-6 >= abs(value(
        #     vap_in.get_enthalpy_flow_terms("Vap") +
        #     liq_in.get_enthalpy_flow_terms("Liq") -
        #     vap_out.get_enthalpy_flow_terms("Vap") -
        #     liq_out.get_enthalpy_flow_terms("Liq")))

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_profiles(self, model):
        x_vap = {0.0: 0.11453,
                 0.1: 0.11300439732982537,
                 0.2: 0.1110705169224961,
                 0.3: 0.10890070870032177,
                 0.4: 0.10633635548565122,
                 0.5: 0.1031044672408336,
                 0.6: 0.09871166679631208,
                 0.7: 0.09215715243761836,
                 0.8: 0.08116495101562057,
                 0.9: 0.06050758465490609,
                 1.0: 0.027515852017664082}
        T_vap = {0.0: 317.88,
                 0.1: 322.5643364536555,
                 0.2: 324.5743696981395,
                 0.3: 326.48541575084346,
                 0.4: 328.4978967317451,
                 0.5: 330.72705826832924,
                 0.6: 333.3121273697413,
                 0.7: 336.42786487515275,
                 0.8: 340.23215075294877,
                 0.9: 344.47699470956695,
                 1.0: 346.1767147775371}
        x_liq = {0.0: 0.062232533591032965,
                 0.1: 0.06147280117214211,
                 0.2: 0.06055384536417231,
                 0.3: 0.05951095522249928,
                 0.4: 0.058258853253299284,
                 0.5: 0.05664899054741099,
                 0.6: 0.054400803666666456,
                 0.7: 0.05091079166159317,
                 0.8: 0.0446919002682971,
                 0.9: 0.03194066113506298,
                 1.0: 0.00963}
        T_liq = {0.0: 322.95307143292825,
                 0.1: 324.8692840916428,
                 0.2: 326.7964684340596,
                 0.3: 328.8489240507993,
                 0.4: 331.14618628284336,
                 0.5: 333.84199540305985,
                 0.6: 337.13340638921903,
                 0.7: 341.1834367288921,
                 0.8: 345.49693261684604,
                 0.9: 345.2935206597128,
                 1.0: 319.87}

        for x in model.fs.unit.vapor_phase.length_domain:
            assert pytest.approx(x_vap[x], rel=1e-5) == value(
                model.fs.unit.vapor_phase.properties[
                    0, x].mole_frac_comp["CO2"])
            assert pytest.approx(T_vap[x], rel=1e-5) == value(
                model.fs.unit.vapor_phase.properties[0, x].temperature)
            assert pytest.approx(x_liq[x], rel=1e-5) == value(
                model.fs.unit.liquid_phase.properties[
                    0, x].mole_frac_comp["CO2"])
            assert pytest.approx(T_liq[x], rel=1e-5) == value(
                model.fs.unit.liquid_phase.properties[0, x].temperature)


class TestStripperSS:
    @pytest.fixture(scope="class")
    def model(self):
        # Spacial domain finite elemets and finite element list
        x_nfe = 10
        x_nfe_list = [i / x_nfe for i in range(x_nfe + 1)]

        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        # Set up property package
        m.fs.vapor_properties = VaporParameterBlock(
            default={'process_type': ProcessType.stripper})
        m.fs.liquid_properties = LiquidParameterBlock(
            default={'process_type': ProcessType.stripper})

        # Create instance of column on flowsheet
        m.fs.unit = PackedColumn(default={
            "process_type": ProcessType.stripper,
            "finite_elements": x_nfe,
            "length_domain_set": x_nfe_list,
            "transformation_method": "dae.finite_difference",
            "vapor_side": {
                "transformation_scheme": "BACKWARD",
                "property_package": m.fs.vapor_properties,
                "has_pressure_change": False},
            "liquid_side":
            {
                "transformation_scheme": "FORWARD",
                "property_package": m.fs.liquid_properties
            }})

        # Fix  input variables
        m.fs.unit.diameter_column.fix(0.64135)
        m.fs.unit.length_column.fix(12.1)

        # Vapor
        m.fs.unit.vapor_inlet.flow_mol.fix(17.496)
        m.fs.unit.vapor_inlet.temperature.fix(396.6)
        m.fs.unit.vapor_inlet.pressure.fix(183430)
        m.fs.unit.vapor_inlet.mole_frac_comp[0, "CO2"].fix(0.0145)
        m.fs.unit.vapor_inlet.mole_frac_comp[0, "H2O"].fix(0.9855)
        # Liquid
        m.fs.unit.liquid_inlet.flow_mol.fix(84.48)
        m.fs.unit.liquid_inlet.temperature.fix(382.15)
        m.fs.unit.liquid_inlet.mole_frac_comp[0, "CO2"].fix(0.0331)
        m.fs.unit.liquid_inlet.mole_frac_comp[0, "H2O"].fix(0.8547)
        m.fs.unit.liquid_inlet.mole_frac_comp[0, "MEA"].fix(0.1122)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.dynamic is False
        assert model.fs.unit.config.liquid_side.transformation_scheme == \
            'FORWARD'
        assert model.fs.unit.config.vapor_side.transformation_scheme == \
            'BACKWARD'

        assert hasattr(model.fs.unit, "vapor_inlet")
        assert hasattr(model.fs.unit, "vapor_outlet")
        assert hasattr(model.fs.unit, "liquid_inlet")
        assert hasattr(model.fs.unit, "liquid_outlet")

    @pytest.mark.unit
    def test_degrees_of_freedom(self, model):
        assert degrees_of_freedom(model) == 0

    @pytest.mark.component
    def test_units(self, model):
        assert_units_consistent(model)

    @pytest.mark.component
    def test_initialization(self, model):
        initialization_tester(
            model,
            homotopy_steps_h=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1])

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, model):
        res = solver.solve(model)

        # Solver status/condition
        assert res.solver.termination_condition == TerminationCondition.optimal
        assert res.solver.status == SolverStatus.ok

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, model):
        assert pytest.approx(9.38012, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.flow_mol[0])
        assert pytest.approx(0.143431, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.mole_frac_comp[0, "CO2"])
        assert pytest.approx(0.856569, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.mole_frac_comp[0, "H2O"])
        assert pytest.approx(183430, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.pressure[0])
        assert pytest.approx(396.516, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.temperature[0])

        assert pytest.approx(92.5959, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.flow_mol[0])
        assert pytest.approx(0.0184088, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.mole_frac_comp[0, "CO2"])
        assert pytest.approx(0.879225, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.mole_frac_comp[0, "H2O"])
        assert pytest.approx(0.102366, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.mole_frac_comp[0, "MEA"])
        assert pytest.approx(183430, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.pressure[0])
        assert pytest.approx(393.866, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.temperature[0])

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, model):
        vap_in = model.fs.unit.vapor_phase.properties[0, 0]
        vap_out = model.fs.unit.vapor_phase.properties[0, 1]
        liq_in = model.fs.unit.liquid_phase.properties[0, 1]
        liq_out = model.fs.unit.liquid_phase.properties[0, 0]

        # Material conservation
        for j in ["CO2", "H2O"]:
            assert 1e-6 >= abs(value(
                vap_in.get_material_flow_terms("Vap", j) +
                liq_in.get_material_flow_terms("Liq", j) -
                vap_out.get_material_flow_terms("Vap", j) -
                liq_out.get_material_flow_terms("Liq", j)))
        for j in ["MEA"]:
            assert 1e-6 >= abs(value(
                liq_in.get_material_flow_terms("Liq", j) -
                liq_out.get_material_flow_terms("Liq", j)))

        # Energy conservation
        # assert 1e-6 >= abs(value(
        #     vap_in.get_enthalpy_flow_terms("Vap") +
        #     liq_in.get_enthalpy_flow_terms("Liq") -
        #     vap_out.get_enthalpy_flow_terms("Vap") -
        #     liq_out.get_enthalpy_flow_terms("Liq")))

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_profiles(self, model):
        x_vap = {0.0: 0.0145,
                 0.1: 0.024829517200206354,
                 0.2: 0.03394523170577613,
                 0.3: 0.042681395439366006,
                 0.4: 0.051580499284370994,
                 0.5: 0.06109516751748907,
                 0.6: 0.07168363867574917,
                 0.7: 0.08386912065257068,
                 0.8: 0.09829528027395286,
                 0.9: 0.11604469797368286,
                 1.0: 0.14343116744332465}
        T_vap = {0.0: 396.6,
                 0.1: 395.0616706558291,
                 0.2: 394.67382091442073,
                 0.3: 394.410790762689,
                 0.4: 394.189465140404,
                 0.5: 393.9844363857082,
                 0.6: 393.77927146867927,
                 0.7: 393.5593472951173,
                 0.8: 393.3200157928309,
                 0.9: 393.25645193200296,
                 1.0: 396.51607033845767}
        x_liq = {0.0: 0.01840878931275917,
                 0.1: 0.020328445578453303,
                 0.2: 0.021922229234957163,
                 0.3: 0.02336458000118307,
                 0.4: 0.024752361861265264,
                 0.5: 0.02614988027446931,
                 0.6: 0.027606312061395594,
                 0.7: 0.029161777338846732,
                 0.8: 0.03084342381065094,
                 0.9: 0.03259048959589782,
                 1.0: 0.0331}
        T_liq = {0.0: 393.8656274735747,
                 0.1: 393.67873057023945,
                 0.2: 393.4811086735324,
                 0.3: 393.26365847546737,
                 0.4: 393.0170549701345,
                 0.5: 392.7288985133659,
                 0.6: 392.3826549487772,
                 0.7: 391.95489848366833,
                 0.8: 391.3876619344523,
                 0.9: 390.1784251859471,
                 1.0: 382.15}

        for x in model.fs.unit.vapor_phase.length_domain:
            assert pytest.approx(x_vap[x], rel=1e-5) == value(
                model.fs.unit.vapor_phase.properties[
                    0, x].mole_frac_comp["CO2"])
            assert pytest.approx(T_vap[x], rel=1e-5) == value(
                model.fs.unit.vapor_phase.properties[0, x].temperature)
            assert pytest.approx(x_liq[x], rel=1e-5) == value(
                model.fs.unit.liquid_phase.properties[
                    0, x].mole_frac_comp["CO2"])
            assert pytest.approx(T_liq[x], rel=1e-5) == value(
                model.fs.unit.liquid_phase.properties[0, x].temperature)
