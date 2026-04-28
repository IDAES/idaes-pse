import pytest

from pyomo.environ import ConcreteModel, Expression, units as pyunits, Var
from idaes.core import (
    Component,
    declare_process_block_class,
    Phase,
    LiquidPhase,
    SolidPhase,
    VaporPhase,
    PhysicalParameterBlock,
    StateBlock,
    StateBlockData,
)
from idaes.core.scaling import get_scaling_factor, set_scaling_factor
from idaes.models.properties.modular_properties.base.generic_property import (
    ModularPropertiesScaler,
)


@declare_process_block_class("DensityParameters")
class _DensityParameters(PhysicalParameterBlock):
    def build(self):
        super(_DensityParameters, self).build()

        self.p1 = Phase()
        self.p2 = SolidPhase()
        self.p3 = LiquidPhase()
        self.p4 = VaporPhase()
        # self.c1 = Component()
        self.phase_list = ["p1", "p2", "p3", "p4"]

    @classmethod
    def define_metadata(cls, obj):
        obj.define_custom_properties(
            {
                "vol_mol_phase": {"method": "_vol_mol_phase_method"},
                "dens_mol_phase": {"method": "_dens_mol_phase_method"},
            }
        )
        obj.add_default_units(
            {
                "time": pyunits.s,
                "mass": pyunits.kg,
                "length": pyunits.m,
                "amount": pyunits.mol,
                "temperature": pyunits.K,
            }
        )


class _DensityStateBlock(StateBlock):
    pass


class DummyVDScaler(ModularPropertiesScaler):
    def variable_scaling_routine(self):
        raise AssertionError()

    def constraint_scaling_routine(self):
        raise AssertionError()

    def _volume_density_scaling(self, model, overwrite, sf_mw_phase):
        return super()._volume_density_scaling(model, overwrite, sf_mw_phase)

    def _estimate_volume_density_scaling_factors(self, model, phase, sf_mw_phase):
        # This function will be tested separately.
        # Right now we will just test that it was run.
        model._sf_mw_phase = sf_mw_phase

        return (4999, 2063)


@declare_process_block_class("DensityState", block_class=_DensityStateBlock)
class _DensityState(StateBlockData):
    def build(self):
        super(_DensityState, self).build()

    def _vol_mol_phase_method(self):
        @self.Expression(self.phase_list)
        def vol_mol_phase(b, p):
            return 137 * pyunits.m**3 / pyunits.mol

    def _dens_mol_phase_method(self):
        @self.Expression(self.phase_list)
        def dens_mol_phase(b, p):
            return 314159 * pyunits.mol / pyunits.m**3


class TestDensityVolumeScaling(object):
    _sf_mw_phase = {"p1": 1 / 2, "p2": 1 / 5, "p3": 1 / 10, "p4": 1 / 17}

    @pytest.fixture(scope="function")
    def model(self):
        m = ConcreteModel()
        m.props = DensityParameters()
        m.state = DensityState(parameters=m.props)
        return m

    @pytest.mark.parametrize(
        "overwrite",
        [False, True],
    )
    @pytest.mark.parametrize("pass_phase_default_sf", [False, True])
    @pytest.mark.parametrize("pass_overall_default_sf", [False, True])
    @pytest.mark.unit
    def test_properties_not_created(
        self, overwrite, pass_phase_default_sf, pass_overall_default_sf, model
    ):
        scaler_obj = DummyVDScaler()
        if pass_phase_default_sf:
            for j, p in enumerate(model.state.phase_list):
                scaler_obj.default_scaling_factors[f"dens_mol_phase[{p}]"] = (
                    j + 3
                ) ** 3
                scaler_obj.default_scaling_factors[f"vol_mol_phase[{p}]"] = (
                    (j + 5) * (j + 6) / 2
                )
        if pass_overall_default_sf:
            scaler_obj.default_scaling_factors["dens_mol_phase"] = 137
            scaler_obj.default_scaling_factors["vol_mol_phase"] = 137

        assert not model.state.is_property_constructed("vol_mol_phase")
        assert not model.state.is_property_constructed("dens_mol_phase")
        assert not model.state.is_property_constructed("scaling_factor")
        assert not model.state.is_property_constructed("scaling_hint")

        scaler_obj._volume_density_scaling(
            model.state, overwrite=overwrite, sf_mw_phase=self._sf_mw_phase
        )
        if pass_phase_default_sf or pass_overall_default_sf:
            assert not model.state.is_property_constructed("_sf_mw_phase")
        else:
            assert model.state._sf_mw_phase is self._sf_mw_phase

        assert not model.state.is_property_constructed("vol_mol_phase")
        assert not model.state.is_property_constructed("dens_mol_phase")
        assert not model.state.is_property_constructed("scaling_factor")
        assert not model.state.is_property_constructed("scaling_hint")

    @pytest.mark.parametrize("pass_phase_default_sf", [False, True])
    @pytest.mark.parametrize("pass_overall_default_sf", [False, True])
    @pytest.mark.parametrize("construct_other", [False, True])
    @pytest.mark.unit
    def test_only_density_sf_set_overwrite_false(
        self, pass_phase_default_sf, pass_overall_default_sf, construct_other, model
    ):
        for j, p in enumerate(model.state.phase_list):
            set_scaling_factor(model.state.dens_mol_phase[p], j + 1)
        if construct_other:
            model.state.vol_mol_phase
        else:
            assert not model.state.is_property_constructed("vol_mol_phase")
        assert not model.state.is_property_constructed("scaling_factor")
        assert model.state.is_property_constructed("scaling_hint")

        scaler_obj = DummyVDScaler()
        if pass_phase_default_sf:
            for j, p in enumerate(model.state.phase_list):
                scaler_obj.default_scaling_factors[f"dens_mol_phase[{p}]"] = (
                    j + 3
                ) ** 3
                scaler_obj.default_scaling_factors[f"vol_mol_phase[{p}]"] = (
                    (j + 5) * (j + 6) / 2
                )
        if pass_overall_default_sf:
            scaler_obj.default_scaling_factors["dens_mol_phase"] = 137
            scaler_obj.default_scaling_factors["vol_mol_phase"] = 314

        scaler_obj._volume_density_scaling(
            model.state, overwrite=False, sf_mw_phase=self._sf_mw_phase
        )

        if construct_other:
            for j, p in enumerate(model.state.phase_list):
                assert get_scaling_factor(model.state.vol_mol_phase[p]) == 1 / (j + 1)
        else:
            assert not model.state.is_property_constructed("vol_mol_phase")
        assert not model.state.is_property_constructed("scaling_factor")
        assert model.state.is_property_constructed("scaling_hint")
        assert not model.state.is_property_constructed("_sf_mw_phase")

        for j, p in enumerate(model.state.phase_list):
            assert get_scaling_factor(model.state.dens_mol_phase[p]) == j + 1

    @pytest.mark.parametrize("pass_phase_default_sf", [False, True])
    @pytest.mark.parametrize("pass_overall_default_sf", [False, True])
    @pytest.mark.parametrize("construct_other", [False, True])
    @pytest.mark.unit
    def test_only_volume_sf_set_overwrite_false(
        self, pass_phase_default_sf, pass_overall_default_sf, construct_other, model
    ):
        for j, p in enumerate(model.state.phase_list):
            set_scaling_factor(model.state.vol_mol_phase[p], j + 1)
        if construct_other:
            model.state.dens_mol_phase
        else:
            assert not model.state.is_property_constructed("dens_mol_phase")
        assert not model.state.is_property_constructed("scaling_factor")
        assert model.state.is_property_constructed("scaling_hint")

        scaler_obj = DummyVDScaler()

        if pass_phase_default_sf:
            for j, p in enumerate(model.state.phase_list):
                scaler_obj.default_scaling_factors[f"dens_mol_phase[{p}]"] = (
                    j + 3
                ) ** 3
                scaler_obj.default_scaling_factors[f"vol_mol_phase[{p}]"] = (
                    (j + 5) * (j + 6) / 2
                )
        if pass_overall_default_sf:
            scaler_obj.default_scaling_factors["dens_mol_phase"] = 137
            scaler_obj.default_scaling_factors["vol_mol_phase"] = 314

        scaler_obj._volume_density_scaling(
            model.state, overwrite=False, sf_mw_phase=self._sf_mw_phase
        )

        if construct_other:
            for j, p in enumerate(model.state.phase_list):
                assert get_scaling_factor(model.state.dens_mol_phase[p]) == 1 / (j + 1)

        else:
            assert not model.state.is_property_constructed("dens_mol_phase")
        assert not model.state.is_property_constructed("scaling_factor")
        assert model.state.is_property_constructed("scaling_hint")
        assert not model.state.is_property_constructed("_sf_mw_phase")

        for j, p in enumerate(model.state.phase_list):
            assert get_scaling_factor(model.state.vol_mol_phase[p]) == j + 1

    @pytest.mark.parametrize("pass_phase_default_sf", [False, True])
    @pytest.mark.parametrize("pass_overall_default_sf", [False, True])
    @pytest.mark.unit
    def test_both_sf_set_overwrite_false(
        self, pass_phase_default_sf, pass_overall_default_sf, model
    ):
        for j, p in enumerate(model.state.phase_list):
            set_scaling_factor(model.state.dens_mol_phase[p], j + 1)
            set_scaling_factor(model.state.vol_mol_phase[p], 2 * j + 5)
        assert not model.state.is_property_constructed("scaling_factor")
        assert model.state.is_property_constructed("scaling_hint")

        scaler_obj = DummyVDScaler()

        if pass_phase_default_sf:
            for j, p in enumerate(model.state.phase_list):
                scaler_obj.default_scaling_factors[f"dens_mol_phase[{p}]"] = (
                    j + 3
                ) ** 3
                scaler_obj.default_scaling_factors[f"vol_mol_phase[{p}]"] = (
                    (j + 5) * (j + 6) / 2
                )
        if pass_overall_default_sf:
            scaler_obj.default_scaling_factors["dens_mol_phase"] = 137
            scaler_obj.default_scaling_factors["vol_mol_phase"] = 314

        scaler_obj._volume_density_scaling(
            model.state, overwrite=False, sf_mw_phase=self._sf_mw_phase
        )

        assert not model.state.is_property_constructed("scaling_factor")
        assert model.state.is_property_constructed("scaling_hint")
        assert not model.state.is_property_constructed("_sf_mw_phase")

        for j, p in enumerate(model.state.phase_list):
            assert get_scaling_factor(model.state.dens_mol_phase[p]) == j + 1
            assert get_scaling_factor(model.state.vol_mol_phase[p]) == 2 * j + 5

    @pytest.mark.parametrize(
        "overwrite",
        [False, True],
    )
    @pytest.mark.parametrize(
        # Parameterize over two values simultaneously so we don't test
        # the (False, False) combination (which was already tested above)
        "construct_density, construct_volume",
        [(False, True), (True, False), (True, True)],
    )
    @pytest.mark.parametrize("default_density", [False, True])
    @pytest.mark.parametrize("default_volume", [False, True])
    @pytest.mark.parametrize("pass_phase_default_sf", [False, True])
    @pytest.mark.unit
    def test_default_sf(
        self,
        overwrite,
        construct_density,
        construct_volume,
        default_density,
        default_volume,
        pass_phase_default_sf,
        model,
    ):
        if construct_density:
            model.state.dens_mol_phase
            if overwrite:
                set_scaling_factor(model.state.dens_mol_phase["p3"], 5303)
        if construct_volume:
            model.state.vol_mol_phase
            if overwrite:
                set_scaling_factor(model.state.vol_mol_phase["p1"], 6277)

        assert not model.state.is_property_constructed("scaling_factor")
        if overwrite:
            assert model.state.is_property_constructed("scaling_hint")
        else:
            assert not model.state.is_property_constructed("scaling_hint")

        scaler_obj = DummyVDScaler()

        if pass_phase_default_sf:
            for j, p in enumerate(model.state.phase_list):
                if default_density:
                    scaler_obj.default_scaling_factors[f"dens_mol_phase[{p}]"] = (
                        j + 3
                    ) ** 3
                if default_volume:
                    scaler_obj.default_scaling_factors[f"vol_mol_phase[{p}]"] = (
                        (j + 5) * (j + 6) / 2
                    )

        if default_density:
            scaler_obj.default_scaling_factors["dens_mol_phase"] = 137
        if default_volume:
            scaler_obj.default_scaling_factors["vol_mol_phase"] = 314

        scaler_obj._volume_density_scaling(
            model.state, overwrite=overwrite, sf_mw_phase=self._sf_mw_phase
        )

        if pass_phase_default_sf:
            if default_density and default_volume:
                assert not model.state.is_property_constructed("_sf_mw_phase")
                for j, p in enumerate(model.state.phase_list):
                    if construct_density:
                        assert (
                            get_scaling_factor(model.state.dens_mol_phase[p])
                            == (j + 3) ** 3
                        )
                    else:
                        assert not model.state.is_property_constructed("dens_mol_phase")
                    if construct_volume:
                        assert (
                            get_scaling_factor(model.state.vol_mol_phase[p])
                            == (j + 5) * (j + 6) / 2
                        )
                    else:
                        assert not model.state.is_property_constructed("vol_mol_phase")
            elif default_density:
                assert not model.state.is_property_constructed("_sf_mw_phase")
                for j, p in enumerate(model.state.phase_list):
                    if construct_density:
                        assert (
                            get_scaling_factor(model.state.dens_mol_phase[p])
                            == (j + 3) ** 3
                        )
                    else:
                        assert not model.state.is_property_constructed("dens_mol_phase")
                    if construct_volume:
                        assert get_scaling_factor(
                            model.state.vol_mol_phase[p]
                        ) == pytest.approx(1 / (j + 3) ** 3, rel=1e-8)
                    else:
                        assert not model.state.is_property_constructed("vol_mol_phase")

            elif default_volume:
                assert not model.state.is_property_constructed("_sf_mw_phase")
                for j, p in enumerate(model.state.phase_list):
                    if construct_density:
                        assert get_scaling_factor(
                            model.state.dens_mol_phase[p]
                        ) == pytest.approx(2 / ((j + 5) * (j + 6)), rel=1e-8)
                    else:
                        assert not model.state.is_property_constructed("dens_mol_phase")
                    if construct_volume:
                        assert (
                            get_scaling_factor(model.state.vol_mol_phase[p])
                            == (j + 5) * (j + 6) / 2
                        )
                    else:
                        assert not model.state.is_property_constructed("vol_mol_phase")
            else:
                assert model.state._sf_mw_phase is self._sf_mw_phase
                if not construct_density:
                    assert not model.state.is_property_constructed("dens_mol_phase")
                if not construct_volume:
                    assert not model.state.is_property_constructed("vol_mol_phase")

        else:
            if default_density and default_volume:
                assert not model.state.is_property_constructed("_sf_mw_phase")
                for j, p in enumerate(model.state.phase_list):
                    if construct_density:
                        assert get_scaling_factor(model.state.dens_mol_phase[p]) == 137
                    else:
                        assert not model.state.is_property_constructed("dens_mol_phase")
                    if construct_volume:
                        assert get_scaling_factor(model.state.vol_mol_phase[p]) == 314
                    else:
                        assert not model.state.is_property_constructed("vol_mol_phase")
            elif default_density:
                assert not model.state.is_property_constructed("_sf_mw_phase")
                for j, p in enumerate(model.state.phase_list):
                    if construct_density:
                        assert get_scaling_factor(model.state.dens_mol_phase[p]) == 137
                    else:
                        assert not model.state.is_property_constructed("dens_mol_phase")
                    if construct_volume:
                        assert get_scaling_factor(
                            model.state.vol_mol_phase[p]
                        ) == pytest.approx(1 / 137, rel=1e-8)
                    else:
                        assert not model.state.is_property_constructed("vol_mol_phase")
            elif default_volume:
                assert not model.state.is_property_constructed("_sf_mw_phase")
                for j, p in enumerate(model.state.phase_list):
                    if construct_density:
                        assert get_scaling_factor(
                            model.state.dens_mol_phase[p]
                        ) == pytest.approx(1 / 314, rel=1e-8)
                    else:
                        assert not model.state.is_property_constructed("dens_mol_phase")
                    if construct_volume:
                        assert get_scaling_factor(model.state.vol_mol_phase[p]) == 314
                    else:
                        assert not model.state.is_property_constructed("vol_mol_phase")
            else:
                assert model.state._sf_mw_phase is self._sf_mw_phase
                if not construct_density:
                    assert not model.state.is_property_constructed("dens_mol_phase")
                if not construct_volume:
                    assert not model.state.is_property_constructed("vol_mol_phase")
