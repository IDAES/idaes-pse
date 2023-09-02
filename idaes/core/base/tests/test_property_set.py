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
import pytest

from pyomo.environ import units

from idaes.core.base.property_set import (
    PropertyMetadata,
    PropertySetBase,
    StandardPropertySet,
    ElectrolytePropertySet,
    _PropertyMetadataIndex,
)
from idaes.core.base.property_meta import UnitSet
from idaes.core.util.exceptions import PropertyPackageError


class DummyMeta:
    pass


class Test_PropertyMetadataIndex:
    @pytest.fixture(scope="function")
    def meta(self):
        p = DummyMeta()

        meta = _PropertyMetadataIndex(parent=p, idx="baz")

        return meta

    @pytest.mark.unit
    def test_all_args(self, meta):
        p = DummyMeta()

        meta = _PropertyMetadataIndex(
            parent=p,
            idx="bar",
            method="foo",
            supported=True,
            required=True,
            valid_range=(0, 1),
        )
        assert meta.method == "foo"
        assert meta._parent is p
        assert meta._idx == "bar"
        assert meta.supported
        assert meta.required
        assert meta.valid_range == (0, 1)

    @pytest.mark.unit
    def test_default_value(self, meta):
        assert meta.method is None
        assert not meta.supported
        assert not meta.required

        assert meta._idx == "baz"
        assert meta.valid_range is None

    @pytest.mark.unit
    def test_set_method(self, meta):
        assert meta.method is None
        meta.set_method("foo")
        assert meta.method == "foo"

    @pytest.mark.unit
    def test_set_required(self, meta):
        assert not meta.required
        meta.set_required(True)
        assert meta.required

    @pytest.mark.unit
    def test_set_supported(self, meta):
        assert not meta.supported
        meta.set_supported(True)
        assert meta.supported

    @pytest.mark.unit
    def test_valid_range_verification_invalid(self):
        p = DummyMeta()

        with pytest.raises(
            ValueError, match=r"valid_range must be a tuple of length 2 \(got foo\)"
        ):
            meta = _PropertyMetadataIndex(
                parent=p,
                idx="bar",
                method="foo",
                supported=True,
                required=True,
                valid_range="foo",
            )

    @pytest.mark.unit
    def test_valid_range_verification_length(self):
        p = DummyMeta()

        with pytest.raises(
            ValueError,
            match=r"valid_range must be a tuple of length 2 \(got \(1, 2, 3\)\)",
        ):
            meta = _PropertyMetadataIndex(
                parent=p,
                idx="bar",
                method="foo",
                supported=True,
                required=True,
                valid_range=(1, 2, 3),
            )

    @pytest.mark.unit
    def test_valid_range_verification_values(self):
        p = DummyMeta()

        with pytest.raises(
            ValueError,
            match=r"valid_range must be a 2-tuple with form \(lower, upper\): first value "
            r"was greater than second value: \(1, 0\)",
        ):
            meta = _PropertyMetadataIndex(
                parent=p,
                idx="bar",
                method="foo",
                supported=True,
                required=True,
                valid_range=(1, 0),
            )

    @pytest.mark.unit
    def test_update_property_args(self, meta):
        assert meta.method is None
        assert not meta.supported
        assert not meta.required
        assert meta.valid_range is None

        meta.update_property(method="foo", required=True, supported=False)

        assert meta.method == "foo"
        assert not meta.supported
        assert meta.required
        assert meta.valid_range is None

        meta.update_property(
            method="foo", required=False, supported=True, valid_range=(0, 10)
        )

        assert meta.method == "foo"
        assert meta.supported
        assert not meta.required
        assert meta.valid_range == (0, 10)


class TestPropertyMetadata:
    @pytest.fixture(scope="function")
    def meta(self):
        meta = PropertyMetadata(
            name="test",
            units=units.dimensionless,
        )

        return meta

    @pytest.mark.unit
    def test_no_name(self):
        with pytest.raises(TypeError, match='"name" is required'):
            meta = PropertyMetadata(
                name=None,
                units=units.dimensionless,
            )

    @pytest.mark.unit
    def test_all_args(self):
        meta = PropertyMetadata(
            name="test",
            units=units.dimensionless,
            indices=["a", "b", "none"],
            initialize={
                None: {"method": "foo", "required": True, "supported": True},
                "a": {"method": "bar", "required": False, "supported": True},
                "b": {"method": "baz", "required": True, "supported": False},
            },
        )

        assert meta.name == "test"
        assert meta.units is units.dimensionless
        assert meta._indices == ["a", "b", "none"]

        assert isinstance(meta._none, _PropertyMetadataIndex)
        assert meta._none.method == "foo"
        assert meta._none.required
        assert meta._none.supported

        assert isinstance(meta._a, _PropertyMetadataIndex)
        assert meta._a.method == "bar"
        assert not meta._a.required
        assert meta._a.supported

        assert isinstance(meta._b, _PropertyMetadataIndex)
        assert meta._b.method == "baz"
        assert meta._b.required
        assert not meta._b.supported

    @pytest.mark.unit
    def test_default_value(self, meta):
        assert meta.name == "test"
        assert meta.units is units.dimensionless
        assert meta._indices is None
        assert isinstance(meta._none, _PropertyMetadataIndex)


class TestPropertySetBase:
    @pytest.fixture(scope="function")
    def pset(self):
        p = DummyMeta()
        p.default_units = "foo"

        pset = PropertySetBase(parent=p)

        return pset

    @pytest.mark.unit
    def test_expected_attrs(self):
        p = object()
        pset = PropertySetBase(parent=p)

        assert pset._parent_block is p
        assert pset._defined_properties == []
        assert pset._defined_indices == ["comp", "phase", "phase_comp"]

    @pytest.mark.unit
    def test_direct_assignment(self):
        p = object()
        pset = PropertySetBase(parent=p)

        with pytest.raises(
            TypeError,
            match="PropertySets do not support direct assignment. Please use define_property",
        ):
            pset.foo = "bar"

    @pytest.mark.unit
    def test_add_property(self, pset):
        pset._add_property(
            name="foo",
            method="baz",
            supported=True,
            required=True,
            units=units.dimensionless,
        )

        assert isinstance(pset.foo, PropertyMetadata)
        assert pset.foo.name == "foo"
        assert pset.foo._none.method == "baz"
        assert pset.foo._none.supported
        assert pset.foo._none.required
        assert pset.foo.units is units.dimensionless
        for i in pset.foo._indices:
            assert i in ["comp", "phase", "phase_comp", "none"]

        assert "foo" in pset._defined_properties

        assert isinstance(pset.foo._none, _PropertyMetadataIndex)
        assert isinstance(pset.foo._comp, _PropertyMetadataIndex)
        assert isinstance(pset.foo._phase, _PropertyMetadataIndex)
        assert isinstance(pset.foo._phase_comp, _PropertyMetadataIndex)

    @pytest.mark.unit
    def test_add_property_unindexed_init(self, pset):
        pset._add_property(
            name="foo",
            initialize={
                "none": {
                    "method": "baz",
                    "supported": True,
                    "required": True,
                },
            },
            units=units.dimensionless,
        )

        assert isinstance(pset.foo, PropertyMetadata)
        assert pset.foo.name == "foo"
        assert pset.foo._none.method == "baz"
        assert pset.foo._none.supported
        assert pset.foo._none.required
        assert pset.foo.units is units.dimensionless
        for i in pset.foo._indices:
            assert i in ["comp", "phase", "phase_comp", "none"]

        assert "foo" in pset._defined_properties

        assert isinstance(pset.foo._none, _PropertyMetadataIndex)
        assert isinstance(pset.foo._comp, _PropertyMetadataIndex)
        assert isinstance(pset.foo._phase, _PropertyMetadataIndex)
        assert isinstance(pset.foo._phase_comp, _PropertyMetadataIndex)

    @pytest.mark.unit
    def test_add_property_w_index(self, pset):
        pset._add_property(
            name="foo",
            units=units.dimensionless,
            indices=["a", "b"],
        )

        assert isinstance(pset.foo, PropertyMetadata)
        assert pset.foo.name == "foo"
        assert pset.foo.units is units.dimensionless

        for i in pset.foo._indices:
            assert i in ["a", "b"]

        assert "foo" in pset._defined_properties

        assert isinstance(pset.foo._a, _PropertyMetadataIndex)
        assert isinstance(pset.foo._b, _PropertyMetadataIndex)
        assert not hasattr(pset.foo, "_none")

    @pytest.mark.unit
    def test_add_property_w_index_other(self, pset):
        with pytest.raises(
            ValueError,
            match="Cannot assign method, required or supported attributes directly to indexed "
            "properties. Please use the initialize argument instead.",
        ):
            pset._add_property(
                name="foo",
                method="baz",
                supported=True,
                required=True,
                units=units.dimensionless,
                indices=["a", "b"],
            )

    @pytest.mark.unit
    def test_add_property_init_and_other(self, pset):
        with pytest.raises(
            ValueError,
            match="Cannot provide values for initialize and any of method, required or supported "
            "arguments. Please use the one approach or the other.",
        ):
            pset._add_property(
                name="foo",
                method="baz",
                supported=True,
                required=True,
                units=units.dimensionless,
                initialize="bar",
            )

    @pytest.mark.unit
    def test_define_property(self, pset):
        pset.define_property(
            name="foo",
            method="baz",
            supported=True,
            required=True,
            units=units.dimensionless,
        )

        assert isinstance(pset.foo, PropertyMetadata)
        assert pset.foo.name == "foo"
        assert pset.foo._none.method == "baz"
        assert pset.foo._none.supported
        assert pset.foo._none.required
        assert pset.foo.units is units.dimensionless

        assert isinstance(pset.foo._none, _PropertyMetadataIndex)
        assert isinstance(pset.foo._comp, _PropertyMetadataIndex)
        assert isinstance(pset.foo._phase, _PropertyMetadataIndex)
        assert isinstance(pset.foo._phase_comp, _PropertyMetadataIndex)

    @pytest.mark.unit
    def test_define_property_duplicate(self, pset):
        pset.define_property(
            name="foo",
            method="baz",
            supported=True,
            required=True,
            units=units.dimensionless,
        )

        with pytest.raises(
            PropertyPackageError,
            match="A property with the name foo already exists. "
            "Please use update_property method if you wish to update "
            "an existing property's metadata.",
        ):
            pset.define_property(
                name="foo",
                method="baz",
                supported=True,
                required=True,
                units=units.dimensionless,
            )

    @pytest.mark.unit
    def test_getitem(self, pset):
        pset.define_property(
            name="FOO",
            method="baz",
            supported=True,
            required=True,
            units=units.dimensionless,
        )

        assert pset["FOO"] is pset.FOO._none
        assert pset["FOO_comp"] is pset.FOO._comp
        assert pset["FOO_phase"] is pset.FOO._phase
        assert pset["FOO_phase_comp"] is pset.FOO._phase_comp

    @pytest.mark.unit
    def test_getitem_no_present(self, pset):
        with pytest.raises(
            ValueError,
            match="Unhandled property: foo. This is mostly likely due "
            "to the property not being defined in this PropertySet.",
        ):
            pset["foo"]

    @pytest.mark.unit
    def test_iter(self, pset):
        pset.define_property(
            name="FOO",
            method="baz",
            supported=True,
            required=True,
            units=units.dimensionless,
        )
        pset.define_property(
            name="_BAR",
            method="baz",
            supported=True,
            required=True,
            units=units.dimensionless,
        )

        for i in pset:
            assert i in [pset.FOO, pset._BAR]

    @pytest.mark.unit
    def test_list_required_properties(self, pset):
        pset.define_property(
            name="FOO",
            units=units.dimensionless,
        )
        pset.FOO["none"].set_required(True)
        pset.FOO["phase"].set_required(True)

        rp = pset.list_required_properties()
        for p in rp:
            assert p in [pset.FOO[None], pset.FOO["phase"]]

    @pytest.mark.unit
    def test_check_required_properties(self, pset):
        # ADd some required properties
        pset.define_property(
            name="FOO",
            units=units.dimensionless,
        )
        pset.FOO["none"].set_required(True)
        pset.FOO["phase"].set_required(True)
        pset.FOO["phase_comp"].set_required(True)

        # Create a second property set, which only supports some required properties
        pset2 = PropertySetBase(parent=pset._parent_block)
        # In metadata but explicitly not supported
        pset2.define_property(
            name="FOO",
            units=units.dimensionless,
        )
        pset2.FOO["none"].set_supported(True)
        pset2.FOO["phase"].set_supported(True)
        pset2.FOO["phase_comp"].set_supported(False)

        unsupported = pset.check_required_properties(pset2)
        assert "FOO_phase_comp" in unsupported
        assert len(unsupported) == 1

    @pytest.mark.unit
    def test_list_supported_properties(self, pset):
        pset.define_property(
            name="FOO",
            units=units.dimensionless,
        )
        pset.FOO["none"].set_supported(True)
        pset.FOO["phase"].set_supported(True)

        sp = pset.list_supported_properties()
        for p in sp:
            assert p in [pset.FOO[None], pset.FOO["phase"]]

    @pytest.mark.unit
    def test_list_unsupported_properties(self, pset):
        pset.define_property(
            name="FOO",
            units=units.dimensionless,
        )
        pset.FOO["none"].set_supported(True)
        pset.FOO["phase"].set_supported(True)

        usp = pset.list_unsupported_properties()
        for p in usp:
            assert p in [pset.FOO["comp"], pset.FOO["phase_comp"]]

    @pytest.mark.unit
    def test_unitset(self):
        p = DummyMeta()
        p.default_units = "foo"

        pset = PropertySetBase(parent=p)

        assert pset.unitset is p.default_units

    @pytest.mark.unit
    def test_get_name_and_index(self):
        p = DummyMeta()
        p.default_units = "foo"

        pset = PropertySetBase(parent=p)
        pset._defined_properties.append("foo_bar")

        n, i = pset.get_name_and_index("foo_bar")
        assert n == "foo_bar"
        assert i is None

        n, i = pset.get_name_and_index("foo_bar_phase")
        assert n == "foo_bar"
        assert i == "phase"

        n, i = pset.get_name_and_index("foo_bar_comp")
        assert n == "foo_bar"
        assert i == "comp"

        n, i = pset.get_name_and_index("foo_bar_phase_comp")
        assert n == "foo_bar"
        assert i == "phase_comp"

        with pytest.raises(
            ValueError,
            match="Unhandled property: baz. This is mostly likely due to the "
            "property not being defined in this PropertySet.",
        ):
            pset.get_name_and_index("baz")

    @pytest.mark.unit
    def test_get_name_and_index_phase_frac(self):
        p = DummyMeta()
        p.default_units = "foo"

        pset = PropertySetBase(parent=p)
        pset._defined_properties.append("phase_frac")

        n, i = pset.get_name_and_index("phase_frac")
        assert n == "phase_frac"
        assert i is None

        n, i = pset.get_name_and_index("phase_frac_comp")
        assert n == "phase_frac"
        assert i == "comp"

        n, i = pset.get_name_and_index("phase_frac_phase")
        assert n == "phase_frac"
        assert i == "phase"

        n, i = pset.get_name_and_index("phase_frac_phase_comp")
        assert n == "phase_frac"
        assert i == "phase_comp"

    @pytest.mark.unit
    def test_get_name_and_index_custom(self):
        p = DummyMeta()
        p.default_units = "foo"

        pset = PropertySetBase(parent=p)
        pset._defined_properties.append("foo")
        pset._defined_indices.append("bar")

        n, i = pset.get_name_and_index("foo_bar")
        assert n == "foo"
        assert i == "bar"


class DerivedPropertySet(PropertySetBase):
    foo = PropertyMetadata(name="foo")
    bar = PropertyMetadata(name="bar")
    baz = "baz"


class TestDerivedPropertySet:
    @pytest.mark.unit
    def test_standard_properties(self):
        p = DummyMeta()
        p.default_units = "foo"

        pset = DerivedPropertySet(parent=p)

        for i in pset._defined_properties:
            assert i in ["foo", "bar"]

        assert isinstance(pset.foo, PropertyMetadata)
        assert isinstance(pset.bar, PropertyMetadata)
        assert not isinstance(pset.baz, PropertyMetadata)

        # Check that we didn't accidentally overwrite anything
        assert pset.foo is not DerivedPropertySet.foo
        assert pset.bar is not DerivedPropertySet.bar
        assert pset.baz is DerivedPropertySet.baz


class TestStandardPropertySet:
    @pytest.mark.unit
    def test_defined_indices(self):
        p = DummyMeta()
        p.default_units = UnitSet()

        pset = StandardPropertySet(parent=p)

        assert pset._defined_indices == ["comp", "phase", "phase_comp"]

    @pytest.mark.unit
    def test_unindexed_properties(self):
        p = DummyMeta()
        p.default_units = UnitSet()

        pset = StandardPropertySet(parent=p)

        assert pset.temperature._indices == ["none"]
        assert pset.temperature_bubble._indices == ["none"]
        assert pset.temperature_dew._indices == ["none"]
        assert pset.temperature_red._indices == ["none"]

    @pytest.mark.unit
    def test_get_name_and_index(self):
        p = DummyMeta()
        p.default_units = UnitSet()

        pset = PropertySetBase(parent=p)
        pset._defined_properties.append("foo_bar")

        n, i = pset.get_name_and_index("foo_bar")
        assert n == "foo_bar"
        assert i is None

        n, i = pset.get_name_and_index("foo_bar_phase")
        assert n == "foo_bar"
        assert i == "phase"

        n, i = pset.get_name_and_index("foo_bar_comp")
        assert n == "foo_bar"
        assert i == "comp"

        n, i = pset.get_name_and_index("foo_bar_phase_comp")
        assert n == "foo_bar"
        assert i == "phase_comp"

        with pytest.raises(
            ValueError,
            match="Unhandled property: baz. This is mostly likely due to the "
            "property not being defined in this PropertySet.",
        ):
            pset.get_name_and_index("baz")


class TestElectrolytePropertySet:
    @pytest.mark.unit
    def test_defined_indices(self):
        p = DummyMeta()
        p.default_units = UnitSet()

        pset = ElectrolytePropertySet(parent=p)

        assert pset._defined_indices == [
            "comp",
            "phase",
            "phase_comp",
            "phase_comp_apparent",
            "phase_comp_true",
        ]

    @pytest.mark.unit
    def test_get_name_and_index(self):
        p = DummyMeta()
        p.default_units = UnitSet()

        pset = ElectrolytePropertySet(parent=p)
        pset._defined_properties.append("foo_bar")

        n, i = pset.get_name_and_index("foo_bar")
        assert n == "foo_bar"
        assert i is None

        n, i = pset.get_name_and_index("foo_bar_phase")
        assert n == "foo_bar"
        assert i == "phase"

        n, i = pset.get_name_and_index("foo_bar_comp")
        assert n == "foo_bar"
        assert i == "comp"

        n, i = pset.get_name_and_index("foo_bar_phase_comp")
        assert n == "foo_bar"
        assert i == "phase_comp"

        n, i = pset.get_name_and_index("foo_bar_phase_comp_apparent")
        assert n == "foo_bar"
        assert i == "phase_comp_apparent"

        n, i = pset.get_name_and_index("foo_bar_phase_comp_true")
        assert n == "foo_bar"
        assert i == "phase_comp_true"

        with pytest.raises(
            ValueError,
            match="Unhandled property: baz. This is mostly likely due to the "
            "property not being defined in this PropertySet.",
        ):
            pset.get_name_and_index("baz")
