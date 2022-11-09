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
import pytest

from pyomo.environ import units

from idaes.core.base.property_set import PropertyMetadata, PropertySetBase
from idaes.core.util.exceptions import PropertyPackageError


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
    def test_all_args(self, meta):
        meta = PropertyMetadata(
            name="test",
            units=units.dimensionless,
            method="foo",
            supported=True,
            required=True,
        )
        assert meta.name == "test"
        assert meta.units is units.dimensionless
        assert meta.method == "foo"
        assert meta.supported
        assert meta.required

    @pytest.mark.unit
    def test_default_value(self, meta):
        assert meta.name == "test"
        assert meta.units is units.dimensionless
        assert meta.method is None
        assert not meta.supported
        assert not meta.required

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
    def test_update_property_args(self, meta):
        assert meta.name == "test"
        assert meta.units is units.dimensionless
        assert meta.method is None
        assert not meta.supported
        assert not meta.required

        meta.update_property(method="foo", required=True, supported=False)

        assert meta.name == "test"
        assert meta.units is units.dimensionless
        assert meta.method == "foo"
        assert not meta.supported
        assert meta.required

        meta.update_property(method="foo", required=False, supported=True)

        assert meta.name == "test"
        assert meta.units is units.dimensionless
        assert meta.method == "foo"
        assert meta.supported
        assert not meta.required


class DummyMeta:
    pass


class TestPropertySetBase:
    @pytest.fixture(scope="function")
    def pset(self):
        p = DummyMeta()
        p.default_units = "foo"

        pset = PropertySetBase(parent=p)

        return pset

    @pytest.mark.unit
    def test_parent(self):
        p = object()
        pset = PropertySetBase(parent=p)

        assert pset._parent_block is p

    @pytest.mark.unit
    def test_add_property_meta(self, pset):
        pset._add_property_meta(
            local_name="foo",
            name="bar",
            method="baz",
            supported=True,
            required=True,
            units=units.dimensionless,
        )

        assert isinstance(pset.foo, PropertyMetadata)
        assert pset.foo.name == "bar"
        assert pset.foo.method == "baz"
        assert pset.foo.supported
        assert pset.foo.required
        assert pset.foo.units is units.dimensionless

    @pytest.mark.unit
    def test_add_property_meta_duplicate(self, pset):
        pset._add_property_meta(
            local_name="foo",
            name="bar",
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
            pset._add_property_meta(
                local_name="foo",
                name="bar",
                method="baz",
                supported=True,
                required=True,
                units=units.dimensionless,
            )

    @pytest.mark.unit
    def test_getitem(self, pset):
        pset._add_property_meta(
            local_name="FOO",
            name="bar",
            method="baz",
            supported=True,
            required=True,
            units=units.dimensionless,
        )

        # Test both capitalized and uncapitalized
        assert pset["foo"] is pset.FOO
        assert pset["FOO"] is pset.FOO

    @pytest.mark.unit
    def test_getitem_no_present(self, pset):
        with pytest.raises(
            KeyError, match="Property foo is not defined in this PropertySet."
        ):
            pset["foo"]

    @pytest.mark.unit
    def test_iter(self, pset):
        pset._add_property_meta(
            local_name="FOO",
            name="bar",
            method="baz",
            supported=True,
            required=True,
            units=units.dimensionless,
        )
        pset._add_property_meta(
            local_name="_BAR",
            name="bar",
            method="baz",
            supported=True,
            required=True,
            units=units.dimensionless,
        )

        # Only FOO should be returned by iter
        for i in pset:
            assert i is pset.FOO

    @pytest.mark.unit
    def test_define_property(self, pset):
        pset.define_property(
            name="bar",
            method="baz",
            supported=True,
            required=True,
            units=units.dimensionless,
        )

        assert isinstance(pset.bar, PropertyMetadata)
        assert pset.bar.name == "bar"
        assert pset.bar.method == "baz"
        assert pset.bar.supported
        assert pset.bar.required
        assert pset.bar.units is units.dimensionless

    @pytest.mark.unit
    def test_define_stanadrd_property(self, pset):
        pset._define_standard_property(
            name="bar",
            method="baz",
            supported=True,
            required=True,
            units=units.dimensionless,
        )

        assert isinstance(pset._bar, PropertyMetadata)
        assert pset._bar.name == "bar"
        assert pset._bar.method == "baz"
        assert pset._bar.supported
        assert pset._bar.required
        assert pset._bar.units is units.dimensionless

    @pytest.mark.unit
    def test_list_required_properties(self, pset):
        pset._add_property_meta(
            local_name="FOO",
            name="bar",
            method="baz",
            supported=True,
            required=True,
            units=units.dimensionless,
        )
        pset._add_property_meta(
            local_name="BAR",
            name="bar",
            method="baz",
            supported=True,
            required=False,
            units=units.dimensionless,
        )

        assert pset.list_required_properties() == [pset.FOO]

    @pytest.mark.unit
    def test_check_required_proeprties(self, pset):
        # ADd some required properties
        pset._add_property_meta(
            local_name="FOO",
            name="foo",
            method="baz",
            supported=True,
            required=True,
            units=units.dimensionless,
        )
        pset._add_property_meta(
            local_name="BAR",
            name="bar",
            method="baz",
            supported=True,
            required=True,
            units=units.dimensionless,
        )
        pset._add_property_meta(
            local_name="BAZ",
            name="baz",
            method="baz",
            supported=True,
            required=True,
            units=units.dimensionless,
        )

        # Create a second property set, which only supports some required properties
        pset2 = PropertySetBase(parent=pset._parent_block)
        # In metadata but explicitly not supported
        pset2._add_property_meta(
            local_name="FOO",
            name="foo",
            method="baz",
            supported=False,
            required=True,
            units=units.dimensionless,
        )
        # In metadata and supported
        pset2._add_property_meta(
            local_name="BAR",
            name="bar",
            method="baz",
            supported=True,
            required=True,
            units=units.dimensionless,
        )
        # BAZ is missing from metadata

        # Check required properties should return FOO and BAZ
        unsupported = pset.check_required_properties(pset2)
        assert "foo" in unsupported
        assert "baz" in unsupported
        assert len(unsupported) == 2

    @pytest.mark.unit
    def test_list_supported_properties(self, pset):
        pset._add_property_meta(
            local_name="FOO",
            name="bar",
            method="baz",
            supported=True,
            required=True,
            units=units.dimensionless,
        )
        pset._add_property_meta(
            local_name="BAR",
            name="bar",
            method="baz",
            supported=False,
            required=False,
            units=units.dimensionless,
        )

        assert pset.list_supported_properties() == [pset.FOO]

    @pytest.mark.unit
    def test_list_unsupported_properties(self, pset):
        pset._add_property_meta(
            local_name="FOO",
            name="bar",
            method="baz",
            supported=True,
            required=True,
            units=units.dimensionless,
        )
        pset._add_property_meta(
            local_name="BAR",
            name="bar",
            method="baz",
            supported=False,
            required=False,
            units=units.dimensionless,
        )

        assert pset.list_unsupported_properties() == [pset.BAR]

    @pytest.mark.unit
    def test_unitset(self):
        p = DummyMeta()
        p.default_units = "foo"

        pset = PropertySetBase(parent=p)

        return pset.unitset is p.default_units


# TODO: How to test standard and electrolyte property sets?
