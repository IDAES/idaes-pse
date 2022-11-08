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

from idaes.core.base.property_set import PropertyMetadata


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
    def test_no_units(self):
        with pytest.raises(TypeError, match='"units" is required'):
            meta = PropertyMetadata(
                name="test",
                units=None,
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

    @pytest.mark.unit
    def test_update_property_dict(self, meta):
        assert meta.name == "test"
        assert meta.units is units.dimensionless
        assert meta.method is None
        assert not meta.supported
        assert not meta.required

        meta.update_property({"method": "foo", "required": True, "supported": False})

        assert meta.name == "test"
        assert meta.units is units.dimensionless
        assert meta.method == "foo"
        assert not meta.supported
        assert meta.required

        meta.update_property({"method": "foo", "required": False, "supported": True})

        assert meta.name == "test"
        assert meta.units is units.dimensionless
        assert meta.method == "foo"
        assert meta.supported
        assert not meta.required

    @pytest.mark.unit
    def test_update_property_dict_and_other(self, meta):
        with pytest.raises(
            ValueError,
            match="If dict is provided, cannot provide values for method, required and supported",
        ):
            meta.update_property(
                dict={"method": "foo", "required": True, "supported": False},
                method="bar",
            )

        with pytest.raises(
            ValueError,
            match="If dict is provided, cannot provide values for method, required and supported",
        ):
            meta.update_property(
                dict={"method": "foo", "required": True, "supported": False},
                required="bar",
            )

        with pytest.raises(
            ValueError,
            match="If dict is provided, cannot provide values for method, required and supported",
        ):
            meta.update_property(
                dict={"method": "foo", "required": True, "supported": False},
                supported="bar",
            )
