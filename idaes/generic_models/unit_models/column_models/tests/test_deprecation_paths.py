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
Tests for deprecation paths for renamed column_models models
"""
import pytest


@pytest.mark.unit
def test_import_condenser_from_init(caplog):
    from idaes.generic_models.unit_models.column_models import Condenser
    assert (
        "DEPRECATED: the 'Condenser' class has been moved to"
        "'idaes.models_extra.column_models.condenser.Condenser'"
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_old_condenser(caplog):
    from idaes.generic_models.unit_models.column_models.condenser import Condenser
    assert (
        "DEPRECATED: the 'Condenser' class has been moved to"
        "'idaes.models_extra.column_models.condenser.Condenser'"
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_old_condenser_type(caplog):
    from idaes.generic_models.unit_models.column_models.condenser import \
        CondenserType
    assert (
        "DEPRECATED: the 'CondenserType' class has been moved to"
        "'idaes.models_extra.column_models.condenser.CondenserType'"
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_old_temperature_spec(caplog):
    from idaes.generic_models.unit_models.column_models.condenser import \
        TemperatureSpec
    assert (
        "DEPRECATED: the 'TemperatureSpec' class has been moved to"
        "'idaes.models_extra.column_models."
        "condenser.TemperatureSpec'."
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_reboiler_from_init(caplog):
    from idaes.generic_models.unit_models.column_models import Reboiler
    assert (
        "DEPRECATED: the 'Reboiler' class has been moved to"
        "'idaes.models_extra.column_models.reboiler.Reboiler'"
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_old_reboiler(caplog):
    from idaes.generic_models.unit_models.column_models.reboiler import Reboiler
    assert (
        "DEPRECATED: the 'Reboiler' class has been moved to"
        "'idaes.models_extra.column_models.reboiler.Reboiler'"
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_tray_from_init(caplog):
    from idaes.generic_models.unit_models.column_models import Tray
    assert (
        "DEPRECATED: the 'Tray' class has been moved to"
        "'idaes.models_extra.column_models.tray.Tray'"
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_old_tray(caplog):
    from idaes.generic_models.unit_models.column_models.tray import Tray
    assert (
        "DEPRECATED: the 'Tray' class has been moved to"
        "'idaes.models_extra.column_models.tray.Tray'"
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_tray_column_from_init(caplog):
    from idaes.generic_models.unit_models.column_models import TrayColumn
    assert (
        "DEPRECATED: the 'TrayColumn' class has been moved to"
        "'idaes.models_extra.column_models."
        "tray_column.TrayColumn'") in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_old_tray_column(caplog):
    from idaes.generic_models.unit_models.column_models.tray_column import \
        TrayColumn
    assert (
        "DEPRECATED: the 'TrayColumn' class has been moved to"
        "'idaes.models_extra.column_models."
        "tray_column.TrayColumn'") in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_old_solvent_column(caplog):
    from idaes.generic_models.unit_models.column_models.solvent_column import \
        PackedColumn
    assert (
        "DEPRECATED: the 'PackedColumn' class has been moved to"
        "'idaes.models_extra.column_models."
        "solvent_column.PackedColumn'") in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_old_solvent_condenser(caplog):
    from idaes.generic_models.unit_models.column_models.solvent_condenser import \
        SolventCondenser
    assert (
        "DEPRECATED: the 'SolventCondenser' class has been moved to"
        "'idaes.models_extra.column_models."
        "solvent_condenser.SolventCondenser'") in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_old_solvent_reboiler(caplog):
    from idaes.generic_models.unit_models.column_models.solvent_reboiler import \
        SolventReboiler
    assert (
        "DEPRECATED: the 'SolventReboiler' class has been moved to"
        "'idaes.models_extra.column_models."
        "solvent_reboiler.SolventReboiler'") in caplog.text.replace('\n', '')
