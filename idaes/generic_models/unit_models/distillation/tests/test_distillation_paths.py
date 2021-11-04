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
Tests for deprecation paths for renamed distillation models
"""
import pytest


@pytest.mark.unit
def test_import_condenser_from_init(caplog):
    from idaes.generic_models.unit_models.distillation import Condenser
    assert (
        "DEPRECATED: the 'Condenser' class has been moved to"
        "'idaes.generic_models.unit_models.column_models.condenser.Condenser'"
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_old_condenser(caplog):
    from idaes.generic_models.unit_models.distillation.condenser import Condenser
    assert (
        "DEPRECATED: the 'Condenser' class has been moved to"
        "'idaes.generic_models.unit_models.column_models.condenser.Condenser'"
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_old_condenser_type(caplog):
    from idaes.generic_models.unit_models.distillation.condenser import \
        CondenserType
    assert (
        "DEPRECATED: the 'CondenserType' class has been moved to "
        "'idaes.generic_models.unit_models.column_models."
        "condenser.CondenserType'"
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_old_temperature_spec(caplog):
    from idaes.generic_models.unit_models.distillation.condenser import \
        TemperatureSpec
    assert (
        "DEPRECATED: the 'TemperatureSpec' class has been moved to "
        "'idaes.generic_models.unit_models.column_models."
        "condenser.TemperatureSpec'."
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_reboiler_from_init(caplog):
    from idaes.generic_models.unit_models.distillation import Reboiler
    assert (
        "DEPRECATED: the 'Reboiler' class has been moved to"
        "'idaes.generic_models.unit_models.column_models.reboiler.Reboiler'"
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_old_reboiler(caplog):
    from idaes.generic_models.unit_models.distillation.reboiler import Reboiler
    assert (
        "DEPRECATED: the 'Reboiler' class has been moved to"
        "'idaes.generic_models.unit_models.column_models.reboiler.Reboiler'"
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_tray_from_init(caplog):
    from idaes.generic_models.unit_models.distillation import Tray
    assert (
        "DEPRECATED: the 'Tray' class has been moved to"
        "'idaes.generic_models.unit_models.column_models.tray.Tray'"
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_old_tray(caplog):
    from idaes.generic_models.unit_models.distillation.tray import Tray
    assert (
        "DEPRECATED: the 'Tray' class has been moved to"
        "'idaes.generic_models.unit_models.column_models.tray.Tray'"
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_tray_column_from_init(caplog):
    from idaes.generic_models.unit_models.distillation import TrayColumn
    assert (
        "DEPRECATED: the 'TrayColumn' class has been moved to "
        "'idaes.generic_models.unit_models.column_models."
        "tray_column.TrayColumn'") in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_old_tray_column(caplog):
    from idaes.generic_models.unit_models.distillation.tray_column import \
        TrayColumn
    assert (
        "DEPRECATED: the 'TrayColumn' class has been moved to "
        "'idaes.generic_models.unit_models.column_models."
        "tray_column.TrayColumn'") in caplog.text.replace('\n', '')
