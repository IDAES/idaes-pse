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
"""Commandline Utilities for Managing the IDAES Config files"""
# TODO: Missing docstrings
# pylint: disable=missing-function-docstring

# TODO: protected access issues
# pylint: disable=protected-access

__author__ = "John Eslick"

import click

from idaes.commands import cb


@cb.command(
    name="config-write",
    help="This command has been removed.",
    context_settings=dict(ignore_unknown_options=True, allow_extra_args=True),
)
def config_write():
    click.echo("This command has been removed.")


@cb.command(
    name="config-file",
    help="This command has been removed.",
    context_settings=dict(ignore_unknown_options=True, allow_extra_args=True),
)
def config_file():
    click.echo("This command has been removed.")


@cb.command(
    name="config-set",
    help="This command has been removed.",
    context_settings=dict(ignore_unknown_options=True, allow_extra_args=True),
)
def config_set():
    click.echo("This command has been removed.")


@cb.command(
    name="config-display",
    help="This command has been removed.",
    context_settings=dict(ignore_unknown_options=True, allow_extra_args=True),
)
def config_display():
    click.echo("This command has been removed.")
