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
"""Commandline Utilities for Managing the IDAES Config files"""

__author__ = "John Eslick"

import idaes
import click
import json
import pyomo.common.config
from idaes.commands import cb
import idaes.config as cfg


@cb.command(name="config-write", help="Write the IDAES config")
@click.option("--file", default="idaes_config.json", help="File name to write.")
@click.option("--default", is_flag=True, help="Write the default config")
def config_write(file, default):
    idaes.write_config(file, default)
    click.echo(f"Wrote config to {file}")


@cb.command(name="config-file", help="Show the config file paths")
@click.option("--global", "gbl", is_flag=True, help="Global config path")
@click.option("--local", is_flag=True, help="Local config path")
def config_file(gbl, local):
    if gbl:
        click.echo(idaes._global_config_file)
    elif local:
        click.echo(idaes._local_config_file)
    else:
        click.echo(f"global: {idaes._global_config_file}")
        click.echo(f"local: {idaes._local_config_file}")


@cb.command(name="config-set", help="Set a set a configuration option")
@click.argument("key")
@click.argument("value")
@click.option("--global", "glb", is_flag=True, help="Write to global config")
@click.option("--local", is_flag=True, help="Write to local config")
@click.option("--file", help="Write to specified file", type=str)
@click.option("--display", is_flag=True, help="Show new config")
@click.option("--add", is_flag=True, help="Add to list or set")
@click.option("--del", "dlt", is_flag=True, help="Delete from to list or set")
@click.option(
    "--file_as_global",
    is_flag=True,
    help="Testing option, alternate global config location",
)
@click.option(
    "--file_as_local",
    is_flag=True,
    help="Testing option, alternate local config file location",
)
def config_set(
    glb, local, key, value, file, file_as_local, file_as_global, display, add, dlt
):
    # get locations for the local and global config files
    global_config_file = idaes._global_config_file
    local_config_file = idaes._local_config_file
    # for testing allow an alternate location to be used for local and global
    # config files, to avoid messing with the real config files
    if file_as_local:
        local_config_file = file
        file = None
    if file_as_global:
        global_config_file = file
        file = None
    # make sure one and only one place is specified to write to
    if not (glb ^ local ^ bool(file)):
        click.echo("Must specify exactly one of --global, --local, --file")
        return
    # if global make sure you don't pick up local config and write it back
    if glb:
        idaes.cfg = cfg._new_idaes_config_block()
        idaes.read_config(global_config_file)
        idaes.cfg.display()
    elif local:
        idaes.cfg = cfg._new_idaes_config_block()
        idaes.read_config(global_config_file)
        idaes.read_config(local_config_file)
    idaes.cfg.display()
    # get the config block entry
    key = key.split(":")  # use : instead of . due to logger names
    value = value.replace("'", '"')  # " expected by json but gets eaten by shell
    value = json.loads(value)
    c = idaes.cfg
    for k in key[:-1]:
        try:
            c = c[k]
        except KeyError:
            if isinstance(c, pyomo.common.config.ConfigBlock):
                c[k] = pyomo.common.config.ConfigBlock(implicit=True)
            else:
                c[k] = {}
            c = c[k]
    # Set the config value
    if add and isinstance(c[key[-1]], list):
        c[key[-1]].append(value)
    elif add and isinstance(c[key[-1]], set):
        c[key[-1]].add(value)
    elif dlt and isinstance(c[key[-1]], (set, list)):
        c[key[-1]].remove(value)
    elif dlt or add:
        click.echo("--add or --del options must be used on a list or set")
    else:
        c[key[-1]] = value
    # Write changes
    if glb:
        idaes.write_config(global_config_file)
    elif local:
        idaes.write_config(local_config_file)
    elif file:
        idaes.write_config(file)
    if display:
        idaes.cfg.display()


@cb.command(name="config-display", help="Show IDAES config")
def config_display():
    idaes.cfg.display()
