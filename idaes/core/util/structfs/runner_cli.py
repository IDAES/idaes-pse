#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Command-line interface to run a module.
"""
import os
import argparse
import importlib
import importlib.util
import json
import logging
import traceback
from io import FileIO
import sys

from idaes.core.util.structfs.runner import Runner

_log = logging.getLogger("idaes.core.util.structfs.runner_cli")


def _error(ofile: FileIO, msg: str, code: int = -1) -> int:
    stack_trace = traceback.format_exc()
    d = {"status": code, "error": msg, "error_detail": stack_trace}
    json.dump(d, ofile)
    _log.error(msg)
    return code


def main():
    """Program entry point."""
    p = argparse.ArgumentParser()
    p.add_argument("module", help="Python module name or path to .py file")
    p.add_argument("output", help="Output file for result JSON")
    p.add_argument("--object", help="Object in module (default=FS)", default="FS")
    p.add_argument(
        "--to", help="Step name to run to (default=all)", default=Runner.STEP_ANY
    )
    p.add_argument(
        "--info",
        action="store_true",
        help="Instead of running module, get information like list of steps",
        default=False,
    )
    args = p.parse_args()

    try:
        ofile = open(args.output, mode="w")
    except IOError as err:
        return _error(
            sys.stderr, f"Cannot open output fNoneile '{args.output}': {err}", -1
        )

    module_name = args.module
    try:
        mod = _load_module(module_name)
    except ValueError as err:
        return _error(ofile, str(err), 1)
    except ModuleNotFoundError:
        return _error(ofile, f"Could not find module '{module_name}'", 2)

    obj_name = args.object
    try:
        obj = getattr(mod, obj_name)
        if not isinstance(obj, Runner):
            return _error(
                ofile,
                f"Object must be an instance of the Runner class, got '{obj.__class__.__name__}'",
                3,
            )
    except AttributeError:
        return _error(
            ofile, f"Could not find object '{obj_name}' in module '{module_name}'", 4
        )

    if args.info:
        report = {"steps": obj.list_steps(), "class_name": obj.__class__.__name__}
    else:
        to_step = args.to
        try:
            obj.run_steps(first=Runner.STEP_ANY, last=to_step)
        except Exception as e:  # pylint: disable=broad-exception-caught
            return _error(ofile, f"While running steps: {e}", 5)

        report = obj.report()
        report["status"] = 0

    json.dump(report, ofile)

    return 0


def _load_module(module_or_path: str):
    """
    Load a module - supports both module names and file paths.

    Args:
        module_or_path: Can be either:
            - Module name: "idaes.models.flash_flowsheet"
            - File path: "/Users/user/Downloads/my_flowsheet.py"
    Returns:
        module: The loaded Python module object.

    Note:
        For file paths, this function sets up a pseudo-package structure to
        support relative imports (e.g., 'from ..sibling import something').
    """
    # Check if input is a file path
    if module_or_path.endswith(".py") or os.path.isfile(module_or_path):
        # This is a file path
        file_path = os.path.abspath(module_or_path)

        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File not found: {file_path}")

        # Get directory structure for package simulation
        dir_path = os.path.dirname(file_path)  # e.g., /Users/user/workspace/subdir
        parent_dir = os.path.dirname(dir_path)  # e.g., /Users/user/workspace
        package_name = os.path.basename(dir_path)  # e.g., "subdir"
        module_basename = os.path.splitext(os.path.basename(file_path))[
            0
        ]  # e.g., "test"
        full_module_name = f"{package_name}.{module_basename}"  # e.g., "subdir.test"

        # Add parent directory to sys.path so Python can find sibling packages
        # This enables imports like "from ..fsrunner import ..."
        if parent_dir not in sys.path:
            sys.path.insert(0, parent_dir)

        # Create module spec with submodule_search_locations for package support
        spec = importlib.util.spec_from_file_location(
            full_module_name, file_path, submodule_search_locations=[dir_path]
        )

        if spec is None or spec.loader is None:
            raise ImportError(f"Cannot create module spec for {file_path}")

        # Create the module object from spec
        module = importlib.util.module_from_spec(spec)

        # KEY: Set __package__ so relative imports know the package context
        module.__package__ = package_name

        # Register in sys.modules so other imports can find it
        sys.modules[full_module_name] = module

        # Execute the module code (this actually loads the content)
        spec.loader.exec_module(module)
        return module
    else:
        if module_or_path.startswith("."):
            raise ValueError("Relative module names not allowed")
        # This is a module name, use the original logic
        return importlib.import_module(module_or_path)


if __name__ == "__main__":
    sys.exit(main())
