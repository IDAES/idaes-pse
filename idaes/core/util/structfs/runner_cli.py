"""
Command-line interface to run a module.
"""

import argparse
import importlib
import json
import logging
import traceback
from io import FileIO
import sys

from idaes.core.util.structfs.runner import Runner

_log = logging.getLogger("idaes.core.util.structfs.runner_cli")


def error(ofile: FileIO, msg: str, code: int = -1) -> int:
    stack_trace = traceback.format_exc()
    d = {"status": code, "error": msg, "error_detail": stack_trace}
    json.dump(d, ofile)
    _log.error(msg)
    return code


def main():
    p = argparse.ArgumentParser()
    p.add_argument("module", help="Python module name", default=None)
    p.add_argument(
        "--output",
        default="result.json",
        help="Output file for result JSON (default=result.json)",
    )
    p.add_argument("--object", help="Object in module (default=FS)", default="FS")
    p.add_argument(
        "--to", help="Step name to run to (default=all)", default=Runner.STEP_ANY
    )
    args = p.parse_args()

    try:
        ofile = open(args.output, mode="w")
    except Exception as err:
        return error(sys.stderr, "Cannot open output file: {err}", -1)

    if args.module is None:
        return error(ofile, f"Argument --module is required")
    if args.to is None:
        return error(ofile, f"Argument --to is required")

    module_name = args.module
    if module_name.startswith("."):
        return error(ofile, "Relative module names not allowed", 1)
    try:
        mod = importlib.import_module(module_name)
    except ModuleNotFoundError:
        return error(ofile, f"Could not find module '{module_name}'", 2)

    obj_name = args.object
    try:
        obj = getattr(mod, obj_name)
        if not isinstance(obj, Runner):
            return error(
                ofile,
                f"Object must be an instance of the Runner class, got '{obj.__class__.__name__}'",
                3,
            )
    except AttributeError:
        return error(
            ofile, f"Could not find object '{obj_name}' in module '{module_name}'", 4
        )

    to_step = args.to
    try:
        obj.run_steps(first=Runner.STEP_ANY, last=to_step)
    except Exception as e:
        return error(ofile, f"While running steps: {e}", 5)

    report = obj.report()
    report["status"] = 0
    json.dump(report, ofile)
    print(f"===> {ofile.name}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
