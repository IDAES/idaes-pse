"""
Run all build steps for IDAES documentation.
"""
# standard libary
import argparse
import glob
import logging
import os
from pathlib import Path
import shutil
import subprocess
import sys

__author__ = "Dan Gunter"


_log = logging.getLogger("build_docs")
__h = logging.StreamHandler()
__h.setFormatter(logging.Formatter("%(levelname)s %(asctime)s - %(message)s"))
_log.addHandler(__h)
_log.setLevel(logging.ERROR)


class CommandError(Exception):
    def __init__(self, cmd, why):
        super().__init__(why)
        self.command = cmd
        self.message = why


def pipeline(*commands, **params):
    first = True
    for cmd in commands:
        func = globals()[f"run_{cmd}"]
        print_header(cmd, first=first)
        func(**params)
        first = False


def run_apidoc(clean=True):
    """Run the sphinx-apidoc extension to build API docs.
    """
    if clean:
        # remove old docs
        print_status("Remove old docs")
        for old_rst in glob.glob("apidoc/*.rst"):
            os.unlink(old_rst)
    # run sphinx-apidoc
    print_status("Run sphinx-apidoc")
    os.environ["SPHINX_APIDOC_OPTIONS"] = "members,ignore-module-all,noindex"
    _run(
        "apidoc",
        ["sphinx-apidoc", "--module-first", "--output-dir", "apidoc", "../idaes", "../idaes/*tests*"],
        60,
    )
    postprocess_apidoc(Path("apidoc"))


def postprocess_apidoc(root):
    """Perform postprocessing on generated apidoc files
    """
    # Remove :noindex: from all entries in given modules
    remove_noindex = ["idaes.dmf"]
    for module in remove_noindex:
        module_path = root / (module + ".rst")
        _log.debug(f"Looking for :noindex: in {module_path}")
        lines, has_noindex = [], False
        with module_path.open("r") as input_file:
            for line in input_file:
                if ":noindex:" in line:
                    has_noindex = True
                else:
                    lines.append(line)
        if has_noindex:
            _log.debug(f"Removing :noindex: from {module_path}")
            with module_path.open("w") as output_file:
                for line in lines:
                    output_file.write(line)


def run_html(clean=True):
    """Run sphinx-build to create HTML.
    """
    build_dir = "build"
    output_file = "sphinx-errors.txt"
    if clean:
        if os.path.exists(build_dir):
            print_status(f"Remove build dir '{build_dir}'")
            shutil.rmtree(build_dir)
    if os.path.exists(output_file):
        os.unlink(output_file)
    # run
    print_status("Run sphinx-build")
    _run(
        "html",
        ["sphinx-build", "-M", "html", ".", build_dir, "-w", output_file, "-q"],
        180,
    )
    # check output file
    print_status("Check output file")
    num_lines = num_warnings = num_errors = 0
    if os.path.exists(output_file):
        # Count lines, warnings and errors
        with open(output_file) as outf:
            for line in outf:
                num_lines += 1
                num_warnings = num_warnings+1 if "WARNING: " in line else num_warnings
                num_errors = num_errors+1 if "ERROR: " in line else num_errors
    if num_lines > 0:
        raise CommandError(
            "html",
            f"sphinx-build had {num_warnings} warnings, {num_errors} errors in {num_lines} lines of output\n"
            f"These will cause tests to fail\n"
            f"See file '{output_file}' for details",
        )


def _run(what, args, timeout):
    """Run a command with a timeout.
    """
    _log.debug(f"command='{' '.join(args)}'")
    try:
        proc = subprocess.Popen(args)
    except Exception as err:
        exe = args[0]
        raise CommandError(what, f"Could not run '{exe}':\n{err}")
    try:
        proc.wait(timeout)
    except subprocess.TimeoutExpired:
        raise CommandError(what, f"Timed out after {timeout} seconds")


def print_header(msg: str, first: bool = False):
    if not first:
        print()
    print(f"=== {msg}")


def print_status(msg: str):
    print(f"--- {msg}")


def print_error(cmd, msg):
    lines = msg.split("\n")
    print()
    print(f"*** ERROR in '{cmd}'")
    print("***")
    for line in lines:
        print(f"*** {line}")
    print("***")
    print()


def main() -> int:
    """Entry point for this module.
    """
    prs = argparse.ArgumentParser(description=__doc__.strip())
    prs.add_argument(
        "--dirty",
        action="store_true",
        help="Do not clean files before running commands",
    )
    prs.add_argument(
        "-v",
        "--verbose",
        action="count",
        dest="vb",
        default=0,
        help="Print some debugging information",
    )
    args = prs.parse_args()
    if args.vb > 0:
        _log.setLevel(logging.DEBUG)
    try:
        pipeline("apidoc", "html", clean=not args.dirty)
    except CommandError as err:
        print_error(err.command, err.message)
        return -1
    print_header("SUCCESS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
