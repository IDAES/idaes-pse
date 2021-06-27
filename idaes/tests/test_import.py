# coding: utf-8
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
Make sure all the files behave well when imported

test_import:
    * imports without errors
    * imports quickly (what code is at top level doesn't take long)

test_toplevel_statements:
    * minimal code at top level
"""
# stdlib
import os
from pathlib import Path
import re
import subprocess
import sys
import time
from typing import List, Union, Dict

# third-party
import pytest

# this package
import idaes


__author__ = "Dan Gunter (LBNL)"


# The root of the source code tree
idaes_path = Path(idaes.__file__).parent

# Constants controlling the expected behavior
IMPORT_TIME_MAX = 5     # Seconds to import a file
TOPLEVEL_LINES_MAX = 5  # Lines of top-level `code` allowed
OUTPUT_OK = True        # Is extra output from the module ok?
BATCH_SIZE = 10         # Size of batch for speeding up imports

####################################################################################


def modules_batched(size=BATCH_SIZE):
    """Get list of modules in batches of provided size.

    This needs to precede `test_import`, since it is called by its pytest
    'parametrize' decorator.
    """
    all_modules = sorted(
        [
            str(p.relative_to(idaes_path)).replace(os.sep, ".")[:-3]
            for p in idaes_path.rglob("*.py")
            if not p.name == "__init__.py"
        ]
    )
    batched, num_modules = [], len(all_modules)
    for i in range(0, num_modules, size):
        batched.append(all_modules[i : min(i + size, num_modules)])
    return batched


@pytest.mark.integration
@pytest.mark.parametrize("module_names", modules_batched())
def test_import(module_names):
    """Run the Python interpreter to import batches of modules.

    Since most of the overhead of running this is in starting the interpreter,
    it is much faster to do, e.g., 10 imports at a time -- as long as failure is
    rare. The drawback is a little extra complexity if there is an error in an
    import, since one needs to go back through the batch and see which module
    actually caused the failure. Note that if failure is *not* rare, we have bigger
    problems than whether this test is slow.
    """
    # Try to import everything in the batch, hope for success
    batch_max = IMPORT_TIME_MAX * len(module_names) + 1
    t0 = time.time()
    results = collect_import_time_batch(module_names, batch_max, allow_output=OUTPUT_OK)
    t1 = time.time()
    # If it succeeded, check that all were under IMPORT_TIME_MAX
    if results:
        for module_name, result in results.items():
            if result.timing > IMPORT_TIME_MAX:
                pytest.fail(
                    f"Import too slow for {module_name}: "
                    f"{result.timing}s > {IMPORT_TIME_MAX}s"
                )
        # success
        print(f"Time for success: {t1 - t0}s = {(t1 - t0) / BATCH_SIZE}")
    # If there was a more serious error, look at the batch one-by-one to
    # figure out which one caused it
    else:
        erred = {}
        for module_name in module_names:
            result = collect_import_time(
                module_name, IMPORT_TIME_MAX, allow_output=OUTPUT_OK
            )
            if result.error is not None:
                erred[module_name] = result
        if erred:  # should always be true
            errlist = [
                f"Import error for {module_name}:\n{result.error}"
                for module_name, result in erred.items()
            ]
            pytest.fail("\n".join(errlist))


class Result:
    def __init__(self):
        self.timing, self.error = -1.0, None

    def add_error(self, e):
        self.error = e if self.error is None else self.error + "\n" + e


def collect_import_time_batch(
    modules: List[str], max_time: float, allow_output: bool = True
) -> Union[Dict, None]:
    """Import a batch of modules by running with 'importtime' option.

    Args:
        modules: List of absolute name of modules to import (idaes.<something>)
        max_time: Timeout for process that imports the modules
        allow_output: If False, consider output during import an error; otherwise
            ignore it.

    Returns:
        Either a map of result objects, if all succeeded, or None
    """
    results, proc = {m: Result() for m in modules}, None

    # Run import in new interpreter with '-X importtime' option
    import_stmt = ";".join([f"import idaes.{m}" for m in modules])
    print(f"import batch: '{import_stmt}'")
    try:
        proc = subprocess.run(
            [sys.executable, "-X", "importtime", "-c", import_stmt],
            stderr=subprocess.STDOUT,
            stdout=subprocess.PIPE,
            timeout=max_time,
        )
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
        return None

    if proc.returncode != 0:
        return None

    # Parse output of '-X importtime'
    output = proc.stdout.decode()
    for line in output.splitlines():
        # For rows with timings
        if line.startswith("import time:"):
            # Save if match with a module in the batch
            for module in modules:
                if line.endswith(module):
                    times = re.findall(r"\d+", line)
                    results[module].timing = float(times[1]) / 1e6
                    break
        elif not allow_output:
            return None

    # Return timing and/or error information
    return results


def collect_import_time(
    module: str, max_time: float, allow_output: bool = True
) -> Result:
    """Import a module by running with 'importtime' option.

    Args:
        modules: List of absolute name of modules to import (idaes.<something>)
        max_time: Timeout for process that imports the module
        allow_output: If False, consider output during import an error; otherwise
        ignore it.

    Returns:
        Result object with an error and, if no error, a timing
    """
    result, proc = Result(), None

    # Run import in new interpreter with '-X importtime' option
    import_stmt = f"import idaes.{module}"
    print(f"import single: '{import_stmt}'")
    try:
        proc = subprocess.run(
            [sys.executable, "-X", "importtime", "-c", import_stmt],
            stderr=subprocess.STDOUT,
            stdout=subprocess.PIPE,
            timeout=max_time,
        )
    except subprocess.CalledProcessError as err:  # total failure!
        result.error = str(err)
    except subprocess.TimeoutExpired:  # timeout
        result.error = f"Timeout (> {max_time}s)"

    if proc is not None:
        output = proc.stdout.decode()
        if proc.returncode != 0:
            # find "traceback" text in output
            result.error = "Non-zero return code:\n"
            for line in output.splitlines():
                if not line.startswith("import time:"):
                    result.add_error(line)
        else:
            # Parse output of '-X importtime'
            prev_import = None
            for line in output.splitlines():
                # For rows with timings
                if line.startswith("import time:"):
                    # Save if match with current module
                    if line.endswith(module):
                        times = re.findall(r"\d+", line)
                        result.timing = float(times[1]) / 1e6
                    # In any case, remember module name
                    prev_import = line[line.rfind("|") + 1:].strip()
                elif not allow_output:
                    # Associate dis-allowed output with last seen module name
                    result.add_error(f"[{prev_import}] {line}")

    # Return timing and/or error information
    return result


####################################################################################


@pytest.mark.unit
@pytest.mark.parametrize("module_path", list(idaes_path.rglob("*.py")))
def test_toplevel_statements(module_path):
    """Look for "too many" top-level statements.

    This test is far quicker than `test_import`, but they test different things.
    It does not have any way of knowing how long it takes to import the module
    (or whether it can be imported).

    It also has a pretty crude idea of what a top-level statement is (and isn't)
    """
    # Skip some paths
    if module_path == idaes_path / "__init__.py":
        return

    with module_path.open(encoding="utf8") as f:
        num_statements = count_statements(f)
        assert num_statements <= TOPLEVEL_LINES_MAX


def count_statements(f):
    statements = 0
    comment_block = False
    lineno = 0

    # Look at file line-by-line
    for line in f:
        lineno += 1
        is_statement = False
        # if in a comment block, wait for it to end
        if comment_block:
            if re.match(r".*('''|\"\"\")", line):
                comment_block = False
        # look for start of comment block
        elif re.match(r"\s*('''|\"\"\")", line):
            comment_block = True
        # skip empy, indented, comment lines
        elif len(line) == 0 or re.match(r"^\s+$", line) or re.match(r"\s|#", line):
            pass
        # skip function, class defs or decorators
        elif line.startswith("def") or line.startswith("class") or line.startswith("@"):
            pass
        # match anything that looks like a function call, except some special ones
        elif re.match(r".*\w+\(", line):
            is_statement = True
            # these things *are* allowed at the top-level
            if re.match(
                ".*(add_command|CONFIG|list\(|dict\(|attempt_import|\.declare\()", line
            ):
                is_statement = False
        # If it looks like a statement, increment the count
        if is_statement:
            statements += 1
            # This makes it easy to see what matched if the file fails the test
            print(f"{lineno:4d}: {line}", end="")

    return statements
