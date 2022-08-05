#!/usr/bin/env python
import argparse
import gc
import logging
import os
import platform
import sys
import time

try:
    import ujson as json
except ImportError:
    import json

from collections import OrderedDict

import pyomo.common.unittest as unittest
import pyomo
from pyomo.version import version_info as pyomo_version
from pyomo.common.timing import TicTocTimer


class TimingHandler(logging.Handler):
    def __init__(self):
        super(TimingHandler, self).__init__()
        self._testRecord = None
        self.enabled = True

    def setTest(self, testRecord):
        self._testRecord = testRecord["timing"] = OrderedDict()

    def clearTest(self):
        self._testRecord = None

    def emit(self, record):
        if self._testRecord is None:
            return
        cat = record.msg.__class__.__name__
        if cat in self._testRecord:
            _cat = self._testRecord[cat]
        else:
            _cat = self._testRecord[cat] = OrderedDict()
        if isinstance(record.msg, str):
            name = record.msg
        else:
            try:
                name = record.msg.obj.name
            except AttributeError:
                name = record.msg.obj.__class__.__name__
        if name in _cat:
            _val = _cat[name]
            if type(_val) is not list:
                _val = [_val]
            _val.append(record.msg.timer)
        else:
            _val = record.msg.timer
        _cat[name] = _val


class DataRecorder(object):
    def __init__(self, data):
        self._data = data
        self._timer = TicTocTimer()
        self._category = {}

        self._timingHandler = TimingHandler()
        timing_logger = logging.getLogger("pyomo.common.timing")
        timing_logger.setLevel(logging.INFO)
        timing_logger.addHandler(self._timingHandler)

    @unittest.pytest.fixture(autouse=True)
    def add_data_attribute(self, request):
        # set a class attribute on the invoking test context
        request.cls.testdata = OrderedDict()
        self._data[request.node.nodeid] = request.cls.testdata
        yield
        request.cls.testdata = None

    @unittest.pytest.hookimpl(hookwrapper=True)
    def pytest_runtest_call(self, item):
        self._timingHandler.setTest(self._data[item.nodeid])
        # Trigger garbage collection (try and get a "clean" environment)
        gc.collect()
        gc.collect()
        self._timer.tic("")
        # Run the test
        yield
        # Collect the timing and clean up
        self._data[item.nodeid]["test_time"] = self._timer.toc("")
        self._timingHandler.clearTest()


def getPyomoInfo():
    cwd = os.getcwd()
    sha = None
    diffs = []
    branch = None
    try:
        os.chdir(os.path.dirname(pyomo.__file__))
        sha = os.popen("git rev-parse HEAD").read().strip()
        diffs = os.popen("git diff-index --name-only HEAD").read()
        diffs = diffs.strip().split()
        branch = os.popen("git symbolic-ref -q --short HEAD").read().strip()
    finally:
        os.chdir(cwd)
    return {
        "branch": branch,
        "sha": sha,
        "diffs": diffs,
        "pyomo_version": pyomo_version,
    }


def getRunInfo(cython):
    info = {
        "time": time.time(),
        "python_implementation": platform.python_implementation(),
        "python_version": tuple(sys.version_info),
        "python_build": platform.python_build(),
        "platform": platform.system(),
        "hostname": platform.node(),
    }
    if info["python_implementation"].lower() == "pypy":
        info["pypy_version"] = tuple(sys.pypy_version_info)
    if cython:
        import Cython

        info["cython"] = tuple(int(x) for x in Cython.__version__.split("."))
    info.update(getPyomoInfo())
    return info


def run_tests(cython, argv):
    gc.collect()
    gc.collect()
    results = (getRunInfo(cython), OrderedDict())
    recorder = DataRecorder(results[1])
    unittest.pytest.main(argv, plugins=[recorder])
    gc.collect()
    gc.collect()
    return results


def main(argv):
    parser = argparse.ArgumentParser(
        epilog="Remaining arguments are passed to nosetests"
    )
    parser.add_argument(
        "-o",
        "--output",
        action="store",
        dest="output",
        default=None,
        help="Store the test results to the specified file.",
    )
    parser.add_argument(
        "-d",
        "--dir",
        action="store",
        dest="output_dir",
        default=None,
        help="Store the test results in the specified directory.  If -o "
        "is not specified, then a file name is automatically generated "
        "based on the git branch and hash.",
    )
    parser.add_argument(
        "-n",
        "--replicates",
        action="store",
        dest="replicates",
        type=int,
        default=1,
        help="Number of replicates to run.",
    )
    parser.add_argument(
        "--with-cython",
        action="store_true",
        dest="cython",
        help="Cythonization enabled.",
    )

    options, argv = parser.parse_known_args(argv)
    argv.append("-W ignore::Warning")
    cython = options.cython
    results = tuple(run_tests(cython, argv) for i in range(options.replicates))
    results = (results[0][0],) + tuple(r[1] for r in results)

    if options.output_dir:
        if not options.output:
            options.output = "perf-%s-%s-%s-%s.json" % (
                results[0]["branch"],
                results[0]["sha"][:7] + ("_mod" if results[0]["diffs"] else ""),
                results[0]["python_implementation"].lower()
                + (".".join(str(i) for i in results[0]["python_version"][:3])),
                time.strftime("%y%m%d_%H%M", time.localtime()),
            )
        options.output = os.path.join(options.output_dir, options.output)
    if options.output:
        print(f"Writing results to {options.output}")
        ostream = open(options.output, "w")
        close_ostream = True
    else:
        ostream = sys.stdout
        close_ostream = False
    try:
        # Note: explicitly specify sort_keys=False so that ujson
        # preserves the OrderedDict keys in the JSON
        json.dump(results, ostream, indent=2, sort_keys=False)
    finally:
        if close_ostream:
            ostream.close()
    print("Performance run complete.")


if __name__ == "__main__":
    main(sys.argv)
