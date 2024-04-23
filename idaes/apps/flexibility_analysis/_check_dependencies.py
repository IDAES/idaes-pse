from pyomo.common.dependencies import attempt_import
import unittest

np, nump_available = attempt_import("numpy")
if not nump_available:
    raise unittest.SkipTest("flexibility_analysis tests require numpy")
