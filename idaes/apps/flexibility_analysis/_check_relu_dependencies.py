from pyomo.common.dependencies import attempt_import
import unittest

tensorflow, tensorflow_available = attempt_import("tensorflow")
omlt, nump_available = attempt_import("omlt")
if not tensorflow_available or not nump_available:
    raise unittest.SkipTest("flexibility_analysis tests require tensorflow and omlt")
