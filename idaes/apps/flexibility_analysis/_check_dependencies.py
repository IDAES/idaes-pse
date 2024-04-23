from pyomo.common.dependencies import attempt_import
import unittest

coramin, coramin_available = attempt_import("coramin")
np, nump_available = attempt_import("numpy")
if not coramin_available or not nump_available:
    raise unittest.SkipTest("flexibility_analysis tests require coramin and numpy")
