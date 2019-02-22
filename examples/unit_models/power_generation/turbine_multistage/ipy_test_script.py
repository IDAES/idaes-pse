"""
Ipython (or jupter notebook, jupyter qtconsole) script that allows interactive
testing and experimenting with the multistage turbine model.  This script
reusues the automated test script to construct, initialize, and solve the
initial model.  From there you can modify input, run the model, and inspect the
model stucture.

To use this script:
    [1]: %gui qt
    [2]: %run ipy_test_script.py #or path to script
    [3]: m = test_initialize() #construct, initialize, and solve pressure-driven
    [4]: ui = get_mainwindow_nb(model=m) #optionally open model viewer
    [5]: #have a lot of fun
"""

from pyomo.environ import value
from pyomo.contrib.viewer.ui import get_mainwindow_nb

from idaes.unit_models.power_generation.tests.test_turbine_multistage import \
    test_initialize, solver
