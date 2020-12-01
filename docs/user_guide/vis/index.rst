.. _IFV:

IDAES Flowsheet Visualizer
===========================

Concepts
--------
The IDAES Flowsheet Visualizer, or IFV for short, is a graphical user interface to help understand an IDAES flowsheet.
It uses the visual language of traditional process engineering diagrams, somewhat simplified, to display
the components of a given flowsheet,
their connections, and values associated with unit models and streams. At this time, it is read-only, i.e. nothing
you can do with the flowsheet will alter the underlying IDAES models. On the other hand, the tool is wired so that
changes in the underlying models can be shown in real-time in the IFV.

Quickstart
----------
This section will walk you through constructing a basic flowsheet and visualizing it with the IFV.

Before starting to use the IFV, you must have installed IDAES. For detailed instructions on how to do this,
see :ref:`IDAES Installation`.

In this quickstart we will use a Jupyter Notebook to compose our flowsheet and view it with the IFV. So, first
you need to open a new, blank, Jupyter Notebook. You can start the Jupyter Notebook (server) from the terminal command-line with
the command ``jupyter notebook``. See the `official Jupyter Notebook documentation <https://jupyter-notebook.readthedocs.io/>`_
for detailed instructions. The Jupyter Notebook is a web application, so you should get a web page (often as a tab
on an open browser window) that shows the contents of the directory from where you issued the command. Navigate to
a directory where you want to save your work, then click on the "New" button in the upper-right corner to create
a new notebook. You may have multiple options under the "Python" sub-menu; make sure you select a Python installation
that has IDAES already installed.

We need to create, in this Jupyter Notebook, an IDAES flowsheet to visualize. Paste the following code (taken from the Flash unit tutorial)
into the first cell in the notebook::

    from pyomo.environ import ConcreteModel, SolverFactory, Constraint, value
    from idaes.core import FlowsheetBlock
    from idaes.generic_models.properties.activity_coeff_models.BTX_activity_coeff_VLE \
        import BTXParameterBlock
    from idaes.generic_models.unit_models import Flash
    # Model and flowsheet
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    # Flash properties
    m.fs.properties = BTXParameterBlock(default={"valid_phase": ('Liq', 'Vap'),
                                                "activity_coeff_model": "Ideal",
                                                "state_vars": "FTPz"})
    # Flash unit
    m.fs.flash = Flash(default={"property_package": m.fs.properties})
    m.fs.flash.inlet.flow_mol.fix(1)
    m.fs.flash.inlet.temperature.fix(368)
    m.fs.flash.inlet.pressure.fix(101325)
    m.fs.flash.inlet.mole_frac_comp[0, "benzene"].fix(0.5)
    m.fs.flash.inlet.mole_frac_comp[0, "toluene"].fix(0.5)
    m.fs.flash.heat_duty.fix(0)
    m.fs.flash.deltaP.fix(0)
    # Initialize
    m.fs.flash.initialize()


Run that cell (hit shift-enter), and you should see some output from the last call to ``initialize()``, ending with something like ::

    2020-11-24 08:41:53 [INFO] idaes.init.fs.flash: Initialization Complete: optimal - Optimal Solution Found

This means the initialization succeeded, and you can move on. If you see instead something like this ::

    ---------------------------------------------------------------------------
    ModuleNotFoundError                       Traceback (most recent call last)
    <ipython-input-2-ad88ea7f7460> in <module>
    ----> 1 from pyomo.environ import ConcreteModel, SolverFactory, Constraint, value
          2 from idaes.core import FlowsheetBlock
          3 from idaes.generic_models.properties.activity_coeff_models.BTX_activity_coeff_VLE \
          4     import BTXParameterBlock
          5 from idaes.generic_models.unit_models import Flash

    ModuleNotFoundError: No module named 'pyomo'

In this case, you probably didn't run the notebook within a Python environment where IDAES is installed. See if
there are other Python environments, or "kernels", you can run under the 'Kernel -> Change kernel' menu. If not,
review the :ref:`installation instructions <IDAES Installation>`, close this notebook and exit the current running
Jupyter, and try again to run the command "jupyter notebook" right after setting up your environment to be the one in which you installed IDAES.

Before continuing, save the notebook ('File -> Save' or Ctrl-S) with an appropriate name, like "Hello World".

When you ran the cell above, it created a new blank cell for you to continue editing. We will use this cell to visualize
our initialized (but not solved) model. In the new cell, type in and run (shift-enter) the following code::

    m.fs.visualize("Hello, World", save_as="hello_world.json")

This will create a new browser tab or window with the IFV displaying the flowsheet:

.. image:: ../../_images/ifv_helloworld_1.png
    :width: 800

For the initial layout, the components
have just been placed in a diagonal. You can rearrange the diagram with the mouse (the components can all
be moved), and for more details on the available functions, see the next section. If you hit "Save", the IFV will save
your changes in the layout to the destination that you passed to "save_as", in this case the file
"hello_world.json", in the current directory.

.. TODO Tell user how to see values on the unit model and streams

User Guide
----------
.. More detailed guided tour of functionality -- leverage notebook from quickstart

Reference
---------
.. Alphabetical reference of functionality



