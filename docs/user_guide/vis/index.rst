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

OK, now that we have a new notebook, we need to create an IDAES flowsheet to visualized. Since the details of how
to do this are not our main concern, just cut and paste the following code (taken from the Flash unit tutorial)
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

This means the initialization succeeded, and you can move on. If you see, something else -- well, many things could
have gone wrong -- but one common problem is you see something like this ::

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
review the :ref:`installation instructions <IDAES Installation>`, close this notebook and exit the current running Jupyter, and try again to run
the command "jupyter notebook" right after setting up your environment to be the one in which you installed IDAES.

Before continuing, save the notebook ('File -> Save' or Ctrl-S) with an appropriate name, like "Hello World".

When you ran the cell above, it created a new blank cell for you to continue editing. We will use this cell to visualize
our initialized (but not solved) model. In the new cell, type in and run (shift-enter) the following code::

    m.fs.visualize("Hello, World")

This will produce some diagnostic output, but more importantly it should create a new browser tab or window with the
IFV displaying the flowsheet. This will look similar to this:

.. image:: ../../_images/ifv_helloworld_1.png
    :width: 800

You'll notice that the layout is not too great, and in fact if you look closely you'll see that the components
have just been placed in a diagonal. You can try rearranging the diagram with the mouse (the components can all
be moved), and for more details on the available functions, see the next section.

.. TODO Tell user how to see values on the unit model and streams

But before that, try one more thing.
Go back to the Jupyter notebook and solve the Flash optimization problem, but adding and running a new cell with
the following code::

    solver = SolverFactory('ipopt')
    status = solver.solve(m, tee=True)

Since this is such a simple optimization problem, the solver should finish very quickly with, at the bottom
of its diagnostic output, the message ``EXIT: Optimal Solution Found``.

Now go back to your IFV window and click on "Refresh Graph".

.. TODO Point out how the values have changed to reflect solved model

User Guide
----------
.. More detailed guided tour of functionality -- leverage notebook from quickstart

Reference
---------
.. Alphabetical reference of functionality

OLD CONTENT:

Overview
--------

The Flowsheet Visualizer is a service that starts a flask server to
display an interactive webpage with the current flowsheet's unit models and
arcs. Users may manipulate the display by clicking and dragging the unit 
models, streams, and stream labels.

Installation instructions
-------------------------

1. :ref:`Ensure that the latest IDAES is installed. <getting_started/index:Installation>` 

.. _usage:

Usage
-----

1. Create a flowsheet in a Jupyter Notebook. For the purpose of these 
   instructions the model will be `m` and the flowsheet will be `m.fs`

2. Call the method `visualize()` from the flowsheet with a model name 
   as a string:
   `m.fs.visualize('model_name')`

.. image:: ../../_images/modelvis/fs_visualize_jupyter_notebook.png

3. A webpage should display:

.. image:: ../../_images/modelvis/initial_layout.png

If a webpage does not display then copy and
paste the URL that outputs from the visualize command:

.. image:: ../../_images/modelvis/circled_url.png

4. Manipulate the layout of the model display as desired:

.. image:: ../../_images/modelvis/modified_layout.png

5. If the flowsheet is later modified, click the Refresh Graph button to
   see the changes.

The displayed layout is preserved as much as possible, with new components
appearing along a diagonal line. 

.. note::
    This feature is still under development. 
    Several types of changes to the flowsheet currently cause the entire user-
    modified layout to be lost. Consider saving the layout often (see below).

.. image:: ../../_images/modelvis/new_unit_model_layout.png

6. Save the displayed layout using the save button on the visualization page. 
   This writes the visualization to a file in the user's home directory under 
   `.idaes/viz` using the model name provided to `visualize()`. 
   In this example the filename would be `model_name.viz`.

.. _streamlabels:

Stream Labels
-------------

The initial layout loads with the stream labels hidden. Show or hide all of 
the stream labels by clicking the button with the speech bubbles, 
on the toolbar.

Show or hide an individual label by right clicking on the stream or its label.

.. _miscfeatures:

Misc. Features
--------------

* Right click on an icon to rotate it by 90 degrees.

* Create anchor points on a stream by left clicking on the stream. The stream 
  will be forced to connect through each anchor point, typically adding right angles.

