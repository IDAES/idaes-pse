.. _modelvis::

JupyterLab Model Visualization Tool
===================================

.. note::
    The IDAES-modelvis viewer requires the use of JupyterLab. 


Overview
--------

The *idaes.dmf.ui* subpackage enables the serialization of flowsheets to "idaes.vis" 
files, which in conjunction with the idaes-modelvis viewer produces interactive visual 
representations of flowsheets. The resultant diagrams can be rearranged and saved.


Instructions
------------

1. :ref:`Ensure that the latest IDAES is installed. <idaes_installation>` 

2. `Install JupyterLab. <https://jupyterlab.readthedocs.io/en/stable/getting_started/installation.html>`_
   If you are not using Conda environments, use the `pip install` instructions.
  
3. Build and install the IDAES-modelvis JupyterLab extension.

.. code-block:: sh

	cd <repository>/ui/modelvis/idaes-model-vis
    npm install
    npm run build
    jupyter labextension link .

4. Launch JupyterLab (run ``jupyter lab`` from a folder you wish to work out of).

5. In Python (in a Notebook tab within JupyterLab), construct a flowsheet as usual 
(add unit models, set connections, etc.).

6. In Python (in the same Notebook), import the flowsheet serializer module, create an instance
of the flowsheet serializer, and run the serialization as shown below:
     
.. code-block:: python

    from idaes.dmf.ui.flowsheet_serializer import FlowsheetSerializer
    flowsheet_serializer = FlowsheetSerializer()
    # Substitute 'm.fs' with the variable holding your flowsheet
    # and 'myflowsheetname' with a file name to save to.
    flowsheet_serializer.save(m.fs, 'myflowsheetname')


7. A ``.idaes.vis`` file should be created with the chosen filename, and
become visible in the JupyterLab file browser. If there is an existing
file with the same name, you must either choose a different filename
or use ``FlowsheetSerializer.save()`` with the additional optional argument
``overwrite=True`` (in which case the file will be overwritten).

8. Open the created ``.idaes.vis`` file in JupyterLab. A tab should open and display
a graph representation of the serialized flowsheet; the components are
tiled diagonally by default, and can be rearranged to your liking.

9. The layout of the graph can be saved into the serialized file by using JupyterLab's
``File->Save`` menu item. Autosaving can also be configured by using JupyterLab's
``Settings->Advanced Settings Editor`` option under ``Document Manager``.


Developer notes
---------------

Rebuilding
^^^^^^^^^^

After making changes to the TypeScript, rebuild the extension and reinstall it into JupyterLab:

.. code-block:: sh

    npm run build
    jupyter lab build
