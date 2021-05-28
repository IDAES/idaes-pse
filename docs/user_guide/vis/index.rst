.. _IFV:

IDAES Flowsheet Visualizer
===========================

.. contents::
    :depth: 2

Introduction
------------
The IDAES Flowsheet Visualizer, or IFV for short, is a web-based user interface (UI) that lets you:

* View any IDAES flowsheet as a process engineering diagram
* Export flowsheet diagrams as images (`SVG <https://www.w3.org/Graphics/SVG/>`_ format)
* View and export the "stream table" for the flowsheet
* Rearrange the flowsheet diagram to your taste and save the arrangement for next time
* Dynamically refresh the displayed values to reflect changes in the underlying IDAES model

To use the IFV, first :ref:`install IDAES <IDAES Installation>`.
The IFV can be invoked from a Jupyter Notebook or a Python script. It does not require that you run any
other application. Currently the IFV is only for viewing the flowsheet on your own computer. :sup:`1`

Starting and stopping the IFV is fast and does not consume many resources.

:sup:`1` *But, since it is a web application, a shared service for viewing flowsheets stored remotely is definitely possible.*

Guide
-----
This guide describes how to invoke (i.e., start) and use the IFV.

Invocation
^^^^^^^^^^
The IFV visualizes *flowsheets*. To get started with creating a flowsheet with IDAES, see the
:ref:`Flowsheet models<user_guide/components/flowsheet/index:Flowsheet>` documentation page.
Once you have created your flowsheet, simply call the :ref:`visualize <visualize-function>` method on that object,
passing some parameters to give it a name and optional file for saving changes::

    # First, create your IDAES model, "m", which has an attribute ".fs" for the flowsheet
    # Then, invoke the `visualize` method
    m.fs.visualize("My Flowsheet")

The invocation of the `visualize` method will pop up a browser tab or window with the UI, displaying the
flowsheet and, if the information is available, the stream table. In the notebook or script, you can continue
to run more code and the UI will continue to work in the background. You can close the UI at any time. If you
exit the script or notebook while the UI is running, you can still manipulate the diagram and stream table,
but you will not be able to save
or refresh, since these require communication with the Python process that no longer exists.

There are three ways to invoke the `visualize` functionality, which in the end do the same thing and
have the same arguments.

1. Use the `visualize` method of a flowsheet (as above)
2. Call the `visualize` function from the package `idaes.ui.fsvis`, passing it a flowsheet object
3. Call the same `visualize` function from the module `idaes.ui.fsvis.fsvis`, passing it a flowsheet object

In all cases, the arguments and behavior are the same.
See the :ref:`visualize function documentation <visualize-function>` for details on parameters to this function.

.. note:: You can continue to modify and use your model and flowsheet after calling `visualize()`. This
  may update the visualization, but nothing you do in the IFV will affect the Python model; it is *read-only*.

User Interface
^^^^^^^^^^^^^^

This section describes how to use the graphical web user interface. We start with a screenshot of the UI, with
the main areas highlighted. Then we zoom in on each area and describe how to use it.

.. _ifv-screenshot:

.. only:: html

    **In the screenshot of the IFV UI below, you can click on one of the indicated areas to jump to the associated
    description.**

    .. raw:: html

        <img src="../../_static/extra/ifv_screenshot_window.png" usemap="#image-map" width="800px">
        <map name="image-map">
            <area target="_self" alt="Top bar" title="Top bar" href="#top-bar" coords="20,10,800,70" shape="rect"></area>
            <area target="_self" alt="Diagram Controls" title="Diagram Controls" href="#diagram-controls" coords="510,80,800,125" shape="rect"></area>
            <area target="_self" alt="Diagram/Flowsheet" title="Diagram/Flowsheet" href="#diagram" coords="25,130,800,570" shape="rect"></area>
            <area target="_self" alt="Stream Table" title="Stream Table" href="#stream-table" coords="25,570,800,890" shape="rect"></area>
        </map>

.. only:: latex

    .. figure:: ../../static/extra/ifv_screenshot_window.png
        :width: 800

        Screenshot of the main window of the IFV UI

Top bar
+++++++

.. figure:: /images/ifv_screenshot_topbar.png
    :width: 600

    Screenshot of the top bar of the IFV UI

The *top bar* has a title bar, which contains the IDAES logo and the name of the flowsheet being visualized,
and a menu. The menu items are:

.. _export-menu:

* **Refresh**: Update the view with any changes made to the flowsheet from the Python side.
  This also has the effect of saving the current layout.
* **Save**: Save the current layout to the data store that was specified with the visualization
  was launched. Note that this does *not* update with any changes made to the flowsheet in Python (use *Refresh* for that).
  Neither does it have any effect on the Python flowsheet values, as the IFV cannot modify the underlying flowsheet.
* **Export**: Save the flowsheet or stream table as a file.

  * **Flowsheet**: Save the flowsheet as a Scalable Vector Graphics (SVG) file, a common format for
    images that consist of "vector" elements like boxes, lines, and text. SVG files can be viewed like images
    by most programs that allow image viewing, and even edited with a program like `Inkscape <https://inkscape.org/>`_.
    You will get a preview of the image and a "Download" button that will save in a file named for the
    flowsheet, with a ".svg" extension.
  * **Stream table**: Save the flowsheet as comma-separated values. The result will be a text file, called "export.csv",
    that contains the data.

* **View**: Toggle the visibility of the flowsheet (diagram) area and the stream table area.
* **Help**: Load this documentation page.

:ref:`Back to main window screenshot <ifv-screenshot>`

Diagram
++++++++

.. figure:: /images/ifv_screenshot_diagram.png
    :width: 800

    Screenshot of the main *diagram* (or flowsheet) area of the IFV UI

The *diagram* (or *flowsheet*) area lets you rearrange the flowsheet as you need and zoom in on particular sections.
You can interact with the components on the diagram:

Shapes
    Geometric shapes on the flowsheet represent unit models, inlets and outlets, and other IDAES components.
    They are connected by lines, and each has a name. All shapes can be moved by clicking and dragging them.
    If you right-click on a shape, it will rotate 90 degrees.

Lines
    The lines connecting units can be manipulated by clicking and dragging.
    You can click on a line to create a new segment that can be used for routing the line
    around objects. You can eliminate a segment by clicking on the dot that appears as you hover over
    the line. There are also pill-shaped handles that appear on the lines for moving them.
    The endpoints of the lines are determined by the flowsheet and cannot be changed. For the same reason, you
    also cannot add or remove lines.

Labels
    Both the shapes and lines have associated values that can be shown, which pop up over the
    lines if you toggle the "Show labels" control. See the :ref:`diagram-controls` section for details.

More details on mouse and keyboard actions for the diagram are available in the `documentation of the
underlying Rappid toolkit <https://google.com/>`_.

:ref:`Back to main window screenshot <ifv-screenshot>`

.. _diagram-controls:

Diagram Controls
++++++++++++++++

.. figure:: /images/ifv_screenshot_diagramcontrols.png
    :width: 400

    Screenshot of the *diagram controls* area of the IFV UI

The *diagram controls* allow you to affect some global properties of the diagram/flowsheet area.

.. |zoomin| image:: /images/ifv_icon_zoomin.png
.. |zoomout| image:: /images/ifv_icon_zoomout.png
.. |zoomfit| image:: /images/ifv_icon_fit.png

View actions
  * Labels: Toggle visibility of the information (*labels*) shown for each stream. This is the same information
    that appears in the :ref:`ifv_streamtable`.
  * Grid: Toggle a background "grid"
  * |zoomin|: Zoom in by 25%
  * |zoomout|: Zoom out by 25%
  * |zoomfit|: Fit the diagram into the current area

:ref:`Back to main window screenshot <ifv-screenshot>`

.. _ifv_streamtable:

Stream Table
++++++++++++

The IFV will show a stream table with variables defined for each stream in the flowsheet, if these
values exist and the flowsheet adheres to the IDAES conventions for naming the inlet and outlet
streams. An example of a stream table is shown below.

.. figure:: /images/ifv_screenshot_streamtable.png
    :width: 800

    Screenshot of an example stream table

There are a number of ways of manipulating this table:

  * The "Hide Fields" pull-down menu provides a list of stream names.
    Select a name to hide/show that column in the table.
  * Click on the column header and drag it left or right to change its order in the table.
  * Resize a column by hovering over a column border until you see the mouse pointer change, then drag it to resize.

You can also export the entire table as a file of comma-separated values. See the :ref:`Export <export-menu>`
documentation for details.

:ref:`Back to main window screenshot <ifv-screenshot>`


Reference
---------

.. currentmodule:: idaes.ui.fsvis

.. _visualize-function:

.. autofunction:: visualize

.. autodata:: VisualizeResult

Software notes
^^^^^^^^^^^^^^
This section provides some additional details for developers or users more interested in the
programming details.

.. _ifv-architecture:

Client/server architecture
++++++++++++++++++++++++++
The ``visualize()`` command works by starting an HTTP server in a separate thread, and serving
requests from the UI (or any other requester). The server only responds to requests from your computer,
not the internet. When you exit the script or Jupyter Notebook that called `visualize` then you will also
stop the server -- and the associated IFV page will no longer be able to save or refresh the flowsheet.
The architecture diagram is shown below.

.. note ; the figure below was generated with asciiflow infinity, but is just text and can
.. be edited in any way. For HTML and PDF these are rendered as nice little diagrams by the
.. Sphinx plugin "sphinxcontrib.aafig" using the Python "aafigure" https://pypi.org/project/aafigure/ package

.. code-block:: text

    +-------------------+                        +--------------------+
    |                   |                        |    Web browser     |
    |  'Python script'  |                        +--------------------+
    |  'or Jupyter'     |    +---------------+   | IFV web interface  |
    |  'Notebook'       |    | 'HTTP server' |   +--------------------+
    |                   |    | 'running in'  |   |    +--+            |
    |                   |    | 'a separate'  |   |    +--+            |
    |                   |    | 'thread'      |   |      |      +--+   |
    |                   |    |               <--->      +----> +--+   |
    |  m.fs.visualize   +---->  Load/Save    |   |                    |
    |                   |    |               |   |                    |
    +-------------------+    +-----^---------+   +--------------------+
                                   |
                                   |
                            +------v--------------+
                            |   Local Storage     |
                            +---------------------+

Persistence architecture
++++++++++++++++++++++++
.. py:currentmodule:: idaes.ui.fsvis.persist

The saving of the model uses the the module :mod:`idaes.ui.fsvis.persist`.
This module implements the well-known "|factory-link|", which makes it easy to extend by adding
a new :class:`~.DataStore` sub-class and updating the logic in the factory method,
:func:`~.DataStore.create`, to create and return instances of that class for a given input type.
The input in this case comes from the ``save_as`` argument to the *visualize()* method.

.. |factory-link| raw:: html

    <a href="https://en.wikipedia.org/wiki/Factory_(object-oriented_programming)" target="_blank" style="text-decoration: none;">factory pattern</a>

.. TODO: add an example of extending it, e.g. to save in an S3 bucket
