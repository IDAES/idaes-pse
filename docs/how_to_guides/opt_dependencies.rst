.. _optional-dependencies:

Installing optional dependencies
================================

Depending on the installation method (Conda, pip) and the IDAES version, not all of IDAES's dependencies will be installed by default.

Installing these *optional dependencies* might require additional steps, described below, in addition to the default installation steps described elsewhere in the documentation.

For pip installations
^^^^^^^^^^^^^^^^^^^^^

.. important:: Users who installed IDAES by running ``pip install idaes-pse`` should follow these steps.

When installing IDAES using pip, optional dependencies can be installed by specifying one or more "``extras_require`` targets":

.. code-block:: bash

   # install the `ui` target
   pip install "idaes-pse[ui]"

   # install the `ui` and `dmf` targets
   pip install "idaes-pse[ui,dmf]"

.. important:: The ``pip install`` argument should be wrapped in double quotes (``"``) when it contains square brackets. Otherwise, it might cause an error to be reported by the shell/command interface used to invoke the command.

Available optional dependencies targets
---------------------------------------

As of IDAES 2.1, the following ``extras_require`` targets are available:

* ``ui``: for the `IDAES User Interface (UI) <https://github.com/IDAES/idaes-ui>`_ (`UI documentation <https://idaes-ui.readthedocs.io/en/latest/>`_)
* ``dmf``: for the :ref:`Data Management Framework <dmf-overview>`
* ``grid``: for the :ref:`IDAES Grid integration <idaes-grid>`
* ``omlt``: for the :ref:`OMLT integration <omlt>`
* ``coolprop``: for the :py:mod:`idaes.models.properties.modular_properties.coolprop` property package

Specifying optional dependencies when installing a prerelease version of IDAES
------------------------------------------------------------------------------

.. code-block:: bash

   pip install "idaes-pse[ui,dmf] @ git+https://github.com/IDAES/idaes-pse@main"

For Conda installations
^^^^^^^^^^^^^^^^^^^^^^^

.. important:: These steps apply to users who installed IDAES by running ``conda install <...>``. They **do not apply** if IDAES was installed using pip.

As of IDAES 2.1, the IDAES Conda package available already includes several dependencies that are optional for pip installations.
This difference is due to a variety of reasons, including the fact that not all optional dependencies are available as Conda packages (e.g. OMLT); or that Conda does not have an equivalent of ``extras_require`` targets to specify optional dependencies at install-time.
The consequence of this is that some optional dependencies might be unavailable if IDAES was installed using Conda.

Therefore, we encourage users who encounter issues after a Conda installation of IDAES to **try installing IDAES again using pip instead** and see if the issues are resolved.
