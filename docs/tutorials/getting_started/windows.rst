Windows Installation Guide
==========================

Quickstart
----------
The following commands should be sufficient to get you started with installing and using IDAES. For more information on these steps, see below.

.. sidebar:: Learn more about the Installation
   :class: quickstart-toc

   * `Install Prerequisites`_
   * `Install IDAES-PSE`_
   * `Install IDAES Extensions`_
   * `Install IDAES Examples`_
   * `Run IDAES Tests`_

.. Import quick start guide, including OS specific lines & skipping non-OS lines & comments
.. literalinclude:: install_templates/quickstart.txt
   :language: bash

------------------------------------------------

Installing IDAES
----------------
To get IDAES fully set up on your machine, we'll go through the steps to get idaes-pse package installed as well as setting up the IDAES extensions, which includes some extra solvers and function libraries, the IDAES example files, and the IDAES tests.

.. include:: install_templates/idaes_data.txt

Install Prerequisites
^^^^^^^^^^^^^^^^^^^^^

**Install Miniconda**

1. Download & install `Miniconda <https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe>`_.
2. Install anaconda from the downloaded & open the Anaconda Prompt (Start -> "Anaconda Prompt").

Install IDAES-PSE
^^^^^^^^^^^^^^^^^^

.. include:: install_templates/conda_idaes_pse.txt

Install IDAES Extensions
^^^^^^^^^^^^^^^^^^^^^^^^

.. include:: install_templates/extensions.txt

Install IDAES Examples
^^^^^^^^^^^^^^^^^^^^^^

.. include:: install_templates/examples.txt

Run IDAES Tests
^^^^^^^^^^^^^^^

.. include:: install_templates/run_tests.txt
