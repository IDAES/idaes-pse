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
.. literalinclude:: install_templates/quickstart.md
   :language: bash
   :linenos:
   :start-after: modification warning
   :lines: 1-6,8,11-

------------------------------------------------

Installing IDAES
----------------
To get IDAES fully set up on your machine, we'll go through the steps to get idaes-pse package installed as well as setting up the IDAES extensions, which includes some extra solvers and function libraries, the IDAES example files, and the IDAES tests.

Install Prerequisites
^^^^^^^^^^^^^^^^^^^^^

**Install Miniconda**

1. Download & install `Miniconda <https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe>`_.
2. Install anaconda from the downloaded & open the Anaconda Prompt (Start -> "Anaconda Prompt").

Install IDAES-PSE
^^^^^^^^^^^^^^^^^^

.. include:: install_templates/conda_idaes_pse.md

Install IDAES Extensions
^^^^^^^^^^^^^^^^^^^^^^^^

.. include:: install_templates/extensions.md

Install IDAES Examples
^^^^^^^^^^^^^^^^^^^^^^

.. |os_specific_fpath| replace:: `C:\\Users\\MyName\\IDAES\\Examples`
.. include:: install_templates/examples.md

Run IDAES Tests
^^^^^^^^^^^^^^^

.. include:: install_templates/run_tests.md
