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
   :lines: 1-5,7,10-

------------------------------------------------

Installing IDAES
----------------

Install Prerequisites
^^^^^^^^^^^^^^^^^^^^^

**Install Miniconda (optional)**

1. Download & install `Anaconda <https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe>`_.
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
