Generic Install
---------------

The remaining steps performed in either the Linux or OSX Terminal or Powershell.
If you installed Miniconda on Windows use the Anaconda Prompt or Anaconda
Powershell Prompt.  Regardless of OS and shell, the following steps are the same.


Install IDAES Base
^^^^^^^^^^^^^^^^^^

Using pip
"""""""""
.. code-block:: none

   # Install latest stable release
   pip install idaes-pse

.. raw:: html

  <details>
   <summary>Installing other versions of IDAES</summary>

a. To get a previous `IDAES release <https://github.com/IDAES/idaes-pse/releases>`_, for example v1.13

.. code-block:: none

    pip install idaes-pse==1.13

b. To get the latest version from the GitHub main branch

.. code-block:: none

    pip install 'idaes-pse[prerelease] @ https://github.com/IDAES/idaes-pse/archive/main.zip'

c. To get a specific fork or branch, for example myfork (of idaes-pse) and mybranch

.. code-block:: none

    pip install 'idaes-pse[prerelease] @ https://github.com/myfork/idaes-pse/archive/mybranch.zip'

d. For developers: follow the :ref:`advanced user installation<tutorials/advanced_install/index:Advanced User Installation>`.

.. raw:: html

  </details>
  </br>


Using conda
"""""""""""
.. code-block:: none
   
   # Create a new Conda environment (Here named my-idaes-env, but this is up to you)
   conda create --yes --name my-idaes-env python=3.8
   conda activate my-idaes-env

   # Install IDAES Conda package
   conda install --yes -c IDAES-PSE -c conda-forge idaes-pse

Install IDAES extensions, examples, & tests
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
1. Run the :doc:`idaes get-extensions command<../../../reference_guides/commands/get_extensions>`
   to install the compiled binaries. These binaries include solvers and function libraries.
   See :ref:`Binary Packages <tutorials/getting_started/binaries:Binary Packages>`
   for more details::

    idaes get-extensions

..

    .. note:: If you are not able to successfully run the ``idaes get-extensions``
              command due to network security settings or another reason, you can
              download binary release files from
              https://github.com/IDAES/idaes-ext/releases, and extract them in the
              directory indicated by the ``idaes bin-directory`` command. You will
              need both the ``idaes-lib-*`` and ``idaes-solvers-*`` files
              appropriate for your operating system.

..

   .. warning:: The IDAES binary extensions are not yet supported on Mac/OSX.

                As a fallback (assuming you are using a conda env) you can install
                the generic ipopt solver with the command ``conda install -c
                conda-forge ipopt`` though this will not have all the features
                of our extensions package.

2. Run the :doc:`idaes get-examples command <../../../reference_guides/commands/get_examples>` to download
   and install the example files::

    idaes get-examples

..

    By default this will install in a folder "examples" in the current directory.
    The command has many options, but an important
    one is ``--dir``, which specifies the folder in which to install.::
    
       idaes get-examples --dir <PATH>

    An example of the path formatting would look like |os_specific_fpath|

3. Run tests::

    pytest --pyargs idaes -W ignore

4. You should see the tests run and all should pass to ensure the installation worked. You
   may see some "Error" level log messages, but they are okay, and produced by tests for
   error handling. The number of tests that failed and succeeded is reported at the end of the pytest
   output. You can report problems on the |github-issues|
   (Please try to be specific about the command and the offending output.)
