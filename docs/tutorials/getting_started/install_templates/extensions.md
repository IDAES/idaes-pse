Run the :doc:`idaes get-extensions command<../../../reference_guides/commands/get_extensions>` to install the compiled binaries. These binaries include solvers and function libraries. The PETSc solver is optional and can be omitted. See :ref:`Binary Packages <tutorials/getting_started/binaries:Binary Packages>` for more details :: 

   idaes get-extensions --extra petsc

.. note:: If you are not able to successfully run the ``idaes get-extensions``
          command due to network security settings or another reason, you can
          download binary release files from
          https://github.com/IDAES/idaes-ext/releases, and extract them in the
          directory indicated by the ``idaes bin-directory`` command. You will
          need both the ``idaes-lib-*`` and ``idaes-solvers-*`` files
          appropriate for your operating system.