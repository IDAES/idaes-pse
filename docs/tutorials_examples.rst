Tutorials and Examples
======================

Viewing online
--------------
Tutorials and examples for IDAES are located on the |examples-site|.
If you are new to IDAES, it is strongly recommended to start with the tutorials.
There are also pre-recorded tutorial videos on the `IDAES page on YouTube <https://www.youtube.com/channel/UCpp3J_990C0Oz_CbxRDfr6g>`_.

Running locally
---------------
To run and use the examples on your own computer, once you have installed IDAES,
you should run ``idaes get-examples`` in a command-line shell.
Please see :ref:`this page<user_guide/commands/get_examples:idaes get-examples: Fetch example scripts and Jupyter Notebooks>` for details on how to use this program.

Once you have installed the examples, change the directory where you downloaded them and
run the following command\ :sup:`1` ::

        jupyter notebook notebook_index.ipynb

This will open a descriptive Jupyter Notebook with links that allow you to open and run any of
the other tutorial and example notebooks. Happy modeling!

:sup:`1` If you have installed `JupyterLab <https://jupyterlab.readthedocs.io/en/stable/index.html>`_,
you can use it by simply substituting "jupyter lab" for "jupyter notebook" in the given command.

Additional information
----------------------
The sources for the tutorials and examples are maintained on the
`IDAES examples repository <https://github.com/IDAES/examples-pse>`_.

If you want to develop custom unit and property models, please refer to the
:ref:`advanced user guide <advanced_user_guide/index:Advanced User Guide>`.


