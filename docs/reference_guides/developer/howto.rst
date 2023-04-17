
.. _idaes-contrib-guide:

IDAES Contributor Guide
========================

About
-----
This page tries to give all the essential information needed
to contribute software to the IDAES project. It is designed
to be useful to both internal and external collaborators.

Code and other file locations
-----------------------------
Source code
    The main Python package is under the `idaes/` directory.
    Sub-directories, aka subpackages, should be documented elsewhere.
    If you add a new directory in this tree, be sure to add a `__init__.py` in that directory
    so Python knows it is a subpackage with Python modules.
    Code that is not part of the core package is under `apps/`. This code can have any
    layout that the creator wants.

Documentation
    The documentation for the core package is under `docs`.

Examples
    Examples are under the `examples/` directory.
    Tutorials from workshops are under the `examples/workshops/` subdirectory.
    

Developer environment
---------------------
Development of IDAES will require an extra set of required package not needed by regular users.
To install those extra developer tools use the command ``pip install -r requirements-dev.txt``
rather than ``pip install -r requirements.txt``


Code style
------------
The code style is not entirely consistent. But some general guidelines are:

* follow the `PEP8`_ style (or variants such as `Black`_)
* use `Google-style`_ docstrings on classes, methods, and functions
* format your docstrings as `reStructuredText`_ so they can be nicely rendered as HTML by Sphinx
* check your spelling using `crate-ci-typos`_
* add logging to your code by creating and using a global log object named
  for the module, which can be created like: ``_log = logging.getLogger(__name__)``
* take credit by adding a global author variable: ``__author__ = 'yourname'``

.. _PEP8: https://www.python.org/dev/peps/pep-0008/
.. _Black: https://github.com/python/black
.. _Google-style: https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html
.. _reStructuredText: http://docutils.sourceforge.net/rst.html
.. _crate-ci-typos: https://github.com/crate-ci/typos

Tests
-----
For general information about writing tests in Python, see :ref:`tst-top`.

There are three types of tests:

Python source code
    The Python tests are integrated into the Python source code directories.
    Every package (directory with `.py` modules and an `__init__.py` file)
    should also have a `tests/` sub-package, in which are test files. These,
    by convention are named `test_<something>.py`.

Doctests
    With some special reStructuredText "directives" (see "Writing tests"), the documentation
    can contain tests. This is particularly useful for making sure examples in the
    documentation still run without errors.

Jupyter notebook tests
    (coming soon)


Writing tests
^^^^^^^^^^^^^
We use `pytest`_ to run our tests. The main advantage of this framework over
the built-in `unittest` that comes with Python is that almost no boilerplate
code is required. You write a function named `test_<something>()` and,
inside it, use the (pytest-modified) `assert` keyword to check that things
are correct.

Writing the Python unit tests in the `tests/` directory is,
hopefully, quite straightforward.
Here is an example (out of context) that tests a couple of 
things related to configuration in the core unit model library::

    def test_config_block():
        m = ConcreteModel()

        m.u = Unit()

        assert len(m.u. config) == 2
        assert m.u.config.dynamic == useDefault

See the existing tests for many more examples.

For tests in the documentation, you need to wrap the test itself
in a directive called `testcode`. Here is an example::

    .. testcode::

        from pyomo.environ import *
        from pyomo.common.config import ConfigValue
        from idaes.core import ProcessBlockData, declare_process_block_class

        @declare_process_block_class("MyBlock")
        class MyBlockData(ProcessBlockData):
            CONFIG = ProcessBlockData.CONFIG()
            CONFIG.declare("xinit", ConfigValue(default=1001, domain=float))
            CONFIG.declare("yinit", ConfigValue(default=1002, domain=float))
            def build(self):
                super(MyBlockData, self).build()
                self.x = Var(initialize=self.config.xinit)
                self.y = Var(initialize=self.config.yinit)

First, note that reStructuredText directive and indented Python code. The indentation of the
Python code is important. You have to write an entire program here, so all the
imports are necessary (unless you use the `testsetup` and `testcleanup` directives,
but honestly this isn't worth it unless you are doing a lot of tests in one file).
Then you write your Python code as usual.

Running tests
^^^^^^^^^^^^^
Running all tests is done by, at the top directory, running the command: ``pytest``.

The documentation test code will actually be run by a special hook in the pytest configuration that
treats the Makefile like a special kind of test.
As a result, *when you run pytest in any way
that includes the "docs/" directory (including the all tests mode), then all the documentation tests will run,
and errors/etc. will be reported through pytest*. A useful corollary is that, to run
documentation tests, do: ``pytest docs/Makefile``

You can run specific tests using the pytest syntax, see its documentation or ``pytest -h`` for details.

.. _pytest: https://docs.pytest.org/en/latest/

Documentation
--------------
The documentation is built from its sources with a tool called Sphinx.
The sources for the documentation are:

* hand-written text files, under `docs/`, with the extension ".rst" for `reStructuredText`_.
* the Python source code
* selected Jupyter Notebooks 

Building documentation
^^^^^^^^^^^^^^^^^^^^^^

.. note:: To build the documentation locally, you will need to have the Sphinx tools installed.
       This will be done for you by running ``pip install requirements-dev.txt`` ("developer" setup)
       as opposed to the regular ``pip install requirements.txt`` ("user" setup).

To build the documentation locally, use our custom `build.py` script.

    cd docs
    python build.py

The above commands will do a completely clean build to create HTML output.

If the command succeeds, the final line will look like::

    === SUCCESS

If it fails, it will instead print something like::

    *** ERROR in 'html'
    ***
    *** message about the command that failed
    *** and any additional info
    ***

If you want to see the commands actually being run, add `-v` to the command line.

By default the build command removes all existing built files before running the
Sphinx commands. To turn this off, and rebuild only "new" things, add `--dirty`
to the command line.

Previewing documentation
^^^^^^^^^^^^^^^^^^^^^^^^
The generated documentation can be previewed locally by opening
the generated HTML files in a web browser. The files are under the `docs/build/`
directory, so you can open the file ``docs/build/index.html`` to get started.
