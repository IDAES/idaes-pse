<!-- BEGIN Status badges -->
[![CircleCI](https://circleci.com/gh/IDAES/idaes.svg?style=svg&circle-token=863b50f5cba60e880ef338300fa545a3a140973c)](https://circleci.com/gh/IDAES/idaes)
[![Coverage Status](https://coveralls.io/repos/github/IDAES/idaes/badge.svg?branch=master&t=9GISoZ)](https://coveralls.io/github/IDAES/idaes?branch=master)
<!-- END Status badges -->

# IDAES PSE Framework

The IDAES PSE Framwork aims to provide multi-scale, simulation-based, open source
computational tools and models to support the design, analysis, optimization,
scale-up, operation and troubleshooting of innovative, advanced energy systems.

## Contributing

**By contributing to this repository, you are agreeing to all the terms set out
in the LICENSE.txt and COPYRIGHT.txt files in this directory.**

## Getting Started
For installation instructions, please refer to the [online documentation](https://idaes.github.io/idaes-pse/).

The documentation for IDAES is built using [Sphinx](http://www.sphinx-doc.org/). To generate the HTML version of the documentation, first make sure Sphinx is installed for your version of Python, then go to the "docs/" subdirectory and run the _Makefile_:

```
cd docs
make html
```

To view the documentation you just built, open the file
`doc/build/html/index.html` in a web browser.


## Running tests

After you install, you can run tests to make sure everything is working. We use [pytest](https://pytest.org/) for testing and generating code coverage reports.  The `pytest` command should be available in the conda environment created by running the `install.sh` script as described in the installation instructions.

To run tests against the core modules and DMF (not `idaes/contrib`), and generate a coverage report:

```
$ source activate <idaes_conda_env>  # If you used "install.sh <idaes_conda_env>" to install
$ pytest
```

To run tests in `idaes/contrib` just add that to the pytest command:

```
$ pytest ideas/contrib  # These are not guarenteed to all succeed...
```

If there are errors, or you are having trouble, you can use our [issue tracker on Github](https://github.com/IDAES/idaes/issues) to look for other users experiencing similar problems, or to report a new bug.


## Running a notebook

There are example [Jupyter](https://jupyter.org) notebook(s) in the `examples/` directory. To run them, you should invoke Jupyter on a Notebook file (these end in the extension `.ipynb`).

```
jupyter notebook examples/run-mea-model.ipynb
```

This should start up a notebook server and then pop up a tab or window in your default web browser showing the Notebook. For more information on how to use Jupyter, see the "Help" menu in the Notebook window itself, and the extensive documentation on the [Jupyter website](https://jupyter.org). 

## Contacts and more information

Please see the
[IDAES main website](https://www.idaes.org) for general information
and people to contact.
