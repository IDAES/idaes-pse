# IDAES Toolkit

The IDAES Toolkit aims to provide multi-scale, simulation-based, open source
computational tools and models to support the design, analysis, optimization,
scale-up, operation and troubleshooting of innovative, advanced energy systems.

<!-- BEGIN Status badges -->
## Build statuses
![GitHub CI (main)](https://github.com/IDAES/idaes-pse/workflows/Pull%20request%20(main)%20CI%20tests/badge.svg)
<!-- END Status badges -->

## System requirements

The code and examples have been tested with the following operating systems:

|Operating system|Supported versions  |
|----------------|--------------------|
| Linux          | Any modern Linux   |
| Windows        | Windows 10         |
| Mac OSX        | Not supported*     |

*For advanced users, Mac OSX installation may be performed with some small changes to the Linux installation instructions.

Most of the functionality is implemented in Python. In accordance with
the end-of-life for many Python 2 libraries, the IDAES Toolkit is written
for Python 3. The following sub-versions are supported:

* Python 3.6
* Python 3.7
* Python 3.7+ (should work, not explicitly tested)

Note that Python 3.5 is *not* supported.

## Contributing

**By contributing to this repository, you are agreeing to all the terms set out
in the LICENSE.txt and COPYRIGHT.txt files in this directory.**

## Quickstart

To install with `pip`:

```bash
pip install idaes_pse
```

To install with Anaconda's `conda`: **coming soon**

## Getting Started
For installation instructions, please refer to the [online documentation](https://idaes-pse.readthedocs.io/en/stable/).

The documentation for IDAES is built using [Sphinx](http://www.sphinx-doc.org/). To generate the HTML version of the documentation, first make sure Sphinx is installed for your version of Python,  
then go to the "docs/" subdirectory and run the `build.py` command:

```
cd docs
python build.py
```

To view the documentation you just built, open the file
`docs/build/index.html` in a web browser.


## Running tests

After you install, you can run tests to make sure everything is working. We use [pytest](https://pytest.org/) for testing and generating code coverage reports.  The `pytest` command should be available in the conda environment created by running the `install.sh` script as described in the installation instructions.

To run tests against the core modules, unit models and DMF, and generate a coverage report, run tests in `idaes/` with the following command:

```
$ pytest  # Please note some tests may be skipped based on solver availability. 
```

If there are errors, or you are having trouble, you can use our [issue tracker on Github](https://github.com/IDAES/idaes/issues) to look for other users experiencing similar problems, or to report a new bug.


## Running a Jupyter notebook

There are example [Jupyter](https://jupyter.org) notebook(s) in the `examples/` 
directory. To access them, you should start up a Jupyter Lab notebook server using the
following command.

```
jupyter lab
```

This should start up a server and then pop up a tab or window in your default 
web browser showing the Jupyter UI. On the left hand side you can browse to 
available notebooks (files ending in ".ipynb"). For more information on how to 
use Jupyter Lab, use the built-in *Help* menu and the extensive documentation 
on the [Jupyter website](https://jupyter.org).
For more details on the examples, please refer to the 
[online documentation](https://idaes-pse.readthedocs.io/en/latest/). 

## Contacts and more information

General, background and overview information is available at the [IDAES main website](https://www.idaes.org).
Framework development happens at our [GitHub repo](https://github.com/IDAES/idaes-pse) where you can [report issues/bugs](https://github.com/IDAES/idaes-pse/issues) or [make contributions](https://github.com/IDAES/idaes-pse/pulls).
For further enquiries, send an email to: <idaes-support@idaes.org>

