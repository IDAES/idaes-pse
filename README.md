# IDAES Toolkit

The IDAES Toolkit aims to provide multi-scale, simulation-based, open source
computational tools and models to support the design, analysis, optimization,
scale-up, operation and troubleshooting of innovative, advanced energy systems.

<!-- BEGIN Status badges -->
## Build statuses
![Tests](https://github.com/IDAES/idaes-pse/workflows/Tests/badge.svg?branch=main)
![Integration](https://github.com/IDAES/idaes-pse/workflows/Integration/badge.svg?branch=main)
[![codecov](https://codecov.io/gh/IDAES/idaes-pse/branch/main/graph/badge.svg?token=1lNQNbSB29)](https://codecov.io/gh/IDAES/idaes-pse)
[![Documentation Status](https://readthedocs.org/projects/idaes-pse/badge/?version=latest)](https://idaes-pse.readthedocs.io/en/latest/?badge=latest)
[![GitHub contributors](https://img.shields.io/github/contributors/IDAES/idaes-pse.svg)](https://github.com/IDAES/idaes-pse/graphs/contributors)
[![Merged PRs](https://img.shields.io/github/issues-pr-closed-raw/IDAES/idaes-pse.svg?label=merged+PRs)](https://github.com/IDAES/idaes-pse/pulls?q=is:pr+is:merged)
[![Issue stats](http://isitmaintained.com/badge/resolution/IDAES/idaes-pse.svg)](http://isitmaintained.com/project/IDAES/idaes-pse)
<!-- END Status badges -->

## Getting Started

Our [complete documentation is online](https://idaes-pse.readthedocs.io/en/stable/) but here is a summarized set of steps to get started using the framework.

While not required, we encourage the installation of [Anaconda](https://www.anaconda.com/products/individual#Downloads) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) and using the `conda` command to create a separate python environment in which to install the IDAES Toolkit.

Use conda to create a new "idaes-pse" (could be any name you like) environment then activate that environment:
```bash
conda create --name idaes-pse python=3
conda activate idaes-pse
```

Now, in that "idaes-pse" environment, install the IDAES Toolkit using either `pip install` or `conda install` (but not both):

```bash
# install latest stable release
pip install idaes_pse
# install latest version from the main branch of this repository
pip install 'idaes-pse[prerelease] @ https://github.com/IDAES/idaes-pse/archive/main.zip'
```

You can check the version installed with the command:

```bash
idaes --version
```

Now install the examples and the pre-build extensions (binary solvers):

```bash
idaes get-examples
idaes get-extensions  # on MacOS use: conda install -c conda-forge ipopt
```

This will install the examples into an `examples` subdirectory which can be opened using a [Jypter](https://jupyter.org) Notebook:

```bash
jupyter notebook examples/notebook_index.ipynb
```
From there you can explore the examples and tutorials.

For more information on how to use Jupyter Lab, use the built-in *Help* menu and the extensive documentation on the [Jupyter website](https://jupyter.org).

Finally, refer to the [complete idaes-pse documentation](https://idaes-pse.readthedocs.io/en/stable) for more detailed [installation instructions](https://idaes-pse.readthedocs.io/en/stable/getting_started/), [user guide](https://idaes-pse.readthedocs.io/en/stable/user_guide/), examples, technical specification, etc.

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
* Python 3.8
* Python 3.9

Note that Python 3.5 is *not* supported.

## Contacts and more information

General, background and overview information is available at the [IDAES main website](https://www.idaes.org).
Framework development happens at our [GitHub repo](https://github.com/IDAES/idaes-pse) where you can ask questions by starting a [discussion](https://github.com/IDAES/idaes-pse/discussions), [report issues/bugs](https://github.com/IDAES/idaes-pse/issues) or [make contributions](https://github.com/IDAES/idaes-pse/pulls).
For further enquiries, send an email to: <idaes-support@idaes.org>

## Contributing

Please see our [Advanced User Guide](https://idaes-pse.readthedocs.io/en/stable/advanced_user_guide/) and [Developer Documentation](https://idaes-pse.readthedocs.io/en/stable/advanced_user_guide/developer/) on how to work with the idaes-pse source code and contirbute changes to the project.

**By contributing to this repository, you are agreeing to all the terms set out in the LICENSE.md and COPYRIGHT.md files in this directory.**
