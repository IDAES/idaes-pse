# IDAES Toolkit

The IDAES Toolkit aims to provide multi-scale, simulation-based, open source
computational tools and models to support the design, analysis, optimization,
scale-up, operation and troubleshooting of innovative, advanced energy systems.

<!-- BEGIN Status badges -->
## Project Build and Download Statuses
[![Tests](https://github.com/IDAES/idaes-pse/actions/workflows/core.yml/badge.svg)](https://github.com/IDAES/idaes-pse/actions/workflows/core.yml)
[![Integration](https://github.com/IDAES/idaes-pse/actions/workflows/integration.yml/badge.svg)](https://github.com/IDAES/idaes-pse/actions/workflows/integration.yml)
[![codecov](https://codecov.io/gh/IDAES/idaes-pse/branch/main/graph/badge.svg?token=1lNQNbSB29)](https://codecov.io/gh/IDAES/idaes-pse)
[![Documentation Status](https://readthedocs.org/projects/idaes-pse/badge/?version=latest)](https://idaes-pse.readthedocs.io/en/latest/?badge=latest)
[![Services](https://github.com/Pyomo/jenkins-status/blob/main/idaes_services.svg)](https://pyomo-jenkins.sandia.gov/)
[![GitHub contributors](https://img.shields.io/github/contributors/IDAES/idaes-pse.svg)](https://github.com/IDAES/idaes-pse/graphs/contributors)
[![Merged PRs](https://img.shields.io/github/issues-pr-closed-raw/IDAES/idaes-pse.svg?label=merged+PRs)](https://github.com/IDAES/idaes-pse/pulls?q=is:pr+is:merged)
[![Issue stats](http://isitmaintained.com/badge/resolution/IDAES/idaes-pse.svg)](http://isitmaintained.com/project/IDAES/idaes-pse)
[![Downloads](https://pepy.tech/badge/idaes-pse)](https://pepy.tech/project/idaes-pse)
<!-- END Status badges -->

## Getting Started

Our [complete documentation is online](https://idaes-pse.readthedocs.io/en/stable/) but here is a summarized set of steps to get started using the framework.

While not required, we encourage the installation of [Anaconda](https://www.anaconda.com/products/individual#Downloads) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) and using the `conda` command to create a separate python environment in which to install the IDAES Toolkit.

Use conda to create a new "idaes-pse" (could be any name you like) environment then activate that environment:
```bash
conda create --name idaes-pse python=3.10
conda activate idaes-pse
```

Now, in that "idaes-pse" environment, install the IDAES Toolkit using either `pip install` or `conda install` (but not both):

```bash
# install latest stable release
pip install idaes-pse
# install latest stable release with one set of optional dependencies, e.g. `ui` for the user interface
pip install "idaes-pse[ui]"
# install latest stable release with multiple sets of optional dependencies
pip install "idaes-pse[ui,dmf,omlt,grid,coolprop]"
# install latest version from the main branch of this repository
pip install "idaes-pse @ git+https://github.com/IDAES/idaes-pse@main"
# install from the `mybranch` branch of the fork belonging to `myuser`
pip install "idaes-pse @ git+https://github.com/myuser/idaes-pse@mybranch"
```

You can check the version installed with the command:

```bash
idaes --version
```

Now install the pre-built extensions (binary solvers):

```bash
idaes get-extensions
```

The IDAES examples can be installed by running:

```bash
pip install idaes-examples
```

For more information, refer to the [IDAES/examples](https://github.com/IDAES/examples) repository, as well as the online static version of the examples available at <https://idaes-examples.readthedocs.org>.

Finally, refer to the [complete idaes-pse documentation](https://idaes-pse.readthedocs.io/en/latest) for detailed [installation instructions](https://idaes-pse.readthedocs.io/en/latest/tutorials/getting_started/index.html), examples, guides, and reference.

## System requirements

The code and examples have been tested with the following operating systems:

|Operating system|Supported versions  |
|----------------|--------------------|
| Linux          | Any modern Linux   |
| Windows        | Windows 10         |
| macOS          | Partly supported*  |

*HSL is not currently provided for macOS on Intel processors, so some features may be limited or not available.

Most of the functionality is implemented in Python. In accordance with
the end-of-life for many Python 2 libraries, the IDAES Toolkit is written
for Python 3. The following sub-versions are supported:

* Python 3.8
* Python 3.9
* Python 3.10
* Python 3.11

Note that Python 3.6 is *not* supported.

## Contacts and more information

General, background and overview information is available at the [IDAES main website](https://www.idaes.org).
Framework development happens at our [GitHub repo](https://github.com/IDAES/idaes-pse) where you can ask questions by starting a [discussion](https://github.com/IDAES/idaes-pse/discussions), [report issues/bugs](https://github.com/IDAES/idaes-pse/issues) or [make contributions](https://github.com/IDAES/idaes-pse/pulls).
For further enquiries, send an email to: <idaes-support@idaes.org>

## Funding acknowledgements

This work was conducted as part of the [Institute for the Design of Advanced Energy Systems (IDAES)](https://idaes.org)
with support through the [Simulation-Based Engineering, Crosscutting Research Program](https://netl.doe.gov/coal/simulation-based-engineering)
within the U.S. Department of Energy’s [Office of Fossil Energy and Carbon Management (FECM)](https://www.energy.gov/fecm/office-fossil-energy-and-carbon-management).
As of 2021, additional support was provided by FECM’s [Solid Oxide Fuel Cell Program](https://www.energy.gov/fecm/science-innovation/clean-coal-research/solid-oxide-fuel-cells),
and [Transformative Power Generation Program](https://www.energy.gov/fecm/science-innovation/office-clean-coal-and-carbon-management/advanced-energy-systems/transformative).

## Contributing

Please see our [Advanced User Installation](https://idaes-pse.readthedocs.io/en/stable/tutorials/advanced_install/) and [How-to Guides](https://idaes-pse.readthedocs.io/en/stable/how_to_guides/) on how to work with the idaes-pse source code and contribute changes to the project.

**By contributing to this repository, you are agreeing to all the terms set out in the LICENSE.md and COPYRIGHT.md files in this directory.**
