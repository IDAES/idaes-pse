# Thermodynamics-and-kinetics

This repository contains ALAMO and RIPE wrappers for the
IDAES project.

## Contributing

**By contributing to this repository, you are agreeing to all the terms set out
in the LICENSE.txt and COPYRIGHT.txt files in this directory.**

## ALAMOPy: Python interface to ALAMO

The ALAMO Python module facilitates easy access to ALAMO and provides additional tools for uncertainty analysis. Detailed documentation of ALAMO Python can be found in [ALAMOPy documentation](./documentation.pdf). The tools in the package require ALAMO and BARON executable files, availible at www.minlp.com

### Installation

Change to the `alamo_python` subdirectory and use the standard
Python install formula:

```
python setup.py install
```


### Usage

The ALAMOpy wrapper can be used to call alamo from Python. ALAMO should be located in the user's path, or the location of the ALAMO executable can be provided to the wrapper as a key word arguement. 

ALAMO can be called through the wrapper using the doalamo subroutine

alamopy.doalamo(xdata,zdata)

### Results

By default, two files will be generated:
* logscratch - A file containing the output ALAMO provide to the terminal
* tempalm.py - A Python usable file named tempalm.py by default, this file contains the resultant ALAMO model. The key word argument almname can be used to control the name of this file.

The python file produced by ALAMOpy can be used by importing it into Python.

import tempalm
z = tempalm.f(x)

The resultant function defined in tempalm.py can be easily incorporated into PYOMO models, as demonstrated in Jupyter notebook ./sixhumpcamel.ipynb included with the examples.

#### Additional results and options

res = alamopy.doalamo(xdata,zdata)

If the results of the call are assigned to a varaible, additional information is provided in a dictionary including:
* res['model']    - a string representation of the resultant ALAMO model
* res['f(model)'] - a lambda function of the resultant model
* res['ssr']      - the sum of squared residuals on the provided training set

The ALAMO wrapper allows for any options defined in the ALAMO user manual to be appended directly to any constructed .alm file by using the keyword argument almopt='/optfile'.

The additional keyword arguemnts may be provided directly in the call to alamopy.doalamo including:
* All basis fucntion transformations availible through ALAMO (monomialpower, multi2power, etc... )
* xlabels = list['str', ...] - labels given to input variables as a list of strings
* zlabels = list['str', ...] - labels given to outputs
* xval = 1/2d list or np.array - validaiton data for alamo
* zval = 1/2d list or np.array - response validation data for alamo
* modeler = int - modeler value used in alamo
* savescratch = T/F - saves .alm and .lst
* savetrace = T/F - saves trace file
* expandoutput = T/F - provides additional information in the output dictionary
* cvfun = T/F - additional file with variable coefficients useful for cross validation
* almname = 'str' - specify a name for the .alm file and resultant .py model file
* almpath = 'str' - path to the alamo executable
* simulator = 'str' - path to a simulator or python function to interface with ALAMO


The additional output information provided by the expandoutput option is also contained in the associated ALAMO .trc and .lst files. A thorough explanation can be obtained in the ALAMO manual.


## ALAMOPy Examples

Three examples are included with the current version of alamopy:
* [camel6.py](./camel6.py) with a [Jupyter notebook](./sixhumpcamel.ipynb)
* [ackley.py](./ackley.py)
* [branin.py](./branin.py)

## RIPE : Reaction Identification and Parameter Estimation

The RIPE module provides tools for reaction network identification. RIPE uses reactor data consisting of concentration, or conversion, values for multiple species that are obtained dynamically, or at multiple process conditions (temperatures, flow rates, working volumes) to identify probable reaction kinetics. The RIPE module also contains tools to facilitate adaptive experimental design. The experimental design tools in RIPE require the use of the python package RBFopt. More information for RBFopt is availible at www.github.com/coin-or/rbfopt

### Installation

Change to the `ripe_python` subdirectory and use the standard
Python install formula:

    python setup.py install

### Usage

RIPE can be used to build models for static datasets through the function ripe.ripemodel

ripe_results = ripe.ripemodel(data, kwargs)

* data is provided to RIPE as one, two, or three dimensional python data structures, where the first axis corresponds to observations at different process conditions, the second axis corresponds to observations of different chemical species, and the third axis corresponds to dynamic observation of a chemical species at a specified process condition.

RIPE adaptive experimental design can be accessed using ripe.ems

[proposed_x, errors] = ripe.ems(ripe_results, simulator, l_bounds, u_bounds, n_species, kwargs)

* ripe_results - The results from ripe.ripemodel, additional information provided in the results section
* simulator - a black-box simulator for the unknown process.
* l_bounds/u_bounds - lower andupper bounds for the input variables in the adaptive design
* nspecies - the number of chemical species present in the black-box system

Reaction stoichiometries and mechanisms are provided explicitly to ripemodel through the keyword arguments mechanisms and stoichiometry. Detailed explanations of the forms of these arguments are provided in the stoiciometry and mechanism specification section. Additional keyword arguments can be found in the additional options section.

### Results

By default, one file will be generated
* riperesults.txt - a file containing the selected reactions and parameter estimates

### Stiochiometry and Mechanism Specification

Considered reaction stiochiometries are provided through keyword arguments.

#### Stoichiometry

Considered reaction stoichiometries are defiend as a list of list, where reactants and products are defined as negative and positive integers , respectively, according to their stoichiometric coefficeints. A set of considered reaction stoichiometries must be provided. If process data consists of species conversion, a positive coefficient should be specified.

#### Mechanisms

Considered reaction mechanisms are provided explicitly to RIPE through q keyword argument. If no kinetic mechanisms are specified, mass action kinetics are ascribed to every considered stoichiometry. RIPE contains kinetic mechanisms defined internally, and called through ripe.mechs.<mechanism>. The availible mechanisms include:

* massact - mass action kinetics, order informed by reaction stoichiometry

19 empirical rate forms included relate specifically to catalyst conversion in chemical looping combustion reactors include:

* Random nucleation
* Power law models
* Avrami-Erofeev models

These internal kinetics can be specified by calling ripe.mechs.massact or ripe.mechs.clcforms respectively. User-defined kinetic mechanisms can also be supplied to RIPE as python functions. An example is provided in the file crac.py.

### Additional Results and Options

In addition to the arguments stoichiometry and mechanism, a number of other optional arguments are availible, including:

Arguments relating to process conditions

* x0 - initial concentration at each process condition for every species
* time - time associated with dynamic samples for every process condition
* temp - temperature associated with every process condition
* flow - flow rate at every process condition for every species
* vol - reactor volume at every process condition

Arguments related to RIPE algorithmic function

* tref - reference termpeature for reformulated Arrhenius models
* ccon - specified cardinality constraint instead of BIC objective
* sigma - expected variance of noise, estimated if not provided
* onemechper - one mechanism per stoichiometry in selected model, true by default 

Additional arguments

* minlp_path - path to baron or other minlp solver, can also be set in shared.py
* alamo_path - path to alamo, can also be set in shared.py
* expand_output - provide estimates for noise variance in model resutls
* zscale - linear scaling of observed responses between -1 and 1
* ascale - linear scaling of activities between -1 and 1
* hide_output - surpress output to terminal
* keepfiles - keep scratch files for debugging
* showpyomo - show pyomo output to terminal, false by default

### RIPE Examples

Three examples are included with RIPE. These examples demonstrate different use cases, and provide a template for utilizing user-defined mechanisms. 

* clc.py - a chemical looping combustion example in which catalyst conversion is observed over time
* isoT.py - an example that utilizes both ripe.ripemodel and ripe.ems
* crac.py - an example that utilizes user-defined reaction mechanisms

All of these examples are built for Linux machines. They can be called from the command line by calling python directly, or can be called from inside a python environment using execfile().
