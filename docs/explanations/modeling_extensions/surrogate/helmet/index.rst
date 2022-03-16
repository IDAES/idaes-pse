.. alamopy documentation master file, created by
   sphinx-quickstart on Wed Mar 21 17:35:59 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. index::
    pair: alamo;alamopy


HELMET: HELMholtz Energy Thermodynamics
========================================

The purpose of HELMET (HELMholtz Energy Thermodynamics) is to provide a framework for regressing multiparameter equations of state that identify an equation for Helmholtz energy and multiple thermodynamic properties simultaneously. HELMET uses best subset selection to simultaneously model various thermodynamic properties based on the properties thermodynamic relation to Helmholtz energy. The generated model is a function of reduced density and inverse reduced temperature and uses partial derivatives to calculate the different properties. Constraints are placed on the regression to maintain thermodynamically feasible values and improve extrapolation and behavior of the model based on physical restrictions.

.. warning::
  This is the first public release of HELMET. Future work will include mixtures, regression using Pyomo models, and increased plotting and preprocessing capabilities.

Basic Usage
-----------

.. warning::
  To use this software, ALAMOPY and the solver BARON are required.

For the basic use of HELMET, the main regression steps can be imported from helmet.HELMET. These functions provide general capabilities of HELMET for new users.

.. code-block:: python

  import helmet.Helmet as Helmet

The methods available in helmet.Helmet peform the necessary steps of the regression properties.

1. **initialize(\*\*kargs)**

   Initializes key thermodynamic constants, the location of data and sampling, properties to be fit, and optimization settings

  * **molecule** - name of the chemical of interest, directs naming of files and where the data should exist
  * **fluid_data** - a tuple containing key thermodynamic constants (critical temperature, critical pressure, critical density, molecular weight, triple point, accentric factor)
  * **filename** - used for location of data
  * **gamsname** - used for naming of files
  * **max_time** - max time used for the solver
  * **props** - list of thermodynamic properties to be fit 

    Supported thermodynamic properties are 

      * Pressure: 'PVT'
      * Isochoric heat capacity: 'CV'
      * Isobaric heat capacity: 'CP'
      * Speed of Sound: 'SND'

  * **sample** - sample ratio, ex. sample = 3 then a third of datapoints will be used 

2. **prepareAncillaryEquations(plot=True)**

  Fits equations to saturated vapor and liquid density and vapor pressure. The keyword argument plot defaults to False

3. **viewPropertyData()**
  
  Plots the different thermodynamic properties available and a way to check that the importing of data is successful

4. **setupRegression(numTerms = 12, gams=True)**

  Writes the optimization program for modelling the thermodynamic properties. Currently this is through GAMS but in the future it can also be solved using Pyomo.

5. **runRegression()**

  Begins the modelling of the multiparameter equation

6. **viewResults(filename)**

  Based on the optimization settings, the solution of the regression is parsed and fitness metrics are calculated. The results can be visualized with different plots.




HELMET Output
-----------------

The output for HELMET is a single equation representing Helmholtz energy. Partial derivatives of this equation will give you the fit thermodynamic properties as well as other properties related to Helmholtz energy.


HELMET Examples
----------------

The provided HELMET example uses data modified for this application and made available by the IAPWS orgnization at http://www.iapws.org/95data.html for IAPWS Formulation 1995 for Thermodynamic Properties of Ordinary Water Substance for General and Scientific Use.


..   import helmet.Helmet as Helmet

..   num_terms = 14
..   max_time = 1000

..   Fluids = {'H2O': (647.096, 22.064, 17.8737279956, 18.015268, 273.16, 0.344 )}

..   molecule = 'H2O'
..   (critT, critP, critD, M, triple, acc) = Fluids[molecule]
..   R = 8.314472; # J mol^-1 K^-1 

..   # Constants for a molecule 
..   Helmet.initialize(molecule=molecule, 
..                     fluid_data = Fluids[molecule], 
..                     filename = os.getcwd() + "/%s"%molecule, 
..                     gamsname= os.getcwd() + "/%s"%molecule, 
..                     max_time = max_time, 
..                     props=['PVT','CV', 'CP','SND'], 
..                     sample =3)
   
..   # Prepare Ancillary Equations of sat liq/vapor density and vapor pressure
..   Helmet.prepareAncillaryEquations(plot = True)  # plot=True

..   # View data used for regression
..   Helmet.viewPropertyData()

..   # Write and runs GAMS data file and regression file
..   Helmet.setupRegression(numTerms = 14, gams=True)
..   Helmet.runRegression(gams=True)
          

..   # View Results by importing the data
..   Helmet.viewResults("H2Omain.lst")

