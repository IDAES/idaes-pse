.. alamopy documentation master file, created by
   sphinx-quickstart on Wed Mar 21 17:35:59 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

ALAMO Python (alamopy)
======================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Basic Usage
-----------

    import alamopy
    results = alamopy.doalamo(xdata, zdata)
    # OR:
    # results = alamopy.doalamo(xdata, zdata, xval=xv, zval=z)
    confidence = alamopy.almconfidence(results, __xdata*__, __zdata*__)
    alamopy.almpick(results, 'function_name')
    results = alamopy.almunpickle('function_name')

Options available for *doalamo*
-------------------------------
Possible arguments to be passed to ALAMO through do alamo and additional arguments that govern the behavior of doalamo.

* xlabels - list of strings to label the input variables
* zlabels - list of strings to label the output variables
* functions - logfcns, expfcns, cosfcns, sinfcns, linfcns, intercept
  * These are '0-1' options to activate these functions
* monomialpower, multi2power, multi3power, ratiopower
  * List of terms to be used in the respective basis functions
* modeler - integer 1-7 determines the choice of fitness metrice
* solvemip - '0-1' option that will force the solving of the .gms file


Options specific to *doalamo*
-----------------------------
These options are specific to alamopy and will not change the behavior of the underlying .alm file.

* expandoutput - '0-1' option that can be used to collect more information from the ALAMO .lst and .trc file
* showalm - '0-1' option that controlif the ALAMO output is printed to screen
* almname - A string that will assign the name of the .alm file
* outkeys - '0-1' option for dictionary indexing according to the output labels
* outkeys - '0-1' option for dictionary indexing according to the output labels
* outkeys - '0-1' option for dictionary indexing according to the output labels
* savetrace - '0-1' option that controls the status of the trace file
* savescratch - '0-1' option to save the .alm and .lst files
* almopt -  A string option that will append a text file of the same name to the end of each .alm fille to faciliate advanced user access in an automated fashion


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
