.. _pecos:

Pecos: Data Quality Control and Fault Detection
===============================================

Before using plant data in process models, quality control and fault detection analysis is recommended to identify 
potential data issues (e.g., missing or corrupt data) and data points that are not suitable for the intended analysis 
(e.g., abnormal plant behavior).
The following documentation describes methods to run data quality control and fault detection analysis using Pecos.

Pecos is an open-source Python package designed to monitor performance of time series data, subject to a series of quality control tests. 
The software includes methods to run quality control tests defined by the user and generate reports which include 
test results and graphics. Results from the quality control analysis can be used to extract "clean data" 
which removes data points that failed quality control inspection.
Pecos was originally developed for the U.S. Department of Energy in 2016 to monitor solar photovoltaic systems and has been used for a wide 
range of applications since.  The software was updated for IDAES to facilitate near real-time analysis using continuous data streams.  

More information on Pecos can be found in the online documentation at https://pecos.readthedocs.io.
The software repository is located at https://github.com/sandialabs/pecos.

The following functionality is available in Pecos:

* Check data for missing, non-monotonic, and duplicate time stamps
* Check for missing data
* Check for corrupt data
* Check for data that are outside the expected range
* Check for stagnant data and/or abrupt changes in the data using the difference between max and min values within a rolling window
* Check for outliers using normalized data within a rolling window

The analysis generates the following information:

* Cleaned data (data that failed a test are removed)
* Boolean mask (indicates which data points failed a test)
* Summary of the quality control test results (includes the variable name, start and end time for each failure, and an error message)

The test results summary and accompanying graphics can be included in HTML or LATEX reports generated using Pecos.

Pecos supports both static and streaming analysis along with custom quality control functions:
 
* Static analysis operates on the entire data set to determine if all data points are normal or anomalous. 
  While this can include operations like moving window statistics, the quality control tests operate on the entire data set at once. 
* Streaming analysis loops through each data point using a quality control test that relies on information from "clean data" in a moving window. 
  If a data point is determined to be anomalous, it is not included in the window for subsequent analysis. 
* The user can define custom quality control functions used to determine if data is anomalous and return custom metadata from the analysis.

Data points that do not pass quality control inspection should be
removed or replaced by various means (interpolation, data from a duplicate sensor, values from a model) before using the data for further analysis.
Data replacement strategies are generally defined on a case-by-case basis. 
If large sections of the data failed quality control tests, the data might not be suitable for use.

The raw data, results from the quality control analysis, and the analysis files used to run Pecos can be stored in the 
Data Management Framework (DMF) to ensure reproducible results.

Installation
------------

To install Pecos using pip::

	pip install pecos 
	
To install Pecos using git::

	git clone https://github.com/sandialabs/pecos
	cd pecos
	python setup.py install
	
Examples
--------
IDAES examples that use Pecos are listed on the |examples-site|. 
Pecos also includes several general examples, located at https://github.com/sandialabs/pecos/tree/master/examples.
