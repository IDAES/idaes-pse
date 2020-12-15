.. _qc:

Data Quality Control and Fault Detection
============================================

Before using plant data in process models, quality control and fault detection analysis is recommended to identify 
potential data issues (i.e., missing or corrupt data) and data points that are not suitable for the intended analysis (i.e., abnormal plant behavior).
This workflow describes methods to run data quality control and fault detection analysis using Pecos within the IDAES framework, 
available in the ``idaes.apps.pecos`` module.  

Pecos is an open-source Python package designed to monitor performance of time series data, subject to a series of quality control tests. 
The software includes methods to run quality control tests defined by the user and generate reports which include performance metrics, 
test results, and graphics. Results from the quality control analysis can be used to extract "clean data" which removes data points that failed quality control inspection.
The software can be customized for specific applications.  More information on Pecos can be found at https://pecos.readthedocs.io.

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

The test results summary and accompanying graphics can then be included in HTML or LATEX reports generated using Pecos.

Pecos supports both static and streaming analysis along with custom quality control functions:
 
* Static analysis operates on the entire data set to determine if all data points are normal or anomalous. 
  While this can include operations like moving window statistics, the quality control tests operates on the entire data set at once. 
* Streaming analysis loops through each data point using a quality control tests that relies on information from "clean data" in a moving window. 
  If a data point is determined to be anomalous, it is not included in the window for subsequent analysis. 
* The user can define custom quality control functions used to determine if data is anomalous and return custom metadata from the analysis.

Data points that do not pass quality control inspection should be
removed or replaced by various means (interpolation, data from a duplicate sensor, values from a model) before using the data for further analysis.
Data replacement strategies are generally defined on a case-by-case basis. 
If large sections of the data failed quality control tests, the data might not be suitable for use.

The raw data, results from the quality control analysis, and the analysis files used to run Pecos can be stored in the Data Management Framework (DMF) to ensure reproducibility.

Example Workflow
-----------------
Pecos includes a ``PerformanceMonitoring`` class which is the base class used to define the quality control analysis.
In the following simple example, a Pandas DataFrame (loaded from a file or extracted from the DMF) is added to a Pecos PerformanceMonitoring object 
using the ``add_dataframe`` method and quality control tests are run. 
Quality control tests can be run using the entire data set or run on specific columns of data.
The example uses the following methods: 

* ``check_timestamp`` to see if data is monotonically increasing and sampled every 10 minutes (600 seconds)
* ``check_missing`` to identify missing data
* ``check_corrupt`` to identify corrupt data values
* ``check_range`` to see if data is between -5 and 5
* ``check_delta`` to identify data that changes by less than 1 within a 20 minute moving window

.. doctest::
    :hide:

    >>> import pandas as pd
    >>> data = pd.DataFrame(data=[[4.9,1.4,-2.5], [3.5,-4.5,-4.0],[3.1,3.5,4.7],[3.0,5,9],[3.2,10.6,-999],[-4.8,0.6,15.4]], 
    ...                     index=['1/1/2020 00:00:00', '1/1/2020 00:10:00', '1/1/2020 00:30:00', '1/1/2020 00:40:00', '1/1/2020 00:20:00', '1/1/2020 00:50:00'],
    ...                     columns=['A', 'B', 'C'])
    >>> data.index = pd.to_datetime(data.index)
	
.. doctest::

    >>> print(data)
                           A     B      C
    2020-01-01 00:00:00  4.9   1.4   -2.5
    2020-01-01 00:10:00  3.5  -4.5   -4.0
    2020-01-01 00:30:00  3.1   3.5    4.7
    2020-01-01 00:40:00  3.0   5.0    9.0
    2020-01-01 00:20:00  3.2  10.6 -999.0
    2020-01-01 00:50:00 -4.8   0.6   15.4

.. doctest::

    >>> import idaes.apps.pecos
    >>> pm = idaes.apps.pecos.monitoring.PerformanceMonitoring()
    >>> pm.add_dataframe(data)
    >>> pm.check_timestamp(600)
    >>> pm.check_missing()
    >>> pm.check_corrupt([-999])
    >>> pm.check_range([-5,5])
    >>> pm.check_delta([1, None], 1200)
	
A summary of test results and cleaned data can then be extracted.  

.. doctest::

    >>> print(pm.test_results)
      Variable Name          Start Time            End Time Timesteps              Error Flag
    0               2020-01-01 00:20:00 2020-01-01 00:20:00         1  Nonmonotonic timestamp
    1             C 2020-01-01 00:20:00 2020-01-01 00:20:00         1            Corrupt data
    2             B 2020-01-01 00:20:00 2020-01-01 00:20:00         1   Data > upper bound, 5
    3             C 2020-01-01 00:40:00 2020-01-01 00:50:00         2   Data > upper bound, 5
    4             A 2020-01-01 00:10:00 2020-01-01 00:40:00         4  Delta < lower bound, 1

.. doctest::

    >>> print(pm.cleaned_data)
                           A    B    C
    2020-01-01 00:00:00  4.9  1.4 -2.5
    2020-01-01 00:10:00  NaN -4.5 -4.0
    2020-01-01 00:20:00  NaN  NaN  NaN
    2020-01-01 00:30:00  NaN  3.5  4.7
    2020-01-01 00:40:00  NaN  5.0  NaN
    2020-01-01 00:50:00 -4.8  0.6  NaN
	
Results can be included in HTML or LATEX formatted reports.  
The ``plot_test_results`` function creates a graphic for each variable that includes a quality control test failure, highlighting data points that failed a test.
The ``write_monitoring_report`` generates an report (HTML by default) that includes the test results summary and graphics.

.. doctest::
    
    >>> test_results_graphics = idaes.apps.pecos.graphics.plot_test_results(pm.data, pm.test_results)
    >>> filename = idaes.apps.pecos.io.write_monitoring_report(pm.data, pm.test_results, test_results_graphics)
    
Additional Examples
----------------------

- Pecos includes several general examples, located at https://github.com/sandialabs/pecos/tree/master/examples
- An IDAES example using plant data will be added to the examples-pse repository