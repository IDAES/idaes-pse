Data Reconciliation and Parameter Estimation
============================================

This workflow generally describes features of the IDAES framework that are useful
for data reconciliation and parameter estimation.  Many of these features can
be used for any task where plant data is to be used in conjunction with
a process model. It is assumed that the user is familiar with the IDAES modeling
platform and Pyomo. See the :ref:`General Workflow <user_guide/workflow/general:General Workflow>`
for more information on how to set up a model.

This provides general information about IDAES functionality laid out in
terms of typical use cases. See
:ref:`﻿Tutorials and Examples<tutorials_examples:Tutorials and Examples>`, for
specific complete examples of data reconciliation and parameter estimation examples.
Relevant tutorials can be found in ``tutorials/advanced/data_recon_and_parameter_estimation``.

Data Management
---------------

.. currentmodule:: idaes.dmf.model_data

IDAES has functions to read and mange process data process data. Data management
functions are contained in the ``idaes.dmf.model_data`` module.


Reading Data
~~~~~~~~~~~~
A set of process data can be stored in two files, a data file and a metadata file.
The data file is a CSV file where the first column contains a data point index,
usually a timestamp.  The first row contains a header with the measurement tag names.
The rest of the files contains measurement data. The data file format is shown in
the table below.

+----------------------+-----------+-----------+------+
| index tag (optional) | tag 1     | tag 2     | ...  |
+======================+===========+===========+======+
| index 1              | data(1,1) | data(1,2) | ...  |
+----------------------+-----------+-----------+------+
| index 2              | data(2,1) | data(2,2) | ...  |
+----------------------+-----------+-----------+------+
| ...                  | ...       | ...       | ...  |
+----------------------+-----------+-----------+------+

Metadata is provided for each tag in a separate CSV file. The metadata file has no
header row, and aside from the tag name, any column my be blank. In the metadata csv
file, the first column is the measurement tag name, the second column is a string
that maps the tag to a specific model quantity, the third column is a description of
the tag, and the fourth column is a for units of measure. Any columns past the fourth
column are ignored and can be used to store any additional information.

+--------+--------------------------+---------------+--------------------+
| tag 1  | model reference string 1 | description 1 | units of measure 1 |
+--------+--------------------------+---------------+--------------------+
| tag 2  | model reference string 2 | description 2 | units of measure 1 |
+--------+--------------------------+---------------+--------------------+
| ...    | ...                      | ...           | ...                |
+--------+--------------------------+---------------+--------------------+

The unit strings should be interpretable by `pint <https://pint.readthedocs.io/en/stable/>`_,
with the additional unit strings given in :ref:`Unit String Information <user_guide/workflow/data_rec_parmest:APPENDIX: Unit String Information>`.
The model reference string is the a string to reverence a model quantity.  In the
reference string the top-level model is always represented by ``m``.  For example,
the reference string for a heater block outlet temperature could be
``m.flowsheet.heater.control_volume.properties_out[:].temperature``. This
reference will be indexed by time.

Reading data can be done with the ``idaes.dmf.model_data.read_data()`` function.

.. autofunction:: read_data
  :noindex:

Binning Data
~~~~~~~~~~~~

Process data can be divided into bins based on some criteria.  This allows for
estimating measurement uncertainty when no better information is available, and
provides a way to look at how different process measurements vary for operating
conditions that should be similar in some way. As an example, power plant data
could be binned by power output, and assuming that operating procedures are
standard, it could be assumed that measurements in each bin should be about the
same.  The variance of a measurement in a bin could be used a an approximation of
uncertainty.  Binning the data by load and time could show how process measurements
change over time and be useful for things like fault detection and equipment
degradation.

Adding bin information to a data frame is done with the
``idaes.dmf.model_data.bin_data()`` function.

.. autofunction:: bin_data
  :noindex:

The ``idaes.dmf.model_data.bin_stdev()`` function can be used to calculate the
standard deviation for each measurement in each bin.

.. autofunction:: bin_stdev
  :noindex:

The function ``idaes.dmf.model_data.data_plot_book()`` creates a multipage PDF
containing box plots for all the measurements based on the bins.

.. autofunction:: data_plot_book
  :noindex:

To compare reconciled data to original measurements the
``idaes.dmf.data_rec_plot_book()`` is used.  For each bin there are two box
plots one for the original data and one for reconciled data.

.. autofunction:: data_rec_plot_book
  :noindex:


Tagging the Model
-----------------

Mapping process data to a model is typically done by creating model tag dictionaries.
Where the dictionary key is a measurement tag and the value is a reference to a model
variable, expression, or parameter. The tags may be process measurement tags, or any
other convenient string. IDAES has some utilities to help facilitate tagging of models.

If model reference strings where provided in the tag metadata file, the tag metadata
from the ``idaes.dmf.model_data.read_data()`` will contain model references. These
references can be accessed in the metadata dictionary as ``metadata[tag]["reference"]``.

The reference strings are optional in the tag metadata file and can be added after reading
the data. To add a reference string to the metadata, you can update the metadata dictionary
like ``metadata[tag]["reference_string"] = reference_string``. If you update the reference
string, the idaes.dmf.model_data.upadate_metadata_model_references()`` function can be used
to update the references in the tag metadata.

.. autofunction:: upadate_metadata_model_references
  :noindex:

An easy way to create a new tag dictionary from tag metadata is to use dictionary
comprehension like so: ``data_tags = {k:v["reference"][0] for k, v in metadata.items() if v["reference"] is not None}``

Often data reconciliation is performed using process data as the first step of
parameter estimation or model validation. Data reconciliation can often fill in
information for many unmeasured quantities.  Most of this data can be associated
with process streams.  To make managing proliferation of data tags easier,
it is often desirable to create a new set of tags based on stream (Arc) names that
can be automatically obtained from the model.

.. currentmodule:: idaes.core.util.tables

The first step to creating a new set of tags based on streams is to get a dictionary
of streams and their associated state blocks, with can be done with the
``idaes.core.util.tables.arcs_to_stream_dict()`` function.

.. autofunction:: arcs_to_stream_dict
  :noindex:

The stream dictionary can be converted to a corresponding dictionary of state at a
specific time with the ``idaes.core.util.tables.stream_states_dict()`` function.

.. autofunction:: stream_states_dict
  :noindex:

With the dictionary of states, a tag dictionary can be created automatically with
the ``idaes.core.util.tables.stream_states_dict()`` function.

.. autofunction:: stream_states_dict
  :noindex:


Objective Function
------------------

For either parameter estimation or data reconciliation the objective function is
often written in the form:


.. math::

  \min \sum_i \frac{(x_{\text{data},i} - x_{\text{model},i})^2}{\sigma_i^2}

To add the data to be used in the objective standard practice has been to add a
mutable parameter for data and standard deviation indexed by measurement tags.  The
parameter values can be set from the measurement data frame to a specific index.

The following code snippet exemplifies the use of data parameters in a model.

.. code-block:: python

  # df is from reading process data into a pandas.DataFrame. bin_stdev comes
  # from binning the data and calculating the standard deviations, as described
  # in the Data Management section.

  # Add data parameters
  m.data = pyo.Param(data_tags, mutable=True)
  m.data_stdev = pyo.Param(data_tags, mutable=True)

  # A function to set the data parameters from measurement data
  def set_data(m, df, data_tags, index=None):
    m.bin_no = df.iloc[index]["bin_no"]
    for t in data_tags:
        m.data[t] = df.iloc[index][t]
        m.data_stdev[t] = bin_stdev[m.bin_no][t]

  # Expressions for error in the objective function
  @m.Expression(data_tags)
  def err(m, i):
      return (m.data[i] - data_tags[i])/m.data_stdev[i]


Parameter Estimation
--------------------

For parameter estimation and data reconciliation, it is recommended to use
`Paramest <https://pyomo.readthedocs.io/en/stable/contributed_packages/parmest/>`_.
If more control is needed a user can also set up their own parameter estimation
problem, by combining multiple process models into a larger parameter estimation
model.


APPENDIX: Unit String Information
---------------------------------

This section provides additional detail about units strings that can be used to read
data with the ``read_data()`` function.

Temperature Differences
~~~~~~~~~~~~~~~~~~~~~~~

The unit conversion for temperatures with offsets (C and F) are deferent depending
on whether a measurement is temperature or temperature difference.  It is important
to ensure the temperature units are correctly specified before reading data. The units
"delta_degC and delta_degF" are defined to handle temperature differences.

Units Not Converted
~~~~~~~~~~~~~~~~~~~

The following units are not affected by unit conversion.

  - percent
  - ppm
  - ppb
  - pH
  - VAR
  - MVAR
  - H2O
  - percent open
  - percent closed

Gauge Pressure
~~~~~~~~~~~~~~

The table below shows a list of unit string that are taken to be gauge pressure
when data is read. Gauge pressures are converted to absolute pressure in the unit
conversion process.

===================   ======================
Gauge Pressure Unit   Absolute Pressure Unit
===================   ======================
psig                  psi
in water gauge        in water
in hg gauge           in hg
===================   ======================


Additional Unit Definitions
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Some units are common enough in process data that the ``read_data()`` function will
recognize them and convert them to standard Pint unit strings. Unit strings are case
sensitive to handle things like milli (m) and mega (M) prefixes.

The table below shows the additional units.

.. note:
    This is not a complete list of units, just strings that are recognized in addition
    to what is recognized by Pint. Additional unit definitions can be passed to the
    ``read_data()`` function if required.

+---------------------+-------------------------+
| Unit String         |  Pint Unit String       |
+=====================+=========================+
| Pressure                                      |
+---------------------+-------------------------+
| PSI                 | psi                     |
+---------------------+-------------------------+
| PSIA                | psi                     |
+---------------------+-------------------------+
| psia                | psi                     |
+---------------------+-------------------------+
| PSIG                | psig                    |
+---------------------+-------------------------+
| INWC                | in water                |
+---------------------+-------------------------+
| IN WC               | in water                |
+---------------------+-------------------------+
| IN/WC               | in water                |
+---------------------+-------------------------+
| " H2O               | in water                |
+---------------------+-------------------------+
| INHG                | in hg                   |
+---------------------+-------------------------+
| IN HG               | in hg                   |
+---------------------+-------------------------+
| IN/HG               | in hg                   |
+---------------------+-------------------------+
| HGA                 | in hg                   |
+---------------------+-------------------------+
| IN HGA              | in hg                   |
+---------------------+-------------------------+
| Fraction                                      |
+---------------------+-------------------------+
| PCT                 | percent                 |
+---------------------+-------------------------+
| pct                 | percent                 |
+---------------------+-------------------------+
| PERCT               | percent                 |
+---------------------+-------------------------+
| PERCT.              | percent                 |
+---------------------+-------------------------+
| PCNT                | percent                 |
+---------------------+-------------------------+
| PPM                 | ppm                     |
+---------------------+-------------------------+
| PPB                 | ppb                     |
+---------------------+-------------------------+
| % OPEN              | percent open            |
+---------------------+-------------------------+
| % CLSD              | percent closed          |
+---------------------+-------------------------+
| % CLOSED            | percent closed          |
+---------------------+-------------------------+
| Length                                        |
+---------------------+-------------------------+
| IN                  | in                      |
+---------------------+-------------------------+
| INS                 | in                      |
+---------------------+-------------------------+
| INCHES              | in                      |
+---------------------+-------------------------+
| Inches              | in                      |
+---------------------+-------------------------+
| FT                  | ft                      |
+---------------------+-------------------------+
| FEET                | ft                      |
+---------------------+-------------------------+
| FOOT                | ft                      |
+---------------------+-------------------------+
| Feet                | ft                      |
+---------------------+-------------------------+
| MILS                | minch                   |
+---------------------+-------------------------+
| Speed                                         |
+---------------------+-------------------------+
| MPH                 | mile/hr                 |
+---------------------+-------------------------+
| IPS                 | in/s                    |
+---------------------+-------------------------+
| Volume                                        |
+---------------------+-------------------------+
| KGAL                | kgal                    |
+---------------------+-------------------------+
| Vol Flow                                      |
+---------------------+-------------------------+
| GPM                 | gal/min                 |
+---------------------+-------------------------+
| gpm                 | gal/min                 |
+---------------------+-------------------------+
| CFM                 | ft^3/min                |
+---------------------+-------------------------+
| KCFM                | ft^3/mmin               |
+---------------------+-------------------------+
| SCFM                | ft^3/min                |
+---------------------+-------------------------+
| KSCFM               | ft^3/mmin               |
+---------------------+-------------------------+
| Angle                                         |
+---------------------+-------------------------+
| DEG                 | deg                     |
+---------------------+-------------------------+
| Angular Speed                                 |
+---------------------+-------------------------+
| RPM                 | rpm                     |
+---------------------+-------------------------+
| Frequency                                     |
+---------------------+-------------------------+
| HZ                  | hz                      |
+---------------------+-------------------------+
| Temperature                                   |
+---------------------+-------------------------+
| DEG F               | degF                    |
+---------------------+-------------------------+
| Deg F               | degF                    |
+---------------------+-------------------------+
| deg F               | degF                    |
+---------------------+-------------------------+
| DEG C               | degC                    |
+---------------------+-------------------------+
| Deg C               | degC                    |
+---------------------+-------------------------+
| deg C               | degC                    |
+---------------------+-------------------------+
| DEGF                | degF                    |
+---------------------+-------------------------+
| DegF                | degF                    |
+---------------------+-------------------------+
| DEGC                | degC                    |
+---------------------+-------------------------+
| DegC                | degC                    |
+---------------------+-------------------------+
| Temperature Difference                        |
+---------------------+-------------------------+
| DELTA DEG F         | delta_degF              |
+---------------------+-------------------------+
| DETLA Deg F         | delta_degF              |
+---------------------+-------------------------+
| DETLA deg F         | delta_degF              |
+---------------------+-------------------------+
| DETLA DEG C         | delta_degC              |
+---------------------+-------------------------+
| DETLA Deg C         | delta_degC              |
+---------------------+-------------------------+
| DELTA deg C         | delta_degC              |
+---------------------+-------------------------+
| DELTA DEGF          | delta_degF              |
+---------------------+-------------------------+
| DELTA DegF          | delta_degF              |
+---------------------+-------------------------+
| DELTA degF          | delta_degF              |
+---------------------+-------------------------+
| DELTA DEGC          | delta_degC              |
+---------------------+-------------------------+
| DELTA DegC          | delta_degC              |
+---------------------+-------------------------+
| DELTA degC          | delta_degC              |
+---------------------+-------------------------+
| Delta DEG F         | delta_degF              |
+---------------------+-------------------------+
| Delta Deg F         | delta_degF              |
+---------------------+-------------------------+
| Delta deg F         | delta_degF              |
+---------------------+-------------------------+
| Delta DEG C         | delta_degC              |
+---------------------+-------------------------+
| Delta Deg C         | delta_degC              |
+---------------------+-------------------------+
| Delta deg C         | delta_degC              |
+---------------------+-------------------------+
| Delta DEGF          | delta_degF              |
+---------------------+-------------------------+
| Delta DegF          | delta_degF              |
+---------------------+-------------------------+
| Delta degF          | delta_degF              |
+---------------------+-------------------------+
| Delta DEGC          | delta_degC              |
+---------------------+-------------------------+
| Delta DegC          | delta_degC              |
+---------------------+-------------------------+
| Delta degC          | delta_degC              |
+---------------------+-------------------------+
| delta DEG F         | delta_degF              |
+---------------------+-------------------------+
| delta Deg F         | delta_degF              |
+---------------------+-------------------------+
| delta deg F         | delta_degF              |
+---------------------+-------------------------+
| delta DEG C         | delta_degC              |
+---------------------+-------------------------+
| delta Deg C         | delta_degC              |
+---------------------+-------------------------+
| delta deg C         | delta_degC              |
+---------------------+-------------------------+
| delta DEGF          | delta_degF              |
+---------------------+-------------------------+
| delta DegF          | delta_degF              |
+---------------------+-------------------------+
| delta degF          | delta_degF              |
+---------------------+-------------------------+
| delta DEGC          | delta_degC              |
+---------------------+-------------------------+
| delta DegC          | delta_degC              |
+---------------------+-------------------------+
| delta degC          | delta_degC              |
+---------------------+-------------------------+
| Energy                                        |
+---------------------+-------------------------+
| MBTU                | kbtu                    |
+---------------------+-------------------------+
| Mass                                          |
+---------------------+-------------------------+
| MLB                 | klb                     |
+---------------------+-------------------------+
| K LB                | klb                     |
+---------------------+-------------------------+
| K LBS               | klb                     |
+---------------------+-------------------------+
| lb.                 | lb                      |
+---------------------+-------------------------+
| Mass flow           |                         |
+---------------------+-------------------------+
| TPH                 | ton/hr                  |
+---------------------+-------------------------+
| tph                 | ton/hr                  |
+---------------------+-------------------------+
| KLB/HR              | klb/hr                  |
+---------------------+-------------------------+
| KPPH                | klb/hr                  |
+---------------------+-------------------------+
| Current                                       |
+---------------------+-------------------------+
| AMP                 | amp                     |
+---------------------+-------------------------+
| AMPS                | amp                     |
+---------------------+-------------------------+
| Amps                | amp                     |
+---------------------+-------------------------+
| Amp                 | amp                     |
+---------------------+-------------------------+
| AMP AC              | amp                     |
+---------------------+-------------------------+
| pH                                            |
+---------------------+-------------------------+
| PH                  | pH                      |
+---------------------+-------------------------+
| VARS (volt-amp reactive)                      |
+---------------------+-------------------------+
| VARS                | VAR                     |
+---------------------+-------------------------+
| MVARS               | MVAR                    |
+---------------------+-------------------------+
