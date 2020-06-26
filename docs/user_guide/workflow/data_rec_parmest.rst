Data Reconciliation and Parameter Estimation
============================================

This workflow generally describes features of the IDAES framework that are useful
for data reconciliation and parameter estimation.  Many of these features can
generally be used for any task where plant data is to be used in conjunction with
a process model.

This workflow provides general information about IDAES functionality laid out in
terms of typical use cases.  To see a complete specific examples of complete
data reconciliation and parameter estimation examples, see
:ref:`﻿Tutorials and Examples<tutorials_examples:Tutorials and Examples>`.
Relevant tutorials can be found in
``tutorials/advanced/data_recon_and_parameter_estimation``.

Data Management
---------------

.. currentmodule:: idaes.dmf.model_data

IDAES has functions to read and mange process data process data. Data management
functions are contained in the ``idaes.dmf.model_data`` module.


Reading Data
~~~~~~~~~~~~
A set of process data can be stored in two files, a data file and a metadata file.
The data file is a CSV file where the first column first column contains a data
point index, usually a timestamp.  The first row contains a header with the
measurement tag names. The rest of the files contains measurement data.  The
data file format is shown in the table below.

+---------+-----------+------+
|         | tag 1     | ...  |
+=========+===========+======+
| index 1 | data(1,1) | ...  |
+---------+-----------+------+
| index 2 | data(2,1) | ...  |
+---------+-----------+------+
| ...     | ...       | ...  |
+---------+-----------+------+

Meta data is provided for each tag in a separate CSV file.  The meta data file
has no header row, and aside from the tag name any column my be blank.  In the meta
data csv file, the first column is the measurement tag name, the second column is
a string that maps the tag to a specific model quantity, the third column is a
description of the tag, and the fourth column is a for units of measure. Any columns
past the fourth column are ignore and can be used to store any additional information.

+--------+--------------------------+---------------+--------------------+
| tag 1  | model reference string 1 | description 1 | units of measure 1 |
+--------+--------------------------+---------------+--------------------+
| tag 2  | model reference string 2 | description 2 | units of measure 1 |
+--------+--------------------------+---------------+--------------------+
| ...    | ...                      | ...           | ...                |
+--------+--------------------------+---------------+--------------------+

The unit strings should be interpretable by `pint <https://pint.readthedocs.io/en/stable/>`_,
with some additional unit strings defined in ``idaes.dmf.model_data._unit_strings``.
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

Process data can be divvied into bins based on some criteria.  This allows for
estimating measurement uncertainty when no better information is available, and
provides a way to look at how different process measurements vary for operating
conditions that should be similar in some way. As an example consider a power
plant data could be binned by power output, and assuming that operating
procedures are standard, it could be assumed that measurements in each bin should
be about the same.  The variance of a measurement in a bin could be used a an
approximation of uncertainty.  Binning the data by load and time could show how
process measurements change over time and be useful for things like fault detection
and equipment degradation.

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
``idaes.dmf.data_rec_plot_book()``

.. autofunction:: data_rec_plot_book
  :noindex:


Tagging the Model
-----------------

Mapping process data to a model is typically done by creating model tag
dictionaries.  Where the dictionary key is a measurement tag and the value is
a reference to a quantity in a model. The tags may be process measurement tags,
or any other convenient string. IDAES has some utilities to help facilitate
tagging of models.

If model reference strings where provided in the tag metadata file, the tag metadata
will from the ``idaes.dmf.model_data.read_data()`` will contain model references.
These references can be accessed in the metadata dictionary as ``metadata[tag]["reference"]``.

The reference strings are optional in the tag metadata file are optional and can be
added after reading the data. to add a reference string to the metadata you can
update the metadata dictionary like ``metadata[tag]["reference_string"] = reference_string``.
If you update the reference string, the idaes.dmf.model_data.upadate_metadata_model_references()``
function can be used to update the references in the tag metadata.

.. autofunction:: upadate_metadata_model_references
  :noindex:

An easy way to create a new tag dictionary from tag metadata is to use dictionary
comprehension like so ``data_tags = {k:v["reference"][0] for k, v in metadata.items() if v["reference"] is not None}``

Often data reconciliation is performed using process data as the first step of
parameter estimation or model validation. Data reconciliation can often fill in
information for many unmeasured quantities.  Most of this data can be associated
with process streams.  To make managing proliferation of data tags easier to manage
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
