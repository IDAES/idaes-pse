Defining Pure Component Properties
==================================

Most methods for calculating the thermophysical properties of materials start from estimating the properties of each component in its pure form, before applying mixing rules todetermine the properties of the mixture. Pure component properties generally take the form of empirical correlations as a fucntion of material state (generally temperature) derived from experimental data. Data and correlations for many component components are readily avaialble in literature. However due to the empirical nature of these correlations and the wide range of data available, different sources use differnt forms for their correlations.

Within the IDAES Generic Property Package Framework, pure component property correlations are provided in the form of Python methods whcih return a Pyomo expression relating the pure component property to the material state (using the :ref:`standard naming conventions<standards:Standard Variable Names>`. IDAES provides a number of libraries of containing common forms for these correaltions, and a list of the libraries currently supported by IDAES is given below.

A list of all the pure component properties currently supported by the IDAES Generic proeprty Package Framework can be found in the :ref:`developers section<property_models/general/developers:Pure Component Properties>` of this documentation.

Pure Component Libraries
------------------------

.. toctree::
    :maxdepth: 1

    pure/NIST
    pure/Perrys
    pure/RPP
