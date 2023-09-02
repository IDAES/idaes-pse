How-to use the Data Management Framework
========================================
The Data Management Framework (DMF) is used to organize data within IDAES. It is often not visible, but does have
an API and command-line tool ('dmf') to view and edit its information.

.. contents::

Configure
---------

Create initial configuration (bootstrap)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The DMF uses a file in the user's home directory called, by default, ".dmf" to give the location of the currently
active workpace (and, potentially, other things in the future). To create this file and initialize it with the
default workspace for the current IDAES installation, you can use either a command-line or Python call. You only
need to do this once per home-directory (which for most people means once, period).

Command-line::

    dmf setup

API call::

    from idaes.core.dmf import create_configuration
    result = create_configuration()
    # result will be the location of the configuration file, as a
    # Python pathlib.Path object.


Tabular Data
-------------

The DMF has built-in support for a set of data tables (in CSV or Excel) associated with a reference
such as a journal publication. This page documents how to :ref:`use existing tables <dmftbl_find>` and how to
:ref:`add your own tables <dmftbl_create>` for others to use.

.. _dmftbl_find:

Find and use existing data
^^^^^^^^^^^^^^^^^^^^^^^^^^

Use the `idaes.core.datasets` module.
You can list all available datasets with `available()`.

.. code-block:: python

    from idaes.core import datasets
    datasets.available()
    # returns a dictionary keyed by name with namedtuple value, like:
    # {'Pitzer': AvailableResult(Class=<class 'idaes.core.datasets.Pitzer'>,
    #                           description='Pitzer(1984) publication and related tables.')}

You can either use the class directly by its name, or pull it out of the `available()` result.
In either case, instantiate the class to pull its data from the DMF:

.. code-block:: python

    # the next 2 lines are equivalent
    pz = datasets.Pitzer()
    pz = datasets.available()["Pitzer"].Class()

Once you have the dataset, you can see which tables are present with the `list_tables()` method.

.. code-block:: python

    pz.list_tables()
    # sample output:
    # ['Standard G', 'Standard S', ..., 'Specific Enthalpy']

Then to retrieve and use a given table, use `get_table()` with the name of the desired table:

.. code-block:: python

    t = pz.get_table('Standard G')
    # t.data is a Pandas dataframe
    # t.description has a (text) description of the table contents
    # t.units and t.units_dict give the units by column index and name

.. _dmftbl_create:

Create and add new data
^^^^^^^^^^^^^^^^^^^^^^^

.. py:currentmodule:: idaes.core.dmf.datasets

To create your own tables, you need to create a new directory and put two things in it:
a configuration file called `dataset.json`, and data files in Excel or CSV format.

An excerpt of a `dataset.json` is shown below:

.. code-block:: json

    {
        "name": "ThermoStuff:1999",
        "text": {
            "file": "Thermodynamic Properties.pdf",
            "title": "Thermodynamic Properties of Some Stuff",
            "date": "1999",
            "authors": "Joe Bazooka, Carl Froffenheffer, Andrew Lee",
            "venue": "Journal of Interesting Observations",
            "doi": "https://doi.org/10.1063/1.5551212"
        },
        "tables": [
          {
            "name": "Standard G",
            "description": "Standard Gibbs energies, according to some guy I know",
            "datafile": "gibbs_std.csv"
          },
        ]
    }

.. note:: Make sure you pick a unique string for the ``name`` field, since this will be the key by which this
    publication and associated data are found. Using the same name for different publications will result in one
    overwriting the other and other bad behavior.

The `dataset.json` above referred to one file that had the text of the publication
and one file with comma-separated values of the data. Copy or move all these files into the same
directory, let's call it `DataDir`, whose contents will now be:

    * dataset.json
    * Thermodynamic Properties.pdf
    * gibbs_std.csv

Then you can load this directory of data into the DMF with the following command-line:

.. code-block:: shell

    dmf load --global DataDir

The ``--global`` option means "use the default global DMF workspace instead of any current workspace". If you
choose to use your own workspace instead, you'll have to pass it in explicitly later, e.g., to the subclass of Publication that you
create below.

The data file format is a header row plus data, with the only "special" aspect being that
if the name of a column in the header ends with some text in square brackets, that text is assumed
to be the units for the values. For example:

::

    Temperature [K], Pressure [Pa], Value
    100, 90, 12.34

In this table, the units "K" and "Pa" would be parsed out of the first two columns, and the units of the
third column would be empty.

Finally, you can make the dataset accessible as a class by subclassing from :class:`Publication`
and invoking the superclass with the name of the dataset.
The key part to get right here is that the name used in the class constructor must match the ``name``
field from the `dataset.json` configuration file. For example, with the configuration given above:

.. code-block:: python

    from idaes.core.dmf.datasets import Publication

    class ThermoStuff(Publication):
        def __init__(self, **kwargs):
            super().__init__("ThermoStuff:1999", **kwargs)

General API usage
-----------------

Create a DMF resource
^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    from idaes.core.dmf.resource import Resource

    r = Resource()
    r.v["version_info"]["version"] = test_version
    r.v["collaborators"] = [
        {"name": "Clark Kent", "email": "ckent@dailyplanet.com"},
        {"name": "Superman", "email": "sman@fortress.solitude.org"},
    ]
    r.v["sources"].append(
        {
            "isbn": "978-0201616224",
            "source": 'Hunt, A. and Thomas, D., "The Pragmatic Programmer", '
            "Addison-Wesley, 1999",
            "date": "1999-01-01",
        }
    )
    r.v["codes"].append(
        {
            "type": "function",
            "name": "test_resource_full",
            "desc": "The test function",
            "location": "test_files.py",
            "version": test_version,
        }
    )
    r.v["datafiles"].append({"path": "/etc/passwd"})
    r.v["aliases"] = ["test_resource_full"]
    r.v["tags"] = ["test", "resource"]
    r.data = {"arbitrary": {"values": [1, 2, 3]}}
    return r

