Data Management Framework Resource Examples
============================================

The ``Resource`` class, in the ``idaes.dmf.resource`` package, is the base unit of information
that is stored in the Data Management Framework. Below are examples of creating and using
``Resource`` objects.

.. contents::
    :local:
    :depth: 1

Create a resource
-----------------
.. code-block:: python

    from idaes.dmf.resource import Resource

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


Create a resource representing a CSV table
------------------------------------------
.. code-block:: python

    from idaes.dmf.resource import Resource
    from idaes.dmf.tables import Table

    rsrc = Resource(type_=ResourceTypes.tabular)
    rsrc.set_field("name", "table_test_resource")
    rsrc.set_field("desc", "Example table")

    t = Table()
    # Assume existing comma-separated values file "data.csv", with a header row, e.g.
    #     Index,Pressure,Temperature
    #     1,256.12,38.7
    #     2,125.61,73.8
    t.read_csv("data.csv")
    t.add_to_resource(rsrc)

Create a resource representing a publication and multiple CSV tables
---------------------------------------------------------------------
In this case, the publication is the source for 10 related CSV tables.
Each table has a header, with optional units in square brackets, and then value rows:

=====  ======= =============  ==============  ===================  ===========  ===========  ===================
T [C]  P [bar] G0/RT H2O [-]  G0/RT NaCl [-]  A phi [(kg/mol^0.5]  B0 [kg/mol]  B1 [kg/mol]  10^3xC [(kg/mol)^2]
=====  ======= =============  ==============  ===================  ===========  ===========  ===================
0      1       -23.4638       -13.836         0.3767               0.0493       0.2462       2.58
10     1       -22.91         -13.871         0.3821               0.0619       0.2612       1.68
...    ...     ...            ...             ...                  ...          ...          ...
=====  ======= =============  ==============  ===================  ===========  ===========  ===================

The source publication and the 10 datasets (CSV files) are each separate `resources` in the DMF, with a
relationship of "derived [from]" linking them. This can be pictured like this::

                  ┌───────────────┐
                  │ pitzer-source │
                  └──────▲────────┘
                         │    derived from
               ┌─────────┴────────┬─────────────────────┐
               │                  │                     │
        ┌──────┴───────┐  ┌───────┴─────┐      ┌────────┴─────┐
        │ dataset1     │  │ dataset2    │      │ dataset10    │
        │ pitzer_1.csv │  │ pitzer_2.csv│ ...  │ pitzer_10.csv│
        └──────────────┘  └─────────────┘      └──────────────┘

And finally, the code:

.. code-block:: python

    # imports
    from idaes.dmf import DMF, resource, tables
    from pathlib import Path

    workspace_name = "pitzer_data"

    # Create DMF workspace to store data
    dmf = DMF(workspace_name, create=True)

    # Create/add resource representing source publication
    src = resource.Resource(name="pitzer-source", type_=resource.ResourceTypes.publication)
    src.v["sources"].append({"date": "1984",
                       "doi": "10.1063/1.555709",
                       "source": 'Pitzer, Kenneth S., J. Christopher Peiper, and R. H. Busey. '
                       '"Thermodynamic properties of aqueous sodium chloride solutions." '
                       'Journal of Physical and Chemical Reference Data 13.1 (1984): 1-102.'})
    dmf.add(src)


    # Find tables to read, i.e. look for CSV files
    data_files = Path(".").glob("*.csv")

    # Create/add resources for each table & link them to source publication
    for df in data_files:
        name = df.stem
        tbl = tables.Table()
        tbl.read_csv(df)
        r = resource.Resource(name=name, type_=resource.ResourceTypes.tabular)
        dmf.add(r)
        # Create the link: <publication source> -- derived --> <this table>
        resource.create_relation(src, resource.Predicates.derived, r)
        print(f"created resource {name}")
    # register all the relations we just created
    dmf.update()

    # Find all the resources derived from the Pitzer publication.
    pub = dmf.find_one(name="pitzer-source")
    for r in dmf.find_related_resources(pub, outgoing=True):
        print(f"Resource: {r.name}")
