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
