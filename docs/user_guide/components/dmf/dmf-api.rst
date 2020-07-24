.. _dmf-api:

DMF Application Programming Interface (API)
===========================================
This page describes how to use the DMF when you create and save your models.
For information on performing some DMF functions from the command-line,
see :ref:`dmf-cli`. All the modules referenced here are in the :mod:`idaes.dmf` subpackage.

.. contents::
    :local:
    :depth: 1

.. currentmodule:: idaes.dmf

Initialization
--------------
You can create a new :class:`dmfbase.DMF` instance quite simply:

.. code-block:: python

    from idaes.dmf import DMF
    dmf = DMF()  # new DMF instance

When initialized this way, the DMF will use the :ref:`configuration <dmf-config>` it finds
in a file called `.dmf` in the user's home directory. You can specify another
configuration to use. The configuration of the DMF specifies where the workspace is located, which can be
retrieved through the attribute `workspace_path`.

Adding data
-----------
The data in the DMF is broken down into "resources". When adding data to the DMF with
the Python API, you create resources and add them to the DMF instance. A resource
describes one dataset, and contains:

  - *metadata* about the creator, creation and modification time, version, names, and description
  - *provenance* about the sources of the information
  - *data*, as a references to files, embedded structured (JSON) data, or both
  - *codes*, as references to code file locations or module paths, and optionally specific sections of that file or module
  - *relations*, i.e. labeled connections to and from other resources. The following

To add a dataset, you first create a "resource",
which is an instance of :class:`resource.Resource` (in module `resource`). This class
provides some convenience methods for manipulating the underlying structure of the resource,
which is contained in a Python dictionary (in an attribute called `v`, for "values")
and described, using JSON Schema syntax. The schema is contained in the module variable
:attr:`resource.RESOURCE_SCHEMA`. An example of creating a new Resource object:

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

You can also create a resource that describes a file by using the
:meth:`resource.Resource.from_file` method. This will fill in the `datafiles` section of the
resource object.

Once you have the resource object populated, you can add it
to the DMF instance (and, thus, its workspace) with the :meth:`dmfbase.DMF.add` method:

.. code-block:: python

    from idaes.dmf import DMF
    from idaes.dmf.resource import Resource

    r = Resource()
    # ... create resource ...
    dmf = DMF()
    dmf.add(r)

You can create a resource and add it to the DMF in a single step with the :meth:`dmfbase.DMF.new`
method:

.. code-block:: python

    from idaes.dmf import DMF

    dmf = DMF()
    r = dmf.new(file="/path/to/breaking_news.doc",
                author={"name": "Clark Kent", "email": "ckent@dailyplanet.com"})

Once a resource is added to a DMF instance, you can still modify its content,
but you need to call :meth:`dmfbase.DMF.update` to synchronize those changes with the stored
values. This is necessary for adding *relations* between two resources,
which you simply cannot do until both of them are created. But it can also be used
to do things like add a description:

.. code-block:: python

    from idaes.dmf import DMF

    dmf = DMF()
    # create and add resource
    r = dmf.new(file="/path/to/breaking_news.doc")
    # add a description to the resource
    r.v["description"] = get_description()
    # sync the description to the stored value
    dmf.update()


Adding relations
----------------
One of the main functions of the DMF is to track the relationships, or relations, between
its resources. In the lingo of graphs of objects, and in particular the Resource Description Framework (RDF)
that is used as the foundation for many provenance systems, these relations are directed edges between
objects, labeled by "predicates". In this terminology, the resource from which the directed edge starts
is called the "subject" of the relation, and the resource from which the directed edge ends is the "object".
The DMF defines the following predicates (associated module string constants are shown in parentheses):

* (:attr:`resource.PR_DERIVED`) derived: object is derived from subject
* (:attr:`resource.PR_CONTAINS`) contains: object contains the subject
* (:attr:`resource.PR_USES`) uses: object uses the subject
* (:attr:`resource.PR_VERSION`) version: object is a (new) version of the subject

Adding a relation between two resources is pretty straightforward. You create both resources and add them to the DMF, then
create a "triple" to describe the connection between them (with the "predicate" that labels that connection),
with the :func:`resource.create_relation` function.
Then you call the :meth:`dmfbase.DMF.update` function on the DMF instance to save the relation:

.. code-block:: python

    from idaes.dmf import DMF
    from idaes.dmf.resource import Triple, PR_DERIVED
    from idaes.dmf.resource import create_relation_args

    dmf = DMF()
    # create and add resources
    r1 = dmf.new(file="/path/to/breaking_news.doc")
    r2 = dmf.new(file="/path/to/interview_notes.txt")
    # create relation (news --was derived from--> notes)
    create_relation(r1, PR_DERIVED, r2)
    # sync the relation to the DMF
    dmf.update()

