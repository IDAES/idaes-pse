Github repository overview
==========================
This section describes the layout of the
`Github repositories <https://help.github.com/articles/about-repositories/>`_.
Later sections will give guidelines for contributing code to these
repositories.

Repositories
------------
+-----------------+----------+-----------------------+
| Repository name | Public?  | Description           |
+=================+==========+=======================+
| idaes-pse       | Yes      | Main public           |
|                 |          | repository, including |
|                 |          | core framework and    |
|                 |          | integrated tools      |
+-----------------+----------+-----------------------+
| idaes-dev       | No       | Main private          |
|                 |          | repository, where     |
|                 |          | code is contributed   |
|                 |          | before being          |
|                 |          | “mirrored” to the     |
|                 |          | public `ideas-pse`    |
|                 |          | repository            |
+-----------------+----------+-----------------------+
| workspace       | No       | Repository for code   |
|                 |          | that does not belong  |
|                 |          | to any particular     |
|                 |          | CRADA or NDA, but     |
|                 |          | also is never         |
|                 |          | intended to be        |
|                 |          | released open-source  |
+-----------------+----------+-----------------------+

The URL for an IDAES repository, e.g. “some-repo”, will be
``https://github.com/IDAES/some-repo``.

Public vs. Private
------------------
All these repositories except for “idaes-pse” will only be visible on
Github, on the web, for people who have been added to the IDAES
developer team in the IDAES “organization” (See `About Github
organizations <https://help.github.com/articles/about-organizations/>`_).
The `idaes-pse` repository will be visible to anyone, even
people without a Github account.
