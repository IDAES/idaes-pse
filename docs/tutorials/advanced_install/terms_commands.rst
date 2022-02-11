Terminology and Commands
========================

This section gives a high-level introduction to Git and GitHub terminology and commands.

More details resources include `Atlassian Github tutorials <https://www.atlassian.com/git/tutorials>`_ , 
`GitHub help <https://help.github.com/>`_, and `Git documentation <https://git-scm.com/doc>`_.

.. contents:: :local:

Terminology
-----------

Summary
^^^^^^^

branch
    A name for a series of commits. See :ref:`devterm-branch`.

fork
    Copy of a repository in GitHub. See :ref:`devterm-fork`.

pull request (PR)
    A request to compare and merge code in a GitHub repository. See :ref:`devterm-pr`.

.. _devterm-branch:

Branches
^^^^^^^^

A *branch* is a series of commits that allows you to separate the code development 
from the main code. There is a good description of what Git branches are 
and `how they work here <https://git-scm.com/book/en/v1/Git-Branching-What-a-Branch-Is>`__.
Understanding this takes a little study, but this pays off by making
Git’s behavior much less mysterious. The short, practical version is
that a branch is a name for a series of commits that you want to group
together, and keep separable from other series of commits. From Git's perspective,
the branch is just a name for the first commit in that series.

It is recommended that you create new branches on which to develop your work,
and reserve the "main" branch for merging in work that has been completed
and approved on GitHub. One way to do this is to create branches that correspond
directly to issues on GitHub, and include the issue number in the branch name.

.. _devterm-fork:

Forks
^^^^^

A *fork* is a copy of a repository, in the GitHub shared space (a copy of
a repository from GitHub down to your local disk is called a “clone”).
In this context, that means a copy of the “idaes-dev” repository from
the IDAES organization (https://github.com/IDAES/idaes-dev) to your
own user space, e.g., https://github.com/myname/idaes-dev). The
mechanics of creating and using forks on GitHub are given 
`here <https://help.github.com/articles/fork-a-repo/>`__.

.. _devterm-pr:

Pull Requests
^^^^^^^^^^^^^

A fundamental procedure in the development lifecycle is what is called a
“pull request”. Understanding what these are, and do, is important for
participating fully in the software development process. First,
understand that pull requests are for collaborative development (GitHub)
and not part of the core revision control functionality that is offered
by Git. The official GitHub description of pull requests is
`here <https://help.github.com/articles/about-pull-requests>`_. However,
it gets technical rather quickly, so a higher-level explanation may be
helpful:

Pull requests are a mechanism that GitHub provides to look at what the
code on some branch from your fork of the repository would be like if it
were merged with the main branch in the main (e.g., idaes-pse/idaes-dev)
repository. You can think of it as a staging area where the code is merged
and all the tests are run, without changing the target repository.
Everyone on the team can see a pull request, comment on it, and review
it.

Git Commands
------------

The Git tool has many different commands, but there are several really
important ones that tend to get used as verbs in software development
conversations, and therefore are good to know:

add
    Put a file onto the list of “things I want to commit” (see "commit"),
    called “staging” the file.

commit
    Save the changes in “staged” files into Git (since the last time you did
    this), along with a user-provided description of what the changes mean
    (called the “commit message”).

push
    Move local committed changes to the GitHub-hosted “remote”
    repository by “pushing” them across the network.

pull
    Update your local files with changes from the GitHub-hosted
    “remote” repository by “pulling” them across the network.

Note that the `push` and `pull` commands require GitHub (or some other service
that can host a remote copy of the repository).

