Developer introductory material
===============================
This section gives terms and information targeted at people who are new to collaborative
software development. It serves as background for understanding the
:ref:`collab_dev` procedures.

Terminology
-----------
`Git <https://git-scm.com/>`__
    A “version control system”, for keeping track of changes in a set of files

branch
    A name for a series of commits that you want to group together.

fork
    Copy of a repository, in Github, typically from its original author or
    organization into your own namespace in Github. For example, you can
    fork the `idaes-dev` repository from its home in the `idaes` organization,
    which is written `idaes/idaes-dev`, to your own namespace
    `myname/idaes-dev`.

`Github <https://github.com>`__
    A hosting service for Git
    repositories that adds many other features that are useful for
    collaborative software development. Two key concepts that are part of
    Github are:

        * Fork: A copy of a Github repository. See “Forks” below.

        * Pull Request: A request to compare and merge code in a Github
                        repository. See “Pull Requests” below.

Git commands
^^^^^^^^^^^^
The Git tool has many different commands, but there are several really
important ones that tend to get used as verbs in software development
conversations, and therefore are good to know:

add
    Put a file onto the list of “things I want to commit” (seecommit),
    called “staging” the file.

commit
    Save the changes in “staged” files (since the last time you did
    this), along with a user-provided description of what the changes mean
    (called the “commit message”).

push
    Move local committed changes to the Github-hosted “remote”
    repository by “pushing” them across the network.

pull
    Update your local files with changes from the Github-hosted
    “remote” repository by “pulling” them across the network.

.. note:: The `push` and `pull` commands require Github (or some other service
          that can host a remote copy of the repository).

Branches
^^^^^^^^
There is a good description of what git branches are and how they work
`here <https://git-scm.com/book/en/v1/Git-Branching-What-a-Branch-Is>`_.
Understanding this takes a little study, but this pays off by making
git’s behavior much less mysterious. The short, practical version is
that a branch is a name for a series of commits that you want to group
together,and keep separable from other series of commits. From git's perspective,
the branch is just a name for the first commit in that series.

It is recommended that you create new branches on which to develop your work,
and reserve the "master" branch for merging in work that has been completed
and approved on Github. One way to do this is to create branches that correspond
directly to issues on Github, and include the issue number in the branch name.

.. _devterm_fork:

Forks
^^^^^
A `fork` is a copy of a repository, in the Github shared space (a copy of
a repository from Github down to your local disk is called a “clone”).
In this context, that means a copy of the “idaes-dev” repository from
the IDAES organization (\ https://github.com/IDAES/idaes-dev\ ) toyour
own user space (e.g.\ https://github.com/dangunter/idaes-dev\ ). The
mechanics of creating and using forks onGithub are given
here:\ https://help.github.com/articles/fork-a-repo/\ . In the Software
Development section, details are given on how to create and use this
fork for collaborative software development.

.. _devterm_pr:

Pull Requests
^^^^^^^^^^^^^
A fundamental procedure in the development lifecycle is what is called a
“pull request”. Understanding what these are, and do, is important for
participating fully in the software development process. First,
understand that pull requests are for collaborative development (Github)
and not part of the core revision control functionality that is offered
by Git. The official Github description of pull requests is
`here <https://help.github.com/articles/about-pull-requests>`_ . However,
it gets technical rather quickly, so a higher-level explanation may be
helpful:

Pull requests are a mechanism that Github provides to look at what the
code on some branch from your fork of the repository would be like if it
were merged with the master branch in the main (e.g., idaes/idaes-dev)
repository. You can think of it as a staging area where the code is merged
and all the tests are run, without changing the target repository.
Everyone on the team can see a pull request, comment on it, and review
it.

Details on how to use pull requests (or “PRs”) are given in the
:ref:`collab_dev` section.
For now, suffice it to say that PRs are central to
collaborative development. You _must_ use them to add code to the project.
