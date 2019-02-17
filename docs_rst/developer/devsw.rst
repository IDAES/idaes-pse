.. _collab_dev:

Collaborative software development
==================================

.. note:: This section of the developer documentation is a work in progress.

This section gives guidance for all developers on the project.

Although the main focus of this project is developing open source software (OSS),
it is also true that some of the software may be developed internally or in
coordination with industry under a :term:`CRADA` or :term:`NDA`.

*It is the developer's responsibility*, for a given development effort,
to keep in mind what role you must assume and thus which set of procedures
must be followed.

CRADA/NDA
  If you are developing software covered by a CRADA, NDA, or other legal
  agreement that does not explicitly allow the data and/or code to be
  released as open-source under the IDAES license, then you must follow
  procedures under :ref:`collab_dev_nda`.

Internal
  If you are developing non-CRADA/NDA software, which is not intended to be
  part of the core framework or (ever) released as open-source then follow procedures
  under :ref:`collab_dev_int`.

Core/open-source
  If you are developing software with no proprietary data or code, which
  is intended to be released as open-source with the core framework, then follow
  procedures under :ref:`collab_dev_oss`.

.. _collab_dev_nda:

Developing Software with Proprietary Content
--------------------------------------------
TBD

.. _collab_dev_int:

Developing Software for Internal Use
------------------------------------
TBD

.. _collab_dev_oss:

Developing software for Open-source Release
-------------------------------------------
In this section we lay out all the essential steps for developing
software that will end up in the open-source releases for IDEAS.
We can break this process into four distinct steps:

1. Setup: Prepare your local system for collaborative development.
2. Initiate: Notify collaborators of intent to make some changes
3. Develop: Make local changes.
4. Collaborate: Push the changes to Github and get feedback from other developers.
5. Merge: Merge the changes into the shared "master" branch.

As illustrated in the diagram below, the first phase only needs to happen once,
whereas the remaining phases
are performed for every new "topic" (e.g. a bugfix or new feature). The develop and
collaborate phases are performed in a loop until changes are approved by the team.

.. figure:: ../_static/sw-dev-workflow.png
    :align: right
    :height: 200px


    Software development workflow

1. Setup
^^^^^^^^
Before you can start developing software collaboratively,
you need to make sure you are set up in Github and set up your local development environment.

Github setup
~~~~~~~~~~~~
To work within the project, you need to create a login on `Github`_. You also
need to make sure that this login has been added to the IDAES organization.

If these steps are successful, you should be able to login to Github, visit the
`IDAES Github organization <https://github.com/IDAES/>`_, and see "Private" repositories
such as `idaes-dev` and `workspace`.

.. _Github: https://github.com/

Fork the repo
~~~~~~~~~~~~~
You use a "fork" of a repository (or "repo" for short) to create a space where you
can save changes without directly affecting the main repository. Then, as we will see,
you _request_ that these changes be incorporated (after review).

This section assumes that the repository in question is ``idaes-dev``,
but the idea is the same for any other repo.

You should first visit the repo on Github
by pointing your browser to https://github.com/IDAES/idaes-dev/. Then you should
fork the repo into a repo of the same name under your name.

.. figure:: ../_static/github-fork-repo.png
    :align: center

    Screenshot showing where to click to fork the Github repo

Clone your fork
~~~~~~~~~~~~~~~
A "clone" is a copy of a Github repository on your local machine. This is what
you need to do in order to actually edit and change the files.
To make a clone of the fork you created in the previous step,
change to a directory where you want to put the source code and run the command::

    git clone git@github.com:MYNAME/idaes-dev.git

Of course, replace MYNAME with your login name. This will download all the files in
the latest version of the repository onto your local disk.

Create the Python environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Once you have the repo cloned, you can change into that directory (by default, it
will be called "idaes-dev" like the repo) and install the Python packages.

But before you do that, you need to get the Python package manager fully up and
running. We use a Python packaging system called Conda_.
Below are instructions for installing a minimal version of Conda, called Miniconda_.
The full version installs a large number of scientific analysis and visualization libraries
that are not required by the IDAES framework.

.. _Conda: https://conda.io/
.. _Miniconda: https://conda.io/en/latest/miniconda.html

.. code-block:: sh

    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh

Create and activate a conda environment (along with its own copy of ``pip``)
for the new IDAES installation **(you will need to** ``conda activate idaes``
**when you open a fresh terminal window and wish to use IDAES)**:

.. code-block:: sh

    conda create -n idaes pip
    conda activate idaes

Now that conda and pip are installed, and you are in the "idaes" conda environment,
you can run the standard steps for installing a Python package in development mode:

.. code-block:: sh

    pip install -r requirements.txt
    python setup.py develop

You can test that everything is installed properly by running the tests with
Pytest_:

.. code-block:: sh

    pytest

.. _Pytest: https://pytest.org/

2. Initiate
^^^^^^^^^^^
We will call a set of changes that belong together, e.g. because they depend on
each other to work, a "topic". This section describes how to start work on a new
topic. The workflow for initiating a topic is shown in the diagram below.

.. figure:: ../_static/sw-init-workflow.png
    :align: right
    :height: 400px

    Initiate topic workflow


Create an issue on Github
~~~~~~~~~~~~~~~~~~~~~~~~~
To create an issue on Github, simply navigate to the repository page and click on
the "Issues" tab. Then click on the "Issues" button and fill in a title and brief
description of the issue. You do not need to list details about sub-steps required
for the issue, as this sort of information is better put in the (related) pull
request that you will create later. Assign the issue to the appropriate people,
which is often yourself.

There is one more important step to take, that will allow the rest of the project
to easily notice your issue: add the issue to the "Priorities" project. The screenshot
below shows where you need to click to do this.

.. image:: ../_static/github-issue-priority.png
    :align: center

Create a branch on your fork
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
It is certainly possible to do your work on your fork in the "master"
branch. The problem that can arise here is if you need to do two unrelated
things at the same time, for example working on a new feature and fixing
a bug in the current code. This can be quite tricky to manage as a single set
of changes, but very easy to handle by putting each new set of changes in
its own branch, which as was mentioned earlier we call a *topic* branch.
When all the changes in the branch are done and merged, you can delete it
both locally and in your fork so you don't end up with a bunch of old branches
cluttering up your git history.

The command for doing this is simple:

.. code-block:: sh

    git co -b <BRANCH-NAME>

The branch name should be one word, with dashes or underscores as needed.
One convention for the name that can be helpful is to include the Issue number
at the end, e.g. ``git co -b mytopic-issue42``. This is especially useful later
when you are cleaning up old branches, and you can quickly see which branches
are related to issues that are completed.

Make local edits and push changes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
A new branch, while it feels like a change, is not really a change in the
eyes of Git or Github, and by itself will not allow you to start a new pull
request (which is the goal of this whole phase). The easiest thing to do do is
a special "empty" commit::

    git commit --allow-empty -m 'Empty commit so I can open a PR'


Since this is your first "push" to this branch, you are going to need to set an upstream
branch on the remote that should receive the changes. If this sounds complicated,
it's OK because git actually gives you cut-and-paste instructions. Just run
the ``git push`` command with no other arguments::

    $ git push
    fatal: The current branch mybranch-issue3000 has no upstream branch.
    To push the current branch and set the remote as upstream, use

        git push --set-upstream origin mybranch-issue3000

Cut and paste the suggested command, and you're ready to go. Subsequent
calls to "push" will not require any additional arguments to work.

Start a new Pull Request on Github
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Finally, you are ready to initiate the pull request. Right after you perform the
``push`` command above, head to the repository
URL in Github (https://github.com/IDAES/idaes-dev) and you should see a highlighted
bar below the tabs, as in the figure below, asking if you want to start a pull-request.

.. image:: ../_static/github-start-pullrequest.png
    :align: center

Click on this and fill in the requested information. Remember to link to the issue
you created earlier.

Depending on the Github plan, there may be a pull-down menu for creating the pull
request that lets you create a "draft" pull request. If that is not present, you
can signal this the old-fashioned way by adding "[WIP]" (for Work-in-Progress) at
the beginning of the pull request title.

Either way, create the pull request. Do *not* assign reviewers until you are done
making your changes (which is probably not now). This way the assigning of reviewers
becomes an unambiguous signal that the PR is actually ready for review.

3. Develop
^^^^^^^^^^
The development process is a loop of adding code, testing and
debugging, and committing and pushing to Github. You may go through many (many!)
iterations of this loop before the code is ready for review.

.. note:: Avoid having pull requests that take months to complete. It is
          better to divide up the work, even artificially, into a piece that
          can be reviewed and merged into the main repository within a week or two.

Run tests
~~~~~~~~~
After significant batches of changes, you should make sure you have tests
for the new/changed functionality. This involves writing :ref:`unit-tests` as
well as running the test suite and examining the results of the :ref:`coverage-tests`.
The automated testing that will occur later will fail if the tests fail, or the
changes result in less overall testing (aka "lower test coverage"). See the
linked sections for detailed instructions.

.. _git-commit:

Commit changes
~~~~~~~~~~~~~~
TBD

.. _git-push:

Push changes to Github
~~~~~~~~~~~~~~~~~~~~~~
TBD

4. Collaborate
^^^^^^^^^^^^^^
TBD

Request review
~~~~~~~~~~~~~~
TBD

Keep your branch up to date
~~~~~~~~~~~~~~~~~~~~~~~~~~~
TBD

5. Merge
^^^^^^^^
TBD

Code Review Procedures
^^^^^^^^^^^^^^^^^^^^^^
.. note:: “It’s a simple 3-step process. Step one: Fix! Step two: It! Step three:
Fix it!” -- Oscar Rogers (Kenan Thompson), Saturday Night Live, 2/2009

Code review is the last line of defense between a mistake that the IDAES
team will see and a mistake the whole world will see. In the case of
that mistake being a leak of proprietary information, the entire project
is jeopardized, so we need to take this process seriously.

Automated Checks
~~~~~~~~~~~~~~~~
The first level of code review is a set of automated checks that *must* pass
before the code is ready for people to review it. These checks will run
on the initiation of a :ref:`pull request <devterm_pr>` and on every new commit to that pull
request that is pushed to Github (thus the name “continuous
integration”).

The “continuous integration” of the code is hosted by an online service
– we use `CircleCI <https://circleci.com>`_ -- that can automatically
rerun the tests after every change (in this case, every new Pull Request
or update to the code in an existing Pull Request) and report the
results back to Github for display in the web pages. This status
information can then be used as an automatic gatekeeper on whether the
code can be merged into the master branch – if tests fail, then no merge
is allowed. Following this procedure, it is not possible for the master
branch to ever be failing its own tests.

Testing
^^^^^^^

.. note:: “Discovering the unexpected is more important than confirming the known.“
          -- George E. P. Box

Testing is essential to the process of creating software. The longer one programs,
the more the idea that "If it isn't tested, it probably doesn't work" makes sense.
This section outlines the different kinds of testing and how to write/use them in
the project.

.. _unit-tests:

Unit tests
~~~~~~~~~~
Testing individual pieces of functionality, including the
ability to report the correct kind of errors from bad inputs. Unit tests
must always run quickly. If it takes more than 5 seconds, it is not a unit
test, and it is expected that most unit tests take well under 1 second.
The reason for this is that the entire unit test suite is run on every
change in a Pull Request, and should also be run relatively frequently
on local developer machines. If this suite of hundreds of tests takes
more than a couple of minutes to run, it will introduce a significant
bottleneck in the development workflow.

.. _coverage-tests:

Code coverage tests
~~~~~~~~~~~~~~~~~~~
The “coverage” of the code refers to what percentage of
the code (“lines covered” divided by total lines) is executed by the
automated tests. This is important because passing automated tests is
only meaningful if the automated tests cover the majority of the code’s
behavior. This is not a perfect measure, of course, since simply
executing a line of code under one condition does not mean it would
execute correctly under all conditions. The code coverage is evaluated
locally and then integrated with Github through a tool called `Coveralls
<https://coveralls.io>`_.

Code Review
^^^^^^^^^^^

Summary
~~~~~~~
Every piece of code must be reviewed by at least two people.

In every case, one of those people will be a designated “gatekeeper” and
the one or more others will be “technical reviewers”.

The technical reviewers are expected to consider various aspects of the
proposed changes (details below), and engage the author in a discussion
on any aspects that are deemed lacking or missing.

The gatekeeper is expected to make sure all criteria have been met, and
actually merge the PR.

Assigning Roles

The gatekeeper is a designated person, who will always be added to
review a Pull Request (PR)

Gatekeeper is a role that will be one (?) person for some period like a
week or two weeks

The role should rotate around the team, it’s expected to be a fair
amount of work and should be aligned with availability and paper
deadlines, etc.

The originator of the PR will add as reviewers the gatekeeper and 1+
technical reviewers.

Originator responsibilities

The originator of the PR should include in the PR itself information
about where to find:

Changes to code/data

Tests of the changes

Documentation of the changes

The originator should be responsive to the reviewers

Technical reviewer responsibilities

The technical reviewer(s) should look at the proposed changes for

Technical correctness (runs properly, good style, internal code
documentation, etc.)

Tests

Documentation

No proprietary / sensitive information

Until they approve, the conversation in the PR is between the technical
reviewers and the originator (the gatekeeper is not required to
participate, assuming they have many PRs to worry about)

Gatekeeper responsibilities

The gatekeeper does not need to engage until there is at least one
approving technical review.

Once there is, they should verify that:

Changes do not contain proprietary data

Tests are adequate and do not fail

Documentation is adequate

Once everything is verified, the gatekeeper merges the PR

 
