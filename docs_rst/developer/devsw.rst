.. _collab_dev:

Collaborative software development
==================================

.. note:: This section of the developer documentation is a work in progress.

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
    :align: center

    Software development workflow

1. Setup
^^^^^^^^
TBD

Fork the repo
~~~~~~~~~~~~~
TBD

Clone your fork
~~~~~~~~~~~~~~~
TBD

2. Initiate
^^^^^^^^^^^
We will call a set of changes that belong together, e.g. because they depend on
each other to work, a "topic". This section describes how to start work on a new
topic.

Create an issue on Github
~~~~~~~~~~~~~~~~~~~~~~~~~

Create a branch on your fork
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
You can have multiple branches at the same time. This is why we are
recommending you always create a new branch for a new issue instead of
just using the “master” branch, which git would definitely let you do.
This type of branch is called a topic branch. The advantage is, if you
have to delay your work and switch to fixing some other issue, you can
keep the work separated by simply switching back and forth between
branches. This is especially useful for coordination inGithub, since
each pull request (PR) is associated with a specific branch. You can
have any number of active PRs in parallel if you use separate branches
for each one.

Tip: Assuming that you have followed the procedure of creating a Github
issue first, you can include a descriptive word and the issue number in
the branch name, e.g. “gibbs-issue163”. This will make it easier to keep
track of branches, and easier to feel confident about deleting them once
changes are merged.

Start a new Pull Request on Github
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TBD

3. Develop
^^^^^^^^^^
TBD

Run tests
~~~~~~~~~
TBD

Commit changes
~~~~~~~~~~~~~~
TBD

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
The first level of code review is a set of automated checks that _must_ pass
before the code is ready for people to review it. These checks will run
on the initiation of a :ref:`pull request <devterm_pr>` and on every new commit to that pull
request that is pushed to Github (thus the name “continuous
integration”).

The “continuous integration” of the code is hosted by an online service
– we use CircleCI (\ https://circleci.com\ )-- that can automatically
rerun the tests after every change (in this case, every new Pull Request
or update to the code in an existing Pull Request) and report the
results back to Github for display in the web pages. This status
information can then be used as an automatic gatekeeper on whether the
code can be merged into the master branch – if tests fail, then no merge
is allowed. Following this procedure, it is not possible for the master
branch to ever be failing its own tests.

Types of tests
~~~~~~~~~~~~~~
Unit tests: Testing individual pieces of functionality, including the
ability to report the correct kind of errors from bad inputs. Unit tests
mustalwaysrun quickly. If it takes more than 5 seconds, it is not a unit
test, and it is expected that most unit tests take well under 1 second.
The reason for this is that the entire unit test suite is run on every
change in a Pull Request, and should also be run relatively frequently
on local developer machines. If this suite of hundreds of tests takes
more than a couple of minutes to run, it will introduce a significant
bottleneck in the development workflow.

Code coverage: The “coverage” of the code refers to what percentage of
the code (“lines covered” divided by total lines) is executed by the
automated tests. This is important because passing automated tests is
only meaningful if the automated tests cover the majority of the code’s
behavior. This is not a perfect measure, of course, since simply
executing a line of code under one condition does not mean it would
execute correctly under all conditions. The code coverage is evaluated
locally and then integrated with Github through a tool called Coveralls
(\ https://coveralls.io/\ ).

Code Review

Summary

Every piece of code must be reviewed by at least two people.

In every case, one of those people will be a designated “gatekeeper” and
the one or more others will be “technical reviewers”.

Thetechnical reviewersare expected to consider various aspects of the
proposed changes (details below), and engage the author in a discussion
on any aspects that are deemed lacking or missing.

Thegatekeeperis expected to make sure all criteria have been met, and
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
approving technical review

Once there is, they should verify that

Changes do not contain proprietary data

Tests are adequate and do not fail

Documentation is adequate

Once everything is verified, the gatekeeper merges the PR

 
