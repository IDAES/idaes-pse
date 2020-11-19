.. _sw-code-review:

Code Review
===========

*“It’s a simple 3-step process. Step one: Fix! Step two: It! Step three:
Fix it!” -- Oscar Rogers (Kenan Thompson), Saturday Night Live, 2/2009*

Code review is the last line of defense between a mistake that the IDAES
team will see and a mistake the whole world will see. In the case of
that mistake being a leak of proprietary information, the entire project
is jeopardized, so we need to take this process seriously.

Summary
-------

.. warning:: This section is an incomplete set of notes

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

Automated Checks
~~~~~~~~~~~~~~~~
The first level of code review is a set of automated checks that *must* pass
before the code is ready for people to review it. These checks will run
on the initiation of a :ref:`pull request <devterm-pr>` and on every new commit to that pull
request that is pushed to Github (thus the name “continuous
integration”).

The “continuous integration” of the code is hosted by an online service
– we use `CircleCI <https://circleci.com>`_ -- that can automatically
rerun the tests after every change (in this case, every new Pull Request
or update to the code in an existing Pull Request) and report the
results back to Github for display in the web pages. This status
information can then be used as an automatic gatekeeper on whether the
code can be merged into the main branch – if tests fail, then no merge
is allowed. Following this procedure, it is not possible for the main
branch to ever be failing its own tests.
