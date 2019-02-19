.. _tst-top:

Software Testing
================

Testing is essential to the process of creating software. The longer one programs,
the more the idea that "If it isn't tested, it probably doesn't work" makes sense.
This section outlines the different kinds of testing and how to write/use them in
the project.

.. _tst-unit:

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

.. _tst-coverage:

Code coverage
~~~~~~~~~~~~~
The “coverage” of the code refers to what percentage of
the code (“lines covered” divided by total lines) is executed by the
automated tests. This is important because passing automated tests is
only meaningful if the automated tests cover the majority of the code’s
behavior. This is not a perfect measure, of course, since simply
executing a line of code under one condition does not mean it would
execute correctly under all conditions. The code coverage is evaluated
locally and then integrated with Github through a tool called `Coveralls
<https://coveralls.io>`_.

