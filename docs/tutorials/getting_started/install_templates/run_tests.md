Run the IDAES test suite to verify that your installation has been successful. :: 

    pytest --pyargs idaes -W ignore

You should see the tests run and all should pass to ensure the installation worked. You
may see some "Error" level log messages, but they are okay, and produced by tests for
error handling. The number of tests that failed and succeeded is reported at the end of the pytest
output. 

You can ask questions using the `Github Discussion Forum <https://github.com/IDAES/idaes-pse/discussions>`_ or report problems on the |github-issues|
(Please try to be specific about the command and the offending output.)
