Scaling Toolbox
===============

Introduction
------------

The numerical scaling of a model is critical to its robustness and tractability, and can mean the difference between finding a solution and a solver failure. Scaling of non-linear models is a key step in the model development and application workflow, and often requires a significant amount of time and effort. The IDAES Scaling Toolbox aims to assist users in this task.

Key Points
----------

* Scaling is critical for model performance
* Scaling is not a one-off process, but something that needs to be reassessed before each solver call
* There are three aspects to scaling which are equally important
* Model scaling is more about avoiding poor scaling than finding “good” scaling
* Order-of-magnitude scaling factors are often better than exact values
* A number of tools are available to assist modelers with this task

Topics
------

.. toctree::
    :maxdepth: 1

    scaling_theory
    diagnosing_scaling_issues
    scaling_toolbox
    applying_scaling_tools
    scaling_workflow

