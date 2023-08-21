Mathematical Utility Methods
============================

.. contents:: Contents
    :depth: 2

There are a number of useful mathematical operations that either exhibit singularities or non-smooth behavior that can cause issues when used in equation-oriented models. This library contains smooth formulations for some of these operations which can be used to approximate these operations in models.

Users should note that these formulations involve smoothing parameters which should be tuned for their particular system. Larger values generally result in smoother behavior but greater error.

.. automodule:: idaes.core.util.math
  :members:

