
Flexibility Analysis Overview
=============================

.. contents::
    :depth: 3
    :local:

Introduction
------------

The flexibility analysis module within IDAES provides a framework for
evaluating how well a given system performs with respect to a set of
uncertain parameters. Two methods are provided. The flexibility (or feasibility) test
(FT) can be used to determine if a set of performance constraints can
be satisfied for any realization of uncertain parameters. The
flexibility index (FI) can be used to quantify the size of the
uncertainty region for which the performance constraints can be
satisfied [Grossmann1987]_.

The FT is given by

.. _FT:

.. math::

   \phi(\underline{\theta}, \overline{\theta}) = &\max_{\theta} \min_{z} u \\
   & s.t. \\
   & g_{j}(x,z,\theta) \leq u \\
   & h(x,z,\theta) = 0 \\
   & \underline{\theta} \leq \theta \leq \overline{\theta}

where the uncertain parameters are given by :math:`\theta`, :math:`z`
are the controls, :math:`u` is the maximum constraint violation,
:math:`g_j` are the performance constraints, and :math:`h` are the
constraints which represent physics (e.g., mass balances). Note that
the dimension of :math:`x` must match the dimension of :math:`h`. In
other words, if :math:`\theta` and :math:`z` are fixed, then :math:`h`
should completely determine :math:`x`. If
:math:`\phi(\underline{\theta}, \overline{\theta})` is less than or
equal to zero, then the FT passes indicating that, for any
:math:`\theta` between :math:`\underline{\theta}` and
:math:`\overline{\theta}`, there exists a :math:`z` such that
:math:`g_j(x, z, \theta) \leq 0` and :math:`h(x, z, \theta) = 0`. If
:math:`\phi(\underline{\theta}, \overline{\theta})` is greater than
zero, then the FT fails indicating that the performance constraints
cannot be satisfied for at least one value of :math:`\theta` between
:math:`\underline{\theta}` and :math:`\overline{\theta}`. Also note
that this formulation assumes that :math:`h(x,z,\theta) = 0` can be
satisfied for any :math:`\theta` between :math:`\underline{\theta}`
and :math:`\overline{\theta}`.

The FI is given by 

.. math::

   \psi(\theta^{N}, \Delta \theta) = &\max \delta \\
   & s.t. \\
   & \phi(\underline{\theta}, \overline{\theta}) \leq 0 \\
   & \underline{\theta} = \theta^{N} - \delta \Delta \theta \\
   & \overline{\theta} = \theta^{N} + \delta \Delta \theta

where :math:`\theta^{N}` is a "nominal" point. The goal of the FI is
to find the largest region around this nominal point for which the
performance constraints can be satisfied. As written, the FI searches
for the largest hyperrectangle, but other shapes (e.g., ellipses) can
be used. The hyperrectangle is all that is currently supported in the
flexibility analysis module in IDAES. Typically, :math:`\delta` is
bounded between 0 and 1.

Flexibility Test Solution Methods
---------------------------------

Vertex Enumeration
^^^^^^^^^^^^^^^^^^

Vertex enumeration solves the inner minimization problem

.. _innerProblem:

.. math::

   & \min_{z} u \\
   & s.t. \\
   & g_{j}(x,z,\theta) \leq u \\
   & h(x,z,\theta) = 0

at each vertex of the hyperrectangle :math:`[\underline{\theta}, \overline{\theta}]`. For certain problem types (e.g., linear), the solution to :ref:`FT<FT>` is guaranteed to be at one of the vertices of this hyperrectangle [Swaney1985]_. For other problem types, vertex enumeration is only a heuristic that may not find the value of :math:`\theta` that maximizes the violation of the performance constraints.

Active Constraint Method
^^^^^^^^^^^^^^^^^^^^^^^^

The active constraint method converts the bilevel problem given by :ref:`FT<FT>` to a single-level problem by formulating the KKT conditions of the :ref:`inner minimization problem<innerProblem>` and embedding them as constraints in the outer problem [Grossmann1987]_. Note that this method assumes the Haar Conditions hold and a constraint is added enforcing that the number of active inequalities (i.e., :math:`g_{j}(x,z,\theta) = u`) is equal to the number of controls plus one. If the resulting single-level problem is solved to global optimality, this method will be conservative because it will find the "worst" local minima of the inner minimization problem.

Sampling
^^^^^^^^

Sampling is similar to `Vertex Enumeration`_ except that the :ref:`inner minimization problem<innerProblem>` is solved at random samples of :math:`\theta` rather than only at the vertices of :math:`[\underline{\theta}, \overline{\theta}]`. This can be useful for nonlinear problems but still may miss the worst-case :math:`\theta`.

Decision Rules
^^^^^^^^^^^^^^

Decision rules can be used to convert the bilevel problem given by :ref:`FT<FT>` to a single-level problem by removing all degrees of freedom of the inner problem with a control policy. Suppose we have a decision rule give by :math:`z = d(\theta)`. Because the only degrees of freedom in the inner problem are :math:`z`, the :ref:`FT<FT>` may be reformulated as

.. math::

   & \max_{\theta} \overline{u} \\
   & s.t. \\
   & g_{j}(x,z,\theta) = u_{j} \\
   & h(x,z,\theta) = 0 \\
   & \overline{u} = \sum u_{j} y_{j} \\
   & \sum y_{j} = 1 \\
   & z = d(\theta) \\
   & \underline{\theta} \leq \theta \leq \overline{\theta}

Currently, the two types of decision rules supported are linear decision rules and neural network decision rules with ReLU activation functions. Because the decision rules result in suboptimal values of :math:`z`, this method is conservative.

Flexibility Index Solution Methods
----------------------------------

Bisection Method
^^^^^^^^^^^^^^^^

The bisection method simply uses bisection to find the :math:`\delta` such that :math:`\phi(\underline{\theta}, \overline{\theta}) = 0` (:math:`\phi(\underline{\theta}, \overline{\theta})` is monotonically increasing with :math:`\delta`). Each subproblem solves the :ref:`FT<FT>` using one of the methods described above.

Usage
-----

The flexibility analysis module within IDAES provides two primary
functions. The first is
:meth:`solve_flextest<idaes.apps.flexibility_analysis.solve_flextest>`. The
:class:`FlexTestConfig<idaes.apps.flexibility_analysis.FlexTestConfig>`
specifies how the flexibility test should be solved.

.. [Grossmann1987] Grossmann, Ignacio E., and
                   Christodoulos A. Floudas. "Active constraint
                   strategy for flexibility analysis in chemical
                   processes." Computers & Chemical Engineering 11.6
                   (1987): 675-693.

.. [Swaney1985] Swaney, Ross Edward, and Ignacio E. Grossmann. 
                "An index for operational flexibility in chemical 
                process design. Part I: Formulation and theory." 
                AIChE Journal 31.4 (1985): 621-630.
