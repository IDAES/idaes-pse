
Flexibility Analysis Overview
=============================

The flexibility analysis module within IDAES provides a framework for
evaluating how well a given system performs with respect to a set of
uncertain parameters. Two methods are provided. The flexibility (or feasibility) test
(FT) can be used to determine if a set of performance constraints can
be satisfied for any realization of uncertain parameters. The
flexibility index (FI) can be used to quantify the size of the
uncertainty region for which the performance constraints can be
satisfied [Grossmann1987]_.

The FT is given by

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
