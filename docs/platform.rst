IDAES Integrated Platform
=========================
The IDAES integrated platform (IDAES-IP) supports the full process modeling lifecycle from conceptual design to dynamic optimization and
control within a single modeling environment. At the center of this platform is the Core Modeling Framework
(:term:`IDAES-CMF`) which leverages the open source, U.S. Department of Energy-funded extensible algebraic modeling
environment, :term:`Pyomo`.

Below is a diagram showing the components of the IDAES Integrated Platform (IDAES-IP).


.. image:: _images/IDAES-IP.png

This diagram shows the relationship between the major IDAES-IP components:

.. glossary::

    IDAES-IP
        IDAES integrated platform, described on this page

    IDAES-Core
        The software package that includes the Core Modeling Framework (IDAES-CMF); process, unit, and property model libraries;
        data management, artificial intelligence and uncertainty quantification tools; and graphical user interfaces.

    IDAES-CMF
        The IDAES-CMF is the center of the IDAES-Core. It extends :term:`Pyomo`'s block-based hierarchical modeling
        constructs to create a library of models for common process unit operations and thermophysical properties,
        along with a framework for the rapid development of process flowsheets.

    IDAES-AI
        Artificial intelligence and machine learning tools

    IDAES-UQ
        Tools supporting rigorous uncertainty quantification and optimization under uncertainty

    Graphical User Interfaces
        Tools for graphical interactive work, such as visualization of IDAES flowsheets

    Python programming environments
        Jupyter Notebook examples and extensions for interactive scripting in Python

    Data Management (IDAES-DMF)
        Data Management Framework (DMF) supporting provenance for IDAES workflows

    Pyomo
        Open source, U.S. Department of Energy-funded extensible algebraic modeling environment (AML).
        For more information, see the `Pyomo website <https://pyomo.org>`_.

    IDAES-Materials
    IDAES-Design
    IDAES-Enterprise
    IDAES-Operations
         Domain-specific tools for materials design, process design, enterprise-wide optimization, and control.

