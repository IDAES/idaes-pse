Phase Equilibria
================

The IDAES toolset contains methods for generating and displaying phase equilibria data.

Available Methods
-----------------

The phase equilibria methods include methods to calculate bubble and dew temperatures ( ``Txy_data()``) and include this data in a readable ``Class`` (``TXYDataClass``). This class can be used in the method ``build_txy_diagrams()`` to create T-x-y diagrams.

.. autofunction:: idaes.core.util.phase_equilibria.Txy_data

.. autofunction:: idaes.core.util.phase_equilibria.TXYDataClass

.. autofunction:: idaes.core.util.phase_equilibria.build_txy_diagrams

The methods also include ``Txy_diagram()`` which is a method that calculates the data and creates the plots automatically.

.. autofunction:: idaes.core.util.phase_equilibria.Txy_diagram

.. note::Currently these methods work with **FTPx** state definitions.
