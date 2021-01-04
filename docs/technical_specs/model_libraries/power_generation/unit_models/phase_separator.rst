HelmPhaseSeparator Model
========================

.. index::
  pair: idaes.power_generation.unit_models.helm.phase_separator;HelmPhaseSeparator


.. currentmodule:: idaes.power_generation.unit_models.helm.phase_separator

Introduction
------------

The HelmPhaseSeparator model consists of a simple phase separator to be used only with the Helmholtz equation of state. 
The two-phase mixture at the inlet is separated into the vapor and liquid streams at the two corresponding outlets. This simple unit includes one state block (mixed_state) for the inlet, and two state blocks, one for liquid (liq_state) and the other for vapor (vap_state). 
Note that this water-specific flash model replaces IDAES' generic flash unit operation model.

Model inputs:

* mixed_state, variables (flow_mol, enth_mol, and pressure), port name = inlet

Model Outputs:

* liq_state, variables (flow_mol, enth_mol, and pressure), port name = liq_outlet
* vap_state, variables (flow_mol, enth_mol, and pressure), port name = vap_outlet


Degrees of Freedom
------------------

The HelmPhaseSeparator model consist of nine variables and six constraints. By fixing the inlet state (self.inlet.flow_mol,self.inlet.flow_mol, and self.inlet.flow_mol) or three degrees
of freedom, the system will be fully specified.


Variables
---------

============= ================== =========== ===============
Variable      Symbol             Index Sets  Doc
============= ================== =========== ===============
flow_mol      :math:`F`          time        molar flowrate
enth_mol      :math:`h`          time        molar enthalpy
pressure      :math:`P`          time        pressure
============= ================== =========== ===============


Constraints
-----------

The phase separator model uses the IAPWS95 property package to calculate the vapor fraction and enthalpies of the vapor and liquid phases at the inlet of the unit.
The flowrates of the vap_outlet and liq_outlet streams are calculate as the products of the inlet flow rate and corresponding phase fractions for vapor and liquid, respectively.
The enthalpies of the vapor and liquid phases in the inlet stream are assigned to the enthalpies of the vap_outlet and the liq_outlet streams, respectively.
The pressure of the two outlet streams are identical to that of the inlet stream.

Material Balances:
Vapor State:

.. math::
  flow\_mol_{mixed\_state}*vapor\_frac_{mixed\_state} = flow\_mol_{vap_outlet}

Liquid State:

.. math::
  flow\_mol_{mixed\_state}*(1 - vapor\_frac_{mixed\_state}) = flow\_mol_{liq_outlet}
  
Energy Balances:

.. math::
  enth\_mol\_phase_{mixed\_state}_[Vap] = enth\_mol_{vap\_state}

.. math::
  enth\_mol\_phase_{mixed\_state}_[Liq] = enth\_mol_{liq\_state}

Momentum Balances:

.. math::
  pressure_{mixed\_state}_[Liq] = pressure_{liq\_state} = pressure_{vap\_state}
