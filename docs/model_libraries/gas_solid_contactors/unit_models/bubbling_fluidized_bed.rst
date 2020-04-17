Bubbling Fluidized Bed Reactor
==================

The IDAES Bubbling Fluidized Bed Reactor (BFB) model represents a unit operation
 where two material streams – a solid phase and a gas phase – pass through a linear
  reactor vessel while undergoing chemical reaction(s).
   
The BFB model is represented as a 1-D axially discretized model
with two phases (gas and solid), and two regions (bubble and emulsion)
resulting in 3 control volume_1D blocks (bubble, gas_emulsion solid_emulsion).
The model captures the gas-solid interaction between both phases and regions
through reaction, mass and heat transfer.

This model is a simplified version of the Kunii and Levenspiel 3-region model with the assumption
that cloud-wake region effects are negligble. Other model assumptions are

Assumptions:
Gas emulsion is at minimum fluidization conditions
Gas feeds into emulsion region before the excess enters into the bubble region

The BFB model equations are derived from:
A. Lee, D.C. Miller. A one-dimensional (1-D) three-region model for a bubbling
fluidized-bed Adsorber, Ind. Eng. Chem. Res. 52 (2013) 469–484.

Degrees of Freedom
------------------

BFBs generally have at least 2 (or more) degrees of freedom. Typically fixed variables are reactor length and diameter.

Model Structure
---------------

The core BFB unit model consists of two inlet ports (named gas_inlet and solid_inlet), two outlet ports
 (named gas_outlet and solid_outlet) and three ControlVolume1DBlock Blocks
 (named bubble_region, gas_emulsion_region and solid_emulsion_region).

Variables
---------

ControlVolume1DBlock Blocks add the following additional Variables:

====================== ======= ===============================================================
Variable               Name    Notes
====================== ======= ===============================================================
:math:`L`              length  Reference to control_volume.length
:math:`A`              area    Reference to control_volume.area
:math:`V`              volume  Reference to control_volume.volume
:math:`Q_{t,x}`        heat    Only if has_heat_transfer = True, reference to holdup.heat
:math:`\Delta P_{t,x}` deltaP  Only if has_pressure_change = True, reference to holdup.deltaP
====================== ======= ===============================================================

Constraints
-----------

Geometric Constraints
---------------------

Area of orifice:

.. math:: A_{o} = \frac{1}{n_{o}}

Bed cross-sectional area:

.. math:: A = \pi \frac{D^{2}}{4}

Area of bubble region:

.. math:: A_{b,t,x} = \delta_{t,x} A

Area of gas emulsion region:

.. math:: A_{ge,t,x} = \delta_{e,t,x} \varepsilon_{e,t,x} A

Area of solid emulsion region:

.. math:: A_{se,t,x} = \delta_{e,t,x} {\left(1 - \varepsilon_{e,t,x} \right)} A

Length of bubble region:

.. math:: L_{b} = L

Length of gas emulsion region:

.. math:: L_{ge} = L

Length of solid emulsion region:

.. math:: L_{se} = L

Hydrodynamic Constraints
------------------------
Emulsion region volume fraction:

.. math:: \delta_{e,t,x} = 1 -\delta_{t,x}

Average cross-sectional voidage:

.. math:: \varepsilon_{t,x} = 1 - \left(1 - \varepsilon_{e,t,x} \right) \left(1 - \delta_{t,x} \right)

Emulsion region voidage:

.. math:: \varepsilon_{e,t,x} = \varepsilon_{mf}

Bubble growth coefficient:

.. math:: \gamma_{t,x} = \frac{0.0256}{v_{mf}} {\left(\frac{D}{g} \right)}^{0.5}

Maximum bubble diameter:

.. math:: d_{b,m,t,x}^{5}g = 2.59^{5} {\left([v_{g,t,x} - v_{ge,t,x}]A \right)}^{2}

Bubble diameter (gas inlet, x = 0):

.. math:: d_{b,t,x} = 1.38g^{-0.2} {\left([v_{g,t,x} - v_{ge,t,x}]A_{o} \right)}^{0.4}

Bubble diameter (x > 0):

.. math:: \frac{d_{b,t,x}}{ dx } = \frac{0.3}{D} L {\left(d_{b,m,t,x} - d_{b,t,x} - \gamma_{t,x}{\left(D_{t} d_{b,t,x}\right)}^{0.5}\right)}

Bubble rise velocity:

.. math:: v_{br,t,x}^{2} = 0.711^{2} g d_{b,t,x}

Bubble velocity:

.. math:: v_{b,t,x} = v_{g,t,x} - v_{mf} + v_{br,t,x}

Emulsion region gas velocity:

.. math:: v_{ge,t,x} = v_{mf}

Superficial gas velocity:

.. math:: v_{g,t,x} = v_{b,t,x} \delta_{t,x} + v_{ge,t,x}\delta_{e,t,x}

Gas emulst,xon pressure drop:

if 'has_pressure_change' is 'True':

.. math:: \Delta P_{t,x} = - U_{1} g {\left(1 - \varepsilon_{e,t,x} \right)} \rho_{mass,se,t,x}

Mass Transfer Constraints
-------------------------

Bubble to emulsion gas mass transfer coefficient:

.. math:: K_{be,t,x,j} d_{b,t,x}^{1.25} = 5.94 v_{mf} d_{b,t,x}^{0.25} + 5.85 D_{vap,t,x,j}^{0.5} g^{0.25}

Bulk gas mass transfer:

.. math:: K_{gbulkc,t,x,j} = 6 K_{d} \delta_{t,x} A {\left(C_{ge,total,t,x} - C_{b,total,t,x} \right)} d_{b,t,x} y_{ge,t,x,j}

Heat Transfer Constraints
-------------------------

Bubble to emulsion gas heat transfer coefficient:

.. math:: H_{be,t,x,j} d_{b,t,x}^{1.25} = 4.5 v_{mf} c_{p\_vap,b,t,x} C_{b,total,t,x} d_{b,t,x}^{0.25} + 5.85 {\left(k_{vap,b,t,x}  C_{b,total,t,x} c_{p\_vap,b,t,x} \right)}^{0.5} g^{0.25}

Convective heat transfer coefficient:

.. math:: h_{tc,t,x} d_{p} = 0.03 k_{vap,e,t,x} {\left(v_{ge,t,x} d_{p} \frac{C_{ge,total,t,x}}{\mu_{vap,t,x}} \right)}^{1.3}

Emulsion region gas-solids convective heat transfer:

.. math:: h_{t\_gs,t,x} d_{p} = 6 \delta_{e,t,x} {\left(1 - \varepsilon_{e,t,x} \right)} h_{tc,t,x} {\left(T_{ge,t,x} - T_{se,t,x} \right)}

Bulk gas heat transfer:

.. math:: H_{gbulk,t,x} =  K_{d} \delta_{t,x} A {\left(C_{ge,total,t,x} - C_{b,total,t,x} \right)} d_{b,t,x} h_{vap,e,t,x}

Mass and heat transfer terms in control volumes
-----------------------------------------------

Bubble mass transfer '(p=vap)':

.. math:: M_{tr,b,t,x,p,j} = K_{gbulkc,t,x,j} - A_{b,t,x} K_{be,t,x,j} {\left(C_{b,total,t,x} - C_{ge,total,t,x} \right)}

Gas emulsion mass transfer '(p=vap)':

.. math:: M_{tr,ge,t,x,p,j} = - K_{gbulkc,t,x,j} + A_{b,t,x} K_{be,t,x,j} {\left(C_{b,total,t,x} - C_{ge,total,t,x} \right)} + r_{hetero,ge,t,x,j}

Solid emulsion mass transfer '(p=sol)':

.. math:: M_{tr,se,t,x,p,j} = 0

if 'energy_balance_type' is 'EnergyBalanceType.none':

Bubble heat transfer:

.. math:: H_{tr, b, t,x} = H_{gbulk,t,x} - A_{b,t,x} H_{be,t,x,j} {\left(T_{b,t,x} - T_{ge,t,x} \right)} 

Gas emulsion heat transfer:

.. math:: H_{tr, ge, t,x} = - H_{gbulk,t,x} + A_{b,t,x} H_{be,t,x,j} {\left(T_{b,t,x} - T_{ge,t,x} \right)} - h_{t\_gs,t,x} A

Solid emulsion heat transfer:

.. math:: H_{tr, se, t,x} =  h_{t\_gs,t,x} A

Reaction constraints
--------------------

if 'homogeneous reaction package' is not 'None':
 
Bubble rate reaction extent:

.. math:: r_{ext,b,t,x,r} = A_{b,t,x} r_{b,t,x,r}

Gas emulsion rate reaction extent:

.. math:: r_{ext,ge,t,x,r} = A_{ge,t,x} r_{ge,t,x,r}

if 'heterogeneous reaction package' is not 'None':

Solid emulsion rate reaction extent:

.. math:: r_{ext,se,t,x,r} = A_{se,t,x} r_{se,t,x,r}

Gas emulsion heterogeneous rate reaction extent:

.. math:: r_{hetero,ge,t,x,j} = A_{se,t,x} \sum_{r}^{reactions} {s_{j,r} r_{se,t,x,r}}

Flowrate constraints
--------------------

Bubble gas flowrate:

.. math:: F_{mol,b,t,x} = A \delta_{t,x} v_{b,t,x} C_{b,total,t,x}

Emulsion gas flowrate:

.. math:: F_{mol,ge,t,x} = A v_{ge,t,x} C_{ge,total,t,x}

Inlet boundary conditions
-------------------------

Gas emulsion pressure at inlet:

.. math:: P_{ge,t,0} = P_{g,t,inlet} - \Delta P_or

Total gas balance at inlet:

.. math:: F_{mol,g,inlet,t} = F_{mol,b,t,0} + F_{mol,ge,t,0}

Superficial gas velocity at inlet:

.. math:: v_{g,t,0} = \frac{F_{mol,g,inlet,t}}{C_{ge,total,t,0} A} 

Bubble mole fraction at inlet:

.. math:: y_{g,inlet,t,j} = y_{b,t,0,j}

Gas emulsion mole fraction at inlet:

.. math:: y_{g,inlet,t,j} = y_{ge,t,0,j}

Bubble mole fraction at inlet:

.. math:: y_{g,inlet,t,j} = y_{b,t,0,j}

Solid emulsion mass flow at inlet:

if 'flow_type' is 'co_current' x = 0 else if 'flow_type' is 'counter_current' x = 1:

.. math:: F_{mass,s,inlet,t} = F_{mass,se,t,x}

Solid emulsion mass fraction at inlet:

if 'flow_type' is 'co_current' x = 0 else if 'flow_type' is 'counter_current' x = 1:

.. math:: x_{s,inlet,t} = x_{se,t,x}

if 'energy_balance_type' is not 'EnergyBalanceType.none':

Gas inlet energy balance:

.. math:: H_{g,inlet,t} = H_{b,t,0} + H_{ge,t,0}

Gas emulsion temperature at inlet:

.. math:: T_{g,inlet,t} = H_{b,t,0} + H_{ge,t,0}

if 'flow_type' is 'co_current' x = 0 else if 'flow_type' is 'counter_current' x = 1:

Solid inlet energy balance:

.. math:: H_{s,inlet,t} = H_{se,t,x}

Outlet boundary conditions
-------------------------

Gas emulsion pressure at outlet:

.. math:: P_{g,inlet,t} = P_{ge,t,1} 

Total gas balance at outlet:

.. math:: F_{mol,g,outlet,t} = F_{mol,b,t,1} + F_{mol,ge,t,1}

Solid outlet material balance:

if 'flow_type' is 'co_current' x = 1 else if 'flow_type' is 'counter_current' x = 0:

.. math:: F_{mass,s,outlet,t} = F_{mass,se,t,x}

if 'energy_balance_type' is not 'EnergyBalanceType.none':

Gas outlet energy balance:

.. math:: H_{g,outlet,t} = H_{b,t,1} + H_{ge,t,1}

Solid outlet energy balance:

if 'flow_type' is 'co_current' x = 1 else if 'flow_type' is 'counter_current' x = 0:

.. math:: H_{g,outlet,t} = H_{ge,t,x}

                         
Initialization
--------------

The initialization method for this model uses a multi-step, sequential, hierarchical initialization approach.

MBR Class
---------

.. module:: idaes.gas_solid_contactors.unit_models.moving_bed

.. autoclass:: MBR
    :members:

MBRData Class
-------------

.. autoclass:: MBRData
    :members:
