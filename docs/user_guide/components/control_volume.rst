Control Volume
==============

A key feature of the IDAES modeling framework is the use of Control Volume Blocks. Control 
Volumes represent a volume of material over which material, energy and/or momentum balances 
can be performed. Control Volume Blocks contain methods to automate the task of writing common 
forms of these balance equations. Control Volume Blocks can also automate the creation of 
StateBlocks and ReactionBlocks associated with the control volume.

3. Control Volume Blocks - material, energy and momentum balances and the associated terms. 
These include:

    - balance equations
    - holdup volume
    - material and energy holdups; both variables and constraints
    - material and energy accumulation terms (Pyomo.dae handles the creation of the associated derivative constraints)
    - material generation terms (kinetic reactions, chemical and phase equilibrium, mass transfer)
    - extent of reaction terms and constraints relating these to the equivalent generation terms
    - phase fraction within the holdup volume and constrain on the sum of phase fractions
    - heat and work transfer terms
    - pressure change term
    - diffusion and conduction terms (where applicable) and associated constraints
    - Mixer and Splitter blocks for handling multiple inlets/outlets


