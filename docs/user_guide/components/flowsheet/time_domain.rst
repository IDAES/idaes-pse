Time Domain
===========

Time domain is an essential component of the IDAES framework. When a user first declares a 
Flowsheet model a time domain is created, the form of which depends on whether the Flowsheet 
is declared to be dynamic or steady-state 
(see :ref:`FlowsheetBlock <technical_specs/core/flowsheet_block:Flowsheet Block>`). 
In situations where the user makes use of nested flowsheets, each sub-flowsheet refers to its 
parent flowsheet for the time domain.

Different models may handle the time domain differently, but in general all IDAES models refer 
to the time domain of their parent flowsheet. The only exception to this are blocks associated 
with Property calculations. PropertyBlocks (i.e. StateBlocks and ReactionBlocks) represent the state of the material at a single point 
in space and time, and thus do not contain the time domain. Instead, PropertyBlocks are indexed 
by time (and space where applicable) - i.e. there is a separate StateBlock for each point in 
time. The user should keep this in mind when working with IDAES models, as it is important for 
understanding where the time index appears within a model.

In order to facilitate referencing of the time domain, all Flowsheet objects have a `time` 
configuration argument which is a reference to the time domain for that flowsheet. All IDAES 
models contain a `flowsheet` method which returns the parent flowsheet object, thus a reference 
to the time domain can always be found using the following code: `flowsheet().config.time`.

Another important thing to note is that steady-state models do contain a time domain. While the
time domain for steady-stage models is a single point at time = 0.0, they still contain a 
reference to the time domain and the components (e.g. StateBlocks) are indexed by time.