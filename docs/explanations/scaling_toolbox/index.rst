Scaling Toolbox
===============

.. contents::
    :depth: 3
    :local:

.. toctree::
    :maxdepth: 1


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
* A number of tools are avaialble to assist modelers with this task


What is Model Scaling?
----------------------

The primary goal of model scaling is to minimize the effects of numerical round-off errors (due to machine precision) and to ensure that the model equations are well-posed and non-singular (to numerical tolerance). In reality, there are actually three inter-related scaling concepts which contribute to the overall performance and tractability of a model.

Types of Scaling
''''''''''''''''

* Variable scaling determines what sort of change in a variable is “significant” or not. This factors into things like, for example, how an interior point method (such as IPOPT) will behave (i.e., how far into the interior of the feasible region a variable should be) and step size determination.

* Constraint residual scaling refers to the magnitude of the residual of all constraints in the model. This is important as this is required to determine whether or not a model has converged, thus it is important that a residual equal to the solver tolerance does not significantly alter the solution to the model.

  * E.g. consider a constraint A=B. If the magnitude of A and B are 1e6 and the solver tolerance is 1e-6, this means that A and B need to be solved to 12 significant figures of precision to converge the constraint (which may be unnecessarily onerous). Similarly, if A and B were of order 1e-6, then they will be considered equal if A=1e-6 and B=0; i.e. ±100% error.
  
* Jacobian scaling refers to the overall scaling and conditioning of the problem Jacobian. This is important as this determines the ability of the solver to find a path from the initial state to the final solution. It is important to ensure that the Jacobian is well-conditioned both in order to get good numerical behavior from the solver and also ensure that floating point round-off error does not result in large changes to the variables at the solution.

These aspects are not always complimentary, and there are often cases where improving one aspect of the model scaling can negatively affect another aspect (e.g., focusing too much on Jacobian condition number can often result in poor constraint residual tolerances). Additionally, each aspect affects different parts of the solver routine; Jacobian scaling is important for determining the step direction and size at each iteration, whilst constraint residual scaling is important for determining if and when a solution is found. The relative importance of each of these depends on the solver being used (e.g., for IPOPT constraint residual scaling is more important than Jacobian scaling as IPOPT does its own internal scaling of the Jacobian), however in general all of these will have an impact upon the solver behavior and users should endeavor to have good scaling for all aspects of the model. 

What Determines Good Scaling?
'''''''''''''''''''''''''''''

Due to the three different aspects of scaling and the fact that they can sometimes compete with each other, it is hard (if not impossible) to define a single metric for what defines “good” scaling. However, common indicators of bad scaling are:

* Variables with extremely large or small scaled magnitudes (>1e6 or <1e-6) – you should generally aim to have most variables have scaled magnitudes between 1 and 10.

  * An important exception is any variable that can pass through zero, like enthalpy or bidirectional fluxes. These variables should be scaled based on “significant variation” (think something like a standard deviation) of the variable throughout different conditions of the model.
  * Zero flux boundary conditions should be scaled like the other fluxes in the model.

* Constraints with 1 or more very large terms – as each constraint must be solved to solver tolerance in order to meet the convergence criteria, it is important that change on the order of magnitude of the solver tolerance is meaningful for each constraint term.
* Constraints in which all terms are very small – having a few small terms combined with moderate terms is often acceptable (e.g., some terms might go to zero under certain circumstances), but constraints where all terms are of small magnitude is indicative of poor scaling.
* Constraint with terms that are of a significantly different orders of magnitude. This means that a change that is significant in the small terms may be negligible for the larger terms (i.e., you can make significant changes in the small term with no effect on the larger term), whilst a change in the larger terms may result in huge changes in the smaller terms (i.e., a change in the final significant figure of the large term might result in orders of magnitude change in the smaller term, thus causing a cascade of changes through your model).
* Any time a difference is taken between several large numbers to obtain a small number, the model is at risk of catastrophic cancellation and ill-conditioning
* Very large Jacobian Condition Number (>1e10). The condition number represents a worst-case error amplification in the problem. For example, with a condition number of 1e10 then a round-off error due to machine precision (~1e-16) could result in a change in model state of 1e-6 (= 1e10 x 1e-16); note that most solvers use a default tolerance of 1e-6 to 1e8, thus with a condition number of 1e10 or greater, there is a risk that your solution may be dominated by round-off errors.

The :ref:`IDAES Diagnostics Toolbox<explanations/model_diagnostics/index:Model Diagnostics Workflow>` contains methods to identify and report these issues, and these should be used to help identify potential scaling issues.

Approaches to Scaling
---------------------

Ultimately, the modeler is the one with the most understanding of the system being modeled and the expected solution. In an ideal world, the modeler would have sufficient time and knowledge to develop a set of scaling factors customized for their specific use case. However, process models are generally very large and modelers often lack the time and understanding to develop tailored scaling routines, thus  there is a need for tools to support the modeler in this.

In general, there are two types of scaling routines that can be developed for a model:
* those that depend on the current state of the model (often referred to as auto-scalers), and thus require an initial solution, and
* those that make use of “best-guess” measures and heuristics for scaling and can thus be used before a model is initialized but require development by an expert modeler.

In either case, scaling depends heavily on input from the end user. Whilst it might be possible to get a good idea of the scaling for the current state of the model (assuming the model has been initialized), this may not be indicative of the scaling at the final solution (generally, we want to solve or optimize for some unknown solution). Thus, the modeler needs to provide as much information as they can of the expected scaling near the solution point as possible (based on engineering knowledge or intuition, or solution for similar problems in the past). It is important to note that scaling does not need to be exact – order of magnitude approximations are often as good (or better) than precise values – and the aim is less about providing “good” scaling as it is about avoiding bad scaling.

Auto-Scalers (Current State Scaling)
''''''''''''''''''''''''''''''''''''

Pros:

* fully automated, can scale an entire model with one command
* useful for getting initial values for scaling factors, or filling-in missing scaling factors

Cons:

* require an initial solution, and thus not useful for pre-initialization scaling
* consider only the current model state, and often overscale the problem

For models with an initial (ideally feasible) solution, full information on the state of the model and the Jacobian is available. A number of approaches are available that can take this information and generate scaling factors for the entire model in order to meet some overall model characteristic (e.g., minimizing the norm of the Jacobian matrix). These have the advantage of requiring minimal input from the modeler and being able to scale an entire model in one go. However, these approaches require an initial solution (and thus are not useful for pre-initialization scaling) and as they consider only a single characteristic metric calculated at the current model state, they can often over scale the model and may not provide the best performance.

A suite of autoscaler methods are available as part of the IDAES Scaling Toolbox through the :ref:`AutoScaler Class<reference_guides/scaling/autoscaler:AutoScaler Class>`.

Custom Scalers (Best-Guess Scaling)
'''''''''''''''''''''''''''''''''''

Pros:

* independent of model state, and can thus be used on uninitialized models

Cons:

* dependent on developer skill and foresight, and may not give good results for all cases

The alternative to model-state based scaling is to develop customized scaling routines for individual models which take into account the model structure and behavior to estimate scaling factors. These routines are generally written by the original model developer, and thus depend heavily on the skill and foresight of the developer. On the other hand, as these routines depend on knowledge of the model structure rather than the model state, these routines can be applied to uninitialized models (given sufficient estimates of a few key variables).

A suite of methods to assist with developing custom scaling routines are available as part of the IDAES Scaling Toolbox through the :ref:`CustomScalerBase Class<reference_guides/scaling/custom_scaler:CustomScalerBase Class>`.

Utility Functions
'''''''''''''''''

The IDAES Scaling Toolbox also contains a number of utility functions for working with and reporting model scaling. More details on these functions can be :ref:`found here<reference_guides/scaling/scaling_utils:Scaling Utility Functions>`.

How to Use the IDAES Scaling Toolbox
------------------------------------
Show general workflow

Mention model scalers
