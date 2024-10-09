Scaling Toolbox
===============

.. contents:: :local:
    :depth: 1

Approaches to Scaling
---------------------

Ultimately, the modeler is the one with the most understanding of the system being modeled and the expected solution. In an ideal world, the modeler would have sufficient time and knowledge to develop a set of scaling factors customized for their specific use case. However, process models are generally very large and modelers often lack the time and understanding to develop tailored scaling routines, thus there is a need for tools to support the modeler in this.

In general, there are two types of scaling routines that can be developed for a model:

* those that depend on the current state of the model (often referred to as auto-scalers), and thus require an initial solution, and
* those that make use of “best-guess” measures and heuristics for scaling and can thus be used before a model is initialized but require development by an expert modeler.

In either case, scaling depends heavily on input from the end user. Whilst it might be possible to get a good idea of the scaling for the current state of the model (assuming the model has been initialized), this may not be indicative of the scaling at the final solution (generally, we want to solve or optimize for some unknown solution). Thus, the modeler needs to provide as much information as they can of the expected scaling near the solution point as possible (based on engineering knowledge or intuition, or solutions for similar problems in the past). It is important to note that scaling does not need to be exact – order of magnitude approximations are often as good (or better) than precise values – and the aim is less about providing “good” scaling as it is about avoiding bad scaling.

Auto-Scalers (Current State Scaling)
------------------------------------

Pros:

* fully automated, can scale an entire model with one command
* applicable to any model
* useful for getting initial values for scaling factors, or filling-in missing scaling factors

Cons:

* require an initial solution, and thus not useful for pre-initialization scaling
* consider only the current model state, and often overscale the problem

For models with an initial (ideally feasible) solution, full information on the state of the model and the Jacobian is available. A number of approaches are available that can take this information and generate scaling factors for the entire model in order to meet some overall model characteristic (e.g., minimizing the norm of the Jacobian matrix). These have the advantage of requiring minimal input from the modeler and being able to scale an entire model in one go. However, these approaches require an initial solution (and thus are not useful for pre-initialization scaling) and as they consider only a single characteristic metric calculated at the current model state, they can often over scale the model and may not provide the best performance.

A suite of autoscaler methods is available as part of the IDAES Scaling Toolbox through the :ref:`AutoScaler Class<reference_guides/scaling/autoscaler:AutoScaler Class>`.

Custom Scalers (Best-Guess Scaling)
-----------------------------------

Pros:

* independent of model state, and can thus be used on uninitialized models

Cons:

* specific to a given model type - Custom Scalers developed for one model are not applicable to other models
* dependent on developer skill and foresight, and may not give good results for all cases

The alternative to model-state based scaling is to develop customized scaling routines for individual models which take into account the model structure and behavior to estimate scaling factors. These routines are generally written by the original model developer, and thus depend heavily on the skill and foresight of the developer. On the other hand, as these routines depend on knowledge of the model structure rather than the model state, these routines can be applied to uninitialized models (given sufficient estimates of a few key variables).

A suite of methods to assist with developing custom scaling routines are available as part of the IDAES Scaling Toolbox through the :ref:`CustomScalerBase Class<reference_guides/scaling/custom_scaler:CustomScalerBase Class>`. Many models in the core IDAES libraries will have custom Scalers available for use - see the documentation of individual models for these.

Utility Functions
-----------------

The IDAES Scaling Toolbox also contains a number of utility functions for working with and reporting model scaling. More details on these functions can be :ref:`found here<reference_guides/scaling/scaling_utils:Scaling Utility Functions>`.

