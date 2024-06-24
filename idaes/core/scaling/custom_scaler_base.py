#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Base class for custom scaling routines.

Author: Andrew Lee
"""
from copy import copy

from pyomo.common.config import document_kwargs_from_configdict
from pyomo.environ import units

from idaes.core.scaling.scaling_base import CONFIG, ScalerBase
from idaes.core.scaling.util import get_scaling_factor
import idaes.logger as idaeslog

# Set up logger
_log = idaeslog.getLogger(__name__)

CSCONFIG = CONFIG()

DEFAULT_UNIT_SCALING = {
    "K": 1e-2,
    "Pa": 1e-5,
}


@document_kwargs_from_configdict(CSCONFIG)
class CustomScalerBase(ScalerBase):
    """
    Base class for custom scaling routines.

    """

    # Common data structures for default scaling
    # DEFAULT_SCALING_FACTORS = {"component_local_name": DEFAULT_SCALING}
    DEFAULT_SCALING_FACTORS = None

    # UNIT_SCALING_FACTORS = {"units": UNIT_BASED_SCALING}
    UNIT_SCALING_FACTORS = copy(DEFAULT_UNIT_SCALING)

    def __init__(self, **kwargs):
        self.config = CSCONFIG(kwargs)

        if self.DEFAULT_SCALING_FACTORS is not None:
            self.default_scaling_factors = copy(self.DEFAULT_SCALING_FACTORS)
        else:
            self.default_scaling_factors = {}

        if self.UNIT_SCALING_FACTORS is not None:
            self.unit_scaling_factors = copy(self.UNIT_SCALING_FACTORS)
        else:
            self.unit_scaling_factors = {}

    def scale_model(
        self,
        model,
        first_stage_fill_in: list = None,
        second_stage_fill_in: list = None,
        submodel_scalers: dict = None,
    ):
        """
        Default model scaling routine.

        This method performs a four-step scaling routine:

        1. Scale variables using variable_scaling_routine
        2. Perform first-stage scaling factor fill in using user provided method(s), called in order declared
        3. Scale constraints using constraint_scaling_routine
        4. Perform second-stage scaling factor fill in using user provided method(s), called in order declared

        Args:
            model - model to be scaled
            first_stage_fill_in - list of methods to use for first-stage scaling factor fill in
            second_stage_fill_in - list of methods to use for second-stage scaling factor fill in
            submodel_scalers - dict of Scalers to use for sub-models, keyed by submodel local name

        Returns:
            dict of additional scaling information
        """
        # Step 0: Verify model
        self.verify_model(model)

        # Step 1: Call variable scaling routine
        var_scaling = self.variable_scaling_routine(
            model, overwrite=self.config.overwrite, submodel_scalers=submodel_scalers
        )
        # Step 2: Call variable fill in
        if first_stage_fill_in is not None:
            for i in first_stage_fill_in:
                i(model)

        # Step 3: Call constraint scaling routine
        cons_scaling = self.constraint_scaling_routine(
            model, overwrite=self.config.overwrite, submodel_scalers=submodel_scalers
        )

        # Step 4: Call constraint fill in
        if second_stage_fill_in is not None:
            for i in second_stage_fill_in:
                i(model)

        # Step 5: Return scaling information for parent model
        scaling_data = {}
        if var_scaling is not None:
            scaling_data.update(var_scaling)
        if cons_scaling is not None:
            scaling_data.update(cons_scaling)
        return scaling_data

    def verify_model(self, model):
        """
        Verify that model to be scaled meets expectations.

        Derived classes must overload this method.

        Args:
            model - model to be scaled

        Raises:
            Exception if model does not match expected form
        """
        raise NotImplementedError(
            "Custom Scaler has not implemented a verify_model method."
        )

    def variable_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        """
        Routine to apply scaling factors to variables in model.

        Derived classes must overload this method.

        Args:
            model - model to be scaled
            overwrite - whether to overwrite existing scaling factors
            submodel_scalers - dict of Scalers to use for sub-models, keyed by submodel local name

        Returns:
            dict of additional scaling information
        """
        raise NotImplementedError(
            "Custom Scaler has not implemented a variable_scaling_routine method."
        )

    def constraint_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        """
        Routine to apply scaling factors to constraints in model.

        Derived classes must overload this method.

        Args:
            model - model to be scaled
            overwrite - whether to overwrite existing scaling factors
            submodel_scalers - dict of Scalers to use for sub-models, keyed by submodel local name

        Returns:
            dict of additional scaling information
        """
        raise NotImplementedError(
            "Custom Scaler has not implemented a constraint_scaling_routine method."
        )

    def get_default_scaling_factor(self, component):
        """
        Get scaling factor for component from dict of default values.

        Args:
            component - component to get default scaling factor for

        Returns:
            default scaling factor if it exists, else None
        """
        try:
            return self.default_scaling_factors[component.local_name]
        except KeyError:
            _log.debug(f"No default scaling factor found for {component.name}")
            return None

    # Common methods for variable scaling
    def scale_variable_by_component(
        self, target_variable, scaling_component, overwrite: bool = False
    ):
        """
        Set scaling factor for target_variable equal to that of scaling_component.

        Args:
            target_variable - variable to set scaling factor for
            scaling_component - component to use for scaling factor
            overwrite - whether to overwrite existing scaling factors

        Returns:
            None
        """
        sf = get_scaling_factor(scaling_component)

        if sf is not None:
            self.set_variable_scaling_factor(
                variable=target_variable, scaling_factor=sf, overwrite=overwrite
            )
        else:
            _log.debug(
                f"Could not set scaling factor for {target_variable.name}, "
                f"no scaling factor set for {scaling_component.name}"
            )

    def scale_variable_by_bounds(self, variable, overwrite: bool = False):
        """
        Set scaling factor for variable based on bounds.

        If variable has both upper and lower bounds, scaling factor will be based on the
        mean of the bounds. If variable has only one bound, scaling factor will be based
        on that bound. If variable has no bounds, scaling factor will not be set.

        Args:
            variable - variable to set scaling factor for
            overwrite - whether to overwrite existing scaling factors

        Returns:
            None
        """
        if variable.lb is not None:
            if variable.ub is not None:
                # Both bounds, use mean
                xmag = 0.5 * (variable.ub + variable.lb)
            else:
                # Only lower bound
                xmag = variable.lb
        elif variable.ub is not None:
            # Only upper bound
            xmag = variable.ub
        else:
            # No bounds
            _log.debug(
                f"No scaling factor set for {variable.name}; variable has no bounds."
            )
            return

        if xmag == 0:
            sf = 1
        else:
            sf = 1 / abs(xmag)

        self.set_variable_scaling_factor(
            variable=variable, scaling_factor=sf, overwrite=overwrite
        )

    def scale_variable_by_default(self, variable, overwrite: bool = False):
        """
        Set scaling factor for variable based on default scaling factor.

        Args:
            variable - variable to set scaling factor for
            overwrite - whether to overwrite existing scaling factors

        Returns:
            None
        """
        sf = self.get_default_scaling_factor(variable)
        if sf is not None:
            self.set_variable_scaling_factor(
                variable=variable, scaling_factor=sf, overwrite=overwrite
            )
        else:
            _log.debug(
                f"Could not set scaling factor for {variable.name}, "
                f"no default scaling factor set."
            )

    def scale_variable_by_units(self, variable, overwrite: bool = False):
        """
        Set scaling factor for variable based on units of measurement.

        Units of measurement for variable are compared to those stored in
        self.unit_scaling_factors, and if a mathc is found the scaling factor set
        using the associated value.

        Args:
            variable - variable to set scaling factor for
            overwrite - whether to overwrite existing scaling factors

        Returns:
            None
        """
        uom = variable.units

        if uom in self.unit_scaling_factors:
            self.set_variable_scaling_factor(
                variable=variable,
                scaling_factor=self.unit_scaling_factors[uom],
                overwrite=overwrite,
            )
        else:
            _log.debug(
                f"No scaling factor set for {variable.name}; no match for units {uom} found"
                "in self.unit_scaling_factors"
            )

    # Common methods for constraint scaling
    def scale_constraint_by_component(
        self, target_constraint, scaling_component, overwrite: bool = False
    ):
        """
        Set scaling factor for target_constraint equal to that of scaling_component.

        Args:
            target_constraint - constraint to set scaling factor for
            scaling_component - component to use for scaling factor
            overwrite - whether to overwrite existing scaling factors

        Returns:
            None
        """
        sf = get_scaling_factor(scaling_component)
        if sf is not None:
            self.set_constraint_scaling_factor(
                constraint=target_constraint, scaling_factor=sf, overwrite=overwrite
            )
        else:
            _log.debug(
                f"Could not set scaling factor for {target_constraint.name}, "
                f"no scaling factor set for {scaling_component.name}"
            )

    def scale_constraint_by_default(self, constraint, overwrite: bool = False):
        """
        Set scaling factor for constraint based on default scaling factor.

        Args:
            constraint - constraint to set scaling factor for
            overwrite - whether to overwrite existing scaling factors

        Returns:
            None
        """
        sf = self.get_default_scaling_factor(constraint)
        if sf is not None:
            self.set_constraint_scaling_factor(
                constraint=constraint, scaling_factor=sf, overwrite=overwrite
            )
        else:
            _log.debug(
                f"Could not set scaling factor for {constraint.name}, "
                f"no default scaling factor set."
            )

    # by magnitude of terms - harmonic mean, max, min, 1-norm, 2-norm, ...
