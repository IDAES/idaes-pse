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
from enum import Enum

from pyomo.environ import ComponentMap, units, value
from pyomo.core.base.units_container import UnitsError
from pyomo.core.expr import identify_variables
from pyomo.core.expr.calculus.derivatives import Modes, differentiate

from idaes.core.scaling.scaling_base import CONFIG, ScalerBase
from idaes.core.scaling.util import get_scaling_factor, NominalValueExtractionVisitor
import idaes.logger as idaeslog

# Set up logger
_log = idaeslog.getLogger(__name__)

CSCONFIG = CONFIG()

DEFAULT_UNIT_SCALING = {
    # "QuantityName: (reference units, scaling factor)
    # Model developers should be careful when using these, especially when
    # dealing with differential measurements (e.g. pressure and temperature differences)
    "Temperature": (units.K, 1e-2),
    "Pressure": (units.Pa, 1e-5),
}


class ConstraintScalingScheme(str, Enum):
    """
    Schemes available for calculating constraint scaling factors.

    * harmonicMean ('harmonic_mean'): sf = sum(1/abs(nominal value))
    * inverseSum ('inverse_sum'): sf = 1 / sum(abs(nominal value))
    * inverseRSS ('inverse_root_sum_squared'): sf = 1 / sqrt(sum(abs(nominal value)**2))
    * inverseMaximum ('inverse_maximum'): sf =  1 / max(abs(nominal value)
    * inverseMinimum ('inverse_minimum'): sf = 1 / min(abs(nominal value)
    """

    harmonicMean = "harmonic_mean"
    inverseSum = "inverse_sum"
    inverseRSS = "inverse_root_sum_squared"
    inverseMaximum = "inverse_maximum"
    inverseMinimum = "inverse_minimum"


class CustomScalerBase(ScalerBase):
    """
    Base class for custom scaling routines.

    """

    CONFIG = ScalerBase.CONFIG()

    # Common data structures for default scaling
    # DEFAULT_SCALING_FACTORS = {"component_local_name": DEFAULT_SCALING}
    DEFAULT_SCALING_FACTORS = None

    # UNIT_SCALING_FACTORS = {"units": UNIT_BASED_SCALING}
    UNIT_SCALING_FACTORS = copy(DEFAULT_UNIT_SCALING)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

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
        submodel_scalers: ComponentMap = None,
    ):
        """
        Default model scaling routine.

        This method performs a four-step scaling routine:

        1. Scale variables using variable_scaling_routine
        2. Perform first-stage scaling factor fill in using user provided method(s), called in order declared
        3. Scale constraints using constraint_scaling_routine
        4. Perform second-stage scaling factor fill in using user provided method(s), called in order declared

        Args:
            model: model to be scaled
            first_stage_fill_in: list of methods to use for first-stage scaling factor fill in
            second_stage_fill_in: list of methods to use for second-stage scaling factor fill in
            submodel_scalers: ComponentMap of Scalers to use for sub-models

        Returns:
            None
        """
        # Step 1: Call variable scaling routine
        self.variable_scaling_routine(
            model, overwrite=self.config.overwrite, submodel_scalers=submodel_scalers
        )
        # Step 2: Call variable fill in
        if first_stage_fill_in is not None:
            for i in first_stage_fill_in:
                i(model)

        # Step 3: Call constraint scaling routine
        self.constraint_scaling_routine(
            model, overwrite=self.config.overwrite, submodel_scalers=submodel_scalers
        )

        # Step 4: Call constraint fill in
        if second_stage_fill_in is not None:
            for i in second_stage_fill_in:
                i(model)

    def variable_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: ComponentMap = None
    ):
        """
        Routine to apply scaling factors to variables in model.

        Derived classes must overload this method.

        Args:
            model: model to be scaled
            overwrite: whether to overwrite existing scaling factors
            submodel_scalers: ComponentMap of Scalers to use for sub-models

        Returns:
            None
        """
        raise NotImplementedError(
            "Custom Scaler has not implemented a variable_scaling_routine method."
        )

    def constraint_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: ComponentMap = None
    ):
        """
        Routine to apply scaling factors to constraints in model.

        Derived classes must overload this method.

        Args:
            model: model to be scaled
            overwrite: whether to overwrite existing scaling factors
            submodel_scalers: ComponentMap of Scalers to use for sub-models

        Returns:
            None
        """
        raise NotImplementedError(
            "Custom Scaler has not implemented a constraint_scaling_routine method."
        )

    def get_default_scaling_factor(self, component):
        """
        Get scaling factor for component from dict of default values.

        Args:
            component: component to get default scaling factor for

        Returns:
            default scaling factor if it exists, else None
        """
        try:
            return self.default_scaling_factors[component.local_name]
        except KeyError:
            # Might be indexed, see if parent component has default scaling
            parent = component.parent_component()
            try:
                return self.default_scaling_factors[parent.local_name]
            except KeyError:
                # No default scaling found, give up
                pass

            # Log a message and return nothing
            _log.debug(f"No default scaling factor found for {component.name}")
            return None

    # Common methods for variable scaling
    def scale_variable_by_component(
        self, target_variable, scaling_component, overwrite: bool = False
    ):
        """
        Set scaling factor for target_variable equal to that of scaling_component.

        Args:
            target_variable: variable to set scaling factor for
            scaling_component: component to use for scaling factor
            overwrite: whether to overwrite existing scaling factors

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
            variable: variable to set scaling factor for
            overwrite: whether to overwrite existing scaling factors

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
            variable: variable to set scaling factor for
            overwrite: whether to overwrite existing scaling factors

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
        self.unit_scaling_factors, and if a match is found the scaling factor set
        using the associated value.

        Args:
            variable: variable to set scaling factor for
            overwrite: whether to overwrite existing scaling factors

        Returns:
            None
        """
        uom = units.get_units(variable)

        sf = None
        # Keys in self.unit_scaling_factors are not used - only required because Pyomo
        # units are non-hashable and thus cannot be keys
        for refunits, unit_scale in self.unit_scaling_factors.values():
            try:
                # Try convert the reference scaling factor to variable units
                # TODO: Have not found a more efficient way to do this, as Pyomo basically
                # TODO: involves a try/except anyway
                # Need to invert reference scaling factor, and then invert result
                sf = 1 / units.convert_value(
                    1 / unit_scale, from_units=refunits, to_units=uom
                )
                # Break once we have a match - no need to continue
                break
            except UnitsError:
                pass

        if sf is not None:
            self.set_variable_scaling_factor(
                variable=variable,
                scaling_factor=sf,
                overwrite=overwrite,
            )
        else:
            _log.debug(
                f"No scaling factor set for {variable.name}; no match for units {uom} found "
                "in self.unit_scaling_factors"
            )

    # Common methods for constraint scaling
    def scale_constraint_by_component(
        self, target_constraint, scaling_component, overwrite: bool = False
    ):
        """
        Set scaling factor for target_constraint equal to that of scaling_component.

        Args:
            target_constraint: constraint to set scaling factor for
            scaling_component: component to use for scaling factor
            overwrite: whether to overwrite existing scaling factors

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
            constraint: constraint to set scaling factor for
            overwrite: whether to overwrite existing scaling factors

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

    def get_expression_nominal_values(self, expression):
        """
        Calculate nominal values for each additive term in a Pyomo expression.

        The nominal value of any Var is defined as the inverse of its scaling factor
        (if assigned, else 1).

        Args:
            expression: Pyomo expression to collect nominal values for

        Returns:
            list of nominal values for each additive term
        """
        # For convenience, if expression is a Pyomo component with an expr attribute,
        # redirect to the expr attribute
        if hasattr(expression, "expr"):
            expression = expression.expr

        return NominalValueExtractionVisitor().walk_expression(expression)

    def scale_constraint_by_nominal_value(
        self,
        constraint,
        scheme: ConstraintScalingScheme = ConstraintScalingScheme.inverseMaximum,
        overwrite: bool = False,
    ):
        """
        Set scaling factor for constraint based on the nominal value(s).

        Terms with expected magnitudes of 0 will be ignored.

        Args:
            constraint: constraint to set scaling factor for
            scheme: ConstraintScalingScheme Enum indicating method to apply
              for determining constraint scaling.
            overwrite: whether to overwrite existing scaling factors

        Returns:
            None
        """
        nominal = self.get_expression_nominal_values(constraint.expr)

        # Remove any 0 terms
        nominal = [j for j in nominal if j != 0]

        if len(nominal) == 0:
            # No non-zero terms...
            sf = 1
        elif scheme == ConstraintScalingScheme.harmonicMean:
            sf = sum(1 / abs(i) for i in nominal)
        elif scheme == ConstraintScalingScheme.inverseSum:
            sf = 1 / sum(abs(i) for i in nominal)
        elif scheme == ConstraintScalingScheme.inverseRSS:
            sf = 1 / sum(abs(i) ** 2 for i in nominal) ** 0.5
        elif scheme == ConstraintScalingScheme.inverseMaximum:
            sf = 1 / max(abs(i) for i in nominal)
        elif scheme == ConstraintScalingScheme.inverseMinimum:
            sf = 1 / min(abs(i) for i in nominal)
        else:
            raise ValueError(
                f"Invalid value for 'scheme' argument ({scheme}) in "
                "scale_constraint_by_nominal_value."
            )

        self.set_constraint_scaling_factor(
            constraint=constraint, scaling_factor=sf, overwrite=overwrite
        )

    def scale_constraint_by_nominal_derivative_norm(
        self, constraint, norm: int = 2, overwrite: bool = False
    ):
        """
        Scale constraint by norm of partial derivatives.

        Calculates partial derivatives of constraint at nominal variable values,
        and then scaled the constraint by the user-selected norm of these derivatives.
        Given perfect variable scaling, this should provide a similar result to
        applying scaling based on the Jacobian norm, however this approach does not
        require an initial solution for the problem (relying on nominal values instead).

        Args:
            constraint: constraint to be scaled.
            norm: type of norm to use for scaling. Must be a positive integer.
            overwrite: whether to overwrite existing scaling factors.

        Returns:
            None
        """
        # Cast norm to int to make sure it is valid
        norm = int(norm)
        if norm < 1:
            raise ValueError(f"norm must be a positive integer (received {norm})")

        var_data = []
        try:
            # Iterate over all variables in constraint
            for v in identify_variables(constraint.body):
                # Store current value for restoration
                ov = v.value  # original value
                sf = self.get_scaling_factor(v)  # scaling factor
                if sf is None:
                    # If no scaling factor set, use nominal value of 1
                    sf = 1

                var_data.append((v, ov, sf))

            # Get partial derivatives
            pjac = []
            for v in var_data:
                # Iterate over all variable and set values
                for w in var_data:
                    if w is not v:
                        # Set all other variables to their nominal magnitude
                        # nominal_value = 1/ scaling_factor
                        w[0].value = 1 / w[2]
                    else:
                        # Set derivative var to scaled value of 1
                        # With perfect scaling, scaling_factor * value = 1
                        w[0].value = 1

                pjac.append(
                    value(
                        differentiate(
                            expr=constraint.body, wrt=v[0], mode=Modes.reverse_symbolic
                        )
                        * (1 / v[2])  # Need to divide by scaling_factor
                    )
                )

        finally:
            # Restore all values for clean up
            for v in var_data:
                v[0].value = v[1]

        # Calculate norm
        cnorm = sum(abs(j) ** norm for j in pjac) ** (1 / norm)
        if cnorm != 0:
            sf = 1 / cnorm
        else:
            sf = 1
        self.set_constraint_scaling_factor(constraint, sf, overwrite=overwrite)

    # Other methods
    def propagate_state_scaling(
        self, target_state, source_state, overwrite: bool = False
    ):
        """
        Propagate scaling of state variables from one StateBlock to another.

        Indexing of target and source StateBlocks must match.

        Args:
            target_state: StateBlock to set scaling factors on
            source_state: StateBlock to use as source for scaling factors
            overwrite: whether to overwrite existing scaling factors

        Returns:
            None
        """
        for bidx, target_data in target_state.items():
            target_vars = target_data.define_state_vars()
            source_vars = source_state[bidx].define_state_vars()

            for state, var in target_vars.items():
                for vidx, vardata in var.items():
                    self.scale_variable_by_component(
                        target_variable=vardata,
                        scaling_component=source_vars[state][vidx],
                        overwrite=overwrite,
                    )

    def call_submodel_scaler_method(
        self,
        model,
        submodel: str,
        method: str,
        submodel_scalers: ComponentMap = None,
        overwrite: bool = False,
    ):
        """
        Call scaling method for submodel.

        Scaler for submodel is taken from submodel_scalers if present, otherwise the
        default scaler for the submodel is used.

        Args:
            model: parent model of submodel
            submodel: local name of submodel to be scaled as str
            submodel_scalers: user provided ComponentMap of Scalers to use for submodels
            method: name of method to call from submodel (as string)

        Returns:
            None
        """
        # Get actual submodel object from name
        # For this method, we have to use the component name as the Scaler is written
        # before the model is constructed.
        sm_obj = model.find_component(submodel)

        if submodel_scalers is None:
            submodel_scalers = {}

        # Iterate over indices of submodel
        for smdata in sm_obj.values():
            # Get Scaler for submodel
            if sm_obj in submodel_scalers:
                scaler = submodel_scalers[sm_obj]
                if callable(scaler):
                    # Check to see if Scaler is callable - this implies it is a class and not an instance
                    # Call the class to create an instance
                    scaler = scaler()
                _log.debug(f"Using user-defined Scaler for {model}.{submodel}.")
            else:
                try:
                    scaler = smdata.default_scaler
                    _log.debug(f"Using default Scaler for {model}.{submodel}.")
                except AttributeError:
                    _log.debug(
                        f"No default Scaler set for {model}.{submodel}. Cannot call {method}."
                    )
                    return
                if scaler is not None:
                    scaler = scaler()
                else:
                    _log.debug(
                        f"No Scaler found for {model}.{submodel}. Cannot call {method}."
                    )

            # If a Scaler is found, call desired method
            if scaler is not None:
                try:
                    smeth = getattr(scaler, method)
                except AttributeError:
                    raise AttributeError(
                        f"Scaler for {model}.{submodel} does not have a method named {method}."
                    )
                smeth(smdata, overwrite=overwrite)
