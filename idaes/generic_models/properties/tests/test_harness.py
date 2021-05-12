##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
import pytest

from pyomo.environ import (ConcreteModel,
                           Expression,
                           Set,
                           Var,
                           Constraint)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import (ControlVolume0DBlock,
                        FlowsheetBlock,
                        PhysicalParameterBlock,
                        StateBlockData,
                        MaterialBalanceType,
                        EnergyBalanceType,
                        MaterialFlowBasis)
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util import get_solver
import idaes.core.util.scaling as iscale


from pyomo.core.base.var import _VarData
from pyomo.core.base.param import _ParamData
from pyomo.core.base.expression import _ExpressionData

_scalable = (_VarData, _ParamData, _ExpressionData)

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


# -----------------------------------------------------------------------------
class PropertyTestHarness(object):
    def configure_class(self, m):
        self.prop_pack = None
        self.param_args = {}
        self.prop_args = {}

        self.has_density_terms = True

        self.configure()

        # Need to attach these to m for later use
        m.has_density_terms = self.has_density_terms
        m.prop_args = self.prop_args

    def configure(self):
        # Placeholder method to allow user to setup test harness
        pass

    @pytest.fixture(scope="class")
    def frame(self):
        m = ConcreteModel()
        self.configure_class(m)

        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.params = self.prop_pack(default=self.param_args)

        return m

    def test_param_block(self, frame):
        if not isinstance(frame.fs.params, PhysicalParameterBlock):
            raise TypeError(
                "Parameter Block does not inherit from PhysicalParameterBlock")

    def test_component_list(self, frame):
        if not isinstance(frame.fs.params.component_list, Set):
            raise TypeError(
                "Parameter block does not have a component_list, or "
                "component_list is not a Pyomo Set.")

    def test_phase_list(self, frame):
        if not isinstance(frame.fs.params.phase_list, Set):
            raise TypeError(
                "Parameter block does not have a phase_list, or "
                "phase_list is not a Pyomo Set.")

    def test_state_block_class(self, frame):
        if not hasattr(frame.fs.params, "state_block_class"):
            raise AttributeError(
                "Parameter block does not specify state_block_class.")

    def test_properties_meta_data(self, frame):
        if frame.fs.params.get_metadata().properties is None:
            raise AttributeError(
                "Parameter block has not specified properties metadata.")

    def test_default_units_meta_data(self, frame):
        if frame.fs.params.get_metadata().default_units is None:
            raise AttributeError(
                "Parameter block has not specified default_units metadata.")

    def test_state_block_construction(self, frame):
        frame.fs.props = frame.fs.params.build_state_block(
                [1], default={"defined_state": True,
                              **frame.prop_args})

        if not isinstance(frame.fs.props[1], StateBlockData):
            raise TypeError(
                "State Block does not inherit from StateBlockData.")

    def test_state_block_temperature(self, frame):
        if not isinstance(frame.fs.props[1].temperature, (Var, Expression)):
            raise AttributeError(
                "State Block does not contain temperature, or temperature is "
                "not a Pyomo Var or Expression.")

    def test_state_block_pressure(self, frame):
        if not isinstance(frame.fs.props[1].pressure, (Var, Expression)):
            raise AttributeError(
                "State Block does not contain pressure, or pressure is "
                "not a Pyomo Var or Expression.")

    def test_state_block_mole_frac_phase_comp(self, frame):
        if not isinstance(frame.fs.props[1].mole_frac_phase_comp,
                          (Var, Expression)):
            raise AttributeError(
                "State Block does not contain mole_frac_phase_comp, or "
                "mole_frac_phase_comp is not a Pyomo Var or Expression.")
        err = False
        if len(frame.fs.props[1].mole_frac_phase_comp) != (
                len(frame.fs.params.component_list) *
                len(frame.fs.params.phase_list)):
            err = True

        for k in frame.fs.props[1].mole_frac_phase_comp:
            if (k[0] not in frame.fs.params.phase_list or
                    k[1] not in frame.fs.params.component_list):
                err = True

        if err:
            raise ValueError(
                "mole_frac_phase_comp is not indexed by phase and component.")

    def test_get_material_flow_terms(self, frame):
        try:
            for p in frame.fs.params.phase_list:
                for j in frame.fs.params.component_list:
                    term = frame.fs.props[1].get_material_flow_terms(p, j)
                    # Assert that the term can be assigned a scale factor
                    assert isinstance(term, _scalable)
        except KeyError:
            raise KeyError(
                "get_material_flow_terms method is not indexed by phase and "
                "component.")
        except AttributeError:
            raise AttributeError(
                "State block has not implemented get_material_flow_terms "
                "method.")

    def test_get_enthalpy_flow_terms(self, frame):
        try:
            for p in frame.fs.params.phase_list:
                term = frame.fs.props[1].get_enthalpy_flow_terms(p)
                # Assert that the term can be assigned a scale factor
                assert isinstance(term, _scalable)
        except KeyError:
            raise KeyError(
                "get_enthalpy_flow_terms method is not indexed by phase.")
        except AttributeError:
            raise AttributeError(
                "State block has not implemented get_enthalpy_flow_terms "
                "method.")

    def test_get_material_density_terms(self, frame):
        if frame.has_density_terms:
            try:
                for p in frame.fs.params.phase_list:
                    for j in frame.fs.params.component_list:
                        term = frame.fs.props[1].get_material_density_terms(p, j)
                        # Assert that the term can be assigned a scale factor
                        assert isinstance(term, _scalable)
            except KeyError:
                raise KeyError(
                    "get_material_density_terms method is not indexed by phase"
                    " and component.")
            except AttributeError:
                raise AttributeError(
                    "State block has not implemented "
                    "get_material_density_terms method.")

    def test_get_energy_density_terms(self, frame):
        if frame.has_density_terms:
            try:
                for p in frame.fs.params.phase_list:
                    term = frame.fs.props[1].get_energy_density_terms(p)
                    assert isinstance(term, _scalable)
            except KeyError:
                raise KeyError(
                    "get_enthalpy_density_terms method is not indexed by "
                    "phase.")
            except AttributeError:
                raise AttributeError(
                    "State block has not implemented "
                    "get_enthalpy_density_terms method.")

    def test_default_material_balance_type(self, frame):
        if frame.fs.props[1].default_material_balance_type() not in \
                MaterialBalanceType:
            raise ValueError(
                "Invalid value for default_material_balance_type.")

    def test_default_energy_balance_type(self, frame):
        if frame.fs.props[1].default_energy_balance_type() not in \
                EnergyBalanceType:
            raise ValueError(
                "Invalid value for default_energy_balance_type.")

    def test_get_material_flow_basis(self, frame):
        if frame.fs.props[1].get_material_flow_basis() not in \
                MaterialFlowBasis:
            raise ValueError(
                "Invalid value for get_material_flow_basis.")

    def test_define_state_vars(self, frame):
        try:
            sv = frame.fs.props[1].define_state_vars()
        except AttributeError:
            raise AttributeError(
                "State block has not implemented define_state_vars method.")

        if not isinstance(sv, dict):
            raise TypeError(
                "define_state_vars must return a dictionary of Vars.")

        for v in sv.values():
            if not isinstance(v, Var):
                raise TypeError(
                    "Invlaid entry in define_state_Vars, {}. All members must "
                    "be Pyomo Vars.".format(v))

    def test_unit_consistency(self, frame):
        assert_units_consistent(frame)

    def test_initialize(self, frame):
        frame._init_dof = degrees_of_freedom(frame.fs.props[1])

        if not hasattr(frame.fs.props, "initialize") or \
                not callable(frame.fs.props.initialize):
            raise AttributeError(
                "State Block has not implemented an initialize method.")

        frame._flags = frame.fs.props.initialize(hold_state=True)

        if degrees_of_freedom(frame.fs.props[1]) != 0:
            raise Exception(
                "initialize did not result in a State BLock with 0 "
                "degrees of freedom.")

    def test_release_state(self, frame):
        if not hasattr(frame.fs.props, "release_state") or \
                not callable(frame.fs.props.release_state):
            raise AttributeError(
                "State Block has not implemented a release_state method.")

        frame.fs.props.release_state(frame._flags)

        if degrees_of_freedom(frame.fs.props[1]) != frame._init_dof:
            raise Exception(
                "release state did not restore State Block to original "
                "degrees of freedom.")

    def test_CV_integration(self, frame):
        frame.fs.cv = ControlVolume0DBlock(default={
                "property_package": frame.fs.params})

        frame.fs.cv.add_geometry()

        frame.fs.cv.add_state_blocks(has_phase_equilibrium=True)

        frame.fs.cv.add_material_balances(has_phase_equilibrium=True)

        frame.fs.cv.add_energy_balances()

        frame.fs.cv.add_momentum_balances()

    def test_default_scaling_factors(self, frame):
        # check that the calculate_scaling_factors method successfully copies
        # the default scaling factors to the scaling suffixes.  If there are
        # no default scaling factors, this should pass
        iscale.calculate_scaling_factors(frame) #this also ensure, it doesn't except 
        for v in frame.fs.props[1].component_data_objects(
            (Constraint, Var, Expression),
            descend_into=False):
            name = v.getname().split("[")[0]
            index = v.index()
            print(v)
            assert (iscale.get_scaling_factor(v) ==
                frame.fs.props[1].config.parameters.get_default_scaling(name, index)) or (
                frame.fs.props[1].config.parameters.get_default_scaling(name, index) is None)
