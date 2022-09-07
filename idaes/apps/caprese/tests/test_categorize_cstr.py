#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
"""
Test for Caprese's module for NMPC.
"""

import pytest
from pyomo.environ import Var
from pyomo.common.collections import ComponentSet
from pyomo.dae.flatten import flatten_dae_components

from idaes.apps.caprese.categorize import (
    categorize_dae_variables,
)
from idaes.apps.caprese.common.config import VariableCategory as VC
from idaes.apps.caprese.examples.cstr_model import make_model

__author__ = "Robert Parker"


@pytest.mark.unit
def test_categorize_1():
    """
    This test categorizes the enzyme cstr "as intended."
    Volume is used as a measurement, and solvent flow rates are
    fixed, but otherwise equations are as expected.
    """
    model = make_model(horizon=1, ntfe=5, ntcp=2)
    time = model.fs.time

    init_input_list = [
        model.fs.mixer.S_inlet.flow_vol[0],
        model.fs.mixer.E_inlet.flow_vol[0],
    ]
    init_input_set = ComponentSet(init_input_list)

    init_deriv_list = [
        model.fs.cstr.control_volume.energy_accumulation[0, "aq"],
        *list(model.fs.cstr.control_volume.material_accumulation[0, "aq", :]),
    ]
    init_deriv_set = ComponentSet(init_deriv_list)

    init_diff_list = [
        model.fs.cstr.control_volume.energy_holdup[0, "aq"],
        *list(model.fs.cstr.control_volume.material_holdup[0, "aq", :]),
    ]
    init_diff_set = ComponentSet(init_diff_list)

    init_fixed_list = [
        model.fs.mixer.E_inlet.temperature[0],
        model.fs.mixer.S_inlet.temperature[0],
        *list(model.fs.mixer.E_inlet.conc_mol[0, :]),
        *list(model.fs.mixer.S_inlet.conc_mol[0, :]),
        model.fs.cstr.outlet.conc_mol[0, "Solvent"],
        model.fs.mixer.outlet.conc_mol[0, "Solvent"],
    ]
    init_fixed_set = ComponentSet(init_fixed_list)

    init_meas_list = [
        model.fs.cstr.control_volume.energy_holdup[0, "aq"],
        model.fs.cstr.control_volume.volume[0],
        *list(model.fs.cstr.control_volume.material_holdup[0, "aq", :]),
    ]
    init_meas_set = ComponentSet(init_meas_list)
    # Solvent holdup is not a measurement; we measure volume instead
    init_meas_set.remove(
        model.fs.cstr.control_volume.material_holdup[0, "aq", "Solvent"]
    )

    init_alg_list = [
        model.fs.cstr.outlet.flow_vol[0],
        model.fs.cstr.outlet.temperature[0],
        model.fs.cstr.inlet.flow_vol[0],
        model.fs.cstr.inlet.temperature[0],
        model.fs.mixer.outlet.flow_vol[0],
        model.fs.mixer.outlet.temperature[0],
        model.fs.cstr.control_volume.volume[0],
        *list(model.fs.cstr.control_volume.properties_out[0].flow_mol_comp[:]),
        *list(model.fs.cstr.inlet.conc_mol[0, :]),
        *list(model.fs.cstr.control_volume.properties_in[0].flow_mol_comp[:]),
        *list(model.fs.cstr.control_volume.rate_reaction_generation[0, "aq", :]),
        *list(model.fs.mixer.mixed_state[0].flow_mol_comp[:]),
        *list(model.fs.mixer.E_inlet_state[0].flow_mol_comp[:]),
        *list(model.fs.mixer.S_inlet_state[0].flow_mol_comp[:]),
        *list(model.fs.cstr.outlet.conc_mol[0, :]),
        *list(model.fs.mixer.outlet.conc_mol[0, :]),
        *list(model.fs.cstr.control_volume.reactions[0].reaction_coef[:]),
        *list(model.fs.cstr.control_volume.reactions[0].reaction_rate[:]),
        *list(model.fs.cstr.control_volume.rate_reaction_extent[0, :]),
    ]
    init_alg_set = ComponentSet(init_alg_list)
    # solvent outlet concentrations are not algebraic variables as
    # it is fixed.
    init_alg_set.remove(model.fs.cstr.outlet.conc_mol[0, "Solvent"])
    init_alg_set.remove(model.fs.mixer.outlet.conc_mol[0, "Solvent"])

    scalar_vars, dae_vars = flatten_dae_components(model, time, ctype=Var)
    category_dict = categorize_dae_variables(dae_vars, time, init_input_list)

    input_vars = category_dict[VC.INPUT]
    diff_vars = category_dict[VC.DIFFERENTIAL]
    deriv_vars = category_dict[VC.DERIVATIVE]
    fixed_vars = category_dict[VC.FIXED]
    alg_vars = category_dict[VC.ALGEBRAIC]
    meas_vars = category_dict[VC.MEASUREMENT]

    assert len(input_vars) == len(init_input_set)
    for v in input_vars:
        assert v[0] in init_input_set

    assert len(deriv_vars) == len(init_deriv_set)
    for v in deriv_vars:
        assert v[0] in init_deriv_set

    assert len(diff_vars) == len(init_deriv_set)
    for v in diff_vars:
        assert v[0] in init_diff_set

    assert len(fixed_vars) == len(init_fixed_set)
    for v in fixed_vars:
        assert v[0] in init_fixed_set

    assert len(alg_vars) == len(init_alg_set)
    for v in alg_vars:
        assert v[0] in init_alg_set

    assert len(meas_vars) == len(init_meas_set)
    for v in meas_vars:
        assert v[0] in init_meas_set

    assert len(scalar_vars) == 0


@pytest.mark.unit
def test_categorize_2():
    """
    In this instance, temperature is "measured" (used as an initial condition)
    instead of energy_holdup, conc[P] is measured instead of holdup[P],
    and accumulation[C] is measured instead of holdup[C].
    """
    model = make_model(horizon=1, ntfe=5, ntcp=2)
    time = model.fs.time
    t0 = time.first()

    material_holdup = model.fs.cstr.control_volume.material_holdup
    material_accumulation = model.fs.cstr.control_volume.material_accumulation
    energy_holdup = model.fs.cstr.control_volume.energy_holdup
    energy_accumulation = model.fs.cstr.control_volume.energy_accumulation
    conc_mol = model.fs.cstr.outlet.conc_mol

    # Specify temperature instead of energy holdup
    energy_holdup[t0, "aq"].unfix()
    model.fs.cstr.outlet.temperature[t0].fix(300)

    # Specify C_P instead of holdup
    material_holdup[t0, "aq", "P"].unfix()
    conc_mol[t0, "P"].fix(0.0)

    # Specify accumulation of C instead of holdup
    # You might want to do this if you want to start at steady state,
    # but don't know the value of every variable at the steady state
    # you want to start at...
    material_holdup[t0, "aq", "C"].unfix()
    material_accumulation[t0, "aq", "C"].fix(0.0)

    init_input_list = [
        model.fs.mixer.S_inlet.flow_vol[0],
        model.fs.mixer.E_inlet.flow_vol[0],
    ]
    init_input_set = ComponentSet(init_input_list)

    init_deriv_list = [
        model.fs.cstr.control_volume.energy_accumulation[0, "aq"],
        *list(model.fs.cstr.control_volume.material_accumulation[0, "aq", :]),
    ]
    init_deriv_set = ComponentSet(init_deriv_list)

    init_diff_list = [
        model.fs.cstr.control_volume.energy_holdup[0, "aq"],
        *list(model.fs.cstr.control_volume.material_holdup[0, "aq", :]),
    ]
    init_diff_set = ComponentSet(init_diff_list)

    init_fixed_list = [
        model.fs.mixer.E_inlet.temperature[0],
        model.fs.mixer.S_inlet.temperature[0],
        *list(model.fs.mixer.E_inlet.conc_mol[0, :]),
        *list(model.fs.mixer.S_inlet.conc_mol[0, :]),
        model.fs.cstr.outlet.conc_mol[0, "Solvent"],
        model.fs.mixer.outlet.conc_mol[0, "Solvent"],
    ]
    init_fixed_set = ComponentSet(init_fixed_list)

    init_meas_list = [
        model.fs.cstr.outlet.temperature[0],
        model.fs.cstr.control_volume.volume[0],
        model.fs.cstr.control_volume.material_holdup[0, "aq", "S"],
        model.fs.cstr.control_volume.material_holdup[0, "aq", "E"],
        model.fs.cstr.control_volume.material_holdup[0, "aq", "S"],
        model.fs.cstr.control_volume.material_accumulation[0, "aq", "C"],
        model.fs.cstr.outlet.conc_mol[0, "P"],
    ]
    init_meas_set = ComponentSet(init_meas_list)
    # No need to remove solvent holdup here as it was not added to this list.

    init_alg_list = [
        model.fs.cstr.outlet.flow_vol[0],
        model.fs.cstr.outlet.temperature[0],
        model.fs.cstr.inlet.flow_vol[0],
        model.fs.cstr.inlet.temperature[0],
        model.fs.mixer.outlet.flow_vol[0],
        model.fs.mixer.outlet.temperature[0],
        model.fs.cstr.control_volume.volume[0],
        *list(model.fs.cstr.control_volume.properties_out[0].flow_mol_comp[:]),
        *list(model.fs.cstr.inlet.conc_mol[0, :]),
        *list(model.fs.cstr.control_volume.properties_in[0].flow_mol_comp[:]),
        *list(model.fs.cstr.control_volume.rate_reaction_generation[0, "aq", :]),
        *list(model.fs.mixer.mixed_state[0].flow_mol_comp[:]),
        *list(model.fs.mixer.E_inlet_state[0].flow_mol_comp[:]),
        *list(model.fs.mixer.S_inlet_state[0].flow_mol_comp[:]),
        *list(model.fs.cstr.outlet.conc_mol[0, :]),
        *list(model.fs.mixer.outlet.conc_mol[0, :]),
        *list(model.fs.cstr.control_volume.reactions[0].reaction_coef[:]),
        *list(model.fs.cstr.control_volume.reactions[0].reaction_rate[:]),
        *list(model.fs.cstr.control_volume.rate_reaction_extent[0, :]),
    ]
    init_alg_set = ComponentSet(init_alg_list)
    # solvent outlet concentrations are not algebraic variables as
    # it is fixed.
    init_alg_set.remove(model.fs.cstr.outlet.conc_mol[0, "Solvent"])
    init_alg_set.remove(model.fs.mixer.outlet.conc_mol[0, "Solvent"])

    scalar_vars, dae_vars = flatten_dae_components(model, time, ctype=Var)
    category_dict = categorize_dae_variables(dae_vars, time, init_input_list)

    input_vars = category_dict[VC.INPUT]
    diff_vars = category_dict[VC.DIFFERENTIAL]
    deriv_vars = category_dict[VC.DERIVATIVE]
    fixed_vars = category_dict[VC.FIXED]
    alg_vars = category_dict[VC.ALGEBRAIC]
    meas_vars = category_dict[VC.MEASUREMENT]

    assert len(input_vars) == len(init_input_set)
    for v in input_vars:
        assert v[0] in init_input_set

    assert len(deriv_vars) == len(init_deriv_set)
    for v in deriv_vars:
        assert v[0] in init_deriv_set

    assert len(diff_vars) == len(init_deriv_set)
    for v in diff_vars:
        assert v[0] in init_diff_set

    assert len(fixed_vars) == len(init_fixed_set)
    for v in fixed_vars:
        assert v[0] in init_fixed_set

    assert len(alg_vars) == len(init_alg_set)
    for v in alg_vars:
        assert v[0] in init_alg_set

    assert len(meas_vars) == len(init_meas_set)
    for v in meas_vars:
        assert v[0] in init_meas_set

    assert len(scalar_vars) == 0


@pytest.mark.unit
def test_categorize_3():
    """
    In this test we fix one of the differential variables.
    This is the case where somebody wants to run an isothermal
    CSTR.
    """
    model = make_model(horizon=1, ntfe=5, ntcp=2)
    time = model.fs.time

    CV = model.fs.cstr.control_volume
    CV.energy_holdup.fix(300)
    CV.energy_accumulation_disc_eq.deactivate()

    init_input_list = [
        model.fs.mixer.S_inlet.flow_vol[0],
        model.fs.mixer.E_inlet.flow_vol[0],
    ]
    init_input_set = ComponentSet(init_input_list)

    init_deriv_list = [
        *list(model.fs.cstr.control_volume.material_accumulation[0, "aq", :]),
    ]
    init_deriv_set = ComponentSet(init_deriv_list)

    init_diff_list = [
        *list(model.fs.cstr.control_volume.material_holdup[0, "aq", :]),
    ]
    init_diff_set = ComponentSet(init_diff_list)

    init_fixed_list = [
        # Energy holdup has been fixed
        model.fs.cstr.control_volume.energy_holdup[0, "aq"],
        model.fs.mixer.E_inlet.temperature[0],
        model.fs.mixer.S_inlet.temperature[0],
        *list(model.fs.mixer.E_inlet.conc_mol[0, :]),
        *list(model.fs.mixer.S_inlet.conc_mol[0, :]),
        model.fs.cstr.outlet.conc_mol[0, "Solvent"],
        model.fs.mixer.outlet.conc_mol[0, "Solvent"],
    ]
    init_fixed_set = ComponentSet(init_fixed_list)

    init_meas_list = [
        model.fs.cstr.control_volume.volume[0],
        *list(model.fs.cstr.control_volume.material_holdup[0, "aq", :]),
    ]
    init_meas_set = ComponentSet(init_meas_list)
    # Solvent holdup is not a measurement; we measure volume instead
    init_meas_set.remove(
        model.fs.cstr.control_volume.material_holdup[0, "aq", "Solvent"]
    )

    init_alg_list = [
        # Since energy_holdup is fixed, energy_accumulation is not
        # considered a derivative. Instead it is algebraic.
        model.fs.cstr.control_volume.energy_accumulation[0, "aq"],
        model.fs.cstr.outlet.flow_vol[0],
        model.fs.cstr.outlet.temperature[0],
        model.fs.cstr.inlet.flow_vol[0],
        model.fs.cstr.inlet.temperature[0],
        model.fs.mixer.outlet.flow_vol[0],
        model.fs.mixer.outlet.temperature[0],
        model.fs.cstr.control_volume.volume[0],
        *list(model.fs.cstr.control_volume.properties_out[0].flow_mol_comp[:]),
        *list(model.fs.cstr.inlet.conc_mol[0, :]),
        *list(model.fs.cstr.control_volume.properties_in[0].flow_mol_comp[:]),
        *list(model.fs.cstr.control_volume.rate_reaction_generation[0, "aq", :]),
        *list(model.fs.mixer.mixed_state[0].flow_mol_comp[:]),
        *list(model.fs.mixer.E_inlet_state[0].flow_mol_comp[:]),
        *list(model.fs.mixer.S_inlet_state[0].flow_mol_comp[:]),
        *list(model.fs.cstr.outlet.conc_mol[0, :]),
        *list(model.fs.mixer.outlet.conc_mol[0, :]),
        *list(model.fs.cstr.control_volume.reactions[0].reaction_coef[:]),
        *list(model.fs.cstr.control_volume.reactions[0].reaction_rate[:]),
        *list(model.fs.cstr.control_volume.rate_reaction_extent[0, :]),
    ]
    init_alg_set = ComponentSet(init_alg_list)
    # solvent outlet concentrations are not algebraic variables as
    # it is fixed.
    init_alg_set.remove(model.fs.cstr.outlet.conc_mol[0, "Solvent"])
    init_alg_set.remove(model.fs.mixer.outlet.conc_mol[0, "Solvent"])

    scalar_vars, dae_vars = flatten_dae_components(model, time, ctype=Var)
    category_dict = categorize_dae_variables(dae_vars, time, init_input_list)

    input_vars = category_dict[VC.INPUT]
    diff_vars = category_dict[VC.DIFFERENTIAL]
    deriv_vars = category_dict[VC.DERIVATIVE]
    fixed_vars = category_dict[VC.FIXED]
    alg_vars = category_dict[VC.ALGEBRAIC]
    meas_vars = category_dict[VC.MEASUREMENT]

    assert len(input_vars) == len(init_input_set)
    for v in input_vars:
        assert v[0] in init_input_set

    assert len(deriv_vars) == len(init_deriv_set)
    for v in deriv_vars:
        assert v[0] in init_deriv_set

    assert len(diff_vars) == len(init_deriv_set)
    for v in diff_vars:
        assert v[0] in init_diff_set

    assert len(fixed_vars) == len(init_fixed_set)
    for v in fixed_vars:
        assert v[0] in init_fixed_set

    assert len(alg_vars) == len(init_alg_set)
    for v in alg_vars:
        assert v[0] in init_alg_set

    assert len(meas_vars) == len(init_meas_set)
    for v in meas_vars:
        assert v[0] in init_meas_set

    assert len(scalar_vars) == 0


@pytest.mark.unit
def test_categorize_4():
    """
    This tests categorization when a psuedo-steady state
    approximation is used. Energy accumulation and accumulation of E
    are fixed, the corresponding initial conditions unfixed, and
    the corresponding discretization equations deactivated.
    """
    model = make_model(horizon=1, ntfe=5, ntcp=2)
    time = model.fs.time

    CV = model.fs.cstr.control_volume
    CV.energy_accumulation[:, "aq"].fix(0.0)
    CV.material_accumulation[:, "aq", "E"].fix(0.0)
    CV.energy_holdup[0, "aq"].unfix()
    CV.material_holdup[0, "aq", "E"].unfix()
    CV.energy_accumulation_disc_eq.deactivate()
    CV.material_accumulation_disc_eq.deactivate()

    init_input_list = [
        model.fs.mixer.S_inlet.flow_vol[0],
        model.fs.mixer.E_inlet.flow_vol[0],
    ]
    init_input_set = ComponentSet(init_input_list)

    init_deriv_list = [
        CV.material_accumulation[0, "aq", "C"],
        CV.material_accumulation[0, "aq", "S"],
        CV.material_accumulation[0, "aq", "P"],
        CV.material_accumulation[0, "aq", "Solvent"],
    ]
    init_deriv_set = ComponentSet(init_deriv_list)

    init_diff_list = [
        CV.material_holdup[0, "aq", "C"],
        CV.material_holdup[0, "aq", "S"],
        CV.material_holdup[0, "aq", "P"],
        CV.material_holdup[0, "aq", "Solvent"],
    ]
    init_diff_set = ComponentSet(init_diff_list)

    init_fixed_list = [
        model.fs.cstr.control_volume.energy_accumulation[0, "aq"],
        CV.material_accumulation[0, "aq", "E"],
        model.fs.mixer.E_inlet.temperature[0],
        model.fs.mixer.S_inlet.temperature[0],
        *list(model.fs.mixer.E_inlet.conc_mol[0, :]),
        *list(model.fs.mixer.S_inlet.conc_mol[0, :]),
        model.fs.cstr.outlet.conc_mol[0, "Solvent"],
        model.fs.mixer.outlet.conc_mol[0, "Solvent"],
    ]
    init_fixed_set = ComponentSet(init_fixed_list)

    init_meas_list = [
        model.fs.cstr.control_volume.volume[0],
        CV.material_holdup[0, "aq", "P"],
        CV.material_holdup[0, "aq", "C"],
        CV.material_holdup[0, "aq", "S"],
    ]
    init_meas_set = ComponentSet(init_meas_list)

    init_alg_list = [
        model.fs.cstr.control_volume.energy_holdup[0, "aq"],
        CV.material_holdup[0, "aq", "E"],
        model.fs.cstr.outlet.flow_vol[0],
        model.fs.cstr.outlet.temperature[0],
        model.fs.cstr.inlet.flow_vol[0],
        model.fs.cstr.inlet.temperature[0],
        model.fs.mixer.outlet.flow_vol[0],
        model.fs.mixer.outlet.temperature[0],
        model.fs.cstr.control_volume.volume[0],
        *list(model.fs.cstr.control_volume.properties_out[0].flow_mol_comp[:]),
        *list(model.fs.cstr.inlet.conc_mol[0, :]),
        *list(model.fs.cstr.control_volume.properties_in[0].flow_mol_comp[:]),
        *list(model.fs.cstr.control_volume.rate_reaction_generation[0, "aq", :]),
        *list(model.fs.mixer.mixed_state[0].flow_mol_comp[:]),
        *list(model.fs.mixer.E_inlet_state[0].flow_mol_comp[:]),
        *list(model.fs.mixer.S_inlet_state[0].flow_mol_comp[:]),
        *list(model.fs.cstr.outlet.conc_mol[0, :]),
        *list(model.fs.mixer.outlet.conc_mol[0, :]),
        *list(model.fs.cstr.control_volume.reactions[0].reaction_coef[:]),
        *list(model.fs.cstr.control_volume.reactions[0].reaction_rate[:]),
        *list(model.fs.cstr.control_volume.rate_reaction_extent[0, :]),
    ]
    init_alg_set = ComponentSet(init_alg_list)
    # solvent outlet concentrations are not algebraic variables as
    # it is fixed.
    init_alg_set.remove(model.fs.cstr.outlet.conc_mol[0, "Solvent"])
    init_alg_set.remove(model.fs.mixer.outlet.conc_mol[0, "Solvent"])

    scalar_vars, dae_vars = flatten_dae_components(model, time, ctype=Var)
    category_dict = categorize_dae_variables(dae_vars, time, init_input_list)

    input_vars = category_dict[VC.INPUT]
    diff_vars = category_dict[VC.DIFFERENTIAL]
    deriv_vars = category_dict[VC.DERIVATIVE]
    fixed_vars = category_dict[VC.FIXED]
    alg_vars = category_dict[VC.ALGEBRAIC]
    meas_vars = category_dict[VC.MEASUREMENT]

    assert len(input_vars) == len(init_input_set)
    for v in input_vars:
        assert v[0] in init_input_set

    assert len(deriv_vars) == len(init_deriv_set)
    for v in deriv_vars:
        assert v[0] in init_deriv_set

    assert len(diff_vars) == len(init_deriv_set)
    for v in diff_vars:
        assert v[0] in init_diff_set

    assert len(fixed_vars) == len(init_fixed_set)
    for v in fixed_vars:
        assert v[0] in init_fixed_set

    assert len(alg_vars) == len(init_alg_set)
    for v in alg_vars:
        assert v[0] in init_alg_set

    assert len(meas_vars) == len(init_meas_set)
    for v in meas_vars:
        assert v[0] in init_meas_set

    assert len(scalar_vars) == 0


@pytest.mark.unit
def test_categorize_error():
    model = make_model(horizon=1, ntfe=5, ntcp=2)
    time = model.fs.time

    scalar_vars, dae_vars = flatten_dae_components(model, time, ctype=Var)

    # Add a dummy var to treat as an input.
    # This var is not in `dae_vars`, so it will not be located during
    # categorization, which should fail.
    model.dummy_var = Var(time)

    init_input_list = [
        model.fs.mixer.S_inlet.flow_vol[0],
        model.fs.mixer.E_inlet.flow_vol[0],
        model.dummy_var[0],
    ]

    with pytest.raises(RuntimeError, match=r"Not all inputs could be found"):
        category_dict = categorize_dae_variables(
            dae_vars,
            time,
            init_input_list,
        )

    # Re-run flattener. Now `dummy_var` should be included in `dae_vars`.
    scalar_vars, dae_vars = flatten_dae_components(model, time, ctype=Var)
    category_dict = categorize_dae_variables(
        dae_vars,
        time,
        init_input_list,
    )


if __name__ == "__main__":
    test_categorize_1()
    test_categorize_2()
    test_categorize_3()
    test_categorize_4()
    test_categorize_error()
