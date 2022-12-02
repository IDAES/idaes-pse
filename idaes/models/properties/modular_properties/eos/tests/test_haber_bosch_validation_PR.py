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
Tests for Cubic equation of state methods

Author: Douglas Allan

Cross reference calculated values of partial molar enthalpies, volumes, and
chemical potentials against values provided in the supplemental material of
(Rahbari, et al., 2018). Because the agreement for those values is looser than
preferable, cross-cross reference values with CoolProp while also validating
other thermodynamic properties, like specific entropy.

Ultimately, this test doesn't give, tight, independent validation of partial
molar enthalpy or entropy, but it does validate partial molar volumes and
chemical potentials. Special thanks to Rahbari et al. for including such
detailed supplementary material to make an apples-to-apples comparison possible

Rahbari, A., Hens, R., Nikolaidis, I. K., Poursaeidesfahani, A., Ramdin, M.,
Economou, I. G., … Vlugt, T. J. H. (2018). Computation of partial molar
properties using continuous fractional component Monte Carlo. Molecular
Physics, 116(21–22), 3331–3344. https://doi.org/10.1080/00268976.2018.1451663

"""

from copy import deepcopy
from io import StringIO
from inspect import cleandoc

import pandas as pd

import pytest
from pytest import approx

from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    Constraint,
    log,
    Reals,
    SolverFactory,
    value,
    Var,
    units as pyunits,
)


from idaes.core import VaporPhase, Component, PhaseType
from idaes.models.properties.modular_properties.eos.ceos import Cubic, CubicType
from idaes.models.properties.modular_properties.eos.ideal import Ideal


from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)


from idaes.models.properties.modular_properties.pure import NIST
from idaes.models.properties.modular_properties.state_definitions import FTPx


try:
    import CoolProp.CoolProp as CP
except:
    # If CoolProp isn't available, the test is skipped later on
    pass

standard_temp = 298.15
standard_pressure = 1e5


def _CoolProp_available():
    try:
        import CoolProp

        return True
    except:
        return False


_phase_dicts = {
    "Vap": {
        "type": VaporPhase,
        "equation_of_state": Cubic,
        "equation_of_state_options": {"type": CubicType.PR},
    }
}

_component_params = {
    "H2": {
        "type": Component,
        "valid_phase_types": [PhaseType.vaporPhase],
        "elemental_composition": {"H": 2},
        "enth_mol_ig_comp": NIST,
        "entr_mol_ig_comp": NIST,
        "cp_mol_ig_comp": NIST,
        "parameter_data": {
            "mw": (0.0020159, pyunits.kg / pyunits.mol),
            "pressure_crit": (1296400, pyunits.Pa),
            "temperature_crit": (33.14, pyunits.K),
            # This isn't the best value of the acentric factor of hydrogen to
            # use, but it's what the reference uses
            "omega": -0.219,
            "kappa1": 0,
            "cp_mol_ig_comp_coeff": {
                "A": 33.066178,
                "B": -11.363417,
                "C": 11.432816,
                "D": -2.772874,
                "E": -0.158558,
                "F": -9.980797,
                "G": 172.707974,
                "H": 0.0,
            },
        },
    },
    "N2": {
        "type": Component,
        "valid_phase_types": [PhaseType.vaporPhase],
        "elemental_composition": {"N": 2},
        "enth_mol_ig_comp": NIST,
        "entr_mol_ig_comp": NIST,
        "cp_mol_ig_comp": NIST,
        "parameter_data": {
            "mw": (0.0280134, pyunits.kg / pyunits.mol),
            "pressure_crit": (3395800, pyunits.Pa),
            "temperature_crit": (126.19, pyunits.K),
            "omega": 0.0372,
            "cp_mol_ig_comp_coeff": {
                "A": 19.50583,
                "B": 19.88705,
                "C": -8.598535,
                "D": 1.369784,
                "E": 0.527601,
                "F": -4.935202,
                "G": 212.39,
                "H": 0.0,
            },
        },
    },
    "NH3": {
        "type": Component,
        "valid_phase_types": [PhaseType.vaporPhase],
        "elemental_composition": {"N": 1, "H": 3},
        "enth_mol_ig_comp": NIST,
        "entr_mol_ig_comp": NIST,
        "cp_mol_ig_comp": NIST,
        "parameter_data": {
            "mw": (0.017031, pyunits.kg / pyunits.mol),
            "pressure_crit": (11333000, pyunits.Pa),
            "temperature_crit": (405.4, pyunits.K),
            "omega": 0.25601,
            "cp_mol_ig_comp_coeff": {
                "A": 19.99563,
                "B": 49.77119,
                "C": -15.37599,
                "D": 1.921168,
                "E": 0.189174,
                "F": -53.30667,
                "G": 203.8591,
                "H": -45.89806,
            },
        },
    },
}

# Shamelessly stolen from natural_gas_PR.py
# returns a configuration dictionary for the list of specified components
def _get_prop(components=None, phases="Vap"):
    if components is None:
        components = list(_component_params.keys())
    configuration = {
        "components": {},  # fill in later based on selected components
        "parameter_data": {},
        "phases": {},
        # Set base units of measurement
        "base_units": {
            "time": pyunits.s,
            "length": pyunits.m,
            "mass": pyunits.kg,
            "amount": pyunits.mol,
            "temperature": pyunits.K,
        },
        # Specifying state definition
        "state_definition": FTPx,
        "state_bounds": {
            "flow_mol": (0, 1, 10, pyunits.mol / pyunits.s),
            "temperature": (253.15, 573, 1000, pyunits.K),
            "pressure": (100, 40e6, 1e9, pyunits.Pa),
        },
        "pressure_ref": (1e5, pyunits.Pa),
        "temperature_ref": (298.15, pyunits.K),
    }

    c = configuration["components"]
    for comp in components:
        c[comp] = deepcopy(_component_params[comp])
    if isinstance(phases, str):
        phases = [phases]
    for k in phases:
        configuration["phases"][k] = deepcopy(_phase_dicts[k])

    # Fill the binary parameters with zeros.
    d = configuration["parameter_data"]
    d["PR_kappa"] = {(a, b): 0 for a in c for b in c}
    return configuration


@pytest.fixture()
def model():
    m = ConcreteModel()

    configuration = _get_prop(["H2", "N2", "NH3"], ["Vap"])

    m.params = GenericParameterBlock(**configuration)
    m.props = m.params.state_block_class(
        defined_state=True, parameters=m.params, has_phase_equilibrium=False
    )
    configuration["phases"]["Vap"]["equation_of_state"] = Ideal
    m.params_IG = GenericParameterBlock(**configuration)
    m.props_IG = m.params.state_block_class(
        defined_state=True, parameters=m.params_IG, has_phase_equilibrium=False
    )

    Ntot = 18 + 53 + 385
    for prop in [m.props, m.props_IG]:
        prop.flow_mol.fix(1)
        prop.temperature.fix(573)
        prop.pressure.fix(40e6)

        prop.mole_frac_comp["N2"].fix(18 / Ntot)
        prop.mole_frac_comp["H2"].fix(53 / Ntot)
        prop.mole_frac_comp["NH3"].fix(385 / Ntot)
        prop.initialize()
    return m


@pytest.fixture()
def dataframe():
    # Need to use cleandoc to get rid of hanging indents inside the string
    string = cleandoc(
        """\
        P,N_NH3,N_N2,N_H2,H_NH3,H_res_NH3,V_NH3,H_N2,H_res_N2,V_N2,H_H2,H_res_H2,V_H2,G_res_NH3,G_res_N2,G_res_H2,dH_rxn
        10,306,57,171,-36.61,-1.65,0.441,8.28,0.2,0.504,8.57,0.56,0.501,-0.38,0.26,0.25,-107.23
        20,349,36,107,-38.43,-3.47,0.205,9.1,1.02,0.274,9.76,1.74,0.27,-0.79,0.68,0.64,-115.26
        30,371,25,74,-40.13,-5.17,0.13,10.45,2.37,0.2,11.36,3.35,0.193,-1.14,1.19,1.11,-124.81
        40,385,18,53,-41.56,-6.6,0.096,12.01,3.93,0.162,13.02,5.01,0.153,-1.43,1.75,1.59,-134.23
        50,394,13,40,-42.65,-7.7,0.078,13.41,5.33,0.136,14.4,6.38,0.126,-1.65,2.3,2.03,-141.94
        60,400,10,31,-43.45,-8.49,0.068,14.46,6.38,0.118,15.35,7.34,0.107,-1.81,2.8,2.42,-147.44
        70,404,8,25,-44.02,-9.06,0.061,15.18,7.1,0.103,15.97,7.95,0.092,-1.91,3.26,2.76,-151.14
        80,407,7,20,-44.42,-9.46,0.056,15.68,7.6,0.093,16.36,8.34,0.081,-1.97,3.67,3.05,-153.61
        """
    )

    return pd.read_csv(StringIO(string))


# TODO: At present, this method is unnecessary because we're comparing residual
# entropies, not full entropies. We can revisit it once CoolProp patches
# current issues it has with its implementation---either use it or delete it
def get_reference_entropy(comp):
    m = ConcreteModel()
    m.params = GenericParameterBlock(**_get_prop([comp]))
    m.props = m.params.state_block_class(
        defined_state=True, parameters=m.params, has_phase_equilibrium=False
    )

    # Want to make sure intermediate quantities are constructed before we
    # solve the block. Therefore introduce this variable and constraint

    m.S_ref = Var(domain=Reals, initialize=1, units=pyunits.J / pyunits.mol)

    def rule_S_ref(blk, j):
        return m.S_ref == m.props.entr_mol

    m.S_ref_eq = Constraint(m.params.component_list, rule=rule_S_ref)

    m.props.mole_frac_comp[comp].fix(1)
    m.props.temperature.fix(standard_temp)
    m.props.pressure.fix(standard_pressure)

    solver = SolverFactory("ipopt")
    results = solver.solve(m)

    assert check_optimal_termination(results)
    assert value(m.S_ref) == approx(
        value(m.props.entr_mol_phase_comp["Vap", comp]), rel=1e-12
    )

    return value(m.S_ref)


# TODO: This class needs to be revisited after CoolProp fixes current
# implementation bugs. Right now it wraps CoolProp's AbsractState class
# in order to perform calculation corrections that are necessary.
class CoolPropTester:
    component_list = ["H2", "N2", "NH3"]
    Href = {}
    Sref = {}

    def __init__(self, mole_frac_comp):
        # You're supposed to be able to set reference states in CoolProp, but
        # that functionality is broken for cubic EoSes at the moment
        for comp in self.component_list:
            self.Href[comp] = CP.PropsSI(
                "HMOLAR", "T", standard_temp, "P", standard_pressure, "PR::" + comp
            )
            self.Sref[comp] = CP.PropsSI(
                "SMOLAR", "T", standard_temp, "P", standard_pressure, "PR::" + comp
            )
        self.Href["NH3"] -= -45949  # Heat of formation in J/mol
        self.AS = CP.AbstractState("PR", "H2&N2&NH3")
        self.set_mole_frac_comp(mole_frac_comp)
        self.AS.update(CP.PT_INPUTS, 4e7, 573)

    def set_temperature(self, temperature):
        self.AS.update(CP.PT_INPUTS, self.AS.p(), temperature)
        self._set_phase()

    def set_pressure(self, pressure):
        self.AS.update(CP.PT_INPUTS, pressure, self.AS.T())
        self._set_phase()

    def set_mole_frac_comp(self, mole_frac_comp):
        self.mole_frac_comp = mole_frac_comp
        self.AS.set_mole_fractions(
            [mole_frac_comp[comp] for comp in self.component_list]
        )

        # Since the phase doesn't appear to matter, all this code does is
        # add extra cycles

        # pts = self.AS.all_critical_points()
        # assert len(pts) < 2
        # if len(pts) == 1:
        #     self.Tcrit = pts[0].T
        #     self.Pcrit = pts[0].p
        # else:
        #     self.Tcrit = float('inf')
        #     self.Pcrit = float('inf')
        #
        # self._set_phase()

    def _set_phase(self):
        self.AS.specify_phase(CP.iphase_gas)

    # I don't think there is any distinction in CoolProp between a gas,
    # a supercritical gas, and a supercritical fluid for a cubic EoS
    # Keep this code around in case I'm wrong---Doug

    # if self.AS.T() < self.Tcrit and self.AS.p() < self.Pcrit:
    #     self.AS.specify_phase(CP.iphase_gas)
    # elif self.AS.T() > self.Tcrit and self.AS.p() < self.Pcrit:
    #     self.AS.specify_phase(CP.iphase_supercritical_gas)
    # else:
    #     self.AS.specify_phase(CP.iphase_supercritical)

    def enth_mol(self):
        return self.AS.hmolar() - sum(
            self.mole_frac_comp[comp] * self.Href[comp] for comp in self.component_list
        )

    def entr_mol(self):
        return self.AS.smolar() - sum(
            self.mole_frac_comp[comp] * self.Sref[comp] for comp in self.component_list
        )

    def gibbs_mol(self):
        return self.AS.gibbsmolar() - sum(
            self.mole_frac_comp[comp]
            * (self.Href[comp] - self.AS.T() * self.Sref[comp])
            for comp in self.component_list
        )

    def enth_res_mol(self):
        return self.AS.hmolar_residual()

    def gibbs_res_mol_comp(self, comp):
        # CoolProp's gmolar_residual uses an ideal gas reference state with
        # the same *density*, whereas we use one with the same *pressure*.
        # Therefore get the residual gibbs energy by summing excess chemical
        # potentials from the fugacity coefficients
        return (
            8.314472
            * self.AS.T()
            * log(self.AS.fugacity_coefficient(self.component_list.index(comp)))
        )

    def entr_res_mol(self):
        # CoolProp's smolar_residual uses an ideal gas reference state with
        # the same *density*, whereas we use one with the same *pressure*.
        # Therefore this backwards calculation
        return (self.enth_res_mol() - self.gibbs_res_mol()) / self.AS.T()

    def gibbs_res_mol(self):
        return sum(
            self.mole_frac_comp[comp] * self.gibbs_res_mol_comp(comp)
            for comp in self.component_list
        )

    def vol_mol(self):
        return 1 / self.AS.rhomolar()


# Subdividing the tests really doesn't make sense in this case
@pytest.mark.integration
@pytest.mark.skipif(not _CoolProp_available(), reason="CoolProp not available")
def test_thermo(model, dataframe):
    m = model
    Ntot = 18 + 53 + 385
    tester = CoolPropTester({"H2": 53 / Ntot, "N2": 18 / Ntot, "NH3": 385 / Ntot})
    tester.set_temperature(573)
    df = dataframe
    solver = SolverFactory("ipopt")

    # Correct CoolProp's reference entropy to match what the Shomate
    # equation gives at 298.15 K
    for comp in tester.component_list:
        tester.Sref[comp] -= get_reference_entropy(comp)

    # Construct these variables and constraints to make sure that all
    # intermediate variables and constraints are created before we solve
    m.H_res = Var(m.params.component_list, initialize=0, units=pyunits.kJ / pyunits.mol)
    m.V = Var(m.params.component_list, initialize=1, units=pyunits.L / pyunits.mol)
    m.G_res = Var(m.params.component_list, initialize=0, units=pyunits.kJ / pyunits.mol)

    def rule_H_res(blk, j):
        return m.H_res[j] == (
            1e-3
            * (pyunits.kJ / pyunits.J)
            * (
                m.props.enth_mol_phase_comp["Vap", j]
                - m.props_IG.enth_mol_phase_comp["Vap", j]
            )
        )

    def rule_G_res(blk, j):
        return m.G_res[j] == (
            1e-3
            * (pyunits.kJ / pyunits.J)
            * (
                m.props.gibbs_mol_phase_comp["Vap", j]
                - m.props_IG.gibbs_mol_phase_comp["Vap", j]
            )
        )

    def rule_V(blk, j):
        return m.V[j] == (
            1e3 * (pyunits.L / pyunits.m**3) * m.props.vol_mol_phase_comp["Vap", j]
        )

    m.H_res_eq = Constraint(m.params.component_list, rule=rule_H_res)
    m.G_res_eq = Constraint(m.params.component_list, rule=rule_G_res)
    m.V_eq = Constraint(m.params.component_list, rule=rule_V)

    mole_frac_comp = {}
    for idx, row in df.iterrows():
        Ntot = row["N_N2"] + row["N_H2"] + row["N_NH3"]
        for comp in tester.component_list:
            mole_frac_comp[comp] = row["N_" + comp] / Ntot
        tester.set_mole_frac_comp(mole_frac_comp)

        for prop in [m.props, m.props_IG]:
            for comp in tester.component_list:
                prop.mole_frac_comp[comp].fix(mole_frac_comp[comp])
            prop.pressure.fix(row["P"] * 1e6)
        tester.set_pressure(row["P"] * 1e6)
        results = solver.solve(m)
        assert check_optimal_termination(results)

        # Make sure that log mole fractions have been created and validated
        for comp in tester.component_list:
            assert value(m.props.log_mole_frac_phase_comp["Vap", comp]) == approx(
                value(log(m.props.mole_frac_phase_comp["Vap", comp])), rel=1e-6
            )

        # Assert thermodynamic consistency---that all the partial molar
        # quantities add up to the specific molar quantity
        # This is nontrivial---the specific molar quantities are calculated
        # with different equations than the partial molar quantities

        assert approx(value(m.props.entr_mol_phase["Vap"]), rel=1e-6) == value(
            sum(
                m.props.mole_frac_comp[comp] * m.props.entr_mol_phase_comp["Vap", comp]
                for comp in tester.component_list
            )
        )

        assert approx(value(m.props.enth_mol_phase["Vap"]), rel=1e-6) == value(
            sum(
                m.props.mole_frac_comp[comp] * m.props.enth_mol_phase_comp["Vap", comp]
                for comp in tester.component_list
            )
        )

        assert approx(value(m.props.gibbs_mol_phase["Vap"]), rel=1e-6) == value(
            sum(
                m.props.mole_frac_comp[comp] * m.props.gibbs_mol_phase_comp["Vap", comp]
                for comp in tester.component_list
            )
        )

        assert approx(value(m.props.vol_mol_phase["Vap"]), rel=1e-6) == value(
            sum(
                m.props.mole_frac_comp[comp] * m.props.vol_mol_phase_comp["Vap", comp]
                for comp in tester.component_list
            )
        )

        for comp in tester.component_list:
            assert approx(
                value(m.props.gibbs_mol_phase_comp["Vap", comp]), rel=1e-6
            ) == value(
                m.props.enth_mol_phase_comp["Vap", comp]
                - m.props.temperature * m.props.entr_mol_phase_comp["Vap", comp]
            )

        assert (
            approx(
                value(m.props.entr_mol_phase["Vap"] - m.props_IG.entr_mol_phase["Vap"]),
                rel=1e-4,
            )
            == tester.entr_res_mol()
        )

        print("============================================")
        print("Temperature: {0:.1f} K".format(value(m.props.temperature)))
        print("Pressure: {0:.1f} MPa".format(1e-6 * value(m.props.pressure)))
        for comp in tester.component_list:
            print(
                "Mole fraction "
                + comp
                + ": {0:.3f}".format(value(m.props.mole_frac_comp[comp]))
            )

        print("")
        print("Excess enthalpy (kJ/mol):")
        print(
            "IDAES:    {0:.4f}".format(
                1e-3
                * value(
                    m.props.enth_mol_phase["Vap"] - m.props_IG.enth_mol_phase["Vap"]
                )
            )
        )
        print("CoolProp: {0:.4f}".format(1e-3 * tester.enth_res_mol()))
        print(
            "Paper:    {0:.4f}".format(
                sum(
                    mole_frac_comp[comp] * row["H_res_" + comp]
                    for comp in tester.component_list
                )
            )
        )

        assert (
            approx(
                value(m.props.enth_mol_phase["Vap"] - m.props_IG.enth_mol_phase["Vap"]),
                rel=1e-4,
            )
            == tester.enth_res_mol()
        )

        print("")
        print("Partial molar excess enthalpy (kJ/mol):")
        print("         IDAES  Paper")
        for comp in tester.component_list:
            # I do not like this loose tolerance, but this paper used a
            # possibly sketchy numerical integration scheme to calculate
            # H_res instead of differentiating the fugacity coefficient like
            # we did. The fact that their residual chemical potentials are also
            # off beyond what we'd expect from sig-figs I think vindicates our
            # calculations in general, but it would be nice to have an independent
            # calculation of H_res that matched a bit more closely
            assert approx(value(m.H_res[comp]), rel=7.5e-2) == row["H_res_" + comp]
            if comp == "NH3":
                print(
                    comp
                    + ":  {0:.2f}  {1:.2f}".format(
                        value(m.H_res[comp]), row["H_res_" + comp]
                    )
                )
            else:
                print(
                    comp
                    + ":    {0:.2f}   {1:.2f}".format(
                        value(m.H_res[comp]), row["H_res_" + comp]
                    )
                )
        print("")

        print("Excess Gibbs energy (kJ/mol):")
        print(
            "IDAES:    {0:.4f}".format(
                1e-3
                * value(
                    m.props.gibbs_mol_phase["Vap"] - m.props_IG.gibbs_mol_phase["Vap"]
                )
            )
        )
        print("CoolProp: {0:.4f}".format(1e-3 * tester.gibbs_res_mol()))
        print(
            "Paper:    {0:.4f}".format(
                sum(
                    mole_frac_comp[comp] * row["G_res_" + comp]
                    for comp in tester.component_list
                )
            )
        )
        assert (
            approx(
                value(
                    m.props.gibbs_mol_phase["Vap"] - m.props_IG.gibbs_mol_phase["Vap"]
                ),
                rel=1e-4,
            )
            == tester.gibbs_res_mol()
        )

        print("")
        print("Excess chemical potential (kJ/mol):")
        print("         IDAES  CoolProp  Paper")
        for comp in tester.component_list:
            assert approx(
                value(m.G_res[comp]), rel=1e-4
            ) == 1e-3 * tester.gibbs_res_mol_comp(comp)
            if comp == "NH3":
                print(
                    comp
                    + ":  {0:.4f}  {1:.4f}  {2:.4f}".format(
                        value(m.G_res[comp]),
                        1e-3 * tester.gibbs_res_mol_comp(comp),
                        row["G_res_" + comp],
                    )
                )
            else:
                print(
                    comp
                    + ":    {0:.4f}   {1:.4f}   {2:.4f}".format(
                        value(m.G_res[comp]),
                        1e-3 * tester.gibbs_res_mol_comp(comp),
                        row["G_res_" + comp],
                    )
                )
        print("")

        print("")
        print("Specific molar volume (L/mol):")
        print("IDAES:    {0:.4f}".format(1e3 * value(m.props.vol_mol_phase["Vap"])))
        print("CoolProp: {0:.4f}".format(1e3 * tester.vol_mol()))
        print(
            "Paper:    {0:.4f}".format(
                sum(
                    mole_frac_comp[comp] * row["V_" + comp]
                    for comp in tester.component_list
                )
            )
        )
        assert approx(value(m.props.vol_mol_phase["Vap"]), rel=1e-4) == tester.vol_mol()
        print("")
        print("Partial molar volume (L/mol):")
        print("         IDAES  Paper")
        for comp in tester.component_list:
            assert approx(value(m.V[comp]), rel=1e-2) == row["V_" + comp]
            if comp == "NH3":
                print(
                    comp
                    + ":   {0:.4f}  {1:.4f}".format(value(m.V[comp]), row["V_" + comp])
                )
            else:
                print(
                    comp
                    + ":    {0:.4f}  {1:.4f}".format(value(m.V[comp]), row["V_" + comp])
                )
        print("")
        print("")
