#######################################################################
# Thermodynamic Framework with application to MEA-H2O-CO2 system
# 1 .Symmetric E-NRTL model
# 2. Excess enthalpy model
# 3. Enthalpy model
# 4  Chemical equilibruim for MEA-CO2-H2O systems
# 5  Physical equlibruim for MEA-CO2-H2O system(CO2 partial pressure)

########################################################################

# import Pyomo libraries
from pyomo.environ import (ConcreteModel,
                           Var, Param,
                           Reals,
                           exp, log,
                           value, Constraint, Expression)

# import idaes libraies
from idaes.core.solvers import get_solver

# import third party libraries
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from idaes.models_extra.column_models.properties.MEA_solvent_eNRTL import (
    configuration as liquid_config,
)
from idaes.models.properties.modular_properties import GenericParameterBlock

#==============================================================================
#NMR (Speciation) Data from Hilliard and Buttinger

#MEACOO lean loading from HILLIARD Data
MEACOO_Hilliard_x = [0.018537701,    0.022108392, 0.025679084, 0.029696112, 0.041747197, 0.051120262, 0.068973721, 0.098431927, 0.119856078, 0.149314284, 0.168953088, 0.197964958, 0.218050099, 0.247508306, 0.268039783, 0.296605316, 0.317136794, 0.346595,    0.366233805, 0.395692011, 0.416223488 ,0.445681695, 0.466213172, 0.495225042, 0.515310183, 0.544768389, 0.564407194, 0.584045998, 0.5938654,   0.613504204, 0.623323607, 0.633143009, 0.642962411, 0.652781813, 0.662601215, 0.672420617, 0.68224002,  0.692059422, 0.701953213, 0.711698226, 0.721517628, 0.741156433, 0.780344774, 0.790253443, 0.799924067, 0.759136213, 0.770653378]

# MEACOO speciation from Hilliard Data
MEACOO_Hilliard_y = [0.002689693,    0.003255221, 0.003923521, 0.004769957, 0.007172606, 0.008563292, 0.011550161, 0.016280013, 0.019928495, 0.024540686, 0.027922143, 0.032878914, 0.036147048, 0.040822002, 0.044247891, 0.049022623, 0.052180641, 0.056728841, 0.060218309, 0.064980415, 0.068108002, 0.072405391, 0.075494994, 0.078224594, 0.077608323, 0.074210408, 0.071492735, 0.067601691, 0.066949109, 0.063779488, 0.061696488, 0.060488513, 0.059149512, 0.057538828, 0.055805196, 0.054164197, 0.053242566, 0.051139896, 0.050658824, 0.047889768, 0.046690095, 0.04313945,  0.036551172 ,0.035005666, 0.033583512, 0.04015814,  0.038870725]

# MEAH lean laoding from Hilliard Data
MEAH_Hilliard_x = [0.0094126,    0.034159477, 0.059154319, 0.083256488, 0.107358657, 0.133246171, 0.15734834,  0.182343182, 0.207338024, 0.232332866, 0.257327708, 0.280537204, 0.305532046, 0.330526888, 0.355521729, 0.380516571, 0.40461874,  0.429613582, 0.457732779, 0.471658477, 0.491208014, 0.521038167, 0.540305025, 0.550124427, 0.559943829, 0.569763231, 0.589402035, 0.618860242, 0.638499046, 0.667957253, 0.687596057, 0.717054264, 0.736693068, 0.766151274, 0.785418132]

# MEAH speciation from Hilliard Data
MEAH_Hilliard_y = [0.001759266,  0.006111763, 0.010498049, 0.014896761, 0.019084427, 0.023641977, 0.027901316, 0.032079201, 0.036663164, 0.040974588, 0.045215804, 0.049267093, 0.054002775, 0.058186453, 0.062274632, 0.066999083, 0.070855443, 0.075213857, 0.079535627, 0.0801089,   0.082259266, 0.086309693, 0.089763469, 0.091658607, 0.092897986, 0.094647195, 0.097952251, 0.102666812, 0.105777281, 0.110951229, 0.113632347, 0.118835298, 0.121072891, 0.126900277, 0.127636246]


#MEACOO LLD from BUTTINGER Data
MEACOO_Buttinger_x = [0.125, 0.205, 0.236, 0.28, 0.343, 0.354, 0.46, 0.466,
                      0.487, 0.508, 0.663, 0.683, 0.762]

#MEACOO speciation from Buttinger Data
MEACOO_Buttinger_y = [0.0144, 0.0241,  0.0271,  0.0327,  0.0396,  0.041,
   0.0492,  0.0502,  0.0511,  0.0483,  0.0347,  0.0323,  0.0188]

#MEAH LLD from BUTTINGER Data
MEAH_Buttinger_x = [0.125, 0.205, 0.236, 0.28, 0.343, 0.354, 0.46, 0.466,
                      0.487, 0.508, 0.663, 0.683, 0.762]

#MEAH speciation from Buttinger Data
MEAH_Buttinger_y = [0.1009,  0.0934,  0.0876,  0.084,   0.0758,  0.0749,
  0.065,   0.0663,  0.0646,  0.0647,  0.0759,  0.0739,  0.0807]


#==============================================================================
# MEA-H2O binary parameters

A_MEA_H2O = {
    "Liu": -1.609693,
    "Zhang": 1.5201,
    "Morgan": 3.25,
    "Hessen": 1.998,
    "Schmidt": 0.0607,
    "Austgen": 0.0}

B_MEA_H2O = {
    "Liu": -265.1962,
    "Zhang": -910.3,
    "Morgan": 0,
    "Hessen": -1030,
    "Schmidt": -887.3,
    "Austgen": -649.75}

A_H2O_MEA = {
    "Liu": 1.648945,
    "Zhang": 0.1559,
    "Morgan": 4.34,
    "Hessen": -0.5897,
    "Schmidt": 0.6011,
    "Austgen": 1.674}

B_H2O_MEA = {
    "Liu": 125.0832,
    "Zhang": 110.8,
    "Morgan": -2200,
    "Hessen": 250.6,
    "Schmidt": 1036,
    "Austgen": 0}
#=============================================================================


# create list and sets
class Enrtl_sets(object):

    '''
    This class contains Lists and Sets that will be used to create
    indexed Parameters,Variables and Expressions for the ENRTL model

    Lists will be first created and then converted to sets to perform
    set operations for creation of other super sets

    ===========================
    LISTS/Sets
    ===========================
    1. mm   -> mm',m'm,mm
    2. cm   -> cm,mc
    3. am   -> am,ma
    4. ca   -> ca,ac
    5. caca -> caca',cac'a,caca
    6 .cam  -> cam,mca
    7  E2   -> union of all sets with 2 indices (mm,cm,am,ca)
    '''

    # sets for CO2-MEA-H2O system application
    Ss = solvents = {'MEA', 'H2O'}
    Sa = anions = {'MEACOO-', 'HCO3-'}
    Sc = cations = {'MEAH+'}
    Sm = molecules = {'CO2', 'H2O', 'MEA'}
    Sma = Sm.union(Sa)
    Smc = Sm.union(Sc)
    Sca = Sc.union(Sa)
    Smac = Sm.union(Sca)
    Sh = Sm - Ss
    solutes = Smac - Ss

    # sets for ENRTL model
    # create a list of tuples and then convert to sets
    # LISTS
    # cam
    salt_molecule_interaction = []
    for i in Sc:
        for j in Sa:
            for k in Sm:
                # ca-m
                salt_molecule_interaction.append((i, j, k))
                # m-ca
                salt_molecule_interaction.append((k, i, j))

    # caca
    salt_salt_interaction = []

    for i in Sc:
        for j in Sa:
            # caca
            salt_salt_interaction.append((i, j, i, j))
            for k in Sa:
                if i + j != i + k:
                    # ca-ca'
                    salt_salt_interaction.append((i, j, i, k))
    for i in Sa:
        for j in Sc:
            for k in Sc:
                if i + j != i + k:
                    # ca-c'a
                    salt_salt_interaction.append((j, i, k, i))

    # mm
    molecule_molecule_interaction = []
    for i in Sm:
        for j in Sm:
            molecule_molecule_interaction.append((i, j))

    # cm
    cation_molecule_interaction = []
    for i in Sc:
        for j in Sm:
            # cm
            cation_molecule_interaction.append((i, j))
            # mc
            cation_molecule_interaction.append((j, i))

    # am
    anion_molecule_interaction = []
    for i in Sa:
        for j in Sm:
            # a-m
            anion_molecule_interaction.append((i, j))
            # m-a
            anion_molecule_interaction.append((j, i))

    # ca
    cation_anion_interaction = []
    for i in Sc:
        for j in Sa:
            # ca
            cation_anion_interaction.append((i, j))
            # ac
            cation_anion_interaction.append((j, i))

    mm = set(molecule_molecule_interaction)
    cm = set(cation_molecule_interaction)
    am = set(anion_molecule_interaction)
    ca = set(cation_anion_interaction)
    cam = set(salt_molecule_interaction)
    caca = set(salt_salt_interaction)
    E2 = mm.union(cm, am, ca)


class PhyCons(object):
    R_gas = 8.314*1e-3
    NA = 6.022137e23         # avogadro number
    kb = 1.381e-23  # Boltzmann constant
    e_charge = 1.6021773e-19  # absolute electron charge
    eo = 8.854e-12  # permitivity of free space
    rB = 3e-10


def thermodynamic_model(data):

    # create the concret model
    model = ConcreteModel()

    #--------------------------------------------------------------------------
    # INPUTS: Independent Varaibles To be Fixed From Data

    # Temperature
    model.T = Param(initialize=float(data['T']),doc='Temperature in K')

    # CO2 lean loading
    xa_id = ['ref','fin']
    lld_input ={
        'ref':float(data['lld']),
        'fin':float(data['lld']+ 0.001)
        }

    model.LLD = Param(xa_id,initialize=lld_input,doc='CO2 lean loading')

    # MEA wt fraction
    model.MEA_wt = Param(initialize=float(data['wt']), doc='MEA percent fraction')

    #--------------------------------------------------------------------------
    #Moles for CO2-MEA-H2O system
    model.n_MEA = Expression(expr= model.MEA_wt/61.08)
    model.n_H2O = Expression(expr= (100-model.MEA_wt)/18.015)

    def rule_nCO2(b,i):
        return b.LLD[i]*b.n_MEA
    model.n_CO2 = Expression(xa_id,rule=rule_nCO2)

    def rule_dnCO2(b):
        return 0.001*b.n_MEA
    model.dn_CO2 = Expression(rule=rule_dnCO2)

    def rule_nT(b,i):
        return b.n_CO2[i] + b.n_MEA + b.n_H2O
    model.nT = Expression(xa_id,rule=rule_nT)

    #--------------------------------------------------------------------------
    #apparent mole fraction
    def rule_xa(b,i,j):
        if j == 'CO2':
            return b.n_CO2[i]/b.nT[i]
        elif j == 'MEA':
            return b.n_MEA/b.nT[i]
        elif j == 'H2O':
            return b.n_H2O/b.nT[i]
    model.xa = Expression(xa_id,Enrtl_sets.Sm,rule=rule_xa)
    #--------------------------------------------------------------------------
    model.P = Param(initialize=107.325, doc='System Pressure in kPa')#107.325

    # create molecule-molecule interaction parameters(No regression)
    model.Amm = Param(Enrtl_sets.mm, mutable=True,
                      doc='Molecule-Molecule interaction parameters')
    model.Bmm = Param(Enrtl_sets.mm, mutable=True,
                      doc='Molecule-Molecule interaction parameters')
    for i in Enrtl_sets.Sm:
        for j in Enrtl_sets.Sm:
            if i == 'H2O' and j == 'MEA':
                model.Amm[i, j] = A_H2O_MEA['Hessen']
                model.Bmm[i, j] = B_H2O_MEA['Hessen']
            elif i == 'MEA' and j == 'H2O':
                model.Amm[i, j] = A_MEA_H2O['Hessen']
                model.Bmm[i, j] = B_MEA_H2O['Hessen']
            else:
                model.Amm[i, j] = 0.0
                model.Bmm[i, j] = 0.0

    # create salt-molecule interaction parameters as variables
    model.Ccam = Var(Enrtl_sets.salt_molecule_interaction, within=Reals,
                    doc='Salt-Molecule interaction parameters')
    model.Dcam = Var(Enrtl_sets.salt_molecule_interaction, within=Reals,
                     doc='Salt-Molecule interaction parameters')
    for i in Enrtl_sets.cam:
        if i == ('H2O', 'MEAH+', 'MEACOO-'):
            model.Ccam[i].fix(18.588108594566588)#6.732)
            model.Dcam[i].fix(-1533.6489230168652)#-9.2407)
        elif i == ('H2O', 'MEAH+', 'HCO3-'):
            model.Ccam[i].fix(8.5721)
            model.Dcam[i].fix(0.0)
        elif i == ('MEAH+', 'MEACOO-', 'H2O'):
            model.Ccam[i].fix(-7.055694714651596)#-3.1628)
            model.Dcam[i].fix(1224.0939968361579)#39.498)
        elif i == ('MEAH+', 'HCO3-', 'H2O'):
            model.Ccam[i].fix(-4.0092)
            model.Dcam[i].fix(0.0)
        else:
            if i[0] in Enrtl_sets.Sm:
                model.Ccam[i].fix(8)
                model.Dcam[i].fix(0.0)
            elif i[2] in Enrtl_sets.Sm:
                model.Ccam[i].fix(-4)
                model.Dcam[i].fix(0.0)

    # Energy interaction parameters with four indices
    model.t4 = Param(Enrtl_sets.caca, default=0.0,
                     doc='Salt-Salt interaction parameters')
    # create true species mole fraction with indices:
    # fin ---> final : used for phase equilibruim
    # ref ---> reference used for heat of absorption calculation
    xt = dict()
    xt[('fin', 'MEA')]=0.0038594305800293276
    xt[('fin', 'HCO3-')]=0.03665389616090885
    xt[('fin', 'MEAH+')]=0.06684528655689173
    xt[('fin', 'CO2')]=0.0340508209760122
    xt[('fin', 'MEACOO-')]=0.030191390395982877
    xt[('fin', 'H2O')]=0.761553888773283
    xt[('ref', 'MEA')]=0.0038643433557289295
    xt[('ref', 'HCO3-')]=0.03663207202113963
    xt[('ref', 'MEAH+')]=0.06683700862502923
    xt[('ref', 'CO2')]=0.03396837367103388
    xt[('ref', 'MEACOO-')]=0.030204936603889593
    xt[('ref', 'H2O')]=0.7616562570981495

    #get initial values for nt
    nt = dict()
    for i in xa_id:
        for j in Enrtl_sets.Smac:
            nt[i,j] = xt[i,j]*value(model.nT[i])

    model.xt = Var(xa_id, Enrtl_sets.Smac, bounds=(0, 1),
                   initialize=xt, doc='True species mole fraction')
    model.nt = Var(xa_id, Enrtl_sets.Smac,
                   initialize=nt, doc='True species mole number')

    def rule_nt_tot(b,p):
        return sum(b.nt[p,i] for i in Enrtl_sets.Smac)
    model.nt_tot = Expression(xa_id,rule=rule_nt_tot)

    def rule_xt(b,p,i):
        return b.xt[p,i]*b.nt_tot[p] == b.nt[p,i]
    model.eq_xt = Constraint(xa_id,Enrtl_sets.Smac,rule=rule_xt)

    # absolute ionic charge used for NRTL model
    # z = |z| for ions
    # z = 1 for molecules
    model.z = Param(Enrtl_sets.Smac, default=1,
                    doc='absolute value of ionic charge:NRTL')

    # Effective Mole fraction
    def rule_effective_mole_fraction(m, i, j):
        return m.xt[i, j] * m.z[j]

    model.X = Expression(xa_id, Enrtl_sets.Smac,
                         rule=rule_effective_mole_fraction)
    #==========================================================================
    # anionic charge fraction
    model.Ya = Var(xa_id, Enrtl_sets.Sa,
                   initialize=0.5, doc='anionic charge fraction')

    def rule_anionic_charge_fraction(m, i, a):
        return m.Ya[i, a] * sum(m.X[i, j] for j in Enrtl_sets.Sa) == m.X[i, a]

    model.eq_Ya = Constraint(xa_id, Enrtl_sets.Sa,
                             rule=rule_anionic_charge_fraction)

    # cationic charge fraction
    model.Yc = Var(xa_id, Enrtl_sets.Sc,
                   initialize=0.5, doc='cationic charge fraction')

    def rule_cationic_charge_fraction(m, i, c):
        return m.Yc[i, c] * sum(m.X[i, j] for j in Enrtl_sets.Sc) == m.X[i, c]

    model.eq_Yc = Constraint(xa_id, Enrtl_sets.Sc,
                             rule=rule_cationic_charge_fraction)

    #==========================================================================
    # NONRANDOMNESS FACTORS
    model.a3 = Param(Enrtl_sets.cam, default=0.2,
                     doc='salt-molecule non-randomness parameter')

    model.a4 = Param(Enrtl_sets.caca, default=0.2,
                     doc='salt-salt non-randomness parameter')

    # add (1-sumYa) or (1-sumYc) so that eNRTL can reduce to NRTL
    model.a2 = Expression(xa_id, Enrtl_sets.E2)

    for j in xa_id:
        for i in Enrtl_sets.E2:
            # mm
            if i[0] in Enrtl_sets.Sm and i[1] in Enrtl_sets.Sm:
                model.a2[j, i] = 0.2
            # cm
            elif (i[0] in Enrtl_sets.Sc and i[1] in Enrtl_sets.Sm):
                model.a2[j, i] = sum(model.Ya[j, a] *
                                     model.a3[i[0], a, i[1]]
                                     for a in Enrtl_sets.Sa) + 1\
                    - sum(model.Ya[j, a] for a in Enrtl_sets.Sa)
            # mc
            elif (i[0] in Enrtl_sets.Sm and i[1] in Enrtl_sets.Sc):
                model.a2[j, i] = sum(model.Ya[j, a] *
                                     model.a3[i[1], a, i[0]]
                                     for a in Enrtl_sets.Sa) + 1\
                    - sum(model.Ya[j, a] for a in Enrtl_sets.Sa)
            # am
            elif (i[0] in Enrtl_sets.Sa and i[1] in Enrtl_sets.Sm):
                model.a2[j, i] = sum(model.Yc[j, c] *
                                     model.a3[c, i[0], i[1]]
                                     for c in Enrtl_sets.Sc) + 1\
                    - sum(model.Yc[j, c] for c in Enrtl_sets.Sc)
            # ma
            elif (i[0] in Enrtl_sets.Sm and i[1] in Enrtl_sets.Sa):
                model.a2[j, i] = sum(model.Yc[j, c] *
                                     model.a3[c, i[1], i[0]]
                                     for c in Enrtl_sets.Sc) + 1\
                    - sum(model.Yc[j, c] for c in Enrtl_sets.Sc)
            # ca
            elif (i[0] in Enrtl_sets.Sc and i[1] in Enrtl_sets.Sa):
                model.a2[j, i] = \
                    sum(model.Yc[j, c] * model.a4[i[0], i[1], c, i[1]]
                        for c in Enrtl_sets.Sc) + 1\
                    - sum(model.Yc[j, c] for c in Enrtl_sets.Sc)
            # ac
            elif (i[0] in Enrtl_sets.Sa and i[1] in Enrtl_sets.Sc):
                model.a2[j, i] = \
                    sum(model.Ya[j, a] * model.a4[i[1], i[0], i[1], a]
                        for a in Enrtl_sets.Sa) + 1\
                    - sum(model.Ya[j, a] for a in Enrtl_sets.Sa)

    #==========================================================================
    # Energy interaction parameters : t3 and (t2 for mm only)

    model.t3 = Expression(Enrtl_sets.cam)
    for i in Enrtl_sets.cam:
        model.t3[i] = model.Ccam[i] + model.Dcam[i] / model.T

    # t2
    model.t2 = Expression(xa_id, Enrtl_sets.E2)

    # t2 : mm
    for j in xa_id:
        for i in Enrtl_sets.E2:
            if i[0] in Enrtl_sets.Sm and i[1] in Enrtl_sets.Sm:
                model.t2[j, i] = model.Amm[i] + model.Bmm[i] / model.T

    #==========================================================================
    # Boltzmann kind factors
    # G3
    model.G3 = Expression(Enrtl_sets.cam)
    for i in Enrtl_sets.cam:
        model.G3[i] = exp(-model.a3[i] * model.t3[i])

    # G4
    model.G4 = Expression(Enrtl_sets.caca)
    for i in Enrtl_sets.caca:
        model.G4[i] = exp(-model.a4[i] * model.t4[i])
    # G2 : # add (1-sumYa) or (1-sumYc) so that eNRTL can reduce to NRTL
    model.G2 = Expression(xa_id, Enrtl_sets.E2)
    for j in xa_id:
        for i in Enrtl_sets.E2:
            # mm
            if i[0] in Enrtl_sets.Sm and i[1] in Enrtl_sets.Sm:
                model.G2[j, i] = exp(-model.a2[j, i] * model.t2[j, i])
            # cm
            elif (i[0] in Enrtl_sets.Sc and i[1] in Enrtl_sets.Sm):
                model.G2[j, i] = sum(model.Ya[j, a] *
                                           model.G3[i[0], a, i[1]]
                                           for a in Enrtl_sets.Sa) + 1\
                    - sum(model.Ya[j, a] for a in Enrtl_sets.Sa)
            # mc
            elif (i[0] in Enrtl_sets.Sm and i[1] in Enrtl_sets.Sc):
                model.G2[j, i] = sum(model.Ya[j, a] *
                                           model.G3[i[0], i[1], a]
                                           for a in Enrtl_sets.Sa) + 1\
                    - sum(model.Ya[j, a] for a in Enrtl_sets.Sa)
            # am
            elif (i[0] in Enrtl_sets.Sa and i[1] in Enrtl_sets.Sm):
                model.G2[j, i] = sum(model.Yc[j, c] *
                                           model.G3[c, i[0], i[1]]
                                           for c in Enrtl_sets.Sc) + 1\
                    - sum(model.Yc[j, c] for c in Enrtl_sets.Sc)
            # ma
            elif (i[0] in Enrtl_sets.Sm and i[1] in Enrtl_sets.Sa):
                model.G2[j, i] = sum(model.Yc[j, c] *
                                           model.G3[i[0], c, i[1]]
                                           for c in Enrtl_sets.Sc) + 1\
                    - sum(model.Yc[j, c] for c in Enrtl_sets.Sc)
            # ca
            elif (i[0] in Enrtl_sets.Sc and i[1] in Enrtl_sets.Sa):
                model.G2[j, i] = \
                    sum(model.Yc[j, c] *
                        model.G4[i[0], i[1], c, i[1]]
                        for c in Enrtl_sets.Sc) + 1\
                    - sum(model.Yc[j, c] for c in Enrtl_sets.Sc)
            # ac
            elif (i[0] in Enrtl_sets.Sa and i[1] in Enrtl_sets.Sc):
                model.G2[j, i] = \
                    sum(model.Ya[j, a] *
                        model.G4[i[1], i[0], i[1], a]
                        for a in Enrtl_sets.Sa) + 1\
                    - sum(model.Ya[j, a] for a in Enrtl_sets.Sa)

    #======================================================================
    # Energy interaction parameters
    # t2 : cm,am,ca
    for j in xa_id:
        for i in Enrtl_sets.E2:
            if not (i[0] in Enrtl_sets.Sm and i[1] in Enrtl_sets.Sm):
                model.t2[j, i] = -1 * log(model.G2[j, i]) / model.a2[j, i]
    #==========================================================================
    #absolute ionic charge used for PDH model
    #Z = |z| for ions
    #Z = 0 for molecules

    model.Z = Param(Enrtl_sets.Smac,mutable = True,
        doc='absolute value of ionic charge:PDH')
    for i in Enrtl_sets.Sca:
        model.Z[i]=1
    for i in Enrtl_sets.Sm:
        model.Z[i]=0

    #Molecular weight g/mol
    mw = {
        "CO2":44.01,
        "H2O":18.02,
        "N2":28.01,
        "MEA":61.08,
        "O2":32.00}
    model.M = Param(mw.keys(),initialize=mw,doc='molecular weight')

    #mass fraction of solvent on solute free basis
    def rule_ws(m,i,j):
        return m.M[j]*m.xt[i,j]/sum(m.M[k]*m.xt[i,k] for k in Enrtl_sets.Ss)
    model.ws = Expression(xa_id,Enrtl_sets.Ss,rule=rule_ws)

    #mole fraction of solvent on solute free basis
    def rule_xs(m,i,j):
        return m.xt[i,j]/sum(m.xt[i,k] for k in Enrtl_sets.Ss)
    model.xs = Expression(xa_id,Enrtl_sets.Ss,rule=rule_xs)

    #ionic strength
    def rule_Ix(m,i):
        return 0.5*sum(m.xt[i,k]*m.Z[k]**2 for k in Enrtl_sets.Sca)
    model.Ix = Expression(xa_id,rule=rule_Ix)

    #molar volume Parameters
    vm_a = {'MEA':-5.35162e-1,
            'H2O':-3.2484}
    vm_b = {'MEA':-4.51417e2,
            'H2O':1.65e3}
    vm_c = {'MEA':1.19451e6,
            'H2O':7.93e5}

    model.Av = Param(Enrtl_sets.Ss,initialize=vm_a,
                     doc='molar volume Parameter A')
    model.Bv = Param(Enrtl_sets.Ss,initialize=vm_b,
                     doc='molar volume Parameter B')
    model.Cv = Param(Enrtl_sets.Ss,initialize=vm_c,
                     doc='molar volume Parameter C')

    #molar volume of pure solvent m3/mol
    def rule_Vs(m,i):
        return m.M[i]/(m.Av[i]*m.T**2 + m.Bv[i]*m.T + m.Cv[i])
    model.Vs = Expression(Enrtl_sets.Ss,rule=rule_Vs)


    #Dielectric constant parameters
    Es_a = {'MEA':31.07,
            'H2O':78.71}
    Es_b = {'MEA':15128,
            'H2O':31989}
    Es_c = {'MEA':298.15,
            'H2O':298.15}

    model.Ae = Param(Enrtl_sets.Ss,initialize=Es_a,
                     doc='Dielectric constant Parameter A')
    model.Be = Param(Enrtl_sets.Ss,initialize=Es_b,
                     doc='Dielectric constant Parameter B')
    model.Ce = Param(Enrtl_sets.Ss,initialize=Es_c,
                     doc='Dielectric constant Parameter C')

    # Dielectric constant of pure solvent
    def rule_Es(m,i):
        return m.Ae[i] + m.Be[i]*(1/m.T - 1/m.Ce[i])
    model.Es = Expression(Enrtl_sets.Ss,rule=rule_Es)

    def rule_EwT(m):
        return m.Es['H2O']*m.T
    model.EwT = Expression(rule=rule_EwT)

    # Dielectric constant of mixed  solvent

    def rule_Emx(m,i):
        return sum(m.ws[i,k]*m.Es[k] for k in Enrtl_sets.Ss)
    model.Emx = Expression(xa_id, rule=rule_Emx)


    def rule_EsT(m,i):
        return m.Emx[i]*m.T
    model.EsT = Expression(xa_id, rule=rule_EsT)


    # molar density of the mixed solvent mol/m3
    def rule_dens(m,i):
        return 1/sum(m.xs[i,k]*m.Vs[k] for k in Enrtl_sets.Ss)
    model.dens = Expression(xa_id, rule=rule_dens)

    #Debye-Huckel parameter on mole fraction basis
    def rule_Ax(m,i):
        NA =PhyCons.NA
        kb=PhyCons.kb
        e_charge=PhyCons.e_charge
        eo=PhyCons.eo
        return 1/3*(2*3.142*NA*m.dens[i])**0.5 * \
               (e_charge**2/(4*3.142*eo*kb*m.EsT[i]))**(3/2)
    model.Ax = Expression(xa_id, rule=rule_Ax)

    #==========================================================================
    # infinite dilution activity coefficients of solutes(ionic & molecular)
    # based on aqueous reference state
    def rule_ln_ac_inf_aq(m,i,j):
        if j in Enrtl_sets.Sh:
            return m.t2[i,'H2O',j] + m.G2[i,j,'H2O']*m.t2[i,j,'H2O']

        elif j in Enrtl_sets.Sc:
            return m.z[j]*(m.t2[i,'H2O',j] + m.G2[i,j,'H2O']*m.t2[i,j,'H2O'])

        elif j in Enrtl_sets.Sa:
            return m.z[j]*(m.t2[i,'H2O',j] + m.G2[i,j,'H2O']*m.t2[i,j,'H2O'])

    model.ln_ac_inf_aq = Expression(xa_id,Enrtl_sets.solutes,
                                    rule=rule_ln_ac_inf_aq)

    # infinite dilution activity coefficients of henry component in MEA
    def rule_ln_ac_inf_mea(m,i,j):
        return m.t2[i,'MEA',j] + m.G2[i,j,'MEA']*m.t2[i,j,'MEA']

    model.ln_ac_inf_mea = Expression(xa_id,Enrtl_sets.Sh,
                                     rule=rule_ln_ac_inf_mea)

    #infinite dilution activity coefficient of solutes based on MIXED
    #solvent solution(MEA & H2O)
    def rule_ln_ac_inf_mx(b,i,m):
        ss = Enrtl_sets.Ss
        xs = b.xs
        t2 = b.t2
        G2 = b.G2
        return sum(xs[i,s1]*G2[i,s1,m]*t2[i,s1,m] for s1 in ss)/\
               sum(xs[i,s2]*G2[i,s2,m] for s2 in ss)+\
               sum(xs[i,s]*G2[i,m,s]/sum(xs[i,s3]*G2[i,s3,s] for s3 in ss)
                *(t2[i,m,s] - sum(xs[i,s4]*G2[i,s4,s]*t2[i,s4,s] for s4 in ss)/sum(xs[i,s5]*G2[i,s5,s] for s5 in ss)) for s in ss)
    model.ln_ac_inf_mx = Expression(xa_id,Enrtl_sets.Sh,
                                    rule=rule_ln_ac_inf_mx)


    def rule_ln_ac_SR(b,i,j):
        z = b.z
        X = b.X
        G2 = b.G2
        t2 = b.t2
        mac = Enrtl_sets.Smac
        ma = Enrtl_sets.Sma
        mc = Enrtl_sets.Smc
        sm = Enrtl_sets.Sm
        sc = Enrtl_sets.Sc
        sa = Enrtl_sets.Sa

        if j in sm:
            SR = sum(X[i,k]*G2[i,k,j]*t2[i,k,j] for k in mac)/\
                 sum(X[i,r]*G2[i,r,j] for r in mac) +\
                 sum(X[i,m]*G2[i,j,m]/sum(X[i,r1]*G2[i,r1,m] for r1 in mac)*
                    (t2[i,j,m] -sum(X[i,r2]*G2[i,r2,m]*t2[i,r2,m] for r2 in mac)/sum(X[i,r3]*G2[i,r3,m] for r3 in mac)) for m in sm)+\
                 sum(X[i,c]*G2[i,j,c]/sum(X[i,p1]*G2[i,p1,c] for p1 in ma)*
                    (t2[i,j,c] - sum(X[i,p2]*G2[i,p2,c]*t2[i,p2,c]
                                     for p2 in ma)/
                                 sum(X[i,p3]*G2[i,p3,c] for p3 in ma))
                     for c in sc) +\
                 sum(X[i,a]*G2[i,j,a]/sum(X[i,k1]*G2[i,k1,a] for k1 in mc)*
                    (t2[i,j,a] - sum(X[i,k2]*G2[i,k2,a]*t2[i,k2,a] for k2 in mc)/sum(X[i,k3]*G2[i,k3,a] for k3 in mc)) for a in sa)
        elif j in sc:
            SR = z[j]*(sum(X[i,k]*G2[i,k,j]*t2[i,k,j] for k in ma)/
                 sum(X[i,r]*G2[i,r,j] for r in ma) +
                 sum(X[i,m]*G2[i,j,m]/sum(X[i,r1]*G2[i,r1,m] for r1 in mac) *
                    (t2[i,j,m] - sum(X[i,r2]*G2[i,r2,m]*t2[i,r2,m] for r2 in mac)/sum(X[i,r3]*G2[i,r3,m] for r3 in mac)) for m in sm) +
                 sum(X[i,a]*G2[i,j,a]/sum(X[i,k1]*G2[i,k1,a] for k1 in mc)*
                    (t2[i,j,a] - sum(X[i,k2]*G2[i,k2,a]*t2[i,k2,a] for k2 in mc)/sum(X[i,k3]*G2[i,k3,a] for k3 in mc)) for a in sa))
        elif j in sa:
            SR = z[j]*(sum(X[i,k]*G2[i,k,j]*t2[i,k,j] for k in mc)/
                 sum(X[i,r]*G2[i,r,j] for r in mc) +
                 sum(X[i,m]*G2[i,j,m]/sum(X[i,r1]*G2[i,r1,m] for r1 in mac) *
                    (t2[i,j,m] - sum(X[i,r2]*G2[i,r2,m]*t2[i,r2,m] for r2 in mac)/sum(X[i,r3]*G2[i,r3,m] for r3 in mac)) for m in sm) +
                 sum(X[i,c]*G2[i,j,c]/sum(X[i,k1]*G2[i,k1,c] for k1 in ma)*
                    (t2[i,j,c] - sum(X[i,k2]*G2[i,k2,c]*t2[i,k2,c] for k2 in ma)/sum(X[i,k3]*G2[i,k3,c] for k3 in ma)) for c in sc))

        return SR

    model.ln_ac_SR = Expression(xa_id,Enrtl_sets.Smac,
                                 rule=rule_ln_ac_SR)

    #Born Term contribution
    def rule_ln_ac_born(b,i,j):
        kb=PhyCons.kb
        e_charge=PhyCons.e_charge
        eo=PhyCons.eo
        rB = PhyCons.rB
        Z = b.Z
        EsT = b.EsT
        EwT = b.EwT
        return (e_charge*Z[j])**2*(1/EsT[i]-1/EwT)/(800*3.142*eo*rB*kb)
    model.ln_ac_born = Expression(xa_id, Enrtl_sets.Smac,
                       rule=rule_ln_ac_born)

    #Long range(LR) contribution to natural log of activity coef PDH model
    def rule_PDH(b,i,j):
        Ax = b.Ax
        Ix = b.Ix
        ee = 14.9 # approach parameter
        Z = b.Z
        PDH = -Ax[i]*(2*Z[j]**2/ee*log(1+ee*Ix[i]**0.5) +
                        (Z[j]**2*Ix[i]**0.5 - 2*Ix[i]**1.5)/(1+ee*Ix[i]**0.5))
        return PDH
    model.ln_ac_pdh = Expression(xa_id,Enrtl_sets.Smac,
                    rule=rule_PDH)

    #  un-symmetric SR activity coefficent based on aqueous phase
    # reference state
    def rule_ln_ac_SR_aq(b,i,j):
        return b.ln_ac_SR[i,j] - b.ln_ac_inf_aq[i,j]
    model.ln_ac_SR_aq = Expression(xa_id,Enrtl_sets.solutes,
                                rule=rule_ln_ac_SR_aq)

    #un-symmetric SR activity coefficent of henry components based on mixed solvent
    #reference state
    def rule_ln_ac_SR_mx(b,i,j):
        return b.ln_ac_SR[i,j] - b.ln_ac_inf_mx[i,j]

    model.ln_ac_SR_mx = Expression(xa_id,Enrtl_sets.Sh,
                                   rule=rule_ln_ac_SR_mx)

    #natural log of activity coefficent of   solvent species
    def rule_ln_ac_solvent(b,i,j):
        return b.ln_ac_SR[i,j] + b. ln_ac_pdh[i,j]
    model.ln_ac_solvent = Expression(xa_id, Enrtl_sets.Ss,
                                     rule=rule_ln_ac_solvent)

    #activity coefficent of  solutes species based on aqueous reference state
    def rule_ln_ac_solute(b,i,j):
        return b.ln_ac_SR_aq[i,j] + b. ln_ac_pdh[i,j] + b.ln_ac_born[i,j]
    model.ln_ac_solute = Expression(xa_id, Enrtl_sets.solutes,
                                     rule=rule_ln_ac_solute)

    # activity coef. of  henry component based on mixed solvent reference state
    def rule_ln_ac_solute_mx(b,i,j):
        return b.ln_ac_SR_mx[i,j] + b. ln_ac_pdh[i,j]
    model.ln_ac_solute_mx = Expression(xa_id, Enrtl_sets.Sh,
                                     rule=rule_ln_ac_solute_mx)

    #==========================================================================
    # Excess enthalpy calculation

    #Born term
    def rule_GB1(b,j):
        kb = PhyCons.kb
        e_charge = PhyCons.e_charge
        eo = PhyCons.eo
        rB = PhyCons.rB
        return e_charge**2 / (800*3.142*eo*rB*kb) *\
               sum(b.xt[j,i]*b.Z[i]**2 for i in Enrtl_sets.Smac)
    model.GB1 = Expression(xa_id, rule=rule_GB1)

    model.GB2 = Expression(expr = (model.Ae['H2O']-1/model.Ce['H2O'])/(model.EwT)**2)

    def rule_GB3(b,j):
        return sum(b.ws[j,i]*(b.Ae[i] -1/b.Ce[i]) for i in Enrtl_sets.Ss)/\
              (b.EsT[j])**2
    model.GB3 = Expression(xa_id, rule=rule_GB3)

    def rule_dGB(b,j):
        return b.GB1[j]*(b.GB2-b.GB3[j])
    model.dGB = Expression(xa_id, rule=rule_dGB)

    # dLR: Hex contribution from Long Range(LR)
    def rule_dLR1(b,j):
        return -4*b.Ix[j]/14.9*log(1+14.9*b.Ix[j]**0.5)
    model.dLR1 = Expression(xa_id, rule=rule_dLR1)

    def rule_dLR2(b):
        NA = PhyCons.NA
        kb = PhyCons.kb
        e_charge = PhyCons.e_charge
        eo = PhyCons.eo
        return 1/3*(2*3.142*NA)**0.5*(e_charge**2/(4*3.142*eo*kb))**(3/2)
    model.dLR2 = Expression(rule=rule_dLR2)

    def rule_dEstdT(b,j):
        return sum(b.ws[j,i]*(b.Ae[i] -1/b.Ce[i]) for i in Enrtl_sets.Ss)
    model.dEstdT = Expression(xa_id, rule=rule_dEstdT)

    def rule_ddensdT(b,j):
        return (b.dens[j])**2*\
            sum(b.xs[j,i]*b.Vs[i]*(2*b.Av[i]*b.T + b.Bv[i])/
            (b.Av[i]*b.T**2 + b.Bv[i]*b.T +b.Cv[i]) for i in Enrtl_sets.Ss)
    model.ddensdT = Expression(xa_id, rule=rule_ddensdT)

    def rule_dLR(b,j):
        return b.dLR1[j]*b.dLR2*(
               (b.dens[j])**0.5*(-1.5)*b.EsT[j]**(-5/2)*b.dEstdT[j] +
                b.EsT[j]**(-1.5)*0.5*(b.dens[j])**(-0.5)*b.ddensdT[j])
    model.dLR = Expression(xa_id, rule=rule_dLR)

    # dSR: Hex contribution from Short Range(SR)
    #dtijdT
    model.dtimdT = Expression(xa_id,Enrtl_sets.Smac,Enrtl_sets.Sm)
    for p in xa_id:
        for j in Enrtl_sets.Sm:
            for i in Enrtl_sets.Smac:
                if i in Enrtl_sets.Sm:
                    model.dtimdT[p,i,j] = -model.Bmm[i,j]/model.T**2
                elif i in Enrtl_sets.Sc:
                    model.dtimdT[p,i,j] = -sum(model.Ya[p,a]*model.Dcam[i,a,j]*
                        model.a3[i,a,j]*model.G3[i,a,j] for a in Enrtl_sets.Sa)/(model.a2[p,i,j]*model.G2[p,i,j]*model.T**2)
                elif i in Enrtl_sets.Sa:
                    model.dtimdT[p,i,j] = -sum(model.Yc[p,c]*model.Dcam[c,i,j]*
                        model.a3[c,i,j]*model.G3[c,i,j] for c in Enrtl_sets.Sc)/(model.a2[p,i,j]*model.G2[p,i,j]*model.T**2)

    model.dticdT = Expression(xa_id,Enrtl_sets.Sma,Enrtl_sets.Sc)
    for p in xa_id:
        for j in Enrtl_sets.Sc:
            for i in Enrtl_sets.Sma:
                if i in Enrtl_sets.Sm:
                    model.dticdT[p,i,j] = -sum(model.Ya[p,a]*model.Dcam[i,j,a]*model.a3[i,j,a]*model.G3[i,j,a] for a in Enrtl_sets.Sa)/(model.a2[p,i,j]*model.G2[p,i,j]*model.T**2)
                elif i in Enrtl_sets.Sa:
                    model.dticdT[p,i,j] = 0

    model.dtiadT = Expression(xa_id,Enrtl_sets.Smc,Enrtl_sets.Sa)
    for p in xa_id:
        for j in Enrtl_sets.Sa:
            for i in Enrtl_sets.Smc:
                if i in Enrtl_sets.Sm:
                    model.dtiadT[p,i,j] = -sum(model.Yc[p,c]*model.Dcam[i,c,j]*model.a3[i,c,j]*model.G3[i,c,j] for c in Enrtl_sets.Sc)/(model.a2[p,i,j]*model.G2[p,i,j]*model.T**2)
                elif i in Enrtl_sets.Sc:
                    model.dtiadT[p,i,j] = 0

    # dGijdT
    model.dGimdT = Expression(xa_id,Enrtl_sets.Smac,Enrtl_sets.Sm)
    for p in xa_id:
        for j in Enrtl_sets.Sm:
            for i in Enrtl_sets.Smac:
                if i in Enrtl_sets.Sm:
                    model.dGimdT[p,i,j] = model.Bmm[i,j]*model.a2[p,i,j]*model.G2[p,i,j]/model.T**2
                elif i in Enrtl_sets.Sc:
                    model.dGimdT[p,i,j] = sum(model.Ya[p,a]*model.Dcam[i,a,j]*model.a3[i,a,j]*
                                       model.G3[i,a,j] for a in Enrtl_sets.Sa)/model.T**2
                elif i in Enrtl_sets.Sa:
                    model.dGimdT[p,i,j] = sum(model.Yc[p,c]*model.Dcam[c,i,j]*model.a3[c,i,j]*
                                       model.G3[c,i,j] for c in Enrtl_sets.Sc)/model.T**2

    model.dGicdT = Expression(xa_id,Enrtl_sets.Sma,Enrtl_sets.Sc)
    for p in xa_id:
        for j in Enrtl_sets.Sc:
            for i in Enrtl_sets.Sma:
                if i in Enrtl_sets.Sm:
                    model.dGicdT[p,i,j] = sum(model.Ya[p,a]*model.Dcam[i,j,a]*model.a3[i,j,a]*
                                       model.G3[i,j,a] for a in Enrtl_sets.Sa)/model.T**2
                elif i in Enrtl_sets.Sa:
                    model.dGicdT[p,i,j] = 0


    model.dGiadT = Expression(xa_id,Enrtl_sets.Smc,Enrtl_sets.Sa)
    for p in xa_id:
        for j in Enrtl_sets.Sa:
            for i in Enrtl_sets.Smc:
                if i in Enrtl_sets.Sm:
                    model.dGiadT[p,i,j] = sum(model.Yc[p,c]*model.Dcam[i,c,j]*model.a3[i,c,j]*
                                       model.G3[i,c,j] for c in Enrtl_sets.Sc)/model.T**2
                elif i in Enrtl_sets.Sc:
                    model.dGiadT[p,i,j] = 0

        # derivative of infinite dilution activity
    model.dlncdT = Expression(xa_id,Enrtl_sets.Sc)
    model.dlnadT = Expression(xa_id,Enrtl_sets.Sa)
    model.dlnhdT = Expression(xa_id,Enrtl_sets.Sh)

    for p in xa_id:
        for i in Enrtl_sets.Sc:
            model.dlncdT[p,i] = model.z[i]*(model.dticdT[p,'H2O',i] +
                        model.G2[p,i,'H2O']*model.dtimdT[p,i,'H2O'] +
                        model.t2[p,i,'H2O']*model.dGimdT[p,i,'H2O'])

        for i in Enrtl_sets.Sa:
            model.dlnadT[p,i] = model.z[i]*(model.dtiadT[p,'H2O',i] +
                        model.G2[p,i,'H2O']*model.dtimdT[p,i,'H2O'] +
                        model.t2[p,i,'H2O']*model.dGimdT[p,i,'H2O'])

        for i in Enrtl_sets.Sh:
            model.dlnhdT[p,i] = (model.dtimdT[p,'H2O',i] +
                        model.G2[p,i,'H2O']*model.dtimdT[p,i,'H2O'] +
                        model.t2[p,i,'H2O']*model.dGimdT[p,i,'H2O'])

    # dGSR/dT = SRm + SRc + SRa
    def rule_dSRm(b,p):
        sm = Enrtl_sets.Sm
        mac = Enrtl_sets.Smac

        return sum(b.X[p,m]*(sum(b.X[p,i]*(b.G2[p,i,m]*b.dtimdT[p,i,m]+
                                 b.t2[p,i,m]*b.dGimdT[p,i,m]) for i in mac)/
                             sum(b.X[p,j]*b.G2[p,j,m] for j in mac)
                             -
                             sum(b.X[p,k]*b.G2[p,k,m]*b.t2[p,k,m] for k in mac)*
                             sum(b.X[p,h]*b.dGimdT[p,h,m] for h in mac)/
                             sum(b.X[p,r]*b.G2[p,r,m] for r in mac)**2)
                   for m in sm)
    model.dSRm = Expression(xa_id, rule=rule_dSRm)

    def rule_dSRc(b,p):
        sc = Enrtl_sets.Sc
        ma = Enrtl_sets.Sma
        return sum(b.X[p,c]*(sum(b.X[p,i]*(b.G2[p,i,c]*b.dticdT[p,i,c]+
                                 b.t2[p,i,c]*b.dGicdT[p,i,c]) for i in ma)/
                             sum(b.X[p,j]*b.G2[p,j,c] for j in ma)
                             -
                             sum(b.X[p,k]*b.G2[p,k,c]*b.t2[p,k,c] for k in ma)*
                             sum(b.X[p,h]*b.dGicdT[p,h,c] for h in ma)/
                             sum(b.X[p,r]*b.G2[p,r,c] for r in ma)**2)
                   for c in sc)

    model.dSRc = Expression(xa_id, rule=rule_dSRc)

    def rule_dSRa(b,p):
        sa = Enrtl_sets.Sa
        mc = Enrtl_sets.Smc
        return sum(b.X[p,a]*(sum(b.X[p,i]*(b.G2[p,i,a]*b.dtiadT[p,i,a]+
                                 b.t2[p,i,a]*b.dGiadT[p,i,a]) for i in mc)/
                             sum(b.X[p,j]*b.G2[p,j,a] for j in mc)
                             -
                             sum(b.X[p,k]*b.G2[p,k,a]*b.t2[p,k,a] for k in mc)*
                             sum(b.X[p,h]*b.dGiadT[p,h,a] for h in mc)/
                             sum(b.X[p,r]*b.G2[p,r,a] for r in mc)**2)
                   for a in sa)

    model.dSRa = Expression(xa_id, rule=rule_dSRa)

    def rule_dSR(b,p):
        sc = Enrtl_sets.Sc
        sa = Enrtl_sets.Sa
        sh = Enrtl_sets.Sh
        return b.dSRm[p] + b.dSRc[p] + b.dSRa[p]  \
               - sum(b.xt[p,c]*b.dlncdT[p,c] for c in sc) \
               - sum(b.xt[p,a]*b.dlnadT[p,a] for a in sa) \
               - sum(b.xt[p,h]*b.dlnhdT[p,h] for h in sh)

    model.dSR = Expression(xa_id, rule=rule_dSR)

    def rule_Hex(b,p):
        R = PhyCons.R_gas
        return -R*b.T**2*(b.dSR[p] + b.dLR[p] + b.dGB[p])
    model.Hex = Expression(xa_id, rule=rule_Hex)

    #==========================================================================
    #GFORM_ig (kJ/mol)
    gform_ig = {
    'MEA':-103.300,
    'H2O':-228.590,
    'CO2':-394.370
    }
    model.GFORM_ig = Param(Enrtl_sets.Sm,initialize=gform_ig)

    #HFORM_ig (kJ/mol)
    hform_ig = {
    'MEA':-206.700,
    'H2O':-241.820,
    'CO2':-393.510
    }
    model.HFORM_ig = Param(Enrtl_sets.Sm,initialize=hform_ig)


    # Cp_ig (kJ/mol K)
    cp_param = {
        (1,'MEA'):1.3207E1*1e-3,
        (1,'H2O'):32.22*1e-3,
        (1,'CO2'):1.9795E1*1e-3,
        (2,'MEA'):2.8158E-1*1e-3,
        (2,'H2O'):1.923E-3*1e-3,
        (2,'CO2'):7.3437E-2*1e-3,
        (3,'MEA'):-0.1513E-3*1e-3,
        (3,'H2O'):10.548E-6*1e-3,
        (3,'CO2'):-5.6019E-5*1e-3,
        (4,'MEA'):3.1287E-8*1e-3,
        (4,'H2O'):-3.594E-9*1e-3,
        (4,'CO2'):1.7153E-8*1e-3
        }
    model.Cp_ig = Param([1,2,3,4],Enrtl_sets.Sm,initialize=cp_param)

    # HCO3- anion
    model.GFORM_HCO3 = Param(initialize=-576.81367)
    model.HFORM_HCO3 = Param(initialize=-691.990)
    model.Cp_HCO3 = Param(initialize=-29.26*1e-3)

    #GFORM_aq (kJ/mol): MEACOO- & MEAH+
    gform_aq = {
        'MEACOO-':-481.9013924509483,#-492.99,
        'MEAH+':-179.33957077685685#-189.62
        }
    model.GFORM_aq = Var(['MEACOO-','MEAH+'],initialize=gform_aq)
    model.GFORM_aq.fix()

    #HFORM_aq (kJ/mol): MEACOO- & MEAH+
    hform_aq = {
        'MEACOO-':-711.9870815027019,#-707.47,
        'MEAH+':-324.00385529510794#-331.64
        }
    model.HFORM_aq = Var(['MEACOO-','MEAH+'],initialize=hform_aq)
    model.HFORM_aq.fix()

    #Cp_aq (J/mol): MEACOO- & MEAH+
    #A
    cp_aq_A = {
        'MEACOO-':0.11354500314776703,
        'MEAH+': 0.2584025436154042
        }
    model.Acp_aq = Var(['MEACOO-','MEAH+'],initialize=cp_aq_A,bounds=(-1,1))
    model.Acp_aq.fix()
    # B
    cp_aq_B = {
        'MEACOO-':0,#4.23e-2*1e-3,
        'MEAH+':0#-1.63*1e-3
        }
    model.Bcp_aq = Var(['MEACOO-','MEAH+'],initialize=cp_aq_B,
                    bounds=(-1,1),within=Reals)
    model.Bcp_aq.fix()

    #==========================================================================
    #Enthalpy model
    R = PhyCons.R_gas
    # critical pressure in kPa
    Pc_param = {'H2O':22064,'MEA':7124 ,'CO2':7383 }
    model.Pc = Param(Enrtl_sets.Sm,initialize=Pc_param)

    # critical temperature in K
    Tc_param = {'H2O':647.096 ,'MEA':678.2 ,'CO2':304.2 }
    model.Tc = Param(Enrtl_sets.Sm,initialize=Tc_param)

    # omega
    omega_param = {'H2O':0.344 ,'MEA':0.446737 ,'CO2':0.224 }
    model.omega = Param(Enrtl_sets.Sm,initialize=omega_param)

    # Heat of vaporization
    hvap_param = {
        ('MEA','a') :0.3288,
        ('MEA','b') :-0.0857,
        ('MEA','t1'):126.67,
        ('MEA','tc'):398.25,
        ('MEA','dH'):54835.8*1e-3,
        ('H2O','a') :0.3106,
        ('H2O','b') :0.0,
        ('H2O','t1'):100,
        ('H2O','tc'):373.946,
        ('H2O','dH'):40655*1e-3
        }
    def rule_hvap(b,i):
        return hvap_param[i,'dH']*((1- (b.T-273.15)/hvap_param[i,'tc'])/
                        (1- hvap_param[i,'t1']/hvap_param[i,'tc']))**(hvap_param[i,'a']+ hvap_param[i,'b']*(1-(b.T-273.15)/hvap_param[i,'tc']))
    model.Hvap = Expression(Enrtl_sets.Ss,rule=rule_hvap)

    psat_par = {
    (1,'MEA'):165.87,
    (1,'H2O'):73.649+ log(1e-3),
    (2,'MEA'):-13492,
    (2,'H2O'):-7258.2,
    (3,'MEA'):-21.9,
    (3,'H2O'):-7.3037,
    (4,'MEA'):1.38e-5,
    (4,'H2O'):4.1653e-6
    }

    # saturation pressure of solvents
    def rule_lnpsat(b,i):
        return psat_par[1,i] + psat_par[2,i]/b.T +\
               psat_par[3,i]*log(b.T) + psat_par[4,i]*b.T**2

    model.lnPsat = Expression(Enrtl_sets.Ss,rule=rule_lnpsat,
                    doc='natural log of solvent saturation pressure')

    def rule_psat(b,i):
        return exp(b.lnPsat[i])
    model.Psat = Expression(Enrtl_sets.Ss,rule=rule_psat,
                    doc='solvent saturation pressure in kPa')

    # second viral cofficient
    def rule_bii(b,i):
        return R*b.Tc[i]/b.Pc[i]*\
               (0.083-0.422/(b.T/b.Tc[i])**1.6 +
                b.omega[i]*(0.139-0.172/(b.T/b.Tc[i])**4.2))
    model.Bii = Expression(Enrtl_sets.Sm,rule=rule_bii)

    def rule_dbiidT(b,i):
        return R/b.Pc[i]*(0.6752/(b.T/b.Tc[i])*2.6 + b.omega[i]*0.7224/(b.T/b.Tc[i])**5.2)
    model.dBiidT = Expression(Enrtl_sets.Sm,rule=rule_dbiidT)

    #Departure terms for molecular solutes : MEA,H2O, CO2
    def rule_dep(b,i):
        if i in Enrtl_sets.Ss:
            return b.Psat[i]*(b.Bii[i]-b.T*b.dBiidT[i])
        elif i == 'CO2':
            return b.P*(b.Bii[i]-b.T*b.dBiidT[i])
    model.dep = Expression(Enrtl_sets.Sm,rule=rule_dep)

    # Vapor Enthalpy of CO2
    def rule_HvCO2(b):
        To = 298.15
        HvCO2 = b.HFORM_ig['CO2'] + b.dep['CO2'] +\
                b.Cp_ig[1,'CO2']*(b.T-To) + \
                b.Cp_ig[2,'CO2']*(b.T**2 - To**2)/2 + \
                b.Cp_ig[3,'CO2']*(b.T**3 - To**3)/3 + \
                b.Cp_ig[4,'CO2']*(b.T**4 - To**4)/4
        return HvCO2
    model.Hv_CO2 = Expression(rule=rule_HvCO2)

    #Pure component liquid enthalpy
    def rule_HLi(b,i):
        To = 298.15
        if i == 'CO2':
            # henry parameters for CO2 in H2O
            #cis = [-8477.7,-21.957, 0.00578]
            cis = [-5876,-8.598, -0.012]
            return b.Hv_CO2 - R*(-cis[0]+cis[1]*b.T+cis[2]*b.T**2)
        elif i in Enrtl_sets.Ss:
            return b.HFORM_ig[i] + b.dep[i] - b.Hvap[i] + \
                      b.Cp_ig[1,i]*(b.T-To) + \
                      b.Cp_ig[2,i]*(b.T**2 -To**2)/2 + \
                      b.Cp_ig[3,i]*(b.T**3 -To**3)/3 + \
                      b.Cp_ig[4,i]*(b.T**4 -To**4)/4
        elif i in ['MEACOO-','MEAH+']:
            return b.HFORM_aq[i] + b.Acp_aq[i]*(b.T-To) + \
                      b.Bcp_aq[i]*(b.T**2-To**2)/2
        elif i == 'HCO3-':
            return b.HFORM_HCO3 + b.Cp_HCO3*(b.T-To)

    model.HLi = Expression(Enrtl_sets.Smac,rule=rule_HLi)

    # # Mixture Enthalpy(Liquid in kJ/mol)
    # def rule_HL(b,p):
    #     return  sum(b.xt[p,i]*b.HLi[i] for i in Enrtl_sets.Smac) + b.Hex[p]
    # model.HL = Expression(xa_id,rule=rule_HL)


    #=================================================================
    # projector term to convert true Hab to Habs_apparent
    def rule_dG_TdT(b,i):
        To = 298.15
        if i in Enrtl_sets.Ss:
            dG_TdT  = (R*(-psat_par[2,i]/b.T**2 + psat_par[3,i]/b.T +
                           2*psat_par[4,i]*b.T +
                           b.Bii[i]/(R*b.T)*(-psat_par[2,i]/b.T**2 + psat_par[3,i]/b.T + 2*psat_par[4,i]*b.T)*b.Psat[i] +
                           b.Psat[i]/(R*b.T)*b.dBiidT[i] -
                           b.Psat[i]*b.Bii[i]/(R*b.T**2) +
                           1/R*(-b.Vs[i]/b.T*b.Psat[i]*(-psat_par[2,i]/b.T**2 + psat_par[3,i]/b.T + 2*psat_par[4,i]*b.T)-
                           (b.P-b.Psat[i])*b.Vs[i]*
                           (3*b.Av[i]*b.T+2*b.Bv[i] + b.Cv[i]/b.T)/
                           (b.Av[i]*b.T**3+ b.Bv[i]*b.T**2 + b.Cv[i]*b.T))
                           )
                       -b.HFORM_ig[i]/b.T**2 +
                       b.Cp_ig[1,i]*(To/b.T**2 - 1/b.T) -
                       b.Cp_ig[2,i]/2*(1 - (To/b.T)**2) -
                       b.Cp_ig[3,i]/3*(b.T - To**3/b.T**2)-
                       b.Cp_ig[4,i]/4*(b.T**2 -To**4/b.T**2))

        elif i == 'CO2':
            # henry parameters for CO2 in H2O
            his = [91.344,-5876,-8.598, -0.012]
            dG_TdT = (R*(-his[1]/b.T**2+ his[2]/b.T+ his[3])
                      -b.HFORM_ig[i]/b.T**2 +
                      b.Cp_ig[1,i]*(To/b.T**2 - 1/b.T) -
                      b.Cp_ig[2,i]/2*(1 - (To/b.T)**2) -
                      b.Cp_ig[3,i]/3*(b.T - To**3/b.T**2)-
                      b.Cp_ig[4,i]/4*(b.T**2 -To**4/b.T**2))
        elif i in ['MEACOO-','MEAH+']:
            dG_TdT = -b.HFORM_aq[i]/b.T**2 + b.Acp_aq[i]*(To/b.T**2 - 1/b.T)
        elif i == 'HCO3-':
            dG_TdT = -b.HFORM_HCO3/b.T**2 + b.Cp_HCO3*(To/b.T**2 - 1/b.T)

        return dG_TdT
    model.dG_TdT = Expression(Enrtl_sets.Smac, rule=rule_dG_TdT)

    #=================================================================
    def rule_lnqsat(b,i):
        return b.Bii[i]*b.Psat[i]/(R*b.T)
    model.lnQsat = Expression(Enrtl_sets.Ss,rule=rule_lnqsat,
                              doc='natural log of saturation fugacity coefficent')
    def rule_lnfug(b,i):
        return b.lnQsat[i] + b.lnPsat[i] + \
               b.Vs[i]*(b.P - b.Psat[i])/(R*b.T)
    model.lnfug = Expression(Enrtl_sets.Ss,rule=rule_lnfug,
                             doc='natural log of liquid fugacity')

    # Gibbs free energy
    def rule_GFE(b,i):
        To = 298.15
        Pref = 101.325
        if i in Enrtl_sets.Ss:
            GFE  = (R*b.T*(b.lnfug[i]-log(Pref)) +
                   b.GFORM_ig[i]*(b.T/To) +
                   b.HFORM_ig[i]*(1-(b.T/To))+
                   b.Cp_ig[1,i]*To*((b.T/To)-1-(b.T/To)*log((b.T/To))) -
                   b.Cp_ig[2,i]/2*To**2*((b.T/To)-1)**2 -
                   b.Cp_ig[3,i]/6*To**3*((b.T/To)-1)**2*((b.T/To)+2)-
                   b.Cp_ig[4,i]/12*To**4*((b.T/To)-1)**2*((b.T/To)**2+2*(b.T/To)+3)
                   )
        elif i == 'CO2':
            # henry parameters for CO2 in H2O
            #his = [170.713,-8477.7,-21.957, 0.00578]
            his = [91.344,-5876,-8.598, -0.012]
            GFE  = (R*b.T*(his[0]+ his[1]/b.T+ his[2]*log(b.T)+ his[3]*b.T
                  -log(Pref) + log(1e-3)) +
                   b.GFORM_ig[i]*(b.T/To) +
                   b.HFORM_ig[i]*(1-(b.T/To))+
                   b.Cp_ig[1,i]*To*((b.T/To)-1-(b.T/To)*log((b.T/To))) -
                   b.Cp_ig[2,i]/2*To**2*((b.T/To)-1)**2 -
                   b.Cp_ig[3,i]/6*To**3*((b.T/To)-1)**2*((b.T/To)+2)-
                   b.Cp_ig[4,i]/12*To**4*((b.T/To)-1)**2*((b.T/To)**2+2*(b.T/To)+3)
                   )
        elif i in ['MEACOO-','MEAH+']:
            GFE = b.GFORM_aq[i]*(b.T/To) + \
                  b.HFORM_aq[i]*(1-(b.T/To)) +\
                  b.Acp_aq[i]*To*((b.T/To)-1-(b.T/To)*log((b.T/To))) - \
                  b.Bcp_aq[i]/2*To**2*((b.T/To)-1)**2
        elif i == 'HCO3-':
            GFE = b.GFORM_HCO3*(b.T/To) + \
                  b.HFORM_HCO3*(1-(b.T/To)) +\
                  b.Cp_HCO3*To*((b.T/To)-1-(b.T/To)*log((b.T/To)))
        return GFE
    model.GFE = Expression(Enrtl_sets.Smac, rule=rule_GFE)

    #=================================================================
    #equilibruim constant
    stoich = dict()
    stoich[1,'CO2'] = -1
    stoich[1,'MEA'] = -2
    stoich[1,'H2O'] = 0
    stoich[1,'MEAH+'] = 1
    stoich[1,'HCO3-'] = 0
    stoich[1,'MEACOO-'] = 1
    stoich[2,'CO2'] = -1
    stoich[2,'MEA'] = -1
    stoich[2,'H2O'] = -1
    stoich[2,'MEAH+'] = 1
    stoich[2,'HCO3-'] = 1
    stoich[2,'MEACOO-'] = 0


    def rule_lnK(b,r):
        return -sum(stoich[r,i]*b.GFE[i] for i in Enrtl_sets.Smac)/(R*b.T)
    model.lnK = Expression([1,2],rule=rule_lnK)

    def rule_dlnKdT(b,r):
        return -b.T**2/R*sum(stoich[r,i]*b.dG_TdT[i] for i in Enrtl_sets.Smac)
    model.dlnKdT = Expression([1,2],rule=rule_dlnKdT)

    def rule_projector(b,p):
        return R*(b.dlnKdT[1]*b.nt[p,'MEACOO-'] +
                  b.dlnKdT[2]*b.nt[p,'HCO3-'])
    model.projector = Expression(xa_id,rule=rule_projector)

    def rule_Hex_app(b,p):
        return (b.Hex[p]*b.nt_tot[p] + b.projector[p])
    model.Hex_app = Expression(xa_id, rule=rule_Hex_app)

    def rule_Hex_app_per_mol(b,p):
        return (b.Hex[p]*b.nt_tot[p] + b.projector[p])/b.nT[p]
    model.Hex_app_per_mol = Expression(xa_id, rule=rule_Hex_app_per_mol)

    # Mixture Enthalpy(Liquid in kJ/mol)
    def rule_HLt(b,p):
        return  sum(b.xt[p,i]*b.HLi[i] for i in Enrtl_sets.Smac)  + b.Hex[p]
    model.HLt = Expression(xa_id,rule=rule_HLt)

    def rule_HL(b,p):
        return  b.n_H2O*b.HLi['H2O'] + b.n_MEA*b.HLi['MEA'] \
              + b.n_CO2[p]*b.HLi['CO2'] + b.Hex_app[p]
    model.HL = Expression(xa_id,rule=rule_HL)

    #molecular weight in kg/mol
    MW = {'MEA':61.08/1000,
          'H2O':18.02/1000,
          'CO2':44.01/1000}
    model.MW = Param(Enrtl_sets.Sm,initialize=MW,
                     doc='Molecular weight')

    def rule_mw_avg(b,p):
        return sum(b.xa[p,i]*b.MW[i] for i in Enrtl_sets.Sm)
    model.mw_avg = Expression(xa_id,rule=rule_mw_avg)

    def rule_mHL(b,p):
        return  (b.HL[p]/b.nT[p])/b.mw_avg[p]
    model.mHL = Expression(xa_id,rule=rule_mHL)

    #==================================================================
    # Speciation model.
    def rule_speciation(b,i,j):
        if i == 1:
            return b.lnK[1] == (log(b.xt[j,'MEACOO-']) +
                                log(b.xt[j,'MEAH+']) -
                                log(b.xt[j,'CO2'])-
                              2*log(b.xt[j,'MEA']) +
                                    b.ln_ac_solute[j,'MEACOO-']+
                                    b.ln_ac_solute[j,'MEAH+'] -
                                    b.ln_ac_solute[j,'CO2'] -
                                  2*b.ln_ac_solvent[j,'MEA'])
        elif i == 2:
            return b.lnK[2] == (log(b.xt[j,'HCO3-']) +
                                log(b.xt[j,'MEAH+']) -
                                log(b.xt[j,'CO2']) -
                                log(b.xt[j,'MEA']) -
                                log(b.xt[j,'H2O']) +
                                    b.ln_ac_solute[j,'MEACOO-']+
                                    b.ln_ac_solute[j,'MEAH+'] -
                                    b.ln_ac_solute[j,'CO2'] -
                                    b.ln_ac_solvent[j,'MEA']-
                                    b.ln_ac_solvent[j,'H2O'])
        elif i ==3:
            return b.n_MEA == \
                   b.nt[j,'MEA'] + b.nt[j,'MEAH+'] + b.nt[j,'MEACOO-']
        elif i ==4:
            return b.n_CO2[j] == \
                   b.nt[j,'CO2'] + b.nt[j,'HCO3-'] + b.nt[j,'MEACOO-']
        elif i ==5:
            return b.n_H2O == b.nt[j,'H2O'] + b.nt[j,'HCO3-']
        elif i ==6:
            return b.nt[j,'MEAH+']== b.nt[j,'HCO3-'] + b.nt[j,'MEACOO-']

    model.speciation = Constraint([1,2,3,4,5,6],xa_id,rule=rule_speciation)

    #================================================================
    #Phase equilibruim
    # interaction paramter for Henry's mixing rule
    phi_set = []
    for i in Enrtl_sets.Ss:
        for j in Enrtl_sets.Ss:
            phi_set.append((i,j))
    model.phi = Expression(phi_set)

    for i in phi_set:
        if i == ('MEA','H2O'):
            model.phi[i]= 8.7105 - 1.324*log(model.T)
        elif i == ('H2O','MEA'):
            model.phi[i]= 1/(8.7105 - 1.324*log(model.T))
        else:
            model.phi[i]= 1

    def rule_wh(b,i,j):
        return b.xs[i,j]/sum(b.xs[i,k]*b.phi[j,k] for k in Enrtl_sets.Ss)
    model.wh = Expression(xa_id,Enrtl_sets.Ss, rule=rule_wh,
                          doc='weighting factorfor Henry mixing rule')

    def rule_lnHis(b,i):
        if i == 'H2O':
            # henry parameters for CO2 in H2O
            #his = [170.713,-8477.7,-21.957, 0.00578]
            his = [91.344,-5876,-8.598, -0.012]
            return his[0] + his[1]/b.T +his[2]*log(b.T) + his[3]*b.T+log(1e-3)
        elif i == 'MEA':
            return log(6.6434e8*1e-3*exp(-896.5/b.T))
    model.lnHis = Expression(Enrtl_sets.Ss, rule=rule_lnHis)


    def rule_lnHmx(b,j):
        return b.ln_ac_inf_mx[j,'CO2'] + \
            b.wh[j,'H2O']*(b.lnHis['H2O'] - b.ln_ac_inf_aq[j,'CO2'])+\
            b.wh[j,'MEA']*(b.lnHis['MEA'] - b.ln_ac_inf_mea[j,'CO2'])

    model.lnHmx = Expression(xa_id,rule=rule_lnHmx,
                            doc='CO2 Henry constant in mixed solvent')

    def rule_lnPCO2(b,j):
        return log(b.xt[j,'CO2']) + b.ln_ac_solute_mx[j,'CO2'] +  b.lnHmx[j]
    model.lnPCO2 = Expression(xa_id,rule=rule_lnPCO2)


    def rule_Habs(b):
        return (b.HL['fin']- b.HL['ref']- b.dn_CO2*b.Hv_CO2)/(b.dn_CO2)
    model.Habs = Expression(rule=rule_Habs)

    return model

def make_idaes_model():
    model = ConcreteModel()
    #model.fs = FlowsheetBlock(time_set=[0])
    model.params = GenericParameterBlock(**liquid_config)
    model.state_block = model.params.build_state_block(
        [0],
        has_phase_equilibrium=False,
        defined_state=True,
    )
    return model

if __name__ == "__main__":
    #VLE data
    #40
    AL40 = [0.102,0.206,0.25,0.337,0.353,0.401,0.417,0.421,
            0.433,0.447,0.464,0.476,0.477,0.485,0.489,0.516,0.524]
    A40 =[0.0016,0.0123,0.0246,0.0603,0.0851,0.1835,0.2928,0.3188,0.3809,
          0.5702, 1.0662,1.8326,1.8278,2.3193,2.8577,8.5583,11.812]

    JL40 =[0.0888,0.203,0.365,0.461,0.513,0.557,0.609,0.646,0.709,
           0.794,0.844,0.965]
    J40 = [0.00147,0.00896,0.0677,0.604,2.57,8.09,36.1,103,293,593,993,
           2992]

    HL40 = [0.153,0.163,0.17,0.191,0.194, 0.232,0.246,0.269,0.272,0.326,0.348,
            0.35,0.36, 0.382,0.386, 0.389,0.4, 0.464,0.466,0.481,0.491, 0.501,
            0.518,0.591]
    H40 = [0.0057, 0.00664, 0.00721, 0.00995, 0.00985, 0.0146,  0.0191,
           0.0231, 0.0224, 0.0485, 0.0662, 0.0721, 0.0966, 0.131, 0.12, 0.113,
           0.128, 0.75, 0.574, 0.883, 1.1, 1.87, 3.03, 28.3]

    #80
    AL80 = [0.017,   0.04,    0.075 ,  0.122,   0.155,   0.216,   0.271,   0.347,
       0.4, 0.398]
    A80 = [0.0056,   0.0219,  0.0557,  0.1406,  0.2485,  0.6137,  1.2538,  3.7522,  7.9387,  8.3031]

    JL80 = [0.0174,  0.0236,  0.0781,  0.118,   0.187,   0.348,   0.46,    0.517,   0.576,0.640,0.680,0.797,0.875,0.953,0.998]
    J80 = [0.00296,  0.00484, 0.0465 , 0.0992,  0.278,   2.67,    16,  56,  235,
           558,958,2951,5943,9925,14902]
    #120
    ML120 = [0.155,  0.1766,  0.1843,  0.2085,  0.2326,  0.2381,  0.256,   0.2901,  0.2967,  0.3004,  0.3125,  0.3191,  0.3298,  0.3424,  0.3424,  0.35,    0.3594,  0.3882,  0.4182]
    M120 = [7.354,   9.314,   9.045,   15.51,   19.62,   25.2,    27.71,   39.18,   40.4 ,   43.49,   51.82,   58.57,   62.88,   77.59,   74.95,   83.61,   92.79 ,  137.9,   191.9]

    JL120 = [0.00333,   0.0112,  0.0247,  0.119,   0.349,   0.403,   0.444 , 0.473,
            0.536,0.644,0.719,0.78,0.829,0.863]
    J120 = [0.00202, 0.0221,  0.0984,  2.29,    46.8 ,   122, 222, 422, 822,2804,
            5809,9770,14741,17723]

    KL120 = [0.07,0.16,0.25,0.33,0.4,0.46,0.48]
    K120 = [1.26,7.65,21.32,54.28,145.30,357.11,430.36]

    #heat of absorption data Kim(2014)
    H_data_40 = [84.03,84.36,85.44,85.36,85.76,81.55,71.30,48.60,40.24,36.67,35.51]
    L_data_40 = [0.06, 0.15,0.23,0.31,0.39,0.44,0.51,0.56,0.61,0.65,0.69]

    # kim 2007
    L_data_40b = [0.041, 0.084,   0.124,   0.167 ,  0.211 ,  0.256,   0.31,    0.352,   0.397,   0.454 ,  0.494 ,  0.535,   0.576 ,  0.608 ,  0.64,    0.662,   0.684 ,  0.704 ,  0.715,   0.037,   0.081,   0.124,   0.165 ,  0.21,    0.25,    0.29 ,   0.337 ,  0.375 ,  0.417,   0.463 ,  0.507 ,  0.543,   0.576,   0.604,   0.634,   0.651 ,  0.678]

    H_data_40b = [83.575, 82.067,  85.574,  85.138,  85.312,  84.531,  84.377,  84.052,  84.494,  84.255,  84.074,  87.098,  52.746,  46.21,   40.548,  37.492,  36.39,   33.897,  36.16,   84.079,  82.123,  81.983,  81.547,  81.533,  83.199,  82.781,  78.333,  86.437,  84.836,  80.703,  73.567,  59.395,  56.821,  51.026,  40.471,  35.974,  32.71]

    H_data_80 = [88.55,88.25,88.52,88.05,88.27,85.55,73.71,58.94,51.05]
    L_data_80 = [0.06,0.14,0.2,0.28,0.36,0.45,0.50,0.55,0.58]

    H_data_120 = [99.08,96.95,100.77,99.28,95.55,86.05,79.74]
    L_data_120 = [0.07,0.16,0.25,0.33,0.4,0.46,0.48]

    #KIm 2007
    L_data_80b = [0.047, 0.09, 0.159,0.203,0.242,0.282,0.333,0.378,0.423,0.464,   0.504,0.545,0.583,0.621,   0.648,   0.043,   0.092,   0.137,   0.181,   0.225 ,  0.276,   0.326,   0.375,   0.424 ,  0.472,   0.518,   0.57,    0.594,   0.605]
    H_data_80b = [93.363,    90.904,  93.74,   93.072 , 90.456 , 92.071,  90.846,  92.954,  94.039,  93.691,  96.132 , 85.178,  67.069 , 54.808,  36.984,  95.258,  91.558,  87.308,  95.685,  90.069,  95.604,  96.949,  94.175,  94.533,  98.246,  86.77,   68.329,  55.471,  42.259]

    #speciation data
    NMR_L = [0.09, 0.21,  0.27,  0.26,  0.39,  0.59,  0.78,  0.95]
    NMR_CO2 = [0,  0, 0, 0, 0, 0, 0.0009,  0.0032]
    NMR_MEA = [0.0985, 0.0715,  0.0563,  0.0559,0.0376,0.0028,0.0012, 0.0014]
    NMR_HCO3 = [0.0001,  0.002 ,0.0003,  0.0004,  0.0012,0.0151,0.0448, 0.0696]
    NMR_CO3 = [0.0001, 0.0027,  0.0004,  0.0004,  0.001, 0.0009,  0.0009,  0.0011]
    NMR_MEACOO = [0.0103,  0.018, 0.0284,  0.0261,  0.0495,  0.0494,0.0388, 0.0312]
    NMR_MEAH = [0.0065 , 0.016, 0.0225,  0.0218,  0.0473,  0.061, 0.0743 , 0.0843]


    #list for CO2 partial pressure
    p40 =[]
    p80 =[]
    p120 =[]
    #list for excess enthalpy
    Hex_40 = []
    Hex_80 = []
    Hex_120 = []
    Hexp_40 = []
    Hexp_80 = []
    Hexp_120 = []
    #list for Heat of absorption
    Habs_40 = []
    Habs_80 = []
    Habs_120 = []
    #list for speciation
    d_CO2 = []
    d_MEA = []
    d_HCO3 = []
    d_MEACOO = []
    d_MEAH = []

    solver = get_solver("ipopt")

    LD = np.linspace(0.01,0.8,100)
    t  = [40,80,120]
    T = [i+273.15 for i in t]
    import pandas as pd

    idaes_model = make_idaes_model()
    comp_dict = {
        "H2O": "H2O",
        "MEA": "MEA",
        "CO2": "CO2",
        "MEA_+": "MEAH+",
        "HCO3_-": "HCO3-",
        "MEACOO_-": "MEACOO-"
    }
    for i in T:
        for j in range(len(LD)):
            data = {'experiment': 277, 'T': i, 'wt': 30.0, 'lld': LD[j],
                   'PCO2': 27.71}
            model = thermodynamic_model(data)
            solver.solve(model, tee=True)

            idaes_model.state_block[0].temperature.value = model.T.value
            idaes_model.state_block[0].mole_frac_phase_comp_true["Liq","H2O"].value = model.xt['ref','H2O'].value
            idaes_model.state_block[0].mole_frac_phase_comp_true["Liq","CO2"].value = model.xt['ref','CO2'].value
            idaes_model.state_block[0].mole_frac_phase_comp_true["Liq","MEA"].value = model.xt['ref','MEA'].value
            idaes_model.state_block[0].mole_frac_phase_comp_true["Liq","MEACOO_-"].value = model.xt['ref','MEACOO-'].value
            idaes_model.state_block[0].mole_frac_phase_comp_true["Liq","HCO3_-"].value = model.xt['ref','HCO3-'].value
            idaes_model.state_block[0].mole_frac_phase_comp_true["Liq","MEA_+"].value = model.xt['ref','MEAH+'].value

            xH2O = model.xt['ref','H2O'].value + model.xt['ref','HCO3-'].value
            xMEA = model.xt['ref','MEA'].value+ model.xt['ref','MEACOO-'].value + model.xt['ref','MEAH+'].value
            xCO2 = model.xt['ref','CO2'].value+ model.xt['ref','MEACOO-'].value + model.xt['ref','HCO3-'].value
            xTot = xH2O + xMEA + xCO2
            idaes_model.state_block[0].mole_frac_phase_comp["Liq","H2O"].value = xH2O/xTot
            idaes_model.state_block[0].mole_frac_phase_comp["Liq","CO2"].value = xCO2/xTot
            idaes_model.state_block[0].mole_frac_phase_comp["Liq","MEA"].value = xMEA/xTot
            idaes_model.state_block[0].mole_frac_comp["H2O"].value = xH2O/xTot
            idaes_model.state_block[0].mole_frac_comp["CO2"].value = xCO2/xTot
            idaes_model.state_block[0].mole_frac_comp["MEA"].value = xMEA/xTot

            for k in ["H2O","MEA"]:
                assert value(abs(idaes_model.state_block[0].Liq_log_gamma[k]-model.ln_ac_solvent["ref", comp_dict[k]]))<1e-2

            for k, k_prime in comp_dict.items():
                if k not in {"H2O", "MEA"}:
                    assert value(abs(idaes_model.state_block[0].Liq_log_gamma[k]-model.ln_ac_solute["ref", comp_dict[k]]))<1e-2
                assert value(abs(idaes_model.state_block[0].Liq_X[k]-model.X["ref", k_prime]))<1e-6
                assert value(abs(idaes_model.state_block[0].Liq_log_gamma_lc_I[k] - model.ln_ac_SR["ref",k_prime])) < 1e-2
                for l, l_prime in comp_dict.items():
                    if k[-1] == "+" and l[-1] == "+":
                        continue
                    elif k[-1] == "-" and l[-1] == "-":
                        continue
                    assert value(abs(idaes_model.state_block[0].Liq_G[k, l] - model.G2["ref", k_prime, l_prime])) < 5e-4

            if i == 40+273.15:
                p40.append(value(exp(model.lnPCO2['ref'])))
                Hex_40.append(value(model.Hex['ref']))
                Habs_40.append(value(-model.Habs))
                #speciation
                d_CO2.append(value(model.xt['ref','CO2']))
                d_MEA.append(value(model.xt['ref','MEA']))
                d_HCO3.append(value(model.xt['ref','HCO3-']))
                d_MEACOO.append(value(model.xt['ref','MEACOO-']))
                d_MEAH.append(value(model.xt['ref','MEAH+']))
            if i == 80+273.15:
                p80.append(value(exp(model.lnPCO2['ref'])))
                Hex_80.append(value(model.Hex['ref']))
                Habs_80.append(value(-model.Habs))
            if i == 120+273.15:
                p120.append(value(exp(model.lnPCO2['ref'])))
                Hex_120.append(value(model.Hex['ref']))
                Habs_120.append(value(-model.Habs))

    #==========================================================================
    # CO2 Partial Pressure
    fig = plt.figure(figsize=(8,7))
    #40
    plt.semilogy(LD,p40,linestyle="-",color='k',label='40: Model')
    plt.semilogy(AL40,A40,marker='o',markersize=8,mec="r",mfc="None",
                 linestyle="",c='r',label='40: Aronu et al.(2011)')
    plt.semilogy(JL40,J40,marker='s',markersize=8,mec="b",mfc="None",
                 linestyle="",c='b',label='40: Jou et al.(1995)')
    plt.semilogy(HL40,H40,marker='^',markersize=8,mec="g",mfc="None",
                 linestyle="",c='g',label='40: Hilliard(2008)')
    #80
    plt.semilogy(LD,p80,linestyle="--",color='k',label='80: Model')
    plt.semilogy(AL80,A80,marker='*',markersize=8,mec="r",mfc="None",
                 linestyle="",c='r',label='80: Aronu et al.(2011)')
    plt.semilogy(JL80,J80,marker='x',markersize=8,mec="b",mfc="None",
                 linestyle="",c='b',label='80: Jou et al.(1995)')
    #120
    plt.semilogy(LD,p120,linestyle=":",color='k',label='120: Model')
    plt.semilogy(JL120,J120,marker='^',markersize=8,mec="b",
                 linestyle="",c='b',label='120: Jou et al.(1995)')
    plt.semilogy(ML120,M120,marker='o',markersize=8,mec="r",
                 linestyle="",c='r',label="120: Ma'mun et al.(2005)")
    plt.semilogy(KL120,K120,marker='d',markersize=8,mec="g",mfc="None",
                 linestyle="",c='g',label="120: Kim et al.(2014)")

    plt.ylabel('CO$_{2}$ pressure, kPa',fontsize=12)
    plt.xlim(0,0.6)
    plt.legend(fontsize=12,loc='center left', bbox_to_anchor=(1, 0.5))
    plt.tick_params(labelsize=12)
    plt.xlabel('Loading, mol CO$_{2}$/mol MEA',fontsize=12)
    plt.tight_layout()
    plt.show()

    # ========================================================================
    # SPECIATION
    sp = plt.figure(figsize=(8, 7))
    ax = sp.add_subplot(111)

    ax.set_title('Speciation at 40$^{\mathrm{o}}$C',fontsize=18,
                 fontweight='bold')
    ax.set_ylabel('Mole fraction',fontsize=12)
    ax.set_xlabel('Loading, mol CO$_{2}$/mol MEA',fontsize=12)

    #model

    lb2=ax.semilogy(LD,d_MEA,linestyle="-",c='b',label='Model: x$_{MEA}$')
    lb3=ax.semilogy(LD,d_MEACOO,linestyle="-.",c='r',label='Model: x$_{MEACOO-}$')
    lb4=ax.semilogy(LD,d_MEAH,linestyle=":",c='g',label='Model: x$_{MEA+}$')
    lb5=ax.semilogy(LD,d_HCO3,linestyle="--",c='k',label='Model: x$_{HCO3-}$')

    #data
    lb7=ax.semilogy(NMR_L,NMR_MEA,marker='o',markersize=10,mec="b",mfc="None",
                 linestyle="",c='b',label='NMR: x$_{MEA}$')
    lb8=ax.semilogy(NMR_L,NMR_MEACOO,marker='s',markersize=10,mec="r",#mfc="None",
                 linestyle="",c='r',label='NMR: x$_{MEACOO-}$')
    lb9=ax.semilogy(NMR_L,NMR_MEAH,marker='*',markersize=10,mec="g",#mfc="None",
                 linestyle="",c='g',label='NMR: x$_{MEA+}$')
    lb10=ax.semilogy(NMR_L,NMR_HCO3,marker='o',markersize=10,mec="k",
                 linestyle="",c='k',label='NMR: x$_{HCO3-}$')


    # HILLIARD
    lb12=ax.semilogy(MEACOO_Hilliard_x,
                MEACOO_Hilliard_y,
                marker='s',markersize=7,mec="r",mfc="None",
                 linestyle="",c='r',label='Hilliard: x$_{MEACOO-}$')

    lb13=ax.semilogy(MEAH_Hilliard_x,
                MEAH_Hilliard_y,
                marker='*',markersize=8,mec="g",mfc="None",
                 linestyle="",c='g',label='Hilliard: x$_{MEA+}$')

    # BUTTINGER
    lb14=ax.semilogy(MEACOO_Buttinger_x,
                MEACOO_Buttinger_y,
                marker='x',markersize=7,mec="r",
                 linestyle="",c='r',
                 label='B${\mathrm{\ddot{o}}}$ttinger: x$_{MEACOO-}$')


    lab  = lb2+lb3+lb4+lb5+lb7+lb8+lb9+lb10+lb12+lb13+lb14
    labels = [l.get_label() for l in lab]
    ax.legend(lab,labels,fontsize=12,
              loc='center left', bbox_to_anchor=(1, 0.5))
    plt.tick_params(labelsize=12)
    plt.xlim(0,0.8)
    plt.xlabel('Loading, mol CO$_{2}$/mol MEA',fontsize=12)
    plt.tight_layout()
    plt.show()

    #==========================================================================
    #HEAT OF ABSORPTION
    fig = plt.figure(figsize=(8,7),)
    #model
    plt.plot(LD,Habs_40,
                 linestyle="-",c='g',label='Model: 40$^{\mathrm{o}}$C')
    plt.plot(LD,Habs_80,
                 linestyle="--",c='b',label='Model: 80$^{\mathrm{o}}$C')
    plt.plot(LD,Habs_120,
                 linestyle=":",c='r',label='Model: 120$^{\mathrm{o}}$C')

    #data
    plt.plot(L_data_40,H_data_40,marker='^',markersize=10,mec="g",mfc="None",
             linestyle="",c='g',label='Dataset: 40A')
    plt.plot(L_data_40b,H_data_40b,marker='^',markersize=10,mec="g",mfc="g",
             linestyle="",c='b',label='Dataset: 40B')
    plt.plot(L_data_80,H_data_80,marker='o',markersize=10,mec="b",mfc="None",
             linestyle="",c='b',label='Dataset: 80A')
    plt.plot(L_data_80b,H_data_80b,marker='o',markersize=10,mec="b",mfc="b",
             linestyle="",c='b',label='Dataset: 80B')
    plt.plot(L_data_120,H_data_120,marker='s',markersize=10,mec="r",mfc="None",
             linestyle="",c='r',label='Dataset: 120A')

    plt.ylabel('Negative H$_{\mathrm{abs}}$, kJ/mol CO$_{2}$',fontsize=12)
    plt.legend(fontsize=12)
    plt.tick_params(labelsize=12)
    plt.xlabel('Loading, mol CO$_{2}$/mol MEA',fontsize=12)
    plt.xlim(0.03,0.6)
    plt.tight_layout()
    plt.show()

    #==========================================================================
    # Excess Enthalpy spanned by true species
    fig = plt.figure(figsize=(13, 9))
    ax1 = fig.add_subplot(111)
    ax1.set_title('True Excess Enthalpy (40-120$^{\mathrm{o}}$C)',
        fontsize=20,
        fontweight='bold')
    ax1.set_ylabel('Excess Enthalpy, kJ/mol',fontsize=20)
    ax1.set_xlabel('Loading, mol CO$_{2}$/mol MEA',fontsize=20)

    lab1=ax1.plot(LD,Hex_40,linestyle="--",c='g',label='Hex at 40 C ')
    lab2=ax1.plot(LD,Hex_80,linestyle="-",c='b',label='Hex at 80 C')
    lab3=ax1.plot(LD,Hex_120,linestyle="-.",c='r',label='Hex at 120 C')

    lab_1 = lab1+lab2+lab3
    labels_1 = [l.get_label() for l in lab_1]
    ax1.legend(lab_1,labels_1,fontsize=20)
    ax1.tick_params(labelsize=20,direction='in')
    plt.tight_layout()

    plt.show()
