import math
import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt

# import DataImport
# import DataManipulation
# import BasisFunctions

# import helmet
from helmet import DataImport, BasisFunctions, AncillaryEquations


np.seterr(all="ignore")
P = []
T = []
D = []
global rho_start, Res, Tc, Pc, Dc, critT, critD, critP
critT, critD = 0, 0
rho_start = []
Res = []

# global critT, critP, critD


def molData(fluidData, Dmolecule, RVal):
    """ Passing of the Data from the main module ::module:: MPEOSDeveloperModule

        :param fluidData: (critT, critP, critD, M, triple, acc).
        :type fluidData: array
        :param Dmolecule: Name of the molecule of interest.
        :type Dmolecule: str.
        :param RVal: Gas Constant.
        :type RVal: int.
    """

    global critT, critP, critD, M, triple, acc, R, Rm, molecule

    global Tc, Pc, Dc
    molecule = Dmolecule
    (critT, critP, critD, M, triple, acc) = fluidData
    Tc = critT
    Pc = critP
    Dc = critD
    R = RVal
    Rm = RVal


def Vapor_Phase(T):
    global critT, molecule
    # print molecule
    T_PV = critT / float(T)
    Theta = 1 - float(T) / critT

    PVf = AncillaryEquations.getPV()
    PV = np.exp(PVf.f(T_PV, Theta)) * critP

    return PV


def Sat_Liq_Density(Tl):
    global critT, critD
    X = 1 - Tl / critT

    DL = AncillaryEquations.getDL()
    SatLiqDens = (DL.f(X) + 1) * critD
    return SatLiqDens


def Sat_Vap_Density(Tv):
    global critT, critD
    Ts = 1 - float(Tv) / float(critT)
    DV = AncillaryEquations.getDV()

    try:
        SatVapDens = np.exp(np.exp(DV.f(Ts)) * float(critD))
    except OverflowError:
        SatVapDens = critD
    return SatVapDens


def TPcoords(t, lnT, lnp, rlnT=0.1, rlnp=0.1):
    return np.exp(lnT + rlnT * np.cos(t)), np.exp(lnp + rlnp * np.sin(t))


def obj_circle(t, lnT, lnp, rho_guess, Ts, Ti, Ps):
    T2, P2 = TPcoords(t, lnT, lnp)
    r = Dens_calc(lnp, T2, Tc / T2, Ps)
    return r


def Dens_calc1(D, Ts, Ti, Ps):
    Tau = Ts
    Delta = D / Dc
    R = 8.314472
    n = [
        0.96464,
        -2.7855,
        0.86712,
        -0.18860,
        0.11804,
        0.00025181,
        0.57196,
        -0.029287,
        -0.43351,
        -0.12540,
        -0.028207,
        0.014076,
    ]

    basis_d = [
        Tau ** 0.25,  # 1
        Tau ** 1.125,  # 2
        Tau ** 1.5,  # 3
        2 * Delta * Tau ** 1.375,  # 4
        3 * Delta ** 2 * Tau ** 0.25,  # 5
        7 * Delta ** 6 * Tau ** 0.875,  # 6
        Tau ** 0.625 * (2 * Delta * np.exp(-Delta) - Delta ** 2 * np.exp(-Delta)),  # 7
        Tau ** 1.75
        * (5 * Delta ** 4 * np.exp(-Delta) - Delta ** 5 * np.exp(-Delta)),  # 8
        Tau ** 3.625
        * (np.exp(-Delta ** 2) - 2 * Delta ** 2 * np.exp(-Delta ** 2)),  # 9
        Tau ** 3.625
        * (
            -2 * Delta ** 5.0 * np.exp(-Delta ** 2)
            + 4 * Delta ** 3 * np.exp(-Delta ** 2)
        ),  # 10
        Tau ** 14.5
        * (
            3 * Delta ** 2 * np.exp(-Delta ** 3) - 3 * Delta ** 5 * np.exp(-Delta ** 3)
        ),  # 11
        Tau ** 12
        * (4 * Delta ** 3 * np.exp(-Delta ** 3) - 3 * Delta ** 6 * np.exp(-Delta ** 3)),
    ]
    # 12

    # cols = len(basis_d)
    deriv = [x * y for x, y in zip(n, basis_d)]
    denom = (1 + Delta * sum(deriv)) * (D * 1000 * R * Ti)
    Z = Ps
    # print denom

    X = Z - denom
    # print Ps, denom, X
    return X


def Dens_calc(D, Ts, Ti, Ps, Y=[], Beta=[]):
    global critT, critP, critD, M, triple, acc, R, Rm, molecule
    # Tc = critT
    Rm = R
    # Pc = critP
    Dc = critD

    Tau = Ts
    Delta = D / Dc
    R = 8.314472

    if Y == []:
        if molecule == "TOL":
            Beta = [
                0.96464,
                -2.7855,
                0.86712,
                -0.18860,
                0.11804,
                0.00025181,
                0.57196,
                -0.029287,
                -0.43351,
                -0.12540,
                -0.028207,
                0.014076,
            ]
            Y = [2, 9, 12, 23, 26, 79, 124, 202, 240, 327, 376, 385]

            # Y = [3,13,23,24, 32, 45, 138, 173, 190, 238, 327, 364];
            # Beta = [-1.575970, 0.912713, 0.704467, -0.959453, -0.264234, 0.071262, 0.676527, 0.031623, -0.097217, -0.14363, -0.277574, -0.017548]

        # Beta = [1.12166, -4.56888, 0.119044, 0.551162, 0.06269140, 0.00019414, -0.251272, -0.254633, -1.6979, 2.15312, 1.79384, -2.06678]
        # Beta = [1.12165698, -4.56888269, 0.11904442, 0.55116207, 0.06269145, 0.00019414, -0.25127216, -0.25463277, -1.69790030, 2.15311549, 1.79383858, -2.06677957]

        # Y = [12, 49, 61, 74, 81, 95, 96, 97, 165, 215, 225, 240];
        # Beta =[-4.46705, 0.06059820,  -0.02425780, 0.00364612, -0.00179127, -0.00109395, 0.00150141, 3.19540, 4.92228, -3.45665, 5, -0.92253]

        # Y= [2,9,12,23,26,79,124,202,240,327,376,385]
        # Beta = [.96464, -2.7855, 0.86712, -0.18860, 0.11804, 0.0025181, 0.57196, -0.029287, -0.43351, -0.12540, -0.028207, 0.014076];
        # Beta = [1.12166, -4.56888, 0.119044, 0.551162, 0.06269140, 0.00019414, -0.251272, -0.254633, -1.6979, 2.15312, 1.79384, -2.06678]

        # # Y = [1, 9, 10, 14, 21, 27, 39, 66, 68, 163, 168, 174, 215, 217];
        # # Beta = [1.70189, 1.53231, -5, -1.2564, -1.23229, 0.587247, -0.03631580, -0.00509255, 0.0079239, 4.37337, 0.314477, 0.553013, -3.39858, 3.74245];

        # # Y = [12, 26, 49, 61, 74, 81, 95, 96, 97, 145, 165, 215, 225, 240, 266];
        # # Beta = [-5, -0.06159220, 0.04166310, -0.01340710, 0.00156098, 0.00080525, -0.00097637, 0.00085627, 0.561646, 1.7803, 3.98788, -0.299768, 4.59962, -3.23868, -0.0894113]

        # # Y = [12, 27, 28, 49, 61, 74, 81, 95, 96, 97, 165, 215, 240];
        # # Beta = [-4.85306, 0.979153, -0.822724, -0.02513, 0.00617703, -0.00075077, 0.00268564, -0.00085908, 0.00043310, 2.1987,5, -1.5095, 3.1139, -0.836486]

        # # Y = [2,11, 14, 60, 69, 81, 82, 174, 185, 214];
        # # Beta = [0.1771139, -2.2471602, 0.39909357, -0.056842, 0.02005314, -0.00937278, 0.00809853, -1.12192773, 1.48328352, -0.02654376];

        # Y = [1, 9, 20, 27, 31, 34, 40, 43, 46, 56, 58, 70]
        # Beta = [0.836033, -1.7237, -0.211423, 0.109762, -0.256459, .536034, -0.0417975, 0.200962, -0.371925, -0.024484, 0.081429, -0.00469394 ]

        # # Y = [11, 12, 23, 24, 37, 43, 64, 68, 78, 80, 86, 212];
        # # Beta = [-3.3947, 1.88573, 1.22755, -0.789829, 0.072962, -0.114456, -0.001639, 0.004292, 0.003248,-0.001788, -0.000337, -0.009488]
        if molecule == "H2O":
            Y = [12, 56, 82, 85, 87]
            Beta = [-0.509759, 0.001015, 0.000068, 0.000153, -0.000155]

            Y = [1, 12, 23, 28, 46, 60, 61, 64, 72, 75, 212, 270]
            Beta = [
                0.730649,
                -2.45798,
                0.608661,
                -0.087838,
                -0.006013,
                -0.001654,
                0.003015,
                -0.001517,
                0.000689,
                -0.000220,
                0.002403,
                0.018922,
            ]

            Y = [12, 29, 56, 72, 73, 84, 142, 241]
            Beta = [
                -0.822598,
                0.060244,
                -0.006249,
                0.001418,
                0.000143,
                -0.000171,
                -0.261943,
                0.044109,
            ]

            Y = [
                2,
                10,
                13,
                22,
                116,
                149,
                192,
                196,
                229,
                327,
                337,
                358,
                366,
                370,
                371,
                381,
                387,
                396,
            ]
            Beta = [
                0.207943,
                -1.901351,
                0.259927,
                0.255975,
                -1.501166,
                0.476397,
                -0.031241,
                -0.122461,
                1.303569,
                -0.086988,
                0.118378,
                -1.794894,
                3.324689,
                -1.977635,
                -0.329108,
                0.400920,
                0.236498,
                -0.241413,
            ]

            Y = [
                4,
                12,
                38,
                51,
                57,
                61,
                126,
                142,
                166,
                189,
                205,
                231,
                240,
                322,
                359,
                370,
                371,
                398,
            ]
            Beta = [
                0.949598,
                -2.0205,
                0.033214,
                -0.007369,
                0.002074,
                0.000265,
                0.167416,
                -0.389556,
                0.369538,
                -0.17986,
                0.054316,
                1.5576,
                -1.70785,
                -0.052242,
                0.196935,
                -0.116019,
                0.060112,
                -0.003816,
            ]

        else:
            Y = [11, 12, 23, 24, 37, 43, 64, 68, 78, 80, 86, 212]
            Beta = [
                -3.3947,
                1.88573,
                1.22755,
                -0.789829,
                0.072962,
                -0.114456,
                -0.001639,
                0.004292,
                0.003248,
                -0.001788,
                -0.000337,
                -0.009488,
            ]

    basis_d = BasisFunctions.drdRes(Delta, Tau, Y, Beta)

    # basis_d = [ Tau**0.25,     #1
    #           Tau**1.125,      #2
    #           Tau**1.5,        #3
    #            2*Delta*Tau**1.375,    #4
    #            3*Delta**2*Tau**0.25,   #5
    #            7*Delta**6*Tau**0.875,  #6
    #            Tau**0.625*(2*Delta*np.exp(-Delta) - Delta**2*np.exp(-Delta)),                             #7
    #            Tau**1.75*(5*Delta**4*np.exp(-Delta) - Delta**5*np.exp(-Delta)),                           #8
    #            Tau**3.625*(np.exp(-Delta**2) - 2*Delta**2*np.exp(-Delta**2)),                            #9
    #            Tau**3.625*(-2*Delta**5.*np.exp(-Delta**2) + 4*Delta**3*np.exp(-Delta**2)),                #10
    #            Tau**14.5*(3*Delta**2*np.exp(-Delta**3) - 3*Delta**5*np.exp(-Delta**3)),                  #11
    #            Tau**12*(4*Delta**3*np.exp(-Delta**3) - 3* Delta**6*np.exp(-Delta**3))];                  #12

    # cols = len(basis_d);
    # deriv = [x*y for x,y in zip(n,basis_d)]

    denom = (1 + basis_d) * (D / 1000 * R * Ti)
    Z = Ps

    # print denom

    X = Z - denom
    # print Ps, denom, X
    return X


def Dens_calc_min(D, Ts, Ti, Ps, solPrint=False, Y=[], Beta=[]):
    global critT, critP, critD, M, triple, acc, R, Rm, molecule
    Dc = critD

    Tau = Ts
    # if D == 0:
    #   return 100
    Delta = D / Dc
    R = 8.314472

    # print Y
    if Y == []:
        Y = [2, 9, 12, 23, 26, 79, 124, 202, 240, 327, 376, 385]
        Beta = [
            0.96464,
            -2.7855,
            0.86712,
            -0.18860,
            0.11804,
            0.0025181,
            0.57196,
            -0.029287,
            -0.43351,
            -0.12540,
            -0.028207,
            0.014076,
        ]
        Beta = [
            1.12165698,
            -4.56888269,
            0.11904442,
            0.55116207,
            0.06269145,
            0.00019414,
            -0.25127216,
            -0.25463277,
            -1.69790030,
            2.15311549,
            1.79383858,
            -2.06677957,
        ]

        Y = [2, 11, 14, 60, 69, 81, 82, 174, 185, 214]
        Beta = [
            0.1771139,
            -2.2471602,
            0.39909357,
            -0.056842,
            0.02005314,
            -0.00937278,
            0.00809853,
            -1.12192773,
            1.48328352,
            -0.02654376,
        ]

        Y = [1, 9, 20, 27, 31, 34, 40, 43, 46, 56, 58, 70]
        Beta = [
            0.836033,
            -1.7237,
            -0.211423,
            0.109762,
            -0.256459,
            0.536034,
            -0.0417975,
            0.200962,
            -0.371925,
            -0.024484,
            0.081429,
            -0.00469394,
        ]

        Y = [12, 35, 38, 39, 66, 70, 141, 189, 190, 207, 235, 239]
        Beta = [
            -0.368208,
            -0.075016,
            -0.140981,
            0.167270,
            -0.000537,
            0.000272,
            -1.211720,
            1.140863,
            -1.150929,
            0.0164,
            0.514016,
            -0.959177,
        ]

        # CO
        Y = [13, 23, 24, 25, 33, 34, 35, 36, 37, 73, 85, 383]
        Beta = [
            1.397840,
            -2.015520,
            1.042380,
            -1.006470,
            1.907410,
            -3.37696,
            2.025670,
            -0.316729,
            0.256581,
            -0.003715,
            0.000624,
            0.004008,
        ]

        # Y = [12, 26, 49, 61, 74, 81, 95, 96, 97, 145, 165, 215, 225, 240, 266];
        # Beta = [-5, -0.06159220, 0.04166310, -0.01340710, 0.00156098, 0.00080525, -0.00097637, 0.00085627, 0.561646, 1.7803, 3.98788, -0.299768, 4.59962, -3.23868, -0.0894113]

        # Y = [12, 27, 28, 49, 61, 74, 81, 95, 96, 97, 165, 215, 240];
        # Beta = [-4.85306, 0.979153, -0.822724, -0.02513, 0.00617703, -0.00075077, 0.00268564, -0.00085908, 0.00043310, 2.1987,5, -1.5095, 3.1139, -0.836486]

        # Y = [1, 9, 10, 14, 21, 27, 39, 66, 68, 163, 168, 174, 215, 217];
        # Beta = [1.70189, 1.53231, -5, -1.2564, -1.23229, 0.587247, -0.03631580, -0.00509255, 0.0079239, 4.37337, 0.314477, 0.553013, -3.39858, 3.74245];

    basis_d = BasisFunctions.drdRes(Delta, Tau, Y, Beta)

    # n = [ 0.96464, -2.7855, 0.86712, -0.18860, 0.11804, 0.00025181, 0.57196, -0.029287, -0.43351, -0.12540, -0.028207, 0.014076];

    # basis_d = [ Tau**0.25,     #1
    #           Tau**1.125,      #2
    #           Tau**1.5,        #3
    #            2*Delta*Tau**1.375,    #4
    #            3*Delta**2*Tau**0.25,   #5
    #            7*Delta**6*Tau**0.875,  #6
    #            Tau**0.625*(2*Delta*np.exp(-Delta) - Delta**2*np.exp(-Delta)),                             #7
    #            Tau**1.75*(5*Delta**4*np.exp(-Delta) - Delta**5*np.exp(-Delta)),                           #8
    #            Tau**3.625*(np.exp(-Delta**2) - 2*Delta**2*np.exp(-Delta**2)),                            #9
    #            Tau**3.625*(-2*Delta**5.*np.exp(-Delta**2) + 4*Delta**3*np.exp(-Delta**2)),                #10
    #            Tau**14.5*(3*Delta**2*np.exp(-Delta**3) - 3*Delta**5*np.exp(-Delta**3)),                  #11
    #            Tau**12*(4*Delta**3*np.exp(-Delta**3) - 3* Delta**6*np.exp(-Delta**3))];                  #12

    # cols = len(basis_d);
    # N = n*ones(1,cols(2));]
    # print N

    # deriv = [x*y for x,y in zip(n,basis_d)]
    # print deriv
    # print sum(deriv)

    # denom =( 1 + Delta*basis_d)*(D/1000*R*Ti);
    denom = (1 + basis_d) * (D / 1000 * R * Ti)
    Z = Ps

    # if solPrint: print "Z: ", Z, "denom :", denom, "Density: ", D
    X = (Z - denom) ** 2
    return X


def predictRho(P, T, D, indexcheck=-1):

    global rho_start
    global critT, critP, critD, M, triple, acc, R, Rm, molecule
    Tc = critT
    Rm = R
    # print R
    Pc = critP
    Dc = critD
    # print Tc, Pc, acc

    # Rm = 8.314472;
    if molecule == "H2O":
        Rm = 8.314472  # *10**(-6)

    b = 0.08664 * Rm * Tc / Pc
    m = 0.480 + 1.574 * acc - 0.176 * acc ** 2
    a = 0.42747 * (Rm ** 2) * (Tc ** 2) / Pc * ((1 + m * (1 - ((T / Tc) ** 0.5))) ** 2)

    A = a * P / (Rm ** 2) / (T ** 2)
    B = b * P / Rm / T
    r = (3.0 * (A - B - (B ** 2)) - 1.0) / 3.0
    q = -2.0 / 27.0 + 1.0 / 3.0 * (A - B - (B ** 2)) - A * B
    Det = (r / 3) ** 3 + (q / 2) ** 2

    rho_start = []
    num = 0
    # num = len(Det)
    try:
        num = len(Det)
        # print Det
        pass
    except Exception:
        num = 1
        P = [P]
        T = [T]
        D = [D]
        Det = [Det]
        r = [r]
        q = [q]

    for i in range(num):
        if i == indexcheck:
            print("Pressure: ", P[i], "Temperature: ", T[i])
        if P[i] == -1:
            SatLiqDens = Sat_Liq_Density(T[i])
            rho_start.append(SatLiqDens)
        elif P[i] == -2:
            SatVapDens = Sat_Vap_Density(T[i])
            rho_start.append(SatVapDens)
        # SCENARIOS
        # Near Critical
        elif (
            abs(Det[i]) <= 10 ** -6 and abs(r[i]) <= 10 ** -6
        ):  # change det rang 10**-6
            if i == indexcheck:
                print("Part 1: Near Critical: ", Det[i], r[i])
            rho_start.append(Dc)
        # One Root
        elif Det[i] >= 0:  # and (-q[i]/2.0 + Det[i]**.5)>=0:
            if i == indexcheck:
                print("Part 2: One Root", Det[i], r[i], q[i])
            u = (-q[i] / 2.0 + Det[i] ** 0.5) ** (1.0 / 3.0)
            if i == indexcheck:
                print("Part 2: NAN? ", q[i], Det[i], r[i], u)
            Y1 = u - r[i] / 3.0 / u
            # print 'RM',Rm
            Root = 1000 * P[i] / Rm / T[i] / (Y1 + 1.0 / 3.0)
            if np.isnan(Root):
                Root = Dc
            try:
                SatVapDens = Sat_Vap_Density(T[i])
                SatLiqDens = Sat_Liq_Density(T[i])
                if molecule == "H2O":
                    # print molecule
                    SatVapDens = SatVapDens
                    SatLiqDens = SatLiqDens
            except RuntimeWarning:
                if i == indexcheck:
                    print("Part 2: Runtime Warning")
                SatVapDens = 1
                SatLiqDens = 9

            if i == indexcheck:
                print(
                    "Part 2: Root: ",
                    Root,
                    " SatLiq: ",
                    SatLiqDens,
                    " SatVap: ",
                    SatVapDens,
                )
            if Root > SatVapDens and Root < SatLiqDens:
                # if i == indexcheck: print np.exp(Vapor_Phase(T[i]))*Pc, P[i]
                # if T[i]>Tc:
                #   rho_start.append(SatVapDens)
                if Root < SatLiqDens:
                    rho_start.append(SatLiqDens)
                    if i == indexcheck:
                        print("Part 2: SatLiq: ", SatLiqDens)
                else:
                    if i == indexcheck:
                        print("Part 2: SatVap: ", SatVapDens)
                    rho_start.append(SatVapDens)
            else:
                # Yval = Root;
                rho_start.append(Root)

        # Evaluate the multiple roots
        elif Det[i] < 0:  # ) or (-q[i]/2.0 + Det[i]**.5)<0):
            # checkedval = (-q[i]/2.0 + Det[i]**.5)
            if i == indexcheck:
                print("Part 3: Multiple Roots", Det[i])
            O_var = (-(r[i] ** 3.0) / 27.0) ** (1.0 / 2.0)
            Phi = np.arccos(-q[i] / 2.0 / O_var)
            # try:
            #   cosval = -q[i]/2.0/O;
            #   while cosval>1:
            #       cosval=cosval-2;
            #   while cosval <-1:
            #       cosval = cosval +2;
            #       # print cosval
            #   # Phi = np.arccos(-q[i]/2.0/O);
            #   Phi = np.arccos(cosval);
            # except Exception, e:
            #   # rho_start.append(Dc);
            #   Phi = np.arccos(-1.0/2.0/O)
            # if np.isnan(Phi) and i == indexcheck:
            #   print "NAN PHI", r[i], q[i], O
            # rho_start.append(Dc)
            #   Phi = np.arccos(1.0/2.0/O)
            #   print q[i],np.nan, Phi, O, r[i]

            Y1 = 2 * (O_var ** (1.0 / 3.0)) * np.cos(Phi / 3.0) + 1.0 / 3.0
            Y2 = (
                2 * (O_var ** (1.0 / 3.0)) * np.cos(Phi / 3.0 + 2 * np.pi / 3.0) + 1.0 / 3.0
            )
            Y3 = (
                2 * (O_var ** (1.0 / 3.0)) * np.cos(Phi / 3.0 + 4 * np.pi / 3.0) + 1.0 / 3.0
            )
            D1 = 1000 * P[i] / Rm / T[i] / (Y1 + 1.0 / 3.0)
            D2 = 1000 * P[i] / Rm / T[i] / (Y2 + 1.0 / 3.0)
            D3 = 1000 * P[i] / Rm / T[i] / (Y3 + 1.0 / 3.0)
            Yval = max([D1, D2, D3])
            if i == indexcheck:
                print("Part 3: Multiple Roots Ds ", [D1, D2, D3])
            if i == indexcheck:
                print("Part 3: Multiple Roots Ys ", [Y1, Y2, Y3])
            # Dlist = [D1, D2, D3]
            # Ylist = [Y1, Y2, Y3]
            Dreal = []
            # Ypos = []
            PV = np.exp(Vapor_Phase(T[i])) * Pc
            # for j in range(3):
            #   if Ylist[j] >0:
            #       Dval = 1000*P[i]/Rm/T[i]/(Ylist[j]+1.0/3.0);
            #       if len(Ypos) == 0:
            #           Ypos = [Dval]
            #       else:
            #           Ypos.append(Dval);
            #   if (np.isreal(Dlist[j]) or abs(np.imag(Dlist[j])) <10**-8):
            #       Dreal = [Dreal, Dlist[j]];
            if i == indexcheck:
                print("Part 3: Liquid", P[i], PV)  # np.exp(PV)*Pc
            if len(Dreal) == 1:
                val = Dreal(1)
                Yval = np.real(val)
            elif P[i] > np.exp(PV) * Pc * 0.97:  # np.exp(PV)*Pc*0.95):# Liquid
                if i == indexcheck:
                    print("Part 3: Liquid", P[i], PV, np.log(PV / Pc))  # np.exp(PV)*Pc
                # if i == indexcheck: print "Part 3: Multiple Roots", Ypos
                # if len(Ypos) >1:
                Yval = max([D1, D2, D3])
                # else:
                #   try:
                #       Y = Ypos[0];
                #   except Exception, e:
                #       print i, Ypos, Ylist, Dlist, P[i], T[i], D[i], Det[i], (-q[i]/2.0 + Det[i]**.5), q, Phi
                #       raise e
            else:  # %Vapor
                if i == indexcheck:
                    print(
                        "Part 3: Vapor", P[i], np.log(PV / Pc), PV, PV * 0.95, Pc
                    )  # PV, np.exp(PV)*Pc*0.95, Pc
                SatVapDens = Sat_Vap_Density(T[i])
                SatLiqDens = Sat_Liq_Density(T[i])
                if i == indexcheck:
                    print(
                        "Part 3: Vapor",
                        "SatVapDens",
                        SatVapDens,
                        "SatLiqDens",
                        SatLiqDens,
                    )
                # if len(Ypos) >1:
                Yval = min([D1, D2, D3])
                # else:
                #   try:
                #       Y = Ypos[0];
                #   except Exception, e:
                #       print i, Ypos, Ylist, Dlist, P[i], T[i], D[i], Det[i], (-q[i]/2.0 + Det[i]**.5), q, Phi
                #       raise e

            if i == indexcheck:
                print("Part 4: Saturation Replace")
            Yval = np.real(Yval)
            SatVapDens = Sat_Vap_Density(T[i])
            # print SatVapDens
            if SatVapDens > 10 ** 6:
                SatVapDens = 0
            SatLiqDens = Sat_Liq_Density(T[i])
            # DiffVap = abs(Yval - SatVapDens)
            # DiffLiq = abs(SatLiqDens - Yval)
            # print Y, SatVapDens, SatLiqDens, DiffVap, DiffLiq
            if i == indexcheck:
                print(
                    "Part 4: Checks", P[i], np.log(PV / Pc), PV * 0.95, T[i]
                )  # PV, np.exp(PV)*Pc*0.95
            if Yval > SatVapDens and Yval < SatLiqDens:
                if (
                    P[i] > (PV * 0.97) and PV > 0
                ):  # and np.exp(PV)*Pc*0.97 > 10**-3): #and (np.exp(PV*Pc*0.97)) < 10**-4 ): #% 0.95?? %|| DiffLiq>DiffVap ):
                    if i == indexcheck:
                        print(
                            "Part 4: SatLiqDens replace: ", Yval, SatLiqDens, SatVapDens
                        )
                    # if(np.exp(PV)*Pc*0.97 > 10**-3):
                    Yval = SatLiqDens
                    # %disp('Replaced Root with Sat Vapor');
                else:  # %if(DiffLiq<DiffVap)
                    if i == indexcheck:
                        print(
                            "Part 4: SatVapDens replace: ", Yval, SatVapDens, SatLiqDens
                        )
                    # if(np.exp(PV)*Pc*0.97 < 10**-3 and SatVapDens > 10**-2):
                    #   Y = SatLiqDens;
                    # else:
                    Yval = SatVapDens
            # else:
            # if(P[i]>(PV) and PV>0):# and np.exp(PV)*Pc*0.97 > 10**-3): #and (np.exp(PV*Pc*0.97)) < 10**-4 ): #% 0.95?? %|| DiffLiq>DiffVap ):
            #   if i == indexcheck: print "Part 4: SatLiqDens replace: ", Yval, SatLiqDens, SatVapDens
            #   # if(np.exp(PV)*Pc*0.97 > 10**-3):
            #   Yval = SatLiqDens;
            #  #%disp('Replaced Root with Sat Vapor');
            # else: # %if(DiffLiq<DiffVap)
            #   if i == indexcheck: print "Part 4: SatVapDens replace: ", Yval, SatVapDens, SatLiqDens
            #   # if(np.exp(PV)*Pc*0.97 < 10**-3 and SatVapDens > 10**-2):
            #   #   Y = SatLiqDens;
            #   # else:
            #   Yval = SatVapDens;
            # %disp('Replaced Root with Sat Liquid');

            rho_start.append(np.real(Yval))
            # Root[0] = np.real(Y);
        else:
            rho_start.append(Dc)

        # if rho_start[-1]<0: print i, rho_start[-1]
    if num == 1:
        # print num
        return rho_start[0]
    return rho_start


def predictRhoSing(P, T, D):
    global rho_start
    # print "ACC", acc;

    Rm = 8.314472
    b = 0.08664 * Rm * Tc / Pc
    m = 0.480 + 1.574 * acc - 0.176 * acc ** 2
    a = 0.42747 * (Rm ** 2) * (Tc ** 2) / Pc * ((1 + m * (1 - ((T / Tc) ** 0.5))) ** 2)

    # Z = p/d/R/T
    A = a * P / (Rm ** 2) / (T ** 2)
    B = b * P / Rm / T
    # Y = Z-1/3
    # Y**3 + r*Y + q = 0
    r = (3.0 * (A - B - B ** 2) - 1.0) / 3.0
    # print 2.0/27.0
    q = -2.0 / 27.0 + 1.0 / 3.0 * (A - B - (B ** 2)) - A * B
    Det = (r / 3) ** 3 + (q / 2) ** 2

    rho_start = 0

    if abs(Det) <= 10 ** -6 and abs(r) <= 10 ** -1:
        rho_start = Dc
    elif Det >= 0 and (-q / 2.0 + Det ** 0.5) >= 0:
        u = (-q / 2 + Det ** 0.5) ** (1.0 / 3.0)
        Y1 = u - r / 3.0 / u
        Root = 1000 * P / Rm / T / (Y1 + 1.0 / 3.0)
        try:
            Sat_VapD = Sat_Vap_Density(T)
            SatVap_L = Sat_Liq_Density(T)
        except RuntimeWarning:
            Sat_VapD = 1
            SatVap_L = 9

        if Root > Sat_VapD or Root < SatVap_L:
            if Root < SatVap_L:
                rho_start = SatVap_L
            else:
                rho_start = Sat_VapD
        else:
            rho_start = Root

    elif Det < 0 or (-q / 2.0 + Det ** 0.5) < 0:
        O_var = (-r ** 3.0 / 27.0) ** (1.0 / 2.0)
        try:
            Phi = np.arccos(-q / 2.0 / O_var)
        except Exception:
            Phi = np.arccos(-1.0 / 2.0 / O_var)

        Y1 = 2 * O_var ** (1.0 / 3.0) * np.cos(Phi / 3.0) + 1.0 / 3.0
        Y2 = 2 * O_var ** (1.0 / 3.0) * np.cos(Phi / 3.0 + 2 * np.pi / 3.0) + 1.0 / 3.0
        Y3 = 2 * O_var ** (1.0 / 3.0) * np.cos(Phi / 3.0 + 4 * np.pi / 3.0) + 1.0 / 3.0
        D1 = 1000 * P / Rm / T / (Y1 + 1.0 / 3.0)
        D2 = 1000 * P / Rm / T / (Y2 + 1.0 / 3.0)
        D3 = 1000 * P / Rm / T / (Y3 + 1.0 / 3.0)
        Y = max([D1, D2, D3])
        Dlist = [D1, D2, D3]
        Dreal = []
        PV = Vapor_Phase(T)
        for j in range(3):
            if np.isreal(Dlist[j]) or abs(np.imag(Dlist[j])) < 10 ** -8:
                Dreal = [Dreal, Dlist[j]]
        if len(Dreal) == 1:
            val = Dreal(1)
            Y = np.real(val)
        elif P > np.exp(PV) * Pc * 0.95:  # np.exp(PV)*Pc*0.95):# Liquid
            Y = max([D1, D2, D3])
        else:  # %Vapor
            Y = min([D1, D2, D3])
        Y = np.real(Y)
        Sat_VapD = Sat_Vap_Density(T)
        if Sat_VapD > 10 ** 6:
            Sat_VapD = 0
        SatVap_L = Sat_Liq_Density(T)
        # DiffVap = abs(Y - Sat_VapD)
        # DiffLiq = abs(SatVap_L - Y)
        if Y > Sat_VapD and Y < SatVap_L:
            if P > (
                np.exp(PV) * Pc * 0.97
            ):  # and (np.exp(PV*Pc*0.97)) < 10**-4 ): #% 0.95?? %|| DiffLiq>DiffVap ):
                Y = SatVap_L
                # %disp('Replaced Root with Sat Vapor');
            else:  # %if(DiffLiq<DiffVap)
                Y = Sat_VapD
                # %disp('Replaced Root with Sat Liquid');
        rho_start = np.real(Y)
        # Root[0] = np.real(Y);
    else:
        rho_start = Dc

    # if rho_start[0] <0:
    #   rho_start[0] = Dc

    # print rho_start

    return rho_start


def findRhoRootSing(rho_guess, Ti, Ps, D, indexcheck=-1):
    global Res
    Res = 0

    solPrint = False
    # appended = -1

    # a = 2369
    # b = 2378
    index = 0

    Ts = Tc / Ti
    Ps = Ps * 10 ** 6
    try:
        i = 1

        while i < 3:
            rho_lower = rho_guess - i
            rho_upper = rho_guess + i
            if rho_guess - i < 0:
                rho_lower = 0.01
            if (
                Dens_calc(rho_upper, Ts, Ti, Ps) > 0
                and Dens_calc(rho_lower, Ts, Ti, Ps) < 0
            ) or (
                Dens_calc(rho_upper, Ts, Ti, Ps) < 0
                and Dens_calc(rho_lower, Ts, Ti, Ps) > 0
            ):
                if solPrint:
                    print("brent btwn", rho_upper, rho_lower)
                Root = scipy.optimize.brentq(
                    Dens_calc, rho_lower, rho_upper, args=(Ts, Ti, Ps)
                )
                if solPrint:
                    print(Ps, Root)  # realD
                # appended = index
                Res = Root
                break
            elif (
                Dens_calc(rho_upper, Ts, Ti, Ps) > 0
                and Dens_calc(rho_guess, Ts, Ti, Ps) < 0
            ) or (
                Dens_calc(rho_upper, Ts, Ti, Ps) < 0
                and Dens_calc(rho_guess, Ts, Ti, Ps) > 0
            ):
                Root = scipy.optimize.brentq(
                    Dens_calc, rho_guess, rho_upper, args=(Ts, Ti, Ps)
                )
                # print Root
                Res = Root
                break
            elif (
                Dens_calc(rho_lower, Ts, Ti, Ps) > 0
                and Dens_calc(rho_guess, Ts, Ti, Ps) < 0
            ) or (
                Dens_calc(rho_lower, Ts, Ti, Ps) < 0
                and Dens_calc(rho_guess, Ts, Ti, Ps) > 0
            ):
                Root = scipy.optimize.brentq(
                    Dens_calc, rho_lower, rho_guess, args=(Ts, Ti, Ps)
                )
                Res = Root
                break
            else:
                i += 0.01

        if i > 3.0:
            # print index
            Root = scipy.optimize.fsolve(Dens_calc_min, rho_guess, args=(Ts, Ti, Ps))
            # Root = scipy.optimize.fmin_powell(Dens_calc_min, rho_guess, args = (Ts, Ti, Ps), xtol=0.0001, ftol=0.0001, maxiter=None, maxfun=None, full_output=0, disp=0, retall=0, callback=None, direc=None)
            Res = Root[0]

    except RuntimeError:
        Res = Root
    except RuntimeWarning:
        print(i, Root)
        Res = Root["x"]

        index += 1
    return Res


def evalRhoPred(rho_guesses, T, P, D):
    global Res

    # a = 2369
    # b = 2378
    index = 0
    # solPrint = False
    # for rho_guess, Ti, Ps, realD in zip(rho_start, T[a:b], P[a:b], D[a:b]):
    for rho_guess, Ti, Ps, realD in zip(rho_guesses, T, P, D):
        # print "Critical Temperature: ", Tc
        # Ts = Tc / Ti
        # print rho_guess
        # Ps = Ps*10**6;
        # Ps = Ps
        # *10**6

        Res.append(rho_guess)

        index += 1


def findRhoRoot(rho_guesses, T, P, D, indexCheck=-1, Y=[], Beta=[]):
    global Res
    global critT, critP, critD, M, triple, acc, R, Rm, molecule
    Tc = critT
    Rm = R
    # Pc = critP
    # Dc = critD

    # print Beta
    # a = 2369
    # b = 2378
    index = 0
    solPrint = False
    appended = -1
    Res = []
    # for rho_guess, Ti, Ps, realD in zip(rho_start, T[a:b], P[a:b], D[a:b]):
    for rho_guess, Ti, Ps, realD in zip(rho_guesses, T, P, D):
        # print "Critical Temperature: ", Tc
        Ts = Tc / Ti
        # print rho_guess
        # Ps = Ps*10**6;
        Ps = Ps
        # *10**6

        if Ps > 0:
            try:
                i = 1

                if index == indexCheck:
                    solPrint = True
                else:
                    solPrint = False

                if solPrint:
                    print(
                        Dens_calc(realD, Ts, Ti, Ps, Y, Beta),
                        Dens_calc(rho_guess, Ts, Ti, Ps, Y, Beta),
                        Dens_calc(rho_guess - 1, Ts, Ti, Ps, Y, Beta),
                        Dens_calc(rho_guess + 1, Ts, Ti, Ps, Y, Beta),
                    )
                while i < 3:
                    rho_lower = rho_guess - i
                    rho_upper = rho_guess + i
                    if rho_guess - i < 0:
                        rho_lower = 0.01
                    # if solPrint: print Dens_calc(rho_guess,Ts,Ti,Ps), Dens_calc(rho_lower,Ts,Ti,Ps), Dens_calc(rho_upper,Ts,Ti,Ps)
                    if (
                        Dens_calc(rho_upper, Ts, Ti, Ps, Y, Beta) > 0
                        and Dens_calc(rho_lower, Ts, Ti, Ps, Y, Beta) < 0
                    ) or (
                        Dens_calc(rho_upper, Ts, Ti, Ps, Y, Beta) < 0
                        and Dens_calc(rho_lower, Ts, Ti, Ps, Y, Beta) > 0
                    ):
                        if solPrint:
                            print("brent btwn", rho_upper, rho_lower)
                        Root = scipy.optimize.brentq(
                            Dens_calc, rho_lower, rho_upper, args=(Ts, Ti, Ps, Y, Beta)
                        )
                        if solPrint:
                            print(Ps, realD, Root)
                        appended = index
                        Res.append(Root)
                        break
                    if (
                        Dens_calc(rho_upper, Ts, Ti, Ps, Y, Beta) > 0
                        and Dens_calc(rho_guess, Ts, Ti, Ps, Y, Beta) < 0
                    ) or (
                        Dens_calc(rho_upper, Ts, Ti, Ps, Y, Beta) < 0
                        and Dens_calc(rho_guess, Ts, Ti, Ps, Y, Beta) > 0
                    ):
                        if solPrint:
                            print("brent up", rho_upper, rho_guess)
                        Root = scipy.optimize.brentq(
                            Dens_calc, rho_guess, rho_upper, args=(Ts, Ti, Ps, Y, Beta)
                        )
                        if solPrint:
                            print(Ps, realD, Root)
                        appended = index
                        Res.append(Root)
                        break
                    elif (
                        Dens_calc(rho_lower, Ts, Ti, Ps, Y, Beta) > 0
                        and Dens_calc(rho_guess, Ts, Ti, Ps, Y, Beta) < 0
                    ) or (
                        Dens_calc(rho_lower, Ts, Ti, Ps, Y, Beta) < 0
                        and Dens_calc(rho_guess, Ts, Ti, Ps, Y, Beta) > 0
                    ):
                        if solPrint:
                            print("brent down", rho_lower, rho_guess)
                        Root = scipy.optimize.brentq(
                            Dens_calc, rho_lower, rho_guess, args=(Ts, Ti, Ps, Y, Beta)
                        )
                        if solPrint:
                            print(Ps, realD, Root)
                        appended = index
                        Res.append(Root)
                        break
                    else:
                        i += 0.01

                if i > 3:
                    # print index
                    fsRoot = scipy.optimize.fsolve(
                        Dens_calc_min, rho_guess, args=(Ts, Ti, Ps, solPrint, Y, Beta)
                    )
                    # lsRoot = scipy.optimize.leastsq(
                    #     Dens_calc_min, rho_guess, args=(Ts, Ti, Ps, solPrint, Y, Beta)
                    # )
                    # powRoot = scipy.optimize.fmin_powell(
                    #     Dens_calc_min,
                    #     rho_guess,
                    #     args=(Ts, Ti, Ps, solPrint, Y, Beta),
                    #     xtol=0.0001,
                    #     ftol=0.0001,
                    #     maxiter=None,
                    #     maxfun=None,
                    #     full_output=0,
                    #     disp=0,
                    #     retall=0,
                    #     callback=None,
                    #     direc=None,
                    # )
                    # if solPrint: print "after 10", Ps, realD, Root, rho_guess

                    # if solPrint:
                    # rhos = np.linspace(rho_guess-5, rho_guess+5);
                    # plt.plot(rhos, Dens_calc_min(rhos, Ts, Ti, Ps));
                    if appended == index:
                        "Double appended - ******", appended
                    # print realD, rho_guess, fsRoot, lsRoot, powRoot
                    # print 'FSROOT: ', fsRoot
                    Res.append(fsRoot)

            # except ValueError as VE:
            # print "ValueError", i, Root,
            # Res.append(Root)
            except RuntimeError:
                print("Runtime Error", i, Root)
                Res.append(Root)
            except RuntimeWarning:
                print("Runtime Warning", i, Root)
                Res.append(Root["x"])
            # print "RES: ", Res
            index += 1
        else:
            Res.append(rho_guess)

    return Res


def runSoave():
    global molecule
    # print molecule;
    P = []
    T = []
    D = []
    # %SOAVE EQUATION
    # p(T,v) = RT/(v-b) - a(T)/v(v+b)
    # DataImport.filename = 'TOLUENE_edit';
    # DataImport.PVT('TOL')
    BasisFunctions.formCustomBasis()
    # molecule = 'TOL'

    DataImport.filename = molecule
    Fluids = {
        "TOL": (591.75, 4.126, 3.169, 92.13842, 178.0, 0.2657),
        "CO": (132.86, 3.494, 10.85, 28.0101, 68.16, 0.0497),
        "H2O": (647.096, 22.064, 17.8737279956, 18.015268, 273.16, 0.344),
    }
    (Tc, Pc, Dc, M, triple, acc) = Fluids[molecule]
    # print Fluids[molecule]
    # print "ACC",acc;

    DataImport.molData(Fluids[molecule], R)

    global critT, critP, critD
    critT = Tc
    critP = Pc
    critD = Dc

    DataImport.PVT(molecule)
    PVTValues = np.array(DataImport.PVTValues)

    # PVTValues = np.array([[27.094, 27.7181608317, 653.182094024],[27.094, 27.7181608317, 653.182094024]])

    for i in PVTValues[:, 0]:
        P.append(float(i))
    P = np.array(P)
    for i in PVTValues[:, 2]:
        T.append(float(i))
    T = np.array(T)
    for i in PVTValues[:, 1]:
        D.append(float(i))
    D = np.array(D)

    # P.append(float(PVTValues[1,0]));
    # T.append(float(PVTValues[1,2]));
    # D.append(float(PVTValues[1,1]));
    # P.append(float(PVTValues[2,0]));
    # T.append(float(PVTValues[2,2]));
    # D.append(float(PVTValues[2,1]));
    # P = np.array(P);
    # T = np.array(T);
    # D = np.array(D);

    # print max(P)

    evalindexcheck = 2698
    # evalindexcheck = 1

    rho_start = predictRho(P, T, D, evalindexcheck)

    # print 'Length', rho_start

    Res = findRhoRoot(np.array(rho_start), np.array(T), P, D, evalindexcheck)
    # evalRhoPred(rho_start, T, P, D);

    # print len(T)

    Tbad = []
    Dbad = []
    ResBad = []
    startBad = []
    indexBad = []

    index = 0
    for t, x, y, z in zip(T, D, Res, rho_start):
        if np.abs(x - y) > 10:
            Tbad.append(t)
            Dbad.append(x)
            ResBad.append(y)
            startBad.append(z)
            indexBad.append(index)
        index += 1

    for i, d, r, p, g in zip(indexBad, Dbad, Tbad, ResBad, startBad):
        if i == evalindexcheck:
            print(
                "Index: ",
                i,
                " Density: ",
                d,
                " Result: ",
                p,
                " Start: ",
                g,
                " Temp: ",
                r,
            )
        print(
            "Index: ", i, " Density: ", d, " Result: ", p, " Start: ", g, " Temp: ", r
        )

    #  PLOTTING

    fig = plt.figure()
    # fig1 = plt.figure()
    fig2 = plt.figure()
    # fig3 = plt.figure()
    fig.suptitle("Density calculated with Soave EOS", fontweight="bold")
    ax = fig.add_subplot(111)
    # ax1 = fig1.add_subplot(111)
    ax2 = fig2.add_subplot(111)
    # ax3 = fig3.add_subplot(111)

    ax.scatter(T, D, s=10, edgecolor="r", facecolor="w", marker="o", label="data")
    ax.scatter(T, Res, s=10, c="g", marker="o", label="calculated")
    # ax1.scatter(Tbad, Dbad, s=15, c='b', marker="o", label='bad Density')
    # ax1.scatter(Tbad, ResBad, s=15, c='k', marker="o", label='bad Density')
    ax2.scatter(D, Res[:], s=8, c="b", marker="o", label="R")
    ax2.plot(np.linspace(0, np.max(D)), np.linspace(0, np.max(D)))

    # diff = D - Res

    # ax3.scatter(T, diff, s=15, c='b', marker="o", label="R")
    # ax.scatter(Ta, Zval, s=15, c='b', marker='o', labe80'fit')
    ax.set_ylabel("Density (dM)", fontsize=20)
    ax.set_xlabel("Temperature (K)", fontsize=20)

    ax2.set_ylabel("Calculated", fontsize=20)
    ax2.set_xlabel("Measured", fontsize=20)

    # ax1.scatter(Ta, delPcalc, s = 25)
    ax.legend(loc="upper right", fontsize=20)
    # ax.set_title('Data Fit',fontsize = 20);
    # ax1.set_title('Residuals',fontsize = 20)
    # ax1.set_ylabel("Percent Deviation (%)",fontsize = 20)
    # ax1.set_xlabel("Temperature (K)",fontsize = 20)
    # plt.show();

    # SNDValues = []
    # DataImport.SND(molecule)
    # Values = DataImport.SNDValues;
    # #importSND(molecule);
    # DT = []
    # WS = []
    # for x in Values:
    #   Theta = float(critT)/float(x[1]);
    #   Psi = float(x[0])/float(critP);
    #   #z = float(x[2])**2//float(x[1]);
    #   #z = (float(x[2])**2)/R/float(x[1]);
    #   z = float(x[2])
    #   #D = DataManipulation.density(float(x[1]), float(x[0]))
    #   if(float(x[0])>0):
    #       De = DataManipulation.calcDensity(float(x[1]), float(x[0]))
    #       De = De/critD
    #       if De != 1:
    #           SNDValues.append((De, float(x[1]), z, float(x[0])))
    #           DT.append((De, Theta))
    #           WS.append(z);
    # SNDnp = np.asarray(SNDValues);
    # T_SND = SNDnp[:,1]
    # D_SND = SNDnp[:,0]
    # w = SNDnp[:,2]
    # P = SNDnp[:,3]

    # w = (w**2)/R/1000*M/T_SND*float(critT)

    # # fig = plt.figure()
    # # ax = fig.add_subplot(111, projection='3d')

    # ax.scatter(T_SND, D_SND*critD, s=25, c ='r', marker="o", label='NIST data');

    plt.show()


# test
# runSoave()


# Before 12/8
def findRhoRoot1(rho_guesses, T, P, D):
    global Res

    # a = 2369
    # b = 2378
    index = 0
    # for rho_guess, Ti, Ps, realD in zip(rho_start, T[a:b], P[a:b], D[a:b]):
    for rho_guess, Ti, Ps, realD in zip(rho_guesses, T, P, D):
        Ts = Tc / Ti
        print(rho_guess)
        # Ps = Ps*10**6;
        Ps = Ps
        # *10**6
        try:
            i = 1

            while i < 10:
                rho_lower = rho_guess - i
                rho_upper = rho_guess + i
                if rho_guess - i < 0:
                    rho_lower = 0.01
                if (
                    Dens_calc(rho_upper, Ts, Ti, Ps) > 0
                    and Dens_calc(rho_guess, Ts, Ti, Ps) < 0
                ) or (
                    Dens_calc(rho_upper, Ts, Ti, Ps) < 0
                    and Dens_calc(rho_guess, Ts, Ti, Ps) > 0
                ):
                    # if(index == 2328):
                    #   print 'brent up', rho_upper, rho_guess
                    Root = scipy.optimize.brentq(
                        Dens_calc, rho_guess, rho_guess + i, args=(Ts, Ti, Ps)
                    )
                    print(Ps, Root)
                    Res.append(Root)
                    break
                elif (
                    Dens_calc(rho_lower, Ts, Ti, Ps) > 0
                    and Dens_calc(rho_guess, Ts, Ti, Ps) < 0
                ) or (
                    Dens_calc(rho_lower, Ts, Ti, Ps) < 0
                    and Dens_calc(rho_guess, Ts, Ti, Ps) > 0
                ):
                    # if(index == 2328):
                    #   print 'brent down', rho_lower, rho_guess
                    Root = scipy.optimize.brentq(
                        Dens_calc, rho_guess - i, rho_guess, args=(Ts, Ti, Ps)
                    )
                    print(Ps, Root)
                    Res.append(Root)
                    break
                else:
                    i += 0.03

            if i > 10.0:
                # print index
                Root = scipy.optimize.fsolve(
                    Dens_calc_min, rho_guess, args=(Ts, Ti, Ps)
                )
                # Root = scipy.optimize.fmin_powell(Dens_calc_min, rho_guess, args = (Ts, Ti, Ps), xtol=0.0001, ftol=0.0001, maxiter=None, maxfun=None, full_output=0, disp=0, retall=0, callback=None, direc=None)
                print("after 10", Ps, Root)
                Res.append(Root)
        except ValueError:
            print("ValueError", i, Root)
            Res.append(Root)
        except RuntimeError:
            print("Runtime Error", i, Root)
            Res.append(Root)
        except RuntimeWarning:
            print("Runtime Warning", i, Root)
            Res.append(Root["x"])

        index += 1

        # removed 12/18


def Dens_calc_min1(D, Ts, Ti, Ps):
    Tau = Ts
    # if D == 0:
    #   return 100
    Delta = D / Dc
    R = 8.314472

    n = [
        0.96464,
        -2.7855,
        0.86712,
        -0.18860,
        0.11804,
        0.00025181,
        0.57196,
        -0.029287,
        -0.43351,
        -0.12540,
        -0.028207,
        0.014076,
    ]

    basis_d = [
        Tau ** 0.25,  # 1
        Tau ** 1.125,  # 2
        Tau ** 1.5,  # 3
        2 * Delta * Tau ** 1.375,  # 4
        3 * Delta ** 2 * Tau ** 0.25,  # 5
        7 * Delta ** 6 * Tau ** 0.875,  # 6
        Tau ** 0.625 * (2 * Delta * np.exp(-Delta) - Delta ** 2 * np.exp(-Delta)),  # 7
        Tau ** 1.75
        * (5 * Delta ** 4 * np.exp(-Delta) - Delta ** 5 * np.exp(-Delta)),  # 8
        Tau ** 3.625
        * (np.exp(-Delta ** 2) - 2 * Delta ** 2 * np.exp(-Delta ** 2)),  # 9
        Tau ** 3.625
        * (
            -2 * Delta ** 5.0 * np.exp(-Delta ** 2)
            + 4 * Delta ** 3 * np.exp(-Delta ** 2)
        ),  # 10
        Tau ** 14.5
        * (
            3 * Delta ** 2 * np.exp(-Delta ** 3) - 3 * Delta ** 5 * np.exp(-Delta ** 3)
        ),  # 11
        Tau ** 12
        * (4 * Delta ** 3 * np.exp(-Delta ** 3) - 3 * Delta ** 6 * np.exp(-Delta ** 3)),
    ]
    # 12

    # cols = len(basis_d)
    # N = n*ones(1,cols(2));]
    # print N

    deriv = [x * y for x, y in zip(n, basis_d)]
    # print deriv
    # print sum(deriv)

    denom = (1 + Delta * sum(deriv, 1)) * (D * 1000 * R * Ti)
    Z = Ps

    X = (Z - denom) ** 2
    return X
