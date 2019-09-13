# module: Data Manipulation

import math

molecule = ""
critT, critD, critP, acc, R, M, Rm = [0, 0, 0, 0, 0, 0, 0]


def molData(fluidData, mol, RVal):
    global critT, critP, critD, M, triple, acc, R, Rm, molecule
    (critT, critP, critD, M, triple, acc) = fluidData
    molecule = mol
    R = RVal
    Rm = RVal


def density(Temp, P):
    T = 1 - Temp / critT
    DL = (
        25.994558623605964697845 * T ** (1.000000)
        - 478.64575811196385757285 * T ** (2.000000)
        + 5515.1342679951330865151 * T ** (3.000000)
        - 37526.587391966946597677 * T ** (4.000000)
        + 158477.39820627312292345 * T ** (5.000000)
        - 427355.45491787482751533 * T ** (6.000000)
        + 736859.28930252173449844 * T ** (7.000000)
        - 786188.12076019600499421 * T ** (8.000000)
        + 472843.31603601045208052 * T ** (9.000000)
        - 122584.49532037785684224 * T ** (10.000000)
        + 1.1780271036143130647389
    )
    DV = (
        -43.623794551394261986843 * T ** (1.000000)
        + 1297.1857120402635246137 * T ** (2.000000)
        - 25235.922066945880942512 * T ** (3.000000)
        + 277481.28943618264747784 * T ** (4.000000)
        - 1852260.5602481754031032 * T ** (5.000000)
        + 7783224.4511505197733641 * T ** (6.000000)
        - 20703758.371924288570881 * T ** (7.000000)
        + 33821499.021056793630123 * T ** (8.000000)
        - 30956786.949815854430199 * T ** (9.000000)
        + 12152491.109968656674027 * T ** (10.000000)
        - 0.23298818266597640103122
    )

    DV = (
        -32.684780103039813070609 * T ** (1.000000)
        + 280.32471802420684525714 * T ** (2.000000)
        - 1510.1997961204849616479 * T ** (3.000000)
        + 3703.8931081660693962476 * T ** (4.000000)
        - 3507.3868095829366211547 * T ** (5.000000)
        + 5664.3175341201904302579 * T ** (10.000000)
    )
    DL = (
        -32.684780103039813070609 * T ** (1.000000)
        + 280.32471802420684525714 * T ** (2.000000)
        - 1510.1997961204849616479 * T ** (3.000000)
        + 3703.8931081660693962476 * T ** (4.000000)
        - 3507.3868095829366211547 * T ** (5.000000)
        + 5664.3175341201904302579 * T ** (10.000000)
    )

    # try:
    # 	DVa = math.exp(DV)*critD;
    # 	DLa= DL*critD; #math.exp(DL)*critD
    # except OverflowError:
    # 	D = DL*critD
    # 	print "Overflow %f %f" % (T,D)
    # 	return D

    # if T > 0 and P > critP:
    # 	return DL * critD;
    # if T < 0:
    # 	return math.exp(DV)*critD

    if DV < 0:
        # try:
        D = math.exp(DV) * critD
        # except OverflowError:
        # 	D = DL*critD
        # print "Overflow %f %f" % (T,D)
        # return D
        # print 'V'
    else:
        # D = math.exp(DL)*critD;
        D = 1 + DL * critD
        # 	#print 'L'
        # print (Temp, DLa, DVa, D)
        # if D > 20 or D <0:
        # 	 print (Temp, DL, DV, D)
    return D


def calcDensity(T, P):
    # D = 0.10698034493832267108337E-001 * 1*P**1 * T**1.375000 - 0.27378491773854597181315E-004 * 2*P**2 * T**1.375000 + 0.31548613059452712368735E-007 * 3*P**3 * T**1.000000 + 3753.9363550676634986303 * (1*P**1 * T**0.125000 * exp(-P**1)+ P**2 * T**0.125000 * (-1)*P**0 * exp(-P**1)) - 4104.1407915656318436959 * (1*P**1 * T**0.250000 * exp(-P**1)+ P**2 * T**0.250000 * (-1)*P**0 * exp(-P**1)) + 703.79073920278040077392 * (1*P**1 * T**1.875000 * exp(-P**1)+ P**2 * T**1.875000 * (-1)*P**0 * exp(-P**1)) - 1245.7011880618581471936 * (1*P**1 * T**2.750000 * exp(-P**1)+ P**2 * T**2.750000 * (-1)*P**0 * exp(-P**1)) + 897.04205653273436382733 * (1*P**1 * T**2.875000 * exp(-P**1)+ P**2 * T**2.875000 * (-1)*P**0 * exp(-P**1)) - 13182.531376504457512056 * (2*P**2 * T**0.125000 * exp(-P**1)+ P**3 * T**0.125000 * (-1)*P**0 * exp(-P**1)) + 36939.124040959948615637 * (2*P**2 * T**0.250000 * exp(-P**1)+ P**3 * T**0.250000 * (-1)*P**0 * exp(-P**1)) - 34471.692633433762239292 * (2*P**2 * T**0.375000 * exp(-P**1)+ P**3 * T**0.375000 * (-1)*P**0 * exp(-P**1)) + 10714.875400222023017704 * (2*P**2 * T**0.500000 * exp(-P**1)+ P**3 * T**0.500000 * (-1)*P**0 * exp(-P**1)) + 421.41552012412671501806 * (4*P**4 * T**0.125000 * exp(-P**1)+ P**5 * T**0.125000 * (-1)*P**0 * exp(-P**1)) - 914.79819770409983448189 * (4*P**4 * T**0.250000 * exp(-P**1)+ P**5 * T**0.250000 * (-1)*P**0 * exp(-P**1)) + 510.02344270157743721938 * (4*P**4 * T**0.375000 * exp(-P**1)+ P**5 * T**0.375000 * (-1)*P**0 * exp(-P**1)) - 16.610600348701765938131 * (4*P**4 * T**1.000000 * exp(-P**1)+ P**5 * T**1.000000 * (-1)*P**0 * exp(-P**1)) + 0.48300787515721871345775 * (2*P**2 * T**2.500000 * exp(-P**2)+ P**3 * T**2.500000 * (-2)*P**1 * exp(-P**2)) - 0.54203731763982276881109 * (3*P**3 * T**0.125000 * exp(-P**2)+ P**4 * T**0.125000 * (-2)*P**1 * exp(-P**2)) + 0.10956318549655791414615E-002 * (3*P**3 * T**4.875000 * exp(-P**3)+ P**4 * T**4.875000 * (-3)*P**2 * exp(-P**3)) + 2.0151758106596644459785
    # D =0.16599700296811911726103 * math.log(P) - 1.8228259746722570433519 * math.log(T) + 4.7722058187082154745440 * T - 0.10694443334317172561443E-002 * P**2 + 0.60078126850290943712662E-005 * P**3 - 0.79850625916240425272719E-001 * T**5 + 0.19160644625617791192429E-001 * T**6 + 4.8629635693401471741026 * 1*P**1 * T**1.125000 - 8.8202421772495878116160 * 1*P**1 * T**1.250000 + 3.9934342113725054268514 * 1*P**1 * T**1.375000 + 0.22049900982351808287341E-002 * 2*P**2 * T**1.250000 - 0.18005997492775232225737E-002 * 2*P**2 * T**1.375000 - 0.74124514746363916474470E-005 * 3*P**3 * T**0.625000 + 0.55851005524067735911456E-005 * 3*P**3 * T**0.750000 - 16020194.437006814405322 * (1*P**1 * T**0.125000 * math.exp(-P**1)+ P**2 * T**0.125000 * (-1)*P**0 * math.exp(-P**1)) + 54638902.801874086260796 * (1*P**1 * T**0.250000 * math.exp(-P**1)+ P**2 * T**0.250000 * (-1)*P**0 * math.exp(-P**1)) - 62743349.601962096989155 * (1*P**1 * T**0.375000 * math.exp(-P**1)+ P**2 * T**0.375000 * (-1)*P**0 * math.exp(-P**1)) + 24287817.546043690294027 * (1*P**1 * T**0.500000 * math.exp(-P**1)+ P**2 * T**0.500000 * (-1)*P**0 * math.exp(-P**1)) - 6942566.9620791021734476 * (1*P**1 * T**2.500000 * math.exp(-P**1)+ P**2 * T**2.500000 * (-1)*P**0 * math.exp(-P**1)) + 15343071.260329255834222 * (1*P**1 * T**2.625000 * math.exp(-P**1)+ P**2 * T**2.625000 * (-1)*P**0 * math.exp(-P**1)) - 11432124.796559454873204 * (1*P**1 * T**2.750000 * math.exp(-P**1)+ P**2 * T**2.750000 * (-1)*P**0 * math.exp(-P**1)) + 2868447.2914309846237302 * (1*P**1 * T**2.875000 * math.exp(-P**1)+ P**2 * T**2.875000 * (-1)*P**0 * math.exp(-P**1)) + 2613454.6495949728414416 * (2*P**2 * T**0.125000 * math.exp(-P**1)+ P**3 * T**0.125000 * (-1)*P**0 * math.exp(-P**1)) - 12185294.806305957958102 * (2*P**2 * T**0.250000 * math.exp(-P**1)+ P**3 * T**0.250000 * (-1)*P**0 * math.exp(-P**1)) + 22684146.332617521286011 * (2*P**2 * T**0.375000 * math.exp(-P**1)+ P**3 * T**0.375000 * (-1)*P**0 * math.exp(-P**1)) - 21075362.925843134522438 * (2*P**2 * T**0.500000 * math.exp(-P**1)+ P**3 * T**0.500000 * (-1)*P**0 * math.exp(-P**1)) + 9771994.8955730032175779 * (2*P**2 * T**0.625000 * math.exp(-P**1)+ P**3 * T**0.625000 * (-1)*P**0 * math.exp(-P**1)) - 1808941.1855218114797026 * (2*P**2 * T**0.750000 * math.exp(-P**1)+ P**3 * T**0.750000 * (-1)*P**0 * math.exp(-P**1)) + 6377.5284854477249609772 * (3*P**3 * T**0.125000 * math.exp(-P**1)+ P**4 * T**0.125000 * (-1)*P**0 * math.exp(-P**1)) - 17950.602755125717521878 * (3*P**3 * T**0.250000 * math.exp(-P**1)+ P**4 * T**0.250000 * (-1)*P**0 * math.exp(-P**1)) + 16836.854826272749050986 * (3*P**3 * T**0.375000 * math.exp(-P**1)+ P**4 * T**0.375000 * (-1)*P**0 * math.exp(-P**1)) - 5262.7033198032941072597 * (3*P**3 * T**0.500000 * math.exp(-P**1)+ P**4 * T**0.500000 * (-1)*P**0 * math.exp(-P**1)) - 3.2695037732013205733494 * (4*P**4 * T**0.125000 * math.exp(-P**1)+ P**5 * T**0.125000 * (-1)*P**0 * math.exp(-P**1)) + 3.0624062857157605677116 * (4*P**4 * T**0.250000 * math.exp(-P**1)+ P**5 * T**0.250000 * (-1)*P**0 * math.exp(-P**1)) + 0.46530967190589622717312 * (3*P**3 * T**0.125000 * math.exp(-P**2)+ P**4 * T**0.125000 * (-2)*P**1 * math.exp(-P**2)) - 0.73945859735560309777824E-001 * (3*P**3 * T**2.875000 * math.exp(-P**2)+ P**4 * T**2.875000 * (-2)*P**1 * math.exp(-P**2)) - 3.9887088278569757804348
    # PV = nRT
    # D = P/R/T/1000
    # D = 0.28127164603432880385370 * math.log(P) + 2.4750111455058325660161 * math.log(T) + 0.48688384348486563046876 * math.exp(T) - 0.31799974046157301887927E-002 * P + 0.51456683316407483692032E-005 * P**2 - 0.36339383311418899102918 * T**3

    D = critD
    # VAPOR PHASE
    T_PV = critT / T
    Theta = 1 - T / critT
    PV = (
        -21.588463092387922159787 * T_PV * Theta ** (1.000000)
        + 832.63804785273657671496 * T_PV * Theta ** (2.000000)
        - 16928.315362090459530009 * T_PV * Theta ** (3.000000)
        + 172321.08345621117041446 * T_PV * Theta ** (4.000000)
        - 1004489.4144252967089415 * T_PV * Theta ** (5.000000)
        + 3566067.1640294389799237 * T_PV * Theta ** (6.000000)
        - 7846778.2891382686793804 * T_PV * Theta ** (7.000000)
        + 10450008.859835047274828 * T_PV * Theta ** (8.000000)
        - 7714735.7110287053510547 * T_PV * Theta ** (9.000000)
        + 2422742.8102968716993928 * T_PV * Theta ** (10.000000)
    )
    PV = -PV

    # Saturated Curve
    Ts = 1 - T / critT
    EQ_DV = (
        -43.623794551394261986843 * Ts ** (1.000000)
        + 1297.1857120402635246137 * Ts ** (2.000000)
        - 25235.922066945880942512 * Ts ** (3.000000)
        + 277481.28943618264747784 * Ts ** (4.000000)
        - 1852260.5602481754031032 * Ts ** (5.000000)
        + 7783224.4511505197733641 * Ts ** (6.000000)
        - 20703758.371924288570881 * Ts ** (7.000000)
        + 33821499.021056793630123 * Ts ** (8.000000)
        - 30956786.949815854430199 * Ts ** (9.000000)
        + 12152491.109968656674027 * Ts ** (10.000000)
        - 0.23298818266597640103122
    )
    # EQ_DV =  - 32.684780103039813070609 .* T.^(1.000000) + 280.32471802420684525714 .* T.^(2.000000) - 1510.1997961204849616479 .* T.^(3.000000) + 3703.8931081660693962476 .* T.^(4.000000) - 3507.3868095829366211547 .* T.^(5.000000) + 5664.3175341201904302579 .* T.^(10.000000);
    EQ_DL = (
        25.994558623605964697845 * Ts ** (1.000000)
        - 478.64575811196385757285 * Ts ** (2.000000)
        + 5515.1342679951330865151 * Ts ** (3.000000)
        - 37526.587391966946597677 * Ts ** (4.000000)
        + 158477.39820627312292345 * Ts ** (5.000000)
        - 427355.45491787482751533 * Ts ** (6.000000)
        + 736859.28930252173449844 * Ts ** (7.000000)
        - 786188.12076019600499421 * Ts ** (8.000000)
        + 472843.31603601045208052 * Ts ** (9.000000)
        - 122584.49532037785684224 * Ts ** (10.000000)
        + 1.1780271036143130647389
    )
    # EQ_DL = 1  - 32.684780103039813070609 .* T.^(1.000000) + 280.32471802420684525714 .* T.^(2.000000) - 1510.1997961204849616479 .* T.^(3.000000) + 3703.8931081660693962476 .* T.^(4.000000) - 3507.3868095829366211547 .* T.^(5.000000) + 5664.3175341201904302579 .* T.^(10.000000);

    # SOAVE EQUATION
    # p(T,v) = RT/(v-b) - a(T)/v(v+b)
    Rm = 8.314472 / 345
    # /345
    b = 0.08664 * Rm * critT / critP
    m = 0.480 + 1.574 * acc - 0.176 * acc ** 2
    a = 0.42747 * Rm ** 2 * critT ** 2 / critP * (1 + m * (1 - (T / critT) ** 0.5)) ** 2

    # Z = p/d/R/T
    A = a * P / (Rm ** 2) / (T ** 2)
    B = b * P / Rm / T
    # Y = Z-1/3
    # Y^3 + r*Y + q = 0
    r = (3 * (A - B - B ** 2) - 1) / 3
    q = -2 / 27 + 1 / 3 * (A - B - B ** 2) - A * B
    Det = (r / 3) ** 3 + (q / 2) ** 2

    # num = length(Det);
    # Root = zeros(num,5);
    # Root(:,3) = Det;
    # Root(:,4) = P;
    # Root(:,5) = T;
    if P < (PV - 0.01) or (P + 0.01) > PV:
        if (abs(Det) < 10 ** -8) and (abs(r) < 10 ** -3):
            print("Critical Zone %f, %f" % (T, P))
            D = critD
        # elif(Det>= 0):
        #     u = (-q/2 + Det**.5)**(1/3);
        #     Y1 = u - r/3/u;
        #     dens = P/Rm/T/(Y1+1/3);
        #     print "One Root %f, %f, %f, %f" % (T,P, dens, Det)
        #     D= dens;
        elif Det >= 0:  # <0
            s = abs(r)
            O_var = ((s ** 3) / 27) ** 0.5
            # s should be -r
            val = math.radians(-1 / 2 / O_var)
            # print val
            try:
                Phi = math.acos(val)
            except ValueError:
                return critD
            Y1 = 2 * O_var ** (3 / 2) * math.cos(Phi / 3) + 1 / 3
            Y2 = 2 * O_var ** (3 / 2) * math.cos(Phi / 3 + 2 * math.pi / 3) + 1 / 3
            Y3 = 2 * O_var ** (3 / 2) * math.cos(Phi / 3 + 4 * math.pi / 3) + 1 / 3
            D1 = P / Rm / T / (Y1)
            D2 = P / Rm / T / (Y2)
            D3 = P / Rm / T / (Y3)
            if PV > P:  # vapor phase
                Y = min([D1, D2, D3])
            # print "Vapor Root %f, %f, %f" % (T,P, Y)
            else:  # Liquid
                Y = max([D1, D2, D3]) / M
                if Y > 50:
                    print(Y / M)
                    Y = EQ_DL * critD
                # print "Liquid Root %f, %f, %f" % (T,P, Y)
            if Y < 0:
                if P == -1:
                    Y = EQ_DL * critD
                else:
                    Y = math.exp(EQ_DV) * critD
                    Y = EQ_DL * critD
                # print "Saturated Root %f, %f, %f" % (T,P, Y)
                # if Y > critD:
                # 	Y = critD
                # 	print "CritD"
            D = Y
            # Root(i,3) = (Y2+1/3)*P(i)/Rm/T(i);
            # Root(i,4) = (Y3+1/3)*P(i)/Rm/T(i);
    return D


def PVT(x):
    """
    Calculate dimensionless compressibility
    Inputs:
        X = [Pressure, Density, Temperature]
    OutputS:
        X = [Delta, Tau, Compressibility]
    """
    Pressure = float(x[0])
    Density = float(x[1])
    Temperature = float(x[2])

    Tau = float(critT) / Temperature
    Delta = Density / float(critD)

    if Density != 0:
        Z = (Pressure * 1000 / R / Temperature / Density) - 1.00
    else:
        Z = 0
    return [Delta, Tau, Z]


def P(x):
    Temperature = float(x[2])
    Density = float(x[1])
    Pressure = float(x[0])
    Tau = float(critT) / Temperature
    Delta = Density / float(critD)
    return [Delta, Tau, Pressure]


def CP(x):
    """
    Calculate dimensionless isobaric heat capacity
    Inputs:
        X = [Density, Temperature, Isobaric Heat Capacity]
    Outputs:
        X = [Delta, Tau, CP]
    """
    Delta = float(x[0]) / float(critD)
    Tau = float(critT) / float(x[1])
    CP = float(x[2]) / R

    return [Delta, Tau, CP]


def CV(x):
    """
    Calculate dimensionless isochoric heat capacity
    Inputs:
        X = [Density, Temperature, Isochoric Heat Capacity]
    Outputs:
        X = [Delta, Tau, CV]
    """
    Delta = float(x[0]) / float(critD)
    Tau = float(critT) / float(x[1])
    CV = float(x[2]) / R  # removed /1000
    return [Delta, Tau, CV]


def SND(x):  # m^2/s^2
    """
    Calculate dimensionless speed of sound
    Inputs:
        X = [Density, Temperature, Speed of Sound]
    Outputs:
        X = [Delta, Tau, W]
    """
    Delta = float(x[0]) / float(critD)
    Tau = float(critT) / float(x[1])
    w = (float(x[2]) ** 2) / R / 1000 * M / float(x[1])
    return [Delta, Tau, w]


def CP0(x):  # CP0 given as J/mol K
    Cp0 = float(x[1]) / Rm / 1000.0
    Temp = float(x[0]) / 1000.0
    return [Temp, Cp0]


def DL(x):
    """
    Calculate Theta and Delta for saturated liquid density
    Inputs:
         X = [Density, Temperature]
    """
    Theta = 1 - float(x[1]) / float(critT)
    Delta = float(x[0]) / float(critD) - 1
    # if molecule == "H2O":
    #     # print M
    #     Theta = 1 - float(x[0]) / float(critT)
    #     Delta = float(x[1]) / M / float(critD) - 1
    return [Theta, Delta]


def DV(x):
    """
    Calculate Theta and Delta for saturated vapor density
    Inputs:
         X = [Density, Temperature]
    """
    Theta = 1 - float(x[1]) / float(critT)
    Delta = math.log(float(x[0]) / float(critD))
    # if molecule == "H2O":
        # Theta = 1 - float(x[0]) / float(critT)
        # Delta = math.log(float(x[1]) / M / float(critD))
        # print Delta
    return [Theta, Delta]


def Dsat(x):
    Theta = float(critT) / float(x[1])
    Delta = float(x[0]) / float(critD)
    return [Delta, Theta]


def PV(x):
    """
    Calculate Theta, Tau, and Psi for saturated liquid density
    Inputs:
         X = [Pressure, Temperature]
    """
    Theta = 1 - float(x[1]) / float(critT)
    Psi = math.log(float(x[0]) / float(critP))
    Tau = float(critT) / float(x[1])
    # if molecule == "H2O":
    #     Theta = 1 - float(x[0]) / float(critT)
    #     Psi = math.log(float(x[1]) / float(critP))
    #     Tau = float(critT) / float(x[0])
    return [Tau, Theta, Psi]
