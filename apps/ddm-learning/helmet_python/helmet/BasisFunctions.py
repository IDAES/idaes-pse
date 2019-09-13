# module GAMS Terms

# import math
# from pylab import *
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import sympy as sy

# y = sy.symbols('x')

drd_vals = []
ar_vals = []
ai_vals = []
itt_val = 0
rtt_vals = []
dtrdt_vals = []
d2rd_vals = []
d3rd_vals = []
d4rd_vals = []
d5rd_Vals = []
critT, critD, critP, acc, R, M, Rm = [0, 0, 0, 0, 0, 0, 0]
molecule = ""
coeffs = []
indexes = []


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
    molecule = Dmolecule
    (critT, critP, critD, M, triple, acc) = fluidData
    R = RVal
    Rm = RVal


# def formCustomMixBasis(LemJac=False):
#     """ Basis Functions developed a bank of terms based on literature (Lemmon, Span, Wagner) """
#     # print "terms"
#     global coeffs, indexes
#     coeffs = []
#     if LemJac:
#         for i in range(1, 6):
#             for j in range(1, 13):  # 12
#                 coeffs.append([i, j / 8.0, 0, 0])
#         indexes.append(len(coeffs))
#         for i in range(1, 6):
#             for j in range(1, 24):  # 24
#                 coeffs.append([i, j / 8.0, 1, 0])
#         for i in range(1, 6):
#             for j in range(1, 30):  # 24
#                 coeffs.append([i, j / 8.0, 2, 0])
#         for i in range(2, 5):
#             for j in range(24, 38):  # 38
#                 coeffs.append([i, j / 2.0, 3, 0])
#         indexes.append(len(coeffs))
#         for i in range(1, 6):
#             for j in range(1, 24):  # 24
#                 for m in range(1, 7):
#                     coeffs.append([i, j / 8.0, 1, m / 2])
#         for i in range(1, 6):
#             for j in range(1, 30):  # 24
#                 for m in range(1, 7):
#                     coeffs.append([i, j / 8.0, 2, m / 2])
#         for i in range(2, 5):
#             for j in range(24, 38):  # 38
#                 for m in range(1, 7):
#                     coeffs.append([i, j / 2.0, 3, m / 2])
#         indexes.append(len(coeffs))
#     else:
#         for i in range(1, 6):
#             for j in range(1, 13):  # 12
#                 coeffs.append([i, j / 8.0, 0, 0])
#         indexes.append(len(coeffs))
#         for i in range(1, 6):
#             for j in range(1, 24):  # 24
#                 coeffs.append([i, j / 8.0, 1, 0])
#         for i in range(1, 6):
#             for j in range(1, 30):  # 24
#                 coeffs.append([i, j / 8.0, 2, 0])
#         for i in range(2, 5):
#             for j in range(24, 38):  # 38
#                 coeffs.append([i, j / 2.0, 3, 0])
#         indexes.append(len(coeffs))


def formCustomBasis(LemJac=False):
    """ Basis Functions developed a bank of terms based on literature (Lemmon, Span, Wagner) """
    # print "terms"
    global coeffs, indexes
    coeffs = []
    if LemJac:
        for i in range(1, 9):
            for j in range(1, 13):  # 12
                coeffs.append([i, j / 8.0, 0, 0])
        indexes.append(len(coeffs))
        for i in range(1, 6):
            for j in range(1, 24):  # 24
                coeffs.append([i, j / 8.0, 1, 0])
        for i in range(1, 6):
            for j in range(1, 30):  # 24
                coeffs.append([i, j / 8.0, 2, 0])
        for i in range(2, 5):
            for j in range(24, 38):  # 38
                coeffs.append([i, j / 2.0, 3, 0])
        indexes.append(len(coeffs))
        for i in range(1, 6):
            for j in range(1, 24):  # 24
                for m in range(1, 7):
                    coeffs.append([i, j / 8.0, 1, m / 2])
        for i in range(1, 6):
            for j in range(1, 30):  # 24
                for m in range(1, 7):
                    coeffs.append([i, j / 8.0, 2, m / 2])
        for i in range(2, 5):
            for j in range(24, 38):  # 38
                for m in range(1, 7):
                    coeffs.append([i, j / 2.0, 3, m / 2])
        indexes.append(len(coeffs))
    else:
        for i in range(1, 9):
            for j in range(1, 13):  # 12
                coeffs.append([i, j / 8.0, 0, 0])
        indexes.append(len(coeffs))
        for i in range(1, 6):
            for j in range(1, 24):  # 24
                coeffs.append([i, j / 8.0, 1, 0])
        for i in range(1, 6):
            for j in range(1, 30):  # 24
                coeffs.append([i, j / 8.0, 2, 0])
        for i in range(2, 5):
            for j in range(24, 38):  # 38
                coeffs.append([i, j / 2.0, 3, 0])
        indexes.append(len(coeffs))


# def formCustomBasesFlex():
#     global coeffs, indexes
#     coeffs = []

#     coeffs.append([5, 1, 3, 0])
#     coeffs.append([5, 2, 1, 0])
#     coeffs.append([5, 1, 1, 0])
#     coeffs.append([0.0001, 1, 3, 0])
#     coeffs.append([4.11737, 1, 0, 0])
#     coeffs.append([2.49949, 3, 0, 0])
#     coeffs.append([2.653320, 2, 0, 0])
#     coeffs.append([0.0001, 3, 0, 0])
#     coeffs.append([2.893860, 1, 0, 0])
#     coeffs.append([4.41217, 2, 1, 0])
#     coeffs.append([1.25253, 3, 0, 0])
#     coeffs.append([0.0001, 2, 0, 0])


# def formSpanWagnerNonPolarBasis():
#     global coeffs, indexes
#     coeffs = []

#     coeffs.append([1, 0.25, 0])
#     coeffs.append([1, 1.125, 0])
#     coeffs.append([1, 1.5, 0])
#     coeffs.append([2, 1.375, 0])
#     coeffs.append([3, 0.25, 0])
#     coeffs.append([7, 0.875, 0])

#     # for i in range(1,9):
#     #   for j in range(1,13): #12
#     #       coeffs.append([i,j/8.0, 0])

#     # if(i==1 and j/8.0 == 0.25):
#     #   print len(coeffs)
#     # if(i==1 and j/8.0 == 1.125):
#     #   print len(coeffs)
#     # if(i==1 and j/8.0 == 1.5):
#     #   print len(coeffs)
#     # if(i==2 and j/8.0 == 1.375):
#     #   print len(coeffs)
#     # if(i==3 and j/8.0 == 0.25):
#     #   print len(coeffs)
#     # if(i==7 and j/8.0 == 0.875):
#     #   print len(coeffs)

#     # print("out")

#     # d.append(i);
#     # t.append(j/8.0);
#     # c.append(0);
#     # x = "\t %s^%d * %s^%f\n" %(CritDensLabel, i, CritTempLabel, j/8.0)
#     # k+=1
#     # textFile.write(x);

#     coeffs.append([2, 0.625, 1])
#     coeffs.append([5, 1.75, 1])
#     # indexes.append(len(coeffs))
#     # for i in range(1,6):
#     #   for j in range(1,24): #24
#     #       coeffs.append([i,j/8.0, 1])
#     # if(i==2 and j/8.0 == 0.625):
#     #   print len(coeffs)
#     # if(i==5 and j/8.0 == 1.75):
#     #   print len(coeffs)
#     # print("out")

#     # d.append(i);
#     # t.append(j/8.0);
#     # c.append(1)
#     # y = "\t %s^%d * %s^%f * np.exp(-%s^%d)\n" %(CritDensLabel,i, CritTempLabel, j/8.0, CritDensLabel,1)
#     # k +=1
#     # textFile.write(y);

#     coeffs.append([1, 3.625, 2])
#     coeffs.append([4, 3.625, 2])

#     # for i in range(1,6):
#     #   for j in range(1,30): #24
#     #       coeffs.append([i,j/8.0, 2])
#     # if(i==1 and j/8.0 == 3.625):
#     #   print len(coeffs)
#     # if(i==4 and j/8.0 == 3.625):
#     #   print len(coeffs)
#     # print("out")

#     # d.append(i);
#     # t.append(j/8.0);
#     # c.append(2)
#     # k+= 1
#     # y = "\t %s^%d * %s^%f * np.exp(-%s^%d)\n" %(CritDensLabel,i, CritTempLabel, j/8.0, CritDensLabel,2)
#     # textFile.write(y);

#     coeffs.append([3, 14.5, 3])
#     coeffs.append([4, 12.0, 3])

#     # for i in range(2,5):
#     #   for j in range(24,38): #38
#     #       coeffs.append([i,j/2.0, 3])
#     # if(i==3 and j/2.0 == 14.5):
#     #   print len(coeffs)
#     # if(i==4 and j/2.0 == 12.0):
#     #   print len(coeffs)
#     # d.append(i);
#     # t.append(j/8.0);
#     # c.append(1)
#     # k+=1
#     # y = "\t %s^%d * %s^%f * np.exp(-%s^%d)\n" %(CritDensLabel,i, CritTempLabel, j/2.0, CritDensLabel,3)
#     # textFile.write(y);
#     indexes.append(len(coeffs))
#     # print len(coeffs)


# def formIdealBasis():
#     # a = [3.5241174832, 1.1360823464]
#     c = 4.0
#     # B = [
#     #     0.96464,
#     #     -2.7855,
#     #     0.86712,
#     #     -0.18860,
#     #     0.11804,
#     #     0.0025181,
#     #     0.57196,
#     #     -0.029287,
#     #     -0.43351,
#     #     -0.12540,
#     #     -0.028207,
#     #     0.014076,
#     # ]
#     v = [1.6994, 8.0577, 17.059, 8.4567, 8.6423]
#     u = [190, 797, 1619, 3072, 7915]

#     T = CVnp[:, 1]
#     iThetaTT = -np.ones_like(T) * (c - 1) / T ** 2
#     for uit, vi in zip(u, v):
#         ui = uit / float(critT)
#         iThetaTT = iThetaTT - vi * (ui) ** 2 * np.exp(-ui * T) * (1 - np.exp(-ui * T)) ** (-2)
#     return iThetaTT


def getTerm(Y):
    for y in Y:
        print(y)
        d, t, c, m = coeffs[y - 1]
        print("Delta^%f Tau^%f np.exp(-Delta^%f) np.exp(-Tau^%f)" % (d, t, c, m))


def ideal(D, T):
    global ai_vals

    # TOLUENE
    a = [3.5241174832, 1.1360823464]
    c = [4.0, 0, 0]
    # B = [
    #     0.96464,
    #     -2.7855,
    #     0.86712,
    #     -0.18860,
    #     0.11804,
    #     0.0025181,
    #     0.57196,
    #     -0.029287,
    #     -0.43351,
    #     -0.12540,
    #     -0.028207,
    #     0.014076,
    # ]
    v = [1.6994, 8.0577, 17.059, 8.4567, 8.6423]
    u = [190, 797, 1619, 3072, 7915]

    # Carbon Monoxide
    a = [-3.3728318564, 3.3683460039]
    c = [3.5, 0.22311e-6, 1.5]
    # B = [
    #     0.90554,
    #     -2.4515,
    #     0.53149,
    #     0.024173,
    #     0.072156,
    #     0.00018818,
    #     0.19405,
    #     -0.043268,
    #     -0.12778,
    #     -0.027896,
    #     -0.03414,
    #     0.016329,
    # ]
    v = [1.0128]
    u = [3089]

    idealA = (
        a[0]
        + a[1] * T
        + np.log(D)
        + (c[0] - 1) * np.log(T)
        - c[1] * (critT) ** c[2] / (c[2] * (c[2] + 1)) * T ** (-c[2])
    )
    for uit, vi in zip(u, v):
        ui = uit / float(critT)
        idealA = idealA + vi * np.log(1 - np.exp(-ui * T))
    ai_vals = idealA


def ar(D, T):
    global ar_vals
    global coeffs
    for d, t, c, m in coeffs:
        val = 0
        if c == 0:
            val = (D ** d) * (T ** t)
        elif m == 0:
            val = np.exp(-D ** c) * (D ** (d)) * (T ** t)
        else:
            val = np.exp(-D ** c) * (D ** (d)) * (T ** t) * np.exp(-T ** m)
        ar_vals.append(val)


def idealBY(D, T, Y, Beta):

    # TOLUENE
    if molecule == "TOL":
        a = [3.5241174832, 1.1360823464]
        c = [4.0, 0, 0]
        # B = [
        #     0.96464,
        #     -2.7855,
        #     0.86712,
        #     -0.18860,
        #     0.11804,
        #     0.0025181,
        #     0.57196,
        #     -0.029287,
        #     -0.43351,
        #     -0.12540,
        #     -0.028207,
        #     0.014076,
        # ]
        v = [1.6994, 8.0577, 17.059, 8.4567, 8.6423]
        u = [190, 797, 1619, 3072, 7915]

    # Carbon Monoxide
    if molecule == "CO":
        a = [-3.3728318564, 3.3683460039]
        c = [3.5, 0.22311e-6, 1.5]
        # B = [
        #     0.90554,
        #     -2.4515,
        #     0.53149,
        #     0.024173,
        #     0.072156,
        #     0.00018818,
        #     0.19405,
        #     -0.043268,
        #     -0.12778,
        #     -0.027896,
        #     -0.03414,
        #     0.016329,
        # ]
        v = [1.0128]
        u = [3089]

    # Water
    if molecule == "H2O":
        a = [-8.32044648201, 6.6832105268]
        c = [3.00632 + 1, 0, 0]
        v = [0.012436, 0.97315, 1.27950, 0.96956, 0.24873]
        u = [
            1.2878967 * float(critT),
            3.53734222 * float(critT),
            7.74073708 * float(critT),
            9.24437796 * float(critT),
            27.5075105 * float(critT),
        ]

    if c[2] > 0:
        idealA = (
            a[0]
            + a[1] * T
            + np.log(D)
            + (c[0] - 1) * np.log(T)
            - c[1] * (critT) ** c[2] / (c[2] * (c[2] + 1)) * T ** (-c[2])
        )
    else:
        idealA = a[0] + a[1] * T + np.log(D) + (c[0] - 1) * np.log(T)
    for uit, vi in zip(u, v):
        ui = uit / float(critT)
        idealA = idealA + vi * np.log(1 - np.exp(-ui * T))
    ai_vals = idealA
    return ai_vals


def arBY(D, T, Y, Beta):
    global drd_vals
    global coeffs
    # drd_vals = [];
    val = 0
    # print len(Y)
    if type(Y) is not int:
        for x, y in zip(Y, Beta):
            x = x - 1
            [d, t, c, m] = coeffs[x]
            if c == 0:
                try:
                    val += y * (D ** d) * (T ** t)
                    pass
                except Exception:
                    print("Term", x, d, t, T, D)
            elif m == 0:
                val = np.exp(-D ** c) * (D ** (d)) * (T ** t)
            else:
                val = np.exp(-D ** c) * (D ** (d)) * (T ** t) * np.exp(-T ** m)

        return val


# @staticmethod
def drd(D, T):
    global drd_vals
    global coeffs

    drd_vals = []
    for d, t, c, m in coeffs:
        val = 0
        if c == 0:
            val = d * (D ** d) * (T ** t)
        elif m == 0:
            # x,y = sy.symbols('x, y')
            # expr = x**d * y**t * sy.exp(-(x**c))
            # f_prime = expr.diff(x,1)
            # checkVal = f_prime.subs([(x,D),(y,T)])

            val = np.exp(-D ** c) * ((D ** (d)) * (T ** t) * (d - c * (D ** c)))
            # print D*checkVal, val
        else:
            val = (
                np.exp(-D ** c) * (D ** (d)) * (T ** t) * np.exp(-T ** m) * (d - c * (D ** c))
            )
        drd_vals.append(val)
        # print val;
    # print len(drd_vals);


# @staticmethod
def d2rd(D, T):
    global d2rd_vals
    global coeffs
    d2rd_vals = []
    for d, t, c, m in coeffs:
        val = 0
        if c == 0:
            val = d * (d - 1) * (D ** d) * (T ** t)
        elif m == 0:  # TODO
            val = np.exp(-D ** c) * (
                (D ** (d))
                * (T ** t)
                * ((d - c * (D ** c)) * (d - 1 - c * (D ** c)) - (c ** 2) * (D ** c))
            )
            # Base = (D ** (d)) * (T ** (t)) * np.exp(-(D ** c))
            # A = (d - c * (D ** c)) * (d - 1 - c * (D ** c)) - (c ** 2) * (D ** c)
            # valCheck = Base * A

            # x,y = sy.symbols('x, y')
            # expr = x**d * y**t * sy.exp(-(x**c))
            # f_prime = expr.diff(x,2)
            # checkVal = f_prime.subs([(x,D),(y,T)])

            # print val, valCheck, checkVal*(D**2)

        else:
            val = np.exp(-D ** c) * (
                (D ** (d))
                * (T ** t)
                * np.exp(-T ** m)
                * ((d - c * (D ** c)) * (d - 1 - c * (D ** c)) - (c ** 2) * (D ** c))
            )
        d2rd_vals.append(val)
        # print val;
    # print len(drd_vals);


def d2rdt(D, T):
    global d2rdt_vals
    global coeffs
    # d2rd_vals = []
    for d, t, c, m in coeffs:
        val = 0
        if c == 0:
            val = d * (d - 1) * t * (D ** d) * (T ** t)
        elif m == 0:  # TODO
            val = np.exp(-D ** c) * (
                (D ** (d))
                * (T ** t)
                * ((d - c * (D ** c)) * (d - 1 - c * (D ** c)) - (c ** 2) * (D ** c))
            )
            # Base = t * (D ** (d)) * (T ** (t)) * np.exp(-(D ** c))
            # A = (d - c * (D ** c)) * (d - 1 - c * (D ** c)) - (c ** 2) * (D ** c)
            # valCheck = Base * A

            # x,y = sy.symbols('x, y')
            # expr = x**d * y**t * sy.exp(-(x**c))
            # f_prime = expr.diff(x,2)
            # checkVal = f_prime.subs([(x,D),(y,T)])

            # print val, valCheck, checkVal*(D**2)

        else:
            val = (
                np.exp(-D ** c)
                * (
                    (D ** (d))
                    * (T ** t)
                    * np.exp(-T ** m)
                    * (
                        (d - c * (D ** c)) * (d - 1 - c * (D ** c))
                        - (c ** 2) * (D ** c)
                    )
                )
                * (t - m * (T ** m))
            )
        d2rdt_vals.append(val)
        # print val;
    # print len(drd_vals);


# @staticmethod
def d3rd(D, T):
    global d3rd_vals
    global coeffs
    d3rd_vals = []
    # i = 1
    for d, t, c, m in coeffs:
        val = 0
        if c == 0:
            val = d * (d - 1) * (d - 2) * (D ** d) * (T ** t)
        elif m == 0:
            # val = (D**(d))*(T**t)*exp(-(D**c))*(d*(d-1)*(d-2) + (D**c)*(-2*c + 6*d*c-3*(d**2)*c - 3*d*(c**2) + 3*(c**2) - (c**3)) + (D**(2*c))*(3*d*(c**2) -3*(c**2) + 3*(c**3)) - (c**3)*(D**(3*c)))
            val = (
                (D ** (d))
                * (T ** (t))
                * np.exp(-(D ** c))
                * (
                    (d - c * (D ** c)) * (d - 1 - c * (D ** c)) * (d - 2 - c * (D ** c))
                    + (3 * d - 3 + c - 3 * c * (D ** c)) * (-(c ** 2) * (D ** (c)))
                )
            )

            # Base = (D ** (d)) * (T ** (t)) * np.exp(-(D ** c))
            # A = (d - c * (D ** c)) * (d - 1 - c * (D ** c)) * (d - 2 - c * (D ** c))
            # B = -(c ** 2) * (D ** (c)) * (3 * d - 3 + c - 3 * c * (D ** c))
            # valCheck = Base * (A + B)
            # fx = (x**d)*(y**(t))*exp(-(x**c))
            # valCheck = sy.diff(fx,(D,T),3)
            # # valCheck = sp.diff(lambda x,y: (x**d)*(y**t)*exp(-(x**c)), (D,T), (3,0))
            # if i < 10:
            #   print val, valCheck
            #   i = i +1;
            # x,y = sy.symbols('x, y')
            # expr = x**d * y**t * sy.exp(-(x**c))
            # f_prime = expr.diff(x,3)
            # checkVal = f_prime.subs([(x,D),(y,T)])

            # print val, valCheck, checkVal*(D**3)

        else:
            val = (
                (D ** (d))
                * (T ** t)
                * np.exp(-(T ** m))
                * np.exp(-(D ** c))
                * (
                    d * (d - 1) * (d - 2)
                    + (D ** c)
                    * (
                        -2 * c
                        + 6 * d * c
                        - 3 * (d ** 2) * c
                        - 3 * d * (c ** 2)
                        + 3 * (c ** 2)
                        - (c ** 3)
                    )
                    + (D ** (2 * c)) * (3 * d * (c ** 2) - 3 * (c ** 2) + 3 * (c ** 3))
                    - (c ** 3) * (D ** (3 * c))
                )
            )
            # print val
        d3rd_vals.append(val)


# TODO
def d4rd(D, T):
    global d4rd_vals
    global coeffs
    d4rd_vals = []
    # i = 1
    for d, t, c, m in coeffs:
        val = 0
        if c == 0:
            val = d * (d - 1) * (d - 2) * (d - 3) * (D ** d) * (T ** t)
        elif m == 0:
            # val = (D**(d))*(T**t)*exp(-(D**c))*(d*(d-1)*(d-2) + (D**c)*(-2*c + 6*d*c-3*(d**2)*c - 3*d*(c**2) + 3*(c**2) - (c**3)) + (D**(2*c))*(3*d*(c**2) -3*(c**2) + 3*(c**3)) - (c**3)*(D**(3*c)))
            val = (
                (D ** (d))
                * (T ** (t))
                * np.exp(-(D ** c))
                * (
                    (d - c * (D ** c))
                    * (d - 1 - c * (D ** c))
                    * (d - 2 - c * (D ** c))
                    * (d - 3 - c * (D ** c))
                    + (
                        6 * (c ** 2) * D ** (2 * c)
                        - 7 * (c ** 2) * (D ** c)
                        + (c ** 2)
                        - 12 * c * d * (D ** c)
                        + 4 * c * d
                        + 18 * c * (D ** c)
                        - 6 * c
                        + 6 * d ** 2
                        - 18 * d
                        + 11
                    )
                    * (-(c ** 2) * (D ** (c)))
                )
            )
            # Base = (D ** (d)) * (T ** (t)) * np.exp(-D ** c)
            # A = (
            #     (d - c * (D ** c))
            #     * (d - 1 - c * (D ** c))
            #     * (d - 2 - c * (D ** c))
            #     * (d - 3 - c * (D ** c))
            # )
            # # B = (3*d-3+c-3*c*(D**c))*(d-3-c*(D**c));
            # # E = c*(3*d-3+c-6*c*(D**c))
            # # C = (d-c*(D**c))*(d-1-c*(D**c)) + (d-1-c*(D**c))*(d-2-c*(D**c)) + (d-c*(D**c))*(d-2-c*(D**c))

            # BCE = (
            #     (c ** 2) * ((D ** c) - 1) * (6 * D ** c - 1)
            #     - 2 * c * (2 * d - 3) * (3 * (D ** c) - 1)
            #     + 6 * d ** 2
            #     - 18 * d
            #     + 11
            # )
            # valCheck = Base*(A-(c**2)*(D**(c))*(B+ C+E));
            # valCheck1 = Base * (A - (c ** 2) * (D ** (c)) * (BCE))
            # if i < 10:
            #   print val, valCheck, valCheck1
            #   i = i +1;

            # x,y = sy.symbols('x, y')
            # expr = x**d * y**t * sy.exp(-(x**c))
            # f_prime = expr.diff(x,4)
            # checkVal = f_prime.subs([(x,D),(y,T)])

            # print val, valCheck1, checkVal*(D**4)

        else:  # TODO(
            val = (
                (D ** (d))
                * (T ** t)
                * np.exp(-(T ** m))
                * np.exp(-(D ** c))
                * (
                    d * (d - 1) * (d - 2)
                    + (D ** c)
                    * (
                        -2 * c
                        + 6 * d * c
                        - 3 * (d ** 2) * c
                        - 3 * d * (c ** 2)
                        + 3 * (c ** 2)
                        - (c ** 3)
                    )
                    + (D ** (2 * c)) * (3 * d * (c ** 2) - 3 * (c ** 2) + 3 * (c ** 3))
                    - (c ** 3) * (D ** (3 * c))
                )
            )
            val = 0
            # print val
        d4rd_vals.append(val)


def d5rd(D, T):
    global d5rd_vals
    global coeffs
    d5rd_vals = []
    for d, t, c, m in coeffs:
        val = 0
        if c == 0:
            val = d * (d - 1) * (d - 2) * (d - 3) * (d - 4) * (D ** d) * (T ** t)
        elif m == 0:
            # Base = (D**(d))*(T**(t))*exp(-(D**c))
            # A = (d-c*(D**c))*(d-1-c*(D**c))*(d-2-c*(D**c))*(d-3-c*(D**c))*(d-4-c*(D**c))

            # # from  d4rd
            # # B = (3*d-3+c-3*c*(D**c))*(d-3-c*(D**c));
            # # E = c*(3*d-3+c-6*c*(D**c))
            # # C = (d-c*(D**c))*(d-1-c*(D**c)) + (d-1-c*(D**c))*(d-2-c*(D**c)) + (d-c*(D**c))*(d-2-c*(D**c))
            # BCE = (c**2)*((D**c)-1)*(6*D**c-1) - 2*c*(2*d-3)*(3*(D**c)-1)+6*d**2-18*d + 11
            # # dr4dvalCheck = Base*(A-(c**2)*(D**(c))*(BCE));

            #                   #3                                          #0                                                  #2                                  #1
            # F = (d-c*(D**c))*(d-1-c*(D**c))*(d-2-c*(D**c)) + (d-1-c*(D**c))*(d-2-c*(D**c))*(d-3-c*(D**c)) + (d-c*(D**c))*(d-1-c*(D**c))*(d-3-c*(D**c)) + (d-c*(D**c))*(d-2-c*(D**c))*(d-3-c*(D**c))
            # G = (-c**2*(D**c))*(c**2*(D**(c)))*(12*c*(D**c) - 7*c -12*d + 18)
            # H = (-c**2*(D**c))*c*(BCE)

            # val = Base*(d-4-c*(D**c)*(A - (c**2)*(D**(c))*(BCE))) + Base*(-c**2*(D**c)*F) + Base*(G+H)

            x, y = sy.symbols("x, y")
            expr = x ** d * y ** t * sy.exp(-(x ** c))
            f_prime = expr.diff(x, 5)
            val = f_prime.subs([(x, D), (y, T)]) * (D ** 5)

            # print 'dr5d', val, checkVal*(D**5)

        else:  # TODO
            val = (
                (D ** (d))
                * (T ** t)
                * np.exp(-(T ** m))
                * np.exp(-(D ** c))
                * (
                    d * (d - 1) * (d - 2)
                    + (D ** c)
                    * (
                        -2 * c
                        + 6 * d * c
                        - 3 * (d ** 2) * c
                        - 3 * d * (c ** 2)
                        + 3 * (c ** 2)
                        - (c ** 3)
                    )
                    + (D ** (2 * c)) * (3 * d * (c ** 2) - 3 * (c ** 2) + 3 * (c ** 3))
                    - (c ** 3) * (D ** (3 * c))
                )
            )
            val = 0
            # print val
        d5rd_vals.append(val)



def dtrdt(D, T):
    global dtrdt_vals
    global coeffs
    dtrdt_vals = []
    for d, t, c, m in coeffs:
        val = 0
        if c == 0:
            val = d * t * (D ** d) * (T ** t)
        elif m == 0:
            val = np.exp(-D ** c) * ((D ** (d)) * (T ** t) * (d - c * (D ** c)))
        else:
            val = (
                np.exp(-D ** c)
                * ((D ** (d)) * (T ** t) * (d - c * (D ** c)))
                * np.exp(-T ** m)
                * (t - m * (T ** m))
            )
        dtrdt_vals.append(val)


def rTT(D, T):
    global rtt_vals
    global coeffs
    rtt_vals = []
    for d, t, c, m in coeffs:
        val = 0
        if c == 0:
            val = t * (t - 1) * (D ** d) * (T ** (t))
        elif m == 0:
            val = np.exp(-D ** c) * ((D ** (d)) * (T ** (t))) * t * (t - 1)
        else:  # NOT DONE
            val = np.exp(-D ** c) * (D ** (d)) * (T ** t) * np.exp(-T ** t)
        rtt_vals.append(val)


def d3rdRes(D, T, Y, Beta):
    global d3rd_vals
    global coeffs
    d3rd_vals = []
    val = 0
    for x, y in zip(Y, Beta):
        x = x - 1
        [d, t, c, m] = coeffs[x]
        if c == 0:
            val += y * d * (d - 1) * (d - 2) * (D ** d) * (T ** t)
        elif m == 0:  # TODO
            val += (
                y
                * np.exp(-D ** c)
                * (D ** (d))
                * (T ** t)
                * (
                    d * (d - 1) * (d - 2)
                    + c
                    * (D ** c)
                    * (-2 + 6 * d - 3 * (d ** 2) - 3 * d * c + 3 * c - (c ** 2))
                    + (D ** (2 * c)) * (3 * d * (c ** 2) - 3 * (c ** 2) + 3 * (c ** 3))
                    - (c ** 3) * (D ** (3 * c))
                )
            )
        else:
            val += (
                y
                * np.exp(-D ** c)
                * (D ** (d))
                * (T ** t)
                * np.exp(-T ** m)
                * (
                    d * (d - 1) * (d - 2)
                    + c
                    * (D ** c)
                    * (-2 + 6 * d - 3 * (d ** 2) - 3 * d * c + 3 * c - (c ** 2))
                    + (D ** (2 * c)) * (3 * d * (c ** 2) - 3 * (c ** 2) + 3 * (c ** 3))
                    - (c ** 3) * (D ** (3 * c))
                )
            )
        # print val
    return val


# USED ENGLE - 9/19

# PVT derivatives
def drdRes(D, T, Y, Beta):
    """
    Calculates the partial derivaties w.r.t. density
    Inputs:
        D - Delta
        T - Tau
        Y - index of basis Function (int or array)
        Beta - weighting (float or array)
    """
    global drd_vals
    global coeffs
    val = 0
    if type(Y) is not int:
        for x, y in zip(Y, Beta):
            x = x - 1
            [d, t, c, m] = coeffs[x]
            if c == 0:
                try:
                    val += y * d * (D ** d) * (T ** t)
                    pass
                except Exception:
                    print("Term", x, d, t, T, D)
            elif m == 0:
                val += y * np.exp(-D ** c) * ((D ** (d)) * (T ** t) * (d - c * (D ** c)))
            else:
                val += (
                    y
                    * np.exp(-D ** c)
                    * (D ** (d))
                    * (T ** t)
                    * np.exp(-T ** m)
                    * (d - c * (D ** c))
                )
        return val
    else:
        x = Y - 1
        y = Beta
        [d, t, c, m] = coeffs[x]
        if c == 0:
            val += y * d * (D ** d) * (T ** t)
        elif m == 0:
            val += y * np.exp(-D ** c) * ((D ** (d)) * (T ** t) * (d - c * (D ** c)))
        else:
            val += (
                y
                * np.exp(-D ** c)
                * (D ** (d))
                * (T ** t)
                * np.exp(-T ** m)
                * (d - c * (D ** c))
            )
        return val

# CV derivatives
def iTT(D, T):
    global itt_val

    # TOLUENE
    if molecule == "TOL":
        # a = [3.5241174832, 1.1360823464]
        c = [4.0, 0, 0]
        # B = [.96464, -2.7855, 0.86712, -0.18860, 0.11804, 0.0025181, 0.57196, -0.029287, -0.43351, -0.12540, -0.028207, 0.014076];
        v = [1.6994, 8.0577, 17.059, 8.4567, 8.6423]
        u = [190, 797, 1619, 3072, 7915]

    # Carbon Monoxide
    if molecule == "CO":
        # a = [-3.3728318564, 3.3683460039]
        c = [3.5, 0.22311e-6, 1.5]
        # B = [0.90554, -2.4515, 0.53149, 0.024173, 0.072156, 0.00018818, 0.19405, -0.043268, -0.12778, -0.027896, -0.03414, 0.016329]
        v = [1.0128]
        u = [3089]

    # Carbon Dioxide
    if molecule == "CO2":
        # a = [8.37304456, -3.70454304, 2.5]
        c = [2.5 + 1, 0, 0]
        v = [1.99427042, 0.62105248, 0.41195293, 1.04028922, 0.08327678]
        u = [
            3.15163 * float(critT),
            6.11190 * float(critT),
            6.77708 * float(critT),
            11.32384 * float(critT),
            27.0 * float(critT),
        ]

    # Water
    if molecule == "H2O":
        # a = [-8.32044648201, 6.6832105268]
        c = [3.00632 + 1, 0, 0]
        v = [0.012436, 0.97315, 1.27950, 0.96956, 0.24873]
        u = [
            1.2878967 * float(critT),
            3.53734222 * float(critT),
            7.74073708 * float(critT),
            9.24437796 * float(critT),
            27.5075105 * float(critT),
        ]

    itt_val = 0
    if c[2] != 0:
        itt_val = -(c[0] - 1) - c[1] * (critT ** c[2]) / (c[2] * (c[2] + 1)) * (
            T ** (-c[2])
        ) * (-c[2]) * (-c[2] - 1)
        # /T**2;
    else:
        itt_val = -(c[0] - 1)
    for uit, vi in zip(u, v):
        ui = uit / float(critT)
        itt_val = itt_val - vi * ((ui) ** 2) * (T ** 2) * np.exp(-ui * T) * (
            (1 - np.exp(-ui * T)) ** (-2)
        )
    return itt_val


def rTTRes(D, T, Y, Beta):
    global rtt_vals
    global coeffs
    val = 0
    if type(Y) is int:
        x = Y - 1
        y = Beta
        [d, t, c, m] = coeffs[x]
        if c == 0:
            value = t * (t - 1) * (D ** d) * (T ** (t))
        elif m == 0:  
            value = np.exp(-D ** c) * ((D ** (d)) * (T ** (t))) * t * (t - 1)
        else:
            value = (
                np.exp(-D ** c)
                * ((D ** (d)) * (T ** (t)))
                * ((t - m * (T ** m)) * (t - 1 - m * (T ** m)) - (m ** 2) * (T ** m))
            )
        val = y * value
    else:
        for x, y in zip(Y, Beta):
            x = x - 1
            [d, t, c, m] = coeffs[x]
            if c == 0:
                value = t * (t - 1) * (D ** d) * (T ** (t))
            elif m == 0: 
                value = np.exp(-D ** c) * ((D ** (d)) * (T ** (t))) * t * (t - 1)
            else:
                value = (
                    np.exp(-D ** c)
                    * ((D ** (d)) * (T ** (t)))
                    * (
                        (t - m * (T ** m)) * (t - 1 - m * (T ** m))
                        - (m ** 2) * (T ** m)
                    )
                )
            val += y * value
    return val


def d2rdRes(D, T, Y, Beta):
    global d2rd_vals
    global coeffs
    val = 0
    d2rd_vals = []
    if type(Y) is int:
        x = Y - 1
        y = Beta
        [d, t, c, m] = coeffs[x]
        if c == 0:
            val = y * d * (d - 1) * (D ** d) * (T ** t)
        elif m != 0:  # TODO
            val = np.exp(-D ** c) * (D ** (d)) * (T ** t) * np.exp(-T ** t)
        else:
            val = (
                y
                * np.exp(-D ** c)
                * (
                    (D ** (d))
                    * (T ** t)
                    * (
                        (d - c * (D ** c)) * (d - 1 - c * (D ** c))
                        - (c ** 2) * (D ** c)
                    )
                )
            )
    else:
        for x, y in zip(Y, Beta):
            x = x - 1
            [d, t, c, m] = coeffs[x]
            if c == 0:
                val += y * d * (d - 1) * (D ** d) * (T ** t)
            elif m == 0:  # TODO
                val += (
                    y
                    * np.exp(-D ** c)
                    * (
                        (D ** (d))
                        * (T ** t)
                        * (
                            (d - c * (D ** c)) * (d - 1 - c * (D ** c))
                            - (c ** 2) * (D ** c)
                        )
                    )
                )
            else:
                val += (
                    y
                    * np.exp(-D ** c)
                    * (
                        (D ** (d))
                        * (T ** t)
                        * np.exp(-T ** m)
                        * (
                            (d - c * (D ** c)) * (d - 1 - c * (D ** c))
                            - (c ** 2) * (D ** c)
                        )
                    )
                )
    return val


def dtrdtRes(D, T, Y, Beta):
    global dtrdt_vals
    global coeffs
    val = 0
    dtrdt_vals = []
    if type(Y) is int:
        x = Y - 1
        y = Beta
        [d, t, c, m] = coeffs[x]
        if c == 0:
            val = y * d * t * (D ** d) * (T ** t)
        elif m == 0:  # TODO
            val = y * t * np.exp(-D ** c) * ((D ** (d)) * (T ** t) * (d - c * (D ** c)))
        else:
            val = (
                y
                * t
                * np.exp(-D ** c)
                * ((D ** (d)) * (T ** t) * np.exp(-T ** m) * (d - c * (D ** c)))
                * (t - m * (T ** m))
            )
    else:
        for x, y in zip(Y, Beta):
            x = x - 1
            [d, t, c, m] = coeffs[x]
            if c == 0:
                val += y * d * t * (D ** d) * (T ** t)
            elif m == 0:
                val += (
                    y * t * np.exp(-D ** c) * ((D ** (d)) * (T ** t) * (d - c * (D ** c)))
                )
            else:
                val += (
                    y
                    * t
                    * np.exp(-D ** c)
                    * ((D ** (d)) * (T ** t) * np.exp(-T ** m) * (d - c * (D ** c)))
                    * (t - m * (T ** m))
                )
    return val


def d2rdrtRes(D, T, Y, Beta):
    global d2rd_vals
    global coeffs
    val = 0
    d2rd_vals = []
    if type(Y) is int:
        x = Y - 1
        y = Beta
        [d, t, c, m] = coeffs[x]
        if c == 0:
            val = y * d * (d - 1) * (D ** d) * t * (T ** t)
        elif m != 0:  # TODO
            val = (
                np.exp(-D ** c)
                * (D ** (d))
                * (T ** t)
                * np.exp(-T ** t)
                * (t - m * (T ** m))
            )
        else:
            val = (
                y
                * np.exp(-D ** c)
                * (
                    (D ** (d))
                    * (T ** t)
                    * (
                        (d - c * (D ** c)) * (d - 1 - c * (D ** c))
                        - (c ** 2) * (D ** c)
                    )
                )
                * t
            )
    else:
        for x, y in zip(Y, Beta):
            x = x - 1
            [d, t, c, m] = coeffs[x]
            if c == 0:
                val += y * d * (d - 1) * (D ** d) * t * (T ** t)
            elif m == 0:  # TODO
                val += (
                    y
                    * np.exp(-D ** c)
                    * (
                        (D ** (d))
                        * (T ** t)
                        * (
                            (d - c * (D ** c)) * (d - 1 - c * (D ** c))
                            - (c ** 2) * (D ** c)
                        )
                    )
                    * t
                )
            else:
                val += (
                    y
                    * np.exp(-D ** c)
                    * (
                        (D ** (d))
                        * (T ** t)
                        * np.exp(-T ** m)
                        * (
                            (d - c * (D ** c)) * (d - 1 - c * (D ** c))
                            - (c ** 2) * (D ** c)
                        )
                    )
                    * (t - m * (T ** m))
                )
    return val