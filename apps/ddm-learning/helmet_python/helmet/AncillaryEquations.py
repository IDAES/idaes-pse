# module for ALAMO Writing.

# import math
import sys
import alamopy
import importlib

from helmet import DataImport
from helmet import DataManipulation


molecule = ""
max_time = 1000



global DLfun, DVfun, PVfun
DLfun, DVfun, PVfun = None, None, None

DataToWrite = []


def DL():
    global DLfun
    DataImport.DL(molecule)
    Values = DataImport.DLValues

    xval, zval = [], []
    for x in Values:
        manVal = DataManipulation.DL(x)
        xval.append(manVal[0])
        zval.append(manVal[1])

    res = alamopy.alamo(
        xval,
        zval,
        zlabels=["%sDL" % molecule],
        almname="%sDL" % molecule,
        monomialpower=(-1, 1, 2, 3),
        expfcns=1,
        savetrace=True,
        savescratch=True,
    )
    # print(res['model'], res['ssr'], res['R2'])

    # conf_inv = alamopy.almconfidence(res)
    # print('Confidence Intervals : {}'.format(conf_inv['conf_inv']))

    # alamopy.almplot(res)
    if res is None:
        raise Exception("Model does not exist for Saturated Liquid Density")

    DLfun = importlib.import_module("%sDL" % molecule, "f")


def getDL():
    global DLfun
    if DLfun is None:
        DLfun = importlib.import_module("%sDL" % molecule, "f")
    return DLfun


def DV():
    global DVfun
    DataImport.DV(molecule)
    Values = DataImport.DVValues

    xval, zval = [], []
    for x in Values:
        manVal = DataManipulation.DV(x)
        xval.append(manVal[0])
        zval.append(manVal[1])

    res = alamopy.alamo(
        xval,
        zval,
        zlabels=["%sDV" % molecule],
        almname="%sDV" % molecule,
        monomialpower=(1, 2, 3, 4, 5, 6),
        savetrace=True,
    )
    # print(res['model'], res['ssr'], res['R2'])

    if res is None:
        raise Exception("Model does not exist for Saturated Vapor Density")

    DVfun = importlib.import_module("%sDV" % molecule, "f")


def getDV():
    global DVfun
    if DVfun is None:
        DVfun = importlib.import_module("%sDV" % molecule, "f")
    return DVfun


def PV():
    global PVfun
    DataImport.PV(molecule)
    Values = DataImport.Values

    xval, zval = [], []
    for x in Values:
        manVal = DataManipulation.PV(x)
        xval.append([manVal[0], manVal[1]])
        zval.append(manVal[2])

    res = alamopy.alamo(
        xval,
        zval,
        zlabels=["%sPV" % molecule],
        almname="%sPV" % molecule,
        monomialpower=(1, 2, 3, 4, 5, 6),
    )
    # print(res['model'], res['ssr'], res['R2'])

    if res is None:
        raise Exception("Model does not exist for Vapor Pressure")

    PVfun = importlib.import_module("%sPV" % molecule, "f")


def getPV():
    global PVfun
    if PVfun is None:
        PVfun = importlib.import_module("%sPV" % molecule, "f")
    return PVfun




# Opening an ALAMO File
# def openFile(data_name):
#     global textFile
#     print("Creating new text file")
#     name = molecule + data_name + ".alm"
#     print(name)
#     try:
#         textFile = open(name, "w")
#     except Exception:
#         print("couldn't open the file Something wrong")
#         sys.exit(0)


# # Closing an ALAMO File
# def closeFile():
#     textFile.close()
#     print("Closing File")


# # Set of Basis Functions used
# def CP0CustomBasisMin(num_custom):
#     textFile.write("ncustombas %d\n" % num_custom)
#     textFile.write("BEGIN_CUSTOMBAS\n")
#     for i in range(1, 8000):
#         textFile.write(
#             "(min(30,%f)/T)^2 * exp(min(30,%f)/T)/(exp(min(30,%f)/T)-1)^2\n"
#             % (i / 1000.0, i / 1000.0, i / 1000.0)
#         )
#     # 		textFile.write("(T%d\n" %(i, i, i) );
#     textFile.write("END_CUSTOMBAS\n\n")


# def CP0CustomBasisArr(num_custom):
#     textFile.write("ncustombas %d\n" % num_custom)
#     textFile.write("BEGIN_CUSTOMBAS\n")
#     for i in range(1, 8000):  # 1-5
#         textFile.write(
#             "exp(%f/T)*(%f/T/(exp(%f/T)-1))^2\n" % (i / 1000.0, i / 1000.0, i / 1000.0)
#         )
#     # 		textFile.write("(T%d\n" %(i, i, i) );
#     textFile.write("END_CUSTOMBAS\n\n")

# def CP0CustomBasis(num_custom):
#     textFile.write("ncustombas %d\n" % num_custom)
#     textFile.write("BEGIN_CUSTOMBAS\n")
#     for i in range(1, 8000):  # 1-5
#         textFile.write(
#             "(%f/T)^2 * exp(%f/T)/(exp(%f/T)-1)^2\n"
#             % (i / 1000.0, i / 1000.0, i / 1000.0)
#         )
#         # 		textFile.write("(T%d\n" %(i, i, i) );
#         # try:
#         #     calc = (i / 5) ** 2 * math.exp(i / 5) / (math.exp(i / 5) - 1) ** 2
#         #     # print l
#         # except OverflowError:
#         #     # print "OverFlow %f" %(i)
#         #     continue
#         # except ZeroDivisionError:
#         #     # print "Zero %f" %(i)
#         #     continue

#     textFile.write("END_CUSTOMBAS\n\n")




# def CP0():
#     DataImport.CP0(molecule)
#     CP0Values = DataImport.CP0Values
#     openFile("CP0")
#     mins = [0, 0]
#     maxs = [0, 0]
#     DataToWrite = []
#     for x in CP0Values:
#         manVal = DataManipulation.CP0(x)
#         DataToWrite.append(manVal)
#     mins[0] = min(y[0] for y in DataToWrite)
#     mins[1] = min(y[1] for y in DataToWrite)
#     maxs[0] = max(y[0] for y in DataToWrite)
#     maxs[1] = max(y[1] for y in DataToWrite)
#     minmax = "xmin %d \nxmax %d \n" % (mins[0], maxs[0])
#     writeHeader(1, 1, len(CP0Values), len(CP0Values), ["T"], ["CP0"], minmax)
#     CP0CustomBasis(7999)
#     # CP0CustomBasisArr(7999);
#     textFile.write("BEGIN_DATA\n")
#     for x in DataToWrite:
#         textFile.write("	%s   %s\n" % (x[0], x[1]))
#     textFile.write("END_DATA\n\n")
#     closeFile()

