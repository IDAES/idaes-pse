import os
from idaes.surrogate.helmet import Helmet


def main():
    # Push information for the molecule

    num_terms = 14

    Fluids = {'H2O': (647.096, 22.064, 17.8737279956, 18.015268, 273.16, 0.344 )}

    molecule = 'H2O'
    (critT, critP, critD, M, triple, acc) = Fluids[molecule]
    R = 8.314472; # J mol^-1 K^-1

    # Constants for a molecule
    Helmet.initialize(molecule=molecule, fluid_data = Fluids[molecule], filename = os.getcwd() + "/%s"%molecule, gamsname= os.getcwd() + "/%s"%molecule, props=['PVT','CV', 'CP','SND'], sample =3)

    # Prepare Ancillary Equations of sat liq/vapor density and vapor pressure
    Helmet.prepareAncillaryEquations()  # plot=True

    # View data used for regression
    Helmet.viewPropertyData()

    # Write GAMS Gdx data file and regression file
    Helmet.setupRegression(numTerms = num_terms, gams=True)

    # runs the GAMS gdx and regression file
    # Helmet.runRegression(gams=True)


    # View Results by importing the data
    #Helmet.viewResults("H2Omain.lst", plot=True)


if __name__ == "__main__":
    main()
