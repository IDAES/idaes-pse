Natural Gas Property Package
============================
The file named "natural_gas_PR" is a parameter database for natural gas properties using the GenericParameterBlock. It contains relevant information and functions for working with natural gas components and reactions.

1. Introduction:
The "natural_gas_PR" file serves as a parameter database specifically designed for natural gas properties. It utilizes the GenericParameterBlock to store and access information related to natural gas components and their respective properties. This documentation provides an overview of the file structure, the available functions, supported components, included gas reactions, and the sources of property data.

2. File Structure:
The "natural_gas_PR" file consists of the following elements:

- get_prop Function: This function generates an instance's configuration to pass to the GenericParameterBlock. It allows users to retrieve the properties of a specific natural gas component.

- get_rxn Function: This function is used to configure the gas reactions using the GenericReactionParameterBlock. It allows users to specify the reactions and their associated parameters.

- Available Components: The file supports various natural gas components that can be selected when using the get_prop function. The available components are as follows:
  - H2
  - CO
  - H2O
  - CO2
  - O2
  - N2
  - Ar
  - CH4
  - C2H6
  - C3H8
  - C4H10
  - H2S
  - SO2
  - C2H4

- Gas Reactions: The "natural_gas_PR" file includes natural gas-related reactions within the get_rxn function. These reactions can be configured using the GenericReactionParameterBlock. The following reactions are included:
  - CH4 combustion
  - C2H6 combustion
  - C3H8 combustion
  - C4H10 combustion
  - CO combustion
  - H2 combustion
  - H2S combustion
  - C2H4 combustion
  - Water gas shift reaction

3. Property Sources:
The properties used in the "natural_gas_PR" file are obtained from the following sources:

- Ideal gas assumptions:
  - Source: NIST webbook
  - Properties: Heat capacity coefficients for all species except ethane, propane, and butane. Reference enthalpies and entropies for all species.

- Peng-Robinson Equation of State assumptions:
  - Source: The Properties of Gases and Liquids (1987) 4th edition, Chemical Engineering Series - Robert C. Reid
  - Properties: Critical temperatures and pressures. Omega. Heat capacity coefficients for ethane, propane, and butane.

4. Conclusion:
The "natural_gas_PR" file provides a parameter database for natural gas properties using the GenericParameterBlock. It offers functions to retrieve properties of specific natural gas components, configure gas reactions using the GenericReactionParameterBlock, and includes predefined gas reactions for modeling purposes. The properties used in the file are based on ideal gas assumptions from the NIST webbook and Peng-Robinson Equation of State assumptions from "The Properties of Gases and Liquids" by Robert C. Reid.
 




