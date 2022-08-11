+---------------------------**** DEVELOPER WARNING ****----------------------------+
| OS installation guides include direct, numbered line references to this document |
| Modifications to this file will require updates to the Quickstart guide for      |
| Mac, Windows, & Linux                                                            |
+--------------------------**** modification warning ****--------------------------+

# Set up & activate Conda new environment with IDAES-PSE
conda create --yes --name my-idaes-env -c conda-forge -c IDAES-PSE python=3.10 idaes-pse
conda activate my-idaes-env

# Install IDAES Extensions
.. Windows&Linux
idaes get-extensions
.. MacOSX
conda install -c conda-forge ipopt

# Install IDAES Examples (optional)
idaes get-examples

# Run Tests
pytest --pyargs idaes -W ignore