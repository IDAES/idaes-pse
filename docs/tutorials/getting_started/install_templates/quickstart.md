# Set up & activate Conda new environment with IDAES-PSE
conda create --yes --name my-idaes-env -c conda-forge -c IDAES-PSE python=3.10 idaes-pse
conda activate my-idaes-env

# Install IDAES Extensions
idaes get-extensions --extra petsc

# Install IDAES Examples (optional)
idaes get-examples

# Run Tests
pytest --pyargs idaes -W ignore