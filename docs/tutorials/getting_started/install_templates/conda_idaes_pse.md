We recommend using Conda to manage your environment & modules.

.. code-block:: console
   
   # Create a new Conda environment (Here named my-idaes-env, but this is up to you)
   conda create --yes --name my-idaes-env python=3.8
   conda activate my-idaes-env

   # Install IDAES Conda package
   conda install --yes -c IDAES-PSE -c conda-forge idaes-pse