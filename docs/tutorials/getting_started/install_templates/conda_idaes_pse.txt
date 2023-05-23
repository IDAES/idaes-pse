We recommend using Conda to manage your environment & modules.

.. code-block:: console
   
   # Create a new Conda environment (Here named my-idaes-env, but this is up to you)
   conda create --yes --name my-idaes-env python=3.10
   conda activate my-idaes-env

   # Install IDAES Conda package
   conda install --yes -c IDAES-PSE -c conda-forge idaes-pse

.. note:: The command above will install the most recent stable (release) version of IDAES.
          To install other versions of IDAES, including pre-release versions, 
          refer to the :doc:`Versioned IDAES Installation <../../how_to_guides/versioned_idaes_install>` guide
