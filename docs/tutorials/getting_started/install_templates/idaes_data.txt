.. note::
    If you want to install multiple versions of IDAES in different environments, you can set the ``IDAES_DATA`` environment variable before installing
    or using IDAES.  The easiest way to do this is to set the environment variable as part of activating your environment before installing IDAES.  
    Setting this environment variable will keep binary extensions and configuration files separate for the different IDAES versions. For more information
    about setting environment variables in conda environments see https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#saving-environment-variables. Using virtual environments (venv), the easiest way to set environment variables is to edit the ``activate`` script.
    You can edit the deactivate section to unset the variable when deactivating the environment as well.