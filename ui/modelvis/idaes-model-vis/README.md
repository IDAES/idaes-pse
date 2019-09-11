This is a ui component to the IDAES project. The FlowsheetSerializer module will serialize a flowsheet into a .idaes.vis file which this extension will render in JupyterLab. 

## This will only work in JupyterLab. It will NOT work in a Jupyter Notebook.

## Install

Activate the environment that you have idaes installed in. 

Run the following commands from the top of your git clone.

```
conda install -c conda-forge nodejs jupyterlab
cd ui/modelvis/idaes-model-vis
npm install
npm run build
jupyter labextension link .
```

## Running

This can be run from any directory. NOTE: JupyterLab will use your current working directory as your home directory.

`jupyter lab`

## To use

Launch a Jupyter notebook or a python file from inside of Jupyter Lab.

First create your model then run the following commands to serialize your model.
```
from idaes.dmf.ui.flowsheet_serializer import FlowsheetSerializer
flowsheet_serializer = FlowsheetSerializer()
flowsheet_serializer.save(flowsheet, "file_prefix")
```
This will create a file called file_prefix.idaes.vis in the directory that you are in in Jupyter Lab.

Double click the created file in the file tree. This should open the file and you will be able to drag and drop elements.

Anytime your mouse leaves the graph area the layout will be saved.

NOTE: If you call flowsheet_serializer.save again, it will not re-save the file unless you specify `overwrite=True`. This is done in an effort to preserve the layout that was saved in the idaes vis tab. If you specify `overwrite=True` then the layout will be lost and the model will be rendered using the default layout.

Example: `flowsheet_serializer.save(flowsheet, "file_prefix", overwrite=True)`
