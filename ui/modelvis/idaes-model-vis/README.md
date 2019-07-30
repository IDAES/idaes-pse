This is a ui component to the IDAES project. The FlowsheetSerializer module will serialize a flowsheet into a .idaes.vis file which this extension will render in JupyterLab.

## Install

From the top of your git clone

Follow all instructions to import and install idaes
```
conda install -c conda-forge nodejs jupyterlab
cd ui/modelvis/idaes-model-vis
npm install
npm run build
jupyter labextension link .
```

## Running

This can be run from whatever directory you wish. JupyterLab will use this as your home directory.

`jupyter lab`

## To use

Inside of a jupyter notebook or a python file

Create your model
```
from idaes.dmf.ui.flowsheet_serializer import FlowsheetSerializer
flowsheet_serializer = FlowsheetSerializer()
flowsheet_serializer.save(flowsheet, "file_prefix")
```
This will create a file called file_prefix.idaes.vis in the directory that you are in.

Double click the created file in the file tree. This should open the file and you will be able to drag and drop elements.

Once you have your model in a layout you want to save then click on the `Save Graph to File` button. It will save according to your browser's download settings, which defaults to your downloads folder as `serialized_graph.idaes.vis`.

You can then open that file by double clicking on it in the file tree.

