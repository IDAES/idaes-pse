import { IRenderMime } from '@jupyterlab/rendermime-interfaces';


import { JSONObject } from '@phosphor/coreutils';


import { Widget } from '@phosphor/widgets';

import 'jquery';
import 'lodash';
import 'backbone';
import {dia} from 'jointjs';

import '../style/index.css';
import '../style/joint.css';

/**
 * The default mime type for the extension.
 */
const MIME_TYPE = 'application/vnd.idaes.model';

/**
 * The class name added to the extension.
 */
const CLASS_NAME = 'mimerenderer-mimerenderer-idaes-model';

/**
 * A widget for rendering mimerenderer-idaes-model.
 */
export class OutputWidget extends Widget implements IRenderMime.IRenderer {
  /**
   * Construct a new output widget.
   */
  constructor(options: IRenderMime.IRendererOptions) {
    super();
    this._mimeType = options.mimeType;
    this.addClass(CLASS_NAME);

    this.save_file_button = document.createElement("button");
    this.save_file_button.id = "save_file_button";
    this.save_file_button.innerText = "Save Graph to File";

    ////this.load_file_button = document.createElement("button");
    //this.load_file_button.id = "load_file_button";
    //this.load_file_button.name = "load file";
    //this.load_file_button.innerText = "Load File";

    this.myholder = document.createElement("div");
    this.myholder.id = "myholder";
  }

  /**
   * Render mimerenderer-idaes-model into this widget's node.
   */
  renderModel(model: IRenderMime.IMimeModel): Promise<void> {
    
    let data = model.data[this._mimeType] as JSONObject;
    //this.node.textContent = JSON.stringify(data);


    this.node.appendChild(this.save_file_button);
    //this.node.appendChild(this.load_file_button);

    this.node.appendChild(this.myholder);

    var somelink = document.createElement('a');
    somelink.href = "http://www.google.com";
    somelink.innerHTML = "something";
    this.myholder.appendChild(somelink);

    var graph = new dia.Graph;

    console.log("before paper");
    //var paper = new dia.Paper({
    new dia.Paper({
                    el: document.getElementById('myholder'),
                    model: graph,
                    width: 1000,
                    height: 1000,
                    gridSize: 1
                });
                
                //console.log(paper);

    console.log("after paper");
    graph.fromJSON(data);

    console.log("after graph load");
    /*
    <body>
        <button onclick="saveFile(this)" id="save_file_button">Save Graph to File</button><br>
            //Load Graph From File: <input type="file" id="load_file_button" name="load file"/><br>
        A sample saved graph is available at idaes/ui/html/demo_graph.json
        <!-- content -->
        <div id="myholder"></div>

        <!-- dependencies -->
        <script type="module">
              
        
    </script>

    </body>
    */ 

    function saveFile(evt: MouseEvent) {
        var jsonstring = JSON.stringify(graph.toJSON());
        var data = "text/json;charset=utf-8," + encodeURIComponent(jsonstring);
        var fakeobj = document.createElement('a');
        fakeobj.href = "data:" + data;
        fakeobj.download = 'serialized_graph.idaes.vis';
        fakeobj.click();
        fakeobj.remove();
    }

    /*
    function loadFile(evt: InputEvent) {
        // getting a hold of the file reference
        var file = evt.target.files[0]; 
        var reader = new FileReader();
        reader.readAsText(file);

        // here we tell the reader what to do when it's done reading...
        reader.onload = readerEvent => {
                  var content = readerEvent.target.result;
                  graph.fromJSON(JSON.parse(content));
               }

    }
    */

    //this.load_file_button.addEventListener('change', loadFile, false);
    this.save_file_button.onclick = saveFile;

    return Promise.resolve();
  }

  private _mimeType: string;
  private save_file_button: HTMLButtonElement;
  //private load_file_button: HTMLButtonElement;
  private myholder: HTMLDivElement;
}

/**
 * A mime renderer factory for mimerenderer-idaes-model data.
 */
export const rendererFactory: IRenderMime.IRendererFactory = {
  safe: true,
  mimeTypes: [MIME_TYPE],
  createRenderer: options => new OutputWidget(options)
};

/**
 * Extension definition.
 */
const extension: IRenderMime.IExtension = {
  id: 'idaes-model-vis:plugin',
  rendererFactory,
  rank: 0,
  dataType: 'json',
  fileTypes: [
    {
      name: 'mimerenderer-idaes-model',
      mimeTypes: [MIME_TYPE],
      extensions: ['.idaes.vis']
    }
  ],
  documentWidgetFactoryOptions: {
    name: 'IDAES Model ',
    primaryFileType: 'mimerenderer-idaes-model',
    fileTypes: ['mimerenderer-idaes-model'],
    defaultFor: ['mimerenderer-idaes-model']
  }
};

export default extension;
