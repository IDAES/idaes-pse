import { IRenderMime } from '@jupyterlab/rendermime-interfaces';


import { JSONObject } from '@phosphor/coreutils';


import { Widget } from '@phosphor/widgets';

import 'jquery';
import 'lodash';
import 'backbone';
import {dia} from 'jointjs';
import {shapes} from 'jointjs';

import '../style/index.css';
import '../style/joint.css';

import '../style/iconmapping.css';

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

    new dia.Paper({
                    el: this.myholder,
                    model: graph,
                    width: 1000,
                    height: 1000,
                    gridSize: 1
                });
                
    console.log("after paper");
    graph.fromJSON(data);

	var rect = new shapes.standard.Rectangle();
    rect.position(100, 30);
    rect.resize(100, 40);
    rect.attr({
        body: {
            fill: 'blue'
        },
        label: {
            text: 'Hello',
            fill: 'white'
        }
    });
	rect.addTo(graph);

    const image = document.createElement('div');
    image.className = 'mixer';

	var testimage = new shapes.standard.Image();
    testimage.position(100, 150);

    testimage.resize(100, 100);
    testimage.attr('image/xlinkHref', image);
    //testimage.markup = [{"tagName": "image", "selector": "image"}];
    testimage.addTo(graph);


    console.log("after graph load");

    function saveFile(evt: MouseEvent) {
        var jsonstring = JSON.stringify(graph.toJSON());
        var data = "text/json;charset=utf-8," + encodeURIComponent(jsonstring);
        var fakeobj = document.createElement('a');
        fakeobj.href = "data:" + data;
        fakeobj.download = 'serialized_graph.idaes.vis';
        fakeobj.click();
        fakeobj.remove();
    }

    this.save_file_button.onclick = saveFile;

    return Promise.resolve();
  }

  private _mimeType: string;
  private save_file_button: HTMLButtonElement;
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
