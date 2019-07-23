import { IRenderMime } from '@jupyterlab/rendermime-interfaces';
import { JSONObject } from '@phosphor/coreutils';
import { Widget } from '@phosphor/widgets';
import * as joint from 'jointjs';
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

namespace utils {
  type map = {
    [key: string]: string;
  }
  const staticURLToFilepath: map = {};
  // Takes filepath like 'icons/compressor_2_flipped.svg'
  // Returns a CSS var that will hold the built static url for the icon
  // Eg., '--idaes-compressor-2-flipped-icon'
  // These are defined in ../style/index.css
  function getCSSVar(filepath: string) {
    const filename = filepath.replace(/^icons\/|.svg$/gi, '');
    return `--idaes-${filename.replace(/_/g, '-')}-icon`;
  }

  // Remaps the image hrefs in a .vis file provided as JSON in 1 of 2 directions
  // forDisplay: filepath => staticURL
  // forStorage: staticURL => filepath
  export function remapIcons(data: JSONObject, mode: 'forDisplay' | 'forStorage') {
    const rootCSS = getComputedStyle(document.documentElement);
    const cells = data.cells as any[];
    cells.filter(c => c.type === 'standard.Image').forEach(cell => {
      const image = cell.attrs.image;
      switch(mode) {
        case 'forDisplay': {
          const filepath = image.xlinkHref
          const cssVar = getCSSVar(filepath);
          // Ref: https://stackoverflow.com/a/19826476
          const staticURL = ((
            rootCSS.getPropertyValue(cssVar).match(/url\([^\)]+\)/gi) || ['']
          )[0].split(/[()'"]+/)[1] || '').replace(/\\/g, '');
          if (staticURL) {
            image.xlinkHref = staticURL;
            // Needed to convert the other way later
            staticURLToFilepath[staticURL] = filepath;
          }
          break;
        }
        case 'forStorage': {
          // Assume: data has been displayed, populating staticURLToFilepath
          const staticURL = image.xlinkHref;
          const filepath = staticURLToFilepath[staticURL];
          image.xlinkHref = filepath;
          break;
        }
        default:
          break;
      }
    });
  }
}

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
  async renderModel(model: IRenderMime.IMimeModel) {

    let data = model.data[this._mimeType] as JSONObject;
    //this.node.textContent = JSON.stringify(data);
    this.node.appendChild(this.save_file_button);
    this.node.appendChild(this.myholder);

    // var somelink = document.createElement('a');
    // somelink.href = "http://www.google.com";
    // somelink.innerHTML = "something";
    // this.myholder.appendChild(somelink);

    var standard = joint.shapes.standard;
    var graph = new joint.dia.Graph([], { cellNamespace: { standard } });
    new joint.dia.Paper({
      el: this.myholder,
      model: graph,
      cellViewNamespace: { standard },
      width: 1000,
      height: 1000,
      gridSize: 1,
    });
    utils.remapIcons(data, 'forDisplay');
    graph.fromJSON(data);

    function saveFile(evt: MouseEvent) {
      let json = graph.toJSON();
      utils.remapIcons(json, 'forStorage');
      var data = "text/json;charset=utf-8," + encodeURIComponent(JSON.stringify(json));
      var fakeobj = document.createElement('a');
      fakeobj.href = "data:" + data;
      fakeobj.download = 'serialized_graph.idaes.vis';
      fakeobj.click();
      fakeobj.remove();
    }
    this.save_file_button.onclick = saveFile;
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
