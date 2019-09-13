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
const CLASS_NAME = 'mimerenderer-idaes-model';

namespace utils {
  type map = {
    [key: string]: string;
  }
  const staticURLToFilepath: map = {};
  // Takes filepath like 'icons/compressor_2.svg'
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

    this.holder = document.createElement("div");
    this.holder.id = "mimerenderer-idaes-model";

    // We need to create the graph and paper in the constructor
    // If you try to create them in renderModel (which is called everytime the user
    // opens an .idaes.vis file or changes tabs) then you get icons that do not drop
    // when the mouseup event is emitted
    var standard = joint.shapes.standard;
    var width = 10000;
    var height = 10000;
    var gridSize = 1;
    this.graph = new joint.dia.Graph([], { cellNamespace: { standard } });
    this.paper = new joint.dia.Paper({
      el: this.holder,
      model: this.graph,
      cellViewNamespace: { standard },
      width: width,
      height: height,
      gridSize: gridSize,
      interactive: true
    });

    // Adds link tools (adding vertices, moving segments) to links when your mouse over
    this.paper.on("cell:mouseover", function(cellView, evt) {
      if (cellView.model.isLink()) {
        var verticesTool = new joint.linkTools.Vertices({
          focusOpacity: 0.5,
          redundancyRemoval: true,
          snapRadius: 20,
          vertexAdding: true,
        });
        var segmentsTool = new joint.linkTools.Segments();

        var toolsView = new joint.dia.ToolsView({
          tools: [
            verticesTool, segmentsTool
          ]
        });
        cellView.addTools(toolsView)
        cellView.showTools()
      }
    })

    // Removes the link tools when you leave the link
    this.paper.on("cell:mouseout", function(cellView, evt) {
      if (cellView.model.isLink()) {
        cellView.hideTools()
      }
    })

    // Icons rotate 90 degrees on right click. Replaces browser context menu
    this.paper.on("element:contextmenu", function(cellView, evt) {
      cellView.model.rotate(90)
    })

    // Constrain the elements to the paper
    this.paper.on('cell:pointermove', function (cellView, evt, x, y) {

      var bbox = cellView.getBBox();
      var constrained = false;

      var constrainedX = x;

      if (bbox.x <= 0) { constrainedX = x + gridSize; constrained = true }
      if (bbox.x + bbox.width >= width) { constrainedX = x - gridSize; constrained = true }

      var constrainedY = y;

      if (bbox.y <= 0) {  constrainedY = y + gridSize; constrained = true }
      if (bbox.y + bbox.height >= height) { constrainedY = y - gridSize; constrained = true }

      //if you fire the event all the time you get a stack overflow
      if (constrained) { cellView.pointermove(evt, constrainedX, constrainedY) }
    });
  }

  /**
   * Render mimerenderer-idaes-model into this widget's node.
   */
  async renderModel(model: IRenderMime.IMimeModel) {

    let data = model.data[this._mimeType] as JSONObject;
    this.node.appendChild(this.holder);

    utils.remapIcons(data, 'forDisplay');
    this.graph.fromJSON(data);

    // We need to remove the mouseleave events from the paper or every time that renderModel is called
    // it will create a new event listener causing the file to be saved multiple times
    this.paper.off('paper:mouseleave');
    // When the mouse leaves the paper save the data and layout.
    // This makes the layout persist when switching tabs and reopening the file
    this.paper.on('paper:mouseleave', evt => {
      let json_data = JSON.stringify(this.graph.toJSON())
      let newData = {'application/vnd.idaes.model': json_data };
      model.setData({ data: newData });
    });
  }

  private _mimeType: string;
  private holder: HTMLDivElement;
  private paper: joint.dia.Paper;
  private graph: joint.dia.Graph;
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
