
export class Paper {
    constructor(app) {
        this._app = app;
        var standard = joint.shapes.standard;
        var width = 800;
        var height = 800;
        var gridSize = 1;

        // Default values for the highlighting events
        this._originalLinkStroke = "#979797";
        this._originalLinkStrokeWidth = 2;
        this._highlightLinkStroke = "#0B79BD";
        this._highlightLinkStrokeWidth = 4;

        this._graph = new joint.dia.Graph([], { cellNamespace: { standard } });
        this._paper = new joint.dia.Paper({
            model: this._graph,
            cellViewNamespace: { standard },
            width: width,
            height: height,
            gridSize: gridSize,
            drawGrid: false,
            interactive: true
        });

        this._paperScroller = new joint.ui.PaperScroller({
            paper: this._paper,
            autoResizePaper: true,
            scrollWhileDragging: true,
            baseWidth: 10,
            baseHeight: 10,
            cursor: 'grab'
        });

        // We need to save this to a variable so that we can access it later
        self = this;

        // The container of the paperscroller needs to be a fixed size or the paperScroller
        // will explode in the y direction when you drag an unit model off of the paper
        // We want all of the elements to be the same width so set the width equal to the 
        // stream table
        let stream_table = document.getElementById("stream-table");
        $('#idaes-canvas').css({ height: stream_table.offsetHeight });
        $("#idaes-canvas")[0].append(self._paperScroller.render().el);

        self.preSetupRegisterEvents();
    }

    get graph() {
        return this._graph
    }

    set graph(data) {
        this._graph.fromJSON(data);
    }

    get paper() {
        return self._paper
    }

    get paperScroller() {
        return self._paperScroller
    }

    translate_for_angle(angle, width, height) {
       // TODO: replace with geometry that considers width and height
       const angle_translation = {0: [0, 5], 90: [38, -35], 180: [0, -72], 270: [-38, -34]};
       return angle_translation[angle];
    }

    /**
     * Register Events before the graph model is loaded
     */
    preSetupRegisterEvents() {
        let model_id = $("#idaes-fs-name").data("flowsheetId");
        let url = "/fs?id=".concat(model_id);

        // Getting the main elements for the idaes canvas and the stream table
        // to be able to dispatch highlighting events to the streams existing
        // on paper and in the stream table
        let idaesCanvas = document.querySelector('#idaes-canvas');
        let streamTable = document.querySelector('#stream-table-data');

        // Setup paper resize on window resize
        window.onresize = function() {
            let stream_table = document.getElementById("stream-table");
            $('#idaes-canvas').css({ height: stream_table.offsetHeight });
        }

        // Registering listeners to idaes-canvas to highlight the correct
        // streams in the paper
        idaesCanvas.addEventListener('HighlightStream', (event) => {
            const relatedLinkElement = idaesCanvas.querySelector(
                `[model-id=${event.detail.streamId}]`
            );
            if (relatedLinkElement) {
                relatedLinkElement.dispatchEvent(new Event('HighlightStream'));
            }
        });
        // Registering listeners to idaes-canvas to remove the highlight from
        // the correct streams in the paper
        idaesCanvas.addEventListener('RemoveHighlightStream', (event) => {
            const relatedLinkElement = idaesCanvas.querySelector(
                `[model-id=${event.detail.streamId}]`
            );
            if (relatedLinkElement) {
                relatedLinkElement.dispatchEvent(new Event('RemoveHighlightStream'));
            }
        });

        // /images/icons rotate 90 degrees on right click. Replaces browser 
        // context menu
        self._paper.on("element:contextmenu", function(cellView, evt) {
            cellView.model.rotate(90)
            // This is needed to keep the text labels for the unit models in the correct orientation
            // x and y were specifically picked to keep the label in the same place 
            // in relation to the unit model (bottom middle)
            // TODO Make this figuring out the x and y positions a function so that we can compute it
            const angle = cellView.model.angle()
            const angle_translation = self.translate_for_angle(angle, 0, 0);
            if (angle_translation === undefined) {
                console.error(`Angle of unit model must be either 0, 90, 180, or 270. Angle is ${angle}`);
            }
            else {
                cellView.model.attr("label/transform", `translate(${angle_translation[0]}, ${angle_translation[1]}) rotate(-${angle})`);
            }
        });

        // Setup event when a link in the paper is hovered upon
        self._paper.on("link:mouseenter", function(linkView) {
            // Adds link tools (adding vertices, moving segments) to links when your
            // mouse over
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
            linkView.addTools(toolsView);
            linkView.showTools();

            // Highlight the corresponding Link and the column in the Stream Table
            const highlightStreamEvent = new CustomEvent(
                'HighlightStream',
                {
                    detail: {
                        streamId: linkView.model.id
                    }
                }
            );
            streamTable.dispatchEvent(highlightStreamEvent);
            idaesCanvas.dispatchEvent(highlightStreamEvent);
        });

        // Setup event when the hovering over link ends
        self._paper.on("link:mouseleave", function(linkView) {
            // Removes the link tools when you leave the link
            linkView.hideTools();

            // Remove the highlight from the link and the column in the
            // Stream Table when the hovering ends
            const removeHighlightStreamEvent = new CustomEvent(
                'RemoveHighlightStream',
                {
                    detail: {
                        streamId: linkView.model.id
                    }
                }
            );
            streamTable.dispatchEvent(removeHighlightStreamEvent);
            idaesCanvas.dispatchEvent(removeHighlightStreamEvent);
        });

        // Send a post request to the server with the new this._graph 
        // This is essentially the saving mechanism (for a server instance) for 
        // right now
        // See the comments above the save button for more saving TODOs
        self._paper.on('paper:mouseleave', () => {this._app.saveModel(url, self._graph)});
        //     $.ajax({
        //         type: 'POST',
        //         contentType: 'application/json',
        //         data: JSON.stringify(self._graph.toJSON()),
        //         dataType: 'json',
        //         url: url,
        //         success: function (data) {
        //         },
        //         error: function(error) {
        //             console.log(error);
        //         }
        //     });
        // });

        // Link labels will appear and disapper on right click. Replaces browser context menu
        self._paper.on("link:contextmenu", function(linkView, evt) {
            if (linkView.model.label(0)["attrs"]["text"]["display"] === 'none') {
                linkView.model.label(0, {
                    attrs: {
                        text: {
                            display: "block",
                        },
                        rect: { "fill-opacity": "1" }
                    } 
                });
            }
            else {
                linkView.model.label(0, {
                    attrs: {
                        text: {
                            display: "none",
                        },
                        rect: { "fill-opacity": "0" }
                    } 
                });
            }
        });
    }

    /**
     * Register Events after the graph model is loaded
     */
    postSetupRegisterEvents() {
        // Setup event listeners for the links in Paper/Graph
        this._graph.getLinks().forEach((link) => {
            let linkView = link.findView(this._paper);
            linkView.el.addEventListener('HighlightStream', () => {
                linkView.model.attr({
                    line: {
                        stroke: this._highlightLinkStroke,
                        'stroke-width': this._highlightLinkStrokeWidth
                    }
                });
            });

            linkView.el.addEventListener('RemoveHighlightStream', () => {
                linkView.model.attr({
                    line: {
                        stroke: this._originalLinkStroke,
                        'stroke-width': this._originalLinkStrokeWidth
                    }
                });
            });
        });
    }

    /**
     * Setup the graph model
     */
    setup(model) {
        this._graph.fromJSON(model);
        this.postSetupRegisterEvents();
    }
    
};
