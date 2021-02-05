
export class Paper {
    constructor(app) {
        this._app = app;
        var standard = joint.shapes.standard;
        var width = 800;
        var height = 800;
        var gridSize = 1;

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
        let stream_table = document.getElementById("stream-table-container");
        $('#idaes-canvas').css({ width: stream_table.offsetWidth, height: 800 });
        $("#idaes-canvas")[0].append(self._paperScroller.render().el);

        self.setupEvents();

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

    setupEvents() {
        let model_id = $("#idaes-fs-name").data("flowsheetId");
        let url = "/fs?id=".concat(model_id);

        // /images/icons rotate 90 degrees on right click. Replaces browser 
        // context menu
        self._paper.on("element:contextmenu", function(cellView, evt) {
            cellView.model.rotate(90)
            // This is needed to keep the text labels for the unit models in the correct orientation
            // x and y were specifically picked to keep the label in the same place 
            // in relation to the unit model (bottom middle)
            // TODO Make this figuring out the x and y positions a function so that we can compute it
            const angle = cellView.model.angle()
            switch (angle) {
                case 0: {
                    const x = 0;
                    const y = 5;
                    cellView.model.attr("label/transform", `translate(${x}, ${y}) rotate(-${angle})`)
                    break;
                }
                case 90: {
                    const x = 38;
                    const y = -35;
                    cellView.model.attr("label/transform", `translate(${x}, ${y}) rotate(-${angle})`)
                    break;
                }
                case 180: {
                    const x = 0;
                    const y = -72;
                    cellView.model.attr("label/transform", `translate(${x}, ${y}) rotate(-${angle})`)
                    break;
                }
                case 270: {
                    const x = -38;
                    const y = -34;
                    cellView.model.attr("label/transform", `translate(${x}, ${y}) rotate(-${angle})`)
                    break;
                }
                default: {
                    const x = 0;
                    const y = 0;
                    cellView.model.attr("label/transform", `translate(${x}, ${y}) rotate(-${angle})`)
                }
            }
        });

        // Adds link tools (adding vertices, moving segments) to links when your 
        // mouse over
        self._paper.on("link:mouseover", function(cellView, evt) {
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
        });

        // Removes the link tools when you leave the link
        self._paper.on("link:mouseout", function(cellView, evt) {
            cellView.hideTools()
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
    
};
