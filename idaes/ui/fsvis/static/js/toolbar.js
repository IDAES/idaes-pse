export class Toolbar { 
    constructor(app, paper) {
        this._app = app;
        this._paper = paper;
        this.setupToolbar();
    }

    setGrid(gridSize, color) {
        // Set grid size on the JointJS paper object (joint.dia.paper instance)
        self._paper.options.gridSize = gridSize;
        // Draw a grid into the HTML 5 canvas and convert it to a data URI image
        let canvas = $('<canvas/>', { width: gridSize, height: gridSize });
        let context = canvas[0].getContext('2d');
        canvas[0].width = gridSize;
        canvas[0].height = gridSize;
        context.beginPath();
        context.rect(1, 1, 1, 1);
        context.fillStyle = color || '#AAAAAA';
        context.fill();
        // Finally, set the grid background image of the paper container element.
        self._paper.$el.css('background-image', 'url("' + canvas[0].toDataURL('image/png') + '")');
    };

    setupToolbar() {
        // Grab the model information from the div tag so that we can use it in our ajax calls
        let model_id = $("#idaes-fs-name").data("flowsheetId");
        let url = `/fs?id=${ model_id }`;
        let model_server_url = $("#model-server-url").data("modelurl");

        let toolbar = new joint.ui.Toolbar({
            autoToggle: true,
            references: {
                paper: this._paper,
                paperScroller: this._paper.paperScroller
            },
            tools: [
                { type: 'button', name: 'refresh', text: 'Refresh Graph'},
                { type: 'separator' },
                { type: 'toggle', name: 'labels', label: 'Labels:' },
                { type: 'separator' },
                { type: 'button', name: 'save', text: 'Save'},
                { type: 'separator' },
                { type: 'button', name: 'svg', text: 'Export SVG'},
                { type: 'separator' },
                { type: 'toggle', name: 'grid', label: 'Grid:'},
                { type: 'separator' },
                { type: 'zoomIn' },
                { type: 'zoomOut' },
                { type: 'zoomToFit' },
                { type: 'separator' },
                { type: 'label', text: 'Canvas size:' },
                { type: 'selectBox', name: 'sizeBox', options: [{content: 800, selected: true}, {content: 1600, selected: true}, {content: 3200}, {content: 6400}] },
                { type: 'separator' },
                { type: 'button', name: 'help', text: 'Help'},
            ],
        });

        toolbar.on('refresh:pointerclick', () => this._app.refreshModel(url, this._paper));

        // Note: it is not good style to embed all the processing here.
        // Instead, call functions as in the refresh:pointerclick event above.

        toolbar.on('labels:change', (value, event) => {
            // Go through all of the links and set the display values
            if (value === true) {
                this._paper._graph.getLinks().forEach(function (link) {
                    link.label(0, {
                        attrs: {
                            text: {
                                display: "block",
                            },
                            rect: { "fill-opacity": "1" }
                        }
                    });
                });
            }
            else {
                this._paper._graph.getLinks().forEach(function (link) {
                    link.label(0, {
                        attrs: {
                            text: {
                                display: "none",
                            },
                            rect: { "fill-opacity": "0" }
                        }
                    });
                });
            }
        });

        toolbar.on('save:pointerclick', () => this._app.saveModel(url, this._paper.graph));
        // {
        //     // Send an ajax POST request to the flask server with the jointjs model
        //     // This request also sets the request header to Source: save_button so
        //     // that the server knows to save the jointjs model to a file
        //     // If this isn't set then the flask server will just save the jointjs
        //     // to the database rather than saving it to the file
        //     $.ajax({
        //         type: 'POST',
        //         contentType: 'application/json',
        //         data: JSON.stringify(self._graph.toJSON()),
        //         dataType: 'json',
        //         url: url,
        //         beforeSend: function(request) {
        //             request.setRequestHeader("Source", "save_button");
        //         },
        //         success: function (e) {
        //         },
        //         error: function(error) {
        //             console.log(error);
        //         }
        //     });
        // });

        toolbar.on('svg:pointerclick', (event) => {
            let p = this._paper.paper;
            // Make sure to hide all of the vertices and bars on the links 
            // so they don't show up in the SVG
            p.hideTools();
            p.toSVG(function(svg) {
                new joint.ui.Lightbox({
                    image: 'data:image/svg+xml,' + encodeURIComponent(svg),
                    downloadable: true,
                    fileName: model_id.concat(".svg")
                }).open();
            }, {
                preserveDimensions: true,
                convertImagesToDataUris: true,
                useComputedStyles: true,
                stylesheet: '.scalable * { vector-effect: non-scaling-stroke }'
            });
        });

        toolbar.on('grid:change', (value, event) => {
            if (value === true) {
                this.setGrid(10, '#BBBBBB');
            }
            else {
                this.setGrid(10, '#FFFFFF');
            }
        });

        toolbar.on('sizeBox:option:select', (value, event) => {
            $('#idaes-canvas').css({ width: value["content"], height: value["content"] });
        });

        toolbar.on('help:pointerclick', (event) => {
            globalThis.open("https://idaes-pse.readthedocs.io/en/stable/user_guide/vis/index.html")
        });

        $('#toolbar-container').append(toolbar.render().el);

        this.setGrid(10, '#FFFFFF');
    }

}
