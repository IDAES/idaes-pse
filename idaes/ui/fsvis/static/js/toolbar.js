export class Toolbar { 
    constructor(graph, paper, paperScroller) {
        this._graph = graph;
        this._paper = paper;
        this._paperScroller = paperScroller;

        // We need to save this to a variable so that we can access it later
        self = this;

        self.setupToolbar();
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

        var toolbar = new joint.ui.Toolbar({
            autoToggle: true,
            references: {
                paper: self._paper,
                paperScroller: self._paperScroller
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

        toolbar.on('refresh:pointerclick', event => this.refreshModel(url, self._paper));

        // Note: it is not good style to embed all the processing here.
        // Instead, call functions as in the refresh:pointerclick event above.

        toolbar.on('labels:change', function(value, event) {
            // Go through all of the links and set the display values
            if (value == true) {
                self._paper.model.getLinks().forEach(function (link) {
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
                self._paper.model.getLinks().forEach(function (link) {
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

        toolbar.on('save:pointerclick', function(event) {
            // Send an ajax POST request to the flask server with the jointjs model
            // This request also sets the request header to Source: save_button so 
            // that the server knows to save the jointjs model to a file
            // If this isn't set then the flask server will just save the jointjs 
            // to the database rather than saving it to the file
            $.ajax({
                type: 'POST',
                contentType: 'application/json',
                data: JSON.stringify(self._graph.toJSON()),
                dataType: 'json',
                url: url,
                beforeSend: function(request) {
                    request.setRequestHeader("Source", "save_button");
                },
                success: function (e) {
                },
                error: function(error) {
                    console.log(error);
                }
            });
        });

        toolbar.on('svg:pointerclick', function(event) {
            // Make sure to hide all of the vertices and bars on the links 
            // so they don't show up in the SVG
            self._paper.hideTools()
            self._paper.toSVG(function(svg) {
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

        toolbar.on('grid:change', function(value, event) {
            if (value == true) {
                self.setGrid(10, '#BBBBBB');
            }
            else {
                self.setGrid(10, '#FFFFFF');
            }
        });

        toolbar.on('sizeBox:option:select', function(value, event) {
            $('#idaes-canvas').css({ width: value["content"], height: value["content"] });
        });

        toolbar.on('help:pointerclick', function(event) {
            window.open("https://idaes-pse.readthedocs.io/en/stable/user_guide/vis/index.html")
        });

        $('#toolbar-container').append(toolbar.render().el);

        self.setGrid(10, '#FFFFFF');
    }

    /**
     * Inform the user of some event.
     *
     * Do NOT use this for internal messages or debugging.
     *
     * @param level The level of 'severity' of the message. 0=info, 1=warning, 2=error
     * @param message The message to show
     * @param duration Duration, in seconds, to show the message. 0=forever
     */
    informUser(level, message, duration) {
        // TODO: Write into a status area
        // Write to console
        switch(level) {
            case 0:
                console.log(message);
                break;
            case 1:
                console.warn(message);
                break;
            case 2:
                console.error(message);
                break;
            default:
                console.log(message);
        }
    }

    /**
     * Save current model value and then update with value in the Python process.
     *
     * This makes two calls to the server: first a PUT to save the current model, and
     * second a GET to retrieve the new values. If either of these fails, the method will fail
     * and not make any changes to its inputs.
     *
     * If this succeeds, the value of the model in the Paper instance will be replaced with the
     * new value sent from the server in the Python process.
     *
     * @param url The HTTP server that is running in the Python process
     * @param paper Instance of Paper that has model in 'model' attribute.
     */
    refreshModel(url, paper) {
        // Inform user of progress (1)
        console.debug("paper.model=", paper.model);
        this.informUser(0, "Refresh: save current values from model");
        // First save our version of the model
        let clientModel = paper.model;
        let clientData = JSON.stringify(clientModel.toJSON());
        console.debug(`Sending to ${url}: ` + clientData);
        $.ajax({url: url, type: 'PUT', contentType: "application/json", data: clientData})
            // On failure inform user and stop
            .fail(error => this.informUser(
                2, "Fatal error: cannot save current model before refresh: " + error))
            // On success, continue on to fetch new model
            .done(data => {
                // Inform user of progress (2)
                this.informUser(0, "Refresh: load new model values from Python program");
                $.ajax({url: url, dataType: "json"})
                    // If we got the model, save it
                    .done(data => {paper.model.fromJSON(data)})
                    // Otherwise fail
                    .fail((jqXHR, textStatus, errorThrown) => {
                        this.informUser(2, "Fatal error: Could not retrieve new model from Python program: " +
                            textStatus + ", error=" + errorThrown);
                    });
            });
    }
}
