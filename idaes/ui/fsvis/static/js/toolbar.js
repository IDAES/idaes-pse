export class Toolbar { 
    constructor(graph, paper, paperScroller, model_id) {
        this._graph = graph;
        this._paper = paper;
        this._paperScroller = paperScroller;
        self = this;
        self.setupToolbar();
    };

    setGrid(gridSize, color) {
        // Set grid size on the JointJS paper object (joint.dia.paper instance)
        self._paper.options.gridSize = gridSize;
        // Draw a grid into the HTML 5 canvas and convert it to a data URI image
        var canvas = $('<canvas/>', { width: gridSize, height: gridSize });
        canvas[0].width = gridSize;
        canvas[0].height = gridSize;
        var context = canvas[0].getContext('2d');
        context.beginPath();
        context.rect(1, 1, 1, 1);
        context.fillStyle = color || '#AAAAAA';
        context.fill();
        // Finally, set the grid background image of the paper container element.
        var gridBackgroundImage = canvas[0].toDataURL('image/png');
        self._paper.$el.css('background-image', 'url("' + gridBackgroundImage + '")');
    };

    setupToolbar() {
        var data_model = $("#model").data("model");
        var model_id = data_model.model.id;
        var url = "/fs?id=".concat(model_id);

        var toolbar = new joint.ui.Toolbar({
            autoToggle: true,
            references: {
                paper: self._paper,
                paperScroller: self._paperScroller
            },
            tools: [
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

        toolbar.on('labels:change', function(value, event) {
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
    };
};
