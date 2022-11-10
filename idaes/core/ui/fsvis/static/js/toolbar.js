/**
 * The Institute for the Design of Advanced Energy Systems Integrated Platform
 * Framework (IDAES IP) was produced under the DOE Institute for the
 * Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
 * by the software owners: The Regents of the University of California, through
 * Lawrence Berkeley National Laboratory,  National Technology & Engineering
 * Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
 * Research Corporation, et al.  All rights reserved.
 *
 * Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
 * license information.
 */


/**
 * This class is responsible for all the two toolbars that exist in the page.
 * 
 * The first toolbar is the one that is seen in the page's header. And the
 * second toolbar is seen on the top right corner of the visualizer paper
 * graph.
 */
export class Toolbar { 
    constructor(app, paper, stream_table) {
        this._app = app;
        this._paper = paper;
        this._stream_table = stream_table;
        this.setupPageToolbar();
        this.setupFlowsheetToolbar();
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

    setContainerView(container) {
        // Set the viewablility of the container
        if (container.style.display === "none") {
                container.style.display = "block";
        }
        else {
            container.style.display = "none";
        };
    };

    addSpinner($this, loadingText, timeout) {
        // Add a spinner to a button for the length of timeout (in milliseconds)
        if ($this.html() !== loadingText) {
            $this.data('original-text', $this.html());
            $this.html(loadingText);
        }
        setTimeout(function() {
            $this.html($this.data('original-text'));
        }, timeout);
    }

    setupPageToolbar() {
        // Grab the model information from the div tag so that we can use it in our ajax calls
        const model_id = $("#idaes-fs-name").data("flowsheetId");
        const url = `/fs?id=${ model_id }`;
        const model_server_url = $("#model-server-url").data("modelurl");

        // Moved spinner adding before other events so the spinner gets added first
        let app = this;
        // Add spinner to refresh button on click
        $('#refresh-btn').on('click', function() {
            var $this = $(this);
            var loadingText = '<i class="fa fa-circle-o-notch fa-spin"></i>Refresh';
            app.addSpinner($this, loadingText, 1000)
         });

        // Add spinner to save button on click
        $('#save-btn').on('click', function() {
            var $this = $(this);
            var loadingText = '<i class="fa fa-circle-o-notch fa-spin"></i>Save';
            app.addSpinner($this, loadingText, 1000)
         });

        // Save event listener
        document.querySelector("#save-btn").addEventListener("click", () => {
            this._app.saveModel(url, this._paper.graph);
        });

        // Refresh event listener
        document.querySelector("#refresh-btn").addEventListener("click", () => {
            this._app.refreshModel(url, this._paper)
        });

        // Flowsheet to PNG export event listener
        document.querySelector("#export-flowsheet-png-btn").addEventListener("click", () => {
            let p = this._paper.paper;
            // Make sure to hide all of the vertices and bars on the links
            // so they don't show up in the PNG
            p.hideTools();
            p.toPNG(function(png) {
                new joint.ui.Lightbox({
                    image: png,
                    downloadable: true,
                    fileName: model_id.concat(".png")
                }).open();
            }, {
                preserveDimensions: true,
                convertImagesToDataUris: true,
                useComputedStyles: true,
                stylesheet: '.scalable * { vector-effect: non-scaling-stroke }'
            });
        });

        // Flowsheet to JPEG export event listener
        document.querySelector("#export-flowsheet-jpg-btn").addEventListener("click", () => {
            let p = this._paper.paper;
            // Make sure to hide all of the vertices and bars on the links
            // so they don't show up in the JPEG
            p.hideTools();
            p.toJPEG(function(jpeg) {
                new joint.ui.Lightbox({
                    image: jpeg,
                    downloadable: true,
                    fileName: model_id.concat(".jpeg")
                }).open();
            }, {
                preserveDimensions: true,
                convertImagesToDataUris: true,
                useComputedStyles: true,
                stylesheet: '.scalable * { vector-effect: non-scaling-stroke }'
            });
        });

        // Flowsheet to SVG export event listener
        document.querySelector("#export-flowsheet-svg-btn").addEventListener("click", () => {
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

        // Stream table to CSV export event listener
        document.querySelector("#export-stream-table-btn").addEventListener("click", () => {
            this._stream_table._gridOptions.api.exportDataAsCsv({suppressQuotes: true});
        });

        // Show/hide flowsheet event listener
        document.querySelector("#view-flowsheet-btn").addEventListener("change", function() {
            let container = document.querySelector("#idaes-canvas");
            if (this.checked) {
                container.style.display = "block";
            }
            else {
                container.style.display = "none";
            };
        });

        // Show/hide stream table event listener
        document.querySelector("#view-stream-table-btn").addEventListener("change", function() {
            let container = document.querySelector("#stream-table");
            if (this.checked) {
                container.style.display = "block";
            }
            else {
                container.style.display = "none";
            };
        });
    };

    setupFlowsheetToolbar() {
        let app = this;

        // Labels toggle event listener
        document.querySelector("#labels-toggle").addEventListener("change", function() {
            if (this.checked) {
                app._paper._graph.getLinks().forEach(function (link) {
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
                app._paper._graph.getLinks().forEach(function (link) {
                    link.label(0, {
                        attrs: {
                            text: {
                                display: "none",
                            },
                            rect: { "fill-opacity": "0" }
                        }
                    });
                });
            };
        });

        // Streams toggle event listener
        document.querySelector("#stream-names-toggle").addEventListener("change", function() {
            if (this.checked) {
                app._paper._graph.getLinks().forEach(function (link) {
                    link.label(1, {
                        attrs: {
                            text: {
                                display: "block",
                            }
                        }
                    });
                });
            }
            else {
                app._paper._graph.getLinks().forEach(function (link) {
                    link.label(1, {
                        attrs: {
                            text: {
                                display: "none",
                            }
                        }
                    });
                });
            };
        });

        // Grid toggle event listener
        document.querySelector("#grid-toggle").addEventListener("change", function() {
            if (this.checked) {
                app.setGrid(10, '#BBBBBB');
            }
            else {
                app.setGrid(10, '#FFFFFF');
            };
        });

        // Zoom in event listener
        document.querySelector("#zoom-in-btn").addEventListener("click", () => {
            this._paper.paperScroller.zoom(0.2, { max: 4 });
        });

        // Zoom out event listener
        document.querySelector("#zoom-out-btn").addEventListener("click", () => {
            this._paper.paperScroller.zoom(-0.2, { min: 0.2 });
        });

        // Zoom to fit event listener
        document.querySelector("#zoom-fit-btn").addEventListener("click", () => {
            this._paper.zoomToFit();
        });
    }

}
