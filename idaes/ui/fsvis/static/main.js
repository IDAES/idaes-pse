joint.setTheme('modern');

// We need to create the graph and paper in the beginning
// If you try to create them in renderModel (which is called for every get request) 
// then you get /images/icons that do not drop
// when the mouseup event is emitted
var standard = joint.shapes.standard
var width = 800;
var height = 800;
var gridSize = 1;
var graph = new joint.dia.Graph([], { cellNamespace: { standard } });
var paper = new joint.dia.Paper({
    model: graph,
    cellViewNamespace: { standard },
    width: width,
    height: height,
    gridSize: gridSize,
    drawGrid: true,
    interactive: true
});

// Get the model from the div tag (see the html file for an explanation)
var data_model = $("#model").data("model");
var model_id = data_model.model.id;
var url = "/fs?id=".concat(model_id)
renderModel(data_model)

// var paperScroller = new joint.ui.PaperScroller({
//     paper: paper,
//     autoResizePaper: true,
//     scrollWhileDragging: true,
//     cursor: 'grab'
// });

// $("#idaes-canvas")[0].append(paperScroller.render().el);

var toolbar = new joint.ui.Toolbar({
    autoToggle: true,
    references: {
        paper: paper,
    },
    tools: [
        { type: 'toggle', name: 'labels', label: 'Labels' },
        { type: 'separator' },
        { type: 'button', name: 'zoomIn', text: 'Zoom In'}
    ]
});

toolbar.on('zoomIn:click', function(value, event) {
    paperScroller.zoom(0.2, { max: 4 });
});

toolbar.on('labels:change', function(value, event) {
    if (value == true) {
        paper.model.getLinks().forEach(function (link) {
            link.label(0, {
                attrs: {
                    text: {
                        text: link.label(0)["attrs"]["text"]["text"],
                        fill: link.label(0)["attrs"]["text"]["fill"],
                        "text-anchor": link.label(0)["attrs"]["text"]["text-anchor"],
                        display: "block",
                    },
                    rect: { fill: '#d7dce0', stroke: 'white', 'stroke-width': 0, "fill-opacity": "1" }
                }  
            });
        });
    }
    else {
        paper.model.getLinks().forEach(function (link) {
            link.label(0, {
                attrs: {
                    text: {
                        text: link.label(0)["attrs"]["text"]["text"],
                        fill: link.label(0)["attrs"]["text"]["fill"],
                        "text-anchor": link.label(0)["attrs"]["text"]["text-anchor"],
                        display: "none",
                    },
                    rect: { fill: '#d7dce0', stroke: 'white', 'stroke-width': 0, "fill-opacity": "0" }
                }
            });
        });
    }
});

$('#toolbar-container').append(toolbar.render().el);

// stencil_obj.render().load();

// If we decide to do this, add a button to the side panel that enables/disables 
// it (just change the grid color to white)

// function setGrid(paper, gridSize, color) {
//     // Set grid size on the JointJS paper object (joint.dia.Paper instance)
//     paper.options.gridSize = gridSize;
//     // Draw a grid into the HTML 5 canvas and convert it to a data URI image
//     var canvas = $('<canvas/>', { width: gridSize, height: gridSize });
//     canvas[0].width = gridSize;
//     canvas[0].height = gridSize;
//     var context = canvas[0].getContext('2d');
//     context.beginPath();
//     context.rect(1, 1, 1, 1);
//     context.fillStyle = color || '#AAAAAA';
//     context.fill();
//     // Finally, set the grid background image of the paper container element.
//     var gridBackgroundImage = canvas[0].toDataURL('image/png');
//     paper.$el.css('background-image', 'url("' + gridBackgroundImage + '")');
// }

// setGrid(paper, 10, '#FF0000');

// /images/icons rotate 90 degrees on right click. Replaces browser context menu
paper.on("element:contextmenu", function(cellView, evt) {
    cellView.model.rotate(90)
});

// Adds link tools (adding vertices, moving segments) to links when your mouse over
paper.on("link:mouseover", function(cellView, evt) {
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
paper.on("link:mouseout", function(cellView, evt) {
    cellView.hideTools()
});



// Send a post request to the server with the new graph 
// This is essentially the saving mechanism (for a server instance) for right now
// See the comments above the save button for more saving TODOs
paper.on('paper:mouseleave', evt => {
    $.ajax({
        type: 'POST',
        contentType: 'application/json',
        data: JSON.stringify(graph.toJSON()),
        dataType: 'json',
        url: url,
        success: function (data) {
        },
        error: function(error) {
            console.log(error);
        }
    });
});

// Take a model and imports with graph.fromJSON
function renderModel(model) {
    $('#idaes-fs-name').text(model.model.id);  // set flowsheet name
    graph.fromJSON(model);
}


// // Set up the help button
// // Not implemented yet
// var help_button = document.getElementById("help_button");

// //help_button.innerText = "Help";
// help_button.onclick = () => {
//     window.open("https://idaes-pse.readthedocs.io/en/stable/user_guide/vis/index.html")
// }

// Link labels will appear and disapper on right click. Replaces browser context menu
paper.on("link:contextmenu", function(linkView, evt) {
    if (linkView.model.label(0)["attrs"]["text"]["display"] == 'none') {
        linkView.model.label(0, {
            attrs: {
                text: {
                    text: linkView.model.label(0)["attrs"]["text"]["text"],
                    fill: linkView.model.label(0)["attrs"]["text"]["fill"],
                    "text-anchor": linkView.model.label(0)["attrs"]["text"]["text-anchor"],
                    display: "block",
                },
                rect: { fill: '#d7dce0', stroke: 'white', 'stroke-width': 0, "fill-opacity": "1" }
            }  
        });
    }
    else {
        linkView.model.label(0, {
            attrs: {
                text: {
                    text: linkView.model.label(0)["attrs"]["text"]["text"],
                    fill: linkView.model.label(0)["attrs"]["text"]["fill"],
                    "text-anchor": linkView.model.label(0)["attrs"]["text"]["text-anchor"],
                    display: "none",
                },
                rect: { fill: '#d7dce0', stroke: 'white', 'stroke-width': 0, "fill-opacity": "0" }
            }
        });
    }
});

// When the save button is clicked send a post request to the server with the layout
// We still need to differentiate between the saving the layout in the server and the
// save button save
// I (Makayla) assume that we'll need to modify this to send a signal or something to the server to
// save the json to a file.
$(document).ready( function() {
    $('#save_button').click(function() {
       $.ajax({
            type: 'POST',
            contentType: 'application/json',
            data: JSON.stringify(graph.toJSON()),
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

});
