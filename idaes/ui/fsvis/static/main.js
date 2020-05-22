// We need to create the graph and paper in the constructor
// If you try to create them in renderModel (which is called for every get request) 
// then you get /images/icons that do not drop
// when the mouseup event is emitted
var standard = joint.shapes.standard
var width = 800;
var height = 800;
var gridSize = 1;
var graph = new joint.dia.Graph([], { cellNamespace: { standard } });
var holder = $("#idaes-canvas")[0];
var paper = new joint.dia.Paper({
  el: holder,
  model: graph,
  cellViewNamespace: { standard },
  width: width,
  height: height,
  gridSize: gridSize,
  interactive: true
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

// /images/icons rotate 90 degrees on right click. Replaces browser context menu
paper.on("element:contextmenu", function(cellView, evt) {
  cellView.model.rotate(90)
});

// Link labels will appear and disapper on right click. Replaces browser context menu
paper.on("link:contextmenu", function(linkView, evt) {
  if (linkView.model.attr('text/display') == 'none') {
    linkView.model.attr({
      'text': {
        display: "block",
      },
      'rect': { fill: 'white', stroke: 'black', 'stroke-width': 1, "fill-opacity": "1" }
    });
  }
  else {
    linkView.model.attr({
      'text': {
        display: "none",
      },
      'rect': { fill: 'white', stroke: 'white', 'stroke-width': 0, "fill-opacity": "0" }
    });
  }
});

// Constrain the elements to the paper on the top and left side as the elements can be lost
// if they go off the paper in those directions. Elements can be recovered from the right 
// and bottom if the user resizes the paper
paper.on('cell:pointermove', function (cellView, evt, x, y) {
  var bbox = cellView.getBBox();
  var constrained = false;

  var constrainedX = x;
  if (bbox.x <= 0) { constrainedX = x + gridSize; constrained = true }

  var constrainedY = y;
  if (bbox.y <= 0) {  constrainedY = y + gridSize; constrained = true }

  //if you fire the event all the time you get a stack overflow
  if (constrained) { cellView.pointermove(evt, constrainedX, constrainedY) }
});

// Get the model from the div tag (see the html file for an explanation)
var data_model = $("#model").data("model");
var model_id = data_model.model.id;
var url = "/fs?id=".concat(model_id)
renderModel(data_model)

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
            console.log(data);
        },
        error: function(error) {
            console.log(error);
        }
    });
});

// Take a model and imports with graph.fromJSON
function renderModel(model) {
    console.log("rendering model..", model.model.id);
    $('#idaes-fs-name').text(model.model.id);  // set flowsheet name
    graph.fromJSON(model);
}


// Set up the help button
// Not implemented yet
var help_button = document.getElementById("help_button");

//help_button.innerText = "Help";
help_button.onclick = () => {
  window.alert("Not implemented yet")
}

// Set up the toggle arc label button
var show_hide_button = $("#show_hide_all_button");
// show_hide_button.innerText = "Show/Hide Arc Labels";
show_hide_button.on('click', function() {
  paper.model.getLinks().forEach(function (link) {
    if (link.attr('text/display') == 'none') {
      link.attr({
        'text': {
          display: "block",
        },
        'rect': { fill: '#d7dce0', stroke: 'white', 'stroke-width': 0, "fill-opacity": "1" }
      });
    }
    else {
      link.attr({
        'text': {
          display: "none",
        },
        'rect': { fill: '#d7dce0', stroke: 'white', 'stroke-width': 0, "fill-opacity": "0" }
      });
    }
  });
});

// When the save button is clicked send a post request to the server with the layout
// We still need to differentiate between the saving the layout in the server and the
// save button save
// I assume that we'll need to modify this to send a signal or something to the server to
// save the json to a file.
$(document).ready( function() {
    $('#save_button').click(function() {
       $.ajax({
            type: 'POST',
            contentType: 'application/json',
            data: JSON.stringify(graph.toJSON()),
            dataType: 'json',
            url: url,
            success: function (e) {
                console.log(e);
            },
            error: function(error) {
                console.log(error);
            }
        });
    });

    $(window).resize(function() {
        var canvas = $("#idaes-canvas")[0].getBoundingClientRect();
        paper.setDimensions(canvas.width, canvas.height);

    });

});
