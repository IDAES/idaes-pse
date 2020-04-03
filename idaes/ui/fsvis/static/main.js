document.getElementById("text").innerHTML = "Javascript run!";

import * as jquery from "/lib/node_modules/jquery/jquery.js"
import * as lodash from "/lib/node_modules/lodash/lodash.js"
import * as backbone from "/lib/node_modules/backbone/backbone.js";
import { dia } from '/lib/jointjs/src/core.mjs';
import * as standard from '/lib/jointjs/src/shapes/standard.mjs';

console.log(JSON.stringify(dia));

var holder = document.getElementById("idaes_flowsheet_visualizer");

// We need to create the graph and paper in the constructor
// If you try to create them in renderModel (which is called everytime the user
// opens an .idaes.vis file or changes tabs) then you get icons that do not drop
// when the mouseup event is emitted
var width = 10000;
var height = 10000;
var gridSize = 1;
var graph = new dia.Graph([], { cellNamespace: { standard } });
var paper = new dia.Paper({
  el: holder,
  model: graph,
  cellViewNamespace: { standard },
  width: width,
  height: height,
  gridSize: gridSize,
  interactive: true
});

var show_hide_button = document.createElement("button");

// show_hide_button.innerText = "Show/Hide Arc Labels";  
// show_hide_button.onclick = () => {
//   paper.model.getLinks().forEach(function (link) {        
//     if (link.attr('text/display') == 'none') {
//       link.attr({
//         'text': {
//           display: "block",
//         },
//         'rect': { fill: '#d7dce0', stroke: 'white', 'stroke-width': 0, "fill-opacity": "1" }
//       });
//     }
//     else {
//       link.attr({
//         'text': {
//           display: "none",
//         },
//         'rect': { fill: '#d7dce0', stroke: 'white', 'stroke-width': 0, "fill-opacity": "0" }
//       });
//     }
//   });
// }

// var save_button = document.createElement("button");

// save_button.innerText = "Save";  
// save_button.onclick = () => {
//   console.log("Not implemented yet")
// }

// var help_button = document.createElement("button");

// help_button.innerText = "Help";  
// help_button.onclick = () => {
//   console.log("Not implemented yet")
// }

// // Adds link tools (adding vertices, moving segments) to links when your mouse over
// paper.on("link:mouseover", function(cellView, evt) {
//   var verticesTool = new joint.linkTools.Vertices({
//     focusOpacity: 0.5,
//     redundancyRemoval: true,
//     snapRadius: 20,
//     vertexAdding: true,
//   });
//   var segmentsTool = new joint.linkTools.Segments();

//   var sourceArrowheadTool = new joint.linkTools.SourceArrowhead();
//   var targetArrowheadTool = new joint.linkTools.TargetArrowhead();

//   var toolsView = new joint.dia.ToolsView({
//     tools: [
//       verticesTool, segmentsTool,
//       sourceArrowheadTool, targetArrowheadTool
//     ]
//   });
//   cellView.addTools(toolsView)
//   cellView.showTools()
// });

// // Removes the link tools when you leave the link
// paper.on("link:mouseout", function(cellView, evt) {
//   cellView.hideTools()
// });

// // Icons rotate 90 degrees on right click. Replaces browser context menu
// paper.on("element:contextmenu", function(cellView, evt) {
//   cellView.model.rotate(90)
// });

// // Link labels will appear and disapper on right click. Replaces browser context menu
// paper.on("link:contextmenu", function(linkView, evt) {
//   if (linkView.model.attr('text/display') == 'none') {
//     linkView.model.attr({
//       'text': {
//         display: "block",
//       },
//       'rect': { fill: 'white', stroke: 'black', 'stroke-width': 1, "fill-opacity": "1" }
//     });
//   }
//   else {
//     linkView.model.attr({
//       'text': {
//         display: "none",
//       },
//       'rect': { fill: 'white', stroke: 'white', 'stroke-width': 0, "fill-opacity": "0" }
//     });
//   }
// });

// // Constrain the elements to the paper
// paper.on('cell:pointermove', function (cellView, evt, x, y) {

//   var bbox = cellView.getBBox();
//   var constrained = false;

//   var constrainedX = x;

//   if (bbox.x <= 0) { constrainedX = x + gridSize; constrained = true }
//   if (bbox.x + bbox.width >= width) { constrainedX = x - gridSize; constrained = true }

//   var constrainedY = y;

//   if (bbox.y <= 0) {  constrainedY = y + gridSize; constrained = true }
//   if (bbox.y + bbox.height >= height) { constrainedY = y - gridSize; constrained = true }

//   //if you fire the event all the time you get a stack overflow
//   if (constrained) { cellView.pointermove(evt, constrainedX, constrainedY) }
// });
