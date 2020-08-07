import { Paper } from './paper.js';

// Take a model and imports with this.graph.fromJSON
function renderModel(model, paper) {
    $('#idaes-fs-name').text(model.model.id);  // set flowsheet name
    paper.graph.fromJSON(model);
}

$( document ).ready(function() {
    // Get the model from the div tag (see the html file for an explanation)
    var data_model = $("#model").data("model");

    var paper = new Paper();
    
    renderModel(data_model, paper);
});
