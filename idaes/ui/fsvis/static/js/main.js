import { Paper } from './paper.js';

// Take a model and imports with this.graph.fromJSON
function renderModel(model, paper) {
    $('#idaes-fs-name').text(model.model.id);  // set flowsheet name
    paper.graph.fromJSON(model);
}

function getFlowsheet(flowsheetId) {
    const url = `/fs?id=${ flowsheetId }`;
    return $.ajax({
        type: 'GET',
        url: url,
        datatype: 'json',
        error: (error) => { console.log(error) }
    })
}

// =====================
//    Main function
// =====================
$( document ).ready(function() {
    // Get the model from the div tag (see the html file for an explanation)
    let flowsheetId = $("#idaes-fs-name").data("flowsheetId");
    let paper = new Paper();
    getFlowsheet(flowsheetId)
        .done( (model) => {
            renderModel(model, paper);
        });
});
