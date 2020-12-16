import { Paper } from './paper.js';
import { Toolbar } from './toolbar.js';

export class App {
    constructor (flowsheetId) {
        this.paper = new Paper(this);
        const url = `/fs?id=${ flowsheetId }`;
        $.ajax({url: url, datatype: 'json'})
            .done((model) => {
                this.renderModel(model);
                this.toolbar = new Toolbar(this, this.paper);
            })
            .fail((error) => { console.log(error) });
    }

    renderModel(model) {
        $('#idaes-fs-name').text(model.model.id);  // set flowsheet name
        this.paper.graph.fromJSON(model);
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
        // console.debug("paper.model=", paper.model);
        this.informUser(0, "Refresh: save current values from model");
        // First save our version of the model
        let clientModel = paper.graph;
        let clientData = JSON.stringify(clientModel.toJSON());
        // console.debug(`Sending to ${url}: ` + clientData);
        $.ajax({url: url, type: 'PUT', contentType: "application/json", data: clientData})
            // On failure inform user and stop
            .fail(error => this.informUser(
                2, "Fatal error: cannot save current model before refresh: " + error))
            // On success, continue on to fetch new model
            .done(() => {
                // Inform user of progress (2)
                this.informUser(0, "Refresh: load new model values from Python program");
                $.ajax({url: url, dataType: "json"})
                    // If we got the model, save it
                    .done(data => {paper.graph = data}) // uses setter
                    // Otherwise fail
                    .fail((jqXHR, textStatus, errorThrown) => {
                        this.informUser(2, "Fatal error: Could not retrieve new model from Python program: " +
                            textStatus + ", error=" + errorThrown);
                    });
            });
    }

    /**
     * Save the model value.
     *
     * This sends a PUT to the server to save the current model value.
     *
     * @param url The HTTP server that is running in the Python process
     * @param model The model to save
     */
    saveModel(url, model) {
        let clientData = JSON.stringify(model.toJSON());
        // console.debug(`Sending to ${url}: ` + clientData);
        this.informUser(0, "Save current values from model");
        $.ajax({url: url, type: 'PUT', contentType: "application/json", data: clientData})
            // On failure inform user and stop
            .fail(error => this.informUser(
                2, "Fatal error: cannot save current model: " + error))
            .done(() => {
                this.informUser(0, "Saved new model values");
            });
    }
}

// =====================
//    Main function
// =====================
$( document ).ready(function() {
    let flowsheetId = $("#idaes-fs-name").data("flowsheetId");
    globalThis.app = new App(flowsheetId);
});
