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

import { Paper } from './paper.js';
import { StreamTable } from './stream_table.js';
import { Toolbar } from './toolbar.js';
import { JointJsCellConfig } from './cell_config.js';


/**
 * The main client app responsible for IDAES related visualizations. Here in
 * the main file, we:
 *     1. Render the model.
 *     2. Display the stream table.
 */
export class App {
    constructor (flowsheetId) {
        this.paper = new Paper(this);
        const url = `/fs?id=${ flowsheetId }`;

        // Adding a special flag to mark that the graph changed
        this._is_graph_changed = false;
        // Setting name (key) that defines the save model time interval
        this._save_time_interval_key = 'save_time_interval';
        this._default_save_time_interval = 5000; // Default time interval
        this._save_time_interval = this.getSaveTimeInterval();

        this.setupGraphChangeChecker(this._save_time_interval);

        $.ajax({url: url, datatype: 'json'})
            .done((model) => {
                this.renderModel(model);
                this.stream_table = new StreamTable(this, model);
                this.toolbar = new Toolbar(this, this.paper, this.stream_table);
            })
            .fail((xhr, status, error) => { 
                console.log(error);
                console.log(status); 
            });
        // Make the dropdowns not disappear on use. They will disappear when the user clicks off the dropdown
        $(document).on('click', '.dropdown-menu', function (e) {
            e.stopPropagation();
        });
    }

    renderModel(model) {
        $('#idaes-fs-name').text(model.model.id);  // set flowsheet name
        var jjCellConfig = new JointJsCellConfig(model);
        var processed_model = jjCellConfig.processRoutingConfig();
        this.paper.setup(processed_model);
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
        this.informUser(0, "Refresh: save current values from model");
        // First save our version of the model
        let clientModel = paper.graph;
        let clientData = JSON.stringify(clientModel.toJSON());
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
                    .done(data => {
                        // Display views before refreshing
                        const viewFlowsheet = document.querySelector("#view-flowsheet-btn");
                        const viewStreamTable = document.querySelector("#view-stream-table-btn");

                        const clickEvent = new MouseEvent('click');

                        if (!viewFlowsheet.checked) {
                            viewFlowsheet.dispatchEvent(clickEvent);
                        }
                        if (!viewStreamTable.checked) {
                            viewStreamTable.dispatchEvent(clickEvent);
                        }

                        // Refresh
                        this.renderModel(data);
                        this.stream_table.initTable(data);
                    })
                    // Otherwise fail
                    .fail((jqXHR, textStatus, errorThrown) => {
                        this.informUser(2, "Fatal error: Could not retrieve new model from Python program: " +
                            textStatus + ", error=" + errorThrown);
                    });
            });
    }

    /**
     * Get the save time interval value from the application's setting block.
     */
    getSaveTimeInterval() {
        let settings_url = "/setting?setting_key=".concat(this._save_time_interval_key);

        let save_time_interval = this._default_save_time_interval;

        $.ajax({url: settings_url, type: 'GET', contentType: "application/json"})
            // On failure inform user and stop
            .fail(error => this.informUser(
                2, "Fatal error: cannot get setting value: " + error))
            .done((response) => {
                if (response.value != 'None') {
                    save_time_interval = response.value;
                }
                else {
                    this.informUser(1, "Warning: save_time_interval was not set correctly. " +
                        "Default time value of " + this._default_save_time_interval.toString() + "will be set.");
                }
            });
        return save_time_interval;
    }

    /**
     * Set `_is_graph_changed` flag to true.
     *
     * An example application for this flag is to save the model whenever the
     * graph is changed.
     */
    graphChanged() {
        this._is_graph_changed = true;
    }

    /**
     * Setup an JS interval that check if the graph has changed and saveModel
     * if it does change.
     *
     * @param wait waiting time before actually saving the model
     */
    setupGraphChangeChecker(wait) {
        let model_id = $("#idaes-fs-name").data("flowsheetId");
        let flowsheet_url = "/fs?id=".concat(model_id);

        var graphChangedChecker = setInterval(() => {
            if (this._is_graph_changed) {
                this.saveModel(flowsheet_url, this.paper.graph);
                // reset flag
                this._is_graph_changed = false;
            }
        }, wait);
        return graphChangedChecker;
    }

    /**
     * Save the model value. Waiting time could be specified to
     * disable multiple redundant saves caused by a stream of events
     *
     * Changing cell positions & link vertices fire multiple events
     * subsequently. That's why we add waiting time before actually
     * saving the model.
     *
     * This sends a PUT to the server to save the current model value.
     *
     * @param url The HTTP server that is running in the Python process
     * @param model The model to save
     */
    saveModel(url, model) {
        let clientData = JSON.stringify(model.toJSON());
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
