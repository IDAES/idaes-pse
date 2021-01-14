export class StreamTable {
    constructor(app, model) {
        this._app = app;
        this.initTable(model);
    };

    initTable(model) {
        // Clear the table first in case this is a call on refresh then call all of the methods to fill the table and setup the events
        this.clearTable();
        this.fillTable(model);
        this.setupEvents();
    }

    clearTable() {
        // Clear the table
        $("#hide-fields-list").empty();
        console.log(document.querySelector("#stream-table-data"))
        $("#stream-table-data").empty();
    }

    fillTable(model) {
        // This method fills the table with the streams and information from the model

        // Get the stream table data from the html
        let stream_table_data = model["model"]["stream_table"];

        // Get the hide fields list
        const hide_fields_list = document.querySelector("#hide-fields-list");

        // Specify the column headers
        let columns = stream_table_data["columns"];
        let column_defs = [];
        for (let col in columns) {
            // There is an empty column because of the way that the pandas dataframe was oriented so 
            // only add the columns that don't have an empty column header
            if (columns[col] !== "") {
                column_defs.push({headerName: columns[col], field: columns[col], filter: 'agTextColumnFilter', sortable: true});
                let list_item = document.createElement("li");
                let checkbox_item = document.createElement("div");
                checkbox_item.class = "checkbox";
                // checkbox_item.id = columns[col] + "-checkbox"
                checkbox_item.innerHTML = '<label><input type="checkbox" value="' + columns[col] + '" id="' + columns[col] + '" checked>' + columns[col] + '</label>';
                list_item.appendChild(checkbox_item);
                hide_fields_list.appendChild(list_item);
            };
        };

        // Set the row data
        let variables = stream_table_data["index"];
        let data_arrays = stream_table_data["data"];
        let row_data = [];
        for (let var_index in variables) {
            let row_object = {};
            let data = data_arrays[var_index];
            for (let col_index in columns) {
                row_object[columns[col_index]] = data[col_index];
            };
            row_data.push(row_object);
        };

        // let the grid know which columns and what data to use
        this._gridOptions = {
            columnDefs: column_defs,
            rowData: row_data,
        };

        // lookup the container we want the Grid to use
        let eGridDiv = document.querySelector('#stream-table-data');

        // create the grid passing in the div to use together with the columns & data we want to use
        new agGrid.Grid(eGridDiv, this._gridOptions);
    };

    setupEvents() {
        // This method sets up the event listeners for the table 

        // Set up the show/hide checkboxes for the Hide Field dropdown in the nav bar
        let checkboxes = document.querySelectorAll("input[type=checkbox]");
        // We need to save this to another variable temporarily to avoid collisions with this
        let app = this
        checkboxes.forEach(function(checkbox) {
            checkbox.addEventListener('change', function() {
                if (this.checked) {
                    app._gridOptions.columnApi.setColumnVisible(this.id, true)
                } 
                else {
                    app._gridOptions.columnApi.setColumnVisible(this.id, false)
                };
            });
        });

        // Set up an event listener to show and hide the table div 
        document.querySelector("#hide-table-btn").addEventListener("click", function() {
            let container = document.querySelector("#stream-table-container");
            let hide_fields_dropdown = document.querySelector("#hide-fields-dropdown");
            let export_csv_button = document.querySelector("#export-table-csv-btn");
            if (container.style.display === "none") {
                container.style.display = "block";
                hide_fields_dropdown.style.display = "block";
                export_csv_button.style.display = "block";
            } 
            else {
                container.style.display = "none";
                hide_fields_dropdown.style.display = "none";
                export_csv_button.style.display = "none";
            };
        });

        // Set up the table csv export
        document.querySelector("#export-table-csv-btn").addEventListener("click", function() {
            app._gridOptions.api.exportDataAsCsv({suppressQuotes: true});
        });
    };
};




