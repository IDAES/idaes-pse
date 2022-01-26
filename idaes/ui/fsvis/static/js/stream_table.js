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
            // Also ignore the "Units" column header
            let column_header = columns[col];
            if (column_header !== "" && column_header !== "Units") {
                // If the column_header is Variable then we don't want the column to be right-aligned and we want the column to be pinned to the left so when the user scrolls the column scrolls with them
                if (column_header === "Variable") {
                    column_defs.push({
                        headerName: column_header,
                        field: column_header,
                        filter: 'agTextColumnFilter',
                        sortable: true,
                        resizable: true,
                        pinned: 'left',
                        cellRenderer: (params) => {
                            return '<span class="streamtable-variable">' + params.value + '</span>';
                        }
                    });
                }
                // If the column header isn't "Variable" then we assume that the contents of the column are numbers so they should be right aligned
                else {
                    column_defs.push({
                        headerName: column_header,
                        field: column_header,
                        filter: 'agTextColumnFilter',
                        sortable: true,
                        resizable: true,
                        cellStyle: {"text-align": "right"}
                    });
                }
                let list_item = document.createElement("li");
                let checkbox_item = document.createElement("div");
                checkbox_item.class = "checkbox";
                // checkbox_item.id = column_header + "-checkbox"
                checkbox_item.innerHTML = '<label class="fancy-checkbox"><input type="checkbox" value="' + column_header + '" id="' + column_header + '" checked><i class="fas fa-check checked"></i><i class="far fa-circle unchecked"></i>' + column_header + '</label>';
                list_item.appendChild(checkbox_item);
                hide_fields_list.appendChild(list_item);
            };
        };

        // Set the row data
        let variables = stream_table_data["index"];
        let data_arrays = stream_table_data["data"];
        let row_data = [];
        let variable_col = "Variable";
        for (let var_index in variables) {
            let row_object = {};
            let data = data_arrays[var_index];
            for (let col_index in columns) {
                if (columns[col_index] === "Units") {
                    if (data[col_index] && data[col_index] !== 'None') {
                        row_object[variable_col] = row_object[variable_col] + '<span class="streamtable-units">' + data[col_index].html + '</span>';
                    }
                    else {
                        row_object[variable_col] = row_object[variable_col] + '<span class="streamtable-units">&ndash;</span>';
                    }
                }
                else {
                    row_object[columns[col_index]] = data[col_index];
                }
            };
            row_data.push(row_object);
        };

        // let the grid know which columns and what data to use
        this._gridOptions = {
            columnDefs: column_defs,
            rowData: row_data,
            suppressColumnVirtualisation: true,
        };

        // Color the even rows grey
        this._gridOptions.getRowStyle = function(params) {
            if (params.node.rowIndex % 2 === 0) {
                return { background: "#f3f3f3" };
            }
        }

        // lookup the container we want the Grid to use
        let eGridDiv = document.querySelector('#stream-table-data');

        // create the grid passing in the div to use together with the columns & data we want to use
        new agGrid.Grid(eGridDiv, this._gridOptions);
        this._gridOptions.columnApi.autoSizeAllColumns();
    };

    setupEvents() {
        // This method sets up the event listeners for the table 

        // Set up the show/hide checkboxes for the Hide Field dropdown in the nav bar
        let hide_fields_list = document.querySelector("#hide-fields-list")
        let checkboxes = hide_fields_list.querySelectorAll("input[type=checkbox]");
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

        // Getting the main elements for the idaes canvas and the stream table
        // to be able to dispatch highlighting events to the streams existing
        // on paper and in the stream table
        let streamTable = document.querySelector('#stream-table-data');
        let idaesCanvas = document.querySelector('#idaes-canvas');

        // Registering listeners to the stream table to highlight the correct
        // streams in the stream table
        streamTable.addEventListener('HighlightStream', (event) => {
            var streamGridCells = streamTable.querySelectorAll(
                `[col-id=${event.detail.streamId}]`
            );
            streamGridCells.forEach((gridCell, index) => {
                if (gridCell.getAttribute('role') == 'columnheader') {
                    gridCell.classList.add('link-streamtable-hover-columnheader');
                }
                else if (index == streamGridCells.length - 1) {
                    gridCell.classList.add('link-streamtable-hover-lastrow');
                }
                else {
                    gridCell.classList.add('link-streamtable-hover');
                }
            });
        });

        // Registering listeners to idaes-canvas to remove the highlight from
        // the correct streams in the stream table
        streamTable.addEventListener('RemoveHighlightStream', (event) => {
            var streamGridCells = streamTable.querySelectorAll(
                `[col-id=${event.detail.streamId}]`
            );
            streamGridCells.forEach((gridCell) => {
                gridCell.classList.remove('link-streamtable-hover-columnheader');
                gridCell.classList.remove('link-streamtable-hover-lastrow');
                gridCell.classList.remove('link-streamtable-hover');
            });
        });

        let streamGridCells = document.querySelectorAll('[col-id]');
        streamGridCells.forEach((gridCell) => {
            // When the mouse hovers over a grid cell, the link as well as the
            // stream column that represents the correct stream will be highlighted.
            gridCell.addEventListener('mouseenter', function(event) {
                const highlightStreamEvent = new CustomEvent(
                    'HighlightStream',
                    {
                        detail: {
                            streamId: event.target.attributes['col-id'].value
                        }
                    }
                );
                streamTable.dispatchEvent(highlightStreamEvent);
                idaesCanvas.dispatchEvent(highlightStreamEvent);
            });

            // When the mouse leaves a grid cell, the link as well as the
            // stream column that represents the correct stream will remove
            // the highlighting feature.
            gridCell.addEventListener('mouseleave', function(event) {
                const removeHighlightStreamEvent = new CustomEvent(
                    'RemoveHighlightStream',
                    {
                        detail: {
                            streamId: event.target.attributes['col-id'].value
                        }
                    }
                );
                streamTable.dispatchEvent(removeHighlightStreamEvent);
                idaesCanvas.dispatchEvent(removeHighlightStreamEvent);
            });
        });
    };
};
