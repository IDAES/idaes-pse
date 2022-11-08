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
 * StreamTable class is responsible to handle creating the Stream Table for
 * the flowsheet that is displayed on the Idaes paper graph. This class does
 * the data filling and data styling for each cell in the grid.
 */
export class StreamTable {
    // Variable Types
    UNFIXED = 'unfixed';
    FIXED = 'fixed';
    PARAMETER = 'parameter';
    EXPRESSION = 'expression';

    // Brushing Event handlers
    highlightFn = null;
    removeHighlightFn = null;
    gridCellMouseEnterFn = null;
    gridCellMouseLeaveFn = null;


    constructor(app, model) {
        this._app = app;
        // Define brushing event handlers
        this.defineTableBrushingFns();

        this.initTable(model);

        // Keeping track of existing variable types. e.g fixed, Parameter, Expression
        this.existing_var_types = new Set();
    };

    initTable(model) {
        // Clear the table first in case this is a call on refresh then call all of the methods to fill the table and setup the events
        this.clearTable();
        this.emptyVarTypesPanel();
        this.fillTable(model);
        this.setupEvents();
    }

    clearTable() {
        // Clear the table
        $("#hide-fields-list").empty();
        $("#stream-table-data").empty();
    }

    /**
     * Clear list of existing variable types
     */
    emptyVarTypesPanel () {
        this.existing_var_types = new Set();
        const var_types_panel = document.querySelector('#existing-variable-types');
        var_types_panel.innerHTML = "";
    }

    /**
     * Create a panel for each variable type that exists in the Stream Table
     */
    fillVarTypesPanel() {
        const var_types_panel = document.querySelector('#existing-variable-types');
        const stream_table_class = 'streamtable-vartype-element';

        // Adding header
        if (this.existing_var_types.has(this.FIXED) ||
            this.existing_var_types.has(this.PARAMETER) ||
            this.existing_var_types.has(this.EXPRESSION)) {
            const header_vartype = document.createElement('p');
            header_vartype.innerHTML = 'Annotated Variable Types:';
            header_vartype.className = stream_table_class;
            var_types_panel.appendChild(header_vartype);
        }
        // Adding each type
        this.existing_var_types.forEach(var_type => {
            switch (var_type) {
                case this.UNFIXED:
                    // This will execute once since this.existing_var_types is a set
                    console.debug(`Unfixed variables don't have a visual indicator`);
                    break;
                case this.FIXED:
                case this.PARAMETER:
                case this.EXPRESSION:
                    const elem_vartype = document.createElement('span'); // Parent node
                    elem_vartype.className = stream_table_class;

                    // Create dot with the right color and the right variable type text
                    const elem_dot = document.createElement('span');
                    const elem_text = document.createElement('span');
                    elem_text.className = 'streamtable-vartype-text';
                    elem_dot.className = `streamtable-vartype-${var_type}`;
                    elem_dot.title = var_type;
                    elem_text.innerHTML = var_type;
                    elem_vartype.appendChild(elem_dot);
                    elem_vartype.appendChild(elem_text);
                    var_types_panel.appendChild(elem_vartype);
                    break;
                default:
                    console.warn(`Couldn't identify Variable type: ${data[col_index]}`);
            };
        });
    }

    /**
     * This method fills the table with the streams and information from the model
     */
    fillTable(model) {
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
            if (column_header !== "" && column_header !== "Units" && !column_header.includes("_vartype")) {
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
                            return '<span class="streamtable-cell">' + params.value + '</span>';
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
                        cellRenderer: (params) => {
                            return '<span class="streamtable-cell">' + params.value + '</span>';
                        }
                    });
                    let list_item = document.createElement("li");
                    let checkbox_item = document.createElement("div");
                    checkbox_item.class = "checkbox";
                    // checkbox_item.id = column_header + "-checkbox"
                    checkbox_item.innerHTML = '<label class="fancy-checkbox"><input type="checkbox" value="' + column_header + '" id="' + column_header + '" checked><i class="fas fa-check checked"></i><i class="far fa-circle unchecked"></i>' + column_header + '</label>';
                    list_item.appendChild(checkbox_item);
                    hide_fields_list.appendChild(list_item);
                }
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
                    if (data[col_index] && data[col_index].html) {
                        row_object[variable_col] = row_object[variable_col] + '<span class="streamtable-units">' + data[col_index].html + '</span>';
                    }
                    else {
                        row_object[variable_col] = row_object[variable_col] + '<span class="streamtable-units">&ndash;</span>';
                    }
                }
                else if (columns[col_index] === "Variable") {
                    row_object[columns[col_index]] = data[col_index];
                }
                else {
                    var [value, type] = data[col_index];
                    let cell_style = "";
                    switch (type) {
                        case this.UNFIXED:
                            this.existing_var_types.add(type);
                            break;
                        case this.FIXED:
                        case this.PARAMETER:
                        case this.EXPRESSION:
                            this.existing_var_types.add(type);
                            cell_style = `<span class="streamtable-vartype-${type}" style="margin-top: 7%;" title="${type}"></span>`;
                            break;
                        default:
                            console.warn(`Couldn't identify Variable type: ${type}`);
                    };
                    row_object[columns[col_index]] = cell_style + '<span class="streamtable-variable-value">' + value + '</span>';
                }
            };
            row_data.push(row_object);
        };

        // Fill the Variable Types panel
        this.fillVarTypesPanel();

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

    /**
     * Define event handlers and save them as objects to be able
     * to remove these events later if a column in the Stream Table is removed
     */
    defineTableBrushingFns() {
        // Getting the main elements for the idaes canvas and the stream table
        // to be able to dispatch highlighting events to the streams existing
        // on paper and in the stream table
        let streamTable = document.querySelector('#stream-table-data');
        let idaesCanvas = document.querySelector('#idaes-canvas');

        // Function to highlight a stream table column
        this.highlightFn = (event) => {
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
        }

        // Function to undo the highlighting of a stream table column
        this.removeHighlightFn = (event) => {
            var streamGridCells = streamTable.querySelectorAll(
                `[col-id=${event.detail.streamId}]`
            );
            streamGridCells.forEach((gridCell) => {
                gridCell.classList.remove('link-streamtable-hover-columnheader');
                gridCell.classList.remove('link-streamtable-hover-lastrow');
                gridCell.classList.remove('link-streamtable-hover');
            });
        }

        // Function to trigger the right highlighting events when the mouse is
        // hovering over a cell in the stream table
        this.gridCellMouseEnterFn = (event) => {
            if (document.querySelector("#view-stream-highlight-btn").checked) {
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
            }
        }

        // Function to trigger undoing the highlighting events when the mouse
        // leaves the hovered upon cell in the stream table
        this.gridCellMouseLeaveFn = (event) => {
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
        }
    }

    /**
     * Register Stream Brushing when hovering over a grid cell
     */
    registerTableBrushing() {
        // Getting the main elements for the idaes canvas and the stream table
        // to be able to dispatch highlighting events to the streams existing
        // on paper and in the stream table
        let streamTable = document.querySelector('#stream-table-data');
        let idaesCanvas = document.querySelector('#idaes-canvas');

        let streamGridCells = document.querySelectorAll('[col-id]');

        // Cleaning up events
        streamTable.removeEventListener('HighlightStream', this.highlightFn);
        streamTable.removeEventListener('RemoveHighlightStream', this.removeHighlightFn);
        streamGridCells.forEach((gridCell) => {
            // When the mouse hovers over a grid cell, the link as well as the
            // stream column that represents the correct stream will be highlighted.
            gridCell.removeEventListener('mouseenter', this.gridCellMouseEnterFn);

            // When the mouse leaves a grid cell, the link as well as the
            // stream column that represents the correct stream will remove
            // the highlighting feature.
            gridCell.removeEventListener('mouseleave', this.gridCellMouseLeaveFn);
        });

        // Registering listeners to the stream table to highlight the correct
        // streams in the stream table
        streamTable.addEventListener('HighlightStream', this.highlightFn);

        // Registering listeners to idaes-canvas to remove the highlight from
        // the correct streams in the stream table
        streamTable.addEventListener('RemoveHighlightStream', this.removeHighlightFn);

        streamGridCells.forEach((gridCell) => {
            // When the mouse hovers over a grid cell, the link as well as the
            // stream column that represents the correct stream will be highlighted.
            gridCell.addEventListener('mouseenter', this.gridCellMouseEnterFn);

            // When the mouse leaves a grid cell, the link as well as the
            // stream column that represents the correct stream will remove
            // the highlighting feature.
            gridCell.addEventListener('mouseleave', this.gridCellMouseLeaveFn);
        });
    }

    /**
     * This method sets up the event listeners for the table
     */
    setupEvents() {
        // Set up the show/hide checkboxes for the Hide Field dropdown in the nav bar
        let hide_fields_list = document.querySelector("#hide-fields-list")
        let checkboxes = hide_fields_list.querySelectorAll("input[type=checkbox]");
        // We need to save this to another variable temporarily to avoid collisions with this
        let app = this
        checkboxes.forEach(function(checkbox) {
            checkbox.addEventListener('change', function() {
                if (this.checked) {
                    app._gridOptions.columnApi.setColumnVisible(this.id, true)
                    app.registerTableBrushing()
                }
                else {
                    app._gridOptions.columnApi.setColumnVisible(this.id, false)
                };
            });
        });

        this.registerTableBrushing()
    };
};
