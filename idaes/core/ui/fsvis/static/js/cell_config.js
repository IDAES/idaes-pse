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
 * 
 * Author: Abdelrahman Elbashandy <aaelbashandy@lbl.gov>
 */

/**
 * Interpret cell configuration specified by each icon chosen in the model.
 * Configuring cell parameters such as the routing paths that are taken by
 * links connecting unit models.
 */
export class JointJsCellConfig {
	constructor(model) {
	    this._model = model;
	}
    
	get model() {
	    return this._model
	}
    
	set model(model) {
	    this._model = model;
	}

    /**
     * Finding the correct cell index based on the given cell name 'cellName'
     */
    findCellIndex(cellName, cellType) {
        for (let i = 0; i < this.model['cells'].length; i++) {
            const cell = this.model['cells'][i];
            if (cell.id == cellName && cell.type == cellType) {
                return i;
            }
        }
        // If an index is not returned, that means the link was not found
        throw new Error(`Link with linkName: ${cellName} was not found`);
    }

    /**
     * Generate a custom function that handles the router 'gap' option.
     * The 'gap' option is specified by the users to choose the two vertices
     * that the link will take to connect the unit models (elements).
     */
    routerGapFnFactory(gap) {
        // Link router functions in JointJS execute automatically when their
        // elements that they are connected to rotate.
        var router_fn = (vertices, opt, linkView) => {
            const src = linkView.getEndAnchor('source');
            const tgt = linkView.getEndAnchor('target');

            // Calculates the original points that the link has to go through
            // to give the effect of having a gap.
            const p1 = new g.Point(
                src.x + gap.source.x,
                src.y + gap.source.y
            );
            const p2 = new g.Point(
                tgt.x + gap.destination.x,
                tgt.y + gap.destination.y
            );

            // Getting the current angles of the icons to calculate the right
            // position for the original points calculated above
            const src_angle = linkView.getEndView('source').model.angle();
            const tgt_angle = linkView.getEndView('target').model.angle();

            // Flip direction if perpendicular
            const src_orth = src_angle % 180 === 0 ? 1 : -1;
            const tgt_orth = tgt_angle % 180 === 0 ? 1 : -1;

            // Rotate the vector based on the angle of the element
            const src_rotated = p1.rotate(src, src_orth * src_angle);
            const tgt_rotated = p2.rotate(tgt, tgt_orth * tgt_angle);

            // TODO: Research if there's a better pre-implemented router
            // function than the manhattan function.
            return joint.routers.manhattan([src_rotated, ...vertices, tgt_rotated], opt, linkView);
        }

        return router_fn;
    }

    /**
     * Read the routing config for each link in jointjs and create a custom
     * routing function for it based on the routing configuration.
     */
	processRoutingConfig() {
        const src = "source";
        const dest = "destination";

        var routing_config = this._model['routing_config'];
        for (var linkName in routing_config) {
            // Create routing function
            var routing_fn = this.routerGapFnFactory(
                routing_config[linkName].cell_config.gap
            );

            // Because the jointjs object expects all the cell objects to exist
            // in an array, then we gotta find the right index for that link.
            const cell_index = this.findCellIndex(linkName, "standard.Link");

            // Assign the routing function to the right cell/link
            this._model['cells'][cell_index].router = routing_fn;
        }
        return this._model;
	}
};
