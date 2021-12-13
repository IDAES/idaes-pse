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

    // TODO: Implement a function that find the orthogonal path from a
    // source to a target.
    // Orthogonal Connector Routing:
    // https://link.springer.com/content/pdf/10.1007%2F978-3-642-11805-0_22.pdf
    //
    // routerGapFnFactory(linkEnds, direction, minGap) {
    //     var router_fn = null;
    //     switch(direction) {
    //         case "left":
    //             router_fn = (vertices, opt, linkView) => {
    //                 const a = linkView.getEndAnchor('source');
    //                 const b = linkView.getEndAnchor('target');
    //                 const minGap = link.destination.gap.distance;
    //                 const x1 = Math.min(b.x - minGap, (a.x + b.x) / 2);
    //                 const p1 = {
    //                     x: x1,
    //                     y: a.y
    //                 };
    //                 const p2 = {
    //                     x: x1,
    //                     y: b.y
    //                 };
    //                 return [p1, ...vertices, p2];
    //             }
    //             break;
    //         default:
    //             throw Error('Unknown direction for Link routing');
    //     }

    //     return router_fn;
    // }

    /**
     * Generate a custom function that handles the router 'gap' option.
     * The 'gap' option is specified by the users to choose the paths that the
     * link will take to connect unit models.
    */
    routerGapFnFactory(direction, gap) {
        var router_fn = null;
        switch(direction) {
            case "left":
                router_fn = (vertices, opt, linkView) => {
                    const a = linkView.getEndAnchor('source');
                    const b = linkView.getEndAnchor('target');
                    const minGap = gap
                    const x1 = Math.min(b.x - minGap, (a.x + b.x) / 2);
                    const p1 = {
                        x: x1,
                        y: a.y
                    };
                    const p2 = {
                        x: x1,
                        y: b.y
                    };
                    return [p1, ...vertices, p2];
                }
                break;
            default:
                throw Error('Unsupported direction for Link routing');
        }

        return router_fn;
    }

	processRoutingConfig() {
        var routing_config = this._model['routing_config'];
        for (var link in routing_config) {
            var routing_fn = null;
            // TODO: Implement for source as well
            if ('destination' in routing_config[link]) {
                routing_fn = this.routerGapFnFactory(
                    routing_config[link].destination.gap.direction,
                    routing_config[link].destination.gap.distance
                );

                this._model['cells'][routing_config[link].cell_index].router = routing_fn;
            }
        }
        return this._model;
	}
	
    };
    