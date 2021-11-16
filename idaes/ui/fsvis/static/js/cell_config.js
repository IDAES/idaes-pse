
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
    