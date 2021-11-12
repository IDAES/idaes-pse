
export class JointJsCellConfig {
	constructor(cells_json) {
	    this._cells_json = cell_json;
	}
    
	get cells_json() {
	    return this._cells_json
	}
    
	set cell_json(cell_json) {
	    this._cells_json = cell_json;
	}
    
	parse() {
        var routing_fns = {}
        for (link in this._cells_json['routing_config']) {
            var routing_fn = null;
            // TODO: Implement for source as well
            if ('destination' in link) {
                routing_fn = function(vertices, opt, linkView) {
                    const sourceBBox = linkView.sourceBBox;
                    const a = linkView.getEndAnchor('source');
                    const b = linkView.getEndAnchor('target');
                    const minGap = link.destination.gap.distance;
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
            }

            routing_fns[link] = routing_fn;
        }


	}
	
    };
    