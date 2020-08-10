##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################

"""
Description
-----------
This mapping defines the ports for a given icon which makes the arcs point 
to the correct place on the icon

When a unit model or icon is added
----------------------------------
If it does not have a section already in this mapping:
Copy "cstr" and paste it at the bottom of link_position_mapping. Then, modify 
the "x" and "y" of each of the groups inside of the "group" dict in order to 
get the arcs to point to the appropriate icon edges.

If it does have a section in the mapping:
Modify the "x" and "y" of each of the groups inside of the "group" dict to 
get the arcs to point to the appropriate icon edges

Dependencies
------------
flowsheet_serializer depends on this to give the arcs a place to point to

If you do this it will break
----------------------------
Removing anything from inside of one of the groups (ie "in" or "out") in 
"groups" you will break the ports that allow for the arcs to point to the 
correct place on the icon

Note
----
The current implementation is simple and only handles one place of input 
and one place of output. I (Makayla) think that adding a new input or output 
will be easy in terms of the mapping (just add another group) but then you 
have to figure out the name of the port to point to/from when serializing 
the model
"""
link_position_mapping = {
    "default": {
        "groups": {
            "in":{
                "position":{
                    "name":"left",
                    "args":{
                       "x":15,
                       "y":0,
                       "dx":1,
                       "dy":1
                    }
                },
                "attrs": {
                    "rect": {
                        "stroke": '#000000',
                        'stroke-width': 0,
                        "width": 0,
                        "height": 0
                    }
                },
                "markup": '<g><rect/></g>'
            },
            "out": {
                "position":{
                    "name":"left",
                    "args":{
                       "x":48,
                       "y":45,
                       "dx":1,
                       "dy":1
                    }
                },
                "attrs": {
                    "rect": {
                        "stroke": '#000000',
                        'stroke-width': 0,
                        "width": 0,
                        "height": 0
                    }
                },
                "markup": '<g><rect/></g>'
            },
        },
        "items":[
           {
              "group":"in",
              "id":"in"
           },
           {
              "group":"out",
              "id":"out"
           }
        ]
    },
    "cstr": {
        "groups": {
            "in":{
                "position":{
                    "name":"left",
                    "args":{
                       "x":15,
                       "y":0,
                       "dx":1,
                       "dy":1
                    }
                },
                "attrs": {
                    "rect": {
                        "stroke": '#000000',
                        'stroke-width': 0,
                        "width": 0,
                        "height": 0
                    }
                },
                "markup": '<g><rect/></g>'
            },
            "out": {
                "position":{
                    "name":"left",
                    "args":{
                       "x":48,
                       "y":45,
                       "dx":1,
                       "dy":1
                    }
                },
                "attrs": {
                    "rect": {
                        "stroke": '#000000',
                        'stroke-width': 0,
                        "width": 0,
                        "height": 0
                    }
                },
                "markup": '<g><rect/></g>'
            },
        },
        "items":[
           {
              "group":"in",
              "id":"in"
           },
           {
              "group":"out",
              "id":"out"
           }
        ]
    },
    "equilibrium_reactor": {
        "groups": {
            "in":{
                "position":{
                    "name":"left",
                    "args":{
                       "x":5,
                       "y":10,
                       "dx":1,
                       "dy":1
                    }
                },
                "attrs": {
                    "rect": {
                        "stroke": '#000000',
                        'stroke-width': 0,
                        "width": 0,
                        "height": 0
                    }
                },
                "markup": '<g><rect/></g>'
            },
            "out": {
                "position":{
                    "name":"left",
                    "args":{
                       "x":45,
                       "y":45,
                       "dx":1,
                       "dy":1
                    }
                },
                "attrs": {
                    "rect": {
                        "stroke": '#000000',
                        'stroke-width': 0,
                        "width": 0,
                        "height": 0
                    }
                },
                "markup": '<g><rect/></g>'
            },
        },
        "items":[
           {
              "group":"in",
              "id":"in"
           },
           {
              "group":"out",
              "id":"out"
           }
        ]
    },
    "gibbs_reactor": {
        "groups": {
            "in":{
                "position":{
                    "name":"left",
                    "args":{
                       "x":5,
                       "y":10,
                       "dx":1,
                       "dy":1
                    }
                },
                "attrs": {
                    "rect": {
                        "stroke": '#000000',
                        'stroke-width': 0,
                        "width": 0,
                        "height": 0
                    }
                },
                "markup": '<g><rect/></g>'
            },
            "out": {
                "position":{
                    "name":"left",
                    "args":{
                       "x":45,
                       "y":45,
                       "dx":1,
                       "dy":1
                    }
                },
                "attrs": {
                    "rect": {
                        "stroke": '#000000',
                        'stroke-width': 0,
                        "width": 0,
                        "height": 0
                    }
                },
                "markup": '<g><rect/></g>'
            },
        },
        "items":[
           {
              "group":"in",
              "id":"in"
           },
           {
              "group":"out",
              "id":"out"
           }
        ]
    },
    "pfr": {
        "groups": {
            "in":{
                "position":{
                    "name":"left",
                    "args":{
                       "x":5,
                       "y":10,
                       "dx":1,
                       "dy":1
                    }
                },
                "attrs": {
                    "rect": {
                        "stroke": '#000000',
                        'stroke-width': 0,
                        "width": 0,
                        "height": 0
                    }
                },
                "markup": '<g><rect/></g>'
            },
            "out": {
                "position":{
                    "name":"left",
                    "args":{
                       "x":45,
                       "y":45,
                       "dx":1,
                       "dy":1
                    }
                },
                "attrs": {
                    "rect": {
                        "stroke": '#000000',
                        'stroke-width': 0,
                        "width": 0,
                        "height": 0
                    }
                },
                "markup": '<g><rect/></g>'
            },
        },
        "items":[
           {
              "group":"in",
              "id":"in"
           },
           {
              "group":"out",
              "id":"out"
           }
        ]
    },
    "stoichiometric_reactor": {
        "groups": {
            "in":{
                "position":{
                    "name":"left",
                    "args":{
                       "x":5,
                       "y":10,
                       "dx":1,
                       "dy":1
                    }
                },
                "attrs": {
                    "rect": {
                        "stroke": '#000000',
                        'stroke-width': 0,
                        "width": 0,
                        "height": 0
                    }
                },
                "markup": '<g><rect/></g>'
            },
            "out": {
                "position":{
                    "name":"left",
                    "args":{
                       "x":45,
                       "y":45,
                       "dx":1,
                       "dy":1
                    }
                },
                "attrs": {
                    "rect": {
                        "stroke": '#000000',
                        'stroke-width': 0,
                        "width": 0,
                        "height": 0
                    }
                },
                "markup": '<g><rect/></g>'
            },
        },
        "items":[
           {
              "group":"in",
              "id":"in"
           },
           {
              "group":"out",
              "id":"out"
           }
        ]
    },
    "flash": {
        "groups": {
            "in":{
                "position":{
                    "name":"left",
                    "args":{
                       "x":8,
                       "y":25,
                       "dx":1,
                       "dy":1
                    }
                },
                "attrs": {
                    "rect": {
                        "stroke": '#000000',
                        'stroke-width': 0,
                        "width": 0,
                        "height": 0
                    }
                },
                "markup": '<g><rect/></g>'
            },
            "top": {
                "position":{
                    "name":"top",
                    "args":{
                       "x":25,
                       "y":0,
                       "dx":1,
                       "dy":1
                    }
                },
                "attrs": {
                    "rect": {
                        "stroke": '#000000',
                        'stroke-width': 0,
                        "width": 0,
                        "height": 0
                    }
                },
                "markup": '<g><rect/></g>'
            },
            "bottom": {
                "position":{
                    "name":"bottom",
                    "args":{
                       "x":25,
                       "y":50,
                       "dx":1,
                       "dy":1
                    }
                },
                "attrs": {
                    "rect": {
                        "stroke": '#000000',
                        'stroke-width': 0,
                        "width": 0,
                        "height": 0
                    }
                },
                "markup": '<g><rect/></g>'
            },
        },
        "items":[
           {
              "group":"in",
              "id":"in"
           },
           {
              "group":"top",
              "id":"top"
           },
           {
              "group":"bottom",
              "id":"bottom"
           }
        ]
    },
    "mixer": {
        "groups": {
            "in":{
                "position":{
                    "name":"left",
                    "args":{
                       "x":2,
                       "y":25,
                       "dx":1,
                       "dy":1
                    }
                },
                "attrs": {
                    "rect": {
                        "stroke": '#000000',
                        'stroke-width': 0,
                        "width": 0,
                        "height": 0
                    }
                },
                "markup": '<g><rect/></g>'
            },

            "out": {
                "position":{
                    "name":"left",
                    "args":{
                       "x":48,
                       "y":25,
                       "dx":1,
                       "dy":1
                    }
                },
                "attrs": {
                    "rect": {
                        "stroke": '#000000',
                        'stroke-width': 0,
                        "width": 0,
                        "height": 0
                    }
                },
                "markup": '<g><rect/></g>'
            },
        },
        "items":[
           {
              "group":"in",
              "id":"in"
           },
           {
              "group":"out",
              "id":"out"
           }
        ]
    },
    "feed": {
        "groups": {
            "out":{
                "position": {
                    "name": "left",
                    "args": {
                       "x": 48,
                       "y": 25,
                       "dx": 1,
                       "dy": 1
                    }
                },
                "attrs": {
                    "rect": {
                        "stroke": '#000000',
                        'stroke-width': 0,
                        "width": 0,
                        "height": 0
                    }
                },
                "markup": '<g><rect/></g>'
            }
        },
        "items":[
           {
              "group":"out",
              "id":"out"
           },
        ]
    },
    "feed_flash": {
        "groups": {
            "in":{
                "position":{
                    "name":"left",
                    "args":{
                       "x":25,
                       "y":0,
                       "dx":1,
                       "dy":1
                    }
                },
                "attrs": {
                    "rect": {
                        "stroke": '#000000',
                        'stroke-width': 0,
                        "width": 0,
                        "height": 0
                    }
                },
                "markup": '<g><rect/></g>'
            },

            "out": {
                "position":{
                    "name":"left",
                    "args":{
                       "x":25,
                       "y":50,
                       "dx":1,
                       "dy":1
                    }
                },
                "attrs": {
                    "rect": {
                        "stroke": '#000000',
                        'stroke-width': 0,
                        "width": 0,
                        "height": 0
                    }
                },
                "markup": '<g><rect/></g>'
            },
        },
        "items":[
           {
              "group":"in",
              "id":"in"
           },
           {
              "group":"out",
              "id":"out"
           }
        ]
    },
    "product": {
        "groups": {
            "in":{
                "position":{
                    "name":"left",
                    "args":{
                       "x":2,
                       "y":25,
                       "dx":1,
                       "dy":1
                    }
                },
                "attrs": {
                    "rect": {
                        "stroke": '#000000',
                        'stroke-width': 0,
                        "width": 0,
                        "height": 0
                    }
                },
                "markup": '<g><rect/></g>'
            },
        },
        "items":[
           {
              "group":"in",
              "id":"in"
           },
        ]
    },
    "separator": {
        "groups": {
            "in": {
                "position": {
                    "name":"left",
                    "args": {
                       "x": 2,
                       "y": 25,
                       "dx": 1,
                       "dy": 1
                    }
                },
                "attrs": {
                    "rect": {
                        "stroke": '#000000',
                        'stroke-width': 0,
                        "width": 0,
                        "height": 0
                    }
                },
                "markup": '<g><rect/></g>'
            },

            "out": {
                "position":{
                    "name":"right",
                    "args":{
                       "x":48,
                       "y":25,
                       "dx":1,
                       "dy":1
                    }
                },
                "attrs": {
                    "rect": {
                        "stroke": '#000000',
                        'stroke-width': 0,
                        "width": 0,
                        "height": 0
                    }
                },
                "markup": '<g><rect/></g>'
            },
        },
        "items":[
           {
              "group":"in",
              "id":"in"
           },
           {
              "group":"out",
              "id":"out"
           }
        ]
    },
    "heater": {
        "groups": {
            "in": {
                "position": {
                    "name": "left",
                    "args": {
                       "x": 6,
                       "y": 25,
                       "dx": 1,
                       "dy": 1
                    }
                },
                "attrs": {
                    "rect": {
                        "stroke": '#000000',
                        'stroke-width': 0,
                        "width": 0,
                        "height": 0
                    }
                },
                "markup": '<g><rect/></g>'
            },

            "out": {
                "position":{
                    "name":"left",
                    "args":{
                       "x":43,
                       "y":25,
                       "dx":1,
                       "dy":1
                    }
                },
                "attrs": {
                    "rect": {
                        "stroke": '#000000',
                        'stroke-width': 0,
                        "width": 0,
                        "height": 0
                    }
                },
                "markup": '<g><rect/></g>'
            },
        },
        "items":[
           {
              "group":"in",
              "id":"in"
           },
           {
              "group":"out",
              "id":"out"
           }
        ]
    },
    "pressure_changer": {
        "groups": {
            "in":{
                "position":{
                    "name":"left",
                    "args":{
                       "x":2,
                       "y":25,
                       "dx":1,
                       "dy":1
                    }
                },
                "attrs": {
                    "rect": {
                        "stroke": '#000000',
                        'stroke-width': 0,
                        "width": 0,
                        "height": 0
                    }
                },
                "markup": '<g><rect/></g>'
            },

            "out": {
                "position":{
                    "name":"left",
                    "args":{
                       "x":48,
                       "y":25,
                       "dx":1,
                       "dy":1
                    }
                },
                "attrs": {
                    "rect": {
                        "stroke": '#000000',
                        'stroke-width': 0,
                        "width": 0,
                        "height": 0
                    }
                },
                "markup": '<g><rect/></g>'
            },
        },
        "items":[
           {
              "group":"in",
              "id":"in"
           },
           {
              "group":"out",
              "id":"out"
           }
        ]
    },
    "heat_exchanger": {
        "groups": {
            "in":{
                "position":{
                    "name":"left",
                    "args":{
                       "x":2,
                       "y":25,
                       "dx":1,
                       "dy":1
                    }
                },
                "attrs": {
                    "rect": {
                        "stroke": '#000000',
                        'stroke-width': 0,
                        "width": 0,
                        "height": 0
                    }
                },
                "markup": '<g><rect/></g>'
            },

            "out": {
                "position":{
                    "name":"left",
                    "args":{
                       "x":48,
                       "y":25,
                       "dx":1,
                       "dy":1
                    }
                },
                "attrs": {
                    "rect": {
                        "stroke": '#000000',
                        'stroke-width': 0,
                        "width": 0,
                        "height": 0
                    }
                },
                "markup": '<g><rect/></g>'
            },
        },
        "items":[
           {
              "group":"in",
              "id":"in"
           },
           {
              "group":"out",
              "id":"out"
           }
        ]
    },
    # TODO When we get to using different icons we need to figure out how
    #  to deal with mutiple inlets and outlets in different places
    "heat_exchanger1d": {
        "groups": {
            "in":{
                "position":{
                    "name":"left",
                    "args":{
                       "x":2,
                       "y":25,
                       "dx":1,
                       "dy":1
                    }
                },
                "attrs": {
                    "rect": {
                        "stroke": '#000000',
                        'stroke-width': 0,
                        "width": 0,
                        "height": 0
                    }
                },
                "markup": '<g><rect/></g>'
            },

            "out": {
                "position":{
                    "name":"left",
                    "args":{
                       "x":48,
                       "y":25,
                       "dx":1,
                       "dy":1
                    }
                },
                "attrs": {
                    "rect": {
                        "stroke": '#000000',
                        'stroke-width': 0,
                        "width": 0,
                        "height": 0
                    }
                },
                "markup": '<g><rect/></g>'
            },
        },
        "items":[
           {
              "group":"in",
              "id":"in"
           },
           {
              "group":"out",
              "id":"out"
           }
        ]
    },
    "packed_column": {
        "groups": {
            "in":{
                "position":{
                    "name":"left",
                    "args":{
                       "x":48,
                       "y":10,
                       "dx":1,
                       "dy":1
                    }
                },
                "attrs": {
                    "rect": {
                        "stroke": '#000000',
                        'stroke-width': 0,
                        "width": 0,
                        "height": 0
                    }
                },
                "markup": '<g><rect/></g>'
            },

            "out": {
                "position":{
                    "name":"left",
                    "args":{
                       "x":48,
                       "y":40,
                       "dx":1,
                       "dy":1
                    }
                },
                "attrs": {
                    "rect": {
                        "stroke": '#000000',
                        'stroke-width': 0,
                        "width": 0,
                        "height": 0
                    }
                },
                "markup": '<g><rect/></g>'
            },
        },
        "items":[
           {
              "group":"in",
              "id":"in"
           },
           {
              "group":"out",
              "id":"out"
           }
        ],
    },
    "tray_column": {
        "groups": {
            "in":{
                "position":{
                    "name":"left",
                    "args":{
                       "x":48,
                       "y":10,
                       "dx":1,
                       "dy":1
                    }
                },
                "attrs": {
                    "rect": {
                        "stroke": '#000000',
                        'stroke-width': 0,
                        "width": 0,
                        "height": 0
                    }
                },
                "markup": '<g><rect/></g>'
            },

            "out": {
                "position":{
                    "name":"left",
                    "args":{
                       "x":48,
                       "y":40,
                       "dx":1,
                       "dy":1
                    }
                }
            },
        },
        "items":[
           {
              "group":"in",
              "id":"in"
           },
           {
              "group":"out",
              "id":"out"
           }
        ],
    },
    "default": {
        "groups": {
            "in":{
                "position":{
                    "name":"left",
                    "args":{
                       "x":2,
                       "y":0,
                       "dx":1,
                       "dy":1
                    }
                },
                "attrs": {
                    "rect": {
                        "stroke": '#000000',
                        'stroke-width': 0,
                        "width": 0,
                        "height": 0
                    }
                },
                "markup": '<g><rect/></g>'
            },

            "out": {
                "position":{
                    "name":"left",
                    "args":{
                       "x":48,
                       "y":50,
                       "dx":1,
                       "dy":1
                    }
                },
                "attrs": {
                    "rect": {
                        "stroke": '#000000',
                        'stroke-width': 0,
                        "width": 0,
                        "height": 0
                    }
                },
                "markup": '<g><rect/></g>'
            },
        },
        "items":[
           {
              "group":"in",
              "id":"in"
           },
           {
              "group":"out",
              "id":"out"
           }
        ]
    }
}
