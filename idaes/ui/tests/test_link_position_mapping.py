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
import pytest

from idaes.ui import link_position_mapping


@pytest.mark.unit
@pytest.mark.parametrize(
    "test_input,expected",
    [
        (
            "cstr",
            {
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
        ),
        (
            "equilibrium_reactor",
            {
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
        ),
        (
            "gibbs_reactor",
            {
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
        ),
        (
            "pfr",
            {
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
        ),
        (
            "stoichiometric_reactor",
            {
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
        ),
        (
            "flash",
            {
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
        ),
        (
            "mixer",
            {
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
        ),
        (
            "feed",
            {
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
        ),
        (
            "feed_flash",
            {
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
        ),
        (
            "product",
            {
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
        ),
        (
            "separator",
            {
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
        ),
        (
            "heater",
            {
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
        ),
        (
            "pressure_changer",
            {
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
        ),
        (
            "heat_exchanger",
            {
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
        ),
        (
            "heat_exchanger1d",
            {
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
        ),
        (
            "packed_column",
            {
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
        ),
        (
            "tray_column",
            {
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
        ),
    ],
)
@pytest.mark.unit
def test_link_position_mapping(test_input, expected):
    assert link_position_mapping.link_position_mapping[test_input] == expected
