#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
from typing import Dict


class UnitModelIcon:
    """Represents icon display information for a given unit model.
    """

    #: Name of default unit_model to use
    DEFAULT = "default"

    def __init__(self, unit_model: str = None, default: str = DEFAULT):
        """Construct with a given unit model type name.

        Args:
            unit_model: Name of the unit model type. If not given, use value in attribute `DEFAULT`
            default: If the given `unit_model` is not found, use this one; but then if this is falsy raise a KeyError

        Raises:
            KeyError if unit_model name is not found and `default` arg is falsy (e.g. None or "")
        """
        if unit_model is None:
            unit_model = self.DEFAULT
        self._model = unit_model
        try:
            self._model_details = self._mapping[unit_model]
        except KeyError:
            if not default:
                raise
            self._model_details = self._mapping[self.DEFAULT]
        self._pos = self._build_link_positions()

    @property
    def icon(self) -> str:
        """Get the name of the icon.
        """
        # return self._info[0]
        return self._model_details['image']

    @property
    def link_positions(self) -> Dict:
        """Get the link positions.

        Example result::

            {
                "port_groups": {
                    "in": {
                        "position": {
                            "name": "left",
                            "args": {"x": 15, "y": 0, "dx": 1, "dy": 1},
                        },
                        "attrs": {
                            "rect": {
                                "stroke": "#000000",
                                "stroke-width": 0,
                                "width": 0,
                                "height": 0,
                            }
                        },
                        "markup": "<g><rect/></g>",
                    },
                    "out": {
                        "position": {
                            "name": "left",
                            "args": {"x": 48, "y": 45, "dx": 1, "dy": 1},
                        },
                        "attrs": {
                            "rect": {
                                "stroke": "#000000",
                                "stroke-width": 0,
                                "width": 0,
                                "height": 0,
                            }
                        },
                        "markup": "<g><rect/></g>",
                    },
                },
                "items": []
            }

        Returns:
            The link position (see example result)
        """
        return self._pos
    
    @property
    def routing_config(self) -> Dict:
        """Get the Unit model routing config to be used to add jointjs vertices
        for layout control within the created graph.

        Example result::

            {
                "in": {
                    "gap": {
                        "direction": "left",
                        "distance": 20
                    }
                }
            }

        Returns:
            The routing configuration (see example result)
        """
        if "routing_config" in self._model_details:
            return self._model_details["routing_config"]
        else:
            return {}

    def _build_link_positions(self) -> Dict:
        """Fill in boilerplate based on raw info and place built value in class cache.
        """
        # build link positions from info
        groups, items = {}, []
        for group_name, group_config in self._model_details["port_groups"].items():
            groups[group_name] = group_config
            groups[group_name].update({
                "attrs": {
                    "rect": {
                        "stroke": "#000000",
                        "stroke-width": 0,
                        "width": 0,
                        "height": 0,
                    }
                },
                "markup": "<g><rect/></g>",
            })

        # set new link positions attr and place in cache
        positions = {"groups": groups, "items": []}
        return positions

    # === Data ===

    # Name is unit name, value is the information for its icon
    # Value is a tuple: icon image, and one or more position tuples: (group [in/out], name [side], (x, y, dx, dy))
    # Notes for updating:
    #  - Use 'cstr' as your template for new entries
    #  - Do not remove in/out entries in existing entries, or arcs won't connect
    # TODO: Move this mapping to its own directory/files.
    _mapping = {
        "cstr": {
            "image": "reactor_c.svg",
            "port_groups": {
                "in": {
                    "position": {
                        "name": "left",
                        "args": {
                            "x": 15,
                            "y": 0,
                            "dx": 1,
                            "dy": 1
                        }
                    }
                },
                "out": {
                    "position": {
                        "name": "left",
                        "args": {
                            "x": 48,
                            "y": 45,
                            "dx": 1,
                            "dy": 1
                        }
                    }
                }
            }
        },
        "flash": {
            "image": "flash.svg",
            "port_groups": {
                "bottom": {
                    "position": {
                        "name": "bottom",
                        "args": {
                            "x": 25,
                            "y": 50,
                            "dx": 1,
                            "dy": 1
                        }
                    }
                },
                "in": {
                    "position": {
                        "name": "left",
                        "args": {
                            "x":  8,
                            "y": 25,
                            "dx": 1,
                            "dy": 1
                        }
                    }
                },
                "top": {
                    "position": {
                        "name": "top",
                        "args": {
                            "x": 25,
                            "y":  0,
                            "dx": 1,
                            "dy": 1
                        }
                    }
                },
                # added by AR
                "out": {
                    "position": {
                        "name": "left",
                        "args": {
                            "x": 45,
                            "y": 45,
                            "dx": 1,
                            "dy": 1
                        }
                    }
                },
                #end
            }
        },
        "gibbs_reactor": {
            "image": "reactor_g.svg",
            "port_groups": {
                "in": {
                    "position": {
                        "name": "left",
                        "args": {
                            "x":  5,
                            "y": 10,
                            "dx": 1,
                            "dy": 1
                        }
                    }
                },
                "out": {
                    "position": {
                        "name": "left",
                        "args": {
                            "x": 45,
                            "y": 45,
                            "dx": 1,
                            "dy": 1
                        }
                    }
                }
            }
        },
        "heat_exchanger": {
            "image": "heat_exchanger_1.svg",
            "port_groups": {
                "in": {
                    "position": {
                        "name": "left",
                        "args": {
                            "x":  2,
                            "y": 25,
                            "dx": 1,
                            "dy": 1
                        }
                    }
                },
                "out": {
                    "position": {
                        "name": "left",
                        "args": {
                            "x": 48,
                            "y": 25,
                            "dx": 1,
                            "dy": 1
                        }
                    }
                }
            }
        },
        "heater": {
            "image": "heater_2.svg",
            "port_groups": {
                "in": {
                    "position": {
                        "name": "left",
                        "args": {
                            "x":  6,
                            "y": 25,
                            "dx": 1,
                            "dy": 1
                        }
                    }
                },
                "out": {
                    "position": {
                        "name": "left",
                        "args": {
                            "x": 43,
                            "y": 25,
                            "dx": 1,
                            "dy": 1
                        }
                    }
                }
            }
        },
        "heat_exchanger_1D": {
            "image": "heat_exchanger_1.svg",
            "port_groups": {
                "in": {
                    "position": {
                        "name": "left",
                        "args": {
                            "x": 15,
                            "y":  0,
                            "dx": 1,
                            "dy": 1
                        }
                    }
                },
                "out": {
                    "position": {
                        "name": "left",
                        "args": {
                            "x": 48,
                            "y": 45,
                            "dx": 1,
                            "dy": 1
                        }
                    }
                }
            }
        },
        "mixer": {
            "image": "mixer.svg",
            "port_groups": {
                "in": {
                    "position": {
                        "name": "left",
                        "args": {}
                    }
                },
                "out": {
                    "position": {
                        "name": "left",
                        "args": {
                            "x": 48,
                            "y": 25,
                            "dx": 1,
                            "dy": 1
                        }
                    }
                }
            },
            "routing_config": {
                "in": {
                    "gap": {
                        "direction": "left",
                        "distance": 30
                    }
                }
            }
        },
        "plug_flow_reactor": {
            "image": "reactor_pfr.svg",
            "port_groups": {
                "in": {
                    "position": {
                        "name": "left",
                        "args": {
                            "x": 15,
                            "y":  0,
                            "dx": 1,
                            "dy": 1
                        }
                    }
                },
                "out": {
                    "position": {
                        "name": "left",
                        "args": {
                            "x": 48,
                            "y": 45,
                            "dx": 1,
                            "dy": 1
                        }
                    }
                }
            }
        },
        "pressure_changer": {
            "image": "compressor.svg",
            "port_groups": {
                "in": {
                    "position": {
                        "name": "left",
                        "args": {
                            "x":  2,
                            "y": 25,
                            "dx": 1,
                            "dy": 1
                        }
                    }
                },
                "out": {
                    "position": {
                        "name": "left",
                        "args": {
                            "x": 48,
                            "y": 25,
                            "dx": 1,
                            "dy": 1
                        }
                    }
                }
            }
        },
        "separator": {
            "image": "splitter.svg",
            "port_groups": {
                "in": {
                    "position": {
                        "name": "left",
                        "args": {
                            "x":  2,
                            "y": 25,
                            "dx": 1,
                            "dy": 1
                        }
                    }
                },
                "out": {
                    "position": {
                        "name": "right",
                        "args": {
                            "x": 48,
                            "y": 25,
                            "dx": 1,
                            "dy": 1
                        }
                    }
                }
            }
        },
        "stoichiometric_reactor": {
            "image": "reactor_s.svg",
            "port_groups": {
                "in": {
                    "position": {
                        "name": "left",
                        "args": {
                            "x":  5,
                            "y": 10,
                            "dx": 1,
                            "dy": 1
                        }
                    }
                },
                "out": {
                    "position": {
                        "name": "left",
                        "args": {
                            "x": 45,
                            "y": 45,
                            "dx": 1,
                            "dy": 1
                        }
                    }
                }
            }
        },
        "equilibrium_reactor": {
            "image": "reactor_e.svg",
            "port_groups": {
                "in": {
                    "position": {
                        "name": "left",
                        "args": {
                            "x":  5,
                            "y": 10,
                            "dx": 1,
                            "dy": 1
                        }
                    }
                },
                "out": {
                    "position": {
                        "name": "left",
                        "args": {
                            "x": 45,
                            "y": 45,
                            "dx": 1,
                            "dy": 1
                        }
                    }
                }
            }
        },
        "feed": {
            "image": "feed.svg",
            "port_groups": {
                "out": {
                    "position": {
                        "name": "left",
                        "args": {
                            "x": 48,
                            "y": 25,
                            "dx": 1,
                            "dy": 1
                        }
                    }
                }
            }
        },
        "product": {
            "image": "product.svg",
            "port_groups": {
                "in": {
                    "position": {
                        "name": "left",
                        "args": {
                            "x":  2,
                            "y": 25,
                            "dx": 1,
                            "dy": 1
                        }
                    }
                }
            }
        },
        "feed_flash": {
            "image": "feed.svg",
            "port_groups": {
                "in": {
                    "position": {
                        "name": "left",
                        "args": {
                            "x": 25,
                            "y":  0,
                            "dx": 1,
                            "dy": 1
                        }
                    }
                },
                "out": {
                    "position": {
                        "name": "left",
                        "args": {
                            "x": 25,
                            "y": 50,
                            "dx": 1,
                            "dy": 1
                        }
                    }
                }
            }
        },
        "statejunction": {
            "image": "NONE",
            "port_groups": {
                "in": {
                    "position": {
                        "name": "left",
                        "args": {
                            "x": 15,
                            "y":  0,
                            "dx": 1,
                            "dy": 1
                        }
                    }
                },
                "out": {
                    "position": {
                        "name": "left",
                        "args": {
                            "x": 48,
                            "y": 45,
                            "dx": 1,
                            "dy": 1
                        }
                    }
                }
            }
        },
        "translator": {
            "image": "NONE",
            "port_groups": {
                "in": {
                    "position": {
                        "name": "left",
                        "args": {
                            "x": 15,
                            "y":  0,
                            "dx": 1,
                            "dy": 1
                        }
                    }
                },
                "out": {
                    "position": {
                        "name": "left",
                        "args": {
                            "x": 48,
                            "y": 45,
                            "dx": 1,
                            "dy": 1
                        }
                    }
                }
            }
        },
        "packed_column": {
            "image": "packed_column_1.svg",
            "port_groups": {
                "in": {
                    "position": {
                        "name": "left",
                        "args": {
                            "x": 48,
                            "y": 10,
                            "dx": 1,
                            "dy": 1
                        }
                    }
                },
                "out": {
                    "position": {
                        "name": "left",
                        "args": {
                            "x": 48,
                            "y": 40,
                            "dx": 1,
                            "dy": 1
                        }
                    }
                }
            }
        },
        "tray_column": {
            "image": "tray_column_1.svg",
            "port_groups": {
                "in": {
                    "position": {
                        "name": "left",
                        "args": {
                            "x": 48,
                            "y": 10,
                            "dx": 1,
                            "dy": 1
                        }
                    }
                },
                "out": {
                    "position": {
                        "name": "left",
                        "args": {
                            "x": 48,
                            "y": 40,
                            "dx": 1,
                            "dy": 1
                        }
                    }
                }
            }
        },
        "default": {
            "image": "default.svg",
            "port_groups": {
                "in": {
                    "position": {
                        "name": "left",
                        "args": {
                            "x":  2,
                            "y":  0,
                            "dx": 1,
                            "dy": 1
                        }
                    }
                },
                "out": {
                    "position": {
                        "name": "left",
                        "args": {
                            "x": 48,
                            "y": 50,
                            "dx": 1,
                            "dy": 1
                        }
                    }
                }
            }
        },
    }
