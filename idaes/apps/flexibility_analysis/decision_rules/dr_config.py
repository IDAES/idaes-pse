from pyomo.common.config import ConfigDict


class DRConfig(ConfigDict):
    r"""
    A base class for specifying options for building
    decision rules.
    """

    def __init__(
        self,
        description=None,
        doc=None,
        implicit=False,
        implicit_domain=None,
        visibility=0,
    ):
        super().__init__(
            description=description,
            doc=doc,
            implicit=implicit,
            implicit_domain=implicit_domain,
            visibility=visibility,
        )
