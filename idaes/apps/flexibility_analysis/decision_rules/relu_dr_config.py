from .dr_config import DRConfig
from pyomo.common.config import ConfigValue


class ReluDRConfig(DRConfig):
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
        self.n_layers: int = self.declare(
            "n_layers", ConfigValue(domain=int, default=4)
        )
        self.n_nodes_per_layer: int = self.declare(
            "n_nodes_per_layer", ConfigValue(domain=int, default=4)
        )
        self.tensorflow_seed: int = self.declare(
            "tensorflow_seed", ConfigValue(domain=int, default=0)
        )
        self.scale_inputs: bool = self.declare(
            "scale_inputs", ConfigValue(domain=bool, default=True)
        )
        self.scale_outputs: bool = self.declare(
            "scale_outputs", ConfigValue(domain=bool, default=True)
        )
        self.epochs: int = self.declare("epochs", ConfigValue(domain=int, default=2000))
        self.batch_size: int = self.declare(
            "batch_size", ConfigValue(domain=int, default=20)
        )
        self.learning_rate = self.declare(
            "learning_rate", ConfigValue(default=None)
        )
        self.plot_history: bool = self.declare(
            "plot_history", ConfigValue(domain=bool, default=False)
        )
