"""
Common Surrogate interface for IDAES.
"""
from pathlib import Path
from typing import Dict
import yaml


class Metrics:
    """Names for known types of metrics.

    Use these as keys in dictionaries, e.g.:

        m = {Metrics.RMSE: self.rmse}

    When adding attributes to this class, please include a comment with the
    prefix "#:" immediately above it, so Sphinx knows that this is documentation
    for the attribute. 
    """

    #: Root mean-squared error
    RMSE = "RMSE"


class ConfigurationError(Exception):
    pass


class Config:
    """Configuration for surrogate modeling tool(s).

    Can load/save itself from a file, and retrieve the settings for a given
    tool, or the global settings, as a dict.

    Typical usage (read-only):

        alamo_settings = Config("settings.json").get_tool("alamo")

    The configuration file format (JSON or YAML) is very simple:

        { "shared": { ...shared settings... },
          "<toolA-name>": { ...toolA-specific settings... },
          "<toolB-name>": { ...toolB-specific settings... },
          ...
        } 
    """

    #: Key for global settings section in configuration
    GLOBAL_SECTION = "shared"

    def __init__(self, input=None, values: Dict = None):
        """Constructor.
        
        Args:
            input: Filename, path, or file-like object (has `.read` attr), to load from file.
            values: If present, initial configuration values. If a file is also
                    provided, duplicate keys from the file will overwrite these values.

        Raises:
            ConfigurationError: if there is a problem reading from the file, or parsing the
                                configuration
        """
        # start with initial values
        if values is None:
            self._cfg = {}
        else:
            self._cfg = values
        # add file values
        if input is not None:
            f = self._get_fileobj(input)
            f_values = yaml.load(f)
            self._cfg.update(f_values)
            # remember file's path
            self._input_path = Path(f.name)
        else:
            self._input_path = None
        # ensure global section exists
        if not self.GLOBAL_SECTION in self._cfg:
            self._cfg[self.GLOBAL_SECTION] = {}

    def get_global(self) -> Dict:
        return self._cfg[self.GLOBAL_SECTION].copy()

    def set_global(self, values: Dict):
        self._cfg[self.GLOBAL_SECTION] = values

    def get_tool(self, tool: str) -> Dict:
        """Get configuration settings for a given tool.

        Args:
            tool: Name of tool

        Returns:
            Configuration value dictionary.

        Raises:
            KeyError: If no section for given tool name is found.
        """
        return self._cfg[tool].copy()  # raises KeyError if not found

    def set_tool(self, tool: str, values: Dict):
        """Change stored configuration settings to given values.

        Args:
            tool: Name of tool
            values: New settings

        Raises:
            KeyError: If no section for given tool name is found.
        """
        self._cfg[tool] = values  # raises KeyError if not found

    def save(self, output=None):
        """Save current values to output file.
        
        Args:
            output: Filename, path, or file-like object (has `.read` attr). If not provided,
                    will try to save to input file.
                    
        Raises:
            ConfigurationError: If no output, but an input file was not given or it is not writable.
        """
        if output is None:
            if self._input_path is None:
                raise ConfigurationError(
                    "No output given, and no input given when object was constructed"
                )
            try:
                ofile = self._input_path.open("w")
            except OSError as err:
                raise ConfigurationError(
                    f"Opening configuration file '{self._input_path}' for writing: {err}"
                )
        else:
            ofile = self._get_fileobj(output, mode="w")
        yaml.dump(self._cfg, ofile)

    @staticmethod
    def _get_fileobj(x, mode="r"):
        """Utility method to get a file object from a file (no-op), Path, or filename.
        """
        if hasattr(x, "read"):
            fileobj = x
        else:
            path = Path(x)
            try:
                fileobj = path.open(mode=mode)
            except OSError as err:
                modeing = "Reading from" if mode == "r" else "Writing to"
                raise ConfigurationError(
                    f"{modeing} configuration file '{path}': {err}"
                )
        return fileobj


class Surrogate:
    """Surrogate model interface definition.

    This is an abstract base class: do not instantiate it directly.
    """

    def __init__(self, model=None, settings=None):
        self.model = model
        self.settings = settings
        self._results = None
        self._vdata = None

    def build_model(self, **kwargs):
        pass

    @property
    def validation_data(self):
        """Get validation data.

        Returns:
             current validation data (may be None).
        """
        return self._vdata

    @validation_data.setter
    def validation_data(self, v):
        """Set validation data
        """
        self._vdata = v
        # TODO: invalidate any cached metrics
