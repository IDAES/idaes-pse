"""
Common Surrogate interface for IDAES.
"""
from pathlib import Path
from typing import Dict
import yaml

from mypy_extensions import TypedDict

class SurrogateModeler:
    """
    Default setup for known modeler setups
    """
    def surrogate_analysis():
        # external utitily
        # run all surrogate tools - paremeter()

        # ex 1: - pysmo and alamo

        # ex 2: range 1-3 monomial powers

        # analyze metrics
        # return best
        print("hello, I fit nothing")

    # MIGUEL
    # translation of similar commands across the two 
    # polynomial_order = 3 => both surrogate tools

    # General Settings that trigger differnt models

    # Type set-up for different tools
    PysmoPolyRegression = TypedDict('PysmoPolyRegression', {'number_of_crossvalidations': int,
                                                            'maximum_polynomial_order': int,
                                                            'training_split': float,
                                                            'solution_method': str,
                                                            'multinomials': bool,
                                                            'additional_features_list': list
                                                            }, total=False
                                    )

    PysmoRBF = TypedDict('PysmoRBF', {'basis_function': str,
                                      'solution_method': str,
                                      'regularization': bool
                                      }, total=False
                         )

    PysmoKriging = TypedDict('PysmoKriging', {'numerical_gradients': bool,
                                              'regularization': bool
                                              }, total=False
                             )

    # ALAMo argument types here

    AlamoDict = TypedDict('Alamo', {'xlabels': list,
                                 'zlabels': list,
                                 'xval': list,
                                 'zval': list,
                                 'xmin': list,
                                 'xmax': list,
                                 'modeler': int,
                                 'linfcns': int,
                                 'expfcns': int,
                                 'logfcns': int,
                                 'sinfcns': int,
                                 'cosfcns': int,
                                 'monomialpower': list,
                                 'multi2power': list,
                                 'multi3power': list,
                                 'ratiopower': list,
                                 'screener': int,
                                 'almname': str,
                                 'savescratch': bool,
                                 'savetrace': bool,
                                 'expandoutput': bool,
                                 'almopt': str,
                                 'loo': bool,
                                 'lmo': bool,
                                 'maxiter': int,
                                 'simulator': str
                                }, total=False
                         )

    # Base set-up for different tools: user cannot access
    pysmo_polyregression_base = PysmoPolyRegression(number_of_crossvalidations=2,
                                                    maximum_polynomial_order=2,
                                                    training_split=0.8,
                                                    solution_method='pyomo',
                                                    multinomials=False,
                                                    additional_features_list=[]
                                                    )

    pysmo_rbf_base = PysmoRBF(basis_function='gaussian',
                              solution_method='pyomo',
                              regularization=True
                              )

    pysmo_kriging_base = PysmoKriging(numerical_gradients=False,
                                      regularization=True
                                      )

    # Base settings for ALAMO modeller here

    # Issues multiple output from sets of inputs
    alamo_base = AlamoDict(monomialpower=(1, 2, 3, 4, 5, 6),
                           multi2power=(1, 2),
                           expandoutput=True
                )


class SurrogateFactory:
    """ 
    Construct default setting of a tool (CONFIG) 
    Changed by user to Surrogate instance Config


    Ex.
        alternative_cases = {'run_1': dict(pysmo_polyregression_base, maximum_polynomial_order=3),
                     'run_2': dict(pysmo_polyregression_base, solution_method='mle'),
                     'run_3': dict(pysmo_rbf_base, basis_function='cubic')
                     }

        tools = [
            ['pysmo_polygression', 'pysmo_rbf', 'pysmo_kriging', 'pysmo_polygression', 'pysmo_polygression', 'pysmo_rbf'],
            [pysmo_polyregression_base, pysmo_rbf_base, pysmo_kriging_base] + list(alternative_cases.values())
            ]  

        user loops:
            Surrogate = new instance()
    """


    def __init__(self, tool=None):
        """
            One tool 
        """
        self.tool = tool
        self.modeler = None

    def initialize_Modeler(self):
        """
        Intialize surrogate modeler
        Return Surrogate object
        """

        return self.modeler


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

    #: Mean-squared error
    MSE = "MSE"

    #: Sum of squared error
    SSE = "SSE"

    #: Time
    Time = "Time"

    #: Order
    Order = "Order"

    #: R-squared
    R2 = "R2"


class ConfigurationError(Exception):
    pass


# ConfigBlock
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


# Single Surrogate Modeler
class Surrogate:

    def __init__(self, settings=None):

        # Config
        self.settings = settings # Config ?
        self.modeler = None

        # Results
        self._results = None
        self._metrics = None
        self._model = None
        self._b_built = False # flag for regression

        # Data
        self._rdata_in = None
        self._rdata_out = None
        self._vdata_in = None
        self._vdata_out = None


    def modify_config(self, config: Config):
        _b_built = False
        self.settings = config

    # Build

    def build_model(self):
        self._b_built = True
        pass


    def update_model(self):
        self.build_model()
        pass

    def generate_expression(self, variable_list):
        pass

    # Get Results

    def get_model(self): # Pyomo Expression
        return self._model

    def get_metrics(self): # Metrics Object
        return self._metrics

    def get_results(self): # Pyomo Expression, Metrics Object
        return self._results

    # Data Handling

    def get_regressed_data(self):
    	return (self._rdata_in, self._rdata_out)

    def regressed_data(self, r_in, r_out): # 2D Numparray
        self._rdata_in = r_in
        self._rdata_out = r_out
     
    def get_validation_data(self):
        return (self._vdata_in, self._vdata.out)

    def validation_data(self, v_in, v_out): # 2D Numparray
        self._vdata_in = v_in
        self._vdata_out = v_out

    # Using regressed model
    def calculate_outputs(inputs): # 2D Numparray, use pyomo expression
        outputs = self._model(inputs)
        return

    # Additional Metrics - MUST NOT OVERWRITE MODELER METRICS

    def get_trained_metrics(): # 2D Nparray, 2D Numparray
        if not _b_built:
            print("Warning: No surrogate model regressed.")
            return
        pass

    def get_validated_metrics(xval, zval): # 2D, 2D Numparray, use pyomo expression
        if not _b_built:
            print("Warning: No surrogate model regressed.")
            return
        pass


