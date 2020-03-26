"""
Common Surrogate interface for IDAES.
"""
from pathlib import Path
from typing import Dict
import yaml
from pyomo.common.config import ConfigBlock

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
                                                            'additional_features_list': list,
                                                            'pyomo_vars': list
                                                            }, total=False
                                    )

    PysmoRBF = TypedDict('PysmoRBF', {'basis_function': str,
                                      'solution_method': str,
                                      'regularization': bool,
                                      'pyomo_vars': list
                                      }, total=False
                         )

    PysmoKriging = TypedDict('PysmoKriging', {'numerical_gradients': bool,
                                              'regularization': bool,
                                              'pyomo_vars': list
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
                                                    additional_features=[]
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


# Single Surrogate Modeler
class Surrogate:

    CONFIG = ConfigBlock()

    # Common Declarations
    # CONFIG.declare()

    def __init__(self, **settings):

        # Config
        self.config = self.CONFIG(settings)
        self.modeler = None

        # Results
        self._results = {}
        self._metrics = None
        self._model = None
        self._b_built = False  # flag for regression

        # Data
        self._rdata_in = None
        self._rdata_out = None
        self._vdata_in = None
        self._vdata_out = None


    def modify_config(self, **kwargs):
        _b_built = False
        self.config = self.CONFIG(kwargs)

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
