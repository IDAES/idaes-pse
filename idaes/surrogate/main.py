from idaes.surrogate import alamopy
from idaes.surrogate.my_surrogate_base import Surrogate
from idaes.surrogate.my_surrogate_base import Metrics

from sympy.parsing.sympy_parser import parse_expr
from sympy import sympify
from sympy import symbols, var
import sympy 

from pyomo.core.expr.sympy_tools import PyomoSympyBimap, sympy_available, Sympy2PyomoVisitor, sympy2pyomo_expression
import pyomo.environ as pyo
import numpy as np


class Alamopy(Surrogate):

	def __init__(self, settings):
		super().__init__(settings)
		self.alamopy_results = None

	def build_model(self):
		super().build_model()
		if self._vdata_in is not None and self._vdata_out is not None:
			self.alamopy_results = alamopy.alamo(self._rdata_in, self._rdata_out,  xval=self._vdata_in, zval=self._vdata_out, **self.settings)
		else:
			self.alamopy_results = alamopy.alamo(self._rdata_in, self._rdata_out, **self.settings)

		self._res = self.alamopy_results
		self._model = self.alamopy_results['model']

		self.handle_results(self._res)

	def handle_results(self, res):
		# sympy to pyomo converter
		# sympyp2pyomo_expression

		self._results[Metrics.RMSE] = self.alamopy_results['rmse']
		self._results[Metrics.SSE]  = self.alamopy_results['ssr']
		self._results[Metrics.Time] = self.alamopy_results['totaltime']
		self._results[Metrics.MSE]  = float(self.alamopy_results['rmse'])**2
		self._results[Metrics.Order] = self.alamopy_results['nbas']
		self._results[Metrics.R2]   = self.alamopy_results['R2']


		# Generate pyomo expression
		m = pyo.ConcreteModel()
		m.x = pyo.Var(range(len(self._rdata_in)))
		
		obj_map = PyomoSympyBimap()
		obj_map.sympy2pyomo = {}
		sympy_locals = {}
		i = 1
		for label in res['xlabels']:
			sympy_locals[label] = sympy.Symbol(label)
			sympy_obj = sympy.Symbol(label)
			obj_map.sympy2pyomo[sympy_obj] = m.x[i]
			i+=1

		model_string = res['model'].split('=')[1]
		model_symp = parse_expr(model_string.replace("^", "**"), local_dict = sympy_locals)
		model_pyomo = sympy2pyomo_expression(model_symp, obj_map)
		self._model = model_pyomo


