##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
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
	Property regression models

	Ability to use DMF csv files to regress important parameters 
	for literature based models and ALAMO models.

	Tracking property metadata and new model metadata
	* See :meth:'PropertyData.from_csv', for expected format for data.
	* See :meth:'PropertyMetadata()' for the expected format for metadata.

	Author: Marissa Engle
"""
##### WARNING WARNING WARNING WARNING WARNING WARNING
##### For now this module is going to get commented out.
##### It has not been tested in too long
##### WARNING WARNING WARNING WARNING WARNING WARNING

# # Local Imports
# from __future__ import division
# from pyomo.environ import *
# from pyomo.opt import SolverFactory
# # import idaes_dmf.propdb.types as prop_types # older version
#
# from idaes_dmf import dmf, resource, propdata
# from idaes_dmf.resource import Resource, FilePath;
# import pendulum
#
# import itertools
# import visPlot
# from idaes.surrogate import alamopy
#
# # Third-party imports
# try:
# 	import pandas as pd
# 	import numpy as np
# 	import sympy
# 	from sympy.parsing.sympy_parser import parse_expr
# except ImportError:
# 	np, pd = None, None
#
# # Change to your location of mydata in DMF
# my_dmf = dmf.DMF("/home/mengle/DMF/mydata/")
#
# class PropertyPackage(object):
# 	"Property package containing multiple Property Models"
#
# 	def __init__(self):
# 		self.textFile = None;
# 		self.build_package();
#
# 	def build_package(self):
# 		self.openFile();
# 		self.importData();
# 		self.buildModels();
# 		self.closeFile();
#
# 	# TODO: Write the property package resource
# 	# def create_dmf_PropertyPackageResource(self):
# 	# 	print 'TODO: Property Package'
#
# 	def openFile(self):
# 		print('Creating new text file');
# 		try:
# 			self.textFile = open('propTrial.py', 'w')
# 		except:
# 			print("couldn\'t open the file. Something wrong")
# 			sys.exit();
#
# 	def closeFile(self):
# 		""" Closes the GAMS file"""
# 		self.textFile.close();
# 		print('Closing the File');
#
# class PropertyModel(object):
#
# 	def __init__(self, data=None, metadata=None, outputFile = None, almFile = None, desc = None, **kwargs):
# 		self.textfile = outputFile;
# 		self.opt = SolverFactory('ipopt')
# 		self.model = ConcreteModel();
#
# 		self.alamo_model = ALAMOModel(data, metadata, almFile)
# 		self.alamo_model.set_name(almFile);
# 		self.table = propdata.PropertyTable(data=data, metadata=metadata)
# 		self.data = self.table.data
#
# 		#adds data resource
# 		# TODO: add if it doesn't exist
# 		self.create_dmf_DataResource(data, metadata,desc);
# 		# d.add(rsrc)
#
# 		names = data.names()
# 		# print list(names)
#
# 		# print self.alamo_model.xdata;
# 		# self.solve_model(self.model)
# 		if 'almopt' in kwargs:
# 			self.build_alamo_model(self.alamo_model, kwargs)
# 		else:
# 			self.build_alamo_model(self.alamo_model)
# 		# self.create_dmf_PropertyModelResoure();
#
# 		# result = self.choose(self.model, self.alamo_model)
#
# 	def create_dmf_DataResource(self, data=None, metadata= None, desc= None):
# 		self.table = propdata.PropertyTable(data=data, metadata=metadata)
# 		rsrc = resource.PropertyDataResource(self.table);
# 		spec = {'desc!':desc}
# 		if desc is not None:
# 			rsrc.desc = desc
#
# 		resources = my_dmf.find(spec)
# 		try:
# 			findResc = next(resources);
# 			print 'Existed Data Resource', desc
# 		except Exception, e:
# 			print "No Resource found end of generator issue"
# 		finally:
# 			print 'New Data Resource', desc
# 			my_dmf.add(rsrc)
# 		# my_dmf.remove(filter_dict={'desc':'test'})
# 		self.dataResource = rsrc
#
#
# 	def create_dmf_PropertyModelResoure(self):
# 		print 'Property Model Resource'
#
# 		rsrc = resource.Resource(type='property_model')
# 		rsrc.created = pendulum.now()
# 		rsrc.desc = 'property model for something'
# 		# rsrc.sources = [self.dataResource.sources];
# 		# rsrc.codes  = [resource.Code(type='python function', desc='ALAMO generated property model', name='SteamEnthalpyalm.py'),
# 				# resource.Code(type='python function', desc ='ALAMO generated parameter property model', name='SteamEnthalpycv.py')]
# 		# rsrc.datafiles = [self.dataResource.file];
# 		my_dmf.add(rsrc)
#
# 	def get_dataframe(self):
# 		df = self.data.values_dataframe();
# 		return df;
#
# 	def get_residuals(self):
# 		residuals = None;
# 		return residuals;
#
# 	def plot_residuals(self):
# 		# def plot_residual(flowsheet, idaes_model_object_name, x_name,
#                   # y_name, plot_title='Some default title',
#                   # save_directory='', file_format='png'):
# 		names = data.names();
#
# 		residuals = self.model.get_residuals();
# 		ALAMO_residuals = self.alamo_model.get_residuals();
#
# 		visPlot.plot_residual(None, self, 'X_label', 'residuals', "Residuals", "ResidualPlot");
#
# 	def build_model(self, model):
# 		self.data_sets()
# 		self.sets(model)
# 		self.parameters(model)
# 		self.variables(model)
# 		self.objective(model)
#
# 	def solve_model(self, model):
# 		self.build_model(model)
# 		self.results = self.opt.solve(model);
# 		# self.model_output(model, self.textfile);
#
# 		SSE_params = model.OBJ();
# 		self.results_str(model);
#
# 	def build_alamo_model(self, alamo_model):
# 		self.alamo_datasets(alamo_model);
# 		# self.alamo_datasets(alamo_model);
#
# 	def choose(self, model, alamo_model):
# 		# output for DMF and output of python file, need to change 1/10/2018
# 		AModels = alamo_model.models;
# 		results = min(AModels, key=lambda p: float(p['ssr']))
#
# 		SSE_params = model.OBJ();
# 		self.results_str(model);
#
# 		if SSE_params < float(results['ssr']):
# 			print "Write param model"
# 			self.model_output(model, self.textfile);
# 		else:
# 			self.alamo_output(results, self.textfile)
#
# class ALAMOModel(object):
#
# 	def __init__ (self,data, metadata, almname):
# 		self.xdata = []
# 		self.zdata = []
# 		self.models = [];
# 		self.almname = 'Property.alm'
#
#
# 	def get_residuals():
# 		residuals = None;
# 		return residuals;
#
# 	def set_name(self,almname):
# 		self.almname=almname;
#
# 	def build_alamo_model(self, xdata, zdata, xlabels = None, zlabels='e', **kwargs):
# 		# results=alamopy.doalamo(xdata,zdata,xlabels=xlbls,  zlabels='expr',  monomialpower='1,2,3', multi3power='1,2,3', multi2power='1,2,3', ratiopower='1 2', logfcns=1, expfcns=1, savescratch=1, expandoutput=1); #logfcns=1, expfcns=1, ratiopower='1 2 ', multi3power='1 2', multi2power='1 2 '
# 		# results=alamopy.doalamo(xdata,zdata,xlabels=xlbls,  zlabels='e',  monomialpower='1 2 3', multi3power='1 2 3', multi2power='1 2 3', savescratch=1, expandoutput=1); #logfcns=1, expfcns=1, ratiopower='1 2 ', multi3power='1 2', multi2power='1 2 '
# 		# results=alamopy.doalamo(xdata,zdata, almname='SteamEnthalpy.alm', almopt='steamcustom.alm', monomialpower='1 2 3', multi3power='1 2 3', multi2power='1 2 3', savescratch=1, expandoutput=1); #logfcns=1, expfcns=1, ratiopower='1 2 ', multi3power='1 2', multi2power='1 2 '
#
# 		if 'almopt' in kwargs:
# 			results=alamopy.doalamo(xdata,zdata, xlabels=xlabels,  zlabels='e', almname=self.almname, almopt=kwargs.get("almopt"), monomialpower='1 2 3', multi3power='1 2 3', multi2power='1 2 3', savescratch=1, expandoutput=1); #logfcns=1, expfcns=1, ratiopower='1 2 ', multi3power='1 2', multi2power='1 2 '
# 		else:
# 			results=alamopy.doalamo(xdata,zdata, almname=self.almname, monomialpower='1 2 3', multi3power='1 2 3', multi2power='1 2 3', savescratch=1, expandoutput=1); #logfcns=1, expfcns=1, ratiopower='1 2 ', multi3power='1 2', multi2power='1 2 '
#
# 		print "Model function: %s" % results['model'];
# 		print "Sum of squared error = %s" % results['ssr']
# 		print "R^2 = %s" % results['R2'];
# 		res = alamopy.almconfidence(results, xdata, zdata);
# 		print "Confidence Interval =%s" % results['conf_inv']
# 		alamopy.almplot(res)
# 		print results['model'].replace("^","**")
# 		alm_model = results['model']
# 		self.models.append(results);
#
# 	def enumerate_datasets(self, xvars, xlabels, xdata, **kwargs):
# 		xsets= range(0,xvars)
#
# 		xlabels = np.array(xlabels)
#
# 		# for i in range(1, xvars+1):
# 		# 	for subset in itertools.combinations(xsets, i):
# 		# 		sub = np.array(subset)
# 		# 		#print self.xdata[:,sub]
# 		# 		self.build_alamo_model(self.xdata[:,sub], self.zdata, tuple(xlabels[sub]),'expr')
# 		self.build_alamo_model(self.xdata[:,:], self.zdata, xlabels= xlabels, zlabels='e', **kwargs)
# 		#xdata.append(self.data.get_column('r').values)
#
# 	# def model_output(self,  model, textfile):
# 	# 	print "Write ALAMO Model"
# 	# 	textfile.write("def _V_liq(self): \n");
# 	# 	textfile.write("\t %s)\n\n" %model['model'])
#
#
# # Example with SpiraxSteamTableEnthalpy.csv
# class SteamPropertyPackage(PropertyPackage):
# 	def __init__(self):
# 		self.textFile = None;
# 		self.build_package();
#
# 	def importData(self):
# 		# self.data1 = propdata.PropertyData.from_csv("/home/mengle/DMF/mydata/SpiraxSteamTableEnthalpy.csv",2);
# 		# self.metadata1 = propdata.PropertyMetadata.from_csv("/home/mengle/DMF/mydata/SpiraxSteamTableEnthalpy-meta.csv")
#
# 		self.data1 = propdata.PropertyData.from_csv("/home/mengle/DMF/mydata/SpiraxEnthalpy.csv",2);
# 		self.metadata1 = propdata.PropertyMetadata.from_csv("/home/mengle/DMF/mydata/SpiraxEnthalpy-meta.csv")
#
# 	def buildModels(self):
# 		print "Enthalpy Model"
# 		Enthalpy = EnthalpyModel(self.data1, self.metadata1, self.textFile, "Enthalpy", "spirax enthalpy");
#
# class EnthalpyModel(PropertyModel):
#
#     def __init__(self, data, metadata, outputFile, almtrace, desc):
#     	PropertyModel.__init__(self,data, metadata, outputFile, almtrace, desc);
#
#     def results_str(self, model):
#     	print "Sum of squared error = %f" % model.OBJ();
#
#     def model_output(self, model, textfile):
#     	textfile.write("def _h_vap(self): \n")
#     	textfile.write("\tTref=PhysPropParams.Tref")
#     	textfile.write("\tPref=PhysPropParams.Pref")
#
#     def data_sets(self):
#     	self.I = np.asarray(range(self.data.num_rows))
#     	self.Tr = dict(zip(self.I,self.data.get_column('Tr').values))
#     	self.Pr = dict(zip(self.I,self.data.get_column('Pr').values))
#     	# self.h_vap = dict(zip(self.I,self.data.get_column('Jkg').values))
#     	self.h_vap = dict(zip(self.I,self.data.get_column('Jmol').values))
#
#
#     def sets(self,model):
#     	model.I = Set(initialize=self.I)
#
#     def parameters(self,model):
#     	# Given parameters of the system
#     	model.Tr = Param(model.I, initialize = self.Tr)
#     	model.Pr = Param(model.I, initialize = self.Pr)
#     	model.h_vap = Param(model.I, intialize  = self.h_vap)
#
#     def alamo_datasets(self,alamo_model):
#     	# zdata = np.asarray(self.data.get_column('Jkg').values)
#     	zdata = np.asarray(self.data.get_column('Jmol').values)
#     	self.alamo_model.zdata = zdata.reshape(len(zdata),1);
#
#     	xdata = []
#     	xdata.append(self.data.get_column('Tr').values);
#     	xdata.append(self.data.get_column('Pr').values);
#
#     	self.alamo_model.xdata = np.transpose(xdata);
#     	self.alamo_model.enumerate_datasets(2, ['Tr', 'Pr'], self.alamo_model.xdata)#, almopt='steamcustom.alm')
#
#     def alamo_output(self,alamo_model, textfile):
#     	textfile.write("def h_vap(self): \n")
#     	textfile.write("\t# ALAMO Model")
#
#
#
# 	# def alamo_output(self, alamo_model,textfile):
# 	# 	textfile.write("def _h_vap(self): \n");
# 	# 	textfile.write("\t# ALAMO Model\n")
# 	# 	textfile.write("\t(T, P, y) = self.deref_tpy()\n");
# 	# 	textfile.write("\tself.var_or_expr(name=\"visc_liq\", doc=\"Liquid molar volume\", expr=%s)\n\n" %alamo_model['model'])
#
# 	#     # def _h_vap(self):
# 	#     #     c = CpIdealGasParams.c
# 	#     #     Tref = PhysPropParams.Tref
# 	#     #     hv_h2o = 43594
# 	#     #     hv_mea = 58000
# 	#     #     (T, P, y) = self.deref_tpy()
# 	#     #     expr_h_vap = sum((
# 	#     #         (c[i][0]*T + c[i][1]/2*T**2 + c[i][2]/3*T*3 + c[i][3]/4*T**4) -
# 	#     #         (c[i][0]*Tref + c[i][1]/2*Tref**2 + c[i][2]/3*Tref*3 + c[i][3]/4*Tref**4)
# 	#     #         )*y[i] for i in y.keys()) + y["H2O"]*hv_h2o + y["MEA"]*hv_mea
# 	#     #     self.var_or_expr(name="h_vap", expr=expr_h_vap, doc="Ideal gas enthalpy")
#
#
#
# class MEAPropertyPackage(PropertyPackage):
#
# 	def __init__(self):
# 		self.textFile = None;
# 		self.build_package();
#
# 	def importData(self):
# 		d = dmf.DMF("/home/mengle/DMF/mydata/");
# 		self.data1 = propdata.PropertyData.from_csv('/home/mengle/DMF/examples/jayarathna-data.csv', 2)
# 		self.metadata1 = propdata.PropertyMetadata.from_csv('/home/mengle/DMF/examples/jayarathna-meta.csv')#, 'MEA')
#
# 		self.data2 = propdata.PropertyData.from_csv("/home/mengle/DMF/mydata/amundsen-data.csv", 2)
# 		self.metadata2 = propdata.PropertyMetadata.from_csv("/home/mengle/DMF/mydata/amundsen-meta.csv")
#
# 		self.data3 = propdata.PropertyData.from_csv("/home/mengle/DMF/mydata/han-data.csv", 2)
# 		num_added = self.data3.add_csv("/home/mengle/DMF/mydata/jayarathna-data.csv")
# 		num_added = self.data3.add_csv("/home/mengle/DMF/mydata/amundsen-data2.csv")
# 		self.metadata3 = propdata.PropertyMetadata.from_csv("/home/mengle/DMF/mydata/han-meta.csv")
#
# 	def buildModels(self):
# 		print "Surface Tension Model"
# 		SurfaceTension = SurfaceTensionModel(self.data1, self.metadata1, self.textFile, "SurfaceTension", "jayrathna surface tension");
# 		# SurfaceTension.get_dataframe();
#
# 		print "Viscosity Model"
# 		Viscosity = ViscosityModel(self.data2,self.metadata2, self.textFile, "Viscosity", "amundsen viscosity")
#
# 		print "Density Model"
# 		Density = DensityModel(self.data3, self.metadata3, self.textFile, "Density", "amundsen jayrathna hans density")
#
# class DensityModel(PropertyModel):
#
# 	# class VLiqParams(object):
# 		#    # Liquid molar volume parameters (Weiland et al. 1998, adjusted for si units)
# 		#    a = {'MEA':-5.35162e-4, 'H2O':-3.2484e-3, 'CO2':0}
# 		#    b = {'MEA':-0.451417, 'H2O':1.65, 'CO2':0}
# 		#    c = {'MEA':1194.51, 'H2O':793.0, 'CO2':927.0}
# 		#    Vst = -1.8218e-6
#
# 		# def _V_liq(self):
# 	    #    (T, P, y) = self.deref_tpy()
# 	    #    a = VLiqParams.a
# 	    #    b = VLiqParams.b
# 	    #    c = VLiqParams.c
# 	    #    Vst = VLiqParams.Vst
# 	    #    mw = PhysPropParams.mw
# 	    #    if self.external:
# 	    #        self.var_or_expr(name="V_liq2",
# 	    #            expr=self.extf_V_liq(T, y['CO2'], y['H2O'], 0.0, y['MEA']),
# 	    #            doc="Liquid molar volume")
# 	    #        return
# 	    #    self.var_or_expr(name="V_liq", doc="Liquid molar volume",
# 	    #        expr=sum(y[i]*mw[i]/(a[i]*T**2 + b[i]*T + c[i]) for i in y.keys()) +
# 	    #        y['MEA']*y['H2O']*Vst)
#
# 	def __init__(self, data, metadata, outputFile, almtrace, description):
# 		PropertyModel.__init__(self, data, metadata, outputFile, almtrace, description);
#
#
# 	def results_str(self, model):
# 		print model.a_mea.value, model.b_mea.value, model.c_mea.value, model.c_CO2.value, model.Vstar.value;
# 		print "Sum of squared error = %f" % model.OBJ();
#
#
# 	def model_output(self, model, textfile):
# 		textfile.write("def _V_liq(self): \n");
# 		textfile.write("\t(T, P, y) = self.deref_tpy()\n");
# 		for v in model.component_objects(Var, active = True):
# 			textfile.write("\t%s = %s\n" %(v.name, v.value))
# 		textfile.write("\tVst = -1.8218e-6\n")
# 		textfile.write("\tmw = PhysPropParams.mw\n")
# 		textfile.write("\tr = self.wt_frac[\"MEA\"]\n\n");
# 		expr = "sum(y[i]*mw[i]/(a[i]*T**2 + b[i]*T + c[i]) for i in y.keys()) + y['MEA']*y['H2O']*Vst)"
# 		textfile.write("\tself.var_or_expr(name=\"V_liq\", doc=\"Liquid molar volume\", expr=%s)\n\n" %expr)
#
# 	def data_sets(self):
# 		self.I = np.asarray(range(self.data.num_rows));
# 		self.Temp = dict(zip(self.I,self.data.get_column('T').values));
# 		self.alpha = dict(zip(self.I,self.data.get_column('CO2 Loading').values));
# 		self.rstar = dict(zip(self.I,self.data.get_column('r').values));
# 		self.rho = dict(zip(self.I,self.data.get_column('Density Data').values));
#
# 	def sets(self, model):
# 		model.I = Set(initialize=self.I)
#
# 	def parameters(self, model):
# 		# Given parameters of the system
# 		model.alpha = Param(model.I, initialize = self.alpha)
# 		model.r = Param(model.I, initialize=self.rstar)
# 		model.T = Param(model.I, initialize=self.Temp)
#
# 		# Measured values of interest
# 		model.rho_sln = Param(model.I, initialize=self.rho)
#
# 	def variables(self, model):
# 		# Regressed parameter values
# 		model.a_mea = Var(domain=Reals, initialize = 1)
# 		model.b_mea = Var(domain= Reals, initialize = 1)
# 		model.c_mea = Var(domain=Reals, initialize = 1)
#
# 		#model.a_CO2 = Var(domain=Reals, initialize = 0)
# 		#model.b_CO2 = Var(domain= Reals, initialize = 0)
# 		model.c_CO2 = Var(domain=Reals, initialize = 1)
#
# 		#model.d = Var(domain=Reals, initialize = 1)
# 		#model.e = Var(domain=Reals, initialize = 1)
# 		#model.V_CO2 = Var(domain=NonNegativeReals, initialize = 1)
# 		model.Vstar = Var(domain=Reals, initialize = -1)
#
# 	def objective(self, model):
#
# 		def obj_expression(model):
# 			expr = 0
#
# 			# Molecular Weight
# 			MW_mea = 61.08
# 			MW_H2O = 18.02
# 			MW_CO2 = 44.01
#
# 			for i in model.I:
#
# 				# Mol Fractions
# 				X_mea = 1/(1+ model.alpha[i] + (MW_mea/MW_H2O)*((1-model.r[i])/model.r[i]))
# 				X_CO2 = X_mea * model.alpha[i]
# 				X_H2O =  1 - X_mea - X_CO2
# 				MW_sln = X_mea*MW_mea + X_CO2*MW_CO2 + X_H2O*MW_H2O
#
# 				# Volume
# 				V_mea = MW_mea/((model.a_mea/1e7)*model.T[i]**2 + (model.b_mea/1e4)*model.T[i] + model.c_mea)
# 				V_CO2 = MW_CO2/model.c_CO2;
# 				#V_CO2 = MW_CO2/((model.a_CO2/1e7)*model.T[i]**2 + (model.b_CO2/1e4)*model.T[i] + model.c_CO2)
# 				V_H2O = MW_H2O/(-(32.484/1e7)*model.T[i]**2+ (16.5/1e4)*model.T[i] + 0.793)
# 				#Vstar2 = model.d + model.e*X_mea
#
# 				#V_sln = X_mea*V_mea + X_H2O*V_H2O + X_CO2*model.V_CO2 + X_mea*X_H2O*model.Vstar # + X_mea*X_CO2*Vstar2
# 				V_sln = X_mea*V_mea + X_H2O*V_H2O + X_CO2*V_CO2 + X_mea*X_H2O*model.Vstar # + X_mea*X_CO2*Vstar2
#
# 				expr+= (model.rho_sln[i] - MW_sln/V_sln)**2
# 			return expr;
#
# 		model.OBJ = Objective(rule=obj_expression)
#
# 	def alamo_datasets(self,alamo_model):
# 		# Measured value of interest
# 		zdata = np.asarray(self.data.get_column('Density Data').values);
# 		self.alamo_model.zdata = zdata.reshape((len(zdata),1))
#
# 		# Molecular Weight
# 		MW_mea = 61.08
# 		MW_H2O = 18.02
# 		MW_CO2 = 44.01
#
# 		xdata = []
# 		xdata.append(self.data.get_column('T').values)
#
# 		alpha = self.data.get_column('CO2 Loading').values;
# 		rstar = self.data.get_column('r').values;
#
# 		X_mea = [1/(1+ a + (MW_mea/MW_H2O)*((1-r)/r)) for a,r in zip(alpha, rstar)];
# 		X_CO2 = [X_m* a for X_m,a in zip(X_mea,alpha)]
# 		X_H2O =  [1 - X_m - X_C for X_m, X_C in zip(X_mea, X_CO2)];
#
# 		xdata.append(X_mea);
# 		xdata.append(X_CO2);
# 		xdata.append(X_H2O);
#
# 		self.alamo_model.xdata = np.transpose(xdata)
# 		# self.alamo_model.enumerate_datasets(4, ['T','y[\'MEA\']','y[\'CO2\']', 'y[\'H2O\']'], self.alamo_model.xdata)
# 		self.alamo_model.enumerate_datasets(4, ['T','yMEA','yCO2', 'yH2O'], self.alamo_model.xdata)
# 		# self.alamo_model.enumerate_datasets(4,'T yMEA yCO2 yH2O', self.alamo_model.xdata)
#
# 	def alamo_output(self, alamo_model,textfile):
# 		textfile.write("def _V_liq(self): \n");
# 		textfile.write("\t# ALAMO Model\n")
# 		textfile.write("\t(T, P, y) = self.deref_tpy()\n");
# 		textfile.write("\tVst = -1.8218e-6\n")
# 		textfile.write("\tmw = PhysPropParams.mw\n")
# 		textfile.write("\tr = self.wt_frac[\"MEA\"]\n\n");
# 		expr = "sum(y[i]*mw[i]/(a[i]*T**2 + b[i]*T + c[i]) for i in y.keys()) + y['MEA']*y['H2O']*Vst)"
# 		textfile.write("\tself.var_or_expr(name=\"visc_liq\", doc=\"Liquid molar volume\", expr=%s)\n\n" %alamo_model['model'])
#
# class SurfaceTensionModel(PropertyModel):
#
# 	# class SurfaceTensionParams(object):
# 		#     ### Surface tension parameters
# 		#     ###
# 		#     # Pure component surface tension for MEA and H2O
# 		#     # (Asprion 2005)
# 		#     C = {'MEA':(0.09945, 1.067,  0,      0    ),
# 		#          'H2O':(0.18548, 2.717, -3.554, 2.047)}
# 		#     # Parameters for CO2 surface tension
# 		#     S = {'CO2':(-5.987, 3.7699, -0.43164, 0.018155, -0.01207, 0.002119)}
# 		#     F = {
# 		#         'a':  2.4558,
# 		#         'b': -1.5311,
# 		#         'c':  3.4994,
# 		#         'd': -5.6398,
# 		#         'e': 10.2109,
# 		#         'f':  2.3122,
# 		#         'g':  4.5608,
# 		#         'h': -2.3924,
# 		#         'i':  5.3324,
# 		#         'j':-12.0494}
#
#
# 	    # def _surf_tens(self):
# 	     #    (T, P, y) = self.deref_tpy()
# 	     #    if self.external:
# 	     #        self.var_or_expr(name="surf_tens",
# 	     #            expr=self.extf_sig(T, y['CO2'], y['H2O'], y['MEA']))
# 	     #        return
# 	     #    F = SurfaceTensionParams.F
# 	     #    Tc = PhysPropParams.Tc
# 	     #    c = SurfaceTensionParams.C
# 	     #    S = SurfaceTensionParams.S
# 	     #    # Variables
# 	     #    r = self.wt_frac['CO2']
# 	     #    alpha = self.expr_loading = Expression(expr=y["CO2"]/y["MEA"])
#
# 	     #    # Constraints
# 	     #    def rule_sigma_i(blk, i):
# 	     #        if i in ['MEA', 'H2O']: #correlation type 1
# 	     #            return c[i][0]*(1-T/Tc[i])**\
# 	     #                (c[i][1] + c[i][2]*T/Tc[i] + c[i][3]*(T/Tc[i])**2)
# 	     #        elif i=='CO2': #correlation type 2
# 	     #            return S['CO2'][0]*r**2 + S['CO2'][1]*r + S['CO2'][2] +\
# 	     #                T*(S['CO2'][3]*r**2 + S['CO2'][4]*r+S['CO2'][5])
# 	     #    self.var_or_expr(["MEA", "CO2", "H2O"], name="sigma_i",
# 	     #        rule=rule_sigma_i)
# 	     #    sigma_i = self.sigma_i
# 	     #    self.var_or_expr(name="surf_tens", expr=sigma_i['H2O'] +
# 	     #        (sigma_i['H2O'] - sigma_i['CO2'])*y['CO2']*
# 	     #        (F['a']+F['b']*alpha+F['c']*alpha**2+F['d']*r+F['e']*r**2) +
# 	     #        (sigma_i['MEA'] - sigma_i['H2O'])*y['MEA']*
# 	     #        (F['f']+F['g']*alpha+F['h']*alpha**2+F['i']*r+F['j']*r**2))
#
#
# 	def __init__(self, data, metadata, outputFile, almtrace, description):
# 		PropertyModel.__init__(self, data, metadata, outputFile, almtrace, description);
#
#
# 	def results_str(self, model):
# 		# print self.results;
# 		print model.a.value, model.b.value, model.c.value, model.d.value, model.e.value, model.f.value, model.g.value, model.h.value, model.i.value, model.j.value
# 		print "Sum of squared error = %f" % model.OBJ()
#
#
# 	def model_output(self, model, textfile):
# 		textfile.write("def _surf_tens(self): \n");
# 		textfile.write("\t# Pure component surface tension for MEA and H2O\n")
# 		textfile.write("\t# (Asprion 2005)\n");
# 		textfile.write("\t(T, P, y) = self.deref_tpy()\n")
# 		textfile.write("\tC = {'MEA':(0.09945, 1.067,  0,      0    ),\n")
# 		textfile.write("\t'H2O':(0.18548, 2.717, -3.554, 2.047)}\n\n");
#
# 		textfile.write("\t# Parameters for CO2 surface tension\n")
# 		textfile.write("\tS = {'CO2':(-5.987, 3.7699, -0.43164, 0.018155, -0.01207, 0.002119)}\n");
# 		textfile.write("\tF = {")
# 		for v in model.component_objects(Var, active = True):
# 			if v.name== 'j':
# 				textfile.write("\n\t'%s' = %s}\n" %(v.name, v.value))
# 			else:
# 				textfile.write("\n\t'%s' = %s," %(v.name, v.value))
# 		textfile.write("\tTc = PhysPropParams.Tc\n")
#
# 		textfile.write("\t# Variables\n")
# 		textfile.write("\tr = self.wt_frac['CO2']\n");
# 		textfile.write("\talpha = self.expr_loading = Expression(expr=y[\"CO2\"]/y[\"MEA\"])\n");
#
#
# 		textfile.write("\t# Constraints\n");
# 		textfile.write("\tdef rule_sigma_i(blk, i):\n");
# 		textfile.write("\t\tif i in ['MEA', 'H2O']: #correlation type 1\n");
# 		textfile.write("\t\t\treturn c[i][0]*(1-T/Tc[i])** (c[i][1] + c[i][2]*T/Tc[i] + c[i][3]*(T/Tc[i])**2)\n");
# 		textfile.write("\t\telif i=='CO2': #correlation type 2\n");
# 		textfile.write("\t\t\treturn S['CO2'][0]*r**2 + S['CO2'][1]*r + S['CO2'][2] + T*(S['CO2'][3]*r**2 + S['CO2'][4]*r+S['CO2'][5])\n");
# 		textfile.write("\tself.var_or_expr([\"MEA\", \"CO2\", \"H2O\"], name=\"sigma_i\", rule=rule_sigma_i)\n");
# 		textfile.write("\tsigma_i = self.sigma_i\n");
# 		textfile.write("\tself.var_or_expr(name=\"surf_tens\", expr=sigma_i['H2O'] +(sigma_i['H2O'] - sigma_i['CO2'])*y['CO2']*(F['a']+F['b']*alpha+F['c']*alpha**2+F['d']*r+F['e']*r**2) +(sigma_i['MEA'] - sigma_i['H2O'])*y['MEA']*(F['f']+F['g']*alpha+F['h']*alpha**2+F['i']*r+F['j']*r**2))\n\n");
#
# 	def data_sets(self):
# 		self.I = np.asarray(range(self.data.num_rows));
#
# 		self.Temp = dict(zip(self.I,self.data.get_column('T').values));
# 		self.alpha = dict(zip(self.I,self.data.get_column('CO2 Loading').values));
# 		self.rstar = dict(zip(self.I,self.data.get_column('r').values));
#
# 		# Measured value of interest
# 		self.yVal = dict(zip(self.I,self.data.get_column('Surface Tension').values));
#
# 	def sets(self, model):
# 		model.I = Set(initialize=self.I)
#
# 	def parameters(self, model):
# 		# Given parameters of the system
# 		model.alpha = Param(model.I, initialize=self.alpha)
# 		model.rstar = Param(model.I, initialize=self.rstar)
# 		model.T = Param(model.I, initialize=self.Temp)
#
# 		# Measured value of interest
# 		model.y = Param(model.I, initialize=self.yVal)
#
# 	def variables(self, model):
# 		model.a = Var(domain=Reals, initialize = 0)
# 		model.b = Var(domain= Reals, initialize = 0)
# 		model.c = Var(domain=Reals, initialize = 0)
# 		model.d = Var(domain=Reals, initialize = 0)
# 		model.e = Var(domain=Reals, initialize = 0)
# 		model.f = Var(domain=Reals, initialize = 0)
# 		model.g = Var(domain=Reals, initialize = 0)
# 		model.h = Var(domain=Reals, initialize = 0)
# 		model.i = Var(domain=Reals, initialize = 0)
# 		model.j = Var(domain=Reals, initialize = 0)
#
# 	def objective(self, model):
#
# 		def obj_expression(model):
# 			expr = 0
# 			MW_mea = 61.08
# 			MW_h2o = 18.02
# 			MW_co2 = 44.01
#
# 			Tc_h2o = 647.13
# 			Tc_mea = 614.45
#
# 			c1_h2o = 0.18548
# 			c2_h2o = 2.717
# 			c3_h2o = -3.554
# 			c4_h2o = 2.047
# 			c1_mea = 0.09945
# 			c2_mea = 1.067
# 			c3_mea = 0
# 			c4_mea = 0
#
# 			S1 = -5.987
# 			S2 = 3.7699
# 			S3 = -0.43164
# 			S4 = 0.018155
# 			S5 = -0.01207
# 			S6 = 0.002119
#
# 			for i in model.I:
# 				X_mea = 1/(1+model.alpha[i]+(MW_mea/MW_h2o)*(1-model.rstar[i])/model.rstar[i])
# 				X_co2 = model.alpha[i]*X_mea
#
# 				sigma_h2o = c1_h2o*(1-model.T[i]/Tc_h2o)**(c2_h2o+c3_h2o*(model.T[i]/Tc_h2o)+c4_h2o*(model.T[i]/Tc_h2o)**2)
# 				sigma_mea = c1_mea*(1-model.T[i]/Tc_mea)**(c2_mea+c3_mea*(model.T[i]/Tc_mea)+c4_mea*(model.T[i]/Tc_mea)**2)
# 				sigma_co2 = S1*model.rstar[i]**2+S2*model.rstar[i]+S3+model.T[i]*(S4*model.rstar[i]**2+S5*model.rstar[i]+S6)
#
# 				f_func = model.a+model.b*model.alpha[i]+model.c*model.alpha[i]**2+model.d*model.rstar[i]+model.e*model.rstar[i]**2
# 				g_func = model.f+model.g*model.alpha[i]+model.h*model.alpha[i]**2+model.i*model.rstar[i]+model.j*model.rstar[i]**2
#
# 				sigma_calc=sigma_h2o+(sigma_co2-sigma_h2o)*f_func*X_co2+(sigma_mea-sigma_h2o)*g_func*X_mea
#
# 				expr+= (model.y[i] - sigma_calc)**2
#
# 			return expr
#
# 		model.OBJ = Objective(rule=obj_expression)
#
# 	def alamo_datasets(self,alamo_model):
# 		# Measured value of interest
# 		zdata = np.asarray(self.data.get_column('Surface Tension').values);
# 		#print zdata
# 		self.alamo_model.zdata = zdata.reshape((len(zdata),1))
#
# 		MW_mea = 61.08
# 		MW_h2o = 18.02
# 		MW_co2 = 44.01
#
# 		Tc_h2o = 647.13
# 		Tc_mea = 614.45
#
# 		c1_h2o = 0.18548
# 		c2_h2o = 2.717
# 		c3_h2o = -3.554
# 		c4_h2o = 2.047
# 		c1_mea = 0.09945
# 		c2_mea = 1.067
# 		c3_mea = 0
# 		c4_mea = 0
#
# 		S1 = -5.987
# 		S2 = 3.7699
# 		S3 = -0.43164
# 		S4 = 0.018155
# 		S5 = -0.01207
# 		S6 = 0.002119
#
# 		xdata = []
# 		T = self.data.get_column('T').values;
# 		alpha = self.data.get_column('CO2 Loading').values;
# 		rstar = self.data.get_column('r').values;
#
# 		X_mea = [1/(1+a+(MW_mea/MW_h2o)*(1-r)/r) for a,r in zip(alpha, rstar)]
# 		X_co2 = [a*x_m for a, x_m in zip(alpha, X_mea)]
#
# 		xdata.append(self.data.get_column('T').values)
# 		xdata.append(X_mea)
# 		xdata.append(X_co2)
#
# 		sigma_h2o = [c1_h2o*(1-t/Tc_h2o)**(c2_h2o+c3_h2o*(t/Tc_h2o)+c4_h2o*(t/Tc_h2o)**2) for t in T]
# 		sigma_mea = [c1_mea*(1-t/Tc_mea)**(c2_mea+c3_mea*(t/Tc_mea)+c4_mea*(t/Tc_mea)**2) for t in T]
# 		sigma_co2 = [S1*r**2+S2*r+S3+t*(S4*r**2+S5*r+S6) for r,t in zip(rstar, T)]
#
# 		xdata.append(sigma_h2o)
# 		xdata.append(sigma_mea)
# 		xdata.append(sigma_co2)
#
# 		xdata.append(rstar)
#
# 		self.alamo_model.xdata = np.transpose(xdata)
# 		#self.alamo_model.enumerate_datasets(4, ['T','y[\'MEA\']','y[\'CO2\']', 'sigma_i[\'H2O\']', 'sigma_i[\'MEA\']', 'sigma_i[\'CO2\']', 'r'], self.alamo_model.xdata)
# 		self.alamo_model.enumerate_datasets(4, ['T','yM','yC', 'sigma_iH2O', 'sigma_iMEA', 'sigma_iCO2', 'r'], self.alamo_model.xdata)
#
# 	def alamo_output(self, alamo_model, textfile):
# 		textfile.write("def _surf_tens(self): \n");
# 		textfile.write("\t# ALAMO Model\n")
# 		textfile.write("\t# Pure component surface tension for MEA and H2O\n")
# 		textfile.write("\t# (Asprion 2005)\n");
# 		textfile.write("\t(T, P, y) = self.deref_tpy()\n")
# 		textfile.write("\tC = {'MEA':(0.09945, 1.067,  0,      0    ),\n")
# 		textfile.write("\t'H2O':(0.18548, 2.717, -3.554, 2.047)}\n\n");
#
# 		textfile.write("\t# Parameters for CO2 surface tension\n")
# 		textfile.write("\tS = {'CO2':(-5.987, 3.7699, -0.43164, 0.018155, -0.01207, 0.002119)}\n");
# 		textfile.write("\tF = {")
# 		# for v in model.component_objects(['T','y[\'MEA\']','y[\'CO2\']', 'sigma_i[\'H2O\']', 'sigma_i[\'MEA\']', 'sigma_i[\'CO2\']'], active = True):
# 		# 	if v.name== 'j':
# 		# 		textfile.write("\n\t'%s' = %s}\n" %(v.name, v.value))
# 		# 	else:
# 		# 		textfile.write("\n\t'%s' = %s," %(v.name, v.value))
# 		textfile.write("}\n")
# 		textfile.write("\tTc = PhysPropParams.Tc\n")
#
# 		textfile.write("\t# Variables\n")
# 		textfile.write("\tr = self.wt_frac['CO2']\n");
# 		textfile.write("\talpha = self.expr_loading = Expression(expr=y[\"CO2\"]/y[\"MEA\"])\n");
#
#
# 		textfile.write("\t# Constraints\n");
# 		textfile.write("\tdef rule_sigma_i(blk, i):\n");
# 		textfile.write("\t\tif i in ['MEA', 'H2O']: #correlation type 1\n");
# 		textfile.write("\t\t\treturn c[i][0]*(1-T/Tc[i])** (c[i][1] + c[i][2]*T/Tc[i] + c[i][3]*(T/Tc[i])**2)\n");
# 		textfile.write("\t\telif i=='CO2': #correlation type 2\n");
# 		textfile.write("\t\t\treturn S['CO2'][0]*r**2 + S['CO2'][1]*r + S['CO2'][2] + T*(S['CO2'][3]*r**2 + S['CO2'][4]*r+S['CO2'][5])\n");
# 		textfile.write("\tself.var_or_expr([\"MEA\", \"CO2\", \"H2O\"], name=\"sigma_i\", rule=rule_sigma_i)\n");
# 		textfile.write("\tsigma_i = self.sigma_i\n");
# 		textfile.write("\tself.var_or_expr(name=\"surf_tens\", expr=%s)\n\n" %alamo_model['model'].replace("^","**"));
#
# class ViscosityModel(PropertyModel):
#
# 	# class ViscLiqParams(object):
# 		#     """
# 		#     Viscosity parameters for an MEA solvent
# 		#     """
# 		#     a = -0.0838
# 		#     b = 2.8817
# 		#     c = 33.651
# 		#     d = 1817.0
# 		#     e = 0.00847
# 		#     f = 0.0103
# 		#     g = -2.3890
#
#
# 	    # def _visc_liq(self):
# 	    #     (T, P, y) = self.deref_tpy()
# 	    #     if self.external:
# 	    #         self.var_or_expr(name="visc_liq",
# 	    #             expr=self.extf_mu_liq(T, y['CO2'], y['MEA']))
# 	    #         return
# 	    #     a = ViscLiqParams.a
# 	    #     b = ViscLiqParams.b
# 	    #     c = ViscLiqParams.c
# 	    #     d = ViscLiqParams.d
# 	    #     e = ViscLiqParams.e
# 	    #     f = ViscLiqParams.f
# 	    #     g = ViscLiqParams.g
# 	    #     r = self.wt_frac["MEA"]
#
# 	    #     self.var_or_expr(name="visc_liq",
# 	    #         expr=1.002e-3*10**((1.3272*(293.15 - T - 0.001053*(T-293.15)**2))/
# 	    #         (T-168.15))*exp(100*r*(T*(a*100*r + b) + c*100*r + d)*
# 	    #         (y["CO2"]/y["MEA"]*(e*100*r + f*T + g) + 1)/T**2))
#
# 	def __init__(self, data, metadata, outputFile, almtrace, description):
# 		PropertyModel.__init__(self, data, metadata, outputFile, almtrace, description);
#
#
# 	def results_str(self, model):
# 		print model.a.value, model.b.value, model.c.value, model.d.value, model.e.value, model.f.value, model.g.value;
# 		print "Sum of squared error = %f" % model.OBJ()
#
# 	def model_output(self, model, textfile):
# 		textfile.write("def _visc_liq(self): \n");
# 		textfile.write("\t(T, P, y) = self.deref_tpy()\n");
# 		for v in model.component_objects(Var, active = True):
# 			textfile.write("\t%s = %s\n" %(v.name, v.value))
# 		textfile.write("\tr = self.wt_frac[\"MEA\"]\n\n");
# 		expr = "1.002e-3*10**((1.3272*(293.15 - T - 0.001053*(T-293.15)**2))/(T-168.15))*exp(100*r*(T*(a*100*r + b) + c*100*r + d)*(y[\"CO2\"]/y[\"MEA\"]*(e*100*r + f*T + g) + 1)/T**2)"
# 		textfile.write("\tself.var_or_expr(name=\"visc_liq\", expr=%s)\n\n" %expr)
#
#
# 	def data_sets(self):
# 		self.I = np.asarray(range(self.data.num_rows));
# 		W_mea_vals = [v*100 for v in self.data.get_column('r').values]
# 		self.W_mea = dict(zip(self.I, W_mea_vals));
# 		self.Temp = dict(zip(self.I,self.data.get_column('T').values));
# 		self.alpha = dict(zip(self.I,self.data.get_column('CO2 Loading').values));
#
# 		# Measured value of interest
# 		self.u_sln = dict(zip(self.I,self.data.get_column('Viscosity Value').values));
#
# 	def sets(self, model):
# 		model.I = Set(initialize=self.I)
#
# 	def parameters(self, model):
# 		# Given parameters of the system
# 		model.alpha = Param(model.I, initialize=self.alpha)
# 		model.W_mea = Param(model.I, initialize=self.W_mea)
# 		model.T = Param(model.I, initialize=self.Temp)
#
# 		# Measured value of interest
# 		model.u_sln = Param(model.I, initialize=self.u_sln)
#
# 	def variables(self, model):
# 		model.a = Var(domain=Reals, initialize = 0)
# 		model.b = Var(domain= Reals, initialize = 0)
# 		model.c = Var(domain=Reals, initialize = 0)
# 		model.d = Var(domain=Reals, initialize = 0)
# 		model.e = Var(domain=Reals, initialize = 0)
# 		model.f = Var(domain=Reals, initialize = 0)
# 		model.g = Var(domain=Reals, initialize = 0)
#
# 	def objective(self, model):
#
# 		def obj_expression(model):
# 				expr = 0;
# 				for i in model.I:
# 					u_H2O = 1.002*10**((1.3272*(293.15 - model.T[i] - 0.001053*(model.T[i]-293.15)**2))/(model.T[i]-168.15))
# 					expr+= (model.u_sln[i] - u_H2O*exp((((model.a*model.W_mea[i]+model.b)*model.T[i] + model.c*model.W_mea[i] + model.d)*(model.alpha[i]*(model.e*model.W_mea[i] + model.f*model.T[i] + model.g) + 1)* model.W_mea[i])/(model.T[i]**2)))**2
# 				return expr
#
# 		model.OBJ = Objective(rule=obj_expression)
#
# 	def alamo_datasets(self,alamo_model):
# 		# Measured value of interest
# 		zdata = np.asarray(self.data.get_column('Viscosity Value').values);
# 		zdata = np.log(zdata)
# 		#print zdata
# 		self.alamo_model.zdata = zdata.reshape((len(zdata),1))
#
# 		xdata = []
# 		xdata.append(self.data.get_column('T').values)
# 		xdata.append(self.data.get_column('CO2 Loading').values)
# 		W_mea_vals = [v*100 for v in self.data.get_column('r').values]
# 		xdata.append(W_mea_vals)
#
# 		u_H2O = [1.002*10**((1.3272*(293.15 - t - 0.001053*(t-293.15)**2))/(t-168.15)) for t in self.data.get_column('T').values]
# 		# xdata.append(u_H2O)
#
# 		self.alamo_model.xdata = np.transpose(xdata)
# 		self.alamo_model.enumerate_datasets(3, ['T','a','W'], self.alamo_model.xdata)
# 		# self.alamo_model.enumerate_datasets(3, ['T','y[\"CO2\"]/y[\"MEA\"]','r'], self.alamo_model.xdata)
#
# 	def alamo_output(self, alamo_model, textfile):
# 		textfile.write("def _visc_liq(self): \n");
# 		textfile.write("\t# ALAMO Model\n")
# 		textfile.write("\t(T, P, y) = self.deref_tpy()\n");
# 		# for v in model.component_objects(Var, active = True):
# 		# 	textfile.write("\t%s = %s\n" %(v.name, v.value))
# 		textfile.write("\tr = self.wt_frac[\"MEA\"]\n\n");
# 		expr = "1.002e-3*10**((1.3272*(293.15 - T - 0.001053*(T-293.15)**2))/(T-168.15))*exp(100*r*(T*(a*100*r + b) + c*100*r + d)*(y[\"CO2\"]/y[\"MEA\"]*(e*100*r + f*T + g) + 1)/T**2)"
# 		textfile.write("\tself.var_or_expr(name=\"visc_liq\", expr=%s)\n\n" %alamo_model['model'])
#
# # for method_name in dir(EnthalpyModel):
# # 	print method_name
# # help(EnthalpyModel.alamo_datasets)
#
# # SteamPropertyPackage()
#
# MEAPropertyPackage()
#
#
#
#
# ########### SNIPPETS
#
# 	# def variables(self, model):
# 	# 	# Regressed parameter values
# 	# 	model.a_mea = Var(domain=Reals, initialize = 1)
# 	# 	model.b_mea = Var(domain= Reals, initialize = 1)
# 	# 	model.c_mea = Var(domain=Reals, initialize = 1)
#
# 	# 	#model.a_CO2 = Var(domain=Reals, initialize = 0)
# 	# 	#model.b_CO2 = Var(domain= Reals, initialize = 0)
# 	# 	model.c_CO2 = Var(domain=Reals, initialize = 1)
#
# 	# 	#model.d = Var(domain=Reals, initialize = 1)
# 	# 	#model.e = Var(domain=Reals, initialize = 1)
# 	# 	#model.V_CO2 = Var(domain=NonNegativeReals, initialize = 1)
# 	# 	model.Vstar = Var(domain=Reals, initialize = -1)
#
# 	# def objective(self, model):
#
# 	# 	def obj_expression(model):
# 	# 		expr = 0
#
# 	# 		# Molecular Weight
# 	# 		MW_mea = 61.08
# 	# 		MW_H2O = 18.02
# 	# 		MW_CO2 = 44.01
#
# 	# 		for i in model.I:
#
# 	# 			# Mol Fractions
# 	# 			X_mea = 1/(1+ model.alpha[i] + (MW_mea/MW_H2O)*((1-model.r[i])/model.r[i]))
# 	# 			X_CO2 = X_mea * model.alpha[i]
# 	# 			X_H2O =  1 - X_mea - X_CO2
# 	# 			MW_sln = X_mea*MW_mea + X_CO2*MW_CO2 + X_H2O*MW_H2O
#
# 	# 			# Volume
# 	# 			V_mea = MW_mea/((model.a_mea/1e7)*model.T[i]**2 + (model.b_mea/1e4)*model.T[i] + model.c_mea)
# 	# 			V_CO2 = MW_CO2/model.c_CO2;
# 	# 			#V_CO2 = MW_CO2/((model.a_CO2/1e7)*model.T[i]**2 + (model.b_CO2/1e4)*model.T[i] + model.c_CO2)
# 	# 			V_H2O = MW_H2O/(-(32.484/1e7)*model.T[i]**2+ (16.5/1e4)*model.T[i] + 0.793)
# 	# 			#Vstar2 = model.d + model.e*X_mea
#
# 	# 			#V_sln = X_mea*V_mea + X_H2O*V_H2O + X_CO2*model.V_CO2 + X_mea*X_H2O*model.Vstar # + X_mea*X_CO2*Vstar2
# 	# 			V_sln = X_mea*V_mea + X_H2O*V_H2O + X_CO2*V_CO2 + X_mea*X_H2O*model.Vstar # + X_mea*X_CO2*Vstar2
#
# 	# 			expr+= (model.rho_sln[i] - MW_sln/V_sln)**2
# 	# 		return expr;
#
# 	# 	model.OBJ = Objective(rule=obj_expression)
#
