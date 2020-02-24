##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
# =============================================================================

from pyomo.environ import (Constraint,
                           Expression,
                           Var,
                           Param,
                           exp,
                           log10)

def global_costing_parameters(self, year=None):
    if year is None:
        year='2018'
    # Cost index $/year (method argument or 2018 default)
    ce_index_dic = {'2019':680, '2018':671.1, '2017':567.5, '2016':541.7,
                '2015':556.8, '2014':576.1, '2013':567.3, '2012':584.6,
                '2011':585.7, '2010':550.8}
    
    self.CE_index = Param(mutable = True, initialize = ce_index_dic[year],
                          doc= 'CE cost index $ year')

def _make_vars(self):
    self.base_cost = Var(initialize = 1e5,
                         bounds=(0, 1e12),
                         doc='Heat Exchanger Base Cost CB cost in $')
    self.purchase_cost = Var(initialize =1e4,
                             bounds = (0,1e12),
                             doc = 'Heat Exchanger Purchase Cost CP in $')
    self.FP = Var(initialize = 100,
                  doc = 'Pressure factor')
    
    self.FM = Param(mutable = True, initialize = 3,
                    doc= 'construction material correction factor')
    
    
def hx_costing(self, hx_type='U-tube', FM='stain_steel', 
                 FL = '12ft'):
    '''
    Heat exchanger costing method
    Source: Process and Product Design Principles: Synthesis, Analysis, 
    and Evaluation
    Seider, Seader, Lewin, Windagdo, 3rd Ed. John Wiley and Sons
    Chapter 22. Cost Accounting and Capital Cost Estimation
    22.2 Cost Indexes and Capital Investment
    
    This method computes the purchase cost (CP) for a shell and tube heat
    exchanger (Eq. 22.43), the model computes the base cost (CB for 4 type 
    of heat exchangers, such as floating head, fixed head, U-tube, and 
    Kettle vaporizer), construction material factor (FM), pressure design 
    factor (FP), and tube length correction factor (FL), using CE base cost 
    index of 500. 
    
    Cp = FP*FM*FL*CB
    
    Args:
        hx_type : 'floating_head', 'fixed_head', 'U-tube'*, 'Kettle_vap'
        
        material factor (FM): 'stain_steel'*, 'carb_steel'
        
        tube length (FL): '8ft'*, '12ft', '16ft', '20ft'
        
        * --> default options
    
    Author: MZamarripa (created 9/17/2019; updated 2/21/2020)
    '''
    # ------------------------ Heat Exchanger cost ------------------------
    # heat exchanger cost

    _make_vars(self)
    self.FL = Param(mutable = True, initialize = 1.12,
                    doc='hx tube length correction factor')
    
    self.hx_os = Param(mutable = True, initialize = 1.1,
                       doc= 'hx oversize factor 1.1 to 1.5')

    # select length correction factor 
    c_fl = {'8ft':1.25, '12ft':1.12, '16ft':1.05, '20ft':1.00}
    self.FL = c_fl[FL]

    # select construction material (stainless steel default)
    c_material = {'stain_steel':3, 'carbon_steel':2}
    self.FM = c_material[FM]
    # ToDo (9/17/2019):FM can be calculated for advanced materials Eq.22.44 
    # i.e. Carbon Steel/Cr-Mo steel, Cr-Mo steel/Cr-Mo steel, Monel/Monel, 
    # Titanium/titanium, etc.
    
    #--------------------------------------------------
    # base cost calculation
#       # select heat exchanger type:
    alf1 = {'floating_head':11.9052,'fixed_head':11.2927,
            'U-tube':11.3852,'Kettle_vap':12.2052}
    alf2 = {'floating_head':0.8709,'fixed_head':0.8228,
            'U-tube':0.9186,'Kettle_vap':0.8709}
    alf3 = {'floating_head':0.09005,'fixed_head':0.09861,
            'U-tube':0.09790,'Kettle_vap':0.09005}
    
    if (self.parent_block().config.tube.property_package.get_metadata().
            default_units['length']) == 'm':
        area = self.parent_block().area*10.7639
    elif (self.parent_block().config.tube.property_package.get_metadata().
            default_units['length']) == 'ft':
        area = self.area
    else:
        raise Exception('area units not suported')

    def hx_cost_rule(self):
        return self.base_cost == self.FL * exp(alf1[hx_type]
                    - alf2[hx_type]*log10(area*self.hx_os)
                    + alf3[hx_type]*log10(area*self.hx_os)**2)
    self.base_cost_eq = Constraint(rule=hx_cost_rule)
    
    #------------------------------------------------------
    # Pressure factor calculation
    # doublecheck units (higher pressure fluid should be tube side)
    if (self.parent_block().config.tube.property_package.get_metadata().
            properties['pressure']['units']) == 'Pa':
        pressure = self.parent_block().tube.properties_in[0].pressure*14.69/1.01325e5
    elif (self.parent_block().config.tube.property_package.get_metadata().
            properties['pressure']['units']) == 'psig':
        pressure = self.tube.properties_in[0].pressure

    #units must be in psig 
    #(ToDo(9/17/2019): check inlet pressure to select correlation)
    def hx_P_factor(self):
        # eq valid from 600 pisg to 3000 psig
        #    return self.hx_FP == 0.8510 + 0.1292*(pressure/600) 
        #                                + 0.0198*(pressure/600)**2
        # eq valid from 100 pisg to 2000 psig
        return self.FP == 0.9803 + (0.018*(pressure*14.69/1.01325e5/100)
                                    + 0.0017*(pressure*14.69/1.01325e5/100)**2)
    self.p_factor_eq = Constraint(rule=hx_P_factor)


    #---------------------------------------------------------
    # purchase cost equation
    def hx_CP_rule(self):
        return self.purchase_cost == (self.FP*self.FM*self.FL*
                                      (self.parent_block().flowsheet().costing.CE_index/500)*self.base_cost)
    self.cp_cost_eq = Constraint(rule=hx_CP_rule)


def pressure_changer_costing(self, FM="stain_steel", # pump, compressor
                             mover_type = "compressor", #fan, blower, compressor
                             compressor_type="centrifugal", # only compressors
                             driver_mover_type: "electrical" # only compressors
                             ):
    '''    
    Pressure changer costing method
    Source: Process and Product Design Principles: Synthesis, Analysis, 
    and Evaluation
    Seider, Seader, Lewin, Windagdo, 3rd Ed. John Wiley and Sons
    Chapter 22. Cost Accounting and Capital Cost Estimation
    22.2 Cost Indexes and Capital Investment
    
    This method computes the purchase cost (CP) for a pressure changer unit, 
    and can be used for costing of pumps, fan, blower, compressor, or Turbine.  
    
    Author: MZamarripa (created 2/21/2020)
    
    Arguments: 

        if config.compressor 
        if config.compressor is True:
        mover_type: Fan, Blower, Compressor*
            Fan_type:
            Blower_type:
            Compressor_type: centrifugal*, reciprocating, screw
                drive_type: electric_motor*, steam_turbine, gas_turbine
                
        FM - construction material: 
            Pump: stain_steel*
            Compressor: stain_steel*, nickel_alloy 
    *default option
    '''
    _make_vars(self)
    
    # if compressor is == False, that means pressure changer in a Turbine
    if self.parent_block().config.compressor == False:
        # turbine purchase cost calculation W / 746 = horsepower
        if(self.parent_block().config.property_package.get_metadata().
           properties['enth_mol']['units']) == 'J/mol':
            def CP_rule(self):
                return self.purchase_cost == 530*(-1*self.parent_block().work_mechanical[0]/746)**0.81
            self.cp_cost_eq = Constraint(rule=CP_rule)
        else:
            raise Exception("units not supported")
    #if compressor is True (the unit could be a Pump, Fan, Blower, or Compressor)
    else:
        #if ThermodynamicAssumption == Pump (the unit is a Pump)
        if self.parent_block().config.ThermodynamicAssumption == ThermodynamicAssumption.pump:
            # ToDo: pump costing code
            print('pump code')
        else:
            if mover_type == "compressor":
                # Compressor Purchase Cost Correlation
                # compressor drive factor
                FD_param = {'electrical_motor':1, 'steam_turbine':1.15, 'gas_turbine':1.25}
                self.FD = Param(mutable = True, 
                                initialize = FD_param[driver_mover_type],
                                doc='drive factor')
                # compressor material factor
                FM_comp = {'stain_steel':2.5, 'nickel_alloy':5.0}
                self.FM = FM_comp[FM]
                c_alf1 = {'centrifugal':7.58,
                          'reciprocating':7.9661, 
                          'screw':8.1238}
                c_alf2 = {'centrifugal':0.8,
                          'reciprocating':0.8, 
                          'screw':0.7243}
                #Purchase cost rule
                def CP_rule(self):
                    return self.purchase_cost == exp(alf1[compressor_type])
                self.cp_cost_eq = Constraint(rule=CP_rule)

            elif mover_type == "Fan":
                # fan cost correlation
            elif mover_type == "Blower":
                #Blower Cost Correlation
            else:
                raise Exception('mover type not supported')
        
    