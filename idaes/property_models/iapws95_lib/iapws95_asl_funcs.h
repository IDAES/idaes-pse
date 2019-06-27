/*------------------------------------------------------------------------------
 Institute for the Design of Advanced Energy Systems Process Systems
 Engineering Framework (IDAES PSE Framework) Copyright (c) 2018, by the
 software owners: The Regents of the University of California, through
 Lawrence Berkeley National Laboratory,  National Technology & Engineering
 Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
 University Research Corporation, et al. All rights reserved.

 Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
 license information, respectively. Both files are also available online
 at the URL "https://github.com/IDAES/idaes".
------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------
 This provides the ASL function interface for IAPWS R6-95(2016) property
 functions.

 Author: John Eslick
 File: iapws95_asl_funcs.h
------------------------------------------------------------------------------*/

#include <funcadd.h>

#ifndef _INCLUDE_IAPWS95_ASL_FUNCS_H_
#define _INCLUDE_IAPWS95_ASL_FUNCS_H_

double p_asl(arglist *al);
double u_asl(arglist *al);
double s_asl(arglist *al);
double h_asl(arglist *al);
double g_asl(arglist *al);
double f_asl(arglist *al);
double cv_asl(arglist *al);
double cp_asl(arglist *al);
double w_asl(arglist *al);

double hvpt_asl(arglist *al);
double hlpt_asl(arglist *al);
double vf_asl(arglist *al);
double tau_asl(arglist *al);

double delta_sat_l_asl(arglist *al);
double delta_sat_v_asl(arglist *al);
double p_sat_asl(arglist *al);
double tau_sat_asl(arglist *al);

double delta_liq_asl(arglist *al);
double delta_vap_asl(arglist *al);

double phi0_asl(arglist *al);
double phi0_delta_asl(arglist *al);
double phi0_delta2_asl(arglist *al);
double phi0_tau_asl(arglist *al);
double phi0_tau2_asl(arglist *al);
double phir_asl(arglist *al);
double phir_delta_asl(arglist *al);
double phir_delta2_asl(arglist *al);
double phir_tau_asl(arglist *al);
double phir_tau2_asl(arglist *al);
double phir_delta_tau_asl(arglist *al);

#endif
