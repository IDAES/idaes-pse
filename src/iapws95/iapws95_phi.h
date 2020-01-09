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
 IAPWS R6-95(2016) dimensionless Helmholtz free enegy calcuations and up to 4th
 derivatives phi0 is the ideal part phir is the residual part. Memoization of all
 phir functions is included. Non-analytic terms (55 and 56) were not included.
 Accuracy is slightly decreased near the critical point, but the function is
 much better behaved, so it will work better for optimization.

 Author: John Eslick
 File: iapws95_phi.h
------------------------------------------------------------------------------*/

#include "iapws95_config.h"

#ifndef _INCLUDE_IAPWS95_PHI_H_
#define _INCLUDE_IAPWS95_PHI_H_


// Ideal part and derivatives of dimensionless Helmholtz
s_real phi0(s_real delta, s_real tau); //Ideal part of dimensionless Helmholtz
s_real phi0_delta(s_real delta); //derivative of phi0 wrt delta, delta=rho/rho_c
s_real phi0_delta2(s_real delta); //2nd derivative of phi0 wrt delta
s_real phi0_tau(s_real tau); //derivative of phi0 wrt tau, tau=T_c/T
s_real phi0_tau2(s_real tau); //2nd derivative of phi0 wrt tau
//phi0_delta_tau is 0 so don't need a function for it

s_real phi0_delta3(s_real delta);
s_real phi0_delta4(s_real delta);
s_real phi0_tau3(s_real tau);
s_real phi0_tau4(s_real tau);

// Residual part and derivatives of dimensionless Helmholtz
s_real phir(s_real delta, s_real tau);
s_real phir_delta(s_real delta, s_real tau);
s_real phir_delta2(s_real delta, s_real tau);
s_real phir_tau(s_real delta, s_real tau);
s_real phir_tau2(s_real delta, s_real tau);
s_real phir_delta_tau(s_real delta, s_real tau);
s_real phir_delta3(s_real delta, s_real tau);
s_real phir_delta4(s_real delta, s_real tau);
s_real phir_delta2_tau(s_real delta, s_real tau);
s_real phir_delta2_tau2(s_real delta, s_real tau);
s_real phir_delta_tau2(s_real delta, s_real tau);
s_real phir_delta_tau3(s_real delta, s_real tau);
s_real phir_tau3(s_real delta, s_real tau);
s_real phir_tau4(s_real delta, s_real tau);
s_real phir_delta3_tau(s_real delta, s_real tau);


/*------------------------------------------------------------------------------
  Some phi0 derivatives that are 0, just to be more explicit in calculations
  phi0 is the ideal part of dimensionless Helmholtz free energy
------------------------------------------------------------------------------*/
static const s_real phi0_delta_tau = 0;  //derivative of phi0 wrt delta and tau
static const s_real phi0_delta2_tau = 0;
static const s_real phi0_delta_tau2 = 0;
static const s_real phi0_delta3_tau = 0;
static const s_real phi0_delta2_tau2 = 0;
static const s_real phi0_delta_tau3 = 0;

#endif
