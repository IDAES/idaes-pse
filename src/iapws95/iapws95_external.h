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
 Main IAPWS R6-95(2016) steam calculations. For now, the non-analytic terms are
 not included, because they cause a singularity at the critical point.  It is
 assumed that we will generally not be operating very near the critical point,
 so well behaved calculations are prefered to high accuracy in near the critical
 point.

 For references see iapws95.h

 Author: John Eslick
 File iapws95.h
------------------------------------------------------------------------------*/

#include "iapws95_config.h"
#include "iapws95_phi.h"

#ifndef _INCLUDE_IAPWS95_H_
#define _INCLUDE_IAPWS95_H_

/*------------------------------------------------------------------------------
  Basic thermodynamic quantities as functions of delta (rho/rho_c) and
  tau (T_c/T)
------------------------------------------------------------------------------*/
s_real p(s_real delta, s_real tau);
s_real u(s_real delta, s_real tau);
s_real s(s_real delta, s_real tau);
s_real h(s_real delta, s_real tau);
s_real f(s_real delta, s_real tau);
s_real cv(s_real delta, s_real tau);
s_real cp(s_real delta, s_real tau);
s_real g(s_real delta, s_real tau);
s_real w(s_real delta, s_real tau);

/*------------------------------------------------------------------------------
  Function to calculate delta (rho/rho_c) from P and tau (T_c/T)
------------------------------------------------------------------------------*/

// Solve for delta using Halley's method delta0 in an inital guess, tol is
// the absolute residual tolerance, nit provides a means to return number of
// iterations for testing there is more than 1 solution maybe even more than one
// physically meaningfull solution, so delta0 is very important
s_real delta_p_tau(s_real p, s_real tau, s_real delta_0, s_real tol=1e-10,
                   int *nit=NULL, s_real *grad=NULL, s_real *hes=NULL);

/*------------------------------------------------------------------------------
  Functions for saturation pressure and density as a function to tau (T_c/T).
  These use the method of Akasaka (2008), which is good at least up to 647.093 K
  In this case the last bit of the density curves from 647.090K to 647.096K was
  filled in with cubics that match density, first and second derivatives at the
  transition point and critical density at the critical point.
------------------------------------------------------------------------------*/
s_real sat_delta_liq(s_real tau); // saturated liquid density at tau
s_real sat_delta_vap(s_real tau); // saturated vapor density at tau
s_real sat_p(s_real tau);         // saturation pressure at tau
int sat(s_real tau, s_real *delta_l_sol, s_real *delta_v_sol); //sat solver

// Initial guesses for phase densities on saturation curve
s_real delta_sat_l_approx(s_real tau);
s_real delta_sat_v_approx(s_real tau);

/*------------------------------------------------------------------------------
  Basic thermodynamic functions with gradient and Hessian as functions of
  delta (rho/rho_c) and tau (T_c/T)
------------------------------------------------------------------------------*/
s_real p_with_derivs(s_real delta, s_real tau, s_real *grad, s_real *hes);
s_real u_with_derivs(s_real delta, s_real tau, s_real *grad, s_real *hes);
s_real h_with_derivs(s_real delta, s_real tau, s_real *grad, s_real *hes);
s_real s_with_derivs(s_real delta, s_real tau, s_real *grad, s_real *hes);
s_real f_with_derivs(s_real delta, s_real tau, s_real *grad, s_real *hes);
s_real g_with_derivs(s_real delta, s_real tau, s_real *grad, s_real *hes);
s_real cp_with_derivs(s_real delta, s_real tau, s_real *grad, s_real *hes);
s_real cv_with_derivs(s_real delta, s_real tau, s_real *grad, s_real *hes);
s_real w_with_derivs(s_real delta, s_real tau, s_real *grad, s_real *hes);

s_real hvpt_with_derivs(s_real pr, s_real tau, s_real *grad, s_real *hes);
s_real hlpt_with_derivs(s_real pr, s_real tau, s_real *grad, s_real *hes);
s_real vf_with_derivs(s_real ht, s_real pr, s_real *grad, s_real *hes);
s_real tau_with_derivs(s_real ht, s_real pr, s_real *grad, s_real *hes);

/*------------------------------------------------------------------------------
  Functions for saturation pressure and reduced density as a function of tau
  (T_c/T) with gradient and Hessian.
------------------------------------------------------------------------------*/
s_real sat_delta_liq_with_derivs(s_real tau, s_real *grad, s_real *hes);
s_real sat_delta_vap_with_derivs(s_real tau, s_real *grad, s_real *hes);
s_real sat_p_with_derivs(s_real tau, s_real *grad, s_real *hes, bool limit=1);
s_real sat_tau_with_derivs(s_real pr, s_real *grad, s_real *hes, int *nit=NULL);

/*------------------------------------------------------------------------------
  Functions for saturation delta (rho/roh_c) as a function of P and tau (T_c/T)
  with gradient and Hessian.  These functions are well tested and just feed the
  right initial guess to delta_p_tau() gauranteeing the correct solution.
------------------------------------------------------------------------------*/
s_real delta_vap(s_real p, s_real tau, s_real *grad=NULL, s_real *hes=NULL, int *nit=NULL);
s_real delta_liq(s_real p, s_real tau, s_real *grad=NULL, s_real *hes=NULL, int *nit=NULL);
#endif
