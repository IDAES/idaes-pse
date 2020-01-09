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
 This file provides some configuration parameters.

 Author: John Eslick
 File: iapws95_config.h
------------------------------------------------------------------------------*/

#ifndef _INCLUDE_IAPWS95_CONFIG_H_
#define _INCLUDE_IAPWS95_CONFIG_H_

#include<cmath>

// max memo table size (0 deactivates memoization)
#define MAX_MEMO 100000000
// Precision: {PRECISION_LONG_DOUBLE, PRECISION_DOUBLE} the exact meaning of
// that depends on the machine and compiler.
//#define PRECISION_LONG_DOUBLE
#define PRECISION_DOUBLE
// Redidual abs tolerance for solving for vapor reduced density from p, tau
#define TOL_DELTA_VAP 1e-11
// Residual abs tolerance for solving for liquid reduced density from p, tau
#define TOL_DELTA_LIQ 1e-11
// Max iterations for solving for delta (reduced density) from p, tau
#define MAX_IT_DELTA 20
// Use bracketing methods in particularly difficult areas when solving for
// density from temperature and pressure.  This lets me get a very accurate
// initial guess that I just feed to the newton type solve.
// Bracketing methods tolerance, for absolute error on delta (reduced density)
#define TOL_BRACKET 1e-9
// Bracketing methods iteration limit
#define MAX_IT_BRACKET 30
// Saturation curve relative tolerances for phase Gibbs free enegy difference
#define TOL_REL_SAT_G 1e-11
// Saturation curve max iterations
#define MAX_IT_SAT 10
// Saturation solver gamma factor Akasaka (2008)
#define SAT_GAMMA 1.0
//Parameters for solving for tau_sat as a function of pressure
#define MAX_IT_SAT_TAU 20
#define TOL_SAT_TAU 1e-11

// Set types and functions that specify precision
#ifdef PRECISION_DOUBLE
  typedef double s_real;
  #define s_exp exp
  #define s_log log
  #define s_pow pow
  #define s_sqrt sqrt
  #define s_fabs fabs
#endif

#ifdef PRECISION_LONG_DOUBLE
  typedef long double s_real;
  #define s_exp expl
  #define s_log logl
  #define s_pow powl
  #define s_sqrt sqrtl
  #define s_fabs fabsl
#endif


#endif
