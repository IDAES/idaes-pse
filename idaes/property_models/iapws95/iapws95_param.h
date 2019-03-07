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
 This file provides the IAPWS R6-95(2016) parameters and a few more parameters
 related to the solution methods. The iX() functions shift the indexes for
 the arrays to match the IAPWS calculations without allocating extra storage.

 References:
   International Association for the Properties of Water and Steam (2016).
       IAPWS R6-95 (2016), "Revised Release on the IAPWS Formulation 1995 for
       the Properties of Ordinary Water Substance for General Scientific Use,"
       URL: http://iapws.org/relguide/IAPWS95-2016.pdf
   Wagner, W.,  A. Pruss (2002). "The IAPWS Formulation 1995 for the
       Thermodynamic Properties of Ordinary Water Substance for General and
       Scientific Use." J. Phys. Chem. Ref. Data, 31, 387-535.
   Wagner, W. et al. (2000). "The IAPWS Industrial Formulation 1997 for the
       Thermodynamic Properties of Water and Steam," ASME J. Eng. Gas Turbines
       and Power, 122, 150-182.
   Akasaka, R. (2008). "A Reliable and Useful Method to Determine the Saturation
       State from Helmholtz Energy Equations of State." Journal of Thermal
       Science and Technology, 3(3), 442-451.

 Author: John Eslick
 File: iapws95_param.h
------------------------------------------------------------------------------*/

#ifndef _INCLUDE_IAPWS95_PARAM_H_
#define _INCLUDE_IAPWS95_PARAM_H_

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

// Critiacal point for water and R
static const s_real T_c = 647.096;   // Critical T (K)
static const s_real rho_c = 322;     // Critical density (kg/m^3)
static const s_real R = 0.46151805;  // Gas constant (kJ/kg/K)
static const s_real P_c = 2.2064e4;  // Critical Pressure (kPa)

// functions to convert indexes from IAPWS equations to actual offsets
// (to save wasted memory and keep the equations understandable.)
inline unsigned char i1(unsigned char i){return i - 1;}
inline unsigned char i4(unsigned char i){return i - 4;}
inline unsigned char i8(unsigned char i){return i - 8;}
inline unsigned char i52(unsigned char i){return i - 52;}
inline unsigned char i55(unsigned char i){return i - 55;}

// Some phi0 derivatives that are 0, just to be more explicit in calculations
// phi0 is the ideal part of dimensionless Helmholtz free energy
static const s_real phi0_delta_tau = 0;  //derivative of phi0 wrt delta and tau
static const s_real phi0_delta2_tau = 0;
static const s_real phi0_delta_tau2 = 0;
static const s_real phi0_delta3_tau = 0;
static const s_real phi0_delta2_tau2 = 0;
static const s_real phi0_delta_tau3 = 0;

// Functions to convert to dimensionless state variables from T and Rho.
// temperature and density
inline s_real f_tau(s_real T){return T_c/T;}
inline s_real f_delta(s_real rho){return rho/rho_c;}


// psat curve parameters from IAPWS-97 (industial formulation)
// Use this as an initial guess when solving for the saturation curve using the
// more consitent set of IAPWS-95 (scientific formulation) equations.
static const s_real n_psat[] = {
   0.11670521452767e4, //1
  -0.72421316703206e6, //2
  -0.17073846940092e2, //3
   0.12020824702470e5, //4
  -0.32325550322333e7, //5
   0.14915108613530e2, //6
  -0.48232657361591e4, //7
   0.40511340542057e6, //8
  -0.23855557567849,   //9
   0.65017534844798e3  //10
};

//
// Constants from IAPWS95 (R6 2016) from here to end
//

static const s_real n0[] = {
    -8.3204464837497,  //1 (starts at 1, i1(1) = 0)
     6.6832105275932,  //2
     3.00632,          //3
     0.012436,         //4
     0.97315,          //5
     1.27950,          //6
     0.96956,          //7
     0.24873};         //8 (starts at 1, i1(8) = 7)

static const s_real gamma0[] = {
    1.28728967,  //4 (starts at 4, i4(4) = 0)
    3.53734222,  //5
    7.74073708,  //6
    9.24437796,  //7
    27.5075105}; //8

static const unsigned char c[] = {
    1,   //8 (starts at 8, i8(8) = 0)
    1,   //9
    1,   //10
    1,   //11
    1,   //12
    1,   //13
    1,   //14
    1,   //15
    1,   //16
    1,   //17
    1,   //18
    1,   //19
    1,   //20
    1,   //21
    1,   //22
    2,   //23
    2,   //24
    2,   //25
    2,   //26
    2,   //27
    2,   //28
    2,   //29
    2,   //30
    2,   //31
    2,   //32
    2,   //33
    2,   //34
    2,   //35
    2,   //36
    2,   //37
    2,   //38
    2,   //39
    2,   //40
    2,   //41
    2,   //42
    3,   //43
    3,   //44
    3,   //45
    3,   //46
    4,   //47
    6,   //48
    6,   //49
    6,   //50
    6};  //51

static const unsigned char d[] = {
     1,  //1
     1,  //2
     1,  //3
     2,  //4
     2,  //5
     3,  //6
     4,  //7
     1,  //8
     1,  //9
     1,  //10
     2,  //11
     2,  //12
     3,  //13
     4,  //14
     4,  //15
     5,  //16
     7,  //17
     9,  //18
    10,  //19
    11,  //20
    13,  //21
    15,  //22
     1,  //23
     2,  //24
     2,  //25
     2,  //26
     3,  //27
     4,  //28
     4,  //29
     4,  //30
     5,  //31
     6,  //32
     6,  //33
     7,  //34
     9,  //35
     9,  //36
     9,  //37
     9,  //38
     9,  //39
     10, //40
     10, //41
     12, //42
      3, //43
      4, //44
      4, //45
      5, //46
     14, //47
      3, //48
      6, //49
      6, //50
      6, //51
      3, //52
      3, //53
      3}; //54

static const s_real t[] = {
    -0.5,    //1
     0.875,  //2
     1,      //3
     0.5,    //4
     0.75,   //5
     0.375,  //6
     1,      //7
     4,      //8
     6,      //9
    12,      //10
     1,      //11
     5,      //12
     4,      //13
     2,      //14
    13,      //15
     9,      //16
     3,      //17
     4,      //18
    11,      //19
     4,      //20
    13,      //21
     1,      //22
     7,      //23
     1,      //24
     9,      //25
    10,      //26
    10,      //27
     3,      //28
     7,      //29
    10,      //30
    10,      //31
     6,      //32
    10,      //33
    10,      //34
     1,      //35
     2,      //36
     3,      //37
     4,      //38
     8,      //39
     6,      //40
     9,      //41
     8,      //42
    16,      //43
    22,      //44
    23,      //45
    23,      //46
    10,      //47
    50,      //48
    44,      //49
    46,      //50
    50,      //51
     0,      //52
     1,      //53
     4};     //54

static const s_real n[] = {
     0.12533547935523e-1,  //1 (starts at 1, i1(1) = 0)
     0.78957634722828e1,   //2
    -0.87803203303561e1,   //3
     0.31802509345418,     //4
    -0.26145533859358,     //5
    -0.78199751687981e-2,  //6
     0.88089493102134e-2,  //7
    -0.66856572307965,     //8
     0.20433810950965,     //9
    -0.66212605039687e-4,  //10
    -0.19232721156002,     //11
    -0.25709043003438,     //12
     0.16074868486251,     //13
    -0.40092828925807e-1,  //14
     0.39343422603254e-6,  //15
    -0.75941377088144e-5,  //16
     0.56250979351888e-3,  //17
    -0.15608652257135e-4,  //18
     0.11537996422951e-8,  //19
     0.36582165144204e-6,  //20
    -0.13251180074668e-11, //21
    -0.62639586912454e-9,  //22
    -0.10793600908932,     //23
     0.17611491008752e-1,  //24
     0.22132295167546,     //25
    -0.40247669763528,     //26
     0.58083399985759,     //27
     0.49969146990806e-2,  //28
    -0.31358700712549e-1,  //29
    -0.74315929710341,     //30
     0.47807329915480,     //31
     0.20527940895948e-1,  //32
    -0.13636435110343,     //33
     0.14180634400617e-1,  //34
     0.83326504880713e-2,  //35
    -0.29052336009585e-1,  //36
     0.38615085574206e-1,  //37
    -0.20393486513704e-1,  //38
    -0.16554050063743e-2,  //39
     0.19955571979541e-2,  //40
     0.15870308324157e-3,  //41
    -0.16388568342530e-4,  //42
     0.43613615723811e-1,  //43
     0.34994005463765e-1,  //44
    -0.76788197844621e-1,  //45
     0.22446277332006e-1,  //46
    -0.62689710414685e-4,  //47
    -0.55711118565645e-9,  //48
    -0.19905718354408,     //49
     0.31777497330738,     //50
    -0.11841182425981,     //51
    -0.31306260323435e2,   //52
     0.31546140237781e2,   //5
    -0.25213154341695e4,   //54
    -0.14874640856724,     //55
     0.31806110878444};    //56

static const s_real alpha[] =      {20,    20,    20   }; //52, 53, 54
static const s_real gamma_[] =     { 1.21,  1.21,  1.25}; //52, 53, 54
static const unsigned char eps[] = { 1,     1,     1   }; //52, 53, 54

static const s_real beta[] = {150, 150, 250, 0.3, 0.3}; //52, 53, 54, 55, 56

static const s_real a[] = {  3.5,    3.5 };  //55, 56
static const s_real b[] = {  0.85,   0.95};  //55, 56
static const s_real A[] = {  0.32,   0.32};  //55, 56
static const s_real B[] = {  0.2,    0.2 };  //55, 56
static const s_real C[] = { 28,     32   };  //55, 56
static const s_real D[] = {700,    800   };  //55, 56

#endif
