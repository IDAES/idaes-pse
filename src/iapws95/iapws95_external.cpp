/*------------------------------------------------------------------------------
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
------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------
 Main IAPWS R6-95(2016) steam calculations. For now, the non-analytic terms are
 not included, because they cause a singularity at the critical point.  It is
 assumed that we will generally not be operating very near the critical point,
 so well behaved calculations are prefered to high accuracy in near the critical
 point.

 For references see iapws95.h

 Author: John Eslick
 File iapws95.cpp
------------------------------------------------------------------------------*/
#include<stdio.h>
#include<cmath>
#include<iostream>

#include"iapws95_memo.h"
#include"iapws95_external.h"
#include"iapws95_config.h"
#include"iapws95_deriv_parts.h"

#define TAU_LOW 0.15
#define TAU_HIGH 4.0

void zero_derivs2(s_real *grad, s_real *hes){
  if(grad!=NULL){
    grad[0] = 0;
    grad[1] = 0;
    if(hes!=NULL){
      hes[0] = 0;
      hes[1] = 0;
      hes[2] = 0;
    }
  }
}

/*------------------------------------------------------------------------------
Basic Property Calculations

Your guide to variables in this section:
  T_c: critical temperature (K) defined in iapws95_param.h
  rho_c: critical density (kg/m3) defined in iapws95_param.h
  R: ideal gas constant (kJ/kg/K) defined in iapws95_param.h
  delta: reduced density rho/rho_c
  tau: 1/reduced temperature T_c/T
*-----------------------------------------------------------------------------*/

s_real p(s_real delta, s_real tau){ //pressure (kPa)
  return R*(delta*rho_c)*(T_c/tau)*(1 + delta*phir_delta(delta, tau));
}

s_real u(s_real delta, s_real tau){ // internal energy (kJ/kg)
  return R*T_c*(phi0_tau(tau) + phir_tau(delta, tau));
}

s_real s(s_real delta, s_real tau){ // entropy (kJ/kg/K)
  return R*(tau*(phi0_tau(tau) + phir_tau(delta, tau)) -
         phi0(delta, tau) - phir(delta, tau));
}

s_real h(s_real delta, s_real tau){ // enthalpy (kJ/kg)
  return R*(T_c/tau)*(1 + tau*(phi0_tau(tau) + phir_tau(delta, tau)) +
         delta*phir_delta(delta, tau));
}

s_real f(s_real delta, s_real tau){ // Helmholtz free energy (kJ/kg)
  return R*(T_c/tau)*(phi0(delta, tau) + phir(delta, tau));
}

s_real g(s_real delta, s_real tau){ // Gibbs free energy (kJ/kg)
  return h(delta, tau) - T_c/tau*s(delta, tau);
}

s_real cv(s_real delta, s_real tau){ // Constant volume heat capacity (kJ/kg/K)
  return -R*tau*tau*(phi0_tau2(tau) + phir_tau2(delta, tau));
}

s_real cp(s_real delta, s_real tau){ // Constant pressure heat capacity (kJ/kg/K)
  return cv(delta, tau) + R*XD/XE;
}

s_real w(s_real delta, s_real tau){ // Speed of sound (m/s)
  return s_sqrt(XB*1000*(XE + R*XD/cv(delta, tau)));
}

/*------------------------------------------------------------------------------
In this section you will find function for calculating density from pressure
and teperature. This lets you calculate properties in the more standard way as
a function of temperature and pressure instead of density and pressure.

Unfortunatly this is difficult and a bit messy.
*-----------------------------------------------------------------------------*/

s_real delta_p_tau_rf(s_real pr, s_real tau, s_real a, s_real b, bool bisect){
  /*----------------------------------------------------------------------------
  Bracketing methods, false position and bisection, for finding a better initial
  guess for density when solving for density from temperature and pressure. This
  is only used in particularly diffucult areas.  At this point it probably
  overused, but there are places where it is probably necessary.

  Args:
    pr: pressure (kPa)
    tau: inverse of reduced pressure Tc/T
    bisect: 1 to use bisection (probably slow), 0 to use false position
            bisection isn't really used, but I kept it around for debugging
    a: first density bound (kg/m3)  | really density, why didn't use delta?
    b: second density bound (kg/m3) | it was easier to think in terms of density
  Returns:
    delta: reduced density at pr and tau (more or less approximate)
  ----------------------------------------------------------------------------*/
  s_real c=a, fa, fb, fc;
  int it = 0;
  // If right by critical point guess critical density. (okay this isn't part
  // of a backeting method but it was conveneint).
  if( fabs(T_c/tau - T_c) < 1e-7 && fabs(pr - P_c) < 1e-4) return 1;
  // solve f(delta, tau) = 0; f(delta, tau) = p(delta, tau) - pr
  a /= rho_c; // convert to reduced density
  b /= rho_c; // convert to reduced density
  fa = p(a, tau) - pr; // initial f(a, tau)
  fb = p(b, tau) - pr; // initial f(b, tau)
  while(it < MAX_IT_BRACKET && (a - b)*(a - b) > TOL_BRACKET){
    if(bisect) c = (a + b)/2.0; // bisection
    else c = b - fb*(b - a)/(fb - fa); //regula falsi
    fc = p(c, tau) - pr; // calcualte f(c)
    if(fc*fa >= 0){a = c; fa = fc;}
    else{b = c; fb = fc;}
    ++it;
  }
  return (a+b)/2.0;
}

s_real delta_p_tau(s_real pr, s_real tau, s_real delta_0, s_real tol, int *nit,
  s_real *grad, s_real *hes){
  /*----------------------------------------------------------------------------
  Halley's method to calculate density from temperature and pressure

  Args:
   pr: pressure (kPa)
   tau: inverse of reduced pressure Tc/T
   delta0: initial guess for delta
   tol: absolute residual tolerance for convergence
   nit: pointer to return number of iterations to, or NULL
   grad: location to return gradient of delta wrt pr and tau or NULL
   hes: location to return hessian (upper triangle) of NULL
  Returns:
   delta: reduced density (accuracy depends on tolarance and function shape)
  ----------------------------------------------------------------------------*/
  s_real delta = delta_0, fun, gradp[2], hesp[3];
  int it = 0; // iteration count
  fun = p_with_derivs(delta, tau, gradp, hesp) - pr;
  while(fabs(fun) > tol && it < MAX_IT_DELTA){
    delta = delta - fun*gradp[0]/(gradp[0]*gradp[0] - 0.5*fun*hesp[0]);
    fun = p_with_derivs(delta, tau, gradp, hesp) - pr;
    ++it;
  }
  if(nit != NULL) *nit = it;
  if(grad != NULL){ // calculate gradient if needed
    grad[0] = 1.0/gradp[0];
    grad[1] = -gradp[1]*grad[0]; //triple product
    if(hes != NULL){ // calculate hession if needed.
      hes[0] = -hesp[0]*grad[0]/gradp[0]/gradp[0];
      hes[1] = -(hesp[1] + hesp[0]*grad[1])/gradp[0]/gradp[0];
      hes[2] = -(grad[0]*(hesp[2] + grad[1]*hesp[1]) + gradp[1]*hes[1]);}}
  return delta;
}

s_real delta_liq(s_real p, s_real tau, s_real *grad, s_real *hes, int *nit){
  /*----------------------------------------------------------------------------
  Get a good liquid or super critical inital density guess then call
  delta_p_tau() to calculate density. In difficult cases an interval with the
  desired soultion is used with a backeting method to produce a better inital
  guess.  There is a area around the critical point and in the super critical
  region where this is hard to solve so there is some pretty extensive code for
  guessing.

  Since solving this takes some time and this function is called a lot with the
  exact same inputs when using the Pyomo wrapper, this function is memoized.

  Args:
   p: pressure (kPa)
   tau: inverse of reduced pressure Tc/T
   grad: location to return gradient of delta wrt pr and tau or NULL
   hes: location to return hessian (upper triangle) of NULL
   nit: pointer to return number of iterations to, or NULL
  Returns:
   delta: reduced density (accuracy depends on tolarance and function shape)
  ----------------------------------------------------------------------------*/
  s_real val = memoize::get_bin(memoize::delta_liq, p, tau, grad, hes);
  if(!std::isnan(val)) return val;
  bool free_grad = 0, free_hes = 0; // grad and/or hes not provided so allocate
  // Since I'm going to cache results, grad and hes will get calculated
  // whether requested or not.  If they were NULL allocate space.
  if(grad==NULL){grad = new s_real[2]; free_grad = 1;}
  if(hes==NULL){hes = new s_real[3]; free_hes = 1;}
  s_real delta, a=265, b=435;
  s_real delta0, tol=TOL_DELTA_LIQ;
  bool use_rf = 0;
  if(tau <= 1.0 || p >= P_c){
    //over critical T or P
    s_real T = T_c/tau;
    if(T <= 700 && p >= 201.130822*T - 108125.908560
      && p <= 513.615856*T - 308340.302065 && tau <= 1.0 && p >= P_c){
       //Supercritical and below 700K and above 230 kg/m3 and 500 kg/m3
      if(p <= 209.854135*T - 113739.536798){use_rf = 1; a = 200; b = 270;} //to 240
      else if(p <= 218.470983*T - 119297.736552){use_rf = 1; a = 210; b = 280;} //to 250
      else if(p <= 227.005303*T - 124813.248970){use_rf = 1; a = 220; b = 290;} //to 260
      else if(p <= 235.485796*T - 130302.357984){use_rf = 1; a = 230; b = 300;} //to 270
      else if(p <= 243.944935*T - 135784.135585){use_rf = 1; a = 240; b = 310;} //to 280
      else if(p <= 252.418201*T - 141279.863409){use_rf = 1; a = 250; b = 320;} //290
      else if(p <= 260.943770*T - 146812.785189){use_rf = 1; a = 260; b = 330;} //300
      else if(p <= 269.562490*T - 152408.083194){use_rf = 1; a = 270; b = 340;} //310
      else if(p <= 278.317658*T - 158092.719009){use_rf = 1; a = 280; b = 350;} //320
      else if(p <= 280.088999*T - 159242.829624){use_rf = 1; a = 290; b = 352;} //322
      else if(p <= 287.254070*T - 163894.756927){use_rf = 1; a = 292; b = 360;} //330
      else if(p <= 296.416346*T - 169842.139858){use_rf = 1; a = 300; b = 370;} //340
      else if(p <= 315.588522*T - 182278.125877){use_rf = 1; a = 310; b = 390;} //360
      else if(p <= 325.682543*T - 188817.815209){use_rf = 1; a = 330; b = 400;} //370
      else if(p <= 336.179290*T - 195610.653555){use_rf = 1; a = 340; b = 410;} //380
      else if(p <= 347.141083*T - 202694.642666){use_rf = 1; a = 350; b = 420;} //390
      else if(p <= 358.643866*T - 210116.489968){use_rf = 1; a = 360; b = 430;} //400
      else if(p <= 370.770948*T - 217927.339647){use_rf = 1; a = 370; b = 440;} //410
      else if(p <= 383.599932*T - 226173.842356){use_rf = 1; a = 380; b = 450;} //420
      else if(p <= 397.187936*T - 234888.065525){use_rf = 1; a = 390; b = 460;} //430
      else if(p <= 411.562304*T - 244081.188734){use_rf = 1; a = 400; b = 470;} //440
      else if(p <= 426.721613*T - 253744.275496){use_rf = 1; a = 410; b = 480;} //450
      else if(p <= 442.646089*T - 263855.481587){use_rf = 1; a = 420; b = 490;} //460
      else if(p <= 459.311492*T - 274389.563829){use_rf = 1; a = 430; b = 500;} //470
      else if(p <= 476.699579*T - 285324.927152){use_rf = 1; a = 440; b = 505;} //480
      else if(p <= 494.801592*T - 296645.796447){use_rf = 1; a = 450; b = 510;} //490
      else {use_rf = 1; a = 460; b = 515;}}
    else if(p >= -0.030275366391367*T*T + 800.362108600627*T - 468414.629148942){//rho >= 600
      delta0 = 1100.0/rho_c;}
    else if(p >= -0.000693273085341*T*T + 539.413129025933*T - 324914.831221387){//rho >= 500
      delta0 = 600.0/rho_c;}
    else if((p >= 0.456576012809866*T*T - 224.822653977863*T - 23457.7703540548)&& (T <= 700)){ //rho >= 425 and T < 700
      delta0 = 430/rho_c;}
    else if(p <= -2.26154966031622E-06*T*T + 0.467571780470989*T - 4.3044421839477){//rho <= 1
      delta0 = 0.5/rho_c;}
    else if(p <= -0.001481257848455*T*T + 15.4473970626314*T - 2739.92167514421){//rho <= 25
      delta0 = 10.0/rho_c;}
    else if(p <= -0.005800701021532*T*T + 38.396281025698*T - 10729.2461152518){//rho <= 50
      delta0 = 35/rho_c;}
    else if(p <= -0.019380208563602*T*T + 99.0627216060896*T - 38041.4518082429){//rho <= 100
      delta0 = 75/rho_c;}
    else if(p <= -0.039367518187518*T*T + 222.056595942056*T - 105356.343016244){//rho <= 200
      delta0 = 150/rho_c;}
    else if(p <= -0.028017508220982*T*T + 277.53686176142*T - 145798.12140543){//rho <= 275
      delta0 = 250/rho_c;}
    else if(T>=700){delta0=400/rho_c;}
    else if(p <= 0.378699785930901*T*T - 204.557196883595*T - 4116.16210692618){
      //rho about between 275 and 350
      use_rf = 1;
      if(T<= 647.096) a=322;
      if(p <= P_c) b = 322;
      else b = 360;
    }
    else{
      // rho about between 350 and 425
      use_rf = 1;
      if(T<= 647.096) a=340;
    }
  }
  else if(p > 21.5 && tau < T_c/645){
    use_rf = 1;
    a = 322;
    b = 480;
    delta0 = 500/rho_c;
  }
  else delta0=1100/rho_c;
  if(use_rf) delta0 = delta_p_tau_rf(p, tau, a, b, 0); //bracket for better i.g.
  delta = delta_p_tau(p, tau, delta0, tol, nit, grad, hes); //solve
  if(std::isnan(delta) || delta < 1e-12 || delta > 5.0){
    // This is just to avoid evaluation errors.  Want to be able to calucalte
    // vapor properties even when vapor doesn't exist.  In the IDAES Framework
    // these properties may be calculated and multipled by a zero liquid fraction,
    // so it doesn't mater that they are wrong.
    delta = 3.1;
    zero_derivs2(grad, hes);
  }
  memoize::add_bin(memoize::delta_liq, p, tau, delta, grad, hes); //store
  if(free_grad) delete[] grad; // free grad and hes if not allocated by calling
  if(free_hes) delete[] hes;   //   function
  return delta;
}

s_real delta_vap(s_real p, s_real tau, s_real *grad, s_real *hes, int *nit){
  /*----------------------------------------------------------------------------
  Get a good vapor or super critical inital density guess then call
  delta_p_tau() to calculate density. In the supercritical region this just
  calls the liquid function. In the rest of the vapor region the inital guess
  is pretty easy.

  Since solving this takes some time and this function is called a lot with the
  exact same inputs when using the Pyomo wrapper, this function is memoized.

  Args:
   p: pressure (kPa)
   tau: inverse of reduced pressure Tc/T
   grad: location to return gradient of delta wrt pr and tau or NULL
   hes: location to return hessian (upper triangle) of NULL
   nit: pointer to return number of iterations to, or NULL
  Returns:
   delta: reduced density (accuracy depends on tolarance and function shape)
  ----------------------------------------------------------------------------*/
  s_real val = memoize::get_bin(memoize::DV_FUNC, p, tau, grad, hes);
  if(!std::isnan(val)) return val; // return stored result if available
  s_real delta, delta0 = 0.01;
  bool free_grad = 0, free_hes = 0; // grad and/or hes not provided so allocate
  // If supercritical use the liquid calculation, which includes sc region
  if(tau <= 1.0 || p >= P_c) return delta_liq(p, tau, grad, hes);
  // Since I'm going to cache results, grad and hes will get calculated
  // whether requested or not.  If they were NULL allocate space.
  if(grad==NULL){grad = new s_real[2]; free_grad = 1;}
  if(hes==NULL){hes = new s_real[3]; free_hes = 1;}
  if(p < 1e3) delta0 = 0.0001; // Low pressure so extra low density guess
  else delta0 = 0.001; // this guess should work for rest of vapor region
  delta = delta_p_tau(p, tau, delta0, TOL_DELTA_VAP, nit, grad, hes);
  if(std::isnan(delta) || delta < 1e-12 || delta > 5.0){
    // This is just to avoid evaluation errors.  Want to be able to calucalte
    // vapor properties even when vapor doesn't exist.  In the IDAES Framework
    // these properties may be calculated and multipled by a zero vapor fraction,
    // so it doesn't mater that they are wrong.
    delta = 0.001;
    zero_derivs2(grad, hes);
  }
  memoize::add_bin(memoize::DV_FUNC, p, tau, delta, grad, hes); // store result
  if(free_grad) delete[] grad; // free grad and hes if not allocated by calling
  if(free_hes) delete[] hes; //   function
  return delta;
}


/*------------------------------------------------------------------------------
In this section you will find functions for calculating the saturation curve.
Staturation pressure and density as a function of temperature.
*-----------------------------------------------------------------------------*/
s_real sat_delta_liq(s_real tau){ //caculate saturated liquid density from tau
  s_real delta_l, delta_v;
  sat(tau, &delta_l, &delta_v);
  return delta_l;
}

s_real sat_delta_vap(s_real tau){ //caculate saturated vapor density from tau
  s_real delta_l, delta_v;
  sat(tau, &delta_l, &delta_v);
  return delta_v;
}

s_real p_sat_iapws97(s_real tau){ //saturation pressure from tau IAPWS-97 eq.
  //the IAPWS-97 isn't as consitent as IAPWS-95, but this provides a good guess
  s_real T = T_c/tau;
  s_real tt = T + n_psat[8]/(T - n_psat[9]);
  s_real A = tt*tt + n_psat[0]*tt + n_psat[1];
  s_real B = n_psat[2]*tt*tt + n_psat[3]*tt + n_psat[4];
  s_real C = n_psat[5]*tt*tt + n_psat[6]*tt + n_psat[7];
  return 1000*pow(2*C/(-B + pow(B*B - 4*A*C, 0.5)), 4);
}

s_real delta_sat_v_approx(s_real tau){ //approximate saturated vapor density
  // This equation is from the original IAPWS-95 paper
  s_real XX = 1 - 1.0/tau;
  s_real delta = exp(-2.03150240*pow(XX,2.0/6.0)
            - 2.68302940*pow(XX,4.0/6.0)
            - 5.38626492*pow(XX,8.0/6.0)
            - 17.2991605*pow(XX,18.0/6.0)
            - 44.7586581*pow(XX,37.0/6.0)
            - 63.9201063*pow(XX,71.0/6.0));
  return delta_p_tau(p_sat_iapws97(tau), tau, delta);
}

s_real delta_sat_l_approx(s_real tau){ //approximate saturated vapor density
  // This equation is from the original IAPWS-95 paper.
  s_real XX = 1 - 1.0/tau;
  s_real delta = 1.001 + 1.99274064*pow(XX,1.0/3.0)
           + 1.09965342*pow(XX,2.0/3.0)
           - 0.510839303*pow(XX,5.0/3.0)
           - 1.75493479*pow(XX,16.0/3.0)
           - 45.5170352*pow(XX,43.0/3.0)
           - 6.74694450e5*pow(XX,110.0/3.0);
  return delta_p_tau(p_sat_iapws97(tau), tau, delta);
}

inline s_real J(s_real delta, s_real tau){
  // Term from Akasaka method for saturation state
  return delta*(1+delta*phir_delta(delta, tau));
}

inline s_real K(s_real delta, s_real tau){
  // Term from Akasaka method for saturation state
  return delta*phir_delta(delta,tau) + phir(delta, tau) + log(delta);
}

inline s_real J_delta(s_real delta, s_real tau){
  return 1.0 + 2.0*delta*phir_delta(delta, tau) + delta*delta*phir_delta2(delta, tau);
}

inline s_real K_delta(s_real delta, s_real tau){
  return 2.0*phir_delta(delta, tau) + delta*phir_delta2(delta, tau) + 1.0/delta;
}

inline s_real Delta_Aka(s_real delta_l, s_real delta_v, s_real tau){
  return J_delta(delta_v, tau)*K_delta(delta_l, tau) -
         J_delta(delta_l, tau)*K_delta(delta_v, tau);
}

int sat(s_real tau, s_real *delta_l_sol, s_real *delta_v_sol){
    //Get stautated phase densities at tau by Akasaka (2008) method
    s_real delta_l, delta_v, fg;
    s_real gradl[1], hesl[1], gradv[1], hesv[1];
    int n = 0, max_it=MAX_IT_SAT;
    s_real agamma = SAT_GAMMA, tol = TOL_REL_SAT_G;

    if(tau < T_c/647.094){
      agamma = 0.40;
      if(tau - 1 < 1e-13){
        delta_l = 1.0;
        delta_v = 1.0;
        max_it=0;
      }
      else{
        delta_l = 1.02;
        delta_v = 0.98;
        tol *= 1e-2;
        max_it *= 4;
      }
    }
    else{
      // okay so you've decided to solve this thing
      delta_l = delta_sat_l_approx(tau); // guess based on IAPWS-97
      delta_v = delta_sat_v_approx(tau); // guess based on IAPWS-97
    }
    // Since the equilibrium conditions are gl = gv and pl = pv, I am using the
    // the relative differnce in g as a convergence criteria, that is easy to
    // understand.  fg < tol for convergence, fg is calucalted upfront in the
    // off chance that the guess is the solution
    *delta_l_sol = delta_l; // just in case we don't do at least 1 iteration
    *delta_v_sol = delta_v; // just in case we don't do at least 1 iteration
    fg = fabs((g(delta_v, tau) - g(delta_l, tau))/g(delta_l, tau));
    while(n<max_it && fg > tol){
      ++n; // Count iterations
      //calculations deltas at next step (Akasaka (2008))
      *delta_l_sol = delta_l + agamma/Delta_Aka(delta_l, delta_v, tau)*(
             (K(delta_v, tau) - K(delta_l, tau))*J_delta(delta_v,tau) -
             (J(delta_v, tau) - J(delta_l, tau))*K_delta(delta_v,tau));
      *delta_v_sol = delta_v + agamma/Delta_Aka(delta_l, delta_v, tau)*(
             (K(delta_v, tau) - K(delta_l, tau))*J_delta(delta_l,tau) -
             (J(delta_v, tau) - J(delta_l, tau))*K_delta(delta_l,tau));
      delta_v = *delta_v_sol; //step
      delta_l = *delta_l_sol;
      //calculate convergence criterium
      fg = fabs((g(delta_v, tau) - g(delta_l, tau))/g(delta_l, tau));
    }

    //Calculate grad and hes for and memoize

    gradv[0] = LHM/LGM;
    gradl[0] = gradv[0]*LBV/LBL + (LCV - LCL)/LBL;
    hesv[0] = LdHdt(delta_l, delta_v, tau, gradl[0], gradv[0])/LGM
             - LHM/LGM/LGM*LdGdt(delta_l, delta_v, tau, gradl[0], gradv[0]);
    hesl[0] = hesv[0]*LBV*LFL + gradv[0]*(LBVt + LBVd*gradv[0])*LFL
              + gradv[0]*LBV*(LFLt + LFLd*gradl[0]) + (LFLt + LFLd*gradl[0])*(LCV - LCL)
              + LFL*(LCVt - LCLt + LCVd*gradv[0] - LCLd*gradl[0]);

    memoize::add_un(memoize::DL_SAT_FUNC, tau, delta_l, gradl, hesl);
    memoize::add_un(memoize::DV_SAT_FUNC, tau, delta_v, gradv, hesv);
    return n;
}

/*------------------------------------------------------------------------------
Basic Properties Calaculations, 1st and 2nd derivatives
*-----------------------------------------------------------------------------*/

s_real p_with_derivs(s_real delta, s_real tau, s_real *grad, s_real *hes){
  //pressure derivatives are extra important in here, so gonna keep s_real, and
  //convert to s_real in ASL function
  s_real val = memoize::get_bin(memoize::P_FUNC, delta, tau, grad, hes);
  if(!std::isnan(val)) return val;
  bool free_grad = 0, free_hes = 0; // grad and/or hes not provided so allocate
  // Since I'm going to cache results, grad and hes will get calculated
  // whether requested or not.  If they were NULL allocate space.
  if(grad==NULL){grad = new s_real[2]; free_grad = 1;}
  if(hes==NULL){hes = new s_real[3]; free_hes = 1;}
  grad[0] = XA_d + XA_d*XC + XA*XC_d;
  grad[1] = XA_t + XA_t*XC + XA*XC_t;
  hes[0] = XA_dd + XA_dd*XC + 2*XA_d*XC_d + XA*XC_dd;
  hes[1] = XA_dt + XA_dt*XC + XA_d*XC_t + XA_t*XC_d + XA*XC_dt;
  hes[2] = XA_tt + XA_tt*XC + 2*XA_t*XC_t + XA*XC_tt;
  s_real pr = p(delta, tau);
  memoize::add_bin(memoize::P_FUNC, delta, tau, pr, grad, hes);
  if(free_grad) delete[] grad; // free grad and hes if not allocated by calling
  if(free_hes) delete[] hes; //   function
  return pr;
}

s_real u_with_derivs(s_real delta, s_real tau, s_real *grad, s_real *hes){
  if(grad!=NULL){
    grad[0] = XB*XH_d;
    grad[1] = XB_t*XH + XB*XH_t;
    if(hes!=NULL){
      hes[0] = XB*XH_dd;
      hes[1] = XB_t*XH_d + XB*XH_dt;
      hes[2] = XB_tt*XH + 2*XB_t*XH_t + XB*XH_tt;}}
  return u(delta, tau);
}

s_real s_with_derivs(s_real delta, s_real tau, s_real *grad, s_real *hes){
  if(grad!=NULL){
    grad[0] = R*(XH_d - phi0_delta(delta) - phir_delta(delta, tau));
    grad[1] = R*(XH_t - phi0_tau(tau) - phir_tau(delta, tau));
    if(hes!=NULL){
      hes[0] = R*(XH_dd - phi0_delta2(delta) - phir_delta2(delta, tau));
      hes[1] = R*(XH_dt - phir_delta_tau(delta, tau));
      hes[2] = R*(XH_tt - phi0_tau2(tau) - phir_tau2(delta, tau));}}
  return s(delta, tau);
}

s_real h_with_derivs(s_real delta, s_real tau, s_real *grad, s_real *hes){
  if(grad!=NULL){
    grad[0] = XB*(XH_d + XC_d);
    grad[1] = XB_t*(1 + XH + XC) + XB*(XH_t + XC_t);
    if(hes!=NULL){
      hes[0] = XB*(XH_dd + XC_dd);
      hes[1] = XB_t*(XH_d + XC_d) + XB*(XH_dt + XC_dt);
      hes[2] = XB_tt*(1 + XH + XC) + 2*XB_t*(XH_t + XC_t) + XB*(XH_tt + XC_tt);}}
  return h(delta, tau);
}

s_real hvpt_with_derivs(s_real pr, s_real tau, s_real *grad, s_real *hes){
  s_real gradh[2], hesh[3]; // derivatives of h(delta, tau)
  s_real gradd[2], hesd[3]; // derivatives of delta(p, tau)
  s_real delta, h;
  delta = delta_vap(pr, tau, gradd, hesd);
  h = h_with_derivs(delta, tau, gradh, hesh);
  if(grad!=NULL){
    grad[0] = gradh[0]*gradd[0];
    grad[1] = gradh[1] + gradh[0]*gradd[1];
    if(hes!=NULL){
      hes[0] = hesh[0]*gradd[0]*gradd[0] + gradh[0]*hesd[0];
      hes[1] = hesh[1]*gradd[0] + hesh[0]*gradd[0]*gradd[1] + gradh[0]*hesd[1];
      hes[2] = hesh[2] + 2*hesh[1]*gradd[1] + hesh[0]*gradd[1]*gradd[1] + gradh[0]*hesd[2];}}
  return h;
}

s_real hlpt_with_derivs(s_real pr, s_real tau, s_real *grad, s_real *hes){
  s_real gradh[2], hesh[3]; // derivatives of h(delta, tau)
  s_real gradd[2], hesd[3]; // derivatives of delta(p, tau)
  s_real delta, h;
  delta = delta_liq(pr, tau, gradd, hesd);
  h = h_with_derivs(delta, tau, gradh, hesh);
  if(grad!=NULL){
    grad[0] = gradh[0]*gradd[0];
    grad[1] = gradh[1] + gradh[0]*gradd[1];
    if(hes!=NULL){
      hes[0] = hesh[0]*gradd[0]*gradd[0] + gradh[0]*hesd[0];
      hes[1] = hesh[1]*gradd[0] + hesh[0]*gradd[0]*gradd[1] + gradh[0]*hesd[1];
      hes[2] = hesh[2] + 2*hesh[1]*gradd[1] + hesh[0]*gradd[1]*gradd[1] + gradh[0]*hesd[2];}}
  return h;
}

s_real vf_with_derivs(s_real ht, s_real pr, s_real *grad, s_real *hes){
    s_real tau, gradt[1], hest[1], gradhv[2], heshv[3], gradhl[2], heshl[3];
    tau = sat_tau_with_derivs(pr, gradt, hest);
    s_real hv = hvpt_with_derivs(pr, tau, gradhv, heshv);
    s_real hl = hlpt_with_derivs(pr, tau, gradhl, heshl);

    if(pr >= P_c){zero_derivs2(grad, hes); return 0.0;} //Classify supercritical as liquid
    else if(hl > ht){ zero_derivs2(grad, hes); return 0.0;}
    else if(hv < ht){ zero_derivs2(grad, hes); return 1.0;}

    if(grad != NULL){
        s_real dhvdp = gradhv[0] + gradhv[1]*gradt[0];
        s_real dhldp = gradhl[0] + gradhl[1]*gradt[0];
        grad[0] = 1.0/(hv - hl);
        grad[1] = -dhldp/(hv - hl) - (ht - hl)/(hv-hl)/(hv-hl)*(dhvdp - dhldp);
        if(hes != NULL){
          s_real d2hvdp2 = heshv[0] + 2*heshv[1]*gradt[0] + heshv[2]*gradt[0]*gradt[0] + gradhv[1]*hest[0];
          s_real d2hldp2 = heshl[0] + 2*heshl[1]*gradt[0] + heshl[2]*gradt[0]*gradt[0] + gradhl[1]*hest[0];
          hes[0] = 0;
          hes[1] = -1.0/(hv-hl)/(hv-hl)*(dhvdp - dhldp);
          hes[2] = -d2hldp2/(hv-hl) + 2*dhldp/(hv-hl)/(hv-hl)*(dhvdp - dhldp) +
                    2*(ht-hl)/(hv-hl)/(hv-hl)/(hv-hl)*(dhvdp - dhldp)*(dhvdp - dhldp) -
                    (ht-hl)/(hv-hl)/(hv-hl)*(d2hvdp2 - d2hldp2);
        }
    }
    return (ht - hl)/(hv - hl);
}

s_real tau_with_derivs(s_real ht, s_real pr, s_real *grad, s_real *hes){
    s_real tau_sat;
    tau_sat = sat_tau_with_derivs(pr, NULL, NULL);
    s_real hv = hvpt_with_derivs(pr, tau_sat, NULL, NULL);
    s_real hl = hlpt_with_derivs(pr, tau_sat, NULL, NULL);
    s_real fun, tau, gradh[2], hesh[3], tol = 1e-11;
    int it = 0, max_it = 20;
    if(hl > ht){
      tau = tau_sat + 0.2;
      fun = hlpt_with_derivs(pr, tau, gradh, hesh) - ht;
      while(fabs(fun) > tol && it < max_it){
        tau = tau - fun*gradh[1]/(gradh[1]*gradh[1] - 0.5*fun*hesh[2]);
        fun = hlpt_with_derivs(pr, tau, gradh, hesh) - ht;
        ++it;
      }
    }
    else if(hv < ht){
      tau = tau_sat - 0.1*(ht - hv)/1000;
      fun = hvpt_with_derivs(pr, tau, gradh, hesh) - ht;
      while(fabs(fun) > tol && it < max_it){
        tau = tau - fun*gradh[1]/(gradh[1]*gradh[1] - 0.5*fun*hesh[2]);
        fun = hvpt_with_derivs(pr, tau, gradh, hesh) - ht;
        ++it;
      }
    }
    else{
      zero_derivs2(grad, hes);
      return tau_sat;
    }
    if(tau < 0.0 || tau > TAU_HIGH){
        std::cerr << "IAPWS LOW T CLIP WARNING: h = " << ht << " P= " << pr << " tau = " << tau << "\n";
        tau = TAU_HIGH;
        fun = hvpt_with_derivs(pr, tau, gradh, hesh) - ht;
    }
    else if(tau < TAU_LOW){
        std::cerr << "IAPWS HIGH T CLIP WARNING: h = " << ht << " P= " << pr << " tau = " << tau << "\n";
        tau = TAU_LOW;
        fun = hvpt_with_derivs(pr, tau, gradh, hesh) - ht;
    }
    if(grad != NULL){
        grad[0] = 1.0/gradh[1];
        grad[1] = -grad[0]*gradh[0];
        if(hes != NULL){
          hes[0] = -grad[0]*grad[0]*grad[0]*hesh[2];
          hes[1] = -grad[0]*grad[0]*(hesh[1] + hesh[2]*grad[1]);
          hes[2] = -hes[1]*gradh[0] - grad[0]*(hesh[0] + hesh[1]*grad[1]);
        }
    }
    return tau;
}

s_real g_with_derivs(s_real delta, s_real tau, s_real *grad, s_real *hes){
  s_real h_grad[2], h_hes[3], s_grad[2], s_hes[3], s_val;
  s_real *h_grad_ptr=NULL, *h_hes_ptr=NULL, *s_grad_ptr=NULL, *s_hes_ptr=NULL;

  if(grad!=NULL){
    h_grad_ptr=h_grad;
    s_grad_ptr=s_grad;
    if(hes!=NULL){
      h_hes_ptr=h_hes;
      s_hes_ptr=s_hes;}}

  if(grad!=NULL){
    h_with_derivs(delta, tau, h_grad_ptr, h_hes_ptr);
    s_val = s_with_derivs(delta, tau, s_grad_ptr, s_hes_ptr);
    grad[0] = h_grad[0] - T_c/tau*s_grad[0];
    grad[1] = h_grad[1] + T_c/tau/tau*s_val - T_c/tau*s_grad[1];
    if(hes!=NULL){
      hes[0] = h_hes[0] - T_c/tau*s_hes[0];
      hes[1] = h_hes[1] + T_c/tau/tau*s_grad[0] - T_c/tau*s_hes[1];
      hes[2] = h_hes[2] - 2.0*T_c/tau/tau/tau*s_val + 2.0*T_c/tau/tau*s_grad[1]
        - T_c/tau*s_hes[2];}}
  return g(delta, tau);
}

s_real f_with_derivs(s_real delta, s_real tau, s_real *grad, s_real *hes){
  if(grad!=NULL){
    grad[0] = XB*(phi0_delta(delta) + phir_delta(delta,tau));
    grad[1] = XB_t*(phi0(delta, tau) + phir(delta,tau)) +
      XB*(phi0_tau(tau) + phir_tau(delta,tau));
    if(hes!=NULL){
      hes[0] = XB*(phi0_delta2(delta) + phir_delta2(delta,tau));
      hes[1] = XB_t*(phi0_delta(delta) + phir_delta(delta,tau)) +
        XB*phir_delta_tau(delta,tau);
      hes[2] = XB_tt*(phi0(delta, tau) + phir(delta,tau)) +
        2.0*XB_t*(phi0_tau(tau) + phir_tau(delta,tau)) +
        XB*(phi0_tau2(tau) + phir_tau2(delta,tau));}}
  return f(delta, tau);
}

s_real cv_with_derivs(s_real delta, s_real tau, s_real *grad, s_real *hes){
  if(grad!=NULL){
    grad[0] = -tau*tau*R*(phi0_delta_tau2 + phir_delta_tau2(delta, tau));
    grad[1] = -2.0*tau*R*(phi0_tau2(tau) + phir_tau2(delta, tau)) -
      tau*tau*R*(phi0_tau3(tau) + phir_tau3(delta, tau));
    if(hes!=NULL){
      hes[0] = -tau*tau*R*(phi0_delta2_tau2 + phir_delta2_tau2(delta, tau));
      hes[1] = -2.0*tau*R*(phi0_delta_tau2 + phir_delta_tau2(delta, tau)) -
        tau*tau*R*(phi0_delta_tau3 + phir_delta_tau3(delta, tau));
      hes[2] = -2.0*R*(phi0_tau2(tau) + phir_tau2(delta, tau)) -
        4.0*tau*R*(phi0_tau3(tau) + phir_tau3(delta, tau)) -
        tau*tau*R*(phi0_tau4(tau) + phir_tau4(delta, tau));}}
  return cv(delta, tau);
}

s_real cp_with_derivs(s_real delta, s_real tau, s_real *grad, s_real *hes){
  s_real cv_grad[2], cv_hes[3];
  s_real *cv_grad_ptr=NULL, *cv_hes_ptr=NULL;
  if(grad!=NULL){
    cv_grad_ptr = cv_grad;
    if(hes!=NULL){
      cv_hes_ptr = cv_hes;}}
  if(grad!=NULL){
    cv_with_derivs(delta, tau, cv_grad_ptr, cv_hes_ptr);
    grad[0] = cv_grad[0] + R*XD_d/XE - R*XD*XE_d/XE/XE;
    grad[1] = cv_grad[1] + R*XD_t/XE - R*XD*XE_t/XE/XE;
    if(hes!=NULL){
      hes[0] = cv_hes[0] + R*XD_dd/XE - 2.0*R*XD_d*XE_d/XE/XE -
        R*XD*XE_dd/XE/XE + 2.0*R*XD*XE_d*XE_d/XE/XE/XE;
      hes[1] = cv_hes[1] + R*XD_dt/XE - R*XD_d*XE_t/XE/XE - R*XD_t*XE_d/XE/XE -
        R*XD*XE_dt/XE/XE + 2.0*R*XD*XE_d*XE_t/XE/XE/XE;
      hes[2] = cv_hes[2] + R*XD_tt/XE - 2.0*R*XD_t*XE_t/XE/XE -
        R*XD*XE_tt/XE/XE + 2.0*R*XD*XE_t*XE_t/XE/XE/XE;}}
  return cp(delta, tau);
}

s_real w_with_derivs(s_real delta, s_real tau, s_real *grad, s_real *hes){
  //TODO<jce> hey you forgot one
  if(grad!=NULL){
    grad[0] = 0;
    grad[1] = 0;
    if(hes!=NULL){
      hes[0] = 0;
      hes[1] = 0;
      hes[2] = 0;}}
  return w(delta, tau);
}

s_real sat_delta_liq_with_derivs(s_real tau, s_real *grad, s_real *hes){
  // For efficency I'm requiring memoization here.  Delta_v and delta_l are
  // calculated together so should just remember both.  This will break if you
  // turn off memoiztion.
  s_real delta_l, delta_v;
  s_real val = memoize::get_un(memoize::DL_SAT_FUNC, tau, grad, hes);
  if(!std::isnan(val)) return val;
  sat(tau, &delta_l, &delta_v);
  // Should be there now
  val = memoize::get_un(memoize::DL_SAT_FUNC, tau, grad, hes);
  if(!std::isnan(val)) return val;
  return (s_real)NAN;
}

s_real sat_delta_vap_with_derivs(s_real tau, s_real *grad, s_real *hes){
  // For efficency I'm requiring memoization here.  Delta_v and delta_l are
  // calculated together so should just remember both.  This will break if you
  // turn off memoiztion.
  s_real delta_l, delta_v;
  s_real val = memoize::get_un(memoize::DV_SAT_FUNC, tau, grad, hes);
  if(!std::isnan(val)) return val;
  sat(tau, &delta_l, &delta_v);
  // Should be there now
  val = memoize::get_un(memoize::DV_SAT_FUNC, tau, grad, hes);
  if(!std::isnan(val)) return val;
  return (s_real)NAN;
}

s_real sat_p_with_derivs(s_real tau, s_real *grad, s_real *hes, bool limit){
  //Before getting into the real calculation, check if outside the allowed range
  //of 240K to T_c, the low end of this range doen't mean anything.  Its too cold
  //to expect liquid, but the calucations hold up there so for numerical reasons
  //I'll allow it.  Above the critical temperture there is only a single phase
  if(tau > 647.096/240.0 && limit){ // below 270 K
    if(grad != NULL) grad[0] = 0;
    if(hes != NULL) hes[0] = 0;
    return +3.76195238e-02;
  }
  else if(tau < 1.0 && limit){ // above critical T
    if(grad != NULL) grad[0] = 0;
    if(hes != NULL) hes[0] = 0;
    return P_c;
  }

  //Now check if there is a stored value
  s_real val = memoize::get_un(memoize::P_SAT_FUNC, tau, grad, hes);
  if(!std::isnan(val)) return val;
  // grad and/or hes not provided so allocate for memoization
  bool free_grad = 0, free_hes = 0;
  if(grad==NULL){grad = new s_real[1]; free_grad = 1;}
  if(hes==NULL){hes = new s_real[1]; free_hes = 1;}
  s_real grad_delta[1], hes_delta[1], grad_p[2], hes_p[3];
  s_real delta = sat_delta_vap_with_derivs(tau, grad_delta, hes_delta);
  s_real Psat = p_with_derivs(delta, tau, grad_p, hes_p);
  grad[0] = grad_p[1] + grad_p[0]*grad_delta[0];
  hes[0] = hes_p[2] + 2*hes_p[1]*grad_delta[0] + grad_p[0]*hes_delta[0]
           + hes_p[0]*grad_delta[0]*grad_delta[0];
  memoize::add_un(memoize::P_SAT_FUNC, tau, Psat, grad, hes);
  if(free_grad) delete[] grad; // free grad and hes if not allocated by calling
  if(free_hes) delete[] hes;   // function
  return Psat;
}

s_real sat_tau_with_derivs(s_real pr, s_real *grad, s_real *hes, int *nit){
  //Before getting into the real calculation, check if outside the allowed range
  //of 240K to T_c, the low end of this range doen't mean anything.  Its too cold
  //to expect liquid, but the calucations hold up there so for numerical reasons
  //I'll allow it.  Above the critical temperture there is only a single phase
  if(pr < +3.76195238e-02){ // below 270 K
    if(grad != NULL) grad[0] = 0;
    if(hes != NULL) hes[0] = 0;
    return 647.096/240.0;
  }
  else if(pr > P_c){ // above critical T
    if(grad != NULL) grad[0] = 0;
    if(hes != NULL) hes[0] = 0;
    return 1.0;
  }

  //Now check if there is a stored value
  s_real val = memoize::get_un(memoize::TAU_SAT_FUNC, pr, grad, hes);
  if(!std::isnan(val)) return val;
  // grad and/or hes not provided so allocate for memoization
  bool free_grad = 0, free_hes = 0;
  if(grad==NULL){grad = new s_real[1]; free_grad = 1;}
  if(hes==NULL){hes = new s_real[1]; free_hes = 1;}
  s_real tau = 1.5, fun, gradp[1], hesp[1], tol=TOL_SAT_TAU;
  if(P_c - pr < 1e-3){
    pr = P_c - 1e-3;
  }
  int it = 0; // iteration count
  fun = sat_p_with_derivs(tau, gradp, hesp, 0) - pr;
  while(fabs(fun) > tol && it < MAX_IT_SAT_TAU){
    tau = tau - fun*gradp[0]/(gradp[0]*gradp[0] - 0.5*fun*hesp[0]);
    fun = sat_p_with_derivs(tau, gradp, hesp, 0) - pr;
    ++it;
  }
  grad[0] = 1.0/gradp[0];
  hes[0] = -grad[0]*grad[0]*grad[0]*hesp[0];
  memoize::add_un(memoize::TAU_SAT_FUNC, pr, tau, grad, hes);
  if(free_grad) delete[] grad; // free grad and hes if not allocated by calling
  if(free_hes) delete[] hes;   // function
  if(nit != NULL) *nit = it;
  return tau;
}
