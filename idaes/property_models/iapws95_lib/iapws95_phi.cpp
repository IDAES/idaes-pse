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
 File: iapws95_phi.cpp
------------------------------------------------------------------------------*/

#include<stdio.h>
#include<cmath>
#include<utility>
#include<iostream>

#include"iapws95_phi.h"
#include"iapws95_deriv_parts.h"
#include"iapws95_memo.h"


/*------------------------------------------------------------------------------
  phi0 and derivatives
------------------------------------------------------------------------------*/

s_real phi0(s_real delta, s_real tau){
  s_real sum4;
  //Calculate sums
  sum4 = 0;
  for(int i=4; i<9; ++i){
    sum4 += n0[i1(i)]*s_log(1 - s_exp(-gamma0[i4(i)]*tau));
  }
  //return result
  return s_log(delta) + n0[i1(1)] + n0[i1(2)]*tau + n0[i1(3)]*s_log(tau) + sum4;
}

s_real phi0_delta(s_real delta){
  return 1.0/delta;
}

s_real phi0_delta2(s_real delta){
  return -1.0/delta/delta;
}

s_real phi0_tau(s_real tau){
  s_real sum4;
  sum4 = 0;
  for(int i=4; i<9; ++i){
    sum4 += n0[i1(i)]*gamma0[i4(i)]*(1.0/(1 - s_exp(-gamma0[i4(i)]*tau)) - 1);
  }
  return n0[i1(2)] + n0[i1(3)]/tau + sum4;
}

s_real phi0_tau2(s_real tau){
  s_real sum4;
  sum4 = 0;
  for(int i=4; i<9; ++i){
    sum4 += n0[i1(i)]*gamma0[i4(i)]*gamma0[i4(i)]*s_exp(-gamma0[i4(i)]*tau)
            *s_pow(1 - s_exp(-gamma0[i4(i)]*tau), -2);
  }
  return -n0[i1(3)]/tau/tau - sum4;
}

s_real phi0_delta3(s_real delta){
  return 2.0/delta/delta/delta;
}

s_real phi0_tau3(s_real tau){
  s_real E, Et, gam, sum4;
  sum4 = 0;
  for(int i=4; i<9; ++i){
    gam = gamma0[i4(i)];
    E = s_exp(-gam*tau);
    Et = -gam*E;
    sum4 += n0[i1(i)]*gam*gam*(Et/(1-E)/(1-E) + 2.0*E*Et/(1-E)/(1-E)/(1-E));
  }
  return 2.0*n0[i1(3)]/tau/tau/tau - sum4;
}

s_real phi0_delta4(s_real delta){
  return -6.0/delta/delta/delta/delta;
}

s_real phi0_tau4(s_real tau){
  s_real E, Et, Ett, gam, sum4;
  sum4 = 0;
  for(int i=4; i<9; ++i){
    gam = gamma0[i4(i)];
    E = s_exp(-gam*tau);
    Et = -gam*E;
    Ett = gam*gam*E;
    sum4 += n0[i1(i)]*gam*gam*(
      Ett/(1-E)/(1-E) + 2.0*Et/(1-E)/(1-E)/(1-E)*Et + 2.0*Et/(1-E)/(1-E)/(1-E)*Et
      + 2.0*E/(1-E)/(1-E)/(1-E)*Ett + 6.0*E*Et*Et/(1-E)/(1-E)/(1-E)/(1-E));
  }
  return -6.0*n0[i1(3)]/tau/tau/tau/tau - sum4;
}

/*------------------------------------------------------------------------------
  phir and derivatives
------------------------------------------------------------------------------*/
s_real phir(s_real delta, s_real tau){
  //Check if stored and return stored value if so
  double val = memoize::get_bin0(memoize::phir, delta, tau);
  if(!std::isnan(val)) return val;
  s_real sum = 0;
  //Calculate sums
  for(int i=1; i<8; ++i){
    sum += n[i1(i)]*s_pow(delta, d[i1(i)])*s_pow(tau, t[i1(i)]);
  }
  for(int i=8; i<52; ++i){
    sum += n[i1(i)]*s_pow(delta, d[i1(i)])*s_pow(tau, t[i1(i)])
            *s_exp(-s_pow(delta,c[i8(i)]));
  }
  for(int i=52; i<55; ++i){
    sum += n[i1(i)]*s_pow(delta, d[i1(i)])*s_pow(tau, t[i1(i)])
             *s_exp(-alpha[i52(i)]*s_pow(delta - eps[i52(i)],2)
                   -beta[i52(i)]*s_pow(tau - gamma_[i52(i)],2));
  }
  //Store
  memoize::add_bin0(memoize::phir, delta, tau, sum);
  return sum;
}

s_real phir_delta(s_real delta, s_real tau){
  //Check if stored and return stored value if so
  double val = memoize::get_bin0(memoize::phir_delta, delta, tau);
  if(!std::isnan(val)) return val;
  s_real sum=0;
  //Calculate sums
  for(int i=1; i<8; ++i){
    sum += n[i1(i)]*d[i1(i)]*s_pow(delta, d[i1(i)]-1)*s_pow(tau, t[i1(i)]);
  }
  for(int i=8; i<52; ++i){
    sum += n[i1(i)]*s_pow(delta, d[i1(i)]-1)*s_pow(tau, t[i1(i)])
            *s_exp(-s_pow(delta,c[i8(i)]))*(d[i1(i)]-c[i8(i)]*s_pow(delta, c[i8(i)]));
  }
  for(int i=52; i<55; ++i){
    sum += n[i1(i)]*s_pow(delta, d[i1(i)])*s_pow(tau, t[i1(i)])
             *s_exp(-alpha[i52(i)]*s_pow(delta - eps[i52(i)],2)
                   -beta[i52(i)]*s_pow(tau - gamma_[i52(i)],2))*
             (d[i1(i)]/delta - 2*alpha[i52(i)]*(delta-eps[i52(i)]));
  }
  memoize::add_bin0(memoize::phir_delta, delta, tau, sum);
  return sum;
}

s_real phir_delta2(s_real delta, s_real tau){
  //Check if stored and return stored value if so
  double val = memoize::get_bin0(memoize::phir_delta2, delta, tau);
  if(!std::isnan(val)) return val;
  s_real sum = 0;
  //Calculate sums
  for(int i=1; i<8; ++i){
    sum += n[i1(i)]*d[i1(i)]*(d[i1(i)]-1)*s_pow(delta, d[i1(i)]-2)*s_pow(tau, t[i1(i)]);
  }
  for(int i=8; i<52; ++i){
    sum += n[i1(i)]*s_exp(-s_pow(delta,c[i8(i)]))*s_pow(delta, d[i1(i)]-2)*s_pow(tau, t[i1(i)])
            *((d[i1(i)]-c[i8(i)]*s_pow(delta, c[i8(i)]))
            *(d[i1(i)]-1-c[i8(i)]*s_pow(delta,c[i8(i)])) -
            c[i8(i)]*c[i8(i)]*s_pow(delta,c[i8(i)]));
  }
  for(int i=52; i<55; ++i){
    sum += n[i1(i)]*s_pow(tau, t[i1(i)])
             *s_exp(-alpha[i52(i)]*(delta - eps[i52(i)])*(delta - eps[i52(i)])
                   -beta[i52(i)]*(tau - gamma_[i52(i)])*(tau - gamma_[i52(i)]))*
             (-2.0*alpha[i52(i)]*s_pow(delta, d[i1(i)]) +
              4.0*alpha[i52(i)]*alpha[i52(i)]*s_pow(delta, d[i1(i)])*
                  (delta-eps[i52(i)])*(delta-eps[i52(i)]) -
              4.0*alpha[i52(i)]*d[i1(i)]*s_pow(delta, d[i1(i)]-1)*(delta-eps[i52(i)]) +
              d[i1(i)]*(d[i1(i)]-1)*s_pow(delta, d[i1(i)]-2));
  }
  memoize::add_bin0(memoize::phir_delta2, delta, tau, sum);
  return sum;
}

s_real phir_tau(s_real delta, s_real tau){
  //Check if stored and return stored value if so
  double val = memoize::get_bin0(memoize::phir_tau, delta, tau);
  if(!std::isnan(val)) return val;
  s_real sum = 0;
  //Calculate sums
  for(int i=1; i<8; ++i){
    sum += n[i1(i)]*t[i1(i)]*s_pow(delta, d[i1(i)])*s_pow(tau, t[i1(i)]-1);
  }
  for(int i=8; i<52; ++i){
    sum += n[i1(i)]*t[i1(i)]*s_pow(delta, d[i1(i)])*s_pow(tau, t[i1(i)]-1)*
            s_exp(-s_pow(delta,c[i8(i)]));
  }
  for(int i=52; i<55; ++i){
    sum += n[i1(i)]*s_pow(delta, d[i1(i)])*s_pow(tau, t[i1(i)])
             *s_exp(-alpha[i52(i)]*s_pow(delta - eps[i52(i)],2)
                   -beta[i52(i)]*s_pow(tau - gamma_[i52(i)],2))
             *(t[i1(i)]/tau - 2*beta[i52(i)]*(tau-gamma_[i52(i)]));
  }
  memoize::add_bin0(memoize::phir_tau, delta, tau, sum);
  return sum;
}

s_real phir_delta_tau(s_real delta, s_real tau){
  //Check if stored and return stored value if so
  double val = memoize::get_bin0(memoize::phir_delta_tau, delta, tau);
  if(!std::isnan(val)) return val;
  s_real sum = 0;
  //Calculate sums
  for(int i=1; i<8; ++i){
    sum += n[i1(i)]*t[i1(i)]*d[i1(i)]*s_pow(delta, d[i1(i)]-1)*s_pow(tau, t[i1(i)]-1);
  }
  for(int i=8; i<52; ++i){
    sum += n[i1(i)]*t[i1(i)]*s_pow(delta, d[i1(i)]-1)*s_pow(tau, t[i1(i)]-1)*
            s_exp(-s_pow(delta,c[i8(i)]))*(d[i1(i)] - c[i8(i)]*s_pow(delta,c[i8(i)]));
  }
  for(int i=52; i<55; ++i){
    sum += n[i1(i)]*s_pow(delta, d[i1(i)])*s_pow(tau, t[i1(i)])
             *s_exp(-alpha[i52(i)]*s_pow(delta - eps[i52(i)],2)
                   -beta[i52(i)]*s_pow(tau - gamma_[i52(i)],2))
             *(d[i1(i)]/delta - 2*alpha[i52(i)]*(delta-eps[i52(i)]))
             *(t[i1(i)]/tau - 2*beta[i52(i)]*(tau-gamma_[i52(i)]));
  }
  memoize::add_bin0(memoize::phir_delta_tau, delta, tau, sum);
  return sum;
}

s_real phir_tau2(s_real delta, s_real tau){
  //Check if stored and return stored value if so
  double val = memoize::get_bin0(memoize::phir_tau2, delta, tau);
  if(!std::isnan(val)) return val;
  s_real sum = 0;
  //Calculate sums
  for(int i=1; i<8; ++i){
    sum += n[i1(i)]*t[i1(i)]*(t[i1(i)]-1)*s_pow(delta, d[i1(i)])*s_pow(tau, t[i1(i)]-2);
  }
  for(int i=8; i<52; ++i){
    sum += n[i1(i)]*t[i1(i)]*(t[i1(i)]-1)*s_pow(delta, d[i1(i)])*s_pow(tau, t[i1(i)]-2)*
            s_exp(-s_pow(delta,c[i8(i)]));
  }
  for(int i=52; i<55; ++i){
    sum += n[i1(i)]*s_pow(delta, d[i1(i)])*s_pow(tau, t[i1(i)])
             *s_exp(-alpha[i52(i)]*s_pow(delta - eps[i52(i)],2)
                   -beta[i52(i)]*s_pow(tau - gamma_[i52(i)],2))
             *((t[i1(i)]/tau - 2*beta[i52(i)]*(tau-gamma_[i52(i)]))
             *(t[i1(i)]/tau - 2*beta[i52(i)]*(tau-gamma_[i52(i)]))
             -t[i1(i)]/tau/tau - 2*beta[i52(i)]);
  }
  memoize::add_bin0(memoize::phir_tau2, delta, tau, sum);
  return sum;
}

s_real phir_delta3(s_real delta, s_real tau){
  double val = memoize::get_bin0(memoize::phir_delta3, delta, tau);
  if(!std::isnan(val)) return val;
  s_real sum = 0;
  //Calculate sums
  for(int i=1; i<8; ++i){
    sum += n[i1(i)]*d[i1(i)]*(d[i1(i)]-1)*(d[i1(i)]-2)*
      s_pow(delta, d[i1(i)]-3)*s_pow(tau, t[i1(i)]);
  }
  for(int i=8; i<52; ++i){
    sum += n[i1(i)]*s_pow(tau, t[i1(i)])*(Ea_d*Eg + Ea*Eg_d);
  }
  for(int i=52; i<55; ++i){
    sum += n[i1(i)]*s_pow(tau, t[i1(i)])*(Fa_d*Fm + Fa*Fm_d);
  }
  memoize::add_bin0(memoize::phir_delta3, delta, tau, sum);
  return sum;
}

s_real phir_delta4(s_real delta, s_real tau){
  double val = memoize::get_bin0(memoize::phir_delta4, delta, tau);
  if(!std::isnan(val)) return val;
  s_real sum = 0;
  //Calculate sums
  for(int i=1; i<8; ++i){
    sum += n[i1(i)]*d[i1(i)]*(d[i1(i)]-1)*(d[i1(i)]-2)*(d[i1(i)]-3)*
      s_pow(delta, d[i1(i)]-4)*s_pow(tau, t[i1(i)]);
  }
  for(int i=8; i<52; ++i){
    sum += n[i1(i)]*s_pow(tau, t[i1(i)])*(Ea_dd*Eg + 2.0*Ea_d*Eg_d + Ea*Eg_dd);
  }
  for(int i=52; i<55; ++i){
    sum += n[i1(i)]*s_pow(tau, t[i1(i)])*(Fa_dd*Fm + 2.0*Fa_d*Fm_d + Fa*Fm_dd);
  }
  memoize::add_bin0(memoize::phir_delta4, delta, tau, sum);
  return sum;
}

s_real phir_delta2_tau(s_real delta, s_real tau){
  double val = memoize::get_bin0(memoize::phir_delta2_tau, delta, tau);
  if(!std::isnan(val)) return val;
  s_real sum = 0;
  //Calculate sums
  for(int i=1; i<8; ++i){
    sum += n[i1(i)]*d[i1(i)]*t[i1(i)]*(d[i1(i)]-1)*
      s_pow(delta, d[i1(i)]-2)*s_pow(tau, t[i1(i)]-1);
  }
  for(int i=8; i<52; ++i){
    sum += n[i1(i)]*t[i1(i)]*s_pow(tau, t[i1(i)] - 1)*Ea*Ed*(Ef - Ee);
  }
  for(int i=52; i<55; ++i){
    sum +=
      n[i1(i)]*(t[i1(i)]*s_pow(tau, t[i1(i)]-1)*Fa*Fm + s_pow(tau, t[i1(i)])*Fa_t*Fm);
  }
  memoize::add_bin0(memoize::phir_delta2_tau, delta, tau, sum);
  return sum;
}

s_real phir_delta2_tau2(s_real delta, s_real tau){
  double val = memoize::get_bin0(memoize::phir_delta2_tau2, delta, tau);
  if(!std::isnan(val)) return val;
  s_real sum = 0;
  //Calculate sums
  for(int i=1; i<8; ++i){
    sum += n[i1(i)]*d[i1(i)]*t[i1(i)]*(d[i1(i)]-1)*(t[i1(i)]-1)*
      s_pow(delta, d[i1(i)]-2)*s_pow(tau, t[i1(i)]-2);
  }
  for(int i=8; i<52; ++i){
    sum += n[i1(i)]*t[i1(i)]*(t[i1(i)]-1)*s_pow(tau, t[i1(i)]-2)*Ea*Ed*(Ef - Ee);
  }
  for(int i=52; i<55; ++i){
    sum += n[i1(i)]*(t[i1(i)]*(t[i1(i)]-1)*s_pow(tau, t[i1(i)]-2)*Fa*Fm +
      2*t[i1(i)]*s_pow(tau, t[i1(i)]-1)*Fa_t*Fm + s_pow(tau, t[i1(i)])*Fa_tt*Fm);
  }
  memoize::add_bin0(memoize::phir_delta2_tau2, delta, tau, sum);
  return sum;
}

s_real phir_delta_tau2(s_real delta, s_real tau){
  double val = memoize::get_bin0(memoize::phir_delta_tau2, delta, tau);
  if(!std::isnan(val)) return val;
  s_real sum = 0;
  //Calculate sums
  for(int i=1; i<8; ++i){
    sum += n[i1(i)]*d[i1(i)]*t[i1(i)]*(t[i1(i)]-1)*
      s_pow(delta, d[i1(i)]-1)*s_pow(tau, t[i1(i)]-2);
  }
  for(int i=8; i<52; ++i){
    sum += n[i1(i)]*t[i1(i)]*(t[i1(i)]-1)*s_pow(delta, d[i1(i)]-1)*
      s_pow(tau, t[i1(i)]-2)*Eb*s_exp(-s_pow(delta, c[i8(i)]));
  }
  for(int i=52; i<55; ++i){
    sum += n[i1(i)]*(t[i1(i)]*s_pow(delta, d[i1(i)])*s_pow(tau, t[i1(i)] -1)*
      Fa*Ga*Gb + s_pow(delta, d[i1(i)])*s_pow(tau, t[i1(i)])*
      Fa_t*Ga*Gb + s_pow(delta, d[i1(i)])*s_pow(tau, t[i1(i)])*Fa*Ga*Gb_t);
  }
  memoize::add_bin0(memoize::phir_delta_tau2, delta, tau, sum);
  return sum;
}

s_real phir_delta_tau3(s_real delta, s_real tau){
  double val = memoize::get_bin0(memoize::phir_delta_tau3, delta, tau);
  if(!std::isnan(val)) return val;
  s_real sum=0;
  //Calculate sums
  for(int i=1; i<8; ++i){
    sum += n[i1(i)]*d[i1(i)]*t[i1(i)]*(t[i1(i)]-1)*(t[i1(i)]-2)*
      s_pow(delta, d[i1(i)]-1)*s_pow(tau, t[i1(i)]-3);
  }
  for(int i=8; i<52; ++i){
    sum += n[i1(i)]*t[i1(i)]*(t[i1(i)] - 1)*(t[i1(i)] - 2)*
      s_pow(delta, d[i1(i)] - 1)*s_pow(tau, t[i1(i)] - 3)*Eb*
      s_exp(-s_pow(delta, c[i8(i)]));;
  }
  for(int i=52; i<55; ++i){
    sum += n[i1(i)]*s_pow(delta, d[i1(i)])*Ga*(
      t[i1(i)]*(t[i1(i)] - 1)*s_pow(tau, t[i1(i)] - 2)*Fa*Gb +
      2.0*t[i1(i)]*s_pow(tau, t[i1(i)] - 1)*Fa_t*Gb +
      2.0*t[i1(i)]*s_pow(tau, t[i1(i)] - 1)*Fa*Gb_t +
      s_pow(tau, t[i1(i)])*Fa_tt*Gb +
      2.0*s_pow(tau, t[i1(i)])*Fa_t*Gb_t +
      s_pow(tau, t[i1(i)])*Fa*Gb_tt);
  }
  memoize::add_bin0(memoize::phir_delta_tau3, delta, tau, sum);
  return sum;
}

s_real phir_tau3(s_real delta, s_real tau){
  double val = memoize::get_bin0(memoize::phir_tau3, delta, tau);
  if(!std::isnan(val)) return val;
  s_real sum = 0;
  //Calculate sums
  for(int i=1; i<8; ++i){
    sum += n[i1(i)]*t[i1(i)]*(t[i1(i)]-1)*(t[i1(i)]-2)*
      s_pow(delta, d[i1(i)])*s_pow(tau, t[i1(i)]-3);
  }
  for(int i=8; i<52; ++i){
    sum += n[i1(i)]*t[i1(i)]*(t[i1(i)]-1)*(t[i1(i)]-2)*
      s_pow(delta, d[i1(i)])*s_pow(tau, t[i1(i)]-3)*s_exp(-s_pow(delta,c[i8(i)]));
  }
  for(int i=52; i<55; ++i){
    sum += n[i1(i)]*s_pow(delta, d[i1(i)])*(
      t[i1(i)]*s_pow(tau, t[i1(i)]-1)*Fa*(Ha + Hb) +
      s_pow(tau, t[i1(i)])*Fa_t*(Ha + Hb) +
      s_pow(tau, t[i1(i)])*Fa*(Ha_t + Hb_t));
  }
  memoize::add_bin0(memoize::phir_tau3, delta, tau, sum);
  return sum;
}

s_real phir_tau4(s_real delta, s_real tau){
  double val = memoize::get_bin0(memoize::phir_tau4, delta, tau);
  if(!std::isnan(val)) return val;
  s_real sum = 0;
  //Calculate sums
  for(int i=1; i<8; ++i){
    sum += n[i1(i)]*t[i1(i)]*(t[i1(i)]-1)*(t[i1(i)]-2)*(t[i1(i)]-3)*
      s_pow(delta, d[i1(i)])*s_pow(tau, t[i1(i)]-4);
  }
  for(int i=8; i<52; ++i){
    sum += n[i1(i)]*t[i1(i)]*(t[i1(i)]-1)*(t[i1(i)]-2)*(t[i1(i)]-3)*
      s_pow(delta, d[i1(i)])*s_pow(tau, t[i1(i)]-4)*s_exp(-s_pow(delta,c[i8(i)]));;
  }
  for(int i=52; i<55; ++i){
    sum += n[i1(i)]*s_pow(delta, d[i1(i)])*(
      t[i1(i)]*(t[i1(i)]-1)*s_pow(tau, t[i1(i)]-2)*Fa*(Ha + Hb) +
      2.0*t[i1(i)]*s_pow(tau, t[i1(i)]-1)*Fa_t*(Ha + Hb) +
      2.0*t[i1(i)]*s_pow(tau, t[i1(i)]-1)*Fa*(Ha_t + Hb_t) +
      2.0*s_pow(tau, t[i1(i)])*Fa_t*(Ha_t + Hb_t) +
      s_pow(tau, t[i1(i)])*Fa_tt*(Ha + Hb) +
      s_pow(tau, t[i1(i)])*Fa*(Ha_tt + Hb_tt));
  }
  memoize::add_bin0(memoize::phir_tau4, delta, tau, sum);
  return sum;
}

s_real phir_delta3_tau(s_real delta, s_real tau){
  double val = memoize::get_bin0(memoize::phir_delta3_tau, delta, tau);
  if(!std::isnan(val)) return val;
  s_real sum = 0;
  //Calculate sums
  for(int i=1; i<8; ++i){
    sum += n[i1(i)]*d[i1(i)]*(d[i1(i)]-1)*(d[i1(i)]-2)*
      s_pow(delta, d[i1(i)]-3)*t[i1(i)]*s_pow(tau, t[i1(i)] - 1);
  }
  for(int i=8; i<52; ++i){
    sum += n[i1(i)]*t[i1(i)]*s_pow(tau,t[i1(i)] - 1)*(Ea_d*Eg + Ea*Eg_d);
  }
  for(int i=52; i<55; ++i){
    sum += n[i1(i)]*t[i1(i)]*s_pow(tau,t[i1(i)] - 1)*(Fa_d*Fm + Fa*Fm_d) +
      n[i1(i)]*s_pow(tau,t[i1(i)])*(Fa_dt*Fm + Fa_t*Fm_d);
  }
  memoize::add_bin0(memoize::phir_delta3_tau, delta, tau, sum);
  return sum;
}
