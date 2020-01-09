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
 Terms for breaking down derivatives into managable parts.  Will add
 documentation to explan calculations at some point.

 Author: John Eslick
 File: iapws95_deriv_parts.h
------------------------------------------------------------------------------*/

#include "iapws95_external.h"
#include "iapws95_param.h"

#ifndef _INCLUDE_IAPWS95_DERIV_PARTS_H_
#define _INCLUDE_IAPWS95_DERIV_PARTS_H_

/*
ok I know this is probably sort of lazy, but I know there are three parameters
I'm not going to use here and since they are declared as static in the
parameters header file, I get local copies of the pointers anyway.  I'm going to
define a function with those parameter just to avoid the unused variable warning.
I guess it's my stupid way of saying yeah I know, don't worry it's cool.
*/
inline s_real _dummy_unused(unsigned char i){
  return n[i] + n0[i] + gamma0[i];
}


/*
This next set of macros lets me use the functions for derivative parts without
typing the arguments.  The goal is to make the derivative equations more
readable. Of course, that means where these macros are used you need to have
variable names i, delta, and tau. Luckly no one should need to use these.

Look below this section for the actual functions
*/
#define Ea      Ea_F(i, delta)
#define Ea_d    Ea_d_F(i, delta)
#define Ea_dd   Ea_dd_F(i, delta)
#define Eb      Eb_F(i, delta)
#define Eb_d    Eb_d_F(i, delta)
#define Eb_dd   Eb_dd_F(i, delta)
#define Ec      Ec_F(i, delta)
#define Ec_d    Ec_d_F(i, delta)
#define Ec_dd   Ec_dd_F(i, delta)
#define Ed      Ed_F(i, delta)
#define Ed_d    Ed_d_F(i, delta)
#define Ed_dd   Ed_dd_F(i, delta)
#define Ee      Ee_F(i, delta)
#define Ee_d    Ee_d_F(i, delta)
#define Ee_dd   Ee_dd_F(i, delta)
#define Ef      Ef_F(i, delta)
#define Ef_d    Ef_d_F(i, delta)
#define Ef_dd   Ef_dd_F(i, delta)
#define Eg      Eg_F(i, delta)
#define Eg_d    Eg_d_F(i, delta)
#define Eg_dd   Eg_dd_F(i, delta)

#define Fa      Fa_F(i, delta, tau)
#define Fa_d    Fa_d_F(i, delta, tau)
#define Fa_t    Fa_t_F(i, delta, tau)
#define Fa_dd   Fa_dd_F(i, delta, tau)
#define Fa_dt   Fa_dt_F(i, delta, tau)
#define Fa_tt   Fa_tt_F(i, delta, tau)
#define Fb      Fb_F(i, delta)
#define Fb_d    Fb_d_F(i, delta)
#define Fb_dd   Fb_dd_F(i, delta)
#define Fc      Fc_F(i, delta)
#define Fc_d    Fc_d_F(i, delta)
#define Fc_dd   Fc_dd_F(i, delta)
#define Fd      Fd_F(i, delta)
#define Fd_d    Fd_d_F(i, delta)
#define Fd_dd   Fd_dd_F(i, delta)
#define Fe      Fe_F(i, delta)
#define Fe_d    Fe_d_F(i, delta)
#define Fe_dd   Fe_dd_F(i, delta)
#define Ff      Ff_F(i, delta)
#define Ff_d    Ff_d_F(i, delta)
#define Ff_dd   Ff_dd_F(i, delta)
#define Fg      Fg_F(i, delta)
#define Fg_d    Fg_d_F(i, delta)
#define Fg_dd   2.0
#define Fh      Fh_F(i, delta)
#define Fh_d    Fh_d_F(i, delta)
#define Fh_dd   Fh_dd_F(i, delta)
#define Fk      Fk_F(i, delta)
#define Fk_d    1.0
#define Fk_dd   0.0
#define Fm      Fm_F(i, delta)
#define Fm_d    Fm_d_F(i, delta)
#define Fm_dd   Fm_dd_F(i, delta)

#define Ga      Ga_F(i, delta)
#define Gb      Gb_F(i, tau)
#define Gb_t    Gb_t_F(i, tau)
#define Gb_tt   Gb_tt_F(i, tau)

#define Ha      Ha_F(i, tau)
#define Ha_t    Ha_t_F(i, tau)
#define Ha_tt   Ha_tt_F(i, tau)
#define Hb      Hb_F(i, tau)
#define Hb_t    Hb_t_F(i, tau)
#define Hb_tt   Hb_tt_F(i, tau)

#define XA      XA_F(delta, tau)
#define XA_d    XA_d_F(tau)
#define XA_t    XA_t_F(delta, tau)
#define XA_dd   0
#define XA_dt   XA_dt_F(tau)
#define XA_tt   XA_tt_F(delta, tau)

#define XB      XB_F(tau)
#define XB_t    XB_t_F(tau)
#define XB_tt   XB_tt_F(tau)

#define XC      XC_F(delta, tau)
#define XC_d    XC_d_F(delta, tau)
#define XC_t    XC_t_F(delta, tau)
#define XC_dd   XC_dd_F(delta, tau)
#define XC_dt   XC_dt_F(delta, tau)
#define XC_tt   XC_tt_F(delta, tau)

#define XD      XD_F(delta, tau)
#define XD_d    XD_d_F(delta, tau)
#define XD_t    XD_t_F(delta, tau)
#define XD_dd   XD_dd_F(delta, tau)
#define XD_dt   XD_dt_F(delta, tau)
#define XD_tt   XD_tt_F(delta, tau)

#define XE      XE_F(delta, tau)
#define XE_d    XE_d_F(delta, tau)
#define XE_t    XE_t_F(delta, tau)
#define XE_dd   XE_dd_F(delta, tau)
#define XE_dt   XE_dt_F(delta, tau)
#define XE_tt   XE_tt_F(delta, tau)

#define XF      XF_F(delta, tau)
#define XF_d    XF_d_F(delta, tau)
#define XF_t    XF_t_F(delta, tau)
#define XF_dd   XF_dd_F(delta, tau)
#define XF_dt   XF_dt_F(delta, tau)
#define XF_tt   XF_tt_F(delta, tau)

#define XG      XG_F(delta, tau)
#define XG_d    XG_d_F(delta, tau)
#define XG_t    XG_t_F(delta, tau)
#define XG_dd   XG_dd_F(delta, tau)
#define XG_dt   XG_dt_F(delta, tau)
#define XG_tt   XG_tt_F(delta, tau)

#define XH      XH_F(delta, tau)
#define XH_d    XH_d_F(delta, tau)
#define XH_t    XH_t_F(delta, tau)
#define XH_dd   XH_dd_F(delta, tau)
#define XH_dt   XH_dt_F(delta, tau)
#define XH_tt   XH_tt_F(delta, tau)

//
#define LAV     LA(delta_v, tau)
#define LAL     LA(delta_l, tau)
#define LAVd    LA_d(delta_v, tau)
#define LAVt    LA_t(delta_v, tau)
#define LALd    LA_d(delta_l, tau)
#define LALt    LA_t(delta_l, tau)
#define LAVdd   LA_dd(delta_v, tau)
#define LAVtt   LA_tt(delta_v, tau)
#define LAVdt   LA_dt(delta_v, tau)
#define LALdd   LA_dd(delta_l, tau)
#define LALtt   LA_tt(delta_l, tau)
#define LALdt   LA_dt(delta_l, tau)

#define LBV     LB(delta_v, tau)
#define LBVt    LB_t(delta_v, tau)
#define LBVd    LB_d(delta_v, tau)
#define LBL     LB(delta_l, tau)
#define LBLt    LB_t(delta_l, tau)
#define LBLd    LB_d(delta_l, tau)

#define LCV     LC(delta_v, tau)
#define LCVt    LC_t(delta_v, tau)
#define LCVd    LC_d(delta_v, tau)
#define LCL     LC(delta_l, tau)
#define LCLt    LC_t(delta_l, tau)
#define LCLd    LC_d(delta_l, tau)

#define LDV     LD(delta_v, tau)
#define LDVt    LD_t(delta_v, tau)
#define LDVd    LD_d(delta_v, tau)
#define LDL     LD(delta_l, tau)
#define LDLt    LD_t(delta_l, tau)
#define LDLd    LD_d(delta_l, tau)

#define LEV     LE(delta_v, tau)
#define LEVt    LE_t(delta_v, tau)
#define LEVd    LE_d(delta_v, tau)
#define LEL     LE(delta_l, tau)
#define LELt    LE_t(delta_l, tau)
#define LELd    LE_d(delta_l, tau)

#define LFV     LF(delta_v, tau)
#define LFVt    LF_t(delta_v, tau)
#define LFVd    LF_d(delta_v, tau)
#define LFL     LF(delta_l, tau)
#define LFLt    LF_t(delta_l, tau)
#define LFLd    LF_d(delta_l, tau)

#define LGM     LG(delta_l, delta_v, tau)
#define LHM     LH(delta_l, delta_v, tau)


// Functions for derivative parts, eventually will document in a paper or report

inline s_real Ea_F(unsigned char i, s_real delta){
  return s_exp(-s_pow(delta,c[i]));}

inline s_real Ea_d_F(unsigned char i, s_real delta){
  return -c[i]*s_pow(delta,c[i] - 1)*Ea;}

inline s_real Ed_F(unsigned char i, s_real delta){
  return s_pow(delta, d[i] - 2);}

inline s_real Ed_d_F(unsigned char i, s_real delta){
  return (d[i] - 2)*s_pow(delta, d[i] - 3);}

inline s_real Eb_F(unsigned char i, s_real delta){
   return d[i] - c[i]*s_pow(delta, c[i]);}

inline s_real Eb_d_F(unsigned char i, s_real delta){
   return -c[i]*c[i]*s_pow(delta, c[i] - 1);}

inline s_real Ec_F(unsigned char i, s_real delta){
   return d[i] - 1 - c[i]*s_pow(delta, c[i]);}

inline s_real Ec_d_F(unsigned char i, s_real delta){
  return -c[i]*c[i]*s_pow(delta, c[i] - 1);}

inline s_real Ee_F(unsigned char i, s_real delta){
  return c[i]*c[i]*s_pow(delta, c[i]);}

inline s_real Ee_d_F(unsigned char i, s_real delta){
  return c[i]*c[i]*c[i]*s_pow(delta, c[i] - 1);}

inline s_real Ef_F(unsigned char i, s_real delta){
  return Eb*Ec;}

inline s_real Ef_d_F(unsigned char i, s_real delta){
  return Eb_d*Ec + Eb*Ec_d;}

inline s_real Eg_F(unsigned char i, s_real delta){
  return Ed*(Ef - Ee);}

inline s_real Eg_d_F(unsigned char i, s_real delta){
  return Ed_d*(Ef - Ee) + Ed*(Ef_d - Ee_d);}

inline s_real Ea_dd_F(unsigned char i, s_real delta){
  return -c[i]*(c[i] - 1)*s_pow(delta,c[i] - 2)*Ea +
         -c[i]*s_pow(delta,c[i] - 1)*Ea_d;}

inline s_real Eb_dd_F(unsigned char i, s_real delta){
  return -c[i]*c[i]*(c[i] - 1)*s_pow(delta, c[i] - 2);}

inline s_real Ec_dd_F(unsigned char i, s_real delta){
  return -c[i]*c[i]*(c[i] - 1)*s_pow(delta, c[i] - 2);}

inline s_real Ed_dd_F(unsigned char i, s_real delta){
  return (d[i] - 2)*(d[i] - 3)*s_pow(delta, d[i] - 4);}

inline s_real Ee_dd_F(unsigned char i, s_real delta){
  return c[i]*c[i]*c[i]*(c[i] - 1)*s_pow(delta, c[i] - 2);}

inline s_real Ef_dd_F(unsigned char i, s_real delta){
  return Eb_dd*Ec + 2.0*Eb_d*Ec_d + Eb*Ec_dd;}

inline s_real Eg_dd_F(unsigned char i, s_real delta){
  return Ed_dd*(Ef - Ee) + 2.0*Ed_d*(Ef_d - Ee_d) + Ed*(Ef_dd - Ee_dd);}

inline s_real Fa_F(unsigned char i, s_real delta, s_real tau){
  return s_exp(-alpha[i]*(delta - eps[i])*(delta - eps[i]) -
    beta[i]*(tau - theta[i])*(tau - theta[i]));}

inline s_real Fa_d_F(unsigned char i, s_real delta, s_real tau){
  return -2.0*alpha[i]*(delta - eps[i])*Fa;}

inline s_real Fb_F(unsigned char i, s_real delta){
  return -2.0*alpha[i]*s_pow(delta, d[i]);}

inline s_real Ff_F(unsigned char i, s_real delta){
  return 4.0*alpha[i]*alpha[i]*s_pow(delta, d[i]);}

inline s_real Fg_F(unsigned char i, s_real delta){
  return (delta - eps[i])*(delta - eps[i]);}

inline s_real Fc_F(unsigned char i, s_real delta){return Ff*Fg;}

inline s_real Fh_F(unsigned char i, s_real delta){
  return -4.0*d[i]*alpha[i]*s_pow(delta, d[i] - 1);}

inline s_real Fk_F(unsigned char i, s_real delta){
  return (delta - eps[i]);}

inline s_real Fd_F(unsigned char i, s_real delta){return Fh*Fk;}

inline s_real Fe_F(unsigned char i, s_real delta){
  return d[i]*(d[i] - 1)*s_pow(delta, d[i] - 2);}

inline s_real Fb_d_F(unsigned char i, s_real delta){
  return -2*alpha[i]*d[i]*s_pow(delta, d[i] - 1);}

inline s_real Ff_d_F(unsigned char i, s_real delta){
  return 4*alpha[i]*alpha[i]*d[i]*s_pow(delta, d[i] - 1);}

inline s_real Fg_d_F(unsigned char i, s_real delta){
  return 2*(delta - eps[i]);}

inline s_real Fc_d_F(unsigned char i, s_real delta){
  return Ff_d*Fg + Ff*Fg_d;}

inline s_real Fh_d_F(unsigned char i, s_real delta){
  return -4.0*d[i]*alpha[i]*(d[i]-1)*s_pow(delta, d[i] - 2);}

inline s_real Fd_d_F(unsigned char i, s_real delta){
  return Fh_d*Fk + Fh*Fk_d;}

inline s_real Fe_d_F(unsigned char i, s_real delta){
  return d[i]*(d[i] - 1)*(d[i] - 2)*s_pow(delta, d[i] - 3);}

inline s_real Fa_dd_F(unsigned char i, s_real delta, s_real tau){
  return -2.0*alpha[i]*(delta - eps[i])*Fa_d - 2.0*alpha[i]*Fa;}

inline s_real Fb_dd_F(unsigned char i, s_real delta){
  return -2.0*alpha[i]*d[i]*(d[i] - 1)*s_pow(delta, d[i] - 2);}

inline s_real Fe_dd_F(unsigned char i, s_real delta){
  return d[i]*(d[i] - 1)*(d[i] - 2)*(d[i] - 3)*s_pow(delta, d[i] - 4);}

inline s_real Ff_dd_F(unsigned char i, s_real delta){
  return 4.0*alpha[i]*alpha[i]*d[i]*(d[i] - 1)*s_pow(delta, d[i] - 2);}

inline s_real Fh_dd_F(unsigned char i, s_real delta){
  return -4.0*d[i]*alpha[i]*(d[i]-1)*(d[i]-2)*s_pow(delta, d[i] - 3);}

inline s_real Fc_dd_F(unsigned char i, s_real delta){
  return Ff_dd*Fg + 2*Ff_d*Fg_d + Ff*Fg_dd;}

inline s_real Fd_dd_F(unsigned char i, s_real delta){
  return Fh_dd*Fk + 2*Fh_d*Fk_d + Fh*Fk_dd;}

inline s_real Fm_F(unsigned char i, s_real delta){
  return Fb + Fc + Fd + Fe;}

inline s_real Fm_d_F(unsigned char i, s_real delta){
  return Fb_d + Fc_d + Fd_d + Fe_d;}

  inline s_real Fm_dd_F(unsigned char i, s_real delta){
    return Fb_dd + Fc_dd + Fd_dd + Fe_dd;}

inline s_real Fa_t_F(unsigned char i, s_real delta, s_real tau){
  return -2.0*beta[i]*(tau - theta[i])*Fa;}

inline s_real Fa_tt_F(unsigned char i, s_real delta, s_real tau){
  return -2.0*beta[i]*Fa - 2.0*beta[i]*(tau - theta[i])*Fa_t;}

inline s_real Fa_dt_F(unsigned char i, s_real delta, s_real tau){
  return -2.0*beta[i]*(tau - theta[i])*Fa_d;}

inline s_real Ga_F(unsigned char i, s_real delta){
  return d[i]/delta - 2.0*alpha[i]*(delta - eps[i]);}

inline s_real Gb_F(unsigned char i, s_real tau){
  return t[i]/tau - 2.0*beta[i]*(tau - theta[i]);}

inline s_real Gb_t_F(unsigned char i, s_real tau){
  return -t[i]/tau/tau - 2.0*beta[i];}

inline s_real Gb_tt_F(unsigned char i, s_real tau){
  return 2.0*t[i]/tau/tau/tau;}

inline s_real Ha_F(unsigned char i, s_real tau){
  return s_pow(t[i]/tau -2.0*beta[i]*(tau - theta[i]), 2);}

inline s_real Hb_F(unsigned char i, s_real tau){
  return -t[i]/tau/tau - 2.0*beta[i];
}

inline s_real Ha_t_F(unsigned char i, s_real tau){
  return 2.0*Hb*(t[i]/tau -2.0*beta[i]*(tau - theta[i]));
}

inline s_real Hb_t_F(unsigned char i, s_real tau){
  return 2.0*t[i]/tau/tau/tau;
}

inline s_real Ha_tt_F(unsigned char i, s_real tau){
  return 2.0*Hb_t*(t[i]/tau -2.0*beta[i]*(tau - theta[i])) + 2.0*Hb*Hb;
}

inline s_real Hb_tt_F(unsigned char i, s_real tau){
  return -6.0*t[i]/tau/tau/tau/tau;
}

inline s_real XA_F(s_real delta, s_real tau){
  return R*delta*rho_c*T_c/tau;
}

inline s_real XA_d_F(s_real tau){
  return R*rho_c*T_c/tau;
}

inline s_real XA_t_F(s_real delta, s_real tau){
  return -R*delta*rho_c*T_c/tau/tau;
}

inline s_real XA_dt_F(s_real tau){
  return -R*rho_c*T_c/tau/tau;
}

inline s_real XA_tt_F(s_real delta, s_real tau){
  return 2.0*R*delta*rho_c*T_c/tau/tau/tau;
}

inline s_real XB_F(s_real tau){
  return R*T_c/tau;
}

inline s_real XB_t_F(s_real tau){
  return -R*T_c/tau/tau;
}

inline s_real XB_tt_F(s_real tau){
  return 2*R*T_c/tau/tau/tau;
}

inline s_real XC_F(s_real delta, s_real tau){
  return delta*phir_delta(delta, tau);
}

inline s_real XC_d_F(s_real delta, s_real tau){
  return phir_delta(delta, tau) + delta*phir_delta2(delta, tau);
}

inline s_real XC_t_F(s_real delta, s_real tau){
  return delta*phir_delta_tau(delta, tau);
}

inline s_real XC_dd_F(s_real delta, s_real tau){
  return 2.0*phir_delta2(delta, tau) + delta*phir_delta3(delta, tau);
}

inline s_real XC_dt_F(s_real delta, s_real tau){
  return phir_delta_tau(delta, tau) + delta*phir_delta2_tau(delta, tau);
}

inline s_real XC_tt_F(s_real delta, s_real tau){
  return delta*phir_delta_tau2(delta, tau);
}

inline s_real XG_F(s_real delta, s_real tau){
  return 1 + delta*phir_delta(delta, tau) - delta*tau*phir_delta_tau(delta, tau);
}

inline s_real XG_d_F(s_real delta, s_real tau){
  return phir_delta(delta, tau) + delta*phir_delta2(delta, tau)
    - tau*phir_delta_tau(delta, tau) - delta*tau*phir_delta2_tau(delta, tau);
}

inline s_real XG_t_F(s_real delta, s_real tau){
  return -delta*tau*phir_delta_tau2(delta, tau);
}

inline s_real XG_dd_F(s_real delta, s_real tau){
  return 2*phir_delta2(delta, tau) + delta*phir_delta3(delta, tau)
    -2*tau*phir_delta2_tau(delta, tau) - delta*tau*phir_delta3_tau(delta, tau);
}

inline s_real XG_dt_F(s_real delta, s_real tau){
  return -tau*phir_delta_tau2(delta, tau) - delta*tau*phir_delta2_tau2(delta, tau);
}

inline s_real XG_tt_F(s_real delta, s_real tau){
  return -delta*phir_delta_tau2(delta, tau) - delta*tau*phir_delta_tau3(delta, tau);
}

inline s_real XD_F(s_real delta, s_real tau){
  return XG*XG;
}

inline s_real XD_d_F(s_real delta, s_real tau){
  return 2.0*XG_d*XG;
}

inline s_real XD_t_F(s_real delta, s_real tau){
  return 2.0*XG_t*XG;
}

inline s_real XD_dd_F(s_real delta, s_real tau){
  return 2.0*XG_dd*XG + 2.0*XG_d*XG_d;
}

inline s_real XD_dt_F(s_real delta, s_real tau){
  return 2.0*XG_dt*XG + 2.0*XG_d*XG_t;
}

inline s_real XD_tt_F(s_real delta, s_real tau){
  return 2.0*XG_tt*XG + 2.0*XG_t*XG_t;
}

inline s_real XE_F(s_real delta, s_real tau){
  return 1 + 2*delta*phir_delta(delta, tau) + delta*delta*phir_delta2(delta, tau);
}

inline s_real XE_d_F(s_real delta, s_real tau){
  return 2*phir_delta(delta, tau) + 4*delta*phir_delta2(delta, tau)
   + delta*delta*phir_delta3(delta, tau);
}

inline s_real XE_t_F(s_real delta, s_real tau){
  return 2*delta*phir_delta_tau(delta, tau) + delta*delta*phir_delta2_tau(delta, tau);
}

inline s_real XE_dd_F(s_real delta, s_real tau){
  return 6*phir_delta2(delta, tau) + 6*delta*phir_delta3(delta, tau) +
    delta*delta*phir_delta4(delta, tau);
}

inline s_real XE_dt_F(s_real delta, s_real tau){
  return 2*phir_delta_tau(delta, tau) + 4*delta*phir_delta2_tau(delta, tau) +
    delta*delta*phir_delta3_tau(delta, tau);
}

inline s_real XE_tt_F(s_real delta, s_real tau){
  return 2*delta*phir_delta_tau2(delta, tau) +
    delta*delta*phir_delta2_tau2(delta, tau);
}

inline s_real XH_F(s_real delta, s_real tau){
  return tau*(phi0_tau(tau) + phir_tau(delta, tau));
}

inline s_real XH_d_F(s_real delta, s_real tau){
  return tau*phir_delta_tau(delta, tau);
}

inline s_real XH_t_F(s_real delta, s_real tau){
  return phi0_tau(tau) + phir_tau(delta, tau) +
    tau*(phi0_tau2(tau) + phir_tau2(delta, tau));
}

inline s_real XH_dd_F(s_real delta, s_real tau){
  return tau*phir_delta2_tau(delta, tau);
}

inline s_real XH_dt_F(s_real delta, s_real tau){
  return phir_delta_tau(delta, tau) + tau*phir_delta_tau2(delta, tau);
}

inline s_real XH_tt_F(s_real delta, s_real tau){
  return 2*phi0_tau2(tau) + 2*phir_tau2(delta, tau) +
    tau*(phi0_tau3(tau) + phir_tau3(delta, tau));
}

inline s_real LA(s_real delta, s_real tau){
  return delta*phir_delta(delta, tau) + 1;
}

inline s_real LA_d(s_real delta, s_real tau){
  return delta*phir_delta2(delta, tau) + phir_delta(delta, tau);
}

inline s_real LA_t(s_real delta, s_real tau){
  return delta*phir_delta_tau(delta, tau);
}

inline s_real LA_dd(s_real delta, s_real tau){
  return delta*phir_delta3(delta, tau) + 2*phir_delta2(delta, tau);
}

inline s_real LA_dt(s_real delta, s_real tau){
  return delta*phir_delta2_tau(delta, tau) + phir_delta_tau(delta, tau);
}

inline s_real LA_tt(s_real delta, s_real tau){
  return delta*phir_delta_tau2(delta, tau);
}

inline s_real LB(s_real delta, s_real tau){
  return LA_d(delta, tau)*delta + LA(delta, tau);
}

inline s_real LB_t(s_real delta, s_real tau){
  return LA_dt(delta, tau)*delta + LA_t(delta, tau);
}

inline s_real LB_d(s_real delta, s_real tau){
  return LA_dd(delta, tau)*delta + 2*LA_d(delta, tau);
}

inline s_real LC(s_real delta, s_real tau){
  return LA_t(delta, tau)*delta;
}

inline s_real LC_t(s_real delta, s_real tau){
  return LA_tt(delta, tau)*delta;
}

inline s_real LC_d(s_real delta, s_real tau){
  return LA_dt(delta, tau)*delta + LA_t(delta, tau);
}

inline s_real LD(s_real delta, s_real tau){
  return phir_tau(delta, tau) + LA_t(delta, tau);
}

inline s_real LD_t(s_real delta, s_real tau){
  return phir_tau2(delta, tau) + LA_tt(delta, tau);
}

inline s_real LD_d(s_real delta, s_real tau){
  return phir_delta_tau(delta, tau) + LA_dt(delta, tau);
}

inline s_real LE(s_real delta, s_real tau){
  return LA_d(delta, tau) + 1.0/delta + phir_delta(delta, tau);
}

inline s_real LE_t(s_real delta, s_real tau){
  return LA_dt(delta, tau) + phir_delta_tau(delta, tau);
}

inline s_real LE_d(s_real delta, s_real tau){
  return LA_dd(delta, tau) - 1.0/delta/delta + phir_delta2(delta, tau);
}

inline s_real LF(s_real delta, s_real tau){
  return 1.0/LB(delta, tau);
}

inline s_real LF_t(s_real delta, s_real tau){
  return -LB_t(delta, tau)/LB(delta, tau)/LB(delta, tau);
}

inline s_real LF_d(s_real delta, s_real tau){
  return -LB_d(delta, tau)/LB(delta, tau)/LB(delta, tau);
}

inline s_real LG(s_real delta_l, s_real delta_v, s_real tau){
  return LEV - LEL*LBV*LFL;
}

inline s_real LdGdt(s_real delta_l, s_real delta_v, s_real tau, s_real ddldt, s_real ddvdt){
  return LEVt + LEVd*ddvdt - (LELt + LELd*ddldt)*LFL*LBV - (LFLt + LFLd*ddldt)*LEL*LBV -
    (LBVt + LBVd*ddvdt)*LEL*LFL;
}

inline s_real LH(s_real delta_l, s_real delta_v, s_real tau){
  return LDL - LDV + LEL*LFL*(LCV - LCL);
}

inline s_real LdHdt(s_real delta_l, s_real delta_v, s_real tau, s_real ddldt, s_real ddvdt){
  return LDLt - LDVt + LDLd*ddldt - LDVd*ddvdt + (LELt + LELd*ddldt)*LFL*(LCV - LCL)
    + (LFLt + LFLd*ddldt)*LEL*(LCV - LCL)
    + LFL*LEL*(LCVt - LCLt + LCVd*ddvdt - LCLd*ddldt);
}

#endif
