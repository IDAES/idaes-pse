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
 File: iapws95_asl_funcs.cpp
------------------------------------------------------------------------------*/

#include<stdio.h>
#include<math.h>
#include"iapws95_external.h"
#include"iapws95_phi.h"
#include"iapws95_asl_funcs.h"

void funcadd(AmplExports *ae){
    /* Arguments for addfunc (this is not fully detailed see funcadd.h)
     * 1) Name of function in AMPL
     * 2) Function pointer to C function
     * 3) see FUNCADD_TYPE enum in funcadd.h
     * 4) Number of arguments (the -1 is variable arg list length)
     * 5) Void pointer to function info */
    int typ = FUNCADD_REAL_VALUED;
    addfunc("p", (rfunc)p_asl, typ, 2, NULL);
    addfunc("u", (rfunc)u_asl, typ, 2, NULL);
    addfunc("s", (rfunc)s_asl, typ, 2, NULL);
    addfunc("h", (rfunc)h_asl, typ, 2, NULL);
    addfunc("g", (rfunc)g_asl, typ, 2, NULL);
    addfunc("f", (rfunc)f_asl, typ, 2, NULL);
    addfunc("cv", (rfunc)cv_asl, typ, 2, NULL);
    addfunc("cp", (rfunc)cp_asl, typ, 2, NULL);
    addfunc("w", (rfunc)w_asl, typ, 2, NULL);
    addfunc("hvpt", (rfunc)hvpt_asl, typ, 2, NULL);
    addfunc("hlpt", (rfunc)hlpt_asl, typ, 2, NULL);
    addfunc("tau", (rfunc)tau_asl, typ, 2, NULL);
    addfunc("vf", (rfunc)vf_asl, typ, 2, NULL);
    addfunc("delta_liq", (rfunc)delta_liq_asl, typ, 2, NULL);
    addfunc("delta_vap", (rfunc)delta_vap_asl, typ, 2, NULL);
    addfunc("delta_sat_l", (rfunc)delta_sat_l_asl, typ, 1, NULL);
    addfunc("delta_sat_v", (rfunc)delta_sat_v_asl, typ, 1, NULL);
    addfunc("p_sat", (rfunc)p_sat_asl, typ, 1, NULL);
    addfunc("tau_sat", (rfunc)tau_sat_asl, typ, 1, NULL);
    addfunc("phi0", (rfunc)phi0_asl, typ, 2, NULL);
    addfunc("phi0_delta", (rfunc)phi0_delta_asl, typ, 1, NULL);
    addfunc("phi0_delta2", (rfunc)phi0_delta2_asl, typ, 1, NULL);
    addfunc("phi0_tau", (rfunc)phi0_tau_asl, typ, 1, NULL);
    addfunc("phi0_tau2", (rfunc)phi0_tau2_asl, typ, 1, NULL);
    addfunc("phir", (rfunc)phir_asl, typ, 2, NULL);
    addfunc("phir_delta", (rfunc)phir_delta_asl, typ, 2, NULL);
    addfunc("phir_delta2", (rfunc)phir_delta2_asl, typ, 2, NULL);
    addfunc("phir_tau", (rfunc)phir_tau_asl, typ, 2, NULL);
    addfunc("phir_tau2", (rfunc)phir_tau2_asl, typ, 2, NULL);
    addfunc("phir_delta_tau", (rfunc)phir_delta_tau_asl, typ, 2, NULL);
}

void cast_deriv2(s_real *g1, double *g2, s_real *h1, double *h2){
  if(g2 != NULL){
    g2[0] = (double)g1[0];
    g2[1] = (double)g1[1];
  }
  if(h2 != NULL){
    h2[0] = (double)h1[0];
    h2[1] = (double)h1[1];
    h2[2] = (double)h1[2];
  }
}

void cast_deriv1(s_real *g1, double *g2, s_real *h1, double *h2){
  if(g2 != NULL) g2[0] = (double)g1[0];
  if(h2 != NULL) h2[0] = (double)h1[0];
}

double p_asl(arglist *al){
  s_real f, grad[2], hes[3];
  if(al->derivs==NULL && al->hes==NULL){
    return p_with_derivs(al->ra[al->at[0]], al->ra[al->at[1]], NULL, NULL);}
  else{
    f = p_with_derivs(al->ra[al->at[0]], al->ra[al->at[1]], grad, hes);
    cast_deriv2(grad, al->derivs, hes, al->hes);
    return f;
  }
}

double u_asl(arglist *al){
  s_real f, grad[2], hes[3];
  if(al->derivs==NULL && al->hes==NULL){
    return u_with_derivs(al->ra[al->at[0]], al->ra[al->at[1]], NULL, NULL);}
  else{
    f = u_with_derivs(al->ra[al->at[0]], al->ra[al->at[1]], grad, hes);
    cast_deriv2(grad, al->derivs, hes, al->hes);
    return f;
  }
}

double s_asl(arglist *al){
  s_real f, grad[2], hes[3];
  if(al->derivs==NULL && al->hes==NULL){
    return s_with_derivs(al->ra[al->at[0]], al->ra[al->at[1]], NULL, NULL);}
  else{
    f = s_with_derivs(al->ra[al->at[0]], al->ra[al->at[1]], grad, hes);
    cast_deriv2(grad, al->derivs, hes, al->hes);
    return f;
  }
}

double h_asl(arglist *al){
  s_real f, grad[2], hes[3];
  if(al->derivs==NULL && al->hes==NULL){
    return h_with_derivs(al->ra[al->at[0]], al->ra[al->at[1]], NULL, NULL);}
  else{
    f = h_with_derivs(al->ra[al->at[0]], al->ra[al->at[1]], grad, hes);
    cast_deriv2(grad, al->derivs, hes, al->hes);
    return f;
  }
}

double g_asl(arglist *al){
  s_real f, grad[2], hes[3];
  if(al->derivs==NULL && al->hes==NULL){
    return g_with_derivs(al->ra[al->at[0]], al->ra[al->at[1]], NULL, NULL);}
  else{
    f = g_with_derivs(al->ra[al->at[0]], al->ra[al->at[1]], grad, hes);
    cast_deriv2(grad, al->derivs, hes, al->hes);
    return f;
  }
}

double f_asl(arglist *al){
  s_real f, grad[2], hes[3];
  if(al->derivs==NULL && al->hes==NULL){
    return f_with_derivs(al->ra[al->at[0]], al->ra[al->at[1]], NULL, NULL);}
  else{
    f = f_with_derivs(al->ra[al->at[0]], al->ra[al->at[1]], grad, hes);
    cast_deriv2(grad, al->derivs, hes, al->hes);
    return f;
  }
}

double cv_asl(arglist *al){
  s_real f, grad[2], hes[3];
  if(al->derivs==NULL && al->hes==NULL){
    return cv_with_derivs(al->ra[al->at[0]], al->ra[al->at[1]], NULL, NULL);}
  else{
    f = cv_with_derivs(al->ra[al->at[0]], al->ra[al->at[1]], grad, hes);
    cast_deriv2(grad, al->derivs, hes, al->hes);
    return f;
  }
}

double cp_asl(arglist *al){
  s_real f, grad[2], hes[3];
  if(al->derivs==NULL && al->hes==NULL){
    return cp_with_derivs(al->ra[al->at[0]], al->ra[al->at[1]], NULL, NULL);}
  else{
    f = cp_with_derivs(al->ra[al->at[0]], al->ra[al->at[1]], grad, hes);
    cast_deriv2(grad, al->derivs, hes, al->hes);
    return f;
  }
}

double w_asl(arglist *al){
  s_real f, grad[2], hes[3];
  if(al->derivs==NULL && al->hes==NULL){
    return w_with_derivs(al->ra[al->at[0]], al->ra[al->at[1]], NULL, NULL);}
  else{
    f = w_with_derivs(al->ra[al->at[0]], al->ra[al->at[1]], grad, hes);
    cast_deriv2(grad, al->derivs, hes, al->hes);
    return f;
  }
}

double hvpt_asl(arglist *al){
  s_real f, grad[2], hes[3];
  if(al->derivs==NULL && al->hes==NULL){
    return hvpt_with_derivs(al->ra[al->at[0]], al->ra[al->at[1]], NULL, NULL);}
  else{
    f = hvpt_with_derivs(al->ra[al->at[0]], al->ra[al->at[1]], grad, hes);
    cast_deriv2(grad, al->derivs, hes, al->hes);
    return f;
  }
}

double hlpt_asl(arglist *al){
  s_real f, grad[2], hes[3];
  if(al->derivs==NULL && al->hes==NULL){
    return hlpt_with_derivs(al->ra[al->at[0]], al->ra[al->at[1]], NULL, NULL);}
  else{
    f = hlpt_with_derivs(al->ra[al->at[0]], al->ra[al->at[1]], grad, hes);
    cast_deriv2(grad, al->derivs, hes, al->hes);
    return f;
  }
}

double tau_asl(arglist *al){
  s_real f, grad[2], hes[3];
  if(al->derivs==NULL && al->hes==NULL){
    return tau_with_derivs(al->ra[al->at[0]], al->ra[al->at[1]], NULL, NULL);}
  else{
    f = tau_with_derivs(al->ra[al->at[0]], al->ra[al->at[1]], grad, hes);
    cast_deriv2(grad, al->derivs, hes, al->hes);
    return f;
  }
}

double vf_asl(arglist *al){
  s_real f, grad[2], hes[3];
  if(al->derivs==NULL && al->hes==NULL){
    return vf_with_derivs(al->ra[al->at[0]], al->ra[al->at[1]], NULL, NULL);}
  else{
    f = vf_with_derivs(al->ra[al->at[0]], al->ra[al->at[1]], grad, hes);
    cast_deriv2(grad, al->derivs, hes, al->hes);
    return f;
  }
}

double delta_sat_l_asl(arglist *al){
  s_real f, grad[1], hes[1];
  if(al->derivs==NULL && al->hes==NULL){
    return sat_delta_liq_with_derivs(al->ra[al->at[0]], NULL, NULL);}
  else{
    f = sat_delta_liq_with_derivs(al->ra[al->at[0]], grad, hes);
    cast_deriv1(grad, al->derivs, hes, al->hes);
    return f;
  }
}

double delta_sat_v_asl(arglist *al){
  s_real f, grad[1], hes[1];
  if(al->derivs==NULL && al->hes==NULL){
    return sat_delta_vap_with_derivs(al->ra[al->at[0]], NULL, NULL);}
  else{
    f = sat_delta_vap_with_derivs(al->ra[al->at[0]], grad, hes);
    cast_deriv1(grad, al->derivs, hes, al->hes);
    return f;
  }
}

double p_sat_asl(arglist *al){
  s_real f, grad[1], hes[1];
  if(al->derivs==NULL && al->hes==NULL){
    return sat_p_with_derivs(al->ra[al->at[0]], NULL, NULL);}
  else{
    f = sat_p_with_derivs(al->ra[al->at[0]], grad, hes);
    cast_deriv1(grad, al->derivs, hes, al->hes);
    return f;
  }
}

double tau_sat_asl(arglist *al){
  s_real f, grad[1], hes[1];
  if(al->derivs==NULL && al->hes==NULL){
    return sat_tau_with_derivs(al->ra[al->at[0]], NULL, NULL);}
  else{
    f = sat_tau_with_derivs(al->ra[al->at[0]], grad, hes);
    cast_deriv1(grad, al->derivs, hes, al->hes);
    return f;
  }
}

double delta_liq_asl(arglist *al){
  s_real f, grad[2], hes[3];
  if(al->derivs==NULL && al->hes==NULL){
    return delta_liq(al->ra[al->at[0]], al->ra[al->at[1]], NULL, NULL);}
  else{
    f = delta_liq(al->ra[al->at[0]], al->ra[al->at[1]], grad, hes);
    cast_deriv2(grad, al->derivs, hes, al->hes);
    return f;
  }
}

double delta_vap_asl(arglist *al){
  s_real f, grad[2], hes[3];
  if(al->derivs==NULL && al->hes==NULL){
    return delta_vap(al->ra[al->at[0]], al->ra[al->at[1]], NULL, NULL);}
  else{
    f = delta_vap(al->ra[al->at[0]], al->ra[al->at[1]], grad, hes);
    cast_deriv2(grad, al->derivs, hes, al->hes);
    return f;
  }
}

//with deriv functions for phi, only provided for testing through ASL interface
double phi0_derivs(double delta, double tau, double *grad, double *hes){
  if(grad != NULL){
    grad[0] = phi0_delta(delta);
    grad[1] = phi0_tau(tau);
    if(hes != NULL){
      hes[0] = phi0_delta2(delta);
      hes[1] = 0;
      hes[2] = phi0_tau2(tau);
    }
  }
  return phi0(delta, tau);
}

double phi0_delta_derivs(double delta, double *grad, double *hes){
  if(grad != NULL){
    grad[0] = phi0_delta2(delta);
    grad[1] = 0;
    if(hes != NULL){
      hes[0] = phi0_delta3(delta);
      hes[1] = 0;
      hes[2] = 0;
    }
  }
  return phi0_delta(delta);
}

double phi0_delta2_derivs(double delta, double *grad, double *hes){
  if(grad != NULL){
    grad[0] = phi0_delta3(delta);
    grad[1] = 0;
    if(hes != NULL){
      hes[0] = phi0_delta4(delta);
      hes[1] = 0;
      hes[2] = 0;
    }
  }
  return phi0_delta2(delta);
}

double phi0_tau_derivs(double tau, double *grad, double *hes){
  if(grad != NULL){
    grad[0] = phi0_tau2(tau);
    grad[1] = 0;
    if(hes != NULL){
      hes[0] = phi0_tau3(tau);
      hes[1] = 0;
      hes[2] = 0;
    }
  }
  return phi0_tau(tau);
}

double phi0_tau2_derivs(double tau, double *grad, double *hes){
  if(grad != NULL){
    grad[0] = phi0_tau3(tau);
    grad[1] = 0;
    if(hes != NULL){
      hes[0] = phi0_tau4(tau);
      hes[1] = 0;
      hes[2] = 0;
    }
  }
  return phi0_tau2(tau);
}

double phir_derivs(double delta, double tau, double *grad, double *hes){
  if(grad != NULL){
    grad[0] = phir_delta(delta, tau);
    grad[1] = phir_tau(delta, tau);
    if(hes != NULL){
      hes[0] = phir_delta2(delta, tau);
      hes[1] = phir_delta_tau(delta, tau);
      hes[2] = phir_tau2(delta, tau);
    }
  }
  return phir(delta, tau);
}

double phir_delta_derivs(double delta, double tau, double *grad, double *hes){
  if(grad != NULL){
    grad[0] = phir_delta2(delta, tau);
    grad[1] = phir_delta_tau(delta, tau);
    if(hes != NULL){
      hes[0] = phir_delta3(delta, tau);
      hes[1] = phir_delta2_tau(delta, tau);
      hes[2] = phir_delta_tau2(delta, tau);
    }
  }
  return phir_delta(delta, tau);
}

double phir_delta2_derivs(double delta, double tau, double *grad, double *hes){
  if(grad != NULL){
    grad[0] = phir_delta3(delta, tau);
    grad[1] = phir_delta2_tau(delta, tau);
    if(hes != NULL){
      hes[0] = phir_delta4(delta, tau);
      hes[1] = phir_delta3_tau(delta, tau);
      hes[2] = phir_delta2_tau2(delta, tau);
    }
  }
  return phir_delta2(delta, tau);
}

double phir_tau_derivs(double delta, double tau, double *grad, double *hes){
  if(grad != NULL){
    grad[0] = phir_delta_tau(delta, tau);
    grad[1] = phir_tau2(delta, tau);
    if(hes != NULL){
      hes[0] = phir_delta2_tau(delta, tau);
      hes[1] = phir_delta_tau2(delta, tau);
      hes[2] = phir_tau3(delta, tau);
    }
  }
  return phir_tau(delta, tau);
}

double phir_tau2_derivs(double delta, double tau, double *grad, double *hes){
  if(grad != NULL){
    grad[0] = phir_delta_tau2(delta, tau);
    grad[1] = phir_tau3(delta, tau);
    if(hes != NULL){
      hes[0] = phir_delta2_tau2(delta, tau);
      hes[1] = phir_delta_tau3(delta, tau);
      hes[2] = phir_tau4(delta, tau);
    }
  }
  return phir_tau2(delta, tau);
}

double phir_delta_tau_derivs(double delta, double tau, double *grad, double *hes){
  if(grad != NULL){
    grad[0] = phir_delta2_tau(delta, tau);
    grad[1] = phir_delta_tau2(delta, tau);
    if(hes != NULL){
      hes[0] = phir_delta3_tau(delta, tau);
      hes[1] = phir_delta2_tau2(delta, tau);
      hes[2] = phir_delta_tau3(delta, tau);
    }
  }
  return phir_delta_tau(delta, tau);
}


double phi0_asl(arglist *al){
  return phi0_derivs(al->ra[al->at[0]], al->ra[al->at[1]], al->derivs, al->hes);}
double phi0_delta_asl(arglist *al){
  return phi0_delta_derivs(al->ra[al->at[0]], al->derivs, al->hes);}
double phi0_delta2_asl(arglist *al){
  return phi0_delta2_derivs(al->ra[al->at[0]], al->derivs, al->hes);}
double phi0_tau_asl(arglist *al){
  return phi0_tau_derivs(al->ra[al->at[0]], al->derivs, al->hes);}
double phi0_tau2_asl(arglist *al){
  return phi0_tau2_derivs(al->ra[al->at[0]], al->derivs, al->hes);}
double phir_asl(arglist *al){
  return phir_derivs(al->ra[al->at[0]], al->ra[al->at[1]], al->derivs, al->hes);}
double phir_delta_asl(arglist *al){
  return phir_delta_derivs(al->ra[al->at[0]], al->ra[al->at[1]], al->derivs, al->hes);}
double phir_delta2_asl(arglist *al){
  return phir_delta2_derivs(al->ra[al->at[0]], al->ra[al->at[1]], al->derivs, al->hes);}
double phir_tau_asl(arglist *al){
  return phir_tau_derivs(al->ra[al->at[0]], al->ra[al->at[1]], al->derivs, al->hes);}
double phir_tau2_asl(arglist *al){
  return phir_tau2_derivs(al->ra[al->at[0]], al->ra[al->at[1]], al->derivs, al->hes);}
double phir_delta_tau_asl(arglist *al){
  return phir_delta_tau_derivs(al->ra[al->at[0]], al->ra[al->at[1]], al->derivs, al->hes);}
