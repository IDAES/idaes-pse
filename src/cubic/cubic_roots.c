/*
  General root calculator for cubic equations of state.
  Author: John Eslick
  Date: November 22, 2016
  Notes:
*/

#include "cubic_roots.h"



/***********************************************************************
 * CUBIC FORMULA FUNCTION
 * The root finding approach given in CRC Standard Mathematical
 * Tables and Formulae 32nd edition pg 68 was used to identify the
 * roots. There maybe a mistake in it depending on your printing,
 * if you want to check on it also check out the errata:
 * http://www.mathtable.com/smtf/
 **********************************************************************/
double cubic_root2(int phase, double b, double c, double d){
  double F, G, H, I, J, K, M, N, P, R, S, T, U;
  double zr[3], z;

  F = (3*c - b*b)/3.0;
  G = (2*b*b*b - 9*b*c + 27*d)/27.0;
  H = G*G/4.0 + F*F*F/27.0;
  P = -b/3.0;
  if(H <= 0.0){
      //Three roots, can include double or triple roots too
      I = sqrt(G*G/4.0 - H);
      J = curoot(I);
      K = acos(-G/2.0/I);
      M = cos(K/3);
      N = sqrt_3*sin(K/3.0);
      zr[0] = P + 2*J*M;
      zr[1] = P - J*(M + N);
      zr[2] = P - J*(M - N);
      //Sort lowest to highest, may not need this, but for saftey...
      if(zr[0] > zr[1]){z = zr[1]; zr[1] = zr[0]; zr[0] = z;}
      if(zr[1] > zr[2]){z = zr[2]; zr[2] = zr[1]; zr[1] = z;}
      if(zr[0] > zr[1]){z = zr[1]; zr[1] = zr[0]; zr[0] = z;}
      //Vapor is the highest and liquid is the lowest
      if(phase == 0) z = zr[0]; //liquid
      else z = zr[2]; //vapor
  }
  else{ //Only one root
      R = -G/2.0 + sqrt(H);
      T = -G/2.0 - sqrt(H);
      S =  curoot(R);
      U = -curoot(-T);
      z = S + U + P;
    }
    return z;
}

double cubic_root(int phase, eos_indx eos, double A, double B,
                  double *derivs, double *hes){
    /*
     * Find solution to the cubic equation:
     * 0 = z^3 - (1+B-uB)z^2 + (A+-uB-uB^2+wB^2)z - AB-wB^2 - wB^3
     *
     * Also using the definition below (no a, a = 1)
     *
     * 0 = z^3 + b*z2 + c*z + d
     * a = 1
     * b = -(1 + B - u*B)
     * c =  (A - u*B - u*B^2 + w*B^2)
     * d = -(A*B + w*B^2 + w*B^3)
     *
     * Return either what should be the liquid root if phase == 0 or the
     * vapor root if phase == 1, and the derivatives of the returned
     * root
     *
     *
     */
    const double u = eos_u[eos];
    const double w = eos_w[eos];
    static double b, c, d, z;

    b = -(1.0 + B - u*B);
    c = A + w*B*B - u*B - u*B*B;
    d = -A*B - w*B*B - w*B*B*B;
    z = cubic_root2(phase, b, c, d);
    if(derivs) cuderiv(eos, 0, A, B, z, derivs, hes);
    return z;
}

double cubic_root_ext(int phase, eos_indx eos, double A, double B,
                      double *derivs, double *hes){

  const double u = eos_u[eos];
  const double w = eos_w[eos];
  double b, c, d, z, a, pm;
  double bext, cext, dext;
  double det;

  b = -(1.0 + B - u*B);
  c = A + w*B*B - u*B - u*B*B;
  d = -A*B - w*B*B - w*B*B*B;

  if(phase == liquid) pm = -1.0; //liquid
  else pm = 1.0;  //vapor

  det = b*b - 3*c;
  if(det > 0){
    // could need the extension on liquid side so check
    a = -1.0/3.0*b + pm*1.0*sqrt(det)/3.0;
    if( ((a*a*a + b*a*a + c*a + d < 0)&&(phase==liquid)) ||
        ((a*a*a + b*a*a + c*a + d > 0)&&(phase==vapor))){
      //use extension liquid
      a = 2*a;
      bext = -b - 3.0*a;
      cext = 3*a*a + 2*b*a + c;
      dext = d - 0.75*a*a*a - 0.5*b*a*a;
      z = cubic_root2(!phase, bext, cext, dext);
      if(derivs) cuderiv(eos, phase+1, A, B, z, derivs, hes);
    }
    else{
      z = cubic_root2(phase, b, c, d);
      if(derivs) cuderiv(eos, 0, A, B, z, derivs, hes);
    }
  }
  else{
    z = cubic_root2(phase, b, c, d);
    if(derivs) cuderiv(eos, 0, A, B, z, derivs, hes);
  }
  return z;
}

/***********************************************************************
 *
 * LIQUID/VAPOR ROOT PROPERTY FUNCTIONS
 *
 **********************************************************************/

real ceos_z_liq(arglist *al){
  /* This finds the liquid root for a cubic EOS
   *
   * It also calculates first and second derivative.  Returns the vapor root
   * if the liquid phase doesn't exist, so it is discontinuous.  The extended
   * version of this may be better, but this has the nice property that where
   * a vapor root doesn't exist any very small vapor component would have the
   * same properties as the liquid. The exisitence of a liquid root doesn't
   * necessarily mean that liquid is present, so the non-existing liquid phase
   * would not always have the same properies as the vapor.
   * The arguments are:
   *     0) eos_index, specifies the specific eos (e.g. Peng-Robinson)
   *     1) A from the general cubic eos
   *     2) B from the general cubic eos
   */
    double A, B;
    eos_indx eos;
    eos = (eos_indx)(al->ra[al->at[0]] + 0.5);  //cast to int is floor
    A = al->ra[al->at[1]];
    B = al->ra[al->at[2]];
    return cubic_root(0, eos, A, B, al->derivs, al->hes);
}

real ceos_z_vap(arglist *al){
    /* This finds the vapor root for a cubic EOS
     *
     * It also calculates first and second derivative.  Returns the liquid root
     * if the vapor phase doesn't exist, so it is discontinuous.  The extended
     * version of this may be better, but this has the nice property that where
     * a vapor root doesn't exist any very small vapor component would have the
     * same properties as the liquid. The exisitence of a vapor root doesn't
     * necessarily mean that vapor is present, so the non-existing vapor phase
     * would not always have the same properies as the vapor.
     *
     * The arguments are:
     *     0) eos_index, specifies the specific eos (e.g. Peng-Robinson)
     *     1) A from the general cubic eos
     *     2) B from the general cubic eos
     */
    double A, B;
    eos_indx eos;
    eos = (eos_indx)(al->ra[al->at[0]] + 0.5);  //cast to int is floor
    A = al->ra[al->at[1]];
    B = al->ra[al->at[2]];
    return cubic_root(1, eos, A, B, al->derivs, al->hes);
}

real ceos_z_liq_extend(arglist *al){
    /* This finds the liquid root for a cubic EOS
     *
     * It also calculates first and second derivative. If the liquid root
     * doesn't exist it uses a smooth extention of the cubic to find a real
     * fake root, which should be nice numerically, and the property
     * calculations. Should be done in such a way that is the fake root is
     * returned the liquid fraction will be zero (or very close to zero).
     *
     * The arguments are:
     *     0) eos_index, specifies the specific eos (e.g. Peng-Robinson)
     *     1) A from the general cubic eos
     *     2) B from the general cubic eos
     */
    double A, B;
    eos_indx eos;
    eos = (eos_indx)(al->ra[al->at[0]] + 0.5);  //cast to int is floor
    A = al->ra[al->at[1]];
    B = al->ra[al->at[2]];
    return cubic_root_ext(0, eos, A, B, al->derivs, al->hes);
}

real ceos_z_vap_extend(arglist *al){
  /* This finds the vapor root for a cubic EOS
   *
   * It also calculates first and second derivative. If the vapor root
   * doesn't exist it uses a smooth extention of the cubic to find a real
   * fake root, which should be nice numerically, and the property
   * calculations. Should be done in such a way that is the fake root is
   * returned the vapor fraction will be zero (or very close to zero).
   *
   * The arguments are:
   *     0) eos_index, specifies the specific eos (e.g. Peng-Robinson)
   *     1) A from the general cubic eos
   *     2) B from the general cubic eos
   */
    double A, B;
    eos_indx eos;
    eos = (eos_indx)(al->ra[al->at[0]] + 0.5);  //cast to int is floor
    A = al->ra[al->at[1]];
    B = al->ra[al->at[2]];
    return cubic_root_ext(1, eos, A, B, al->derivs, al->hes);
}

void funcadd(AmplExports *ae){
    /* Arguments for addfunc (this is not fully detailed see funcadd.h)
     * 1) Name of function in AMPL
     * 2) Function pointer to C function
     * 3) see FUNCADD_TYPE enum in funcadd.h
     * 4) Number of arguments (the -1 is variable arg list length)
     * 5) Void pointer to function info
     */
    //Real value function, and suppress anoying warings that, I think
    //happen when this gets called on an already loaded library.
    int t = FUNCADD_REAL_VALUED;
    addfunc("ceos_z_vap", (rfunc)ceos_z_vap, t, -1, NULL);
    addfunc("ceos_z_liq", (rfunc)ceos_z_liq, t, -1, NULL);
    addfunc("ceos_z_vap_extend", (rfunc)ceos_z_vap_extend, t, -1, NULL);
    addfunc("ceos_z_liq_extend", (rfunc)ceos_z_liq_extend, t, -1, NULL);
}

/***********************************************************************
 *
 * helpful little functions
 *
 **********************************************************************/
double curoot(double x){
    //cube root function that can take negatives
    if(x < 0) return -pow(-x, 1.0/3.0);
    else return pow(x, 1.0/3.0);
}

/***********************************************************************
 * CUBIC ROOT FINDER DERIVATIVES FUNCTION
 **********************************************************************/
int cubic_derivs(double b, double c, double z, double *grad, double *hes){
      // 0 = z^3 + b*z^2 + c*z + d
      // grad[0] = dz/db, grad[1] = dz/dc, grad[2] = dz/dd
      // hes[0] = d2z/db2, hes[1] = d2z/db/dc, hes[3] = d2z/db/dd
      // hes[2] = d2z/dc2, hes[4] = d2z/dc/dd, hes[5] = d2z/dd2
      // derivatives don't depend on d, so don't need a d argument
      //   (but derivatives with respect to d are not 0)
      double Y = 3.0*z*z + 2.0*b*z + c;
      double Y2 = Y*Y;

      grad[0] = -z*z/Y;
      grad[1] = -z/Y;
      grad[2] = -1.0/Y;

      hes[0] = -2*z*grad[0]/Y + z*z/Y2*(6.0*z*grad[0] + 2.0*b*grad[0] + 2.0*z);
      hes[1] = -2*z*grad[1]/Y + z*z/Y2*(6.0*z*grad[1] + 2.0*b*grad[1] + 1.0);
      hes[3] = -2*z*grad[2]/Y + z*z/Y2*(6.0*z*grad[2] + 2.0*b*grad[2]);

      hes[2] = -grad[1]/Y + z/Y2*(6.0*z*grad[1] + 2.0*b*grad[1] + 1.0);
      hes[4] = -grad[2]/Y + z/Y2*(6.0*z*grad[2] + 2.0*b*grad[2]);

      hes[5] = 1.0/Y2*(6.0*z*grad[2] + 2.0*b*grad[2]);
      return 0;
}

int ext_cubic_derivs(int phase, double b, double c, double z, double *grad, double *hes){
      // derivatives with respect to coefficents for an extended cubic
      // that produces real root answers where real roots shouldn't exist
      // this should help numerically when solving problems with disappearing
      // phases.
      // grad[0] = dz/db, grad[1] = dz/dc, grad[2] = dz/dd
      // hes[0] = d2z/db2, hes[1] = d2z/db/dc, hes[3] = d2z/db/dd
      // hes[2] = d2z/dc2, hes[4] = d2z/dc/dd, hes[5] = d2z/dd2
      double a, bext, cext, pm, det;
      double dadb, dadc;
      double dbextdb, dbextdc, dbextdd;
      double dcextdb, dcextdc, dcextdd;
      double ddextdb, ddextdc, ddextdd;

      double d2bextdb2, d2bextdbdc, d2bextdbdd;
      double d2bextdc2, d2bextdcdd;
      double d2bextdd2;

      double d2cextdb2, d2cextdbdc, d2cextdbdd;
      double d2cextdc2, d2cextdcdd;
      double d2cextdd2;

      double d2dextdb2, d2dextdbdc, d2dextdbdd;
      double d2dextdc2, d2dextdcdd;
      double d2dextdd2;

      double d2adb2, d2adbdc; // d2adbdd;
      double d2adc2; //d2adcdd;
      //double d2add2;

      double ddb_dzdbext, ddc_dzdbext, ddd_dzdbext;
      double ddb_dzdcext, ddc_dzdcext, ddd_dzdcext;
      double ddb_dzddext, ddc_dzddext, ddd_dzddext;

      double grad1[3];
      double hes1[6];

      if(phase == liquid) pm = -1.0; //liquid
      else pm = 1.0;  //vapor

      det = b*b-3*c;
      a = -2.0/3.0*b + pm*2.0*sqrt(det)/3.0;
      bext = -b - 3.0*a;
      cext = 3*a*a + 2*b*a + c;
      //dext = d - 0.75*a*a*a - 0.5*b*a*a;

      cubic_derivs(bext, cext, z, grad1, hes1);

      dadb = -2.0/3.0 + pm*2.0*b/3.0/sqrt(det);
      dadc = -pm/sqrt(det);

      dbextdb = -1.0 - 3.0*dadb;
      dbextdc = -3.0*dadc;
      dbextdd = 0;

      dcextdb = 6.0*a*dadb + 2.0*b*dadb + 2.0*a;
      dcextdc = 6.0*a*dadc + 2.0*b*dadc + 1.0;
      dcextdd = 0;

      ddextdb = -2.25*a*a*dadb - b*a*dadb - 0.5*a*a;
      ddextdc = -2.25*a*a*dadc - b*a*dadc;
      ddextdd = 1;

      grad[0] = grad1[0]*dbextdb + grad1[1]*dcextdb + grad1[2]*ddextdb;
      grad[1] = grad1[0]*dbextdc + grad1[1]*dcextdc + grad1[2]*ddextdc;
      grad[2] = grad1[0]*dbextdd + grad1[1]*dcextdd + grad1[2]*ddextdd;

      d2adb2 = pm*2.0/3.0/sqrt(det) - 2*pm*b*b/3.0/pow(det, 1.5);
      d2adbdc = pm*b/pow(det,1.5);
      //d2adbdd = 0;
      d2adc2 = -pm*1.5/pow(det, 1.5);
      //d2adcdd = 0;
      //d2add2 = 0;

      d2bextdb2 = -3.0*d2adb2;
      d2bextdbdc = -3.0*d2adbdc;
      d2bextdbdd = 0;
      d2bextdc2 = -3.0*d2adc2;
      d2bextdcdd = 0;
      d2bextdd2 = 0;

      d2cextdb2 = 6.0*dadb*dadb + 6.0*a*d2adb2 + 4.0*dadb + 2*b*d2adb2;
      d2cextdbdc = 6.0*dadc*dadb + 6.0*a*d2adbdc + 2.0*b*d2adbdc + 2.0*dadc;
      d2cextdbdd = 0;

      d2cextdc2 = 6.0*dadc*dadc + 6.0*a*d2adc2 + 2.0*b*d2adc2;
      d2cextdcdd = 0;
      d2cextdd2 = 0;

      d2dextdb2 = -4.5*a*dadb*dadb - 2.25*a*a*d2adb2
                  - 2.0*a*dadb - b*dadb*dadb - b*a*d2adb2;
                  // wrong? ddextdb = -2.25*a*a*dadb - b*a*dadb - 0.5*a*a;
      d2dextdbdc = -4.5*a*dadc*dadb - 2.25*a*a*d2adbdc
                   - a*dadc - b*dadc*dadb - b*a*d2adbdc;
      d2dextdbdd = 0;
      d2dextdc2 = -4.5*a*dadc*dadc - 2.25*a*a*d2adc2 - b*dadc*dadc - b*a*d2adc2;
      d2dextdcdd = 0;
      d2dextdd2 = 0;

      ddb_dzdbext = hes1[0]*dbextdb + hes1[1]*dcextdb + hes1[3]*ddextdb;
      ddc_dzdbext = hes1[0]*dbextdc + hes1[1]*dcextdc + hes1[3]*ddextdc;
      ddd_dzdbext = hes1[0]*dbextdd + hes1[1]*dcextdd + hes1[3]*ddextdd;

      ddb_dzdcext = hes1[1]*dbextdb + hes1[2]*dcextdb + hes1[4]*ddextdb;
      ddc_dzdcext = hes1[1]*dbextdc + hes1[2]*dcextdc + hes1[4]*ddextdc;
      ddd_dzdcext = hes1[1]*dbextdd + hes1[2]*dcextdd + hes1[4]*ddextdd;

      ddb_dzddext = hes1[3]*dbextdb + hes1[4]*dcextdb + hes1[5]*ddextdb;
      ddc_dzddext = hes1[3]*dbextdc + hes1[4]*dcextdc + hes1[5]*ddextdc;
      ddd_dzddext = hes1[3]*dbextdd + hes1[4]*dcextdd + hes1[5]*ddextdd;

      hes[0] = ddb_dzdbext*dbextdb + grad1[0]*d2bextdb2
             + ddb_dzdcext*dcextdb + grad1[1]*d2cextdb2
             + ddb_dzddext*ddextdb + grad1[2]*d2dextdb2;  // wrong ?

      hes[1] = ddc_dzdbext*dbextdb + grad1[0]*d2bextdbdc
             + ddc_dzdcext*dcextdb + grad1[1]*d2cextdbdc
             + ddc_dzddext*ddextdb + grad1[2]*d2dextdbdc;

      hes[3] = ddd_dzdbext*dbextdb + grad1[0]*d2bextdbdd
             + ddd_dzdcext*dcextdb + grad1[1]*d2cextdbdd
             + ddd_dzddext*ddextdb + grad1[2]*d2dextdbdd;

      hes[2] = ddc_dzdbext*dbextdc + grad1[0]*d2bextdc2
             + ddc_dzdcext*dcextdc + grad1[1]*d2cextdc2
             + ddc_dzddext*ddextdc + grad1[2]*d2dextdc2;

      hes[4] = ddd_dzdbext*dbextdc + grad1[0]*d2bextdcdd
             + ddd_dzdcext*dcextdc + grad1[1]*d2cextdcdd
             + ddd_dzddext*ddextdc + grad1[2]*d2dextdcdd;

      hes[5] = ddd_dzdbext*dbextdd + grad1[0]*d2bextdd2
             + ddd_dzdcext*dcextdd + grad1[1]*d2cextdd2
             + ddd_dzddext*ddextdd + grad1[2]*d2dextdd2;

      return 0;
}


int AB_derivs(eos_indx eos, char ext, double A, double B, double z, double *grad, double *hes){
    // Calculate derivatives with respect to A and B from the general cubic
    // equation of state.

    double b,c;
    double dbdA, dbdB, dcdA, dcdB, dddA, dddB;
    double d2bdA2, d2bdAdB, d2bdB2;
    double d2cdA2, d2cdAdB, d2cdB2;
    double d2ddA2, d2ddAdB, d2ddB2;
    double ddA_dzdb, ddA_dzdc, ddA_dzdd;
    double ddB_dzdb, ddB_dzdc, ddB_dzdd;
    double grad1[3];
    double hes1[6];

    const double u = eos_u[eos];
    const double w = eos_w[eos];

    b = -(1.0 + B - u*B);
    c = A + w*B*B - u*B - u*B*B;
    //d = -A*B - w*B*B - w*B*B*B;

    if(ext==0) cubic_derivs(b, c, z, grad1, hes1);
    else if(ext==1) ext_cubic_derivs(liquid, b, c, z, grad1, hes1);
    else if(ext==2) ext_cubic_derivs(vapor, b, c, z, grad1, hes1);

    dbdA = 0.0;
    dbdB = u - 1.0;
    dcdA = 1.0;
    dcdB = 2*w*B - u - 2.0*u*B;
    dddA = -B;
    dddB = -A - 2.0*w*B - 3.0*w*B*B;

    d2bdA2 = 0.0;
    d2bdAdB = 0.0;
    d2bdB2 = 0.0;
    d2cdA2 = 0.0;
    d2cdAdB = 0.0;
    d2cdB2 = 2.0*w - 2.0*u;
    d2ddA2 = 0.0;
    d2ddAdB = -1.0;
    d2ddB2 = -2.0*w - 6.0*w*B;

    //grad1[0] = dz/db, grad[1] = dz/dc, grad[2] = dz/dd;
    grad[0] = dbdA*grad1[0] + dcdA*grad1[1] + dddA*grad1[2]; //dzdA
    grad[1] = dbdB*grad1[0] + dcdB*grad1[1] + dddB*grad1[2]; //dzdB

    // hes1[0] = d2z/db2, hes1[1] = d2z/db/dc, hes1[3] = d2z/db/dd
    // hes1[2] = d2z/dc2, hes1[4] = d2z/dc/dd
    // hes1[5] = d2z/dd2
    ddA_dzdb = dbdA*hes1[0] + dcdA*hes1[1] + dddA*hes1[3];
    ddB_dzdb = dbdB*hes1[0] + dcdB*hes1[1] + dddB*hes1[3];
    ddA_dzdc = dbdA*hes1[1] + dcdA*hes1[2] + dddA*hes1[4];
    ddB_dzdc = dbdB*hes1[1] + dcdB*hes1[2] + dddB*hes1[4];
    ddA_dzdd = dbdA*hes1[3] + dcdA*hes1[4] + dddA*hes1[5];
    ddB_dzdd = dbdB*hes1[3] + dcdB*hes1[4] + dddB*hes1[5];

    // hes[0] = d2z/dA2, hes[1] = d2z/dA/dB, hes[2] = d2z/dB2;
    hes[0] = d2bdA2*grad1[0] + dbdA*ddA_dzdb
           + d2cdA2*grad1[1] + dcdA*ddA_dzdc
           + d2ddA2*grad1[2] + dddA*ddA_dzdd;

    hes[1] = d2bdAdB*grad1[0] + dbdA*ddB_dzdb
           + d2cdAdB*grad1[1] + dcdA*ddB_dzdc
           + d2ddAdB*grad1[2] + dddA*ddB_dzdd;

    hes[2] = d2bdB2*grad1[0] + dbdB*ddB_dzdb
           + d2cdB2*grad1[1] + dcdB*ddB_dzdc
           + d2ddB2*grad1[2] + dddB*ddB_dzdd;

    return 0;
}

int cuderiv(eos_indx eos, char ext, double A, double B, double z, double *derivs, double *hes){
    // implicit calculation of first and second partial derivatives of z with
    // respect to A and B for the roots found in the cubic_root function
    double grad1[2], hes1[3];
    AB_derivs(eos, ext, A, B, z, grad1, hes1);
    derivs[0] = 0; //wrt the eos parameter
    derivs[1] = grad1[0]; // dz/dA
    derivs[2] = grad1[1]; // dz/dB

    hes[0] = 0; //wrt eos parameter
    hes[1] = 0; //wrt eos parameter
    hes[3] = 0; //wrt eos parameter

    hes[2] = hes1[0];
    hes[4] = hes1[1];
    hes[5] = hes1[2];

    return 0;
}