
#include <stdio.h>
#include <math.h>

#include "iapws95.h"
#include "iapws95_param.h"


int test_point(s_real T, s_real rho){
  int err;
  s_real tau = f_tau(T);
  s_real delta = f_delta(rho);
  s_real delta_l, delta_v;
  err = sat(tau, &delta_l, &delta_v);
  printf("T = %f, rho = %f\n", (double)T, (double)rho);
  printf("              |       0 (Ideal) |    r (Residual) \n");
  printf("--------------+-----------------+-----------------\n");
  printf("phi           | %+15.8e | %+15.8e\n", (double)phi0(delta, tau), (double)phir(delta, tau));
  printf("phi_delta     | %+15.8e | %+15.8e\n", (double)phi0_delta(delta), (double)phir_delta(delta, tau));
  printf("phi_delta2    | %+15.8e | %+15.8e\n", (double)phi0_delta2(delta), (double)phir_delta2(delta, tau));
  printf("phi_tau       | %+15.8e | %+15.8e\n", (double)phi0_tau(tau), (double)phir_tau(delta, tau));
  printf("phi_tau2      | %+15.8e | %+15.8e\n", (double)phi0_tau2(tau), (double)phir_tau2(delta, tau));
  printf("phi_delta_tau | %+15.8e | %+15.8e\n", (double)phi0_delta_tau, (double)phir_delta_tau(delta, tau));
  printf("\nSome other properties\n");
  printf("p = %+15.8e kPa\n", (double)p(delta, tau));
  printf("h = %+15.8e kJ/kg\n", (double)h(delta, tau));
  printf("s = %+15.8e kJ/kg/K\n", (double)s(delta, tau));
  printf("u = %+15.8e kJ/kg\n", (double)u(delta, tau));
  printf("f = %+15.8e kJ/kg\n", (double)f(delta, tau));
  printf("g = %+15.8e kJ/kg\n", (double)g(delta, tau));
  printf("v = %+15.8e m3/kg\n", (double)(1.0/rho));
  printf("\nSaturated properties at T = %f K\n", (double)T);
  printf("rho_l_sat (aux. estimate) = %+15.8e kg/m3\n", (double)(rho_c*delta_sat_l_approx(tau)));
  printf("rho_v_sat (aux. estimate) = %+15.8e kg/m3\n", (double)(rho_c*delta_sat_v_approx(tau)));
  printf("rho_l_sat = %+15.8e kg/m3\n", (double)(delta_l*rho_c));
  printf("rho_v_sat = %+15.8e kg/m3\n", (double)(delta_v*rho_c));
  printf("psat = %+15.8e kPa\n", (double)p(delta_v, tau));
  return err;
}

int test_15_6(){
  printf("\n\nTable 6 Conditions\n");
  return test_point(500, 838.025);
}

int test_15_7(){
  s_real delta, tau;
  s_real T[4] = {300, 500, 645.5, 900};
  s_real rho[4][4] = {
    {0.9965560e3, 0.1005308e4, 0.1188202e4, 0.0},
    {0.4350000,   0.4532000e1, 0.8380250e3, 0.1084564e4},
    {0.3580000e3, 0.0, 0.0, 0.0},
    {0.2410000, 0.5261500e2, 0.8707690e3}};
  printf("\nTable 7 Conditions\n\n");
  printf("    T (K) |     rho (kg/m3) |         p (kPA) |    cv (kJ/kg/K) |         w (m/s) |     s (kJ/kg/K)\n");
  printf("----------+-----------------+-----------------+-----------------+-----------------+-----------------\n");
  for(int i=0; i<4; ++i){
    for(int j=0; j<4; ++j){
      if(rho[i][j] < 1e-10) continue; //skip zeros, they are not part of the table
      delta = rho[i][j]/rho_c;
      tau = T_c/T[i];
      printf(" %8.2e | %+15.8e | %+15.8e | %+15.8e | %+15.8e | %+15.8e \n",
        (double)T[i], (double)rho[i][j], (double)p(delta, tau), (double)cv(delta, tau), (double)w(delta, tau), (double)s(delta, tau));
    }
  }
  return 0;
}

int test_15_8(){
  int err;
  s_real T[3] = {275, 450, 625};
  s_real rho_l[5], rho_v[5], h_l[5], h_v[5], s_l[5], s_v[5], g_l[5], g_v[5];
  s_real tau, delta_v, delta_l;

  printf("\n\nTable 8 Conditions\n");
  printf("                 |       T = 275 K |       T = 450 K |       T = 625 K\n");
  printf("-----------------+-----------------+-----------------+-----------------\n");
  for(int i=0; i<3; ++i){
    tau = f_tau(T[i]);
    err = sat(tau, &delta_l, &delta_v);
    rho_l[i] = delta_l*rho_c;
    rho_v[i] = delta_v*rho_c;
    h_l[i] = h(delta_l, tau);
    h_v[i] = h(delta_v, tau);
    s_l[i] = s(delta_l, tau);
    s_v[i] = s(delta_v, tau);
    g_l[i] = g(delta_l, tau);
    g_v[i] = g(delta_v, tau);
  }
  printf("     p_sat (kPa) | %+15.8e | %+15.8e | %+15.8e\n",
    (double)p(rho_v[0]/rho_c, T_c/T[0]), (double)p(rho_v[1]/rho_c, T_c/T[1]), (double)p(rho_v[2]/rho_c, T_c/T[2]));
  printf("   rho_l (kg/m3) | %+15.8e | %+15.8e | %+15.8e\n",
    (double)(rho_l[0]), (double)(rho_l[1]), (double)(rho_l[2]));
  printf("   rho_v (kg/m3) | %+15.8e | %+15.8e | %+15.8e\n",
    (double)(rho_v[0]), (double)(rho_v[1]), (double)(rho_v[2]));
  printf("     h_l (kJ/kg) | %+15.8e | %+15.8e | %+15.8e\n",
    (double)h_l[0], (double)h_l[1], (double)h_l[2]);
  printf("     h_v (kJ/kg) | %+15.8e | %+15.8e | %+15.8e\n",
    (double)h_v[0], (double)h_v[1], (double)h_v[2]);
  printf("   s_l (kJ/kg/K) | %+15.8e | %+15.8e | %+15.8e\n",
    (double)s_l[0], (double)s_l[1], (double)s_l[2]);
  printf("   s_v (kJ/kg/K) | %+15.8e | %+15.8e | %+15.8e\n",
    (double)s_v[0], (double)s_v[1], (double)s_v[2]);
  printf("     g_l (kJ/kg) | %+15.8e | %+15.8e | %+15.8e\n",
    (double)g_l[0], (double)g_l[1], (double)g_l[2]);
  printf("     g_v (kJ/kg) | %+15.8e | %+15.8e | %+15.8e\n",
    (double)g_v[0], (double)g_v[1], (double)g_v[2]);
  return err;
}

int test_sat(s_real T1, s_real T2, int n){
  int m;
  s_real delta_l, delta_v, tau, T;
  s_real inc = (T2 - T1)/(n - 1);
  printf("\n\nSaturation curve\n");
  printf("       T (K) |   rho_l (kg/m3) |   rho_v (kg/m3) |     g_l (kJ/kg) |     g_v (kJ/kg) |       p_l (kPa) |       p_v (kPa) | iter\n");
  printf("-------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+------\n");
  for(int i=0; i<n; ++i){
    T = T1 + i*inc;
    tau = f_tau(T);
    m = sat(tau, &delta_l, &delta_v);
    printf(" %11.5f | %+15.8e | %+15.8e | %+15.8e | %+15.8e | %+15.8e | %+15.8e | %d\n",
      (double)T, (double)(delta_l*rho_c), (double)(delta_v*rho_c),
      (double)g(delta_l, tau), (double)g(delta_v, tau),
      (double)p(delta_l, tau), (double)p(delta_v, tau), m);
  }
  return 0;
}

int test_sat_tau(s_real P1, s_real P2, int n){
  s_real tau, P;
  s_real inc = (P2 - P1)/(n - 1);
  int nit=0;
  printf("\n\nSaturation curve T as function of P\n");
  printf("         P (kPa) |       T (K) | it \n");
  printf("-----------------+-------------+----\n");
  for(int i=0; i<n; ++i){
    P = P1 + i*inc;
    tau = sat_tau_with_derivs(P, NULL, NULL, &nit);
    printf(" %+15.8e | %11.5f | %d\n", (double)P, (double)(T_c/tau), nit);
  }
  return 0;
}

int test_delta_from_p_tau(){
  s_real pr = 101.325;
  s_real T = 500;
  s_real tau = T_c/T;
  s_real delta, rho, Psat, res, P_c = 2.2064e4;
  s_real grad[2], hes[3];
  int i, j, nit=0;

  int np = 31;
  s_real Plist[32] = {
    50, 100, 101.325, 250, 500, 750, 1000, 2000, 3000, 4000, 5000, 7000, 8000,
    10e3, 12.5e3, 15e3, 17.5e3, 20e3, 22.064e3, 22.5e3, 24.223e3, 25e3,
    30e3, 40e3, 50e3, 75e3, 100e3, 200e3, 400e3, 600e3, 800e3, 1000e3};
  int nt = 24;
  s_real Tlist[25] = {
    271.235, 275, 280, 285, 290, 300, 315, 340, 365, 390, 430, 480, 530, 580,
    610, 630, 646.5, 647.096, 650, 684.218, 725, 850, 1100, 1250, 1273};

  printf("\n\n");
  printf("Vapor density as a function of T and P\n");
  printf("          T (K) |         p (kPa) | rho_vap (kg/m3) |       Residual  | it \n");
  printf("----------------+-----------------+-----------------+-----------------+----\n");
  for(j=0; j<np; ++j){
    for(i=0; i<nt; ++i){
      T = Tlist[i];
      pr = Plist[j];
      tau = T_c/T;
      if(T <= T_c) Psat = sat_p_with_derivs(tau, NULL, NULL);
      else Psat = 1e10;
      if( (pr > Psat && T < T_c) || (T > T_c && pr>P_c) ){
        delta = 4.0/0.0;
      }
      else{
        delta = delta_vap(pr, tau, grad, hes, &nit);
        res = p(delta, tau) - pr;
        rho = rho_c*delta;
        printf("%+15.8e | %+15.8e | %+15.8e | %+15.8e | %d  \n",
          (double)T,
          (double)pr,
          (double)rho,
          (double)res, nit);
      }
    }
  }
  printf("\n\n");
  printf("Liquid density as a function of T and P\n");
  printf("          T (K) |         p (kPa) | rho_liq (kg/m3) |       Residual  | it \n");
  printf("----------------+-----------------+-----------------+-----------------+----\n");
  for(j=0; j<np; ++j){
    for(i=0; i<nt; ++i){
      T = Tlist[i];
      pr = Plist[j];
      tau = T_c/T;
      if(T < T_c) Psat = sat_p_with_derivs(tau, NULL, NULL);
      else Psat = P_c;
      if(pr < Psat){
        delta = 4.0/0.0;
      }
      else{
        delta = delta_liq(pr, tau, grad, hes, &nit);
        res = p(delta, tau) - pr;
        rho = rho_c*delta;
        printf("%+15.8e | %+15.8e | %+15.8e | %+15.8e | %d  \n",
          (double)T,
          (double)pr,
          (double)rho,
          (double)res, nit);
      }
    }
  }
  return 0;
}

int test_spline(){
  s_real grad[1], hes[1];
  s_real delta, tau=T_c/647.0;
  s_real A, B, C, D;
  delta = sat_delta_vap_with_derivs(tau, grad, hes);
  printf("\n\nVap:\n");
  A = (1 - hes[0]/2.0 - grad[0] + hes[0]*tau - delta - hes[0]/2.0*tau*tau + grad[0]*tau)/
      (1 - 3*tau + 3*tau*tau - tau*tau*tau);
  B = hes[0]/2.0 - 3*A*tau;
  C = grad[0] - hes[0]*tau + 3*A*tau*tau;
  D = 1 - A - B - C;
  printf("delta = %+22.15e\n", (double)delta);
  printf("delta' = %+22.15e\n", (double)grad[0]);
  printf("delta'' = %+22.15e\n", (double)hes[0]);
  printf("A = %+22.15e\n", (double)A);
  printf("B = %+22.15e\n", (double)B);
  printf("C = %+22.15e\n", (double)C);
  printf("D = %+22.15e\n", (double)D);
  printf("delta_spline(tau) = %+22.15e\n", (double)(A*tau*tau*tau + B*tau*tau + C*tau + D));
  printf("error(tau) = %+22.15e\n", (double)(A*tau*tau*tau + B*tau*tau + C*tau + D - delta));
  printf("delta_spline'(tau) = %+22.15e\n", (double)(3*A*tau*tau + 2*B*tau + C));
  printf("delta_spline''(tau) = %+22.15e\n", (double)(6*A*tau + 2*B));
  printf("delta_spline(1) = %+22.15e\n", (double)(A + B + C + D));


  delta = sat_delta_liq_with_derivs(tau, grad, hes);
  printf("\n\nLiq:\n");
  A = (1 - hes[0]/2.0 - grad[0] + hes[0]*tau - delta - hes[0]/2.0*tau*tau + grad[0]*tau)/
      (1 - 3*tau + 3*tau*tau - tau*tau*tau);
  B = hes[0]/2.0 - 3*A*tau;
  C = grad[0] - hes[0]*tau + 3*A*tau*tau;
  D = 1 - A - B - C;
  printf("delta = %+22.15e\n", (double)delta);
  printf("delta' = %+22.15e\n", (double)grad[0]);
  printf("delta'' = %+22.15e\n", (double)hes[0]);
  printf("A = %+22.15e\n", (double)A);
  printf("B = %+22.15e\n", (double)B);
  printf("C = %+22.15e\n", (double)C);
  printf("D = %+22.15e\n", (double)D);
  printf("delta_spline(tau) = %+22.15e\n", (double)(A*tau*tau*tau + B*tau*tau + C*tau + D));
  printf("error(tau) = %+22.15e\n", (double)(A*tau*tau*tau + B*tau*tau + C*tau + D - delta));
  printf("delta_spline'(tau) = %+22.15e\n", (double)(3*A*tau*tau + 2*B*tau + C));
  printf("delta_spline''(tau) = %+22.15e\n", (double)(6*A*tau + 2*B));
  printf("delta_spline(1) = %+22.15e\n", (double)(A + B + C + D));
  return 0;
}

int main(){
  test_15_6();
  test_15_7();
  test_15_8();
  test_sat(240, 645, 50);
  test_sat(645, 647.096, 50);
  test_sat(647.093, 647.096, 50);
  test_sat(647.095, 647.096, 100);
  test_sat(647.09599, 647.096, 100);
  test_delta_from_p_tau();
  test_spline();
  test_sat_tau(0.5, P_c, 100);
  printf("\n");
}
