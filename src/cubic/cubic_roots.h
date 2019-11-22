/*
  General root calculator for cubic equations of state.

  Author: John Eslick
  Date: November 22, 2016

  Notes:
*/
#include <time.h>
#include <math.h>
#include "funcadd.h"
#undef fprintf
#include <stdio.h>

/***********************************************************************
 *
 * LOGGING (complie with logging if file is defined)
 *   this keeps detaled information about every function evaluation
 *   and is just for serious debugging.
 *
 **********************************************************************/
//#define LOG1 "c.log"

/***********************************************************************
 *
 * CONSTANTS
 *
 **********************************************************************/
static const double deriv_cap=1e9; // Largest derivative magnitude allowed
static const double sqrt_3=1.73205080756888; //square root of 3

typedef enum{ // Equations of state
    EOS_PR = 0,
    EOS_SRK = 1,
    EOS_END = 2
}eos_indx;

static const double eos_u[EOS_END] = {+2.0, +1.0};
static const double eos_w[EOS_END] = {-1.0, +0.0};

static const char liquid = 0;
static const char vapor = 1;

double curoot(double x);
int cubic_derivs(double b, double c, double z, double *grad, double *hes);
int ext_cubic_derivs(int phase, double b, double c, double z, double *grad, double *hes);
int AB_derivs(eos_indx eos, char ext, double A, double B, double z, double *grad, double *hes);
int cuderiv(eos_indx eos, char ext, double A, double B, double z, double *derivs, double *hes);
