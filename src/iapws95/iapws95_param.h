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
 This file provides the IAPWS R6-95(2016) parameters.

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

#include "iapws95_config.h"

static const s_real R = 0.46151805;  // Gas constant (kJ/kg/K)

// Critiacal point for water and R
static const s_real T_c = 647.096;   // Critical T (K)
static const s_real rho_c = 322;     // Critical density (kg/m^3)
static const s_real P_c = 2.2064e4;  // Critical Pressure (kPa)

// To generalize the equation of state there are parmeters to set the number
// of terms in each sumation.  So far we are looking at IAPWS95 and Span-Wagner
// They have the same types of terms but different numbers of term.  For both,
// we are currently droping the nonanalytic terms, which we've labed S4 there
// See docs for the form of the terms in S1, S2, S3, and S4.

const unsigned char S1_set[2] = {1, 7};
const unsigned char S2_set[2] = {8, 51};
const unsigned char S3_set[2] = {52, 54};
const unsigned char S4_set[2] = {55, 56};  // we don't currently use these

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

static const s_real param0[] = {  //ideal gas parameters n0 and gamma0
     0.0,              //pad,    0   (nonexistant n0)
    -8.3204464837497,  //n1,     1
     6.6832105275932,  //n2,     2
     3.00632,          //n3,     3
     0.012436,         //n4,     4
     0.97315,          //n5,     5  (nonexistant gamma0)
     1.27950,          //n6,     6  (nonexistant gamma1)
     0.96956,          //n7,     7  (nonexistant gamma2)
     0.24873,          //n8,     8  (nonexistant gamma3)
     1.28728967,       //gamma4, 9
     3.53734222,       //gamma5, 10
     7.74073708,       //gamma6, 11
     9.24437796,       //gamma7, 12
     27.5075105        //gamma8, 13
};

static const s_real *n0 = param0;
static const s_real *gamma0 = param0 + 5;

static const s_real param[] = {
    0,   //0,    c0
    0,   //1,    c1
    0,   //2,    c2
    0,   //3,    c3
    0,   //4,    c4
    0,   //5,    c5
    0,   //6,    c6
    0,   //7,    c7
    1,   //8,    c8
    1,   //9,    c9
    1,   //10,   c10
    1,   //11,
    1,   //12,
    1,   //13,
    1,   //14,
    1,   //15,
    1,   //16,
    1,   //17,
    1,   //18,
    1,   //19,
    1,   //20,
    1,   //21,
    1,   //22,
    2,   //23,
    2,   //24,
    2,   //25,
    2,   //26,
    2,   //27,
    2,   //28,
    2,   //29,
    2,   //30,
    2,   //31,
    2,   //32,
    2,   //33,
    2,   //34,
    2,   //35,
    2,   //36,
    2,   //37,
    2,   //38,
    2,   //39,
    2,   //40,
    2,   //41,
    2,   //42,
    3,   //43,
    3,   //44,
    3,   //45,
    3,   //46,
    4,   //47,
    6,   //48,
    6,   //49,
    6,   //50,
    6,   //51,   c51   d0
    1,   //52,   d1
    1,   //53,   d2
    1,   //54,   d3
    2,   //55,   d4
    2,   //56
    3,   //57
    4,   //58
    1,   //59
    1,   //60
    1,   //61
    2,   //62
    2,   //63
    3,   //64
    4,   //65
    4,   //66
    5,   //67
    7,   //68
    9,   //69
    10,  //70
    11,  //71
    13,  //72
    15,  //73
    1,   //74
    2,   //75
    2,   //76
    2,   //77
    3,   //78
    4,   //79
    4,   //80
    4,   //81
    5,   //82
    6,   //83
    6,   //84
    7,   //85
    9,   //86
    9,   //87
    9,   //88
    9,   //89
    9,   //90
    10,  //91
    10,  //92
    12,  //93
    3,   //94
    4,   //95
    4,   //96
    5,   //97
    14,  //98
    3,   //99
    6,   //100
    6,   //101
    6,    //102
    3,     //103
    3,      //104
    3,       //105,  d54,  t0
    -0.5,    //106,  t1
     0.875,  //107,  t2
     1,      //108,  t3
     0.5,    //109
     0.75,   //110
     0.375,  //111
     1,      //112
     4,      //113
     6,      //114
    12,      //115
     1,      //116
     5,      //117
     4,      //118
     2,      //119
    13,      //120
     9,      //121
     3,      //122
     4,      //123
    11,      //124
     4,      //125
    13,      //126
     1,      //127
     7,      //128
     1,      //129
     9,      //130
    10,      //131
    10,      //132
     3,      //133
     7,      //134
    10,      //135
    10,      //136
     6,      //137
    10,      //138
    10,      //139
     1,      //140
     2,      //141
     3,      //142
     4,      //143
     8,      //144
     6,      //145
     9,       //146
     8,        //147
    16,         //148
    22,          //149
    23,           //150
    23,            //151
    10,             //152
    50,              //153
    44,               //154
    46,                //155
    50,                 //156
     0,                  //157
     1,                   //158
     4,                    //159,  t54,  n0
     0.12533547935523e-1,  //160,  n1
     0.78957634722828e1,   //161,  n2
    -0.87803203303561e1,   //162,  n3
     0.31802509345418,     //163
    -0.26145533859358,     //164
    -0.78199751687981e-2,  //165
     0.88089493102134e-2,  //166
    -0.66856572307965,     //167
     0.20433810950965,     //168
    -0.66212605039687e-4,  //169
    -0.19232721156002,     //170
    -0.25709043003438,     //171
     0.16074868486251,     //172
    -0.40092828925807e-1,  //173
     0.39343422603254e-6,  //174
    -0.75941377088144e-5,  //175
     0.56250979351888e-3,  //176
    -0.15608652257135e-4,  //177
     0.11537996422951e-8,  //178
     0.36582165144204e-6,  //179
    -0.13251180074668e-11, //180
    -0.62639586912454e-9,  //181
    -0.10793600908932,     //182
     0.17611491008752e-1,  //183
     0.22132295167546,     //184
    -0.40247669763528,     //185
     0.58083399985759,     //186
     0.49969146990806e-2,  //187
    -0.31358700712549e-1,  //188
    -0.74315929710341,     //189
     0.47807329915480,     //190
     0.20527940895948e-1,  //191
    -0.13636435110343,     //192
     0.14180634400617e-1,  //193
     0.83326504880713e-2,  //194
    -0.29052336009585e-1,  //195
     0.38615085574206e-1,  //196
    -0.20393486513704e-1,  //197
    -0.16554050063743e-2,  //198
     0.19955571979541e-2,  //199
     0.15870308324157e-3,  //200
    -0.16388568342530e-4,  //201
     0.43613615723811e-1,  //202
     0.34994005463765e-1,  //203
    -0.76788197844621e-1,  //204
     0.22446277332006e-1,  //205
    -0.62689710414685e-4,  //206
    -0.55711118565645e-9,  //207
    -0.19905718354408,     //208
     0.31777497330738,     //209
    -0.11841182425981,     //210
    -0.31306260323435e2,   //211
     0.31546140237781e2,   //212
    -0.25213154341695e4,   //213
    -0.14874640856724,     //214
     0.31806110878444,     //215,  n56,
     20,                   //216,  alpha52
     20,                   //217,  alpha53
     20,                   //218,  alpha54
     1.21,                 //219,  theta52
     1.21,                 //220,  theta53
     1.25,                 //221,  theta54
     1,                    //222,  eps52
     1,                    //223,  eps53
     1,                    //224,  eps54
     150,                  //225,  beta52
     150,                  //226,  beta53
     250,                  //227,  beta54
     0.3,                  //228,  beta55
     0.3,                  //229,  beta56
};

static const s_real *c = param + 0;
static const s_real *d = param + 51;
static const s_real *t = param + 105;
static const s_real *n = param + 159;
static const s_real *alpha = param + 216 - 52;
static const s_real *theta = param + 219 - 52;
static const s_real *eps = param + 222 - 52;
static const s_real *beta = param + 225 -52;

// Functions to convert to dimensionless state variables from T and Rho.
// temperature and density
inline s_real f_tau(s_real T){return T_c/T;}
inline s_real f_delta(s_real rho){return rho/rho_c;}

#endif
