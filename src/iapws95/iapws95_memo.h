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

/*-------------------------------------------------
 Simple memoization for IAPWS95 calcualtions

 Author: John Eslick
-------------------------------------------------*/

#include<unordered_map>
#include<boost/functional/hash.hpp>
#include"iapws95_config.h"

#ifndef _INCLUDE_IAPWS95_MEMO_H_
#define _INCLUDE_IAPWS95_MEMO_H_

namespace memoize{
  static const unsigned int max_memo=MAX_MEMO;

  typedef struct{  // storage type for no derivatives
    s_real val = (s_real)NAN;
  } memo0;

  typedef struct{ // storage to derivatives w.r.t. 1 var
    s_real val = (s_real)NAN;
    s_real grad[1] = {(s_real)NAN};
    s_real hes[1] = {(s_real)NAN};
  } memo1;

  typedef struct{ // storage to derivatives w.r.t. 2 vars
    s_real val = (s_real)NAN;
    s_real grad[2] = {(s_real)NAN, (s_real)NAN};
    s_real hes[3] = {(s_real)NAN, (s_real)NAN, (s_real)NAN};
  } memo2;

  typedef std::tuple<unsigned char, s_real, s_real> args_bin;
  typedef std::tuple<unsigned char, s_real> args_un;

  //binary functions without derivatives
  static const unsigned char
    phir = 1,
    phir_delta = 2,
    phir_tau = 3,
    phir_delta2 = 4,
    phir_delta_tau = 5,
    phir_tau2 = 6,
    phir_delta3 = 7,
    phir_delta2_tau = 8,
    phir_delta_tau2 = 9,
    phir_delta4 = 10,
    phir_delta2_tau2 = 11,
    phir_delta3_tau = 12,
    phir_delta_tau3 = 13,
    phir_tau3 = 14,
    phir_tau4 = 15;
  // binary functions with derivatives
  static const unsigned char
    P_FUNC = 1,
    H_FUNC = 2,
    S_FUNC = 3,
    delta_liq = 3,
    DV_FUNC = 4;
  // unary functions with derivatives
  static const unsigned char
    P_SAT_FUNC = 1,
    DL_SAT_FUNC = 2,
    DV_SAT_FUNC = 3,
    TAU_SAT_FUNC = 4;

  unsigned int add_bin(unsigned char f, s_real x, s_real y, s_real val, s_real *grad, s_real *hes);
  unsigned int add_un(unsigned char f, s_real x, s_real val, s_real *grad,  s_real *hes);
  s_real get_bin(unsigned char f, s_real x, s_real y, s_real *grad, s_real *hes);
  s_real get_un(unsigned char f, s_real x, s_real *grad, s_real *hes);

  unsigned int add_bin0(unsigned char f, s_real x, s_real y, s_real val);
  unsigned int add_un0(unsigned char f, s_real x, s_real val);
  s_real get_bin0(unsigned char f, s_real x, s_real y);
  s_real get_un0(unsigned char f, s_real x);
}

#endif
