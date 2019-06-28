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

#include"iapws95_memo.h"

using namespace memoize;

std::unordered_map<args_bin, memo2, boost::hash<args_bin>> table_bin;
std::unordered_map<args_un, memo1, boost::hash<args_un>> table_un;
std::unordered_map<args_bin, memo0, boost::hash<args_bin>> table_bin0;
std::unordered_map<args_un, memo0, boost::hash<args_un>> table_un0;

unsigned int memoize::add_bin0(unsigned char f, s_real x, s_real y, s_real val){
  if(max_memo == 0) return 0;
  if(table_bin0.size() > max_memo) table_bin0.clear(); //wipe out if it goes oversize
  memo0 *data = &table_bin0[std::make_tuple(f, x, y)];
  data->val = val;
  return table_bin0.size();
}

unsigned int memoize::add_un0(unsigned char f, s_real x, s_real val){
  if(max_memo == 0) return 0;
  if(table_un0.size() > max_memo) table_un0.clear(); //wipe out if it goes oversize
  memo0 *data = &table_un0[std::make_tuple(f, x)];
  data->val = val;
  return table_un0.size();
}

s_real memoize::get_bin0(unsigned char f, s_real x, s_real y){
  if(max_memo == 0) return (s_real)NAN;
  memo0 *data = &table_bin0[std::make_tuple(f, x, y)];
  return data->val;
}

s_real memoize::get_un0(unsigned char f, s_real x){
  if(max_memo == 0) return (s_real)NAN;
  memo0 *data = &table_un0[std::make_tuple(f, x)];
  return data->val;
}

unsigned int memoize::add_bin(unsigned char f, s_real x, s_real y, s_real val,
                     s_real *grad, s_real *hes){
  if(max_memo == 0) return 0;
  if(table_bin.size() > max_memo) table_bin.clear(); //wipe out if it goes oversize
  memo2 *data = &table_bin[std::make_tuple(f, x, y)];
  data->val = val;
  data->grad[0] = grad[0];
  data->grad[1] = grad[1];
  data->hes[0] = hes[0];
  data->hes[1] = hes[1];
  data->hes[2] = hes[2];
  return table_bin.size();
}

unsigned int memoize::add_un(unsigned char f, s_real x,
                    s_real val, s_real *grad, s_real *hes){
  if(max_memo == 0) return 0;
  if(table_un.size() > max_memo) table_un.clear(); //wipe out if it goes oversize
  memo1 *data = &table_un[std::make_tuple(f, x)];
  data->val = val;
  data->grad[0] = grad[0];
  data->hes[0] = hes[0];
  return table_un.size();
}

s_real memoize::get_bin(unsigned char f, s_real x, s_real y, s_real *grad, s_real *hes){
  if(max_memo == 0) return (s_real)NAN;
  memo2 *data = &table_bin[std::make_tuple(f, x, y)];
  if(!std::isnan(data->val)){
    if(grad!=NULL){
      grad[0] = data->grad[0];
      grad[1] = data->grad[1];
    }
    if(hes!=NULL){
      hes[0] = data->hes[0];
      hes[1] = data->hes[1];
      hes[2] = data->hes[2];
    }
  }
  return data->val;
}

s_real memoize::get_un(unsigned char f, s_real x, s_real *grad, s_real *hes){
  if(max_memo == 0) return (s_real)NAN;
  memo1 *data = &table_un[std::make_tuple(f, x)];
  if(!std::isnan(data->val)){
    if(grad!=NULL)grad[0] = data->grad[0];
    if(hes!=NULL) hes[0] = data->hes[0];
  }
  return data->val;
}
