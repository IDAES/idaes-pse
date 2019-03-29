#include "functions.h"
#include <math.h>

void funcadd(AmplExports *ae){
    /* Arguments for addfunc (this is not fully detailed see funcadd.h)
     * 1) Name of function in AMPL
     * 2) Function pointer to C function
     * 3) see FUNCADD_TYPE enum in funcadd.h
     * 4) Number of arguments (the -1 is variable arg list length)
     * 5) Void pointer to function info */
    int typ = FUNCADD_REAL_VALUED;
    addfunc("cbrt", (rfunc)scbrt, typ, 1, NULL);
}

extern real scbrt(arglist *al){
    real x = al->ra[al->at[0]];
    if(al->derivs!=NULL){
      if(fabs(x) < 6e-9) al->derivs[0] = 1e5;
      else al->derivs[0] = pow(cbrt(x), -2.0)/3.0;
      if(al->hes!=NULL){
        al->hes[0] = -2.0*pow(cbrt(x), -5.0)/9.0;
      }
    }
    return cbrt(x);
}
