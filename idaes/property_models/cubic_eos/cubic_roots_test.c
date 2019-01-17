#define No_AE_redefs
#include "funcadd.h"
#define Printf printf
#undef fprintf

#include <dlfcn.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void *dlhandle; // dynamic library handle

/* Function pointer type for AMPL user function */
typedef real (*auf_t)(arglist *al);
typedef real (*croot_t)(int phase, int eos, double A, double B,
                        double *derivs, double *hes);

double residual(int eos, double A, double B, double z){
    double b,c,d,u,w;
    if(eos==0){
        u = 2;
        w = -1;}
    else{
        u = 1;
        w = 0;}
    b = -(1 + B - u*B);
    c = A + w*B*B - u*B - u*B*B;
    d = -A*B - w*B*B - w*B*B*B;
    return z*z*z + b*z*z + c*z + d;
}

int main(int argc, char **argv){
    char *filename; // .so file name
    croot_t f;
    double zv, zl, z2, z3;
    double derivs[3], derivs2[3], derivs3[3];
    double fd_dzdAv, fd_dzdBv;
    double fd_dzdAl, fd_dzdBl;
    double fd_d2zdA2v, fd_d2zdB2v, fd_d2zdAdBv;
    double fd_d2zdA2l, fd_d2zdB2l, fd_d2zdAdBl;
    double dzdAv, dzdBv, d2zdA2v, d2zdB2v, d2zdAdBv;
    double dzdAl, dzdBl, d2zdA2l, d2zdB2l, d2zdAdBl;
    double A, B, T, P, res1, res2;
    int eos=0;
    double h=1e-6;
    double hes[6], hes2[6];
    FILE *file_ptr;
    char *line = NULL;
    size_t llen = 1000;

    /* Print some heading */
    printf("\n\n");
    printf("+-----------------------------+\n");
    printf("|   Property library tests    |\n");
    printf("+-----------------------------+\n\n");
    /* Check that file may have been specified */
    if(argc < 2){
        printf("Must specify library file name.\n");
        return 1;
    }
    /* Load so file */
    filename = argv[1];
    printf("Library file: %s\n", filename);
    dlhandle = dlopen(filename, RTLD_NOW|RTLD_GLOBAL);
    if(!dlhandle){ // if didn't load print error
        printf("%s\n", dlerror());
        return 2;
    }
    printf("Loaded library\n");
    f = (croot_t)dlsym(dlhandle, "cubic_root_ext");

    file_ptr = fopen("eos_test.csv", "r");
    getline(&line, &llen, file_ptr);
    while(getline(&line, &llen, file_ptr)!=-1){
      sscanf(line, "%lf, %lf, %lf, %lf", &T, &P, &A, &B);
      printf("\nT=%f, P=%f, A=%f, B=%f\n", T, P, A, B);
      zv = f(1, eos, A, B, derivs, hes);
      dzdAv = derivs[1];
      dzdBv = derivs[2];
      d2zdA2v = hes[2];
      d2zdB2v = hes[5];
      d2zdAdBv = hes[4];
      //First derivative FD
      z2 = f(1, eos, A+h/2.0, B, derivs, hes2);
      z3 = f(1, eos, A-h/2.0, B, derivs, hes2);
      fd_dzdAv = (z2 - z3)/h;
      z2 = f(1, eos, A, B+h/2.0, derivs, hes2);
      z3 = f(1, eos, A, B-h/2.0, derivs, hes2);
      fd_dzdBv = (z2 - z3)/h;
      //Second derivative FD on first derivatives
      z2 = f(1, eos, A+h/2.0, B, derivs2, hes2);
      z3 = f(1, eos, A-h/2.0, B, derivs3, hes2);
      fd_d2zdA2v = (derivs2[1] - derivs3[1])/h;
      fd_d2zdAdBv =(derivs2[2] - derivs3[2])/h;
      z2 = f(1, eos, A, B+h/2.0, derivs2, hes2);
      z3 = f(1, eos, A, B-h/2.0, derivs3, hes2);
      fd_d2zdB2v = (derivs2[2] - derivs3[2])/h;

      zl = f(0, eos, A, B, derivs, hes);
      dzdAl = derivs[1];
      dzdBl = derivs[2];
      d2zdA2l = hes[2];
      d2zdB2l = hes[5];
      d2zdAdBl = hes[4];
      //First derivative FD
      z2 = f(0, eos, A+h/2.0, B, derivs, hes2);
      z3 = f(0, eos, A-h/2.0, B, derivs, hes2);
      fd_dzdAl = (z2 - z3)/h;
      z2 = f(0, eos, A, B+h/2.0, derivs, hes2);
      z3 = f(0, eos, A, B-h/2.0, derivs, hes2);
      fd_dzdBl = (z2 - z3)/h;
      //Second derivative FD on first derivatives
      z2 = f(0, eos, A+h/2.0, B, derivs2, hes2);
      z3 = f(0, eos, A-h/2.0, B, derivs3, hes2);
      fd_d2zdA2l = (derivs2[1] - derivs3[1])/h;
      fd_d2zdAdBl =(derivs2[2] - derivs3[2])/h;
      z2 = f(0, eos, A, B+h/2.0, derivs2, hes2);
      z3 = f(0, eos, A, B-h/2.0, derivs3, hes2);
      fd_d2zdB2l = (derivs2[2] - derivs3[2])/h;

      res1 = residual(eos, A, B, zl);
      res2 = residual(eos, A, B, zv);
      printf("-------------------------------------------------------------\n");
      printf("zl=%f, resl=%f, zv=%f, resv=%f\n", zl, res1, zv, res2);
      printf("Exact: dzv/dA=%f, dzv/dB=%f, d2zv/dA2=%f, d2zv/dAdB=%f, d2zv/dB2=%f\n", dzdAv, dzdBv, d2zdA2v, d2zdAdBv, d2zdB2v);
      printf("F.D.:  dzv/dA=%f, dzv/dB=%f, d2zv/dA2=%f, d2zv/dAdB=%f, d2zv/dB2=%f\n", fd_dzdAv, fd_dzdBv, fd_d2zdA2v, fd_d2zdAdBv, fd_d2zdB2v);
      printf("Exact: dzl/dA=%f, dzl/dB=%f, d2zl/dA2=%f, d2zl/dAdB=%f, d2zl/dB2=%f\n", dzdAl, dzdBl, d2zdA2l, d2zdAdBl, d2zdB2l);
      printf("F.D.:  dzl/dA=%f, dzl/dB=%f, d2zl/dA2=%f, d2zl/dAdB=%f, d2zl/dB2=%f\n", fd_dzdAl, fd_dzdBl, fd_d2zdA2l, fd_d2zdAdBl, fd_d2zdB2l);
    }
    free(line);
    fclose(file_ptr);

    dlclose(dlhandle);
    return 0;
}
