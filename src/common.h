#ifndef _COMMON_
#define _COMMON_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_roots.h>
#include "params.h"
#include "dam_utils.h"
#include "cosmo_mad.h"

typedef struct {
  gsl_interp_accel *intacc;
  gsl_spline *spline;
  double x0,xf;
  double y0,yf;
} SplPar;

typedef struct {
  double om,ol,ob;
  double w0,wa,h0;
  double ns,s8;
  char fname_window[256];
  char fname_bias[256];
  char fname_sbias[256];
  char fname_pk_l[256];
  char fname_pk_nl[256];
  char prefix_out[256];
  int lmax;
  Csm_params *cpar;
  double chi_horizon;
  double chi_LSS;
  double prefac_lensing;
  double dchi;
  int do_nc;
  int do_shear;
  int do_cmblens;
  int has_bg;
  int has_dens;
  int has_rsd;
  int has_lensing;
  SplPar *aofchi;
  SplPar *zofchi;
  SplPar *hofchi;
  SplPar *gfofchi;
  SplPar *fgofchi;
  SplPar *wind_0;
  SplPar *wind_M;
  SplPar *wind_L;
  SplPar *bias;
  SplPar *sbias;
  double *cl_dd,*cl_dl,*cl_dc,*cl_ll,*cl_lc,*cl_cc;
} RunParams;

//Defined in common.c
SplPar *spline_init(int n,double *x,double *y,double y0,double yf);
double spline_eval(double x,SplPar *spl);
void spline_free(SplPar *spl);
RunParams *param_new(void);
void param_free(RunParams *par);

//Defined in cosmo.c
RunParams *init_params(char *fname_ini);

//Defined in transfers.c
double transfer_wrap(int l,double k,RunParams *par,char *trtype);

//Defined in spectra.c
void compute_spectra(RunParams *par);

//Defined in io.c
int read_parameter_file(char *fname,RunParams *par);
void write_output(RunParams *par);

#endif //_COMMON_
