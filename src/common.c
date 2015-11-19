#include "common.h"

SplPar *spline_init(int n,double *x,double *y,double y0,double yf)
{
  SplPar *spl=(SplPar *)dam_malloc(sizeof(SplPar));
  spl->intacc=gsl_interp_accel_alloc();
  spl->spline=gsl_spline_alloc(gsl_interp_cspline,n);
  gsl_spline_init(spl->spline,x,y,n);
  spl->x0=x[0];
  spl->xf=x[n-1];
  spl->y0=y0;
  spl->yf=yf;

  return spl;
}

double spline_eval(double x,SplPar *spl)
{
  if(x<=spl->x0)
    return spl->y0;
  else if(x>=spl->xf) 
    return spl->yf;
  else
    return gsl_spline_eval(spl->spline,x,spl->intacc);
}

void spline_free(SplPar *spl)
{
  gsl_spline_free(spl->spline);
  gsl_interp_accel_free(spl->intacc);
  free(spl);
}

RunParams *param_new(void)
{
  RunParams *par=(RunParams *)dam_malloc(sizeof(RunParams));
  par->om=0.3;
  par->ol=0.7;
  par->ob=0.05;
  par->w0=-1.;
  par->wa=0.;
  par->ns=0.96;
  par->s8=0.8;
  sprintf(par->fname_window,"default");
  sprintf(par->fname_bias,"default");
  sprintf(par->fname_sbias,"default");
  sprintf(par->fname_pk,"default");
  sprintf(par->prefix_out,"default");
  par->lmax=100;
  par->cpar=NULL;
  par->chi_horizon=-1.;
  par->chi_LSS=-1.;
  par->prefac_lensing=-1.;
  par->dchi=-1.;
  par->aofchi=NULL;
  par->zofchi=NULL;
  par->hofchi=NULL;
  par->gfofchi=NULL;
  par->fgofchi=NULL;
  par->wind_0=NULL;
  par->wind_M=NULL;
  par->wind_L=NULL;
  par->bias=NULL;
  par->sbias=NULL;
  par->do_nc=0;
  par->do_shear=0;
  par->do_cmblens=0;
  par->has_bg=0;
  par->has_dens=0;
  par->has_rsd=0;
  par->has_lensing=0;
  par->cl_dd=NULL;
  par->cl_dl=NULL;
  par->cl_dc=NULL;
  par->cl_ll=NULL;
  par->cl_lc=NULL;
  par->cl_cc=NULL;
  return par;
}

void param_free(RunParams *par)
{
  csm_params_free(par->cpar);
  if(par->has_bg) {
    spline_free(par->aofchi);
    spline_free(par->zofchi);
    spline_free(par->hofchi);
    spline_free(par->gfofchi);
    spline_free(par->fgofchi);
  }
  if(par->do_nc || par->do_shear)
    spline_free(par->wind_0);
  if(par->do_nc) {
    free(par->cl_dd);
    if(par->do_shear)
      free(par->cl_dl);
    if(par->do_cmblens)
      free(par->cl_dc);
    if(par->has_dens)
      spline_free(par->bias);
    if(par->has_lensing) {
      spline_free(par->sbias);
      spline_free(par->wind_M);
    }
  }
  if(par->do_shear) {
    spline_free(par->wind_L);
    free(par->cl_ll);
    if(par->do_cmblens)
      free(par->cl_lc);
  }
  if(par->do_cmblens)
    free(par->cl_cc);
  free(par);
}
