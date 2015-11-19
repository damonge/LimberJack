#include "common.h"

typedef struct {
  int l;
  RunParams *par;
  char *tr1;
  char *tr2;
} IntPar;

static double cl_integrand(double lk,void *params)
{
  double d1,d2;
  IntPar *p=(IntPar *)params;
  double k=pow(10.,lk);
  double pk=csm_Pk_linear_0(p->par->cpar,k);
  d1=transfer_wrap(p->l,k,p->par,p->tr1);
  if(!strcmp(p->tr1,p->tr2))
    d2=d1;
  else
    d2=transfer_wrap(p->l,k,p->par,p->tr2);

  return k*d1*d2*pk;
}

static double spectra(char *tr1,char *tr2,int l,RunParams *par)
{
  IntPar ipar;
  double result=0,eresult;
  gsl_function F;
  gsl_integration_workspace *w=gsl_integration_workspace_alloc(1000);
  ipar.l=l;
  ipar.par=par;
  ipar.tr1=tr1;
  ipar.tr2=tr2;
  F.function=&cl_integrand;
  F.params=&ipar;
  gsl_integration_qag(&F,D_LKMIN,D_LKMAX,0,1E-4,1000,GSL_INTEG_GAUSS41,w,&result,&eresult);
  gsl_integration_workspace_free(w);

  return M_LN10*2*result/(2*l+1.);
}

void compute_spectra(RunParams *par)
{
  printf("Computing power spectra\n");
#ifdef _HAS_OMP
#pragma omp parallel default(none) shared(par)
  {
#endif //_HAS_OMP
    int l;
#ifdef _HAS_OMP
#pragma omp for
#endif //_HAS_OMP
    for(l=0;l<=par->lmax;l++) {
      printf("%d \n",l);
      if(par->do_nc) {
	par->cl_dd[l]=spectra("nc","nc",l,par);
	if(par->do_shear)
	  par->cl_dl[l]=spectra("nc","shear",l,par);
	if(par->do_cmblens)
	  par->cl_dc[l]=spectra("nc","cmblens",l,par);
      }
      if(par->do_shear) {
	par->cl_ll[l]=spectra("shear","shear",l,par);
	if(par->do_cmblens)
	  par->cl_lc[l]=spectra("shear","cmblens",l,par);
      }
      if(par->do_cmblens) {
	par->cl_cc[l]=spectra("cmblens","cmblens",l,par);
      }
    } //end omp for
#ifdef _HAS_OMP
  } //end omp parallel
#endif //_HAS_OMP
}
