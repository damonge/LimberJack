#include "common.h"

typedef struct {
  int l;
  RunParams *par;
  char *tr1;
  char *tr2;
} IntClPar;

static double cl_integrand(double lk,void *params)
{
  double d1,d2;
  IntClPar *p=(IntClPar *)params;
  double k=pow(10.,lk);
  double pk=csm_Pk_linear_0(p->par->cpar,k);
  if(p->par->r_smooth>0)
    pk*=exp(-k*k*p->par->r_smooth*p->par->r_smooth);
  d1=transfer_wrap(p->l,k,p->par,p->tr1,0);
  d2=transfer_wrap(p->l,k,p->par,p->tr2,1);

  return k*k*k*d1*d2*pk;
}

static void get_k_interval(RunParams *par,char *tr1,char *tr2,int l,double *lkmin,double *lkmax)
{
  double chimin,chimax;
  if(l<par->l_limber_min) {
    chimin=2*(l+0.5)/pow(10.,D_LKMAX);
    chimax=0.5*(l+0.5)/pow(10.,D_LKMIN);
  }
  else {
    if(!strcmp(tr1,"nc")) {
      if(!strcmp(tr2,"nc")) {
	chimin=fmax(par->chimin_nc[0],par->chimin_nc[1]);
	chimax=fmin(par->chimax_nc[0],par->chimax_nc[1]);
      }
      else {
	chimin=par->chimin_nc[0];
	chimax=par->chimax_nc[0];
      }
    }
    else if(!strcmp(tr2,"nc")) {
      chimin=par->chimin_nc[1];
      chimax=par->chimax_nc[1];
    }
    else {
      chimin=0.5*(l+0.5)/pow(10.,D_LKMAX);
      chimax=2*(l+0.5)/pow(10.,D_LKMIN);
    }
  }    

  if(chimin<=0)
    chimin=0.5*(l+0.5)/pow(10.,D_LKMAX);

  *lkmax=fmin(D_LKMAX,log10(2*(l+0.5)/chimin));
  *lkmin=fmax(D_LKMIN,log10(0.5*(l+0.5)/chimax));
}

static double spectra(char *tr1,char *tr2,int l,RunParams *par)
{
  IntClPar ipar;
  double result=0,eresult;
  double lkmax,lkmin;
  gsl_function F;
  gsl_integration_workspace *w=gsl_integration_workspace_alloc(1000);
  ipar.l=l;
  ipar.par=par;
  ipar.tr1=tr1;
  ipar.tr2=tr2;
  F.function=&cl_integrand;
  F.params=&ipar;
  get_k_interval(par,tr1,tr2,l,&lkmin,&lkmax);
  gsl_integration_qag(&F,lkmin,lkmax,0,1E-4,1000,GSL_INTEG_GAUSS41,w,&result,&eresult);
  gsl_integration_workspace_free(w);

  return result*M_LN10*2./M_PI;
}

void compute_spectra(RunParams *par)
{
  double *l_sampling=my_malloc(par->n_ell*sizeof(double));

  printf("Computing power spectra\n");
#ifdef _HAS_OMP
#pragma omp parallel default(none) shared(par,l_sampling)
  {
#endif //_HAS_OMP
    int ii;
#ifdef _HAS_OMP
#pragma omp for schedule(dynamic)
#endif //_HAS_OMP
    for(ii=0;ii<par->n_ell;ii++) {
      int l=par->l_sampling_arr[ii];
#ifdef _DEBUG
      printf("%d \n",l);
#endif //_DEBUG
      l_sampling[ii]=l;
      if(par->do_nc) {
	par->cl_dd[ii]=spectra("nc","nc",l,par);
	if(par->do_shear) {
	  par->cl_d1l2[ii]=spectra("nc","shear",l,par);
	  par->cl_d2l1[ii]=spectra("shear","nc",l,par);
	}
	if(par->do_cmblens)
	  par->cl_dc[ii]=spectra("nc","cmblens",l,par);
	if(par->do_isw)
	  par->cl_di[ii]=spectra("nc","isw",l,par);
      }
      if(par->do_shear) {
	par->cl_ll[ii]=spectra("shear","shear",l,par);
	if(par->do_cmblens)
	  par->cl_lc[ii]=spectra("shear","cmblens",l,par);
	if(par->do_isw)
	  par->cl_li[ii]=spectra("shear","isw",l,par);
      }
      if(par->do_cmblens) {
	par->cl_cc[ii]=spectra("cmblens","cmblens",l,par);
	if(par->do_isw)
	  par->cl_ci[ii]=spectra("cmblens","isw",l,par);
      }
      if(par->do_isw)
	par->cl_ii[ii]=spectra("isw","isw",l,par);
    } //end omp for
#ifdef _HAS_OMP
  } //end omp parallel
#endif //_HAS_OMP

  //Initialize splines
  if(par->do_nc) {
    par->cl_dd_s=spline_init(par->n_ell,l_sampling,par->cl_dd,0,0);
    if(par->do_shear) {
      par->cl_d1l2_s=spline_init(par->n_ell,l_sampling,par->cl_d1l2,0,0);
      par->cl_d2l1_s=spline_init(par->n_ell,l_sampling,par->cl_d2l1,0,0);
    }
    if(par->do_cmblens)
      par->cl_dc_s=spline_init(par->n_ell,l_sampling,par->cl_dc,0,0);
    if(par->do_isw)
      par->cl_di_s=spline_init(par->n_ell,l_sampling,par->cl_di,0,0);
  }
  if(par->do_shear) {
    par->cl_ll_s=spline_init(par->n_ell,l_sampling,par->cl_ll,0,0);
    if(par->do_cmblens)
      par->cl_lc_s=spline_init(par->n_ell,l_sampling,par->cl_lc,0,0);
    if(par->do_isw)
      par->cl_li_s=spline_init(par->n_ell,l_sampling,par->cl_li,0,0);
  }
  if(par->do_cmblens) {
    par->cl_cc_s=spline_init(par->n_ell,l_sampling,par->cl_cc,0,0);
    if(par->do_isw)
      par->cl_ci_s=spline_init(par->n_ell,l_sampling,par->cl_ci,0,0);
  }
  if(par->do_isw)
    par->cl_ii_s=spline_init(par->n_ell,l_sampling,par->cl_ii,0,0);
  
  //Sample into arrays of delta_ell=1
  int ell;
  for(ell=0;ell<=par->lmax;ell++) {
    if(par->do_nc) {
      par->cl_dd[ell]=spline_eval((double)ell,par->cl_dd_s);
      if(par->do_shear) {
	par->cl_d1l2[ell]=spline_eval((double)ell,par->cl_d1l2_s);
	par->cl_d2l1[ell]=spline_eval((double)ell,par->cl_d2l1_s);
      }
      if(par->do_cmblens)
	par->cl_dc[ell]=spline_eval((double)ell,par->cl_dc_s);
      if(par->do_isw)
	par->cl_di[ell]=spline_eval((double)ell,par->cl_di_s);
    }
    if(par->do_shear) {
      par->cl_ll[ell]=spline_eval((double)ell,par->cl_ll_s);
      if(par->do_cmblens)
	par->cl_lc[ell]=spline_eval((double)ell,par->cl_lc_s);
      if(par->do_isw)
	par->cl_li[ell]=spline_eval((double)ell,par->cl_li_s);
    }
    if(par->do_cmblens) {
      par->cl_cc[ell]=spline_eval((double)ell,par->cl_cc_s);
      if(par->do_isw)
	par->cl_ci[ell]=spline_eval((double)ell,par->cl_ci_s);
    }
    if(par->do_isw)
      par->cl_ii[ell]=spline_eval((double)ell,par->cl_ii_s);
  }
}

typedef struct {
  int i_bessel;
  RunParams *par;
  double th;
  SplPar *clsp;
  double *cl;
} IntWtPar;

static double wt_integrand(double l,void *params)
{
  IntWtPar *p=(IntWtPar *)params;
  double x=l*p->th;
  //  double cl=spline_eval(l,p->clsp);
  double cl=p->cl[(int)l];
  double jbes;

  if(p->i_bessel)
    jbes=gsl_sf_bessel_Jn(p->i_bessel,x);
  else
    jbes=gsl_sf_bessel_J0(x);

  return l*jbes*cl;
}

static void compute_wt_single(RunParams *par,double *cl,double *wt,double *llist,int bessel_order)
{
#ifdef _HAS_OMP
#pragma omp parallel default(none)		\
  shared(par,cl,wt,llist,bessel_order)
  {
#endif //_HAS_OMP
    int ith;
    double result,eresult;
    gsl_function F;
    gsl_integration_workspace *w=gsl_integration_workspace_alloc(1000);
    SplPar *clsp=spline_init((par->lmax+1),llist,cl,0,0);
    IntWtPar ipar;

    ipar.i_bessel=bessel_order;
    ipar.par=par;
    ipar.clsp=clsp;
    ipar.cl=cl;
    
#ifdef _HAS_OMP
#pragma omp for
#endif //_HAS_OMP
    for(ith=0;ith<par->n_th;ith++) {
      if(par->do_w_theta_logbin)
	ipar.th=DTOR*par->th_max*pow(10.,(ith+0.5-par->n_th)/par->n_th_logint);
      else
	ipar.th=DTOR*(par->th_min+(par->th_max-par->th_min)*(ith+0.5)/par->n_th);
      F.function=&wt_integrand;
      F.params=&ipar;
      gsl_integration_qag(&F,llist[0],llist[par->lmax],0,1E-4,1000,GSL_INTEG_GAUSS41,w,&result,&eresult);
      wt[ith]=result/(2*M_PI);
    }//end omp for
    gsl_integration_workspace_free(w);
    spline_free(clsp);
#ifdef _HAS_OMP
  } //end omp parallel
#endif //_HAS_OMP
}

void compute_w_theta(RunParams *par)
{
  if(par->do_w_theta) {
    int l;
    double *llist;

    printf("Computing correlation functions\n");

    llist=my_malloc((par->lmax+1)*sizeof(double));
    for(l=0;l<=par->lmax;l++)
      llist[l]=(float)l;

    if(par->do_nc) {
      compute_wt_single(par,par->cl_dd,par->wt_dd,llist,0);
      if(par->do_shear) {
	compute_wt_single(par,par->cl_d1l2,par->wt_d1l2,llist,0);
	compute_wt_single(par,par->cl_d2l1,par->wt_d2l1,llist,0);
      }
      if(par->do_cmblens)
	compute_wt_single(par,par->cl_dc,par->wt_dc,llist,0);
      if(par->do_isw)
	compute_wt_single(par,par->cl_di,par->wt_di,llist,0);
    }
    if(par->do_shear) {
      compute_wt_single(par,par->cl_ll,par->wt_ll_pp,llist,0);
      compute_wt_single(par,par->cl_ll,par->wt_ll_mm,llist,4);
      if(par->do_cmblens)
	compute_wt_single(par,par->cl_lc,par->wt_lc,llist,0);
      if(par->do_isw)
	compute_wt_single(par,par->cl_li,par->wt_li,llist,0);
    }
    if(par->do_cmblens) {
      compute_wt_single(par,par->cl_cc,par->wt_cc,llist,0);
      if(par->do_isw)
	compute_wt_single(par,par->cl_ci,par->wt_ci,llist,0);
    }
    if(par->do_isw)
      compute_wt_single(par,par->cl_ii,par->wt_ii,llist,0);
    free(llist);
  }
  else {
    printf("Skipping correlation functions\n");
  }
}
