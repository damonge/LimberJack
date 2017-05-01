#include "common.h"

//TODO: check non-limber rsd
//TODO: check spin  prefactors (sqrt?)
//TODO: check lensing signs

#define N_CHI 1000

static double f_dens(double chi,RunParams *par,int ibin)
{
  double gf=spline_eval(chi,par->gfofchi);
  double h=spline_eval(chi,par->hofchi);
  double z=spline_eval(chi,par->zofchi);
  double pz=spline_eval(z,par->wind_0[ibin]);
  double bz=spline_eval(z,par->bias);

  return pz*bz*gf*h;
}

static double transfer_dens_nolim(int l,double k,RunParams *par,int ibin)
{
  double ret;
  if(l>par->l_limber_min)
    ret=f_dens((l+0.5)/k,par,ibin)*sqrt(M_PI/(2*l+1.))/k;
  else {
    int i;
    double dchi=(par->chimax_nc[ibin]-par->chimin_nc[ibin])/N_CHI;

    ret=0;
    for(i=0;i<N_CHI;i++) {
      double chi=par->chimin_nc[ibin]+(i+0.5)*dchi;
      double jl=csm_j_bessel(l,k*chi);
      
      ret+=f_dens(chi,par,ibin)*jl;
    }
    ret*=dchi;
  }

  return ret;
}

static double f_rsd(double chi,RunParams *par,int ibin)
{
  double z =spline_eval(chi,par->zofchi);
  double gf=spline_eval(chi,par->gfofchi);
  double fg=spline_eval(chi,par->fgofchi);
  double h =spline_eval(chi,par->hofchi);
  double pz=spline_eval(z  ,par->wind_0[ibin]);

  return pz*fg*gf*h;
}

static double transfer_rsd_nolim(int l,double k,RunParams *par,int ibin)
{
  double ret;

  if(l>par->l_limber_min) {
    double x0=(l+0.5),x1=(l+1.5);
    double chi0=x0/k,chi1=x1/k;
    double f0=f_rsd(chi0,par,ibin),f1=f_rsd(chi1,par,ibin);
    ret=(f0*(1.-l*(l-1.)/(x0*x0))*sqrt(M_PI/(2*l+1.))-f1*2*sqrt(M_PI/(2*l+3.))/x1)/k;
  }
  else {
    int i;
    double dchi=(par->chimax_nc[ibin]-par->chimin_nc[ibin])/N_CHI;

    ret=0;
    for(i=0;i<N_CHI;i++) {
      double chi=par->chimin_nc[ibin]+(i+0.5)*dchi;
      double x=k*chi;
      double jl=csm_j_bessel(l,x);
      double jlp1=csm_j_bessel(l+1,x);

      ret+=f_rsd(chi,par,ibin)*((x*x-l*(l-1))*jl-2*x*jlp1)/(x*x);
    }
    ret*=dchi;
  }

  return ret;
}

static double f_magnification(double chi,RunParams *par,int ibin)
{
  double gf=spline_eval(chi,par->gfofchi);
  double a=spline_eval(chi,par->aofchi);
  double w=spline_eval(chi,par->wind_M[ibin]);

  if(w<=0)
    return 0;
  else
    return gf*w/(a*chi);
}

static double transfer_magnification_nolim(int l,double k,RunParams *par,int ibin)
{
  double ret;
  if(l>par->l_limber_min)
    ret=f_magnification((l+0.5)/k,par,ibin)*sqrt(M_PI/(2*l+1.))/k;
  else {
    int i;
    double dchi=par->chimax_nc[ibin]/N_CHI;

    ret=0;
    for(i=0;i<N_CHI;i++) {
      double chi=(i+0.5)*dchi;
      double jl=csm_j_bessel(l,k*chi);
    
      ret+=f_magnification(chi,par,ibin)*jl;
    }
    ret*=dchi;
  }

  return -2*par->prefac_lensing*l*(l+1)*ret/(k*k);
}

static double f_lensing(double chi,RunParams *par,int ibin)
{
  double gf=spline_eval(chi,par->gfofchi);
  double a=spline_eval(chi,par->aofchi);
  double w=spline_eval(chi,par->wind_L[ibin]);

  if(w<=0)
    return 0;
  else
    return gf*w/(a*chi);
}

static double transfer_lensing_nolim(int l,double k,RunParams *par,int ibin)
{
  double ret;
  if(l>par->l_limber_min)
    ret=f_lensing((l+0.5)/k,par,ibin)*sqrt(M_PI/(2*l+1.))/k;
  else {
    int i;
    double dchi=par->chimax_nc[ibin]/N_CHI;
    
    ret=0;
    for(i=0;i<N_CHI;i++) {
      double chi=(i+0.5)*dchi;
      double jl=csm_j_bessel(l,k*chi);
    
      ret+=f_lensing(chi,par,ibin)*jl;
    }
    ret*=dchi;
  }

  return par->prefac_lensing*sqrt((l+2.)*(l+1.)*l*(l-1.))*ret/(k*k);
}

static double f_IA_NLA(double chi,RunParams *par,int ibin)
{
  double gf=spline_eval(chi,par->gfofchi);
  double h=spline_eval(chi,par->hofchi);
  double z=spline_eval(chi,par->zofchi);
  double pz=spline_eval(z,par->wind_0[ibin]);
  double az=spline_eval(z,par->abias);

  return pz*az*gf*h/(chi*chi);
}

static double transfer_IA_NLA_nolim(int l,double k,RunParams *par,int ibin)
{
  double ret;
  if(l>par->l_limber_min)
    ret=f_IA_NLA((l+0.5)/k,par,ibin)*sqrt(M_PI/(2*l+1.))/k;
  else {
    int i;
    double dchi=(par->chimax_nc[ibin]-par->chimin_nc[ibin])/N_CHI;

    ret=0;
    for(i=0;i<N_CHI;i++) {
      double chi=par->chimin_nc[ibin]+(i+0.5)*dchi;
      double jl=csm_j_bessel(l,k*chi);

      ret+=f_IA_NLA(chi,par,ibin)*jl;
    }
    ret*=dchi;
  }

  return sqrt((l+2.)*(l+1.)*l*(l-1.))*ret/(k*k);
}

static double f_cmblens(double chi,RunParams *par)
{
  if(chi>=par->chi_kappa)
    return 0;
  else {
    double gf=spline_eval(chi,par->gfofchi);
    double a=spline_eval(chi,par->aofchi);
    double w=1-chi/par->chi_kappa;

    return gf*w/(a*chi);
  }
}

static double transfer_cmblens_nolim(int l,double k,RunParams *par)
{
  double ret;
  if(l>par->l_limber_min)
    ret=f_cmblens((l+0.5)/k,par)*sqrt(M_PI/(2*l+1.))/k;
  else {
    int i;
    double dchi=par->chi_kappa/N_CHI;

    ret=0;
    for(i=0;i<N_CHI;i++) {
      double chi=(i+0.5)*dchi;
      double jl=csm_j_bessel(l,k*chi);

      ret+=f_cmblens(chi,par)*jl;
    }
    ret*=dchi;
  }

  return par->prefac_lensing*l*(l+1)*ret/(k*k);
}

static double f_isw(double chi,RunParams *par)
{
  if(chi>=par->chi_isw)
    return 0;
  else {
    double gf=spline_eval(chi,par->gfofchi);
    double h=spline_eval(chi,par->hofchi);
    double fg=spline_eval(chi,par->fgofchi);

    return h*gf*(1-fg);
  }
}

static double transfer_isw_nolim(int l,double k,RunParams *par)
{
  double ret;
  if(l>par->l_limber_min)
    ret=f_isw((l+0.5)/k,par)*sqrt(M_PI/(2*l+1.))/k;
  else {
    int i;
    double dchi=par->chi_isw/N_CHI;
    
    ret=0;
    for(i=0;i<N_CHI;i++) {
      double chi=(i+0.5)*dchi;
      double jl=csm_j_bessel(l,k*chi);

      ret+=f_isw(chi,par)*jl;
    }
    ret*=dchi;
  }

  return 2*par->prefac_lensing*ret/(k*k);
}

double transfer_wrap(int l,double k,RunParams *par,char *trtype,int ibin)
{
  double tr=-1;
  if(!strcmp(trtype,"nc")) {
    if(par->do_nc!=1)
      report_error(1,"Asked to calculate NC transfer function, but can't!\n");
    else {
      tr=0;
      if(par->has_dens)
	tr+=transfer_dens_nolim(l,k,par,ibin);
      if(par->has_rsd)
	tr+=transfer_rsd_nolim(l,k,par,ibin);
      if(par->has_lensing)
	tr+=transfer_magnification_nolim(l,k,par,ibin);
    }
  }
  else if(!strcmp(trtype,"shear")) {
    if(par->do_shear!=1)
      report_error(1,"Asked to calculate shear transfer function, but can't!\n");
    else {
      tr=transfer_lensing_nolim(l,k,par,ibin);
      if(par->has_intrinsic_alignment)
	tr+=transfer_IA_NLA_nolim(l,k,par,ibin);
    }
  }
  else if(!strcmp(trtype,"cmblens")) {
    if(par->do_cmblens!=1)
      report_error(1,"Asked to calculate shear transfer function, but can't!\n");
    else
      tr=transfer_cmblens_nolim(l,k,par);
  }
  else if(!strcmp(trtype,"isw")) {
    if(par->do_isw!=1)
      report_error(1,"Asked to calculate ISW transfer function, but can't!\n");
    else
      tr=transfer_isw_nolim(l,k,par);
  }
  else {
    report_error(1,"Unknown transfer type %s\n",trtype);
    tr=-1.;
  }

  return tr;
}
