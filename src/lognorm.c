#include "common.h"
#include "fftlog.h"

#define NZ_LGN 200
#define NK_LGN 4096
void compute_lognorm_bias(RunParams* par) {
  printf ("Computing lognormal transformation bias\n");
  double zmax=par->wind_0[0]->xf;
  double zmin=0;
  double lkmin=log(1e-5);
  double lkmax=log(1000.);
  double *za = malloc(NZ_LGN * sizeof(double));
  double *lka = malloc(NK_LGN * sizeof(double));
  double *ka = malloc(NK_LGN * sizeof(double));
  double *ra = malloc(NK_LGN * sizeof(double));
  double *Pk = malloc(NK_LGN * sizeof(double));
  double *Pkt = malloc(NK_LGN * sizeof(double));
  double *xi = malloc(NK_LGN * sizeof(double));
  double *ba = malloc(NZ_LGN*NK_LGN* sizeof(double));
    
  for (int i=0;i<NZ_LGN; i++)
    za[i]=zmin+1.0*i/(NZ_LGN-1)*(zmax-zmin);
  for (int i=0;i<NK_LGN; i++) {
    lka[i]=lkmin+(i+0.5)/NK_LGN*(lkmax-lkmin);
    ka[i]=exp(lka[i]);
  }

  double g0=csm_growth_factor(par->cpar,1.0);
  double rsm2=par->r_smooth*par->r_smooth;
  if (rsm2==0) {
    printf ("Cannot do log predictions withmout smoothing!");
    exit(1);
  }
  // now do the actual calculation
  for (int i=0; i<NZ_LGN; i++) {
    double curz=za[i];
    double gf=csm_growth_factor(par->cpar,curz)/g0;
    double bias=spline_eval(curz,par->bias);
    for (int j=0; j<NK_LGN; j++) 
      Pk[j]=csm_Pk_linear_0(par->cpar,ka[j])*gf*gf*exp(-ka[j]*ka[j]*rsm2);
    pk2xi(NK_LGN, ka, Pk, ra, xi);
    for (int j=0; j<NK_LGN; j++) 
      xi[j]=exp(xi[j]*bias*bias)-1.0;
    xi2pk(NK_LGN, ra,xi,ka,Pkt);
    // now we have the effective extra bias
    for (int j=0; j<NK_LGN; j++) {
      double val=Pkt[j]/Pk[j]/bias/bias;
      if ((val>1000) || (val<0) || isnan(val)) val=1.0;
      ba[j*NZ_LGN+i]=sqrt(val);
    }
  }

  /*  for (int i=0;i<NZ_LGN; i++)
    printf ("%g ",za[i]);
  printf ("aaa\n");
  for (int i=0;i<NK_LGN; i++) 
    printf ("%g ",lka[i]);
   printf ("aaa\n");
  for (int i=0; i<NZ_LGN*NK_LGN; i++)
  if (ba[i]==0) printf ("^^^^^^^^^^^^^^^^^^^^^ %i",i); */
   par->lognorm_bias=spline2D_init(NZ_LGN, NK_LGN, za, lka, ba, 0., 100.);
}
