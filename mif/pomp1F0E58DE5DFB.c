/* pomp model file: pomp1F0E58DE5DFB */

#include <pomp.h>
#include <R_ext/Rdynload.h>

 
#define Beta	(__p[__parindex[0]])
#define mu_I	(__p[__parindex[1]])
#define rho	(__p[__parindex[2]])
#define mu_R1	(__p[__parindex[3]])
#define mu_R2	(__p[__parindex[4]])
#define S	(__x[__stateindex[0]])
#define I	(__x[__stateindex[1]])
#define R1	(__x[__stateindex[2]])
#define R2	(__x[__stateindex[3]])
#define B	(__y[__obsindex[0]])
#define C	(__y[__obsindex[1]])
#define DS	(__f[__stateindex[0]])
#define DI	(__f[__stateindex[1]])
#define DR1	(__f[__stateindex[2]])
#define DR2	(__f[__stateindex[3]])
#define TBeta	(__pt[__parindex[0]])
#define Tmu_I	(__pt[__parindex[1]])
#define Trho	(__pt[__parindex[2]])
#define Tmu_R1	(__pt[__parindex[3]])
#define Tmu_R2	(__pt[__parindex[4]])
#define lik	(__lik[0])

void __pomp_initializer (double *__x, const double *__p, double t, const int *__stateindex, const int *__parindex, const int *__covindex, const double *__covars)
{

 S=762;
 I=1;
 R1=0;
 R2=0;
 
}


void __pomp_par_trans (double *__pt, const double *__p, const int *__parindex)
{

 TBeta = exp(Beta);
 Tmu_I = exp(mu_I);
 Trho = expit(rho);
 
}


void __pomp_par_untrans (double *__pt, const double *__p, const int *__parindex)
{

 TBeta = log(Beta);
 Tmu_I = log(mu_I);
 Trho = logit(rho);
 
}


void __pomp_rmeasure (double *__y, const double *__x, const double *__p, const int *__obsindex, const int *__stateindex, const int *__parindex, const int *__covindex, int __ncovars, const double *__covars, double t)
{
 
  B = rpois(rho*R1+1e-6);
  C = rpois(rho*R2);
 
}


void __pomp_dmeasure (double *__lik, const double *__y, const double *__x, const double *__p, int give_log, const int *__obsindex, const int *__stateindex, const int *__parindex, const int *__covindex, int __ncovars, const double *__covars, double t)
{
 
  lik = dpois(B,rho*R1+1e-6,give_log);
 
}


void __pomp_stepfn (double *__x, const double *__p, const int *__stateindex, const int *__parindex, const int *__covindex, int __covdim, const double *__covars, double t, double dt)
{

  double t1 = rbinom(S,1-exp(-Beta*I*dt));
  double t2 = rbinom(I,1-exp(-dt*mu_I));
  double t3 = rbinom(R1,1-exp(-dt*mu_R1));
  double t4 = rbinom(R2,1-exp(-dt*mu_R2));
  S -= t1;
  I += t1 - t2;
  R1 += t2 - t3;
  R2 += t3 - t4;
 
}

#undef Beta
#undef mu_I
#undef rho
#undef mu_R1
#undef mu_R2
#undef S
#undef I
#undef R1
#undef R2
#undef B
#undef C
#undef DS
#undef DI
#undef DR1
#undef DR2
#undef TBeta
#undef Tmu_I
#undef Trho
#undef Tmu_R1
#undef Tmu_R2

static int __pomp_load_stack = 0;

void __pomp_load_stack_incr (void) {
  ++__pomp_load_stack;
}

void __pomp_load_stack_decr (int *val) {
  *val = --__pomp_load_stack;
}

void R_init_pomp1F0E58DE5DFB (DllInfo *info)
{
R_RegisterCCallable("pomp1F0E58DE5DFB", "__pomp_load_stack_incr", (DL_FUNC) __pomp_load_stack_incr);
R_RegisterCCallable("pomp1F0E58DE5DFB", "__pomp_load_stack_decr", (DL_FUNC) __pomp_load_stack_decr);
R_RegisterCCallable("pomp1F0E58DE5DFB", "__pomp_initializer", (DL_FUNC) __pomp_initializer);
R_RegisterCCallable("pomp1F0E58DE5DFB", "__pomp_par_trans", (DL_FUNC) __pomp_par_trans);
R_RegisterCCallable("pomp1F0E58DE5DFB", "__pomp_par_untrans", (DL_FUNC) __pomp_par_untrans);
R_RegisterCCallable("pomp1F0E58DE5DFB", "__pomp_rmeasure", (DL_FUNC) __pomp_rmeasure);
R_RegisterCCallable("pomp1F0E58DE5DFB", "__pomp_dmeasure", (DL_FUNC) __pomp_dmeasure);
R_RegisterCCallable("pomp1F0E58DE5DFB", "__pomp_stepfn", (DL_FUNC) __pomp_stepfn);
}

