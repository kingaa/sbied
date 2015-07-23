/* pomp model file: pomp2E3B53FC7CC8 */

#include <pomp.h>
#include <R_ext/Rdynload.h>

 
#define b1	(__p[__parindex[0]])
#define b2	(__p[__parindex[1]])
#define b3	(__p[__parindex[2]])
#define b4	(__p[__parindex[3]])
#define b5	(__p[__parindex[4]])
#define b6	(__p[__parindex[5]])
#define psi	(__p[__parindex[6]])
#define rho	(__p[__parindex[7]])
#define tau	(__p[__parindex[8]])
#define sigma_dem	(__p[__parindex[9]])
#define sigma_env	(__p[__parindex[10]])
#define SO_0	(__p[__parindex[11]])
#define IO_0	(__p[__parindex[12]])
#define delta	(__p[__parindex[13]])
#define K	(__p[__parindex[14]])
#define SB1_0	(__p[__parindex[15]])
#define SB2_0	(__p[__parindex[16]])
#define SB3_0	(__p[__parindex[17]])
#define SB4_0	(__p[__parindex[18]])
#define SB5_0	(__p[__parindex[19]])
#define SB6_0	(__p[__parindex[20]])
#define SB1	(__x[__stateindex[0]])
#define SB2	(__x[__stateindex[1]])
#define SB3	(__x[__stateindex[2]])
#define SB4	(__x[__stateindex[3]])
#define SB5	(__x[__stateindex[4]])
#define SB6	(__x[__stateindex[5]])
#define IB	(__x[__stateindex[6]])
#define SO	(__x[__stateindex[7]])
#define IO	(__x[__stateindex[8]])
#define xi1	(__covars[__covindex[0]])
#define B	(__covars[__covindex[1]])
#define P	(__covars[__covindex[2]])
#define cases	(__y[__obsindex[0]])
#define DSB1	(__f[__stateindex[0]])
#define DSB2	(__f[__stateindex[1]])
#define DSB3	(__f[__stateindex[2]])
#define DSB4	(__f[__stateindex[3]])
#define DSB5	(__f[__stateindex[4]])
#define DSB6	(__f[__stateindex[5]])
#define DIB	(__f[__stateindex[6]])
#define DSO	(__f[__stateindex[7]])
#define DIO	(__f[__stateindex[8]])
#define Tb1	(__pt[__parindex[0]])
#define Tb2	(__pt[__parindex[1]])
#define Tb3	(__pt[__parindex[2]])
#define Tb4	(__pt[__parindex[3]])
#define Tb5	(__pt[__parindex[4]])
#define Tb6	(__pt[__parindex[5]])
#define Tpsi	(__pt[__parindex[6]])
#define Trho	(__pt[__parindex[7]])
#define Ttau	(__pt[__parindex[8]])
#define Tsigma_dem	(__pt[__parindex[9]])
#define Tsigma_env	(__pt[__parindex[10]])
#define TSO_0	(__pt[__parindex[11]])
#define TIO_0	(__pt[__parindex[12]])
#define Tdelta	(__pt[__parindex[13]])
#define TK	(__pt[__parindex[14]])
#define TSB1_0	(__pt[__parindex[15]])
#define TSB2_0	(__pt[__parindex[16]])
#define TSB3_0	(__pt[__parindex[17]])
#define TSB4_0	(__pt[__parindex[18]])
#define TSB5_0	(__pt[__parindex[19]])
#define TSB6_0	(__pt[__parindex[20]])
#define lik	(__lik[0])

void __pomp_initializer (double *__x, const double *__p, double t, const int *__stateindex, const int *__parindex, const int *__covindex, const double *__covars)
{

  SB1 = SB1_0;
  SB2 = SB2_0;
  SB3 = SB3_0;
  SB4 = SB4_0;
  SB5 = SB5_0;
  SB6 = SB6_0;
  IB = 0;
  IO = IO_0 * P;
  SO = SO_0 * P;
 
}


void __pomp_par_trans (double *__pt, const double *__p, const int *__parindex)
{

 Tpsi = exp(psi);
 Trho = expit(rho);
 Ttau = exp(tau);
 Tsigma_dem = exp(sigma_dem);
 Tsigma_env = exp(sigma_env);
 TSO_0 =  expit(SO_0);
 TIO_0 = expit(IO_0);
 
}


void __pomp_par_untrans (double *__pt, const double *__p, const int *__parindex)
{

 Tpsi = log(psi);
 Trho = logit(rho);
 Ttau = log(tau);
 Tsigma_dem = log(sigma_dem);
 Tsigma_env = log(sigma_env);
 TSO_0 =  logit(SO_0);
 TIO_0 = logit(IO_0);
 
}


void __pomp_rmeasure (double *__y, const double *__x, const double *__p, const int *__obsindex, const int *__stateindex, const int *__parindex, const int *__covindex, int __ncovars, const double *__covars, double t)
{
 
  cases = rnorm(rho*IO, sqrt( pow(tau*IO,2) + rho*IO ) );
  if (cases > 0.0) {
    cases = nearbyint(cases);
  } else {
    cases = 0.0;
  }
 
}


void __pomp_dmeasure (double *__lik, const double *__y, const double *__x, const double *__p, int give_log, const int *__obsindex, const int *__stateindex, const int *__parindex, const int *__covindex, int __ncovars, const double *__covars, double t)
{
 
  double tol = 1.0e-25;
  double mean_cases = rho*IO;
  double sd_cases = sqrt(pow(tau*IO,2) + mean_cases);
  if(cases > 0.0){
    lik = pnorm(cases+0.5,mean_cases,sd_cases,1,0) - pnorm(cases-0.5,mean_cases,sd_cases,1,0) + tol; 
  } else{
    lik = pnorm(cases+0.5,mean_cases,sd_cases,1,0) + tol;
  }
  if (give_log) lik = log(lik);
 
}


void __pomp_stepfn (double *__x, const double *__p, const int *__stateindex, const int *__parindex, const int *__covindex, int __covdim, const double *__covars, double t, double dt)
{

  double lambda, beta, var_epsilon, p, q;
 
  beta = exp(dot_product( (int) K, &xi1, &b1));
  lambda = (beta * (IO+IB) / P + psi);
  var_epsilon = pow(sigma_dem,2)/ lambda +  pow(sigma_env,2);
  lambda *= (var_epsilon < 1.0e-6) ? 1 : rgamma(1/var_epsilon,var_epsilon);
  p = exp(- (delta+lambda)/12);
  q = (1-p)*lambda/(delta+lambda);
  SB1 = B;
  SB2= SB1*p;
  SB3=SB2*p;
  SB4=SB3*p;
  SB5=SB4*p;
  SB6=SB5*p;
  SO= (SB6+SO)*p;
  IB=(SB1+SB2+SB3+SB4+SB5+SB6)*q;
  IO=SO*q;
 
}

#undef b1
#undef b2
#undef b3
#undef b4
#undef b5
#undef b6
#undef psi
#undef rho
#undef tau
#undef sigma_dem
#undef sigma_env
#undef SO_0
#undef IO_0
#undef delta
#undef K
#undef SB1_0
#undef SB2_0
#undef SB3_0
#undef SB4_0
#undef SB5_0
#undef SB6_0
#undef SB1
#undef SB2
#undef SB3
#undef SB4
#undef SB5
#undef SB6
#undef IB
#undef SO
#undef IO
#undef xi1
#undef B
#undef P
#undef cases
#undef DSB1
#undef DSB2
#undef DSB3
#undef DSB4
#undef DSB5
#undef DSB6
#undef DIB
#undef DSO
#undef DIO
#undef Tb1
#undef Tb2
#undef Tb3
#undef Tb4
#undef Tb5
#undef Tb6
#undef Tpsi
#undef Trho
#undef Ttau
#undef Tsigma_dem
#undef Tsigma_env
#undef TSO_0
#undef TIO_0
#undef Tdelta
#undef TK
#undef TSB1_0
#undef TSB2_0
#undef TSB3_0
#undef TSB4_0
#undef TSB5_0
#undef TSB6_0

static int __pomp_load_stack = 0;

void __pomp_load_stack_incr (void) {
  ++__pomp_load_stack;
}

void __pomp_load_stack_decr (int *val) {
  *val = --__pomp_load_stack;
}

void R_init_pomp2E3B53FC7CC8 (DllInfo *info)
{
R_RegisterCCallable("pomp2E3B53FC7CC8", "__pomp_load_stack_incr", (DL_FUNC) __pomp_load_stack_incr);
R_RegisterCCallable("pomp2E3B53FC7CC8", "__pomp_load_stack_decr", (DL_FUNC) __pomp_load_stack_decr);
R_RegisterCCallable("pomp2E3B53FC7CC8", "__pomp_initializer", (DL_FUNC) __pomp_initializer);
R_RegisterCCallable("pomp2E3B53FC7CC8", "__pomp_par_trans", (DL_FUNC) __pomp_par_trans);
R_RegisterCCallable("pomp2E3B53FC7CC8", "__pomp_par_untrans", (DL_FUNC) __pomp_par_untrans);
R_RegisterCCallable("pomp2E3B53FC7CC8", "__pomp_rmeasure", (DL_FUNC) __pomp_rmeasure);
R_RegisterCCallable("pomp2E3B53FC7CC8", "__pomp_dmeasure", (DL_FUNC) __pomp_dmeasure);
R_RegisterCCallable("pomp2E3B53FC7CC8", "__pomp_stepfn", (DL_FUNC) __pomp_stepfn);
}

