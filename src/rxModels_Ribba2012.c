#include <RxODE.h>
#include <RxODE_model.h>
#define __MAX_PROD__ 0
#define _CMT CMT
extern void  rxModels_Ribba2012_ode_solver_solvedata (rx_solve *solve){
  _solveData = solve;
}
extern rx_solve *rxModels_Ribba2012_ode_solver_get_solvedata(){
  return _solveData;
}
SEXP rxModels_Ribba2012_model_vars();
double _theta[17];
extern double* rxModels_Ribba2012_theta(double *theta){
  _theta[0] = 100.0000000000000000; _theta[1] = 0.2400000000000000; _theta[2] = 0.0000000000000000; _theta[3] = 0.0295000000000000; _theta[4] = 0.0000000000000000; _theta[5] = 0.0031000000000000; _theta[6] = 0.0000000000000000; _theta[7] = 0.1210000000000000; _theta[8] = 0.0000000000000000; _theta[9] = 0.7290000000000000; _theta[10] = 0.0000000000000000; _theta[11] = 0.0086700000000000; _theta[12] = 0.0000000000000000; _theta[13] = 7.1299999999999999; _theta[14] = 0.0000000000000000; _theta[15] = 41.2000000000000028; _theta[16] = 0.0000000000000000;
  return _theta;
}


// prj-specific differential eqns
void rxModels_Ribba2012_dydt(int *_neq, double t, double *__zzStateVar__, double *__DDtStateVar__)
{
  int _cSub = _neq[1];
  double   k,
  tkde,
  eta_DoT_tkde,
  kde,
  tkpq,
  eta_DoT_kpq,
  kpq,
  tkqpp,
  eta_DoT_kqpp,
  kqpp,
  tlambdap,
  eta_DoT_lambdap,
  lambdap,
  tgamma,
  eta_DoT_gamma,
  gamma,
  tdeltaqp,
  eta_DoT_deltaqp,
  deltaqp,
  pstar,
  pt,
  q,
  qp,
  c,
  tpt0,
  eta_DoT_pt0,
  pt0,
  tq0,
  eta_DoT_q0,
  q0;

  (void)t;
  (void)k;
  (void)tkde;
  (void)eta_DoT_tkde;
  (void)kde;
  (void)tkpq;
  (void)eta_DoT_kpq;
  (void)kpq;
  (void)tkqpp;
  (void)eta_DoT_kqpp;
  (void)kqpp;
  (void)tlambdap;
  (void)eta_DoT_lambdap;
  (void)lambdap;
  (void)tgamma;
  (void)eta_DoT_gamma;
  (void)gamma;
  (void)tdeltaqp;
  (void)eta_DoT_deltaqp;
  (void)deltaqp;
  (void)pstar;
  (void)pt;
  (void)q;
  (void)qp;
  (void)c;
  (void)tpt0;
  (void)eta_DoT_pt0;
  (void)pt0;
  (void)tq0;
  (void)eta_DoT_q0;
  (void)q0;

  _update_par_ptr(t, _cSub, _solveData, _idx);
  k = _PP[0];
  tkde = _PP[1];
  eta_DoT_tkde = _PP[2];
  tkpq = _PP[3];
  eta_DoT_kpq = _PP[4];
  tkqpp = _PP[5];
  eta_DoT_kqpp = _PP[6];
  tlambdap = _PP[7];
  eta_DoT_lambdap = _PP[8];
  tgamma = _PP[9];
  eta_DoT_gamma = _PP[10];
  tdeltaqp = _PP[11];
  eta_DoT_deltaqp = _PP[12];
  tpt0 = _PP[13];
  eta_DoT_pt0 = _PP[14];
  tq0 = _PP[15];
  eta_DoT_q0 = _PP[16];

  c = __zzStateVar__[0]*((double)(_ON[0]));
  pt = __zzStateVar__[1]*((double)(_ON[1]));
  q = __zzStateVar__[2]*((double)(_ON[2]));
  qp = __zzStateVar__[3]*((double)(_ON[3]));

  kde=tkde*exp(eta_DoT_tkde);
  kpq=tkpq*exp(eta_DoT_kpq);
  kqpp=tkqpp*exp(eta_DoT_kqpp);
  lambdap=tlambdap*exp(eta_DoT_lambdap);
  gamma=tgamma*exp(eta_DoT_gamma);
  deltaqp=tdeltaqp*exp(eta_DoT_deltaqp);
  pstar=pt+q+qp;
  __DDtStateVar__[0] = ((double)(_ON[0]))*(_IR[0] -kde*c);
  __DDtStateVar__[1] = ((double)(_ON[1]))*(_IR[1] + lambdap*pt*(1-pstar/safe_zero(k))+kqpp*qp-kpq*pt-gamma*c*kde*pt);
  __DDtStateVar__[2] = ((double)(_ON[2]))*(_IR[2] + kpq*pt-gamma*c*kde*q);
  __DDtStateVar__[3] = ((double)(_ON[3]))*(_IR[3] + gamma*c*kde*q-kqpp*qp-deltaqp*qp);
  pt0=tpt0*exp(eta_DoT_pt0);
  q0=tq0*exp(eta_DoT_q0);
  (&_solveData->subjects[_cSub])->dadt_counter[0]++;
}

// Jacobian derived vars
void rxModels_Ribba2012_calc_jac(int *_neq, double t, double *__zzStateVar__, double *__PDStateVar__, unsigned int __NROWPD__) {
  int _cSub=_neq[1];
  (&_solveData->subjects[_cSub])->jac_counter[0]++;
}
// Functional based initial conditions.
void rxModels_Ribba2012_inis(int _cSub, double *__zzStateVar__){
  double t=0;
  double   k,
  tkde,
  eta_DoT_tkde,
  kde,
  tkpq,
  eta_DoT_kpq,
  kpq,
  tkqpp,
  eta_DoT_kqpp,
  kqpp,
  tlambdap,
  eta_DoT_lambdap,
  lambdap,
  tgamma,
  eta_DoT_gamma,
  gamma,
  tdeltaqp,
  eta_DoT_deltaqp,
  deltaqp,
  pstar,
  pt,
  q,
  qp,
  c,
  tpt0,
  eta_DoT_pt0,
  pt0,
  tq0,
  eta_DoT_q0,
  q0;

  (void)t;
  (void)k;
  (void)tkde;
  (void)eta_DoT_tkde;
  (void)kde;
  (void)tkpq;
  (void)eta_DoT_kpq;
  (void)kpq;
  (void)tkqpp;
  (void)eta_DoT_kqpp;
  (void)kqpp;
  (void)tlambdap;
  (void)eta_DoT_lambdap;
  (void)lambdap;
  (void)tgamma;
  (void)eta_DoT_gamma;
  (void)gamma;
  (void)tdeltaqp;
  (void)eta_DoT_deltaqp;
  (void)deltaqp;
  (void)pstar;
  (void)pt;
  (void)q;
  (void)qp;
  (void)c;
  (void)tpt0;
  (void)eta_DoT_pt0;
  (void)pt0;
  (void)tq0;
  (void)eta_DoT_q0;
  (void)q0;

  _update_par_ptr(0.0, _cSub, _solveData, _idx);
  k = _PP[0];
  tkde = _PP[1];
  eta_DoT_tkde = _PP[2];
  tkpq = _PP[3];
  eta_DoT_kpq = _PP[4];
  tkqpp = _PP[5];
  eta_DoT_kqpp = _PP[6];
  tlambdap = _PP[7];
  eta_DoT_lambdap = _PP[8];
  tgamma = _PP[9];
  eta_DoT_gamma = _PP[10];
  tdeltaqp = _PP[11];
  eta_DoT_deltaqp = _PP[12];
  tpt0 = _PP[13];
  eta_DoT_pt0 = _PP[14];
  tq0 = _PP[15];
  eta_DoT_q0 = _PP[16];

  c = __zzStateVar__[0]*((double)(_ON[0]));
  pt = __zzStateVar__[1]*((double)(_ON[1]));
  q = __zzStateVar__[2]*((double)(_ON[2]));
  qp = __zzStateVar__[3]*((double)(_ON[3]));

  kde=tkde*exp(eta_DoT_tkde);
  kpq=tkpq*exp(eta_DoT_kpq);
  kqpp=tkqpp*exp(eta_DoT_kqpp);
  lambdap=tlambdap*exp(eta_DoT_lambdap);
  gamma=tgamma*exp(eta_DoT_gamma);
  deltaqp=tdeltaqp*exp(eta_DoT_deltaqp);
  pstar=pt+q+qp;
  pt0=tpt0*exp(eta_DoT_pt0);
  q0=tq0*exp(eta_DoT_q0);
  pt=pt0;
  q=q0;
  __zzStateVar__[0]=((double)(_ON[0]))*(c);
  __zzStateVar__[1]=((double)(_ON[1]))*(pt);
  __zzStateVar__[2]=((double)(_ON[2]))*(q);
  __zzStateVar__[3]=((double)(_ON[3]))*(qp);
}
// prj-specific derived vars
void rxModels_Ribba2012_calc_lhs(int _cSub, double t, double *__zzStateVar__, double *_lhs) {
  double   __DDtStateVar_0__,
  __DDtStateVar_1__,
  __DDtStateVar_2__,
  __DDtStateVar_3__,
  k,
  tkde,
  eta_DoT_tkde,
  kde,
  tkpq,
  eta_DoT_kpq,
  kpq,
  tkqpp,
  eta_DoT_kqpp,
  kqpp,
  tlambdap,
  eta_DoT_lambdap,
  lambdap,
  tgamma,
  eta_DoT_gamma,
  gamma,
  tdeltaqp,
  eta_DoT_deltaqp,
  deltaqp,
  pstar,
  pt,
  q,
  qp,
  c,
  tpt0,
  eta_DoT_pt0,
  pt0,
  tq0,
  eta_DoT_q0,
  q0;

  (void)t;
  (void)__DDtStateVar_0__;
  (void)__DDtStateVar_1__;
  (void)__DDtStateVar_2__;
  (void)__DDtStateVar_3__;
  (void)k;
  (void)tkde;
  (void)eta_DoT_tkde;
  (void)kde;
  (void)tkpq;
  (void)eta_DoT_kpq;
  (void)kpq;
  (void)tkqpp;
  (void)eta_DoT_kqpp;
  (void)kqpp;
  (void)tlambdap;
  (void)eta_DoT_lambdap;
  (void)lambdap;
  (void)tgamma;
  (void)eta_DoT_gamma;
  (void)gamma;
  (void)tdeltaqp;
  (void)eta_DoT_deltaqp;
  (void)deltaqp;
  (void)pstar;
  (void)pt;
  (void)q;
  (void)qp;
  (void)c;
  (void)tpt0;
  (void)eta_DoT_pt0;
  (void)pt0;
  (void)tq0;
  (void)eta_DoT_q0;
  (void)q0;

  _update_par_ptr(t, _cSub, _solveData, _idx);
  k = _PP[0];
  tkde = _PP[1];
  eta_DoT_tkde = _PP[2];
  tkpq = _PP[3];
  eta_DoT_kpq = _PP[4];
  tkqpp = _PP[5];
  eta_DoT_kqpp = _PP[6];
  tlambdap = _PP[7];
  eta_DoT_lambdap = _PP[8];
  tgamma = _PP[9];
  eta_DoT_gamma = _PP[10];
  tdeltaqp = _PP[11];
  eta_DoT_deltaqp = _PP[12];
  tpt0 = _PP[13];
  eta_DoT_pt0 = _PP[14];
  tq0 = _PP[15];
  eta_DoT_q0 = _PP[16];

  c = __zzStateVar__[0]*((double)(_ON[0]));
  pt = __zzStateVar__[1]*((double)(_ON[1]));
  q = __zzStateVar__[2]*((double)(_ON[2]));
  qp = __zzStateVar__[3]*((double)(_ON[3]));

  kde=tkde*exp(eta_DoT_tkde);
  kpq=tkpq*exp(eta_DoT_kpq);
  kqpp=tkqpp*exp(eta_DoT_kqpp);
  lambdap=tlambdap*exp(eta_DoT_lambdap);
  gamma=tgamma*exp(eta_DoT_gamma);
  deltaqp=tdeltaqp*exp(eta_DoT_deltaqp);
  pstar=pt+q+qp;
  __DDtStateVar_0__ = ((double)(_ON[0]))*(_IR[0] -kde*c);
  __DDtStateVar_1__ = ((double)(_ON[1]))*(_IR[1] + lambdap*pt*(1-pstar/safe_zero(k))+kqpp*qp-kpq*pt-gamma*c*kde*pt);
  __DDtStateVar_2__ = ((double)(_ON[2]))*(_IR[2] + kpq*pt-gamma*c*kde*q);
  __DDtStateVar_3__ = ((double)(_ON[3]))*(_IR[3] + gamma*c*kde*q-kqpp*qp-deltaqp*qp);
  pt0=tpt0*exp(eta_DoT_pt0);
  q0=tq0*exp(eta_DoT_q0);

  _lhs[0]=kde;
  _lhs[1]=kpq;
  _lhs[2]=kqpp;
  _lhs[3]=lambdap;
  _lhs[4]=gamma;
  _lhs[5]=deltaqp;
  _lhs[6]=pstar;
  _lhs[7]=pt0;
  _lhs[8]=q0;
}
// Functional based bioavailability
double rxModels_Ribba2012_F(int _cSub,  int _cmt, double _amt, double t){
 return _amt;
}
// Functional based absorption lag
double rxModels_Ribba2012_Lag(int _cSub,  int _cmt, double t){
 return t;
  double   k,
  tkde,
  eta_DoT_tkde,
  kde,
  tkpq,
  eta_DoT_kpq,
  kpq,
  tkqpp,
  eta_DoT_kqpp,
  kqpp,
  tlambdap,
  eta_DoT_lambdap,
  lambdap,
  tgamma,
  eta_DoT_gamma,
  gamma,
  tdeltaqp,
  eta_DoT_deltaqp,
  deltaqp,
  pstar,
  pt,
  q,
  qp,
  c,
  tpt0,
  eta_DoT_pt0,
  pt0,
  tq0,
  eta_DoT_q0,
  q0;

  (void)t;
  (void)k;
  (void)tkde;
  (void)eta_DoT_tkde;
  (void)kde;
  (void)tkpq;
  (void)eta_DoT_kpq;
  (void)kpq;
  (void)tkqpp;
  (void)eta_DoT_kqpp;
  (void)kqpp;
  (void)tlambdap;
  (void)eta_DoT_lambdap;
  (void)lambdap;
  (void)tgamma;
  (void)eta_DoT_gamma;
  (void)gamma;
  (void)tdeltaqp;
  (void)eta_DoT_deltaqp;
  (void)deltaqp;
  (void)pstar;
  (void)pt;
  (void)q;
  (void)qp;
  (void)c;
  (void)tpt0;
  (void)eta_DoT_pt0;
  (void)pt0;
  (void)tq0;
  (void)eta_DoT_q0;
  (void)q0;

  _update_par_ptr(NA_REAL, _cSub, _solveData, _idx);
  k = _PP[0];
  tkde = _PP[1];
  eta_DoT_tkde = _PP[2];
  tkpq = _PP[3];
  eta_DoT_kpq = _PP[4];
  tkqpp = _PP[5];
  eta_DoT_kqpp = _PP[6];
  tlambdap = _PP[7];
  eta_DoT_lambdap = _PP[8];
  tgamma = _PP[9];
  eta_DoT_gamma = _PP[10];
  tdeltaqp = _PP[11];
  eta_DoT_deltaqp = _PP[12];
  tpt0 = _PP[13];
  eta_DoT_pt0 = _PP[14];
  tq0 = _PP[15];
  eta_DoT_q0 = _PP[16];

}
// Modeled zero-order rate
double rxModels_Ribba2012_Rate(int _cSub,  int _cmt, double _amt, double t){
 return 0.0;
  double   k,
  tkde,
  eta_DoT_tkde,
  kde,
  tkpq,
  eta_DoT_kpq,
  kpq,
  tkqpp,
  eta_DoT_kqpp,
  kqpp,
  tlambdap,
  eta_DoT_lambdap,
  lambdap,
  tgamma,
  eta_DoT_gamma,
  gamma,
  tdeltaqp,
  eta_DoT_deltaqp,
  deltaqp,
  pstar,
  pt,
  q,
  qp,
  c,
  tpt0,
  eta_DoT_pt0,
  pt0,
  tq0,
  eta_DoT_q0,
  q0;

  (void)t;
  (void)k;
  (void)tkde;
  (void)eta_DoT_tkde;
  (void)kde;
  (void)tkpq;
  (void)eta_DoT_kpq;
  (void)kpq;
  (void)tkqpp;
  (void)eta_DoT_kqpp;
  (void)kqpp;
  (void)tlambdap;
  (void)eta_DoT_lambdap;
  (void)lambdap;
  (void)tgamma;
  (void)eta_DoT_gamma;
  (void)gamma;
  (void)tdeltaqp;
  (void)eta_DoT_deltaqp;
  (void)deltaqp;
  (void)pstar;
  (void)pt;
  (void)q;
  (void)qp;
  (void)c;
  (void)tpt0;
  (void)eta_DoT_pt0;
  (void)pt0;
  (void)tq0;
  (void)eta_DoT_q0;
  (void)q0;

  _update_par_ptr(NA_REAL, _cSub, _solveData, _idx);
  k = _PP[0];
  tkde = _PP[1];
  eta_DoT_tkde = _PP[2];
  tkpq = _PP[3];
  eta_DoT_kpq = _PP[4];
  tkqpp = _PP[5];
  eta_DoT_kqpp = _PP[6];
  tlambdap = _PP[7];
  eta_DoT_lambdap = _PP[8];
  tgamma = _PP[9];
  eta_DoT_gamma = _PP[10];
  tdeltaqp = _PP[11];
  eta_DoT_deltaqp = _PP[12];
  tpt0 = _PP[13];
  eta_DoT_pt0 = _PP[14];
  tq0 = _PP[15];
  eta_DoT_q0 = _PP[16];

}
// Modeled zero-order duration
double rxModels_Ribba2012_Dur(int _cSub,  int _cmt, double _amt, double t){
 return 0.0;
}
// Model Times
void rxModels_Ribba2012_mtime(int _cSub, double *_mtime){
}
extern SEXP rxModels_Ribba2012_model_vars(){
  int pro=0;
  SEXP _mv = PROTECT(_rxGetModelLib("rxModels_Ribba2012_model_vars"));pro++;
  if (!_rxIsCurrentC(_mv)){
    SEXP lst      = PROTECT(allocVector(VECSXP, 20));pro++;
    SEXP names    = PROTECT(allocVector(STRSXP, 20));pro++;
    SEXP sNeedSort = PROTECT(allocVector(INTSXP,1));pro++;
    int *iNeedSort  = INTEGER(sNeedSort);
    iNeedSort[0] = 84;
    SEXP sMtime = PROTECT(allocVector(INTSXP,1));pro++;
    int *iMtime  = INTEGER(sMtime);
    iMtime[0] = 0;
    SEXP sExtraCmt = PROTECT(allocVector(INTSXP,1));pro++;
    int *iExtraCmt  = INTEGER(sExtraCmt);
    iExtraCmt[0] = 0;
    SEXP params   = PROTECT(allocVector(STRSXP, 17));pro++;
    SEXP lhs      = PROTECT(allocVector(STRSXP, 9));pro++;
    SEXP state    = PROTECT(allocVector(STRSXP, 4));pro++;
  SEXP extraState = PROTECT(allocVector(STRSXP, 0));pro++;
    SEXP stateRmS = PROTECT(allocVector(INTSXP, 4));pro++;
    SEXP timeInt = PROTECT(allocVector(INTSXP, 1));pro++;
    INTEGER(timeInt)[0] = 1559282355;
    SEXP sens     = PROTECT(allocVector(STRSXP, 0));pro++;
    SEXP normState= PROTECT(allocVector(STRSXP, 4));pro++;
    SEXP fn_ini   = PROTECT(allocVector(STRSXP, 0));pro++;
    SEXP dfdy     = PROTECT(allocVector(STRSXP, 0));pro++;
    SEXP tran     = PROTECT(allocVector(STRSXP, 20));pro++;
    SEXP trann    = PROTECT(allocVector(STRSXP, 20));pro++;
    SEXP mmd5     = PROTECT(allocVector(STRSXP, 2));pro++;
    SEXP mmd5n    = PROTECT(allocVector(STRSXP, 2));pro++;
    SEXP model    = PROTECT(allocVector(STRSXP, 1));pro++;
    SEXP modeln   = PROTECT(allocVector(STRSXP, 1));pro++;
    SEXP version    = PROTECT(allocVector(STRSXP, 3));pro++;
    SEXP versionn   = PROTECT(allocVector(STRSXP, 3));pro++;
    SET_STRING_ELT(version,0,mkChar("0.8.1-0"));
    SET_STRING_ELT(version,1,mkChar("https://github.com/nlmixrdevelopment/RxODE"));
    SET_STRING_ELT(version,2,mkChar("2ef740cbdf83ec58b612d7c6d08130e4"));
    SET_STRING_ELT(versionn,0,mkChar("version"));
    SET_STRING_ELT(versionn,1,mkChar("repo"));
    SET_STRING_ELT(versionn,2,mkChar("md5"));
    SET_STRING_ELT(params,0,mkChar("k"));
    SET_STRING_ELT(params,1,mkChar("tkde"));
    SET_STRING_ELT(params,2,mkChar("eta.tkde"));
  SET_STRING_ELT(lhs,0,mkChar("kde"));
    SET_STRING_ELT(params,3,mkChar("tkpq"));
    SET_STRING_ELT(params,4,mkChar("eta.kpq"));
  SET_STRING_ELT(lhs,1,mkChar("kpq"));
    SET_STRING_ELT(params,5,mkChar("tkqpp"));
    SET_STRING_ELT(params,6,mkChar("eta.kqpp"));
  SET_STRING_ELT(lhs,2,mkChar("kqpp"));
    SET_STRING_ELT(params,7,mkChar("tlambdap"));
    SET_STRING_ELT(params,8,mkChar("eta.lambdap"));
  SET_STRING_ELT(lhs,3,mkChar("lambdap"));
    SET_STRING_ELT(params,9,mkChar("tgamma"));
    SET_STRING_ELT(params,10,mkChar("eta.gamma"));
  SET_STRING_ELT(lhs,4,mkChar("gamma"));
    SET_STRING_ELT(params,11,mkChar("tdeltaqp"));
    SET_STRING_ELT(params,12,mkChar("eta.deltaqp"));
  SET_STRING_ELT(lhs,5,mkChar("deltaqp"));
  SET_STRING_ELT(lhs,6,mkChar("pstar"));
    SET_STRING_ELT(params,13,mkChar("tpt0"));
    SET_STRING_ELT(params,14,mkChar("eta.pt0"));
  SET_STRING_ELT(lhs,7,mkChar("pt0"));
    SET_STRING_ELT(params,15,mkChar("tq0"));
    SET_STRING_ELT(params,16,mkChar("eta.q0"));
  SET_STRING_ELT(lhs,8,mkChar("q0"));
    SET_STRING_ELT(state,0,mkChar("c"));
    SET_STRING_ELT(normState,0,mkChar("c"));
    _SR[0] = 0;
    SET_STRING_ELT(state,1,mkChar("pt"));
    SET_STRING_ELT(normState,1,mkChar("pt"));
    _SR[1] = 0;
    SET_STRING_ELT(state,2,mkChar("q"));
    SET_STRING_ELT(normState,2,mkChar("q"));
    _SR[2] = 0;
    SET_STRING_ELT(state,3,mkChar("qp"));
    SET_STRING_ELT(normState,3,mkChar("qp"));
    _SR[3] = 0;
    SET_STRING_ELT(modeln,0,mkChar("normModel"));
    SET_STRING_ELT(model,0,mkChar("k=100;\ntkde=0.24;\neta.tkde=0;\nkde=tkde*exp(eta.tkde);\ntkpq=0.0295;\neta.kpq=0;\nkpq=tkpq*exp(eta.kpq);\ntkqpp=0.0031;\neta.kqpp=0;\nkqpp=tkqpp*exp(eta.kqpp);\ntlambdap=0.121;\neta.lambdap=0;\nlambdap=tlambdap*exp(eta.lambdap);\ntgamma=0.729;\neta.gamma=0;\ngamma=tgamma*exp(eta.gamma);\ntdeltaqp=0.00867;\neta.deltaqp=0;\ndeltaqp=tdeltaqp*exp(eta.deltaqp);\npstar=pt+q+qp;\nd/dt(c)=-kde*c;\nd/dt(pt)=lambdap*pt*(1-pstar/k)+kqpp*qp-kpq*pt-gamma*c*kde*pt;\nd/dt(q)=kpq*pt-gamma*c*kde*q;\nd/dt(qp)=gamma*c*kde*q-kqpp*qp-deltaqp*qp;\ntpt0=7.13;\neta.pt0=0;\npt0=tpt0*exp(eta.pt0);\ntq0=41.2;\neta.q0=0;\nq0=tq0*exp(eta.q0);\npt(0)=pt0;\nq(0)=q0;\n"));
    SEXP ini    = PROTECT(allocVector(REALSXP,17));pro++;
    SEXP inin   = PROTECT(allocVector(STRSXP, 17));pro++;
    SET_STRING_ELT(inin,0,mkChar("k"));
    REAL(ini)[0] = 100.0000000000000000;
    SET_STRING_ELT(inin,1,mkChar("tkde"));
    REAL(ini)[1] = 0.2400000000000000;
    SET_STRING_ELT(inin,2,mkChar("eta.tkde"));
    REAL(ini)[2] = 0.0000000000000000;
    SET_STRING_ELT(inin,3,mkChar("tkpq"));
    REAL(ini)[3] = 0.0295000000000000;
    SET_STRING_ELT(inin,4,mkChar("eta.kpq"));
    REAL(ini)[4] = 0.0000000000000000;
    SET_STRING_ELT(inin,5,mkChar("tkqpp"));
    REAL(ini)[5] = 0.0031000000000000;
    SET_STRING_ELT(inin,6,mkChar("eta.kqpp"));
    REAL(ini)[6] = 0.0000000000000000;
    SET_STRING_ELT(inin,7,mkChar("tlambdap"));
    REAL(ini)[7] = 0.1210000000000000;
    SET_STRING_ELT(inin,8,mkChar("eta.lambdap"));
    REAL(ini)[8] = 0.0000000000000000;
    SET_STRING_ELT(inin,9,mkChar("tgamma"));
    REAL(ini)[9] = 0.7290000000000000;
    SET_STRING_ELT(inin,10,mkChar("eta.gamma"));
    REAL(ini)[10] = 0.0000000000000000;
    SET_STRING_ELT(inin,11,mkChar("tdeltaqp"));
    REAL(ini)[11] = 0.0086700000000000;
    SET_STRING_ELT(inin,12,mkChar("eta.deltaqp"));
    REAL(ini)[12] = 0.0000000000000000;
    SET_STRING_ELT(inin,13,mkChar("tpt0"));
    REAL(ini)[13] = 7.1299999999999999;
    SET_STRING_ELT(inin,14,mkChar("eta.pt0"));
    REAL(ini)[14] = 0.0000000000000000;
    SET_STRING_ELT(inin,15,mkChar("tq0"));
    REAL(ini)[15] = 41.2000000000000030;
    SET_STRING_ELT(inin,16,mkChar("eta.q0"));
    REAL(ini)[16] = 0.0000000000000000;
    SET_STRING_ELT(names,0,mkChar("params"));
    SET_VECTOR_ELT(lst,  0,params);
    SET_STRING_ELT(names,1,mkChar("lhs"));
    SET_VECTOR_ELT(lst,  1,lhs);
    SET_STRING_ELT(names,2,mkChar("state"));
    SET_VECTOR_ELT(lst,  2,state);
    SET_STRING_ELT(names,3,mkChar("trans"));
    SET_VECTOR_ELT(lst,  3,tran);
    SET_STRING_ELT(names,4,mkChar("model"));
    SET_VECTOR_ELT(lst,  4,model);
    SET_STRING_ELT(names,5,mkChar("ini"));
    SET_VECTOR_ELT(lst,  5,ini);
    SET_STRING_ELT(names,6,mkChar("podo"));
    SET_VECTOR_ELT(lst,   6,ScalarLogical(0));
    SET_STRING_ELT(names,7,mkChar("dfdy"));
    SET_VECTOR_ELT(lst,  7,dfdy);
    SET_STRING_ELT(names,8,mkChar("sens"));
    SET_VECTOR_ELT(lst,  8,sens);
    SET_STRING_ELT(names,9,mkChar("fn.ini"));
    SET_VECTOR_ELT(lst,  9,fn_ini);
    SET_STRING_ELT(names,10,mkChar("state.ignore"));
    SET_VECTOR_ELT(lst,  10,stateRmS);
    SET_STRING_ELT(names,11,mkChar("version"));
    SET_VECTOR_ELT(lst,  11,version);
    SET_STRING_ELT(names,12,mkChar("normal.state"));
    SET_VECTOR_ELT(lst,  12,normState);
    SET_STRING_ELT(names,13,mkChar("needSort"));
    SET_VECTOR_ELT(lst,  13,sNeedSort);
    SET_STRING_ELT(names,14,mkChar("nMtime"));
    SET_VECTOR_ELT(lst,  14,sMtime);
    SET_STRING_ELT(names,15,mkChar("extraCmt"));
    SET_VECTOR_ELT(lst,  15,sExtraCmt);
    SET_STRING_ELT(names, 16, mkChar("stateExtra"));
    SET_VECTOR_ELT(lst,  16, extraState);
    SET_STRING_ELT(names, 17, mkChar("dvid"));
    SEXP sDvid = PROTECT(allocVector(INTSXP,0));pro++;
    SET_VECTOR_ELT(lst, 17, sDvid);
    SET_STRING_ELT(names,18,mkChar("timeId"));
    SET_VECTOR_ELT(lst,  18,timeInt);
    SET_STRING_ELT(names,19,mkChar("md5"));    SET_VECTOR_ELT(lst,  19,mmd5);    SET_STRING_ELT(mmd5n,0,mkChar("file_md5"));
    SET_STRING_ELT(mmd5,0,mkChar("8d2c0b98fcae5c10350179b868ceb202"));
    SET_STRING_ELT(mmd5n,1,mkChar("parsed_md5"));
    SET_STRING_ELT(mmd5,1,mkChar("8d2c0b98fcae5c10350179b868ceb202"));
    SET_STRING_ELT(trann,0,mkChar("lib.name"));
    SET_STRING_ELT(tran, 0,mkChar("rxModels"));
    SET_STRING_ELT(trann,1,mkChar("jac"));
    SET_STRING_ELT(tran,1,mkChar("fullint"));
    SET_STRING_ELT(trann,2,mkChar("prefix"));
    SET_STRING_ELT(tran, 2,mkChar("rxModels_Ribba2012_"));
    SET_STRING_ELT(trann,3,mkChar("dydt"));
    SET_STRING_ELT(tran, 3,mkChar("rxModels_Ribba2012_dydt"));
    SET_STRING_ELT(trann,4,mkChar("calc_jac"));
    SET_STRING_ELT(tran, 4,mkChar("rxModels_Ribba2012_calc_jac"));
    SET_STRING_ELT(trann,5,mkChar("calc_lhs"));
    SET_STRING_ELT(tran, 5,mkChar("rxModels_Ribba2012_calc_lhs"));
    SET_STRING_ELT(trann,6,mkChar("model_vars"));
    SET_STRING_ELT(tran, 6,mkChar("rxModels_Ribba2012_model_vars"));
    SET_STRING_ELT(trann,7,mkChar("theta"));
    SET_STRING_ELT(tran, 7,mkChar("rxModels_Ribba2012_theta"));
    SET_STRING_ELT(trann,8,mkChar("inis"));
    SET_STRING_ELT(tran, 8,mkChar("rxModels_Ribba2012_inis"));
    SET_STRING_ELT(trann,  9,mkChar("dydt_lsoda"));
    SET_STRING_ELT(tran,   9,mkChar("rxModels_Ribba2012_dydt_lsoda"));
    SET_STRING_ELT(trann,10,mkChar("calc_jac_lsoda"));
    SET_STRING_ELT(tran, 10,mkChar("rxModels_Ribba2012_calc_jac_lsoda"));
    SET_STRING_ELT(trann,11,mkChar("ode_solver_solvedata"));
    SET_STRING_ELT(tran, 11,mkChar("rxModels_Ribba2012_ode_solver_solvedata"));
    SET_STRING_ELT(trann,12,mkChar("ode_solver_get_solvedata"));
    SET_STRING_ELT(tran, 12,mkChar("rxModels_Ribba2012_ode_solver_get_solvedata"));
    SET_STRING_ELT(trann,13,mkChar("dydt_liblsoda"));
    SET_STRING_ELT(tran, 13,mkChar("rxModels_Ribba2012_dydt_liblsoda"));
    SET_STRING_ELT(trann,14,mkChar("F"));
    SET_STRING_ELT(tran, 14,mkChar("rxModels_Ribba2012_F"));
    SET_STRING_ELT(trann,15,mkChar("Lag"));
    SET_STRING_ELT(tran, 15,mkChar("rxModels_Ribba2012_Lag"));
    SET_STRING_ELT(trann,16,mkChar("Rate"));
    SET_STRING_ELT(tran, 16,mkChar("rxModels_Ribba2012_Rate"));
    SET_STRING_ELT(trann,17,mkChar("Dur"));
    SET_STRING_ELT(tran, 17,mkChar("rxModels_Ribba2012_Dur"));
    SET_STRING_ELT(trann,18,mkChar("mtime"));
    SET_STRING_ELT(tran, 18,mkChar("rxModels_Ribba2012_mtime"));
    SET_STRING_ELT(trann,19,mkChar("assignFuns"));
    SET_STRING_ELT(tran, 19,mkChar("rxModels_Ribba2012_assignFuns"));
    setAttrib(tran, R_NamesSymbol, trann);
    setAttrib(mmd5, R_NamesSymbol, mmd5n);
    setAttrib(model, R_NamesSymbol, modeln);
    setAttrib(ini, R_NamesSymbol, inin);
    setAttrib(version, R_NamesSymbol, versionn);
    setAttrib(lst, R_NamesSymbol, names);
    SEXP cls = PROTECT(allocVector(STRSXP, 1));pro++;
    SET_STRING_ELT(cls, 0, mkChar("rxModelVars"));
    classgets(lst, cls);
    _assign_ptr(lst);
    UNPROTECT(pro);
    return lst;
  } else {
    UNPROTECT(pro);
    return _mv;
  }
}
extern void rxModels_Ribba2012_dydt_lsoda(int *neq, double *t, double *A, double *DADT)
{
  rxModels_Ribba2012_dydt(neq, *t, A, DADT);
}
extern int rxModels_Ribba2012_dydt_liblsoda(double t, double *y, double *ydot, void *data)
{
  int *neq = (int*)(data);
  rxModels_Ribba2012_dydt(neq, t, y, ydot);
  return(0);
}
extern void rxModels_Ribba2012_calc_jac_lsoda(int *neq, double *t, double *A,int *ml, int *mu, double *JAC, int *nrowpd){
  // Update all covariate parameters
  rxModels_Ribba2012_calc_jac(neq, *t, A, JAC, *nrowpd);
}

//Create function to call from R's main thread that assigns the required functions. Sometimes they don't get assigned.
extern void rxModels_Ribba2012_assignFuns(){
  _assignFuns();
}

//Initialize the dll to match RxODE's calls
void R_init0_rxModels_Ribba2012(){
  // Get C callables on load; Otherwise it isn't thread safe
  _assignFuns();
  R_RegisterCCallable("rxModels","rxModels_Ribba2012_assignFuns", (DL_FUNC) rxModels_Ribba2012_assignFuns);
  R_RegisterCCallable("rxModels","rxModels_Ribba2012_theta", (DL_FUNC) rxModels_Ribba2012_theta);
  R_RegisterCCallable("rxModels","rxModels_Ribba2012_inis",(DL_FUNC) rxModels_Ribba2012_inis);
  R_RegisterCCallable("rxModels","rxModels_Ribba2012_dydt",(DL_FUNC) rxModels_Ribba2012_dydt);
  R_RegisterCCallable("rxModels","rxModels_Ribba2012_calc_lhs",(DL_FUNC) rxModels_Ribba2012_calc_lhs);
  R_RegisterCCallable("rxModels","rxModels_Ribba2012_calc_jac",(DL_FUNC) rxModels_Ribba2012_calc_jac);
  R_RegisterCCallable("rxModels","rxModels_Ribba2012_dydt_lsoda", (DL_FUNC) rxModels_Ribba2012_dydt_lsoda);
  R_RegisterCCallable("rxModels","rxModels_Ribba2012_calc_jac_lsoda", (DL_FUNC) rxModels_Ribba2012_calc_jac_lsoda);
  R_RegisterCCallable("rxModels","rxModels_Ribba2012_ode_solver_solvedata", (DL_FUNC) rxModels_Ribba2012_ode_solver_solvedata);
  R_RegisterCCallable("rxModels","rxModels_Ribba2012_ode_solver_get_solvedata", (DL_FUNC) rxModels_Ribba2012_ode_solver_get_solvedata);
  R_RegisterCCallable("rxModels","rxModels_Ribba2012_F", (DL_FUNC) rxModels_Ribba2012_F);
  R_RegisterCCallable("rxModels","rxModels_Ribba2012_Lag", (DL_FUNC) rxModels_Ribba2012_Lag);
  R_RegisterCCallable("rxModels","rxModels_Ribba2012_Rate", (DL_FUNC) rxModels_Ribba2012_Rate);
  R_RegisterCCallable("rxModels","rxModels_Ribba2012_Dur", (DL_FUNC) rxModels_Ribba2012_Dur);
  R_RegisterCCallable("rxModels","rxModels_Ribba2012_mtime", (DL_FUNC) rxModels_Ribba2012_mtime);
  R_RegisterCCallable("rxModels","rxModels_Ribba2012_dydt_liblsoda", (DL_FUNC) rxModels_Ribba2012_dydt_liblsoda);
}
//Initialize the dll to match RxODE's calls
void R_init_rxModels_Ribba2012(DllInfo *info){
  // Get C callables on load; Otherwise it isn't thread safe
  R_init0_rxModels_Ribba2012();
  static const R_CallMethodDef callMethods[]  = {
    {"rxModels_Ribba2012_model_vars", (DL_FUNC) &rxModels_Ribba2012_model_vars, 0},
    {NULL, NULL, 0}
  };

  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info,FALSE);
}

void R_unload_rxModels_Ribba2012 (DllInfo *info){
  // Free resources required for single subject solve.
  SEXP _mv = PROTECT(_rxGetModelLib("rxModels_Ribba2012_model_vars"));
  if (!isNull(_mv)){
    _rxRmModelLib("rxModels_Ribba2012_model_vars");
  }
  UNPROTECT(1);
}
