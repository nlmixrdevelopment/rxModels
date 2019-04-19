#include <RxODE.h>
#include <RxODE_model.h>
#define __MAX_PROD__ 0
extern void  rxModels_pk2cmt_ode_solver_solvedata (rx_solve *solve){
  _solveData = solve;
}
extern rx_solve *rxModels_pk2cmt_ode_solver_get_solvedata(){
  return _solveData;
}
SEXP rxModels_pk2cmt_model_vars();
double _theta[18];
extern double* rxModels_pk2cmt_theta(double *theta){
  _theta[0] = 1.0000000000000000; _theta[1] = 20.0000000000000000; _theta[2] = 1.0000000000000000; _theta[3] = 10.0000000000000000; _theta[4] = 2.0000000000000000; _theta[5] = 0.0000000000000000; _theta[6] = 0.0000000000000000; _theta[7] = 0.0000000000000000; _theta[8] = 0.0000000000000000; _theta[9] = 0.0000000000000000; _theta[10] = 0.0000000000000000; _theta[11] = 0.0000000000000000; _theta[12] = 0.0000000000000000; _theta[13] = 0.0000000000000000; _theta[14] = 0.0000000000000000; _theta[15] = 0.0000000000000000; _theta[16] = 0.0000000000000000; _theta[17] = 0.0000000000000000;
  return _theta;
}


// prj-specific differential eqns
void rxModels_pk2cmt_dydt(int *_neq, double t, double *__zzStateVar__, double *__DDtStateVar__)
{
  int _cSub = _neq[1];
  double   popCl,
  popV,
  popKa,
  popVp,
  popQ,
  bsvCl,
  bsvV,
  bsvKa,
  bsvVp,
  bsvQ,
  cl,
  v,
  ka,
  q,
  vp,
  popLagDepot,
  popLagCentral,
  popRateCentral,
  popDurCentral,
  bsvLagDepot,
  bsvLagCentral,
  bsvRateCentral,
  bsvDurCentral,
  rx_ka,
  rx_rate,
  rx_dur,
  rx_tlag,
  rx_tlag2,
  rx_F,
  rx_F2,
  rx_v,
  rx_k,
  rx_k12,
  rx_k21,
  rx_beta,
  rx_alpha,
  rx_A,
  rx_B,
  rx_A2,
  rx_B2,
  rx_gamma,
  rx_C,
  rx_C2,
  cp;

  (void)t;
  (void)popCl;
  (void)popV;
  (void)popKa;
  (void)popVp;
  (void)popQ;
  (void)bsvCl;
  (void)bsvV;
  (void)bsvKa;
  (void)bsvVp;
  (void)bsvQ;
  (void)cl;
  (void)v;
  (void)ka;
  (void)q;
  (void)vp;
  (void)popLagDepot;
  (void)popLagCentral;
  (void)popRateCentral;
  (void)popDurCentral;
  (void)bsvLagDepot;
  (void)bsvLagCentral;
  (void)bsvRateCentral;
  (void)bsvDurCentral;
  (void)rx_ka;
  (void)rx_rate;
  (void)rx_dur;
  (void)rx_tlag;
  (void)rx_tlag2;
  (void)rx_F;
  (void)rx_F2;
  (void)rx_v;
  (void)rx_k;
  (void)rx_k12;
  (void)rx_k21;
  (void)rx_beta;
  (void)rx_alpha;
  (void)rx_A;
  (void)rx_B;
  (void)rx_A2;
  (void)rx_B2;
  (void)rx_gamma;
  (void)rx_C;
  (void)rx_C2;
  (void)cp;

  _update_par_ptr(t, _cSub, _solveData, _idx);
  popCl = _PP[0];
  popV = _PP[1];
  popKa = _PP[2];
  popVp = _PP[3];
  popQ = _PP[4];
  bsvCl = _PP[5];
  bsvV = _PP[6];
  bsvKa = _PP[7];
  bsvVp = _PP[8];
  bsvQ = _PP[9];
  popLagDepot = _PP[10];
  popLagCentral = _PP[11];
  popRateCentral = _PP[12];
  popDurCentral = _PP[13];
  bsvLagDepot = _PP[14];
  bsvLagCentral = _PP[15];
  bsvRateCentral = _PP[16];
  bsvDurCentral = _PP[17];


  cl =popCl*exp(bsvCl);
  v =popV*exp(bsvV);
  ka =popKa*exp(bsvKa);
  q =popQ*exp(bsvQ);
  vp =popVp*exp(bsvVp);
  rx_ka =ka;
  rx_rate =popRateCentral*exp(bsvRateCentral);
  rx_dur =popDurCentral*exp(bsvDurCentral);
  rx_tlag =popLagDepot*exp(bsvLagDepot);
  rx_tlag2 =popLagCentral*exp(bsvLagCentral);
  rx_F =1;
  rx_F2 =1;
  rx_v =v;
  rx_k =cl/safe_zero(v);
  rx_k12 =q/safe_zero(v);
  rx_k21 =q/safe_zero(vp);
  rx_beta =0.5*(rx_k12+rx_k21+rx_k-sqrt((rx_k12+rx_k21+rx_k)*(rx_k12+rx_k21+rx_k)-4.0*rx_k21*rx_k));
  rx_alpha =rx_k21*rx_k/safe_zero(rx_beta);
  rx_A =rx_ka/safe_zero((rx_ka-rx_alpha))*(rx_alpha-rx_k21)/safe_zero((rx_alpha-rx_beta))/safe_zero(rx_v);
  rx_B =rx_ka/safe_zero((rx_ka-rx_beta))*(rx_beta-rx_k21)/safe_zero((rx_beta-rx_alpha))/safe_zero(rx_v);
  rx_A2 =(rx_alpha-rx_k21)/safe_zero((rx_alpha-rx_beta))/safe_zero(rx_v);
  rx_B2 =(rx_beta-rx_k21)/safe_zero((rx_beta-rx_alpha))/safe_zero(rx_v);
  rx_gamma =0;
  rx_C =0;
  rx_C2 =0;
  cp=solveLinB(_solveData, _cSub,t,0,rx_A,rx_A2,rx_alpha,rx_B,rx_B2,rx_beta,rx_C,rx_C2,rx_gamma,rx_ka,rx_tlag,rx_tlag2,rx_F,rx_F2,rx_rate,rx_dur);
  (&_solveData->subjects[_cSub])->dadt_counter[0]++;
}

// Jacobian derived vars
void rxModels_pk2cmt_calc_jac(int *_neq, double t, double *__zzStateVar__, double *__PDStateVar__, unsigned int __NROWPD__) {
  int _cSub=_neq[1];
  (&_solveData->subjects[_cSub])->jac_counter[0]++;
}
// Functional based initial conditions.
void rxModels_pk2cmt_inis(int _cSub, double *__zzStateVar__){
}
// prj-specific derived vars
void rxModels_pk2cmt_calc_lhs(int _cSub, double t, double *__zzStateVar__, double *_lhs) {
  double   popCl,
  popV,
  popKa,
  popVp,
  popQ,
  bsvCl,
  bsvV,
  bsvKa,
  bsvVp,
  bsvQ,
  cl,
  v,
  ka,
  q,
  vp,
  popLagDepot,
  popLagCentral,
  popRateCentral,
  popDurCentral,
  bsvLagDepot,
  bsvLagCentral,
  bsvRateCentral,
  bsvDurCentral,
  rx_ka,
  rx_rate,
  rx_dur,
  rx_tlag,
  rx_tlag2,
  rx_F,
  rx_F2,
  rx_v,
  rx_k,
  rx_k12,
  rx_k21,
  rx_beta,
  rx_alpha,
  rx_A,
  rx_B,
  rx_A2,
  rx_B2,
  rx_gamma,
  rx_C,
  rx_C2,
  cp;

  (void)t;
  (void)popCl;
  (void)popV;
  (void)popKa;
  (void)popVp;
  (void)popQ;
  (void)bsvCl;
  (void)bsvV;
  (void)bsvKa;
  (void)bsvVp;
  (void)bsvQ;
  (void)cl;
  (void)v;
  (void)ka;
  (void)q;
  (void)vp;
  (void)popLagDepot;
  (void)popLagCentral;
  (void)popRateCentral;
  (void)popDurCentral;
  (void)bsvLagDepot;
  (void)bsvLagCentral;
  (void)bsvRateCentral;
  (void)bsvDurCentral;
  (void)rx_ka;
  (void)rx_rate;
  (void)rx_dur;
  (void)rx_tlag;
  (void)rx_tlag2;
  (void)rx_F;
  (void)rx_F2;
  (void)rx_v;
  (void)rx_k;
  (void)rx_k12;
  (void)rx_k21;
  (void)rx_beta;
  (void)rx_alpha;
  (void)rx_A;
  (void)rx_B;
  (void)rx_A2;
  (void)rx_B2;
  (void)rx_gamma;
  (void)rx_C;
  (void)rx_C2;
  (void)cp;

  _update_par_ptr(t, _cSub, _solveData, _idx);
  popCl = _PP[0];
  popV = _PP[1];
  popKa = _PP[2];
  popVp = _PP[3];
  popQ = _PP[4];
  bsvCl = _PP[5];
  bsvV = _PP[6];
  bsvKa = _PP[7];
  bsvVp = _PP[8];
  bsvQ = _PP[9];
  popLagDepot = _PP[10];
  popLagCentral = _PP[11];
  popRateCentral = _PP[12];
  popDurCentral = _PP[13];
  bsvLagDepot = _PP[14];
  bsvLagCentral = _PP[15];
  bsvRateCentral = _PP[16];
  bsvDurCentral = _PP[17];


  cl =popCl*exp(bsvCl);
  v =popV*exp(bsvV);
  ka =popKa*exp(bsvKa);
  q =popQ*exp(bsvQ);
  vp =popVp*exp(bsvVp);
  rx_ka =ka;
  rx_rate =popRateCentral*exp(bsvRateCentral);
  rx_dur =popDurCentral*exp(bsvDurCentral);
  rx_tlag =popLagDepot*exp(bsvLagDepot);
  rx_tlag2 =popLagCentral*exp(bsvLagCentral);
  rx_F =1;
  rx_F2 =1;
  rx_v =v;
  rx_k =cl/safe_zero(v);
  rx_k12 =q/safe_zero(v);
  rx_k21 =q/safe_zero(vp);
  rx_beta =0.5*(rx_k12+rx_k21+rx_k-sqrt((rx_k12+rx_k21+rx_k)*(rx_k12+rx_k21+rx_k)-4.0*rx_k21*rx_k));
  rx_alpha =rx_k21*rx_k/safe_zero(rx_beta);
  rx_A =rx_ka/safe_zero((rx_ka-rx_alpha))*(rx_alpha-rx_k21)/safe_zero((rx_alpha-rx_beta))/safe_zero(rx_v);
  rx_B =rx_ka/safe_zero((rx_ka-rx_beta))*(rx_beta-rx_k21)/safe_zero((rx_beta-rx_alpha))/safe_zero(rx_v);
  rx_A2 =(rx_alpha-rx_k21)/safe_zero((rx_alpha-rx_beta))/safe_zero(rx_v);
  rx_B2 =(rx_beta-rx_k21)/safe_zero((rx_beta-rx_alpha))/safe_zero(rx_v);
  rx_gamma =0;
  rx_C =0;
  rx_C2 =0;
  cp=solveLinB(_solveData, _cSub,t,0,rx_A,rx_A2,rx_alpha,rx_B,rx_B2,rx_beta,rx_C,rx_C2,rx_gamma,rx_ka,rx_tlag,rx_tlag2,rx_F,rx_F2,rx_rate,rx_dur);

  _lhs[0]=cp;
}
// Functional based bioavailability
double rxModels_pk2cmt_F(int _cSub,  int _cmt, double _amt, double t){
 return _amt;
}
// Functional based absorption lag
double rxModels_pk2cmt_Lag(int _cSub,  int _cmt, double t){
 return t;
  double   popCl,
  popV,
  popKa,
  popVp,
  popQ,
  bsvCl,
  bsvV,
  bsvKa,
  bsvVp,
  bsvQ,
  cl,
  v,
  ka,
  q,
  vp,
  popLagDepot,
  popLagCentral,
  popRateCentral,
  popDurCentral,
  bsvLagDepot,
  bsvLagCentral,
  bsvRateCentral,
  bsvDurCentral,
  rx_ka,
  rx_rate,
  rx_dur,
  rx_tlag,
  rx_tlag2,
  rx_F,
  rx_F2,
  rx_v,
  rx_k,
  rx_k12,
  rx_k21,
  rx_beta,
  rx_alpha,
  rx_A,
  rx_B,
  rx_A2,
  rx_B2,
  rx_gamma,
  rx_C,
  rx_C2,
  cp;

  (void)t;
  (void)popCl;
  (void)popV;
  (void)popKa;
  (void)popVp;
  (void)popQ;
  (void)bsvCl;
  (void)bsvV;
  (void)bsvKa;
  (void)bsvVp;
  (void)bsvQ;
  (void)cl;
  (void)v;
  (void)ka;
  (void)q;
  (void)vp;
  (void)popLagDepot;
  (void)popLagCentral;
  (void)popRateCentral;
  (void)popDurCentral;
  (void)bsvLagDepot;
  (void)bsvLagCentral;
  (void)bsvRateCentral;
  (void)bsvDurCentral;
  (void)rx_ka;
  (void)rx_rate;
  (void)rx_dur;
  (void)rx_tlag;
  (void)rx_tlag2;
  (void)rx_F;
  (void)rx_F2;
  (void)rx_v;
  (void)rx_k;
  (void)rx_k12;
  (void)rx_k21;
  (void)rx_beta;
  (void)rx_alpha;
  (void)rx_A;
  (void)rx_B;
  (void)rx_A2;
  (void)rx_B2;
  (void)rx_gamma;
  (void)rx_C;
  (void)rx_C2;
  (void)cp;

  _update_par_ptr(NA_REAL, _cSub, _solveData, _idx);
  popCl = _PP[0];
  popV = _PP[1];
  popKa = _PP[2];
  popVp = _PP[3];
  popQ = _PP[4];
  bsvCl = _PP[5];
  bsvV = _PP[6];
  bsvKa = _PP[7];
  bsvVp = _PP[8];
  bsvQ = _PP[9];
  popLagDepot = _PP[10];
  popLagCentral = _PP[11];
  popRateCentral = _PP[12];
  popDurCentral = _PP[13];
  bsvLagDepot = _PP[14];
  bsvLagCentral = _PP[15];
  bsvRateCentral = _PP[16];
  bsvDurCentral = _PP[17];

}
// Modeled zero-order rate
double rxModels_pk2cmt_Rate(int _cSub,  int _cmt, double _amt, double t){
 return 0.0;
  double   popCl,
  popV,
  popKa,
  popVp,
  popQ,
  bsvCl,
  bsvV,
  bsvKa,
  bsvVp,
  bsvQ,
  cl,
  v,
  ka,
  q,
  vp,
  popLagDepot,
  popLagCentral,
  popRateCentral,
  popDurCentral,
  bsvLagDepot,
  bsvLagCentral,
  bsvRateCentral,
  bsvDurCentral,
  rx_ka,
  rx_rate,
  rx_dur,
  rx_tlag,
  rx_tlag2,
  rx_F,
  rx_F2,
  rx_v,
  rx_k,
  rx_k12,
  rx_k21,
  rx_beta,
  rx_alpha,
  rx_A,
  rx_B,
  rx_A2,
  rx_B2,
  rx_gamma,
  rx_C,
  rx_C2,
  cp;

  (void)t;
  (void)popCl;
  (void)popV;
  (void)popKa;
  (void)popVp;
  (void)popQ;
  (void)bsvCl;
  (void)bsvV;
  (void)bsvKa;
  (void)bsvVp;
  (void)bsvQ;
  (void)cl;
  (void)v;
  (void)ka;
  (void)q;
  (void)vp;
  (void)popLagDepot;
  (void)popLagCentral;
  (void)popRateCentral;
  (void)popDurCentral;
  (void)bsvLagDepot;
  (void)bsvLagCentral;
  (void)bsvRateCentral;
  (void)bsvDurCentral;
  (void)rx_ka;
  (void)rx_rate;
  (void)rx_dur;
  (void)rx_tlag;
  (void)rx_tlag2;
  (void)rx_F;
  (void)rx_F2;
  (void)rx_v;
  (void)rx_k;
  (void)rx_k12;
  (void)rx_k21;
  (void)rx_beta;
  (void)rx_alpha;
  (void)rx_A;
  (void)rx_B;
  (void)rx_A2;
  (void)rx_B2;
  (void)rx_gamma;
  (void)rx_C;
  (void)rx_C2;
  (void)cp;

  _update_par_ptr(NA_REAL, _cSub, _solveData, _idx);
  popCl = _PP[0];
  popV = _PP[1];
  popKa = _PP[2];
  popVp = _PP[3];
  popQ = _PP[4];
  bsvCl = _PP[5];
  bsvV = _PP[6];
  bsvKa = _PP[7];
  bsvVp = _PP[8];
  bsvQ = _PP[9];
  popLagDepot = _PP[10];
  popLagCentral = _PP[11];
  popRateCentral = _PP[12];
  popDurCentral = _PP[13];
  bsvLagDepot = _PP[14];
  bsvLagCentral = _PP[15];
  bsvRateCentral = _PP[16];
  bsvDurCentral = _PP[17];

}
// Modeled zero-order duration
double rxModels_pk2cmt_Dur(int _cSub,  int _cmt, double _amt, double t){
 return 0.0;
}
// Model Times
void rxModels_pk2cmt_mtime(int _cSub, double *_mtime){
}
extern SEXP rxModels_pk2cmt_model_vars(){
  int pro=0;
  SEXP _mv = PROTECT(_rxGetModelLib("rxModels_pk2cmt_model_vars"));pro++;
  if (!_rxIsCurrentC(_mv)){
    SEXP lst      = PROTECT(allocVector(VECSXP, 18));pro++;
    SEXP names    = PROTECT(allocVector(STRSXP, 18));pro++;
    SEXP sNeedSort = PROTECT(allocVector(INTSXP,1));pro++;
    int *iNeedSort  = INTEGER(sNeedSort);
    iNeedSort[0] = 0;
    SEXP sMtime = PROTECT(allocVector(INTSXP,1));pro++;
    int *iMtime  = INTEGER(sMtime);
    iMtime[0] = 0;
    SEXP sExtraCmt = PROTECT(allocVector(INTSXP,1));pro++;
    int *iExtraCmt  = INTEGER(sExtraCmt);
    iExtraCmt[0] = 2;
    SEXP params   = PROTECT(allocVector(STRSXP, 18));pro++;
    SEXP lhs      = PROTECT(allocVector(STRSXP, 1));pro++;
    SEXP state    = PROTECT(allocVector(STRSXP, 0));pro++;
    SEXP stateRmS = PROTECT(allocVector(INTSXP, 0));pro++;
    SEXP timeInt = PROTECT(allocVector(INTSXP, 1));pro++;
    INTEGER(timeInt)[0] = 1555680036;
    SEXP sens     = PROTECT(allocVector(STRSXP, 0));pro++;
    SEXP normState= PROTECT(allocVector(STRSXP, 0));pro++;
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
    SET_STRING_ELT(params,0,mkChar("popCl"));
    SET_STRING_ELT(params,1,mkChar("popV"));
    SET_STRING_ELT(params,2,mkChar("popKa"));
    SET_STRING_ELT(params,3,mkChar("popVp"));
    SET_STRING_ELT(params,4,mkChar("popQ"));
    SET_STRING_ELT(params,5,mkChar("bsvCl"));
    SET_STRING_ELT(params,6,mkChar("bsvV"));
    SET_STRING_ELT(params,7,mkChar("bsvKa"));
    SET_STRING_ELT(params,8,mkChar("bsvVp"));
    SET_STRING_ELT(params,9,mkChar("bsvQ"));
    SET_STRING_ELT(params,10,mkChar("popLagDepot"));
    SET_STRING_ELT(params,11,mkChar("popLagCentral"));
    SET_STRING_ELT(params,12,mkChar("popRateCentral"));
    SET_STRING_ELT(params,13,mkChar("popDurCentral"));
    SET_STRING_ELT(params,14,mkChar("bsvLagDepot"));
    SET_STRING_ELT(params,15,mkChar("bsvLagCentral"));
    SET_STRING_ELT(params,16,mkChar("bsvRateCentral"));
    SET_STRING_ELT(params,17,mkChar("bsvDurCentral"));
  SET_STRING_ELT(lhs,0,mkChar("cp"));
    SET_STRING_ELT(modeln,0,mkChar("normModel"));
    SET_STRING_ELT(model,0,mkChar("popCl=1;\npopV=20;\npopKa=1;\npopVp=10;\npopQ=2;\nbsvCl=0;\nbsvV=0;\nbsvKa=0;\nbsvVp=0;\nbsvQ=0;\ncl~popCl*exp(bsvCl);\nv~popV*exp(bsvV);\nka~popKa*exp(bsvKa);\nq~popQ*exp(bsvQ);\nvp~popVp*exp(bsvVp);\npopLagDepot=0;\npopLagCentral=0;\npopRateCentral=0;\npopDurCentral=0;\nbsvLagDepot=0;\nbsvLagCentral=0;\nbsvRateCentral=0;\nbsvDurCentral=0;\nrx_ka~ka;\nrx_rate~popRateCentral*exp(bsvRateCentral);\nrx_dur~popDurCentral*exp(bsvDurCentral);\nrx_tlag~popLagDepot*exp(bsvLagDepot);\nrx_tlag2~popLagCentral*exp(bsvLagCentral);\nrx_F~1;\nrx_F2~1;\nrx_v~v;\nrx_k~cl/v;\nrx_k12~q/v;\nrx_k21~q/vp;\nrx_beta~0.5*(rx_k12+rx_k21+rx_k-sqrt((rx_k12+rx_k21+rx_k)*(rx_k12+rx_k21+rx_k)-4.0*rx_k21*rx_k));\nrx_alpha~rx_k21*rx_k/rx_beta;\nrx_A~rx_ka/(rx_ka-rx_alpha)*(rx_alpha-rx_k21)/(rx_alpha-rx_beta)/rx_v;\nrx_B~rx_ka/(rx_ka-rx_beta)*(rx_beta-rx_k21)/(rx_beta-rx_alpha)/rx_v;\nrx_A2~(rx_alpha-rx_k21)/(rx_alpha-rx_beta)/rx_v;\nrx_B2~(rx_beta-rx_k21)/(rx_beta-rx_alpha)/rx_v;\nrx_gamma~0;\nrx_C~0;\nrx_C2~0;\ncp=solveLinB(rx__PTR__,t,0,rx_A,rx_A2,rx_alpha,rx_B,rx_B2,rx_beta,rx_C,rx_C2,rx_gamma,rx_ka,rx_tlag,rx_tlag2,rx_F,rx_F2,rx_rate,rx_dur);\n"));
    SEXP ini    = PROTECT(allocVector(REALSXP,18));pro++;
    SEXP inin   = PROTECT(allocVector(STRSXP, 18));pro++;
    SET_STRING_ELT(inin,0,mkChar("popCl"));
    REAL(ini)[0] = 1.0000000000000000;
    SET_STRING_ELT(inin,1,mkChar("popV"));
    REAL(ini)[1] = 20.0000000000000000;
    SET_STRING_ELT(inin,2,mkChar("popKa"));
    REAL(ini)[2] = 1.0000000000000000;
    SET_STRING_ELT(inin,3,mkChar("popVp"));
    REAL(ini)[3] = 10.0000000000000000;
    SET_STRING_ELT(inin,4,mkChar("popQ"));
    REAL(ini)[4] = 2.0000000000000000;
    SET_STRING_ELT(inin,5,mkChar("bsvCl"));
    REAL(ini)[5] = 0.0000000000000000;
    SET_STRING_ELT(inin,6,mkChar("bsvV"));
    REAL(ini)[6] = 0.0000000000000000;
    SET_STRING_ELT(inin,7,mkChar("bsvKa"));
    REAL(ini)[7] = 0.0000000000000000;
    SET_STRING_ELT(inin,8,mkChar("bsvVp"));
    REAL(ini)[8] = 0.0000000000000000;
    SET_STRING_ELT(inin,9,mkChar("bsvQ"));
    REAL(ini)[9] = 0.0000000000000000;
    SET_STRING_ELT(inin,10,mkChar("popLagDepot"));
    REAL(ini)[10] = 0.0000000000000000;
    SET_STRING_ELT(inin,11,mkChar("popLagCentral"));
    REAL(ini)[11] = 0.0000000000000000;
    SET_STRING_ELT(inin,12,mkChar("popRateCentral"));
    REAL(ini)[12] = 0.0000000000000000;
    SET_STRING_ELT(inin,13,mkChar("popDurCentral"));
    REAL(ini)[13] = 0.0000000000000000;
    SET_STRING_ELT(inin,14,mkChar("bsvLagDepot"));
    REAL(ini)[14] = 0.0000000000000000;
    SET_STRING_ELT(inin,15,mkChar("bsvLagCentral"));
    REAL(ini)[15] = 0.0000000000000000;
    SET_STRING_ELT(inin,16,mkChar("bsvRateCentral"));
    REAL(ini)[16] = 0.0000000000000000;
    SET_STRING_ELT(inin,17,mkChar("bsvDurCentral"));
    REAL(ini)[17] = 0.0000000000000000;
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
    SET_STRING_ELT(names,16,mkChar("timeId"));
    SET_VECTOR_ELT(lst,  16,timeInt);
    SET_STRING_ELT(names,17,mkChar("md5"));    SET_VECTOR_ELT(lst,  17,mmd5);    SET_STRING_ELT(mmd5n,0,mkChar("file_md5"));
    SET_STRING_ELT(mmd5,0,mkChar("f546b88ef8df2bef4f58f48b4c03d3cd"));
    SET_STRING_ELT(mmd5n,1,mkChar("parsed_md5"));
    SET_STRING_ELT(mmd5,1,mkChar("f546b88ef8df2bef4f58f48b4c03d3cd"));
    SET_STRING_ELT(trann,0,mkChar("lib.name"));
    SET_STRING_ELT(tran, 0,mkChar("rxModels"));
    SET_STRING_ELT(trann,1,mkChar("jac"));
    SET_STRING_ELT(tran,1,mkChar("fullint"));
    SET_STRING_ELT(trann,2,mkChar("prefix"));
    SET_STRING_ELT(tran, 2,mkChar("rxModels_pk2cmt_"));
    SET_STRING_ELT(trann,3,mkChar("dydt"));
    SET_STRING_ELT(tran, 3,mkChar("rxModels_pk2cmt_dydt"));
    SET_STRING_ELT(trann,4,mkChar("calc_jac"));
    SET_STRING_ELT(tran, 4,mkChar("rxModels_pk2cmt_calc_jac"));
    SET_STRING_ELT(trann,5,mkChar("calc_lhs"));
    SET_STRING_ELT(tran, 5,mkChar("rxModels_pk2cmt_calc_lhs"));
    SET_STRING_ELT(trann,6,mkChar("model_vars"));
    SET_STRING_ELT(tran, 6,mkChar("rxModels_pk2cmt_model_vars"));
    SET_STRING_ELT(trann,7,mkChar("theta"));
    SET_STRING_ELT(tran, 7,mkChar("rxModels_pk2cmt_theta"));
    SET_STRING_ELT(trann,8,mkChar("inis"));
    SET_STRING_ELT(tran, 8,mkChar("rxModels_pk2cmt_inis"));
    SET_STRING_ELT(trann,  9,mkChar("dydt_lsoda"));
    SET_STRING_ELT(tran,   9,mkChar("rxModels_pk2cmt_dydt_lsoda"));
    SET_STRING_ELT(trann,10,mkChar("calc_jac_lsoda"));
    SET_STRING_ELT(tran, 10,mkChar("rxModels_pk2cmt_calc_jac_lsoda"));
    SET_STRING_ELT(trann,11,mkChar("ode_solver_solvedata"));
    SET_STRING_ELT(tran, 11,mkChar("rxModels_pk2cmt_ode_solver_solvedata"));
    SET_STRING_ELT(trann,12,mkChar("ode_solver_get_solvedata"));
    SET_STRING_ELT(tran, 12,mkChar("rxModels_pk2cmt_ode_solver_get_solvedata"));
    SET_STRING_ELT(trann,13,mkChar("dydt_liblsoda"));
    SET_STRING_ELT(tran, 13,mkChar("rxModels_pk2cmt_dydt_liblsoda"));
    SET_STRING_ELT(trann,14,mkChar("F"));
    SET_STRING_ELT(tran, 14,mkChar("rxModels_pk2cmt_F"));
    SET_STRING_ELT(trann,15,mkChar("Lag"));
    SET_STRING_ELT(tran, 15,mkChar("rxModels_pk2cmt_Lag"));
    SET_STRING_ELT(trann,16,mkChar("Rate"));
    SET_STRING_ELT(tran, 16,mkChar("rxModels_pk2cmt_Rate"));
    SET_STRING_ELT(trann,17,mkChar("Dur"));
    SET_STRING_ELT(tran, 17,mkChar("rxModels_pk2cmt_Dur"));
    SET_STRING_ELT(trann,18,mkChar("mtime"));
    SET_STRING_ELT(tran, 18,mkChar("rxModels_pk2cmt_mtime"));
    SET_STRING_ELT(trann,19,mkChar("assignFuns"));
    SET_STRING_ELT(tran, 19,mkChar("rxModels_pk2cmt_assignFuns"));
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
extern void rxModels_pk2cmt_dydt_lsoda(int *neq, double *t, double *A, double *DADT)
{
  rxModels_pk2cmt_dydt(neq, *t, A, DADT);
}
extern int rxModels_pk2cmt_dydt_liblsoda(double t, double *y, double *ydot, void *data)
{
  int *neq = (int*)(data);
  rxModels_pk2cmt_dydt(neq, t, y, ydot);
  return(0);
}
extern void rxModels_pk2cmt_calc_jac_lsoda(int *neq, double *t, double *A,int *ml, int *mu, double *JAC, int *nrowpd){
  // Update all covariate parameters
  rxModels_pk2cmt_calc_jac(neq, *t, A, JAC, *nrowpd);
}

//Create function to call from R's main thread that assigns the required functions. Sometimes they don't get assigned.
extern void rxModels_pk2cmt_assignFuns(){
  _assignFuns();
}

//Initialize the dll to match RxODE's calls
void R_init0_rxModels_pk2cmt(){
  // Get C callables on load; Otherwise it isn't thread safe
  _assignFuns();
  R_RegisterCCallable("rxModels","rxModels_pk2cmt_assignFuns", (DL_FUNC) rxModels_pk2cmt_assignFuns);
  R_RegisterCCallable("rxModels","rxModels_pk2cmt_theta", (DL_FUNC) rxModels_pk2cmt_theta);
  R_RegisterCCallable("rxModels","rxModels_pk2cmt_inis",(DL_FUNC) rxModels_pk2cmt_inis);
  R_RegisterCCallable("rxModels","rxModels_pk2cmt_dydt",(DL_FUNC) rxModels_pk2cmt_dydt);
  R_RegisterCCallable("rxModels","rxModels_pk2cmt_calc_lhs",(DL_FUNC) rxModels_pk2cmt_calc_lhs);
  R_RegisterCCallable("rxModels","rxModels_pk2cmt_calc_jac",(DL_FUNC) rxModels_pk2cmt_calc_jac);
  R_RegisterCCallable("rxModels","rxModels_pk2cmt_dydt_lsoda", (DL_FUNC) rxModels_pk2cmt_dydt_lsoda);
  R_RegisterCCallable("rxModels","rxModels_pk2cmt_calc_jac_lsoda", (DL_FUNC) rxModels_pk2cmt_calc_jac_lsoda);
  R_RegisterCCallable("rxModels","rxModels_pk2cmt_ode_solver_solvedata", (DL_FUNC) rxModels_pk2cmt_ode_solver_solvedata);
  R_RegisterCCallable("rxModels","rxModels_pk2cmt_ode_solver_get_solvedata", (DL_FUNC) rxModels_pk2cmt_ode_solver_get_solvedata);
  R_RegisterCCallable("rxModels","rxModels_pk2cmt_F", (DL_FUNC) rxModels_pk2cmt_F);
  R_RegisterCCallable("rxModels","rxModels_pk2cmt_Lag", (DL_FUNC) rxModels_pk2cmt_Lag);
  R_RegisterCCallable("rxModels","rxModels_pk2cmt_Rate", (DL_FUNC) rxModels_pk2cmt_Rate);
  R_RegisterCCallable("rxModels","rxModels_pk2cmt_Dur", (DL_FUNC) rxModels_pk2cmt_Dur);
  R_RegisterCCallable("rxModels","rxModels_pk2cmt_mtime", (DL_FUNC) rxModels_pk2cmt_mtime);
  R_RegisterCCallable("rxModels","rxModels_pk2cmt_dydt_liblsoda", (DL_FUNC) rxModels_pk2cmt_dydt_liblsoda);
}
//Initialize the dll to match RxODE's calls
void R_init_rxModels_pk2cmt(DllInfo *info){
  // Get C callables on load; Otherwise it isn't thread safe
  R_init0_rxModels_pk2cmt();
  static const R_CallMethodDef callMethods[]  = {
    {"rxModels_pk2cmt_model_vars", (DL_FUNC) &rxModels_pk2cmt_model_vars, 0},
    {NULL, NULL, 0}
  };

  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info,FALSE);
}

void R_unload_rxModels_pk2cmt (DllInfo *info){
  // Free resources required for single subject solve.
  SEXP _mv = PROTECT(_rxGetModelLib("rxModels_pk2cmt_model_vars"));
  if (!isNull(_mv)){
    _rxRmModelLib("rxModels_pk2cmt_model_vars");
  }
  UNPROTECT(1);
}
