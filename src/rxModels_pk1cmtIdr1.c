#include <RxODE.h>
#include <RxODE_model_shared.h>
#define _CENTRAL_ 2
#include <RxODE_model_shared.h>
#define __MAX_PROD__ 0
#define _CMT CMT
#define _SYNC_simeps_ for (int _svari=_solveData->neps; _svari--;){  if (_solveData->op->svar[_svari] == 0) {popCl = _PP[0];};   if (_solveData->op->svar[_svari] == 1) {popV = _PP[1];};   if (_solveData->op->svar[_svari] == 2) {popKa = _PP[2];};   if (_solveData->op->svar[_svari] == 3) {bsvCl = _PP[3];};   if (_solveData->op->svar[_svari] == 4) {bsvV = _PP[4];};   if (_solveData->op->svar[_svari] == 5) {bsvKa = _PP[5];};   if (_solveData->op->svar[_svari] == 6) {popLagDepot = _PP[6];};   if (_solveData->op->svar[_svari] == 7) {popLagCentral = _PP[7];};   if (_solveData->op->svar[_svari] == 8) {popRateCentral = _PP[8];};   if (_solveData->op->svar[_svari] == 9) {popDurCentral = _PP[9];};   if (_solveData->op->svar[_svari] == 10) {bsvLagDepot = _PP[10];};   if (_solveData->op->svar[_svari] == 11) {bsvLagCentral = _PP[11];};   if (_solveData->op->svar[_svari] == 12) {bsvRateCentral = _PP[12];};   if (_solveData->op->svar[_svari] == 13) {bsvDurCentral = _PP[13];};   if (_solveData->op->svar[_svari] == 14) {bsvImax = _PP[14];};   if (_solveData->op->svar[_svari] == 15) {popImax = _PP[15];};   if (_solveData->op->svar[_svari] == 16) {popIc50 = _PP[16];};   if (_solveData->op->svar[_svari] == 17) {bsvIc50 = _PP[17];};   if (_solveData->op->svar[_svari] == 18) {popKin = _PP[18];};   if (_solveData->op->svar[_svari] == 19) {bsvKin = _PP[19];};   if (_solveData->op->svar[_svari] == 20) {popKout = _PP[20];};   if (_solveData->op->svar[_svari] == 21) {bsvKout = _PP[21];};   if (_solveData->op->svar[_svari] == 22) {gamma = _PP[22];}; }
#define _SYNC_simeta_ for (int _ovari=_solveData->neta; _ovari--;){  if (_solveData->op->ovar[_ovari] == 0) {popCl = _PP[0];};   if (_solveData->op->ovar[_ovari] == 1) {popV = _PP[1];};   if (_solveData->op->ovar[_ovari] == 2) {popKa = _PP[2];};   if (_solveData->op->ovar[_ovari] == 3) {bsvCl = _PP[3];};   if (_solveData->op->ovar[_ovari] == 4) {bsvV = _PP[4];};   if (_solveData->op->ovar[_ovari] == 5) {bsvKa = _PP[5];};   if (_solveData->op->ovar[_ovari] == 6) {popLagDepot = _PP[6];};   if (_solveData->op->ovar[_ovari] == 7) {popLagCentral = _PP[7];};   if (_solveData->op->ovar[_ovari] == 8) {popRateCentral = _PP[8];};   if (_solveData->op->ovar[_ovari] == 9) {popDurCentral = _PP[9];};   if (_solveData->op->ovar[_ovari] == 10) {bsvLagDepot = _PP[10];};   if (_solveData->op->ovar[_ovari] == 11) {bsvLagCentral = _PP[11];};   if (_solveData->op->ovar[_ovari] == 12) {bsvRateCentral = _PP[12];};   if (_solveData->op->ovar[_ovari] == 13) {bsvDurCentral = _PP[13];};   if (_solveData->op->ovar[_ovari] == 14) {bsvImax = _PP[14];};   if (_solveData->op->ovar[_ovari] == 15) {popImax = _PP[15];};   if (_solveData->op->ovar[_ovari] == 16) {popIc50 = _PP[16];};   if (_solveData->op->ovar[_ovari] == 17) {bsvIc50 = _PP[17];};   if (_solveData->op->ovar[_ovari] == 18) {popKin = _PP[18];};   if (_solveData->op->ovar[_ovari] == 19) {bsvKin = _PP[19];};   if (_solveData->op->ovar[_ovari] == 20) {popKout = _PP[20];};   if (_solveData->op->ovar[_ovari] == 21) {bsvKout = _PP[21];};   if (_solveData->op->ovar[_ovari] == 22) {gamma = _PP[22];}; }

extern void  rxModels_pk1cmtIdr1_ode_solver_solvedata (rx_solve *solve){
  _solveData = solve;
}
extern rx_solve *rxModels_pk1cmtIdr1_ode_solver_get_solvedata(){
  return _solveData;
}
SEXP rxModels_pk1cmtIdr1_model_vars();


// prj-specific differential eqns
void rxModels_pk1cmtIdr1_dydt(int *_neq, double __t, double *__zzStateVar__, double *__DDtStateVar__)
{
  int _cSub = _neq[1];
  double t = __t + _solveData->subjects[_neq[1]].curShift;
  (void)t;
    double   popCl,
  popV,
  popKa,
  bsvCl,
  bsvV,
  bsvKa,
  cl,
  v,
  ka,
  popLagDepot,
  popLagCentral,
  popRateCentral,
  popDurCentral,
  bsvLagDepot,
  bsvLagCentral,
  bsvRateCentral,
  bsvDurCentral,
  cp,
  bsvImax,
  popImax,
  logitImax,
  Imax,
  popIc50,
  bsvIc50,
  ic50,
  popKin,
  bsvKin,
  kin,
  popKout,
  bsvKout,
  kout,
  gamma,
  R;

  (void)t;
  (void)popCl;
  (void)popV;
  (void)popKa;
  (void)bsvCl;
  (void)bsvV;
  (void)bsvKa;
  (void)cl;
  (void)v;
  (void)ka;
  (void)popLagDepot;
  (void)popLagCentral;
  (void)popRateCentral;
  (void)popDurCentral;
  (void)bsvLagDepot;
  (void)bsvLagCentral;
  (void)bsvRateCentral;
  (void)bsvDurCentral;
  (void)cp;
  (void)bsvImax;
  (void)popImax;
  (void)logitImax;
  (void)Imax;
  (void)popIc50;
  (void)bsvIc50;
  (void)ic50;
  (void)popKin;
  (void)bsvKin;
  (void)kin;
  (void)popKout;
  (void)bsvKout;
  (void)kout;
  (void)gamma;
  (void)R;

  cp = _PL[0];
  Imax = _PL[1];
  ic50 = _PL[2];
  kin = _PL[3];
  kout = _PL[4];

  _update_par_ptr(__t, _cSub, _solveData, _idx);
  popCl = _PP[0];
  popV = _PP[1];
  popKa = _PP[2];
  bsvCl = _PP[3];
  bsvV = _PP[4];
  bsvKa = _PP[5];
  popLagDepot = _PP[6];
  popLagCentral = _PP[7];
  popRateCentral = _PP[8];
  popDurCentral = _PP[9];
  bsvLagDepot = _PP[10];
  bsvLagCentral = _PP[11];
  bsvRateCentral = _PP[12];
  bsvDurCentral = _PP[13];
  bsvImax = _PP[14];
  popImax = _PP[15];
  popIc50 = _PP[16];
  bsvIc50 = _PP[17];
  popKin = _PP[18];
  bsvKin = _PP[19];
  popKout = _PP[20];
  bsvKout = _PP[21];
  gamma = _PP[22];

  R = __zzStateVar__[0]*((double)(_ON[0]));

  cl =popCl*exp(bsvCl);
  v =popV*exp(bsvV);
  ka =popKa*exp(bsvKa);
  cp=linCmtA(_solveData, _cSub,t,1,1,1,cl,v,0.0,0.0,0.0,0.0,popLagDepot*exp(bsvLagDepot),1.0,0.0,0.0,ka,popLagCentral*exp(bsvLagCentral),1.0,popRateCentral*exp(bsvRateCentral),popDurCentral*exp(bsvDurCentral));
  logitImax =-_safe_log(1/safe_zero(popImax)-1)+bsvImax;
  Imax=1/safe_zero((1+exp(-logitImax)));
  ic50=popIc50*exp(bsvIc50);
  kin=popKin*exp(bsvKin);
  kout=popKout*exp(bsvKout);
  __DDtStateVar__[0] = ((double)(_ON[0]))*(_IR[0] + kin*(1-Imax*R_pow(_as_dbleps(cp),gamma)/safe_zero((R_pow(_as_dbleps(ic50),gamma)+R_pow(_as_dbleps(cp),gamma))))-kout*R);
  (&_solveData->subjects[_cSub])->dadt_counter[0]++;
}

// Jacobian derived vars
void rxModels_pk1cmtIdr1_calc_jac(int *_neq, double __t, double *__zzStateVar__, double *__PDStateVar__, unsigned int __NROWPD__) {
  int _cSub=_neq[1];
  double t = __t + _solveData->subjects[_neq[1]].curShift;
  (void)t;
    (&_solveData->subjects[_cSub])->jac_counter[0]++;
}
// Functional based initial conditions.
void rxModels_pk1cmtIdr1_inis(int _cSub, double *__zzStateVar__){
  double t=0;
  double   popCl,
  popV,
  popKa,
  bsvCl,
  bsvV,
  bsvKa,
  cl,
  v,
  ka,
  popLagDepot,
  popLagCentral,
  popRateCentral,
  popDurCentral,
  bsvLagDepot,
  bsvLagCentral,
  bsvRateCentral,
  bsvDurCentral,
  cp,
  bsvImax,
  popImax,
  logitImax,
  Imax,
  popIc50,
  bsvIc50,
  ic50,
  popKin,
  bsvKin,
  kin,
  popKout,
  bsvKout,
  kout,
  gamma,
  R;

  (void)t;
  (void)popCl;
  (void)popV;
  (void)popKa;
  (void)bsvCl;
  (void)bsvV;
  (void)bsvKa;
  (void)cl;
  (void)v;
  (void)ka;
  (void)popLagDepot;
  (void)popLagCentral;
  (void)popRateCentral;
  (void)popDurCentral;
  (void)bsvLagDepot;
  (void)bsvLagCentral;
  (void)bsvRateCentral;
  (void)bsvDurCentral;
  (void)cp;
  (void)bsvImax;
  (void)popImax;
  (void)logitImax;
  (void)Imax;
  (void)popIc50;
  (void)bsvIc50;
  (void)ic50;
  (void)popKin;
  (void)bsvKin;
  (void)kin;
  (void)popKout;
  (void)bsvKout;
  (void)kout;
  (void)gamma;
  (void)R;

  cp = _PL[0];
  Imax = _PL[1];
  ic50 = _PL[2];
  kin = _PL[3];
  kout = _PL[4];

  _update_par_ptr(0.0, _cSub, _solveData, _idx);
  popCl = _PP[0];
  popV = _PP[1];
  popKa = _PP[2];
  bsvCl = _PP[3];
  bsvV = _PP[4];
  bsvKa = _PP[5];
  popLagDepot = _PP[6];
  popLagCentral = _PP[7];
  popRateCentral = _PP[8];
  popDurCentral = _PP[9];
  bsvLagDepot = _PP[10];
  bsvLagCentral = _PP[11];
  bsvRateCentral = _PP[12];
  bsvDurCentral = _PP[13];
  bsvImax = _PP[14];
  popImax = _PP[15];
  popIc50 = _PP[16];
  bsvIc50 = _PP[17];
  popKin = _PP[18];
  bsvKin = _PP[19];
  popKout = _PP[20];
  bsvKout = _PP[21];
  gamma = _PP[22];

  R = __zzStateVar__[0]*((double)(_ON[0]));

  cl =popCl*exp(bsvCl);
  v =popV*exp(bsvV);
  ka =popKa*exp(bsvKa);
  cp=linCmtA(_solveData, _cSub,t,1,1,1,cl,v,0.0,0.0,0.0,0.0,popLagDepot*exp(bsvLagDepot),1.0,0.0,0.0,ka,popLagCentral*exp(bsvLagCentral),1.0,popRateCentral*exp(bsvRateCentral),popDurCentral*exp(bsvDurCentral));
  logitImax =-_safe_log(1/safe_zero(popImax)-1)+bsvImax;
  Imax=1/safe_zero((1+exp(-logitImax)));
  ic50=popIc50*exp(bsvIc50);
  kin=popKin*exp(bsvKin);
  kout=popKout*exp(bsvKout);
  R=kin/safe_zero(kout);
  __zzStateVar__[0]=((double)(_ON[0]))*(R);
}
// prj-specific derived vars
void rxModels_pk1cmtIdr1_calc_lhs(int _cSub, double __t, double *__zzStateVar__, double *_lhs) {
   double t = __t + _solveData->subjects[_cSub].curShift;
  (void)t;
    double   __DDtStateVar_0__,
  popCl,
  popV,
  popKa,
  bsvCl,
  bsvV,
  bsvKa,
  cl,
  v,
  ka,
  popLagDepot,
  popLagCentral,
  popRateCentral,
  popDurCentral,
  bsvLagDepot,
  bsvLagCentral,
  bsvRateCentral,
  bsvDurCentral,
  cp,
  bsvImax,
  popImax,
  logitImax,
  Imax,
  popIc50,
  bsvIc50,
  ic50,
  popKin,
  bsvKin,
  kin,
  popKout,
  bsvKout,
  kout,
  gamma,
  R;

  (void)t;
  (void)__DDtStateVar_0__;
  (void)popCl;
  (void)popV;
  (void)popKa;
  (void)bsvCl;
  (void)bsvV;
  (void)bsvKa;
  (void)cl;
  (void)v;
  (void)ka;
  (void)popLagDepot;
  (void)popLagCentral;
  (void)popRateCentral;
  (void)popDurCentral;
  (void)bsvLagDepot;
  (void)bsvLagCentral;
  (void)bsvRateCentral;
  (void)bsvDurCentral;
  (void)cp;
  (void)bsvImax;
  (void)popImax;
  (void)logitImax;
  (void)Imax;
  (void)popIc50;
  (void)bsvIc50;
  (void)ic50;
  (void)popKin;
  (void)bsvKin;
  (void)kin;
  (void)popKout;
  (void)bsvKout;
  (void)kout;
  (void)gamma;
  (void)R;

  cp = _PL[0];
  Imax = _PL[1];
  ic50 = _PL[2];
  kin = _PL[3];
  kout = _PL[4];

  _update_par_ptr(__t, _cSub, _solveData, _idx);
  popCl = _PP[0];
  popV = _PP[1];
  popKa = _PP[2];
  bsvCl = _PP[3];
  bsvV = _PP[4];
  bsvKa = _PP[5];
  popLagDepot = _PP[6];
  popLagCentral = _PP[7];
  popRateCentral = _PP[8];
  popDurCentral = _PP[9];
  bsvLagDepot = _PP[10];
  bsvLagCentral = _PP[11];
  bsvRateCentral = _PP[12];
  bsvDurCentral = _PP[13];
  bsvImax = _PP[14];
  popImax = _PP[15];
  popIc50 = _PP[16];
  bsvIc50 = _PP[17];
  popKin = _PP[18];
  bsvKin = _PP[19];
  popKout = _PP[20];
  bsvKout = _PP[21];
  gamma = _PP[22];

  R = __zzStateVar__[0]*((double)(_ON[0]));

  cl =popCl*exp(bsvCl);
  v =popV*exp(bsvV);
  ka =popKa*exp(bsvKa);
  cp=linCmtA(_solveData, _cSub,t,1,1,1,cl,v,0.0,0.0,0.0,0.0,popLagDepot*exp(bsvLagDepot),1.0,0.0,0.0,ka,popLagCentral*exp(bsvLagCentral),1.0,popRateCentral*exp(bsvRateCentral),popDurCentral*exp(bsvDurCentral));
  logitImax =-_safe_log(1/safe_zero(popImax)-1)+bsvImax;
  Imax=1/safe_zero((1+exp(-logitImax)));
  ic50=popIc50*exp(bsvIc50);
  kin=popKin*exp(bsvKin);
  kout=popKout*exp(bsvKout);
  __DDtStateVar_0__ = ((double)(_ON[0]))*(_IR[0] + kin*(1-Imax*R_pow(_as_dbleps(cp),gamma)/safe_zero((R_pow(_as_dbleps(ic50),gamma)+R_pow(_as_dbleps(cp),gamma))))-kout*R);

  _lhs[0]=cp;
  _lhs[1]=Imax;
  _lhs[2]=ic50;
  _lhs[3]=kin;
  _lhs[4]=kout;
}
// Functional based bioavailability
double rxModels_pk1cmtIdr1_F(int _cSub,  int _cmt, double _amt, double __t, double *__zzStateVar__){
 return _amt;
}
// Functional based absorption lag
double rxModels_pk1cmtIdr1_Lag(int _cSub,  int _cmt, double __t){
  double *restrict _alag = _solveData->subjects[_cSub].alag;
  (void)_alag; 
  double t = __t + _solveData->subjects[_cSub].curShift;
  (void)t;
    _alag[2]=0.0;
  _alag[1]=0.0;
  _alag[0]=0.0;
  double   popCl,
  popV,
  popKa,
  bsvCl,
  bsvV,
  bsvKa,
  cl,
  v,
  ka,
  popLagDepot,
  popLagCentral,
  popRateCentral,
  popDurCentral,
  bsvLagDepot,
  bsvLagCentral,
  bsvRateCentral,
  bsvDurCentral,
  cp,
  bsvImax,
  popImax,
  logitImax,
  Imax,
  popIc50,
  bsvIc50,
  ic50,
  popKin,
  bsvKin,
  kin,
  popKout,
  bsvKout,
  kout,
  gamma,
  R;

  (void)t;
  (void)popCl;
  (void)popV;
  (void)popKa;
  (void)bsvCl;
  (void)bsvV;
  (void)bsvKa;
  (void)cl;
  (void)v;
  (void)ka;
  (void)popLagDepot;
  (void)popLagCentral;
  (void)popRateCentral;
  (void)popDurCentral;
  (void)bsvLagDepot;
  (void)bsvLagCentral;
  (void)bsvRateCentral;
  (void)bsvDurCentral;
  (void)cp;
  (void)bsvImax;
  (void)popImax;
  (void)logitImax;
  (void)Imax;
  (void)popIc50;
  (void)bsvIc50;
  (void)ic50;
  (void)popKin;
  (void)bsvKin;
  (void)kin;
  (void)popKout;
  (void)bsvKout;
  (void)kout;
  (void)gamma;
  (void)R;

  cp = _PL[0];
  Imax = _PL[1];
  ic50 = _PL[2];
  kin = _PL[3];
  kout = _PL[4];

  _update_par_ptr(NA_REAL, _cSub, _solveData, _idx);
  popCl = _PP[0];
  popV = _PP[1];
  popKa = _PP[2];
  bsvCl = _PP[3];
  bsvV = _PP[4];
  bsvKa = _PP[5];
  popLagDepot = _PP[6];
  popLagCentral = _PP[7];
  popRateCentral = _PP[8];
  popDurCentral = _PP[9];
  bsvLagDepot = _PP[10];
  bsvLagCentral = _PP[11];
  bsvRateCentral = _PP[12];
  bsvDurCentral = _PP[13];
  bsvImax = _PP[14];
  popImax = _PP[15];
  popIc50 = _PP[16];
  bsvIc50 = _PP[17];
  popKin = _PP[18];
  bsvKin = _PP[19];
  popKout = _PP[20];
  bsvKout = _PP[21];
  gamma = _PP[22];

  R = NA_REAL;

  cl =popCl*exp(bsvCl);
  v =popV*exp(bsvV);
  ka =popKa*exp(bsvKa);
  _alag[(&_solveData->subjects[_cSub])->linCmt] = popLagDepot*exp(bsvLagDepot);
  _alag[(&_solveData->subjects[_cSub])->linCmt+1] = popLagCentral*exp(bsvLagCentral);
  logitImax =-_safe_log(1/safe_zero(popImax)-1)+bsvImax;
  Imax=1/safe_zero((1+exp(-logitImax)));
  ic50=popIc50*exp(bsvIc50);
  kin=popKin*exp(bsvKin);
  kout=popKout*exp(bsvKout);

  return t + _alag[_cmt] - _solveData->subjects[_cSub].curShift;
}
// Modeled zero-order rate
double rxModels_pk1cmtIdr1_Rate(int _cSub,  int _cmt, double _amt, double __t){
  double *restrict _rate= _solveData->subjects[_cSub].cRate;
  (void)_rate;
   double t = __t + _solveData->subjects[_cSub].curShift;
  (void)t;
    _rate[2]=0.0;
  _rate[1]=0.0;
  _rate[0]=0.0;
  double   popCl,
  popV,
  popKa,
  bsvCl,
  bsvV,
  bsvKa,
  cl,
  v,
  ka,
  popLagDepot,
  popLagCentral,
  popRateCentral,
  popDurCentral,
  bsvLagDepot,
  bsvLagCentral,
  bsvRateCentral,
  bsvDurCentral,
  cp,
  bsvImax,
  popImax,
  logitImax,
  Imax,
  popIc50,
  bsvIc50,
  ic50,
  popKin,
  bsvKin,
  kin,
  popKout,
  bsvKout,
  kout,
  gamma,
  R;

  (void)t;
  (void)popCl;
  (void)popV;
  (void)popKa;
  (void)bsvCl;
  (void)bsvV;
  (void)bsvKa;
  (void)cl;
  (void)v;
  (void)ka;
  (void)popLagDepot;
  (void)popLagCentral;
  (void)popRateCentral;
  (void)popDurCentral;
  (void)bsvLagDepot;
  (void)bsvLagCentral;
  (void)bsvRateCentral;
  (void)bsvDurCentral;
  (void)cp;
  (void)bsvImax;
  (void)popImax;
  (void)logitImax;
  (void)Imax;
  (void)popIc50;
  (void)bsvIc50;
  (void)ic50;
  (void)popKin;
  (void)bsvKin;
  (void)kin;
  (void)popKout;
  (void)bsvKout;
  (void)kout;
  (void)gamma;
  (void)R;

  cp = _PL[0];
  Imax = _PL[1];
  ic50 = _PL[2];
  kin = _PL[3];
  kout = _PL[4];

  _update_par_ptr(NA_REAL, _cSub, _solveData, _idx);
  popCl = _PP[0];
  popV = _PP[1];
  popKa = _PP[2];
  bsvCl = _PP[3];
  bsvV = _PP[4];
  bsvKa = _PP[5];
  popLagDepot = _PP[6];
  popLagCentral = _PP[7];
  popRateCentral = _PP[8];
  popDurCentral = _PP[9];
  bsvLagDepot = _PP[10];
  bsvLagCentral = _PP[11];
  bsvRateCentral = _PP[12];
  bsvDurCentral = _PP[13];
  bsvImax = _PP[14];
  popImax = _PP[15];
  popIc50 = _PP[16];
  bsvIc50 = _PP[17];
  popKin = _PP[18];
  bsvKin = _PP[19];
  popKout = _PP[20];
  bsvKout = _PP[21];
  gamma = _PP[22];

  R = NA_REAL;

  cl =popCl*exp(bsvCl);
  v =popV*exp(bsvV);
  ka =popKa*exp(bsvKa);
  _rate[(&_solveData->subjects[_cSub])->linCmt+1] = popRateCentral*exp(bsvRateCentral);
  logitImax =-_safe_log(1/safe_zero(popImax)-1)+bsvImax;
  Imax=1/safe_zero((1+exp(-logitImax)));
  ic50=popIc50*exp(bsvIc50);
  kin=popKin*exp(bsvKin);
  kout=popKout*exp(bsvKout);

  return _rate[_cmt];
}
// Modeled zero-order duration
double rxModels_pk1cmtIdr1_Dur(int _cSub,  int _cmt, double _amt, double __t){
  double *restrict _dur = _solveData->subjects[_cSub].cDur;
  (void)_dur;
    double t = __t + _solveData->subjects[_cSub].curShift;
  (void)t;
    _dur[2]=0.0;
  _dur[1]=0.0;
  _dur[0]=0.0;
  double   popCl,
  popV,
  popKa,
  bsvCl,
  bsvV,
  bsvKa,
  cl,
  v,
  ka,
  popLagDepot,
  popLagCentral,
  popRateCentral,
  popDurCentral,
  bsvLagDepot,
  bsvLagCentral,
  bsvRateCentral,
  bsvDurCentral,
  cp,
  bsvImax,
  popImax,
  logitImax,
  Imax,
  popIc50,
  bsvIc50,
  ic50,
  popKin,
  bsvKin,
  kin,
  popKout,
  bsvKout,
  kout,
  gamma,
  R;

  (void)t;
  (void)popCl;
  (void)popV;
  (void)popKa;
  (void)bsvCl;
  (void)bsvV;
  (void)bsvKa;
  (void)cl;
  (void)v;
  (void)ka;
  (void)popLagDepot;
  (void)popLagCentral;
  (void)popRateCentral;
  (void)popDurCentral;
  (void)bsvLagDepot;
  (void)bsvLagCentral;
  (void)bsvRateCentral;
  (void)bsvDurCentral;
  (void)cp;
  (void)bsvImax;
  (void)popImax;
  (void)logitImax;
  (void)Imax;
  (void)popIc50;
  (void)bsvIc50;
  (void)ic50;
  (void)popKin;
  (void)bsvKin;
  (void)kin;
  (void)popKout;
  (void)bsvKout;
  (void)kout;
  (void)gamma;
  (void)R;

  cp = _PL[0];
  Imax = _PL[1];
  ic50 = _PL[2];
  kin = _PL[3];
  kout = _PL[4];

  _update_par_ptr(NA_REAL, _cSub, _solveData, _idx);
  popCl = _PP[0];
  popV = _PP[1];
  popKa = _PP[2];
  bsvCl = _PP[3];
  bsvV = _PP[4];
  bsvKa = _PP[5];
  popLagDepot = _PP[6];
  popLagCentral = _PP[7];
  popRateCentral = _PP[8];
  popDurCentral = _PP[9];
  bsvLagDepot = _PP[10];
  bsvLagCentral = _PP[11];
  bsvRateCentral = _PP[12];
  bsvDurCentral = _PP[13];
  bsvImax = _PP[14];
  popImax = _PP[15];
  popIc50 = _PP[16];
  bsvIc50 = _PP[17];
  popKin = _PP[18];
  bsvKin = _PP[19];
  popKout = _PP[20];
  bsvKout = _PP[21];
  gamma = _PP[22];

  R = NA_REAL;

  cl =popCl*exp(bsvCl);
  v =popV*exp(bsvV);
  ka =popKa*exp(bsvKa);
  _dur[(&_solveData->subjects[_cSub])->linCmt+1] = popDurCentral*exp(bsvDurCentral);
  logitImax =-_safe_log(1/safe_zero(popImax)-1)+bsvImax;
  Imax=1/safe_zero((1+exp(-logitImax)));
  ic50=popIc50*exp(bsvIc50);
  kin=popKin*exp(bsvKin);
  kout=popKout*exp(bsvKout);

  return _dur[_cmt];
}
// Model Times
void rxModels_pk1cmtIdr1_mtime(int _cSub, double *_mtime){
}
// Matrix Exponential (0)
void rxModels_pk1cmtIdr1_ME(int _cSub, double _t, double __t, double *_mat, const double *__zzStateVar__){
   double t = __t + _solveData->subjects[_cSub].curShift;
  (void)t;
  }
// Inductive linearization Matf
void rxModels_pk1cmtIdr1_IndF(int _cSub, double _t, double __t, double *_matf){
   double t = __t + _solveData->subjects[_cSub].curShift;
  (void)t;
  }
extern SEXP rxModels_pk1cmtIdr1_model_vars(){
  int pro=0;
  SEXP _mv = PROTECT(_rxGetModelLib("rxModels_pk1cmtIdr1_model_vars"));pro++;
  if (!_rxIsCurrentC(_mv)){
    SEXP lst      = PROTECT(allocVector(VECSXP, 22));pro++;
    SEXP names    = PROTECT(allocVector(STRSXP, 22));pro++;
    SEXP sNeedSort = PROTECT(allocVector(INTSXP,1));pro++;
    SEXP sLinCmt = PROTECT(allocVector(INTSXP,10));pro++;    INTEGER(sLinCmt)[0]= 1;
    INTEGER(sLinCmt)[1]= 1;
    INTEGER(sLinCmt)[2]= 0;
    INTEGER(sLinCmt)[3]= 0;
    INTEGER(sLinCmt)[4]= 0;
    INTEGER(sLinCmt)[5]= 0;
    INTEGER(sLinCmt)[6]= 1;
    INTEGER(sLinCmt)[7]= 0;
    INTEGER(sLinCmt)[8]= 0;
    INTEGER(sLinCmt)[9]= 0;
    SEXP sLinCmtN = PROTECT(allocVector(STRSXP, 10));pro++;    SET_STRING_ELT(sLinCmtN, 0, mkChar("ncmt"));
    SET_STRING_ELT(sLinCmtN, 1, mkChar("ka"));
    SET_STRING_ELT(sLinCmtN, 2, mkChar("linB"));
    SET_STRING_ELT(sLinCmtN, 3, mkChar("maxeta"));
    SET_STRING_ELT(sLinCmtN, 4, mkChar("maxtheta"));
    SET_STRING_ELT(sLinCmtN, 5, mkChar("hasCmt"));
    SET_STRING_ELT(sLinCmtN, 6, mkChar("linCmt"));
    SET_STRING_ELT(sLinCmtN, 7, mkChar("linCmtFlg"));
    SET_STRING_ELT(sLinCmtN, 8, mkChar("nIndSim"));
    SET_STRING_ELT(sLinCmtN, 9, mkChar("simflg"));
   setAttrib(sLinCmt,   R_NamesSymbol, sLinCmtN);
    int *iNeedSort  = INTEGER(sNeedSort);
    iNeedSort[0] = 14;
    SEXP sMtime = PROTECT(allocVector(INTSXP,1));pro++;
    int *iMtime  = INTEGER(sMtime);
    iMtime[0] = 0;
    SEXP sExtraCmt = PROTECT(allocVector(INTSXP,1));pro++;
    int *iExtraCmt  = INTEGER(sExtraCmt);
    iExtraCmt[0] = 2;
    SEXP params   = PROTECT(allocVector(STRSXP, 23));pro++;
    SEXP lhs      = PROTECT(allocVector(STRSXP, 5));pro++;
    SEXP slhs      = PROTECT(allocVector(STRSXP, 4));pro++;
    SEXP state    = PROTECT(allocVector(STRSXP, 1));pro++;
  SEXP extraState = PROTECT(allocVector(STRSXP, 0));pro++;
    SEXP stateRmS = PROTECT(allocVector(INTSXP, 1));pro++;
    SEXP timeInt = PROTECT(allocVector(INTSXP, 1));pro++;
    INTEGER(timeInt)[0] = 1606108256;
    SEXP sens     = PROTECT(allocVector(STRSXP, 0));pro++;
    SEXP normState= PROTECT(allocVector(STRSXP, 1));pro++;
    SEXP dfdy     = PROTECT(allocVector(STRSXP, 0));pro++;
    SEXP tran     = PROTECT(allocVector(STRSXP, 22));pro++;
    SEXP trann    = PROTECT(allocVector(STRSXP, 22));pro++;
    SEXP mmd5     = PROTECT(allocVector(STRSXP, 2));pro++;
    SEXP mmd5n    = PROTECT(allocVector(STRSXP, 2));pro++;
    SEXP model    = PROTECT(allocVector(STRSXP, 2));pro++;
    SEXP modeln   = PROTECT(allocVector(STRSXP, 2));pro++;
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
    SET_STRING_ELT(params,3,mkChar("bsvCl"));
    SET_STRING_ELT(params,4,mkChar("bsvV"));
    SET_STRING_ELT(params,5,mkChar("bsvKa"));
  SET_STRING_ELT(slhs,0,mkChar("cl"));
  SET_STRING_ELT(slhs,1,mkChar("v"));
  SET_STRING_ELT(slhs,2,mkChar("ka"));
    SET_STRING_ELT(params,6,mkChar("popLagDepot"));
    SET_STRING_ELT(params,7,mkChar("popLagCentral"));
    SET_STRING_ELT(params,8,mkChar("popRateCentral"));
    SET_STRING_ELT(params,9,mkChar("popDurCentral"));
    SET_STRING_ELT(params,10,mkChar("bsvLagDepot"));
    SET_STRING_ELT(params,11,mkChar("bsvLagCentral"));
    SET_STRING_ELT(params,12,mkChar("bsvRateCentral"));
    SET_STRING_ELT(params,13,mkChar("bsvDurCentral"));
  SET_STRING_ELT(lhs,0,mkChar("cp"));
    SET_STRING_ELT(params,14,mkChar("bsvImax"));
    SET_STRING_ELT(params,15,mkChar("popImax"));
  SET_STRING_ELT(slhs,3,mkChar("logitImax"));
  SET_STRING_ELT(lhs,1,mkChar("Imax"));
    SET_STRING_ELT(params,16,mkChar("popIc50"));
    SET_STRING_ELT(params,17,mkChar("bsvIc50"));
  SET_STRING_ELT(lhs,2,mkChar("ic50"));
    SET_STRING_ELT(params,18,mkChar("popKin"));
    SET_STRING_ELT(params,19,mkChar("bsvKin"));
  SET_STRING_ELT(lhs,3,mkChar("kin"));
    SET_STRING_ELT(params,20,mkChar("popKout"));
    SET_STRING_ELT(params,21,mkChar("bsvKout"));
  SET_STRING_ELT(lhs,4,mkChar("kout"));
    SET_STRING_ELT(params,22,mkChar("gamma"));
    SET_STRING_ELT(state,0,mkChar("R"));
    SET_STRING_ELT(normState,0,mkChar("R"));
    _SR[0] = 0;
    SET_STRING_ELT(modeln,0,mkChar("normModel"));
    SET_STRING_ELT(model,0,mkChar("popCl=1;\npopV=20;\npopKa=1;\nbsvCl=0;\nbsvV=0;\nbsvKa=0;\ncl~popCl*exp(bsvCl);\nv~popV*exp(bsvV);\nka~popKa*exp(bsvKa);\npopLagDepot=0;\npopLagCentral=0;\npopRateCentral=0;\npopDurCentral=0;\nbsvLagDepot=0;\nbsvLagCentral=0;\nbsvRateCentral=0;\nbsvDurCentral=0;\ncp=linCmtA(rx__PTR__,t,1,1,1,cl,v,0.0,0.0,0.0,0.0,popLagDepot*exp(bsvLagDepot),1.0,0.0,0.0,ka,popLagCentral*exp(bsvLagCentral),1.0,popRateCentral*exp(bsvRateCentral),popDurCentral*exp(bsvDurCentral));\nbsvImax=0;\npopImax=0.9999;\nlogitImax~-log(1/popImax-1)+bsvImax;\nImax=1/(1+exp(-logitImax));\npopIc50=100;\nbsvIc50=0;\nic50=popIc50*exp(bsvIc50);\npopKin=9;\nbsvKin=0;\nkin=popKin*exp(bsvKin);\npopKout=0.3;\nbsvKout=0;\nkout=popKout*exp(bsvKout);\ngamma=1;\nd/dt(R)=kin*(1-Imax*cp^gamma/(ic50^gamma+cp^gamma))-kout*R;\nR(0)=kin/kout;\n"));
    SET_STRING_ELT(modeln,1,mkChar("indLin"));
    SET_STRING_ELT(model,1,mkChar(""));
    SEXP ini    = PROTECT(allocVector(REALSXP,23));pro++;
    SEXP inin   = PROTECT(allocVector(STRSXP, 23));pro++;
    SET_STRING_ELT(inin,0,mkChar("popCl"));
    REAL(ini)[0] = 1.0000000000000000;
    SET_STRING_ELT(inin,1,mkChar("popV"));
    REAL(ini)[1] = 20.0000000000000000;
    SET_STRING_ELT(inin,2,mkChar("popKa"));
    REAL(ini)[2] = 1.0000000000000000;
    SET_STRING_ELT(inin,3,mkChar("bsvCl"));
    REAL(ini)[3] = 0.0000000000000000;
    SET_STRING_ELT(inin,4,mkChar("bsvV"));
    REAL(ini)[4] = 0.0000000000000000;
    SET_STRING_ELT(inin,5,mkChar("bsvKa"));
    REAL(ini)[5] = 0.0000000000000000;
    SET_STRING_ELT(inin,6,mkChar("popLagDepot"));
    REAL(ini)[6] = 0.0000000000000000;
    SET_STRING_ELT(inin,7,mkChar("popLagCentral"));
    REAL(ini)[7] = 0.0000000000000000;
    SET_STRING_ELT(inin,8,mkChar("popRateCentral"));
    REAL(ini)[8] = 0.0000000000000000;
    SET_STRING_ELT(inin,9,mkChar("popDurCentral"));
    REAL(ini)[9] = 0.0000000000000000;
    SET_STRING_ELT(inin,10,mkChar("bsvLagDepot"));
    REAL(ini)[10] = 0.0000000000000000;
    SET_STRING_ELT(inin,11,mkChar("bsvLagCentral"));
    REAL(ini)[11] = 0.0000000000000000;
    SET_STRING_ELT(inin,12,mkChar("bsvRateCentral"));
    REAL(ini)[12] = 0.0000000000000000;
    SET_STRING_ELT(inin,13,mkChar("bsvDurCentral"));
    REAL(ini)[13] = 0.0000000000000000;
    SET_STRING_ELT(inin,14,mkChar("bsvImax"));
    REAL(ini)[14] = 0.0000000000000000;
    SET_STRING_ELT(inin,15,mkChar("popImax"));
    REAL(ini)[15] = 0.9999000000000000;
    SET_STRING_ELT(inin,16,mkChar("popIc50"));
    REAL(ini)[16] = 100.0000000000000000;
    SET_STRING_ELT(inin,17,mkChar("bsvIc50"));
    REAL(ini)[17] = 0.0000000000000000;
    SET_STRING_ELT(inin,18,mkChar("popKin"));
    REAL(ini)[18] = 9.0000000000000000;
    SET_STRING_ELT(inin,19,mkChar("bsvKin"));
    REAL(ini)[19] = 0.0000000000000000;
    SET_STRING_ELT(inin,20,mkChar("popKout"));
    REAL(ini)[20] = 0.3000000000000000;
    SET_STRING_ELT(inin,21,mkChar("bsvKout"));
    REAL(ini)[21] = 0.0000000000000000;
    SET_STRING_ELT(inin,22,mkChar("gamma"));
    REAL(ini)[22] = 1.0000000000000000;
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
    SET_STRING_ELT(names,9,mkChar("state.ignore"));
    SET_VECTOR_ELT(lst,  9,stateRmS);
    SET_STRING_ELT(names,10,mkChar("version"));
    SET_VECTOR_ELT(lst,  10,version);
    SET_STRING_ELT(names,11,mkChar("normal.state"));
    SET_VECTOR_ELT(lst,  11,normState);
    SET_STRING_ELT(names,12,mkChar("needSort"));
    SET_VECTOR_ELT(lst,  12,sNeedSort);
    SET_STRING_ELT(names,13,mkChar("nMtime"));
    SET_VECTOR_ELT(lst,  13,sMtime);
    SET_STRING_ELT(names,14,mkChar("extraCmt"));
    SET_VECTOR_ELT(lst,  14,sExtraCmt);
    SET_STRING_ELT(names, 15, mkChar("stateExtra"));
    SET_VECTOR_ELT(lst,  15, extraState);
    SET_STRING_ELT(names, 16, mkChar("dvid"));
    SEXP sDvid = PROTECT(allocVector(INTSXP,0));pro++;
    SET_VECTOR_ELT(lst, 16, sDvid);
    SET_STRING_ELT(names,20,mkChar("timeId"));
    SET_VECTOR_ELT(lst,  20,timeInt);
    SET_STRING_ELT(names,17,mkChar("indLin"));
SEXP matLst = PROTECT(allocVector(VECSXP, 0));pro++;
 SET_VECTOR_ELT(lst,  17, matLst);
    SET_STRING_ELT(names,21,mkChar("md5"));    SET_VECTOR_ELT(lst,  21,mmd5);    SET_STRING_ELT(names,18,mkChar("flags"));
    SET_VECTOR_ELT(lst,  18,sLinCmt);
    SET_STRING_ELT(names,19,mkChar("slhs"));
    SET_VECTOR_ELT(lst,  19,slhs);
    SET_STRING_ELT(mmd5n,0,mkChar("file_md5"));
    SET_STRING_ELT(mmd5,0,mkChar(""));
    SET_STRING_ELT(mmd5n,1,mkChar("parsed_md5"));
    SET_STRING_ELT(mmd5,1,mkChar("6a324e3940d4737686fdbda1953ad18c"));
    SET_STRING_ELT(trann,0,mkChar("lib.name"));
    SET_STRING_ELT(tran, 0,mkChar("rxModels"));
    SET_STRING_ELT(trann,1,mkChar("jac"));
    SET_STRING_ELT(tran,1,mkChar("fullint"));
    SET_STRING_ELT(trann,2,mkChar("prefix"));
    SET_STRING_ELT(tran, 2,mkChar("rxModels_pk1cmtIdr1_"));
    SET_STRING_ELT(trann,3,mkChar("dydt"));
    SET_STRING_ELT(tran, 3,mkChar("rxModels_pk1cmtIdr1_dydt"));
    SET_STRING_ELT(trann,4,mkChar("calc_jac"));
    SET_STRING_ELT(tran, 4,mkChar("rxModels_pk1cmtIdr1_calc_jac"));
    SET_STRING_ELT(trann,5,mkChar("calc_lhs"));
    SET_STRING_ELT(tran, 5,mkChar("rxModels_pk1cmtIdr1_calc_lhs"));
    SET_STRING_ELT(trann,6,mkChar("model_vars"));
    SET_STRING_ELT(tran, 6,mkChar("rxModels_pk1cmtIdr1_model_vars"));
    SET_STRING_ELT(trann,7,mkChar("theta"));
    SET_STRING_ELT(tran, 7,mkChar("rxModels_pk1cmtIdr1_theta"));
    SET_STRING_ELT(trann,8,mkChar("inis"));
    SET_STRING_ELT(tran, 8,mkChar("rxModels_pk1cmtIdr1_inis"));
    SET_STRING_ELT(trann,  9,mkChar("dydt_lsoda"));
    SET_STRING_ELT(tran,   9,mkChar("rxModels_pk1cmtIdr1_dydt_lsoda"));
    SET_STRING_ELT(trann,10,mkChar("calc_jac_lsoda"));
    SET_STRING_ELT(tran, 10,mkChar("rxModels_pk1cmtIdr1_calc_jac_lsoda"));
    SET_STRING_ELT(trann,11,mkChar("ode_solver_solvedata"));
    SET_STRING_ELT(tran, 11,mkChar("rxModels_pk1cmtIdr1_ode_solver_solvedata"));
    SET_STRING_ELT(trann,12,mkChar("ode_solver_get_solvedata"));
    SET_STRING_ELT(tran, 12,mkChar("rxModels_pk1cmtIdr1_ode_solver_get_solvedata"));
    SET_STRING_ELT(trann,13,mkChar("dydt_liblsoda"));
    SET_STRING_ELT(tran, 13,mkChar("rxModels_pk1cmtIdr1_dydt_liblsoda"));
    SET_STRING_ELT(trann,14,mkChar("F"));
    SET_STRING_ELT(tran, 14,mkChar("rxModels_pk1cmtIdr1_F"));
    SET_STRING_ELT(trann,15,mkChar("Lag"));
    SET_STRING_ELT(tran, 15,mkChar("rxModels_pk1cmtIdr1_Lag"));
    SET_STRING_ELT(trann,16,mkChar("Rate"));
    SET_STRING_ELT(tran, 16,mkChar("rxModels_pk1cmtIdr1_Rate"));
    SET_STRING_ELT(trann,17,mkChar("Dur"));
    SET_STRING_ELT(tran, 17,mkChar("rxModels_pk1cmtIdr1_Dur"));
    SET_STRING_ELT(trann,18,mkChar("mtime"));
    SET_STRING_ELT(tran, 18,mkChar("rxModels_pk1cmtIdr1_mtime"));
    SET_STRING_ELT(trann,19,mkChar("assignFuns"));
    SET_STRING_ELT(tran, 19,mkChar("rxModels_pk1cmtIdr1_assignFuns"));
    SET_STRING_ELT(trann,20,mkChar("ME"));
    SET_STRING_ELT(tran, 20,mkChar("rxModels_pk1cmtIdr1_ME"));
    SET_STRING_ELT(trann,21,mkChar("IndF"));
    SET_STRING_ELT(tran, 21,mkChar("rxModels_pk1cmtIdr1_IndF"));
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
extern void rxModels_pk1cmtIdr1_dydt_lsoda(int *neq, double *t, double *A, double *DADT)
{
  rxModels_pk1cmtIdr1_dydt(neq, *t, A, DADT);
}
extern int rxModels_pk1cmtIdr1_dydt_liblsoda(double __t, double *y, double *ydot, void *data)
{
  int *neq = (int*)(data);
  rxModels_pk1cmtIdr1_dydt(neq, __t, y, ydot);
  return(0);
}
extern void rxModels_pk1cmtIdr1_calc_jac_lsoda(int *neq, double *t, double *A,int *ml, int *mu, double *JAC, int *nrowpd){
  // Update all covariate parameters
  rxModels_pk1cmtIdr1_calc_jac(neq, *t, A, JAC, *nrowpd);
}

//Create function to call from R's main thread that assigns the required functions. Sometimes they don't get assigned.
extern void rxModels_pk1cmtIdr1_assignFuns(){
  _assignFuns();
}

//Initialize the dll to match RxODE's calls
void R_init0_rxModels_pk1cmtIdr1(){
  // Get C callables on load; Otherwise it isn't thread safe
  _assignFuns();
  R_RegisterCCallable("rxModels","rxModels_pk1cmtIdr1_assignFuns", (DL_FUNC) rxModels_pk1cmtIdr1_assignFuns);
  R_RegisterCCallable("rxModels","rxModels_pk1cmtIdr1_inis",(DL_FUNC) rxModels_pk1cmtIdr1_inis);
  R_RegisterCCallable("rxModels","rxModels_pk1cmtIdr1_dydt",(DL_FUNC) rxModels_pk1cmtIdr1_dydt);
  R_RegisterCCallable("rxModels","rxModels_pk1cmtIdr1_calc_lhs",(DL_FUNC) rxModels_pk1cmtIdr1_calc_lhs);
  R_RegisterCCallable("rxModels","rxModels_pk1cmtIdr1_calc_jac",(DL_FUNC) rxModels_pk1cmtIdr1_calc_jac);
  R_RegisterCCallable("rxModels","rxModels_pk1cmtIdr1_dydt_lsoda", (DL_FUNC) rxModels_pk1cmtIdr1_dydt_lsoda);
  R_RegisterCCallable("rxModels","rxModels_pk1cmtIdr1_calc_jac_lsoda", (DL_FUNC) rxModels_pk1cmtIdr1_calc_jac_lsoda);
  R_RegisterCCallable("rxModels","rxModels_pk1cmtIdr1_ode_solver_solvedata", (DL_FUNC) rxModels_pk1cmtIdr1_ode_solver_solvedata);
  R_RegisterCCallable("rxModels","rxModels_pk1cmtIdr1_ode_solver_get_solvedata", (DL_FUNC) rxModels_pk1cmtIdr1_ode_solver_get_solvedata);
  R_RegisterCCallable("rxModels","rxModels_pk1cmtIdr1_F", (DL_FUNC) rxModels_pk1cmtIdr1_F);
  R_RegisterCCallable("rxModels","rxModels_pk1cmtIdr1_Lag", (DL_FUNC) rxModels_pk1cmtIdr1_Lag);
  R_RegisterCCallable("rxModels","rxModels_pk1cmtIdr1_Rate", (DL_FUNC) rxModels_pk1cmtIdr1_Rate);
  R_RegisterCCallable("rxModels","rxModels_pk1cmtIdr1_Dur", (DL_FUNC) rxModels_pk1cmtIdr1_Dur);
  R_RegisterCCallable("rxModels","rxModels_pk1cmtIdr1_mtime", (DL_FUNC) rxModels_pk1cmtIdr1_mtime);
  R_RegisterCCallable("rxModels","rxModels_pk1cmtIdr1_ME", (DL_FUNC) rxModels_pk1cmtIdr1_ME);
  R_RegisterCCallable("rxModels","rxModels_pk1cmtIdr1_IndF", (DL_FUNC) rxModels_pk1cmtIdr1_IndF);
  R_RegisterCCallable("rxModels","rxModels_pk1cmtIdr1_dydt_liblsoda", (DL_FUNC) rxModels_pk1cmtIdr1_dydt_liblsoda);
}
//Initialize the dll to match RxODE's calls
void R_init_rxModels_pk1cmtIdr1(DllInfo *info){
  // Get C callables on load; Otherwise it isn't thread safe
  R_init0_rxModels_pk1cmtIdr1();
  static const R_CallMethodDef callMethods[]  = {
    {"rxModels_pk1cmtIdr1_model_vars", (DL_FUNC) &rxModels_pk1cmtIdr1_model_vars, 0},
    {NULL, NULL, 0}
  };

  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info,FALSE);
}

void R_unload_rxModels_pk1cmtIdr1 (DllInfo *info){
  // Free resources required for single subject solve.
  SEXP _mv = PROTECT(_rxGetModelLib("rxModels_pk1cmtIdr1_model_vars"));
  if (!isNull(_mv)){
    _rxRmModelLib("rxModels_pk1cmtIdr1_model_vars");
  }
  UNPROTECT(1);
}
