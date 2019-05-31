#include <RxODE.h>
#include <RxODE_model.h>
#define __MAX_PROD__ 0
#define _CMT CMT
extern void  rxModels_Jones2013_ode_solver_solvedata (rx_solve *solve){
  _solveData = solve;
}
extern rx_solve *rxModels_Jones2013_ode_solver_get_solvedata(){
  return _solveData;
}
SEXP rxModels_Jones2013_model_vars();
double _theta[53];
extern double* rxModels_Jones2013_theta(double *theta){
  _theta[0] = 70.0000000000000000; _theta[1] = 0.2130000000000000; _theta[2] = 0.0856290000000000; _theta[3] = 0.0200000000000000; _theta[4] = 0.0171000000000000; _theta[5] = 0.0047000000000000; _theta[6] = 0.0044000000000000; _theta[7] = 0.0210000000000000; _theta[8] = 0.0076000000000000; _theta[9] = 0.4000000000000000; _theta[10] = 0.0371000000000000; _theta[11] = 0.0026000000000000; _theta[12] = 0.0100000000000000; _theta[13] = 0.0514000000000000; _theta[14] = 0.0257000000000000; _theta[15] = 0.0424000000000000; _theta[16] = 0.0347000000000000; _theta[17] = 0.0997710000000000; _theta[18] = 0.0500000000000000; _theta[19] = 0.0500000000000000; _theta[20] = 0.1200000000000000; _theta[21] = 0.1464620000000000; _theta[22] = 0.0400000000000000; _theta[23] = 0.1900000000000000; _theta[24] = 0.2153850000000000; _theta[25] = 1.0000000000000000; _theta[26] = 0.1700000000000000; _theta[27] = 0.0500000000000000; _theta[28] = 0.0172310000000000; _theta[29] = 0.0107600000000000; _theta[30] = 0.1038550000000000; _theta[31] = 0.1910000000000000; _theta[32] = 0.3740000000000000; _theta[33] = 0.6060000000000000; _theta[34] = 0.5780000000000000; _theta[35] = 0.5830000000000000; _theta[36] = 0.5970000000000000; _theta[37] = 0.5700000000000000; _theta[38] = 0.6200000000000000; _theta[39] = 0.6220000000000000; _theta[40] = 0.6000000000000000; _theta[41] = 0.5910000000000000; _theta[42] = 0.6000000000000000; _theta[43] = 0.6000000000000000; _theta[44] = 0.6810000000000000; _theta[45] = 0.9800000000000000; _theta[46] = 1.0000000000000000; _theta[47] = 8.0000000000000000; _theta[48] = 0.0000000000000000; _theta[49] = 2.1800000000000002; _theta[50] = 1.0000000000000000; _theta[51] = 108.3299999999999983; _theta[52] = 45.0000000000000000;
  return _theta;
}


// prj-specific differential eqns
void rxModels_Jones2013_dydt(int *_neq, double t, double *__zzStateVar__, double *__DDtStateVar__)
{
  int _cSub = _neq[1];
  double   BW,
  FVad,
  FVbo,
  FVbr,
  FVgu,
  FVhe,
  FVki,
  FVli,
  FVlu,
  FVmu,
  FVsk,
  FVsp,
  FVte,
  FVve,
  FVar,
  FVpl,
  FVrb,
  FVre,
  FQad,
  FQbo,
  FQbr,
  FQgu,
  FQhe,
  FQki,
  FQh,
  FQlu,
  FQmu,
  FQsk,
  FQsp,
  FQte,
  FQre,
  Kpad,
  Kpbo,
  Kpbr,
  Kpgu,
  Kphe,
  Kpki,
  Kpli,
  Kplu,
  Kpmu,
  Kpsk,
  Kpsp,
  Kpte,
  Kpre,
  fup,
  BP,
  fumic,
  HLM_CLint,
  CLrenal,
  Ka,
  F,
  CO,
  Vad,
  Vbo,
  Vbr,
  Vgu,
  Vhe,
  Vki,
  Vli,
  Vlu,
  Vmu,
  Vsk,
  Vsp,
  Vte,
  Vve,
  Var,
  Vpl,
  Vrb,
  Vre,
  Vplas_ven,
  Vplas_art,
  QC,
  Qad,
  Qbo,
  Qbr,
  Qgu,
  Qhe,
  Qki,
  Qh,
  Qha,
  Qsp,
  Qlu,
  Qmu,
  Qsk,
  Qte,
  Qre,
  Cadipose,
  Aad,
  Cbone,
  Abo,
  Cbrain,
  Abr,
  Cgut,
  Agu,
  Cheart,
  Ahe,
  Ckidney,
  Aki,
  Cliver,
  Ali,
  Clung,
  Alu,
  Cmuscle,
  Amu,
  Cskin,
  Ask,
  Cspleen,
  Asp,
  Ctestes,
  Ate,
  Cvenous,
  Ave,
  Carterial,
  Aar,
  Crest,
  Are,
  Absorption,
  D,
  Cliverfree,
  Ckidneyfree,
  MPPGL,
  CLmet,
  Venous,
  Cp;

  (void)t;
  (void)BW;
  (void)FVad;
  (void)FVbo;
  (void)FVbr;
  (void)FVgu;
  (void)FVhe;
  (void)FVki;
  (void)FVli;
  (void)FVlu;
  (void)FVmu;
  (void)FVsk;
  (void)FVsp;
  (void)FVte;
  (void)FVve;
  (void)FVar;
  (void)FVpl;
  (void)FVrb;
  (void)FVre;
  (void)FQad;
  (void)FQbo;
  (void)FQbr;
  (void)FQgu;
  (void)FQhe;
  (void)FQki;
  (void)FQh;
  (void)FQlu;
  (void)FQmu;
  (void)FQsk;
  (void)FQsp;
  (void)FQte;
  (void)FQre;
  (void)Kpad;
  (void)Kpbo;
  (void)Kpbr;
  (void)Kpgu;
  (void)Kphe;
  (void)Kpki;
  (void)Kpli;
  (void)Kplu;
  (void)Kpmu;
  (void)Kpsk;
  (void)Kpsp;
  (void)Kpte;
  (void)Kpre;
  (void)fup;
  (void)BP;
  (void)fumic;
  (void)HLM_CLint;
  (void)CLrenal;
  (void)Ka;
  (void)F;
  (void)CO;
  (void)Vad;
  (void)Vbo;
  (void)Vbr;
  (void)Vgu;
  (void)Vhe;
  (void)Vki;
  (void)Vli;
  (void)Vlu;
  (void)Vmu;
  (void)Vsk;
  (void)Vsp;
  (void)Vte;
  (void)Vve;
  (void)Var;
  (void)Vpl;
  (void)Vrb;
  (void)Vre;
  (void)Vplas_ven;
  (void)Vplas_art;
  (void)QC;
  (void)Qad;
  (void)Qbo;
  (void)Qbr;
  (void)Qgu;
  (void)Qhe;
  (void)Qki;
  (void)Qh;
  (void)Qha;
  (void)Qsp;
  (void)Qlu;
  (void)Qmu;
  (void)Qsk;
  (void)Qte;
  (void)Qre;
  (void)Cadipose;
  (void)Aad;
  (void)Cbone;
  (void)Abo;
  (void)Cbrain;
  (void)Abr;
  (void)Cgut;
  (void)Agu;
  (void)Cheart;
  (void)Ahe;
  (void)Ckidney;
  (void)Aki;
  (void)Cliver;
  (void)Ali;
  (void)Clung;
  (void)Alu;
  (void)Cmuscle;
  (void)Amu;
  (void)Cskin;
  (void)Ask;
  (void)Cspleen;
  (void)Asp;
  (void)Ctestes;
  (void)Ate;
  (void)Cvenous;
  (void)Ave;
  (void)Carterial;
  (void)Aar;
  (void)Crest;
  (void)Are;
  (void)Absorption;
  (void)D;
  (void)Cliverfree;
  (void)Ckidneyfree;
  (void)MPPGL;
  (void)CLmet;
  (void)Venous;
  (void)Cp;

  _update_par_ptr(t, _cSub, _solveData, _idx);
  BW = _PP[0];
  FVad = _PP[1];
  FVbo = _PP[2];
  FVbr = _PP[3];
  FVgu = _PP[4];
  FVhe = _PP[5];
  FVki = _PP[6];
  FVli = _PP[7];
  FVlu = _PP[8];
  FVmu = _PP[9];
  FVsk = _PP[10];
  FVsp = _PP[11];
  FVte = _PP[12];
  FVve = _PP[13];
  FVar = _PP[14];
  FVpl = _PP[15];
  FVrb = _PP[16];
  FVre = _PP[17];
  FQad = _PP[18];
  FQbo = _PP[19];
  FQbr = _PP[20];
  FQgu = _PP[21];
  FQhe = _PP[22];
  FQki = _PP[23];
  FQh = _PP[24];
  FQlu = _PP[25];
  FQmu = _PP[26];
  FQsk = _PP[27];
  FQsp = _PP[28];
  FQte = _PP[29];
  FQre = _PP[30];
  Kpad = _PP[31];
  Kpbo = _PP[32];
  Kpbr = _PP[33];
  Kpgu = _PP[34];
  Kphe = _PP[35];
  Kpki = _PP[36];
  Kpli = _PP[37];
  Kplu = _PP[38];
  Kpmu = _PP[39];
  Kpsk = _PP[40];
  Kpsp = _PP[41];
  Kpte = _PP[42];
  Kpre = _PP[43];
  fup = _PP[44];
  BP = _PP[45];
  fumic = _PP[46];
  HLM_CLint = _PP[47];
  CLrenal = _PP[48];
  Ka = _PP[49];
  F = _PP[50];
  CO = _PP[51];
  MPPGL = _PP[52];

  Aad = __zzStateVar__[0]*((double)(_ON[0]));
  Abo = __zzStateVar__[1]*((double)(_ON[1]));
  Abr = __zzStateVar__[2]*((double)(_ON[2]));
  Agu = __zzStateVar__[3]*((double)(_ON[3]));
  Ahe = __zzStateVar__[4]*((double)(_ON[4]));
  Aki = __zzStateVar__[5]*((double)(_ON[5]));
  Ali = __zzStateVar__[6]*((double)(_ON[6]));
  Alu = __zzStateVar__[7]*((double)(_ON[7]));
  Amu = __zzStateVar__[8]*((double)(_ON[8]));
  Ask = __zzStateVar__[9]*((double)(_ON[9]));
  Asp = __zzStateVar__[10]*((double)(_ON[10]));
  Ate = __zzStateVar__[11]*((double)(_ON[11]));
  Ave = __zzStateVar__[12]*((double)(_ON[12]));
  Aar = __zzStateVar__[13]*((double)(_ON[13]));
  Are = __zzStateVar__[14]*((double)(_ON[14]));
  D = __zzStateVar__[15]*((double)(_ON[15]));

  Vad =BW*FVad;
  Vbo =BW*FVbo;
  Vbr =BW*FVbr;
  Vgu =BW*FVgu;
  Vhe =BW*FVhe;
  Vki =BW*FVki;
  Vli =BW*FVli;
  Vlu =BW*FVlu;
  Vmu =BW*FVmu;
  Vsk =BW*FVsk;
  Vsp =BW*FVsp;
  Vte =BW*FVte;
  Vve =BW*FVve;
  Var =BW*FVar;
  Vpl =BW*FVpl;
  Vrb =BW*FVrb;
  Vre =BW*FVre;
  Vplas_ven =Vpl*Vve/safe_zero((Vve+Var));
  Vplas_art =Vpl*Var/safe_zero((Vve+Var));
  QC =CO/safe_zero(1000)*60*60;
  Qad =QC*FQad;
  Qbo =QC*FQbo;
  Qbr =QC*FQbr;
  Qgu =QC*FQgu;
  Qhe =QC*FQhe;
  Qki =QC*FQki;
  Qh =QC*FQh;
  Qha =Qh-Qgu-Qsp;
  Qlu =QC*FQlu;
  Qmu =QC*FQmu;
  Qsk =QC*FQsk;
  Qsp =QC*FQsp;
  Qte =QC*FQte;
  Qre =QC*FQre;
  Cadipose =Aad/safe_zero(Vad);
  Cbone =Abo/safe_zero(Vbo);
  Cbrain =Abr/safe_zero(Vbr);
  Cgut =Agu/safe_zero(Vgu);
  Cheart =Ahe/safe_zero(Vhe);
  Ckidney =Aki/safe_zero(Vki);
  Cliver =Ali/safe_zero(Vli);
  Clung =Alu/safe_zero(Vlu);
  Cmuscle =Amu/safe_zero(Vmu);
  Cskin =Ask/safe_zero(Vsk);
  Cspleen =Asp/safe_zero(Vsp);
  Ctestes =Ate/safe_zero(Vte);
  Cvenous =Ave/safe_zero(Vve);
  Carterial =Aar/safe_zero(Var);
  Crest =Are/safe_zero(Vre);
  Absorption =Ka*D*F;
  Cliverfree =Cliver*fup;
  Ckidneyfree =Ckidney*fup;
  CLmet =(HLM_CLint/safe_zero(fumic))*MPPGL*Vli*60/safe_zero(1000);
  Venous =Qad*(Cadipose/safe_zero(Kpad)*BP)+Qbo*(Cbone/safe_zero(Kpbo)*BP)+Qbr*(Cbrain/safe_zero(Kpbr)*BP)+Qhe*(Cheart/safe_zero(Kphe)*BP)+Qki*(Ckidney/safe_zero(Kpki)*BP)+Qh*(Cliver/safe_zero(Kpli)*BP)+Qmu*(Cmuscle/safe_zero(Kpmu)*BP)+Qsk*(Cskin/safe_zero(Kpsk)*BP)+Qte*(Ctestes/safe_zero(Kpte)*BP)+Qre*(Crest/safe_zero(Kpre)*BP);
  __DDtStateVar__[0] = ((double)(_ON[0]))*(_IR[0] + Qad*(Carterial-Cadipose/safe_zero(Kpad)*BP));
  __DDtStateVar__[1] = ((double)(_ON[1]))*(_IR[1] + Qbo*(Carterial-Cbone/safe_zero(Kpbo)*BP));
  __DDtStateVar__[2] = ((double)(_ON[2]))*(_IR[2] + Qbr*(Carterial-Cbrain/safe_zero(Kpbr)*BP));
  __DDtStateVar__[3] = ((double)(_ON[3]))*(_IR[3] + Absorption+Qgu*(Carterial-Cgut/safe_zero(Kpgu)*BP));
  __DDtStateVar__[4] = ((double)(_ON[4]))*(_IR[4] + Qhe*(Carterial-Cheart/safe_zero(Kphe)*BP));
  __DDtStateVar__[5] = ((double)(_ON[5]))*(_IR[5] + Qki*(Carterial-Ckidney/safe_zero(Kpki)*BP)-CLrenal*Ckidneyfree);
  __DDtStateVar__[6] = ((double)(_ON[6]))*(_IR[6] + Qha*Carterial+Qgu*(Cgut/safe_zero(Kpgu)*BP)+Qsp*(Cspleen/safe_zero(Kpsp)*BP)-Qh*(Cliver/safe_zero(Kpli)*BP)-Cliverfree*CLmet);
  __DDtStateVar__[7] = ((double)(_ON[7]))*(_IR[7] + Qlu*Cvenous-Qlu*(Clung/safe_zero(Kplu)*BP));
  __DDtStateVar__[8] = ((double)(_ON[8]))*(_IR[8] + Qmu*(Carterial-Cmuscle/safe_zero(Kpmu)*BP));
  __DDtStateVar__[9] = ((double)(_ON[9]))*(_IR[9] + Qsk*(Carterial-Cskin/safe_zero(Kpsk)*BP));
  __DDtStateVar__[10] = ((double)(_ON[10]))*(_IR[10] + Qsp*(Carterial-Cspleen/safe_zero(Kpsp)*BP));
  __DDtStateVar__[11] = ((double)(_ON[11]))*(_IR[11] + Qte*(Carterial-Ctestes/safe_zero(Kpte)*BP));
  __DDtStateVar__[12] = ((double)(_ON[12]))*(_IR[12] + Venous-Qlu*Cvenous);
  __DDtStateVar__[13] = ((double)(_ON[13]))*(_IR[13] + Qlu*(Clung/safe_zero(Kplu)*BP)-Qlu*Carterial);
  __DDtStateVar__[14] = ((double)(_ON[14]))*(_IR[14] + Qre*(Carterial-Crest/safe_zero(Kpre)*BP));
  __DDtStateVar__[15] = ((double)(_ON[15]))*(_IR[15] -Absorption);
  Cvenous =Ave/safe_zero(Vve);
  Cp=Cvenous/safe_zero(BP);
  (&_solveData->subjects[_cSub])->dadt_counter[0]++;
}

// Jacobian derived vars
void rxModels_Jones2013_calc_jac(int *_neq, double t, double *__zzStateVar__, double *__PDStateVar__, unsigned int __NROWPD__) {
  int _cSub=_neq[1];
  (&_solveData->subjects[_cSub])->jac_counter[0]++;
}
// Functional based initial conditions.
void rxModels_Jones2013_inis(int _cSub, double *__zzStateVar__){
}
// prj-specific derived vars
void rxModels_Jones2013_calc_lhs(int _cSub, double t, double *__zzStateVar__, double *_lhs) {
  double   __DDtStateVar_0__,
  __DDtStateVar_1__,
  __DDtStateVar_2__,
  __DDtStateVar_3__,
  __DDtStateVar_4__,
  __DDtStateVar_5__,
  __DDtStateVar_6__,
  __DDtStateVar_7__,
  __DDtStateVar_8__,
  __DDtStateVar_9__,
  __DDtStateVar_10__,
  __DDtStateVar_11__,
  __DDtStateVar_12__,
  __DDtStateVar_13__,
  __DDtStateVar_14__,
  __DDtStateVar_15__,
  BW,
  FVad,
  FVbo,
  FVbr,
  FVgu,
  FVhe,
  FVki,
  FVli,
  FVlu,
  FVmu,
  FVsk,
  FVsp,
  FVte,
  FVve,
  FVar,
  FVpl,
  FVrb,
  FVre,
  FQad,
  FQbo,
  FQbr,
  FQgu,
  FQhe,
  FQki,
  FQh,
  FQlu,
  FQmu,
  FQsk,
  FQsp,
  FQte,
  FQre,
  Kpad,
  Kpbo,
  Kpbr,
  Kpgu,
  Kphe,
  Kpki,
  Kpli,
  Kplu,
  Kpmu,
  Kpsk,
  Kpsp,
  Kpte,
  Kpre,
  fup,
  BP,
  fumic,
  HLM_CLint,
  CLrenal,
  Ka,
  F,
  CO,
  Vad,
  Vbo,
  Vbr,
  Vgu,
  Vhe,
  Vki,
  Vli,
  Vlu,
  Vmu,
  Vsk,
  Vsp,
  Vte,
  Vve,
  Var,
  Vpl,
  Vrb,
  Vre,
  Vplas_ven,
  Vplas_art,
  QC,
  Qad,
  Qbo,
  Qbr,
  Qgu,
  Qhe,
  Qki,
  Qh,
  Qha,
  Qsp,
  Qlu,
  Qmu,
  Qsk,
  Qte,
  Qre,
  Cadipose,
  Aad,
  Cbone,
  Abo,
  Cbrain,
  Abr,
  Cgut,
  Agu,
  Cheart,
  Ahe,
  Ckidney,
  Aki,
  Cliver,
  Ali,
  Clung,
  Alu,
  Cmuscle,
  Amu,
  Cskin,
  Ask,
  Cspleen,
  Asp,
  Ctestes,
  Ate,
  Cvenous,
  Ave,
  Carterial,
  Aar,
  Crest,
  Are,
  Absorption,
  D,
  Cliverfree,
  Ckidneyfree,
  MPPGL,
  CLmet,
  Venous,
  Cp;

  (void)t;
  (void)__DDtStateVar_0__;
  (void)__DDtStateVar_1__;
  (void)__DDtStateVar_2__;
  (void)__DDtStateVar_3__;
  (void)__DDtStateVar_4__;
  (void)__DDtStateVar_5__;
  (void)__DDtStateVar_6__;
  (void)__DDtStateVar_7__;
  (void)__DDtStateVar_8__;
  (void)__DDtStateVar_9__;
  (void)__DDtStateVar_10__;
  (void)__DDtStateVar_11__;
  (void)__DDtStateVar_12__;
  (void)__DDtStateVar_13__;
  (void)__DDtStateVar_14__;
  (void)__DDtStateVar_15__;
  (void)BW;
  (void)FVad;
  (void)FVbo;
  (void)FVbr;
  (void)FVgu;
  (void)FVhe;
  (void)FVki;
  (void)FVli;
  (void)FVlu;
  (void)FVmu;
  (void)FVsk;
  (void)FVsp;
  (void)FVte;
  (void)FVve;
  (void)FVar;
  (void)FVpl;
  (void)FVrb;
  (void)FVre;
  (void)FQad;
  (void)FQbo;
  (void)FQbr;
  (void)FQgu;
  (void)FQhe;
  (void)FQki;
  (void)FQh;
  (void)FQlu;
  (void)FQmu;
  (void)FQsk;
  (void)FQsp;
  (void)FQte;
  (void)FQre;
  (void)Kpad;
  (void)Kpbo;
  (void)Kpbr;
  (void)Kpgu;
  (void)Kphe;
  (void)Kpki;
  (void)Kpli;
  (void)Kplu;
  (void)Kpmu;
  (void)Kpsk;
  (void)Kpsp;
  (void)Kpte;
  (void)Kpre;
  (void)fup;
  (void)BP;
  (void)fumic;
  (void)HLM_CLint;
  (void)CLrenal;
  (void)Ka;
  (void)F;
  (void)CO;
  (void)Vad;
  (void)Vbo;
  (void)Vbr;
  (void)Vgu;
  (void)Vhe;
  (void)Vki;
  (void)Vli;
  (void)Vlu;
  (void)Vmu;
  (void)Vsk;
  (void)Vsp;
  (void)Vte;
  (void)Vve;
  (void)Var;
  (void)Vpl;
  (void)Vrb;
  (void)Vre;
  (void)Vplas_ven;
  (void)Vplas_art;
  (void)QC;
  (void)Qad;
  (void)Qbo;
  (void)Qbr;
  (void)Qgu;
  (void)Qhe;
  (void)Qki;
  (void)Qh;
  (void)Qha;
  (void)Qsp;
  (void)Qlu;
  (void)Qmu;
  (void)Qsk;
  (void)Qte;
  (void)Qre;
  (void)Cadipose;
  (void)Aad;
  (void)Cbone;
  (void)Abo;
  (void)Cbrain;
  (void)Abr;
  (void)Cgut;
  (void)Agu;
  (void)Cheart;
  (void)Ahe;
  (void)Ckidney;
  (void)Aki;
  (void)Cliver;
  (void)Ali;
  (void)Clung;
  (void)Alu;
  (void)Cmuscle;
  (void)Amu;
  (void)Cskin;
  (void)Ask;
  (void)Cspleen;
  (void)Asp;
  (void)Ctestes;
  (void)Ate;
  (void)Cvenous;
  (void)Ave;
  (void)Carterial;
  (void)Aar;
  (void)Crest;
  (void)Are;
  (void)Absorption;
  (void)D;
  (void)Cliverfree;
  (void)Ckidneyfree;
  (void)MPPGL;
  (void)CLmet;
  (void)Venous;
  (void)Cp;

  _update_par_ptr(t, _cSub, _solveData, _idx);
  BW = _PP[0];
  FVad = _PP[1];
  FVbo = _PP[2];
  FVbr = _PP[3];
  FVgu = _PP[4];
  FVhe = _PP[5];
  FVki = _PP[6];
  FVli = _PP[7];
  FVlu = _PP[8];
  FVmu = _PP[9];
  FVsk = _PP[10];
  FVsp = _PP[11];
  FVte = _PP[12];
  FVve = _PP[13];
  FVar = _PP[14];
  FVpl = _PP[15];
  FVrb = _PP[16];
  FVre = _PP[17];
  FQad = _PP[18];
  FQbo = _PP[19];
  FQbr = _PP[20];
  FQgu = _PP[21];
  FQhe = _PP[22];
  FQki = _PP[23];
  FQh = _PP[24];
  FQlu = _PP[25];
  FQmu = _PP[26];
  FQsk = _PP[27];
  FQsp = _PP[28];
  FQte = _PP[29];
  FQre = _PP[30];
  Kpad = _PP[31];
  Kpbo = _PP[32];
  Kpbr = _PP[33];
  Kpgu = _PP[34];
  Kphe = _PP[35];
  Kpki = _PP[36];
  Kpli = _PP[37];
  Kplu = _PP[38];
  Kpmu = _PP[39];
  Kpsk = _PP[40];
  Kpsp = _PP[41];
  Kpte = _PP[42];
  Kpre = _PP[43];
  fup = _PP[44];
  BP = _PP[45];
  fumic = _PP[46];
  HLM_CLint = _PP[47];
  CLrenal = _PP[48];
  Ka = _PP[49];
  F = _PP[50];
  CO = _PP[51];
  MPPGL = _PP[52];

  Aad = __zzStateVar__[0]*((double)(_ON[0]));
  Abo = __zzStateVar__[1]*((double)(_ON[1]));
  Abr = __zzStateVar__[2]*((double)(_ON[2]));
  Agu = __zzStateVar__[3]*((double)(_ON[3]));
  Ahe = __zzStateVar__[4]*((double)(_ON[4]));
  Aki = __zzStateVar__[5]*((double)(_ON[5]));
  Ali = __zzStateVar__[6]*((double)(_ON[6]));
  Alu = __zzStateVar__[7]*((double)(_ON[7]));
  Amu = __zzStateVar__[8]*((double)(_ON[8]));
  Ask = __zzStateVar__[9]*((double)(_ON[9]));
  Asp = __zzStateVar__[10]*((double)(_ON[10]));
  Ate = __zzStateVar__[11]*((double)(_ON[11]));
  Ave = __zzStateVar__[12]*((double)(_ON[12]));
  Aar = __zzStateVar__[13]*((double)(_ON[13]));
  Are = __zzStateVar__[14]*((double)(_ON[14]));
  D = __zzStateVar__[15]*((double)(_ON[15]));

  Vad =BW*FVad;
  Vbo =BW*FVbo;
  Vbr =BW*FVbr;
  Vgu =BW*FVgu;
  Vhe =BW*FVhe;
  Vki =BW*FVki;
  Vli =BW*FVli;
  Vlu =BW*FVlu;
  Vmu =BW*FVmu;
  Vsk =BW*FVsk;
  Vsp =BW*FVsp;
  Vte =BW*FVte;
  Vve =BW*FVve;
  Var =BW*FVar;
  Vpl =BW*FVpl;
  Vrb =BW*FVrb;
  Vre =BW*FVre;
  Vplas_ven =Vpl*Vve/safe_zero((Vve+Var));
  Vplas_art =Vpl*Var/safe_zero((Vve+Var));
  QC =CO/safe_zero(1000)*60*60;
  Qad =QC*FQad;
  Qbo =QC*FQbo;
  Qbr =QC*FQbr;
  Qgu =QC*FQgu;
  Qhe =QC*FQhe;
  Qki =QC*FQki;
  Qh =QC*FQh;
  Qha =Qh-Qgu-Qsp;
  Qlu =QC*FQlu;
  Qmu =QC*FQmu;
  Qsk =QC*FQsk;
  Qsp =QC*FQsp;
  Qte =QC*FQte;
  Qre =QC*FQre;
  Cadipose =Aad/safe_zero(Vad);
  Cbone =Abo/safe_zero(Vbo);
  Cbrain =Abr/safe_zero(Vbr);
  Cgut =Agu/safe_zero(Vgu);
  Cheart =Ahe/safe_zero(Vhe);
  Ckidney =Aki/safe_zero(Vki);
  Cliver =Ali/safe_zero(Vli);
  Clung =Alu/safe_zero(Vlu);
  Cmuscle =Amu/safe_zero(Vmu);
  Cskin =Ask/safe_zero(Vsk);
  Cspleen =Asp/safe_zero(Vsp);
  Ctestes =Ate/safe_zero(Vte);
  Cvenous =Ave/safe_zero(Vve);
  Carterial =Aar/safe_zero(Var);
  Crest =Are/safe_zero(Vre);
  Absorption =Ka*D*F;
  Cliverfree =Cliver*fup;
  Ckidneyfree =Ckidney*fup;
  CLmet =(HLM_CLint/safe_zero(fumic))*MPPGL*Vli*60/safe_zero(1000);
  Venous =Qad*(Cadipose/safe_zero(Kpad)*BP)+Qbo*(Cbone/safe_zero(Kpbo)*BP)+Qbr*(Cbrain/safe_zero(Kpbr)*BP)+Qhe*(Cheart/safe_zero(Kphe)*BP)+Qki*(Ckidney/safe_zero(Kpki)*BP)+Qh*(Cliver/safe_zero(Kpli)*BP)+Qmu*(Cmuscle/safe_zero(Kpmu)*BP)+Qsk*(Cskin/safe_zero(Kpsk)*BP)+Qte*(Ctestes/safe_zero(Kpte)*BP)+Qre*(Crest/safe_zero(Kpre)*BP);
  __DDtStateVar_0__ = ((double)(_ON[0]))*(_IR[0] + Qad*(Carterial-Cadipose/safe_zero(Kpad)*BP));
  __DDtStateVar_1__ = ((double)(_ON[1]))*(_IR[1] + Qbo*(Carterial-Cbone/safe_zero(Kpbo)*BP));
  __DDtStateVar_2__ = ((double)(_ON[2]))*(_IR[2] + Qbr*(Carterial-Cbrain/safe_zero(Kpbr)*BP));
  __DDtStateVar_3__ = ((double)(_ON[3]))*(_IR[3] + Absorption+Qgu*(Carterial-Cgut/safe_zero(Kpgu)*BP));
  __DDtStateVar_4__ = ((double)(_ON[4]))*(_IR[4] + Qhe*(Carterial-Cheart/safe_zero(Kphe)*BP));
  __DDtStateVar_5__ = ((double)(_ON[5]))*(_IR[5] + Qki*(Carterial-Ckidney/safe_zero(Kpki)*BP)-CLrenal*Ckidneyfree);
  __DDtStateVar_6__ = ((double)(_ON[6]))*(_IR[6] + Qha*Carterial+Qgu*(Cgut/safe_zero(Kpgu)*BP)+Qsp*(Cspleen/safe_zero(Kpsp)*BP)-Qh*(Cliver/safe_zero(Kpli)*BP)-Cliverfree*CLmet);
  __DDtStateVar_7__ = ((double)(_ON[7]))*(_IR[7] + Qlu*Cvenous-Qlu*(Clung/safe_zero(Kplu)*BP));
  __DDtStateVar_8__ = ((double)(_ON[8]))*(_IR[8] + Qmu*(Carterial-Cmuscle/safe_zero(Kpmu)*BP));
  __DDtStateVar_9__ = ((double)(_ON[9]))*(_IR[9] + Qsk*(Carterial-Cskin/safe_zero(Kpsk)*BP));
  __DDtStateVar_10__ = ((double)(_ON[10]))*(_IR[10] + Qsp*(Carterial-Cspleen/safe_zero(Kpsp)*BP));
  __DDtStateVar_11__ = ((double)(_ON[11]))*(_IR[11] + Qte*(Carterial-Ctestes/safe_zero(Kpte)*BP));
  __DDtStateVar_12__ = ((double)(_ON[12]))*(_IR[12] + Venous-Qlu*Cvenous);
  __DDtStateVar_13__ = ((double)(_ON[13]))*(_IR[13] + Qlu*(Clung/safe_zero(Kplu)*BP)-Qlu*Carterial);
  __DDtStateVar_14__ = ((double)(_ON[14]))*(_IR[14] + Qre*(Carterial-Crest/safe_zero(Kpre)*BP));
  __DDtStateVar_15__ = ((double)(_ON[15]))*(_IR[15] -Absorption);
  Cvenous =Ave/safe_zero(Vve);
  Cp=Cvenous/safe_zero(BP);

  _lhs[0]=Cp;
}
// Functional based bioavailability
double rxModels_Jones2013_F(int _cSub,  int _cmt, double _amt, double t){
 return _amt;
}
// Functional based absorption lag
double rxModels_Jones2013_Lag(int _cSub,  int _cmt, double t){
 return t;
  double   BW,
  FVad,
  FVbo,
  FVbr,
  FVgu,
  FVhe,
  FVki,
  FVli,
  FVlu,
  FVmu,
  FVsk,
  FVsp,
  FVte,
  FVve,
  FVar,
  FVpl,
  FVrb,
  FVre,
  FQad,
  FQbo,
  FQbr,
  FQgu,
  FQhe,
  FQki,
  FQh,
  FQlu,
  FQmu,
  FQsk,
  FQsp,
  FQte,
  FQre,
  Kpad,
  Kpbo,
  Kpbr,
  Kpgu,
  Kphe,
  Kpki,
  Kpli,
  Kplu,
  Kpmu,
  Kpsk,
  Kpsp,
  Kpte,
  Kpre,
  fup,
  BP,
  fumic,
  HLM_CLint,
  CLrenal,
  Ka,
  F,
  CO,
  Vad,
  Vbo,
  Vbr,
  Vgu,
  Vhe,
  Vki,
  Vli,
  Vlu,
  Vmu,
  Vsk,
  Vsp,
  Vte,
  Vve,
  Var,
  Vpl,
  Vrb,
  Vre,
  Vplas_ven,
  Vplas_art,
  QC,
  Qad,
  Qbo,
  Qbr,
  Qgu,
  Qhe,
  Qki,
  Qh,
  Qha,
  Qsp,
  Qlu,
  Qmu,
  Qsk,
  Qte,
  Qre,
  Cadipose,
  Aad,
  Cbone,
  Abo,
  Cbrain,
  Abr,
  Cgut,
  Agu,
  Cheart,
  Ahe,
  Ckidney,
  Aki,
  Cliver,
  Ali,
  Clung,
  Alu,
  Cmuscle,
  Amu,
  Cskin,
  Ask,
  Cspleen,
  Asp,
  Ctestes,
  Ate,
  Cvenous,
  Ave,
  Carterial,
  Aar,
  Crest,
  Are,
  Absorption,
  D,
  Cliverfree,
  Ckidneyfree,
  MPPGL,
  CLmet,
  Venous,
  Cp;

  (void)t;
  (void)BW;
  (void)FVad;
  (void)FVbo;
  (void)FVbr;
  (void)FVgu;
  (void)FVhe;
  (void)FVki;
  (void)FVli;
  (void)FVlu;
  (void)FVmu;
  (void)FVsk;
  (void)FVsp;
  (void)FVte;
  (void)FVve;
  (void)FVar;
  (void)FVpl;
  (void)FVrb;
  (void)FVre;
  (void)FQad;
  (void)FQbo;
  (void)FQbr;
  (void)FQgu;
  (void)FQhe;
  (void)FQki;
  (void)FQh;
  (void)FQlu;
  (void)FQmu;
  (void)FQsk;
  (void)FQsp;
  (void)FQte;
  (void)FQre;
  (void)Kpad;
  (void)Kpbo;
  (void)Kpbr;
  (void)Kpgu;
  (void)Kphe;
  (void)Kpki;
  (void)Kpli;
  (void)Kplu;
  (void)Kpmu;
  (void)Kpsk;
  (void)Kpsp;
  (void)Kpte;
  (void)Kpre;
  (void)fup;
  (void)BP;
  (void)fumic;
  (void)HLM_CLint;
  (void)CLrenal;
  (void)Ka;
  (void)F;
  (void)CO;
  (void)Vad;
  (void)Vbo;
  (void)Vbr;
  (void)Vgu;
  (void)Vhe;
  (void)Vki;
  (void)Vli;
  (void)Vlu;
  (void)Vmu;
  (void)Vsk;
  (void)Vsp;
  (void)Vte;
  (void)Vve;
  (void)Var;
  (void)Vpl;
  (void)Vrb;
  (void)Vre;
  (void)Vplas_ven;
  (void)Vplas_art;
  (void)QC;
  (void)Qad;
  (void)Qbo;
  (void)Qbr;
  (void)Qgu;
  (void)Qhe;
  (void)Qki;
  (void)Qh;
  (void)Qha;
  (void)Qsp;
  (void)Qlu;
  (void)Qmu;
  (void)Qsk;
  (void)Qte;
  (void)Qre;
  (void)Cadipose;
  (void)Aad;
  (void)Cbone;
  (void)Abo;
  (void)Cbrain;
  (void)Abr;
  (void)Cgut;
  (void)Agu;
  (void)Cheart;
  (void)Ahe;
  (void)Ckidney;
  (void)Aki;
  (void)Cliver;
  (void)Ali;
  (void)Clung;
  (void)Alu;
  (void)Cmuscle;
  (void)Amu;
  (void)Cskin;
  (void)Ask;
  (void)Cspleen;
  (void)Asp;
  (void)Ctestes;
  (void)Ate;
  (void)Cvenous;
  (void)Ave;
  (void)Carterial;
  (void)Aar;
  (void)Crest;
  (void)Are;
  (void)Absorption;
  (void)D;
  (void)Cliverfree;
  (void)Ckidneyfree;
  (void)MPPGL;
  (void)CLmet;
  (void)Venous;
  (void)Cp;

  _update_par_ptr(NA_REAL, _cSub, _solveData, _idx);
  BW = _PP[0];
  FVad = _PP[1];
  FVbo = _PP[2];
  FVbr = _PP[3];
  FVgu = _PP[4];
  FVhe = _PP[5];
  FVki = _PP[6];
  FVli = _PP[7];
  FVlu = _PP[8];
  FVmu = _PP[9];
  FVsk = _PP[10];
  FVsp = _PP[11];
  FVte = _PP[12];
  FVve = _PP[13];
  FVar = _PP[14];
  FVpl = _PP[15];
  FVrb = _PP[16];
  FVre = _PP[17];
  FQad = _PP[18];
  FQbo = _PP[19];
  FQbr = _PP[20];
  FQgu = _PP[21];
  FQhe = _PP[22];
  FQki = _PP[23];
  FQh = _PP[24];
  FQlu = _PP[25];
  FQmu = _PP[26];
  FQsk = _PP[27];
  FQsp = _PP[28];
  FQte = _PP[29];
  FQre = _PP[30];
  Kpad = _PP[31];
  Kpbo = _PP[32];
  Kpbr = _PP[33];
  Kpgu = _PP[34];
  Kphe = _PP[35];
  Kpki = _PP[36];
  Kpli = _PP[37];
  Kplu = _PP[38];
  Kpmu = _PP[39];
  Kpsk = _PP[40];
  Kpsp = _PP[41];
  Kpte = _PP[42];
  Kpre = _PP[43];
  fup = _PP[44];
  BP = _PP[45];
  fumic = _PP[46];
  HLM_CLint = _PP[47];
  CLrenal = _PP[48];
  Ka = _PP[49];
  F = _PP[50];
  CO = _PP[51];
  MPPGL = _PP[52];

}
// Modeled zero-order rate
double rxModels_Jones2013_Rate(int _cSub,  int _cmt, double _amt, double t){
 return 0.0;
  double   BW,
  FVad,
  FVbo,
  FVbr,
  FVgu,
  FVhe,
  FVki,
  FVli,
  FVlu,
  FVmu,
  FVsk,
  FVsp,
  FVte,
  FVve,
  FVar,
  FVpl,
  FVrb,
  FVre,
  FQad,
  FQbo,
  FQbr,
  FQgu,
  FQhe,
  FQki,
  FQh,
  FQlu,
  FQmu,
  FQsk,
  FQsp,
  FQte,
  FQre,
  Kpad,
  Kpbo,
  Kpbr,
  Kpgu,
  Kphe,
  Kpki,
  Kpli,
  Kplu,
  Kpmu,
  Kpsk,
  Kpsp,
  Kpte,
  Kpre,
  fup,
  BP,
  fumic,
  HLM_CLint,
  CLrenal,
  Ka,
  F,
  CO,
  Vad,
  Vbo,
  Vbr,
  Vgu,
  Vhe,
  Vki,
  Vli,
  Vlu,
  Vmu,
  Vsk,
  Vsp,
  Vte,
  Vve,
  Var,
  Vpl,
  Vrb,
  Vre,
  Vplas_ven,
  Vplas_art,
  QC,
  Qad,
  Qbo,
  Qbr,
  Qgu,
  Qhe,
  Qki,
  Qh,
  Qha,
  Qsp,
  Qlu,
  Qmu,
  Qsk,
  Qte,
  Qre,
  Cadipose,
  Aad,
  Cbone,
  Abo,
  Cbrain,
  Abr,
  Cgut,
  Agu,
  Cheart,
  Ahe,
  Ckidney,
  Aki,
  Cliver,
  Ali,
  Clung,
  Alu,
  Cmuscle,
  Amu,
  Cskin,
  Ask,
  Cspleen,
  Asp,
  Ctestes,
  Ate,
  Cvenous,
  Ave,
  Carterial,
  Aar,
  Crest,
  Are,
  Absorption,
  D,
  Cliverfree,
  Ckidneyfree,
  MPPGL,
  CLmet,
  Venous,
  Cp;

  (void)t;
  (void)BW;
  (void)FVad;
  (void)FVbo;
  (void)FVbr;
  (void)FVgu;
  (void)FVhe;
  (void)FVki;
  (void)FVli;
  (void)FVlu;
  (void)FVmu;
  (void)FVsk;
  (void)FVsp;
  (void)FVte;
  (void)FVve;
  (void)FVar;
  (void)FVpl;
  (void)FVrb;
  (void)FVre;
  (void)FQad;
  (void)FQbo;
  (void)FQbr;
  (void)FQgu;
  (void)FQhe;
  (void)FQki;
  (void)FQh;
  (void)FQlu;
  (void)FQmu;
  (void)FQsk;
  (void)FQsp;
  (void)FQte;
  (void)FQre;
  (void)Kpad;
  (void)Kpbo;
  (void)Kpbr;
  (void)Kpgu;
  (void)Kphe;
  (void)Kpki;
  (void)Kpli;
  (void)Kplu;
  (void)Kpmu;
  (void)Kpsk;
  (void)Kpsp;
  (void)Kpte;
  (void)Kpre;
  (void)fup;
  (void)BP;
  (void)fumic;
  (void)HLM_CLint;
  (void)CLrenal;
  (void)Ka;
  (void)F;
  (void)CO;
  (void)Vad;
  (void)Vbo;
  (void)Vbr;
  (void)Vgu;
  (void)Vhe;
  (void)Vki;
  (void)Vli;
  (void)Vlu;
  (void)Vmu;
  (void)Vsk;
  (void)Vsp;
  (void)Vte;
  (void)Vve;
  (void)Var;
  (void)Vpl;
  (void)Vrb;
  (void)Vre;
  (void)Vplas_ven;
  (void)Vplas_art;
  (void)QC;
  (void)Qad;
  (void)Qbo;
  (void)Qbr;
  (void)Qgu;
  (void)Qhe;
  (void)Qki;
  (void)Qh;
  (void)Qha;
  (void)Qsp;
  (void)Qlu;
  (void)Qmu;
  (void)Qsk;
  (void)Qte;
  (void)Qre;
  (void)Cadipose;
  (void)Aad;
  (void)Cbone;
  (void)Abo;
  (void)Cbrain;
  (void)Abr;
  (void)Cgut;
  (void)Agu;
  (void)Cheart;
  (void)Ahe;
  (void)Ckidney;
  (void)Aki;
  (void)Cliver;
  (void)Ali;
  (void)Clung;
  (void)Alu;
  (void)Cmuscle;
  (void)Amu;
  (void)Cskin;
  (void)Ask;
  (void)Cspleen;
  (void)Asp;
  (void)Ctestes;
  (void)Ate;
  (void)Cvenous;
  (void)Ave;
  (void)Carterial;
  (void)Aar;
  (void)Crest;
  (void)Are;
  (void)Absorption;
  (void)D;
  (void)Cliverfree;
  (void)Ckidneyfree;
  (void)MPPGL;
  (void)CLmet;
  (void)Venous;
  (void)Cp;

  _update_par_ptr(NA_REAL, _cSub, _solveData, _idx);
  BW = _PP[0];
  FVad = _PP[1];
  FVbo = _PP[2];
  FVbr = _PP[3];
  FVgu = _PP[4];
  FVhe = _PP[5];
  FVki = _PP[6];
  FVli = _PP[7];
  FVlu = _PP[8];
  FVmu = _PP[9];
  FVsk = _PP[10];
  FVsp = _PP[11];
  FVte = _PP[12];
  FVve = _PP[13];
  FVar = _PP[14];
  FVpl = _PP[15];
  FVrb = _PP[16];
  FVre = _PP[17];
  FQad = _PP[18];
  FQbo = _PP[19];
  FQbr = _PP[20];
  FQgu = _PP[21];
  FQhe = _PP[22];
  FQki = _PP[23];
  FQh = _PP[24];
  FQlu = _PP[25];
  FQmu = _PP[26];
  FQsk = _PP[27];
  FQsp = _PP[28];
  FQte = _PP[29];
  FQre = _PP[30];
  Kpad = _PP[31];
  Kpbo = _PP[32];
  Kpbr = _PP[33];
  Kpgu = _PP[34];
  Kphe = _PP[35];
  Kpki = _PP[36];
  Kpli = _PP[37];
  Kplu = _PP[38];
  Kpmu = _PP[39];
  Kpsk = _PP[40];
  Kpsp = _PP[41];
  Kpte = _PP[42];
  Kpre = _PP[43];
  fup = _PP[44];
  BP = _PP[45];
  fumic = _PP[46];
  HLM_CLint = _PP[47];
  CLrenal = _PP[48];
  Ka = _PP[49];
  F = _PP[50];
  CO = _PP[51];
  MPPGL = _PP[52];

}
// Modeled zero-order duration
double rxModels_Jones2013_Dur(int _cSub,  int _cmt, double _amt, double t){
 return 0.0;
}
// Model Times
void rxModels_Jones2013_mtime(int _cSub, double *_mtime){
}
extern SEXP rxModels_Jones2013_model_vars(){
  int pro=0;
  SEXP _mv = PROTECT(_rxGetModelLib("rxModels_Jones2013_model_vars"));pro++;
  if (!_rxIsCurrentC(_mv)){
    SEXP lst      = PROTECT(allocVector(VECSXP, 20));pro++;
    SEXP names    = PROTECT(allocVector(STRSXP, 20));pro++;
    SEXP sNeedSort = PROTECT(allocVector(INTSXP,1));pro++;
    int *iNeedSort  = INTEGER(sNeedSort);
    iNeedSort[0] = 0;
    SEXP sMtime = PROTECT(allocVector(INTSXP,1));pro++;
    int *iMtime  = INTEGER(sMtime);
    iMtime[0] = 0;
    SEXP sExtraCmt = PROTECT(allocVector(INTSXP,1));pro++;
    int *iExtraCmt  = INTEGER(sExtraCmt);
    iExtraCmt[0] = 0;
    SEXP params   = PROTECT(allocVector(STRSXP, 53));pro++;
    SEXP lhs      = PROTECT(allocVector(STRSXP, 1));pro++;
    SEXP state    = PROTECT(allocVector(STRSXP, 16));pro++;
  SEXP extraState = PROTECT(allocVector(STRSXP, 0));pro++;
    SEXP stateRmS = PROTECT(allocVector(INTSXP, 16));pro++;
    SEXP timeInt = PROTECT(allocVector(INTSXP, 1));pro++;
    INTEGER(timeInt)[0] = 1559283223;
    SEXP sens     = PROTECT(allocVector(STRSXP, 0));pro++;
    SEXP normState= PROTECT(allocVector(STRSXP, 16));pro++;
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
    SET_STRING_ELT(params,0,mkChar("BW"));
    SET_STRING_ELT(params,1,mkChar("FVad"));
    SET_STRING_ELT(params,2,mkChar("FVbo"));
    SET_STRING_ELT(params,3,mkChar("FVbr"));
    SET_STRING_ELT(params,4,mkChar("FVgu"));
    SET_STRING_ELT(params,5,mkChar("FVhe"));
    SET_STRING_ELT(params,6,mkChar("FVki"));
    SET_STRING_ELT(params,7,mkChar("FVli"));
    SET_STRING_ELT(params,8,mkChar("FVlu"));
    SET_STRING_ELT(params,9,mkChar("FVmu"));
    SET_STRING_ELT(params,10,mkChar("FVsk"));
    SET_STRING_ELT(params,11,mkChar("FVsp"));
    SET_STRING_ELT(params,12,mkChar("FVte"));
    SET_STRING_ELT(params,13,mkChar("FVve"));
    SET_STRING_ELT(params,14,mkChar("FVar"));
    SET_STRING_ELT(params,15,mkChar("FVpl"));
    SET_STRING_ELT(params,16,mkChar("FVrb"));
    SET_STRING_ELT(params,17,mkChar("FVre"));
    SET_STRING_ELT(params,18,mkChar("FQad"));
    SET_STRING_ELT(params,19,mkChar("FQbo"));
    SET_STRING_ELT(params,20,mkChar("FQbr"));
    SET_STRING_ELT(params,21,mkChar("FQgu"));
    SET_STRING_ELT(params,22,mkChar("FQhe"));
    SET_STRING_ELT(params,23,mkChar("FQki"));
    SET_STRING_ELT(params,24,mkChar("FQh"));
    SET_STRING_ELT(params,25,mkChar("FQlu"));
    SET_STRING_ELT(params,26,mkChar("FQmu"));
    SET_STRING_ELT(params,27,mkChar("FQsk"));
    SET_STRING_ELT(params,28,mkChar("FQsp"));
    SET_STRING_ELT(params,29,mkChar("FQte"));
    SET_STRING_ELT(params,30,mkChar("FQre"));
    SET_STRING_ELT(params,31,mkChar("Kpad"));
    SET_STRING_ELT(params,32,mkChar("Kpbo"));
    SET_STRING_ELT(params,33,mkChar("Kpbr"));
    SET_STRING_ELT(params,34,mkChar("Kpgu"));
    SET_STRING_ELT(params,35,mkChar("Kphe"));
    SET_STRING_ELT(params,36,mkChar("Kpki"));
    SET_STRING_ELT(params,37,mkChar("Kpli"));
    SET_STRING_ELT(params,38,mkChar("Kplu"));
    SET_STRING_ELT(params,39,mkChar("Kpmu"));
    SET_STRING_ELT(params,40,mkChar("Kpsk"));
    SET_STRING_ELT(params,41,mkChar("Kpsp"));
    SET_STRING_ELT(params,42,mkChar("Kpte"));
    SET_STRING_ELT(params,43,mkChar("Kpre"));
    SET_STRING_ELT(params,44,mkChar("fup"));
    SET_STRING_ELT(params,45,mkChar("BP"));
    SET_STRING_ELT(params,46,mkChar("fumic"));
    SET_STRING_ELT(params,47,mkChar("HLM_CLint"));
    SET_STRING_ELT(params,48,mkChar("CLrenal"));
    SET_STRING_ELT(params,49,mkChar("Ka"));
    SET_STRING_ELT(params,50,mkChar("F"));
    SET_STRING_ELT(params,51,mkChar("CO"));
    SET_STRING_ELT(params,52,mkChar("MPPGL"));
  SET_STRING_ELT(lhs,0,mkChar("Cp"));
    SET_STRING_ELT(state,0,mkChar("Aad"));
    SET_STRING_ELT(normState,0,mkChar("Aad"));
    _SR[0] = 1;
    SET_STRING_ELT(state,1,mkChar("Abo"));
    SET_STRING_ELT(normState,1,mkChar("Abo"));
    _SR[1] = 1;
    SET_STRING_ELT(state,2,mkChar("Abr"));
    SET_STRING_ELT(normState,2,mkChar("Abr"));
    _SR[2] = 1;
    SET_STRING_ELT(state,3,mkChar("Agu"));
    SET_STRING_ELT(normState,3,mkChar("Agu"));
    _SR[3] = 1;
    SET_STRING_ELT(state,4,mkChar("Ahe"));
    SET_STRING_ELT(normState,4,mkChar("Ahe"));
    _SR[4] = 1;
    SET_STRING_ELT(state,5,mkChar("Aki"));
    SET_STRING_ELT(normState,5,mkChar("Aki"));
    _SR[5] = 1;
    SET_STRING_ELT(state,6,mkChar("Ali"));
    SET_STRING_ELT(normState,6,mkChar("Ali"));
    _SR[6] = 1;
    SET_STRING_ELT(state,7,mkChar("Alu"));
    SET_STRING_ELT(normState,7,mkChar("Alu"));
    _SR[7] = 1;
    SET_STRING_ELT(state,8,mkChar("Amu"));
    SET_STRING_ELT(normState,8,mkChar("Amu"));
    _SR[8] = 1;
    SET_STRING_ELT(state,9,mkChar("Ask"));
    SET_STRING_ELT(normState,9,mkChar("Ask"));
    _SR[9] = 1;
    SET_STRING_ELT(state,10,mkChar("Asp"));
    SET_STRING_ELT(normState,10,mkChar("Asp"));
    _SR[10] = 1;
    SET_STRING_ELT(state,11,mkChar("Ate"));
    SET_STRING_ELT(normState,11,mkChar("Ate"));
    _SR[11] = 1;
    SET_STRING_ELT(state,12,mkChar("Ave"));
    SET_STRING_ELT(normState,12,mkChar("Ave"));
    _SR[12] = 1;
    SET_STRING_ELT(state,13,mkChar("Aar"));
    SET_STRING_ELT(normState,13,mkChar("Aar"));
    _SR[13] = 1;
    SET_STRING_ELT(state,14,mkChar("Are"));
    SET_STRING_ELT(normState,14,mkChar("Are"));
    _SR[14] = 1;
    SET_STRING_ELT(state,15,mkChar("D"));
    SET_STRING_ELT(normState,15,mkChar("D"));
    _SR[15] = 1;
    SET_STRING_ELT(modeln,0,mkChar("normModel"));
    SET_STRING_ELT(model,0,mkChar("BW=70;\nFVad=0.213;\nFVbo=0.085629;\nFVbr=0.02;\nFVgu=0.0171;\nFVhe=0.0047;\nFVki=0.0044;\nFVli=0.021;\nFVlu=0.0076;\nFVmu=0.4;\nFVsk=0.0371;\nFVsp=0.0026;\nFVte=0.01;\nFVve=0.0514;\nFVar=0.0257;\nFVpl=0.0424;\nFVrb=0.0347;\nFVre=0.099771;\nFQad=0.05;\nFQbo=0.05;\nFQbr=0.12;\nFQgu=0.146462;\nFQhe=0.04;\nFQki=0.19;\nFQh=0.215385;\nFQlu=1;\nFQmu=0.17;\nFQsk=0.05;\nFQsp=0.017231;\nFQte=0.01076;\nFQre=0.103855;\nKpad=0.191;\nKpbo=0.374;\nKpbr=0.606;\nKpgu=0.578;\nKphe=0.583;\nKpki=0.597;\nKpli=0.57;\nKplu=0.62;\nKpmu=0.622;\nKpsk=0.6;\nKpsp=0.591;\nKpte=0.6;\nKpre=0.6;\nfup=0.681;\nBP=0.98;\nfumic=1;\nHLM_CLint=8;\nCLrenal=0;\nKa=2.18;\nF=1;\nCO=108.33;\nVad~BW*FVad;\nVbo~BW*FVbo;\nVbr~BW*FVbr;\nVgu~BW*FVgu;\nVhe~BW*FVhe;\nVki~BW*FVki;\nVli~BW*FVli;\nVlu~BW*FVlu;\nVmu~BW*FVmu;\nVsk~BW*FVsk;\nVsp~BW*FVsp;\nVte~BW*FVte;\nVve~BW*FVve;\nVar~BW*FVar;\nVpl~BW*FVpl;\nVrb~BW*FVrb;\nVre~BW*FVre;\nVplas_ven~Vpl*Vve/(Vve+Var);\nVplas_art~Vpl*Var/(Vve+Var);\nQC~CO/1000*60*60;\nQad~QC*FQad;\nQbo~QC*FQbo;\nQbr~QC*FQbr;\nQgu~QC*FQgu;\nQhe~QC*FQhe;\nQki~QC*FQki;\nQh~QC*FQh;\nQha~Qh-Qgu-Qsp;\nQlu~QC*FQlu;\nQmu~QC*FQmu;\nQsk~QC*FQsk;\nQsp~QC*FQsp;\nQte~QC*FQte;\nQre~QC*FQre;\nCadipose~Aad/Vad;\nCbone~Abo/Vbo;\nCbrain~Abr/Vbr;\nCgut~Agu/Vgu;\nCheart~Ahe/Vhe;\nCkidney~Aki/Vki;\nCliver~Ali/Vli;\nClung~Alu/Vlu;\nCmuscle~Amu/Vmu;\nCskin~Ask/Vsk;\nCspleen~Asp/Vsp;\nCtestes~Ate/Vte;\nCvenous~Ave/Vve;\nCarterial~Aar/Var;\nCrest~Are/Vre;\nAbsorption~Ka*D*F;\nCliverfree~Cliver*fup;\nCkidneyfree~Ckidney*fup;\nMPPGL=45;\nCLmet~(HLM_CLint/fumic)*MPPGL*Vli*60/1000;\nVenous~Qad*(Cadipose/Kpad*BP)+Qbo*(Cbone/Kpbo*BP)+Qbr*(Cbrain/Kpbr*BP)+Qhe*(Cheart/Kphe*BP)+Qki*(Ckidney/Kpki*BP)+Qh*(Cliver/Kpli*BP)+Qmu*(Cmuscle/Kpmu*BP)+Qsk*(Cskin/Kpsk*BP)+Qte*(Ctestes/Kpte*BP)+Qre*(Crest/Kpre*BP);\nd/dt(Aad)~Qad*(Carterial-Cadipose/Kpad*BP);\nd/dt(Abo)~Qbo*(Carterial-Cbone/Kpbo*BP);\nd/dt(Abr)~Qbr*(Carterial-Cbrain/Kpbr*BP);\nd/dt(Agu)~Absorption+Qgu*(Carterial-Cgut/Kpgu*BP);\nd/dt(Ahe)~Qhe*(Carterial-Cheart/Kphe*BP);\nd/dt(Aki)~Qki*(Carterial-Ckidney/Kpki*BP)-CLrenal*Ckidneyfree;\nd/dt(Ali)~Qha*Carterial+Qgu*(Cgut/Kpgu*BP)+Qsp*(Cspleen/Kpsp*BP)-Qh*(Cliver/Kpli*BP)-Cliverfree*CLmet;\nd/dt(Alu)~Qlu*Cvenous-Qlu*(Clung/Kplu*BP);\nd/dt(Amu)~Qmu*(Carterial-Cmuscle/Kpmu*BP);\nd/dt(Ask)~Qsk*(Carterial-Cskin/Kpsk*BP);\nd/dt(Asp)~Qsp*(Carterial-Cspleen/Kpsp*BP);\nd/dt(Ate)~Qte*(Carterial-Ctestes/Kpte*BP);\nd/dt(Ave)~Venous-Qlu*Cvenous;\nd/dt(Aar)~Qlu*(Clung/Kplu*BP)-Qlu*Carterial;\nd/dt(Are)~Qre*(Carterial-Crest/Kpre*BP);\nd/dt(D)~-Absorption;\nCvenous~Ave/Vve;\nCp=Cvenous/BP;\n"));
    SEXP ini    = PROTECT(allocVector(REALSXP,53));pro++;
    SEXP inin   = PROTECT(allocVector(STRSXP, 53));pro++;
    SET_STRING_ELT(inin,0,mkChar("BW"));
    REAL(ini)[0] = 70.0000000000000000;
    SET_STRING_ELT(inin,1,mkChar("FVad"));
    REAL(ini)[1] = 0.2130000000000000;
    SET_STRING_ELT(inin,2,mkChar("FVbo"));
    REAL(ini)[2] = 0.0856290000000000;
    SET_STRING_ELT(inin,3,mkChar("FVbr"));
    REAL(ini)[3] = 0.0200000000000000;
    SET_STRING_ELT(inin,4,mkChar("FVgu"));
    REAL(ini)[4] = 0.0171000000000000;
    SET_STRING_ELT(inin,5,mkChar("FVhe"));
    REAL(ini)[5] = 0.0047000000000000;
    SET_STRING_ELT(inin,6,mkChar("FVki"));
    REAL(ini)[6] = 0.0044000000000000;
    SET_STRING_ELT(inin,7,mkChar("FVli"));
    REAL(ini)[7] = 0.0210000000000000;
    SET_STRING_ELT(inin,8,mkChar("FVlu"));
    REAL(ini)[8] = 0.0076000000000000;
    SET_STRING_ELT(inin,9,mkChar("FVmu"));
    REAL(ini)[9] = 0.4000000000000000;
    SET_STRING_ELT(inin,10,mkChar("FVsk"));
    REAL(ini)[10] = 0.0371000000000000;
    SET_STRING_ELT(inin,11,mkChar("FVsp"));
    REAL(ini)[11] = 0.0026000000000000;
    SET_STRING_ELT(inin,12,mkChar("FVte"));
    REAL(ini)[12] = 0.0100000000000000;
    SET_STRING_ELT(inin,13,mkChar("FVve"));
    REAL(ini)[13] = 0.0514000000000000;
    SET_STRING_ELT(inin,14,mkChar("FVar"));
    REAL(ini)[14] = 0.0257000000000000;
    SET_STRING_ELT(inin,15,mkChar("FVpl"));
    REAL(ini)[15] = 0.0424000000000000;
    SET_STRING_ELT(inin,16,mkChar("FVrb"));
    REAL(ini)[16] = 0.0347000000000000;
    SET_STRING_ELT(inin,17,mkChar("FVre"));
    REAL(ini)[17] = 0.0997710000000000;
    SET_STRING_ELT(inin,18,mkChar("FQad"));
    REAL(ini)[18] = 0.0500000000000000;
    SET_STRING_ELT(inin,19,mkChar("FQbo"));
    REAL(ini)[19] = 0.0500000000000000;
    SET_STRING_ELT(inin,20,mkChar("FQbr"));
    REAL(ini)[20] = 0.1200000000000000;
    SET_STRING_ELT(inin,21,mkChar("FQgu"));
    REAL(ini)[21] = 0.1464620000000000;
    SET_STRING_ELT(inin,22,mkChar("FQhe"));
    REAL(ini)[22] = 0.0400000000000000;
    SET_STRING_ELT(inin,23,mkChar("FQki"));
    REAL(ini)[23] = 0.1900000000000000;
    SET_STRING_ELT(inin,24,mkChar("FQh"));
    REAL(ini)[24] = 0.2153850000000000;
    SET_STRING_ELT(inin,25,mkChar("FQlu"));
    REAL(ini)[25] = 1.0000000000000000;
    SET_STRING_ELT(inin,26,mkChar("FQmu"));
    REAL(ini)[26] = 0.1700000000000000;
    SET_STRING_ELT(inin,27,mkChar("FQsk"));
    REAL(ini)[27] = 0.0500000000000000;
    SET_STRING_ELT(inin,28,mkChar("FQsp"));
    REAL(ini)[28] = 0.0172310000000000;
    SET_STRING_ELT(inin,29,mkChar("FQte"));
    REAL(ini)[29] = 0.0107600000000000;
    SET_STRING_ELT(inin,30,mkChar("FQre"));
    REAL(ini)[30] = 0.1038550000000000;
    SET_STRING_ELT(inin,31,mkChar("Kpad"));
    REAL(ini)[31] = 0.1910000000000000;
    SET_STRING_ELT(inin,32,mkChar("Kpbo"));
    REAL(ini)[32] = 0.3740000000000000;
    SET_STRING_ELT(inin,33,mkChar("Kpbr"));
    REAL(ini)[33] = 0.6060000000000000;
    SET_STRING_ELT(inin,34,mkChar("Kpgu"));
    REAL(ini)[34] = 0.5780000000000000;
    SET_STRING_ELT(inin,35,mkChar("Kphe"));
    REAL(ini)[35] = 0.5830000000000000;
    SET_STRING_ELT(inin,36,mkChar("Kpki"));
    REAL(ini)[36] = 0.5970000000000000;
    SET_STRING_ELT(inin,37,mkChar("Kpli"));
    REAL(ini)[37] = 0.5700000000000000;
    SET_STRING_ELT(inin,38,mkChar("Kplu"));
    REAL(ini)[38] = 0.6200000000000000;
    SET_STRING_ELT(inin,39,mkChar("Kpmu"));
    REAL(ini)[39] = 0.6220000000000000;
    SET_STRING_ELT(inin,40,mkChar("Kpsk"));
    REAL(ini)[40] = 0.6000000000000000;
    SET_STRING_ELT(inin,41,mkChar("Kpsp"));
    REAL(ini)[41] = 0.5910000000000000;
    SET_STRING_ELT(inin,42,mkChar("Kpte"));
    REAL(ini)[42] = 0.6000000000000000;
    SET_STRING_ELT(inin,43,mkChar("Kpre"));
    REAL(ini)[43] = 0.6000000000000000;
    SET_STRING_ELT(inin,44,mkChar("fup"));
    REAL(ini)[44] = 0.6810000000000001;
    SET_STRING_ELT(inin,45,mkChar("BP"));
    REAL(ini)[45] = 0.9800000000000000;
    SET_STRING_ELT(inin,46,mkChar("fumic"));
    REAL(ini)[46] = 1.0000000000000000;
    SET_STRING_ELT(inin,47,mkChar("HLM_CLint"));
    REAL(ini)[47] = 8.0000000000000000;
    SET_STRING_ELT(inin,48,mkChar("CLrenal"));
    REAL(ini)[48] = 0.0000000000000000;
    SET_STRING_ELT(inin,49,mkChar("Ka"));
    REAL(ini)[49] = 2.1800000000000002;
    SET_STRING_ELT(inin,50,mkChar("F"));
    REAL(ini)[50] = 1.0000000000000000;
    SET_STRING_ELT(inin,51,mkChar("CO"));
    REAL(ini)[51] = 108.3300000000000000;
    SET_STRING_ELT(inin,52,mkChar("MPPGL"));
    REAL(ini)[52] = 45.0000000000000000;
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
    SET_STRING_ELT(mmd5,0,mkChar("21cda00d89b1b5cc6fc1c519135f6d7a"));
    SET_STRING_ELT(mmd5n,1,mkChar("parsed_md5"));
    SET_STRING_ELT(mmd5,1,mkChar("21cda00d89b1b5cc6fc1c519135f6d7a"));
    SET_STRING_ELT(trann,0,mkChar("lib.name"));
    SET_STRING_ELT(tran, 0,mkChar("rxModels"));
    SET_STRING_ELT(trann,1,mkChar("jac"));
    SET_STRING_ELT(tran,1,mkChar("fullint"));
    SET_STRING_ELT(trann,2,mkChar("prefix"));
    SET_STRING_ELT(tran, 2,mkChar("rxModels_Jones2013_"));
    SET_STRING_ELT(trann,3,mkChar("dydt"));
    SET_STRING_ELT(tran, 3,mkChar("rxModels_Jones2013_dydt"));
    SET_STRING_ELT(trann,4,mkChar("calc_jac"));
    SET_STRING_ELT(tran, 4,mkChar("rxModels_Jones2013_calc_jac"));
    SET_STRING_ELT(trann,5,mkChar("calc_lhs"));
    SET_STRING_ELT(tran, 5,mkChar("rxModels_Jones2013_calc_lhs"));
    SET_STRING_ELT(trann,6,mkChar("model_vars"));
    SET_STRING_ELT(tran, 6,mkChar("rxModels_Jones2013_model_vars"));
    SET_STRING_ELT(trann,7,mkChar("theta"));
    SET_STRING_ELT(tran, 7,mkChar("rxModels_Jones2013_theta"));
    SET_STRING_ELT(trann,8,mkChar("inis"));
    SET_STRING_ELT(tran, 8,mkChar("rxModels_Jones2013_inis"));
    SET_STRING_ELT(trann,  9,mkChar("dydt_lsoda"));
    SET_STRING_ELT(tran,   9,mkChar("rxModels_Jones2013_dydt_lsoda"));
    SET_STRING_ELT(trann,10,mkChar("calc_jac_lsoda"));
    SET_STRING_ELT(tran, 10,mkChar("rxModels_Jones2013_calc_jac_lsoda"));
    SET_STRING_ELT(trann,11,mkChar("ode_solver_solvedata"));
    SET_STRING_ELT(tran, 11,mkChar("rxModels_Jones2013_ode_solver_solvedata"));
    SET_STRING_ELT(trann,12,mkChar("ode_solver_get_solvedata"));
    SET_STRING_ELT(tran, 12,mkChar("rxModels_Jones2013_ode_solver_get_solvedata"));
    SET_STRING_ELT(trann,13,mkChar("dydt_liblsoda"));
    SET_STRING_ELT(tran, 13,mkChar("rxModels_Jones2013_dydt_liblsoda"));
    SET_STRING_ELT(trann,14,mkChar("F"));
    SET_STRING_ELT(tran, 14,mkChar("rxModels_Jones2013_F"));
    SET_STRING_ELT(trann,15,mkChar("Lag"));
    SET_STRING_ELT(tran, 15,mkChar("rxModels_Jones2013_Lag"));
    SET_STRING_ELT(trann,16,mkChar("Rate"));
    SET_STRING_ELT(tran, 16,mkChar("rxModels_Jones2013_Rate"));
    SET_STRING_ELT(trann,17,mkChar("Dur"));
    SET_STRING_ELT(tran, 17,mkChar("rxModels_Jones2013_Dur"));
    SET_STRING_ELT(trann,18,mkChar("mtime"));
    SET_STRING_ELT(tran, 18,mkChar("rxModels_Jones2013_mtime"));
    SET_STRING_ELT(trann,19,mkChar("assignFuns"));
    SET_STRING_ELT(tran, 19,mkChar("rxModels_Jones2013_assignFuns"));
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
extern void rxModels_Jones2013_dydt_lsoda(int *neq, double *t, double *A, double *DADT)
{
  rxModels_Jones2013_dydt(neq, *t, A, DADT);
}
extern int rxModels_Jones2013_dydt_liblsoda(double t, double *y, double *ydot, void *data)
{
  int *neq = (int*)(data);
  rxModels_Jones2013_dydt(neq, t, y, ydot);
  return(0);
}
extern void rxModels_Jones2013_calc_jac_lsoda(int *neq, double *t, double *A,int *ml, int *mu, double *JAC, int *nrowpd){
  // Update all covariate parameters
  rxModels_Jones2013_calc_jac(neq, *t, A, JAC, *nrowpd);
}

//Create function to call from R's main thread that assigns the required functions. Sometimes they don't get assigned.
extern void rxModels_Jones2013_assignFuns(){
  _assignFuns();
}

//Initialize the dll to match RxODE's calls
void R_init0_rxModels_Jones2013(){
  // Get C callables on load; Otherwise it isn't thread safe
  _assignFuns();
  R_RegisterCCallable("rxModels","rxModels_Jones2013_assignFuns", (DL_FUNC) rxModels_Jones2013_assignFuns);
  R_RegisterCCallable("rxModels","rxModels_Jones2013_theta", (DL_FUNC) rxModels_Jones2013_theta);
  R_RegisterCCallable("rxModels","rxModels_Jones2013_inis",(DL_FUNC) rxModels_Jones2013_inis);
  R_RegisterCCallable("rxModels","rxModels_Jones2013_dydt",(DL_FUNC) rxModels_Jones2013_dydt);
  R_RegisterCCallable("rxModels","rxModels_Jones2013_calc_lhs",(DL_FUNC) rxModels_Jones2013_calc_lhs);
  R_RegisterCCallable("rxModels","rxModels_Jones2013_calc_jac",(DL_FUNC) rxModels_Jones2013_calc_jac);
  R_RegisterCCallable("rxModels","rxModels_Jones2013_dydt_lsoda", (DL_FUNC) rxModels_Jones2013_dydt_lsoda);
  R_RegisterCCallable("rxModels","rxModels_Jones2013_calc_jac_lsoda", (DL_FUNC) rxModels_Jones2013_calc_jac_lsoda);
  R_RegisterCCallable("rxModels","rxModels_Jones2013_ode_solver_solvedata", (DL_FUNC) rxModels_Jones2013_ode_solver_solvedata);
  R_RegisterCCallable("rxModels","rxModels_Jones2013_ode_solver_get_solvedata", (DL_FUNC) rxModels_Jones2013_ode_solver_get_solvedata);
  R_RegisterCCallable("rxModels","rxModels_Jones2013_F", (DL_FUNC) rxModels_Jones2013_F);
  R_RegisterCCallable("rxModels","rxModels_Jones2013_Lag", (DL_FUNC) rxModels_Jones2013_Lag);
  R_RegisterCCallable("rxModels","rxModels_Jones2013_Rate", (DL_FUNC) rxModels_Jones2013_Rate);
  R_RegisterCCallable("rxModels","rxModels_Jones2013_Dur", (DL_FUNC) rxModels_Jones2013_Dur);
  R_RegisterCCallable("rxModels","rxModels_Jones2013_mtime", (DL_FUNC) rxModels_Jones2013_mtime);
  R_RegisterCCallable("rxModels","rxModels_Jones2013_dydt_liblsoda", (DL_FUNC) rxModels_Jones2013_dydt_liblsoda);
}
//Initialize the dll to match RxODE's calls
void R_init_rxModels_Jones2013(DllInfo *info){
  // Get C callables on load; Otherwise it isn't thread safe
  R_init0_rxModels_Jones2013();
  static const R_CallMethodDef callMethods[]  = {
    {"rxModels_Jones2013_model_vars", (DL_FUNC) &rxModels_Jones2013_model_vars, 0},
    {NULL, NULL, 0}
  };

  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info,FALSE);
}

void R_unload_rxModels_Jones2013 (DllInfo *info){
  // Free resources required for single subject solve.
  SEXP _mv = PROTECT(_rxGetModelLib("rxModels_Jones2013_model_vars"));
  if (!isNull(_mv)){
    _rxRmModelLib("rxModels_Jones2013_model_vars");
  }
  UNPROTECT(1);
}
