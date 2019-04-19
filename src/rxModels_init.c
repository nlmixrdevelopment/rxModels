#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <RxODE.h>
#include <RxODE_model.h>
#include "rxModels_compiled.h"
void R_init_rxModels(DllInfo *info){
  R_init0_rxModels_RxODE_models();
  static const R_CallMethodDef callMethods[]  = {
  compiledModelCall
  {NULL, NULL, 0}
  };
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info,FALSE);
}
void R_unload_rxModels_Jones2013(DllInfo *info);
void R_unload_rxModels_pk1cmt(DllInfo *info);
void R_unload_rxModels_pk2cmt(DllInfo *info);
void R_unload_rxModels_pk3cmt(DllInfo *info);
void R_unload_rxModels(DllInfo *info){
  R_unload_rxModels_Jones2013(info);
  R_unload_rxModels_pk1cmt(info);
  R_unload_rxModels_pk2cmt(info);
  R_unload_rxModels_pk3cmt(info);
}
