#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <RxODE.h>
#include <RxODE_model_shared.h>
#define compiledModelCall \
{"rxModels_Jones2013_model_vars", (DL_FUNC) &rxModels_Jones2013_model_vars, 0},\
{"rxModels_Ribba2012_model_vars", (DL_FUNC) &rxModels_Ribba2012_model_vars, 0},\
{"rxModels_pk1cmt_model_vars", (DL_FUNC) &rxModels_pk1cmt_model_vars, 0},\
{"rxModels_pk2cmt_model_vars", (DL_FUNC) &rxModels_pk2cmt_model_vars, 0},\
{"rxModels_pk3cmt_model_vars", (DL_FUNC) &rxModels_pk3cmt_model_vars, 0},
SEXP rxModels_Jones2013_model_vars();
SEXP rxModels_Ribba2012_model_vars();
SEXP rxModels_pk1cmt_model_vars();
SEXP rxModels_pk2cmt_model_vars();
SEXP rxModels_pk3cmt_model_vars();
void R_init0_rxModels_Jones2013();
void R_init0_rxModels_Ribba2012();
void R_init0_rxModels_pk1cmt();
void R_init0_rxModels_pk2cmt();
void R_init0_rxModels_pk3cmt();
void R_init0_rxModels_RxODE_models(){
  R_init0_rxModels_Jones2013();
  R_init0_rxModels_Ribba2012();
  R_init0_rxModels_pk1cmt();
  R_init0_rxModels_pk2cmt();
  R_init0_rxModels_pk3cmt();
}
