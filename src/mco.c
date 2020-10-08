#include "extern.h"

static const R_CallMethodDef callMethods[]  = {
  {"do_eps_ind", (DL_FUNC) &do_eps_ind, 2},
  {"do_hv", (DL_FUNC) &do_hv, 2},
  {"do_nsga2", (DL_FUNC)&do_nsga2, 14},
  {NULL, NULL, 0}
};

void R_init_mco(DllInfo *info) {
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
  R_forceSymbols(info, TRUE);
}
