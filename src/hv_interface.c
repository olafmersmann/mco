#define STRICT_R_HEADERS 1

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Applic.h>

#include "hv.h"

char * program_invocation_short_name;

#define UNPACK_REAL_VECTOR(S, D, N)             \
  double *D = REAL(S);                          \
  const R_len_t N = length(S);

#define UNPACK_REAL_MATRIX(S, D, N, K)		\
  double *D = REAL(S);                          \
  const R_len_t N = nrows(S);			\
  const R_len_t K = ncols(S);


SEXP do_hv(SEXP s_data, SEXP s_ref) {
  SEXP s_res;

  if (!isReal(s_data))  error("Argument 's_data' is not a real matrix.");
  if (!isReal(s_ref))   error("Argument 's_ref' is not a real vector.");

  /* Unpack arguments */
  UNPACK_REAL_MATRIX(s_data, data, k_data, n_data);
  UNPACK_REAL_VECTOR(s_ref, ref, n_ref);

  if (n_ref != k_data)
    error("ref and data must be same dimension.");

  /* Allocate result */
  PROTECT(s_res = allocVector(REALSXP, 1));
  double *res = REAL(s_res);

  res[0] = fpli_hv(data, k_data, n_data, ref);
  UNPROTECT(1); /* s_res */
  return s_res;
}
