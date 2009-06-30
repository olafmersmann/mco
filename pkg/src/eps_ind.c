/*===========================================================================*
 * eps_ind.c: implements the unary epsilon indicator as proposed in
 *            Zitzler, E., Thiele, L., Laumanns, M., Fonseca, C., and
 *            Grunert da Fonseca, V (2003): Performance Assessment of
 *            Multiobjective Optimizers: An Analysis and Review. IEEE
 *            Transactions on Evolutionary Computation, 7(2), 117-132.
 *
 * IMPORTANT:
 *   To make the epsilon indicator work for mixed optimization problems
 *   where some objectives are to be maximized while others are to be
 *   minimized, in the case of minimization the value -epsilon (for the
 *   additive case) resp. 1/epsilon (for the multiplicative version) is
 *   considered and returned. Example: suppose f1 is to be minimized and
 *   f2 to be maximized, and the multiplicative epsilon value computed by
 *   this program is 3.0; this means that the considered nondominated front
 *   needs to be multiplied by 1/3 for all f1 values and by 3 for all
 *   f2 values. Thus, independently of which type of problem one considers
 *   (minimization, maximization, mixed minimization/maximization), a lower
 *   indicator value corresponds to a better approximation set.
 *
 * Author:
 *   Eckart Zitzler, February 3, 2005 / last update August 9, 2005
 *
 * Adapted for use in R by:
 *   Olaf Mersmann <olafm@statistik.tu-dortmund.de>
 */

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Applic.h>

typedef enum  {additive, multiplicative} method_t;

/*
 * calc_eps_ind:
 *  a - reference front
 *  b - current front
 */
static double calc_eps_ind(double  *a, size_t size_a, 
						   double  *b, size_t size_b, 
						   size_t dim, method_t method) {
    size_t i, j, k;
    double  eps, eps_j, eps_k, eps_temp;
	
	eps = additive == method ? DBL_MIN : 0.0;

    for (i = 0; i < size_a; i++) {
		for (j = 0; j < size_b; j++) {
			for (k = 0; k < dim; k++) {
				switch (method) {
				case additive:
					// eps_temp = b[j * dim + k] - a[i * dim + k];
					eps_temp = b[j + size_b * k] - a[i + size_a * k];
					break;
				case multiplicative:
					eps_temp = b[j * dim + k] / a[i * dim + k];
				}
				if ((0 == k) || (eps_k < eps_temp))
					eps_k = eps_temp;
			}
			if ((0 == j) || (eps_j > eps_k))
				eps_j = eps_k;
		}
		if ((0 == i) || (eps < eps_j))
			eps = eps_j;
    }    
    return eps;
}


#define UNPACK_REAL_VECTOR(S, D, N)             \
  double *D = REAL(S);                          \
  const R_len_t N = length(S);

#define UNPACK_REAL_MATRIX(S, D, N, K)          \
  double *D = REAL(S);                          \
  const R_len_t N = nrows(S);                   \
  const R_len_t K = ncols(S);


SEXP do_eps_ind(SEXP s_data, SEXP s_ref) {
  SEXP s_res;

  if (!isReal(s_data)) error("Argument 's_data' is not a real matrix.");
  if (!isReal(s_ref))  error("Argument 's_ref' is not a real matrix.");

  /* Unpack arguments */
  UNPACK_REAL_MATRIX(s_data, data, n_data, k_data);
  UNPACK_REAL_MATRIX(s_ref, ref, n_ref, k_ref);
  
  if (k_ref != k_data)
	  error("Reference and current front must have the same dimension.");

  /* Allocate result */
  PROTECT(s_res = allocVector(REALSXP, 1));
  double *res = REAL(s_res);
  /* Calculate criterion */
  res[0] = calc_eps_ind(ref, n_ref, data, n_data, k_data, additive);  
  UNPROTECT(1); /* s_res */
  return s_res;
}
